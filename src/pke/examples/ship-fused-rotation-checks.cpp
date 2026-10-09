#include "ship/fused-rotation.h"
#include <iostream>
#include <set>

using namespace ship;

static double Error(const CC& cc, const SK& sk, const CT& ct, const std::vector<C>& expected) {
    Plaintext decoded;
    cc->Decrypt(sk,ct,&decoded);
    decoded->SetLength(expected.size());
    double worst = 0;
    for (size_t i = 0; i < expected.size(); ++i) {
        const double error = std::abs(decoded->GetCKKSPackedValue()[i]-expected[i]);
        if (!std::isfinite(error)) throw std::runtime_error("non-finite rotation output");
        worst = std::max(worst,error);
    }
    return worst;
}

static std::vector<C> Oracle(const std::vector<C>& x, int shift, bool beta = true) {
    std::vector<C> result(x.size());
    for (int i = 0; i < static_cast<int>(x.size()); ++i)
        if (beta) result[i] = x[(i+shift+2*static_cast<int>(x.size()))%x.size()];
    return result;
}

static void CheckMetadata(const CT& input, const CT& output) {
    if (input->GetLevel() != output->GetLevel() || input->GetScalingFactor() != output->GetScalingFactor() ||
        input->GetNoiseScaleDeg() != output->GetNoiseScaleDeg() || input->GetSlots() != output->GetSlots() ||
        input->GetKeyTag() != output->GetKeyTag() ||
        input->GetElements()[0].GetNumOfElements() != output->GetElements()[0].GetNumOfElements())
        throw std::runtime_error("fused rotation changed ciphertext metadata/budget");
}

int main() {
    try {
        const uint32_t slots = 512;
        CCParams<CryptoContextCKKSRNS> p;
        p.SetSecurityLevel(HEStd_NotSet); p.SetRingDim(2*slots); p.SetBatchSize(slots);
        p.SetMultiplicativeDepth(9); p.SetScalingModSize(50); p.SetFirstModSize(60);
        p.SetScalingTechnique(FIXEDMANUAL); p.SetKeySwitchTechnique(HYBRID); p.SetCKKSDataType(COMPLEX);
        auto cc = GenCryptoContext(p);
        cc->Enable(PKE); cc->Enable(KEYSWITCH); cc->Enable(LEVELEDSHE); cc->Enable(ADVANCEDSHE);
        auto keys = cc->KeyGen();
        std::vector<C> x(slots);
        for (uint32_t i = 0; i < slots; ++i) x[i] = C((int(i%17)-8)/16.0,(int(i%19)-9)/16.0);
        auto first = cc->Encrypt(keys.publicKey,cc->MakeCKKSPackedPlaintext(x,1,0,nullptr,slots));
        const auto probeKey = MakeFusedMuxKey(cc,keys.secretKey,1,3);
        if (Error(cc,keys.secretKey,FusedMuxRotate(cc,probeKey,first),Oracle(x,3)) >= 1e-7)
            throw std::runtime_error("fused rotation needs an external rotation key");
        std::cout << "PASS fused rotation without EvalRotateKeyGen\n";

        std::set<int32_t> indices{3,-3};
        for (uint32_t bit = 1; bit < slots; bit <<= 1) { indices.insert(bit); indices.insert(-int32_t(bit)); }
        cc->EvalRotateKeyGen(keys.secretKey,std::vector<int32_t>(indices.begin(),indices.end()));
        uint32_t cases = 0;
        double worstError = 0, worstDifference = 0;
        auto compare = [&](const CT& input, const CT& fused, const CT& reference, const std::vector<C>& expected) {
            CheckMetadata(input,fused);
            const double error = Error(cc,keys.secretKey,fused,expected);
            const double referenceError = Error(cc,keys.secretKey,reference,expected);
            const double difference = Error(cc,keys.secretKey,cc->EvalSub(fused,reference),std::vector<C>(slots));
            if (error >= 1e-7 || referenceError >= 1e-7 || difference >= 1e-7)
                throw std::runtime_error("fused/reference/plaintext rotation mismatch");
            worstError = std::max(worstError,error); worstDifference = std::max(worstDifference,difference);
            ++cases;
        };
        for (uint32_t beta : {0u,1u}) for (int32_t rotation : {0,3,-3}) {
            auto fused = MakeFusedMuxKey(cc,keys.secretKey,beta,rotation);
            auto reference = MakeMuxKey(cc,keys.secretKey,beta,rotation);
            for (uint32_t level : {0u,3u,7u}) {
                auto input = cc->Encrypt(keys.publicKey,cc->MakeCKKSPackedPlaintext(x,1,level,nullptr,slots));
                const auto before = input->GetElements();
                compare(input,FusedMuxRotate(cc,fused,input),MuxRotate(cc,reference,input),Oracle(x,rotation,beta));
                if (before != input->GetElements()) throw std::runtime_error("conditional rotation mutated input");
            }
        }
        for (int direction : {-1,1}) for (uint32_t shift : {0u,3u,128u,511u}) {
            auto fused = MakeFusedBlindKey(cc,keys.secretKey,shift,slots,direction);
            auto reference = MakeBlindKey(cc,keys.secretKey,shift,slots,direction);
            for (uint32_t level : {0u,3u,7u}) {
                auto input = cc->Encrypt(keys.publicKey,cc->MakeCKKSPackedPlaintext(x,1,level,nullptr,slots));
                const auto before = input->GetElements();
                compare(input,FusedBlindRotate(cc,fused,input),BlindRotate(cc,reference,input),Oracle(x,direction*int(shift)));
                if (before != input->GetElements()) throw std::runtime_error("blind rotation mutated input");
            }
        }
        std::cout << "PASS fused/reference/oracle cases=" << cases << " max_error=" << worstError
                  << " max_difference=" << worstDifference << '\n';
        auto wrong = first->Clone(); wrong->SetKeyTag("wrong-key");
        bool rejected = false;
        try { FusedMuxRotate(cc,probeKey,wrong); } catch (const std::invalid_argument&) { rejected = true; }
        if (!rejected) throw std::runtime_error("wrong key accepted");
        std::cout << "ALL FUSED ROTATION CHECKS PASSED (toy parameters; security unverified)\n";
        return 0;
    } catch (const std::exception& e) { std::cerr << "FAIL " << e.what() << '\n'; return 1; }
}
