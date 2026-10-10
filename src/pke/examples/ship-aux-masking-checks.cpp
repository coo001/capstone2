// Auxiliary-modulus masking (paper Sec. 4.1/4.4) versus the rescaling path of Alg. 1.
// Toy parameters only: no security claim, timings are single-machine references.
#include "ship/full-bootstrap.h"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <random>

using namespace lbcrypto;
using namespace ship;
using Clock = std::chrono::steady_clock;

static double Error(const CC& cc, const SK& secret, const CT& ciphertext, const std::vector<C>& expected) {
    Plaintext decoded;
    cc->Decrypt(secret, ciphertext, &decoded);
    decoded->SetLength(expected.size());
    double error = 0;
    for (size_t i = 0; i < expected.size(); ++i) {
        const auto value = decoded->GetCKKSPackedValue()[i];
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag())) throw std::runtime_error("non-finite output");
        error = std::max(error, std::abs(value-expected[i]));
    }
    return error;
}

struct Setup {
    CC cc;
    KeyPair<DCRTPoly> dense;
    SK sparse;
};

static Setup MakeSetup(uint32_t ringDim, uint32_t depth, uint32_t hamming) {
    CCParams<CryptoContextCKKSRNS> p;
    p.SetSecurityLevel(HEStd_NotSet);  // Functional checks only.
    p.SetRingDim(ringDim); p.SetBatchSize(ringDim/2); p.SetMultiplicativeDepth(depth);
    p.SetScalingModSize(50); p.SetFirstModSize(60);
    p.SetScalingTechnique(FIXEDMANUAL); p.SetKeySwitchTechnique(HYBRID); p.SetCKKSDataType(COMPLEX);
    Setup s;
    s.cc = GenCryptoContext(p);
    s.cc->Enable(PKE); s.cc->Enable(KEYSWITCH); s.cc->Enable(LEVELEDSHE); s.cc->Enable(ADVANCEDSHE);
    s.dense = s.cc->KeyGen();
    s.cc->EvalMultKeyGen(s.dense.secretKey);
    DCRTPoly::TugType ternary;
    s.sparse = std::make_shared<PrivateKeyImpl<DCRTPoly>>(s.cc);
    s.sparse->SetPrivateElement(DCRTPoly(ternary, s.dense.secretKey->GetPrivateElement().GetParams(),
                                         Format::EVALUATION, hamming));
    return s;
}

static uint32_t Limbs(const CT& ct) { return ct->GetElements()[0].GetNumOfElements(); }

static bool Rejected(const Setup& s, bool auxMasking) {
    try { MakeFullBootstrapKey(s.cc, s.dense, s.sparse, true, auxMasking); }
    catch (const std::invalid_argument&) { return true; }
    return false;
}

int main() {
    try {
        const uint32_t ringDim = 1024, slots = ringDim/2, h = 8, depth = 9;
        auto s = MakeSetup(ringDim, depth, h);
        const auto P = CkksParams(s.cc)->GetParamsP()->GetModulus();
        std::cout << "N=" << ringDim << " h=" << h << " depth=" << depth << " log2P=" << P.GetMSB() << std::endl;
        auto aux = MakeFullBootstrapKey(s.cc, s.dense, s.sparse, true, true);
        auto rescale = MakeFullBootstrapKey(s.cc, s.dense, s.sparse, true, false);
        std::mt19937 rng(20261010);
        std::vector<double> auxMs, rescaleMs;
        for (uint32_t trial = 0; trial < 4; ++trial) {
            std::vector<C> x(slots);
            for (uint32_t i = 0; i < slots; ++i) {
                if (trial == 1) x[i] = C(0.25, -0.125);
                if (trial == 2) x[i] = C((int(i%9)-4)/8.0, (int(i%7)-3)/8.0);
                if (trial == 3) x[i] = C((int(rng()%1001)-500)/1000.0, (int(rng()%1001)-500)/1000.0);
            }
            auto input = s.cc->Encrypt(s.dense.publicKey, s.cc->MakeCKKSPackedPlaintext(x, 1, depth, nullptr, slots));
            const auto before = input->GetElements();
            auto start = Clock::now();
            auto a = FullBootstrap(s.cc, input, aux);
            auxMs.push_back(std::chrono::duration<double,std::milli>(Clock::now()-start).count());
            start = Clock::now();
            auto r = FullBootstrap(s.cc, input, rescale);
            rescaleMs.push_back(std::chrono::duration<double,std::milli>(Clock::now()-start).count());
            if (before != input->GetElements()) throw std::runtime_error("input mutated");
            const double auxError = Error(s.cc, s.dense.secretKey, a, x);
            const double rescaleError = Error(s.cc, s.dense.secretKey, r, x);
            if (a->GetLevel()+1 != r->GetLevel() || Limbs(a) != Limbs(r)+1)
                throw std::runtime_error("auxiliary masking did not save exactly one level");
            if (a->GetScalingFactor() != r->GetScalingFactor() || a->GetNoiseScaleDeg() != 1 ||
                a->GetKeyTag() != s.dense.secretKey->GetKeyTag())
                throw std::runtime_error("auxiliary masking changed output metadata");
            auto lowered = a->Clone();
            s.cc->LevelReduceInPlace(lowered, nullptr, 1);
            const double difference = Error(s.cc, s.dense.secretKey, s.cc->EvalSub(lowered, r), std::vector<C>(slots));
            std::cout << "trial=" << trial << " aux_error=" << auxError << " rescale_error=" << rescaleError
                      << " path_difference=" << difference << " aux_towers=" << Limbs(a)
                      << " rescale_towers=" << Limbs(r) << " aux_ms=" << auxMs.back()
                      << " rescale_ms=" << rescaleMs.back() << std::endl;
            if (auxError >= 1e-4 || difference >= 5e-6) throw std::runtime_error("auxiliary masking accuracy failure");
            if (trial == 3) {
                // The extra limb must be usable: square, rescale, and bootstrap again.
                auto squared = s.cc->EvalMult(a, a);
                s.cc->RescaleInPlace(squared);
                auto expected = x;
                for (auto& v : expected) v *= v;
                const double squareError = Error(s.cc, s.dense.secretKey, squared, expected);
                s.cc->LevelReduceInPlace(squared, nullptr, Limbs(squared)-1);
                const double refreshError = Error(s.cc, s.dense.secretKey, FullBootstrap(s.cc, squared, aux), expected);
                std::cout << "square_error=" << squareError << " second_refresh_error=" << refreshError << std::endl;
                if (squareError >= 2e-4 || refreshError >= 3e-4) throw std::runtime_error("refresh after aux masking failed");
            }
        }
        std::sort(auxMs.begin(), auxMs.end());
        std::sort(rescaleMs.begin(), rescaleMs.end());
        std::cout << "median_ms aux=" << auxMs[auxMs.size()/2] << " rescale=" << rescaleMs[rescaleMs.size()/2]
                  << " (4 single runs; not a benchmark)" << std::endl;

        // Minimum depth: only the auxiliary path fits one level lower.
        // h=8 -> tree 4; consumed = masking + 4 + 1 (FFT); one spare limb is required.
        {
            auto tight = MakeSetup(ringDim, 6, h);
            const bool auxRejected = Rejected(tight, true), rescaleRejected = Rejected(tight, false);
            std::cout << "depth=6 h=8 aux_rejected=" << auxRejected << " rescale_rejected=" << rescaleRejected << std::endl;
            if (auxRejected || !rescaleRejected) throw std::runtime_error("unexpected budget boundary at h=8");
            auto keys = MakeFullBootstrapKey(tight.cc, tight.dense, tight.sparse, true, true);
            std::vector<C> x(slots, C(0.3, -0.2));
            auto input = tight.cc->Encrypt(tight.dense.publicKey, tight.cc->MakeCKKSPackedPlaintext(x, 1, 6, nullptr, slots));
            auto out = FullBootstrap(tight.cc, input, keys);
            const double error = Error(tight.cc, tight.dense.secretKey, out, x);
            std::cout << "depth=6 h=8 aux_error=" << error << " output_towers=" << Limbs(out) << std::endl;
            if (error >= 1e-4 || Limbs(out) != 2) throw std::runtime_error("tight-depth aux bootstrap failed");
        }
        // Paper's Hamming weight h=31: tree ceil(log2(32))=5, so 5 tree + 1 FFT = 6 levels as in the paper.
        {
            const uint32_t h31 = 31, tightDepth = 7;
            auto tight = MakeSetup(ringDim, tightDepth, h31);
            const bool auxRejected = Rejected(tight, true), rescaleRejected = Rejected(tight, false);
            std::cout << "depth=7 h=31 aux_rejected=" << auxRejected << " rescale_rejected=" << rescaleRejected << std::endl;
            if (auxRejected || !rescaleRejected) throw std::runtime_error("unexpected budget boundary at h=31");
            auto keys = MakeFullBootstrapKey(tight.cc, tight.dense, tight.sparse, true, true);
            std::vector<C> x(slots);
            for (uint32_t i = 0; i < slots; ++i) x[i] = C((int(i%11)-5)/10.0, (int(i%5)-2)/5.0);
            auto input = tight.cc->Encrypt(tight.dense.publicKey, tight.cc->MakeCKKSPackedPlaintext(x, 1, tightDepth, nullptr, slots));
            auto start = Clock::now();
            auto out = FullBootstrap(tight.cc, input, keys);
            const double ms = std::chrono::duration<double,std::milli>(Clock::now()-start).count();
            const double error = Error(tight.cc, tight.dense.secretKey, out, x);
            const uint32_t consumed = (tightDepth+1)-Limbs(out);
            std::cout << "depth=7 h=31 aux_error=" << error << " levels_consumed=" << consumed
                      << " output_towers=" << Limbs(out) << " ms=" << ms << std::endl;
            if (error >= 1e-4 || consumed != 6) throw std::runtime_error("h=31 aux bootstrap did not use 6 levels");
        }
        std::cout << "ALL AUX-MASKING CHECKS PASSED (toy parameters; security unverified)\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL " << e.what() << '\n';
        return 1;
    }
}
