#include "ship/full-bootstrap.h"
#include <chrono>
#include <iostream>
#include <random>

using namespace lbcrypto;
using namespace ship;

static double Error(const CC& cc, const SK& secret, const CT& ciphertext, const std::vector<C>& expected) {
    Plaintext decoded;
    cc->Decrypt(secret,ciphertext,&decoded);
    decoded->SetLength(expected.size());
    double error = 0;
    for (size_t i = 0; i < expected.size(); ++i) {
        const auto value = decoded->GetCKKSPackedValue()[i];
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag())) throw std::runtime_error("non-finite output");
        error = std::max(error,std::abs(value-expected[i]));
    }
    return error;
}

// Test-only oracle: decrypt the bottom polynomial, apply the exact sine target,
// and evaluate the canonical embedding directly (no production packing code).
static std::vector<C> SineOracle(const SK& sparse, const CT& input) {
    auto secret = sparse->GetPrivateElement().GetElementAtIndex(0);
    secret.SetFormat(Format::EVALUATION);
    auto phase = input->GetElements()[0].GetElementAtIndex(0) +
                 input->GetElements()[1].GetElementAtIndex(0) * secret;
    phase.SetFormat(Format::COEFFICIENT);
    const uint32_t slots = phase.GetLength()/2, cycl = 4*slots;
    const long double q = phase.GetModulus().ConvertToInt();
    const long double tau = 2*std::acos(-1.0L), gamma = q/input->GetScalingFactor();
    auto coefficient = [&](uint32_t j) {
        long double centered = phase[j].ConvertToInt();
        if (centered > q/2) centered -= q;
        return static_cast<double>((gamma/tau)*std::sin(tau*centered/q));
    };
    std::vector<C> result(slots);
    uint64_t root = 1;
    for (uint32_t k = 0; k < slots; ++k) {
        for (uint32_t j = 0; j < slots; ++j) {
            double angle = static_cast<double>(tau*((root*j)%cycl)/cycl);
            result[k] += C(coefficient(j),coefficient(j+slots))*std::polar(1.0,angle);
        }
        root = (5*root)%cycl;
    }
    return result;
}

template<class Function> static void RequireRejected(Function&& operation) {
    try { operation(); }
    catch (const std::invalid_argument&) { return; }
    throw std::runtime_error("unsupported input was accepted");
}

int main() {
    try {
        for (const auto& [ringDim,hamming] : std::vector<std::pair<uint32_t,uint32_t>>{{128,4},{1024,8}}) {
            const uint32_t slots = ringDim/2, depth = 9;
            CCParams<CryptoContextCKKSRNS> p;
            p.SetSecurityLevel(HEStd_NotSet); // Functional tests only, not a security parameter set.
            p.SetRingDim(ringDim); p.SetBatchSize(slots); p.SetMultiplicativeDepth(depth);
            p.SetScalingModSize(50); p.SetFirstModSize(60);
            p.SetScalingTechnique(FIXEDMANUAL); p.SetKeySwitchTechnique(HYBRID);
            p.SetCKKSDataType(COMPLEX);
            auto cc = GenCryptoContext(p);
            cc->Enable(PKE); cc->Enable(KEYSWITCH); cc->Enable(LEVELEDSHE); cc->Enable(ADVANCEDSHE);
            auto dense = cc->KeyGen();
            cc->EvalMultKeyGen(dense.secretKey);
            DCRTPoly::TugType ternary;
            auto sparse = std::make_shared<PrivateKeyImpl<DCRTPoly>>(cc);
            sparse->SetPrivateElement(DCRTPoly(ternary,dense.secretKey->GetPrivateElement().GetParams(),
                                               Format::EVALUATION,hamming));
            auto keys = MakeFullBootstrapKey(cc,dense,sparse);
            std::mt19937 rng(20261009);
            for (uint32_t trial = 0; trial < 5; ++trial) {
                std::vector<C> x(slots);
                for (uint32_t i = 0; i < slots; ++i) {
                    if (trial == 1) x[i] = C(0.25,-0.125);
                    if (trial == 2 && i == slots-1) x[i] = C(-0.5,0.375);
                    if (trial == 3) x[i] = C((int(i%9)-4)/8.0,(int(i%7)-3)/8.0);
                    if (trial == 4) x[i] = C((int(rng()%1001)-500)/1000.0,(int(rng()%1001)-500)/1000.0);
                }
                // Ordinary OpenFHE packed encryption at the last level, with no coefficient edits.
                auto input = cc->Encrypt(dense.publicKey,cc->MakeCKKSPackedPlaintext(x,1,depth,nullptr,slots));
                const auto before = input->GetElements();
                auto encapsulated = ApplyBottomSwitch(input,keys.encapsulation);
                double switchError = Error(cc,sparse,encapsulated,x);
                if (switchError >= 1e-7) throw std::runtime_error("bottom key-switch error");
                auto start = std::chrono::steady_clock::now();
                auto result = FullBootstrap(cc,input,keys);
                double ms = std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
                double error = Error(cc,dense.secretKey,result,x);
                const double phaseError = Error(cc,dense.secretKey,result,SineOracle(sparse,encapsulated));
                std::cout << "N=" << ringDim << " h=" << hamming << " trial=" << trial
                          << " switch_error=" << switchError << " slot_error=" << error
                          << " sine_oracle_error=" << phaseError
                          << " output_level=" << result->GetLevel()
                          << " input_towers=1 output_towers=" << result->GetElements()[0].GetNumOfElements()
                          << " elapsed_ms=" << ms << std::endl;
                if (error >= 1e-4) throw std::runtime_error("slot round-trip error");
                if (phaseError >= 5e-6) throw std::runtime_error("independent sine/packing oracle mismatch");
                if (result->GetElements()[0].GetNumOfElements() < 3 || result->GetKeyTag() != dense.secretKey->GetKeyTag())
                    throw std::runtime_error("output budget or output key mismatch");
                if (before != input->GetElements()) throw std::runtime_error("input was mutated");
                if (trial == 0) {
                    auto invalid = input->Clone();
                    invalid->SetKeyTag("different-key");
                    RequireRejected([&] { FullBootstrap(cc,invalid,keys); });
                    invalid = input->Clone(); invalid->SetNoiseScaleDeg(2);
                    RequireRejected([&] { FullBootstrap(cc,invalid,keys); });
                    invalid = input->Clone(); invalid->SetSlots(slots/2);
                    RequireRejected([&] { FullBootstrap(cc,invalid,keys); });
                    invalid = input->Clone(); invalid->SetScalingFactor(0);
                    RequireRejected([&] { FullBootstrap(cc,invalid,keys); });
                    auto tooHigh = cc->Encrypt(dense.publicKey,
                        cc->MakeCKKSPackedPlaintext(x,1,depth-1,nullptr,slots));
                    RequireRejected([&] { FullBootstrap(cc,tooHigh,keys); });
                }
                if (trial == 4) {
                    auto squared = cc->EvalMult(result,result);
                    cc->RescaleInPlace(squared);
                    auto expected = x;
                    for (auto& v : expected) v *= v;
                    auto squareError = Error(cc,dense.secretKey,squared,expected);
                    if (squareError >= 2e-4) throw std::runtime_error("restored budget is not usable");
                    cc->LevelReduceInPlace(squared,nullptr,squared->GetElements()[0].GetNumOfElements()-1);
                    auto refreshedAgain = FullBootstrap(cc,squared,keys);
                    double secondError = Error(cc,dense.secretKey,refreshedAgain,expected);
                    std::cout << "N=" << ringDim << " square_error=" << squareError
                              << " second_refresh_error=" << secondError << std::endl;
                    if (secondError >= 3e-4) throw std::runtime_error("second refresh failed");
                }
            }
            // Deliberately leave the small-angle regime. A successful sine
            // evaluation must NOT be misreported as accurate message recovery.
            std::vector<C> large(slots,C(64,-32));
            auto largeInput = cc->Encrypt(dense.publicKey,
                cc->MakeCKKSPackedPlaintext(large,1,depth,nullptr,slots));
            auto largeSparse = ApplyBottomSwitch(largeInput,keys.encapsulation);
            auto largeResult = FullBootstrap(cc,largeInput,keys);
            const auto largeError = Error(cc,dense.secretKey,largeResult,large);
            const auto largePhaseError = Error(cc,dense.secretKey,largeResult,SineOracle(sparse,largeSparse));
            std::cout << "N=" << ringDim << " large_message_error=" << largeError
                      << " large_sine_oracle_error=" << largePhaseError << std::endl;
            if (largeError < 1 || largePhaseError >= 5e-6)
                throw std::runtime_error("sine approximation limit was not reproduced");
        }
        std::cout << "ALL FULL BOOTSTRAP REFERENCE CHECKS PASSED (toy parameters; not optimized or security validated)\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL " << e.what() << '\n';
        return 1;
    }
}
