// Executable assertions for SHIP's encrypted selection and mask algebra.
#include "ship/reference.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <random>
#include <string>

using namespace lbcrypto;
using namespace ship;

static void Require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

static double Check(const CC& cc, const SK& sk, const CT& ct,
                    const std::vector<C>& expected, double tolerance, const std::string& name) {
    Plaintext pt;
    cc->Decrypt(sk, ct, &pt);
    pt->SetLength(expected.size());
    const auto& values = pt->GetCKKSPackedValue();
    double error = 0;
    for (size_t i = 0; i < expected.size(); ++i) {
        Require(std::isfinite(values[i].real()) && std::isfinite(values[i].imag()), name + ": non-finite");
        error = std::max(error, std::abs(values[i] - expected[i]));
    }
    Require(error < tolerance, name + ": error=" + std::to_string(error));
    std::cout << "PASS " << name << " max_error=" << error << " level=" << ct->GetLevel() << '\n';
    return error;
}

static std::vector<C> Rotate(const std::vector<C>& x, int offset) {
    std::vector<C> y(x.size());
    int n = x.size();
    for (int i = 0; i < n; ++i) y[i] = x[(i + offset % n + n) % n];
    return y;
}

static void CheckMasks() {
    size_t cases = 0;
    std::mt19937 rng(20251022);
    for (uint32_t slots : {4u, 8u, 16u, 32u}) {
        const uint32_t n = 2 * slots;
        std::vector<C> a(n);
        for (auto& x : a) x = std::polar(1.0, (rng() % 10000) / 1591.0);
        for (uint32_t j = 0; j < n; ++j) for (int sign : {-1, 1}) {
            auto masks = Masks(j, sign, slots);
            std::vector<C> selected(slots);
            for (uint32_t i = 0; i < slots; ++i)
                selected[i] = masks[0][i]*a[i] + masks[1][i]*std::conj(a[i]) +
                              masks[2][i]*a[i+slots] + masks[3][i]*std::conj(a[i+slots]);
            auto got = Rotate(selected, -static_cast<int>(j % slots));
            // Independent negacyclic monomial multiplication, including both wrap boundaries.
            std::vector<C> want(n);
            for (uint32_t k = 0; k < n; ++k) {
                int s = (k + j >= n ? -1 : 1) * sign;
                want[(k+j)%n] = s > 0 ? a[k] : std::conj(a[k]);
            }
            for (uint32_t i = 0; i < slots; ++i)
                Require(std::abs(got[i]-want[i]) < 1e-12, "negacyclic mask mismatch");
            ++cases;
        }
    }
    std::cout << "PASS mask_oracle cases=" << cases << '\n';
}

int main() {
    try {
        CheckMasks();
        CCParams<CryptoContextCKKSRNS> params;
        params.SetSecurityLevel(HEStd_NotSet);  // Small, explicitly non-production correctness fixture.
        params.SetRingDim(1024);
        params.SetBatchSize(16);
        params.SetMultiplicativeDepth(8);
        params.SetScalingModSize(50);
        params.SetFirstModSize(55);
        params.SetScalingTechnique(FIXEDMANUAL);
        params.SetCKKSDataType(COMPLEX);
        params.SetKeySwitchTechnique(HYBRID);
        auto cc = GenCryptoContext(params);
        cc->Enable(PKE); cc->Enable(KEYSWITCH); cc->Enable(LEVELEDSHE); cc->Enable(ADVANCEDSHE);
        auto kp = cc->KeyGen();
        cc->EvalMultKeyGen(kp.secretKey);
        cc->EvalRotateKeyGen(kp.secretKey, {1, 2, 4, 8, -1, -2, -4, -8, 3, -3});
        std::vector<C> x(16);
        for (size_t i = 0; i < x.size(); ++i) x[i] = C((double(i)-7)/16.0, (int(i%5)-2)/20.0);
        for (uint32_t level : {0u, 3u, 7u}) {
            auto pt = cc->MakeCKKSPackedPlaintext(x, 1, level, nullptr, x.size());
            auto ct = cc->Encrypt(kp.publicKey, pt);
            auto originalElements = ct->GetElements();
            for (uint32_t beta : {0u, 1u}) for (int rot : {0, 3, -3}) {
                auto key = MakeMuxKey(cc, kp.secretKey, beta, rot);
                auto out = MuxRotate(cc, key, ct);
                auto expected = Rotate(x, rot);
                for (auto& v : expected) v *= beta;
                Require(out->GetLevel() == ct->GetLevel() &&
                        out->GetNoiseScaleDeg() == ct->GetNoiseScaleDeg() &&
                        out->GetScalingFactor() == ct->GetScalingFactor() &&
                        out->GetElements()[0].GetNumOfElements() == ct->GetElements()[0].GetNumOfElements(),
                        "HMuxRot consumed a level or changed scale metadata");
                Check(cc, kp.secretKey, out, expected, 1e-7,
                      "hmux beta="+std::to_string(beta)+" rot="+std::to_string(rot));
            }
            for (uint32_t shift : {0u, 3u, 8u, 15u}) {
                auto keys = MakeBlindKey(cc, kp.secretKey, shift, x.size());
                auto out = BlindRotate(cc, keys, ct);
                Require(out->GetLevel() == ct->GetLevel(), "blind rotation consumed levels");
                Check(cc, kp.secretKey, out, Rotate(x, shift), 1e-7,
                      "blind shift="+std::to_string(shift));
            }
            Require(ct->GetElements() == originalElements, "input ciphertext was modified");
        }
        std::cout << "ALL COMPONENT CHECKS PASSED (toy parameters; not a security/performance claim)\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL " << e.what() << '\n';
        return 1;
    }
}
