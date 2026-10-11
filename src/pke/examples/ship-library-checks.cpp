// Functional checks of the library SHIP bootstrapping (toy, non-secure parameters).
// Secure parameter sets are exercised by ship-paper-bench.
#include "openfhe.h"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <random>

using namespace lbcrypto;
using C     = std::complex<double>;
using Clock = std::chrono::steady_clock;

static double MaxError(const CryptoContext<DCRTPoly>& cc, const PrivateKey<DCRTPoly>& sk, const Ciphertext<DCRTPoly>& ct,
                       const std::vector<C>& expected) {
    Plaintext pt;
    cc->Decrypt(sk, ct, &pt);
    pt->SetLength(expected.size());
    double error = 0;
    for (size_t i = 0; i < expected.size(); ++i) {
        const auto v = pt->GetCKKSPackedValue()[i];
        if (!std::isfinite(v.real()) || !std::isfinite(v.imag()))
            throw std::runtime_error("non-finite output");
        error = std::max(error, std::abs(v - expected[i]));
    }
    return error;
}

static uint32_t Limbs(const Ciphertext<DCRTPoly>& ct) {
    return ct->GetElements()[0].GetNumOfElements();
}

static SHIPContextSpec ToySpec(uint32_t h) {
    SHIPContextSpec s;
    s.ringDim        = 1024;
    s.firstModBits   = 60;
    s.scalingModBits = 50;  // gamma = 2^10
    s.multLevels     = 1;
    s.hammingWeight  = h;
    s.numLargeDigits = 3;
    s.auxModBits     = 60;
    s.securityLevel  = HEStd_NotSet;
    return s;
}

template <class F>
static void RequireThrow(F&& f, const char* what) {
    try {
        f();
    }
    catch (const std::exception&) {
        return;
    }
    throw std::runtime_error(std::string("accepted invalid input: ") + what);
}

static void Run(const char* name, const SHIPParams& sp, bool viaContext) {
    auto cc          = GenSHIPCryptoContext(ToySpec(sp.hammingWeight));
    const uint32_t S = cc->GetRingDimension() / 2;
    const uint32_t Q = cc->GetElementParams()->GetParams().size();
    auto kp          = cc->KeyGen();
    cc->EvalMultKeyGen(kp.secretKey);
    auto start = Clock::now();
    std::shared_ptr<SHIPBootstrapKey> key;
    if (viaContext)
        cc->EvalSHIPBootstrapKeyGen(kp.secretKey, sp);
    else
        key = SHIPKeyGen(kp.secretKey, sp);
    const double keygenMs = std::chrono::duration<double, std::milli>(Clock::now() - start).count();
    auto boot = [&](ConstCiphertext<DCRTPoly>& ct) { return viaContext ? cc->EvalSHIPBootstrap(ct) : SHIPBootstrap(ct, *key); };

    std::mt19937 rng(20261011);
    std::uniform_real_distribution<double> u(-1, 1);
    std::vector<C> x(S);
    for (auto& v : x)
        v = sp.realOnly ? C(u(rng), 0) : C(u(rng), u(rng));
    // Input at the S2C level (two limbs) and, separately, at a higher level (extra limbs are dropped).
    for (uint32_t inputLimbs : {2u, 4u}) {
        auto ct = cc->Encrypt(kp.publicKey, cc->MakeCKKSPackedPlaintext(x, 1, Q - inputLimbs, nullptr, S));
        start   = Clock::now();
        auto out = boot(ct);
        const double ms   = std::chrono::duration<double, std::milli>(Clock::now() - start).count();
        const double err  = MaxError(cc, kp.secretKey, out, x);
        std::cout << name << " input_limbs=" << inputLimbs << " max_error=" << err
                  << " precision_bits=" << -std::log2(err) << " output_limbs=" << Limbs(out)
                  << " ms=" << ms << " keygen_ms=" << keygenMs << std::endl;
        if (err >= 1e-3 || Limbs(out) != Q - SHIPProductTreeDepth(sp.hammingWeight))
            throw std::runtime_error("SHIP bootstrap failure");
        if (inputLimbs == 2) {
            auto sq = cc->EvalMult(out, out);
            cc->RescaleInPlace(sq);
            auto expected = x;
            for (auto& v : expected)
                v *= v;
            auto again        = boot(sq);
            const double err2 = MaxError(cc, kp.secretKey, again, expected);
            std::cout << name << " square_then_bootstrap max_error=" << err2 << std::endl;
            if (err2 >= 2e-3)
                throw std::runtime_error("second bootstrap failure");
            auto bad = ct->Clone();
            bad->SetKeyTag("other");
            RequireThrow([&] { boot(bad); }, "key tag");
            auto low = ct->Clone();
            cc->LevelReduceInPlace(low, nullptr, 1);
            RequireThrow([&] { boot(low); }, "one limb");
            auto unrescaled = cc->EvalMult(ct, ct);
            RequireThrow([&] { boot(unrescaled); }, "noise scale degree 2");
        }
    }
    if (!viaContext)
        std::cout << name << " key_payload_bytes=" << SHIPKeyStoredBytes(*key) << std::endl;
}

int main() {
    try {
        SHIPParams windowed;  // paper structure: windowed sparse key, column + base-4 mux
        windowed.hammingWeight        = 8;
        windowed.window               = 20;
        windowed.columnSize           = 6;
        windowed.muxBase              = 4;
        windowed.encapsulationModBits = 60;
        windowed.realOnly             = true;
        Run("real_windowed_theta6_B4", windowed, false);

        SHIPParams complexMux = windowed;  // Algorithm 3 only: unrestricted key, binary mux
        complexMux.window     = 0;
        complexMux.columnSize = 1;
        complexMux.muxBase    = 2;
        complexMux.realOnly   = false;
        Run("complex_unrestricted_binary_mux", complexMux, false);

        SHIPParams api = windowed;
        api.realOnly   = false;
        Run("complex_windowed_via_CryptoContext", api, true);

        // Missing relinearization key must be reported.
        auto cc = GenSHIPCryptoContext(ToySpec(8));
        auto kp = cc->KeyGen();
        RequireThrow([&] { SHIPKeyGen(kp.secretKey, windowed); }, "missing EvalMultKeyGen");
        std::cout << "ALL SHIP LIBRARY CHECKS PASSED (toy parameters, HEStd_NotSet)" << std::endl;
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "FAIL " << e.what() << std::endl;
        return 1;
    }
}
