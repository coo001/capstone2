// SHIP with 128-bit HE-standard parameter sets corresponding to the paper's LL13 / LL14 (Table 2-3).
// usage: ship-paper-bench [LL13|LL14] [trials]
#include "openfhe.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#ifdef _OPENMP
    #include <omp.h>
#endif

using namespace lbcrypto;
using Clock = std::chrono::steady_clock;

static double Seconds(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

int main(int argc, char** argv) {
    try {
        const std::string set = argc > 1 ? argv[1] : "LL13";
        const uint32_t trials = argc > 2 ? std::stoul(argv[2]) : 5;
        SHIPContextSpec spec;
        SHIPParams sp;
        uint32_t sparseBound = 0;  // log2(q0 p') with 128-bit security for h = 31 (lattice estimator, research/ship/security)
        if (set == "LL13" || set == "LL14" || set == "HT14" || set == "HT15") {
            spec        = set == "LL13" ? SHIPContextSpec::LL13() :
                          set == "LL14" ? SHIPContextSpec::LL14() :
                          set == "HT14" ? SHIPContextSpec::HT14() :
                                          SHIPContextSpec::HT15();
            sp          = SHIPParams::Recommended(spec);
            sparseBound = spec.ringDim == (1 << 13) ? 42 : spec.ringDim == (1 << 14) ? 88 : 105;
        }
        else if (set == "custom") {
            // custom trials N q0 scale aux dnum h w theta pbits [boot]  (HEStd_NotSet: diagnostics only)
            if (argc < 12)
                throw std::invalid_argument("custom trials N q0 scale aux dnum h w theta pbits");
            spec.ringDim            = std::stoul(argv[3]);
            spec.firstModBits       = std::stoul(argv[4]);
            spec.scalingModBits     = std::stoul(argv[5]);
            spec.auxModBits         = std::stoul(argv[6]);
            spec.numLargeDigits     = std::stoul(argv[7]);
            spec.hammingWeight      = sp.hammingWeight = std::stoul(argv[8]);
            sp.window               = std::stoul(argv[9]);
            sp.columnSize           = std::stoul(argv[10]);
            sp.encapsulationModBits = std::stoul(argv[11]);
            spec.multLevels         = 1;
            spec.bootModBits        = argc > 12 ? std::stoul(argv[12]) : 0;
            spec.securityLevel      = HEStd_NotSet;
            sparseBound             = 1000;
        }
        else {
            throw std::invalid_argument("usage: ship-paper-bench [LL13|LL14|HT14|HT15|custom] [trials] ...");
        }
        // SHIP_FACTOR_DIR: keep the per-coefficient keys on disk (needed for HT14/HT15 on 48 GB machines).
        const std::string factorDir = std::getenv("SHIP_FACTOR_DIR") ? std::getenv("SHIP_FACTOR_DIR") : "";
        auto start = Clock::now();
        auto cc    = GenSHIPCryptoContext(spec);  // throws unless log2(QP) meets HEStd_128_classic
        const double contextSeconds = Seconds(start);
        const auto params = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc->GetCryptoParameters());
        const uint32_t N  = cc->GetRingDimension();
        const uint32_t S  = N / 2;
        std::cout << "set=" << set << " N=" << N << " slots=" << S << " log2QP=" << SHIPLogQP(cc)
                  << " HEStd_128_classic_max_log2QP=" << StdLatticeParm::FindMaxQ(HEStd_ternary, HEStd_128_classic, N)
                  << " dnum=" << params->GetNumPartQ() << std::endl;
        std::cout << "Q_bits=";
        for (const auto& q : params->GetElementParams()->GetParams())
            std::cout << q->GetModulus().GetMSB() << ' ';
        std::cout << " P_bits=";
        for (const auto& p : params->GetParamsP()->GetParams())
            std::cout << p->GetModulus().GetMSB() << ' ';
        std::cout << std::endl;
        const uint32_t q0Bits = params->GetElementParams()->GetParams()[0]->GetModulus().GetMSB();
        // SHIP_ENCAPS_BITS: size of p' (for a smaller log2(q0 p') than the published bound).
        if (const char* env = std::getenv("SHIP_ENCAPS_BITS"))
            sp.encapsulationModBits = std::stoul(env);
        if (sp.encapsulationModBits == 0)
            sp.encapsulationModBits = std::min<uint32_t>(60, sparseBound - q0Bits);  // what SHIPKeyGen chooses
        const double windowLog = 0.4 * sp.hammingWeight * std::log2(4.0 * sp.window);
        std::cout << "h=" << sp.hammingWeight << " w=" << sp.window << " theta=" << sp.columnSize
                  << " B=" << sp.muxBase << " log2(q0*p')<=" << q0Bits + sp.encapsulationModBits
                  << " (estimator bound " << sparseBound << ")"
                  << " May21_window_cost_log2=" << windowLog << " (>=115 required)" << std::endl;
        if (set != "custom" && (q0Bits + sp.encapsulationModBits > sparseBound || windowLog < 115))
            throw std::runtime_error("sparse-secret parameters exceed the paper's 128-bit bounds");

        auto kp = cc->KeyGen();
        cc->EvalMultKeyGen(kp.secretKey);
        std::cout << "estimated_ship_key_GiB=" << SHIPEstimateKeyBytes(cc, sp) / 1073741824.0
                  << " factor_dir=" << (factorDir.empty() ? "(memory)" : factorDir) << std::endl;
        if (trials == 0)
            return 0;
        start = Clock::now();
        cc->EvalSHIPBootstrapKeyGen(kp.secretKey, sp, factorDir);
        const double keygenSeconds = Seconds(start);
        const auto key            = SHIPGetBootstrapKey(kp.secretKey->GetKeyTag());
        std::cout << "context_s=" << contextSeconds << " keygen_s=" << keygenSeconds
                  << " ship_key_in_memory_GiB=" << SHIPKeyStoredBytes(*key) / 1073741824.0 << std::endl;

        const uint32_t sizeQ = params->GetElementParams()->GetParams().size();
        std::vector<double> latencies;
        double worst = 0;
        std::vector<int> threadList{0};  // 0: library default
        if (const char* env = std::getenv("SHIP_THREADS")) {
            threadList.clear();
            std::stringstream ss(env);
            for (std::string item; std::getline(ss, item, ',');)
                threadList.push_back(std::stoi(item));
        }
#ifdef _OPENMP
        std::cout << "openmp=on max_threads=" << omp_get_max_threads() << std::endl;
#else
        std::cout << "openmp=off" << std::endl;
#endif
        for (int threads : threadList) {
        (void)threads;
#ifdef _OPENMP
        if (threads > 0)
            omp_set_num_threads(threads);
        std::cout << "threads=" << (threads > 0 ? threads : omp_get_max_threads()) << std::endl;
#endif
        latencies.clear();
        worst = 0;
        std::mt19937 rng(20261011);
        std::uniform_real_distribution<double> u(-1, 1);
        for (uint32_t t = 0; t < trials; ++t) {
            std::vector<double> x(S);
            for (auto& v : x)
                v = u(rng);
            auto ct = cc->Encrypt(kp.publicKey, cc->MakeCKKSPackedPlaintext(x, 1, sizeQ - 2, nullptr, S));
            start          = Clock::now();
            auto out       = cc->EvalSHIPBootstrap(ct);
            const double s = Seconds(start);
            Plaintext pt;
            cc->Decrypt(kp.secretKey, out, &pt);
            pt->SetLength(S);
            double eps = 0;
            for (uint32_t i = 0; i < S; ++i)
                eps = std::max(eps, std::abs(pt->GetCKKSPackedValue()[i].real() - x[i]));
            worst = std::max(worst, eps);
            if (set == "custom") {
                // systematic (scale) vs. random error: least-squares ratio and residual
                double xy = 0, xx = 0;
                for (uint32_t i = 0; i < S; ++i) {
                    xy += pt->GetCKKSPackedValue()[i].real() * x[i];
                    xx += x[i] * x[i];
                }
                const double ratio = xy / xx;
                double residual = 0, imag = 0;
                for (uint32_t i = 0; i < S; ++i) {
                    residual = std::max(residual, std::abs(pt->GetCKKSPackedValue()[i].real() - ratio * x[i]));
                    imag     = std::max(imag, std::abs(pt->GetCKKSPackedValue()[i].imag()));
                }
                size_t at = 0; double top = 0, rest = 0;
                for (uint32_t i = 0; i < S; ++i) {
                    const double e = std::abs(pt->GetCKKSPackedValue()[i].real() - x[i]);
                    if (e > top) { top = e; at = i; }
                    if (i != 0) rest = std::max(rest, e);
                }
                std::cout << "diag worst_slot=" << at << " error_slot0=" << std::abs(pt->GetCKKSPackedValue()[0].real() - x[0])
                          << " max_error_excluding_slot0=" << rest << std::endl;
                std::cout << "diag ratio=" << ratio << " residual_after_ratio=" << residual << " max_imag=" << imag
                          << " slot0=" << pt->GetCKKSPackedValue()[0] << " x0=" << x[0] << std::endl;
            }
            latencies.push_back(s);
            const uint32_t limbs = out->GetElements()[0].GetNumOfElements();
            std::cout << "trial=" << t << " latency_s=" << s << " max_real_error=" << eps
                      << " precision_bits=" << -std::log2(eps) << " output_limbs=" << limbs
                      << " mult_levels_after=" << limbs - 2 << std::endl;
            if (t == 0 && set != "custom" && threads == threadList.front()) {
                // IND-CPA^D decryption with noise flooding, emulated exactly as OpenFHE's NOISE_FLOODING_DECRYPT:
                // noise estimate from the imaginary parts (EXEC_NOISE_ESTIMATION), then Gaussian noise of
                // standard deviation 2^(stat/2 + log2(sqrt(12 q)) + estimate) added to the decryption (here to c0).
                double sumSq = 0;
                for (uint32_t i = 0; i < S; ++i)
                    sumSq += std::pow(pt->GetCKKSPackedValue()[i].imag(), 2);
                const double sigmaSlot  = std::sqrt(sumSq / S);
                const double delta      = std::pow(2.0, params->GetPlaintextModulus());
                const double estimate   = std::log2(sigmaSlot * delta / std::sqrt(N / 2.0));
                const double statSec    = 30, queries = 1;  // OpenFHE defaults
                const double logStd     = statSec / 2 + std::log2(std::sqrt(12 * queries));
                DiscreteGaussianGeneratorImpl<NativeVector> flood(std::pow(2.0, logStd + estimate));
                auto flooded  = out->Clone();
                auto elements = flooded->GetElements();
                // The flooding samples can exceed the small SHIP primes: reduce them modulo every prime explicitly.
                auto samples = flood.GenerateIntVector(N);
                DCRTPoly noise(elements[0].GetParams(), Format::COEFFICIENT, true);
                for (size_t k = 0; k < noise.GetNumOfElements(); ++k) {
                    auto limb        = noise.GetElementAtIndex(k);
                    const int64_t q  = static_cast<int64_t>(limb.GetModulus().ConvertToInt<uint64_t>());
                    for (uint32_t i = 0; i < N; ++i)
                        limb[i] = NativeInteger(static_cast<uint64_t>(((samples.get()[i] % q) + q) % q));
                    noise.SetElementAtIndex(k, std::move(limb));
                }
                noise.SetFormat(Format::EVALUATION);
                elements[0] += noise;
                flooded->SetElements(std::move(elements));
                Plaintext pf;
                cc->Decrypt(kp.secretKey, flooded, &pf);
                pf->SetLength(S);
                double ef = 0;
                for (uint32_t i = 0; i < S; ++i)
                    ef = std::max(ef, std::abs(pf->GetCKKSPackedValue()[i].real() - x[i]));
                std::cout << "noise_flooding: estimate_log2=" << estimate << " flooding_std_log2=" << logStd + estimate
                          << " scale_log2=" << std::log2(delta) << " max_real_error=" << ef
                          << " precision_bits=" << -std::log2(ef) << std::endl;
                // The restored level must be usable: one multiplication, then bootstrap again.
                auto sq = cc->EvalMult(out, out);
                cc->RescaleInPlace(sq);
                auto again = cc->EvalSHIPBootstrap(sq);
                Plaintext p2;
                cc->Decrypt(kp.secretKey, again, &p2);
                p2->SetLength(S);
                double e2 = 0;
                for (uint32_t i = 0; i < S; ++i)
                    e2 = std::max(e2, std::abs(p2->GetCKKSPackedValue()[i].real() - x[i] * x[i]));
                std::cout << "square_then_bootstrap max_real_error=" << e2 << " precision_bits=" << -std::log2(e2)
                          << std::endl;
            }
        }
        std::sort(latencies.begin(), latencies.end());
        std::cout << "summary set=" << set << " trials=" << trials << " median_latency_s=" << latencies[trials / 2]
                  << " worst_precision_bits=" << -std::log2(worst) << std::endl;
        }
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "FAIL " << e.what() << std::endl;
        return 1;
    }
}
