// OpenFHE's standard CKKS bootstrapping (EvalBootstrap) with 128-bit HE-standard parameters, as the baseline
// for ship-paper-bench (paper Section 5.3 compares SHIP with the conventional bootstrapping of HEaaN).
// usage: ship-baseline-bench [scalingModBits] [firstModBits] [trials] [ringDim(0 = chosen by security)] [dnum] [levelBudget]
//                            [uniform|sparse_encapsulated]
// sparse_encapsulated: OpenFHE's sparse-secret encapsulation (h = 32 key modulo q0 * p0), the same idea as SHIP's.
#include "openfhe.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <random>
#include <sstream>
#ifdef _OPENMP
    #include <omp.h>
#endif

using namespace lbcrypto;
using Clock = std::chrono::steady_clock;

int main(int argc, char** argv) {
    try {
        const uint32_t scaleBits = argc > 1 ? std::stoul(argv[1]) : 50;
        const uint32_t firstBits = argc > 2 ? std::stoul(argv[2]) : 60;
        const uint32_t trials    = argc > 3 ? std::stoul(argv[3]) : 3;
        const uint32_t ringDim   = argc > 4 ? std::stoul(argv[4]) : 0;
        const uint32_t dnum      = argc > 5 ? std::stoul(argv[5]) : 0;
        const uint32_t budget    = argc > 6 ? std::stoul(argv[6]) : 3;
        const std::vector<uint32_t> levelBudget{budget, budget};
        const std::string distName = argc > 7 ? argv[7] : "uniform";
        const SecretKeyDist dist   = distName == "sparse_encapsulated" ? SPARSE_ENCAPSULATED : UNIFORM_TERNARY;
        const uint32_t levelsAfter = 1;  // same number of multiplicative levels as SHIP LL13 / LL14
        CCParams<CryptoContextCKKSRNS> p;
        p.SetSecretKeyDist(dist);
        p.SetSecurityLevel(HEStd_128_classic);
        if (ringDim)
            p.SetRingDim(ringDim);
        if (dnum)
            p.SetNumLargeDigits(dnum);
        p.SetScalingModSize(scaleBits);
        p.SetFirstModSize(firstBits);
        p.SetScalingTechnique(FLEXIBLEAUTO);
        p.SetKeySwitchTechnique(HYBRID);
        p.SetCKKSDataType(REAL);
        const uint32_t depth = levelsAfter + FHECKKSRNS::GetBootstrapDepth(levelBudget, dist);
        p.SetMultiplicativeDepth(depth);
        std::cout << "requested_depth=" << depth << std::endl;
        auto cc = GenCryptoContext(p);
        cc->Enable(PKE);
        cc->Enable(KEYSWITCH);
        cc->Enable(LEVELEDSHE);
        cc->Enable(ADVANCEDSHE);
        cc->Enable(FHE);
        const uint32_t N     = cc->GetRingDimension();
        const uint32_t slots = N / 2;
        auto params          = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc->GetCryptoParameters());
        std::cout << "baseline N=" << N << " slots=" << slots << " depth=" << depth
                  << " log2QP=" << params->GetParamsQP()->GetModulus().GetMSB()
                  << " HEStd_128_classic_max_log2QP=" << StdLatticeParm::FindMaxQ(HEStd_ternary, HEStd_128_classic, N)
                  << " secret=" << distName
                  << " log2(q0*p0)=" << params->GetElementParams()->GetParams()[0]->GetModulus().GetMSB() +
                                            params->GetParamsP()->GetParams()[0]->GetModulus().GetMSB()
                  << " scale_bits=" << scaleBits << " first_bits=" << firstBits << " levelBudget=" << budget << "," << budget << std::endl;
        auto start = Clock::now();
        cc->EvalBootstrapSetup(levelBudget, {0, 0}, slots);
        auto kp = cc->KeyGen();
        cc->EvalMultKeyGen(kp.secretKey);
        cc->EvalBootstrapKeyGen(kp.secretKey, slots);
        std::cout << "setup_keygen_s=" << std::chrono::duration<double>(Clock::now() - start).count() << std::endl;
        std::vector<int> threadList{0};
        if (const char* env = std::getenv("SHIP_THREADS")) {
            threadList.clear();
            std::stringstream ss(env);
            for (std::string item; std::getline(ss, item, ',');)
                threadList.push_back(std::stoi(item));
        }
        for (int threads : threadList) {
        (void)threads;
#ifdef _OPENMP
        if (threads > 0)
            omp_set_num_threads(threads);
        std::cout << "threads=" << (threads > 0 ? threads : omp_get_max_threads()) << std::endl;
#endif
        std::mt19937 rng(20261011);
        std::uniform_real_distribution<double> u(-1, 1);
        std::vector<double> latencies;
        double worst = 0;
        for (uint32_t t = 0; t < trials; ++t) {
            std::vector<double> x(slots);
            for (auto& v : x)
                v = u(rng);
            auto ct  = cc->Encrypt(kp.publicKey, cc->MakeCKKSPackedPlaintext(x, 1, depth - 1, nullptr, slots));
            start    = Clock::now();
            auto out = cc->EvalBootstrap(ct);
            const double s = std::chrono::duration<double>(Clock::now() - start).count();
            Plaintext pt;
            cc->Decrypt(kp.secretKey, out, &pt);
            pt->SetLength(slots);
            double eps = 0;
            for (uint32_t i = 0; i < slots; ++i)
                eps = std::max(eps, std::abs(pt->GetCKKSPackedValue()[i].real() - x[i]));
            worst = std::max(worst, eps);
            latencies.push_back(s);
            const uint32_t levels = depth - out->GetLevel() - (out->GetNoiseScaleDeg() - 1);
            std::cout << "trial=" << t << " latency_s=" << s << " max_real_error=" << eps
                      << " precision_bits=" << -std::log2(eps) << " levels_after=" << levels << std::endl;
        }
        std::sort(latencies.begin(), latencies.end());
        std::cout << "summary baseline N=" << N << " median_latency_s=" << latencies[trials / 2]
                  << " worst_precision_bits=" << -std::log2(worst) << std::endl;
        }
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "FAIL " << e.what() << std::endl;
        return 1;
    }
}
