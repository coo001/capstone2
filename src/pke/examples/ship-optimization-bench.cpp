// Controlled comparison of our two implementations, NOT SHIP vs standard CKKS.
#include "ship/full-bootstrap.h"
#include <chrono>
#include <iostream>
#include <set>
#include <string>

using namespace ship;
using Clock = std::chrono::steady_clock;

static double Error(const CC& cc, const SK& sk, const CT& ct, const std::vector<C>& expected) {
    Plaintext p;
    cc->Decrypt(sk,ct,&p); p->SetLength(expected.size());
    double result = 0;
    for (size_t i = 0; i < expected.size(); ++i) {
        const auto error = std::abs(p->GetCKKSPackedValue()[i]-expected[i]);
        if (!std::isfinite(error)) throw std::runtime_error("non-finite benchmark result");
        result = std::max(result,error);
    }
    return result;
}

// Coefficient payload only: excludes allocator/container overhead and context tables.
static uint64_t PolyBytes(const DCRTPoly& p) {
    return uint64_t(p.GetRingDimension())*p.GetNumOfElements()*sizeof(NativeInteger);
}
static uint64_t StoredBytes(const CC& cc, const FullBootstrapKey& k, const std::string& tag) {
    uint64_t total = 0;
    std::set<const void*> seen;
    auto addEval = [&](const EvalKey<DCRTPoly>& key) {
        if (!seen.insert(key.get()).second) return;
        for (const auto& p : key->GetAVector()) total += PolyBytes(p);
        for (const auto& p : key->GetBVector()) total += PolyBytes(p);
    };
    for (const auto& digit : k.encapsulation.digits)
        for (const auto& p : digit) total += uint64_t(p.GetLength())*sizeof(NativeInteger);
    for (const auto& f : k.half.factors) {
        for (const auto& ct : f.selectors) if (ct) for (const auto& p : ct->GetElements()) total += PolyBytes(p);
        if (k.half.auxMasking) for (const auto& pair : f.auxSelectors) for (const auto& p : pair) total += PolyBytes(p);
        for (const auto& pair : f.rotation) for (const auto& key : pair) { addEval(key.body); addEval(key.mask); }
        for (const auto& pair : f.fusedRotation) for (const auto& key : pair) { addEval(key.body); addEval(key.mask); }
    }
    for (const auto& [index,key] : *k.half.conjugation) addEval(key);
    for (const auto& [index,key] : cc->GetEvalAutomorphismKeyMap(tag)) addEval(key);
    for (const auto& p : k.packing.diagonals) total += PolyBytes(p->GetElement<DCRTPoly>());
    return total;
}

struct Run {
    std::string mode;
    FullBootstrapKey keys;
    std::vector<double> times;
    double worst = 0;
    CT last;
};

int main(int argc, char** argv) {
    try {
        const std::string mode = argc > 1 ? argv[1] : "compare";
        if (mode != "compare" && mode != "reference" && mode != "fused")
            throw std::invalid_argument("usage: ship-optimization-bench [compare|reference|fused]");
        const uint32_t n = 1024, slots = n/2, depth = 9, h = 8, repeats = 7;
        CCParams<CryptoContextCKKSRNS> p;
        p.SetSecurityLevel(HEStd_NotSet); p.SetRingDim(n); p.SetBatchSize(slots);
        p.SetMultiplicativeDepth(depth); p.SetScalingModSize(50); p.SetFirstModSize(60);
        p.SetScalingTechnique(FIXEDMANUAL); p.SetKeySwitchTechnique(HYBRID); p.SetCKKSDataType(COMPLEX);
        auto cc = GenCryptoContext(p);
        cc->Enable(PKE); cc->Enable(KEYSWITCH); cc->Enable(LEVELEDSHE); cc->Enable(ADVANCEDSHE);
        auto dense = cc->KeyGen(); cc->EvalMultKeyGen(dense.secretKey);
        DCRTPoly::TugType ternary;
        auto sparse = std::make_shared<PrivateKeyImpl<DCRTPoly>>(cc);
        sparse->SetPrivateElement(DCRTPoly(ternary,dense.secretKey->GetPrivateElement().GetParams(),Format::EVALUATION,h));
        std::vector<C> x(slots);
        for (uint32_t i = 0; i < slots; ++i) x[i] = C((int((37*i)%101)-50)/100.0,(int((61*i)%103)-51)/102.0);
        auto input = cc->Encrypt(dense.publicKey,cc->MakeCKKSPackedPlaintext(x,1,depth,nullptr,slots));
        std::vector<Run> runs;
        for (const auto& name : {std::string("reference"),std::string("fused")}) {
            if (mode != "compare" && mode != name) continue;
            const auto start = Clock::now();
            auto keys = MakeFullBootstrapKey(cc,dense,sparse,name == "fused");
            const double setup = std::chrono::duration<double,std::milli>(Clock::now()-start).count();
            // In compare mode both variants share the context rotation cache, so
            // report memory only for standalone runs with an uncontaminated cache.
            if (mode != "compare")
                std::cout << "mode=" << name << " stored_payload_bytes=" << StoredBytes(cc,keys,dense.secretKey->GetKeyTag())
                          << " setup_ms=" << setup << '\n';
            runs.push_back({name,std::move(keys),{},0,nullptr});
        }
        // One warm-up per variant; alternate evaluation order on measured pairs.
        for (auto& run : runs) {
            auto warmup = FullBootstrap(cc,input,run.keys);
            if (Error(cc,dense.secretKey,warmup,x) >= 1e-4) throw std::runtime_error("warm-up accuracy failure");
        }
        for (uint32_t r = 0; r < repeats; ++r) {
            for (size_t j = 0; j < runs.size(); ++j) {
                auto& run = runs[(r%2) ? runs.size()-1-j : j];
                const auto start = Clock::now();
                run.last = FullBootstrap(cc,input,run.keys);
                const double elapsed = std::chrono::duration<double,std::milli>(Clock::now()-start).count();
                run.times.push_back(elapsed);
                run.worst = std::max(run.worst,Error(cc,dense.secretKey,run.last,x));
                if (run.worst >= 1e-4) throw std::runtime_error("benchmark accuracy failure");
                std::cout << "sample mode=" << run.mode << " repeat=" << r << " ms=" << elapsed << '\n';
            }
        }
        for (auto& run : runs) {
            std::sort(run.times.begin(),run.times.end());
            std::cout << "summary mode=" << run.mode << " N=" << n << " h=" << h << " slots=" << slots
                      << " median_ms=" << run.times[repeats/2] << " min_ms=" << run.times.front()
                      << " max_ms=" << run.times.back() << " max_error=" << run.worst
                      << " output_towers=" << run.last->GetElements()[0].GetNumOfElements() << '\n';
        }
        if (runs.size() == 2) {
            const double delta = Error(cc,dense.secretKey,cc->EvalSub(runs[0].last,runs[1].last),std::vector<C>(slots));
            if (delta >= 5e-6) throw std::runtime_error("full-bootstrap differential mismatch");
            std::cout << "comparison speedup=" << runs[0].times[repeats/2]/runs[1].times[repeats/2]
                      << " output_difference=" << delta << '\n';
        }
        std::cout << "BENCHMARK CHECKS PASSED (toy parameters; not a paper or baseline performance claim)\n";
        return 0;
    } catch (const std::exception& e) { std::cerr << "FAIL " << e.what() << '\n'; return 1; }
}
