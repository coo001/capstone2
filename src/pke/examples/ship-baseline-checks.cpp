// Standard OpenFHE bootstrap regression; separate from the SHIP prototype.
#include "openfhe.h"
#include <algorithm>
#include <cmath>
#include <iostream>

using namespace lbcrypto;
int main() {
    try {
        CCParams<CryptoContextCKKSRNS> p;
        p.SetSecretKeyDist(UNIFORM_TERNARY);
        p.SetSecurityLevel(HEStd_NotSet);
        p.SetRingDim(4096);
        p.SetScalingModSize(59);
        p.SetFirstModSize(60);
        p.SetScalingTechnique(FLEXIBLEAUTO);
        std::vector<uint32_t> budget{4,4};
        uint32_t depth = 10 + FHECKKSRNS::GetBootstrapDepth(budget, UNIFORM_TERNARY);
        p.SetMultiplicativeDepth(depth);
        auto cc = GenCryptoContext(p);
        cc->Enable(PKE); cc->Enable(KEYSWITCH); cc->Enable(LEVELEDSHE);
        cc->Enable(ADVANCEDSHE); cc->Enable(FHE);
        uint32_t slots = cc->GetRingDimension()/2;
        cc->EvalBootstrapSetup(budget);
        auto kp = cc->KeyGen();
        cc->EvalMultKeyGen(kp.secretKey);
        cc->EvalBootstrapKeyGen(kp.secretKey, slots);
        std::vector<double> x(slots);
        for (size_t i = 0; i < x.size(); ++i) x[i] = (int(i%17)-8)/16.0;
        auto input = cc->Encrypt(kp.publicKey, cc->MakeCKKSPackedPlaintext(x,1,depth-1));
        auto output = cc->EvalBootstrap(input);
        Plaintext pt;
        cc->Decrypt(kp.secretKey,output,&pt);
        pt->SetLength(slots);
        double error = 0;
        for (size_t i = 0; i < x.size(); ++i) {
            auto v = pt->GetCKKSPackedValue()[i];
            if (!std::isfinite(v.real()) || !std::isfinite(v.imag())) throw std::runtime_error("non-finite baseline");
            error = std::max(error,std::abs(v-std::complex<double>(x[i],0)));
        }
        int before = depth-input->GetLevel()-(input->GetNoiseScaleDeg()-1);
        int after = depth-output->GetLevel()-(output->GetNoiseScaleDeg()-1);
        std::cout << "standard_bootstrap slots=" << slots << " max_error=" << error
                  << " usable_levels_before=" << before << " usable_levels_after=" << after << std::endl;
        if (error >= 1e-3 || after <= before) throw std::runtime_error("standard bootstrap regression");
        std::cout << "ALL BASELINE CHECKS PASSED (upstream toy parameters)\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL " << e.what() << '\n';
        return 1;
    }
}
