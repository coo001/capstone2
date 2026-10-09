#include "ship/half-bootstrap.h"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <random>

using namespace lbcrypto;
using namespace ship;

int main() {
    try {
        CCParams<CryptoContextCKKSRNS> parameters;
        parameters.SetSecurityLevel(HEStd_NotSet);  // Correctness fixture only.
        parameters.SetRingDim(1024);
        parameters.SetBatchSize(512);
        parameters.SetMultiplicativeDepth(8);
        parameters.SetScalingModSize(55);
        parameters.SetFirstModSize(60);
        parameters.SetScalingTechnique(FIXEDMANUAL);
        parameters.SetCKKSDataType(COMPLEX);
        parameters.SetKeySwitchTechnique(HYBRID);
        auto cc = GenCryptoContext(parameters);
        cc->Enable(PKE); cc->Enable(KEYSWITCH); cc->Enable(LEVELEDSHE); cc->Enable(ADVANCEDSHE);
        auto dense = cc->KeyGen();
        cc->EvalMultKeyGen(dense.secretKey);
        const uint32_t slots = cc->GetRingDimension()/2;
        // Includes both signs and both negacyclic/half-vector wrap boundaries.
        const std::vector<std::pair<uint32_t,int>> support = {
            {0, 1}, {slots-1, -1}, {slots, 1}, {2*slots-1, -1}};
        auto sparse = std::make_shared<PrivateKeyImpl<DCRTPoly>>(cc);
        DCRTPoly sparsePoly(dense.secretKey->GetPrivateElement().GetParams(), Format::COEFFICIENT, true);
        for (size_t limb = 0; limb < sparsePoly.GetNumOfElements(); ++limb) {
            auto p = sparsePoly.GetElementAtIndex(limb);
            for (const auto& [j,s] : support) p[j] = s > 0 ? NativeInteger(1) : p.GetModulus()-NativeInteger(1);
            sparsePoly.SetElementAtIndex(limb, std::move(p));
        }
        sparsePoly.SetFormat(Format::EVALUATION);
        sparse->SetPrivateElement(std::move(sparsePoly));
        const double gamma = 1024;
        auto keys = MakeHalfBootstrapKey(cc, dense, support, sparse->GetKeyTag(), gamma);
        std::mt19937 rng(20251022);
        for (uint32_t trial = 0; trial < 3; ++trial) {
            std::vector<double> message(2*slots);
            for (size_t i = 0; i < message.size(); ++i)
                message[i] = trial == 0 ? 0.0 : trial == 1 ? (int(i%9)-4)/8.0 : (int(rng()%1001)-500)/1000.0;
            // A genuine OpenFHE symmetric RLWE encryption of coefficient-domain data.
            // This fixture starts under the sparse key; dense->sparse encapsulation is not implemented here.
            auto pt = cc->MakeCKKSPackedPlaintext(std::vector<double>(slots,0), 1, 8, nullptr, slots);
            auto poly = pt->GetElement<DCRTPoly>();
            poly.SetFormat(Format::COEFFICIENT);
            auto limb = poly.GetElementAtIndex(0);
            const uint64_t q0 = limb.GetModulus().ConvertToInt();
            for (size_t i = 0; i < message.size(); ++i) {
                int64_t value = std::llround(static_cast<long double>(q0)*message[i]/gamma);
                limb[i] = NativeInteger(value < 0 ? q0-static_cast<uint64_t>(-value) : static_cast<uint64_t>(value));
            }
            poly.SetElementAtIndex(0, std::move(limb));
            poly.SetFormat(Format::EVALUATION);
            pt->GetElement<DCRTPoly>() = std::move(poly);
            pt->SetScalingFactor(static_cast<double>(q0)/gamma);
            auto input = cc->Encrypt(sparse, pt);
            auto before = input->GetElements();
            auto start = std::chrono::steady_clock::now();
            auto output = HalfBootstrap(cc, input, keys);
            auto ms = std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
            Plaintext decoded;
            cc->Decrypt(dense.secretKey, output, &decoded);
            decoded->SetLength(slots);
            double messageError = 0, phaseError = 0, imaginaryError = 0;
            const double pi = std::acos(-1.0);
            for (size_t i = 0; i < slots; ++i) {
                auto z = decoded->GetCKKSPackedValue()[i];
                if (!std::isfinite(z.real()) || !std::isfinite(z.imag())) throw std::runtime_error("non-finite output");
                double sineTarget = gamma/(2*pi)*std::sin(2*pi*message[i]/gamma);
                messageError = std::max(messageError,std::abs(z.real()-message[i]));
                phaseError = std::max(phaseError,std::abs(z.real()-sineTarget));
                imaginaryError = std::max(imaginaryError,std::abs(z.imag()));
            }
            std::cout << "trial=" << trial << " message_error=" << messageError
                      << " phase_error=" << phaseError << " imaginary_error=" << imaginaryError
                      << " input_towers=" << input->GetElements()[0].GetNumOfElements()
                      << " output_towers=" << output->GetElements()[0].GetNumOfElements()
                      << " output_level=" << output->GetLevel() << " elapsed_ms=" << ms << std::endl;
            if (messageError >= 2e-6 || phaseError >= 2e-8 || imaginaryError >= 2e-8)
                throw std::runtime_error("half-bootstrap exceeded error threshold");
            if (before != input->GetElements()) throw std::runtime_error("input mutated");
            if (output->GetElements()[0].GetNumOfElements() <= 1) throw std::runtime_error("no budget restored");
            // Reject mismatched scale and key tags; coefficient layout remains a caller contract.
            auto invalid = input->Clone();
            invalid->SetScalingFactor(2*input->GetScalingFactor());
            bool rejected = false;
            try { HalfBootstrap(cc, invalid, keys); } catch (const std::invalid_argument&) { rejected = true; }
            if (!rejected) throw std::runtime_error("invalid input scale was accepted");
            invalid = input->Clone();
            invalid->SetKeyTag(dense.secretKey->GetKeyTag());
            rejected = false;
            try { HalfBootstrap(cc, invalid, keys); } catch (const std::invalid_argument&) { rejected = true; }
            if (!rejected) throw std::runtime_error("invalid input key was accepted");
        }
        std::cout << "ALL HALF-BOOTSTRAP CHECKS PASSED (toy sparse key; no security or latency claim)\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL " << e.what() << '\n';
        return 1;
    }
}
