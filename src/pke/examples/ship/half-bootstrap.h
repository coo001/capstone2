// Correctness prototype for coefficient -> slot SHIP half bootstrapping.
// Requires FIXEDMANUAL, COMPLEX, full packing and a one-limb sparse-key input.
// The fused path is OpenFHE-specific; paper equivalence/security remain unverified.
#pragma once
#include "fused-rotation.h"
#include <cmath>
#include <set>

namespace ship {
using QPPair = std::array<DCRTPoly, 2>;  // (b,a) over the extended basis QP

struct FactorKey {
    std::array<CT, 4> selectors;          // rescaling path: Enc_Q(Delta * mask)
    std::array<QPPair, 4> auxSelectors;   // auxiliary-modulus path: Enc_QP(P * mask)
    BlindKey rotation;
    FusedBlindKey fusedRotation;
};
struct HalfBootstrapKey {
    std::vector<FactorKey> factors;
    PublicKey<DCRTPoly> outputPublicKey;
    std::shared_ptr<std::map<uint32_t, EvalKey<DCRTPoly>>> conjugation;
    std::string inputKeyTag;
    uint32_t slots;
    double gamma;
    bool useFused;
    bool auxMasking;
};

// Number of HE levels consumed before the product tree (paper Alg. 1 vs Sec. 4.4).
inline uint32_t MaskingLevels(bool auxMasking) { return auxMasking ? 0 : 1; }

inline std::shared_ptr<CryptoParametersCKKSRNS> CkksParams(const CC& cc) {
    auto params = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc->GetCryptoParameters());
    if (!params) throw std::invalid_argument("CKKS RNS parameters required");
    return params;
}

// Signed lift of the bottom limb of a small integer polynomial to every limb of QP.
// Exact only while all centered coefficients stay below q0/4, which is checked.
inline DCRTPoly LiftSmallToQP(const DCRTPoly& poly, const std::shared_ptr<DCRTPoly::Params>& paramsQP) {
    auto limb = poly.GetElementAtIndex(0);
    limb.SetFormat(Format::COEFFICIENT);
    const auto q0 = limb.GetModulus();
    const auto quarter = q0 >> 2;
    for (size_t i = 0; i < limb.GetLength(); ++i)
        if (limb[i] >= quarter && limb[i] <= q0-quarter)
            throw std::range_error("plaintext coefficient too large for an exact QP lift");
    DCRTPoly result(paramsQP, Format::COEFFICIENT, true);
    for (size_t k = 0; k < result.GetNumOfElements(); ++k) {
        auto copy = limb;
        const auto& target = paramsQP->GetParams()[k];
        if (target->GetModulus() != q0)
            copy.SwitchModulus(target->GetModulus(), target->GetRootOfUnity(), 0, 0);
        result.SetElementAtIndex(k, std::move(copy));
    }
    result.SetFormat(Format::EVALUATION);
    return result;
}

// Same secret as in KeySwitchHYBRID::KeySwitchGenInternal, extended from Q to QP.
inline DCRTPoly SecretOverQP(const SK& sk, const std::shared_ptr<DCRTPoly::Params>& paramsQP) {
    auto secret = sk->GetPrivateElement().Clone();
    secret.SetFormat(Format::COEFFICIENT);
    DCRTPoly result(paramsQP, Format::COEFFICIENT, true);
    const size_t sizeQ = secret.GetNumOfElements();
    for (size_t i = 0; i < sizeQ; ++i) result.SetElementAtIndex(i, secret.GetElementAtIndex(i));
    for (size_t j = sizeQ; j < result.GetNumOfElements(); ++j) {
        auto limb = secret.GetElementAtIndex(0);
        const auto& target = paramsQP->GetParams()[j];
        limb.SwitchModulus(target->GetModulus(), target->GetRootOfUnity(), 0, 0);
        result.SetElementAtIndex(j, std::move(limb));
    }
    result.SetFormat(Format::EVALUATION);
    return result;
}

// Enc_QP(P * mask) under the output secret (given over QP). The mask is CKKS slot-encoded with scale P instead of Delta,
// so that a later product with a Delta-scaled plaintext and ModDown by P keeps scale Delta.
// Coefficients: round(m_i * P / Delta) from the exact Delta-encoding m_i (relative error <= 2^-log2(Delta)).
inline QPPair EncryptMaskOverQP(const CC& cc, const DCRTPoly& secretQP,
                                const std::vector<double>& mask, uint32_t slots) {
    const auto params = CkksParams(cc);
    const auto paramsQP = params->GetParamsQP();
    const uint32_t bottom = params->GetElementParams()->GetParams().size()-1;
    auto pt = cc->MakeCKKSPackedPlaintext(mask, 1, bottom, nullptr, slots);
    int exponent = 0;
    const double mantissa = std::frexp(pt->GetScalingFactor(), &exponent);
    if (mantissa != 0.5 || exponent < 2)
        throw std::invalid_argument("auxiliary masking requires a power-of-two scaling factor");
    const uint32_t deltaBits = exponent-1;
    auto encoded = pt->GetElement<DCRTPoly>();
    encoded.SetFormat(Format::COEFFICIENT);
    const auto limb = encoded.GetElementAtIndex(0);
    const auto q0 = limb.GetModulus();
    const BigInteger P = params->GetParamsP()->GetModulus();
    const BigInteger half = BigInteger(1).LShift(deltaBits-1);
    DCRTPoly plain(paramsQP, Format::COEFFICIENT, true);
    std::vector<NativePoly> limbs;
    for (size_t k = 0; k < plain.GetNumOfElements(); ++k) limbs.push_back(plain.GetElementAtIndex(k));
    for (size_t i = 0; i < limb.GetLength(); ++i) {
        const bool negative = limb[i] > (q0 >> 1);
        const uint64_t magnitude = (negative ? q0-limb[i] : limb[i]).ConvertToInt<uint64_t>();
        const BigInteger scaled = (BigInteger(magnitude)*P + half).RShift(deltaBits);
        for (auto& out : limbs) {
            const auto r = out.GetModulus();
            NativeInteger value(scaled.Mod(BigInteger(r.ConvertToInt<uint64_t>())).ConvertToInt<uint64_t>());
            out[i] = (negative && value != NativeInteger(0)) ? r-value : value;
        }
    }
    for (size_t k = 0; k < limbs.size(); ++k) plain.SetElementAtIndex(k, std::move(limbs[k]));
    plain.SetFormat(Format::EVALUATION);
    DCRTPoly::DugType uniform;
    DCRTPoly a(uniform, paramsQP, Format::EVALUATION);
    DCRTPoly e(params->GetDiscreteGaussianGenerator(), paramsQP, Format::EVALUATION);
    DCRTPoly b = e - a*secretQP + plain;
    return {std::move(b), std::move(a)};
}

// Only setup has access to sparse support/signs and the dense output secret.
inline HalfBootstrapKey MakeHalfBootstrapKey(
    const CC& cc, const KeyPair<DCRTPoly>& outputKey,
    const std::vector<std::pair<uint32_t, int>>& sparseSupport,
    const std::string& inputKeyTag, double gamma, bool useFused = true, bool auxMasking = true) {
    const uint32_t slots = cc->GetRingDimension() / 2;
    if (gamma <= 0 || !std::isfinite(gamma) || sparseSupport.empty())
        throw std::invalid_argument("invalid half-bootstrap parameters");
    const auto parameters = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc->GetCryptoParameters());
    if (!parameters || parameters->GetScalingTechnique() != FIXEDMANUAL ||
        parameters->GetKeySwitchTechnique() != HYBRID)
        throw std::invalid_argument("reference half-bootstrap requires FIXEDMANUAL and HYBRID");
    if (outputKey.secretKey->GetCryptoContext() != cc || outputKey.publicKey->GetCryptoContext() != cc)
        throw std::invalid_argument("output key context mismatch");
    std::set<uint32_t> distinct;
    for (const auto& [position, sign] : sparseSupport) {
        if (position >= 2*slots || (sign != -1 && sign != 1) || !distinct.insert(position).second)
            throw std::invalid_argument("invalid or duplicate sparse support");
    }
    HalfBootstrapKey result;
    result.outputPublicKey = outputKey.publicKey;
    result.inputKeyTag = inputKeyTag;
    result.slots = slots;
    result.gamma = gamma;
    result.useFused = useFused;
    result.auxMasking = auxMasking;
    const auto paramsQP = parameters->GetParamsQP();
    const DCRTPoly secretQP = auxMasking ? SecretOverQP(outputKey.secretKey, paramsQP) : DCRTPoly();
    std::vector<int32_t> rotations;
    for (uint32_t s = 1; s < slots; s <<= 1) rotations.push_back(-static_cast<int32_t>(s));
    if (!useFused) cc->EvalRotateKeyGen(outputKey.secretKey, rotations);
    const uint32_t conjugationIndex = 4*slots-1;
    cc->EvalAutomorphismKeyGen(outputKey.secretKey, {conjugationIndex});
    // Repeated setup can return no newly generated keys when the index is cached.
    // Keep our own one-entry map instead of relying on that incremental result.
    result.conjugation = std::make_shared<std::map<uint32_t,EvalKey<DCRTPoly>>>();
    result.conjugation->emplace(conjugationIndex,
        cc->GetEvalAutomorphismKeyMap(outputKey.secretKey->GetKeyTag()).at(conjugationIndex));
    for (const auto& [position, sign] : sparseSupport) {
        FactorKey factor;
        auto masks = Masks(position, sign, slots);
        for (size_t band = 0; band < 4; ++band) {
            if (auxMasking)
                factor.auxSelectors[band] = EncryptMaskOverQP(cc, secretQP, masks[band], slots);
            else
                factor.selectors[band] = cc->Encrypt(outputKey.publicKey,
                    cc->MakeCKKSPackedPlaintext(masks[band], 1, 0, nullptr, slots));
        }
        if (useFused)
            factor.fusedRotation = MakeFusedBlindKey(cc, outputKey.secretKey, position % slots, slots, -1);
        else
            factor.rotation = MakeBlindKey(cc, outputKey.secretKey, position % slots, slots, -1);
        result.factors.push_back(std::move(factor));
    }
    return result;
}

inline CT HalfBootstrap(const CC& cc, const CT& input, const HalfBootstrapKey& key) {
    if (!input || input->GetElements().size() != 2 || input->GetKeyTag() != key.inputKeyTag)
        throw std::invalid_argument("half-bootstrap input/key mismatch");
    if (input->GetElements()[0].GetNumOfElements() != 1 ||
        input->GetElements()[1].GetNumOfElements() != 1 ||
        input->GetElements()[0].GetRingDimension() != 2 * key.slots)
        throw std::invalid_argument("half-bootstrap requires a one-limb, full-ring input");
    if (input->GetCryptoContext() != cc || key.outputPublicKey->GetCryptoContext() != cc)
        throw std::invalid_argument("half-bootstrap context mismatch");
    // OpenFHE stores (b,a), with decryption b+a*s. Read only the bottom limb;
    // CRT interpolation and signed big-integer conversion are unnecessary.
    auto b = input->GetElements()[0].GetElementAtIndex(0);
    auto a = input->GetElements()[1].GetElementAtIndex(0);
    b.SetFormat(Format::COEFFICIENT); a.SetFormat(Format::COEFFICIENT);
    const long double pi = std::acos(-1.0L);
    const long double q0 = b.GetModulus().ConvertToInt();
    if (std::abs(input->GetScalingFactor()/(static_cast<double>(q0)/key.gamma)-1.0) > 1e-12)
        throw std::invalid_argument("input scale must be q0/gamma for coefficient-domain messages");
    const uint32_t slots = key.slots;
    auto omega = [&](const NativeInteger& v) {
        long double angle = 2*pi*static_cast<long double>(v.ConvertToInt())/q0;
        return C(static_cast<double>(std::cos(angle)), static_cast<double>(std::sin(angle)));
    };
    std::array<std::vector<C>, 4> phases;
    for (auto& v : phases) v.resize(slots);
    std::vector<C> initial(slots);
    for (uint32_t i = 0; i < slots; ++i) {
        phases[0][i] = omega(a[i]); phases[1][i] = std::conj(phases[0][i]);
        phases[2][i] = omega(a[i+slots]); phases[3][i] = std::conj(phases[2][i]);
        // -i is essential: adding the conjugate recovers sine, not cosine.
        initial[i] = C(0, -key.gamma/(4*static_cast<double>(pi))) * omega(b[i]);
    }
    std::vector<CT> factors;
    // Every factor starts at the level reached after masking: 0 on the auxiliary path.
    factors.push_back(cc->Encrypt(key.outputPublicKey,
        cc->MakeCKKSPackedPlaintext(initial, 1, MaskingLevels(key.auxMasking), nullptr, slots)));
    std::array<Plaintext, 4> phasePlaintexts;
    std::array<DCRTPoly, 4> phaseQP;
    const auto params = CkksParams(cc);
    const auto paramsQP = params->GetParamsQP();
    const uint32_t bottom = params->GetElementParams()->GetParams().size()-1;
    for (size_t band = 0; band < 4; ++band) {
        if (key.auxMasking) {
            // Public Delta-scaled plaintext, lifted exactly from q0 to QP.
            auto pt = cc->MakeCKKSPackedPlaintext(phases[band], 1, bottom, nullptr, slots);
            if (pt->GetScalingFactor() != factors.front()->GetScalingFactor())
                throw std::logic_error("auxiliary masking scale mismatch");
            phaseQP[band] = LiftSmallToQP(pt->GetElement<DCRTPoly>(), paramsQP);
        }
        else {
            phasePlaintexts[band] = cc->MakeCKKSPackedPlaintext(phases[band], 1, 0, nullptr, slots);
        }
    }
    const auto scheme = cc->GetScheme();
    for (const auto& f : key.factors) {
        CT selected;
        if (key.auxMasking) {
            // Paper Sec. 4.1/4.4: PCMult over PQ, then Rescale_P (OpenFHE ApproxModDown).
            // (P*mask) * (Delta*phase) / P keeps scale Delta and level 0: no Q-level is consumed.
            DCRTPoly b(paramsQP, Format::EVALUATION, true), a(paramsQP, Format::EVALUATION, true);
            for (size_t band = 0; band < 4; ++band) {
                b += f.auxSelectors[band][0]*phaseQP[band];
                a += f.auxSelectors[band][1]*phaseQP[band];
            }
            auto extended = factors.front()->CloneEmpty();
            extended->SetElements({std::move(b), std::move(a)});
            selected = scheme->KeySwitchDown(extended);
        }
        else {
            for (size_t band = 0; band < 4; ++band) {
                auto term = cc->EvalMult(f.selectors[band], phasePlaintexts[band]);
                selected = selected ? cc->EvalAdd(selected, term) : term;
            }
            cc->RescaleInPlace(selected);
        }
        factors.push_back(key.useFused ? FusedBlindRotate(cc, f.fusedRotation, selected) :
                                       BlindRotate(cc, f.rotation, selected));
    }
    auto root = ProductTree(cc, std::move(factors));
    auto conjugate = cc->EvalAutomorphism(root, 4 * slots - 1, *key.conjugation);
    return cc->EvalAdd(root, conjugate);
}
}  // namespace ship
