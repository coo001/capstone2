// Correctness prototype for coefficient -> slot SHIP half bootstrapping.
// Requires FIXEDMANUAL, COMPLEX, full packing and a one-limb sparse-key input.
// The fused path is OpenFHE-specific; paper equivalence/security remain unverified.
#pragma once
#include "fused-rotation.h"
#include <cmath>
#include <set>

namespace ship {
struct FactorKey {
    std::array<CT, 4> selectors;
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
};

// Only setup has access to sparse support/signs and the dense output secret.
inline HalfBootstrapKey MakeHalfBootstrapKey(
    const CC& cc, const KeyPair<DCRTPoly>& outputKey,
    const std::vector<std::pair<uint32_t, int>>& sparseSupport,
    const std::string& inputKeyTag, double gamma, bool useFused = true) {
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
        for (size_t band = 0; band < 4; ++band)
            factor.selectors[band] = cc->Encrypt(outputKey.publicKey,
                cc->MakeCKKSPackedPlaintext(masks[band], 1, 0, nullptr, slots));
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
    std::array<Plaintext, 4> phasePlaintexts;
    for (size_t band = 0; band < 4; ++band)
        phasePlaintexts[band] = cc->MakeCKKSPackedPlaintext(phases[band], 1, 0, nullptr, slots);
    std::vector<CT> factors;
    factors.push_back(cc->Encrypt(key.outputPublicKey,
        cc->MakeCKKSPackedPlaintext(initial, 1, 1, nullptr, slots)));
    for (const auto& f : key.factors) {
        CT selected;
        for (size_t band = 0; band < 4; ++band) {
            auto term = cc->EvalMult(f.selectors[band], phasePlaintexts[band]);
            selected = selected ? cc->EvalAdd(selected, term) : term;
        }
        cc->RescaleInPlace(selected);
        factors.push_back(key.useFused ? FusedBlindRotate(cc, f.fusedRotation, selected) :
                                       BlindRotate(cc, f.rotation, selected));
    }
    auto root = ProductTree(cc, std::move(factors));
    auto conjugate = cc->EvalAutomorphism(root, 4 * slots - 1, *key.conjugation);
    return cc->EvalAdd(root, conjugate);
}
}  // namespace ship
