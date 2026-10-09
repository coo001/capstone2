// Correctness reference for a full-packing, FIXEDMANUAL CKKS slot round trip.
// Research prototype: no security estimate and no optimized SHIP latency claim.
#pragma once
#include "half-bootstrap.h"

namespace ship {

// Publish only bottom-modulus encryptions under the sparse key. Generating a
// full-chain HYBRID key under that sparse key would expose a different problem.
struct BottomSwitchKey {
    std::vector<std::array<NativePoly, 2>> digits;  // (b,a) encrypting B^d*s_dense
    std::string inputTag;
    std::string outputTag;
    uint32_t digitBits;
    CC context;
};

inline BottomSwitchKey MakeBottomSwitchKey(const CC& cc, const SK& dense, const SK& sparse,
                                          uint32_t digitBits = 8) {
    if (!dense || !sparse || dense->GetCryptoContext() != cc || sparse->GetCryptoContext() != cc ||
        digitBits == 0 || digitBits > 16)
        throw std::invalid_argument("invalid bottom-switch setup");
    auto oldSecret = dense->GetPrivateElement().GetElementAtIndex(0);
    auto newSecret = sparse->GetPrivateElement().GetElementAtIndex(0);
    if (*oldSecret.GetParams() != *newSecret.GetParams())
        throw std::invalid_argument("bottom moduli differ");
    oldSecret.SetFormat(Format::EVALUATION);
    newSecret.SetFormat(Format::EVALUATION);
    const auto params = oldSecret.GetParams();
    const auto q = oldSecret.GetModulus();
    const uint32_t count = (q.GetMSB() + digitBits - 1) / digitBits;
    auto crypto = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc->GetCryptoParameters());
    if (!crypto || crypto->GetNoiseScale() != 1)
        throw std::invalid_argument("bottom-switch reference requires CKKS noise scale 1");
    BottomSwitchKey result{{}, dense->GetKeyTag(), sparse->GetKeyTag(), digitBits, cc};
    NativePoly::DugType uniform;
    NativeInteger power(1);
    for (uint32_t d = 0; d < count; ++d) {
        NativePoly a(uniform, params, Format::EVALUATION);
        NativePoly e(crypto->GetDiscreteGaussianGenerator(), params, Format::EVALUATION);
        NativePoly b = e - a * newSecret + oldSecret * power;
        result.digits.push_back({std::move(b), std::move(a)});
        power = power.ModMul(NativeInteger(uint64_t{1} << digitBits), q);
    }
    return result;
}

inline CT ApplyBottomSwitch(const CT& input, const BottomSwitchKey& key) {
    if (!input || input->GetCryptoContext() != key.context || input->GetKeyTag() != key.inputTag ||
        input->GetElements().size() != 2 || key.digits.empty() ||
        input->GetElements()[0].GetNumOfElements() != 1 ||
        input->GetElements()[1].GetNumOfElements() != 1)
        throw std::invalid_argument("bottom-switch input must be a one-limb ciphertext under the dense key");
    auto b = input->GetElements()[0].GetElementAtIndex(0);
    auto a = input->GetElements()[1].GetElementAtIndex(0);
    if (*b.GetParams() != *key.digits[0][0].GetParams())
        throw std::invalid_argument("bottom-switch modulus mismatch");
    b.SetFormat(Format::EVALUATION);
    a.SetFormat(Format::COEFFICIENT);
    NativePoly outA(a.GetParams(), Format::EVALUATION, true);
    const uint64_t mask = (uint64_t{1} << key.digitBits) - 1;
    for (size_t d = 0; d < key.digits.size(); ++d) {
        NativePoly digit(a.GetParams(), Format::COEFFICIENT, true);
        const uint32_t shift = d * key.digitBits;
        for (size_t i = 0; i < a.GetLength(); ++i)
            digit[i] = NativeInteger((a[i].ConvertToInt() >> shift) & mask);
        digit.SetFormat(Format::EVALUATION);
        b += digit * key.digits[d][0];
        outA += digit * key.digits[d][1];
    }
    auto elements = input->GetElements();
    elements[0].SetElementAtIndex(0, std::move(b));
    elements[1].SetElementAtIndex(0, std::move(outA));
    auto result = input->Clone();
    result->SetElements(std::move(elements));
    result->SetKeyTag(key.outputTag);
    return result;
}

// Public multiplication by X^(N/2) or X^(-N/2), without a scale change.
// The positive direction multiplies every CKKS slot by i since 5^k == 1 mod 4.
// The negative direction moves the upper coefficient half to the lower half.
inline CT MultiplyByHalfMonomial(const CT& input, bool inverse) {
    if (!input) throw std::invalid_argument("null monomial input");
    auto elements = input->GetElements();
    for (auto& poly : elements) {
        const auto originalFormat = poly.GetFormat();
        poly.SetFormat(Format::COEFFICIENT);
        for (size_t k = 0; k < poly.GetNumOfElements(); ++k) {
            auto limb = poly.GetElementAtIndex(k);
            auto old = limb;
            const auto q = limb.GetModulus();
            const size_t half = limb.GetLength()/2;
            auto negate = [&](const NativeInteger& value) { return value == NativeInteger(0) ? value : q-value; };
            for (size_t i = 0; i < half; ++i) {
                limb[i] = inverse ? old[i+half] : negate(old[i+half]);
                limb[i+half] = inverse ? negate(old[i]) : old[i];
            }
            poly.SetElementAtIndex(k, std::move(limb));
        }
        poly.SetFormat(originalFormat);
    }
    auto result = input->Clone();
    result->SetElements(std::move(elements));
    return result;
}

struct PackingPlan {
    uint32_t slots;
    uint32_t babyStep;
    uint32_t inputLevel;
    std::vector<Plaintext> diagonals;
    bool hoist;
};

// V[k,j] = exp(2*pi*i*(5^k)*j/(2N)). This is OpenFHE's FFTSpecial,
// applied homomorphically to the slots holding coefficient pairs x_j+i*y_j.
inline PackingPlan MakePackingPlan(const CC& cc, const SK& dense, uint32_t level, bool hoist = true) {
    const uint32_t slots = cc->GetRingDimension()/2;
    const uint32_t cycl = 4*slots;
    uint32_t baby = 1;
    while (baby*baby < slots) baby <<= 1;
    PackingPlan result{slots,baby,level,{},hoist};
    std::vector<int32_t> rotations;
    for (uint32_t j = 1; j < baby; ++j) rotations.push_back(j);
    for (uint32_t g = baby; g < slots; g += baby) rotations.push_back(g);
    cc->EvalRotateKeyGen(dense,rotations);
    std::vector<uint32_t> exponents(slots);
    uint64_t e = 1;
    for (auto& x : exponents) { x = e; e = (5*e)%cycl; }
    const double tau = 2*std::acos(-1.0);
    for (uint32_t r = 0; r < slots; ++r) {
        const uint32_t giant = (r/baby)*baby;
        std::vector<C> d(slots);
        for (uint32_t k = 0; k < slots; ++k) {
            const uint32_t row = (k+slots-giant)%slots;
            const uint32_t column = (row+r)%slots;
            const uint32_t exponent = (uint64_t(exponents[row])*column)%cycl;
            d[k] = std::polar(1.0,tau*exponent/cycl);
        }
        result.diagonals.push_back(cc->MakeCKKSPackedPlaintext(d,1,level,nullptr,slots));
    }
    return result;
}

inline CT ApplyPacking(const CC& cc, const CT& input, const PackingPlan& plan) {
    if (!input || input->GetLevel() != plan.inputLevel || input->GetSlots() != plan.slots)
        throw std::invalid_argument("packing level/slot mismatch");
    std::vector<CT> baby(plan.babyStep);
    baby[0] = input;
    auto digits = plan.hoist ? cc->EvalFastRotationPrecompute(input) : nullptr;
    for (uint32_t j = 1; j < plan.babyStep; ++j)
        baby[j] = plan.hoist ? cc->EvalFastRotation(input,j,4*plan.slots,digits) : cc->EvalRotate(input,j);
    CT result;
    for (uint32_t g = 0; g < plan.slots; g += plan.babyStep) {
        CT group;
        for (uint32_t j = 0; j < plan.babyStep && g+j < plan.slots; ++j) {
            auto term = cc->EvalMult(baby[j],plan.diagonals[g+j]);
            group = group ? cc->EvalAdd(group,term) : term;
        }
        cc->RescaleInPlace(group);
        if (g) group = cc->EvalRotate(group,g);
        result = result ? cc->EvalAdd(result,group) : group;
    }
    return result;
}

struct FullBootstrapKey {
    BottomSwitchKey encapsulation;
    HalfBootstrapKey half;
    PackingPlan packing;
};

inline FullBootstrapKey MakeFullBootstrapKey(const CC& cc, const KeyPair<DCRTPoly>& dense,
                                           const SK& sparse, bool useFused = true) {
    if (!cc || !dense.publicKey || !dense.secretKey || !sparse ||
        dense.publicKey->GetCryptoContext() != cc || dense.secretKey->GetCryptoContext() != cc ||
        sparse->GetCryptoContext() != cc || dense.publicKey->GetKeyTag() != dense.secretKey->GetKeyTag())
        throw std::invalid_argument("invalid bootstrap key context or key pair");
    auto crypto = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc->GetCryptoParameters());
    if (!crypto || crypto->GetScalingTechnique() != FIXEDMANUAL ||
        crypto->GetKeySwitchTechnique() != HYBRID || crypto->GetCKKSDataType() != COMPLEX)
        throw std::invalid_argument("full bootstrap reference requires FIXEDMANUAL/HYBRID/COMPLEX CKKS");
    auto coefficients = sparse->GetPrivateElement().GetElementAtIndex(0);
    coefficients.SetFormat(Format::COEFFICIENT);
    const auto q = coefficients.GetModulus();
    std::vector<std::pair<uint32_t,int>> support;
    for (uint32_t i = 0; i < coefficients.GetLength(); ++i) {
        if (coefficients[i] == NativeInteger(0)) continue;
        if (coefficients[i] == NativeInteger(1)) support.emplace_back(i,1);
        else if (coefficients[i] == q-NativeInteger(1)) support.emplace_back(i,-1);
        else throw std::invalid_argument("sparse secret must be ternary");
    }
    uint32_t treeDepth = 0;
    for (size_t width = 1; width < support.size()+1; width <<= 1) ++treeDepth;
    // One masking level, then product tree, then one homomorphic FFT level.
    if (dense.secretKey->GetPrivateElement().GetNumOfElements() <= treeDepth+3)
        throw std::invalid_argument("insufficient output modulus budget");
    return {MakeBottomSwitchKey(cc,dense.secretKey,sparse),
            MakeHalfBootstrapKey(cc,dense,support,sparse->GetKeyTag(),1.0,useFused),
            MakePackingPlan(cc,dense.secretKey,1+treeDepth,useFused)};
}

// This uses the SHIP small-angle sine approximation; it is not an exact identity
// for arbitrary messages. Keep coefficient magnitudes small relative to gamma.
// No secret key is available to this evaluator.
inline CT FullBootstrap(const CC& cc, const CT& input, const FullBootstrapKey& key) {
    if (!cc || cc != key.encapsulation.context || !input ||
        input->GetSlots() != cc->GetRingDimension()/2 || input->GetNoiseScaleDeg() != 1 ||
        !std::isfinite(input->GetScalingFactor()) || input->GetScalingFactor() <= 0)
        throw std::invalid_argument("full bootstrap requires full packing and a rescaled input");
    auto sparseInput = ApplyBottomSwitch(input,key.encapsulation);
    auto half = key.half;
    const auto q0 = sparseInput->GetElements()[0].GetElementAtIndex(0).GetModulus().ConvertToDouble();
    half.gamma = q0/input->GetScalingFactor();
    auto first = HalfBootstrap(cc,sparseInput,half);
    auto second = HalfBootstrap(cc,MultiplyByHalfMonomial(sparseInput,true),half);
    auto paired = cc->EvalAdd(first,MultiplyByHalfMonomial(second,false));
    return ApplyPacking(cc,paired,key.packing);
}
}  // namespace ship
