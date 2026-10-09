// Fused conditional rotation with shared HYBRID decomposition and delayed ModDown.
// An OpenFHE-specific construction; not a claim of paper-identical HMuxRot.
#pragma once
#include "reference.h"

namespace ship {
struct FusedMuxKey {
    EvalKey<DCRTPoly> body;
    EvalKey<DCRTPoly> mask;
    uint32_t automorphism;
    std::vector<uint32_t> permutation;
};
using FusedBlindKey = std::vector<std::array<FusedMuxKey, 2>>;

inline FusedMuxKey MakeFusedMuxKey(const CC& cc, const SK& sk, uint32_t beta, int32_t rotation) {
    if (!cc || !sk || sk->GetCryptoContext() != cc || beta > 1)
        throw std::invalid_argument("invalid fused conditional-rotation key");
    auto crypto = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc->GetCryptoParameters());
    if (!crypto || crypto->GetKeySwitchTechnique() != HYBRID)
        throw std::invalid_argument("fused conditional rotation requires HYBRID CKKS");
    const auto& secret = sk->GetPrivateElement();
    const auto params = secret.GetParams();
    const uint32_t n = secret.GetRingDimension(), m = 2*n;
    const uint32_t index = FindAutomorphismIndex2nComplex(rotation,m);
    std::vector<uint32_t> permutation(n), inversePermutation(n);
    PrecomputeAutoMap(n,index,&permutation);
    const uint32_t inverse = NativeInteger(index).ModInverse(m).ConvertToInt();
    PrecomputeAutoMap(n,inverse,&inversePermutation);

    // Key switch BEFORE applying sigma. The temporary destination is sigma^-1(s),
    // so sigma maps the resulting ciphertext back to the original output secret s.
    auto destination = std::make_shared<PrivateKeyImpl<DCRTPoly>>(cc);
    destination->SetPrivateElement(secret.AutomorphismTransform(inverse,inversePermutation));
    DCRTPoly constant(params,Format::COEFFICIENT,true);
    for (size_t i = 0; i < constant.GetNumOfElements(); ++i) {
        auto limb = constant.GetElementAtIndex(i);
        limb[0] = NativeInteger(beta);
        constant.SetElementAtIndex(i,std::move(limb));
    }
    constant.SetFormat(Format::EVALUATION);
    auto bodySecret = std::make_shared<PrivateKeyImpl<DCRTPoly>>(cc);
    bodySecret->SetPrivateElement(std::move(constant));
    auto maskSecret = std::make_shared<PrivateKeyImpl<DCRTPoly>>(cc);
    maskSecret->SetPrivateElement(beta ? secret : DCRTPoly(params,Format::EVALUATION,true));
    auto body = cc->KeySwitchGen(bodySecret,destination);
    auto mask = cc->KeySwitchGen(maskSecret,destination);
    // These tags identify the FINAL key after sigma, not the temporary destination.
    body->SetKeyTag(sk->GetKeyTag());
    mask->SetKeyTag(sk->GetKeyTag());
    return {std::move(body),std::move(mask),index,std::move(permutation)};
}

inline void ValidateFusedInput(const CC& cc, const CT& input, const FusedMuxKey& key) {
    if (!cc || !input || input->GetCryptoContext() != cc || input->GetElements().size() != 2 ||
        !key.body || !key.mask || key.body->GetCryptoContext() != cc || key.mask->GetCryptoContext() != cc ||
        input->GetKeyTag() != key.body->GetKeyTag() || input->GetKeyTag() != key.mask->GetKeyTag() ||
        key.permutation.size() != input->GetElements()[0].GetRingDimension())
        throw std::invalid_argument("fused conditional-rotation input/key mismatch");
}

// Internal PQ result. Both contributions must be accumulated before ModDown.
inline std::vector<DCRTPoly> FusedMuxExt(
    const CC& cc, const FusedMuxKey& key,
    const std::shared_ptr<std::vector<DCRTPoly>>& bodyDigits,
    const std::shared_ptr<std::vector<DCRTPoly>>& maskDigits,
    const std::shared_ptr<DCRTPoly::Params>& paramsQl) {
    const auto scheme = cc->GetScheme();
    auto body = scheme->EvalFastKeySwitchCoreExt(bodyDigits,key.body,paramsQl);
    auto mask = scheme->EvalFastKeySwitchCoreExt(maskDigits,key.mask,paramsQl);
    std::vector<DCRTPoly> result{(*body)[0]+(*mask)[0],(*body)[1]+(*mask)[1]};
    if (key.automorphism != 1)
        for (auto& element : result)
            element = element.AutomorphismTransform(key.automorphism,key.permutation);
    return result;
}

inline CT FusedMuxRotate(const CC& cc, const FusedMuxKey& key, const CT& input) {
    ValidateFusedInput(cc,input,key);
    const auto scheme = cc->GetScheme();
    auto body = scheme->EvalKeySwitchPrecomputeCore(input->GetElements()[0],input->GetCryptoParameters());
    auto mask = scheme->EvalKeySwitchPrecomputeCore(input->GetElements()[1],input->GetCryptoParameters());
    auto extended = input->Clone();
    extended->SetElements(FusedMuxExt(cc,key,body,mask,input->GetElements()[0].GetParams()));
    return scheme->KeySwitchDown(extended);
}

inline FusedBlindKey MakeFusedBlindKey(const CC& cc, const SK& sk, uint32_t shift,
                                      uint32_t slots, int direction = 1) {
    if (!slots || (slots & (slots-1)) || shift >= slots || (direction != 1 && direction != -1))
        throw std::invalid_argument("invalid fused blind-rotation dimensions");
    FusedBlindKey keys;
    for (uint32_t step = 1; step < slots; step <<= 1) {
        const uint32_t bit = (shift & step) != 0;
        keys.push_back({MakeFusedMuxKey(cc,sk,bit,direction*static_cast<int32_t>(step)),
                        MakeFusedMuxKey(cc,sk,1-bit,0)});
    }
    return keys;
}

inline CT FusedBlindRotate(const CC& cc, const FusedBlindKey& keys, const CT& input) {
    if (!input) throw std::invalid_argument("null fused blind-rotation input");
    auto output = input->Clone();
    const auto scheme = cc->GetScheme();
    for (const auto& level : keys) {
        ValidateFusedInput(cc,output,level[0]);
        ValidateFusedInput(cc,output,level[1]);
        // Both branches consume exactly the same two decompositions (hoisting).
        auto body = scheme->EvalKeySwitchPrecomputeCore(output->GetElements()[0],output->GetCryptoParameters());
        auto mask = scheme->EvalKeySwitchPrecomputeCore(output->GetElements()[1],output->GetCryptoParameters());
        const auto paramsQl = output->GetElements()[0].GetParams();
        auto rotated = FusedMuxExt(cc,level[0],body,mask,paramsQl);
        auto unchanged = FusedMuxExt(cc,level[1],body,mask,paramsQl);
        output->SetElements({rotated[0]+unchanged[0],rotated[1]+unchanged[1]});
        // One PQ -> Q reduction for the entire bit step; no manual 1/P factor.
        output = scheme->KeySwitchDown(output);
    }
    return output;
}
}  // namespace ship
