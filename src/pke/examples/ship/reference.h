// Experimental SHIP building blocks. See research/ship/README.md.
#pragma once
#include "openfhe.h"
#include <algorithm>
#include <array>
#include <complex>
#include <stdexcept>
#include <vector>

namespace ship {
using namespace lbcrypto;
using CT = Ciphertext<DCRTPoly>;
using CC = CryptoContext<DCRTPoly>;
using SK = PrivateKey<DCRTPoly>;
using C = std::complex<double>;

// A correctness reference, not the paper's fused/hoisted HMuxRot.
// Each key encrypts a gadget multiple of beta or beta*s under s.
// The evaluator receives neither beta nor the secret key.
struct MuxKey {
    EvalKey<DCRTPoly> body;
    EvalKey<DCRTPoly> mask;
    int32_t rotation;
};

inline MuxKey MakeMuxKey(const CC& cc, const SK& sk, uint32_t beta, int32_t rotation) {
    if (beta > 1)
        throw std::invalid_argument("beta must be a bit");
    auto params = sk->GetPrivateElement().GetParams();
    DCRTPoly constant(params, Format::COEFFICIENT, true);
    for (size_t i = 0; i < constant.GetNumOfElements(); ++i) {
        auto limb = constant.GetElementAtIndex(i);
        limb[0] = NativeInteger(beta);
        constant.SetElementAtIndex(i, std::move(limb));
    }
    constant.SetFormat(Format::EVALUATION);
    auto bodySecret = std::make_shared<PrivateKeyImpl<DCRTPoly>>(cc);
    bodySecret->SetPrivateElement(std::move(constant));
    auto maskSecret = std::make_shared<PrivateKeyImpl<DCRTPoly>>(cc);
    maskSecret->SetPrivateElement(beta ? sk->GetPrivateElement() :
        DCRTPoly(params, Format::EVALUATION, true));
    return {cc->KeySwitchGen(bodySecret, sk), cc->KeySwitchGen(maskSecret, sk), rotation};
}

inline CT MuxRotate(const CC& cc, const MuxKey& key, const CT& input) {
    if (!input || input->GetElements().size() != 2 || !key.body || !key.mask)
        throw std::invalid_argument("MuxRotate requires a degree-one ciphertext and two keys");
    if (input->GetKeyTag() != key.mask->GetKeyTag() ||
        input->GetKeyTag() != key.body->GetKeyTag())
        throw std::invalid_argument("MuxRotate key tag mismatch");
    if (input->GetCryptoContext() != cc || key.body->GetCryptoContext() != cc ||
        key.mask->GetCryptoContext() != cc)
        throw std::invalid_argument("MuxRotate context mismatch");
    auto scheme = cc->GetScheme();
    auto body = scheme->KeySwitchCore(input->GetElements()[0], key.body);
    auto mask = scheme->KeySwitchCore(input->GetElements()[1], key.mask);
    auto output = input->Clone();
    output->SetElements({(*body)[0] + (*mask)[0], (*body)[1] + (*mask)[1]});
    // KeySwitchCore already performs PQ -> Q, including division by P.
    // No extra multiplication by 1/P or manual scale/level changes belong here.
    return key.rotation ? cc->EvalRotate(output, key.rotation) : output;
}

using BlindKey = std::vector<std::array<MuxKey, 2>>;

inline BlindKey MakeBlindKey(const CC& cc, const SK& sk, uint32_t shift,
                            uint32_t slots, int direction = 1) {
    if (!slots || (slots & (slots - 1)) || shift >= slots ||
        (direction != 1 && direction != -1))
        throw std::invalid_argument("invalid blind rotation dimensions");
    BlindKey keys;
    for (uint32_t step = 1; step < slots; step <<= 1) {
        uint32_t bit = (shift & step) != 0;
        keys.push_back({MakeMuxKey(cc, sk, bit, direction * static_cast<int32_t>(step)),
                        MakeMuxKey(cc, sk, 1 - bit, 0)});
    }
    return keys;
}

inline CT BlindRotate(const CC& cc, const BlindKey& keys, const CT& input) {
    if (!input) throw std::invalid_argument("null blind-rotation input");
    auto output = input->Clone();
    for (const auto& level : keys)
        output = cc->EvalAdd(MuxRotate(cc, level[0], output), MuxRotate(cc, level[1], output));
    return output;
}

// Pre-rotation masks for the factor exp(2*pi*i*(a*s_j*X^j)_i/q).
// Input order: first half, conjugate first half, second half, conjugate second half.
inline std::array<std::vector<double>, 4> Masks(uint32_t j, int sign, uint32_t slots) {
    if (!slots || j >= 2 * slots || (sign != 1 && sign != -1))
        throw std::invalid_argument("invalid sparse support");
    std::array<std::vector<double>, 4> result;
    for (auto& v : result) v.resize(slots);
    // Derive the mask directly from X^j in Z[X]/(X^N+1).
    // After a right rotation by j % slots, output i reads input k=(i-j)%slots.
    for (uint32_t i = 0; i < slots; ++i) {
        const uint32_t n = 2 * slots;
        const uint32_t source = (i + n - j) % n;
        const int effectiveSign = i < j ? -sign : sign;
        const uint32_t band = (source >= slots ? 2 : 0) + (effectiveSign < 0 ? 1 : 0);
        result[band][source % slots] = 1.0;
    }
    return result;
}

// FIXEDMANUAL only: rescale each product and align odd leaves without changing scale.
inline CT ProductTree(const CC& cc, std::vector<CT> terms) {
    if (terms.empty()) throw std::invalid_argument("empty product tree");
    while (terms.size() > 1) {
        std::vector<CT> next;
        for (size_t i = 0; i < terms.size(); i += 2) {
            if (i + 1 == terms.size()) {
                next.push_back(terms[i]);
                continue;
            }
            auto a = terms[i]->Clone();
            auto b = terms[i+1]->Clone();
            const auto target = std::max(a->GetLevel(), b->GetLevel());
            if (a->GetLevel() < target) cc->LevelReduceInPlace(a, nullptr, target - a->GetLevel());
            if (b->GetLevel() < target) cc->LevelReduceInPlace(b, nullptr, target - b->GetLevel());
            auto product = cc->EvalMult(a, b);
            cc->RescaleInPlace(product);
            next.push_back(product);
        }
        terms = std::move(next);
    }
    return terms.front();
}
}  // namespace ship
