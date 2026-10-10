//==================================================================================
// BSD 2-Clause License
//
// Copyright (c) 2014-2025, NJIT, Duality Technologies Inc. and other contributors
//
// All rights reserved.
//
// Author TPOC: contact@openfhe.org
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//==================================================================================

#include "scheme/ckksrns/ckksrns-ship.h"

#include "cryptocontext.h"
#include "gen-cryptocontext.h"
#include "math/distributiongenerator.h"
#include "math/nbtheory.h"
#include "scheme/ckksrns/ckksrns-cryptoparameters.h"
#include "scheme/ckksrns/gen-cryptocontext-ckksrns.h"
#include "utils/exception.h"

#include <array>
#include <cmath>
#include <complex>
#include <map>
#include <mutex>
#include <random>
#include <set>
#include <vector>

namespace lbcrypto {

namespace {

using CT     = Ciphertext<DCRTPoly>;
using CC     = CryptoContext<DCRTPoly>;
using C      = std::complex<double>;
using QPPair = std::array<DCRTPoly, 2>;  // (b, a) over the extended basis QP

std::shared_ptr<CryptoParametersCKKSRNS> CkksParams(const CC& cc) {
    auto params = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc->GetCryptoParameters());
    if (!params)
        OPENFHE_THROW("SHIP requires CKKS RNS parameters");
    return params;
}

uint32_t CeilLog(uint64_t value, uint64_t base) {
    uint32_t k = 0;
    for (uint64_t power = 1; power < value; power *= base)
        ++k;
    return k;
}

// Rot_t(v)_i = v_{i - t mod S}: right rotation, the paper's convention (Section 1.2, Notations).
std::vector<double> RotateRight(const std::vector<double>& v, uint32_t t) {
    const uint32_t S = v.size();
    std::vector<double> out(S);
    for (uint32_t i = 0; i < S; ++i)
        out[i] = v[(i + S - (t % S)) % S];
    return out;
}

// Pre-rotation masks for the factor exp(2 pi i (a * s_j X^j)_i / q0), Lemma 1 rewritten as Equation (5).
// Bands: omega^{a} first half, omega^{-a} first half, omega^{a} second half, omega^{-a} second half.
std::array<std::vector<double>, 4> PreRotationMasks(uint32_t j, int sign, uint32_t S) {
    std::array<std::vector<double>, 4> result;
    for (auto& v : result)
        v.assign(S, 0.0);
    const uint32_t n = 2 * S;
    for (uint32_t i = 0; i < S; ++i) {
        const uint32_t source      = (i + n - j) % n;
        const int effectiveSign    = i < j ? -sign : sign;
        const uint32_t band        = (source >= S ? 2 : 0) + (effectiveSign < 0 ? 1 : 0);
        result[band][source % S] = 1.0;
    }
    return result;
}

// Signed lift of the bottom limb of a small integer polynomial to all limbs of QP (exact below q0/4).
DCRTPoly LiftSmallToQP(const DCRTPoly& poly, const std::shared_ptr<DCRTPoly::Params>& paramsQP) {
    auto limb = poly.GetElementAtIndex(0);
    limb.SetFormat(Format::COEFFICIENT);
    const auto q0      = limb.GetModulus();
    const auto quarter = q0 >> 2;
    for (size_t i = 0; i < limb.GetLength(); ++i)
        if (limb[i] >= quarter && limb[i] <= q0 - quarter)
            OPENFHE_THROW("SHIP: plaintext coefficient too large for an exact QP lift");
    DCRTPoly result(paramsQP, Format::COEFFICIENT, true);
    for (size_t k = 0; k < result.GetNumOfElements(); ++k) {
        auto copy          = limb;
        const auto& target = paramsQP->GetParams()[k];
        if (target->GetModulus() != q0)
            copy.SwitchModulus(target->GetModulus(), target->GetRootOfUnity(), 0, 0);
        result.SetElementAtIndex(k, std::move(copy));
    }
    result.SetFormat(Format::EVALUATION);
    return result;
}

// The secret extended from Q to QP, as in KeySwitchHYBRID::KeySwitchGenInternal.
DCRTPoly SecretOverQP(const PrivateKey<DCRTPoly>& sk, const std::shared_ptr<DCRTPoly::Params>& paramsQP) {
    auto secret = sk->GetPrivateElement().Clone();
    secret.SetFormat(Format::COEFFICIENT);
    DCRTPoly result(paramsQP, Format::COEFFICIENT, true);
    const size_t sizeQ = secret.GetNumOfElements();
    for (size_t i = 0; i < sizeQ; ++i)
        result.SetElementAtIndex(i, secret.GetElementAtIndex(i));
    for (size_t j = sizeQ; j < result.GetNumOfElements(); ++j) {
        auto limb          = secret.GetElementAtIndex(0);
        const auto& target = paramsQP->GetParams()[j];
        limb.SwitchModulus(target->GetModulus(), target->GetRootOfUnity(), 0, 0);
        result.SetElementAtIndex(j, std::move(limb));
    }
    result.SetFormat(Format::EVALUATION);
    return result;
}

// Enc_QP(P * 2^extraBits * v): v slot-encoded with scaling factor P * 2^extraBits under the output secret.
// round(m_i * P * 2^extraBits / Delta) from the exact Delta-encoding m_i; relative error <= 2^-log2(Delta).
// After the product with a Delta-scaled phase and Rescale_P, the factor has scale 2^extraBits * Delta.
QPPair EncryptScaledByP(const CC& cc, const DCRTPoly& secretQP, const std::vector<double>& values, uint32_t S,
                        uint32_t extraBits) {
    const auto params   = CkksParams(cc);
    const auto paramsQP = params->GetParamsQP();
    const uint32_t bottom = params->GetElementParams()->GetParams().size() - 1;
    auto pt               = cc->MakeCKKSPackedPlaintext(values, 1, bottom, nullptr, S);
    int exponent          = 0;
    if (std::frexp(pt->GetScalingFactor(), &exponent) != 0.5 || exponent < 2)
        OPENFHE_THROW("SHIP requires a power-of-two CKKS scaling factor");
    const uint32_t deltaBits = exponent - 1;
    auto encoded             = pt->GetElement<DCRTPoly>();
    encoded.SetFormat(Format::COEFFICIENT);
    const auto limb       = encoded.GetElementAtIndex(0);
    const auto q0         = limb.GetModulus();
    const BigInteger P    = params->GetParamsP()->GetModulus();
    const BigInteger half = BigInteger(1).LShift(deltaBits - 1);
    DCRTPoly plain(paramsQP, Format::COEFFICIENT, true);
    std::vector<NativePoly> limbs;
    for (size_t k = 0; k < plain.GetNumOfElements(); ++k)
        limbs.push_back(plain.GetElementAtIndex(k));
    for (size_t i = 0; i < limb.GetLength(); ++i) {
        if (limb[i] == NativeInteger(0))
            continue;
        const bool negative       = limb[i] > (q0 >> 1);
        const uint64_t magnitude  = (negative ? q0 - limb[i] : limb[i]).ConvertToInt<uint64_t>();
        const BigInteger scaled   = ((BigInteger(magnitude) * P).LShift(extraBits) + half).RShift(deltaBits);
        for (auto& out : limbs) {
            const auto r = out.GetModulus();
            NativeInteger value(scaled.Mod(BigInteger(r.ConvertToInt<uint64_t>())).ConvertToInt<uint64_t>());
            out[i] = (negative && value != NativeInteger(0)) ? r - value : value;
        }
    }
    for (size_t k = 0; k < limbs.size(); ++k)
        plain.SetElementAtIndex(k, std::move(limbs[k]));
    plain.SetFormat(Format::EVALUATION);
    DCRTPoly::DugType uniform;
    DCRTPoly a(uniform, paramsQP, Format::EVALUATION);
    DCRTPoly e(params->GetDiscreteGaussianGenerator(), paramsQP, Format::EVALUATION);
    DCRTPoly b = e - a * secretQP + plain;
    return {std::move(b), std::move(a)};
}

}  // namespace

// ---------------------------------------------------------------------------------------------
// Key material
// ---------------------------------------------------------------------------------------------

// HMuxRot key (Definition 1 with gadget decomposition, Algorithm 5): the key switch targets
// sigma^{-1}(s) and sigma is applied afterwards, so the decomposition is shared (hoisted) across
// the B branches of one B-to-1 mux-rotate and ModDown is applied once per mux step.
struct SHIPMuxKey {
    EvalKey<DCRTPoly> body;
    EvalKey<DCRTPoly> mask;
    uint32_t automorphism = 1;
    std::vector<uint32_t> permutation;
};

struct SHIPFactorKey {
    uint32_t offset = 0;  // public rotation offset o (mod N/2)
    // column[i][band] = Enc_QP(P * 1_{i = r0} * Rot_{o + r0}(M'_band)), i < theta
    std::vector<std::array<QPPair, 4>> column;
    // mux[t][d] selects a right rotation by theta * d * B^t if the t-th base-B digit of r1 equals d
    std::vector<std::vector<SHIPMuxKey>> mux;
};

struct SHIPEncapsulationKey {
    std::shared_ptr<ILNativeParams> paramsQ0;
    std::shared_ptr<ILNativeParams> paramsP;
    NativePoly bQ0, aQ0, bP, aP;  // Enc_{s_sparse}(p' * s_dense) mod q0 * p', EVALUATION format
    NativeInteger pInvModQ0;
};

struct SHIPPackingPlan {
    uint32_t babyStep   = 0;
    uint32_t inputLevel = 0;
    std::vector<Plaintext> diagonals;
};

class SHIPBootstrapKey {
public:
    SHIPParams params;
    std::string denseTag;
    std::string sparseTag;
    uint32_t ringDim   = 0;
    uint32_t slots     = 0;
    uint32_t treeDepth = 0;
    uint32_t leafExtraBits = 0;  // product-tree leaves have scale 2^(p + leafExtraBits)
    SHIPEncapsulationKey encapsulation;
    std::vector<SHIPFactorKey> factors;
    std::map<uint32_t, std::vector<uint32_t>> columnPermutations;  // right rotation t -> automorphism map
    std::shared_ptr<std::map<uint32_t, EvalKey<DCRTPoly>>> conjugation;
    SHIPPackingPlan packing;
};

namespace {

SHIPMuxKey MakeMuxKey(const CC& cc, const PrivateKey<DCRTPoly>& sk, uint32_t beta, int32_t rotation) {
    const auto& secret  = sk->GetPrivateElement();
    const auto params   = secret.GetParams();
    const uint32_t n    = secret.GetRingDimension();
    const uint32_t m    = 2 * n;
    const uint32_t index = FindAutomorphismIndex2nComplex(rotation, m);
    SHIPMuxKey key;
    key.automorphism = index;
    std::vector<uint32_t> inversePermutation(n);
    key.permutation.resize(n);
    PrecomputeAutoMap(n, index, &key.permutation);
    const uint32_t inverse = NativeInteger(index).ModInverse(m).ConvertToInt();
    PrecomputeAutoMap(n, inverse, &inversePermutation);
    auto destination = std::make_shared<PrivateKeyImpl<DCRTPoly>>(cc);
    destination->SetPrivateElement(secret.AutomorphismTransform(inverse, inversePermutation));
    DCRTPoly constant(params, Format::COEFFICIENT, true);
    for (size_t i = 0; i < constant.GetNumOfElements(); ++i) {
        auto limb = constant.GetElementAtIndex(i);
        limb[0]   = NativeInteger(beta);
        constant.SetElementAtIndex(i, std::move(limb));
    }
    constant.SetFormat(Format::EVALUATION);
    auto bodySecret = std::make_shared<PrivateKeyImpl<DCRTPoly>>(cc);
    bodySecret->SetPrivateElement(std::move(constant));
    auto maskSecret = std::make_shared<PrivateKeyImpl<DCRTPoly>>(cc);
    maskSecret->SetPrivateElement(beta ? secret : DCRTPoly(params, Format::EVALUATION, true));
    key.body = cc->KeySwitchGen(bodySecret, destination);
    key.mask = cc->KeySwitchGen(maskSecret, destination);
    key.body->SetKeyTag(sk->GetKeyTag());  // tag of the final key (after sigma)
    key.mask->SetKeyTag(sk->GetKeyTag());
    return key;
}

// sum over the B branches of one mux step, sharing the decompositions; one ModDown.
CT MuxRotate(const CC& cc, const std::vector<std::vector<SHIPMuxKey>>& steps, CT input) {
    const auto scheme = cc->GetScheme();
    for (const auto& branches : steps) {
        auto body           = scheme->EvalKeySwitchPrecomputeCore(input->GetElements()[0], input->GetCryptoParameters());
        auto mask           = scheme->EvalKeySwitchPrecomputeCore(input->GetElements()[1], input->GetCryptoParameters());
        const auto paramsQl = input->GetElements()[0].GetParams();
        std::vector<DCRTPoly> sum;
        for (const auto& key : branches) {
            auto b = scheme->EvalFastKeySwitchCoreExt(body, key.body, paramsQl);
            auto a = scheme->EvalFastKeySwitchCoreExt(mask, key.mask, paramsQl);
            std::vector<DCRTPoly> term{(*b)[0] + (*a)[0], (*b)[1] + (*a)[1]};
            if (key.automorphism != 1)
                for (auto& element : term)
                    element = element.AutomorphismTransform(key.automorphism, key.permutation);
            if (sum.empty()) {
                sum = std::move(term);
            }
            else {
                sum[0] += term[0];
                sum[1] += term[1];
            }
        }
        auto extended = input->CloneEmpty();
        extended->SetElements(std::move(sum));
        input = scheme->KeySwitchDown(extended);
    }
    return input;
}

// FIXEDMANUAL assumes every rescaling prime equals Delta. With small NTT-friendly primes the
// deviation |q/Delta - 1| reaches a few percent and compounds through the product tree. All leaves
// enter at level 0 with scale exactly Delta, so the true root scale is deterministic: simulate the
// pairing of ProductTree and return (true root scale) / Delta.
long double ProductTreeScaleRatio(size_t leaves, long double leafScale, long double delta,
                                  const std::vector<long double>& moduli) {
    std::vector<std::pair<long double, uint32_t>> nodes(leaves, {leafScale, 0});  // (scale, level)
    const size_t sizeQ = moduli.size();
    while (nodes.size() > 1) {
        std::vector<std::pair<long double, uint32_t>> next;
        for (size_t i = 0; i + 1 < nodes.size(); i += 2) {
            const uint32_t level = std::max(nodes[i].second, nodes[i + 1].second);
            next.push_back({nodes[i].first * nodes[i + 1].first / moduli[sizeQ - 1 - level], level + 1});
        }
        if (nodes.size() % 2)
            next.push_back(nodes.back());
        nodes = std::move(next);
    }
    return nodes.front().first / delta;
}

CT ProductTree(const CC& cc, std::vector<CT> terms) {
    while (terms.size() > 1) {
        std::vector<CT> next((terms.size() + 1) / 2);
#pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < terms.size() / 2; ++i) {
            auto a             = terms[2 * i]->Clone();
            auto b             = terms[2 * i + 1]->Clone();
            const auto target  = std::max(a->GetLevel(), b->GetLevel());
            if (a->GetLevel() < target)
                cc->LevelReduceInPlace(a, nullptr, target - a->GetLevel());
            if (b->GetLevel() < target)
                cc->LevelReduceInPlace(b, nullptr, target - b->GetLevel());
            auto product = cc->EvalMult(a, b);
            cc->RescaleInPlace(product);
            next[i] = product;
        }
        if (terms.size() % 2)
            next.back() = terms.back();
        terms = std::move(next);
    }
    return terms.front();
}

// Public multiplication by X^(N/2) (inverse = false: every slot times i) or X^(-N/2).
CT MultiplyByHalfMonomial(const CT& input, bool inverse) {
    auto elements = input->GetElements();
    for (auto& poly : elements) {
        const auto originalFormat = poly.GetFormat();
        poly.SetFormat(Format::COEFFICIENT);
        for (size_t k = 0; k < poly.GetNumOfElements(); ++k) {
            auto limb         = poly.GetElementAtIndex(k);
            const auto old    = limb;
            const auto q      = limb.GetModulus();
            const size_t half = limb.GetLength() / 2;
            auto negate       = [&](const NativeInteger& v) { return v == NativeInteger(0) ? v : q - v; };
            for (size_t i = 0; i < half; ++i) {
                limb[i]        = inverse ? old[i + half] : negate(old[i + half]);
                limb[i + half] = inverse ? negate(old[i]) : old[i];
            }
            poly.SetElementAtIndex(k, std::move(limb));
        }
        poly.SetFormat(originalFormat);
    }
    auto result = input->Clone();
    result->SetElements(std::move(elements));
    return result;
}

// S2C: slots mu -> coefficients (Re mu | Im mu), i.e. multiplication of the slot vector by
// V[k][j] = exp(2 pi i 5^k j / (2N)) (paper Section 3.2, dense matrix with BSGS, one level).
SHIPPackingPlan MakePackingPlan(const CC& cc, const PrivateKey<DCRTPoly>& sk, uint32_t level) {
    const uint32_t S    = cc->GetRingDimension() / 2;
    const uint32_t cycl = 4 * S;
    uint32_t baby       = 1;
    while (baby * baby < S)
        baby <<= 1;
    SHIPPackingPlan plan;
    plan.babyStep   = baby;
    plan.inputLevel = level;
    std::vector<int32_t> rotations;
    for (uint32_t j = 1; j < baby; ++j)
        rotations.push_back(j);
    for (uint32_t g = baby; g < S; g += baby)
        rotations.push_back(g);
    cc->EvalRotateKeyGen(sk, rotations);
    std::vector<uint32_t> exponents(S);
    uint64_t e = 1;
    for (auto& x : exponents) {
        x = e;
        e = (5 * e) % cycl;
    }
    const double tau = 2 * std::acos(-1.0);
    // The rescale after the diagonal products removes the top prime at this level; scale the
    // diagonals by q/Delta so that the output scale equals the input scale exactly.
    const auto params      = CkksParams(cc);
    const auto& moduli     = params->GetElementParams()->GetParams();
    const double rescaleBy = moduli[moduli.size() - 1 - level]->GetModulus().ConvertToDouble();
    const double ratio     = rescaleBy / params->GetScalingFactorReal(level);
    plan.diagonals.resize(S);
    for (uint32_t r = 0; r < S; ++r) {
        const uint32_t giant = (r / baby) * baby;
        std::vector<C> d(S);
        for (uint32_t k = 0; k < S; ++k) {
            const uint32_t row      = (k + S - giant) % S;
            const uint32_t column   = (row + r) % S;
            const uint32_t exponent = (uint64_t(exponents[row]) * column) % cycl;
            d[k]                    = std::polar(ratio, tau * exponent / cycl);
        }
        plan.diagonals[r] = cc->MakeCKKSPackedPlaintext(d, 1, level, nullptr, S);
    }
    return plan;
}

CT ApplyPacking(const CC& cc, const CT& input, const SHIPPackingPlan& plan, uint32_t S) {
    std::vector<CT> baby(plan.babyStep);
    baby[0]     = input;
    auto digits = cc->EvalFastRotationPrecompute(input);
    for (uint32_t j = 1; j < plan.babyStep; ++j)
        baby[j] = cc->EvalFastRotation(input, j, 4 * S, digits);
    const uint32_t groups = (S + plan.babyStep - 1) / plan.babyStep;
    std::vector<CT> partial(groups);
#pragma omp parallel for schedule(dynamic)
    for (uint32_t gi = 0; gi < groups; ++gi) {
        const uint32_t g = gi * plan.babyStep;
        CT group;
        for (uint32_t j = 0; j < plan.babyStep && g + j < S; ++j) {
            auto term = cc->EvalMult(baby[j], plan.diagonals[g + j]);
            group     = group ? cc->EvalAdd(group, term) : term;
        }
        cc->RescaleInPlace(group);
        partial[gi] = g ? cc->EvalRotate(group, g) : group;
    }
    CT result = partial[0];
    for (uint32_t gi = 1; gi < groups; ++gi)
        result = cc->EvalAdd(result, partial[gi]);
    return result;
}

NativePoly ModDownByP(NativePoly xQ0, NativePoly xP, const SHIPEncapsulationKey& key) {
    xP.SetFormat(Format::COEFFICIENT);
    xP.SwitchModulus(key.paramsQ0->GetModulus(), key.paramsQ0->GetRootOfUnity(), 0, 0);  // centered lift
    xP.SetFormat(Format::EVALUATION);
    return (xQ0 - xP) * key.pInvModQ0;
}

// Sparse secret encapsulation at modulus q0 * p' (paper Section 2.2 and 5.2, [BTH22]).
CT Encapsulate(const CT& input, const SHIPBootstrapKey& key) {
    const auto& k = key.encapsulation;
    auto b        = input->GetElements()[0].GetElementAtIndex(0);
    auto aQ0      = input->GetElements()[1].GetElementAtIndex(0);
    b.SetFormat(Format::EVALUATION);
    aQ0.SetFormat(Format::COEFFICIENT);
    auto aP = aQ0;
    aP.SwitchModulus(k.paramsP->GetModulus(), k.paramsP->GetRootOfUnity(), 0, 0);
    aQ0.SetFormat(Format::EVALUATION);
    aP.SetFormat(Format::EVALUATION);
    auto newB = b + ModDownByP(aQ0 * k.bQ0, aP * k.bP, k);
    auto newA = ModDownByP(aQ0 * k.aQ0, aP * k.aP, k);
    auto elements = input->GetElements();
    elements[0].SetElementAtIndex(0, std::move(newB));
    elements[1].SetElementAtIndex(0, std::move(newA));
    auto result = input->Clone();
    result->SetElements(std::move(elements));
    result->SetKeyTag(key.sparseTag);
    return result;
}

// Algorithm 1 with Algorithm 4 (column + mux) and masks over PQ (Sections 4.1, 4.4): no masking level.
CT HalfBootstrap(const CC& cc, const CT& input, const SHIPBootstrapKey& key, double gamma) {
    const auto params     = CkksParams(cc);
    const auto paramsQP   = params->GetParamsQP();
    const uint32_t S      = key.slots;
    const uint32_t bottom = params->GetElementParams()->GetParams().size() - 1;
    auto b                = input->GetElements()[0].GetElementAtIndex(0);
    auto a                = input->GetElements()[1].GetElementAtIndex(0);
    b.SetFormat(Format::COEFFICIENT);
    a.SetFormat(Format::COEFFICIENT);
    const long double pi = std::acos(-1.0L);
    const long double q0 = b.GetModulus().ConvertToInt();
    auto omega           = [&](const NativeInteger& v) {
        const long double angle = 2 * pi * static_cast<long double>(v.ConvertToInt()) / q0;
        return C(static_cast<double>(std::cos(angle)), static_cast<double>(std::sin(angle)));
    };
    std::array<std::vector<C>, 4> phases;
    for (auto& v : phases)
        v.resize(S);
    std::vector<C> initial(S);
    for (uint32_t i = 0; i < S; ++i) {
        phases[0][i] = omega(a[i]);
        phases[1][i] = std::conj(phases[0][i]);
        phases[2][i] = omega(a[i + S]);
        phases[3][i] = std::conj(phases[2][i]);
        initial[i]   = C(0, -gamma / (4 * static_cast<double>(pi))) * omega(b[i]);  // gamma / (4 i pi) * omega^b
    }
    // Exact compensation of the product-tree scale drift through the public factor pt_0.
    {
        std::vector<long double> moduli;
        for (const auto& q : params->GetElementParams()->GetParams())
            moduli.push_back(q->GetModulus().ConvertToDouble());
        const long double delta     = params->GetScalingFactorReal(0);
        const long double leafScale = std::ldexp(delta, key.leafExtraBits);
        const long double kappa     = ProductTreeScaleRatio(key.factors.size() + 1, leafScale, delta, moduli);
        // pt_0 enters with scale 2^extra * Delta like the other leaves; the root lands exactly at Delta.
        for (auto& v : initial)
            v *= static_cast<double>(std::ldexp(1.0L, key.leafExtraBits) / kappa);
    }
    // pt_0 as a trivial ciphertext at the top level (no public key needed).
    auto pt0     = cc->MakeCKKSPackedPlaintext(initial, 1, 0, nullptr, S);
    auto pt0Poly = pt0->GetElement<DCRTPoly>();
    pt0Poly.SetFormat(Format::EVALUATION);
    auto trivial = std::make_shared<CiphertextImpl<DCRTPoly>>(cc, key.denseTag, CKKS_PACKED_ENCODING);
    trivial->SetElements({pt0Poly, DCRTPoly(pt0Poly.GetParams(), Format::EVALUATION, true)});
    trivial->SetLevel(0);
    trivial->SetNoiseScaleDeg(1);
    trivial->SetScalingFactor(pt0->GetScalingFactor());
    trivial->SetSlots(S);
    std::array<DCRTPoly, 4> phaseQP;
    for (size_t band = 0; band < 4; ++band) {
        auto pt = cc->MakeCKKSPackedPlaintext(phases[band], 1, bottom, nullptr, S);
        if (pt->GetScalingFactor() != pt0->GetScalingFactor())
            OPENFHE_THROW("SHIP: scaling factor mismatch");
        phaseQP[band] = LiftSmallToQP(pt->GetElement<DCRTPoly>(), paramsQP);
    }
    const auto scheme = cc->GetScheme();
    const uint32_t m  = 2 * cc->GetRingDimension();
    std::vector<CT> factors(key.factors.size() + 1);
    factors[0] = trivial;
#pragma omp parallel for schedule(dynamic)
    for (size_t f = 0; f < key.factors.size(); ++f) {
        const auto& fk = key.factors[f];
        DCRTPoly sb(paramsQP, Format::EVALUATION, true), sa(paramsQP, Format::EVALUATION, true);
        for (uint32_t i = 0; i < fk.column.size(); ++i) {
            const uint32_t t     = (fk.offset + i) % S;
            const uint32_t index = FindAutomorphismIndex2nComplex(-static_cast<int32_t>(t), m);
            const auto& perm     = key.columnPermutations.at(t);
            for (size_t band = 0; band < 4; ++band) {
                const auto rotated = t ? phaseQP[band].AutomorphismTransform(index, perm) : phaseQP[band];
                sb += fk.column[i][band][0] * rotated;
                sa += fk.column[i][band][1] * rotated;
            }
        }
        auto extended = trivial->CloneEmpty();
        extended->SetElements({std::move(sb), std::move(sa)});
        auto selected = scheme->KeySwitchDown(extended);  // Rescale_P: level stays 0
        factors[f + 1] = MuxRotate(cc, fk.mux, selected);
    }
    auto root      = ProductTree(cc, std::move(factors));
    auto conjugate = cc->EvalAutomorphism(root, m - 1, *key.conjugation);
    return cc->EvalAdd(root, conjugate);
}

std::mutex& KeyStoreMutex() {
    static std::mutex mutex;
    return mutex;
}
std::map<std::string, std::shared_ptr<SHIPBootstrapKey>>& KeyStore() {
    static std::map<std::string, std::shared_ptr<SHIPBootstrapKey>> store;
    return store;
}

}  // namespace

uint32_t SHIPProductTreeDepth(uint32_t hammingWeight) {
    return CeilLog(uint64_t(hammingWeight) + 1, 2);
}

SHIPContextSpec SHIPContextSpec::LL13() {
    // Paper LL13: Base 25, S2C 20, Mult 20, Boot 25 x 5, Aux 28, dnum 8, log PQ = 218.
    // Here: Base 25, S2C/Mult 20, Boot 24 x 4 + 28 (output back at scale 2^20), Aux 28, dnum 8.
    SHIPContextSpec s;
    s.ringDim        = 1 << 13;
    s.firstModBits   = 25;
    s.scalingModBits = 20;
    s.bootModBits    = 24;
    s.multLevels     = 1;
    s.hammingWeight  = 31;
    s.numLargeDigits = 8;
    s.auxModBits     = 28;
    return s;
}

SHIPContextSpec SHIPContextSpec::LL14() {
    // Paper LL14: Base 48, S2C 37, Mult 37, Boot 43 x 5, Aux 50 x 2, dnum 4, log PQ = 437.
    // Here: Base 48, S2C/Mult 37, Boot 42 x 4 + 47 (output back at scale 2^37), Aux 50 x 2, dnum 4.
    SHIPContextSpec s;
    s.ringDim        = 1 << 14;
    s.firstModBits   = 48;
    s.scalingModBits = 37;
    s.bootModBits    = 42;
    s.multLevels     = 1;
    s.hammingWeight  = 31;
    s.numLargeDigits = 4;
    s.auxModBits     = 50;
    return s;
}

uint32_t SHIPLogQP(const CryptoContext<DCRTPoly>& cc) {
    return CkksParams(cc)->GetParamsQP()->GetModulus().GetMSB();
}

CryptoContext<DCRTPoly> GenSHIPCryptoContext(const SHIPContextSpec& spec) {
    if (spec.ringDim < 16 || (spec.ringDim & (spec.ringDim - 1)) || spec.firstModBits <= spec.scalingModBits ||
        spec.numLargeDigits == 0 || spec.auxModBits == 0)
        OPENFHE_THROW("Invalid SHIP context specification");
    const uint32_t depth = 1 + spec.multLevels + SHIPProductTreeDepth(spec.hammingWeight);
    CCParams<CryptoContextCKKSRNS> p;
    p.SetSecurityLevel(HEStd_NotSet);  // checked below with the actual auxiliary moduli
    p.SetRingDim(spec.ringDim);
    p.SetBatchSize(spec.ringDim / 2);
    p.SetMultiplicativeDepth(depth);
    p.SetScalingModSize(spec.scalingModBits);
    p.SetFirstModSize(spec.firstModBits);
    p.SetScalingTechnique(FIXEDMANUAL);
    p.SetKeySwitchTechnique(HYBRID);
    p.SetNumLargeDigits(spec.numLargeDigits);
    p.SetCKKSDataType(COMPLEX);
    p.SetSecretKeyDist(UNIFORM_TERNARY);
    auto cc     = GenCryptoContext(p);
    auto params = CkksParams(cc);
    if (spec.bootModBits > spec.scalingModBits) {
        // Role-specific primes (paper Table 2). Chain order from the bottom: q0, S2C, Mult..., then the
        // boot primes; the lowest boot prime is removed by the last product-tree rescale and has
        // 2b - p bits so that the bootstrapped output has scale exactly 2^p again.
        const uint64_t m = 2 * uint64_t(spec.ringDim);
        std::set<NativeInteger> used;
        auto closest = [&](uint32_t bits) {
            NativeInteger up   = FirstPrime<NativeInteger>(bits, m);
            NativeInteger down = LastPrime<NativeInteger>(bits, m);
            while (used.count(up))
                up = NextPrime<NativeInteger>(up, m);
            while (used.count(down))
                down = PreviousPrime<NativeInteger>(down, m);
            const NativeInteger target(uint64_t(1) << bits);
            const NativeInteger chosen = (up - target) < (target - down) ? up : down;
            used.insert(chosen);
            return chosen;
        };
        const uint32_t tree = SHIPProductTreeDepth(spec.hammingWeight);
        std::vector<NativeInteger> moduli;
        moduli.push_back(closest(spec.firstModBits));
        for (uint32_t i = 0; i < 1 + spec.multLevels; ++i)
            moduli.push_back(closest(spec.scalingModBits));
        moduli.push_back(closest(2 * spec.bootModBits - spec.scalingModBits));
        for (uint32_t i = 1; i < tree; ++i)
            moduli.push_back(closest(spec.bootModBits));
        std::vector<NativeInteger> roots;
        for (const auto& q : moduli)
            roots.push_back(RootOfUnity<NativeInteger>(m, q));
        params->SetElementParams(std::make_shared<ILDCRTParams<BigInteger>>(m, moduli, roots));
    }
    // OpenFHE fixes 60-bit auxiliary primes; SHIP's budget needs smaller ones (paper Table 2).
    params->PrecomputeCRTTables(params->GetKeySwitchTechnique(), params->GetScalingTechnique(),
                                params->GetEncryptionTechnique(), params->GetMultiplicationTechnique(),
                                params->GetNumPartQ(), spec.auxModBits, 0);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    if (spec.securityLevel != HEStd_NotSet) {
        const uint32_t logQP = SHIPLogQP(cc);
        const uint32_t nMin  = StdLatticeParm::FindRingDim(HEStd_ternary, spec.securityLevel, logQP);
        if (nMin == 0 || nMin > spec.ringDim)
            OPENFHE_THROW("SHIP context does not meet the requested HE-standard security: log2(QP) = " +
                          std::to_string(logQP) + ", ring dimension " + std::to_string(spec.ringDim));
    }
    return cc;
}

std::shared_ptr<SHIPBootstrapKey> SHIPKeyGen(const PrivateKey<DCRTPoly>& privateKey, const SHIPParams& sp) {
    if (!privateKey)
        OPENFHE_THROW("SHIPKeyGen: null private key");
    const CC cc       = privateKey->GetCryptoContext();
    const auto params = CkksParams(cc);
    if (params->GetScalingTechnique() != FIXEDMANUAL || params->GetKeySwitchTechnique() != HYBRID)
        OPENFHE_THROW("SHIP requires FIXEDMANUAL scaling and HYBRID key switching");
    const uint32_t N = cc->GetRingDimension();
    const uint32_t S = N / 2;
    if (cc->GetEncodingParams()->GetBatchSize() != S)
        OPENFHE_THROW("SHIP requires full packing (batch size N/2)");
    if (sp.hammingWeight == 0 || sp.hammingWeight > N || sp.columnSize == 0 || sp.muxBase < 2 ||
        2 * sp.window > N || sp.encapsulationModBits < 20 || sp.encapsulationModBits > 60)
        OPENFHE_THROW("Invalid SHIP parameters");
    const uint32_t sizeQ     = params->GetElementParams()->GetParams().size();
    const uint32_t treeDepth = SHIPProductTreeDepth(sp.hammingWeight);
    if (sizeQ < treeDepth + 2)
        OPENFHE_THROW("SHIP: modulus chain too short: need at least " + std::to_string(treeDepth + 2) + " primes");
    try {
        cc->GetEvalMultKeyVector(privateKey->GetKeyTag());
    }
    catch (...) {
        OPENFHE_THROW("SHIPKeyGen: call EvalMultKeyGen for this key first");
    }

    auto key        = std::make_shared<SHIPBootstrapKey>();
    key->params     = sp;
    key->denseTag   = privateKey->GetKeyTag();
    key->sparseTag  = privateKey->GetKeyTag() + "#ship-sparse";
    key->ringDim    = N;
    key->slots      = S;
    key->treeDepth  = treeDepth;
    auto& prng      = PseudoRandomNumberGenerator::GetPRNG();

    // Sparse secret with regularly spaced nonzero coefficients (Section 5.1).
    // position_k = (o_k + r_k) mod N with public o_k and secret r_k < R.
    const uint32_t h = sp.hammingWeight;
    const uint32_t R = sp.window ? 2 * sp.window : S;
    std::vector<uint32_t> offsets(h), secrets(h), positions(h);
    std::vector<int> signs(h);
    std::set<uint32_t> used;
    std::uniform_int_distribution<uint32_t> signDist(0, 1);
    for (uint32_t k = 0; k < h; ++k) {
        for (uint32_t attempt = 0;; ++attempt) {
            if (attempt > 100000)
                OPENFHE_THROW("SHIPKeyGen: cannot sample distinct sparse positions");
            uint32_t position;
            if (sp.window) {
                const uint64_t center = uint64_t(k) * N / h;
                const uint32_t origin = static_cast<uint32_t>((center + N - sp.window) % N);
                std::uniform_int_distribution<uint32_t> dist(0, R - 1);
                secrets[k] = dist(prng);
                offsets[k] = origin;
                position   = (origin + secrets[k]) % N;
            }
            else {
                std::uniform_int_distribution<uint32_t> dist(0, N - 1);
                position   = dist(prng);
                offsets[k] = 0;
                secrets[k] = position % S;
            }
            if (used.insert(position).second) {
                positions[k] = position;
                break;
            }
        }
        signs[k] = signDist(prng) ? 1 : -1;
    }

    // Dense-to-sparse encapsulation key modulo q0 * p'.
    {
        auto denseLimb   = privateKey->GetPrivateElement().GetElementAtIndex(0);
        auto paramsQ0    = denseLimb.GetParams();
        const auto q0    = paramsQ0->GetModulus();
        std::set<NativeInteger> taken;
        for (const auto& pr : params->GetParamsQP()->GetParams())
            taken.insert(pr->GetModulus());
        NativeInteger pPrime = LastPrime<NativeInteger>(sp.encapsulationModBits, 2 * N);
        while (taken.count(pPrime))
            pPrime = PreviousPrime<NativeInteger>(pPrime, 2 * N);
        auto paramsP = std::make_shared<ILNativeParams>(2 * N, pPrime, RootOfUnity<NativeInteger>(2 * N, pPrime));
        NativePoly sparseQ0(paramsQ0, Format::COEFFICIENT, true);
        for (uint32_t k = 0; k < h; ++k)
            sparseQ0[positions[k]] = signs[k] > 0 ? NativeInteger(1) : q0 - NativeInteger(1);
        auto sparseP = sparseQ0;
        sparseP.SwitchModulus(pPrime, paramsP->GetRootOfUnity(), 0, 0);
        denseLimb.SetFormat(Format::COEFFICIENT);
        auto denseP = denseLimb;
        denseP.SwitchModulus(pPrime, paramsP->GetRootOfUnity(), 0, 0);
        NativePoly eQ0(params->GetDiscreteGaussianGenerator(), paramsQ0, Format::COEFFICIENT);
        auto eP = eQ0;
        eP.SwitchModulus(pPrime, paramsP->GetRootOfUnity(), 0, 0);
        for (auto* poly : {&sparseQ0, &sparseP, &denseLimb, &eQ0, &eP})
            poly->SetFormat(Format::EVALUATION);
        NativePoly::DugType uniform;
        NativePoly aQ0(uniform, paramsQ0, Format::EVALUATION);
        NativePoly aP(uniform, paramsP, Format::EVALUATION);
        const NativeInteger pModQ0 = pPrime.Mod(q0);
        auto& enc     = key->encapsulation;
        enc.paramsQ0  = paramsQ0;
        enc.paramsP   = paramsP;
        enc.bQ0       = eQ0 - aQ0 * sparseQ0 + denseLimb * pModQ0;
        enc.aQ0       = std::move(aQ0);
        enc.bP        = eP - aP * sparseP;
        enc.aP        = std::move(aP);
        enc.pInvModQ0 = pModQ0.ModInverse(q0);
    }

    // Product-tree leaves at the scale of the top (boot) primes (paper Table 2: Boot primes > S2C/Mult).
    {
        const double topBits   = std::log2(params->GetElementParams()->GetParams().back()->GetModulus().ConvertToDouble());
        const double scaleBits = std::log2(params->GetScalingFactorReal(0));
        key->leafExtraBits     = topBits > scaleBits ? static_cast<uint32_t>(std::lround(topBits - scaleBits)) : 0;
    }
    // Column (with fused masks over PQ) and mux keys per nonzero coefficient.
    const auto paramsQP     = params->GetParamsQP();
    const DCRTPoly secretQP = SecretOverQP(privateKey, paramsQP);
    const uint32_t theta    = std::min(sp.columnSize, R);
    const uint64_t R1       = (R + theta - 1) / theta;
    const uint32_t digits   = CeilLog(R1, sp.muxBase);
    const uint32_t m        = 2 * N;
    const std::vector<double> zero(S, 0.0);
    key->factors.resize(h);
    for (uint32_t k = 0; k < h; ++k) {
        auto& fk         = key->factors[k];
        fk.offset        = offsets[k] % S;
        const uint32_t r = secrets[k];
        const uint32_t r0 = r % theta;
        uint64_t r1       = r / theta;
        const auto masks  = PreRotationMasks(positions[k], signs[k], S);
        fk.column.resize(theta);
        for (uint32_t i = 0; i < theta; ++i) {
            const uint32_t t = (fk.offset + i) % S;
            if (!key->columnPermutations.count(t)) {
                std::vector<uint32_t> perm(N);
                PrecomputeAutoMap(N, FindAutomorphismIndex2nComplex(-static_cast<int32_t>(t), m), &perm);
                key->columnPermutations.emplace(t, std::move(perm));
            }
            for (size_t band = 0; band < 4; ++band) {
                const auto& values = (i == r0) ? RotateRight(masks[band], (fk.offset + r0) % S) : zero;
                fk.column[i][band] = EncryptScaledByP(cc, secretQP, values, S, key->leafExtraBits);
            }
        }
        uint64_t step = theta;
        for (uint32_t t = 0; t < digits; ++t) {
            const uint32_t digit = static_cast<uint32_t>(r1 % sp.muxBase);
            r1 /= sp.muxBase;
            std::vector<SHIPMuxKey> branches;
            for (uint32_t d = 0; d < sp.muxBase; ++d) {
                const int32_t rotation = -static_cast<int32_t>((step * d) % S);  // right rotation
                branches.push_back(MakeMuxKey(cc, privateKey, d == digit ? 1 : 0, rotation));
            }
            fk.mux.push_back(std::move(branches));
            step *= sp.muxBase;
        }
    }

    // Conjugation (Algorithm 1, Step 22) and the S2C plan at the two lowest moduli.
    const uint32_t conjugationIndex = m - 1;
    cc->EvalAutomorphismKeyGen(privateKey, {conjugationIndex});
    key->conjugation = std::make_shared<std::map<uint32_t, EvalKey<DCRTPoly>>>();
    key->conjugation->emplace(conjugationIndex,
                              cc->GetEvalAutomorphismKeyMap(privateKey->GetKeyTag()).at(conjugationIndex));
    key->packing = MakePackingPlan(cc, privateKey, sizeQ - 2);
    return key;
}

Ciphertext<DCRTPoly> SHIPBootstrap(ConstCiphertext<DCRTPoly>& ciphertext, const SHIPBootstrapKey& key) {
    if (!ciphertext)
        OPENFHE_THROW("SHIPBootstrap: null ciphertext");
    const CC cc       = ciphertext->GetCryptoContext();
    const auto params = CkksParams(cc);
    const uint32_t S  = key.slots;
    if (cc->GetRingDimension() != key.ringDim || ciphertext->GetKeyTag() != key.denseTag)
        OPENFHE_THROW("SHIPBootstrap: ciphertext does not match the bootstrapping key");
    if (ciphertext->GetElements().size() != 2 || ciphertext->GetNoiseScaleDeg() != 1 ||
        ciphertext->GetSlots() != S || !std::isfinite(ciphertext->GetScalingFactor()) ||
        ciphertext->GetScalingFactor() <= 0)
        OPENFHE_THROW("SHIPBootstrap requires a fully packed, rescaled, degree-one ciphertext");
    const uint32_t towers = ciphertext->GetElements()[0].GetNumOfElements();
    if (towers < 2)
        OPENFHE_THROW("SHIPBootstrap requires at least two RNS limbs (the S2C level)");
    auto ct = ciphertext->Clone();
    if (towers > 2)
        cc->LevelReduceInPlace(ct, nullptr, towers - 2);
    auto coefficients = ApplyPacking(cc, ct, key.packing, S);  // S2C: one level, now at q0
    auto sparse       = Encapsulate(coefficients, key);
    const double q0   = sparse->GetElements()[0].GetElementAtIndex(0).GetModulus().ConvertToDouble();
    const double gamma = q0 / sparse->GetScalingFactor();
    auto real = HalfBootstrap(cc, sparse, key, gamma);
    if (key.params.realOnly)
        return real;
    auto imaginary = HalfBootstrap(cc, MultiplyByHalfMonomial(sparse, true), key, gamma);
    return cc->EvalAdd(real, MultiplyByHalfMonomial(imaginary, false));
}

uint64_t SHIPKeyStoredBytes(const SHIPBootstrapKey& key) {
    auto polyBytes = [](const DCRTPoly& p) {
        return uint64_t(p.GetRingDimension()) * p.GetNumOfElements() * sizeof(NativeInteger);
    };
    uint64_t total = 0;
    for (const auto& f : key.factors) {
        for (const auto& bands : f.column)
            for (const auto& pair : bands)
                total += polyBytes(pair[0]) + polyBytes(pair[1]);
        for (const auto& step : f.mux)
            for (const auto& mk : step)
                for (const auto& ek : {mk.body, mk.mask}) {
                    for (const auto& p : ek->GetAVector())
                        total += polyBytes(p);
                    for (const auto& p : ek->GetBVector())
                        total += polyBytes(p);
                }
    }
    for (const auto& pt : key.packing.diagonals)
        total += polyBytes(pt->GetElement<DCRTPoly>());
    total += 4 * uint64_t(key.ringDim) * sizeof(NativeInteger);
    return total;
}

void SHIPInsertBootstrapKey(const std::string& keyTag, std::shared_ptr<SHIPBootstrapKey> key) {
    std::lock_guard<std::mutex> lock(KeyStoreMutex());
    KeyStore()[keyTag] = std::move(key);
}

std::shared_ptr<SHIPBootstrapKey> SHIPGetBootstrapKey(const std::string& keyTag) {
    std::lock_guard<std::mutex> lock(KeyStoreMutex());
    auto it = KeyStore().find(keyTag);
    if (it == KeyStore().end())
        OPENFHE_THROW("No SHIP bootstrapping key for this key tag: call EvalSHIPBootstrapKeyGen first");
    return it->second;
}

void SHIPClearBootstrapKeys() {
    std::lock_guard<std::mutex> lock(KeyStoreMutex());
    KeyStore().clear();
}

}  // namespace lbcrypto
