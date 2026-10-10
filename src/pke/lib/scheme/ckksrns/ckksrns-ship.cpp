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
#include "utils/prng/blake2engine.h"
#include "utils/serial.h"

#include "cereal/types/array.hpp"
#include "cereal/types/string.hpp"

#include <array>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <map>
#include <mutex>
#include <random>
#include <set>
#include <sstream>
#include <vector>

namespace lbcrypto {

using SHIPSeed = default_prng::Blake2Engine::blake2_seed_array_t;

// Key-switching key over Q_{qLimbs} * P with digits covering the first qLimbs primes of Q.
// Digit j encrypts P * [Q/D_j]^{-1}... in the usual hybrid form (as KeySwitchHYBRID::KeySwitchGenInternal);
// a[j] is uniform and generated from (seed, id + j).
struct SHIPSwitchKey {
    uint32_t qLimbs = 0;
    uint64_t id     = 0;
    std::vector<DCRTPoly> b;
    std::vector<DCRTPoly> a;
};

// HMuxRot key (Definition 1 with gadget decomposition, Algorithm 5): the key switch targets sigma^{-1}(s)
// and sigma is applied afterwards, so the decomposition is shared (hoisted) across the B branches of one
// B-to-1 mux-rotate and ModDown is applied once per mux step.
struct SHIPMuxKey {
    SHIPSwitchKey body;  // encrypts beta
    SHIPSwitchKey mask;  // encrypts beta * s
    uint32_t automorphism = 1;
    std::vector<uint32_t> permutation;
};

// Automorphism (rotation or conjugation) key: key switch s -> sigma^{-1}(s), then sigma.
struct SHIPRotationKey {
    int32_t rotation      = 0;  // OpenFHE EvalRotate convention (left rotation); 0 for conjugation
    uint32_t automorphism = 1;
    SHIPSwitchKey key;
    std::vector<uint32_t> permutation;
};

// Enc_QP(P * 2^extra * v) with uniform component from (seed, id).
struct SHIPColumnKey {
    uint64_t id = 0;
    DCRTPoly b;
    DCRTPoly a;
};

struct SHIPFactorKey {
    uint32_t offset = 0;  // public rotation offset o (mod N/2)
    // column[i][band] = Enc_QP(P * 2^extra * 1_{i = r0} * Rot_{o + r0}(M'_band)), i < theta
    std::vector<std::array<SHIPColumnKey, 4>> column;
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
    std::map<int32_t, SHIPRotationKey> rotations;
};

class SHIPBootstrapKey {
public:
    SHIPParams params;
    CryptoContext<DCRTPoly> context;
    std::string denseTag;
    std::string sparseTag;
    uint32_t ringDim       = 0;
    uint32_t slots         = 0;
    uint32_t treeDepth     = 0;
    uint32_t leafExtraBits = 0;  // product-tree leaves have scale 2^(p + leafExtraBits)
    SHIPSeed seed{};
    SHIPEncapsulationKey encapsulation;
    std::vector<SHIPFactorKey> factors;  // in memory, or empty when factorDirectory is used
    std::string factorDirectory;
    uint32_t numFactors = 0;
    std::map<uint32_t, std::vector<uint32_t>> columnPermutations;  // right rotation t -> automorphism map
    SHIPRotationKey conjugation;
    SHIPPackingPlan packing;
};

namespace {

using CT     = Ciphertext<DCRTPoly>;
using CC     = CryptoContext<DCRTPoly>;
using C      = std::complex<double>;
using InArchive  = cereal::PortableBinaryInputArchive;
using OutArchive = cereal::PortableBinaryOutputArchive;

constexpr uint64_t kColumnId   = uint64_t(1) << 60;
constexpr uint64_t kMuxId      = uint64_t(2) << 60;
constexpr uint64_t kRotationId = uint64_t(3) << 60;
constexpr uint64_t kConjId     = uint64_t(4) << 60;
constexpr uint32_t kFormatVersion = 1;

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
        const uint32_t source   = (i + n - j) % n;
        const int effectiveSign = i < j ? -sign : sign;
        const uint32_t band     = (source >= S ? 2 : 0) + (effectiveSign < 0 ? 1 : 0);
        result[band][source % S] = 1.0;
    }
    return result;
}

// Q_{qLimbs} * P basis (first qLimbs primes of Q followed by all of P).
std::shared_ptr<DCRTPoly::Params> SubParams(const CC& cc, uint32_t qLimbs) {
    const auto params   = CkksParams(cc);
    const auto& primesQ = params->GetElementParams()->GetParams();
    const auto& primesP = params->GetParamsP()->GetParams();
    if (qLimbs == primesQ.size())
        return params->GetParamsQP();
    std::vector<NativeInteger> moduli, roots;
    for (uint32_t i = 0; i < qLimbs; ++i) {
        moduli.push_back(primesQ[i]->GetModulus());
        roots.push_back(primesQ[i]->GetRootOfUnity());
    }
    for (const auto& p : primesP) {
        moduli.push_back(p->GetModulus());
        roots.push_back(p->GetRootOfUnity());
    }
    return std::make_shared<DCRTPoly::Params>(2 * cc->GetRingDimension(), moduli, roots);
}

// Restriction of a polynomial over Q * P to the basis Q_{qLimbs} * P.
DCRTPoly Restrict(const DCRTPoly& fullQP, const std::shared_ptr<DCRTPoly::Params>& sub, uint32_t sizeQ,
                  uint32_t qLimbs) {
    DCRTPoly out(sub, fullQP.GetFormat(), true);
    const uint32_t sizeP = sub->GetParams().size() - qLimbs;
    for (uint32_t i = 0; i < qLimbs; ++i)
        out.SetElementAtIndex(i, fullQP.GetElementAtIndex(i));
    for (uint32_t j = 0; j < sizeP; ++j)
        out.SetElementAtIndex(qLimbs + j, fullQP.GetElementAtIndex(sizeQ + j));
    return out;
}

// Uniform polynomial (evaluation representation) derived from the public seed and an identifier.
DCRTPoly SeededUniform(const std::shared_ptr<DCRTPoly::Params>& params, const SHIPSeed& seed, uint64_t id) {
    SHIPSeed s = seed;
    s[14]      = static_cast<uint32_t>(id);
    s[15]      = static_cast<uint32_t>(id >> 32);
    default_prng::Blake2Engine engine(s, 0);
    DCRTPoly out(params, Format::EVALUATION, true);
    for (size_t k = 0; k < out.GetNumOfElements(); ++k) {
        auto limb           = out.GetElementAtIndex(k);
        const uint64_t q    = limb.GetModulus().ConvertToInt<uint64_t>();
        const uint32_t bits = limb.GetModulus().GetMSB();
        const uint64_t mask = bits >= 64 ? ~uint64_t(0) : ((uint64_t(1) << bits) - 1);
        for (size_t i = 0; i < limb.GetLength(); ++i) {
            uint64_t x;
            do {
                x = ((uint64_t(engine()) << 32) | engine()) & mask;
            } while (x >= q);
            limb[i] = NativeInteger(x);
        }
        out.SetElementAtIndex(k, std::move(limb));
    }
    return out;
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

// A secret over Q extended to QP, as in KeySwitchHYBRID::KeySwitchGenInternal.
DCRTPoly ExtendToQP(const DCRTPoly& secretQ, const std::shared_ptr<DCRTPoly::Params>& paramsQP) {
    auto secret = secretQ;
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

// Hybrid key-switching key from oldSecret (EVALUATION, over Q) to newSecretQP (EVALUATION, over QP),
// restricted to the digits and limbs of Q_{qLimbs}.
SHIPSwitchKey GenSwitchKey(const CC& cc, const DCRTPoly& oldSecret, const DCRTPoly& newSecretQP, uint32_t qLimbs,
                           const SHIPSeed& seed, uint64_t id) {
    const auto params    = CkksParams(cc);
    const uint32_t sizeQ = params->GetElementParams()->GetParams().size();
    const uint32_t alpha = params->GetNumPerPartQ();
    const auto& PModq    = params->GetPModq();
    const auto sub       = SubParams(cc, qLimbs);
    const DCRTPoly newSub = Restrict(newSecretQP, sub, sizeQ, qLimbs);
    SHIPSwitchKey key;
    key.qLimbs            = qLimbs;
    key.id                = id;
    const uint32_t digits = (qLimbs + alpha - 1) / alpha;
    for (uint32_t part = 0; part < digits; ++part) {
        DCRTPoly a = SeededUniform(sub, seed, id + part);
        DCRTPoly e(params->GetDiscreteGaussianGenerator(), sub, Format::EVALUATION);
        DCRTPoly b = e - a * newSub;
        const uint32_t start = alpha * part;
        const uint32_t end   = std::min(qLimbs, start + alpha);
        for (uint32_t i = start; i < end; ++i)
            b.SetElementAtIndex(i, b.GetElementAtIndex(i) + oldSecret.GetElementAtIndex(i) * PModq[i]);
        key.b.push_back(std::move(b));
        key.a.push_back(std::move(a));
    }
    return key;
}

// sigma^{-1}(s) for sigma: X -> X^k.
DCRTPoly InverseAutomorphism(const DCRTPoly& secret, uint32_t k) {
    const uint32_t n = secret.GetRingDimension();
    const uint32_t inverse = NativeInteger(k).ModInverse(NativeInteger(2 * n)).ConvertToInt();
    std::vector<uint32_t> perm(n);
    PrecomputeAutoMap(n, inverse, &perm);
    return secret.AutomorphismTransform(inverse, perm);
}

SHIPMuxKey MakeMuxKey(const CC& cc, const DCRTPoly& secret, const std::shared_ptr<DCRTPoly::Params>& paramsQP,
                      uint32_t beta, int32_t rotation, const SHIPSeed& seed, uint64_t id) {
    const uint32_t n     = secret.GetRingDimension();
    const uint32_t sizeQ = secret.GetNumOfElements();
    SHIPMuxKey key;
    key.automorphism = FindAutomorphismIndex2nComplex(rotation, 2 * n);
    key.permutation.resize(n);
    PrecomputeAutoMap(n, key.automorphism, &key.permutation);
    const DCRTPoly destination = ExtendToQP(InverseAutomorphism(secret, key.automorphism), paramsQP);
    DCRTPoly constant(secret.GetParams(), Format::COEFFICIENT, true);
    for (size_t i = 0; i < constant.GetNumOfElements(); ++i) {
        auto limb = constant.GetElementAtIndex(i);
        limb[0]   = NativeInteger(beta);
        constant.SetElementAtIndex(i, std::move(limb));
    }
    constant.SetFormat(Format::EVALUATION);
    const DCRTPoly zero(secret.GetParams(), Format::EVALUATION, true);
    key.body = GenSwitchKey(cc, constant, destination, sizeQ, seed, id);
    key.mask = GenSwitchKey(cc, beta ? secret : zero, destination, sizeQ, seed, id + 64);
    return key;
}

SHIPRotationKey MakeAutomorphismKey(const CC& cc, const DCRTPoly& secret,
                                    const std::shared_ptr<DCRTPoly::Params>& paramsQP, uint32_t automorphism,
                                    int32_t rotation, uint32_t qLimbs, const SHIPSeed& seed, uint64_t id) {
    const uint32_t n = secret.GetRingDimension();
    SHIPRotationKey key;
    key.rotation     = rotation;
    key.automorphism = automorphism;
    key.permutation.resize(n);
    PrecomputeAutoMap(n, automorphism, &key.permutation);
    const DCRTPoly destination = ExtendToQP(InverseAutomorphism(secret, automorphism), paramsQP);
    key.key                    = GenSwitchKey(cc, secret, destination, qLimbs, seed, id);
    return key;
}

// Enc_QP(P * 2^extraBits * v): v slot-encoded with scaling factor P * 2^extraBits under the output secret.
// round(m_i * P * 2^extraBits / Delta) from the exact Delta-encoding m_i; relative error <= 2^-log2(Delta).
// After the product with a Delta-scaled phase and Rescale_P, the factor has scale 2^extraBits * Delta.
SHIPColumnKey EncryptScaledByP(const CC& cc, const DCRTPoly& secretQP, const std::vector<double>& values,
                               uint32_t S, uint32_t extraBits, const SHIPSeed& seed, uint64_t id) {
    const auto params     = CkksParams(cc);
    const auto paramsQP   = params->GetParamsQP();
    const uint32_t bottom = params->GetElementParams()->GetParams().size() - 1;
    DCRTPoly plain(paramsQP, Format::COEFFICIENT, true);
    bool nonzero = false;
    for (double v : values)
        nonzero = nonzero || v != 0.0;
    if (nonzero) {
        auto pt      = cc->MakeCKKSPackedPlaintext(values, 1, bottom, nullptr, S);
        int exponent = 0;
        if (std::frexp(pt->GetScalingFactor(), &exponent) != 0.5 || exponent < 2)
            OPENFHE_THROW("SHIP requires a power-of-two CKKS scaling factor");
        const uint32_t deltaBits = exponent - 1;
        auto encoded             = pt->GetElement<DCRTPoly>();
        encoded.SetFormat(Format::COEFFICIENT);
        const auto limb       = encoded.GetElementAtIndex(0);
        const auto q0         = limb.GetModulus();
        const BigInteger P    = params->GetParamsP()->GetModulus();
        const BigInteger half = BigInteger(1).LShift(deltaBits - 1);
        std::vector<NativePoly> limbs;
        for (size_t k = 0; k < plain.GetNumOfElements(); ++k)
            limbs.push_back(plain.GetElementAtIndex(k));
        for (size_t i = 0; i < limb.GetLength(); ++i) {
            if (limb[i] == NativeInteger(0))
                continue;
            const bool negative      = limb[i] > (q0 >> 1);
            const uint64_t magnitude = (negative ? q0 - limb[i] : limb[i]).ConvertToInt<uint64_t>();
            const BigInteger scaled  = ((BigInteger(magnitude) * P).LShift(extraBits) + half).RShift(deltaBits);
            for (auto& out : limbs) {
                const auto r = out.GetModulus();
                NativeInteger value(scaled.Mod(BigInteger(r.ConvertToInt<uint64_t>())).ConvertToInt<uint64_t>());
                out[i] = (negative && value != NativeInteger(0)) ? r - value : value;
            }
        }
        for (size_t k = 0; k < limbs.size(); ++k)
            plain.SetElementAtIndex(k, std::move(limbs[k]));
    }
    plain.SetFormat(Format::EVALUATION);
    SHIPColumnKey key;
    key.id = id;
    key.a  = SeededUniform(paramsQP, seed, id);
    DCRTPoly e(params->GetDiscreteGaussianGenerator(), paramsQP, Format::EVALUATION);
    key.b = e - key.a * secretQP + plain;
    return key;
}

// Inner product of hoisted digits with a switching key, result over Q_l P (no ModDown).
std::array<DCRTPoly, 2> InnerExt(const std::vector<DCRTPoly>& digits, const SHIPSwitchKey& key) {
    const auto paramsQlP   = digits[0].GetParams();
    const uint32_t sizeQlP = paramsQlP->GetParams().size();
    const uint32_t sizeP   = key.b[0].GetNumOfElements() - key.qLimbs;
    const uint32_t sizeQl  = sizeQlP - sizeP;
    if (digits.size() > key.b.size() || sizeQl > key.qLimbs)
        OPENFHE_THROW("SHIP: switching key does not cover this level");
    DCRTPoly c0(paramsQlP, Format::EVALUATION, true), c1(paramsQlP, Format::EVALUATION, true);
    for (size_t j = 0; j < digits.size(); ++j) {
        const auto& cj = digits[j];
        for (uint32_t i = 0; i < sizeQlP; ++i) {
            const uint32_t idx = i < sizeQl ? i : key.qLimbs + (i - sizeQl);
            const auto& cji    = cj.GetElementAtIndex(i);
            c0.SetElementAtIndex(i, c0.GetElementAtIndex(i) + cji * key.b[j].GetElementAtIndex(idx));
            c1.SetElementAtIndex(i, c1.GetElementAtIndex(i) + cji * key.a[j].GetElementAtIndex(idx));
        }
    }
    return {std::move(c0), std::move(c1)};
}

// sum over the B branches of one mux step, sharing the decompositions; one ModDown.
CT MuxRotate(const CC& cc, const std::vector<std::vector<SHIPMuxKey>>& steps, CT input) {
    const auto scheme = cc->GetScheme();
    for (const auto& branches : steps) {
        auto body = scheme->EvalKeySwitchPrecomputeCore(input->GetElements()[0], input->GetCryptoParameters());
        auto mask = scheme->EvalKeySwitchPrecomputeCore(input->GetElements()[1], input->GetCryptoParameters());
        std::vector<DCRTPoly> sum;
        for (const auto& key : branches) {
            auto b = InnerExt(*body, key.body);
            auto a = InnerExt(*mask, key.mask);
            std::vector<DCRTPoly> term{b[0] + a[0], b[1] + a[1]};
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

// Automorphism with a SHIP key, given the hoisted decomposition of the second component.
CT ApplyAutomorphism(const CC& cc, const CT& input, const std::vector<DCRTPoly>& digitsA,
                     const SHIPRotationKey& key) {
    const auto scheme = cc->GetScheme();
    auto ext          = InnerExt(digitsA, key.key);
    auto extended     = input->CloneEmpty();
    extended->SetElements({std::move(ext[0]), std::move(ext[1])});
    auto down  = scheme->KeySwitchDown(extended);
    DCRTPoly b = input->GetElements()[0] + down->GetElements()[0];
    DCRTPoly a = down->GetElements()[1];
    auto out   = input->CloneEmpty();
    out->SetElements({b.AutomorphismTransform(key.automorphism, key.permutation),
                      a.AutomorphismTransform(key.automorphism, key.permutation)});
    return out;
}

CT ApplyAutomorphism(const CC& cc, const CT& input, const SHIPRotationKey& key) {
    auto digits = cc->GetScheme()->EvalKeySwitchPrecomputeCore(input->GetElements()[1], input->GetCryptoParameters());
    return ApplyAutomorphism(cc, input, *digits, key);
}

// FIXEDMANUAL assumes every rescaling prime equals Delta. With small NTT-friendly primes the
// deviation |q/Delta - 1| reaches a few percent and compounds through the product tree. All leaves
// enter at level 0 with a known scale, so the true root scale is deterministic: simulate the
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
            auto a            = terms[2 * i]->Clone();
            auto b            = terms[2 * i + 1]->Clone();
            const auto target = std::max(a->GetLevel(), b->GetLevel());
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

uint32_t BabyStep(uint32_t S) {
    uint32_t baby = 1;
    while (baby * baby < S)
        baby <<= 1;
    return baby;
}

std::vector<int32_t> PackingRotations(uint32_t S) {
    const uint32_t baby = BabyStep(S);
    std::vector<int32_t> rotations;
    for (uint32_t j = 1; j < baby; ++j)
        rotations.push_back(j);
    for (uint32_t g = baby; g < S; g += baby)
        rotations.push_back(g);
    return rotations;
}

// S2C: slots mu -> coefficients (Re mu | Im mu), i.e. multiplication of the slot vector by
// V[k][j] = exp(2 pi i 5^k j / (2N)) (paper Section 3.2, dense matrix with BSGS, one level).
// The diagonals are scaled by q/Delta (exact output scale) and by 1/messageBound.
std::vector<Plaintext> PackingDiagonals(const CC& cc, uint32_t level, double messageBound) {
    const uint32_t S    = cc->GetRingDimension() / 2;
    const uint32_t cycl = 4 * S;
    const uint32_t baby = BabyStep(S);
    std::vector<uint32_t> exponents(S);
    uint64_t e = 1;
    for (auto& x : exponents) {
        x = e;
        e = (5 * e) % cycl;
    }
    const double tau       = 2 * std::acos(-1.0);
    const auto params      = CkksParams(cc);
    const auto& moduli     = params->GetElementParams()->GetParams();
    const double rescaleBy = moduli[moduli.size() - 1 - level]->GetModulus().ConvertToDouble();
    const double ratio     = rescaleBy / params->GetScalingFactorReal(level) / messageBound;
    std::vector<Plaintext> diagonals(S);
    for (uint32_t r = 0; r < S; ++r) {
        const uint32_t giant = (r / baby) * baby;
        std::vector<C> d(S);
        for (uint32_t k = 0; k < S; ++k) {
            const uint32_t row      = (k + S - giant) % S;
            const uint32_t column   = (row + r) % S;
            const uint32_t exponent = (uint64_t(exponents[row]) * column) % cycl;
            d[k]                    = std::polar(ratio, tau * exponent / cycl);
        }
        diagonals[r] = cc->MakeCKKSPackedPlaintext(d, 1, level, nullptr, S);
    }
    return diagonals;
}

CT ApplyPacking(const CC& cc, const CT& input, const SHIPPackingPlan& plan, uint32_t S) {
    std::vector<CT> baby(plan.babyStep);
    baby[0]     = input;
    auto digits = cc->GetScheme()->EvalKeySwitchPrecomputeCore(input->GetElements()[1], input->GetCryptoParameters());
#pragma omp parallel for schedule(dynamic)
    for (uint32_t j = 1; j < plan.babyStep; ++j)
        baby[j] = ApplyAutomorphism(cc, input, *digits, plan.rotations.at(j));
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
        partial[gi] = g ? ApplyAutomorphism(cc, group, plan.rotations.at(g)) : group;
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
    auto newB     = b + ModDownByP(aQ0 * k.bQ0, aP * k.bP, k);
    auto newA     = ModDownByP(aQ0 * k.aQ0, aP * k.aP, k);
    auto elements = input->GetElements();
    elements[0].SetElementAtIndex(0, std::move(newB));
    elements[1].SetElementAtIndex(0, std::move(newA));
    auto result = input->Clone();
    result->SetElements(std::move(elements));
    result->SetKeyTag(key.sparseTag);
    return result;
}

// ---------------------------------------------------------------------------------------------
// Serialization helpers
// ---------------------------------------------------------------------------------------------

void SaveSwitchKey(OutArchive& ar, const SHIPSwitchKey& k, bool compact) {
    ar(k.qLimbs, k.id, static_cast<uint32_t>(k.b.size()));
    for (const auto& p : k.b)
        ar(p);
    if (!compact)
        for (const auto& p : k.a)
            ar(p);
}

SHIPSwitchKey LoadSwitchKey(InArchive& ar, bool compact, const CC& cc, const SHIPSeed& seed) {
    SHIPSwitchKey k;
    uint32_t count = 0;
    ar(k.qLimbs, k.id, count);
    k.b.resize(count);
    for (auto& p : k.b)
        ar(p);
    if (compact) {
        const auto sub = SubParams(cc, k.qLimbs);
        for (uint32_t j = 0; j < count; ++j)
            k.a.push_back(SeededUniform(sub, seed, k.id + j));
    }
    else {
        k.a.resize(count);
        for (auto& p : k.a)
            ar(p);
    }
    return k;
}

void SaveRotationKey(OutArchive& ar, const SHIPRotationKey& k, bool compact) {
    ar(k.rotation, k.automorphism);
    SaveSwitchKey(ar, k.key, compact);
}

SHIPRotationKey LoadRotationKey(InArchive& ar, bool compact, const CC& cc, const SHIPSeed& seed) {
    SHIPRotationKey k;
    ar(k.rotation, k.automorphism);
    k.key                = LoadSwitchKey(ar, compact, cc, seed);
    const uint32_t n     = cc->GetRingDimension();
    k.permutation.resize(n);
    PrecomputeAutoMap(n, k.automorphism, &k.permutation);
    return k;
}

void SaveFactor(OutArchive& ar, const SHIPFactorKey& f, bool compact) {
    ar(f.offset, static_cast<uint32_t>(f.column.size()), static_cast<uint32_t>(f.mux.size()));
    for (const auto& bands : f.column)
        for (const auto& c : bands) {
            ar(c.id, c.b);
            if (!compact)
                ar(c.a);
        }
    for (const auto& step : f.mux) {
        ar(static_cast<uint32_t>(step.size()));
        for (const auto& mk : step) {
            ar(mk.automorphism);
            SaveSwitchKey(ar, mk.body, compact);
            SaveSwitchKey(ar, mk.mask, compact);
        }
    }
}

SHIPFactorKey LoadFactor(InArchive& ar, bool compact, const CC& cc, const SHIPSeed& seed) {
    SHIPFactorKey f;
    uint32_t columns = 0, steps = 0;
    ar(f.offset, columns, steps);
    const auto paramsQP = CkksParams(cc)->GetParamsQP();
    f.column.resize(columns);
    for (auto& bands : f.column)
        for (auto& c : bands) {
            ar(c.id, c.b);
            if (compact)
                c.a = SeededUniform(paramsQP, seed, c.id);
            else
                ar(c.a);
        }
    const uint32_t n = cc->GetRingDimension();
    f.mux.resize(steps);
    for (auto& step : f.mux) {
        uint32_t branches = 0;
        ar(branches);
        step.resize(branches);
        for (auto& mk : step) {
            ar(mk.automorphism);
            mk.body = LoadSwitchKey(ar, compact, cc, seed);
            mk.mask = LoadSwitchKey(ar, compact, cc, seed);
            mk.permutation.resize(n);
            PrecomputeAutoMap(n, mk.automorphism, &mk.permutation);
        }
    }
    return f;
}

std::string FactorPath(const std::string& directory, uint32_t index) {
    std::ostringstream name;
    name << directory << "/factor-" << index << ".ship";
    return name.str();
}

std::shared_ptr<const SHIPFactorKey> GetFactor(const CC& cc, const SHIPBootstrapKey& key, uint32_t index) {
    if (key.factorDirectory.empty())
        return std::shared_ptr<const SHIPFactorKey>(&key.factors[index], [](const SHIPFactorKey*) {});
    std::ifstream file(FactorPath(key.factorDirectory, index), std::ios::binary);
    if (!file)
        OPENFHE_THROW("SHIP: cannot open factor key file " + FactorPath(key.factorDirectory, index));
    InArchive ar(file);
    bool compact = false;
    ar(compact);
    return std::make_shared<SHIPFactorKey>(LoadFactor(ar, compact, cc, key.seed));
}

void ComputeColumnPermutations(SHIPBootstrapKey& key, uint32_t offset, uint32_t theta) {
    const uint32_t S = key.slots, N = key.ringDim;
    for (uint32_t i = 0; i < theta; ++i) {
        const uint32_t t = (offset + i) % S;
        if (!key.columnPermutations.count(t)) {
            std::vector<uint32_t> perm(N);
            PrecomputeAutoMap(N, FindAutomorphismIndex2nComplex(-static_cast<int32_t>(t), 2 * N), &perm);
            key.columnPermutations.emplace(t, std::move(perm));
        }
    }
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
    // Exact compensation of the product-tree scale drift and of messageBound through pt_0.
    {
        std::vector<long double> moduli;
        for (const auto& q : params->GetElementParams()->GetParams())
            moduli.push_back(q->GetModulus().ConvertToDouble());
        const long double delta     = params->GetScalingFactorReal(0);
        const long double leafScale = std::ldexp(delta, key.leafExtraBits);
        const long double kappa     = ProductTreeScaleRatio(key.numFactors + 1, leafScale, delta, moduli);
        const long double factor    = std::ldexp(1.0L, key.leafExtraBits) / kappa * key.params.messageBound;
        for (auto& v : initial)
            v *= static_cast<double>(factor);
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
    std::vector<CT> factors(key.numFactors + 1);
    factors[0] = trivial;
#pragma omp parallel for schedule(dynamic)
    for (uint32_t f = 0; f < key.numFactors; ++f) {
        const auto fk = GetFactor(cc, key, f);
        DCRTPoly sb(paramsQP, Format::EVALUATION, true), sa(paramsQP, Format::EVALUATION, true);
        for (uint32_t i = 0; i < fk->column.size(); ++i) {
            const uint32_t t     = (fk->offset + i) % S;
            const uint32_t index = FindAutomorphismIndex2nComplex(-static_cast<int32_t>(t), m);
            const auto& perm     = key.columnPermutations.at(t);
            for (size_t band = 0; band < 4; ++band) {
                const auto rotated = t ? phaseQP[band].AutomorphismTransform(index, perm) : phaseQP[band];
                sb += fk->column[i][band].b * rotated;
                sa += fk->column[i][band].a * rotated;
            }
        }
        auto extended = trivial->CloneEmpty();
        extended->SetElements({std::move(sb), std::move(sa)});
        auto selected  = scheme->KeySwitchDown(extended);  // Rescale_P: level stays 0
        factors[f + 1] = MuxRotate(cc, fk->mux, selected);
    }
    auto root      = ProductTree(cc, std::move(factors));
    auto conjugate = ApplyAutomorphism(cc, root, key.conjugation);
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

uint32_t SparseSecretBound(uint32_t ringDim) {
    // Largest log2(q0 p') for which the dense-to-sparse switching key (h = 31) keeps 128-bit security.
    // The paper (Section 5.2) gives 55, 100, 105 for N = 2^13, 2^14, 2^15. The lattice estimator
    // (commit d8c00b48, 2026-08-19, LWE.estimate attacks) gives 2^121.4 and 2^126.2 for the first two
    // (bdd_mitm_hybrid), so the bounds below are the estimator values: 42 -> 2^129.8, 88 -> 2^130.2,
    // 105 -> 2^148.7 (research/ship/security).
    switch (ringDim) {
        case 1 << 13:
            return 42;
        case 1 << 14:
            return 88;
        case 1 << 15:
            return 105;
        default:
            return 0;
    }
}

}  // namespace

// ---------------------------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------------------------

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
    s.numLargeDigits = 4;
    s.auxModBits     = 50;
    return s;
}

SHIPContextSpec SHIPContextSpec::HT14() {
    // Paper HT14: Base 28, S2C 23, Mult 23 x 9, Boot 25 x 5, Aux 29 x 2, dnum 8, log PQ = 438.
    // Here: Base 28, S2C/Mult 23, Boot 24 x 4 + 25, Aux 55 x 1, dnum 8 (one auxiliary prime suffices).
    SHIPContextSpec s;
    s.ringDim        = 1 << 14;
    s.firstModBits   = 28;
    s.scalingModBits = 23;
    s.bootModBits    = 24;
    s.multLevels     = 9;
    s.numLargeDigits = 8;
    s.auxModBits     = 55;
    return s;
}

SHIPContextSpec SHIPContextSpec::HT15() {
    // Paper HT15: Base 50, S2C 39, Mult 39 x 9, Boot 46 x 5, Aux 53 x 4, dnum 4, log PQ = 881.
    // Here: Base 50, S2C/Mult 39, Boot 46 x 4 + 53, Aux 50 x 3, dnum 6 (fits the 881-bit budget).
    SHIPContextSpec s;
    s.ringDim        = 1 << 15;
    s.firstModBits   = 50;
    s.scalingModBits = 39;
    s.bootModBits    = 46;
    s.multLevels     = 9;
    s.numLargeDigits = 6;
    s.auxModBits     = 50;
    return s;
}

SHIPParams SHIPParams::Recommended(const SHIPContextSpec& spec) {
    const uint32_t bound = SparseSecretBound(spec.ringDim);
    if (!bound)
        OPENFHE_THROW("SHIPParams::Recommended: no published sparse-secret bound for this ring dimension");
    SHIPParams p;
    p.hammingWeight        = spec.hammingWeight;
    p.window               = std::max<uint32_t>(175, spec.ringDim / (2 * spec.hammingWeight));
    p.columnSize           = spec.ringDim == (1 << 13) ? 6 : spec.ringDim == (1 << 14) ? 9 : 17;
    p.muxBase              = 4;
    p.encapsulationModBits = 0;  // chosen from the actual q0 by SHIPKeyGen
    p.realOnly             = true;
    return p;
}

uint32_t SHIPLogQP(const CryptoContext<DCRTPoly>& cc) {
    return CkksParams(cc)->GetParamsQP()->GetModulus().GetMSB();
}

CryptoContext<DCRTPoly> GenSHIPCryptoContext(const SHIPContextSpec& spec) {
    if (spec.ringDim < 16 || (spec.ringDim & (spec.ringDim - 1)) || spec.firstModBits <= spec.scalingModBits ||
        spec.numLargeDigits == 0 || spec.auxModBits == 0 ||
        (spec.bootModBits > spec.scalingModBits && 2 * spec.bootModBits - spec.scalingModBits > 60))
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

uint64_t SHIPEstimateKeyBytes(const CryptoContext<DCRTPoly>& cc, const SHIPParams& sp) {
    const auto params     = CkksParams(cc);
    const uint64_t N      = cc->GetRingDimension();
    const uint64_t S      = N / 2;
    const uint64_t sizeQ  = params->GetElementParams()->GetParams().size();
    const uint64_t sizeP  = params->GetParamsP()->GetParams().size();
    const uint64_t alpha  = params->GetNumPerPartQ();
    const uint64_t dnum   = (sizeQ + alpha - 1) / alpha;
    const uint64_t word   = sizeof(NativeInteger);
    const uint64_t R      = sp.window ? 2 * sp.window : S;
    const uint64_t theta  = std::min<uint64_t>(sp.columnSize, R);
    const uint64_t digits = CeilLog((R + theta - 1) / theta, sp.muxBase);
    const uint64_t qp     = N * (sizeQ + sizeP) * word;
    const uint64_t perFactor = theta * 4 * 2 * qp + digits * sp.muxBase * 2 * dnum * 2 * qp;
    const uint64_t tree      = SHIPProductTreeDepth(sp.hammingWeight);
    const uint64_t rotations = PackingRotations(S).size() * ((2 + alpha - 1) / alpha) * 2 * N * (2 + sizeP) * word;
    const uint64_t conjLimbs = sizeQ - tree;
    const uint64_t conj      = ((conjLimbs + alpha - 1) / alpha) * 2 * N * (conjLimbs + sizeP) * word;
    const uint64_t diagonals = S * 2 * N * word;
    return sp.hammingWeight * perFactor + rotations + conj + diagonals;
}

std::shared_ptr<SHIPBootstrapKey> SHIPKeyGen(const PrivateKey<DCRTPoly>& privateKey, const SHIPParams& sp,
                                             const std::string& factorDirectory) {
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
        2 * sp.window > N || (sp.encapsulationModBits != 0 && (sp.encapsulationModBits < 17 || sp.encapsulationModBits > 60)) ||
        !(sp.messageBound >= 1.0) || !std::isfinite(sp.messageBound))
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

    auto key       = std::make_shared<SHIPBootstrapKey>();
    key->params    = sp;
    key->context   = cc;
    key->denseTag  = privateKey->GetKeyTag();
    key->sparseTag = privateKey->GetKeyTag() + "#ship-sparse";
    key->ringDim   = N;
    key->slots     = S;
    key->treeDepth = treeDepth;
    auto& prng     = PseudoRandomNumberGenerator::GetPRNG();
    for (auto& w : key->seed)
        w = prng();

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
        auto denseLimb = privateKey->GetPrivateElement().GetElementAtIndex(0);
        auto paramsQ0  = denseLimb.GetParams();
        const auto q0  = paramsQ0->GetModulus();
        std::set<NativeInteger> taken;
        for (const auto& pr : params->GetParamsQP()->GetParams())
            taken.insert(pr->GetModulus());
        uint32_t pBits = sp.encapsulationModBits;
        if (pBits == 0) {
            const uint32_t bound = SparseSecretBound(N);
            if (!bound || bound < q0.GetMSB() + 17)
                OPENFHE_THROW("SHIPKeyGen: no published sparse-secret bound to choose p' for this ring dimension");
            pBits = std::min<uint32_t>(60, bound - q0.GetMSB());
        }
        key->params.encapsulationModBits = pBits;
        NativeInteger pPrime = LastPrime<NativeInteger>(pBits, 2 * N);
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
        auto& enc                  = key->encapsulation;
        enc.paramsQ0               = paramsQ0;
        enc.paramsP                = paramsP;
        enc.bQ0                    = eQ0 - aQ0 * sparseQ0 + denseLimb * pModQ0;
        enc.aQ0                    = std::move(aQ0);
        enc.bP                     = eP - aP * sparseP;
        enc.aP                     = std::move(aP);
        enc.pInvModQ0              = pModQ0.ModInverse(q0);
    }

    // Product-tree leaves at the scale of the top (boot) primes (paper Table 2: Boot primes > S2C/Mult).
    {
        const double topBits   = std::log2(params->GetElementParams()->GetParams().back()->GetModulus().ConvertToDouble());
        const double scaleBits = std::log2(params->GetScalingFactorReal(0));
        key->leafExtraBits     = topBits > scaleBits ? static_cast<uint32_t>(std::lround(topBits - scaleBits)) : 0;
    }

    // Column (with fused masks over PQ) and mux keys per nonzero coefficient.
    const auto paramsQP     = params->GetParamsQP();
    const DCRTPoly& secret  = privateKey->GetPrivateElement();
    const DCRTPoly secretQP = ExtendToQP(secret, paramsQP);
    const uint32_t theta    = std::min(sp.columnSize, R);
    const uint64_t R1       = (R + theta - 1) / theta;
    const uint32_t digits   = CeilLog(R1, sp.muxBase);
    const std::vector<double> zero(S, 0.0);
    key->numFactors      = h;
    key->factorDirectory = factorDirectory;
    if (!factorDirectory.empty())
        std::filesystem::create_directories(factorDirectory);
    else
        key->factors.resize(h);
    for (uint32_t k = 0; k < h; ++k) {
        SHIPFactorKey fk;
        fk.offset         = offsets[k] % S;
        const uint32_t r  = secrets[k];
        const uint32_t r0 = r % theta;
        uint64_t r1       = r / theta;
        const auto masks  = PreRotationMasks(positions[k], signs[k], S);
        ComputeColumnPermutations(*key, fk.offset, theta);
        fk.column.resize(theta);
        for (uint32_t i = 0; i < theta; ++i)
            for (size_t band = 0; band < 4; ++band) {
                const auto& values = (i == r0) ? RotateRight(masks[band], (fk.offset + r0) % S) : zero;
                const uint64_t id  = kColumnId | (uint64_t(k) << 24) | (uint64_t(i) << 8) | band;
                fk.column[i][band] = EncryptScaledByP(cc, secretQP, values, S, key->leafExtraBits, key->seed, id);
            }
        uint64_t step = theta;
        for (uint32_t t = 0; t < digits; ++t) {
            const uint32_t digit = static_cast<uint32_t>(r1 % sp.muxBase);
            r1 /= sp.muxBase;
            std::vector<SHIPMuxKey> branches;
            for (uint32_t d = 0; d < sp.muxBase; ++d) {
                const int32_t rotation = -static_cast<int32_t>((step * d) % S);  // right rotation
                const uint64_t id      = kMuxId | (uint64_t(k) << 32) | (uint64_t(t) << 24) | (uint64_t(d) << 16);
                branches.push_back(MakeMuxKey(cc, secret, paramsQP, d == digit ? 1 : 0, rotation, key->seed, id));
            }
            fk.mux.push_back(std::move(branches));
            step *= sp.muxBase;
        }
        if (factorDirectory.empty()) {
            key->factors[k] = std::move(fk);
        }
        else {
            std::ofstream file(FactorPath(factorDirectory, k), std::ios::binary | std::ios::trunc);
            if (!file)
                OPENFHE_THROW("SHIPKeyGen: cannot write " + FactorPath(factorDirectory, k));
            OutArchive ar(file);
            const bool compact = false;  // streamed at every bootstrap: store both components
            ar(compact);
            SaveFactor(ar, fk, compact);
        }
    }

    // Conjugation (Algorithm 1, Step 22) at the product-tree output level, and S2C at the two lowest moduli.
    const uint32_t m = 2 * N;
    key->conjugation = MakeAutomorphismKey(cc, secret, paramsQP, m - 1, 0, sizeQ - treeDepth, key->seed, kConjId);
    key->packing.babyStep   = BabyStep(S);
    key->packing.inputLevel = sizeQ - 2;
    for (int32_t r : PackingRotations(S)) {
        const uint64_t id = kRotationId | (uint64_t(r) << 8);
        key->packing.rotations.emplace(
            r, MakeAutomorphismKey(cc, secret, paramsQP, FindAutomorphismIndex2nComplex(r, m), r, 2, key->seed, id));
    }
    key->packing.diagonals = PackingDiagonals(cc, key->packing.inputLevel, sp.messageBound);
    return key;
}

Ciphertext<DCRTPoly> SHIPBootstrap(ConstCiphertext<DCRTPoly>& ciphertext, const SHIPBootstrapKey& key) {
    if (!ciphertext)
        OPENFHE_THROW("SHIPBootstrap: null ciphertext");
    const CC cc      = ciphertext->GetCryptoContext();
    const uint32_t S = key.slots;
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
    auto coefficients  = ApplyPacking(cc, ct, key.packing, S);  // S2C: one level, now at q0
    auto sparse        = Encapsulate(coefficients, key);
    const double q0    = sparse->GetElements()[0].GetElementAtIndex(0).GetModulus().ConvertToDouble();
    const double gamma = q0 / sparse->GetScalingFactor();
    auto real          = HalfBootstrap(cc, sparse, key, gamma);
    if (key.params.realOnly)
        return real;
    auto imaginary = HalfBootstrap(cc, MultiplyByHalfMonomial(sparse, true), key, gamma);
    return cc->EvalAdd(real, MultiplyByHalfMonomial(imaginary, false));
}

uint64_t SHIPKeyStoredBytes(const SHIPBootstrapKey& key) {
    auto polyBytes = [](const DCRTPoly& p) {
        return uint64_t(p.GetRingDimension()) * p.GetNumOfElements() * sizeof(NativeInteger);
    };
    auto switchBytes = [&](const SHIPSwitchKey& k) {
        uint64_t t = 0;
        for (const auto& p : k.b)
            t += polyBytes(p);
        for (const auto& p : k.a)
            t += polyBytes(p);
        return t;
    };
    uint64_t total = 0;
    for (const auto& f : key.factors) {
        for (const auto& bands : f.column)
            for (const auto& c : bands)
                total += polyBytes(c.b) + polyBytes(c.a);
        for (const auto& step : f.mux)
            for (const auto& mk : step)
                total += switchBytes(mk.body) + switchBytes(mk.mask);
    }
    for (const auto& [r, rk] : key.packing.rotations)
        total += switchBytes(rk.key);
    total += switchBytes(key.conjugation.key);
    for (const auto& pt : key.packing.diagonals)
        total += polyBytes(pt->GetElement<DCRTPoly>());
    total += 4 * uint64_t(key.ringDim) * sizeof(NativeInteger);
    return total;
}

void SHIPSerializeBootstrapKey(std::ostream& os, const SHIPBootstrapKey& key, bool compact) {
    OutArchive ar(os);
    const std::string magic = "OpenFHE-SHIP";
    ar(magic, kFormatVersion, compact);
    const auto& p = key.params;
    ar(p.hammingWeight, p.window, p.columnSize, p.muxBase, p.encapsulationModBits, p.realOnly, p.messageBound);
    ar(key.denseTag, key.sparseTag, key.ringDim, key.slots, key.treeDepth, key.leafExtraBits, key.seed);
    const auto& enc      = key.encapsulation;
    const uint64_t prime = enc.paramsP->GetModulus().ConvertToInt<uint64_t>();
    const uint64_t root  = enc.paramsP->GetRootOfUnity().ConvertToInt<uint64_t>();
    ar(prime, root);
    ar(enc.bQ0, enc.aQ0, enc.bP, enc.aP);
    ar(key.packing.babyStep, key.packing.inputLevel, static_cast<uint32_t>(key.packing.rotations.size()));
    for (const auto& [r, rk] : key.packing.rotations)
        SaveRotationKey(ar, rk, compact);
    SaveRotationKey(ar, key.conjugation, compact);
    ar(key.numFactors);
    for (uint32_t f = 0; f < key.numFactors; ++f) {
        if (key.factorDirectory.empty())
            SaveFactor(ar, key.factors[f], compact);
        else
            SaveFactor(ar, *GetFactor(key.context, key, f), compact);
    }
}

std::shared_ptr<SHIPBootstrapKey> SHIPDeserializeBootstrapKey(std::istream& is, const CryptoContext<DCRTPoly>& cc) {
    InArchive ar(is);
    std::string magic;
    uint32_t version = 0;
    bool compact     = false;
    ar(magic, version, compact);
    if (magic != "OpenFHE-SHIP" || version != kFormatVersion)
        OPENFHE_THROW("SHIPDeserializeBootstrapKey: not a SHIP key or unsupported version");
    auto key     = std::make_shared<SHIPBootstrapKey>();
    key->context = cc;
    auto& p      = key->params;
    ar(p.hammingWeight, p.window, p.columnSize, p.muxBase, p.encapsulationModBits, p.realOnly, p.messageBound);
    ar(key->denseTag, key->sparseTag, key->ringDim, key->slots, key->treeDepth, key->leafExtraBits, key->seed);
    if (key->ringDim != cc->GetRingDimension())
        OPENFHE_THROW("SHIPDeserializeBootstrapKey: ring dimension does not match the context");
    const auto params = CkksParams(cc);
    auto& enc         = key->encapsulation;
    uint64_t pPrime = 0, root = 0;
    ar(pPrime, root);
    enc.paramsQ0 = params->GetElementParams()->GetParams()[0];
    enc.paramsP  = std::make_shared<ILNativeParams>(2 * key->ringDim, NativeInteger(pPrime), NativeInteger(root));
    ar(enc.bQ0, enc.aQ0, enc.bP, enc.aP);
    enc.pInvModQ0 = NativeInteger(pPrime).Mod(enc.paramsQ0->GetModulus()).ModInverse(enc.paramsQ0->GetModulus());
    uint32_t rotations = 0;
    ar(key->packing.babyStep, key->packing.inputLevel, rotations);
    for (uint32_t i = 0; i < rotations; ++i) {
        auto rk = LoadRotationKey(ar, compact, cc, key->seed);
        key->packing.rotations.emplace(rk.rotation, std::move(rk));
    }
    key->conjugation = LoadRotationKey(ar, compact, cc, key->seed);
    ar(key->numFactors);
    key->factors.resize(key->numFactors);
    for (uint32_t f = 0; f < key->numFactors; ++f) {
        key->factors[f] = LoadFactor(ar, compact, cc, key->seed);
        ComputeColumnPermutations(*key, key->factors[f].offset, key->factors[f].column.size());
    }
    key->packing.diagonals = PackingDiagonals(cc, key->packing.inputLevel, p.messageBound);
    return key;
}

const std::string& SHIPKeyTag(const SHIPBootstrapKey& key) {
    return key.denseTag;
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
