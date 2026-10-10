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

/*
 * Unit tests for SHIP bootstrapping (scheme/ckksrns/ckksrns-ship.h). Toy, non-secure ring dimension
 * for speed; the 128-bit parameter sets are covered by examples/pke/ship-paper-bench.cpp.
 */

#include "cryptocontext.h"
#include "gen-cryptocontext.h"
#include "scheme/ckksrns/ckksrns-ship.h"
#include "scheme/ckksrns/gen-cryptocontext-ckksrns.h"

#include "gtest/gtest.h"

#include <cmath>
#include <filesystem>
#include <sstream>
#include <complex>
#include <random>
#include <vector>

using namespace lbcrypto;

namespace {
class UTCKKSRNS_SHIP : public ::testing::Test {
protected:
    void SetUp() {}

    void TearDown() {
        CryptoContextImpl<DCRTPoly>::ClearEvalSHIPBootstrapKeys();
        CryptoContextImpl<DCRTPoly>::ClearEvalMultKeys();
        CryptoContextImpl<DCRTPoly>::ClearEvalAutomorphismKeys();
        CryptoContextFactory<DCRTPoly>::ReleaseAllContexts();
    }
};

SHIPContextSpec ToySpec(uint32_t h, uint32_t bootBits = 0) {
    SHIPContextSpec s;
    s.ringDim        = 1024;
    s.firstModBits   = 60;
    s.scalingModBits = 50;  // gamma = 2^10: sine approximation error about 6e-6 for |m| <= 1
    s.bootModBits    = bootBits;
    s.multLevels     = 1;
    s.hammingWeight  = h;
    s.numLargeDigits = 3;
    s.auxModBits     = 60;
    s.securityLevel  = HEStd_NotSet;
    return s;
}

SHIPParams ToyParams() {
    SHIPParams p;
    p.hammingWeight        = 8;
    p.window               = 20;
    p.columnSize           = 6;
    p.muxBase              = 4;
    p.encapsulationModBits = 60;
    p.realOnly             = true;
    return p;
}

double MaxError(const CryptoContext<DCRTPoly>& cc, const PrivateKey<DCRTPoly>& sk, const Ciphertext<DCRTPoly>& ct,
                const std::vector<std::complex<double>>& expected) {
    Plaintext pt;
    cc->Decrypt(sk, ct, &pt);
    pt->SetLength(expected.size());
    double error = 0;
    for (size_t i = 0; i < expected.size(); ++i)
        error = std::max(error, std::abs(pt->GetCKKSPackedValue()[i] - expected[i]));
    return error;
}

std::vector<std::complex<double>> Random(uint32_t n, bool real, uint32_t seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(-1, 1);
    std::vector<std::complex<double>> x(n);
    for (auto& v : x)
        v = real ? std::complex<double>(u(rng), 0) : std::complex<double>(u(rng), u(rng));
    return x;
}

uint32_t Limbs(const Ciphertext<DCRTPoly>& ct) {
    return ct->GetElements()[0].GetNumOfElements();
}

void CheckRoundTrip(const SHIPContextSpec& spec, const SHIPParams& params) {
    auto cc          = GenSHIPCryptoContext(spec);
    const uint32_t S = cc->GetRingDimension() / 2;
    const uint32_t Q = cc->GetElementParams()->GetParams().size();
    auto kp          = cc->KeyGen();
    cc->EvalMultKeyGen(kp.secretKey);
    cc->EvalSHIPBootstrapKeyGen(kp.secretKey, params);
    auto x  = Random(S, params.realOnly, 7);
    auto ct = cc->Encrypt(kp.publicKey, cc->MakeCKKSPackedPlaintext(x, 1, Q - 2, nullptr, S));
    auto out = cc->EvalSHIPBootstrap(ct);
    EXPECT_LT(MaxError(cc, kp.secretKey, out, x), 1e-4);
    EXPECT_EQ(Limbs(out), Q - SHIPProductTreeDepth(params.hammingWeight));
    EXPECT_EQ(out->GetNoiseScaleDeg(), 1u);
    // The restored level is usable: square, rescale, bootstrap again.
    auto sq = cc->EvalMult(out, out);
    cc->RescaleInPlace(sq);
    for (auto& v : x)
        v *= v;
    EXPECT_LT(MaxError(cc, kp.secretKey, cc->EvalSHIPBootstrap(sq), x), 2e-4);
}
}  // anonymous namespace

TEST_F(UTCKKSRNS_SHIP, RealWindowedColumnAndBase4Mux) {
    CheckRoundTrip(ToySpec(8), ToyParams());
}

TEST_F(UTCKKSRNS_SHIP, ComplexUnrestrictedBinaryMux) {
    auto p       = ToyParams();
    p.window     = 0;
    p.columnSize = 1;
    p.muxBase    = 2;
    p.realOnly   = false;
    CheckRoundTrip(ToySpec(8), p);
}

TEST_F(UTCKKSRNS_SHIP, RoleSpecificBootPrimes) {
    // Larger product-tree primes; the output scale must return to 2^scalingModBits.
    auto spec           = ToySpec(8, 52);  // last boot prime: 2 * 52 - 45 = 59 bits
    spec.scalingModBits = 45;
    spec.firstModBits   = 55;
    CheckRoundTrip(spec, ToyParams());
}

TEST_F(UTCKKSRNS_SHIP, HigherLevelInputIsLevelReduced) {
    auto cc          = GenSHIPCryptoContext(ToySpec(8));
    const uint32_t S = cc->GetRingDimension() / 2;
    auto kp          = cc->KeyGen();
    cc->EvalMultKeyGen(kp.secretKey);
    cc->EvalSHIPBootstrapKeyGen(kp.secretKey, ToyParams());
    auto x  = Random(S, true, 11);
    auto ct = cc->Encrypt(kp.publicKey, cc->MakeCKKSPackedPlaintext(x, 1, 0, nullptr, S));
    EXPECT_LT(MaxError(cc, kp.secretKey, cc->EvalSHIPBootstrap(ct), x), 1e-4);
}

TEST_F(UTCKKSRNS_SHIP, RejectsInvalidInputsAndKeys) {
    auto cc          = GenSHIPCryptoContext(ToySpec(8));
    const uint32_t S = cc->GetRingDimension() / 2;
    const uint32_t Q = cc->GetElementParams()->GetParams().size();
    auto kp          = cc->KeyGen();
    EXPECT_THROW(cc->EvalSHIPBootstrapKeyGen(kp.secretKey, ToyParams()), OpenFHEException);  // no EvalMultKeyGen
    cc->EvalMultKeyGen(kp.secretKey);
    auto bad = ToyParams();
    bad.muxBase = 1;
    EXPECT_THROW(cc->EvalSHIPBootstrapKeyGen(kp.secretKey, bad), OpenFHEException);
    auto x  = Random(S, true, 3);
    auto ct = cc->Encrypt(kp.publicKey, cc->MakeCKKSPackedPlaintext(x, 1, Q - 2, nullptr, S));
    EXPECT_THROW(cc->EvalSHIPBootstrap(ct), OpenFHEException);  // no SHIP key yet
    cc->EvalSHIPBootstrapKeyGen(kp.secretKey, ToyParams());
    auto oneLimb = ct->Clone();
    cc->LevelReduceInPlace(oneLimb, nullptr, 1);
    EXPECT_THROW(cc->EvalSHIPBootstrap(oneLimb), OpenFHEException);
    auto degreeTwo = cc->EvalMult(ct, ct);
    EXPECT_THROW(cc->EvalSHIPBootstrap(degreeTwo), OpenFHEException);
}

TEST_F(UTCKKSRNS_SHIP, CompactSerializationRoundTrip) {
    auto cc          = GenSHIPCryptoContext(ToySpec(8));
    const uint32_t S = cc->GetRingDimension() / 2;
    const uint32_t Q = cc->GetElementParams()->GetParams().size();
    auto kp          = cc->KeyGen();
    cc->EvalMultKeyGen(kp.secretKey);
    cc->EvalSHIPBootstrapKeyGen(kp.secretKey, ToyParams());
    std::stringstream compact, full;
    CryptoContextImpl<DCRTPoly>::SerializeEvalSHIPBootstrapKey(compact, kp.secretKey->GetKeyTag(), true);
    CryptoContextImpl<DCRTPoly>::SerializeEvalSHIPBootstrapKey(full, kp.secretKey->GetKeyTag(), false);
    EXPECT_LT(compact.str().size(), 0.6 * full.str().size());
    CryptoContextImpl<DCRTPoly>::ClearEvalSHIPBootstrapKeys();
    cc->DeserializeEvalSHIPBootstrapKey(compact);
    auto x  = Random(S, true, 5);
    auto ct = cc->Encrypt(kp.publicKey, cc->MakeCKKSPackedPlaintext(x, 1, Q - 2, nullptr, S));
    EXPECT_LT(MaxError(cc, kp.secretKey, cc->EvalSHIPBootstrap(ct), x), 1e-4);
}

TEST_F(UTCKKSRNS_SHIP, DiskBackedFactorKeys) {
    auto cc          = GenSHIPCryptoContext(ToySpec(8));
    const uint32_t S = cc->GetRingDimension() / 2;
    const uint32_t Q = cc->GetElementParams()->GetParams().size();
    auto kp          = cc->KeyGen();
    cc->EvalMultKeyGen(kp.secretKey);
    const auto dir = (std::filesystem::temp_directory_path() / "ship-unittest-factors").string();
    std::filesystem::remove_all(dir);
    cc->EvalSHIPBootstrapKeyGen(kp.secretKey, ToyParams(), dir);
    EXPECT_TRUE(std::filesystem::exists(dir + "/factor-0.ship"));
    auto x  = Random(S, true, 9);
    auto ct = cc->Encrypt(kp.publicKey, cc->MakeCKKSPackedPlaintext(x, 1, Q - 2, nullptr, S));
    EXPECT_LT(MaxError(cc, kp.secretKey, cc->EvalSHIPBootstrap(ct), x), 1e-4);
    std::filesystem::remove_all(dir);
}

TEST_F(UTCKKSRNS_SHIP, MessageBoundWidensInputRange) {
    auto cc          = GenSHIPCryptoContext(ToySpec(8));
    const uint32_t S = cc->GetRingDimension() / 2;
    const uint32_t Q = cc->GetElementParams()->GetParams().size();
    auto kp          = cc->KeyGen();
    cc->EvalMultKeyGen(kp.secretKey);
    auto p         = ToyParams();
    p.messageBound = 64;
    cc->EvalSHIPBootstrapKeyGen(kp.secretKey, p);
    auto x = Random(S, true, 13);
    for (auto& v : x)
        v *= 64.0;
    auto ct = cc->Encrypt(kp.publicKey, cc->MakeCKKSPackedPlaintext(x, 1, Q - 2, nullptr, S));
    EXPECT_LT(MaxError(cc, kp.secretKey, cc->EvalSHIPBootstrap(ct), x), 64 * 1e-4);
}

TEST_F(UTCKKSRNS_SHIP, PresetContextsMeetHEStandard) {
    for (const auto& spec : {SHIPContextSpec::LL13(), SHIPContextSpec::LL14(), SHIPContextSpec::HT14(),
                             SHIPContextSpec::HT15()}) {
        auto cc = GenSHIPCryptoContext(spec);
        EXPECT_LE(SHIPLogQP(cc), StdLatticeParm::FindMaxQ(HEStd_ternary, HEStd_128_classic, spec.ringDim));
        EXPECT_EQ(cc->GetRingDimension(), spec.ringDim);
        EXPECT_EQ(cc->GetElementParams()->GetParams().size(),
                  2 + spec.multLevels + SHIPProductTreeDepth(spec.hammingWeight));
    }
    auto tooLarge           = SHIPContextSpec::LL13();
    tooLarge.scalingModBits = 30;  // exceeds the 128-bit budget at N = 2^13
    tooLarge.bootModBits    = 34;
    tooLarge.firstModBits   = 35;
    EXPECT_THROW(GenSHIPCryptoContext(tooLarge), OpenFHEException);
}
