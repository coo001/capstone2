//==================================================================================
// BSD 2-Clause License
//
// Copyright (c) 2014-2022, NJIT, Duality Technologies Inc. and other contributors
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

Example for CKKS bootstrapping with full packing

*/

#define PROFILE

#include "openfhe.h"

using namespace lbcrypto;

struct ShipTestPrecom {
    uint32_t shipTreeHeight;
    std::vector<int32_t> shipHMuxOffsets;
    // [level][jIdx][0=keep,1=rot]
    std::vector<std::vector<std::array<Plaintext, 2>>> shipHMuxMasks;
};

void SimpleBootstrapExample();

int main(int argc, char* argv[]) {
    SimpleBootstrapExample();
}

struct HMRKeyTest {
    Ciphertext<DCRTPoly> ctBeta;  // Enc(beta)  (beta in {0,1})
    int32_t rot;                  // rotation amount
};
struct HMRKey {
    Ciphertext<DCRTPoly> k0;  // over RPQ
    Ciphertext<DCRTPoly> k1;  // over RPQ
    int32_t rot;              // rotation amount
};

// HMuxRot_test(key, ct) = Enc(beta) * Rot(ct, rot)
// (beta=0이면 거의 0, beta=1이면 거의 Rot(ct))
static Ciphertext<DCRTPoly> HMuxRot_Test(CryptoContext<DCRTPoly> cc, const HMRKeyTest& key,
                                         ConstCiphertext<DCRTPoly>& ct) {
    auto ctRot = (key.rot == 0) ? ct->Clone() : cc->EvalAtIndex(ct, key.rot);

    // ciphertext-ciphertext mult (level consumed)
    auto out = cc->EvalMult(ctRot, key.ctBeta);

    // 보통 CKKS에서는 mult 후 rescale이 자동/반자동인데,
    // 환경에 따라 필요할 수도 있음 (FLEXIBLEAUTO면 대체로 내부에서 처리됨).
    // out = cc->Rescale(out); // 필요하면 주석 해제

    return out;
}

static Ciphertext<DCRTPoly> HMuxRot(CryptoContext<DCRTPoly> cc, const HMRKey& key, ConstCiphertext<DCRTPoly>& ctQ) {
    auto scheme = cc->GetScheme();
    // ctQ: R_Q^2

    // 1) hoisting digits: ctQ를 여러 rotation/keyswitch에 재사용할 수 있게 준비
    auto digits = cc->EvalFastRotationPrecompute(ctQ);

    // 2) (a', b') = Rot_j(a, b) 를 "ext form"으로 얻기
    //    true => result is in "extended" representation (RPQ)
    Ciphertext<DCRTPoly> ctRotExt;
    if (key.rot != 0) {
        ctRotExt = cc->EvalFastRotationExt(ctQ, key.rot, digits, true);  // RPQ ext
    }
    else {
        ctRotExt = cc->KeySwitchExt(ctQ, true);  // RPQ ext (no rotation)
    }

    // 3) outExt = a' * k0 + b' * k1 (mod PQ), in ext form
    //    여기서 "a'와 b'"는 ctRotExt의 두 요소에 들어있고,
    //    scheme->EvalMultExt는 (extCiphertext, ciphertextOverPQ) 곱을 처리한다고 보면 됨.
    auto t0 = scheme->EvalMultExt(ctRotExt, key.k0);
    auto t1 = scheme->EvalMultExt(ctRotExt, key.k1);
    scheme->EvalAddExtInPlace(t0, t1);

    // 4) PQ -> Q로 다운 (keyswitch-down)
    auto outQ = cc->KeySwitchDown(t0);

    // 5) 마지막에 P^{-1} 스케일링(논문에서 ⌊P^{-1}·⌉ )
    //    OpenFHE CKKS에서는 상수 곱으로 근사 처리
    //    (정확한 P는 key 생성에서 사용한 aux prime product와 일치해야 함)
    //    여기서는 outQ의 스케일이 이미 "P 배"로 들어온 상태라고 가정.
    //    만약 OpenFHE 내부 스케일링 정책에 따라 자동 rescale이 걸리면 이 줄은 조정.
    // outQ = cc->EvalMult(outQ, invP);

    return outQ;
}


// Algorithm 3의 BRotMux를 "Enc(bits)"로 구현한 테스트 버전
// mmkey[j0] = for k: (Enc(j_k), rot=2^k) 와 (Enc(1-j_k), rot=0)
static Ciphertext<DCRTPoly> BRotMux(CryptoContext<DCRTPoly> cc,
                                         const std::vector<std::pair<HMRKeyTest, HMRKeyTest>>& mmkey,
                                         ConstCiphertext<DCRTPoly>& ctIn) {
    Ciphertext<DCRTPoly> ctOut = ctIn->Clone();
    for (size_t k = 0; k < mmkey.size(); ++k) {
        auto ct0 = HMuxRot_Test(cc, mmkey[k].first, ctOut);   // rot by 2^k if bit=1
        auto ct1 = HMuxRot_Test(cc, mmkey[k].second, ctOut);  // rot by 0 if bit=0
        ctOut    = cc->EvalAdd(ct0, ct1);
    }
    return ctOut;
}

void SimpleBootstrapExample() {
    CCParams<CryptoContextCKKSRNS> parameters;
    // A. Specify main parameters
    /*  A1) Secret key distribution
    * The secret key distribution for CKKS should either be SPARSE_TERNARY or UNIFORM_TERNARY.
    * The SPARSE_TERNARY distribution was used in the original CKKS paper,
    * but in this example, we use UNIFORM_TERNARY because this is included in the homomorphic
    * encryption standard.
    */
    SecretKeyDist secretKeyDist = UNIFORM_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);

    /*  A2) Desired security level based on FHE standards.
    * In this example, we use the "NotSet" option, so the example can run more quickly with
    * a smaller ring dimension. Note that this should be used only in
    * non-production environments, or by experts who understand the security
    * implications of their choices. In production-like environments, we recommend using
    * HEStd_128_classic, HEStd_192_classic, or HEStd_256_classic for 128-bit, 192-bit,
    * or 256-bit security, respectively. If you choose one of these as your security level,
    * you do not need to set the ring dimension.
    */
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetRingDim(1 << 12);

    /*  A3) Scaling parameters.
    * By default, we set the modulus sizes and rescaling technique to the following values
    * to obtain a good precision and performance tradeoff. We recommend keeping the parameters
    * below unless you are an FHE expert.
    */
#if NATIVEINT == 128
    ScalingTechnique rescaleTech = FIXEDAUTO;
    usint dcrtBits               = 78;
    usint firstMod               = 89;
#else
    ScalingTechnique rescaleTech = FLEXIBLEAUTO;
    usint dcrtBits               = 59;
    usint firstMod               = 60;
#endif

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetFirstModSize(firstMod);

    /*  A4) Multiplicative depth.
    * The goal of bootstrapping is to increase the number of available levels we have, or in other words,
    * to dynamically increase the multiplicative depth. However, the bootstrapping procedure itself
    * needs to consume a few levels to run. We compute the number of bootstrapping levels required
    * using GetBootstrapDepth, and add it to levelsAvailableAfterBootstrap to set our initial multiplicative
    * depth. We recommend using the input parameters below to get started.
    */
    std::vector<uint32_t> levelBudget = {4, 4};

    // Note that the actual number of levels avalailable after bootstrapping before next bootstrapping
    // will be levelsAvailableAfterBootstrap - 1 because an additional level
    // is used for scaling the ciphertext before next bootstrapping (in 64-bit CKKS bootstrapping)
    uint32_t levelsAvailableAfterBootstrap = 40;//10
    usint depth = levelsAvailableAfterBootstrap + FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDist);
    parameters.SetMultiplicativeDepth(depth);

    CryptoContext<DCRTPoly> cryptoContext = GenCryptoContext(parameters);

    cryptoContext->Enable(PKE);
    cryptoContext->Enable(KEYSWITCH);
    cryptoContext->Enable(LEVELEDSHE);
    cryptoContext->Enable(ADVANCEDSHE);
    cryptoContext->Enable(FHE);

    usint ringDim = cryptoContext->GetRingDimension();
    // This is the maximum number of slots that can be used for full packing.
    usint numSlots = ringDim / 2;
    std::cout << "CKKS scheme is using ring dimension " << ringDim << std::endl << std::endl;

    auto keyPair = cryptoContext->KeyGen();
    cryptoContext->EvalBootstrapSetup(levelBudget);

    cryptoContext->EvalMultKeyGen(keyPair.secretKey);
    cryptoContext->EvalBootstrapKeyGen(keyPair.secretKey, numSlots);
    // 회전에 필요한 인덱스들 수집
    std::vector<int32_t> rotIndices;

    // BRotMux에서 쓰는 2^k 회전들
    for (uint32_t k = 0; (1u << k) < numSlots; ++k) {
        rotIndices.push_back(1 << k);
    }

    // 참값 비교용 j0 회전도 추가
    rotIndices.push_back(3);

    // 회전 키 생성
    cryptoContext->EvalAtIndexKeyGen(keyPair.secretKey, rotIndices);


    //여기부터
    std::vector<double> x            = {0.25, 0.5, 0.75, 1.0, 2.0, 3.0, 4.0, 5.0};
    size_t encodedLength             = x.size();

    // We start with a depleted ciphertext that has used up all of its levels.
    Plaintext ptxt = cryptoContext->MakeCKKSPackedPlaintext(x);

    ptxt->SetLength(encodedLength);
    std::cout << "Input: " << ptxt << std::endl;

    Ciphertext<DCRTPoly> ciph = cryptoContext->Encrypt(keyPair.publicKey, ptxt);

    // -------------------------
    // BRotMux 테스트 파라미터
    // -------------------------
    uint32_t j0    = 3;  // 테스트할 rotation index (0<=j0<slots)
    uint32_t nBits = 0;
    while ((1u << nBits) < numSlots)
        nBits++;  // log2(slots) 정도

    // Enc(1) 준비: (모든 슬롯에 1) -> 암호화
    std::vector<std::complex<double>> ones(numSlots, 1.0);
    Plaintext ptOnes = cryptoContext->MakeCKKSPackedPlaintext(ones);
    auto ctOnes      = cryptoContext->Encrypt(keyPair.publicKey, ptOnes);

    // mmkey: k마다 (Enc(j_k), rot=2^k) 와 (Enc(1-j_k), rot=0)
    std::vector<std::pair<HMRKeyTest, HMRKeyTest>> mmkey;
    mmkey.reserve(nBits);

    for (uint32_t k = 0; k < nBits; ++k) {
        uint32_t bit = (j0 >> k) & 1u;

        // Enc(bit)
        std::vector<std::complex<double>> vb(numSlots, (double)bit);
        Plaintext ptBit = cryptoContext->MakeCKKSPackedPlaintext(vb);
        auto ctBit      = cryptoContext->Encrypt(keyPair.publicKey, ptBit);

        // Enc(1-bit) = Enc(1) - Enc(bit)
        auto ctOneMinusBit = cryptoContext->EvalSub(ctOnes, ctBit);

        HMRKeyTest keyIf1{ctBit, (int32_t)(1u << k)};  // rotate by 2^k if bit=1
        HMRKeyTest keyIf0{ctOneMinusBit, 0};           // rotate by 0 if bit=0

        mmkey.push_back({keyIf1, keyIf0});
    }

    auto ctRef = cryptoContext->EvalAtIndex(ciph, (int32_t)j0);
    auto ctMux = BRotMux_Test(cryptoContext, mmkey, ciph);

    // decrypt해서 비교
    Plaintext pRef, pMux;
    cryptoContext->Decrypt(keyPair.secretKey, ctRef, &pRef);
    cryptoContext->Decrypt(keyPair.secretKey, ctMux, &pMux);
    pRef->SetLength(encodedLength);
    pMux->SetLength(encodedLength);

    auto vRef = pRef->GetCKKSPackedValue();
    auto vMux = pMux->GetCKKSPackedValue();

    std::cout << "=== BRotMux_Test j0=" << j0 << " compare ===\n";
    for (size_t i = 0; i < encodedLength; ++i) {
        std::cout << i << " ref=" << vRef[i] << " mux=" << vMux[i] << " diff=" << (vRef[i] - vMux[i]) << "\n";
    }

    // 여기서 return; 해버리면 부트스트랩 전체는 안 돌고 ShipBlindRot만 검증 가능
    return;

    /*
    std::vector<double> x = {0.25, 0.5, 0.75, 1.0, 2.0, 3.0, 4.0, 5.0};
    size_t encodedLength  = x.size();

    // We start with a depleted ciphertext that has used up all of its levels.
    Plaintext ptxt = cryptoContext->MakeCKKSPackedPlaintext(x, 1, depth - 1);

    ptxt->SetLength(encodedLength);
    std::cout << "Input: " << ptxt << std::endl;

    Ciphertext<DCRTPoly> ciph = cryptoContext->Encrypt(keyPair.publicKey, ptxt);

    std::cout << "Initial number of levels remaining: " << depth - ciph->GetLevel() << std::endl;

    */
    
    // Perform the bootstrapping operation. The goal is to increase the number of levels remaining
    // for HE computation.
    //auto inScale         = ciph->GetScalingFactor();
    auto ciphertextAfter = cryptoContext->EvalBootstrap(ciph, 1, 0, keyPair.publicKey);

    std::cout << "Number of levels remaining after bootstrapping: "
              << depth - ciphertextAfter->GetLevel() - (ciphertextAfter->GetNoiseScaleDeg() - 1) << std::endl
              << std::endl;
    
    Plaintext result;
    std::cout << "level" << ciphertextAfter->GetLevel() << std::endl;
    std::cout << "tower" << ciphertextAfter->GetElements()[0].GetNumOfElements() << std::endl;
    std::cout << "ct keytag: " << ciphertextAfter->GetKeyTag() << "\n";
    std::cout << "sk keytag: " << keyPair.secretKey->GetKeyTag() << "\n";
    std::cout << "scale: " << ciphertextAfter->GetScalingFactor() << "\n";
    std::cout << "noiseScaleDeg: " << ciphertextAfter->GetNoiseScaleDeg() << "\n";
    std::cout << "level/towers: " << ciphertextAfter->GetLevel() << "/"
              << ciphertextAfter->GetElements()[0].GetNumOfElements() << "\n";
    auto mod = ciphertextAfter->GetElements()[0].GetModulus();
    std::cout << "log2(Q_remaining) ~= " << mod.GetMSB() << "\n";
    std::cout << "noise" <<  ciphertextAfter->GetNoiseScaleDeg() << "\n";

    cryptoContext->Decrypt(keyPair.secretKey, ciphertextAfter, &result);
    result->SetLength(encodedLength);
    auto v = result->GetCKKSPackedValue();
    for (size_t i = 0; i < 8; ++i) {
        std::cout << i << " : " << v[i] << " (isnan=" << std::isnan(v[i].real())
                  << ", isinf=" << std::isinf(v[i].real()) << ")\n";
    }
    std::cout << "scale = " << result->GetScalingFactor() << std::endl;
    std::cout << "Output after bootstrapping \n\t" << result << std::endl;
}