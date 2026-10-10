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
 * SHIP bootstrapping for CKKS (Cheon, Hanrot, Kim, Stehle, EUROCRYPT 2025,
 * "SHIP: A Shallow and Highly Parallelizable CKKS Bootstrapping Algorithm", ePrint 2025/784).
 *
 * Pipeline (paper Sections 3-5):
 *   S2C (dense matrix, BSGS, one level) at the two lowest moduli
 *   -> dense-to-sparse key switch at modulus q0 * p' (sparse secret encapsulation)
 *   -> Algorithm 1 (half bootstrap) with masks fused into column-method blind rotation over PQ
 *      (Sections 4.1 and 4.4, Algorithm 4) and base-B mux blind rotation (Algorithm 3/5, B-to-1 mux)
 *   -> product tree (ceil(log2(h+1)) levels) and conjugation.
 * Requires CKKS with FIXEDMANUAL scaling, HYBRID key switching and full packing.
 */

#ifndef LBCRYPTO_CRYPTO_CKKSRNS_SHIP_H
#define LBCRYPTO_CRYPTO_CKKSRNS_SHIP_H

#include "ciphertext-fwd.h"
#include "cryptocontext-fwd.h"
#include "key/privatekey-fwd.h"
#include "lattice/lat-hal.h"
#include "lattice/stdlatticeparms.h"

#include <cstdint>
#include <memory>
#include <string>

namespace lbcrypto {

/**
 * @brief Bootstrapping-key parameters of SHIP.
 */
struct SHIPParams {
    /// Hamming weight h of the bottom-modulus sparse ternary secret.
    uint32_t hammingWeight = 31;
    /// The k-th nonzero coefficient index lies in [k*N/h - w, k*N/h + w) (paper Section 5.1).
    /// 0 means unrestricted positions (full-range blind rotation).
    uint32_t window = 175;
    /// theta: size of the column-method part of each blind rotation (paper Table 2).
    uint32_t columnSize = 6;
    /// B: base of the mux-method part of each blind rotation (B-to-1 mux-rotate, paper Section 5.1).
    uint32_t muxBase = 4;
    /// Bit size of the special prime p' of the dense-to-sparse switching key modulo q0 * p'.
    uint32_t encapsulationModBits = 28;
    /// Real-valued messages need one half bootstrap; complex messages need two.
    bool realOnly = true;
};

/**
 * @brief Ring and modulus layout of a SHIP-ready CKKS context (paper Table 2: Base / S2C / Mult / Boot / Aux).
 *
 * Q = q0 (firstModBits) * q_S2C * q_mult^(multLevels) (both ~2^scalingModBits, the CKKS scale)
 *     * boot primes for the ceil(log2(h+1)) product-tree levels, P = auxiliary primes of auxModBits.
 * With bootModBits > scalingModBits the product tree runs at scale 2^bootModBits (better precision);
 * its last (lowest) prime has 2*bootModBits - scalingModBits bits so the output scale is again 2^scalingModBits.
 * gamma = q0 / 2^scalingModBits sets the sine-approximation accuracy.
 */
struct SHIPContextSpec {
    uint32_t ringDim        = 0;
    uint32_t firstModBits   = 0;
    uint32_t scalingModBits = 0;
    uint32_t bootModBits    = 0;  // 0: same as scalingModBits (uniform FIXEDMANUAL chain)
    uint32_t multLevels     = 1;
    uint32_t hammingWeight  = 31;
    uint32_t numLargeDigits = 0;
    uint32_t auxModBits     = 0;
    SecurityLevel securityLevel = HEStd_128_classic;

    /// Counterpart of the paper's LL13 (N = 2^13, about 4.5 bits of precision, one multiplicative level).
    static SHIPContextSpec LL13();
    /// Counterpart of the paper's LL14 (N = 2^14, about 17 bits of precision, one multiplicative level).
    static SHIPContextSpec LL14();
};

class SHIPBootstrapKey;

/// Levels used by the product tree of Algorithm 1: ceil(log2(h + 1)).
uint32_t SHIPProductTreeDepth(uint32_t hammingWeight);

/**
 * @brief Generates a CKKS context for SHIP and checks log2(QP) against the HE standard for a ternary secret.
 * Throws if the requested security level is not met.
 */
CryptoContext<DCRTPoly> GenSHIPCryptoContext(const SHIPContextSpec& spec);

/// log2 of Q*P of a context (for reports).
uint32_t SHIPLogQP(const CryptoContext<DCRTPoly>& cc);

/**
 * @brief Generates the SHIP bootstrapping key for the given (dense) secret key.
 * Samples a fresh bottom sparse secret, which never leaves this function.
 * Also inserts the S2C rotation keys and the conjugation key into the context.
 * EvalMultKeyGen must have been called for the secret key.
 */
std::shared_ptr<SHIPBootstrapKey> SHIPKeyGen(const PrivateKey<DCRTPoly>& privateKey, const SHIPParams& params);

/**
 * @brief SHIP bootstrapping. The input must be fully packed, rescaled (noise scale degree 1) and have at least
 * two RNS limbs; extra limbs are dropped. The output has (number of Q primes - ceil(log2(h+1))) limbs.
 * Message magnitudes must stay small relative to gamma (|m| <= 1 for the preset contexts): the
 * algorithm computes gamma/(2 pi) sin(2 pi m / gamma) (paper Section 1.2).
 */
Ciphertext<DCRTPoly> SHIPBootstrap(ConstCiphertext<DCRTPoly>& ciphertext, const SHIPBootstrapKey& key);

/// Bytes of stored polynomial coefficients in the key (excluding context-level rotation keys).
uint64_t SHIPKeyStoredBytes(const SHIPBootstrapKey& key);

/// Key store used by CryptoContextImpl::EvalSHIPBootstrapKeyGen / EvalSHIPBootstrap (indexed by key tag).
void SHIPInsertBootstrapKey(const std::string& keyTag, std::shared_ptr<SHIPBootstrapKey> key);
std::shared_ptr<SHIPBootstrapKey> SHIPGetBootstrapKey(const std::string& keyTag);
void SHIPClearBootstrapKeys();

}  // namespace lbcrypto

#endif
