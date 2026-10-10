// Client/server deployment of SHIP bootstrapping (toy, non-secure parameters):
//  - the client serializes the context, keys and a ciphertext; all in-memory state is then released,
//  - the server deserializes them, bootstraps and computes, and returns the ciphertext,
//  - the client decrypts.
// Also checks compact (seeded) versus full key sizes, disk-backed factor keys and the message bound.
#include "openfhe.h"

#include "ciphertext-ser.h"
#include "cryptocontext-ser.h"
#include "key/key-ser.h"
#include "scheme/ckksrns/ckksrns-ser.h"

#include <filesystem>
#include <iostream>
#include <random>
#include <sstream>

using namespace lbcrypto;
using C = std::complex<double>;

static SHIPContextSpec ToySpec() {
    SHIPContextSpec s;
    s.ringDim        = 1024;
    s.firstModBits   = 60;
    s.scalingModBits = 50;
    s.multLevels     = 1;
    s.hammingWeight  = 8;
    s.numLargeDigits = 3;
    s.auxModBits     = 60;
    s.securityLevel  = HEStd_NotSet;
    return s;
}

static SHIPParams ToyParams() {
    SHIPParams p;
    p.hammingWeight        = 8;
    p.window               = 20;
    p.columnSize           = 6;
    p.muxBase              = 4;
    p.encapsulationModBits = 60;
    return p;
}

static double MaxError(const CryptoContext<DCRTPoly>& cc, const PrivateKey<DCRTPoly>& sk, const Ciphertext<DCRTPoly>& ct,
                       const std::vector<double>& expected) {
    Plaintext pt;
    cc->Decrypt(sk, ct, &pt);
    pt->SetLength(expected.size());
    double e = 0;
    for (size_t i = 0; i < expected.size(); ++i)
        e = std::max(e, std::abs(pt->GetCKKSPackedValue()[i].real() - expected[i]));
    return e;
}

static void Require(bool ok, const std::string& what) {
    if (!ok)
        throw std::runtime_error(what);
}

static void ReleaseEverything() {
    CryptoContextImpl<DCRTPoly>::ClearEvalMultKeys();
    CryptoContextImpl<DCRTPoly>::ClearEvalAutomorphismKeys();
    CryptoContextImpl<DCRTPoly>::ClearEvalSHIPBootstrapKeys();
    CryptoContextFactory<DCRTPoly>::ReleaseAllContexts();
}

int main() {
    try {
        std::mt19937 rng(20261011);
        // ---------------- client ----------------
        std::stringstream ccBlob, pkBlob, skBlob, multBlob, shipBlob, shipFullBlob, ctBlob;
        std::vector<double> x;
        {
            auto cc          = GenSHIPCryptoContext(ToySpec());
            const uint32_t S = cc->GetRingDimension() / 2;
            const uint32_t Q = cc->GetElementParams()->GetParams().size();
            auto kp          = cc->KeyGen();
            cc->EvalMultKeyGen(kp.secretKey);
            cc->EvalSHIPBootstrapKeyGen(kp.secretKey, ToyParams());
            std::uniform_real_distribution<double> u(-1, 1);
            x.resize(S);
            for (auto& v : x)
                v = u(rng);
            auto ct = cc->Encrypt(kp.publicKey, cc->MakeCKKSPackedPlaintext(x, 1, Q - 2, nullptr, S));
            Serial::Serialize(cc, ccBlob, SerType::BINARY);
            Serial::Serialize(kp.publicKey, pkBlob, SerType::BINARY);
            Serial::Serialize(kp.secretKey, skBlob, SerType::BINARY);
            CryptoContextImpl<DCRTPoly>::SerializeEvalMultKey(multBlob, SerType::BINARY, kp.secretKey->GetKeyTag());
            CryptoContextImpl<DCRTPoly>::SerializeEvalSHIPBootstrapKey(shipBlob, kp.secretKey->GetKeyTag(), true);
            CryptoContextImpl<DCRTPoly>::SerializeEvalSHIPBootstrapKey(shipFullBlob, kp.secretKey->GetKeyTag(), false);
            Serial::Serialize(ct, ctBlob, SerType::BINARY);
            const auto key = SHIPGetBootstrapKey(kp.secretKey->GetKeyTag());
            std::cout << "ship_key_in_memory_bytes=" << SHIPKeyStoredBytes(*key)
                      << " serialized_compact_bytes=" << shipBlob.str().size()
                      << " serialized_full_bytes=" << shipFullBlob.str().size()
                      << " ratio=" << double(shipBlob.str().size()) / shipFullBlob.str().size() << std::endl;
            Require(shipBlob.str().size() < 0.6 * shipFullBlob.str().size(), "compact serialization is not smaller");
        }
        ReleaseEverything();

        // ---------------- server ----------------
        std::stringstream resultBlob;
        {
            CryptoContext<DCRTPoly> cc;
            Serial::Deserialize(cc, ccBlob, SerType::BINARY);
            Require(SHIPLogQP(cc) > 0 && cc->GetRingDimension() == 1024, "context deserialization");
            CryptoContextImpl<DCRTPoly>::DeserializeEvalMultKey(multBlob, SerType::BINARY);
            cc->DeserializeEvalSHIPBootstrapKey(shipBlob);
            Ciphertext<DCRTPoly> ct;
            Serial::Deserialize(ct, ctBlob, SerType::BINARY);
            auto refreshed = cc->EvalSHIPBootstrap(ct);
            auto squared   = cc->EvalMult(refreshed, refreshed);
            cc->RescaleInPlace(squared);
            auto again = cc->EvalSHIPBootstrap(squared);
            Serial::Serialize(again, resultBlob, SerType::BINARY);
            std::cout << "server: bootstrapped, squared, bootstrapped again; output limbs="
                      << again->GetElements()[0].GetNumOfElements() << std::endl;
        }

        // ---------------- client ----------------
        {
            CryptoContext<DCRTPoly> cc;
            std::stringstream ccAgain(ccBlob.str());
            Serial::Deserialize(cc, ccAgain, SerType::BINARY);
            PrivateKey<DCRTPoly> sk;
            Serial::Deserialize(sk, skBlob, SerType::BINARY);
            Ciphertext<DCRTPoly> result;
            Serial::Deserialize(result, resultBlob, SerType::BINARY);
            auto expected = x;
            for (auto& v : expected)
                v *= v;
            const double error = MaxError(cc, sk, result, expected);
            std::cout << "client: max_error=" << error << " precision_bits=" << -std::log2(error) << std::endl;
            Require(error < 2e-4, "deployment round trip failed");
        }
        ReleaseEverything();

        // ---------------- disk-backed factor keys and full-key serialization ----------------
        {
            auto cc          = GenSHIPCryptoContext(ToySpec());
            const uint32_t S = cc->GetRingDimension() / 2;
            const uint32_t Q = cc->GetElementParams()->GetParams().size();
            auto kp          = cc->KeyGen();
            cc->EvalMultKeyGen(kp.secretKey);
            const auto dir = (std::filesystem::temp_directory_path() / "ship-factor-check").string();
            std::filesystem::remove_all(dir);
            auto disk = SHIPKeyGen(kp.secretKey, ToyParams(), dir);
            size_t files = 0;
            for (const auto& entry : std::filesystem::directory_iterator(dir))
                files += entry.is_regular_file();
            auto ct          = cc->Encrypt(kp.publicKey, cc->MakeCKKSPackedPlaintext(x, 1, Q - 2, nullptr, S));
            const double e1  = MaxError(cc, kp.secretKey, SHIPBootstrap(ct, *disk), x);
            std::stringstream blob;
            SHIPSerializeBootstrapKey(blob, *disk, false);
            auto loaded      = SHIPDeserializeBootstrapKey(blob, cc);
            const double e2  = MaxError(cc, kp.secretKey, SHIPBootstrap(ct, *loaded), x);
            std::cout << "disk_factor_files=" << files << " disk_key_in_memory_bytes=" << SHIPKeyStoredBytes(*disk)
                      << " disk_error=" << e1 << " reloaded_full_error=" << e2 << std::endl;
            Require(files == 8 && e1 < 1e-4 && e2 < 1e-4, "disk-backed key failed");
            std::filesystem::remove_all(dir);
        }
        ReleaseEverything();

        // ---------------- message bound ----------------
        {
            auto cc          = GenSHIPCryptoContext(ToySpec());
            const uint32_t S = cc->GetRingDimension() / 2;
            const uint32_t Q = cc->GetElementParams()->GetParams().size();
            auto kp          = cc->KeyGen();
            cc->EvalMultKeyGen(kp.secretKey);
            std::uniform_real_distribution<double> u(-64, 64);
            std::vector<double> big(S);
            for (auto& v : big)
                v = u(rng);
            auto ct = cc->Encrypt(kp.publicKey, cc->MakeCKKSPackedPlaintext(big, 1, Q - 2, nullptr, S));
            auto p1 = ToyParams();
            auto p64 = ToyParams();
            p64.messageBound = 64;
            auto k1  = SHIPKeyGen(kp.secretKey, p1);
            auto k64 = SHIPKeyGen(kp.secretKey, p64);
            const double e1  = MaxError(cc, kp.secretKey, SHIPBootstrap(ct, *k1), big);
            const double e64 = MaxError(cc, kp.secretKey, SHIPBootstrap(ct, *k64), big);
            std::cout << "inputs in [-64,64]: messageBound=1 max_error=" << e1 << ", messageBound=64 max_error=" << e64
                      << " (relative " << e64 / 64 << ")" << std::endl;
            Require(e1 > 0.5 && e64 < 0.01, "message bound did not widen the input range");
        }
        ReleaseEverything();
        std::cout << "ALL SHIP DEPLOYMENT CHECKS PASSED (toy parameters, HEStd_NotSet)" << std::endl;
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "FAIL " << e.what() << std::endl;
        return 1;
    }
}
