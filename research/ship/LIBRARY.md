# SHIP 부트스트래핑의 OpenFHE 라이브러리 구현

2026-10-11. 원논문은 Cheon, Hanrot, Kim, Stehlé, "SHIP: A Shallow and Highly Parallelizable CKKS Bootstrapping Algorithm", EUROCRYPT 2025([ePrint 2025/784](https://eprint.iacr.org/2025/784.pdf), 33쪽 판)이다.
페이지 번호는 ePrint 판 기준이다.

## 요약

- SHIP을 OpenFHE 라이브러리 기능으로 구현했다. 예제 헤더가 아니라 `src/pke/lib`에 들어가며 `CryptoContext` API로 호출한다.
- 논문의 알고리즘 구성(Algorithm 1, Lemma 1 마스크, column과 mux를 결합한 블라인드 회전, base-4 B-to-1 mux, 균등 간격 sparse 키, `P·q0` sparse-secret encapsulation, 대각선 S2C, Table 2의 역할별 소수)을 구현했다.
- 논문의 네 파라미터 세트(LL13, LL14, HT14, HT15)를 모두 128비트 HE 표준 안에서 실행했다. HT14·HT15는 키가 각각 약 31GB, 70GB라 키를 디스크에 두고 인자별로 읽는 방식으로 48GB Mac에서 돌렸다.
- 정밀도는 LL14 **16.0–16.5비트**(논문 16.9), HT15 **18.3–18.4비트**(논문 18.1), HT14 **4.0비트**(논문 4.55), LL13 **1.9–2.5비트**(논문 4.45)였다. 부트스트래핑 뒤 남는 곱셈 레벨은 논문과 같다(LL 1, HT 9).
- 배포에 필요한 기능을 추가했다. 키 직렬화(균등 난수 성분을 시드로 저장해 크기 1/2), 디스크 기반 키, 입력 범위 지정(`messageBound`)이다.
- lattice estimator를 직접 실행했다. dense 키는 네 세트 모두 128비트 이상이었다. **sparse 키(h = 31)는 현재 estimator의 혼합 공격 추정에서 128비트에 못 미쳤다**(아래 "보안 검증").
- OpenFHE 기존 `EvalBootstrap`을 SHIP과 같은 sparse-secret encapsulation(`SPARSE_ENCAPSULATED`)으로 조정해 조용한 상태에서 같은 조건(128비트, 남는 곱셈 레벨 1)으로 비교했다. SHIP LL14는 지연 시간이 기존 방식 `N = 2^16`보다 1스레드 2.0배, 15스레드 4.2배 짧고 메모리는 약 절반이었다. 정밀도는 기존 방식이 약 2비트 높았고, 슬롯 수까지 반영한 처리량도 기존 방식이 1.1–2.3배 높았다. SHIP의 이점은 **지연 시간, 작은 N, 메모리**이다.
- 작은 CKKS 스케일에서 정밀도를 위해 OpenFHE의 hybrid key switching ModDown에 정확한 반올림을 추가했다. `WITH_REDUCED_NOISE` 빌드 옵션에서만 동작하며, OpenFHE 전체 단위 테스트(core, binfhe, pke 1,882개)가 통과했다.

## 사용법

```cpp
#include "openfhe.h"
using namespace lbcrypto;

// 클라이언트: 컨텍스트, 키, SHIP 부트스트래핑 키
auto spec = SHIPContextSpec::LL14();        // LL13, LL14, HT14, HT15
auto cc   = GenSHIPCryptoContext(spec);     // HE 표준 128비트 검사 포함
auto kp   = cc->KeyGen();
cc->EvalMultKeyGen(kp.secretKey);
auto params = SHIPParams::Recommended(spec); // w, θ, p'를 논문 방식으로 선택
params.messageBound = 1.0;                   // 입력 크기 상한 K (아래 참고)
cc->EvalSHIPBootstrapKeyGen(kp.secretKey, params);         // 메모리에 키 보관
// cc->EvalSHIPBootstrapKeyGen(kp.secretKey, params, dir); // 인자별 키를 dir에 저장(HT14·HT15)

std::ofstream os("ship.key", std::ios::binary);
CryptoContextImpl<DCRTPoly>::SerializeEvalSHIPBootstrapKey(os, kp.secretKey->GetKeyTag());  // compact

// 서버: 같은 컨텍스트에서 키를 읽고 부트스트래핑
std::ifstream is("ship.key", std::ios::binary);
cc->DeserializeEvalSHIPBootstrapKey(is);
auto out = cc->EvalSHIPBootstrap(ct);       // ceil(log2(h + 1))레벨을 체인 위쪽에서 사용
```

- 입력은 full packing, rescale을 마친 암호문이며 limb가 2개 이상이어야 한다. 더 많으면 S2C 레벨까지 내린다.
- **입력 범위.** SHIP은 `γ/(2π)·sin(2πm/γ)`를 계산하므로(논문 §1.2) `|m|`이 작아야 한다. `messageBound = K`이면 S2C에서 메시지를 K로 나누고, 곱셈 트리의 공개 평문 `pt0`에 K를 곱해 되돌린다. 추가 레벨은 없고 정밀도가 약 `log2 K`비트 줄어든다.
- 출력은 곱셈 트리 깊이만큼 위쪽 레벨을 쓴 뒤의 암호문이다. LL13·LL14에서는 곱셈 레벨 1개, HT14·HT15에서는 9개가 남는다.
- `params.realOnly = false`이면 복소수 메시지를 위해 half-bootstrap을 두 번 수행한다.
- **직렬화.** `compact = true`(기본)는 키 전환 키의 균등 난수 성분 `a`를 저장하지 않고 시드만 저장해 읽을 때 다시 만든다. S2C 대각선은 읽을 때 다시 계산한다. 디스크 기반 키는 인자별 키 파일까지 포함해 쓴다.
- **디스크 기반 키.** `factorDirectory`를 주면 키의 대부분(계수별 column·mux 키)을 생성 중에 파일로 쓰고, 부트스트래핑 때 인자 단위로 읽는다. 메모리에는 S2C·켤레·encapsulation 키만 남는다.
- `ClearEvalSHIPBootstrapKeys()`로 키를 지운다.

## 파일

| 파일 | 내용 |
| --- | --- |
| `src/pke/include/scheme/ckksrns/ckksrns-ship.h` | 공개 API: `SHIPParams`, `SHIPContextSpec`, `GenSHIPCryptoContext`, `SHIPKeyGen`, `SHIPBootstrap`, 직렬화, `SHIPEstimateKeyBytes` |
| `src/pke/lib/scheme/ckksrns/ckksrns-ship.cpp` | 구현 |
| `src/pke/include/cryptocontext.h` | `EvalSHIPBootstrapKeyGen`, `EvalSHIPBootstrap`, `SerializeEvalSHIPBootstrapKey`, `DeserializeEvalSHIPBootstrapKey`, `ClearEvalSHIPBootstrapKeys` |
| `src/core/include/lattice/hal/default/dcrtpoly-impl.h` | `ApproxModDown`의 정확한 반올림 (`WITH_REDUCED_NOISE`일 때) |
| `src/pke/unittest/utckksrns/UnitTestSHIP.cpp` | gtest 9개 (직렬화 왕복, 디스크 키, 입력 범위, 프리셋 보안 검사 포함) |
| `src/pke/examples/ship-library-checks.cpp` | 기능 검사 (toy 파라미터) |
| `src/pke/examples/ship-deployment-checks.cpp` | 클라이언트→서버→클라이언트 왕복, compact 직렬화 크기, 디스크 키, 입력 범위 |
| `src/pke/examples/ship-paper-bench.cpp` | LL13·LL14·HT14·HT15 측정, `trials = 0`이면 보안 예산과 키 크기만 출력, `SHIP_FACTOR_DIR`로 디스크 키 |
| `src/pke/examples/ship-baseline-bench.cpp` | OpenFHE 기존 `EvalBootstrap` 측정 |
| `research/ship/security/ship-estimator.sage` | lattice estimator 실행 스크립트 |

기존 연구 프로토타입(`src/pke/examples/ship/*.h`)은 비교·기록용으로 남겨 두었다.

## 논문과 구현의 대응

| 논문 | 구현 |
| --- | --- |
| Algorithm 1 (p.14): `pt0`, 4개 위상 평문, product tree, `ct + Conj(ct)` | `HalfBootstrap` |
| Lemma 1, Eq. (4)–(5) (p.16–17): 마스크 | `PreRotationMasks`, column 키에서 `Rot_{o+r0}`로 회전 |
| §4.1 (p.18): 마스크 곱셈을 `PQ`에서 하고 `Rescale_P` | `EncryptScaledByP`(스케일 `P` 마스크), `KeySwitchDown` |
| §4.4, Eq. (8), Algorithm 4 (p.21–22): column 방식과 마스크 결합, 이어서 mux 방식 | 인자마다 `θ×4`개 `PQ` 곱셈 후 한 번의 ModDown, 이어서 `MuxRotate` |
| Definition 1, Algorithm 3·5 (p.19–23): HMuxRot, gadget decomposition | `MakeMuxKey`(임시 목적 키 `σ⁻¹(s)` 후 `σ` 적용), HYBRID 분해 |
| §5.1 (p.23–24): B-to-1 mux와 hoisting | mux 단계마다 분해 1회, B개 분기, ModDown 1회 |
| §5.1 (p.24): 균등 간격 sparse 키, `w = max(175, N/(2h))` | `SHIPKeyGen`의 위치 샘플링, 공개 offset `o_k`와 비밀 `r_k < 2w` |
| §3.2, §5.2 (p.15, p.24): S2C는 대각선 방식과 BSGS, 낮은 모듈러스에서 1레벨 | `MakePackingPlan`, `ApplyPacking` |
| §2.2, §5.2: sparse secret encapsulation, 키 전환 모듈러스 `P·q0` | `q0·p'`의 special prime 키 전환 (`Encapsulate`) |
| Table 2 (p.25): Base/S2C/Mult/Boot/Aux 소수 크기와 dnum | `SHIPContextSpec::LL13()`, `LL14()`, `HT14()`, `HT15()`와 `GenSHIPCryptoContext` |

### OpenFHE에 맞춘 구현 선택

- **출력 스케일.** OpenFHE `FIXEDMANUAL`은 모든 rescale 소수를 `2^p`로 가정한다. 곱셈 트리를 `2^b`(b > p)에서 계산하고 마지막 boot 소수를 `2b−p`비트로 두어 출력 스케일을 다시 `2^p`로 맞췄다. 논문 LL13은 boot 25비트 5개이고, 이 구현은 24비트 4개와 28비트 1개이다.
- **소수 편차 보정.** 작은 NTT 소수는 `2^p`에서 수 % 벗어난다. S2C 대각선에 `q/Δ`를 곱하고, 곱셈 트리 전체의 스케일 오차 κ를 미리 계산해 공개 평문 `pt0`에 `1/κ`를 곱한다.
- **`pt0`는 자명한 암호문**(`(pt0, 0)`)으로 넣는다. 평가에 공개키가 필요 없다.
- **보조 모듈러스 크기.** OpenFHE는 `P` 소수를 60비트로 고정한다. LL13의 218비트 예산을 맞추기 위해 SHIP 컨텍스트 생성 시 `P` 소수 크기를 지정해 CRT 표를 다시 계산한다.

## 보안 검증

lattice estimator([malb/lattice-estimator](https://github.com/malb/lattice-estimator) 커밋 `d8c00b48`, 2026-08-19, SageMath 10.5)를 직접 실행했다.
스크립트는 [security/ship-estimator.sage](security/ship-estimator.sage), 원시 결과는 [security/estimator-results.txt](security/estimator-results.txt)이다.
오차 분포는 `σ = 3.19` 이산 가우시안, 수치는 공격 비용의 `log2`이다.

### 1. dense 키 (평가키 전체, uniform ternary, 모듈러스 `QP`)

`GenSHIPCryptoContext`는 OpenFHE의 HE 표준 표와 `log2(QP)`를 비교해 넘으면 예외를 던진다. estimator로는 다음과 같았다.

| 설정 | `log2(QP)` | uSVP | dual hybrid | 참고: HE 표준 경계값(같은 N) |
| --- | --- | --- | --- | --- |
| LL13 | 218 | 126.0 | 126.5 | 218 → 126.0 |
| LL14 | 437 | 126.2 | 126.9 | 438 → 125.9 |
| HT14 | 434 | 127.3 | 128.1 | 438 → 125.9 |
| HT15 | 827 | 134.7 | 135.7 | 881 → 126.1 |

현재 estimator는 HE 표준 표의 경계값 자체를 126비트 정도로 평가한다. SHIP의 dense 키는 HE 표준 경계값과 같거나 더 안전하다.
"128비트"는 이 문서에서 HE 표준 표 기준이라는 뜻이다.

### 2. sparse 키 (h = 31, encapsulation 키 하나, 모듈러스 `q0·p'`)

sparse 키는 dense 키를 sparse 키로 바꾸는 키 전환 키 하나에만 쓰인다. 논문은 `(log N, log(P·q0)) = (13, 55), (14, 100), (15, 105)`이 128비트라고 적었다(§5.2, p.24).
같은 값을 현재 estimator로 계산하면 혼합 공격(`bdd_mitm_hybrid`, `LWE.estimate`가 sparse 키에 쓰는 공격)에서 두 값이 128비트에 못 미쳤다.

| N | `log2(q0·p')` | 혼합 공격 (`bdd_mitm_hybrid`) | 판단 |
| --- | --- | --- | --- |
| 2^13 | 55 (논문) | **121.4** | 부족 |
| 2^13 | 45 | 127.6 | 부족 |
| 2^13 | **42 (이 구현)** | **129.8** | 128 이상 |
| 2^14 | 100 (논문) | **126.2** | 부족 |
| 2^14 | 92 | 128.9 | 128 이상 |
| 2^14 | **88 (이 구현, LL14·HT14)** | **130.2** | 128 이상 |
| 2^15 | **105 (논문, 이 구현)** | **148.7** | 128 이상 |

- uSVP와 dual hybrid는 이 크기에서 적용되지 않거나(무한대) 300비트 이상이었다. 가장 강한 공격은 sparse 키에 특화된 혼합 공격이다.
- 표본 수를 `m = n`(RLWE 표본 하나)으로 제한해도 2^13·55비트 결과는 121.4로 같았다.
- 논문 값과의 차이는 estimator 버전 차이로 확인했다. 논문이 인용한 `14a3625`에서는 primal 혼합 공격이 2^143.6을 주지만, 현재 버전에서는 2^121.4이다(아래 2-1).

그래서 기본값을 estimator 기준으로 바꿨다. `SHIPParams::Recommended`와 `encapsulationModBits = 0`은 `log2(q0·p')`를 2^13에서 42, 2^14에서 88, 2^15에서 105 이하로 고른다.
`p'`의 최소 크기는 17비트로 낮췄다. LL13은 `q0`가 25비트라 `p'`가 17비트가 된다.
LL14는 정밀도가 바뀌지 않았다(16.0–16.5비트, 15회). LL13은 15회 중 최저 1.9비트로, 이전 기본값(`p'` 30비트)의 최저 2.2비트보다 조금 낮은 시도가 있었다. 작은 `p'`의 키 전환 잡음 때문일 수 있다. HT14는 원래 88이었고 HT15는 105 그대로이다.

### 2-1. dual 중간자 혼합 공격과 estimator 버전 (2026-10-11 추가 검토)

논문이 인용한 estimator 버전(`14a3625`, 2024-09-30)과 현재 버전(`d8c00b48`)으로 기본 공격 목록에 없는 dual 중간자 혼합 공격(`dual_mitm_hybrid`)까지 계산했다. 원시 결과: [security/estimator-versions.txt](security/estimator-versions.txt).

| N | `log2(q0·p')` | 14a3625 primal | 14a3625 dual MITM | 현재 primal | 현재 dual MITM | 정밀도 |
| --- | --- | --- | --- | --- | --- | --- |
| 2^13 | 55 (논문) | 143.6 | **116.8** | 121.4 | **117.7** | – |
| 2^13 | 42 (기본값) | 오류 | 126.8 | 129.8 | 127.6 | LL13 1.9–2.5비트 |
| 2^14 | 100 (논문) | 156.1 | **119.1** | 126.2 | **119.9** | LL14 15.9–16.6비트 |
| 2^14 | 88 (기본값) | 오류 | 123.5 | 130.2 | 125.2 | LL14 16.0–16.5비트 |
| 2^14 | 80 | – | 129.7 | 133.1 | 127.3 | LL14 12.0–12.2비트 |
| 2^14 | 76 | – | 128.1 | 134.6 | 129.0 | LL14 8.0–8.3비트 |
| 2^15 | 105 (논문, 기본값) | 오류 | 140.8 | 148.7 | 141.3 | – |

- 논문이 인용한 버전의 기본 공격(primal 혼합)으로는 `(13,55)`가 2^143.6이라 논문 주장과 맞는다. 그러나 dual 중간자 혼합 공격으로는 2^116.8이다.
- dual 공격까지 2^128 이상을 요구하면 `N = 2^14`는 76비트가 필요하고 LL14 정밀도가 8비트로 떨어진다. LL13은 `q0` 25비트와 `p'` 최소 17비트 때문에 도달할 수 없다.
- 기본값은 42·88·105로 유지했다. 현재 estimator의 primal 혼합으로는 2^128 이상이고, dual 중간자 공격으로는 2.5–4.5비트 부족하다(이 공격의 `rop`은 2^77–2^89 메모리를 반영하지 않는다). 더 엄격한 기준이 필요하면 `encapsulationModBits`로 `p'`를 직접 정한다.

### 2-1-1. 자릿수 분해 encapsulation 키 (2026-10-11, 새 기본값)

`p'`를 줄이면 sparse 키 보안은 올라가지만 encapsulation 잡음이 `q0/p'`만큼 커져 정밀도가 떨어진다(LL14: 88비트 16.2, 80비트 12.1, 76비트 8.1비트).
그래서 encapsulation 키에 gadget 분해를 넣었다. 입력의 `q0` 성분을 `d`비트 부호 있는 자릿수 `k = ceil(log2 q0 / d)`개로 나누고, 자릿수마다 `Enc(p' · 2^(jd) · s)` 키를 둔다. 잡음은 `k · 2^(d-1) · σ√N / p'`로 `q0`와 무관해진다.

| LL14 설정 | `log2(q0·p')` | 정밀도 |
| --- | --- | --- |
| 분해 없음, `p'` 28비트 | 76 | 8.1–8.2비트 |
| 16비트 자릿수 | 76 | 16.0–16.3비트 |
| 12비트 자릿수 | 76 | 16.1–16.4비트 |
| **8비트 자릿수, `p'` 17비트** | **65** | **16.2–16.5비트** |

- 새 기본값(`SHIPParams::Recommended`): `encapsulationDigitBits = 8`, `encapsulationModBits = 17`. `log2(q0·p')`는 LL13 42, LL14 65, HT14 45, HT15 68이다.
- 새 기본값 정밀도: LL13 2.0–2.5, LL14 15.9–16.4, HT14 4.0–4.2, HT15 18.3–18.4비트로 이전과 같다.
- encapsulation 시간은 1–1.5ms에서 3.4–5ms로 늘었다(부트스트래핑 1초 중 0.5% 미만). 키는 자릿수 수(최대 7)만큼 다항식 쌍이 늘어 수백 KB이다.
- 단위 테스트 `DigitDecomposedEncapsulationAllowsSmallSpecialPrime`를 추가했다(SHIP 테스트 10개). 직렬화 형식 버전은 2로 올렸다.

| 새 기본값 (현재 estimator) | 일반 h=31 primal | 일반 dual MITM | 같은 엔트로피 primal | 같은 엔트로피 dual MITM |
| --- | --- | --- | --- | --- |
| LL13 `(13,42)` | 129.8 | 127.6 | 122.0 | 119.2 |
| LL14 `(14,65)` | 139.4 | 135.6 | 127.9 | 122.5 |
| HT14 `(14,45)` | 150.3 | 147.4 | 137.2 | 133.6 |
| HT15 `(15,68)` | 161.3 | 155.7 | 147.0 | 140.7 |

일반 모델로는 LL14·HT14·HT15가 모든 공격에서 2^135 이상이다. 균등 간격 구조를 같은 엔트로피로 근사하면 LL13·LL14는 2^128에 못 미친다. 남은 제약은 encapsulation 모듈러스가 아니라 균등 간격 키 구조다.

### 2-2. 균등 간격 구조를 같은 엔트로피의 일반 sparse 키로 근사 (휴리스틱)

SHIP 키는 구간마다 0이 아닌 계수가 하나라 엔트로피가 일반 sparse 키보다 작다(LL13 293.0비트 대 321.3비트, LL14 311.4비트 대 352.3비트).
같은 엔트로피의 일반 sparse 키(`N = 2^13`은 h = 28, `2^14`는 h = 27)로 현재 estimator를 돌렸다. 원시 결과: [security/estimator-windowed.txt](security/estimator-windowed.txt).

| 지점 | primal 혼합 | dual 중간자 혼합 |
| --- | --- | --- |
| LL13 기본값 `(13,42)` | 122.0 | 119.2 |
| 논문 `(13,55)` | 114.5 | 110.8 |
| LL14 기본값 `(14,88)` | 120.0 | 113.5 |

일반 h = 31 모델보다 7–12비트 낮다. 엔트로피를 맞춘 근사일 뿐이며 실제 공격이 구조를 더 잘 쓸 수도, 못 쓸 수도 있다. `(14,100)`, `(14,76)`, `(15,105)`는 Binder 세션이 끊겨 아직 계산하지 못했다.

### 3. 균등 간격 위치

SHIP의 sparse 키는 0이 아닌 계수가 구간마다 하나씩 놓인다. estimator는 이 구조를 모형화하지 않는다.
논문의 [May21] 공격 비용 추정 `(4w)^(0.4·h)`를 계산해 확인한다(§5.1, p.24). 논문의 기준은 `2^115` 이상이며, LL13 `2^117.2`, LL14·HT14 `2^124.6`, HT15 `2^137.0`이다.
LL13의 값은 128보다 작다. 이 기준은 논문의 선택이며 이 구현이 따로 검증하지 않았다.

평가키가 키에 의존하는 메시지(`β`, `β·s`, 마스크)를 암호화한다는 점은 일반 CKKS key switching 키와 같은 종류의 가정이며, 논문도 같은 구성을 쓴다.
이 문서는 새 보안 증명을 주장하지 않는다.

## 결과

환경: Apple M5 Pro(15코어: 성능·효율 혼합), 48GB RAM, macOS arm64, AppleClang, OpenFHE 1.4.2 기반, `NATIVE_SIZE=64`,
`WITH_REDUCED_NOISE=ON`. 입력은 `[-1, 1]` 균등 난수 실수 벡터이며 정밀도는 논문과 같이 `log2(1/ε)`, ε는 슬롯별 실수부 최대 오차다.
지연 시간은 키 생성을 제외한 `EvalSHIPBootstrap` 한 번이며 중앙값과 [최소–최대]를 적었다.
원시 기록은 [results/library](results/library/)에 있다.

### SHIP (128비트)

| 항목 | LL13 대응 | LL14 대응 |
| --- | --- | --- |
| N / 슬롯 | 2^13 / 4,096 | 2^14 / 8,192 |
| `log2(QP)` / HE 표준 상한 | 218 / 218 | 437 / 438 |
| Q 소수(비트) | 25, 20, 21, 28, 24, 25, 24, 25 | 48, 37, 37, 48, 42, 43, 42, 42 |
| P 소수(비트), dnum | 28, 8 | 50·50, 4 |
| `log2(q0·p')` (논문 / 이 구현의 현재 기본값) | 55 / 42 | 100 / 88 |
| h, w, θ, B | 31, 175, 6, 4 | 31, 264, 9, 4 |
| 부트스트래핑 깊이 / 남는 곱셈 레벨 | 6 / 1 | 6 / 1 |
| 정밀도(시도별 범위, 15회) | 1.9–2.5비트 | 16.0–16.5비트 |
| 제곱 후 재부트스트래핑 정밀도 | 1.5비트 | 15.4비트 |
| 1스레드 | 2.98초 [2.97–3.14] | 5.28초 [5.24–5.74] |
| 4스레드 | 0.90초 [0.87–0.92] | 1.45초 [1.41–1.48] |
| 15스레드 | **0.58초** [0.58–0.59] (5.1배) | **0.83초** [0.81–0.88] (6.4배) |
| 키 생성 | 17초 | 24–26초 |
| SHIP 키 / 최대 RSS | 8.0 GiB / 8.4–8.6 GiB | 12.2 GiB / 13.5–13.8 GiB |
| 논문 (HEaaN, AMD EPYC, 1 / 4 / 32코어) | 4.45비트, 3.02초 / 803ms / 215ms | 16.9비트, 4.91초 / 1.29초 / 333ms |

위 표는 최종 코드(시드 저장 키, 레벨별로 줄인 S2C·켤레 키, `log2(q0·p')` 42·88)로 다른 무거운 작업이 없을 때 잰 값이다(2026-10-11, 설정마다 5회).
OpenMP 빌드에서 프로세스마다 `OMP_NUM_THREADS`를 지정했다. 원시 기록은 [results/library/quiet](results/library/quiet/)에 있다.
초판 코드(2026-10-10)보다 15스레드 지연 시간이 줄었다(LL13 0.64 → 0.58초, LL14 0.95 → 0.83초).

### HT14·HT15 (128비트, 디스크 기반 키, 4스레드)

HT14·HT15의 키는 메모리 48GB에 들어가지 않아 `SHIP_FACTOR_DIR`로 인자별 키를 디스크(내장 SSD)에 두고 실행했다.
측정은 2회이며 다른 프로그램(게임)이 함께 돌았다. 지연 시간은 디스크 읽기를 포함하므로 논문과 직접 비교하지 않는다.

| 항목 | HT14 대응 | HT15 대응 |
| --- | --- | --- |
| N / 슬롯 | 2^14 / 8,192 | 2^15 / 16,384 |
| `log2(QP)` / HE 표준 상한 | 434 / 438 | 827 / 881 |
| Q 소수(비트): Base, S2C, Mult ×9, Boot ×5 | 28, 23, 23–24, 24–25 | 51, 39, 39–40, 46–54 |
| P 소수(비트), dnum | 55, 8 | 50·50·50, 6 (논문 53×4, 4) |
| `log2(q0·p')` | 88 | 105 |
| h, w, θ, B | 31, 264, 9, 4 | 31, 528, 17, 4 |
| 남는 곱셈 레벨 | **9** (논문 9) | **9** (논문 9) |
| 정밀도 | **4.0비트** (논문 4.55) | **18.3–18.4비트** (논문 18.1) |
| 제곱 후 재부트스트래핑 정밀도 | 2.2비트 | 16.0비트 |
| 지연 시간 (4스레드, 디스크 읽기 포함) | 8.1–19.2초 | 49.1–50.5초 |
| 논문 4코어 / 1코어 | 3.62초 / 13.2초 | 6.09초 / 29.6초 |
| 키 생성 | 196초 | 354초 |
| 예상 키 크기 / 디스크 사용 | 31.5 GiB / 29 GB | 69.6 GiB / 약 60 GB |
| 메모리에 남는 키 / 최대 RSS | 2.2 GiB / 14.4 GiB | 8.6 GiB / 27.3 GiB |

HT15는 논문보다 정밀도가 약간 높았다. HT14는 LL13과 같은 이유(작은 스케일에서 OpenFHE dense 키 잡음)로 논문보다 0.5비트 낮다.
지연 시간은 디스크에서 키 수십 GB를 읽는 시간이 대부분이다. 메모리가 충분한 장비에서 다시 재야 한다.

### 배포 기능 검증

| 확인 | 결과 |
| --- | --- |
| 클라이언트에서 키 생성·직렬화 → 서버에서 읽어 부트스트래핑 → 클라이언트 복호화 (`ship-deployment-checks`) | 통과, 15.7비트 |
| compact 직렬화 크기 (toy 파라미터) | 49.4 MB, 전체 저장 98.7 MB의 0.500배 |
| 디스크 기반 키 (toy) | 메모리 12 MB (메모리 모드 107 MB), 결과 오차 동일 |
| `messageBound` (입력 `[-64, 64]`) | K = 1: 최대 오차 1.62 (사인 근사 실패) / K = 64: 4×10⁻⁴ |
| OpenFHE 전체 단위 테스트 (`WITH_REDUCED_NOISE=ON`, 변경 후 재실행) | core 156개, binfhe 78개, pke 1,885개(SHIP 9개 포함) 모두 통과 |

### 복호화 결과 공개 시 noise flooding

CKKS 복호화 결과를 다른 사람에게 공개하면 IND-CPA^D 공격([LM21])을 막기 위해 noise flooding이 필요하다.
OpenFHE의 `NOISE_FLOODING_DECRYPT`는 `EXEC_EVALUATION` 모드에서 스케일을 약 2^55 이상으로 강제하므로 SHIP의 작은 스케일(2^20–2^39) 컨텍스트에는 그대로 쓸 수 없다.
그래서 같은 계산을 벤치마크에서 재현했다. 출력 허수부로 잡음을 추정하고, OpenFHE 기본값(통계적 보안 30비트, 질의 1회)으로 표준편차 `2^(15 + log2√12 + 추정값)`의 가우시안 잡음을 더해 복호화했다.

| 세트 | 스케일 | flooding 표준편차 | 부트스트래핑 정밀도 → flooding 후 |
| --- | --- | --- | --- |
| LL13 | 2^20 | 2^21.4 | 2.4비트 → **복호화 의미 없음** (오차 574) |
| LL14 | 2^37 | 2^21.7 | 16.0비트 → 6.9비트 |
| HT15 | 2^39 | 2^22.2 | 18.3비트 → 7.4비트 |

flooding 잡음은 출력 잡음의 약 2^17배라 스케일이 작은 SHIP에서는 정밀도 손실이 크다. 복호화 결과를 공개하는 응용에서는 flooding 전에 스케일을 키우거나 정밀도 요구를 낮춰야 한다.
HT14 기록의 flooding 줄은 잡음을 작은 소수로 줄이기 전의 코드로 잰 값이라 무효이다.

### OpenFHE 기존 부트스트래핑 (`EvalBootstrap`, 128비트, 남는 곱셈 레벨 1)

| 설정 | N | `log2(QP)` | 정밀도 | 1스레드 | 15스레드 | 최대 RSS |
| --- | --- | --- | --- | --- | --- | --- |
| 스케일 34 / 첫 소수 40, budget 3,3, dnum 11 | 2^15 | 875 / 881 | **복호화 실패** (근사 오차 과대) | – | – | 13.8 GiB |
| 스케일 54 / 56, budget 3,3 | 2^16 | 1671 / 1747 | 5.3–5.6비트 | 11.5초 [10.7–13.8] | 9.4초 [4.5–10.0] | 10.0–10.4 GiB |
| 스케일 59 / 60, budget 2,2 (OpenFHE 권장 스케일) | 2^16 | 1602 / 1747 | 10.3–10.6비트 | 30.3초 [26.8–30.4] | 8.8초 [5.4–12.1] | 24.3–29.9 GiB |

uniform ternary 키의 기존 방식은 부트스트래핑 깊이가 19–21레벨이라 같은 보안에서 `N = 2^16`이 필요했다. 스케일 32/38·34/40으로 `N = 2^15`에 넣으면 복호화가 실패했다.
이 세 설정에서는 SHIP LL14(`N = 2^14`, 16비트, 1스레드 5.3–5.6초)가 정밀도도 높고 빨랐다. 그러나 아래처럼 기존 방식에 sparse-secret encapsulation을 쓰면 정밀도는 기존 방식이 더 높다.

#### sparse-secret encapsulation으로 조정한 기존 방식과 같은 조건 비교

위 세 설정은 uniform ternary 키 그대로였다. OpenFHE에는 SHIP과 같은 발상의 `SPARSE_ENCAPSULATED`(h = 32 sparse 키로 바꿔 ModRaise, 키 전환 모듈러스 `q0·p0`)가 있어 이것으로 다시 쟀다.
`ship-baseline-bench`의 7번째 인자 `sparse_encapsulated`로 켠다. 기존 방식의 sparse 키도 estimator로 확인했다(`N = 2^15`, `log2(q0·p0) = 112`: 2^149.7, `N = 2^16`, 120: 2^172.1).

다른 무거운 작업이 없을 때 스레드 수마다 SHIP LL13, LL14, 기존 방식 두 설정을 연달아 돌렸다(SHIP 5회, 기존 방식 3회, 지연 시간은 중앙값).
처리량은 논문 식 (9)의 `슬롯 × 남는 레벨 × 정밀도 / 지연 시간`이며 정밀도는 최저값을 썼다.

| 설정 | N / 슬롯 | 정밀도 | 1스레드 | 4스레드 | 15스레드 | 처리량 1 / 4 / 15스레드 | 키 생성(1스레드) / 최대 RSS |
| --- | --- | --- | --- | --- | --- | --- | --- |
| SHIP LL13 | 2^13 / 4,096 | 1.9–2.5비트 | 2.98초 | 0.90초 | 0.58초 | 0.3만 / 1.0만 / 1.3만 | 17초 / 8.6 GiB |
| **SHIP LL14** | 2^14 / 8,192 | 16.0–16.5비트 | **5.28초** | **1.45초** | **0.83초** | 2.5만 / 9.0만 / 15.8만 | 26초 / 13.8 GiB |
| 기존 방식, 스케일 59/60, budget 2,2 | 2^16 / 32,768 | **18.3–18.7비트** | 10.48초 | 5.00초 | 3.48초 | **5.8만 / 12.1만 / 17.3만** | 52초 / 25.3 GiB |
| 기존 방식, 스케일 45/52, budget 3,3, dnum 18 | 2^15 / 16,384 | 6.7–7.2비트 | 4.89초 | 2.44초 | 1.85초 | 2.3만 / 4.7만 / 6.1만 | 50초 / 17.3 GiB |

- **지연 시간.** SHIP LL14가 기존 방식 `N = 2^16`보다 1스레드 2.0배, 4스레드 3.4배, 15스레드 4.2배 짧다. 스레드가 늘수록 차이가 커지는 것은 31개 인자를 독립적으로 계산하는 SHIP 구조 때문이다.
- **정밀도.** 기존 방식 `N = 2^16`이 2비트 정도 높다.
- **처리량.** 기존 방식 `N = 2^16`이 1스레드 2.3배, 4스레드 1.3배, 15스레드 1.1배 높다. 슬롯이 4배 많기 때문이다.
- **메모리.** SHIP LL14가 약 절반이다.
- **`N = 2^15`의 기존 방식**은 OpenFHE의 보정 조건(`log2(q0/Δ) ≤ 7`)과 881비트 예산 때문에 스케일을 45비트 이하로 둬야 했고, 정밀도가 7비트에 그쳤다. 스케일 40/47에서는 1비트였다. 같은 `N = 2^15`의 SHIP HT15는 18.3비트이다.
- 남는 곱셈 레벨을 1로 고정한 비교이다. 기존 방식은 레벨을 더 남기도록 설정할 수 있으며, 그때 처리량이 늘어난다(논문 Table 3의 HEaaN Param16은 9레벨).
- 게임이 함께 돌 때 잰 이전 기록(`ab-*.txt`)은 편차가 커서 이 표로 대체한다.

비교의 한계:
- 기존 방식은 sparse-secret encapsulation까지 시험했다. composite scaling, 반복 부트스트래핑(meta-BTS), 남는 레벨을 늘린 설정은 시험하지 않았다.
- 슬롯 수가 다르다(기존 `N = 2^16` 32,768, LL14 8,192). uniform 키 설정을 1스레드 처리량(슬롯 × 레벨 × 정밀도 / 지연)으로 비교하면 LL14 약 2.3만–2.5만, 기존 54/56 약 1.5만, 59/60 약 1.1만이다.
- uniform 키 기존 방식 표는 Spotlight 색인 등 다른 프로세스가 함께 돌 때 쟀다. 15스레드와 24–30GB를 쓰는 설정에서는 시도별 편차가 컸다(4.5–12초).
- 이 Mac의 15코어는 성능·효율 코어가 섞여 있다. 논문의 32코어 확장성(14–18배)과 직접 비교하지 않는다.

## OpenFHE 핵심 코드 변경: 정확한 반올림 ModDown

`DCRTPolyImpl::ApproxModDown`은 `[x]_P`를 `[0, P)` 범위의 근사 기저 변환으로 구한다.
그 결과 나눗셈이 사실상 내림이 되어 반올림 오차의 평균이 1/2이 된다. 이 오차에 dense 키가 곱해지면 계수마다 수십 단위의 **편향된 저주파 오차**가 생긴다.
`P`가 소수 여러 개이면 근사 변환의 `αP` 항 때문에 계수마다 ±1의 오차가 더해진다.

2^40 이상의 스케일에서는 드러나지 않지만, SHIP의 LL13(2^20~2^24)에서는 키 스위칭 수백 번에 걸쳐 누적되어 결과를 망가뜨렸다.
진단에서는 인자마다 0번 슬롯 오차가 0.02~0.25였고, 최종 오차가 0.8~5였다.

`WITH_REDUCED_NOISE`가 켜져 있고 `t = 0`(CKKS, BFV)일 때, `r = Σ y_j·(P/p_j) − round(Σ y_j/p_j)·P`로 중앙값 기준 나머지를 정확히 계산해 빼도록 바꿨다.
결과는 `round(x/P)`이다. BGV(`t > 0`)와 이 옵션이 꺼진 빌드는 기존과 같다.

| 확인 | 결과 |
| --- | --- |
| OpenFHE `core_tests` / `binfhe_tests` / `pke_tests` (`WITH_REDUCED_NOISE=ON`) | 모두 통과 (pke 1,882개) |
| LL13 (균일 22비트 소수, `P` 28비트) | 옵션 없음: 최대 오차 1.5~5 / `WITH_REDUCED_NOISE`: 1.3~1.6비트 |
| LL13 (같은 소수, `P` 60비트, 진단용) | 옵션 없음: 최대 오차 0.55~2.5 (0번 슬롯 집중) / `WITH_REDUCED_NOISE`: 3.3비트 |
| LL14 (균일 40비트 소수) | `WITH_REDUCED_NOISE`: 13.8~13.9비트 / 정확한 반올림 추가: 14.4비트 |
| 역할별 소수(Table 2) 적용 후 | LL13 2.3~2.5비트, LL14 16.1~16.5비트 |

### LL13 정밀도 원인 분석 (2026-10-11)

`N = 2^13`에서 한 번에 하나씩만 바꿔(128비트 예산 밖, 비밀키 암호화, 4회) 입력 복호값 대비 정밀도를 쟀다. 원시 기록: [results/library/ll13-analysis](results/library/ll13-analysis/).

| 변경 | 정밀도 |
| --- | --- |
| 기본(LL13과 같은 크기) | 2.2–2.4비트 |
| `p'` 17 → 30비트 | 2.2–2.7비트 |
| P 28 → 40비트 | 3.2–3.6비트 |
| P 28 → 60비트 | 3.4–3.6비트 |
| P 28비트 2개(dnum 4) | 3.4–3.5비트 |
| 곱셈 트리 소수 26–30비트 (P 28비트 2개로 늘어남) | 3.3–3.7비트 |
| P 60비트 + 트리 소수 30비트 | **5.5–5.7비트** (논문 4.45보다 높음) |

LL13 정밀도는 218비트 예산 때문에 작게 잡을 수밖에 없는 보조 소수 P(28비트)와 트리 소수(24비트)에서 생기는 키 전환 잡음에 주로 묶여 있다. HEaaN은 같은 크기에서 4.45비트를 얻으므로 이 크기에서 OpenFHE 쪽 키 전환 잡음이 더 크다. 이유는 아직 확인하지 못했다.
`q0`를 25비트에서 26–30비트로 늘리면(스케일 2^20 고정) 1비트마다 정밀도가 약 1비트씩 떨어졌다(26비트 1.2, 27비트 0.3, 30비트 −1.6비트). 결함이 아니라 출력 `γ/(2π)·sin(2πμ/γ)`가 곱셈 트리의 상대 오차를 간격 `γ = q0/Δ`만큼 키우기 때문이다.

### 기준선 dnum·스케일 조정 (4스레드, `N = 2^16`, 8,192 슬롯, level budget (3,3))

dnum 2·3·4·6은 2.30·2.20·2.34·2.71초, 스케일 50·55비트는 2.20초(12.0비트)·2.33초(18.1비트)였다. 기본 설정(dnum 3, 59/60비트)이 최선과 3% 이내이므로 비교 결론은 바뀌지 않는다.

## 논문과의 차이와 한계

- **LL13 정밀도.** 논문 4.45비트에 못 미친다. OpenFHE의 uniform ternary dense 키에서는 `N = 2^13`, 스케일 2^20의 공개키 암호화·복호화만으로 슬롯 최대 오차가 0.10–0.11이었다(4,096개 단위 복소수, 평균 0.024). 논문 LL13의 전체 오차(0.046)보다 크다. HEaaN과 OpenFHE의 키 분포와 잡음 특성 차이로 보이며, 같은 키 분포로 직접 비교하지는 않았다.
- **LL14 정밀도.** 16.0–16.5비트로 논문 16.9비트보다 약간 낮다. boot 소수가 42비트(논문 43비트)이고 OpenFHE의 dense 키 잡음이 남는다.
- **HT14 정밀도.** 4.0비트로 논문 4.55비트보다 낮다. LL13과 같은 원인으로 보인다.
- **메모리.** 키는 LL13 약 8GB, LL14 약 12GB, HT14 약 31GB, HT15 약 70GB이다. HT14·HT15는 디스크 기반 키로만 48GB 장비에서 돌릴 수 있고, 그때 지연 시간은 디스크 읽기에 좌우된다. 논문 장비는 512GB RAM이다.
- **병렬성.** 31개 인자 계산과 곱셈 트리, S2C 그룹을 OpenMP로 병렬화했다. 이 Mac은 15코어(성능·효율 혼합)이며 논문은 32~48코어이다.
- **보안.** sparse 키의 기본값을 estimator 기준으로 바꿨지만 dual 중간자 혼합 공격으로는 2^128에 2.5–4.5비트 못 미친다(위 2-1). 또한 균등 간격 구조에 대한 [May21] 기준(LL13 2^117.2)은 논문의 선택을 따랐다. 독립적인 보안 검토는 받지 않았다.
- **복호화 결과 공개.** OpenFHE의 noise flooding 복호화는 SHIP의 작은 스케일과 맞지 않는다. 같은 계산을 적용하면 LL13은 결과가 무너지고 LL14·HT15는 7비트 안팎으로 줄었다.
- **compact 직렬화.** 균등 난수 성분은 키별로 새로 뽑은 시드에서 다시 만든다. 시드는 공개되어도 되는 값(공개 성분 `a`를 정하는 값)이지만, 읽는 쪽과 쓰는 쪽의 OpenFHE 판과 PRNG가 같아야 한다.
- **`FIXEDMANUAL`, `HYBRID`, full packing만 지원한다.**
- 사용자 연산의 곱셈 레벨 소수는 `2^p`에 가장 가까운 NTT 소수이다. LL13처럼 작은 소수에서는 사용자의 rescale에 `|q/2^p − 1|` 크기의 스케일 오차가 생긴다(OpenFHE `FIXEDMANUAL`의 일반 특성).

## 재현

```sh
git submodule update --init third-party/cereal third-party/google-test
bash research/ship/run-checks.sh   # WITH_REDUCED_NOISE=ON 빌드, 연구 검사, 라이브러리 검사, SHIP gtest
build-ship/bin/examples/pke/ship-paper-bench LL13 5
build-ship/bin/examples/pke/ship-paper-bench LL14 5
build-ship/bin/examples/pke/ship-paper-bench HT15 0                       # 보안 예산과 키 크기만 출력
SHIP_FACTOR_DIR=/path/to/disk build-ship/bin/examples/pke/ship-paper-bench HT15 2   # 디스크 기반 키 (약 60GB)
SHIP_ENCAPS_BITS=20 build-ship/bin/examples/pke/ship-paper-bench LL13 2   # p' 크기 직접 지정
build-ship/bin/examples/pke/ship-deployment-checks
sage -python research/ship/security/ship-estimator.sage                   # lattice estimator (수십 분)
```

다중 코어 측정은 `-DWITH_OPENMP=ON`으로 빌드하고 프로세스마다 `OMP_NUM_THREADS`를 지정한다.
OpenFHE의 내부 병렬 루프는 시작 시점의 `OMP_NUM_THREADS`로 스레드 수를 고정하므로 실행 중 `omp_set_num_threads`로는 바뀌지 않는다.
macOS AppleClang에서는 Homebrew `libomp`를 `OpenMP_*` CMake 변수로 지정했다.
