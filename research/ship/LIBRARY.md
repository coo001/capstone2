# SHIP 부트스트래핑의 OpenFHE 라이브러리 구현

2026-10-11. 원논문은 Cheon, Hanrot, Kim, Stehlé, "SHIP: A Shallow and Highly Parallelizable CKKS Bootstrapping Algorithm", EUROCRYPT 2025([ePrint 2025/784](https://eprint.iacr.org/2025/784.pdf), 33쪽 판)이다.
페이지 번호는 ePrint 판 기준이다.

## 요약

- SHIP을 OpenFHE 라이브러리 기능으로 구현했다. 예제 헤더가 아니라 `src/pke/lib`에 들어가며 `CryptoContext` API로 호출한다.
- 논문의 알고리즘 구성(Algorithm 1, Lemma 1 마스크, column과 mux를 결합한 블라인드 회전, base-4 B-to-1 mux, 균등 간격 sparse 키, `P·q0` sparse-secret encapsulation, 대각선 S2C, Table 2의 역할별 소수)을 구현했다.
- 128비트 HE 표준을 만족하는 LL13·LL14 대응 파라미터로 실행했다. `log2(QP)`는 각각 218, 437로 논문 Table 2와 같다.
- LL14는 논문의 16.9비트에 가까운 **16.0–16.6비트**, 1스레드 5.3–5.6초, 15스레드 0.95초였다. LL13은 논문 4.45비트보다 낮은 **2.2–2.6비트**, 1스레드 3.2초, 15스레드 0.64초였다. 차이의 원인은 아래 "논문과의 차이"에 정리했다.
- 같은 128비트·곱셈 레벨 1 조건에서 OpenFHE 기존 `EvalBootstrap`은 `N = 2^16`이 필요했고, 시험한 설정에서 5.3–10.6비트, 1스레드 11.5–30.3초였다.
- 작은 CKKS 스케일에서 정밀도를 위해 OpenFHE의 hybrid key switching ModDown에 정확한 반올림을 추가했다. `WITH_REDUCED_NOISE` 빌드 옵션에서만 동작하며, OpenFHE 전체 단위 테스트(core, binfhe, pke 1,882개)가 통과했다.

## 사용법

```cpp
#include "openfhe.h"
using namespace lbcrypto;

auto cc = GenSHIPCryptoContext(SHIPContextSpec::LL14());  // 128-bit check included
auto kp = cc->KeyGen();
cc->EvalMultKeyGen(kp.secretKey);
SHIPParams params;           // h = 31, w = 175, theta = 6, B = 4, real numbers
params.window = 264;         // LL14: w = max(175, N / (2h))
params.columnSize = 9;
params.encapsulationModBits = 52;
cc->EvalSHIPBootstrapKeyGen(kp.secretKey, params);

auto ct  = cc->Encrypt(kp.publicKey, cc->MakeCKKSPackedPlaintext(x, 1, numPrimes - 2, nullptr, N / 2));
auto out = cc->EvalSHIPBootstrap(ct);  // ceil(log2(h + 1)) levels from the top of the chain
```

- 입력은 full packing, rescale을 마친 암호문이며 limb가 2개 이상이어야 한다. 더 많으면 S2C 레벨까지 내린다.
- 메시지는 `|m| <= 1` 정도로 작아야 한다. SHIP은 `γ/(2π)·sin(2πm/γ)`를 계산하기 때문이다(논문 §1.2).
- 출력은 곱셈 트리 깊이만큼 위쪽 레벨을 쓴 뒤의 암호문이다. LL13·LL14에서는 곱셈 레벨 1개가 남는다.
- `params.realOnly = false`이면 복소수 메시지를 위해 half-bootstrap을 두 번 수행한다.
- `ClearEvalSHIPBootstrapKeys()`로 키를 지운다. 직렬화는 아직 지원하지 않는다.

## 파일

| 파일 | 내용 |
| --- | --- |
| `src/pke/include/scheme/ckksrns/ckksrns-ship.h` | 공개 API: `SHIPParams`, `SHIPContextSpec`, `GenSHIPCryptoContext`, `SHIPKeyGen`, `SHIPBootstrap` |
| `src/pke/lib/scheme/ckksrns/ckksrns-ship.cpp` | 구현 |
| `src/pke/include/cryptocontext.h` | `EvalSHIPBootstrapKeyGen`, `EvalSHIPBootstrap`, `ClearEvalSHIPBootstrapKeys` |
| `src/core/include/lattice/hal/default/dcrtpoly-impl.h` | `ApproxModDown`의 정확한 반올림 (`WITH_REDUCED_NOISE`일 때) |
| `src/pke/unittest/utckksrns/UnitTestSHIP.cpp` | gtest 6개 |
| `src/pke/examples/ship-library-checks.cpp` | 기능 검사 (toy 파라미터) |
| `src/pke/examples/ship-paper-bench.cpp` | LL13·LL14 측정 |
| `src/pke/examples/ship-baseline-bench.cpp` | OpenFHE 기존 `EvalBootstrap` 측정 |

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
| Table 2 (p.25): Base/S2C/Mult/Boot/Aux 소수 크기와 dnum | `SHIPContextSpec::LL13()`, `LL14()`와 `GenSHIPCryptoContext` |

### OpenFHE에 맞춘 구현 선택

- **출력 스케일.** OpenFHE `FIXEDMANUAL`은 모든 rescale 소수를 `2^p`로 가정한다. 곱셈 트리를 `2^b`(b > p)에서 계산하고 마지막 boot 소수를 `2b−p`비트로 두어 출력 스케일을 다시 `2^p`로 맞췄다. 논문 LL13은 boot 25비트 5개이고, 이 구현은 24비트 4개와 28비트 1개이다.
- **소수 편차 보정.** 작은 NTT 소수는 `2^p`에서 수 % 벗어난다. S2C 대각선에 `q/Δ`를 곱하고, 곱셈 트리 전체의 스케일 오차 κ를 미리 계산해 공개 평문 `pt0`에 `1/κ`를 곱한다.
- **`pt0`는 자명한 암호문**(`(pt0, 0)`)으로 넣는다. 평가에 공개키가 필요 없다.
- **보조 모듈러스 크기.** OpenFHE는 `P` 소수를 60비트로 고정한다. LL13의 218비트 예산을 맞추기 위해 SHIP 컨텍스트 생성 시 `P` 소수 크기를 지정해 CRT 표를 다시 계산한다.

## 보안 검증

검증은 세 부분이다. 직접 lattice estimator를 돌리지는 않았고, OpenFHE의 HE 표준 표와 논문이 estimator로 제시한 상한을 따랐다.

1. **dense 키 (평가키 전체).** OpenFHE의 HE 표준 표(uniform ternary, 128-bit classic)와 실제 `log2(QP)`를 비교하며, 넘으면 `GenSHIPCryptoContext`가 예외를 던진다. 평가키·마스크 암호문은 모두 dense 키 아래 `QP`에서 만들어진다.
2. **sparse 키.** sparse 키는 `q0·p'` 모듈러스의 키 전환 키에만 쓰인다. 논문은 lattice estimator로 `h = 31`에서 `(log N, log(P·q0)) = (13, 55), (14, 100)`이 128비트임을 제시했다(§5.2, p.24). 이 구현은 `log2(q0·p')`를 55와 100 이하로 둔다.
3. **균등 간격 위치.** 논문의 [May21] 공격 비용 추정 `(4w)^(0.4·h) > 2^115`를 계산해 확인한다(§5.1, p.24). LL13은 `2^117.2`, LL14는 `2^124.6`이다.

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
| `log2(q0·p')` / 논문 상한 | 55 / 55 | 100 / 100 |
| h, w, θ, B | 31, 175, 6, 4 | 31, 264, 9, 4 |
| 부트스트래핑 깊이 / 남는 곱셈 레벨 | 6 / 1 | 6 / 1 |
| 정밀도(시도별 범위) | 2.2–2.6비트 | 15.9–16.6비트 |
| 제곱 후 재부트스트래핑 정밀도 | 1.5–1.6비트 | 15.1–15.4비트 |
| 1스레드 (OpenMP 빌드) | 3.22초 [3.15–3.50] | 5.59초 [5.58–6.08] |
| 1코어 (OpenMP 없는 빌드) | 3.17초 [3.12–3.43] | 5.30초 [5.12–10.08] |
| 2 / 4 / 8스레드 | 1.77 / 1.03 / 0.69초 | – / 1.82 / 1.22초 |
| 15스레드 | **0.64초** (5.0배) | **0.95초** (5.9배) |
| 키 생성 | 17–27초 | 32–44초 |
| SHIP 키 payload / 최대 RSS | 7.9 GiB / 9.4–9.6 GiB | 12.0 GiB / 15.1–15.5 GiB |
| 논문 (HEaaN, AMD EPYC, 1코어 / 32코어) | 4.45비트, 3.02초 / 215ms | 16.9비트, 4.91초 / 333ms |

### OpenFHE 기존 부트스트래핑 (`EvalBootstrap`, 128비트, 남는 곱셈 레벨 1)

| 설정 | N | `log2(QP)` | 정밀도 | 1스레드 | 15스레드 | 최대 RSS |
| --- | --- | --- | --- | --- | --- | --- |
| 스케일 34 / 첫 소수 40, budget 3,3, dnum 11 | 2^15 | 875 / 881 | **복호화 실패** (근사 오차 과대) | – | – | 13.8 GiB |
| 스케일 54 / 56, budget 3,3 | 2^16 | 1671 / 1747 | 5.3–5.6비트 | 11.5초 [10.7–13.8] | 9.4초 [4.5–10.0] | 10.0–10.4 GiB |
| 스케일 59 / 60, budget 2,2 (OpenFHE 권장 스케일) | 2^16 | 1602 / 1747 | 10.3–10.6비트 | 30.3초 [26.8–30.4] | 8.8초 [5.4–12.1] | 24.3–29.9 GiB |

OpenFHE 기존 방식은 부트스트래핑 깊이가 19–21레벨이라 같은 보안에서 `N = 2^16`이 필요했다. 스케일 32/38·34/40으로 `N = 2^15`에 넣으면 복호화가 실패했다.
SHIP LL14는 `N = 2^14`에서 16비트, 1스레드 5.3–5.6초, 15스레드 0.95초였다. 시험한 기존 방식 설정보다 정밀도가 높고 빠르며 메모리도 비슷하거나 적었다.

비교의 한계:
- 기존 방식의 설정은 위 세 가지만 시험했다. sparse 키, composite scaling, 반복 부트스트래핑(meta-BTS) 같은 OpenFHE의 다른 기법으로 조정하지 않았다.
- 슬롯 수가 다르다(기존 32,768, LL14 8,192). 논문 식 (9)의 처리량(슬롯 × 레벨 × 정밀도 / 지연)으로 1스레드를 비교하면 LL14 약 2.3만–2.5만, 기존 54/56 약 1.5만, 59/60 약 1.1만이다.
- 측정 중 Spotlight 색인 등 다른 프로세스가 함께 돌았다. 기존 방식 15스레드와 24–30GB를 쓰는 설정에서는 시도별 편차가 컸다(4.5–12초).
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

## 논문과의 차이와 한계

- **LL13 정밀도.** 논문 4.45비트에 못 미친다. OpenFHE의 uniform ternary dense 키에서는 `N = 2^13`, 스케일 2^20의 공개키 암호화·복호화만으로 슬롯 최대 오차가 0.10–0.11이었다(4,096개 단위 복소수, 평균 0.024). 논문 LL13의 전체 오차(0.046)보다 크다. HEaaN과 OpenFHE의 키 분포와 잡음 특성 차이로 보이며, 같은 키 분포로 직접 비교하지는 않았다.
- **LL14 정밀도.** 16.0–16.6비트로 논문 16.9비트보다 약간 낮다. boot 소수가 42비트(논문 43비트)이고 OpenFHE의 dense 키 잡음이 남는다.
- **메모리.** 논문 파라미터(dnum 8, mux 12개/인자)의 키는 LL13 약 7.9GB, LL14 약 12GB이다. 최대 RSS는 10GB, 16GB 수준이다. 논문 장비는 512GB RAM이다.
- **병렬성.** 31개 인자 계산과 곱셈 트리, S2C 그룹을 OpenMP로 병렬화했다. 이 Mac은 15코어(성능·효율 혼합)이며 논문은 32~48코어이다.
- **직렬화 미지원.** SHIP 키는 아직 저장·불러오기를 지원하지 않는다.
- **`FIXEDMANUAL`, `HYBRID`, full packing만 지원한다.**
- 사용자 연산의 곱셈 레벨 소수는 `2^p`에 가장 가까운 NTT 소수이다. LL13처럼 작은 소수에서는 사용자의 rescale에 `|q/2^p − 1|` 크기의 스케일 오차가 생긴다(OpenFHE `FIXEDMANUAL`의 일반 특성).

## 재현

```sh
git submodule update --init third-party/cereal third-party/google-test
bash research/ship/run-checks.sh   # WITH_REDUCED_NOISE=ON 빌드, 연구 검사, 라이브러리 검사, SHIP gtest
build-ship/bin/examples/pke/ship-paper-bench LL13 5
build-ship/bin/examples/pke/ship-paper-bench LL14 5
```

다중 코어 측정은 `-DWITH_OPENMP=ON`으로 빌드하고 프로세스마다 `OMP_NUM_THREADS`를 지정한다.
OpenFHE의 내부 병렬 루프는 시작 시점의 `OMP_NUM_THREADS`로 스레드 수를 고정하므로 실행 중 `omp_set_num_threads`로는 바뀌지 않는다.
macOS AppleClang에서는 Homebrew `libomp`를 `OpenMP_*` CMake 변수로 지정했다.
