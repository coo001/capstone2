# SHIP → OpenFHE 연구 기록

**2026-10-11 현재 SHIP은 OpenFHE 라이브러리 기능이다.** `cc->EvalSHIPBootstrapKeyGen`, `cc->EvalSHIPBootstrap`으로 호출하며
구현은 `src/pke/lib/scheme/ckksrns/ckksrns-ship.cpp`에 있다. 논문의 블라인드 회전(column+mux, base-4), 균등 간격 sparse 키,
`P·q0` encapsulation, 대각선 S2C, Table 2의 역할별 소수를 구현했고 논문의 네 파라미터 세트(LL13, LL14, HT14, HT15)를 128비트 HE 표준 안에서 실행했다.
키 직렬화(시드 저장으로 크기 1/2), 디스크 기반 키, 입력 범위 지정(`messageBound`)을 지원한다.
lattice estimator 직접 실행 결과 논문의 sparse 키 경계값 두 개(2^13·55비트, 2^14·100비트)가 128비트에 못 미쳐, 기본값을 estimator 기준(42, 88, 105)으로 바꿨다.
설계·보안 검증·결과는 [LIBRARY.md](LIBRARY.md), 원논문 대조는 [PAPER-COMPARISON.md](PAPER-COMPARISON.md)에 있다.

아래는 그 이전 단계의 연구 프로토타입(`src/pke/examples/ship/*.h`, `ship::FullBootstrap`) 기록이다.
라이브러리 구현과 별도로 비교·기록용으로 남겨 두었다.

이전 프로토타입은 **일반 CKKS 슬롯 암호문을 복원하는 전체 경로의 정확성 참조 구현**이었다.
half-bootstrap에서 dense→sparse 키 전환과 양쪽 계수 복구, 슬롯 변환까지 연결했다.

프로토타입에는 **회전 결합·분해 재사용·지연 ModDown** 최적화가 적용되어 있다.
같은 키·입력의 512-slot 실험에서 기존 참조 경로 대비 중앙값 기준 약 2.13배 빨라졌다.
메모리는 거의 줄지 않았다. 자세한 구성과 측정 한계는 [최적화 기록](OPTIMIZATION.md)에 있다.

2026-10-10에 **마스크 곱셈을 보조 모듈러스 `P`에서 수행하는 경로**를 기본으로 바꿨다(논문 §4.4).
부트스트래핑 깊이가 1 줄어 `h=31`에서 논문과 같은 6레벨을 쓴다. 자세한 내용은 [AUX-MASKING.md](AUX-MASKING.md)에 있다.

## 재현

OpenFHE 의존성이 준비된 저장소에서 실행한다. 새 clone에는 submodule도 받아야 한다.
스크립트는 `WITH_REDUCED_NOISE=ON`으로 빌드하며 연구 프로토타입 검사, 라이브러리 검사, SHIP gtest를 실행한다.

```sh
git submodule update --init --recursive
bash research/ship/run-checks.sh
```

`CMAKE_BIN`, `SHIP_BUILD_DIR`, `SHIP_JOBS`로 CMake 실행 파일, 빌드 경로, 병렬 작업 수를 지정할 수 있다.
검증 환경은 macOS arm64 / AppleClang 21 / CMake 4.4.4 / 64-bit backend / OpenMP OFF다.
스크립트는 아래 다섯 검증 실행 파일을 빌드하고, 오류 발생 시 0이 아닌 상태로 종료한다.
외부의 전체 OpenFHE unit test suite를 실행했다는 의미는 아니다.

| 검증 실행 파일 | 확인하는 내용 |
| --- | --- |
| `ship-baseline-checks` | 표준 OpenFHE 부트스트래핑의 정확성과 연산 여유 복구 |
| `ship-component-checks` | 독립 negacyclic 마스크 oracle 240개, 조건부/블라인드 회전 30개 |
| `ship-fused-rotation-checks` | 별도 회전 키 없는 평가, 기존 경로 및 평문 oracle과 42개 비교, metadata/입력 보존 |
| `ship-half-bootstrap-checks` | 실제 sparse-key 암호문의 첫 계수 절반을 슬롯으로 복구, sine 목표 및 메시지 오차 |
| `ship-full-bootstrap-checks` | 일반 복소 CKKS 입력, dense→sparse 전환, 전체 슬롯 복원, 제곱 및 재부트스트래핑 |

## 전체 경로

핵심 코드는 `src/pke/examples/ship/full-bootstrap.h`에 있다.

1. OpenFHE에서 보통 방법으로 인코딩·암호화한, 마지막 모듈러스의 dense-key 암호문을 받는다.
2. 바닥 모듈러스 `q0`에서만 BV 방식의 자릿수 분해 키 전환을 수행한다. 공개 키 전환 재료는 sparse key 아래의 `B^d * s_dense` 암호화이며, 전체 `Q` 또는 `QP`로 확장하지 않는다.
3. sparse-key 암호문에 half-bootstrap을 수행해 앞쪽 계수 절반을 복구한다.
4. `X^(-N/2)`를 곱한 암호문에도 half-bootstrap을 수행해 뒤쪽 절반을 복구한다.
5. 두 결과를 `x + i*y`로 결합하고, OpenFHE의 슬롯 순서에 맞는 동형 선형 변환을 수행한다. full packing에서 `X^(N/2)`는 모든 슬롯에 `i`를 곱한다.
6. 입력과 같은 dense key로 복호화 가능한 CKKS 슬롯 암호문을 반환한다. 남은 모듈러스로 후속 곱셈을 수행할 수 있다.

슬롯 변환 행렬은 `V[k,j] = exp(2πi * 5^k * j / (2N))`이며 baby-step/giant-step으로 적용한다.
현재는 행렬 대각선을 모두 저장하는 참조 구현이다. 저장량과 plaintext 곱셈 비용이 커서 안전한 큰 차수에서의 성능을 대표하지 않는다.

비밀키는 **키 생성과 테스트의 복호화/oracle에만** 사용한다.
평가 함수 `FullBootstrap(cc, input, keys)`는 비밀키를 인자로 받지 않는다.
테스트의 `SineOracle`은 복호화한 계수와 직접 계산한 DFT를 사용해 생산 코드의 슬롯 변환과 독립적으로 결과를 검증한다.

## 지원 조건과 수치적 한계

- `FIXEDMANUAL`, `HYBRID`, `COMPLEX`, full packing (`N/2` slots).
- 입력은 두 성분, 한 RNS limb이며 rescale을 마친 상태(`noiseScaleDeg=1`)여야 한다.
- 출력 키는 보통의 OpenFHE dense key다. 중간 sparse key는 ternary이며 테스트에서는 OpenFHE의 무작위 sparse sampler를 사용한다.
- half-bootstrap의 계수별 목표는 `f(μ) = γ/(2π) * sin(2πμ/γ)`이고 `γ = q0 / inputScale`이다.
- 따라서 입력 계수 크기가 `γ`에 비해 작아야 원래 메시지에 가까워진다. 모든 입력 범위에서 항등 복원을 보장하지 않는다.
- 사인 근사만의 계수 오차는 `|f(μ)-μ| ≤ 2π²|μ|³/(3γ²)`로 상계할 수 있다. 실제 슬롯 오차에는 키 전환 잡음, CKKS 계산 오차와 슬롯 변환에 따른 누적도 포함된다.

잘못된 key tag, 부분 packing, rescale 전 입력, 0 scale, 바닥보다 높은 level은 테스트에서 거부되는지 확인한다.
메시지 범위는 암호화되어 있으므로 평가 함수가 자동으로 검사하지 않는다.

## 실험 결과

아래는 2026-10-09 **최적화 전** 전체 경로의 검증 관측값이다. 암호 키와 잡음은 매번 무작위이므로 실행마다 변한다.
당시 기록은 `results/full-checks.txt`, 최적화 후 통합 검사 기록은 `results/fused-checks.txt`를 참조한다.

전체 경로는 depth 9, scale 50-bit, 첫 모듈러스 60-bit (`γ ≈ 1024`)에서 검사했다.
각 설정에서 영벡터·복소 상수·마지막 슬롯 impulse·주기 패턴·무작위 복소 입력을 사용했다.
일반 입력은 실수부와 허수부 각각 `[-0.5, 0.5]` 범위다.

| 설정 | 작은 입력 최대 슬롯 오차 | 독립 sine/DFT oracle 최대 오차 | 모듈러스 수 | 제곱 후 재부트스트래핑 오차 |
| --- | ---: | ---: | ---: | ---: |
| `N=128`, 64 slots, `h=4` | `1.06548e-7` | `1.53048e-8` | 1 → 5 | `3.86186e-8` |
| `N=1024`, 512 slots, `h=8` | `3.44761e-7` | `3.33802e-7` | 1 → 4 | `6.07464e-7` |

모듈러스 수는 사용 가능한 곱셈 횟수와 동일한 단위가 아니다. 실제 후속 제곱과 rescale, 다시 바닥으로 낮춘 뒤의 재부트스트래핑을 별도로 검사했다.

추가로 모든 슬롯에 `64-32i`를 넣으면 메시지 복원 오차가 약 **1.64515**로 커진다.
동시에 독립 sine/DFT 목표에 대한 오차는 `3.51e-7` 미만이다.
이는 사인 근사의 입력 범위 제한을 드러내는 의도적인 검사이며, 큰 메시지도 정확히 복원된다는 주장에 대한 반례다.

최종 실행의 전체 평가 시간은 `N=128`에서 약 23–26ms, `N=1024`에서 약 475–497ms였다.
키 생성·선계산을 제외한 단일 장비의 참고 수치다. 안전성이 없는 작은 설정이고 baseline과 조건도 달라서, 논문 또는 표준 부트스트래핑에 대한 속도 우위로 해석할 수 없다.

## 원래 실패와 첫 복구

원래 연구 코드는 커밋 `7b6fe71b5fae904f82f063eb9e8846591f1186fe`에 보존된다.
기존 수정에는 정의되지 않은 함수와 OpenFHE API 불일치가 있어 빌드되지 않았다.
라이브러리/표준 예제 파일 여섯 개를 부모 커밋 `aa391988d354d4360f390f223a90e0d1b98839d7`의 OpenFHE 1.4.2 기준으로 복구했다.

이후 참조 구현에서 확인·수정한 문제는 다음과 같다.

- `MakeCKKSPackedPlaintext(values, ...)`의 두 번째 인자는 level이 아니라 noise scale degree다.
- 실수 인코딩으로는 필요한 복소 위상을 보존할 수 없어 `COMPLEX`를 명시했다.
- half-bootstrap의 초기값에 `-i`가 없으면 목표인 사인 대신 코사인 성분을 복구한다.
- `KeySwitchCore`는 이미 `PQ→Q` 처리를 하므로 `1/P`를 다시 적용하면 안 된다.
- 조건부 회전을 단순 ciphertext 곱셈으로 대체하면 추가 깊이가 필요하다. 초기 참조 경로는 두 번의 key switch와 별도 회전을 사용하며, 현재는 이를 결합한 경로도 제공한다.
- product tree의 피연산자 level을 명시적으로 정렬하고, `FIXEDMANUAL`에서 곱셈 뒤 rescale한다.

첫 복구 실행에서 표준 부트스트래핑은 2048 slots / usable levels 1→10 / 최대 오차 `1.60169e-5`를 보였다.
half-bootstrap은 512 slots × 3종 입력 / 1→5 limbs / 최대 메시지 오차 `7.84491e-7`를 보였다.
최종 스크립트는 이 경로들도 함께 회귀 검사한다.

## 아직 완료하지 않은 연구

1. **원논문과 알고리즘별 대조:** 2026-10-10에 ePrint 전문으로 대조했다([PAPER-COMPARISON.md](PAPER-COMPARISON.md)). Algorithm 1·Lemma 1·HMuxRot의 평문 의미와 §4.4의 레벨 절약은 일치한다. column/mux 혼합 블라인드 회전(Algorithm 4), base-4 mux, 균등 간격 sparse 키, 실수 전용 경로, 논문 순서의 S2C는 아직 구현하지 않았다. 현재 코드를 논문의 완전한 재현이라고 부르지 않는다.
2. **잡음과 키 전환의 동일성:** OpenFHE의 근사 기저 변환·ModDown은 논문의 이상적 연산과 다르다. bottom key switch도 논문이 보안을 평가한 `P·q0` 구성과 다르다.
3. **보안 검증:** `HEStd_NotSet`, `N=128/1024`, `h=4/8`은 기능 검사 전용이다. sparse secret 분포, 공개 평가키와 키 전환의 보안 가정, modulus/잡음 예산을 포함한 분석이 필요하다.
4. **확장 가능한 슬롯 변환과 벤치마크:** 현재 대각선 방식의 큰 메모리/연산 비용을 줄이고, 같은 보안 수준·정밀도·입력 범위에서 OpenFHE baseline과 비교해야 한다.
5. **일반화:** 부분 packing, 자동 scale 기법, 더 넓은 메시지 범위, 대규모 반복 실행에 대한 검증이 남아 있다.

다음 연구는 여전히 큰 평가키 및 슬롯 변환 메모리를 줄이고, 논문 대조와 보안 파라미터 분석을 거쳐 동일 조건 비교를 수행하는 순서로 진행한다.

## 참고 자료

- [SHIP 논문 서지정보 — IACR ePrint 2025/784](https://eprint.iacr.org/2025/784)
- [EUROCRYPT 2025 저자 발표 슬라이드](https://iacr.org/submit/files/slides/2025/eurocrypt/eurocrypt2025/447/447_slides.pdf)
- [Poulpy의 독립 SHIP 구현 설명](https://github.com/poulpy-fhe/poulpy/blob/main/docs/ship.md): 원논문 구현과 동일한 출처로 취급하지 않았다.
- OpenFHE의 `ckkspackedencoding.cpp`, `dftransform.cpp`, `base-pke.cpp`: 인코딩/슬롯 순서와 암호문 부호 관례 확인.
