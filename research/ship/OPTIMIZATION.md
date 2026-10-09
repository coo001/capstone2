# OpenFHE 위의 회전 결합과 hoisting

2026-10-09. 이 문서는 **이 저장소의 두 경로 사이**의 최적화 결과를 기록한다.
논문의 SHIP 또는 표준 OpenFHE 부트스트래핑과의 성능 비교가 아니다.
원논문의 전문·알고리즘·보안 증명과의 대조는 아직 완료하지 못했다.

## 변경한 연산

`src/pke/examples/ship/fused-rotation.h`는 다음 연산을 구현한다.

```text
입력: s 아래의 c=(b,a), 숨겨진 비트 β에 대한 평가키, 공개 회전 σ
목표: β·σ(b+a·s) + (1-β)·(b+a·s)를 s 아래에서 평가
```

키 생성 시 회전 분기의 임시 목적 비밀키를 `t = σ⁻¹(s)`로 둔다.
`β`와 `βs`에 대한 두 HYBRID 평가키를 `t` 아래에 생성한다.
평가 단계에서 입력의 `b`와 `a`를 각각 한 번 분해하고, 확장 기저 `QP`에서 두 평가키와 곱해 합친다.
이 결과에 `σ`를 적용하면 비밀키는 `s`로 돌아오고 메시지는 `P·β·σ(b+a·s)`가 된다.
별도로 만든 `(1-β)` 분기의 `QP` 결과를 더한 뒤 마지막에 한 번만 `QP→Q`를 수행한다.

이는 다음 OpenFHE 코드 관례를 따른다.

- `base-leveledshe.cpp::EvalAutomorphismKeyGen`의 역 automorphism 목적 키.
- `keyswitch-hybrid.cpp::EvalKeySwitchPrecomputeCore`의 분해 및 확장 기저 준비.
- `EvalFastKeySwitchCoreExt`의 확장 기저 누적과 `KeySwitchDown`의 `QP→Q` 처리.

키 태그는 **회전 완료 후의 출력 키**를 나타낸다. 평가자는 임시 비밀키나 `β`를 받지 않는다.
중간 `QP` 결과에는 수동으로 `1/P`를 곱하지 않으며, `KeySwitchDown`에서만 축소한다.
CKKS의 근사 기저 변환/반올림 오차 때문에 이전 경로와 ciphertext가 비트 단위로 같을 필요는 없다.
복호화 결과, 독립 평문 목표, level·scale·limb 수를 비교한다.

블라인드 회전의 **비트 한 단계**당 호출 수는 아래와 같다. 0이 아닌 공개 회전 단계 기준이다.

| 연산 | 기존 참조 경로 | 결합 경로 |
| --- | ---: | ---: |
| 입력 성분 분해/확장 준비 | 5 | 2 |
| 확장 기저의 평가키 내적 | 5 | 4 |
| 두 성분에 대한 `QP→Q` 쌍 | 5 | 1 |
| 별도 `EvalRotate` | 1 | 0 |

두 조건부 분기가 **동일한 입력**을 소비하므로 분해 결과 두 개를 함께 사용한다.
다음 비트 단계에서는 입력이 바뀌기 때문에 다시 분해한다.
또한 슬롯 변환의 baby-step 회전들도 동일 입력의 분해를 한 번만 준비하도록 변경했다.
측정한 전체 속도 개선에는 이 두 변경이 함께 포함되어 있다.

## API 및 참조 경로 유지

기존 호출 `MakeFullBootstrapKey(cc, dense, sparse)`는 결합 경로를 사용한다.
네 번째 인자에 `false`를 주면 초기 참조 경로를 사용할 수 있다.
`MakeHalfBootstrapKey`도 마지막 선택 인자로 이를 제어한다.
선택한 한 경로의 평가키만 보관하며 평가 단계에 비밀키를 전달하지 않는다.

두 경로를 같은 컨텍스트에서 생성하는 비교 과정에서 반복 키 생성 문제도 수정했다.
`EvalAutomorphismKeyGen`은 이미 캐시된 인덱스에 대해 새로 생성한 키가 없을 수 있다.
이 반환값을 바로 보관하면 켤레 키 맵이 비어 있을 수 있으므로, 캐시에서 켤레 키를 가져와 별도 한 항목 맵에 보관한다.
전체 검사에서는 결합 경로 설정 후 같은 키로 참조 경로를 생성하며, 벤치마크는 반대 순서도 검사한다.

## 정확성 검사

`research/ship/run-checks.sh`의 다섯 실행 파일이 모두 통과했다.
새 회전 검사는 별도 `EvalRotateKeyGen` 없이 동작하는 사례와 함께 다음 42개 비교를 수행한다.

- 조건부 회전: `β∈{0,1}`, 회전 `0,+3,-3`, level `0,3,7`의 18개 조합.
- 블라인드 회전: 양방향, shift `0,3,128,511`, 같은 세 level의 24개 조합.
- 각 사례에서 이전 경로 및 직접 계산한 평문 목표와 비교하고, 입력 불변성과 metadata 보존을 확인한다.

최종 실행의 회전 최대 오차는 `4.72178e-10`, 두 경로의 최대 차이는 `4.81911e-10`이었다.
512-slot 전체 경로의 최대 메시지 오차는 `2.71444e-7`, 독립 sine/DFT 목표 대비 최대 오차는 `2.22891e-7`이었다.
제곱 후 다시 부트스트래핑한 오차는 `2.66505e-7`이었다. 출력은 여전히 4 limbs다.
큰 입력 `64-32i`에서 약 `1.64515`의 사인 근사 오차가 발생하는 제한도 그대로 확인했다.
정확한 원시 기록은 [fused-checks.txt](results/fused-checks.txt)에 있다.

## 시간과 메모리 측정

환경: Apple M5 Pro, 물리/논리 코어 15개, RAM 48 GiB, macOS arm64,
AppleClang 21, CMake 4.4.4, Release, 64-bit backend, OpenMP OFF.
설정: `N=1024`, 512 complex slots, sparse `h=8`, depth 9, scale 50-bit, 첫 모듈러스 60-bit, `HEStd_NotSet`.
**기능 검사 설정이며 보안 수준이 검증된 파라미터가 아니다.**

같은 프로세스에서 같은 dense/sparse 키·입력 암호문을 사용하고, 각 경로를 1회 준비 실행한 후 7회씩 측정했다.
측정 순서를 매 회차 번갈아 바꿨고 복호화 검사는 타이머 밖에서 수행했다.
키 생성과 선계산 시간은 평가 시간에 포함하지 않았다.

| 측정 항목 | 참조 경로 | 결합 경로 |
| --- | ---: | ---: |
| 전체 평가 중앙값 | 472.839 ms | 221.610 ms |
| 관측 범위 | 467.609–604.333 ms | 219.086–276.477 ms |
| 7회 평가의 최대 메시지 오차 | 5.66118e-7 | 2.07378e-7 |
| 출력 모듈러스 수 | 4 | 4 |

중앙값 비율은 **2.13365배**다. 두 경로 최종 결과의 최대 차이는 `6.45905e-7`이었다.
일부 실행 시간에 변동이 있으므로 이 비율을 모든 장비·크기에서 보장하지 않는다.
키/암호화 잡음도 무작위여서 수치 오차는 실행마다 변한다.

메모리는 각 경로를 **독립 프로세스**에서 실행해 macOS `/usr/bin/time -l`로 측정했다.
두 경로의 키와 회전 캐시를 동시에 보관하는 비교 프로세스의 메모리는 비교에 사용하지 않았다.

| 측정 항목 | 참조 경로 | 결합 경로 |
| --- | ---: | ---: |
| 최대 resident set size | 280,952,832 bytes | 278,790,144 bytes |
| 선택한 저장 계수 배열의 payload | 260,308,992 bytes | 256,868,352 bytes |

최대 RSS는 약 268→266 MiB로 감소했지만 차이는 **약 0.8%**에 불과하다.
payload는 encapsulation 키, selector 암호문, 조건부 회전/켤레/회전 평가키, 슬롯 변환 plaintext의 계수 배열만 센다.
컨테이너·allocator·permutation map·컨텍스트 표·비밀키·공개키·relinearization key는 이 payload에서 제외되며 RSS에는 포함될 수 있다.
시간을 줄였어도 많은 평가키와 모든 슬롯 변환 대각선의 저장 비용은 여전히 남는다.
원시 표본과 독립 프로세스 측정값은 [optimization-benchmark.txt](results/optimization-benchmark.txt)에 있다.

## 재현 명령

```sh
bash research/ship/run-checks.sh
cmake --build build-ship --target ship-optimization-bench --parallel 4
build-ship/bin/examples/pke/ship-optimization-bench compare
# macOS의 독립 프로세스 최대 RSS 측정
/usr/bin/time -l build-ship/bin/examples/pke/ship-optimization-bench reference
/usr/bin/time -l build-ship/bin/examples/pke/ship-optimization-bench fused
```

다른 CMake 실행 파일이나 빌드 디렉터리를 사용했다면 해당 경로로 바꾼다.
Linux의 `/usr/bin/time`은 옵션과 RSS 단위가 다르므로 macOS의 bytes 값과 직접 혼용하지 않는다.

## 다음 연구의 범위

이 변경으로 OpenFHE 위의 회전 결합과 hoisting의 기능·속도 효과는 확인했다.
그러나 논문 HMuxRot의 정확한 재현, 평가키/순환 보안 가정 분석, sparse-secret 보안 추정은 별개로 남아 있다.
큰 차수에서의 메모리 측정과 슬롯 변환의 구조적 최적화가 필요하며,
안전한 파라미터를 정하기 전에는 표준 OpenFHE나 논문 수치에 대한 속도 우위를 주장하지 않는다.
