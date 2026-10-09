# SHIP / OpenFHE 재시작

목표는 SHIP을 OpenFHE에 적용하는 것이다. 현재 결과는 **실행 가능한 정확성 검증용 half-bootstrap**이며, 논문의 최적화 구현이나 안전한 파라미터에서의 전체 부트스트래핑 재현은 아직 아니다.

## 기준점과 원본 보존

- 원래 연구 코드: `7b6fe71b5fae904f82f063eb9e8846591f1186fe`.
- 기반 OpenFHE: `aa391988d354d4360f390f223a90e0d1b98839d7` (1.4.2).
- 작업 브랜치: `research/ship-restart`.
- 원래 연구 변경은 위 커밋과 `main`에 보존되어 있다. 이 브랜치에서는 변경했던 라이브러리 파일 6개를 기반 버전으로 복구하고, 새 실험을 `src/pke/examples/ship/`에 분리했다.
- 기존 `EvalBootstrap` API는 그대로 사용한다. 아직 검증하지 않은 SHIP 경로를 기본 부트스트래핑에 삽입하지 않는다.

## 이번에 재현한 문제

1. 기존 커밋은 AppleClang에서 빌드되지 않는다. `GetNoiseGenerator`, `GetCryptoContext`, `EvalMultExt` 등을 존재하지 않는 객체 API로 호출하고, `HMuxRotKey`, `LiftToPQAndScale`, `cryptoContext`, `keyPair` 등의 정의가 누락되어 있다. 컴파일러는 오류 20개에서 중단했다. 이것이 당시 실패의 전부였다는 뜻은 아니다.
2. `MakeCKKSPackedPlaintext(values, ptLevel)`의 두 번째 인수는 level이 아니라 noise scale degree다. 새 코드는 level을 세 번째 인수로 명시한다.
3. 복소수 위상을 보존하려면 `SetCKKSDataType(COMPLEX)`가 필요하다. 이 설정이 없는 구성에서 허수부가 사라져 정확성 검증이 실패하는 것을 직접 재현했다.
4. `gamma/(4*pi) * exp(i*phase)`에 켤레를 더하면 cosine이 된다. 현재 실수 메시지 복구에는 **`-i * gamma/(4*pi) * exp(i*phase)`**를 사용해 sine을 얻는다. 0 입력 테스트가 이 차이를 검출한다.
5. `KeySwitchCore`는 이미 QP에서 Q로 내릴 때 P로 나누는 처리를 한다. 이후 CKKS 상수 `1/P`를 다시 곱하지 않는다.
6. 일반 암호문 곱셈으로 선택 비트를 곱하면 회전의 각 비트마다 연산 깊이를 소모한다. 새 조건부 회전은 일반 키 스위칭을 조합해 레벨을 유지한다.

## 구현 범위

### 조건부 회전과 블라인드 회전

OpenFHE의 암호문은 `(b,a)`이고 복호화식은 `b+a*s`이다. 설정 단계에서 dense 출력 비밀키 `s` 아래에 다음 두 키 스위칭 키를 만든다.

- 가상 입력 비밀값 `beta`를 위한 키.
- 가상 입력 비밀값 `beta*s`를 위한 키.

평가 단계에서 각각 `b`, `a`에 `KeySwitchCore`를 적용해 더하면 `beta*(b+a*s)`의 암호문을 얻는다. 이어 공개된 후보 회전을 수행한다. 평가 함수는 `beta`나 비밀키를 받지 않는다.

이는 **두 번의 키 스위칭과 별도 회전을 사용하는 정확성 기준 구현**이다. 논문의 fused/hoisted HMuxRot과 비용이 같다는 주장은 하지 않는다. 비트별 두 후보를 합해 blind rotation을 구성한다. 스케일, 레벨, RNS limb 수를 임의로 덮어쓰지 않는다.

### 계수에서 슬롯으로의 half-bootstrap

입력 계약:

- CKKS, 64-bit native backend, `FIXEDMANUAL`, `HYBRID`, `COMPLEX`.
- 한 RNS limb만 남은, sparse 비밀키 아래의 실제 OpenFHE 암호문.
- 계수 메시지가 `round(q0*mu_i/gamma)`로 들어 있으며 스케일은 `q0/gamma`.
- full packing, 실수 메시지, 충분한 출력 깊이.
- 키 생성에 전달한 sparse support/signs와 input key tag는 실제 입력 비밀키에 대응해야 한다.

평가 순서:

1. 입력의 공개된 `a,b` 계수를 바닥 limb에서 읽고 복소수 위상을 인코딩한다.
2. dense 출력키 아래의 암호화된 4종 마스크로 부호와 negacyclic wrap을 선택한다.
3. 암호화된 선택 키로 각 sparse 항을 blind rotate한다.
4. `h+1`개 항을 이진 곱셈 트리로 곱한다. 홀수 개의 마지막 항은 필요한 레벨에 맞춘다.
5. 켤레를 더해 `gamma/(2*pi)*sin(2*pi*mu/gamma)`에 가까운 실수 슬롯을 얻는다.

평가 함수에는 입력/출력 비밀키가 전달되지 않는다. 테스트에서만 마지막 결과를 복호화한다. `ModRaise`, 기존 `EvalMod`, 기존 `CoeffsToSlots`에 출력을 다시 넣지 않는다.

## 검증 결과 (2026-10-09)

환경: macOS arm64, AppleClang 21, OpenFHE 1.4.2, 64-bit backend, OpenMP OFF.

| 검증 | 결과 |
| --- | --- |
| 원본 OpenFHE 부트스트래핑, N=4096, 2048 slots | usable levels 1 -> 10, 최대 오차 `1.60169e-5` |
| 독립적인 negacyclic oracle | 240개 마스크 사례 통과. 모든 위치, 양/음 부호, 두 wrap 경계 포함 |
| HMuxRot 및 blind rotation | 30개 암호문 사례 통과. level 0/3/7, beta 0/1, 좌/우 회전, 0/마지막 인덱스 포함 |
| 조건부 회전 메타데이터 | level, scale, noise scale degree, limb 수 유지. 입력 불변 확인 |
| half-bootstrap, N=1024, h=4, gamma=1024 | 0/교대 부호/난수 입력의 512개 슬롯씩 검사, 입력 1 limb -> 출력 5 limbs |
| half-bootstrap 최대 메시지 오차 | `7.84491e-7` |
| half-bootstrap 최대 sine 목표 오차 | `6.84029e-10` |
| 잘못된 scale/key tag | 예외로 거부 |

수치는 한 실행의 관측치이며 암호화 잡음에 따라 달라진다. 테스트가 허용하는 메시지 오차는 `2e-6`, sine 목표 및 허수부 오차는 `2e-8`이다. 메시지 오차 대부분은 sine 선형 근사의 오차다.

**N=1024, h=4, HEStd_NotSet은 기능 검증용이다. 보안 수준을 제공하지 않는다.** 기존 부트스트래핑은 N=4096이고 입출력 의미도 다르므로, 이 둘의 실행 시간을 속도 향상 수치로 비교하면 안 된다. 이 작업은 새 키 구성의 보안 증명이나 안전한 SHIP 파라미터 선정까지 검증한 것이 아니다.

## 재실행

Git, C++17 compiler, CMake 3.16+와 빌드 도구가 필요하다. 처음 clone할 때 submodule도 가져온다.

```sh
git clone --recurse-submodules https://github.com/coo001/capstone2.git
cd capstone2
git switch -c research/ship-restart 7b6fe71b5fae904f82f063eb9e8846591f1186fe
git apply /path/to/capstone2-restart.patch
bash research/ship/run-checks.sh
```

이미 이 브랜치의 작업 폴더에 있다면 마지막 명령만 실행하면 된다. `CMAKE_BIN`, `SHIP_BUILD_DIR`, `SHIP_JOBS`로 실행 파일/빌드 폴더/동시 컴파일 수를 지정할 수 있다. 시스템 디렉터리에 라이브러리를 설치하지 않는다.

## 다음 연구 단계와 통과 기준

1. **원논문 전문과 키 설계 대조.** 이번에는 ePrint의 논문 정보, 저자 발표 슬라이드, OpenFHE 소스, Poulpy의 독립 구현 설명을 대조했다. ePrint PDF 전문은 이 환경에서 HTTP 403으로 내려받지 못했다. 논문 알고리즘 번호별 검증이 완료되었다고 간주하지 않는다.
2. **fused/hoisted HMuxRot.** 현재의 두 키 스위칭+별도 회전을 대체한다. 현재 회귀 테스트를 동일하게 통과하고, 깊이 유지/키 크기/잡음을 함께 측정해야 한다.
3. **dense -> sparse encapsulation.** 현재 테스트는 sparse 비밀키로 암호화한 계수 입력에서 시작한다. 바닥 modulus에 한정된 encapsulation과 논문에 맞는 sparse 분포를 설계해야 한다.
4. **전체 슬롯 부트스트래핑 연결.** 실제 CKKS 슬롯 입력의 계수 배치, 필요한 두 half의 처리와 후속 선형 변환을 명시해 같은 메시지 의미로 돌아오는지 검증한다. 현재 출력은 first-half coefficient -> real slots이며 일반 `EvalBootstrap`의 대체물이 아니다.
5. **보안/정밀도/성능 비교.** 안전한 ring/modulus/sparse 분포, 동일 입력 범위, 동일 잔여 깊이·오차, 동일 하드웨어·스레드 수에서 키 생성/평가/메모리를 각각 비교한다.

## 근거 자료

- [SHIP, Cheon–Hanrot–Kim–Stehle, EUROCRYPT 2025 / ePrint 2025/784](https://eprint.iacr.org/2025/784)
- [저자 발표 슬라이드](https://iacr.org/submit/files/slides/2025/eurocrypt/eurocrypt2025/447/447_slides.pdf)
- [OpenFHE HYBRID key switching, 기반 커밋](https://github.com/coo001/capstone2/blob/aa391988d354d4360f390f223a90e0d1b98839d7/src/pke/lib/keyswitch/keyswitch-hybrid.cpp)
- [Poulpy SHIP 구현 설명](https://github.com/poulpy-fhe/poulpy/blob/main/docs/ship.md) 및 해당 저장소의 coefficient encoding / mask / bootstrap 구현. 설계와 부호·위상 해석을 교차 확인한 독립 구현이며, 원논문 저자의 구현으로 간주하지 않는다.
