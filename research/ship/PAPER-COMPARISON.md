# 원논문과 구현의 대조

2026-10-10. 대조 기준은 [ePrint 2025/784](https://eprint.iacr.org/2025/784.pdf) 전문(33쪽, EUROCRYPT 2025 출판본의 minor revision)이다.
페이지 번호는 이 ePrint PDF 기준이며 Springer 출판본(31쪽)과 다르다.
처음 대조(2026-10-10)의 코드 기준은 `research/ship-aux-masking` 브랜치의 연구 프로토타입이다.
**2026-10-11 라이브러리 구현에서 아래 "논문과 다른 부분" 1~6 중 1~5와 6의 B-to-1 mux를 구현했다.** 대응표는 [LIBRARY.md](LIBRARY.md)에 있다.
**2026-10-11 추가.** HT14·HT15 프리셋을 디스크 기반 키로 실행했다(정밀도 HT14 4.0비트 / 논문 4.55, HT15 18.3–18.4비트 / 논문 18.1, 남는 레벨 9로 동일).
논문 §5.2(p.24)의 sparse 키 128비트 경계값 `(13, 55)`, `(14, 100)`은 현재 lattice estimator의 혼합 공격으로 2^121.4, 2^126.2였다. `(15, 105)`는 2^148.7이었다.
라이브러리 기본값은 `(13, 42)`, `(14, 88)`, `(15, 105)`로 바꿨고 정밀도 변화는 없었다. 자세한 내용은 [LIBRARY.md](LIBRARY.md#보안-검증)에 있다.

## 논문과 일치하는 부분

| 논문 | 코드 | 판정 |
| --- | --- | --- |
| Algorithm 1 (p.14): 마스크 → BlindRot → product tree → `ct + Conj(ct)` | `half-bootstrap.h`의 `HalfBootstrap` | 일치 |
| `pt0 = Ecd(γ/(4ιπ)·ω^b)`, `pt1..pt4 = ω^{±a}`의 앞·뒤 절반 | 같은 함수의 `initial`, `phases` | 일치. 코드의 `-i·γ/(4π)`는 `γ/(4ιπ)`와 같다 |
| Lemma 1 (p.16)의 마스크와 회전 전 형태 Eq. (5) | `reference.h`의 `Masks` | 일치. 독립 negacyclic oracle 240건으로 검사한다 |
| Definition 1·Lemma 2 (p.19), Algorithm 5 (p.23)의 HMuxRot | `fused-rotation.h` | 평문 의미는 일치한다. 임시 목적 키 `σ⁻¹(s)`로 만든 키에 `σ`를 적용하면 p.23의 키 관계식과 같아진다 |
| §4.1·§4.4 (p.18, p.21–22): 마스크 곱셈을 보조 모듈러스 `P`에서 수행하고 `Rescale_P` | `EncryptMaskOverQP`, `HalfBootstrap`의 `auxMasking` 경로 | 일치. 2026-10-10에 구현했다. [AUX-MASKING.md](AUX-MASKING.md) 참조 |
| §3.2 (p.15)·§5.2 (p.24): S2C는 대각선 방식과 BSGS로 1레벨 | `full-bootstrap.h`의 `MakePackingPlan`, `ApplyPacking` | 방식은 일치. 위치는 아래 차이 4 참조 |

HMuxRot의 잡음까지 같다고 주장하지 않는다.
OpenFHE의 `ApproxSwitchCRTBasis`·`ApproxModDown`은 논문의 이상적인 lift·반올림과 같은 구현이 아니다.

## 논문과 다른 부분 (연구 프로토타입 기준, 라이브러리에서 해소 여부 표시)

라이브러리 구현에서의 상태: 1 구현(θ, base B), 2 구현(window w), 3 구현(`realOnly`), 4 구현(S2C를 낮은 모듈러스에서 먼저),
5 구현(`q0·p'` special prime), 6 B-to-1 mux와 단계별 hoisting 구현.

1. **블라인드 회전 방식.** 논문은 column 방식과 mux 방식을 θ로 섞고(Algorithm 4, p.22), mux는 base `B=4`로 분해한다(Table 2, p.25).
   코드는 mux 방식만 쓰며 base 2로 `log2(N/2)`비트를 모두 회전한다.
   OpenFHE의 `HYBRID`는 키 스위칭 방식 이름이며 논문의 `BRotHybrid`와 다르다.
2. **비밀키 분포.** 논문은 0이 아닌 계수의 위치가 `[kN/h−w, kN/h+w)`에 있다고 가정해 회전 범위를 `2w`로 줄인다(§5.1, p.24).
   보안은 이 분포에 대해 [May21] 공격 비용으로 따로 제시한다.
   코드는 무작위 sparse ternary 키를 쓰므로 이 최적화와 그 보안 근거가 모두 없다.
3. **실수와 복소수.** 논문의 주요 실험은 실수 벡터이다(Table 3).
   코드는 복소 슬롯을 위해 half-bootstrap을 두 번 수행하므로 블라인드 회전 수가 두 배이다.
4. **슬롯 변환의 위치.** 논문은 낮은 모듈러스에서 S2C를 먼저 하고 sparse 키로 전환한 뒤 SHIP을 호출한다(§3.2, p.15; §5.2, p.24).
   코드는 일반 슬롯 암호문을 받아 SHIP **뒤에** 높은 모듈러스에서 DFT를 수행한다.
   수학적으로는 둘 다 맞지만 큰 파라미터에서는 비용이 달라진다.
5. **dense→sparse 키 전환.** 논문은 전환 키의 보안을 `P·q0`에서 평가한다(p.24).
   코드는 `q0`에서 8비트 digit의 BV 방식을 쓴다. 보안 분석은 코드의 구성에 맞춰 다시 해야 한다.
6. **B-to-1 mux-rotate와 hoisting(§5.1, p.23–24)은 구현하지 않았다.**
   코드의 hoisting은 같은 비트 단계에서 두 분기가 분해를 공유하는 범위에 한정된다.

## 깊이

| 항목 | 논문 (h=31) | 코드, rescale 경로 | 코드, 보조 모듈러스 경로 |
| --- | ---: | ---: | ---: |
| 마스크 | 0 | 1 | 0 |
| product tree `⌈log2(h+1)⌉` | 5 | 5 | 5 |
| S2C/DFT | 1 | 1 | 1 |
| 합계 | 6 | 7 | **6** |

`h=31` 기준으로 코드의 보조 모듈러스 경로가 6레벨을 쓰는 것을 `ship-aux-masking-checks`에서 확인했다.
위치는 다르므로(차이 4) 레벨이 놓이는 모듈러스 크기까지 같다는 뜻은 아니다.

## 논문 실험 조건 (비교할 때 맞출 기준)

- LL13: `N=2^13`, `log PQ=218`, `h=31`, θ=6, `ℓ0/θ=64`, `B=4`, 실수 4.45비트, 곱셈 레벨 1 (Table 2–3, p.25–26).
- 1코어 3.02s, 32코어 215ms, 128비트 보안, 실패 확률 `≤2^-128`.
- HEaaN, AVX2, OpenMP, AMD EPYC 7473X 두 개, 100회 중앙값 (§5.4, p.25).
- 정밀도가 다른 설정끼리 지연 시간만 비교하면 안 된다(예: 4.45비트 LL13 대 16.7비트 Param15).
