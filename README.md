# crsless-go

[SNUCP/crsless](https://github.com/SNUCP/crsless) (Multi-key FHE with Non-Interactive Setup in the Plain Model,
https://eprint.iacr.org/2026/322) 의 Julia 구현을 Go + [Lattigo v6.2.0](https://github.com/tuneinsight/lattigo) 의
`ring` 패키지 위로 옮긴 것입니다. HIENAA 대신 Lattigo 의 RNS 다항식·NTT·샘플러를 씁니다.

## 실행 (Go 1.24 이상)

```
go mod tidy
go run ./cmd/crsless -test 2              # main.jl (K=2)
go run ./cmd/crsless -test 4 -distinct    # K=4, 내부곱을 서로 다른 당사자끼리
go run ./cmd/paramfinder                  # paramfinder.jl
go test ./...                             # K=2 테스트 4개
```

## 파일 대응

| 원본 (Julia) | Go |
|---|---|
| `src/parameters.jl` | `mkfhe/parameters.go` (`MKParamSketch`, `MKParameters`, `FindPrimePow2`) |
| `src/entity.jl` | `mkfhe/entity.go` (`Tensor`, `CGSW`, `EvalKey`, `MLEV`, `MGSW`, `Shrink`) |
| `src/key.jl` | `mkfhe/key.go` (`CGSWKey`, `TernaryCGSWKey`) |
| `src/operator.jl` | `mkfhe/operator.go` (`MKOperator`, `DecomposeTo`, `VectorExtProd`, `HoistedVectorExtProd`, `MGSWExtend`, `HoistedExternalProduct`, `InternalProduct`) |
| `src/encryptor.jl` | `mkfhe/encryptor.go` (`MKEncryptor`, `SalweSampleTo`, `Phase`, `CGSWEncrypt*`, `EvkEncrypt`, `JointDecryptor`) |
| `main.jl` | `mkfhe/tests.go` + `cmd/crsless/main.go` |
| `paramfinder.jl` | `cmd/paramfinder/main.go` |
| HIENAA `Decomposer(Q, 1)` | `MKOperator.Gvec` + `DecomposeTo` (RNS 성분별 중심 표현 확장) |
| HIENAA `SimpleScaler` | `mkfhe/scale.go` (big.Int 로 단순 구현, 검증용) |

## 원본과 다른 점

- 모든 다항식을 항상 NTT 형식으로 보관합니다 (원본은 `isntt` 플래그로 형식을 추적). 분해할 때만 INTT 합니다.
- 당사자 번호 `idx` 는 원본처럼 1부터 시작합니다 (`MGSWExtend(res, ct, idx, evks)`).
- 원본 `rgsw_internal_product_test` 는 `while idx1 != idx2` 때문에 항상 같은 당사자끼리 곱합니다.
  기본값은 원본 동작을 따르고, `-distinct` 를 주면 서로 다른 당사자끼리 곱합니다.
- 소수 선택은 HIENAA `find_prime` 과 같은 알고리즘이라 같은 Q 가 나옵니다 (K=2: 114689, 147457, 163841, 188417, 270337, 286721).

## 메모리

평가키가 당사자당 K 개의 cGSW 이고 각 cGSW 가 glen·(K+1)² 개 다항식이라, K 가 커지면 메모리가 급격히 늘어납니다
(원본도 같은 구조). 대략 K=4 는 1GB 미만, K=8 은 10GB 이상이 필요합니다.

## 평가키 압축 (신규, 원본에 없음)

기존 KeyGen 이 만든 평가키 `G_h = cGSW_t(t_h)` (당사자당 K(K+1)²d 환 원소) 에서
`ek_comp = (U_0, H_⋆, H_Δ, V_1..V_{K-1})` ((5K+1)d 환 원소) 를 추출하고, 서버가 이것만으로 `Ĝ_h` 전체를 복원합니다.
기존 코드는 수정하지 않았습니다. 이번 단계는 통신량 압축의 프로토타입이며 (full 키 생성 후 압축), 클라이언트 KeyGen 비용은 줄이지 않습니다.

| 파일 | 내용 |
|---|---|
| `mkfhe/compress.go` | `CompressEvalKey`, `UpKeyStep` (= R), `UpKeys`, `ExtendCLEV`, `DecompressEvalKey`, 직렬화 |
| `mkfhe/compress_check.go` | 잔차 측정, σ=0 암호화기, 통합 시나리오 `EvalKeyCompressionTest` |
| `mkfhe/compress_test.go` | 단계별 테스트 (σ=0 정확 등식 → σ=3.2 잔차 한계), 벤치마크 |
| `cmd/evkcompress` | `go run ./cmd/evkcompress -test 4 -depth 12` |

잔차 판정 기준은 가젯 자리 복호 정확성 조건 `|e| < Q/(2·max q_j)` 입니다. 정확성 검증일 뿐 보안 증명이 아닙니다.

**알려진 한계**: 복원 키 Ĝ_h 의 블록 1..K−1 은 가젯 곱이 두 번 (R 단계의 `H_Δ D_g(A)`, Extension 의 `U_i D_g(a)`) 쌓여
노이즈가 신선한 키보다 훨씬 큽니다 (K=2: 2^4 → 2^51, K=4: 2^4 → 2^103). 기존 파라미터 (dlen=1, 자릿수 ≈ q_i/2) 에서
MGSWExtend 는 통과하지만 (여유 K=2 약 7비트, K=4 약 5비트) 그 위의 InternalProduct 는 노이즈 예산을 넘어 복호에 실패합니다.
σ=0 에서는 InternalProduct 까지 정확하므로 대수적 구성은 맞고, 해결하려면 더 작은 가젯 자릿수 또는 더 큰 Q 가 필요합니다.
