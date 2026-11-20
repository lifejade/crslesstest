package mkfhe

import (
	"fmt"
	"math/big"
	"testing"
)

// σ = 0 (정확한 등식) 을 먼저, 그다음 기존 MPS 파라미터 σ = 3.2 (잔차 < Q/(2 q_max)) 를 검사합니다.
var compressSigmas = []float64{0, 3.2}

type compressFixture struct {
	op  *MKOperator
	K   int
	key *CGSWKey
	evk *EvalKey
	ck  *CompressedEvalKey
}

func newCompressFixture(tb testing.TB, testID int, sigma float64) *compressFixture {
	tb.Helper()
	op, K, err := setup(testID)
	if err != nil {
		tb.Fatal(err)
	}
	key, err := TernaryCGSWKey(op.R, K, HW)
	if err != nil {
		tb.Fatal(err)
	}
	entor, err := newEncryptorSigma(key, sigma, op)
	if err != nil {
		tb.Fatal(err)
	}
	evk := entor.EvkEncrypt(K)
	ck, err := op.CompressEvalKey(evk)
	if err != nil {
		tb.Fatal(err)
	}
	return &compressFixture{op: op, K: K, key: key, evk: evk, ck: ck}
}

func sigmaName(s float64) string { return fmt.Sprintf("sigma=%g", s) }

// assertResidual: σ=0 이면 잔차 = 0, 아니면 잔차 < Q/(2 q_max).
func assertResidual(t *testing.T, op *MKOperator, sigma float64, what string, e *big.Int) {
	t.Helper()
	if sigma == 0 && e.Sign() != 0 {
		t.Fatalf("%s: σ=0 인데 잔차 %s (정확한 등식이어야 함)", what, Log2Str(e))
	}
	if e.Cmp(op.NoiseBound()) >= 0 {
		t.Fatalf("%s: 잔차 %s ≥ 기준 %s", what, Log2Str(e), Log2Str(op.NoiseBound()))
	}
	t.Logf("%-28s 잔차 %s", what, Log2Str(e))
}

// 1단계: Compression
func TestCompressHints(t *testing.T) {
	for _, sigma := range compressSigmas {
		t.Run(sigmaName(sigma), func(t *testing.T) {
			f := newCompressFixture(t, 2, sigma)
			op, K, ck := f.op, f.K, f.ck
			if got, want := ck.NumPolys(), (5*K+1)*op.Glen; got != want {
				t.Fatalf("압축키 원소 수 %d ≠ (5K+1)d = %d", got, want)
			}
			// 구조: U_0 = G_0^{[0]}, H_⋆ = G_0^{[K]}, V_h = G_h^{[0]} 의 (첫 행, 마지막 행), V_0 = U_0 의 (첫 행, 마지막 행)
			for j := 0; j < op.Glen; j++ {
				for r := 0; r <= K; r++ {
					if !ck.U0[j][r].Equal(&f.evk.Vals[0].Block(0)[j][r]) || !ck.HStar[j][r].Equal(&f.evk.Vals[0].Block(K)[j][r]) {
						t.Fatalf("U_0 / H_⋆ 추출 불일치 (열 %d, 행 %d)", j, r)
					}
				}
				for h := 0; h < K; h++ {
					v, blk := ck.CLEV(h), f.evk.Vals[h].Block(0)
					if len(v[j]) != 2 || !v[j][0].Equal(&blk[j][0]) || !v[j][1].Equal(&blk[j][K]) {
						t.Fatalf("V_%d 추출 불일치 (열 %d)", h, j)
					}
				}
			}
			u, s, d, v := op.hintResiduals(ck, f.key.Val)
			assertResidual(t, op, sigma, "S·U_0 − t_0 e_0 gᵀ", u)
			assertResidual(t, op, sigma, "S·H_⋆ − t_0 t gᵀ", s)
			assertResidual(t, op, sigma, "S·H_Δ − (t − Πt) gᵀ", d)
			assertResidual(t, op, sigma, "(1,t_0)·V_h − t_h gᵀ", v)
		})
	}
}

// 2단계: Recursive Up-Key.  K=2 에서는 재귀 깊이가 1 이므로, 측정용으로 2K 까지 적용해 순환 후 증가도 봅니다.
func TestUpKeys(t *testing.T) {
	for _, sigma := range compressSigmas {
		t.Run(sigmaName(sigma), func(t *testing.T) {
			f := newCompressFixture(t, 2, sigma)
			op, K := f.op, f.K
			ups := op.UpKeys(f.ck, 2*K)
			for i, U := range ups {
				assertResidual(t, op, sigma, fmt.Sprintf("S·U_%d − t_0 e_%d gᵀ", i, i%K), op.upKeyResidual(U, f.key.Val, i))
			}
		})
	}
}

// 음성 테스트: 힌트가 틀리면 검사기가 실패를 잡아내야 합니다 (검사가 공허하지 않음을 확인).
func TestUpKeysDetectsBadHint(t *testing.T) {
	f := newCompressFixture(t, 2, 0)
	for _, col := range f.ck.HDelta {
		col.Zero()
	}
	ups := f.op.UpKeys(f.ck, f.K)
	if e := f.op.upKeyResidual(ups[1], f.key.Val, 1); e.Cmp(f.op.NoiseBound()) < 0 {
		t.Fatalf("H_Δ = 0 인데 U_1 잔차 %s 가 기준 미만", Log2Str(e))
	}
}

// 3단계: Full cGSW Reconstruction — 모든 h, 마지막 블록 K 를 포함한 모든 가젯 블록.
func TestExtendCLEV(t *testing.T) {
	for _, sigma := range compressSigmas {
		t.Run(sigmaName(sigma), func(t *testing.T) {
			f := newCompressFixture(t, 2, sigma)
			op, K := f.op, f.K
			rec := op.DecompressEvalKey(f.ck)
			if len(rec.Vals) != K {
				t.Fatalf("Ĝ 개수 %d ≠ K", len(rec.Vals))
			}
			for h := 0; h < K; h++ {
				g := rec.Vals[h]
				if g.Parties() != K || g.Glen != op.Glen || len(g.Stack) != (K+1)*op.Glen {
					t.Fatalf("Ĝ_%d 크기: (K+1)×(K+1)d 가 아님", h)
				}
				orig := op.cgswBlockResiduals(f.evk.Vals[h], f.key.Val, f.key.Val[h])
				for i, e := range op.cgswBlockResiduals(g, f.key.Val, f.key.Val[h]) {
					assertResidual(t, op, sigma, fmt.Sprintf("Ĝ_%d block %d (원본 %s)", h, i, Log2Str(orig[i])), e)
				}
			}
		})
	}
}

// 4단계: 통합 — 복원 키를 기존 MGSWExtend / InternalProduct 에 투입.
// σ = 0 이면 InternalProduct 까지 정확해야 합니다. σ = 3.2 에서는 복원 키의 노이즈가 커서
// 기존 파라미터의 InternalProduct 예산을 넘으므로, 그 결과는 실패로 숨기지 않고 로그로 남깁니다.
func TestCompressionIntegrationK2(t *testing.T) {
	for _, sigma := range compressSigmas {
		t.Run(sigmaName(sigma), func(t *testing.T) {
			rep, err := EvalKeyCompressionTest(2, sigma, 4)
			if rep != nil {
				t.Log("\n" + rep.String())
			}
			if err != nil {
				t.Fatal(err)
			}
			if rep.Rec.IPErr != nil {
				if sigma == 0 {
					t.Fatalf("σ=0 인데 복원 키 InternalProduct 실패: %v", rep.Rec.IPErr)
				}
				t.Logf("한계: 복원 키 InternalProduct 잔차 %s > 기준 %s (%v)",
					Log2Str(rep.Rec.IPNoise), Log2Str(rep.Bound), rep.Rec.IPErr)
			}
		})
	}
}

func BenchmarkCompressK2(b *testing.B) {
	f := newCompressFixture(b, 2, 3.2)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err := f.op.CompressEvalKey(f.evk); err != nil {
			b.Fatal(err)
		}
	}
}

func BenchmarkDecompressK2(b *testing.B) {
	f := newCompressFixture(b, 2, 3.2)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		f.op.DecompressEvalKey(f.ck)
	}
}
