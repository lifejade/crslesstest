package mkfhe

// 평가키 압축 / 서버 측 복원 (신규, 원본 Julia 에 대응 함수 없음).
//
// 표기 (명세와 동일, 차원 인덱스는 0부터):
//   t = (t_0..t_{K-1}) = CGSWKey.Val,  S = [I_K | t]
//   열 벡터 (B; A) = Tensor (b_0..b_{K-1}, a),  위상 S·x = b_r + a·t_r
//   G_h = EvalKey.Vals[h] = cGSW_t(t_h),  블록 G_h^{[i]} = Stack[i*Glen : (i+1)*Glen]
//
// 압축키 ek_comp = (U_0, H_⋆, H_Δ, V_1..V_{K-1}),  환 원소 (5K+1)·d 개.

import (
	"bytes"
	"fmt"

	"github.com/tuneinsight/lattigo/v6/ring"
)

// GadgetBlock 은 가젯 블록 하나 ((rows)×d 행렬) 를 d 개의 열 Tensor 로 보관합니다.
// cGSW 블록이면 각 열이 K+1 개, 1차원 cLEV 이면 2 개 (b, a) 다항식입니다.
type GadgetBlock []Tensor

func NewGadgetBlock(r *ring.Ring, rows, glen int) GadgetBlock {
	b := make(GadgetBlock, glen)
	for j := range b {
		b[j] = NewTensor(r, rows)
	}
	return b
}

// Block 은 G^{[i]} (i = 0..K) 의 뷰를 돌려줍니다 (데이터 공유).
func (c *CGSW) Block(i int) GadgetBlock {
	return GadgetBlock(c.Stack[i*c.Glen : (i+1)*c.Glen])
}

func (b GadgetBlock) copyNew() GadgetBlock {
	out := make(GadgetBlock, len(b))
	for j, col := range b {
		out[j] = make(Tensor, len(col))
		for r := range col {
			out[j][r] = *col[r].CopyNew()
		}
	}
	return out
}

// CompressedEvalKey = ek_comp.
type CompressedEvalKey struct {
	K, Glen int
	U0      GadgetBlock   // G_0^{[0]}
	HStar   GadgetBlock   // G_0^{[K]}
	HDelta  GadgetBlock   // Σ_i (G_i^{[i]} − G_{i−1}^{[i]})
	V       []GadgetBlock // V[h] = cLEV_{t_0}(t_h) (2 행), h = 1..K−1. V[0] 은 nil (U0 에서 추출)
}

// CLEV 는 V_h 를 돌려줍니다. h = 0 이면 U_0 의 첫 행과 마지막 행을 공유하는 뷰입니다.
func (ck *CompressedEvalKey) CLEV(h int) GadgetBlock {
	if h != 0 {
		return ck.V[h]
	}
	v := make(GadgetBlock, ck.Glen)
	for j := range v {
		v[j] = Tensor{ck.U0[j][0], ck.U0[j][ck.K]}
	}
	return v
}

// NumPolys 는 압축키의 환 원소 개수 ((5K+1)·d).
func (ck *CompressedEvalKey) NumPolys() int {
	n := 0
	for _, b := range append([]GadgetBlock{ck.U0, ck.HStar, ck.HDelta}, ck.V[1:]...) {
		for _, col := range b {
			n += len(col)
		}
	}
	return n
}

// MarshalBinary 는 압축키의 모든 다항식을 Lattigo Poly 직렬화로 이어 붙입니다 (실측 통신량용).
func (ck *CompressedEvalKey) MarshalBinary() ([]byte, error) {
	var buf bytes.Buffer
	for _, b := range append([]GadgetBlock{ck.U0, ck.HStar, ck.HDelta}, ck.V[1:]...) {
		for _, col := range b {
			for _, p := range col {
				if _, err := p.WriteTo(&buf); err != nil {
					return nil, err
				}
			}
		}
	}
	return buf.Bytes(), nil
}

// EvalKeyMarshalBinary 는 기존 EvalKey 전체를 같은 방식으로 직렬화합니다 (비교용).
func EvalKeyMarshalBinary(evk *EvalKey) ([]byte, error) {
	var buf bytes.Buffer
	for _, c := range evk.Vals {
		for _, col := range c.Stack {
			for _, p := range col {
				if _, err := p.WriteTo(&buf); err != nil {
					return nil, err
				}
			}
		}
	}
	return buf.Bytes(), nil
}

// CompressEvalKey 는 기존 KeyGen 결과 (Shrink 되지 않은 K 당사자 EvalKey) 에서 ek_comp 를 추출합니다.
func (op *MKOperator) CompressEvalKey(evk *EvalKey) (*CompressedEvalKey, error) {
	K := len(evk.Vals)
	for h, c := range evk.Vals {
		if c.Parties() != K || c.Glen != op.Glen || len(c.Stack) != (K+1)*op.Glen {
			return nil, fmt.Errorf("CompressEvalKey: G_%d 차원 불일치 (Shrink 전 키가 필요)", h)
		}
	}
	r := op.R
	ck := &CompressedEvalKey{
		K: K, Glen: op.Glen,
		U0:     evk.Vals[0].Block(0).copyNew(),
		HStar:  evk.Vals[0].Block(K).copyNew(),
		HDelta: NewGadgetBlock(r, K+1, op.Glen),
		V:      make([]GadgetBlock, K),
	}
	// H_Δ = Σ_i (G_i^{[i]} − G_{(i−1) mod K}^{[i]})
	for i := 0; i < K; i++ {
		plus, minus := evk.Vals[i].Block(i), evk.Vals[(i-1+K)%K].Block(i)
		for j := 0; j < op.Glen; j++ {
			for row := 0; row <= K; row++ {
				r.Add(ck.HDelta[j][row], plus[j][row], ck.HDelta[j][row])
				r.Sub(ck.HDelta[j][row], minus[j][row], ck.HDelta[j][row])
			}
		}
	}
	// V_h = G_h^{[0]} 의 첫 행 (b_0) 과 마지막 행 (a)
	for h := 1; h < K; h++ {
		blk := evk.Vals[h].Block(0)
		ck.V[h] = make(GadgetBlock, op.Glen)
		for j := range ck.V[h] {
			ck.V[h][j] = Tensor{*blk[j][0].CopyNew(), *blk[j][K].CopyNew()}
		}
	}
	return ck, nil
}

// blockMulDecThenAdd: acc += X · D_g(·)  (dec 는 이미 분해된 d 개 다항식).
func (op *MKOperator) blockMulDecThenAdd(acc Tensor, X GadgetBlock, dec []ring.Poly) {
	for j := range X {
		for row := range acc {
			op.R.MulCoeffsBarrettThenAdd(dec[j], X[j][row], acc[row])
		}
	}
}

// UpKeyStep 은 R(X) = (Π B_X; A_X) − H_Δ D_g(A_X) 를 dst 에 씁니다. (Π b)_i = b_{(i−1) mod K}.
// dst 와 X 는 겹치면 안 됩니다.
func (op *MKOperator) UpKeyStep(dst, X, hDelta GadgetBlock) {
	K := len(X[0]) - 1
	dec := NewTensor(op.R, op.Glen)
	tmp := NewTensor(op.R, K+1)
	for j := range X {
		op.DecomposeTo(dec, X[j][K])
		tmp.Zero()
		op.blockMulDecThenAdd(tmp, hDelta, dec)
		for i := 0; i < K; i++ {
			op.R.Sub(X[j][(i-1+K)%K], tmp[i], dst[j][i])
		}
		op.R.Sub(X[j][K], tmp[K], dst[j][K])
	}
}

// UpKeys 는 U_i = R^i(U_0), i = 0..depth−1 을 만듭니다 (보통 depth = K). U_0 은 복사본.
func (op *MKOperator) UpKeys(ck *CompressedEvalKey, depth int) []GadgetBlock {
	ups := make([]GadgetBlock, depth)
	ups[0] = ck.U0.copyNew()
	for i := 1; i < depth; i++ {
		ups[i] = NewGadgetBlock(op.R, ck.K+1, op.Glen)
		op.UpKeyStep(ups[i], ups[i-1], ck.HDelta)
	}
	return ups
}

// ExtendCLEV 는 1차원 cLEV_{t_0}(μ) = V 와 Up key U_0..U_{K−1}, H_⋆ 로 완전한 K차원 cGSW_t(μ) 를 dst 에 만듭니다.
//
//	C_i = (e_i b; 0) + U_i D_g(a),  0 ≤ i < K
//	C_K = (0_K; b)   + H_⋆ D_g(a)
func (op *MKOperator) ExtendCLEV(dst *CGSW, V GadgetBlock, ups []GadgetBlock, hStar GadgetBlock) {
	K := dst.Parties()
	if len(ups) < K || len(dst.Stack) != (K+1)*op.Glen || len(V) != op.Glen || len(V[0]) != 2 {
		panic("ExtendCLEV: 차원이 맞지 않습니다")
	}
	dec := NewTensor(op.R, op.Glen)
	for j := 0; j < op.Glen; j++ {
		b, a := V[j][0], V[j][1]
		op.DecomposeTo(dec, a) // 모든 블록에서 재사용 (hoisting)
		for i := 0; i < K; i++ {
			col := dst.Stack[i*op.Glen+j]
			col.Zero()
			op.blockMulDecThenAdd(col, ups[i], dec)
			op.R.Add(col[i], b, col[i])
		}
		col := dst.Stack[K*op.Glen+j]
		col.Zero()
		op.blockMulDecThenAdd(col, hStar, dec)
		op.R.Add(col[K], b, col[K])
	}
}

// DecompressEvalKey 는 ek_comp 만으로 Ĝ_0..Ĝ_{K−1} 전체를 복원합니다 (Shrink 전 형태).
func (op *MKOperator) DecompressEvalKey(ck *CompressedEvalKey) *EvalKey {
	ups := op.UpKeys(ck, ck.K)
	evk := NewEvalKey(op.R, ck.K, op.Glen)
	for h := 0; h < ck.K; h++ {
		op.ExtendCLEV(evk.Vals[h], ck.CLEV(h), ups, ck.HStar)
	}
	return evk
}
