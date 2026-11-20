package mkfhe

// 원본 main.jl 의 검증에 쓰인 SimpleScaler 대응 함수들 (big.Int 로 단순 구현).

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v6/ring"
)

// roundDiv 는 round(a/b) (a, b ≥ 0).
func roundDiv(a, b *big.Int) *big.Int {
	t := new(big.Int).Lsh(a, 1)
	t.Add(t, b)
	return t.Quo(t, new(big.Int).Lsh(b, 1))
}

// EncodeConstScaled 는 모든 계수가 round(c·Q/t) 인 다항식을 NTT 형식으로 만듭니다.
// (원본: @. v.coeffs[1] = c; simple_scale_to!(…, SimpleScaler([t], Q)))
func (op *MKOperator) EncodeConstScaled(c, t uint64) ring.Poly {
	r := op.R
	Q := r.ModulusAtLevel[r.Level()]
	x := roundDiv(new(big.Int).Mul(new(big.Int).SetUint64(c), Q), new(big.Int).SetUint64(t))
	p := r.NewPoly()
	for l, s := range r.SubRings {
		v := new(big.Int).Mod(x, new(big.Int).SetUint64(s.Modulus)).Uint64()
		for n := range p.Coeffs[l] {
			p.Coeffs[l][n] = v
		}
	}
	r.NTT(p, p)
	return p
}

// ScaleDown 은 계수 형식 p 의 각 계수에 대해 round(x·t/Q) mod t 를 계산합니다.
func (op *MKOperator) ScaleDown(p ring.Poly, t uint64) []uint64 {
	r := op.R
	Q := r.ModulusAtLevel[r.Level()]
	N := r.N()
	coeffs := make([]*big.Int, N)
	for i := range coeffs {
		coeffs[i] = new(big.Int)
	}
	r.PolyToBigint(p, 1, coeffs)
	T := new(big.Int).SetUint64(t)
	out := make([]uint64, N)
	for i, x := range coeffs {
		x.Mod(x, Q)
		y := roundDiv(new(big.Int).Mul(x, T), Q)
		out[i] = y.Mod(y, T).Uint64()
	}
	return out
}

// GadgetDigitDecode 는 원본 검증 루틴: Q→q_j 로 스케일한 뒤 Π_{k≠j} q_k 를 곱합니다 (mod q_j).
// 위상이 g_j·m + e 이면 m (mod q_j) 이 나옵니다.
func (op *MKOperator) GadgetDigitDecode(p ring.Poly, j int) []uint64 {
	r := op.R
	qj := r.SubRings[j].Modulus
	vals := op.ScaleDown(p, qj)
	prod := uint64(1)
	for k, s := range r.SubRings {
		if k != j {
			prod = ring.BRed(prod, s.Modulus%qj, qj, r.SubRings[j].BRedConstant)
		}
	}
	for i := range vals {
		vals[i] = ring.BRed(vals[i], prod, qj, r.SubRings[j].BRedConstant)
	}
	return vals
}
