package mkfhe

// 원본 src/operator.jl 에 대응합니다.

import (
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/ring"
)

// MKOperator 는 링, 가젯 벡터, 작업 버퍼를 묶은 것입니다.
type MKOperator struct {
	Params MKParameters
	R      *ring.Ring
	Glen   int
	Gvec   [][]uint64 // Gvec[i][l] = g_i mod q_l,  g_i = (Q/q_i)·[(Q/q_i)^{-1}]_{q_i}

	ctBuff  []Tensor      // 원본 ct_buff
	decBuff [][]ring.Poly // 원본 tensor_buff (분해 결과)
	inttBuf ring.Poly
}

func NewMKOperator(p MKParameters) (*MKOperator, error) {
	r, err := ring.NewRing(p.N, p.Q)
	if err != nil {
		return nil, err
	}
	K, qlen := p.MaxParties, len(p.Q)
	glen := int(math.Ceil(float64(qlen) / float64(p.Dlen)))

	// 가젯 벡터 (Decomposer(Q, dlen=1))
	Q := r.ModulusAtLevel[r.Level()]
	gvec := make([][]uint64, glen)
	for i := 0; i < glen; i++ {
		qi := new(big.Int).SetUint64(p.Q[i])
		tilde := new(big.Int).Quo(Q, qi)
		inv := new(big.Int).ModInverse(new(big.Int).Mod(tilde, qi), qi)
		g := new(big.Int).Mul(tilde, inv)
		gvec[i] = make([]uint64, qlen)
		for l := 0; l < qlen; l++ {
			gvec[i][l] = new(big.Int).Mod(g, new(big.Int).SetUint64(p.Q[l])).Uint64()
		}
	}

	op := &MKOperator{Params: p, R: r, Glen: glen, Gvec: gvec, inttBuf: r.NewPoly()}
	op.ctBuff = make([]Tensor, 4)
	for i := range op.ctBuff {
		op.ctBuff[i] = NewTensor(r, K+1)
	}
	op.decBuff = make([][]ring.Poly, K+2)
	for i := range op.decBuff {
		op.decBuff[i] = NewTensor(r, glen)
	}
	return op, nil
}

// ---------------- 다항식 보조 연산 ----------------

func (op *MKOperator) addTensor(res, a, b Tensor) {
	for i := range res {
		op.R.Add(a[i], b[i], res[i])
	}
}

// addRNSScalar: p += c (상수 다항식, NTT 형식이므로 모든 슬롯에 더함)
func (op *MKOperator) addRNSScalar(p ring.Poly, c []uint64) {
	for l, s := range op.R.SubRings {
		s.AddScalar(p.Coeffs[l], c[l], p.Coeffs[l])
	}
}

// mulRNSScalarThenAdd: out += p ⊙ c
func (op *MKOperator) mulRNSScalarThenAdd(p ring.Poly, c []uint64, out ring.Poly) {
	for l, s := range op.R.SubRings {
		s.MulScalarMontgomeryThenAdd(p.Coeffs[l], ring.MForm(c[l], s.Modulus, s.BRedConstant), out.Coeffs[l])
	}
}

// DecomposeTo 는 decompose_to! (dlen=1): x(NTT) 의 각 RNS 성분을 중심 표현으로 전체 Q 에 확장한 뒤 NTT 합니다.
func (op *MKOperator) DecomposeTo(dst []ring.Poly, x ring.Poly) {
	r := op.R
	r.INTT(x, op.inttBuf)
	N := r.N()
	for i := 0; i < op.Glen; i++ {
		qi := r.SubRings[i].Modulus
		half := qi >> 1
		src := op.inttBuf.Coeffs[i]
		d := dst[i]
		for l, s := range r.SubRings {
			ql := s.Modulus
			out := d.Coeffs[l]
			if l == i {
				copy(out, src)
				continue
			}
			for n := 0; n < N; n++ {
				v := src[n]
				if v <= half {
					out[n] = v % ql
				} else {
					neg := (qi - v) % ql
					if neg == 0 {
						out[n] = 0
					} else {
						out[n] = ql - neg
					}
				}
			}
		}
		r.NTT(d, d)
	}
}

// ---------------- 원본 operator.jl 의 함수들 ----------------

// HoistedVectorExtProd = hoisted_vector_extprod_to!(res, decv, ct; zeroidx)
func (op *MKOperator) HoistedVectorExtProd(res Tensor, decv [][]ring.Poly, ct *CGSW, zeroidx int) {
	k := ct.Parties()
	if len(res) != k+1 || (len(decv) != k && len(decv) != k+1) || len(ct.Stack) != (k+1)*ct.Glen {
		panic("HoistedVectorExtProd: 차원이 맞지 않습니다")
	}
	acc := op.ctBuff[0][:k+1]
	acc.Zero()
	for i := 0; i < min(len(decv), zeroidx); i++ {
		for j := 0; j < op.Glen; j++ {
			row := ct.Stack[i*op.Glen+j]
			for idx := 0; idx <= k; idx++ {
				op.R.MulCoeffsBarrettThenAdd(decv[i][j], row[idx], acc[idx])
			}
		}
	}
	res.CopyFrom(acc)
}

// VectorExtProd = vector_extprod_to!(res, v, ct; zeroidx).  zeroidx 를 쓰지 않으려면 math.MaxInt.
func (op *MKOperator) VectorExtProd(res, v Tensor, ct *CGSW, zeroidx int) {
	vlen := len(v)
	for i := 0; i < min(vlen, zeroidx); i++ {
		op.DecomposeTo(op.decBuff[i], v[i])
	}
	op.HoistedVectorExtProd(res, op.decBuff[:vlen], ct, zeroidx)
}

// MGSWExtend = mgsw_extend!(res, ct, idx, evk). idx 는 원본과 같이 1부터 시작하는 당사자 번호입니다.
// ct 는 idx 로, evk[i] 는 i+1 로 Shrink 되어 있어야 합니다.
func (op *MKOperator) MGSWExtend(res *MGSW, ct *CGSW, idx int, evk []*EvalKey) {
	n, glen := len(evk), op.Glen

	for i := 1; i <= glen; i++ {
		res.Basket[0].Stack[i-1].Zero()
		op.addRNSScalar(res.Basket[0].Stack[i-1][0], op.Gvec[i-1])
	}

	for i := 1; i <= n; i++ {
		for j := 1; j <= glen; j++ {
			target := res.Basket[i].Stack[j-1]
			target.Zero()
			acci := target[:i+1]
			tmp := op.ctBuff[1][:i+1]

			if i != idx {
				for k := 1; k <= i; k++ {
					acctmp := res.Basket[k-1].Stack[j-1][:i]
					z := k
					if i > idx {
						z = max(k, idx+1)
					}
					op.VectorExtProd(tmp, acctmp, evk[i-1].Vals[k-1], z)
					op.addTensor(acci, acci, tmp)
				}
			} else {
				decbuff := op.decBuff[:i]
				for k := 1; k <= i; k++ {
					acckprev := res.Basket[k-1].Stack[j-1][:i]
					accknext := res.Basket[k-1].Stack[j-1][:i+1]
					for l := 1; l <= k; l++ {
						op.DecomposeTo(decbuff[l-1], acckprev[l-1])
					}
					op.HoistedVectorExtProd(accknext, decbuff, ct, k)
					op.HoistedVectorExtProd(tmp, decbuff, evk[i-1].Vals[k-1], k)
					op.addTensor(acci, acci, tmp)
				}
				op.VectorExtProd(acci, acci, ct, math.MaxInt)
			}
		}
	}
}

// HoistedExternalProduct = hoisted_external_product_to!(res, decxi, y)
func (op *MKOperator) HoistedExternalProduct(res Tensor, decxi [][]ring.Poly, y *MGSW) {
	k := len(y.Basket[0].Stack[0]) - 1
	if len(res) != k+1 || len(decxi) != k+1 {
		panic("HoistedExternalProduct: 차원이 맞지 않습니다")
	}
	acc := op.ctBuff[0][:k+1]
	acc.Zero()
	for i := 0; i <= k; i++ {
		for j := 0; j < op.Glen; j++ {
			row := y.Basket[i].Stack[j]
			for idx := 0; idx <= k; idx++ {
				op.R.MulCoeffsBarrettThenAdd(decxi[i][j], row[idx], acc[idx])
			}
		}
	}
	res.CopyFrom(acc)
}

// InternalProduct = internal_product_to!(res, x, y): h(x)^T · y
func (op *MKOperator) InternalProduct(res, x, y *MGSW) {
	k := len(x.Basket[0].Stack[0]) - 1
	if len(res.Basket) != k+1 || len(x.Basket) != k+1 || len(y.Basket) != k+1 {
		panic("InternalProduct: 차원이 맞지 않습니다")
	}
	dec := op.decBuff[:k+1]
	for i := 0; i <= k; i++ {
		for j := 0; j < op.Glen; j++ {
			for idx := 0; idx <= k; idx++ {
				op.DecomposeTo(dec[idx], x.Basket[i].Stack[j][idx])
			}
			op.HoistedExternalProduct(res.Basket[i].Stack[j], dec, y)
		}
	}
}
