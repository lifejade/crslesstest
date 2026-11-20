// 원본 paramfinder.jl 의 Go 포팅: 잡음 추정치를 출력합니다.
package main

import (
	"fmt"
	"math"
	"math/big"

	"github.com/yourname/crsless-go/mkfhe"
)

func bf(x float64) *big.Float { return new(big.Float).SetPrec(256).SetFloat64(x) }

func bpow(b float64, e int) *big.Float {
	r := bf(1)
	for i := 0; i < e; i++ {
		r.Mul(r, bf(b))
	}
	return r
}

// bigLn 은 자연로그 ln(x) (x > 0).
func bigLn(x *big.Float) float64 {
	m := new(big.Float)
	exp := x.MantExp(m) // x = m · 2^exp, m ∈ [0.5, 1)
	mf, _ := m.Float64()
	return math.Log(mf) + float64(exp)*math.Ln2
}

// circuitEstim = circuit_estim(params, depth). 원본과 같이 log 는 자연로그입니다.
func circuitEstim(p mkfhe.MKParameters, depth int) float64 {
	K, N, d := p.MaxParties, float64(p.N), float64(len(p.Q))
	hw := float64(mkfhe.HW)
	const Bsigma = 42.0
	var Bh uint64
	for _, q := range p.Q {
		Bh = max(Bh, q)
	}

	ext := bf(1.0 / 3.0)
	ext.Mul(ext, bf(d*N*Bsigma*float64(Bh)))
	ext.Mul(ext, bpow(hw+1, K-2))
	ext.Mul(ext, new(big.Float).Add(bpow(float64(K), 3), new(big.Float).Mul(bpow(float64(K+1), 3), bf(hw+1))))

	intNoise := N * float64(Bh) * float64(K+1) * d
	return bigLn(ext) + float64(depth)*math.Log(intNoise)
}

func main() {
	params, err := mkfhe.NewMKParameters(mkfhe.RecommendedSketch[16])
	if err != nil {
		panic(err)
	}
	logNoise := circuitEstim(params, 2)

	logQ := make([]float64, len(params.Q))
	sum, maxLog := 0.0, 0.0
	for i, q := range params.Q {
		logQ[i] = math.Log2(float64(q))
		sum += logQ[i]
		maxLog = max(maxLog, logQ[i])
	}
	fmt.Println(logQ, ",", sum)
	fmt.Println(logNoise, ",", logNoise+maxLog)
}
