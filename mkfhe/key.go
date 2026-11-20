package mkfhe

// 원본 src/key.jl 에 대응합니다.

import (
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

// CGSWKey 는 k 개의 3진 비밀키 (NTT 형식으로 보관).
type CGSWKey struct {
	Val []ring.Poly
}

// TernaryCGSWKey = ternary_cgswkey(us, N, k, hw). hw=0 이면 {-1,0,1} 균등 분포.
func TernaryCGSWKey(r *ring.Ring, k, hw int) (*CGSWKey, error) {
	prng, err := sampling.NewPRNG()
	if err != nil {
		return nil, err
	}
	dist := ring.Ternary{H: hw}
	if hw == 0 {
		dist = ring.Ternary{P: 2.0 / 3.0}
	}
	ts, err := ring.NewSampler(prng, r, dist, false)
	if err != nil {
		return nil, err
	}
	key := &CGSWKey{Val: make([]ring.Poly, k)}
	for i := range key.Val {
		key.Val[i] = ts.ReadNew()
		r.NTT(key.Val[i], key.Val[i])
	}
	return key, nil
}
