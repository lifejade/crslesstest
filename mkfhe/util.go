package mkfhe

import "math/bits"

func mulmod(a, b, q uint64) uint64 {
	hi, lo := bits.Mul64(a, b)
	_, r := bits.Div64(hi%q, lo, q)
	return r
}
