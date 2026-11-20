// Package mkfhe는 SNUCP/crsless (Multi-key FHE with Non-Interactive Setup in the Plain Model,
// https://eprint.iacr.org/2026/322) 의 Julia 구현을 Lattigo v6 의 ring 패키지 위로 옮긴 것입니다.
//
// 원본 src/parameters.jl 에 대응합니다.
package mkfhe

import (
	"fmt"
	"math"
	"math/big"
)

// MKParamSketch = MKParamSketch(max_parties, CyclotomicParam(M), logQ, Qlen)
type MKParamSketch struct {
	MaxParties int
	M          int     // 원분 차수 (2의 거듭제곱). 링 차수 N = M/2
	LogQ       float64 // Q 전체 비트 수
	Qlen       int     // Q 소수 개수. 0이면 ceil(LogQ/62)
}

// MKParameters = MKParameters(sketch). dlen 은 원본과 같이 항상 1 (RNS 소수 하나당 가젯 원소 하나).
type MKParameters struct {
	MaxParties int
	M, N       int
	Q          []uint64
	Dlen       int
}

func NewMKParameters(s MKParamSketch) (MKParameters, error) {
	if s.M < 4 || s.M&(s.M-1) != 0 {
		return MKParameters{}, fmt.Errorf("M=%d: Lattigo는 2의 거듭제곱 원분환만 지원합니다", s.M)
	}
	qlen := s.Qlen
	if qlen == 0 {
		qlen = int(math.Ceil(s.LogQ / 62))
	}
	q, err := FindPrimePow2(s.M, s.LogQ/float64(qlen), qlen)
	if err != nil {
		return MKParameters{}, err
	}
	return MKParameters{MaxParties: s.MaxParties, M: s.M, N: s.M / 2, Q: q, Dlen: 1}, nil
}

// FindPrimePow2 는 HIENAA find_prime_cyclotomic (m이 2의 거듭제곱인 경우)을 그대로 옮긴 것입니다.
// 2^b 근처에서 시작해 p ≡ 1 (mod m) 인 소수를 위로 찾고, 62비트를 넘으면 아래로 찾습니다.
func FindPrimePow2(m int, b float64, n int) ([]uint64, error) {
	if b > 62 {
		return nil, fmt.Errorf("62비트보다 큰 소수는 지원하지 않습니다")
	}
	step := uint64(m)
	initial := uint64(math.Floor(math.Pow(2, b)/float64(step)))*step + 1
	current := nextPrime(initial, step)
	if initial == current {
		initial -= step
	}

	primes := make([]uint64, 0, n)
	for math.Log2(float64(current)) < 62 && len(primes) < n {
		primes = append(primes, current)
		current = nextPrime(current+step, step)
	}
	if len(primes) < n {
		current = prevPrime(initial, step)
		for len(primes) < n {
			if current <= step {
				return nil, fmt.Errorf("조건을 만족하는 소수가 부족합니다")
			}
			primes = append(primes, current)
			current = prevPrime(current-step, step)
		}
	}
	return primes, nil
}

func isPrime(x uint64) bool { return new(big.Int).SetUint64(x).ProbablyPrime(30) }

func nextPrime(x, step uint64) uint64 {
	for !isPrime(x) {
		x += step
	}
	return x
}

func prevPrime(x, step uint64) uint64 {
	for x > step && !isPrime(x) {
		x -= step
	}
	return x
}
