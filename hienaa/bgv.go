// Package hienaa는 HIENAA.jl의 고수준 API(Scheme / encode / encrypt / evaluate ...)를
// Lattigo v6 위에 얇게 재현한 래퍼입니다.
//
// 주의: Lattigo는 2의 거듭제곱 원분환 Z[X]/(X^N+1)만 지원하므로,
// HIENAA의 임의 원분환(m=4369, 21845 등)·subring 구조는 이 래퍼로 옮겨지지 않습니다.
package hienaa

import (
	"github.com/tuneinsight/lattigo/v6/circuits/bgv/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/bgv"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

// BGVParamSketch는 HIENAA의 BGVParamSketch / BFVParamSketch에 대응합니다.
// HIENAA의 (m, logP, logQ, t) 대신 Lattigo는 (LogN, 소수별 비트 길이 목록, t)를 받습니다.
type BGVParamSketch struct {
	LogN             int    // N = 2^LogN  (HIENAA: m = 2N)
	LogQ             []int  // 각 Q 소수의 비트 길이 (HIENAA: logQ 총합)
	LogP             []int  // 키 스위칭용 P 소수 (HIENAA: logP)
	PlaintextModulus uint64 // t. 슬롯 패킹을 쓰려면 t ≡ 1 (mod 2N)인 소수여야 함
	HammingWeight    int    // 희소 3진 비밀키의 해밍 웨이트 (HIENAA: hw). 0이면 기본 분포
	ScaleInvariant   bool   // true면 BFV 동작 (HIENAA의 BFVScheme)
}

// BGVScheme은 HIENAA의 BGVScheme/BFVScheme에 대응합니다.
type BGVScheme struct {
	Params bgv.Parameters
	Ecd    *bgv.Encoder
	Enc    *rlwe.Encryptor
	Dec    *rlwe.Decryptor
	Eval   *bgv.Evaluator

	kgen           *rlwe.KeyGenerator
	sk             *rlwe.SecretKey
	rlk            *rlwe.RelinearizationKey
	gks            []*rlwe.GaloisKey
	scaleInvariant bool
}

// NewBGVScheme = BGVScheme(sketch) + ternary_ringkey + set_encryptor!
func NewBGVScheme(s BGVParamSketch) (*BGVScheme, error) {
	params, err := bgv.NewParametersFromLiteral(bgv.ParametersLiteral{
		LogN: s.LogN, LogQ: s.LogQ, LogP: s.LogP, PlaintextModulus: s.PlaintextModulus,
	})
	if err != nil {
		return nil, err
	}
	kgen := rlwe.NewKeyGenerator(params)
	var sk *rlwe.SecretKey
	if s.HammingWeight > 0 {
		sk = kgen.GenSecretKeyWithHammingWeightNew(s.HammingWeight)
	} else {
		sk = kgen.GenSecretKeyNew()
	}
	sc := &BGVScheme{
		Params: params, Ecd: bgv.NewEncoder(params),
		Enc: rlwe.NewEncryptor(params, sk), Dec: rlwe.NewDecryptor(params, sk),
		kgen: kgen, sk: sk, scaleInvariant: s.ScaleInvariant,
	}
	sc.refreshEvaluator()
	return sc, nil
}

func (s *BGVScheme) refreshEvaluator() {
	evk := rlwe.NewMemEvaluationKeySet(s.rlk, s.gks...)
	s.Eval = bgv.NewEvaluator(s.Params, evk, s.scaleInvariant)
}

// Slots = scheme.oper.packer.k
func (s *BGVScheme) Slots() int { return s.Params.MaxSlots() }

// RelinKeygenAndSet = relin_keygen + set_relinkey!
func (s *BGVScheme) RelinKeygenAndSet() {
	s.rlk = s.kgen.GenRelinearizationKeyNew(s.sk)
	s.refreshEvaluator()
}

// RotateKeygenAndSet = rotate_keygen + set_rotate_key! (열 회전 k 목록)
func (s *BGVScheme) RotateKeygenAndSet(rots []int) {
	gal := make([]uint64, len(rots))
	for i, k := range rots {
		gal[i] = s.Params.GaloisElementForColRotation(k)
	}
	s.gks = append(s.gks, s.kgen.GenGaloisKeysNew(gal, s.sk)...)
	s.refreshEvaluator()
}

// Encode = encode(msg, scheme)
func (s *BGVScheme) Encode(msg []uint64) (*rlwe.Plaintext, error) {
	pt := bgv.NewPlaintext(s.Params, s.Params.MaxLevel())
	return pt, s.Ecd.Encode(msg, pt)
}

// Encrypt = encrypt(pt, scheme)
func (s *BGVScheme) Encrypt(pt *rlwe.Plaintext) (*rlwe.Ciphertext, error) {
	return s.Enc.EncryptNew(pt)
}

// DecryptDecode = decrypt_to! + decode
func (s *BGVScheme) DecryptDecode(ct *rlwe.Ciphertext) ([]uint64, error) {
	out := make([]uint64, s.Slots())
	return out, s.Ecd.Decode(s.Dec.DecryptNew(ct), out)
}

// EvaluatePoly = evaluate(PolyCoeffs(coeffs, scheme), ct, scheme)
// coeffs[i]는 x^i의 계수 (mod t). 재선형화 키가 필요합니다.
func (s *BGVScheme) EvaluatePoly(coeffs []uint64, ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	poly := bignum.NewPolynomial(bignum.Monomial, coeffs, nil)
	return polynomial.NewEvaluator(s.Params, s.Eval).Evaluate(ct, poly, s.Params.DefaultScale())
}
