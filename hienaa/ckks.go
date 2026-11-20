package hienaa

import (
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/lintrans"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// CKKSParamSketch는 HIENAA의 CKKSParamSketch에 대응합니다.
type CKKSParamSketch struct {
	LogN            int
	LogQ            []int
	LogP            []int
	LogDefaultScale int // HIENAA: scale = 2^40 → 40
	HammingWeight   int
	LogSlots        int // 0이면 N/2 (풀 패킹). 작게 주면 희소 패킹
}

// CKKSScheme은 HIENAA의 CKKSScheme에 대응합니다.
type CKKSScheme struct {
	Params   ckks.Parameters
	Ecd      *ckks.Encoder
	Enc      *rlwe.Encryptor
	Dec      *rlwe.Decryptor
	Eval     *ckks.Evaluator
	LogSlots int

	kgen *rlwe.KeyGenerator
	sk   *rlwe.SecretKey
	rlk  *rlwe.RelinearizationKey
	gks  []*rlwe.GaloisKey
}

func NewCKKSScheme(s CKKSParamSketch) (*CKKSScheme, error) {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN: s.LogN, LogQ: s.LogQ, LogP: s.LogP, LogDefaultScale: s.LogDefaultScale,
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
	logSlots := s.LogSlots
	if logSlots == 0 {
		logSlots = params.LogMaxSlots()
	}
	sc := &CKKSScheme{
		Params: params, Ecd: ckks.NewEncoder(params),
		Enc: rlwe.NewEncryptor(params, sk), Dec: rlwe.NewDecryptor(params, sk),
		LogSlots: logSlots, kgen: kgen, sk: sk,
	}
	sc.refreshEvaluator()
	return sc, nil
}

func (s *CKKSScheme) refreshEvaluator() {
	s.Eval = ckks.NewEvaluator(s.Params, rlwe.NewMemEvaluationKeySet(s.rlk, s.gks...))
}

func (s *CKKSScheme) Slots() int { return 1 << s.LogSlots }

func (s *CKKSScheme) RelinKeygenAndSet() {
	s.rlk = s.kgen.GenRelinearizationKeyNew(s.sk)
	s.refreshEvaluator()
}

// GaloisKeygenAndSet = rotate_keygen(list) + set_rotate_key!(list, ...)
func (s *CKKSScheme) GaloisKeygenAndSet(galEls []uint64) {
	s.gks = append(s.gks, s.kgen.GenGaloisKeysNew(galEls, s.sk)...)
	s.refreshEvaluator()
}

func (s *CKKSScheme) Encode(msg []complex128) (*rlwe.Plaintext, error) {
	pt := ckks.NewPlaintext(s.Params, s.Params.MaxLevel())
	pt.LogDimensions.Cols = s.LogSlots
	return pt, s.Ecd.Encode(msg, pt)
}

func (s *CKKSScheme) Encrypt(pt *rlwe.Plaintext) (*rlwe.Ciphertext, error) {
	return s.Enc.EncryptNew(pt)
}

func (s *CKKSScheme) DecryptDecode(ct *rlwe.Ciphertext) ([]complex128, error) {
	out := make([]complex128, s.Slots())
	return out, s.Ecd.Decode(s.Dec.DecryptNew(ct), out)
}

// PlainMatrix는 HIENAA의 PlainMatrix에 대응합니다.
// 내부적으로 Halevi–Shoup 대각 표현 + BSGS (Lattigo lintrans) 를 사용합니다.
type PlainMatrix struct {
	lt lintrans.LinearTransformation
}

// NewPlainMatrix = PlainMatrix(matrix, scheme). matrix는 n×n (n = Slots()), level은 입력 암호문 레벨.
func (s *CKKSScheme) NewPlainMatrix(matrix [][]complex128, level int) (*PlainMatrix, error) {
	n := s.Slots()
	diags := make(lintrans.Diagonals[complex128])
	for k := 0; k < n; k++ {
		d := make([]complex128, n)
		nonzero := false
		for j := 0; j < n; j++ {
			d[j] = matrix[j][(j+k)%n] // k번째 일반화 대각
			if d[j] != 0 {
				nonzero = true
			}
		}
		if nonzero {
			diags[k] = d
		}
	}
	ltp := lintrans.Parameters{
		DiagonalsIndexList:        diags.DiagonalsIndexList(),
		LevelQ:                    level,
		LevelP:                    s.Params.MaxLevelP(),
		Scale:                     s.Params.GetOptimalScalingFactor(s.Params.DefaultScale(), s.Params.DefaultScale(), level),
		LogDimensions:             ring2D(s.LogSlots),
		LogBabyStepGiantStepRatio: 1,
	}
	lt := lintrans.NewTransformation(s.Params, ltp)
	if err := lintrans.Encode(s.Ecd, diags, lt); err != nil {
		return nil, err
	}
	return &PlainMatrix{lt: lt}, nil
}

// RequiredGaloisElements = get_required_key_list(M)
func (m *PlainMatrix) RequiredGaloisElements(s *CKKSScheme) []uint64 {
	return m.lt.GaloisElements(s.Params)
}

// MatMul = mul(M, ct, scheme)  (곱한 뒤 rescale까지 수행)
func (s *CKKSScheme) MatMul(m *PlainMatrix, ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	res, err := lintrans.NewEvaluator(s.Eval).EvaluateNew(ct, m.lt)
	if err != nil {
		return nil, err
	}
	return res, s.Eval.Rescale(res, res)
}
