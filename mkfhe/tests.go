package mkfhe

// 원본 main.jl 의 테스트 4개. cmd/crsless 와 go test 에서 함께 사용합니다.

import (
	"fmt"
	"math"
	"math/rand"
	"time"
)

const ptxtModulus = 256

func setup(testID int) (*MKOperator, int, error) {
	sketch, ok := RecommendedSketch[testID]
	if !ok {
		return nil, 0, fmt.Errorf("test_id %d 없음", testID)
	}
	params, err := NewMKParameters(sketch)
	if err != nil {
		return nil, 0, err
	}
	op, err := NewMKOperator(params)
	return op, params.MaxParties, err
}

func newParties(op *MKOperator, K int) ([]*CGSWKey, []*MKEncryptor, error) {
	keys := make([]*CGSWKey, K)
	entors := make([]*MKEncryptor, K)
	for i := range keys {
		var err error
		if keys[i], err = TernaryCGSWKey(op.R, K, HW); err != nil {
			return nil, nil, err
		}
		if entors[i], err = NewMKEncryptor(keys[i], 3.2, op); err != nil {
			return nil, nil, err
		}
	}
	return keys, entors, nil
}

// VectorGadgetProdTest = vector_gadgetprod_test
func VectorGadgetProdTest(testID int) error {
	op, K, err := setup(testID)
	if err != nil {
		return err
	}
	// 원본은 K 성분짜리 키 하나로 암호화합니다.
	key, err := TernaryCGSWKey(op.R, K, HW)
	if err != nil {
		return err
	}
	entor, err := NewMKEncryptor(key, 3.2, op)
	if err != nil {
		return err
	}

	msg := uint64(rand.Intn(ptxtModulus))
	ct := entor.CGSWEncryptConst(msg, K)

	v := NewTensor(op.R, K)
	for i := 1; i <= K; i++ {
		v[i-1] = op.EncodeConstScaled(uint64(i), ptxtModulus)
	}
	res := NewTensor(op.R, K+1)

	t0 := time.Now()
	op.VectorExtProd(res, v, ct, math.MaxInt)
	fmt.Printf("Vector Gadget Product Took :           %v\n", time.Since(t0))

	p := entor.Phase(res)
	for i := 1; i <= K; i++ {
		got := op.ScaleDown(p[i-1], ptxtModulus)
		want := (msg * uint64(i)) % ptxtModulus
		for _, g := range got {
			if g != want {
				return fmt.Errorf("FAIL: VECTOR GADGET PRODUCT // TESTID %d", testID)
			}
		}
	}
	return nil
}

// MKVectorGadgetProdTest = mk_vector_gadgetprod_test
func MKVectorGadgetProdTest(testID int) error {
	op, K, err := setup(testID)
	if err != nil {
		return err
	}
	keys, entors, err := newParties(op, K)
	if err != nil {
		return err
	}

	msg := make([]uint64, K)
	cts := make([]*CGSW, K)
	for i := range msg {
		msg[i] = uint64(rand.Intn(ptxtModulus))
		cts[i] = entors[i].CGSWEncryptConst(msg[i], K)
	}

	v := make([]Tensor, K+1)
	for i := range v {
		v[i] = NewTensor(op.R, i+1)
	}
	v[0][0] = op.EncodeConstScaled(1, ptxtModulus)

	t0 := time.Now()
	for i := 1; i <= K; i++ {
		cts[i-1].Shrink(i)
		op.VectorExtProd(v[i], v[i-1], cts[i-1], math.MaxInt)
	}
	fmt.Printf("Multi-Key Vector Gadget Product Took : %v\n", time.Since(t0))

	jd := NewJointDecryptor(keys, op)
	res := jd.JointDecrypt(v[K])
	op.R.INTT(res, res)
	want := uint64(1)
	for _, m := range msg {
		want = want * m % ptxtModulus
	}
	for _, g := range op.ScaleDown(res, ptxtModulus) {
		if g != want {
			return fmt.Errorf("FAIL: MULTI-KEY VECTOR GADGET PRODUCT // TESTID %d", testID)
		}
	}
	return nil
}

// checkMGSW 는 원본의 RGSW 검증 루프: basket[1] 은 g_j·m, basket[i+1] 은 g_j·m·jointkey_i 로 복호화되어야 합니다.
func checkMGSW(op *MKOperator, jd *JointDecryptor, m *MGSW, want uint64, K int, name string) error {
	r := op.R
	for j := 0; j < op.Glen; j++ {
		p := jd.JointDecrypt(m.Basket[0].Stack[j])
		r.INTT(p, p)
		vals := op.GadgetDigitDecode(p, j)
		if vals[0] != want%ptxtModulus {
			return fmt.Errorf("FAIL: %s (상수항)", name)
		}
		for _, x := range vals[1:] {
			if x != 0 {
				return fmt.Errorf("FAIL: %s (비상수항)", name)
			}
		}
		msg := vals[0]
		qj := r.SubRings[j].Modulus

		for i := 1; i <= K; i++ {
			p := jd.JointDecrypt(m.Basket[i].Stack[j])
			r.INTT(p, p)
			vals := op.GadgetDigitDecode(p, j)
			jk := *jd.JointKey[i-1].CopyNew()
			r.INTT(jk, jk)
			for n, x := range vals {
				e := mulmod(msg, jk.Coeffs[j][n], qj)
				if x != e {
					return fmt.Errorf("FAIL: %s (basket %d, digit %d, coeff %d)", name, i+1, j+1, n)
				}
			}
		}
	}
	return nil
}

// RGSWExtendTest = rgsw_extend_test
func RGSWExtendTest(testID int) error {
	op, K, err := setup(testID)
	if err != nil {
		return err
	}
	keys, entors, err := newParties(op, K)
	if err != nil {
		return err
	}
	msgs := make([]uint64, K)
	cts := make([]*CGSW, K)
	evks := make([]*EvalKey, K)
	for i := range msgs {
		msgs[i] = uint64(rand.Intn(ptxtModulus))
		cts[i] = entors[i].CGSWEncryptConst(msgs[i], K)
		evks[i] = entors[i].EvkEncrypt(K)
	}

	mgsw := make([]*MGSW, K)
	for i := range mgsw {
		mgsw[i] = NewMGSW(op.R, K, op.Glen)
	}
	for i := 1; i <= K; i++ {
		cts[i-1].Shrink(i)
		evks[i-1].Shrink(i)
	}
	for i := 1; i <= K-1; i++ {
		op.MGSWExtend(mgsw[i-1], cts[i-1], i, evks)
	}
	t0 := time.Now()
	op.MGSWExtend(mgsw[K-1], cts[K-1], K, evks)
	fmt.Printf("RGSW Extension Took :                  %v\n", time.Since(t0))

	jd := NewJointDecryptor(keys, op)
	for idx := 1; idx <= K; idx++ {
		if err := checkMGSW(op, jd, mgsw[idx-1], msgs[idx-1], K, fmt.Sprintf("RGSW EXTENSION // TESTID %d", testID)); err != nil {
			return err
		}
	}
	return nil
}

// RGSWInternalProductTest = rgsw_internal_product_test
//
// 원본은 `while idx1 != idx2` 로 되어 있어 항상 idx1 == idx2 (같은 암호문의 제곱)가 됩니다.
// distinct=true 로 주면 서로 다른 두 당사자의 MGSW 를 곱합니다.
func RGSWInternalProductTest(testID int, distinct bool) error {
	op, K, err := setup(testID)
	if err != nil {
		return err
	}
	keys, entors, err := newParties(op, K)
	if err != nil {
		return err
	}
	msgs := make([]uint64, K)
	cts := make([]*CGSW, K)
	evks := make([]*EvalKey, K)
	for i := range msgs {
		msgs[i] = uint64(rand.Intn(2))
		cts[i] = entors[i].CGSWEncryptConst(msgs[i], K)
		evks[i] = entors[i].EvkEncrypt(K)
	}

	mgsw := make([]*MGSW, K)
	for i := range mgsw {
		mgsw[i] = NewMGSW(op.R, K, op.Glen)
	}
	for i := 1; i <= K; i++ {
		cts[i-1].Shrink(i)
		evks[i-1].Shrink(i)
	}
	res := NewMGSW(op.R, K, op.Glen)
	for i := 1; i <= K; i++ {
		op.MGSWExtend(mgsw[i-1], cts[i-1], i, evks)
	}

	idx1, idx2 := rand.Intn(K)+1, rand.Intn(K)+1
	if distinct {
		for idx1 == idx2 {
			idx1, idx2 = rand.Intn(K)+1, rand.Intn(K)+1
		}
	} else {
		for idx1 != idx2 { // 원본 그대로
			idx1, idx2 = rand.Intn(K)+1, rand.Intn(K)+1
		}
	}

	t0 := time.Now()
	op.InternalProduct(res, mgsw[idx1-1], mgsw[idx2-1])
	fmt.Printf("RGSW Internal Product Took :           %v\n", time.Since(t0))

	jd := NewJointDecryptor(keys, op)
	return checkMGSW(op, jd, res, msgs[idx1-1]*msgs[idx2-1], K, fmt.Sprintf("RGSW INTERNAL PRODUCT // TESTID %d", testID))
}
