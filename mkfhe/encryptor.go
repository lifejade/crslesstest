package mkfhe

// 원본 src/encryptor.jl 에 대응합니다.

import (
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

type MKEncryptor struct {
	Key *CGSWKey
	Op  *MKOperator
	us  ring.Sampler
	gs  ring.Sampler
	buf ring.Poly
}

// NewMKEncryptor = MKEncryptor(key, σ, oper)
func NewMKEncryptor(key *CGSWKey, sigma float64, op *MKOperator) (*MKEncryptor, error) {
	prng, err := sampling.NewPRNG()
	if err != nil {
		return nil, err
	}
	us, err := ring.NewSampler(prng, op.R, ring.Uniform{}, false)
	if err != nil {
		return nil, err
	}
	gs, err := ring.NewSampler(prng, op.R, ring.DiscreteGaussian{Sigma: sigma, Bound: 6 * sigma}, false)
	if err != nil {
		return nil, err
	}
	return &MKEncryptor{Key: key, Op: op, us: us, gs: gs, buf: op.R.NewPoly()}, nil
}

// SalweSampleTo = salwe_sample_to!: res = (e_1 - a·s_1, ..., e_K - a·s_K, a)
func (e *MKEncryptor) SalweSampleTo(res Tensor) {
	r := e.Op.R
	K := len(res) - 1
	e.us.Read(res[K]) // a: NTT 영역에서 바로 균등 샘플
	for i := 0; i < K; i++ {
		e.gs.Read(res[i])
		r.NTT(res[i], res[i])
		r.MulCoeffsBarrett(res[K], e.Key.Val[i], e.buf)
		r.Sub(res[i], e.buf, res[i])
	}
}

func (e *MKEncryptor) SalweSample(K int) Tensor {
	res := NewTensor(e.Op.R, K+1)
	e.SalweSampleTo(res)
	return res
}

// Phase = phase(ct, entor): res_i = ct_i + ct_{K+1}·s_i (계수 형식으로 반환)
func (e *MKEncryptor) Phase(ct Tensor) Tensor {
	r := e.Op.R
	K := len(ct) - 1
	res := NewTensor(r, K)
	for i := 0; i < K; i++ {
		res[i].Copy(ct[i])
		r.MulCoeffsBarrettThenAdd(ct[K], e.Key.Val[i], res[i])
		r.INTT(res[i], res[i])
	}
	return res
}

// CGSWEncryptConst = cgsw_encrypt(PlainConst(m), k, entor)
func (e *MKEncryptor) CGSWEncryptConst(m uint64, k int) *CGSW {
	res := NewCGSW(e.Op.R, k, e.Op.Glen)
	e.CGSWEncryptConstTo(res, m)
	return res
}

func (e *MKEncryptor) CGSWEncryptConstTo(res *CGSW, m uint64) {
	op := e.Op
	glen := op.Glen
	K := len(res.Stack)/glen - 1
	c := make([]uint64, len(op.R.SubRings))
	for i := 1; i <= K+1; i++ {
		for j := 1; j <= glen; j++ {
			row := res.Stack[(i-1)*glen+j-1]
			e.SalweSampleTo(row)
			for l, s := range op.R.SubRings {
				c[l] = ring.BRed(op.Gvec[j-1][l], m%s.Modulus, s.Modulus, s.BRedConstant)
			}
			op.addRNSScalar(row[i-1], c)
		}
	}
}

// CGSWEncryptPolyTo = cgsw_encrypt_to!(res, PlainPoly(m), entor). m 은 NTT 형식.
func (e *MKEncryptor) CGSWEncryptPolyTo(res *CGSW, m ring.Poly) {
	op := e.Op
	glen := op.Glen
	K := len(res.Stack)/glen - 1
	for i := 1; i <= K+1; i++ {
		for j := 1; j <= glen; j++ {
			row := res.Stack[(i-1)*glen+j-1]
			e.SalweSampleTo(row)
			op.mulRNSScalarThenAdd(m, op.Gvec[j-1], row[i-1])
		}
	}
}

// EvkEncrypt = evk_encrypt(K, entor): 자기 비밀키 각 성분 s_i 의 cGSW 암호문.
func (e *MKEncryptor) EvkEncrypt(K int) *EvalKey {
	res := NewEvalKey(e.Op.R, K, e.Op.Glen)
	for i := 0; i < K; i++ {
		e.CGSWEncryptPolyTo(res.Vals[i], e.Key.Val[i])
	}
	return res
}

// JointDecryptor = JointDecryptor(keys, oper)
type JointDecryptor struct {
	Keys     []*CGSWKey
	JointKey []ring.Poly // NTT 형식
	Op       *MKOperator
}

func NewJointDecryptor(keys []*CGSWKey, op *MKOperator) *JointDecryptor {
	r := op.R
	n := len(keys)
	jk := make([]ring.Poly, n)
	for i := 1; i <= n; i++ {
		// jointkey[i] = s_{i,1} + Σ_{j=2}^{i} jointkey[j-1] · s_{i,j}
		jk[i-1] = *keys[i-1].Val[0].CopyNew()
		for j := 2; j <= i; j++ {
			r.MulCoeffsBarrettThenAdd(jk[j-2], keys[i-1].Val[j-1], jk[i-1])
		}
	}
	return &JointDecryptor{Keys: keys, JointKey: jk, Op: op}
}

// JointDecrypt = joint_decrypt(ct, jd): ct_1 + Σ_i ct_{i+1}·jointkey_i  (NTT 형식으로 반환)
func (jd *JointDecryptor) JointDecrypt(ct Tensor) ring.Poly {
	r := jd.Op.R
	res := *ct[0].CopyNew()
	for i := 1; i <= len(jd.Keys); i++ {
		r.MulCoeffsBarrettThenAdd(ct[i], jd.JointKey[i-1], res)
	}
	return res
}
