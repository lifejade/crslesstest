package mkfhe

// 평가키 압축/복원의 정확성·노이즈·성능 측정 (신규). tests.go 와 같이 go test 와 cmd/evkcompress 에서 함께 씁니다.
//
// 노이즈 잔차는 각 열의 위상에서 기대 메시지를 뺀 값의 중심 표현 무한대 노름입니다.
//   잔차(X; msgs) = max_{j, r} ‖ (X[j][r] + X[j][last]·t_r) − g_j·msgs[r] ‖∞
// 가젯 자리 j 의 복호가 정확하려면 |e| < Q / (2 q_j) 이어야 하므로, 판정 기준은 Q / (2 max_j q_j) 입니다.
// 이것은 정확성 검증일 뿐이며 보안 증명이 아닙니다.

import (
	"fmt"
	"math"
	"math/big"
	"time"

	"github.com/tuneinsight/lattigo/v6/ring"
)

// ---------------- 노이즈 0 암호화기 ----------------

// zeroSampler 는 오차를 항상 0 으로 만드는 ring.Sampler 입니다 (σ = 0 수학적 검증용).
type zeroSampler struct{ r *ring.Ring }

func (z zeroSampler) Read(p ring.Poly)               { p.Zero() }
func (z zeroSampler) ReadNew() ring.Poly             { return z.r.NewPoly() }
func (z zeroSampler) ReadAndAdd(ring.Poly)           {}
func (z zeroSampler) AtLevel(level int) ring.Sampler { return z }

// newEncryptorSigma 는 sigma = 0 이면 오차 없는 암호화기를, 아니면 기존 NewMKEncryptor 를 돌려줍니다.
func newEncryptorSigma(key *CGSWKey, sigma float64, op *MKOperator) (*MKEncryptor, error) {
	if sigma == 0 {
		e, err := NewMKEncryptor(key, 3.2, op)
		if err != nil {
			return nil, err
		}
		e.gs = zeroSampler{op.R}
		return e, nil
	}
	return NewMKEncryptor(key, sigma, op)
}

// ---------------- 잔차 측정 ----------------

// NoiseBound 는 가젯 복호 정확성 기준 Q / (2 max_j q_j).
func (op *MKOperator) NoiseBound() *big.Int {
	Q := op.R.ModulusAtLevel[op.R.Level()]
	qmax := uint64(0)
	for _, s := range op.R.SubRings {
		qmax = max(qmax, s.Modulus)
	}
	return new(big.Int).Quo(Q, new(big.Int).SetUint64(2*qmax))
}

// gadgetResidual 은 블록 X 의 잔차를 계산합니다. keys 는 행 r (마지막 마스크 행 제외) 의 비밀키,
// msgs[r] 는 행 r 의 메시지 (NTT, 영값 ring.Poly{} 이면 0) 입니다.
func (op *MKOperator) gadgetResidual(X GadgetBlock, keys, msgs []ring.Poly) *big.Int {
	r := op.R
	last := len(keys)
	ph, gm := r.NewPoly(), r.NewPoly()
	coeffs := make([]*big.Int, r.N())
	for i := range coeffs {
		coeffs[i] = new(big.Int)
	}
	worst := new(big.Int)
	for j := range X {
		for row := 0; row < last; row++ {
			ph.Copy(X[j][row])
			r.MulCoeffsBarrettThenAdd(X[j][last], keys[row], ph)
			if msgs[row].Coeffs != nil { // 영값 Poly = 메시지 0
				gm.Zero()
				op.mulRNSScalarThenAdd(msgs[row], op.Gvec[j], gm)
				r.Sub(ph, gm, ph)
			}
			r.INTT(ph, ph)
			r.PolyToBigintCentered(ph, 1, coeffs)
			for _, c := range coeffs {
				if c.CmpAbs(worst) > 0 {
					worst.Abs(c)
				}
			}
		}
	}
	return worst
}

// cgswBlockResiduals 는 cGSW_t(μ) 의 블록별 잔차 (길이 K+1) 를 돌려줍니다.
// 기대값: 블록 i<K 는 μ·g·e_i, 블록 K 는 μ·g·t.
func (op *MKOperator) cgswBlockResiduals(c *CGSW, t []ring.Poly, mu ring.Poly) []*big.Int {
	K := c.Parties()
	res := make([]*big.Int, K+1)
	for i := 0; i <= K; i++ {
		msgs := make([]ring.Poly, K)
		if i < K {
			msgs[i] = mu
		} else {
			for row := 0; row < K; row++ {
				msgs[row] = op.R.NewPoly()
				op.R.MulCoeffsBarrett(mu, t[row], msgs[row])
			}
		}
		res[i] = op.gadgetResidual(c.Block(i), t[:K], msgs)
	}
	return res
}

// upKeyResidual: S U_i − t_0 e_{i mod K} gᵀ (i ≥ K 이면 Π 가 한 바퀴 돌아 e_{i mod K}).
func (op *MKOperator) upKeyResidual(U GadgetBlock, t []ring.Poly, i int) *big.Int {
	K := len(t)
	msgs := make([]ring.Poly, K)
	msgs[i%K] = t[0]
	return op.gadgetResidual(U, t, msgs)
}

// hintResiduals 는 압축 힌트 각각의 잔차를 계산합니다.
//
//	S U_0 = t_0 e_0 gᵀ,  S H_⋆ = t_0 t gᵀ,  S H_Δ = (t − Π t) gᵀ,  (1, t_0)·V_h = t_h gᵀ
func (op *MKOperator) hintResiduals(ck *CompressedEvalKey, t []ring.Poly) (u0, hStar, hDelta, vmax *big.Int) {
	r, K := op.R, ck.K
	u0 = op.upKeyResidual(ck.U0, t, 0)

	ms := make([]ring.Poly, K)
	md := make([]ring.Poly, K)
	for row := 0; row < K; row++ {
		ms[row] = r.NewPoly()
		r.MulCoeffsBarrett(t[0], t[row], ms[row])
		md[row] = r.NewPoly()
		r.Sub(t[row], t[(row-1+K)%K], md[row])
	}
	hStar = op.gadgetResidual(ck.HStar, t, ms)
	hDelta = op.gadgetResidual(ck.HDelta, t, md)

	vmax = new(big.Int)
	for h := 0; h < K; h++ {
		e := op.gadgetResidual(ck.CLEV(h), t[:1], []ring.Poly{t[h]})
		if e.Cmp(vmax) > 0 {
			vmax = e
		}
	}
	return
}

// Log2Str 는 잔차를 "0" 또는 "2^x.xx" 로 표시합니다.
func Log2Str(x *big.Int) string {
	if x.Sign() == 0 {
		return "0"
	}
	f, _ := new(big.Float).SetInt(x).Float64()
	return fmt.Sprintf("2^%.2f", math.Log2(f))
}

func maxBig(xs ...*big.Int) *big.Int {
	m := new(big.Int)
	for _, x := range xs {
		if x.Cmp(m) > 0 {
			m = x
		}
	}
	return m
}

// mgswResidual 은 MGSW 의 결합 복호 잔차입니다 (checkMGSW 와 같은 기대값).
//
//	basket 0, 자리 j: JointDecrypt = g_j·m,   basket i ≥ 1: JointDecrypt = g_j·m·jointkey_i
func (op *MKOperator) mgswResidual(jd *JointDecryptor, m *MGSW, mu uint64) *big.Int {
	r := op.R
	coeffs := make([]*big.Int, r.N())
	for i := range coeffs {
		coeffs[i] = new(big.Int)
	}
	c := make([]uint64, len(r.SubRings))
	exp := r.NewPoly()
	worst := new(big.Int)
	for i := range m.Basket {
		for j := 0; j < op.Glen; j++ {
			for l, s := range r.SubRings {
				c[l] = mulmod(op.Gvec[j][l], mu%s.Modulus, s.Modulus)
			}
			exp.Zero()
			if i == 0 {
				op.addRNSScalar(exp, c)
			} else {
				op.mulRNSScalarThenAdd(jd.JointKey[i-1], c, exp)
			}
			p := jd.JointDecrypt(m.Basket[i].Stack[j])
			r.Sub(p, exp, p)
			r.INTT(p, p)
			r.PolyToBigintCentered(p, 1, coeffs)
			for _, x := range coeffs {
				if x.CmpAbs(worst) > 0 {
					worst.Abs(x)
				}
			}
		}
	}
	return worst
}

// ---------------- 통합 시나리오 ----------------

// CompressionReport 는 EvalKeyCompressionTest 의 측정 결과입니다. 잔차는 모든 당사자에 대한 최댓값입니다.
type CompressionReport struct {
	K, Glen, N int
	Sigma      float64
	LogQ       float64
	Bound      *big.Int

	FullPolys, CompPolys int // 당사자 1명 기준 환 원소 수
	FullBytes, CompBytes int // 당사자 1명 기준 직렬화 바이트

	KeyGenTime, CompressTime, ReconstructTime time.Duration // 당사자 1명 평균
	UpKeyTime                                 time.Duration // 당사자 1명 평균 (K−1 회 R 적용)

	HintU0, HintHStar, HintHDelta, HintV *big.Int
	UpNoise                              []*big.Int // 재귀 깊이 i = 0..depth−1 의 S U_i 잔차
	OrigBlock, RecBlock                  []*big.Int // 블록 i = 0..K 별 잔차 (모든 h, 모든 당사자 최댓값)

	Orig, Rec EvalResult // 원본 키 / 복원 키로 돌린 기존 Eval
}

// EvalResult 는 기존 Eval (MGSWExtend 전 당사자 + InternalProduct 1 회) 의 시간과 MGSW 복호 잔차입니다.
type EvalResult struct {
	ExtendTime, IPTime   time.Duration
	ExtendNoise, IPNoise *big.Int
	IPErr                error // InternalProduct 결과의 checkMGSW 실패 (nil 이면 정상 복호)
}

// EvalKeyCompressionTest 는 기존 KeyGen → 압축 → 서버 복원 → cGSW 잔차 비교 → 기존 Eval 투입을 수행합니다.
// sigma = 0 이면 오차 없는 암호화로 정확한 등식 (잔차 0) 을 검사하고, 아니면 잔차 < Q/(2 q_max) 를 검사합니다.
// depth 는 Up key 재귀 노이즈 측정 깊이 (K 미만이면 K 로 맞춤).
func EvalKeyCompressionTest(testID int, sigma float64, depth int) (*CompressionReport, error) {
	op, K, err := setup(testID)
	if err != nil {
		return nil, err
	}
	depth = max(depth, K)
	rep := &CompressionReport{K: K, Glen: op.Glen, N: op.R.N(), Sigma: sigma, Bound: op.NoiseBound()}
	f, _ := new(big.Float).SetInt(op.R.ModulusAtLevel[op.R.Level()]).Float64()
	rep.LogQ = math.Log2(f)
	exact := sigma == 0
	ok := func(x *big.Int) bool {
		if exact {
			return x.Sign() == 0
		}
		return x.Cmp(rep.Bound) < 0
	}

	keys := make([]*CGSWKey, K)
	entors := make([]*MKEncryptor, K)
	for i := range keys {
		if keys[i], err = TernaryCGSWKey(op.R, K, HW); err != nil {
			return nil, err
		}
		if entors[i], err = newEncryptorSigma(keys[i], sigma, op); err != nil {
			return nil, err
		}
	}

	// 1. 기존 KeyGen
	evks := make([]*EvalKey, K)
	t0 := time.Now()
	for p := range evks {
		evks[p] = entors[p].EvkEncrypt(K)
	}
	rep.KeyGenTime = time.Since(t0) / time.Duration(K)

	// 2. Compression
	cks := make([]*CompressedEvalKey, K)
	t0 = time.Now()
	for p := range cks {
		if cks[p], err = op.CompressEvalKey(evks[p]); err != nil {
			return nil, err
		}
	}
	rep.CompressTime = time.Since(t0) / time.Duration(K)

	rep.FullPolys = K * (K + 1) * (K + 1) * op.Glen
	rep.CompPolys = cks[0].NumPolys()
	if rep.CompPolys != (5*K+1)*op.Glen {
		return nil, fmt.Errorf("압축키 원소 수 %d ≠ (5K+1)d = %d", rep.CompPolys, (5*K+1)*op.Glen)
	}
	fb, err := EvalKeyMarshalBinary(evks[0])
	if err != nil {
		return nil, err
	}
	cb, err := cks[0].MarshalBinary()
	if err != nil {
		return nil, err
	}
	rep.FullBytes, rep.CompBytes = len(fb), len(cb)

	// 3. 힌트 잔차, 원본 G_h 잔차
	rep.HintU0, rep.HintHStar, rep.HintHDelta, rep.HintV = new(big.Int), new(big.Int), new(big.Int), new(big.Int)
	rep.OrigBlock = make([]*big.Int, K+1)
	rep.RecBlock = make([]*big.Int, K+1)
	for i := range rep.OrigBlock {
		rep.OrigBlock[i], rep.RecBlock[i] = new(big.Int), new(big.Int)
	}
	for p := range cks {
		t := keys[p].Val
		u, s, d, v := op.hintResiduals(cks[p], t)
		rep.HintU0 = maxBig(rep.HintU0, u)
		rep.HintHStar = maxBig(rep.HintHStar, s)
		rep.HintHDelta = maxBig(rep.HintHDelta, d)
		rep.HintV = maxBig(rep.HintV, v)
		for h := 0; h < K; h++ {
			for i, e := range op.cgswBlockResiduals(evks[p].Vals[h], t, t[h]) {
				rep.OrigBlock[i] = maxBig(rep.OrigBlock[i], e)
			}
		}
	}
	for name, x := range map[string]*big.Int{"U_0": rep.HintU0, "H_*": rep.HintHStar, "H_Δ": rep.HintHDelta, "V_h": rep.HintV} {
		if !ok(x) {
			return rep, fmt.Errorf("FAIL: 압축 힌트 %s 잔차 %s", name, Log2Str(x))
		}
	}

	// 4. 재귀 Up key 와 깊이별 노이즈
	rep.UpNoise = make([]*big.Int, depth)
	for i := range rep.UpNoise {
		rep.UpNoise[i] = new(big.Int)
	}
	var upTime time.Duration
	for p := range cks {
		t0 = time.Now()
		ups := op.UpKeys(cks[p], K)
		upTime += time.Since(t0)
		// 측정용으로 깊이 depth 까지 계속 적용 (i ≥ K 는 순환 후 노이즈 증가 관찰용)
		for i := K; i < depth; i++ {
			next := NewGadgetBlock(op.R, K+1, op.Glen)
			op.UpKeyStep(next, ups[i-1], cks[p].HDelta)
			ups = append(ups, next)
		}
		for i, U := range ups {
			rep.UpNoise[i] = maxBig(rep.UpNoise[i], op.upKeyResidual(U, keys[p].Val, i))
		}
	}
	rep.UpKeyTime = upTime / time.Duration(K)
	for i := 0; i < K; i++ {
		if !ok(rep.UpNoise[i]) {
			return rep, fmt.Errorf("FAIL: Up key U_%d 잔차 %s", i, Log2Str(rep.UpNoise[i]))
		}
	}

	// 5. 기존 Eval (원본 키): MGSWExtend 전 당사자 + InternalProduct.
	// 메시지는 모두 1 로 둡니다 (m·e 항이 사라지지 않는 최악의 경우라 노이즈 측정이 의미 있음).
	msgs := make([]uint64, K)
	cts := make([]*CGSW, K)
	for i := range msgs {
		msgs[i] = 1
		cts[i] = entors[i].CGSWEncryptConst(msgs[i], K)
		cts[i].Shrink(i + 1)
	}
	jd := NewJointDecryptor(keys, op)
	// runEval 은 실패해도 끝까지 측정합니다. MGSWExtend 실패는 err 로, InternalProduct 실패는 r.IPErr 로 돌려줍니다.
	// (복원 키는 노이즈가 커서 기존 파라미터에서는 InternalProduct 예산을 넘을 수 있으므로, 그 결과는
	//  오류로 숨기지 않고 리포트에 기록합니다. 원본 키의 InternalProduct 실패는 오류로 취급합니다.)
	runEval := func(ev []*EvalKey, tag string) (r EvalResult, err error) {
		for i := 1; i <= K; i++ {
			ev[i-1].Shrink(i)
		}
		mg := make([]*MGSW, K)
		for i := range mg {
			mg[i] = NewMGSW(op.R, K, op.Glen)
		}
		t0 := time.Now()
		for i := 1; i <= K; i++ {
			op.MGSWExtend(mg[i-1], cts[i-1], i, ev)
		}
		r.ExtendTime = time.Since(t0)
		r.ExtendNoise = new(big.Int)
		for i := 1; i <= K; i++ {
			r.ExtendNoise = maxBig(r.ExtendNoise, op.mgswResidual(jd, mg[i-1], msgs[i-1]))
			if e := checkMGSW(op, jd, mg[i-1], msgs[i-1], K, fmt.Sprintf("%s RGSW EXTENSION (party %d)", tag, i)); e != nil && err == nil {
				err = e
			}
		}
		i1, i2 := 1, min(2, K)
		res := NewMGSW(op.R, K, op.Glen)
		t0 = time.Now()
		op.InternalProduct(res, mg[i1-1], mg[i2-1])
		r.IPTime = time.Since(t0)
		r.IPNoise = op.mgswResidual(jd, res, msgs[i1-1]*msgs[i2-1])
		r.IPErr = checkMGSW(op, jd, res, msgs[i1-1]*msgs[i2-1], K, tag+" RGSW INTERNAL PRODUCT")
		return
	}
	var evalErr error
	rep.Orig, evalErr = runEval(evks, "ORIGINAL")
	if evalErr == nil {
		evalErr = rep.Orig.IPErr
	}
	evks = nil // 피크 메모리를 줄이기 위해 원본 키를 버림 (이후로는 압축키만 사용)

	// 6. 서버 복원: 압축키만으로 Ĝ_h 생성
	recs := make([]*EvalKey, K)
	t0 = time.Now()
	for p := range recs {
		recs[p] = op.DecompressEvalKey(cks[p])
	}
	rep.ReconstructTime = time.Since(t0) / time.Duration(K)

	// 7. 복원 키의 모든 가젯 블록 검사
	for p := range recs {
		t := keys[p].Val
		for h := 0; h < K; h++ {
			if got := recs[p].Vals[h]; got.Parties() != K || len(got.Stack) != (K+1)*op.Glen {
				return rep, fmt.Errorf("FAIL: Ĝ_%d 크기 불일치", h)
			}
			for i, e := range op.cgswBlockResiduals(recs[p].Vals[h], t, t[h]) {
				rep.RecBlock[i] = maxBig(rep.RecBlock[i], e)
				if !ok(e) {
					return rep, fmt.Errorf("FAIL: party %d Ĝ_%d 블록 %d 잔차 %s", p+1, h, i, Log2Str(e))
				}
			}
		}
	}

	// 8. 복원 키로 기존 Eval
	rep.Rec, err = runEval(recs, "RECONSTRUCTED")
	if evalErr != nil {
		return rep, evalErr
	}
	return rep, err
}

func (rep *CompressionReport) String() string {
	s := fmt.Sprintf("K=%d  N=%d  d=%d  logQ=%.1f  σ=%g  판정기준 Q/(2q_max)=%s\n",
		rep.K, rep.N, rep.Glen, rep.LogQ, rep.Sigma, Log2Str(rep.Bound))
	s += fmt.Sprintf("[통신량/당사자] 원소: full K(K+1)²d=%d  comp (5K+1)d=%d  (×%.2f)\n",
		rep.FullPolys, rep.CompPolys, float64(rep.FullPolys)/float64(rep.CompPolys))
	s += fmt.Sprintf("               직렬화: full %.2f MiB  comp %.2f MiB  (×%.2f)\n",
		float64(rep.FullBytes)/(1<<20), float64(rep.CompBytes)/(1<<20), float64(rep.FullBytes)/float64(rep.CompBytes))
	s += fmt.Sprintf("[시간/당사자] KeyGen %v  Compress %v  UpKeys %v  Reconstruct(전체) %v\n",
		rep.KeyGenTime, rep.CompressTime, rep.UpKeyTime, rep.ReconstructTime)
	s += fmt.Sprintf("[힌트 잔차] U_0 %s  H_* %s  H_Δ %s  V_h %s\n",
		Log2Str(rep.HintU0), Log2Str(rep.HintHStar), Log2Str(rep.HintHDelta), Log2Str(rep.HintV))
	s += "[재귀 깊이별 S·U_i 잔차]"
	for i, e := range rep.UpNoise {
		mark := ""
		if i >= rep.K {
			mark = "*"
		}
		s += fmt.Sprintf("  i=%d%s:%s", i, mark, Log2Str(e))
	}
	s += "   (* = 순환 이후, 측정용)\n[블록별 cGSW 잔차] 원본 G_h → 복원 Ĝ_h\n"
	for i := range rep.OrigBlock {
		s += fmt.Sprintf("  block %d: %s → %s\n", i, Log2Str(rep.OrigBlock[i]), Log2Str(rep.RecBlock[i]))
	}
	str := func(x *big.Int) string {
		if x == nil {
			return "-"
		}
		return Log2Str(x)
	}
	s += "[Eval]                원본 키                 복원 키\n"
	s += fmt.Sprintf("  MGSWExtend 시간     %-23v %v\n", rep.Orig.ExtendTime, rep.Rec.ExtendTime)
	s += fmt.Sprintf("  MGSWExtend 잔차     %-23s %s\n", str(rep.Orig.ExtendNoise), str(rep.Rec.ExtendNoise))
	s += fmt.Sprintf("  InternalProduct 시간 %-23v %v\n", rep.Orig.IPTime, rep.Rec.IPTime)
	s += fmt.Sprintf("  InternalProduct 잔차 %-23s %s\n", str(rep.Orig.IPNoise), str(rep.Rec.IPNoise))
	stat := func(e error) string {
		if e != nil {
			return "복호 실패 (잔차 > 기준)"
		}
		return "정상"
	}
	s += fmt.Sprintf("  InternalProduct 복호 %-23s %s\n", stat(rep.Orig.IPErr), stat(rep.Rec.IPErr))
	return s
}
