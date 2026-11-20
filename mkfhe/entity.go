package mkfhe

// 원본 src/entity.jl 에 대응합니다.
// 모든 다항식은 Lattigo ring.Poly 로, 항상 NTT(평가) 형식으로 보관합니다.

import "github.com/tuneinsight/lattigo/v6/ring"

// Tensor 는 HIENAA 의 Tensor (다항식 벡터) 입니다. t[:i] 슬라이스는 Julia 의 뷰와 같이 데이터를 공유합니다.
type Tensor []ring.Poly

func NewTensor(r *ring.Ring, k int) Tensor {
	t := make(Tensor, k)
	for i := range t {
		t[i] = r.NewPoly()
	}
	return t
}

func (t Tensor) Zero() {
	for i := range t {
		t[i].Zero()
	}
}

// CopyFrom 은 Julia 의 copy!(dst, src) 입니다.
func (t Tensor) CopyFrom(src Tensor) {
	for i := range t {
		t[i].Copy(src[i])
	}
}

// CGSW: glen*(k+1) 개의 Tensor, 각 Tensor 는 k+1 개 다항식 (b_1..b_k, a).
type CGSW struct {
	Glen  int
	Stack []Tensor
}

func NewCGSW(r *ring.Ring, k, glen int) *CGSW {
	c := &CGSW{Glen: glen, Stack: make([]Tensor, glen*(k+1))}
	for i := range c.Stack {
		c.Stack[i] = NewTensor(r, k+1)
	}
	return c
}

// Parties 는 현재 cGSW 가 다루는 당사자 수 k 입니다.
func (c *CGSW) Parties() int { return len(c.Stack[0]) - 1 }

// Shrink 는 shrink!(c, k): K 당사자용 cGSW 를 k 당사자용으로 줄입니다.
func (c *CGSW) Shrink(k int) {
	K, glen := c.Parties(), c.Glen
	for i := 0; i < glen*k; i++ {
		c.Stack[i][k].Copy(c.Stack[i][K])
		c.Stack[i] = c.Stack[i][:k+1]
	}
	for i := 0; i < glen; i++ {
		src := c.Stack[glen*K+i]
		src[k].Copy(src[K])
		c.Stack[glen*k+i] = src[:k+1]
	}
	c.Stack = c.Stack[:glen*(k+1)]
}

// EvalKey: k 개의 cGSW.
type EvalKey struct {
	Vals []*CGSW
}

func NewEvalKey(r *ring.Ring, k, glen int) *EvalKey {
	e := &EvalKey{Vals: make([]*CGSW, k)}
	for i := range e.Vals {
		e.Vals[i] = NewCGSW(r, k, glen)
	}
	return e
}

func (e *EvalKey) Shrink(k int) {
	for _, v := range e.Vals {
		v.Shrink(k)
	}
}

// MLEV: glen 개의 Tensor (각 k+1 개 다항식).
type MLEV struct {
	Stack []Tensor
}

func NewMLEV(r *ring.Ring, k, glen int) *MLEV {
	m := &MLEV{Stack: make([]Tensor, glen)}
	for i := range m.Stack {
		m.Stack[i] = NewTensor(r, k+1)
	}
	return m
}

// MGSW: k+1 개의 MLEV.
type MGSW struct {
	Basket []*MLEV
}

func NewMGSW(r *ring.Ring, k, glen int) *MGSW {
	m := &MGSW{Basket: make([]*MLEV, k+1)}
	for i := range m.Basket {
		m.Basket[i] = NewMLEV(r, k, glen)
	}
	return m
}
