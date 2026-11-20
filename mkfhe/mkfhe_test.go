package mkfhe

import "testing"

// go test ./... 는 K=2 로 원본 main.jl 의 테스트 4개를 실행합니다.
func TestK2(t *testing.T) {
	for name, f := range map[string]func(int) error{
		"VectorGadgetProd":   VectorGadgetProdTest,
		"MKVectorGadgetProd": MKVectorGadgetProdTest,
		"RGSWExtend":         RGSWExtendTest,
	} {
		if err := f(2); err != nil {
			t.Fatalf("%s: %v", name, err)
		}
	}
	for _, distinct := range []bool{false, true} {
		if err := RGSWInternalProductTest(2, distinct); err != nil {
			t.Fatalf("InternalProduct(distinct=%v): %v", distinct, err)
		}
	}
}

func TestPrimesK2(t *testing.T) {
	p, err := NewMKParameters(RecommendedSketch[2])
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("Q = %v", p.Q)
}
