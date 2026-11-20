// 평가키 압축/복원 측정 (신규). 사용법: go run ./cmd/evkcompress -test 4 -depth 12
package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/yourname/crsless-go/mkfhe"
)

func main() {
	testID := flag.Int("test", 2, "당사자 수 (2, 4, 8)")
	depth := flag.Int("depth", 0, "Up key 재귀 노이즈 측정 깊이 (기본 K, K 이상은 순환 후 측정용)")
	flag.Parse()

	ipFail := false
	for _, sigma := range []float64{0, 3.2} {
		rep, err := mkfhe.EvalKeyCompressionTest(*testID, sigma, *depth)
		if rep != nil {
			fmt.Print(rep.String(), "\n")
		}
		if err != nil {
			fmt.Println(err)
			os.Exit(1)
		}
		if rep.Rec.IPErr != nil {
			ipFail = true
			fmt.Printf("한계 (σ=%g): 복원 키 InternalProduct 복호 실패 — %v\n\n", sigma, rep.Rec.IPErr)
		}
	}
	if ipFail {
		fmt.Println("압축/복원/cGSW/MGSWExtend 검증 통과, 복원 키 InternalProduct 는 노이즈 예산 초과")
		os.Exit(2)
	}
	fmt.Println("ALL PASSED")
}
