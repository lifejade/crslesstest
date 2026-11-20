// 원본 main.jl 의 Go 포팅. 사용법: go run ./cmd/crsless -test 2
package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/yourname/crsless-go/mkfhe"
)

func main() {
	testID := flag.Int("test", 2, "당사자 수 (2, 4, 8)")
	distinct := flag.Bool("distinct", false, "내부곱 테스트에서 서로 다른 두 당사자를 곱함 (원본은 같은 당사자)")
	flag.Parse()

	run := func(err error) {
		if err != nil {
			fmt.Println(err)
			os.Exit(1)
		}
	}
	run(mkfhe.VectorGadgetProdTest(*testID))
	run(mkfhe.MKVectorGadgetProdTest(*testID))
	run(mkfhe.RGSWExtendTest(*testID))
	run(mkfhe.RGSWInternalProductTest(*testID, *distinct))
	fmt.Println("ALL PASSED")
}
