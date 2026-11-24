// 평가키 압축/복원 측정 (신규). 사용법: go run ./cmd/evkcompress -test 4 -depth 12
//
// 파라미터 스윕: -test 의 프리셋을 기본으로 두고 -K, -d (= Qlen, RNS 소수 개수 = 가젯 길이), -logq, -logn 으로 덮어씁니다.
//
//	go run ./cmd/evkcompress -K 3 -d 8 -logq 214 -logn 13
//	go run ./cmd/evkcompress -K 2 -d 7 -sigma 3.2 -csv   # 결과를 CSV 한 줄로 (헤더는 -csvheader)
package main

import (
	"flag"
	"fmt"
	"os"
	"strconv"
	"strings"

	"github.com/yourname/crsless-go/mkfhe"
)

func main() {
	testID := flag.Int("test", 2, "기본 프리셋 (2, 4, 8)")
	depth := flag.Int("depth", 0, "Up key 재귀 노이즈 측정 깊이 (기본 K, K 이상은 순환 후 측정용)")
	K := flag.Int("K", 0, "당사자 수 (0 이면 프리셋 값)")
	d := flag.Int("d", 0, "가젯 길이 d = RNS 소수 개수 Qlen (0 이면 프리셋 값). 소수 하나 ≈ logq/d 비트")
	logQ := flag.Float64("logq", 0, "Q 전체 비트 수 목표 (0 이면 프리셋 값)")
	logN := flag.Int("logn", 0, "링 차수 log2 N (0 이면 프리셋 값). M = 2N")
	sigmas := flag.String("sigma", "0,3.2", "쉼표로 구분한 오차 표준편차 목록 (0 = 정확한 등식 검사)")
	csv := flag.Bool("csv", false, "사람용 리포트 대신 CSV 한 줄씩 출력")
	csvHeader := flag.Bool("csvheader", false, "CSV 헤더만 출력하고 종료")
	flag.Parse()

	if *csvHeader {
		fmt.Println(mkfhe.CSVHeader)
		return
	}

	sketch, ok := mkfhe.RecommendedSketch[*testID]
	if !ok {
		fmt.Printf("프리셋 -test %d 없음\n", *testID)
		os.Exit(1)
	}
	if *K > 0 {
		sketch.MaxParties = *K
	}
	if *d > 0 {
		sketch.Qlen = *d
	}
	if *logQ > 0 {
		sketch.LogQ = *logQ
	}
	if *logN > 0 {
		sketch.M = 2 << *logN
	}

	ipFail := false
	for _, f := range strings.Split(*sigmas, ",") {
		sigma, err := strconv.ParseFloat(strings.TrimSpace(f), 64)
		if err != nil {
			fmt.Printf("-sigma 값 %q 를 읽을 수 없음\n", f)
			os.Exit(1)
		}
		rep, err := mkfhe.EvalKeyCompressionTestSketch(sketch, sigma, *depth)
		if *csv {
			fmt.Println(rep.CSVRow(sketch, sigma, err))
			continue
		}
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
	if *csv {
		return
	}
	if ipFail {
		fmt.Println("압축/복원/cGSW/MGSWExtend 검증 통과, 복원 키 InternalProduct 는 노이즈 예산 초과")
		os.Exit(2)
	}
	fmt.Println("ALL PASSED")
}
