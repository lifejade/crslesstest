# crsless-go

SNUCP/crsless (Julia, HIENAA 기반 다중키 FHE) 를 Go + Lattigo v6.2.0 `ring` 으로 옮긴 프로젝트. 자세한 대응은 README.md.

## 명령
- 테스트: `go test ./...` (K=2)
- 실행: `go run ./cmd/crsless -test 2` / `-test 4 -distinct`
- 커밋 전: `gofmt -w . && go vet ./...`

## 규칙
- 모든 ring.Poly 는 NTT 형식으로 유지. 계수 형식이 필요하면 복사본에 INTT.
- 원본 Julia 와 1:1 대응을 유지 (함수 주석에 원본 이름 표기). 당사자 인덱스는 1부터.
- K=8 은 메모리 10GB 이상 필요.

## 원본 참조
`git clone https://github.com/SNUCP/crsless` 와 `https://github.com/snu-lukemin/HIENAA.jl` (src/rlwe/decomposition.jl, src/math/basischange.jl).
