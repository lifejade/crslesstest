package mkfhe

// HW 는 원본 main.jl / paramfinder.jl 의 global hw = 192.
const HW = 192

// RecommendedSketch 는 원본의 recommended_sketch.
// (원본 주석: 각 Qi 가 작아야 하므로 logQ = 100 으로 두며, 실행하면 logQ ≈ 105 가 됩니다.)
var RecommendedSketch = map[int]MKParamSketch{
	2:  {MaxParties: 2, M: 1 << 13, LogQ: 100, Qlen: 6},
	4:  {MaxParties: 4, M: 1 << 14, LogQ: 214, Qlen: 5},
	8:  {MaxParties: 8, M: 1 << 14, LogQ: 214, Qlen: 6},
	16: {MaxParties: 16, M: 1 << 14, LogQ: 214, Qlen: 8}, // paramfinder.jl 전용
}
