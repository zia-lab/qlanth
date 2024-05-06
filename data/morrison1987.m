morrison1987 = <|
	"source" -> "Morrison 1987, Angular Momentum Theory Applied to Interactions in Solids",
	"groupAliases" -> <|
		"Ci" -> {"S2"},
		"Cs" -> {"C1h"},
		"D2" -> {"V"},
		"D2h" -> {"Vh"},
		"D2d" -> {"Vd"},
		"C3i" -> {"S6"}
	|>,
	"tab8p1" ->
		"C2, x00xxxx00xx00xxxxx
Cs, 0xxxxx0xx00xxxxxxx
C2h,000xxx0000000xxxxx
D2, 000xx00000x00xx0x0
C2v,x00xx0x00x000xx0x0
D2h,000xx00000000xx0x0
C4, x00x00x000000x00xx
S4, 000x00000xx00x00xx
C4h,000x000000000x00xx
D4, 000x000000000x00x0
C4v,x00x00x000000x00x0
D2d,000x00000x000x00x0
D4h,000x000000000x00x0",
	"tab8p1params"-> {B01, B11, S11, B02, B22, S22, B03, B13, S13, B23, S23, B33, S33, B04, B24, S24, B44, S44},
	"tab8p2" ->
		"C3, xxxxxxxx
C3i,0x000xxx
D3, 0x00xxx0
C3v,xxxx0xx0
D3d,0x000xx0
C6, xxx00x00
C3h,0xxx0x00
C6h,0x000x00
D6, 0x000x00
C6v,xxx00x00
D3h,0xx00x00
D6h,0x000x00",
	"tab8p2params" -> {B01, B02, B03, B33, S33, B04, B34, S34},
	"T" -> <|
		"BqkSqk" -> {B23, B04, B44, B06, B46},
		"constraints" -> {{B44 -> Sqrt[5/14] B04, B46 -> -Sqrt[(7/2)] B06},
		{B44 -> -Sqrt[(5/14)] B04, B46 -> Sqrt[7/2] B06}}
	|>,
	"C1" -> <|
		"BqkSqk" -> {B01, B11, S11, B02, B22, S22, B03, B13, S13, B23, S23, B33, S33, B04, B24, S24, B44, S44},
		"constraints" -> {}
	|>,
	"Ci" -> <|
		"BqkSqk" -> {B02, B22, S22, B04, B24, S24, B44, S44},
		"constraints" -> {}
	|>
|>