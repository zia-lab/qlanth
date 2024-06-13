Options[ParseLiYF4] = {"Export" -> True};
ParseLiYF4[sheetFname_]:=Module[
    {sheetNames, liyf4ExpData, sheetData, header,
    expCol, excludeCol, stateLabels, expEnergies,
    lastRow, shN},
    (
        sheetNames   = Select[Import[sheetFname, "Sheets"], StringContainsQ[#, "3+"] &];
        liyf4ExpData = Table[
            (
                sheetData  = Import[sheetFname, {"Sheets", shN}];
                header     = sheetData[[1]];
                sheetData  = sheetData[[2 ;;]];
                expCol     = Position[header, "Exp (1/cm)"][[1, 1]];
                excludeCol = Position[header, "Exp:Excluded"][[1, 1]];
                (*first column is always the SLJ coarse label of the state, 
                keep it*)
                sheetData   = Transpose[sheetData];
                stateLabels = sheetData[[1]];
                expEnergies = sheetData[[expCol]];
                excludeCol  = sheetData[[excludeCol]];
                lastRow     = Position[stateLabels, ""];
                lastRow = If[Length[lastRow] > 0,
                    lastRow[[1, 1]] - 1,
                    -1
                ];
                StringTake[shN, {1, 2}] -> Transpose[{expEnergies, stateLabels, excludeCol}][[;; lastRow]]
            ),
            {shN, sheetNames}
        ];
        liyf4ExpData    = Association[liyf4ExpData];
        If[OptionValue["Export"],
            Export[FileNameJoin[{NotebookDirectory[], "data", "LiYF4_exp_data.m"}], liyf4ExpData]
        ];
        Return[liyf4ExpData];
    );
]