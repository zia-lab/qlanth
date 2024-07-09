(* This script has two function, used to fit for the data of LaF3 and LiYF4 *)

FitLaF3[]:=(
    If[
        Not[ValueQ[PrintToOutputNb]],
        PrintToOutputNb = HelperNotebook["Lanthanides : OUT"];
    ];
    PrintAndDefine[message_] := (
        progressMessage = ToString[ln] <> ": " <> ToString[message];
        PrintToOutputNb[message];
    );
    Off[N::preclg];
    withSpinSin = 1;
    restart = False;
    saveEigenvectors = False;
    filePrefix = If[withSpinSin == 1,
        "pathfinder-constraints-with-spinspin-and-no-truncation",
        "pathfinder-constraints-with-spinspin-and-no-truncation"
    ];
    \[Sigma]exp = 1.0;
    trunEnergy = Infinity;
    fitOrder = "Pr Nd Dy Ce Sm Ho Er Tm Yb Tb Eu Gd Pm";
    fitOrder = Position[theLanthanides, #][[1, 1]] & /@ StringSplit[fitOrder, " "];
    doThese = fitOrder;
    
    linearParams       = Join[casimirSymbols, cfSymbols, TSymbols, pseudoMagneticSymbols, marvinSymbols];
    quadraticParams    = {};
    interpolatedParams = Join[slaterSymbols,{\[Zeta]}];
    
    If[Not[restart],
    varModel = AssociationThread[{B02, B04, B06, B22, B24, B26, B44, B46, B66, F2,
        F4, F6, M0, P2, T2, T3, T4, T6, T7, 
        T8, \[Alpha], \[Beta], \[Gamma], \[Zeta]}, 
        ConstantArray[<|"data" -> {}, "function" -> Null|>, 24]];
    truncatedSols = <||>;
    ];
    
    Do[
    (
        numE = numE0;
        ln = theLanthanides[[numE]];
        PrintAndDefine[Style["Working on "<>ln<>" ...", Red]];
        If[ln == "Pm",
            (
                pmSol = <||>;
                pmSol["hamDim"] = Binomial[14, numE];
                pmSol["numE"] = numE;
                pmSol["problemVars"] = variedSymbols["Nd"];
                pmSol["bestParamsWithConstraints"] = <||>;
                pmSol["bestRMS"] = 0;
                pmSol["bestRMSwithTruncation"] = 0;
                pmSol["fittedLevels"] = 0;
                pmSol["solWithUncertainty"] = {};
                Do[
                (
                    bestP = varModel[pVar]["function"][numE];
                    pmSol["bestParamsWithConstraints"][pVar] = bestP;
                ),
                {pVar, variedSymbols["Nd"]}];
                pmSol["bestParamsWithConstraints"][\[Epsilon]] = 0;
                constrainedMsandPs = ("MagneticSimplifier" /. 
                    Options[ClassicalFit]) /. 
                pmSol["bestParamsWithConstraints"];
                pmSol["bestParamsWithConstraints"] = Join[pmSol["bestParamsWithConstraints"], 
                    Association@constrainedMsandPs];
                
                pmSol["bestParams"] = Normal@pmSol["bestParamsWithConstraints"];
                pmSol["bestRMSwithFullDiagonalization"] = 0;
                pmSol["constraints"] = {};
                LogSol[pmSol, filePrefix <> "/" <> "Pm" <> "-" <> ToString[trunEnergy]];
                truncatedSols[numE] = pmSol;
                Continue[];
            )
            ];
        pVars = variedSymbols[ln];
        constraints = Association@caseConstraints[ln];
        fixedConstraints = Keys@Select[constraints, NumericQ];
        
        Do[constraints[fixedCon] = varModel[fixedCon]["function"][numE], 
        {fixedCon, fixedConstraints}];
        constraints = Normal[constraints];
        
        carnalKey = "appendix:" <> ln <> ":RawTable";
        expData = {#[[2]], #[[1]], #[[3]]} & /@ Carnall[carnalKey];
        If[OddQ[numE],
            expData = Flatten[{#, #} & /@ expData, 1]
            ];
        presentDataIndices = Flatten[Position[expData, {_?(NumericQ[#] &), _, _}]];
        excludedDataIndices = If[OddQ[numE],
            Flatten[
            Position[Flatten[{#, #} & /@ (Last /@ Carnall[carnalKey]), 1], 
            "not included"]],
            Flatten[Position[Last /@ Carnall[carnalKey], "not included"]]
            ];
        independentVars = Variables[pVars /. constraints];
        startValues = AssociationThread[independentVars, 
            independentVars /. LoadLaF3Parameters[ln]];
        thisSolution = ClassicalFit[numE,
            expData,
            excludedDataIndices,
            pVars,
            startValues,
            constraints,
            "SaveEigenvectors" -> saveEigenvectors,
            "AddConstantShift" -> True,
            "SignatureCheck" -> False,
            "Energy Uncertainty in K" -> \[Sigma]exp,
            "OtherSimplifier" -> (("OtherSimplifier" /. 
                Options[ClassicalFit]) /. (\[Sigma]SS -> _) :> (\[Sigma]SS ->
                    withSpinSin)),
            "MaxIterations" -> 1000,
            "SaveToLog" -> False,
            "MaxHistory" -> 1000,
            "FilePrefix" -> 
            filePrefix <> "/" <> ln <> "-" <> ToString[trunEnergy],
            "TruncationEnergy" -> trunEnergy,
            "PrintFun" -> PrintAndDefine];
        thisSolution["bestParams"] = Append[thisSolution["bestParams"], (nE -> numE)];
        thisSolution["bestParamsWithConstraints"][nE] = numE;
        If[(2 <= numE <= 12),
            (constrainedMsandPs = ("MagneticSimplifier" /. 
                Options[ClassicalFit]) /. thisSolution["bestParamsWithConstraints"];
            thisSolution["bestParamsWithConstraints"] = 
            Join[thisSolution["bestParamsWithConstraints"], 
            Association@constrainedMsandPs];
            )
            ];
        
        PrintAndDefine[">> Computing energies with full diagonalization ..."];
        params = thisSolution["bestParamsWithConstraints"];
        eigenEnergies = FastIonSolver[params,
            "EnergiesOnly" -> True,
            "MakeNotebook" -> False,
            "PrintFun" -> PrintAndDefine,
            "Append to Filename" -> "fitting-" <> ln,
            "Remove Kramers" -> False,
            "Energy Uncertainty in K" -> \[Sigma]exp,
            "Explorer" -> False
            ];
        
        PrintAndDefine[">> Calculating full diagonalization RMS ..."];
        fullDiffs = (eigenEnergies[[;; Length[expData]]] - 
            First /@ expData)[[
            Complement[presentDataIndices, excludedDataIndices]]];
        fullRMSsquared = Sqrt[Total[fullDiffs^2]/thisSolution["degreesOfFreedom"]];
        thisSolution["bestRMSwithTruncation"] = thisSolution["bestRMS"];
        thisSolution["energiesWithTruncation"] = thisSolution["energies"];
        thisSolution["excludedDataIndices"] = excludedDataIndices;
        thisSolution["bestRMSwithFullDiagonalization"] = fullRMSsquared;
        thisSolution["energiesWithFullDiagonalization"] = eigenEnergies;
        thisSolution = KeyDrop[thisSolution, {"bestRMS", "energies"}];
        truncatedSols[numE0] = thisSolution;
        LogSol[thisSolution, filePrefix <> "/" <> ln <> "-" <> ToString[trunEnergy]];
        
        PrintAndDefine["> Updating param model functions ..."];
        bestParams = Association@thisSolution["bestParams"];
        Do[
        (
            If[Not[MemberQ[Keys[bestParams], pVar]],
            Continue[];
            ];
            varModel[pVar]["data"] = 
            DeleteDuplicates@
            Append[varModel[pVar]["data"], {numE, bestParams[pVar]}];
            numData = Length[varModel[pVar]["data"]];
            If[numData > 0,
            (
                order = Which[
                    numData == 1,
                    0,
                    numData == 2,
                    1,
                    numData >= 3,
                    2];
                varModel[pVar]["function"] = Which[
                    MemberQ[linearParams, pVar],
                        NonlinearModelFit[varModel[pVar]["data"], a0 + a1 x, {a0, a1}, x],
                    MemberQ[quadraticParams, pVar],
                        NonlinearModelFit[varModel[pVar]["data"], a0 + a1 x + a2 x^2, {a0, a1, a2}, x],
                    MemberQ[interpolatedParams, pVar],
                        Interpolation[varModel[pVar]["data"], InterpolationOrder -> order]
                    ]
            )
            ];
        ),
        {pVar, pVars}
        ];
        
        Run["~/Scripts/pushover \"finished fitting " <> ln <> "\""];
    ),
    {numE0, doThese}];
    Return[{varModel, truncatedSols}];
)

FitLiYF4[]:=(
    If[
        Not[ValueQ[PrintToOutputNb]],
        PrintToOutputNb = HelperNotebook["Lanthanides : OUT"];
    ];
    PrintAndDefine[message_] := (
        progressMessage = ToString[ln] <> ": " <> ToString[message];
        PrintToOutputNb[message];
    );
    Off[N::preclg];
    withSpinSin = 1;
    \[Sigma]exp = 1.0;
    restart = False;
    maxData = Infinity;
    saveEigenvectors = False;
    truncationEnergy = Infinity;
    filePrefix = StringJoin[{
                "liyf4-pathfinder-with-",
                ToLowerCase[ToString[truncationEnergy]],
                "-truncation-nominal-sigma-",
                If[withSpinSin == 1, "with-spinspin", "without-spinspin"]
                }];
    (* This according to Cheng probably following Gorller-Walrand *)
    magSimplifier = {
        M2 -> 56/100 M0,
        M4 -> 38/100 M0,
        P4 -> 3/4 P2,
        P6 -> 1/2 P2};
    fitOrder = "Er Nd Eu Ho Sm Pr Tm Yb Ce Tb Dy Gd Pm";
    fitOrder = Position[theLanthanides, #][[1, 1]] & /@ 
    StringSplit[fitOrder, " "];
    doThese = fitOrder;
    
    linearParams       = Join[casimirSymbols, cfSymbols, TSymbols, pseudoMagneticSymbols, marvinSymbols];
    quadraticParams    = {};
    interpolatedParams = Join[slaterSymbols,{\[Zeta]}];

    If[Not[restart],
        (
            varModel = AssociationThread[{B02, B04, B06, B44, B46, F2, F4, F6, M0, P2, T2,
                T3, T4, T6, T7, T8, \[Alpha], \[Beta], \[Gamma], \[Zeta]}, 
                ConstantArray[<|"data" -> {}, "function" -> Null|>,
                    20]];
            truncatedSols = <||>;
        )
    ];
    
    Do[(
        numE = numE0;
        ln   = theLanthanides[[numE]];
        PrintAndDefine[Style["Working on "<>ln<>" ...", Red]];
        (* Cases in which we dont' fit data but interpolate for filling gaps in columns *)
        If[Or[ln == "Pm", ln == "Gd"],
            (
                pmSol           = <||>;
                pmSol["hamDim"] = Binomial[14, numE];
                pmSol["numE"]   = numE;
                pmSol["problemVars"] = variedSymbols["Nd"];
                pmSol["bestParamsWithConstraints"] = <||>;
                pmSol["bestRMS"]               = 0;
                pmSol["bestRMSwithTruncation"] = 0;
                pmSol["fittedLevels"]          = 0;
                pmSol["solWithUncertainty"]    = {};
                Do[
                    (
                        bestP         = varModel[pVar]["function"][numE];
                        pmSol["bestParamsWithConstraints"][pVar] = bestP;
                    ),
                    {pVar, variedSymbolsLiYF4["Er"]}
                ];
                pmSol["bestParamsWithConstraints"][\[Epsilon]] = 0;
                constrainedMsandPs = (magneticSimplifier) /. pmSol["bestParamsWithConstraints"];
                pmSol["bestParamsWithConstraints"] = Join[
                    pmSol["bestParamsWithConstraints"], 
                    Association@constrainedMsandPs];
                
                pmSol["bestParams"] = Normal@pmSol["bestParamsWithConstraints"];
                pmSol["bestRMSwithFullDiagonalization"] = 0;
                pmSol["constraints"] = {};
                LogSol[pmSol, 
                    filePrefix <> "/" <> ln <> "-" <> ToString[truncationEnergy]];
                truncatedSols[numE] = pmSol;
                Continue[];
            )
        ];
        
        pVars       = variedSymbolsLiYF4[ln];
        constraints = Association@caseConstraintsLiYF4[ln];
        
        fixedConstraints = Keys@Select[constraints, NumericQ];
        Do[
            fun = varModel[fixedCon]["function"];
            If[fun === Null,
                PrintAndDefine["No function for ", fixedCon, " in ", ln, ", exiting ..."];
                Return[Null]
            ]; 
            constraints[fixedCon] = varModel[fixedCon]["function"][numE],
            {fixedCon, fixedConstraints}
        ];
        constraints = Normal[constraints];
        
        expData    = liyoriteExpData[ln];
        uniqueData = DeleteDuplicates[First /@ expData];
        lastEnergy = uniqueData[[Min[Length[uniqueData], maxData]]];
        If[OddQ[numE],
            expData = Flatten[{#, #} & /@ expData, 1]
        ];
        presentDataIndices  = Flatten[Position[expData, {_?(NumericQ[#] &), _, _}]];
        excludedDataIndices = Pick[Range[1, Length[expData]], Not[#[[3]] == ""] & /@ expData];
        independentVars = Variables[pVars /. constraints];
        startValues = AssociationThread[
            independentVars, 
            independentVars /. paramsChengLiYF4[ln]];
        thisSolution = ClassicalFit[numE,
            expData,
            excludedDataIndices,
            pVars,
            startValues,
            constraints,
            "SaveEigenvectors" -> saveEigenvectors,
            "AddConstantShift" -> True,
            "SignatureCheck" -> False,
            "RefParamsVintage" -> "LiYF4",
            "Energy Uncertainty in K" -> \[Sigma]exp,
            "MaxIterations" -> 100,
            "MaxHistory" -> 1000,
            "FilePrefix" -> 
            filePrefix <> "/" <> ln <> "-" <> ToString[truncationEnergy],
            "TruncationEnergy" -> truncationEnergy,
            "PrintFun" -> PrintAndDefine,
            "SaveToLog" -> False,
            
            "OtherSimplifier" -> (("OtherSimplifier" /. 
                Options[ClassicalFit]) /. (\[Sigma]SS -> _) :> (\[Sigma]SS ->
                    withSpinSin)),
            
            "MagneticSimplifier" -> magSimplifier,
            "SymmetrySimplifier" -> symmetrySimplifier
            ];
        thisSolution["bestParams"] = Append[thisSolution["bestParams"], (nE -> numE)];
        thisSolution["bestParamsWithConstraints"][nE] = numE;
        
        If[(2 <= numE <= 12),
            (
                constrainedMsandPs = (magSimplifier) /. thisSolution["bestParamsWithConstraints"];
                thisSolution["bestParamsWithConstraints"] = Join[
                    thisSolution["bestParamsWithConstraints"], 
                    Association@constrainedMsandPs];
            )
        ];
        
        PrintAndDefine[">> Computing energies with full diagonalization ..."];
        params = thisSolution["bestParamsWithConstraints"];
        eigenEnergies = FastIonSolver[params,
            "EnergiesOnly" -> True,
            "MakeNotebook" -> False,
            "PrintFun" -> PrintAndDefine,
            "Append to Filename" -> "fitting-" <> ln,
            "Remove Kramers" -> False,
            "Host" -> "LiYF4",
            "Energy Uncertainty in K" -> \[Sigma]exp,
            "Explorer" -> False,
            "SymmetrySimplifier" -> {symGroup, symmetrySimplifier}
            ];
        
        PrintAndDefine[">> Calculating full diagonalization RMS ..."];
        fullDiffs = (eigenEnergies[[;; Length[expData]]] - 
            First /@ expData)[[Complement[presentDataIndices, 
            excludedDataIndices]]];
        fullRMSsquared = Sqrt[Total[fullDiffs^2]/thisSolution["degreesOfFreedom"]];
        thisSolution["bestRMSwithTruncation"] = thisSolution["bestRMS"];
        thisSolution["energiesWithTruncation"] = thisSolution["energies"];
        thisSolution["excludedDataIndices"] = excludedDataIndices;
        thisSolution["bestRMSwithFullDiagonalization"] = fullRMSsquared;
        thisSolution["energiesWithFullDiagonalization"] = eigenEnergies;
        thisSolution = KeyDrop[thisSolution, {"bestRMS", "energies"}];
        
        truncatedSols[numE] = thisSolution;
        LogSol[thisSolution, filePrefix <> "/" <> ln <> "-" <> ToString[truncationEnergy]];
        
        PrintAndDefine["> Updating param model functions ..."];
        bestParams = Association@thisSolution["bestParams"];
        Do[
        (
            If[Not[MemberQ[Keys[bestParams], pVar]],
            Continue[];
            ];
            varModel[pVar]["data"] = 
            DeleteDuplicates@
            Append[varModel[pVar]["data"], {numE, bestParams[pVar]}];
            numData = Length[varModel[pVar]["data"]];
            If[numData > 0,
            (
                order = Which[
                    numData == 1,
                    0,
                    numData == 2,
                    1,
                    numData >= 3,
                    2];
                varModel[pVar]["function"] = Which[
                    MemberQ[linearParams, pVar],
                        NonlinearModelFit[varModel[pVar]["data"], a0 + a1 x, {a0, a1}, x],
                    MemberQ[quadraticParams, pVar],
                        NonlinearModelFit[varModel[pVar]["data"], a0 + a1 x + a2 x^2, {a0, a1, a2}, x],
                    MemberQ[interpolatedParams, pVar],
                        Interpolation[varModel[pVar]["data"], InterpolationOrder -> order]
                    ]
            )
            ];
        ),
        {pVar, pVars}
        ];
        
        Run["~/Scripts/pushover \"finished fitting " <> ln <> "\""];
    ),
    {numE0, doThese}];
    Return[{varModel, truncatedSols}];
    );