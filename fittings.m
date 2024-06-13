(*-----------------------------------------------------------------------+
|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
|~~~+----------------------------------------------------------------+~~~|
|~~~|                  _____ __  __  _                               |~~~|
|~~~|                 / __(_) /_/ /_(_)___  ____ ______              |~~~|
|~~~|                / /_/ / __/ __/ / __ \/ __ `/ ___/              |~~~|
|~~~|               / __/ / /_/ /_/ / / / / /_/ (__  )               |~~~|
|~~~|              /_/ /_/\__/\__/_/_/ /_/\__, /____/                |~~~|
|~~~|                                    /____/                      |~~~|
|~~~|                                                                |~~~|
|~~~+----------------------------------------------------------------+~~~|
|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
|~~~~+--------------------------------------------------------------+~~~~|
|~~~~|                                                              |~~~~|
|~~~~|                                                              |~~~~|
|~~~~|  This script puts together some code useful for fitting the  |~~~~|
|~~~~|                  model Hamiltonian to data.                  |~~~~|
|~~~~|                                                              |~~~~|
|~~~~|                                                              |~~~~|
|~~~~+--------------------------------------------------------------+~~~~|
|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
+-----------------------------------------------------------------------*)

Get["qlanth.m"]
Get["qonstants.m"];
Get["misc.m"];
LoadCarnall[];
LoadFreeIon[];

Jiggle::usage = "Jiggle[num, wiggleRoom] takes a number and randomizes it a little by adding or subtracting a random fraction of itself. The fraction is controlled by wiggleRoom.";
Jiggle[num_, wiggleRoom_ : 0.1] := RandomReal[{1 - wiggleRoom, 1 + wiggleRoom}] * num;

AddToList::usage = "AddToList[list, element, maxSize, addOnlyNew] prepends the element to list and returns the list. If maxSize is reached, the last element is dropped. If addOnlyNew is True (the default), the element is only added if it is different from the last element.";
AddToList[list_, element_, maxSize_, addOnlyNew_ : True] := Module[{
    tempList = If[
      addOnlyNew,
      If[
       Length[list] == 0,
       {element},
       If[
        element != list[[-1]],
        Append[list, element],
        list
        ]
      ],
      Append[list, element]
      ]},
   If[Length[tempList] > maxSize,
    Drop[tempList, Length[tempList] - maxSize],
    tempList]
   ];

ProgressNotebook::usage="ProgressNotebook[] creates a progress notebook for the solver. This notebook includes a plot of the RMS history and the current parameter values. The notebook is returned. The RMS history and the parameter values are updated by setting the variables rmsHistory and paramSols. The variables stringPartialVars and paramSols are used to display the parameter values in the notebook. The notebook is created with the title \"Solver Progress\". The notebook is created with the option WindowSelected->True. The notebook is created with the option TextAlignment->Center. The notebook is created with the option WindowTitle->\"Solver Progress\".";
Options[ProgressNotebook] = {"Basic" -> True};
ProgressNotebook[OptionsPattern[]] := (
   nb = Which[
   OptionValue["Basic"],
   CreateDocument[(
      {
       Dynamic[
        TextCell[
         If[
          Length[paramSols] > 0,
          TableForm[
            Prepend[
              Transpose[{stringPartialVars,
              paramSols[[-1]]}],
              {"RMS",rmsHistory[[-1]]}]
              ],
          ""
          ],
         "Output"
         ],
        TrackedSymbols :> {paramSols, stringPartialVars}
        ]
       }
      ),
     WindowSize -> {600, 1000},
     WindowSelected -> True,
     TextAlignment -> Center,
     WindowTitle -> "Solver Progress"
     ],
   True,
   CreateDocument[(
      {
        "",
        Dynamic[Framed[progressMessage]],
        Dynamic[
        GraphicsColumn[
         {ListPlot[rmsHistory,
           PlotMarkers -> "OpenMarkers",
           Frame -> True,
           FrameLabel -> {"Iteration", "RMS"},
           ImageSize -> 800,
           AspectRatio -> 1/3,
           FrameStyle -> Directive[Thick, 15],
           PlotLabel -> If[Length[rmsHistory] != 0, rmsHistory[[-1]], ""]
           ],
          ListPlot[(#/#[[1]]) & /@ Transpose[paramSols],
           Joined -> True,
           PlotRange -> {All, {-5, 5}},
           Frame -> True,
           ImageSize -> 800,
           AspectRatio -> 1,
           FrameStyle -> Directive[Thick, 15],
           FrameLabel -> {"Iteration", "Params"}
           ]
          }
         ],
        TrackedSymbols :> {rmsHistory, paramSols}],
       Dynamic[
        TextCell[
         If[
          Length[paramSols] > 0,
          TableForm[Transpose[{stringPartialVars, paramSols[[-1]]}]],
          ""
          ],
         "Output"
         ],
        TrackedSymbols :> {paramSols, stringPartialVars}
        ]
       }
      ),
     WindowSize -> {600, 1000},
     WindowSelected -> True,
     TextAlignment -> Center,
     WindowTitle -> "Solver Progress"
     ]
   ];
   Return[nb];
   );

energyCostFunTemplate::usage="energyCostFunTemplate is template used to define the cost function for the energy matching. The template is used to define a function TheRightEnergyPath that takes a list of variables and returns the RMS of the energy differences between the computed and the experimental energies. The template requires the values to the following keys to be provided: `vars` and `varPatterns`";
energyCostFunTemplate = StringTemplate["
TheRightEnergyPath[`varPatterns`]:= (
	{eigenEnergies, eigenVecs} = Eigensystem[compHam[`vars`]];
	ordering      = Ordering[eigenEnergies];
	eigenEnergies = eigenEnergies - Min[eigenEnergies];
	states        = Transpose[{eigenEnergies, eigenVecs}];
	states        = states[[ordering]];
	coarseStates  = ParseStates[states, basis];
  coarseStates  = {#[[1]],#[[-1]]}& /@ coarseStates;
  (* The eigenvectors need to be simplified in order to compare labels to labels *)
	missingLevels = Length[coarseStates]-Length[expData];
  (* The energies are in the first element of the tuples. *)
	energyDiffFun = (Abs[#1[[1]]-#2[[1]]])&;
	(* match disregarding labels *)
	energyFlow    = FlowMatching[coarseStates,
	                expData,
	                \"notMatched\" -> missingLevels,
	                \"CostFun\"    -> energyDiffFun
	                ];
	energyPairs   = {#[[1]][[1]], #[[2]][[1]]} & /@ energyFlow[[1]];
	energyRms     = Sqrt[Total[(Abs[#[[2]]-#[[1]]])^2 & /@ energyPairs] / Length[energyPairs]];
	Return[energyRms];
)"];

AppendToLog[message_, file_String] := Module[
  {timestamp = DateString["ISODateTime"], msgString},
  (
    msgString = ToString[message, InputForm]; (* Convert any expression to a string *)
    OpenAppend[file];
    WriteString[file, timestamp, " - ", msgString, "\n"];
    Close[file];
  )
  ];

energyAndLabelCostFunTemplate::usage="energyAndLabelCostFunTemplate is a template used to define the cost function that includes both the differences between energies and the differences between labels. The template is used to define a function TheRightSignedPath that takes a list of variables and returns the RMS of the energy differences between the computed and the experimental energies together with a term that depends on the differences between the labels. The template requires the values to the following keys to be provided: `vars` and `varPatterns`";
energyAndLabelCostFunTemplate = StringTemplate["
TheRightSignedPath[`varPatterns`] := Module[
  {energyRms, eigenEnergies, eigenVecs, ordering, states, coarseStates, missingLevels, energyDiffFun, energyFlow, energyPairs, energyAndLabelFun, energyAndLabelFlow, totalAvgCost},
  (
    {eigenEnergies, eigenVecs} = Eigensystem[compHam[`vars`]];
    ordering      = Ordering[eigenEnergies];
    eigenEnergies = eigenEnergies - Min[eigenEnergies];
    states        = Transpose[{eigenEnergies, eigenVecs}];
    states        = states[[ordering]]; 
    coarseStates  = ParseStates[states, basis];
    
    (* The eigenvectors need to be simplified in order to compare labels to labels *)
    coarseStates  = {#[[1]],#[[-1]]}& /@ coarseStates;
    missingLevels = Length[coarseStates]-Length[expData];
    
    (* The energies are in the first element of the tuples. *)
    energyDiffFun = ( Abs[#1[[1]]-#2[[1]]] ) &;
    
    (* matching disregarding labels to get overall scale for scaling differences in labels *)
    energyFlow    = FlowMatching[coarseStates,
                    expData,
                    \"notMatched\" -> missingLevels,
                    \"CostFun\"    -> energyDiffFun
                    ];
    energyPairs   = {#[[1]][[1]], #[[2]][[1]]}&/@energyFlow[[1]];
    energyRms     = Sqrt[Total[(Abs[#[[2]]-#[[1]]])^2 & /@ energyPairs]/Length[energyPairs]];
    
    (* matching using both labels and energies  *)
    energyAndLabelFun  = With[{del=energyRms},
                (Abs[#1[[1]]-#2[[1]]] + 
                If[#1[[2]]==#2[[2]],
                  0.,
                  del])&];

    (* energyAndLabelFun  = With[{del=energyRms},
                (Abs[#1[[1]]-#2[[1]]] + 
                del*EditDistance[#1[[2]],#2[[2]]])&]; *)
    energyAndLabelFun  = ( Abs[#1[[1]] - #2[[1]]] + EditDistance[#1[[2]],#2[[2]]] )&;
    energyAndLabelFlow = FlowMatching[coarseStates,
                        expData,
                        \"notMatched\" -> missingLevels,
                        \"CostFun\"    -> energyAndLabelFun
                        ];
    totalAvgCost       = Total[energyAndLabelFun@@# & /@ energyAndLabelFlow[[1]]]/Length[energyAndLabelFlow[[1]]];
    Return[totalAvgCost];
  )
]"];

truncatedEnergyCostTemplate = StringTemplate["
TheTruncatedAndSignedPath[`varsWithNumericQ`] :=
(
  (* Calculate the truncated Hamiltonian *)
  numericalFreeIonHam = compileIntermediateTruncatedHam[`varsMixedWithFixedVals`];
  
  (* Diagonalize it *)
  {truncatedEigenvalues, truncatedEigenVectors} = Eigensystem[numericalFreeIonHam];
  
  (* Using the truncated eigenvectors push them up to the full state space *)
  pulledTruncatedEigenVectors = truncatedEigenVectors.Transpose[truncatedIntermediateBasis];
  states = Transpose[{truncatedEigenvalues, pulledTruncatedEigenVectors}];
  states = SortBy[states, First];
  states = ShiftedLevels[states];

  (* Coarsen the resulting eigenstates *)
  coarseStates = ParseStates[states, basis];

  (* Grab the parts that are needed for fitting *)
  coarseStates = {#[[1]], #[[-1]]} & /@ coarseStates;
  
  (* This cost function takes into account both labels and energies a random factor is added for the sake of stability of the solver*)
  energyAndLabelFun = (
        (
          Abs[#1[[1]] - #2[[1]]] + 
          EditDistance[#1[[2]], #2[[2]]]
        ) * 
        (1 + RandomReal[{-10^-6, 10^-6}])) &;
  
  (* This one only takes into account the energies *)
  energyFun = (Abs[#1[[1]] - #2[[1]]]*(1 + RandomReal[{0, 10^-6}]) &);
  
  (* Choose which cost function to use *)
  costFun = energyAndLabelFun;
  
  (* Not all states are to be matched to the experimental data *)
  missingLevels = Length[coarseStates] - Length[expData];
  
  (* If there are more experimental data than calculated ones, don't leave any state unmatched to those*)
  missingLevels = If[missingLevels < 0, 0, missingLevels];
  
  (* Apply the Hungarian algorithm to match the two sets of data *)
  energyAndLabelFlow = FlowMatching[coarseStates,
    expData,
    \"notMatched\" -> missingLevels,
    \"CostFun\" -> costFun];
  totalCosts   = (costFun @@ #)& /@ energyAndLabelFlow[[1]];
  totalAvgCost = Total[totalCosts] / Length[energyAndLabelFlow[[1]]];
  Return[totalAvgCost]
)"
];

Constrainer::usage = "Constrainer[problemVars, ln] returns a list of constraints for the variables in problemVars for trivalent lanthanide ion ln. problemVars are standard model symbols (F2, F4, ...). The ranges returned are based in the fitted parameters for LaF3 as found in Carnall et al. They could probably be more fine grained, but these ranges are seen to describe all the ions in that case.";
Constrainer[problemVars_, ln_] := (
  slater = Which[
    MemberQ[{"Ce", "Yb"}, ln],
    {},
    True,
    {#, (20000. < # < 120000.)} & /@ {F2, F4, F6}
    ];
  alpha = Which[
    MemberQ[{"Ce", "Yb"}, ln],
    {},
    True,
    {{\[Alpha], 14. < \[Alpha] < 22.}}
    ];
  zeta = {{\[Zeta], 500. < \[Zeta] < 3200.}};
  beta = Which[
    MemberQ[{"Ce", "Yb"}, ln],
    {},
    True,
    {{\[Beta], -1000. < \[Beta] < -400.}}
    ];
  gamma = Which[
    MemberQ[{"Ce", "Yb"}, ln],
    {},
    True,
    {{\[Gamma], 1000. < \[Gamma] < 2000.}}
    ];
  tees = Which[
    ln == "Tm",
    {100. < T2 < 500.},
    MemberQ[{"Ce", "Pr", "Yb"}, ln],
    {},
    True,
    {#, -500. < # < 500.} & /@ {T2, T3, T4, T6, T7, T8}];
  marvins = Which[
    MemberQ[{"Ce", "Yb"}, ln],
    {},
    True,
    {{M0, 1.0 < M0 < 5.0}}
    ];
  peas = Which[
    MemberQ[{"Ce", "Yb"}, ln],
    {},
    True,
    {{P2, -200. < P2 < 1200.}}
    ];
  crystalRanges = {#, (-2000. < # < 2000.)} & /@ (Intersection[
      cfSymbols, problemVars]);
  allCons = 
   Join[slater, zeta, alpha, beta, gamma, tees, marvins, peas, 
    crystalRanges];
  allCons = Select[allCons, MemberQ[problemVars, #[[1]]] &];
  Return[Flatten[Rest /@ allCons]]
  )

LogSol::usage = "LogSol[expr, solHistory, prefix] saves the given expression to a file. The file is named with the given prefix and a created UUID. The file is saved in the \"log\" directory under the current directory. The file is saved in the format of a .m file. The function returns the name of the file.";
LogSol[theSolution_, prefix_] := (
   fname = prefix <> "-sols-" <> CreateUUID[] <> ".m";
   fname = FileNameJoin[{".", "log", fname}];
   Print["Saving solution to: ", fname];
   Export[fname, theSolution];
   Return[fname];
   );
  

FitToHam::usage = "FitToHam[numE, expData, fitToSymbols, simplifier, OptionsPattern[]] fits the model Hamiltonian to the experimental data for the trivalent lanthanide ion with number numE. The experimental data is given in the form of a list of tuples. The first element of the tuple is the energy and the second element is the label. The function saves the results to a file, with the string filePrefix prepended to it, by default this is an empty string, in which case the filePrefix is modified to be the name of the lanthanide.
The fitToSymbols is a list of the symbols to be fit. The simplifier is a list of rules that simplify the Hamiltonian. 
The options and their defaults are:
\"PrintFun\"->PrintTemporary,
\"FilePrefix\"->\"\",
\"SlackChannel\"->None,
\"MaxHistory\"->100,
\"MaxIters\"->100,
\"NumCycles\"->10,
\"ProgressWindow\"->True
The PrintFun option is the function used to print progress messages.
The FilePrefix option is the prefix to use for the file name, by default this is the symbol for the lanthanide.
The SlackChannel option is the channel to post progress messages to.
The MaxHistory option is the maximum number of iterations to keep in the history.
The MaxIters option is the maximum number of iterations for the solver.
The NumCycles option is the number of cycles to run the solver for.
The function returns a list of solutions. The solutions are the results of the NMinimize function. The solutions are a list of tuples. The first element of the tuple is the RMS error and the second element is the parameter values
The function also saves the solutions to a file. The file is named with a prefix and a UUID. The file is saved in the current directory. The file is saved in the format of a .m file.";
Options[FitToHam] = {
      "PrintFun"-> PrintTemporary,
      "FilePrefix" -> "",
      "SlackChannel" -> None,
      "MaxHistory" -> 100,
      "ProgressWindow"-> True,
      "MaxIters" -> 100,
      "NumCycles" -> 10};
FitToHam[numE_Integer, expData_List, fitToSymbols_List, simplifier_List, OptionsPattern[]] :=
  (
    PrintFun        = OptionValue["PrintFun"];
    fitToVars       = ToExpression[ToString[#] <> "v"] & /@ fitToSymbols;
    stringfitToVars = ToString /@ fitToVars;
    slackChan       = OptionValue["SlackChannel"];
    maxHistory      = OptionValue["MaxHistory"];
    maxIters        = OptionValue["MaxIters"];
    numCycles       = OptionValue["NumCycles"];
    ln              = theLanthanides[[numE]];
    logFilePrefix   = If[OptionValue["FilePrefix"] == "", 
                        ToString[theLanthanides[[numE]]],
                        OptionValue["FilePrefix"]];
    PrintFun["Assembling the Hamiltonian for f^", numE, " ..."];
    ham = HamMatrixAssembly[numE];
    PrintFun["Simplifying the symbolic expression for the Hamiltonian in terms of the given simplifier..."];
    ham = ReplaceInSparseArray[ham, simplifier];
    PrintFun["Determining the variables to be fit for ..."];
    (* as they remain after simplifying *)
    fitVars = Variables[Normal[ham]]; 
    (* append v to symbols *)
    varVars = ToExpression[ToString[#] <> "v"] & /@ fitVars;
  
    PrintFun[
    "Compiling a function for efficient evaluation of the Hamiltonian matrix ..."];
    compHam = Compile[Evaluate[fitVars], Evaluate[N[Normal[ham]]]];

    PrintFun[
    "Defining the cost function according to given energies and state labels ..."];

    varPatterns = StringJoin[{ToString[#], "_?NumericQ"}] & /@ fitVars;
    varPatterns = Riffle[varPatterns, ", "];
    varPatterns = StringJoin[varPatterns];
    vars = ToString[#] & /@ fitVars;
    vars = Riffle[vars, ", "];
    vars = StringJoin[vars];

    basis = BasisLSJMJ[numE];

    (* define the cost functions given the problem variables *)
    energyCostFunString = 
    energyCostFunTemplate[<|
      "varPatterns" -> varPatterns, 
      "vars" -> vars|>];
    ToExpression[energyCostFunString];
    energyAndLabelCostFunString = energyAndLabelCostFunTemplate[<|"varPatterns" -> varPatterns, "vars" -> vars|>];
    ToExpression[energyAndLabelCostFunString];

    PrintFun["getting starting values from LaF3..."];
    lnParams = LoadLaF3Parameters[ln];
    bills = Table[lnParams[ToExpression[StringTake[ToString[varvar], {1, -2}]]], {varvar, varVars}];

    (* define the function arguments with the frozen args in place *)
    activeArgs = Table[
      If[MemberQ[fitToVars, varvar], 
      varvar, 
      lnParams[ToExpression[StringTake[ToString[varvar], {1, -2}]]]], 
    {varvar, varVars}
    ];
    activeArgs = StringJoin[Riffle[ToString /@ activeArgs, ", "]];
    (* the constraints, very important *)
    constraints = N[Constrainer[fitToVars, ln]];
    complementaryArgs = Table[
      If[MemberQ[fitToVars, varvar],
      varvar, 
      lnParams[ToExpression[StringTake[ToString[varvar], {1, -2}]]]],
    {varvar, varVars}
    ];
    
    fromBill = {B02v -> B02, B04v -> B04, B06v -> B06, B22v -> B22, 
    B24v -> B24, B26v -> B26, B44v -> B44, B46v -> B46, B66v -> B66, 
    M0v -> M0, P2v -> P2} /. lnParams;
    
    If[Not[ValueQ[noteboo]] && OptionValue["ProgressWindow"],
    noteboo = ProgressNotebook["Basic" -> False];
    ];

    threadHeaderTemplate = StringTemplate[ 
    "(`idx`/`reps`) Fitting data for `ln` using `freeVars`."
    ];
    solutions = {};
    Do[
    (
      (* Remove the downvalues of the cost function *)
      (* DownValues[TheRightSignedPath] = {DownValues[TheRightSignedPath][[-1]]}; *)
      (* start history anew *)
      rmsHistory = {};
      paramSols  = {};
      startTime  = Now;
      threadMessage = threadHeaderTemplate[
          <|"reps" -> numCycles,
          "idx" -> rep, 
          "ln" -> ln,
          "freeVars" -> ToString[fitToVars]|>];
      If[slackChan =!= None,
        threadTS = PostMessageToSlack[threadMessage, slackChan]["ts"]
      ];
      solverTemplateNMini = StringTemplate["
        numIter = 0;
        sol = NMinimize[
          Evaluate[
            Join[{TheRightSignedPath[`activeArgs`]},
              constraints
            ]
          ],
          fitToVars,
          MaxIterations->`maxIterations`,
          Method->`Method`,
          `Monitor`:>(
            currentErr = TheRightSignedPath[`activeArgs`];
            numIter   += 1;
            rmsHistory = AddToList[rmsHistory, currentErr, maxHistory, False];
            paramSols  = AddToList[paramSols, fitToVars, maxHistory, False];
          )
        ]
      "];
      solverCode = solverTemplateNMini[<|
        "maxIterations" -> maxIters,
        "Method" -> "{\"DifferentialEvolution\", 
            \"PostProcess\"   -> False,
            \"ScalingFactor\" -> 0.9,               
            \"RandomSeed\"    -> RandomInteger[{0,1000000}],
            \"SearchPoints\"  -> 10}",
        "Monitor" -> "StepMonitor",
        "activeArgs" -> activeArgs|>];
      ToExpression[solverCode];
      timeTaken = QuantityMagnitude[Now - startTime, "Seconds"];
      Print["Took " <> ToString[timeTaken] <> "s"];
      Print[sol];
      {bestError, bestParams} = sol;
      resultMessage = "sigma=" <> ToString[Round[bestError, 0.1]];
      logFname = LogSol[sol, logFilePrefix];
      If[slackChan =!= None,
      (
        PostMessageToSlack[resultMessage,   slackChan, "threadTS" -> threadTS];
        PostFileToSlack[logFname, logFname, slackChan, "threadTS" -> threadTS];
      )
      ];
     
      vsBill = TableForm[
      Transpose[{
        First /@ fromBill,
        Last /@ fromBill, 
        Round[Last /@ bestParams, 1.]}], 
      TableHeadings -> {None, {"Param", "Bill Bkq", "ql Bkq"}}
      ];
      If[slackChan =!= None,
        PostPdfToSlack[logFname, vsBill, slackChan, "threadTS" -> threadTS]
      ];

      (* analysis code *)

      finalHam = compHam @@ (complementaryArgs /. bestParams);
      {eigenEnergies, eigenVecs} = Eigensystem[finalHam];
      ordering = Ordering[eigenEnergies];
      eigenEnergies = eigenEnergies - Min[eigenEnergies];
      states = Transpose[{eigenEnergies, eigenVecs}];
      states = states[[ordering]];
      coarseStates = ParseStates[states, basis];
      
      (* The eigenvectors need to be simplified in order to compare labels to labels *)
      coarseStates = {#[[1]], #[[-1]]} & /@ coarseStates;
      missingLevels = Length[coarseStates] - Length[expData];
      (* The energies are in the first element of the tuples. *)
      energyDiffFun = (Abs[#1[[1]] - #2[[1]]]) &;
      (* matching disregarding labels to get overall scale for scaling differences in labels *)
      energyFlow = FlowMatching[coarseStates,
              expData, 
              "notMatched" -> missingLevels,
              "CostFun" -> energyDiffFun];
      energyPairs = {#[[1]][[1]], #[[2]][[1]]} & /@ energyFlow[[1]];
      energyRms = Sqrt[Total[(Abs[#[[2]] - #[[1]]])^2 & /@ energyPairs] / Length[energyPairs]];
      (* matching using both labels and energies *)
      energyAndLabelFun = (Abs[#1[[1]] - #2[[1]]] + EditDistance[#1[[2]],#2[[2]]]) &;
      energyAndLabelFlow = FlowMatching[coarseStates,
            expData, 
            "notMatched" -> (Length[coarseStates] - Length[expData]), 
            "CostFun" -> energyAndLabelFun];
      totalAvgCost = Total[energyAndLabelFun @@ # & /@ energyAndLabelFlow[[1]]] / Length[energyAndLabelFlow[[1]]];

      compa = (Flatten /@ energyAndLabelFlow[[1]]);
      compa = Join[
          #,
          {
            #[[2]] == #[[4]], 
            If[NumberQ[#[[1]]],
              Round[#[[1]] - #[[3]],1], 
              ""
            ],
          #[[5]] - #[[3]], 
          Which[
            Round[Abs[#[[1]] - #[[3]]]] < Round[Abs[#[[5]] - #[[3]]]],
              "Better",
            Round[Abs[#[[1]] - #[[3]]]] == Round[Abs[#[[5]] - #[[3]]]],
              "Equal",
            True,
              "Worse"
          ]
          }
      ] & /@ compa;
      atable = TableForm[compa, 
        TableHeadings -> {None, 
          {"ql", "ql", "Bill (exp)", "Bill (exp)",
          "Bill (calc)", "labels=", "ql - exp", "bill - exp"}}
      ];
      atable = Framed[atable, FrameMargins -> 20];
      upsAndDowns = {
        {"Better", Length[Select[compa, #[[-1]] == "Better" &]]},
        {"Equal",  Length[Select[compa, #[[-1]] == "Equal" &]]},
        {"Worse",  Length[Select[compa, #[[-1]] == "Worse" &]]}
      };
      upsAndDowns = TableForm[upsAndDowns];
      If[slackChan =!= None,
        PostPdfToSlack["table", atable, slackChan, "threadTS" -> threadTS];
      ];
      solutions = Append[solutions, sol];
    ),
    {rep, 1, numCycles}
    ];
  )

TruncationFit::usage="TruncationFit[numE, expData, numReps, activeVars, startingValues, Options] fits the given expData in an f^numE configuration, generating numReps different solutions, and varying the symbols in activeVars. The list startingValues is a list with all of the parameters needed to define the Hamiltonian (including values for activeVars, which will be disregarded but are required as position placeholders). The function returns a list of solutions. The solutions are the results of the NMinimize function using the Differential Evolution method. The solutions are a list of tuples. The first element of the tuple is the RMS error and the second element is a list of replacement rules for the fitted parameters. Once each NMinimize is done, the function saves the solutions to a file. The file is named with a prefix and a UUID. The file is saved in the log sub-directory as a .m file. The solver is always constrained by the relevant subsets of constraints for the parameters as provided by the Constrainer function. By default the Differential Evolution method starts with a generation of points within the given constraints, however it is also possible here to have a different region from which the initial points are chosen with the option \"StartingForVars\". 

The following options can be used:
  \"SignatureCheck\" : if True then then the function ends prematurely, printing a list with the symbols that would have defined the Hamiltonian after all simplifications have been applied. Useful to check the entire parameter set that the Hamiltonian has, which has to match one-to-one what is provided by startingValues.
  \"FilePrefix\" : the prefix to use for the file name, by default this is the symbol for the lanthanide.
  \"AccuracyGoal\": sets the accuracy goal for NMinimize, the default is 3.
  \"MaxHistory\" : determines how long the logs for the solver can be.
  \"MaxIterations\": determines the maximum number of iterations used by NMinimize.
  
  \"AccuracyGoal\": the accuracy goal used by NMinimize, default of 3.
  \"TrucationEnergy\": if Automatic then the maximum energy in expData is taken, else it takes the value set by this option. In all cases the energies in expData are truncated to this value.
  \"PrintFun\": the function used to print progress messages, the default is PrintTemporary.
  \"SlackChannel\": name of the Slack channel to which to dump progress messaages, the default is None which disables this option entirely.
  \"ProgressView\": whether or not a progress window will be opened to show the progress of the solver, the default is True.
  
  \"ReturnHashFileNameAndExit\": if True then the function returns the name of the file with the solutions and exits, the default is False.
  \"StartingForVars\": if different from {} then it has to be a list with two elements. The first element being a number that determines the fraction half-width of the interval used for choosing the initial generation of points. The second element being a list with as many elements as activeVars corresponding to the midpoints from which the intial generation points are chosen. The default is {}. 
  \"DE:CrossProbability\": the cross probability used by the Differential Evolution method, the default is 0.5.
  \"DE:ScalingFactor\": the scaling factor used by the Differential Evolution method, the default is 0.6.
  \"DE:SearchPoints\": the number of search points used by the Differential Evolution method, the default is Automatic.
  
  \"MagneticSimplifier\": a list of replacement rules to simplify the Marvin and pesudo-magnetic paramters.
  \"MagFieldSimplifier\": a list of replacement rules to specify a magnetic field (in T), if set to {}, then {Bx, By, Bz} can also then be used as variables to be fitted for.
  \"SymmetrySimplifier\": a list of replacements rules to simplify the crystal field.
  \"OtherSimplifier\": an additiona list of replacement rules that are applied to the Hamiltonian before computing with it.
  \"ThreeBodySimplifier\": the default is an Association that simply states which three body parameters Tk are zero in different configurations, if a list of replacement rules is used then that is used instead for the given problem.
  
  \"FreeIonSymbols\": a list with the symbols to be included in the intermediate coupling basis.
  \"AppendToLogFile\": an association appended to the log file under the key \"Appendix\".
  ";
Options[TruncationFit]={
  "MaxHistory"       -> 200,
  "MaxIterations"    -> 100,
  "FilePrefix"       -> "",
  "AccuracyGoal"     -> 3,
  "TruncationEnergy" -> Automatic,
  "PrintFun"         -> PrintTemporary,
  "SlackChannel"     -> None,
  "ProgressView"     -> True,
  "SignatureCheck"   -> False,
  "AppendToLogFile"  -> <||>,
  "StartingForVars"  -> {},
  "ReturnHashFileNameAndExit" -> False,
  "DE:CrossProbability" -> 0.5,
  "DE:ScalingFactor"    -> 0.6,
  "DE:SearchPoints"     -> Automatic,
  "MagneticSimplifier"  -> {
    M2 -> 56/100 M0,
    M4 -> 31/100 M0,
    P4 -> 1/2 P2,
    P6 -> 1/10 P2},
  "MagFieldSimplifier" -> {
    Bx->0,By->0,Bz->0
    },
  "SymmetrySimplifier" -> {
    B12->0,B14->0,B16->0,B34->0,B36->0,B56->0,
    S12->0,S14->0,S16->0,S22->0,S24->0,S26->0,S34->0,S36->0,
    S44->0,S46->0,S56->0,S66->0
    },
"OtherSimplifier" -> {
    F0->0,
    P0->0,
    \[Sigma]SS->0,
    T11p->0,T11->0,T12->0,T14->0,T15->0,
    T16->0,T18->0,T17->0,T19->0,T2p->0
    },
"ThreeBodySimplifier" -> <|
    1  -> {
      T2->0,T3->0,T4->0,T6->0,T7->0,T8->0,T11p->0,T11->0,T12->0,T14->0,T15->0,T16->0,T18->0,T17->0,T19->0,T2p->0},2->{T2->0,T3->0,T4->0,T6->0,T7->0,T8->0,T11p->0,T11->0,T12->0,T14->0,T15->0,T16->0,T18->0,T17->0,T19->0,T2p->0
      },
    3  -> {},
    4  -> {},
    5  -> {},
    6  -> {},
    7  -> {},
    8  -> {},
    9  -> {},
    10 -> {},
    11 -> {},
    12 -> {
      T3->0,T4->0,T6->0,T7->0,T8->0,T11p->0,T11->0,T12->0,T14->0,T15->0,T16->0,T18->0,T17->0,T19->0,T2p->0
      },
    13->{
      T2->0,T3->0,T4->0,T6->0,T7->0,T8->0,T11p->0,T11->0,T12->0,T14->0,T15->0,T16->0,T18->0,T17->0,T19->0,T2p->0
      }
    |>,
"FreeIonSymbols" -> {F0, F2, F4, F6, M0, P2, \[Alpha], \[Beta], \[Gamma], \[Zeta], T2, T3, T4, T6, T7, T8}
};
TruncationFit[numE_Integer, expData0_List, numReps_Integer, activeVars_List, startingValues_List, OptionsPattern[]]:=(
  ln               = theLanthanides[[numE]];
  expData          = expData0;
  PrintFun         = OptionValue["PrintFun"];
  truncationEnergy = If[OptionValue["TruncationEnergy"]===Automatic,
    Max[First/@expData],
    OptionValue["TruncationEnergy"]
  ];
  oddsAndEnds   = <||>;
  expData       = Select[expData, #[[1]] <= truncationEnergy &];
  maxIterations = OptionValue["MaxIterations"];
  maxHistory    = OptionValue["MaxHistory"];
  slackChan     = OptionValue["SlackChannel"];
  accuracyGoal  = OptionValue["AccuracyGoal"];
  logFilePrefix = If[OptionValue["FilePrefix"] == "", 
                      ToString[theLanthanides[[numE]]],
                      OptionValue["FilePrefix"]];
  
  usingInitialRange  = Not[OptionValue["StartingForVars"] === {}];
  If[usingInitialRange,
    (
      PrintFun["Using the solver for initial values in range ..."];
      {fractionalWidth, startVarValues} = OptionValue["StartingForVars"];
    )
  ];
  
  magneticSimplifier  = OptionValue["MagneticSimplifier"];
  magFieldSimplifier  = OptionValue["MagFieldSimplifier"];
  symmetrySimplifier  = OptionValue["SymmetrySimplifier"];
  otherSimplifier     = OptionValue["OtherSimplifier"];
  threeBodySimplifier = If[Head[OptionValue["ThreeBodySimplifier"]] == Association,
    OptionValue["ThreeBodySimplifier"][numE],
    OptionValue["ThreeBodySimplifier"]
  ]; 
  simplifier          = Join[magneticSimplifier,
                        magFieldSimplifier,
                        symmetrySimplifier,
                        threeBodySimplifier,
                        otherSimplifier];
  freeIonSymbols = OptionValue["FreeIonSymbols"];
  runningInteractive = (Head[$ParentLink] === LinkObject);
  
  oddsAndEnds["simplifier"] = simplifier;
  oddsAndEnds["freeIonSymbols"] = freeIonSymbols;
  oddsAndEnds["truncationEnergy"] = truncationEnergy;
  oddsAndEnds["numE"] = numE;
  oddsAndEnds["expData"] = expData;
  oddsAndEnds["numReps"] = numReps;
  oddsAndEnds["activeVars"] = activeVars;
  oddsAndEnds["startingValues"] = startingValues;
  oddsAndEnds["maxIterations"] = maxIterations;
  oddsAndEnds["PrintFun"] = PrintFun;
  oddsAndEnds["ln"] = ln;
  oddsAndEnds["numE"] = numE;
  oddsAndEnds["accuracyGoal"] = accuracyGoal;
  oddsAndEnds["Appendix"] = OptionValue["AppendToLogFile"];
  
  hamDim        = Binomial[14, numE];
  modelSymbols  = Sort[Select[paramSymbols, Not[MemberQ[Join[racahSymbols, chenSymbols,{t2Switch, \[Epsilon],gs}],#]]&]];
  (* Remove the symbols that will be removed by the simplifier, no symbol should remain here that is not in the symbolic hamiltonian *)
  reducedModelSymbols = Select[modelSymbols, Not[MemberQ[Keys[simplifier],#]]&];
  If[OptionValue["SignatureCheck"],
    (
    Print["Given the model parameters and the simplifying assumptions, the resultant model parameters are:"];
    Print[{reducedModelSymbols}];
    Print["The ordering in these needs to be respected in the startValues parameter ..."];
    Print["Exiting ..."];
    Return[""];
    )
  ];
  
  (*calculate the basis*)
  basis        = BasisLSJMJ[numE];
  (* grab the Hamiltonian preserving its block structure *)
  PrintFun["Assembling the Hamiltonian for f^",numE," keeping the block structure ..."];
  ham          = HamMatrixAssembly[numE, "ReturnInBlocks"->True];
  (* apply the simplifier *)
  PrintFun["Simplifying using the given aggregate set of simplification rules ..."];
  ham          = Map[ReplaceInSparseArray[#, simplifier]&, ham, {2}];
  
  (* Get the reference parameters from LaF3 *)
  PrintFun["Getting reference parameters for ",ln," using LaF3 ..."];
  lnParams = LoadLaF3Parameters[ln];
  freeBies = Prepend[Values[(#->(#/.lnParams))&/@freeIonSymbols],numE];
  (* a more explicit alias *)
  allVars           = reducedModelSymbols;
  
  oddsAndEnds["allVars"]  = allVars;
  oddsAndEnds["freeBies"] = freeBies;
  
  (* reload compiled version if found *)
  varHash                  = Hash[{numE, allVars, freeBies, truncationEnergy}];
  compileIntermediateFname = "compileIntermediateTruncatedHam-"<>ToString[varHash]<>".mx";
  truncadedFname           = "TheTruncatedAndSignedPath-"<>ToString[varHash]<>".mx";
  If[OptionValue["ReturnHashFileNameAndExit"],
    (
      Print[varHash];
      Return[truncadedFname];
    )
  ];
  If[FileExistsQ[compileIntermediateFname],
    PrintFun["This ion and free-ion params have been compiled before (as determined by {numE, allVars, freeBies, truncationEnergy}). Loading the previously saved function and intermediate coupling basis ..."];
    {compileIntermediateTruncatedHam, truncatedIntermediateBasis} = Import[compileIntermediateFname];,
  (
    PrintFun["Zeroing out every symbol in the Hamiltonian that is not a free-ion parameter ..."];
    (* Get the free ion symbols *)
    freeIonSimplifier = (#->0) & /@ Complement[reducedModelSymbols,freeIonSymbols];
    (* Take the diagonal blocks for the intermediate analysis *)
    PrintFun["Grabbing the diagonal blocks of the Hamiltonian ..."];
    diagonalBlocks       = Diagonal[ham];
    (* simplify them to only keep the free ion symbols *)
    PrintFun["Simplifying the diagonal blocks to only keep the free ion symbols ..."];
    diagonalScalarBlocks = ReplaceInSparseArray[#,freeIonSimplifier]&/@diagonalBlocks;
    (* these include the MJ quantum numbers, remove that *)
    PrintFun["Contracting the basis vectors by removing the MJ quantum numbers from the diagonal blocks ..."];
    diagonalScalarBlocks = FreeHam[diagonalScalarBlocks, numE];

    argsOfTheIntermediateEigensystems         = StringJoin[Riffle[Prepend[(ToString[#]<>"v_") & /@ freeIonSymbols,"numE_"],", "]];
    argsForEvalInsideOfTheIntermediateSystems = StringJoin[Riffle[(ToString[#]<>"v") & /@ freeIonSymbols,", "]];
    PrintFun["argsOfTheIntermediateEigensystems = ",argsOfTheIntermediateEigensystems];
    PrintFun["argsForEvalInsideOfTheIntermediateSystems = ",argsForEvalInsideOfTheIntermediateSystems];
    PrintFun["If the following fails, make sure to modify the arguments of TheIntermediateEigensystems to match the ones above ..."];

    (* Compile a function that will effectively calculate the spectrum of all of the scalar blocks given the parameters of the free-ion part of the Hamiltonian *)
    (* Compile one function for each of the blocks *)
    PrintFun["Compiling functions for the diagonal blocks of the Hamiltonian ..."];
    compiledDiagonal = Compile[Evaluate[freeIonSymbols], Evaluate[N[Normal[#]]]]&/@diagonalScalarBlocks;
    (* Use that to create a function that will calculate the free-ion eigensystem *)
    TheIntermediateEigensystems[numEv_, F0v_, F2v_, F4v_, F6v_, M0v_, P2v_, \[Alpha]v_, \[Beta]v_, \[Gamma]v_, \[Zeta]v_, T2v_, T3v_, T4v_, T6v_, T7v_, T8v_] := (
      theNumericBlocks = (#[F0v, F2v, F4v, F6v, M0v, P2v, \[Alpha]v, \[Beta]v, \[Gamma]v, \[Zeta]v, T2v, T3v, T4v, T6v, T7v, T8v]&)/@compiledDiagonal;
      theIntermediateEigensystems = Eigensystem/@theNumericBlocks;
      Js     = AllowedJ[numEv];
      basisJ = BasisLSJMJ[numEv,"AsAssociation"->True];
      (* Having calculated the eigensystems with the removed degeneracies, put the degeneracies back in explicitly *)
      elevatedIntermediateEigensystems = MapIndexed[EigenLever[#1,2Js[[#2[[1]]]]+1]&, theIntermediateEigensystems];
      pivot = If[EvenQ[numEv],0,-1/2];
      LSJmultiplets = (#[[1]]<>ToString[InputForm[#[[2]]]])&/@Select[BasisLSJMJ[numEv],#[[-1]]== pivot &];
      (* Calculate the multiplet assignments that the intermediate basis eigenvectors have *)
      multipletAssignments = Table[
        (
          J         = Js[[idx]];
          eigenVecs = theIntermediateEigensystems[[idx]][[2]];
          majorComponentIndices        = Ordering[Abs[#]][[-1]]&/@eigenVecs;
          majorComponentAssignments    = LSJmultiplets[[#]]&/@majorComponentIndices;
          (* All of the degenerate eigenvectors belong to the same multiplet*)
          elevatedMultipletAssignments = ListRepeater[majorComponentAssignments,2J+1];
          elevatedMultipletAssignments
        ),
      {idx, 1, Length[Js]}
      ];
      (* Put together the multiplet assignments and the energies *)
      freeIenergiesAndMultiplets = Transpose/@Transpose[{First/@elevatedIntermediateEigensystems, multipletAssignments}];
      freeIenergiesAndMultiplets = Flatten[freeIenergiesAndMultiplets, 1];
      (* Calculate the change of basis matrix using the intermediate coupling eigenvectors *)
      basisChanger = BlockDiagonalMatrix[Transpose/@Last/@elevatedIntermediateEigensystems];
      basisChanger = SparseArray[basisChanger];
      Return[{theIntermediateEigensystems, multipletAssignments, elevatedIntermediateEigensystems, freeIenergiesAndMultiplets, basisChanger}]
    );
    
    PrintFun["Calculating the intermediate eigensystems for ",ln," using free-ion params from LaF3 ..."];
    (* Calculate intermediate coupling basis using the free-ion params for LaF3 *)
    {theIntermediateEigensystems, multipletAssignments, elevatedIntermediateEigensystems, freeIenergiesAndMultiplets, basisChanger} = TheIntermediateEigensystems@@freeBies;
    
    (* Use that intermediate coupling basis to compile a function for the full Hamiltonian *)
    allFreeEnergies = Flatten[First/@elevatedIntermediateEigensystems];
    (* Important that the intermediate coupling basis have attached energies, which make possible the truncation *)
    ordering = Ordering[allFreeEnergies];
    (* Sort the free ion energies and determine which indices should be included in the truncation *)
    allFreeEnergiesSorted          = Sort[allFreeEnergies];
    {minFreeEnergy, maxFreeEnergy} = MinMax[allFreeEnergies];
    (* Determine the index at which the energy is equal or larger than the truncation energy *)
    sortedTruncationIndex = Which[
      truncationEnergy > (maxFreeEnergy-minFreeEnergy),
      hamDim,
      True,
      FirstPosition[allFreeEnergiesSorted-Min[allFreeEnergiesSorted],x_/;x>truncationEnergy,{0},1][[1]]
    ];
    (* The actual energy at which the truncation is made *)
    roundedTruncationEnergy = allFreeEnergiesSorted[[sortedTruncationIndex]];

    (* The indices that enact the truncation *)
    truncationIndices = ordering[[;;sortedTruncationIndex]];

    (* Using the ham (with all the symbols) change the basis to the computed one *)
    PrintFun["Changing the basis of the Hamiltonian to the intermediate coupling basis ..."];
    intermediateHam            = Transpose[basisChanger].ArrayFlatten[ham].basisChanger;
    (* Using the truncation indices truncate that one *)
    PrintFun["Truncating the Hamiltonian ..."];
    truncatedIntermediateHam   = intermediateHam[[truncationIndices,truncationIndices]];
    (* These are the basis vectors for the truncated hamiltonian *)
    PrintFun["Saving the truncated intermediate basis ..."];
    truncatedIntermediateBasis = basisChanger[[All,truncationIndices]];

    PrintFun["Compiling a function for the truncated Hamiltonian ..."];
    (* Compile a function that will calculate the truncated Hamiltonian given the parameters in allVars, this is the function to be use in fitting *)
    compileIntermediateTruncatedHam = Compile[Evaluate[allVars], 
                                      Evaluate[N[Normal[truncatedIntermediateHam]]]];
    (* Save the compiled function *)
    PrintFun["Saving the compiled function for the truncated Hamiltonian and the truncatedIntermediateBasis..."];
    Export[compileIntermediateFname, {compileIntermediateTruncatedHam, truncatedIntermediateBasis}];
  )
  ];
  
  TheTruncatedAndSignedPathGenerator::usage = "This function puts together the necessary expression for defining a function which has as arguments all the symbolic values in varsMixedWithVals and which feeds to compileIntermediateTruncatedHam the arguments as given in varsMixedWithVals. varsMixedWithVals needs to respect the order of aruments expected by compileIntermediateTruncatedHam. Once the necessary template has been used this function then results in the definition of the function TheTruncatedAndSignedPath.";
  TheTruncatedAndSignedPathGenerator[varsMixedWithVals_List]:=(
    variableVars  = Select[varsMixedWithVals,Not[NumericQ[#]]&];
    numQSignature = StringJoin[Riffle[(ToString[#]<>"_?NumericQ")&/@variableVars,", "]];
    varWithValsSignature = StringJoin[Riffle[(ToString[#]<>"")&/@varsMixedWithVals,", "]];
    funcString    = truncatedEnergyCostTemplate[<|"varsWithNumericQ"->numQSignature,"varsMixedWithFixedVals" -> varWithValsSignature|>];
    ClearAll[TheTruncatedAndSignedPath];
    ToExpression[funcString]
  );

  (* We need to create a function call that has all the frozen parameters in place and all the active symbols unevaluated *)
  (* find the indices of the activeVars to create the function signature *)
  activeVarIndices = Flatten[Position[allVars,#]&/@activeVars];
  (* we start from the numerical values in the current best*)
  jobVars          = startingValues;
  (* we then put back the symbols that should be unevaluated *)
  jobVars[[activeVarIndices]] = activeVars;
  
  oddsAndEnds["jobVars"] = jobVars;
  (* calculate the constraints *)
  constraints                = N[Constrainer[activeVars,ln]];
  oddsAndEnds["constraints"] = constraints;
  (* This is useful for the progress window *)
  activeVarsString  = StringJoin[Riffle[ToString/@activeVars,", "]];
  TheTruncatedAndSignedPathGenerator[jobVars];
  stringPartialVars = ToString/@activeVars;
  
  activeVarsWithRange = If[usingInitialRange,
    MapIndexed[Flatten[{#1, 
        (1-Sign[startVarValues[[#2]]]*fractionalWidth) * startVarValues[[#2]],
        (1+Sign[startVarValues[[#2]]]*fractionalWidth) * startVarValues[[#2]]
        }]&, activeVars],
    activeVars
  ];
  
  (* this is the template for the minimizer *)
  solverTemplateNMini = StringTemplate["
  numIter = 0;
  sol = NMinimize[
    Evaluate[
      Join[{TheTruncatedAndSignedPath[`activeVarsString`]},
        constraints
      ]
    ],
    activeVarsWithRange,
    AccuracyGoal   -> `accuracyGoal`,
    MaxIterations  -> `maxIterations`,
    Method->`Method`,
    `Monitor`:>(
      currentErr = TheTruncatedAndSignedPath[`activeVarsString`];
      currentParams = activeVars;
      numIter   += 1;
      rmsHistory = AddToList[rmsHistory, currentErr, maxHistory, False];
      paramSols  = AddToList[paramSols, activeVars, maxHistory, False];
      If[Not[runningInteractive],(
        Print[numIter,\"/\",`maxIterations`];
        Print[\"err    = \", ToString[NumberForm[Round[currentErr,0.001],{Infinity,3}]]];
        Print[\"params = \", ToString[NumberForm[Round[#,0.0001],{Infinity,4}]] &/@ currentParams];
        )
      ];
    )
  ]
  "];
  methodStringTemplate = StringTemplate["
      {\"DifferentialEvolution\", 
        \"PostProcess\"   -> False,
        \"ScalingFactor\" -> `DE:ScalingFactor`,
        \"CrossProbability\" -> `DE:CrossProbability`,
        \"RandomSeed\"    -> RandomInteger[{0,1000000}],
        \"SearchPoints\"  -> `DE:SearchPoints`}"];
  methodString = methodStringTemplate[<|
    "DE:ScalingFactor" -> OptionValue["DE:ScalingFactor"],
    "DE:CrossProbability" -> OptionValue["DE:CrossProbability"],
    "DE:SearchPoints" -> OptionValue["DE:SearchPoints"]|>];
  (* Evaluate the template *)
  solverCode = solverTemplateNMini[<|
    "accuracyGoal"  -> accuracyGoal,
    "maxIterations" -> maxIterations,
    "Method"->"{\"DifferentialEvolution\", 
                \"PostProcess\"   -> False,
                \"ScalingFactor\" -> 0.6,
                \"CrossProbability\" -> 0.25,               
                \"RandomSeed\"    -> RandomInteger[{0,1000000}],
                \"SearchPoints\"  -> Automatic}",
    "Monitor"->"StepMonitor",
    "activeVarsString"->activeVarsString|>
  ];
  threadHeaderTemplate = StringTemplate[ "(`idx`/`reps`) Fitting data for `ln` using `freeVars`."];
  (* Find as many solutions as numReps *)
  sols = Table[(
    rmsHistory        = {};
    paramSols         = {};
    openNotebooks     = If[runningInteractive,
              ("WindowTitle"/.NotebookInformation[#]) & /@ Notebooks[],
              {}];
    If[Not[MemberQ[openNotebooks,"Solver Progress"]] && OptionValue["ProgressView"],
      ProgressNotebook["Basic"->False]
    ];
    If[Not[slackChan === None],
    (
      threadMessage = threadHeaderTemplate[<|"reps" -> numReps, "idx" -> rep, "ln" -> ln, 
        "freeVars" -> ToString[activeVars]|>];
      threadTS = PostMessageToSlack[threadMessage, slackChan]["ts"];
    )
    ];
    startTime = Now;
    ToExpression[solverCode];
    
    timeTaken = QuantityMagnitude[Now - startTime, "Seconds"];
    Print["Took " <> ToString[timeTaken] <> "s"];
    Print[sol];
    bestError     = sol[[1]];
    bestParams    = sol[[2]];
    resultMessage = "sigma=" <> ToString[Round[bestError, 0.1]];
    solAssoc = <|
      "bestRMS"      -> bestError,
      "solHistory"   -> rmsHistory,
      "bestParams"   -> bestParams,
      "paramHistory" -> paramSols,
      "timeTaken/s"    -> timeTaken
      |>;
    solAssoc = Join[solAssoc,   oddsAndEnds];
    logFname = LogSol[solAssoc, logFilePrefix];
    
    If[Not[slackChan===None],(
      PostMessageToSlack[resultMessage,   slackChan, "threadTS" -> threadTS];
      PostFileToSlack[StringSplit[logFname,"/"][[-1]], logFname, slackChan, "threadTS" -> threadTS]
    )
    ];
    solAssoc
    ),
    {rep,1,numReps}
  ];
  Return[sols];
);

ClassicalFit::usage="ClassicalFit[numE, expData, excludeDataIndices, problemVars, startValues, \[Sigma]exp, constraints_List, Options] fits the given expData in an f^numE configuration, by using the symbols in problemVars. The symbols given in problemVars may be constrained or held constant, this being controlled by constraints list which is a list of replacement rules expressing desired constraints. The constraints list additional constraints imposed upon the model parameters that remain once other simplifications have been \"baked\" into the compiled Hamiltonians that are used to increase the speed of the calculation. 
 
Important, note that in the case of odd number of electrons the given data must explicitly include the Kramers degeneracy; excludeDataIndices must be compatible with this.
 
The list expData needs to be a list of lists with the only restriction that the first element of them corresponds to energies of levels. In this list, an empty value can be used to indicate known gaps in the data. Even if the energy value for a level is known (and given in expData) certain values can be omitted from the fitting procedure through the list excludeDataIndices, which correspond to indices in expData that should be skipped over.
 
The Hamiltonian used for fitting is version that has been truncated either by using the maximum energy given in expData or by manually setting a truncation energy using the option \"TruncationEnergy\".
 
The argument \[Sigma]exp is the estimated uncertainty in the differences between the calculated and the experimental energy levels. This is used to estimate the uncertainty in the fitted parameters. Admittedly this will be a rough estimate (at least on the contribution of the calculated uncertainty), but it is better than nothing and may at least provide a lower bound to the uncertainty in the fitted parameters. It is assumed that the uncertainty in the differences between the calculated and the experimental energy levels is the same for all of them.
 
The list startValues is a list with all of the parameters needed to define the Hamiltonian (including the initial values for problemVars). 
 
The function saves the solution to a file. The file is named with a prefix (controlled by the option \"FilePrefix\") and a UUID. The file is saved in the log sub-directory as a .m file.
 
Here's a description of the different parts of this function: first the Hamiltonian is assembled and simplified using the given simplifications. Then the intermediate coupling basis is calculated using the free-ion parameters for the given lanthanide. The Hamiltonian is then changed to the intermediate coupling basis and truncated. The truncated Hamiltonian is then compiled into a function that can be used to calculate the energy levels of the truncated Hamiltonian. The function that calculates the energy levels is then used to fit the experimental data. The fitting is done using FindMinimum with the Levenberg-Marquardt method.
 
The function returns an association with the following keys: 
 
- \"bestRMS\" which is the best \[Sigma] value found.
- \"bestParams\" which is the best set of parameters found for the variables that were not constrained.
- \"bestParamsWithConstraints\" which has the best set of parameters (from - \"bestParams\") together with the used constraints. These include all the parameters in the model, even those that were not fitted for.
- \"paramSols\" which is a list of the parameters trajectories during the stepping of the fitting algorithm.
- \"timeTaken/s\" which is the time taken to find the best fit.
- \"simplifier\" which is the simplifier used to simplify the Hamiltonian.
- \"excludeDataIndices\" as given in the input.
- \"startValues\" as given in the input.
 
- \"freeIonSymbols\" which are the symbols used in the intermediate coupling basis.
- \"truncationEnergy\" which is the energy used to truncate the Hamiltonian, if it was set to Automatic, the value here is the actual energy used.
- \"numE\" which is the number of electrons in the f^numE configuration.
- \"expData\" which is the experimental data used for fitting.
- \"problemVars\" which are the symbols considered for fitting
 
- \"maxIterations\" which is the maximum number of iterations used by NMinimize.
- \"hamDim\" which is the dimension of the full Hamiltonian.
- \"allVars\" which are all the symbols defining the Hamiltonian under the aggregate simplifications.
- \"freeBies\" which are the free-ion parameters used to define the intermediate coupling basis.
- \"truncatedDim\" which is the dimension of the truncated Hamiltonian.
- \"compiledIntermediateFname\" the file name of the compiled function used for the truncated Hamiltonian.
 
- \"fittedLevels\" which is the number of levels fitted for.
- \"actualSteps\" the number of steps that FindMiniminum actually took.
- \"solWithUncertainty\" which is a list of replacement rules of the form (paramSymbol -> {bestEstimate, uncertainty}).
- \"rmsHistory\" which is a list of the \[Sigma] values found during the fitting.
- \"Appendix\" which is an association appended to the log file under the key \"Appendix\".
- \"presentDataIndices\" which is the list of indices in expData that were used for fitting, this takes into account both the empty indices in expData and also the indices in excludeDataIndices.
 
- \"states\" which contains a list of eigenvalues and eigenvectors for the fitted model, this is only available if the option \"SaveEigenvectors\" is set to True; if a general shift of energy was allowed for in the fitting, then the energies are shifted accordingly.
- \"energies\" which is a list of the energies of the fitted levels, this is only available if the option \"SaveEigenvectors\" is set to False. If a general shift of energy was allowed for in the fitting, then the energies are shifted accordingly.
- \"degreesOfFreedom\" which is equal to the nunber of fitted state energies minus the number of parameters used in fitting.

The function admits the following options with corresponding default values:
- \"MaxHistory\" : determines how long the logs for the solver can be.
- \"MaxIterations\": determines the maximum number of iterations used by NMinimize.
- \"FilePrefix\" : the prefix to use for the subfolder in the log filder, in which the solution files are saved, by default this is \"calcs\" so that the calculation files are saved under the directory \"log/calcs\".
- \"AddConstantShift\" : if True then a constant shift is allowed in the fitting, default is False. If this is the case the variable \"\[Epsilon]\" is added to the list of variables to be fitted for, it must not be included in problemVars.
  
- \"AccuracyGoal\": the accuracy goal used by NMinimize, default of 5.
- \"TrucationEnergy\": if Automatic then the maximum energy in expData is taken, else it takes the value set by this option. In all cases the energies in expData are only considered up to this value.
- \"PrintFun\": the function used to print progress messages, the default is PrintTemporary.
- \"RefParamsVintage\": the vintage of the reference parameters to use. The reference parameters are both used to determine the truncated Hamiltonian, and also as starting values for the solver. It may be \"LaF3\", in which case reference parameters from Carnall are used. It may also be \"LiYF4\", in which case the reference parameters from the LiYF4 paper are used. It may also be Automatic, in which case the given experimental data is used to determine starting values for F^k and \[Zeta]. It may also be a list or association that provides values for the Slater integrals and spin-orbit coupling, the remaining necessary parameters complemented by using \"LaF3\".
  
- \"SlackChannel\": name of the Slack channel to which to dump progress messaages, the default is None which disables this option.
- \"ProgressView\": whether or not a progress window will be opened to show the progress of the solver, the default is True.
- \"SignatureCheck\": if True then then the function returns prematurely, returning a list with the symbols that would have defined the Hamiltonian after all simplifications have been applied. Useful to check the entire parameter set that the Hamiltonian has, which has to match one-to-one what is provided by startingValues.
- \"SaveEigenvectors\": if True then the both the eigenvectors and eigenvalues are saved under the \"states\" key of the returned association. If False then only the energies are saved, the default is False.
  
- \"AppendToLogFile\": an association appended to the log file under the key \"Appendix\".
- \"MagneticSimplifier\": a list of replacement rules to simplify the Marvin and pesudo-magnetic paramters. Here the ratios of the Marvin parameters and the pseudo-magnetic parameters are defined to simplify the magnetic part of the Hamiltonian.
- \"MagFieldSimplifier\": a list of replacement rules to specify a magnetic field (in T), if set to {}, then {Bx, By, Bz} can also be used as variables to be fitted for.
  
- \"SymmetrySimplifier\": a list of replacements rules to simplify the crystal field.
- \"OtherSimplifier\": an additional list of replacement rules that are applied to the Hamiltonian before computing with it. Here the spin-spin contribution can be turned off by setting \[Sigma]SS->0, which is the default.
";
Options[ClassicalFit] = {
  "MaxHistory"       -> 200,
  "MaxIterations"    -> 100,
  "FilePrefix"       -> "calcs",
  "ProgressView"     -> True,
  "TruncationEnergy" -> Automatic,
  "AccuracyGoal"     -> 5,
  "PrintFun"         -> PrintTemporary,
  "SlackChannel"     -> None,
  "RefParamsVintage" -> "LaF3",
  "ProgressView"     -> True,
  "SignatureCheck"   -> False,
  "AddConstantShift" -> False,
  "SaveEigenvectors" -> False,
  "AppendToLogFile"  -> <||>,
  "Energy Uncertainty in K" -> Automatic,
  "MagneticSimplifier" -> {
    M2 -> 56/100 M0,
    M4 -> 31/100 M0,
    P4 -> 1/2 P2,
    P6 -> 1/10 P2
    },
  "MagFieldSimplifier" -> {
    Bx -> 0,
    By -> 0,
    Bz -> 0
    },
  "SymmetrySimplifier" -> {
    B12->0, B14->0, B16->0, B34->0, B36->0, B56->0, 
    S12->0 ,S14->0, S16->0, S22->0, S24->0, S26->0,
    S34->0 ,S36->0, S44->0, S46->0, S56->0, S66->0
    },
  "OtherSimplifier" -> {
    F0->0,
    P0->0,
    \[Sigma]SS->0,
    T11p->0, T11->0, T12->0, T14->0, T15->0, 
    T16->0, T18->0, T17->0, T19->0, T2p->0,
    wChErrA ->0, wChErrB->0
    },
  "ThreeBodySimplifier" -> <|
    1  -> {
      T2->0, T3->0, T4->0, T6->0, T7->0, T8->0,
      T11p->0, T11->0, T12->0, T14->0, T15->0, T16->0, T18->0, T17->0, T19->0,
      T2p->0}, 
    2 -> {
      T2->0, T3->0, T4->0, T6->0, T7->0, T8->0,
      T11p->0, T11->0, T12->0, T14->0, T15->0, T16->0, T18->0, T17->0, T19->0,
      T2p->0
      },
    3  -> {},
    4  -> {},
    5  -> {},
    6  -> {},
    7  -> {},
    8  -> {},
    9  -> {},
    10 -> {},
    11 -> {},
    12 -> {
      T3->0, T4->0, T6->0, T7->0, T8->0,
      T11p->0, T11->0, T12->0, T14->0, T15->0, T16->0, T18->0, T17->0, T19->0,
      T2p->0
      },
    13->{
      T2->0, T3->0, T4->0, T6->0, T7->0, T8->0,
      T11p->0, T11->0, T12->0, T14->0, T15->0, T16->0, T18->0, T17->0, T19->0,
      T2p->0
      }
    |>,
    "FreeIonSymbols" -> {F0, F2, F4, F6, \[Zeta]}
};
ClassicalFit[numE_Integer, expData_List, excludeDataIndices_List, problemVars_List, startValues_Association, constraints_List, OptionsPattern[]]:=Module[
  {accuracyGoal, activeVarIndices, activeVars, activeVarsString, activeVarsWithRange, allFreeEnergies, allFreeEnergiesSorted, allVars, allVarsVec, argsForEvalInsideOfTheIntermediateSystems, argsOfTheIntermediateEigensystems, aVar, aVarPosition, basis, basisChanger, basisChangerBlocks, bestError, bestParams, bestRMS, blockShifts, blockSizes, colIdx, compiledDiagonal, compiledIntermediateFname, constrainedProblemVars, constrainedProblemVarsList, covMat, currentRMS, degressOfFreedom, dependentVars, diagonalBlocks, diagonalScalarBlocks, diff, eigenEnergies, eigenvalueDispenserTemplate, eigenVectors, elevatedIntermediateEigensystems, endTime, fmSol, fmSolAssoc, fractionalWidth, freeBies, freeIenergiesAndMultiplets, freeionSymbols, fullHam, fullSolVec, funcString, ham, hamDim, hamEigenvaluesTemplate, hamString, hess, indepSolVecVec, indepVars, intermediateHam, isolationValues, jobVars, lin, linMat, ln, lnParams, logFilePrefix, logFname, magneticSimplifier, maxFreeEnergy, maxHistory, maxIterations, methodString, methodStringTemplate, minFreeEnergy, minpoly, modelSymbols, multipletAssignments, needlePosition, numBlocks, numQSignature, numReps, solCompendium, openNotebooks, ordering, othersFixed, otherSimplifier, p0, paramBest, paramSigma, perHam, polySols, presentDataIndices, PrintFun, problemVarsPositions, problemVarsQ, problemVarsQString, problemVarsVec, problemVarsWithStartValues, reducedModelSymbols, resultMessage, roundedTruncationEnergy, rowIdx, runningInteractive, shiftToggle, simplifier, slackChan, sol, solAssoc, sols, solWithUncertainty, sortedTruncationIndex, sqdiff, standardValues, starTime, startingValues, startTime, startVarValues, states, steps, symmetrySimplifier, theIntermediateEigensystems, TheIntermediateEigensystems, TheTruncatedAndSignedPathGenerator, thisPoly, threadHeaderTemplate, threadMessage, threadTS, timeTaken, totalVariance, truncadedFname, truncatedIntermediateBasis, truncatedIntermediateHam, truncationEnergy, truncationIndices, RefParams, truncationUmbral, usingInitialRange, varHash, varIdx, varsWithConstants, varWithValsSignature, \[Lambda]0Vec, \[Lambda]exp},
  (
    \[Sigma]exp   = OptionValue["Energy Uncertainty in K"];
    solCompendium = <||>;
    refParamsVintage = OptionValue["RefParamsVintage"];
    RefParams = Which[
      refParamsVintage === "LaF3",
      LoadLaF3Parameters,
      refParamsVintage === "LiYF4",
      LoadLiYF4Parameters,
      True,
      refParamsVintage 
    ];
    hamDim        = Binomial[14, numE];
    addShift      = OptionValue["AddConstantShift"];
    ln            = theLanthanides[[numE]];
    maxHistory    = OptionValue["MaxHistory"];
    maxIterations = OptionValue["MaxIterations"];
    logFilePrefix = If[OptionValue["FilePrefix"] == "", 
                      ToString[theLanthanides[[numE]]],
                      OptionValue["FilePrefix"]
                    ];
    accuracyGoal = OptionValue["AccuracyGoal"];
    slackChan    = OptionValue["SlackChannel"];
    PrintFun     = OptionValue["PrintFun"];
    freeIonSymbols      = OptionValue["FreeIonSymbols"];
    runningInteractive  = (Head[$ParentLink] === LinkObject);
    magneticSimplifier  = OptionValue["MagneticSimplifier"];
    magFieldSimplifier  = OptionValue["MagFieldSimplifier"];
    symmetrySimplifier  = OptionValue["SymmetrySimplifier"];
    otherSimplifier     = OptionValue["OtherSimplifier"];
    threeBodySimplifier = If[Head[OptionValue["ThreeBodySimplifier"]] == Association,
                            OptionValue["ThreeBodySimplifier"][numE],
                            OptionValue["ThreeBodySimplifier"]
                          ]; 
    
    truncationEnergy = If[OptionValue["TruncationEnergy"]===Automatic,
      PrintFun["Truncation energy set to Automatic, using the maximum energy (+20%) in the data ..."];
      Round[1.2 * Max[Select[First /@ expData, NumericQ[#] &]]],
      OptionValue["TruncationEnergy"]
      ];
    truncationEnergy = Max[50000, truncationEnergy];
    PrintFun["Using a truncation energy of ", truncationEnergy, " K"];
    
    simplifier = Join[magneticSimplifier,
                      magFieldSimplifier,
                      symmetrySimplifier,
                      threeBodySimplifier,
                      otherSimplifier];
    
    PrintFun["Determining gaps in the data ..."];
    (* whatever is non-numeric is assumed as a known gap *)
    presentDataIndices = Flatten[Position[expData, {_?(NumericQ[#] &), ___}]];
    (* some indices omitted here based on the excludeDataIndices argument *)
    presentDataIndices = Complement[presentDataIndices, excludeDataIndices];
    
    solCompendium["simplifier"]         = simplifier;
    solCompendium["excludeDataIndices"] = excludeDataIndices;
    solCompendium["startValues"]        = startValues;
    solCompendium["freeIonSymbols"]     = freeIonSymbols;
    solCompendium["truncationEnergy"]   = truncationEnergy;
    solCompendium["numE"]               = numE;
    solCompendium["expData"]            = expData;
    solCompendium["problemVars"]        = problemVars;
    solCompendium["maxIterations"]      = maxIterations;
    solCompendium["hamDim"]             = hamDim;
    solCompendium["constraints"]        = constraints;
    modelSymbols  = Sort[Select[paramSymbols, Not[MemberQ[Join[racahSymbols, juddOfeltIntensitySymbols, chenSymbols,{t2Switch, \[Epsilon], gs, nE}],#]]&]];
    (* remove the symbols that will be removed by the simplifier, no symbol should remain here that is not in the symbolic Hamiltonian *)
    reducedModelSymbols = Select[modelSymbols, Not[MemberQ[Keys[simplifier],#]]&];
    
    (* this is useful to understand what are the arguments of the truncated compiled Hamiltonian *)
    If[OptionValue["SignatureCheck"],
      (
        Print["Given the model parameters and the simplifying assumptions, the resultant model parameters are:"];
        Print[{reducedModelSymbols}];
        Print["Exiting ..."];
        Return[""];
      )
    ];
    
    (* calculate the basis *)
    PrintFun["Retrieving the LSJMJ basis for f^", numE, " ..."];
    basis    = BasisLSJMJ[numE];
    
    Which[refParamsVintage === Automatic,
      (
        PrintFun["Using the automatic vintage with freshly fitted free-ion parameters and others as in LaF3 ..."];
        lnParams = LoadLaF3Parameters[ln];
        freeIonSol = FreeIonSolver[expData, numE];
        freeIonParams = freeIonSol["bestParams"];
        lnParams = Join[lnParams, freeIonParams];
      ),
      MemberQ[{List, Association}, Head[RefParams]],
      (
        RefParams = Association[RefParams];
        PrintFun["Using the given parameters as a starting point ..."];
        lnParams = RefParams;
        extraParams = LoadLaF3Parameters[ln];
        lnParams = Join[extraParams, lnParams];
      ),
      True,
      (
        (* get the reference parameters from the given vintage *)
        PrintFun["Getting reference free-ion parameters for ", ln, " using ", refParamsVintage, " ..."];
        lnParams = ParamPad[RefParams[ln]];
      )
    ];
    freeBies = Prepend[Values[(#->(#/.lnParams)) &/@ freeIonSymbols],  numE];
    (* a more explicit alias *)
    allVars            = reducedModelSymbols;
    numericConstraints = Association@Select[constraints, NumericQ[#[[2]]] &];
    standardValues     = allVars /. Join[lnParams, numericConstraints];
    solCompendium["allVars"]  = allVars;
    solCompendium["freeBies"] = freeBies;
    
    (* reload compiled version if found *)
    varHash                   = Hash[{numE, allVars, freeBies, truncationEnergy, simplifier}];
    compiledIntermediateFname = ln <> "-compiled-intermediate-truncated-ham-" <> ToString[varHash] <> ".mx";
    compiledIntermediateFname = FileNameJoin[{moduleDir, "compiled", compiledIntermediateFname}];
    solCompendium["compiledIntermediateFname"] = compiledIntermediateFname;
    
    If[FileExistsQ[compiledIntermediateFname],
      PrintFun["This ion, free-ion params, and full set of variables have been used before (as determined by {numE, allVars, freeBies, truncationEnergy,simplifier}). Loading the previously saved compiled function and intermediate coupling basis ..."];
      PrintFun["Using : ", compiledIntermediateFname];
      {compileIntermediateTruncatedHam, truncatedIntermediateBasis} = Import[compiledIntermediateFname];,
    (
      (* grab the Hamiltonian preserving the block structure *)
      PrintFun["Assembling the Hamiltonian for f^",numE," keeping the block structure ..."];
      ham          = HamMatrixAssembly[numE, "ReturnInBlocks"->True];
      (* apply the simplifier *)
      PrintFun["Simplifying using the aggregate set of simplification rules ..."];
      ham          = Map[ReplaceInSparseArray[#, simplifier]&, ham, {2}];
      PrintFun["Zeroing out every symbol in the Hamiltonian that is not a free-ion parameter ..."];
      (* Get the free ion symbols *)
      freeIonSimplifier = (#->0) & /@ Complement[reducedModelSymbols, freeIonSymbols];
      (* Take the diagonal blocks for the intermediate analysis *)
      PrintFun["Grabbing the diagonal blocks of the Hamiltonian ..."];
      diagonalBlocks       = Diagonal[ham];
      (* simplify them to only keep the free ion symbols *)
      PrintFun["Simplifying the diagonal blocks to only keep the free ion symbols ..."];
      diagonalScalarBlocks = ReplaceInSparseArray[#,freeIonSimplifier]&/@diagonalBlocks;
      (* these include the MJ quantum numbers, remove that *)
      PrintFun["Contracting the basis vectors by removing the MJ quantum numbers from the diagonal blocks ..."];
      diagonalScalarBlocks = FreeHam[diagonalScalarBlocks, numE];
      
      argsOfTheIntermediateEigensystems         = StringJoin[Riffle[Prepend[(ToString[#]<>"v_") & /@ freeIonSymbols,"numE_"],", "]];
      argsForEvalInsideOfTheIntermediateSystems = StringJoin[Riffle[(ToString[#]<>"v") & /@ freeIonSymbols,", "]];
      PrintFun["argsOfTheIntermediateEigensystems = ", argsOfTheIntermediateEigensystems];
      PrintFun["argsForEvalInsideOfTheIntermediateSystems = ", argsForEvalInsideOfTheIntermediateSystems];
      PrintFun["(if the following fails, it might help to see if the arguments of TheIntermediateEigensystems match the ones shown above)"];
      
      (* compile a function that will effectively calculate the spectrum of all of the scalar blocks given the parameters of the free-ion part of the Hamiltonian *)
      (* compile one function for each of the blocks *)
      PrintFun["Compiling functions for the diagonal blocks of the Hamiltonian ..."];
      compiledDiagonal = Compile[Evaluate[freeIonSymbols], Evaluate[N[Normal[#]]]]&/@diagonalScalarBlocks;
      (* use that to create a function that will calculate the free-ion eigensystem *)
      TheIntermediateEigensystems[numEv_, F0v_, F2v_, F4v_, F6v_, \[Zeta]v_] := (
        theNumericBlocks = (#[F0v, F2v, F4v, F6v, \[Zeta]v]&) /@ compiledDiagonal;
        theIntermediateEigensystems = Eigensystem /@ theNumericBlocks;
        Js     = AllowedJ[numEv];
        basisJ = BasisLSJMJ[numEv,"AsAssociation"->True];
        (* having calculated the eigensystems with the removed degeneracies, put the degeneracies back in explicitly *) 
        elevatedIntermediateEigensystems = MapIndexed[EigenLever[#1, 2Js[[#2[[1]]]]+1 ]&, theIntermediateEigensystems];
        (* Identify a single MJ to keep *)
        pivot         = If[EvenQ[numEv],0,-1/2];
        LSJmultiplets = (#[[1]]<>ToString[InputForm[#[[2]]]])&/@Select[BasisLSJMJ[numEv],#[[-1]]== pivot &];
        (* calculate the multiplet assignments that the intermediate basis eigenvectors have *)
        needlePosition = 0;
        multipletAssignments = Table[
          (
            J         = Js[[idx]];
            eigenVecs = theIntermediateEigensystems[[idx]][[2]];
            majorComponentIndices        = Ordering[Abs[#]][[-1]]&/@eigenVecs;
            majorComponentIndices       += needlePosition;
            needlePosition              += Length[majorComponentIndices];
            majorComponentAssignments    = LSJmultiplets[[#]]&/@majorComponentIndices;
            (* All of the degenerate eigenvectors belong to the same multiplet*)
            elevatedMultipletAssignments = ListRepeater[majorComponentAssignments, 2J+1];
            elevatedMultipletAssignments
          ),
        {idx, 1, Length[Js]}
        ];
        (* put together the multiplet assignments and the energies *)
        freeIenergiesAndMultiplets = Transpose/@Transpose[{First/@elevatedIntermediateEigensystems, multipletAssignments}];
        freeIenergiesAndMultiplets = Flatten[freeIenergiesAndMultiplets, 1];
        (* calculate the change of basis matrix using the intermediate coupling eigenvectors *)
        basisChanger = BlockDiagonalMatrix[Transpose/@Last/@elevatedIntermediateEigensystems];
        basisChanger = SparseArray[basisChanger];
        Return[{theIntermediateEigensystems, multipletAssignments, elevatedIntermediateEigensystems, freeIenergiesAndMultiplets, basisChanger}]
      );
      
      PrintFun["Calculating the intermediate eigensystems for ",ln," using free-ion params from LaF3 ..."];
      (* calculate intermediate coupling basis using the free-ion params for LaF3 *)
      {theIntermediateEigensystems, multipletAssignments, elevatedIntermediateEigensystems, freeIenergiesAndMultiplets, basisChanger} = TheIntermediateEigensystems@@freeBies;
      
      (* use that intermediate coupling basis to compile a function for the full Hamiltonian *)
      allFreeEnergies = Flatten[First/@elevatedIntermediateEigensystems];
      (* important that the intermediate coupling basis have attached energies, which make possible the truncation *)
      ordering = Ordering[allFreeEnergies];
      (* sort the free ion energies and determine which indices should be included in the truncation *)
      allFreeEnergiesSorted          = Sort[allFreeEnergies];
      {minFreeEnergy, maxFreeEnergy} = MinMax[allFreeEnergies];
      (* determine the index at which the energy is equal or larger than the truncation energy *)
      sortedTruncationIndex = Which[
        truncationEnergy > (maxFreeEnergy-minFreeEnergy),
        hamDim,
        True,
        FirstPosition[allFreeEnergiesSorted - Min[allFreeEnergiesSorted],x_/;x>truncationEnergy,{0},1][[1]]
      ];
      (* the actual energy at which the truncation is made *)
      roundedTruncationEnergy = allFreeEnergiesSorted[[sortedTruncationIndex]];
      
      (* the indices that participate in the truncation *)
      truncationIndices = ordering[[;;sortedTruncationIndex]];
      PrintFun["Computing the block structure of the change of basis array ..."];
      blockSizes         = BlockArrayDimensionsArray[ham];
      basisChangerBlocks = ArrayBlocker[basisChanger, blockSizes];
      blockShifts        = First /@ Diagonal[blockSizes];
      numBlocks          = Length[blockSizes];
      (* using the ham (with all the symbols) change the basis to the computed one *)
      PrintFun["Changing the basis of the Hamiltonian to the intermediate coupling basis ..."];
      (* intermediateHam            = Transpose[basisChanger].ham.basisChanger; *)
      (* Return[{basisChangerBlocks, ham}]; *)
      intermediateHam = BlockMatrixMultiply[ham, basisChangerBlocks];
      PrintFun["Distributing products inside of symbolic matrix elements to keep complexity in check ..."];
      Do[
        intermediateHam[[rowIdx, colIdx]] = MapToSparseArray[intermediateHam[[rowIdx, colIdx]], Distribute /@ # &];,
        {rowIdx, 1, numBlocks},
        {colIdx, 1, numBlocks}
      ];
      intermediateHam = BlockMatrixMultiply[BlockTranspose[basisChangerBlocks],intermediateHam];
      PrintFun["Distributing products inside of symbolic matrix elements to keep complexity in check ..."];
      Do[
        intermediateHam[[rowIdx, colIdx]] = MapToSparseArray[intermediateHam[[rowIdx, colIdx]], Distribute /@ # &];,
        {rowIdx, 1, numBlocks},
        {colIdx, 1, numBlocks}
      ];
      (* using the truncation indices truncate that one *)
      PrintFun["Truncating the Hamiltonian ..."];
      truncatedIntermediateHam = TruncateBlockArray[intermediateHam, truncationIndices, blockShifts];
      (* these are the basis vectors for the truncated hamiltonian *)
      PrintFun["Saving the truncated intermediate basis ..."];
      truncatedIntermediateBasis = basisChanger[[All,truncationIndices]];
      
      PrintFun["Compiling a function for the truncated Hamiltonian ..."];
      (* compile a function that will calculate the truncated Hamiltonian given the parameters in allVars, this is the function to be use in fitting *)
      compileIntermediateTruncatedHam = Compile[Evaluate[allVars], Evaluate[truncatedIntermediateHam]];
      (* save the compiled function *)
      PrintFun["Saving the compiled function for the truncated Hamiltonian and the truncated intermediate basis ..."];
      Export[compiledIntermediateFname, {compileIntermediateTruncatedHam, truncatedIntermediateBasis}];
    )
    ];
    
    truncationUmbral = Dimensions[truncatedIntermediateBasis][[2]];
    PrintFun["The truncated Hamiltonian has a dimension of ", truncationUmbral, "x", truncationUmbral, " ..."];
    presentDataIndices = Select[presentDataIndices, # <= truncationUmbral &];
    solCompendium["presentDataIndices"] = presentDataIndices;
    
    (* the problemVars are the symbols that will be fitted for *)
    
    PrintFun["Starting up the fitting process using the Levenberg-Marquardt method ..."];
    (* using the problemVars I need to create the argument list including _?NumericQ *)
    problemVarsQ       = (ToString[#] <> "_?NumericQ") & /@ problemVars;
    problemVarsQString = StringJoin[Riffle[problemVarsQ, ", "]];
    (* we also need to have the padded arguments with the variables in the right position and the fixed values in the remaining ones *)
    problemVarsPositions = Position[allVars, #][[1, 1]] & /@ problemVars;
    problemVarsString    = StringJoin[Riffle[ToString /@ problemVars, ", "]];
    (* to feed parameters to the Hamiltonian, which includes all parameters, we need to form the set of arguments, with fixed values where needed, and the variables in the right position *)
    varsWithConstants                         = standardValues;
    varsWithConstants[[problemVarsPositions]] = problemVars;
    varsWithConstantsString                   = ToString[varsWithConstants];
    
    (* this following function serves eigenvalues from the Hamiltonian, has memoization so it might grow to use a lot of RAM *)
    Clear[HamSortedEigenvalues];
    hamEigenvaluesTemplate = StringTemplate["
    HamSortedEigenvalues[`problemVarsQ`]:=(
      ham         = compileIntermediateTruncatedHam@@`varsWithConstants`;
      eigenValues = Chop[Sort@Eigenvalues@ham];
      eigenValues = eigenValues - Min[eigenValues];
      HamSortedEigenvalues[`problemVarsString`] = eigenValues;
      Return[eigenValues]
    )"];
    hamString = hamEigenvaluesTemplate[<|
        "problemVarsQ"      -> problemVarsQString,
        "varsWithConstants" -> varsWithConstantsString,
        "problemVarsString" -> problemVarsString
        |>];
    ToExpression[hamString];
    
    (* we also need a function that will pick the i-th eigenvalue, this seems unnecessary but it's needed to form the right functional form expected by the Levenberg-Marquardt method *)
    eigenvalueDispenserTemplate = StringTemplate["
    PartialHamEigenvalues[`problemVarsQ`][i_]:=(
      eigenVals = HamSortedEigenvalues[`problemVarsString`];
      eigenVals[[i]]
    )
    "];
    eigenValueDispenserString = eigenvalueDispenserTemplate[<|
        "problemVarsQ"      -> problemVarsQString,
        "problemVarsString" -> problemVarsString
        |>];
    ToExpression[eigenValueDispenserString];  
    
    PrintFun["Determining the free variables after constraints ..."];
    constrainedProblemVars     = (problemVars /. constraints);
    constrainedProblemVarsList = Variables[constrainedProblemVars];
    If[addShift,
      PrintFun["Adding a constant shift to the fitting parameters ..."];
      constrainedProblemVarsList = Append[constrainedProblemVarsList, \[Epsilon]]
    ];
    
    indepVars = Complement[problemVars, #[[1]] & /@ constraints]; 
    stringPartialVars = ToString/@constrainedProblemVarsList;
    
    paramSols  = {};
    rmsHistory = {};
    steps      = 0;
    problemVarsWithStartValues = KeyValueMap[{#1,#2} &, startValues];
    If[addShift,
      problemVarsWithStartValues = Append[problemVarsWithStartValues, {\[Epsilon] ,0}];
    ];
    openNotebooks     = If[runningInteractive,
              ("WindowTitle"/.NotebookInformation[#]) & /@ Notebooks[],
              {}];
    If[Not[MemberQ[openNotebooks,"Solver Progress"]] && OptionValue["ProgressView"],
      ProgressNotebook["Basic"->False]
    ];
    degressOfFreedom = Length[presentDataIndices] - Length[problemVars] - 1;
    PrintFun["Fitting for ", Length[presentDataIndices], " data points with ", Length[problemVars], " free parameters.", " The effective degrees of freedom are ", degressOfFreedom, " ..."];
    
    PrintFun["Starting the fitting process ..."];
    startTime   = Now;
    shiftToggle = If[addShift, 1, 0];
    sol = FindMinimum[
      Sum[(expData[[j]][[1]] - (PartialHamEigenvalues @@ constrainedProblemVars)[j] - shiftToggle * \[Epsilon])^2,
      {j, presentDataIndices}],
      problemVarsWithStartValues,
      Method        -> "LevenbergMarquardt",
      MaxIterations -> OptionValue["MaxIterations"],
      AccuracyGoal  -> OptionValue["AccuracyGoal"],
      StepMonitor :> (
        steps      += 1;
        currentSqSum = Sum[(expData[[j]][[1]] - (PartialHamEigenvalues @@ constrainedProblemVars)[j] - shiftToggle * \[Epsilon])^2, {j, presentDataIndices}];
        currentRMS = Sqrt[currentSqSum / degressOfFreedom];
        paramSols  = AddToList[paramSols, constrainedProblemVarsList, maxHistory];
        rmsHistory = AddToList[rmsHistory, currentRMS, maxHistory];
        )
      ];
    endTime = Now;
    timeTaken = QuantityMagnitude[endTime - startTime, "Seconds"];
    PrintFun["Solution found in ", timeTaken, "s"];
    
    solVec = constrainedProblemVars /. sol[[-1]];
    indepSolVec = indepVars /. sol[[-1]];
    If[addShift,
      \[Epsilon]Best = \[Epsilon]/. sol[[-1]],
      \[Epsilon]Best = 0
    ];
    fullSolVec = standardValues;
    fullSolVec[[problemVarsPositions]] = solVec;
    PrintFun["Calculating the truncated numerical Hamiltonian corresponding to the solution ..."];
    fullHam = compileIntermediateTruncatedHam @@ fullSolVec;
    PrintFun["Calculating energies and eigenvectors ..."];
    {eigenEnergies, eigenVectors} = Eigensystem[fullHam];
    states = Transpose[{eigenEnergies, eigenVectors}];
    states = SortBy[states, First];
    eigenEnergies = First /@ states;
    PrintFun["Shifting energies to make ground state zero of energy ..."];
    eigenEnergies = eigenEnergies - eigenEnergies[[1]];
    PrintFun["Calculating the linear approximant to each eigenvalue ..."];
    allVarsVec = Transpose[{allVars}];
    p0 = Transpose[{fullSolVec}];
    linMat = {};
    If[addShift,
      tail = -2,
      tail = -1];
    Do[
      (
        aVarPosition = Position[allVars, aVar][[1, 1]];
        isolationValues = ConstantArray[0, Length[allVars]];
        isolationValues[[aVarPosition]] = 1;
        dependentVars = KeyValueMap[{#1, D[#2, aVar]} &, Association[constraints]];
        Do[
          isolationValues[[Position[allVars, dVar[[1]]][[1, 1]]]] = dVar[[2]],
          {dVar, dependentVars}
        ];
        perHam  = compileIntermediateTruncatedHam @@ isolationValues;
        lin     = FirstOrderPerturbation[Last /@ states, perHam];
        linMat  = Append[linMat, lin];
      ),
      {aVar, constrainedProblemVarsList[[;;tail]]}
    ];
    PrintFun["Removing the gradient of the ground state ..."];
    linMat = (# - #[[1]] & /@ linMat);
    PrintFun["Transposing derivative matrices into columns ..."];
    linMat  = Transpose[linMat];
    
    PrintFun["Calculating the eigenvalue vector at solution ..."];
    \[Lambda]0Vec = Transpose[{eigenEnergies[[presentDataIndices]]}];
    PrintFun["Putting together the experimental vector ..."];
    \[Lambda]exp = Transpose[{First /@ expData[[presentDataIndices]]}];
    problemVarsVec = If[addShift,
      Transpose[{constrainedProblemVarsList[[;;-2]]}],
      Transpose[{constrainedProblemVarsList}] 
    ];
    indepSolVecVec = Transpose[{indepSolVec}];
    PrintFun["Calculating the difference between eigenvalues at solution ..."];
    diff = If[linMat=={},
      (\[Lambda]0Vec - \[Lambda]exp)  + \[Epsilon]Best,
      (\[Lambda]0Vec - \[Lambda]exp)  + \[Epsilon]Best + linMat[[presentDataIndices]].(problemVarsVec - indepSolVecVec)
    ];
    PrintFun["Calculating the sum of squares of differences around solution ... "];
    sqdiff        = Expand[(Transpose[diff] . diff)[[1, 1]]];
    PrintFun["Calculating the minimum (which should coincide with sol) ..."];
    minpoly       = sqdiff /. AssociationThread[problemVars -> solVec];
    fmSolAssoc    = Association[sol[[2]]];
    If[\[Sigma]exp == Automatic,
      \[Sigma]exp = Sqrt[minpoly / degressOfFreedom];
    ];
    \[CapitalDelta]\[Chi]2 = Sqrt[degressOfFreedom];
    Amat = (1/\[Sigma]exp^2) * Transpose[linMat[[presentDataIndices]]].linMat[[presentDataIndices]];
    paramIntervals = EllipsoidBoundingBox[Amat, \[CapitalDelta]\[Chi]2];
    PrintFun["Calculating the uncertainty in the parameters ..."];
    solWithUncertainty = Table[
      (
        aVar        = constrainedProblemVarsList[[varIdx]];
        paramBest   = aVar /. fmSolAssoc;
        (aVar -> {paramBest, paramIntervals[[varIdx, 2]]})
      ),
    {varIdx, 1, Length[constrainedProblemVarsList]-shiftToggle}
    ];
    PrintFun["Calculating the covariance matrix ..."];
    hess = If[linmat=={},
      {{Infinity}},
      2 * Transpose[linMat[[presentDataIndices]]] . linMat[[presentDataIndices]]
    ];
    covMat  = If[linmat=={},
      {{0}},
      \[Sigma]exp^2 * Inverse[hess]
    ];
    bestRMS    = Sqrt[minpoly / degressOfFreedom];
    bestParams = sol[[2]];
    bestWithConstraints = Association@Join[constraints, bestParams];
    bestWithConstraints = bestWithConstraints /. bestWithConstraints;
    bestWithConstraints = (# + 0.) & /@ bestWithConstraints;
    
    solCompendium["degreesOfFreedom"]    = degressOfFreedom;
    solCompendium["solWithUncertainty"]  = solWithUncertainty;
    solCompendium["truncatedDim"]        = truncationUmbral;
    solCompendium["fittedLevels"]        = Length[presentDataIndices];
    solCompendium["actualSteps"]         = steps;
    solCompendium["bestRMS"]             = bestRMS;
    solCompendium["problemVars"]         = problemVars;
    solCompendium["paramSols"]           = paramSols;
    solCompendium["rmsHistory"]          = rmsHistory;
    solCompendium["Appendix"]            = OptionValue["AppendToLogFile"];
    solCompendium["timeTaken/s"]         = timeTaken;
    solCompendium["bestParams"]          = bestParams;
    solCompendium["bestParamsWithConstraints"] = bestWithConstraints;
    
    If[OptionValue["SaveEigenvectors"],
        solCompendium["states"] = {#[[1]] + \[Epsilon]Best, #[[2]]} &/@ (Chop /@ ShiftedLevels[states]),
        (
          finalEnergies = Sort[First /@ states];
          finalEnergies = finalEnergies - finalEnergies[[1]];
          finalEnergies = finalEnergies + \[Epsilon]Best;
          finalEnergies  = Chop /@ finalEnergies;
          solCompendium["energies"] = finalEnergies;
      )
    ];
    logFname = LogSol[solCompendium, logFilePrefix];
    PrintFun["Finished ..."];
    Return[solCompendium]; 
  )
];


caseConstraints::usage="This Association contains the constraints that are not the same across all of the lanthanides. For instance, since the ratio between M2 and M0 is assumed the same for all the trivalent lanthanides, that one is not included here. 
This association has keys equal to symbols of lanthanides and values equal to lists of rules that express either a parameter being held fixed or made proportional to another.
In Table I of Carnall 1989 these correspond to cases were values are given in square brackets.";
caseConstraints = <|
"Ce" -> {
    B02 -> -218.,
    B04 -> 738.,
    B06 -> 679.,
    B22 -> -50.,
    B24 -> 431.,
    B26 -> -921.,
    B44 -> 616.,
    B46 -> -348.,
    B66 -> -788.
    },
"Pr" -> {},
"Nd" -> {},
"Pm" -> {},
"Sm" -> {
    B22 -> -50.,
    T2 -> 300.,
    T3 -> 36.,
    T4 -> 56.,
    \[Gamma] -> 1500.
    },
"Eu" -> {
    F4 -> 0.713 F2,
    F6 -> 0.512 F2,
    B22 -> -50.,
    B24 -> 597.,
    B26 -> -706.,
    B44 -> 408.,
    B46 -> -508.,
    B66 -> -692.,
    M0 -> 2.1,
    P2 -> 360.,
    T2 -> 300.,
    T3 -> 40.,
    T4 -> 60.,
    T6 -> -300.,
    T7 -> 370.,
    T8 -> 320.,
    \[Alpha] -> 20.16,
    \[Beta]  -> -566.9,
    \[Gamma] -> 1500.
    },
"Pm" -> {
    B02 -> -245.,
    B04 -> 470.,
    B06 -> 640.,
    B22 -> -50.,
    B24 -> 525.,
    B26 -> -750.,
    B44 -> 490.,
    B46 -> -450.,
    B66 -> -760.,
    F2 -> 76400.,
    F4 -> 54900.,
    F6 -> 37700.,
    M0 -> 2.4,
    P2 -> 275.,
    T2 -> 300.,
    T3 -> 35.,
    T4 -> 58.,
    T6 -> -310.,
    T7 -> 350.,
    T8 -> 320.,
    \[Alpha] -> 20.5, 
    \[Beta] -> -560.,
    \[Gamma] -> 1475.,
    \[Zeta] -> 1025.},
"Gd" -> {
    F4 -> 0.710 F2,
    B02 -> -231.,
    B04 -> 604.,
    B06 -> 280.,
    B22 -> -99.,
    B24 -> 340.,
    B26 -> -721.,
    B44 -> 452.,
    B46 -> -204.,
    B66 -> -509.,
    T2 -> 300.,
    T3 -> 42.,
    T4 -> 62.,
    T6 -> -295.,
    T7 -> 350.,
    T8 -> 310.,
    \[Beta] -> -600.,
    \[Gamma] -> 1575.
    },
"Tb" -> {
    F4 -> 0.707 F2,
    T2 -> 320.,
    T3 -> 40.,
    T4 -> 50.,
    \[Gamma] -> 1650.
    },
"Dy" -> {},
"Ho" -> {
    B02 -> -240.,
    T2 -> 400.,
    \[Gamma] -> 1800.
    },
"Er" -> {
    T2 -> 400.,
    \[Gamma] -> 1800.
    },
"Tm" -> {
    T2 -> 400.,
    \[Gamma] -> 1820.
    },
"Yb" -> {
    B02 -> -249.,
    B04 -> 457.,
    B06 -> 282.,
    B22 -> -105.,
    B24 -> 320.,
    B26 -> -482.,
    B44 -> 428.,
    B46 -> -234.,
    B66 -> -492.
}
|>;

variedSymbols =<|
    "Ce" -> {\[Zeta]}, 
    "Pr" -> {B02, B04, B06, B22, B24, B26, B44, B46, B66,
            F2, F4, F6, 
            M0, P2, 
            \[Alpha], \[Beta], \[Gamma],
            \[Zeta]}, 
    "Nd" -> {B02, B04, B06, B22, B24, B26, B44, B46, B66, 
            F2, F4, F6, 
            M0, P2,
            T2, T3, T4, T6, T7, T8,
            \[Alpha], \[Beta], \[Gamma],
            \[Zeta]},
    "Pm" -> {}, 
    "Sm" -> {B02, B04, B06, B24, B26, B44, B46, B66,
            F2, F4, F6, M0, P2, 
            T6, T7, T8,
            \[Alpha], \[Beta], \[Zeta]}, 
    "Eu" -> {B02, B04, B06,
            F2, F4, F6, \[Zeta]}, 
    "Gd" -> {F2, F4, F6,
            M0, P2,
            \[Alpha], \[Zeta]}, 
    "Tb" -> {B02, B04, B06, B22, B24, B26, B44, B46, B66,
            F2, F4, F6,
            M0, P2,
            T6, T7, T8, 
            \[Alpha], \[Beta], \[Zeta]}, 
    "Dy" -> {B02, B04, B06, B22, B24, B26, B44, B46, B66,
            F2, F4, F6,
            M0, P2,
            T2, T3, T4, T6, T7, T8, 
            \[Alpha], \[Beta], \[Gamma], \[Zeta]},
    "Ho" -> {B04, B06, B22, B24, B26, B44, B46, B66,
            F2, F4, F6,
            M0, P2,
            T3, T4, T6, T7, T8,
            \[Alpha], \[Beta], \[Zeta]}, 
    "Er" -> {B02, B04, B06, B22, B24, B26, B44, B46, B66,
            F2, F4, F6,
            M0, P2,
            T3, T4, T6, T7, T8, \[Alpha], \[Beta], \[Zeta]}, 
    "Tm" -> {B02, B04, B06, B22, B24, B26, B44, B46, B66, 
            F2, F4, F6, 
            M0, P2, 
            \[Alpha], \[Beta], \[Zeta]}, 
    "Yb" -> {\[Zeta]}
    |>;
  
caseConstraintsLiYF4 = <|
  "Ce" -> {
    B04 -> -1043.,
    B44 -> -1249.,
    B06 -> -65.,
    B46 -> -1069.
    }, 
  "Pr" -> {
    \[Beta] -> -644.,
    \[Gamma] -> 1413., 
    M0 -> 1.88,
    P2 -> 244.
    },
  "Nd" -> {
    M0 -> 1.85
    }, 
  "Sm" -> {
    \[Alpha] -> 20.5,
    \[Beta] -> -616.,
    \[Gamma] -> 1565., 
    T2 -> 282.,
    T3 -> 26.,
    T4 -> 71.,
    T6 -> -257.,
    T7 -> 314., 
    T8 -> 328.,
    M0 -> 2.38,
    P2 -> 336.
    }, 
  "Eu" -> {
    T2 -> 370.,
    T3 -> 40.,
    T4 -> 40.,
    T6 -> -300.,
    T7 -> 380., 
    T8 -> 370.
    }, 
  "Tb" -> {
    F4 -> 0.709 F2, 
    F6 -> 0.503 F2,
    \[Alpha] -> 17.6,
    \[Beta] -> -581.,
    \[Gamma] -> 1792.,
    T2 -> 330.,
    T3 -> 40.,
    T4 -> 45.,
    T6 -> -365.,
    T7 -> 320., 
    T8 -> 349.,
    M0 -> 2.7,
    P2 -> 482.
    }, 
  "Dy" -> {
    F4 -> 0.707 F2, 
    F6 -> 0.516 F2,
    \[Alpha] -> 17.9,
    \[Beta] -> -628.,
    \[Gamma] -> 1790.,
    T2 -> 326.,
    T3 -> 23.,
    T4 -> 83.,
    T6 -> -294.,
    T7 -> 403., 
    T8 -> 340.,
    M0 -> 4.46,
    P2 -> 610.,
    B46 -> -700.
    },
  "Ho" -> {
    \[Alpha] -> 17.2,
    \[Beta] -> -596.,
    \[Gamma] -> 1839., 
    T2 -> 365.,
    T3 -> 37.,
    T4 -> 95.,
    T6 -> -274.,
    T7 -> 331., 
    T8 -> 343., 
    P2 -> 582.
    },
  "Er" -> {},
  "Tm" -> {
    \[Alpha] -> 17.3,
    \[Beta] -> -665.,
    \[Gamma] -> 1936.,
    M0 -> 4.93, 
    P2 -> 730.
    },
  "Yb" -> {
    B06 -> -23.,
    B46 -> -512.
    }
  |>;
  
variedSymbolsLiYF4 = <|
  "Ce" -> {
    B02, \[Zeta]
    }, 
  "Pr" -> {
    B02, B04, B06, B44, B46,
    F2, F4, F6,
    \[Alpha], \[Zeta]
    }, 
  "Nd" -> {
    B02, B04, B06, B44, B46,
    F2, F4, F6, 
    P2, 
    T2, T3, T4, T6, T7, T8,
    \[Alpha], \[Beta], \[Gamma], \[Zeta]
    }, 
  "Sm" -> {
    B02, B04, B06, B44, B46,
    F2, F4, F6,
    \[Zeta]
    }, 
  "Eu" -> {
    B02, B04, B06, B44, B46,
    F2, F4, F6,
    M0, P2,
    \[Alpha], \[Beta], \[Gamma], \[Zeta]
    }, 
  "Tb" -> {
    B02, B04, B06, B44, B46,
    F2,
    \[Zeta]
    }, 
  "Dy" -> {
    B02, B04, B06, B44,
    F2,
    \[Zeta]
    },
  "Ho" -> {
    B02, B04, B06, B44, B46,
    F2, F4, F6,
    M0,
    \[Zeta]
    },
  "Er" -> {
    B02, B04, B06, B44, B46,
    F2, F4, F6,
    M0, P2,
    T2, T3, T4, T6, T7,  T8,
    \[Alpha], \[Beta], \[Gamma],
    \[Zeta]
    },
  "Tm" -> {
    B02, B04, B06, B44, B46,
    F2, F4, F6,
    \[Zeta]
    }, 
  "Yb" -> {
    B02, B04, B44,
    \[Zeta]
    }
  |>

paramsChengLiYF4::usage="This association has the model parameters as fitted by Cheng et. al \"Crystal-field analyses for trivalent lanthanide ions in LiYF4\".";
paramsChengLiYF4 = <|
  "Ce" -> {
    \[Zeta] -> 630.,
    B02 -> 354., B04 -> -1043.,
    B44 -> -1249., B06 -> -65.,
    B46 -> -1069.
    }, 
  "Pr" -> {
    F2 -> 68955., F4 -> 50505., F6 -> 33098.,
    \[Zeta] -> 748.,
    \[Alpha] -> 23.3, \[Beta] -> -644., \[Gamma] -> 1413.,
    M0 -> 1.88, P2 -> 244.,
    B02 -> 512.,   B04 -> -1127.,
    B44 -> -1239., B06 -> -85., 
    B46 -> -1205.
    }, 
 "Nd" -> {
    F2 -> 72952., F4 -> 52681., F6 -> 35476.,
    \[Zeta] -> 877.,
    \[Alpha] -> 21., \[Beta] -> -579., \[Gamma] -> 1446.,
    T2 -> 210., T3 -> 41., T4 -> 74., T6 -> -293., T7 -> 321., T8 -> 205., 
    M0 -> 1.85, P2 -> 304.,
    B02 -> 391.,   B04 -> -1031.,
    B44 -> -1271., B06 -> -28.,
    B46 -> -1046.
    }, 
 "Sm" -> {
    F2 -> 79515., F4 -> 56766., F6 -> 40078.,
    \[Zeta] -> 1168.,
    \[Alpha] -> 20.5, \[Beta] -> -616., \[Gamma] -> 1565.,
    T2 -> 282., T3 -> 26., T4 -> 71., T6 -> -257., T7 -> 314., T8 -> 328.,
    M0 -> 2.38, P2 -> 336.,
    B02 -> 370.,  B04 -> -757.,
    B44 -> -941., B06 -> -67., 
    B46 -> -895.
    }, 
 "Eu" -> {
    F2 -> 82573., F4 -> 59646., F6 -> 43203.,
    \[Zeta] -> 1329.,
    \[Alpha] -> 21.6, \[Beta] -> -482., \[Gamma] -> 1140.,
    T2 -> 370., T3 -> 40., T4 -> 40., T6 -> -300., T7 -> 380., T8 -> 370.,
    M0 -> 2.41, P2 -> 332.,
    B02 -> 339.,   B04 -> -733.,
    B44 -> -1067., B06 -> -36., 
    B46 -> -764.
    }, 
 "Tb" -> {
    F2 -> 90972., F4 -> 64499., F6 -> 45759.,
    \[Zeta] -> 1702.,
    \[Alpha] -> 17.6, \[Beta] -> -581., \[Gamma] -> 1792.,
    T2 -> 330., T3 -> 40., T4 -> 45., T6 -> -365., T7 -> 320., T8 -> 349.,
    M0 -> 2.7, P2 -> 482.,
    B02 -> 413.,   B04 -> -867.,
    B44 -> -1114., B06 -> -41., 
    B46 -> -736.
    }, 
 "Dy" -> {
    F0 -> 0, 
    F2 -> 90421., F4 -> 63928., F6 -> 46657.,
    \[Zeta] -> 1895.,
    \[Alpha] -> 17.9, \[Beta] -> -628., \[Gamma] -> 1790.,
    T2 -> 326., T3 -> 23., T4 -> 83., T6 -> -294., T7 -> 403., T8 -> 340.,
    M0 -> 4.46, P2 -> 610.,
    B02 -> 360.,  B04 -> -737.,
    B44 -> -943., B06 -> -35.,
    B46 -> -700.
    }, 
 "Ho" -> {
    F2 -> 93512., F4 -> 66084., F6 -> 49765.,
    \[Zeta] -> 2126.,
    \[Alpha] -> 17.2, \[Beta] -> -596., \[Gamma] -> 1839.,
    T2 -> 365., T3 -> 37., T4 -> 95., T6 -> -274., T7 -> 331., T8 -> 343.,
    M0 -> 3.92, P2 -> 582.,
    B02 -> 386.,  B04 -> -629.,
    B44 -> -841., B06 -> -33., 
    B46 -> -687.
    }, 
 "Er" -> {
    F2 -> 97326., F4 -> 67987., F6 -> 53651.,
    \[Zeta] -> 2377.,
    \[Alpha] -> 18.1, \[Beta] -> -599., \[Gamma] -> 1870.,
    T2 -> 380., T3 -> 41., T4 -> 69., T6 -> -356., T7 -> 239., T8 -> 390.,
    M0 -> 4.41, P2 -> 795.,
    B02 -> 325.,   B04 -> -749.,
    B44 -> -1014., B06 -> -19., 
    B46 -> -635.
    }, 
 "Tm" -> {
    F0 -> 0.,
    T2 -> 0.,
    F2 -> 101938., F4 -> 71553., F6 -> 51359.,
    \[Zeta] -> 2632.,
    \[Alpha] -> 17.3, \[Beta] -> -665., \[Gamma] -> 1936., 
    M0 -> 4.93, P2 -> 730.,
    B02 -> 339.,  B04 -> -627.,
    B44 -> -913., B06 -> -39., 
    B46 -> -584.
    }, 
 "Yb" -> {
    \[Zeta] -> 2916.,
    B02 -> 446.,  B04 -> -560.,
    B44 -> -843., B06 -> -23.,
    B46 -> -512.
    }
  |>

StringToSLJ[string_] := Module[
  {stringed = string, LS, J, LSindex},
  (
  If[StringContainsQ[stringed, "+"],
   Return["mixed"]
   ];
  LS = StringTake[stringed, {1, 2}];
  If[StringContainsQ[stringed, "("],
   (
     LSindex = 
     StringCases[stringed, RegularExpression["\\(([^)]*)\\)"] :> "$1"];
    LS = LS <> LSindex;
    stringed = StringSplit[stringed, ")"][[-1]];
    J = ToExpression[stringed];
    ),
   (
    J = ToExpression@StringTake[stringed, {3, -1}];
    )
   ];
  {LS, J}
  )
];

FreeIonSolver::usage="This function takes a list of experimental data and the number of electrons in the lanthanide ion and returns the free-ion parameters that best fit the data. The options are:
  - F4F6_SlaterRatios: a list of two numbers that represent the ratio of F4 to F2 and F6 to F2, respectively.
  - PrintFun: a function that will be used to print the progress of the fitting process.
  - MaxIterations: the maximum number of iterations that the fitting process will run.
  - MaxMultiplets: the maximum number of multiplets that will be used in the fitting process.
  - MaxPercent: the maximum percentage of the data that can be off by the fitting.
  - SubSetBounds: a list of two numbers that represent the minimum and maximum number of multiplets that will be used in the fitting process.
The function returns an association with the following keys:
  - bestParams: the best parameters found in the fitting.
  - worstRelativeError: the worst relative error in the fitting.
  - SlaterRatios: the Slater ratios used in the fitting.
  - usedBaricenters: the baricenters used in the fitting.
If no acceptable solution is found, the function will return all solutions that are not worse than 10*MaxPercent. A solution is acceptable if the worst relative error is less than the MaxPercent option.
";
Options[FreeIonSolver] = {
  "F4F6_SlaterRatios" -> {0.707, 0.516},
  "PrintFun" -> PrintTemporary,
  "MaxIterations" -> 10000,
  "MaxMultiplets" -> 12,
  "MaxPercent" -> 3.,
  "SubSetBounds" -> {5, 12}
  };
FreeIonSolver[expData_, numE_, OptionsPattern[]] := Module[
  (* {maxMultiplets, maxPercent, F4overF2, F6overF2, PrintFun, minSubSetSize, maxSubSetSize, multipletEnergies, numMultiplets, allEqns, subsetSizes, ln, solutions, subsets, subset, eqns, m, b, meritFun, sol, goodThings, bestThings, bestOfAll, finalSol, usedMultiplets, usedBaricenters}, *)
  {},
  (
    maxMultiplets = OptionValue["MaxMultiplets"];
    maxIterations = OptionValue["MaxIterations"];
    maxPercent    = OptionValue["MaxPercent"];
    F4overF2      = OptionValue["F4F6_SlaterRatios"][[1]];
    F6overF2      = OptionValue["F4F6_SlaterRatios"][[2]];
    PrintFun      = OptionValue["PrintFun"];
    minSubSetSize = OptionValue["SubSetBounds"][[1]];
    maxSubSetSize = OptionValue["SubSetBounds"][[2]];
    freeIonParams = {F0, F2, F4, F6, \[Zeta]};
    ln            = theLanthanides[[numE]];
    
    PrintFun["Parsing the barycenters of the different multiplets ..."];
    multipletEnergies = Map[First, #] & /@ GroupBy[expData, #[[2]] &];
    multipletEnergies = Mean[Select[#, NumberQ]] & /@ multipletEnergies;
    multipletEnergies = Select[multipletEnergies, FreeQ[#, Mean] &];
    multipletEnergies = KeySelect[KeyMap[StringToSLJ, multipletEnergies], # =!= "mixed" &];
    multipletEnergies = KeyMap[Prepend[#, numE] &, multipletEnergies];
    numMultiplets     = Length[multipletEnergies];

    PrintFun["Composing the system of equations for the free-ion energies ..."];
    allEqns     = KeyValueMap[FreeIonTable[#1] == #2 &, multipletEnergies];
    allEqns     = Append[Coefficient[#[[1]], {F0, F2, F4, F6, \[Zeta]}], #[[2]]] & /@ allEqns;
    allEqns     = allEqns[[;; Min[Length[allEqns], maxMultiplets]]];
    subsetSizes = Range[1, numMultiplets];
    numSubsets = {#, Binomial[numMultiplets, #]} & /@ subsetSizes;
    numSubsets = Transpose@SortBy[numSubsets, Last];
    accSizes = Accumulate[numSubsets[[2]]];
    numSubsets = Transpose@Append[numSubsets, accSizes];
    lastSub = SelectFirst[numSubsets, #[[3]] > 1000 &, Last[numSubsets]];
    lastPosition = Position[numSubsets, lastSub][[1, 1]];
    chosenSubsetSizes = #[[1]] & /@ numSubsets[[;; lastPosition]];
    solutions   = <||>;
    
    PrintFun["Selecting subsets of different lengths and fitting with ratio-constraints ..."];
    Do[
      subsets = Subsets[Range[1, Length[allEqns]], {subsetSize}];
      PrintFun["Considering ", Length[subsets], " barycenter subsets of size ", subsetSize, " ..."];
      Do[
      (
        subset = subsets[[subsetIndex]];
        eqns   = allEqns[[subset]];
        m      = #[[;; 5]] & /@ eqns;
        b      = #[[6]] & /@ eqns;
        meritFun = Max[100 * Expand[Abs[(m . freeIonParams - b)]/b]];
        sol = NMinimize[{meritFun,
          F0 > 0,
          F2 > 0,
          F4 == F4overF2 * F2, 
          F6 == F6overF2 * F2,
          \[Zeta] > 0},
          freeIonParams, 
          MaxIterations -> maxIterations,
          Method -> "Convex"];
        solutions[{subsetSize, subset}] = sol;
        )
      , {subsetIndex, 1, Length[subsets]}
      ], 
      {subsetSize, chosenSubsetSizes}
    ];

    PrintFun["Collecting solutions of different subset size ..."];
    goodThings = Table[Normal[Sort[KeySelect[#[[1]] & /@ 
            solutions, #[[1]] == subSize &]]][[1]], 
            {subSize, subsetSizes}];
    
    PrintFun["Picking the solutions that are not worse than ", maxPercent, "% ..."];
    bestThings = Select[goodThings, #[[2]] <= maxPercent &];
    If[bestThings == {},
      Print["No acceptable solution found, consider increasing maxPercent or inspecting the given data ..."];
      Return[goodThings];
    ];

    PrintFun["Keeping the solution with the largest number of used barycenters ..."];
    bestOfAll = bestThings[[-1]];
    sol       = solutions[bestOfAll[[1]]];
    subset    = bestOfAll[[1, 2]];
    eqns      = allEqns[[subset]];
    m         = #[[;; 5]] & /@ eqns;
    b         = #[[6]] & /@ eqns;
    usedMultiplets  = Keys[multipletEnergies][[subset]];
    usedBaricenters = {#, multipletEnergies[#]} & /@ usedMultiplets;
    uniqeLS = DeleteDuplicates[#[[2]] &/@ Keys[multipletEnergies]];
    solAssoc = Association[sol[[2]]];
    usedLaF3 = False;
    If[Length[uniqeLS] == 1,
      (
        Print["There is too little data to find Slater parameters, using the ones for LaF3, and keeping the fitted spin-orbit zeta ..."];
        laf3params = LoadLaF3Parameters[ln];
        usedLaF3   = True;
        solAssoc[F0] = laf3params[F0];
        solAssoc[F2] = laf3params[F2];
        solAssoc[F4] = laf3params[F4];
        solAssoc[F6] = laf3params[F6];
      )
    ];
    finalSol = <|
      "bestParams" -> solAssoc,
      "usedLaF3" -> usedLaF3,
      "worstRelativeError" -> sol[[1]],
      "SlaterRatios" -> {F4overF2, F6overF2},
      "usedBaricenters" -> usedBaricenters|>;
    Return[finalSol];
)
];

