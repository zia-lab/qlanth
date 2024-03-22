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

ProgressNotebook::usage="ProgressNotebook[] creates a progress notebook for the solver. This notebook inclu;des a plot of the RMS history and the current parameter values. The notebook is returned. The RMS history and the parameter values are updated by setting the variables rmsHistory and paramSols. The variables stringPartialVars and paramSols are used to display the parameter values in the notebook. The notebook is created with the title \"Solver Progress\". The notebook is created with the option WindowSelected->True. The notebook is created with the option TextAlignment->Center. The notebook is created with the option WindowTitle->\"Solver Progress\".";
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
          TableForm[Prepend[Transpose[{stringPartialVars, paramSols[[-1]]}],{"RMS",rmsHistory[[-1]]}]],
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

AppendToLog[message_, file_String] := 
  Module[{timestamp = DateString["ISODateTime"], msgString},
    msgString = ToString[message, InputForm]; (* Convert any expression to a string *)
    OpenAppend[file];
    WriteString[file, timestamp, " - ", msgString, "\n"];
    Close[file];
  ];

energyAndLabelCostFunTemplate::usage="energyAndLabelCostFunTemplate is a template used to define the cost function that includes both the differences between energies and the differences between labels. The template is used to define a function TheRightSignedPath that takes a list of variables and returns the RMS of the energy differences between the computed and the experimental energies together with a term that depends on the differences between the labels. The template requires the values to the following keys to be provided: `vars` and `varPatterns`";
energyAndLabelCostFunTemplate = StringTemplate["
TheRightSignedPath[`varPatterns`]:= Module[
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
    lnParams = LoadParameters[ln];
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

TruncationFit::usage="TruncaationFit[numE, expData, numReps, activeVars, startingValues, Options] fits the given expData in an f^numE configuration, generating numReps different solutions, and varying the symbols in activeVars. The list startingValues is a list with all of the parameters needed to define the Hamiltonian (including values for activeVars, which will be disregarded but are required as position placeholders). The function returns a list of solutions. The solutions are the results of the NMinimize function using the Differential Evolution method. The solutions are a list of tuples. The first element of the tuple is the RMS error and the second element is a list of replacement rules for the fitted parameters. Once each NMinimize is done, the function saves the solutions to a file. The file is named with a prefix and a UUID. The file is saved in the log sub-directory as a .m file. The solver is always constrained by the relevant subsets of constraints for the parameters as provided by the Constrainer function. By default the Differential Evolution method starts with a generation of points within the given constraints, however it is also possible here to have a different region from which the initial points are chosen with the option \"StartingForVars\". 

The following options can be used:
  \"SignatureCheck\" : if True then then the function returns prematurely, returning a list with the symbols that would have defined the Hamiltonian after all simplifications have been applied. Useful to check the entire parameter set that the Hamiltonian has, which has to match one-to-one what is provided by startingValues.
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
  lnParams = LoadParameters[ln];
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
      multipletAssingments = Table[
        (
          J         = Js[[idx]];
          eigenVecs = theIntermediateEigensystems[[idx]][[2]];
          majorComponentIndices        = Ordering[Abs[#]][[1]]&/@eigenVecs;
          majorComponentAssignments    = LSJmultiplets[[#]]&/@majorComponentIndices;
          (* All of the degenerate eigenvectors belong to the same multiplet*)
          elevatedMultipletAssignments = ListRepeater[majorComponentAssignments,2J+1];
          elevatedMultipletAssignments
        ),
      {idx, 1, Length[Js]}
      ];
      (* Put together the multiplet assignments and the energies *)
      freeIenergiesAndMultiplets = Transpose/@Transpose[{First/@elevatedIntermediateEigensystems, multipletAssingments}];
      freeIenergiesAndMultiplets = Flatten[freeIenergiesAndMultiplets, 1];
      (* Calculate the change of basis matrix using the intermediate coupling eigenvectors *)
      basisChanger = BlockDiagonalMatrix[Transpose/@Last/@elevatedIntermediateEigensystems];
      basisChanger = SparseArray[basisChanger];
      Return[{theIntermediateEigensystems, multipletAssingments, elevatedIntermediateEigensystems, freeIenergiesAndMultiplets, basisChanger}]
    );
    
    PrintFun["Calculating the intermediate eigensystems for ",ln," using free-ion params from LaF3 ..."];
    (* Calculate intermediate coupling basis using the free-ion params for LaF3 *)
    {theIntermediateEigensystems, multipletAssingments, elevatedIntermediateEigensystems, freeIenergiesAndMultiplets, basisChanger} = TheIntermediateEigensystems@@freeBies;
    
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


ClassicalFit::usage="ClassicalFit[numE, expData, excludeDataIndices, problemVars, startValues, \[Sigma]exp, Options] fits the given expData in an f^numE configuration, by varying the symbols in problemVars. The list expData needs to be a list of lists with the only restriction that the first element of them corresponds to energies of levels. In this list, an empty value can be used to indicate known gaps in the data. Even if the energy value for a level is known (and given in expData) certain values can be omitted from the fitting procedure through the list excludeDataIndices, which correspond to indices in expData that should be skipped over.

The Hamiltonian used for fitting is version that has been truncated either by using the maximum energy given in expData or by manually setting a truncation energy using the option \"TruncationEnergy\".

The argument \[Sigma]exp is the estimated uncertainty in the differences between the calculated and the experimental energy levels. This is used to estimate the uncertainty in the fitted parameters. Admittedly this will be a rough estimate (at least on the contribution of the calculated uncertainty), but it is better than nothing and may at least provide a lower bound to the uncertainty in the fitted parameters. It is assumed that the uncertainty in the differences between the calculated and the experimental energy levels is the same for all of them.

The list startValues is a list with all of the parameters needed to define the Hamiltonian (including the initial values for problemVars). 

The function saves the solution to a file. The file is named with a prefix (controlled by the option \"FilePrefix\") and a UUID. The file is saved in the log sub-directory as a .m file.

The function returns a list of four elements. The first element is the best-fit root mean squared difference. The second element is the filename to which the log file was saved to. The third element is a list of lists of the form {paramSymbol -> {fittedVal, fittedUncertainty}}. The fourth one is the covariance matrix of the fitted parameters.

The function admits the following options with default values:
  \"MaxHistory\" : determines how long the logs for the solver can be.
  \"MaxIterations\": determines the maximum number of iterations used by NMinimize.
  \"FilePrefix\" : the prefix to use for the file name, by default this is the symbol for the lanthanide.
  
  \"AccuracyGoal\": the accuracy goal used by NMinimize, default of 5.
  \"TrucationEnergy\": if Automatic then the maximum energy in expData is taken, else it takes the value set by this option. In all cases the energies in expData are truncated to this value.
  \"PrintFun\": the function used to print progress messages, the default is PrintTemporary.
  
  \"SlackChannel\": name of the Slack channel to which to dump progress messaages, the default is None which disables this option entirely.
  \"ProgressView\": whether or not a progress window will be opened to show the progress of the solver, the default is True.
  \"SignatureCheck\": if True then then the function returns prematurely, returning a list with the symbols that would have defined the Hamiltonian after all simplifications have been applied. Useful to check the entire parameter set that the Hamiltonian has, which has to match one-to-one what is provided by startingValues.
  
  \"AppendToLogFile\": an association appended to the log file under the key \"Appendix\".
  \"MagneticSimplifier\": a list of replacement rules to simplify the Marvin and pesudo-magnetic paramters. Here the ratios of the Marvin parameters and the pseudo-magnetic parameters are defined to simplify the magnetic part of the Hamiltonian.
  \"MagFieldSimplifier\": a list of replacement rules to specify a magnetic field (in T), if set to {}, then {Bx, By, Bz} can also be used as variables to be fitted for.
  
  \"SymmetrySimplifier\": a list of replacements rules to simplify the crystal field.
  \"OtherSimplifier\": an additional list of replacement rules that are applied to the Hamiltonian before computing with it. Here the spin-spin contribution can be turned off by setting \[Sigma]SS->0, which is the default.
";
Options[ClassicalFit] = {
  "MaxHistory"       -> 200,
  "MaxIterations"    -> 100,
  "FilePrefix"       -> "",
  "ProgressView"     -> True,
  "TruncationEnergy" -> Automatic,
  "AccuracyGoal"     -> 5,
  "PrintFun"         -> PrintTemporary,
  "SlackChannel"     -> None,
  "ProgressView"     -> True,
  "SignatureCheck"   -> False,
  "AppendToLogFile"  -> <||>,
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
    B12->0,B14->0,B16->0,B34->0,B36->0,B56->0,
    S12->0,S14->0,S16->0,S22->0,S24->0,S26->0,
    S34->0,S36->0,S44->0,S46->0,S56->0,S66->0
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
    "FreeIonSymbols" -> {F0, F2, F4, F6, \[Zeta]}
};
ClassicalFit[numE_Integer, expData0_List, excludeDataIndices_List, problemVars_List, startValues_List, \[Sigma]exp_Real, OptionsPattern[]]:=
  (
    oddsAndEnds   = <||>;
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
      PrintFun["Truncation energy set to Automatic, using the maximum energy in the data ..."];
      Max[Select[First /@ expData, NumericQ[#] &]],
      OptionValue["TruncationEnergy"]
      ];
    PrintFun["Using a truncation energy of ", truncationEnergy, " K"];

    simplifier          = Join[magneticSimplifier,
                          magFieldSimplifier,
                          symmetrySimplifier,
                          threeBodySimplifier,
                          otherSimplifier];
    
    stringPartialVars = ToString/@problemVars;

    PrintFun["Determining gaps in the data ..."];
    (* the indices that are numeric in expData whatever is non-numeric is assumed as a known gap *)
    presentDataIndices = Flatten[Position[expData, {_?(NumericQ[#] &), _, _}]];
    (* some indices omitted here based on the excludeDataIndices argument *)
    presentDataIndices = Complement[presentDataIndices, excludeDataIndices];
    
    hamDim                          = Binomial[14, numE];
    oddsAndEnds["simplifier"]       = simplifier;
    oddsAndEnds["freeIonSymbols"]   = freeIonSymbols;
    oddsAndEnds["truncationEnergy"] = truncationEnergy;
    oddsAndEnds["numE"]             = numE;
    oddsAndEnds["expData"]          = expData;
    oddsAndEnds["problemVars"]      = problemVars;
    oddsAndEnds["maxIterations"]    = maxIterations;
    oddsAndEnds["hamDim"]           = hamDim;
    modelSymbols  = Sort[Select[paramSymbols, Not[MemberQ[Join[racahSymbols, chenSymbols,{t2Switch, \[Epsilon],gs}],#]]&]];
    (* remove the symbols that will be removed by the simplifier, no symbol should remain here that is not in the symbolic Hamiltonian *)
    reducedModelSymbols = Select[modelSymbols, Not[MemberQ[Keys[simplifier],#]]&];
    
    (* this is useful to understand what are the arguments of the truncated compiled Hamiltonian *)
    If[OptionValue["SignatureCheck"],
      (
        Print["Given the model parameters and the simplifying assumptions, the resultant model parameters are:"];
        Print[{reducedModelSymbols}];
        Print["The ordering in these needs to be respected in the startValues parameter ..."];
        Print["Exiting ..."];
        Return[""];
      )
    ];

    (* calculate the basis *)
    PrintFun["Retrieving the LSJMJ basis for f^", numE, " ..."];
    basis    = BasisLSJMJ[numE];

    (* get the reference parameters from LaF3 *)
    PrintFun["Getting reference free-ion parameters for ",ln," using LaF3 ..."];
    lnParams = LoadParameters[ln];
    freeBies = Prepend[Values[(#->(#/.lnParams)) &/@ freeIonSymbols], numE];
    (* a more explicit alias *)
    allVars           = reducedModelSymbols;
    standardValues    = allVars /. lnParams;
    oddsAndEnds["allVars"]  = allVars;
    oddsAndEnds["freeBies"] = freeBies;
    
    (* reload compiled version if found *)
    varHash                  = Hash[{numE, allVars, freeBies, truncationEnergy}];
    compileIntermediateFname = ln<>"-compiled-intermediate-truncated-ham-"<>ToString[varHash]<>".mx";
    compileIntermediateFname = FileNameJoin[{moduleDir, "compiled", compileIntermediateFname}];

    If[FileExistsQ[compileIntermediateFname],
      PrintFun["This ion, free-ion params, and full set of variables have been used before (as determined by {numE, allVars, freeBies, truncationEnergy}). Loading the previously saved compiled function and intermediate coupling basis ..."];
      PrintFun["Using : ", compileIntermediateFname];
      {compileIntermediateTruncatedHam, truncatedIntermediateBasis} = Import[compileIntermediateFname];,
    (
      (* grab the Hamiltonian preserving its block structure *)
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
      PrintFun["argsOfTheIntermediateEigensystems = ",argsOfTheIntermediateEigensystems];
      PrintFun["argsForEvalInsideOfTheIntermediateSystems = ",argsForEvalInsideOfTheIntermediateSystems];
      PrintFun["(if the following fails, it might help to see if the arguments of TheIntermediateEigensystems match the ones shown above)"];

      (* compile a function that will effectively calculate the spectrum of all of the scalar blocks given the parameters of the free-ion part of the Hamiltonian *)
      (* compile one function for each of the blocks *)
      PrintFun["Compiling functions for the diagonal blocks of the Hamiltonian ..."];
      compiledDiagonal = Compile[Evaluate[freeIonSymbols], Evaluate[N[Normal[#]]]]&/@diagonalScalarBlocks;
      (* use that to create a function that will calculate the free-ion eigensystem *)
      TheIntermediateEigensystems[numEv_, F0v_, F2v_, F4v_, F6v_, \[Zeta]v_] := (
        theNumericBlocks = (#[F0v, F2v, F4v, F6v, \[Zeta]v]&)/@compiledDiagonal;
        theIntermediateEigensystems = Eigensystem/@theNumericBlocks;
        Js     = AllowedJ[numEv];
        basisJ = BasisLSJMJ[numEv,"AsAssociation"->True];
        (* having calculated the eigensystems with the removed degeneracies, put the degeneracies back in explicitly *)
        elevatedIntermediateEigensystems = MapIndexed[EigenLever[#1,2Js[[#2[[1]]]]+1]&, theIntermediateEigensystems];
        (* Identify a single MJ to keep *)
        pivot = If[EvenQ[numEv],0,-1/2];
        LSJmultiplets = (#[[1]]<>ToString[InputForm[#[[2]]]])&/@Select[BasisLSJMJ[numEv],#[[-1]]== pivot &];
        (* calculate the multiplet assignments that the intermediate basis eigenvectors have *)
        multipletAssingments = Table[
          (
            J         = Js[[idx]];
            eigenVecs = theIntermediateEigensystems[[idx]][[2]];
            majorComponentIndices        = Ordering[Abs[#]][[1]]&/@eigenVecs;
            majorComponentAssignments    = LSJmultiplets[[#]]&/@majorComponentIndices;
            (* All of the degenerate eigenvectors belong to the same multiplet*)
            elevatedMultipletAssignments = ListRepeater[majorComponentAssignments,2J+1];
            elevatedMultipletAssignments
          ),
        {idx, 1, Length[Js]}
        ];
        (* put together the multiplet assignments and the energies *)
        freeIenergiesAndMultiplets = Transpose/@Transpose[{First/@elevatedIntermediateEigensystems, multipletAssingments}];
        freeIenergiesAndMultiplets = Flatten[freeIenergiesAndMultiplets, 1];
        (* calculate the change of basis matrix using the intermediate coupling eigenvectors *)
        basisChanger = BlockDiagonalMatrix[Transpose/@Last/@elevatedIntermediateEigensystems];
        basisChanger = SparseArray[basisChanger];
        Return[{theIntermediateEigensystems, multipletAssingments, elevatedIntermediateEigensystems, freeIenergiesAndMultiplets, basisChanger}]
      );

      PrintFun["Calculating the intermediate eigensystems for ",ln," using free-ion params from LaF3 ..."];
      (* calculate intermediate coupling basis using the free-ion params for LaF3 *)
      {theIntermediateEigensystems, multipletAssingments, elevatedIntermediateEigensystems, freeIenergiesAndMultiplets, basisChanger} = TheIntermediateEigensystems@@freeBies;

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
        FirstPosition[allFreeEnergiesSorted-Min[allFreeEnergiesSorted],x_/;x>truncationEnergy,{0},1][[1]]
      ];
      (* the actual energy at which the truncation is made *)
      roundedTruncationEnergy = allFreeEnergiesSorted[[sortedTruncationIndex]];

      (* the indices that enact the truncation *)
      truncationIndices = ordering[[;;sortedTruncationIndex]];
      (* Return[{basisChanger, ham, truncationIndices}]; *)
      (* Computing the block structure of the change of basis array *)
      blockSizes = BlockArrayDimensionsArray[ham];
      basisChangerBlocks = ArrayBlocker[basisChanger, blockSizes];
      blockShifts = First /@ Diagonal[blockSizes];
      numBlocks = Length[blockSizes];
      (* using the ham (with all the symbols) change the basis to the computed one *)
      PrintFun["Changing the basis of the Hamiltonian to the intermediate coupling basis ..."];
      (* intermediateHam            = Transpose[basisChanger].ham.basisChanger; *)
      (* Return[{basisChangerBlocks, ham}]; *)
      intermediateHam = BlockMatrixMultiply[ham, basisChangerBlocks];
      PrintFun["Distributing products inside of symbolic matrix elements to keep complexity in check ..."];
      Do[
        intermediateHam[[i, j]] = MapToSparseArray[intermediateHam[[i, j]], Distribute /@ # &];,
        {i, 1, numBlocks},
        {j, 1, numBlocks}
      ];
      intermediateHam = BlockMatrixMultiply[BlockTranspose[basisChangerBlocks],intermediateHam];
      PrintFun["Distributing products inside of symbolic matrix elements to keep complexity in check ..."];
      Do[
        intermediateHam[[i, j]] = MapToSparseArray[intermediateHam[[i, j]], Distribute /@ # &];,
        {i, 1, numBlocks},
        {j, 1, numBlocks}
      ];
      (* using the truncation indices truncate that one *)
      PrintFun["Truncating the Hamiltonian ..."];
      (* truncatedIntermediateHam   = intermediateHam[[truncationIndices,truncationIndices]]; *)
      truncatedIntermediateHam = TruncateBlockArray[intermediateHam, truncationIndices, blockShifts];
      (* these are the basis vectors for the truncated hamiltonian *)
      PrintFun["Saving the truncated intermediate basis ..."];
      truncatedIntermediateBasis = basisChanger[[All,truncationIndices]];
      
      PrintFun["Compiling a function for the truncated Hamiltonian ..."];
      (* compile a function that will calculate the truncated Hamiltonian given the parameters in allVars, this is the function to be use in fitting *)
      compileIntermediateTruncatedHam = Compile[Evaluate[allVars], Evaluate[truncatedIntermediateHam]];
      (* save the compiled function *)
      PrintFun["Saving the compiled function for the truncated Hamiltonian and the truncatedIntermediateBasis..."];
      Export[compileIntermediateFname, {compileIntermediateTruncatedHam, truncatedIntermediateBasis}];
    )
    ];
    
    truncationUmbral = Dimensions[truncatedIntermediateBasis][[2]];
    PrintFun["The truncated Hamiltonian has a dimension of ", truncationUmbral, "x", truncationUmbral, " ..."];
    presentDataIndices = Select[presentDataIndices, # <= truncationUmbral &];

    PrintFun["Starting up the fitting process using the Levenberg-Marquardt method ..."];
    (* using the problemVars I need to create the argument list including _?NumericQ *)
    problemVarsQ = (ToString[#] <> "_?NumericQ") & /@ problemVars;
    problemVarsQString = StringJoin[Riffle[problemVarsQ, ", "]];
    (* we also need to have the padded arguments with the variables in the right position and the fixed values in the remaining ones *)
    problemVarsPositions = Position[allVars, #][[1, 1]] & /@ problemVars;
    problemVarsString = StringJoin[Riffle[ToString /@ problemVars, ", "]];
    (* to feed parameters to the Hamiltonian, which includes all parameters, we need to form the rist set of arguments, with fixed values where needed, and the variables in the right position *)
    varsWithConstants = standardValues;
    varsWithConstants[[problemVarsPositions]] = problemVars;
    varsWithConstantsString = ToString[varsWithConstants];
    
    (* this following function serves eigenvalues from the Hamiltonian, has memoization so it might grow to use a lot of RAM *)
    Clear[HamSortedEigenvalues];
    hamEigenvaluesTemplate = StringTemplate["
    HamSortedEigenvalues[`problemVarsQ`]:=(
      ham=compileIntermediateTruncatedHam@@`varsWithConstants`;
      eigenValues = Sort@Eigenvalues@ham;
      eigenValues = eigenValues - Min[eigenValues];
      HamSortedEigenvalues[`problemVarsString`] = eigenValues;
      Return[eigenValues]
    )"];
    hamString = hamEigenvaluesTemplate[<|
        "problemVarsQ" -> problemVarsQString,
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
    eigenValueDispenserString = 
      eigenvalueDispenserTemplate[<|
        "problemVarsQ"      -> problemVarsQString,
        "problemVarsString" -> problemVarsString
        |>];
    ToExpression[eigenValueDispenserString];  

    paramSols  = {};
    rmsHistory = {};
    steps      = 0;
    problemVarsWithStartValues = Transpose[{problemVars, startValues}];
    openNotebooks     = If[runningInteractive,
              ("WindowTitle"/.NotebookInformation[#]) & /@ Notebooks[],
              {}];
    If[Not[MemberQ[openNotebooks,"Solver Progress"]] && OptionValue["ProgressView"],
      ProgressNotebook["Basic"->False]
    ];
    degressOfFreedom = Length[presentDataIndices] - Length[problemVars];
    PrintFun["Fitting for ", Length[presentDataIndices], " data points with ", Length[problemVars], " free parameters.", " The effective degrees of freedom are ", degressOfFreedom, " ..."];
    PrintFun["Starting the fitting process ..."];
    startTime=Now;
    fmSol = FindMinimum[
      Sum[(expData[[j]][[1]] - (PartialHamEigenvalues @@ problemVars)[j])^2,
      {j, presentDataIndices}],
      problemVarsWithStartValues,
      Method        -> "LevenbergMarquardt",
      MaxIterations -> OptionValue["MaxIterations"],
      AccuracyGoal  -> OptionValue["AccuracyGoal"],
      StepMonitor :> (
        steps      += 1;
        currentRMS = Sum[(expData[[j]][[1]] - (PartialHamEigenvalues @@ problemVars)[j])^2, {j, presentDataIndices}];
        currentRMS = Sqrt[currentRMS / degressOfFreedom];
        paramSols  = AddToList[paramSols, problemVars, maxHistory];
        rmsHistory = AddToList[rmsHistory, currentRMS, maxHistory];
        )
      ];
    endTime = Now;
    timeTaken = QuantityMagnitude[endTime - startTime, "Seconds"];
    Print["Solution found in ", timeTaken, "s"];
    
    sol = fmSol;
    solVec = Last /@ fmSol[[-1]];
    fullSolVec = standardValues;
    fullSolVec[[problemVarsPositions]] = solVec;
    PrintFun["Calculating the numerical Hamiltonian corresponding to the solution ..."];
    fullHam = compileIntermediateTruncatedHam @@ fullSolVec;
    PrintFun["Calculating energies and eigenvectors ..."];
    {eigenEnergies, eigenVectors} = Eigensystem[fullHam];
    states = Transpose[{eigenEnergies, eigenVectors}];
    states = SortBy[states, First];
    eigenEnergies = First /@ states;
    PrintFun[
      "Shifting energies to make ground state zero of energy ..."];
    eigenEnergies = eigenEnergies - eigenEnergies[[1]];
    PrintFun["Calculating the linear approximant to each eigenvalue ..."];
    allVarsVec = Transpose[{allVars}];
    p0 = Transpose[{fullSolVec}];
    linMat = {};
    Do[
      (
        aVarPosition = Position[allVars, aVar][[1, 1]];
        isolationValues = ConstantArray[0, Length[allVars]];
        isolationValues[[aVarPosition]] = 1;
        perHam  = compileIntermediateTruncatedHam @@ isolationValues;
        lin     = FirstOrderPerturbation[Last /@ states, perHam];
        linMat  = Append[linMat, lin];
      ),
      {aVar, problemVars}
    ];
    PrintFun["Removing the gradient of the ground state ..."];
    linMat = (# - #[[1]] & /@ linMat);
    PrintFun["Transposing derivative matrices into columns ..."];
    linMat  = Transpose[linMat];
    
    PrintFun["Calculating the eigenvalue vector at solution ..."];
    \[Lambda]0Vec = Transpose[{eigenEnergies[[presentDataIndices]]}];
    PrintFun["Putting together the experimental vector ..."];
    \[Lambda]exp = Transpose[{First /@ expData[[presentDataIndices]]}];

    problemVarsVec = Transpose[{problemVars}];
    solVecVec = Transpose[{solVec}];
    PrintFun["Calculating the difference between eigenvalues at solution ..."];
    diff = (\[Lambda]0Vec - \[Lambda]exp) + linMat[[presentDataIndices]].(problemVarsVec - solVecVec);
    PrintFun["Calculating the sum of squares of differences around solution ... "];
    sqdiff = Expand[(Transpose[diff] . diff)[[1, 1]]];
    PrintFun["Calculating the minimum (which should coincide to sol) ..."];
    minpoly = sqdiff /. AssociationThread[cfVars -> solVec];
    fmSolAssoc = Association[sol[[2]]];
    totalVariance = Length[presentDataIndices]*\[Sigma]exp^2;
    PrintFun["Calculating the uncertainty in the parameters ..."];
    solWithUncertainty = Table[
      (
        aVar        = problemVars[[varIdx]];
        paramBest   = aVar /. fmSolAssoc;
        othersFixed = AssociationThread[Delete[problemVars, varIdx] -> Delete[solVec, varIdx]];
        thisPoly    = sqdiff /. othersFixed;
        polySols    = Last /@ Last /@ Solve[thisPoly == minpoly + 1*totalVariance];
        polySols    = Select[polySols, Im[#] == 0 &];
        paramSigma  = Max[polySols] - Min[polySols];
        (aVar -> {paramBest, paramSigma})
      ),
    {varIdx, 1, Length[problemVars]}
    ];
    PrintFun["Calculating the covariance matrix ..."];
    hess = 2 * Transpose[linMat[[presentDataIndices]]] . linMat[[presentDataIndices]];
    covMat  = \[Sigma]exp^2 * Inverse[hess];
    bestRMS = Sqrt[minpoly / degressOfFreedom];
    oddsAndEnds["actualSteps"] = steps;
    oddsAndEnds["bestRMS"] = bestRMS;
    oddsAndEnds["solWithUncertainty"] = solWithUncertainty;
    oddsAndEnds["problemVars"] = problemVars;
    oddsAndEnds["paramSols"] = paramSols;
    oddsAndEnds["rmsHistory"] = rmsHistory;
    oddsAndEnds["Appendix"] = OptionValue["AppendToLogFile"];
    oddsAndEnds["timeTaken/s"] = timeTaken;
    logFname = LogSol[oddsAndEnds, logFilePrefix];
    Return[{bestRMS, logFname, solWithUncertainty, covMat}];
  );
