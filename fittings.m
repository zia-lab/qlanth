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
     WindowSize -> {590, 820},
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
           
           PlotLabel -> 
            If[Length[rmsHistory] != 0, rmsHistory[[-1]], ""]
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
     WindowSize -> {590, 820},
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
    (* Pause[0.010]; *)
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

Constrainer::usage = "Constrainer[problemVars, ln] returns a list of constraints for the variables in problemVars for trivalent lanthanide ion ln. problemVars are standard model symbols with a v appended to them.";
Constrainer[problemVars_, ln_] := (
  varcfSymbols = ToExpression[ToString[#] <> "v"] & /@ cfSymbols;
  slater = Which[
    MemberQ[{"Ce", "Yb"}, ln],
    {},
    True,
    {#, (10000. < # < 100000.)} & /@ {F2v, F4v, F6v}
    ];
  alpha = Which[
    MemberQ[{"Ce", "Yb"}, ln],
    {},
    True,
    {{\[Alpha]v, 0. < \[Alpha]v < 20.}}
    ];
  zeta = {{\[Zeta]v, 0. < \[Zeta]v < 3200.}};
  beta = Which[
    MemberQ[{"Ce", "Yb"}, ln],
    {},
    True,
    {{\[Beta]v, -1000. < \[Beta]v < 0.}}
    ];
  gamma = Which[
    MemberQ[{"Ce", "Yb"}, ln],
    {},
    True,
    {{\[Gamma]v, 1000. < \[Gamma]v < 2000.}}
    ];
  tees = Which[
    ln == "Tm",
    {-500. < T2v < 500.},
    MemberQ[{"Ce", "Pr", "Yb"}, ln],
    {},
    True,
    {#, -500. < # < 500.} & /@ {T2v, T3v, T4v, T6v, T7v, T8v}];
  marvins = Which[
    MemberQ[{"Ce", "Yb"}, ln],
    {},
    True,
    {{M0v, 0 < M0v < 5}}
    ];
  peas = Which[
    MemberQ[{"Ce", "Yb"}, ln],
    {},
    True,
    {{P2v, -200. < P2v < 1200.}}
    ];
  crystalRanges = {#, (-2000. < # < 2000.)} & /@ (Intersection[
      varcfSymbols, problemVars]);
  allCons = 
   Join[slater, zeta, alpha, beta, gamma, tees, marvins, peas, 
    crystalRanges];
  allCons = Select[allCons, MemberQ[problemVars, #[[1]]] &];
  Return[Flatten[Rest /@ allCons]]
  )

LogSol::usage = "LogSol[final, solHistory, prefix] saves the final solution and the solution history to a file. The file is named with a prefix and a UUID. The file is saved in the current directory. The file is saved in the format of a .m file. The function returns the name of the file.";
LogSol[theSolution_, solHistory_, prefix_] := (
   fname = prefix <> "-sols-" <> CreateUUID[] <> ".m";
   {optimum, fitParams} = theSolution;
   Print["Saving solution to: ", fname];
   exporter = Association[{"bestRMS" -> optimum,
      "solHistory" -> solHistory,
      "prefix"     -> prefix,
      "fitParams"  -> fitParams}];
   Export[fname, exporter];
   Return[fname];
   );
  

FitToHam::usage = "FitToHam[numE, expData, fitToSymbols, simplifier, OptionsPattern[]] fits the model Hamiltonian to the experimental data for the trivalent lanthanide ion with number numE. The experimental data is given in the form of a list of tuples. The first element of the tuple is the energy and the second element is the label. 
The fitToSymbols is a list of the symbols to be fit. The simplifier is a list of rules that simplify the Hamiltonian. 
The options and their defaults are:
\"PrintFun\"->PrintTemporary,
\"SlackChannel\"->None,
\"MaxHistory\"->100,
\"MaxIters\"->100,
\"NumCycles\"->10,
\"ProgressWindow\"->True
The PrintFun option is the function used to print progress messages.
The SlackChannel option is the channel to post progress messages to.
The MaxHistory option is the maximum number of iterations to keep in the history.
The MaxIters option is the maximum number of iterations for the solver.
The NumCycles option is the number of cycles to run the solver for.
The function returns a list of solutions. The solutions are the results of the NMinimize function. The solutions are a list of tuples. The first element of the tuple is the RMS error and the second element is the parameter values
The function also saves the solutions to a file. The file is named with a prefix and a UUID. The file is saved in the current directory. The file is saved in the format of a .m file.";
Options[FitToHam] = {
      "PrintFun"-> PrintTemporary,
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
    
    ln = theLanthanides[[numE]];
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
      logFname = LogSol[sol, rmsHistory, "label-and-energy"];
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
