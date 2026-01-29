Needs["qlanth`"];
Needs["fittings`"];
Needs["qonstants`"];
Needs["misc`"];

LoadCarnall[];

workDir = DirectoryName[$InputFileName];

FastIonSolverLaF3Carnall::usage = "FastIonSolverLaF3Carnall[numE] solves the energy levels of the given trivalent lanthanide in LaF3 using the parameter values for the  Hamiltonian as taken from the values quoted by Carnall. It can use precomputed symbolic matrices for the Hamiltonian if they have been loaded already and defined as symbols of the form symbolicHamiltonians[__].

The function returns a list with nine elements {rmsDifference, carnallEnergies, eigenEnergies, ln, carnallAssignments, simplerStateLabels, eigensys, basis, truncatedStates}. The elements of the list are as follows:

1. rmsDifference is the root mean squared difference between the calculated values and those quoted by Carnall;
2. carnallEnergies are the quoted calculated energies from Carnall used for comparison;
3. eigenEnergies are the calculated energies (in the case of an odd number of electrons the Kramers degeneracy may have been removed from this list according to the option \"Remove Kramers\");
4. ln is simply a string labelling the corresponding lanthanide;
5. carnallAssignments is a list of strings providing the multiplet assignments that Carnall assumed;
6. simplerStateLabels is a list of strings providing the multiplet assignments that this function assumes;
7. eigensys is a list of tuples where the first element is the energy corresponding to the eigenvector given as the second element (in the case of an odd number of electrons the Kramers degeneracy may have been removed from this list according to the option \"Remove Kramers\");
8. basis is a list that specifies the basis in which the Hamiltonian was constructed and diagonalized, equal to BasisLSJMJ[numE];
9. truncatedStates is the same as eigensys but with the truncated eigenvectors so that the total probability adds up to at least eigenstateTruncationProbability (determined by the option \"eigenstateTruncationProbability\").

This function admits the following options:
  - \"MakeNotebook\" -> True or False. If True, a notebook with a summary of the data is created. Default is True.
  - \"NotebookSave\" -> True or False. If True, the results notebook is saved automatically. Default is True.
  - \"eigenstateTruncationProbability\" -> 0.9. The probability sum of the truncated eigenvectors. Default is 0.9.
  - \"Include spin-spin\" -> True or False. If True, the spin-spin contribution to the magnetic interactions is included. Default is True.
  - \"Max Eigenstates in Table\" -> 100. The maximum number of eigenstates to be shown in the table shown in the results notebook. Default is 100.
  - \"Sparse\" -> True or False. If True, the numerical Hamiltonian is kept in sparse form. Default is True.
  - \"PrintFun\" -> Print, PrintTemporary, or other to serve as a printer for progress messages. Default is Print.
  - \"SaveData\" -> True or False. If True, the resulting data is saved to disk. Default is True.
  - \"ParamOverride\". An association that can override parameters in the Hamiltonian. Default is <||>. This override cannot change the inclusion or exclusion of the spin-spin contribution to the magnetic interactions, for this purpose use the option \"Include spin-spin\".
  - \"Append to Filename\" -> \"\". A string to append to the filename of the saved notebook and data files. Default is \"\".
  - \"Remove Kramers\" -> True or False. If True, the Kramers degeneracy is removed from the eigenstates. Default is True.
  - \"OutputDirectory\" -> \"calcs\". The directory where the output files are saved. Default is \"calcs\".
  - \"Explorer\" -> True or False. If True, the energy level diagram is interactive. Default is False.
";
Options[FastIonSolverLaF3Carnall] = {
  "Host" -> "LaF3",
  "Sparse" -> True,
  "SaveData" -> True,
  "PrintFun" -> Print,
  "Explorer" -> False,
  "MakeNotebook" -> True,
  "NotebookSave" -> True,
  "ParamOverride" -> <||>,
  "Remove Kramers" -> True,
  "Append to Filename" -> "",
  "Include spin-spin" -> True,
  "OutputDirectory" -> "calcs",
  "Max Eigenstates in Table" -> 100,
  "eigenstateTruncationProbability" -> 0.9
}; 
FastIonSolverLaF3Carnall[numE_, OptionsPattern[]] := Module[
  {makeNotebook, eigenstateTruncationProbability, host,
  ln, terms, termNames, carnallEnergies, eigenEnergies, simplerStateLabels,
  eigensys, basis, assignmentMatches, stateLabels, carnallAssignments, params},
  (
    PrintFun     = OptionValue["PrintFun"];
    makeNotebook = OptionValue["MakeNotebook"];
    eigenstateTruncationProbability = OptionValue["eigenstateTruncationProbability"];
    maxStatesInTable                = OptionValue["Max Eigenstates in Table"];
    Duplicator[aList_] := Flatten[{#, #} & /@ aList];
    host          = OptionValue["Host"];
    ParamOverride = OptionValue["ParamOverride"];
    ln        = theLanthanides[[numE]];
    terms     = AllowedNKSLJTerms[Min[numE, 14 - numE]];
    termNames = First /@ terms;
    (* For labeling the states, the degeneracy in some of the terms is elided *)
    PrintFun["> Calculating simpler term labels ..."];
    termSimplifier = Table[termN -> If[StringLength[termN] == 3,
      StringTake[termN, {1, 2}],
      termN
      ],
      {termN, termNames}
    ];
    
    (*Load the parameters from Carnall*)
    PrintFun["> Loading the fit parameters from Carnall ..."];
    params = LoadLaF3Parameters[ln, "Vintage" -> "Carnall", "With Uncertainties" -> False];
    If[numE>7,
      (
        PrintFun["> Conjugating the parameters accounting for the hole-particle equivalence ..."];
        params = HoleElectronConjugation[params];
        params[t2Switch] = 0;
      ),
      params[t2Switch] = 1;
    ];
    
    (* Apply the parameter override *)
    Do[params[key] = ParamOverride[key], 
      {key, Keys[ParamOverride]}
    ];
    
    (* Import the symbolic Hamiltonian *)
    PrintFun["> Loading the symbolic Hamiltonian for this configuration ..."];
    startTime = Now;
    numH = 14 - numE;
    numEH = Min[numE, numH];
    C2vsimplifier = {
      B12 -> 0, B14 -> 0, B16 -> 0, B34 -> 0, B36 -> 0, B56 -> 0,
      S12 -> 0, S14 -> 0, S16 -> 0, S22 -> 0, S24 -> 0, S26 -> 0, S34 -> 0, S36 -> 0,S44 -> 0, S46 -> 0, S56 -> 0, S66 -> 0,
      T11p -> 0, T12 -> 0, T14 -> 0, T15 -> 0, T16 -> 0, T18 -> 0, T17 -> 0, T19 -> 0};
    (* If the necessary symbolicHamiltonian is defined load it, if not, then calculate it *)
    simpleHam = If[
      ValueQ[symbolicHamiltonians[numEH]],
      symbolicHamiltonians[numEH],
      SimplerEffectiveHamiltonian[numE,
        C2vsimplifier,
        "PrependToFilename" -> "C2v-",
        "Overwrite" -> False]
    ];
    endTime  = Now; 
    loadTime = QuantityMagnitude[endTime - startTime, "Seconds"];
    PrintFun[">> Loading the symbolic Hamiltonian took ", loadTime, " seconds."];
    
    (*Enforce the override to the spin-spin contribution to the magnetic interactions*)
    params[\[Sigma]SS] = If[OptionValue["Include spin-spin"], 1, 0];
    
    (*Everything that is not given is set to zero*)
    params = ParamPad[params, "PrintFun" -> PrintFun];
    PrintFun[params];
    numHam = ReplaceInSparseArray[simpleHam, params];
    If[Not[OptionValue["Sparse"]],
        numHam = Normal[numHam]
    ];
    PrintFun["> Calculating the SLJ basis ..."];
    basis = BasisLSJMJ[numE];
    
    (* Eigensolver *)
    PrintFun["> Diagonalizing the numerical Hamiltonian ..."];
    startTime    = Now;
    eigensys     = Eigensystem[numHam];
    endTime      = Now;
    diagonalTime = QuantityMagnitude[endTime - startTime, "Seconds"];
    PrintFun[">> Diagonalization took ", diagonalTime, " seconds."];
    eigensys     = Chop[eigensys];
    eigensys     = Transpose[eigensys];
    
    (* Shift the baseline energy *)
    eigensys = ShiftedLevels[eigensys];
    (* Sort according to energy *)
    eigensys = SortBy[eigensys, First];
    (* Grab just the energies *)
    eigenEnergies = First /@ eigensys;
    
    (* Energies are doubly degenerate in the case of odd number of electrons, keep only one *)
    If[And[OddQ[numE], OptionValue["Remove Kramers"]], 
        (
        PrintFun["> Since there's an odd number of electrons energies come in pairs, taking just one for each pair ..."];
        eigenEnergies = eigenEnergies[[;; ;; 2]];
        )
    ];
    
    (* Compare against the data quoted by Bill Carnall *)
    PrintFun["> Comparing against the data from Carnall ..."];
    mainKey            = StringTemplate["appendix:`Ln`:Association"][<|"Ln" -> ln|>];
    lnData             = Carnall[mainKey];
    carnalKeys         = lnData // Keys;
    repetitions        = Length[lnData[#]["Calc (1/cm)"]] & /@ carnalKeys;
    carnallAssignments = First /@ Carnall["appendix:" <> ln <> ":RawTable"];
    carnallAssignments = Select[carnallAssignments, Not[# === ""] &];
    carnalKey          = StringTemplate["appendix:`Ln`:Calculated"][<|"Ln" -> ln|>];
    carnallEnergies    = Carnall[carnalKey];

    If[And[OddQ[numE], Not[OptionValue["Remove Kramers"]]],
    (
      PrintFun[">> Duplicating the data from Carnall since the Kramers degenerate states where not removed ..."];
      carnallAssignments = Duplicator[carnallAssignments];
      carnallEnergies    = Duplicator[carnallEnergies];
    )
    ];
    
    (* For the difference take as many energies as quoted by Bill *)
    eigenEnergies = eigenEnergies + carnallEnergies[[1]];
    diffs = Sort[eigenEnergies][[;; Length[carnallEnergies]]] - carnallEnergies;
    (* Remove the differences where the appendix tables have elided values*)
    rmsDifference = Sqrt[Mean[(Select[diffs, FreeQ[#, Missing[]] &])^2]];
    titleTemplate = StringTemplate[
      "Energy Level Diagram of \!\(\*SuperscriptBox[\(`ion`\), \(\(3\)\(+\)\)]\)"];
    title = titleTemplate[<|"ion" -> ln|>];
    parsedStates = ParseStates[eigensys, basis];
    If[And[OddQ[numE],OptionValue["Remove Kramers"]], 
      parsedStates = parsedStates[[;; ;; 2]]
    ];
    
    stateLabels = #[[-1]] & /@ parsedStates;
    simplerStateLabels = ((#[[2]] /. termSimplifier) <> ToString[#[[3]], InputForm]) & /@ parsedStates;
    
    PrintFun[">> Truncating eigenvectors to given probability ..."];
    startTime = Now;
    truncatedStates = ParseStatesByProbabilitySum[eigensys, basis, 
      eigenstateTruncationProbability, 
      0.01];
    endTime = Now;
    truncationTime = QuantityMagnitude[endTime - startTime, "Seconds"];
    PrintFun[">>> Truncation took ", truncationTime, " seconds."];
    
    If[makeNotebook,
    (
      PrintFun["> Putting together results in a notebook ..."];
      energyDiagram = Framed[
        EnergyLevelDiagram[eigensys, "Title" -> title, 
        "Explorer" -> OptionValue["Explorer"],
        "Background" -> White]
        , Background -> White, FrameMargins -> 50];
      appToFname = OptionValue["Append to Filename"];
      PrintFun[">> Comparing the term assignments between qlanth and Carnall ..."];
      AssignmentMatchFunc = Which[
        StringContainsQ[#[[1]], #[[2]]],
        "\[Checkmark]",
        True,
        "X"] &;
      assignmentMatches = AssignmentMatchFunc /@ Transpose[{carnallAssignments, simplerStateLabels[[;; Length[carnallAssignments]]]}];
      assignmentMatches = {{"\[Checkmark]", 
        Count[assignmentMatches, "\[Checkmark]"]}, {"X", 
        Count[assignmentMatches, "X"]}};
      labelComparison = (AssignmentMatchFunc /@ Transpose[{carnallAssignments, simplerStateLabels[[;; Length[carnallAssignments]]]}]);
      labelComparison = PadRight[labelComparison, Length[simplerStateLabels], "-"];
      
      statesTable = Grid[Prepend[{Round[#[[1]]], #[[2]]} & /@ 
          truncatedStates[[;;Min[Length[eigensys],maxStatesInTable]]], {"Energy/\!\(\*SuperscriptBox[\(cm\), \(-1\)]\)", 
          "\[Psi]"}], Frame -> All, Spacings -> {2, 2}, 
          FrameStyle -> Blue, 
          Dividers -> {{False, True, False}, {True, True}}];
      DefaultIfMissing[expr_]:= If[FreeQ[expr, Missing[]], expr,"NA"];
      PrintFun[">> Rounding the energy differences for table presentation ..."];
      roundedDiffs = Round[diffs, 0.1];
      roundedDiffs = PadRight[roundedDiffs, Length[simplerStateLabels], "-"];
      roundedDiffs = DefaultIfMissing /@ roundedDiffs;
      diffs = PadRight[diffs, Length[simplerStateLabels], "-"];
      diffs = DefaultIfMissing /@ diffs;
      diffTableData = Transpose[{simplerStateLabels, eigenEnergies,
        labelComparison,
        PadRight[carnallAssignments, Length[simplerStateLabels], "-"],
        DefaultIfMissing/@PadRight[carnallEnergies, Length[simplerStateLabels], "-"], 
        roundedDiffs}
      ];
      diffTable = TableForm[diffTableData, 
        TableHeadings -> {None, {"qlanth", 
        "E/\!\(\*SuperscriptBox[\(cm\), \(-1\)]\)", "", "Carnall", 
        "E/\!\(\*SuperscriptBox[\(cm\), \(-1\)]\)", 
        "\[CapitalDelta]E/\!\(\*SuperscriptBox[\(cm\), \(-1\)]\)"}}
      ];
      
      diffs = Sort[eigenEnergies][[;; Length[carnallEnergies]]] - carnallEnergies;
      notBad = FreeQ[#,Missing[]]&/@diffs;
      diffs = Pick[diffs, notBad];
      (* diffHistogram = Histogram[diffs,
        Frame -> True,
        ImageSize -> 800, 
        AspectRatio -> 1/3, FrameStyle -> Directive[16], 
        FrameLabel -> {"(qlanth-carnall)/Ky", "Freq"}
      ]; *)
      
      rmsDifference = Sqrt[Total[diffs^2/Length[diffs]]];
      labelTempate  = StringTemplate["\!\(\*SuperscriptBox[\(`ln`\), \(\(3\)\(+\)\)]\)"];
      diffData   = diffs;
      diffLabels = simplerStateLabels[[;;Length[notBad]]];
      diffLabels = Pick[diffLabels, notBad];
      diffPlot   = Framed[
          qplotter`ListLabelPlot[
            diffData[[;;;;If[OddQ[numE],2,1]]],
            diffLabels[[;;;;If[OddQ[numE],2,1]]],
            Frame -> True,
            PlotRange -> All,
            ImageSize -> 1200,
            AspectRatio -> 1/3,
            Filling -> Axis,
            FrameLabel -> {"", 
            "(qlanth-carnall) / \!\(\*SuperscriptBox[\(cm\), \(-1\)]\)"},
            PlotMarkers -> "OpenMarkers",
            PlotLabel -> 
            Style[labelTempate[<|"ln" -> ln|>] <> " | " <> "\[Sigma]=" <> 
              ToString[Round[rmsDifference, 0.01]] <> 
              " \!\(\*SuperscriptBox[\(cm\), \(-1\)]\)\n", 20],
            Background -> White
          ],
          Background -> White,
          FrameMargins -> 50
        ];
      (* now place all of this in a new notebook *)
      nb = CreateDocument[
      {
        TextCell[Style[
          DisplayForm[RowBox[{SuperscriptBox[host <> ":" <> ln, "3+"], "(", SuperscriptBox["f", numE], ")"}]]
          ], "Title", TextAlignment -> Center
        ],
        TextCell["Calculated Energy Diagram",
          "Section",
          TextAlignment -> Center
        ],
        TextCell[energyDiagram,
          TextAlignment -> Center
        ],
        TextCell["Multiplet Assignments & Energy Levels",
          "Section",
          TextAlignment -> Center
        ],
        (* TextCell[diffHistogram, TextAlignment -> Center], *)
        TextCell[diffPlot, "Output", TextAlignment -> Center],
        TextCell[assignmentMatches, "Output", TextAlignment -> Center],
        TextCell[diffTable, "Output", TextAlignment -> Center],
        TextCell["Truncated Eigenstates", "Section", TextAlignment -> Center],
        TextCell["These are some of the resultant eigenstates which add up to at least a total probability of " <> ToString[eigenstateTruncationProbability] <> ".", "Text", TextAlignment -> Center],
        TextCell[statesTable, "Output", TextAlignment -> Center]
      },
      WindowSelected -> True,
      WindowTitle -> ln <> " in " <> "LaF3" <> appToFname,
      WindowSize -> {1600, 800}];
      If[OptionValue["SaveData"],
      (
        exportFname = FileNameJoin[{workDir, OptionValue["OutputDirectory"], ln <> " in " <> "LaF3" <> appToFname <> ".mx"}];
        SelectionMove[nb, After, Notebook];
        NotebookWrite[nb, Cell["Reload Data", "Section", TextAlignment -> Center]];
        NotebookWrite[nb,
          Cell[(
            "{rmsDifference, carnallEnergies, eigenEnergies, ln, carnallAssignments, simplerStateLabels, eigensys, basis, truncatedStates} = Import[FileNameJoin[{NotebookDirectory[],\"" <> StringSplit[exportFname,"/"][[-1]] <> "\"}]];"
            ),"Input"
          ]
        ];
        NotebookWrite[nb,
          Cell[(
            "Manipulate[First[MinimalBy[truncatedStates, Abs[First[#] - energy] &]], {energy,0}]"
            ),"Input"]
        ];
        (* Move the cursor to the top of the notebook *)
        SelectionMove[nb, Before, Notebook];
        Export[exportFname,
          {rmsDifference, carnallEnergies, eigenEnergies, ln, carnallAssignments, simplerStateLabels, eigensys, basis, truncatedStates}
        ];
        tinyexportFname = FileNameJoin[
          {workDir, OptionValue["OutputDirectory"], ln <> " in " <> "LaF3" <> appToFname <> " - tiny.m"}
        ];
        tinyExport = <|"ln"->ln,
                        "carnallEnergies"->carnallEnergies,
                        "rmsDifference"-> rmsDifference,
                        "eigenEnergies"-> eigenEnergies,
                        "carnallAssignments"-> carnallAssignments,
                        "simplerStateLabels" -> simplerStateLabels|>;
        Export[tinyexportFname, tinyExport];
      )
      ];
      If[OptionValue["NotebookSave"],
        (
          nbFname = FileNameJoin[{workDir,OptionValue["OutputDirectory"], ln <> " in " <> "LaF3" <> appToFname <> ".nb"}];
          PrintFun[">> Saving notebook to ", nbFname, " ..."];
          NotebookSave[nb, nbFname];
        )
      ];
    )
    ];

    Return[{rmsDifference, carnallEnergies, 
            eigenEnergies, ln,
            carnallAssignments, simplerStateLabels,
            eigensys, basis,
            truncatedStates}];
  )
];

FlowMatching::usage="FlowMatching[aList, bList] returns a list of pairs of elements from aList and bList that are closest to each other, this is returned in a list together with a mapping of indices from the aList to those in bList to which they were matched. The option \"alistLabels\" can be used to specify labels for the elements in aList. The option \"blistLabels\" can be used to specify labels for the elements in bList. If these options are used, the function returns a list with three elements the pairs of matched elements, the pairs of corresponding matched labels, and the mapping of indices. This is basically a wrapper around Mathematica's FindMinimumCostFlow function. By default the option \"notMatched\" is zero, and this means that all elements of aList must be matched to elements of bList. If this is not the case, the option \"notMatched\" can be used to specify how many elements of aList can be left unmatched. By default the cost function is Abs[#1-#2]&, but this can be changed with the option \"CostFun\", this function needs to take two arguments.";
Options[FlowMatching] = {"alistLabels" -> {}, "blistLabels" -> {}, "notMatched" -> 0, "CostFun"-> (Abs[#1-#2] &)};
FlowMatching[aValues0_, bValues0_, OptionsPattern[]] := Module[{
aValues = aValues0, bValues = bValues0, edgesSourceToA, capacitySourceToA, nA, nB,
costSourceToA, midLayer, midLayerEdges, midCapacities,
midCosts, edgesBtoSink, capacityBtoSink, costBtoSink,
allCapacities, allCosts, allEdges, graph,
flow, bestValues, bestLabels, cFun,
aLabels, bLabels, pairedIndices, matchingLabels},
(
matchingLabels = (Length[OptionValue["alistLabels"]] > 0);
aLabels = OptionValue["alistLabels"];
bLabels = OptionValue["blistLabels"];
cFun    = OptionValue["CostFun"];
nA      = Length[aValues];
nB      = Length[bValues];
(*Build up the edges costs and capacities*)
(*From source to the nodes representing the values of the first \
list*)
edgesSourceToA    = ("source" \[DirectedEdge] {"A", #}) & /@ Range[1, nA];
capacitySourceToA = ConstantArray[1, nA];
costSourceToA     = ConstantArray[0, nA];

(*From all the elements of A to all the elements of B*)
midLayer = Table[{{"A", i} \[DirectedEdge] ({"B", j}), 1, cFun[aValues[[i]], bValues[[j]]]}, {i, 1, nA}, {j, 1, nB}];
midLayer = Flatten[midLayer, 1];
{midLayerEdges, midCapacities, midCosts} = Transpose[midLayer];

(*From the elements of B to the sink*)
edgesBtoSink    = ({"B", #} \[DirectedEdge] "sink") & /@ Range[1, nB];
capacityBtoSink = ConstantArray[1, nB];
costBtoSink     = ConstantArray[0, nB];

(*Put it all together*)
allCapacities = Join[capacitySourceToA, midCapacities, capacityBtoSink];
allCosts      = Join[costSourceToA, midCosts, costBtoSink];
allEdges      = Join[edgesSourceToA, midLayerEdges, edgesBtoSink];
graph         = Graph[allEdges, EdgeCapacity -> allCapacities, 
  EdgeCost -> allCosts];

(*Solve it*)
flow          = FindMinimumCostFlow[graph, "source", "sink", nA - OptionValue["notMatched"], "OptimumFlowData"];
(*Collect the pairs of matched indices*)
pairedIndices = Select[flow["EdgeList"], And[Not[#[[1]] === "source"], Not[#[[2]] === "sink"]] &];
pairedIndices = {#[[1, 2]], #[[2, 2]]} & /@ pairedIndices;
(*Collect the pairs of matched values*)
bestValues    = {aValues[[#[[1]]]], bValues[[#[[2]]]]} & /@ pairedIndices;
(*Account for having been given labels*)
If[matchingLabels,
  (
  bestLabels = {aLabels[[#[[1]]]], bLabels[[#[[2]]]]} & /@ pairedIndices;
  Return[{bestValues, bestLabels, pairedIndices}]
  ),
  (
  Return[{bestValues, pairedIndices}]
  )
  ];
)
];

GreedyMatching::usage="GreedyMatching[aList, bList] returns a list of pairs of elements from aList and bList that are closest to each other, this is returned in a list together with a mapping of indices from the aList to those in bList to which they were matched. The option \"alistLabels\" can be used to specify labels for the elements in aList. The option \"blistLabels\" can be used to specify labels for the elements in bList. If these options are used, the function returns a list with three elements the pairs of matched elements, the pairs of corresponding matched labels, and the mapping of indices.";
Options[GreedyMatching] = {
    "alistLabels" -> {},
    "blistLabels" -> {}};
GreedyMatching[aValues0_, bValues0_, OptionsPattern[]] := Module[{
  aValues = aValues0,
  bValues = bValues0,
  bValuesOriginal = bValues0,
  bestLabels, bestMatches,
  bestLabel, aElement, givenLabels,
  aLabels, aLabel,
  diffs, minDiff,
  bLabels,
  minDiffPosition, bestMatch},
  (
  aLabels     = OptionValue["alistLabels"];
  bLabels     = OptionValue["blistLabels"];
  bestMatches = {};
  bestLabels  = {};
  givenLabels = (Length[aLabels] > 0);
  Do[
    (
    aElement        = aValues[[idx]];
    diffs           = Abs[bValues - aElement];
    minDiff         = Min[diffs];
    minDiffPosition = Position[diffs, minDiff][[1, 1]];
    bestMatch       = bValues[[minDiffPosition]];
    bestMatches     = Append[bestMatches, {aElement, bestMatch}];
    If[givenLabels,
      (
      aLabel     = aLabels[[idx]];
      bestLabel  = bLabels[[minDiffPosition]];
      bestLabels = Append[bestLabels, {aLabel, bestLabel}];
      bLabels    = Drop[bLabels, {minDiffPosition}];
      )
      ];
    bValues = Drop[bValues, {minDiffPosition}];
    If[Length[bValues] == 0, Break[]];
    ),
    {idx, 1, Length[aValues]}
    ];
  pairedIndices = MapIndexed[{#2[[1]], Position[bValuesOriginal, #1[[2]]][[1, 1]]} &, bestMatches];
  If[givenLabels,
    Return[{bestMatches, bestLabels, pairedIndices}],
    Return[{bestMatches, pairedIndices}]
    ]
  )
  ];

StochasticMatching::usage="StochasticMatching[aValues, bValues] finds a better assignment by randomly shuffling the elements of aValues and then applying the greedy assignment algorithm. The function prints what is the range of total absolute differences found during shuffling, the standard deviation of all of them, and the number of shuffles that were attempted. The option \"alistLabels\" can be used to specify labels for the elements in aValues. The option \"blistLabels\" can be used to specify labels for the elements in bValues. If these options are used, the function returns a list with three elements the pairs of matched elements, the pairs of corresponding matched labels, and the mapping of indices.";
Options[StochasticMatching] = {"alistLabels" -> {}, 
  "blistLabels" -> {}};
StochasticMatching[aValues0_, bValues0_, numShuffles_ : 200, OptionsPattern[]] := Module[{
  aValues = aValues0,
  bValues = bValues0,
  matchingLabels, ranger, matches, noShuff, bestMatch, highestCost, lowestCost, dev, sorter, bestValues,
  pairedIndices, bestLabels, matchedIndices, shuffler
  },
  (
  matchingLabels = (Length[OptionValue["alistLabels"]] > 0);
  ranger = Range[1, Length[aValues]];
  matches = If[Not[matchingLabels], (
      Table[(
        shuffler = If[i == 1, ranger, RandomSample[ranger]];
        {bestValues, matchedIndices} = 
        GreedyMatching[aValues[[shuffler]], bValues];
        cost = Total[Abs[#[[1]] - #[[2]]] & /@ bestValues];
        {cost, {bestValues, matchedIndices}}
        ), {i, 1, numShuffles}]
      ),
    Table[(
      shuffler = If[i == 1, ranger, RandomSample[ranger]];
      {bestValues, bestLabels, matchedIndices} = 
        GreedyMatching[aValues[[shuffler]], bValues, 
        "alistLabels" -> OptionValue["alistLabels"][[shuffler]], 
        "blistLabels" -> OptionValue["blistLabels"]];
      cost = Total[Abs[#[[1]] - #[[2]]] & /@ bestValues];
      {cost, {bestValues, bestLabels, matchedIndices}}
      ), {i, 1, numShuffles}]
    ];
  noShuff = matches[[1, 1]];
  matches = SortBy[matches, First];
  bestMatch = matches[[1, 2]];
  highestCost = matches[[-1, 1]];
  lowestCost = matches[[1, 1]];
  dev = StandardDeviation[First /@ matches];
  Print[lowestCost, " <-> ", highestCost, " | \[Sigma]=", dev, 
    " | N=", numShuffles, " | null=", noShuff];
  If[matchingLabels,
    (
    {bestValues, bestLabels, matchedIndices} = bestMatch;
    sorter = Ordering[First /@ bestValues];
    bestValues = bestValues[[sorter]];
    bestLabels = bestLabels[[sorter]];
    pairedIndices = 
      MapIndexed[{#2[[1]], Position[bValues, #1[[2]]][[1, 1]]} &, 
      bestValues];
    Return[{bestValues, bestLabels, pairedIndices}]
    ),
    (
    {bestValues, matchedIndices} = bestMatch;
    sorter = Ordering[First /@ bestValues];
    bestValues = bestValues[[sorter]];
    pairedIndices = 
      MapIndexed[{#2[[1]], Position[bValues, #1[[2]]][[1, 1]]} &, 
      bestValues];
    Return[{bestValues, pairedIndices}]
    )
    ];
  )
  ];

MagneticDipoleTransitionsLaF3Carnall::usage = "MagneticDipoleTransitionsLaF3Carnall[numE] calculates the magnetic dipole transitions for the lanthanide ion numE in LaF3, using the output files from FastIonSolverLaF3Carnall stored in the examples folder. The output is a tabular file, a raw data file, and a CSV file. The tabular file contains the following columns: 
    \[Psi]i:simple, (* main contribution to the wavefuction |i>*)
    \[Psi]f:simple, (* main contribution to the wavefuction |j>*)
    \[Psi]i:idx,    (* index of the wavefuction |i>*)
    \[Psi]f:idx,    (* index of the wavefuction |j>*)
    Ei/K,           (* energy of the initial state in K *)
    Ef/K,           (* energy of the final state in K *)
    \[Lambda]/nm,   (* transition wavelength in nm *)
    \[CapitalDelta]\[Lambda]/nm, (* uncertainty in the transition wavelength in nm *)
    \[Tau]/s,       (* radiative lifetime in s *)
    AMD/s^-1        (* magnetic dipole transition rate in s^-1 *)

The raw data file contains the following keys: 
    - Line Strength, (* Line strength array *)
    - AMD, (* Magnetic dipole transition rates in 1/s *)
    - fMD, (* Oscillator strengths from ground to excited states *)
    - Radiative lifetimes, (* Radiative lifetimes in s *)
    - Transition Energies / K, (* Transition energies in K *)
    - Transition Wavelengths in nm. (* Transition wavelengths in nm *)

The CSV file contains the same information as the tabular file. 

The function also creates a notebook with a Manipulate that allows the user to select a wavelength interval and a lifetime power of ten. The results notebook is saved in the examples directory.

The function takes the following options: 
    - \"Make Notebook\" -> True or False. If True, a notebook with a Manipulate is created. Default is True.
    - \"Print Function\" -> PrintTemporary or Print. The function used to print the progress of the calculation. Default is PrintTemporary.
    - \"Host\" -> \"HostName\". Name for the host material. Default is \"LaF3\".
    - \"Wavelength Range\" -> {50,2000}. The range of wavelengths in nm for the Manipulate object in the created notebook. Default is {50,2000}.

The function returns an association containing the following keys: Line Strength, AMD, fMD, Radiative lifetimes, Transition Energies / K, Transition Wavelengths in nm.";
Options[MagneticDipoleTransitionsLaF3Carnall] = {
  "Host" -> "LaF3",
  "Make Notebook" -> True, 
  "Close Notebook" -> True,
  "Wavelength Range" -> {50, 2000},
  "Print Function" -> PrintTemporary,
  "OutputDirectory" -> "calcs"
  };
MagneticDipoleTransitionsLaF3Carnall[numE_Integer, OptionsPattern[]]:= (
  host           = OptionValue["Host"];
  \[Lambda]Range = OptionValue["Wavelength Range"];
  PrintFun       = OptionValue["Print Function"];
  {\[Lambda]min, \[Lambda]max} = OptionValue["Wavelength Range"];
  
  header    = {"\[Psi]i:simple","\[Psi]f:simple","\[Psi]i:idx","\[Psi]f:idx","Ei/K","Ef/K","\[Lambda]/nm","\[CapitalDelta]\[Lambda]/nm","\[Tau]/s","AMD/s^-1"};
  ln        = {"Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb"}[[numE]];
  fname = FileNameJoin[{workDir, OptionValue["OutputDirectory"],ln<>" in " <> host <> " - example.mx"}];
  {rmsDifference, carnallEnergies, eigenEnergies, ln, carnallAssignments, simplerStateLabels, eigensys, basis, truncatedStates} = Import[fname];
  
  (* Some of the above are not needed here *)
  Clear[truncatedStates];
  Clear[basis];
  Clear[rmsDifference];
  Clear[carnallEnergies];
  Clear[carnallAssignments];
  If[OddQ[numE],
      eigenEnergies      = eigenEnergies[[;;;;2]];
      simplerStateLabels = simplerStateLabels[[;;;;2]];
      eigensys           = eigensys[[;;;;2]];
  ];
  eigenEnergies = eigenEnergies - eigenEnergies[[1]];
  
  magIon = <||>;
  (* PrintFun["Calculating the magnetic dipole line strength array..."];
  magIon["Line Strength"] = MagDipLineStrength[eigensys, numE, "Reload MagOp" -> False, "Units" -> "SI"]; *)
  
  PrintFun["Calculating the M1 spontaneous transition rates ..."];
  magIon["AMD"] = MagDipoleRates[eigensys, numE, "Units"->"SI","Lifetime"->False];
  magIon["AMD"] = magIon["AMD"]/.{0.->Indeterminate};
  
  PrintFun["Calculating the oscillator strengths for transitions from the ground state ..."];
  magIon["fMD"] = GroundMagDipoleOscillatorStrength[eigensys, numE];
  
  PrintFun["Calculating the natural radiative lifetimes ..."];
  magIon["Radiative lifetimes"] = 1/magIon["AMD"];
  
  PrintFun["Calculating the transition energies in K ..."];
  transitionEnergies = Outer[Subtract, First/@eigensys,First/@eigensys];
  magIon["Transition Energies / K"] = ReplaceDiagonal[transitionEnergies,Indeterminate];
  
  PrintFun["Calculating the transition wavelengths in nm ..."];
  magIon["Transition Wavelengths in nm"] = 10^7/magIon["Transition Energies / K"];
  
  PrintFun["Estimating the uncertainties in \[Lambda]/nm assuming a 1 K uncertainty in energies."];
  (*Assuming an uncertainty of 1 K in both energies used to calculate the wavelength*)
  \[Lambda]uncertainty = Sqrt[2]*magIon["Transition Wavelengths in nm"]^2*10^-7;
  
  PrintFun["Formatting a tabular output file ..."];
  numEigenvecs    = Length[eigensys];
  roundedEnergies = Round[eigenEnergies, 1.];
  simpleFromTo    = Outer[{#1,#2}&, simplerStateLabels,  simplerStateLabels];
  fromTo          = Outer[{#1,#2}&, Range[numEigenvecs], Range[numEigenvecs]];
  energyPairs     = Outer[{#1,#2}&, roundedEnergies,     roundedEnergies];
  allTransitions  = {simpleFromTo,
      fromTo,
      energyPairs,
      magIon["Transition Wavelengths in nm"],
      \[Lambda]uncertainty,
      magIon["AMD"],
      magIon["Radiative lifetimes"]
  };
  allTransitions = (Flatten/@Transpose[Flatten[#,1]&/@allTransitions]);
  allTransitions = Select[allTransitions, #[[3]]!=#[[4]]&];
  allTransitions = Select[allTransitions, #[[10]]>0&];
  allTransitions = Transpose[allTransitions];
  
  (*round things up*)
  PrintFun["Rounding wavelengths according to estimated uncertainties ..."];
  {roundedWaves,roundedDeltas} = Transpose[MapThread[RoundValueWithUncertainty,{allTransitions[[7]],allTransitions[[8]]}]];
  allTransitions[[7]]          = roundedWaves;
  allTransitions[[8]]          = roundedDeltas;
  
  PrintFun["Rounding lifetimes and transition rates to three significant figures ..."];
  allTransitions[[9]]          = RoundToSignificantFigures[#,3]&/@(allTransitions[[9]]);
  allTransitions[[10]]         = RoundToSignificantFigures[#,3]&/@(allTransitions[[10]]);
  finalTable                   = Transpose[allTransitions];
  finalTable                   = Prepend[finalTable,header];
  
  (* tabular output *)
  basename     = ln <> " in " <> host <> " - example - " <> "M1 - tabular.zip";
  exportFname  = FileNameJoin[{"./examples",basename}];
  PrintFun["Exporting tabular data to "<>exportFname<>" ..."];
  exportKey    = StringReplace[basename,".zip"->".m"];
  Export[exportFname, <|exportKey->finalTable|>];
  
  (* raw data output *)
  basename       = ln <> " in " <> host <> " - example - " <> "M1 - raw.zip";
  rawexportFname = FileNameJoin[{"../laf3-examples",basename}];
  PrintFun["Exporting raw data as an association to "<>exportFname<>" ..."];
  rawexportKey   = StringReplace[basename,".zip"->".m"];
  Export[rawexportFname, <|rawexportKey->magIon|>];
  
  (* csv output *)
  PrintFun["Formatting and exporting a CSV output..."];
  csvOut = Table[
      StringJoin[Riffle[ToString[#,CForm]&/@finalTable[[i]],","]],
  {i,1,Length[finalTable]}
  ];
  csvOut      = StringJoin[Riffle[csvOut, "\n"]];
  basename    = ln <> " in " <> host <> " - example - " <> "M1.csv";
  exportFname = FileNameJoin[{workDir, OptionValue["OutputDirectory"], basename}];
  PrintFun["Exporting csv data to "<>exportFname<>" ..."];
  Export[exportFname, csvOut, "Text"];
  
  If[OptionValue["Make Notebook"],
  (
      PrintFun["Creating a notebook with a Manipulate to select a wavelength interval and a lifetime power of ten ..."];
      finalTable     = Rest[finalTable];
      finalTable     = SortBy[finalTable,#[[7]]&];
      opticalTable   = Select[finalTable,\[Lambda]min<=#[[7]]<=\[Lambda]max&];
      pows           = Sort[DeleteDuplicates[(MantissaExponent[#[[9]]][[2]]-1)&/@opticalTable]];
      
      man            = Manipulate[
      (
          {\[Lambda]min, \[Lambda]max} = \[Lambda]int;
          table = Select[opticalTable,And[(\[Lambda]min<=#[[7]]<=\[Lambda]max),
                      (MantissaExponent[#[[9]]][[2]]-1)==log10\[Tau]]&];
          tab   = TableForm[table,TableHeadings->{None,header}];
          Column[{{"\[Lambda]min="<>ToString[\[Lambda]min]<>" nm","\[Lambda]max="<>ToString[\[Lambda]max]<>" nm",log10\[Tau]},tab}]
      ),
      {{\[Lambda]int,\[Lambda]Range,"\[Lambda] interval"},
          \[Lambda]Range[[1]],
          \[Lambda]Range[[2]],
          50,
          ControlType->IntervalSlider
      },
      {{log10\[Tau],pows[[-1]]},
          pows
      },
      TrackedSymbols  :> {\[Lambda]int,log10\[Tau]},
      SaveDefinitions -> True
      ];

      nb = CreateDocument[{
          TextCell[Style[DisplayForm[RowBox[{"Magnetic Dipole Transitions", "\n", SuperscriptBox[host<>":"<>ln,"3+"],"(",SuperscriptBox["f",numE],")"}]]],"Title",TextAlignment->Center],
          (* TextCell["Magnetic Dipole Transition Lifetimes","Section",TextAlignment->Center], *)
          TextCell[man,"Output",TextAlignment->Center]
        },
        WindowSelected -> True,
        WindowTitle    -> "M1 - "<>ln<>" in "<>host,
        WindowSize     -> {1600,800}
      ];
      SelectionMove[nb, After, Notebook];
      NotebookWrite[nb, Cell["Reload Data", "Section", TextAlignment -> Center]];
      NotebookWrite[nb, Cell[(
            "magTransitions = Import[FileNameJoin[{NotebookDirectory[],\"" <> StringSplit[rawexportFname,"/"][[-1]] <> "\"}],\""<>rawexportKey<>"\"];"
            ),"Input"]];
      SelectionMove[nb, Before, Notebook];
      nbFname = FileNameJoin[{workDir,OptionValue["OutputDirectory"], "M1 - "<>ln<>" in "<>"LaF3"<>".nb"}];
      PrintFun[">> Saving notebook to ",nbFname," ..."];
      NotebookSave[nb, nbFname];
      If[OptionValue["Close Notebook"],
          NotebookClose[nb];
      ];
  )
  ];

  Return[magIon];
  );

FastIonSolver::usage = "FastIonSolver[params] calculates energy levels and eigenvectors for the given parameters describing a trivalent lanthanide. It can use precomputed symbolic matrices for the Hamiltonian if they have been loaded already and defined as symbols of the form symbolicHamiltonians[{symmetry, numE}].

The function returns a list with six elements {eigenEnergies, ln, simplerStateLabels, eigensys, basis, truncatedStates}. The elements of the list are as follows:

1. eigenEnergies are the calculated energies (in the case of an odd number of electrons the Kramers degeneracy may have been removed from this list according to the option \"Remove Kramers\");
2. ln is simply a string labelling the corresponding lanthanide;
3. simplerStateLabels is a list of strings providing the multiplet assignments that this function assumes;
3. eigensys is a list of tuples where the first element is the energy corresponding to the eigenvector given as the second element (in the case of an odd number of electrons the Kramers degeneracy may have been removed from this list according to the option \"Remove Kramers\");
5. basis is a list that specifies the basis in which the Hamiltonian was constructed and diagonalized, equal to BasisLSJMJ[numE];
6. truncatedStates is the same as eigensys but with the truncated eigenvectors so that the total probability adds up to at least eigenstateTruncationProbability (determined by the option \"eigenstateTruncationProbability\").

This function admits the following options:
  - \"MakeNotebook\" -> True or False. If True, a notebook with a summary of the data is created. Default is True.
  - \"NotebookSave\" -> True or False. If True, the results notebook is saved automatically. Default is True.
  - \"eigenstateTruncationProbability\" -> 0.9. The probability sum of the truncated eigenvectors. Default is 0.9.
  - \"Max Eigenstates in Table\" -> 100. The maximum number of eigenstates to be shown in the table shown in the results notebook. Default is 100.
  - \"Sparse\" -> True or False. If True, the numerical Hamiltonian is kept in sparse form. Default is True.
  - \"PrintFun\" -> Print, PrintTemporary, or other to serve as a printer for progress messages. Default is Print.
  - \"SaveData\" -> True or False. If True, the resulting data is saved to disk. Default is True.
  - \"Append to Filename\" -> \"\". A string to append to the filename of the saved notebook and data files. Default is \"\".
  - \"Remove Kramers\" -> True or False. If True, the Kramers degeneracy is removed from the eigenstates for odd electron cases. Default is True.
  - \"OutputDirectory\" -> \"calcs\". The directory where the output files are saved. Default is \"calcs\", may be a folder not yet created.
  - \"Explorer\" -> True or False. If True, the energy level diagram is interactive. Default is False.
  - \"Host\" -> \"crystalHost\". The host material. Default is \"LaF3\".
  - \"SymmetrySimplifier\" a list of two elements, the first one being the name of the assumed symmetry for the crystal field and the second being a list of replacement rules applied to simplify the Hamiltonian used to make calculations. Default is {\"C2v\", {B12 -> 0, ...}}
";
Options[FastIonSolver] = {
  "Sparse"   -> True,
  "Host"     -> "LaF3",
  "SaveData" -> True,
  "PrintFun" -> Print,
  "Explorer" -> False,
  "MakeNotebook" -> True,
  "NotebookSave" -> True,
  "EnergiesOnly" -> False,
  "Remove Kramers"     -> True,
  "Append to Filename" -> "",
  "Energy Uncertainty in K" -> 1,
  "OutputDirectory"    -> "calcs",
  "Max Eigenstates in Table" -> 100,
  "eigenstateTruncationProbability" -> 0.9,
  "AddToExport" -> <||>,
  "SymmetrySimplifier" -> {"C2v",
    {
      B12 -> 0,  B14 -> 0, B16 -> 0, B34 -> 0, B36 -> 0, B56 -> 0,
      S12 -> 0,  S14 -> 0, S16 -> 0, S22 -> 0, S24 -> 0, S26 -> 0, S34 -> 0, S36 -> 0, S44 -> 0,  S46 -> 0, S56 -> 0, S66 -> 0,
      T11p -> 0, T12 -> 0, T14 -> 0, T15 -> 0, T16 -> 0, T18 -> 0, T17 -> 0, T19 -> 0
    }
    }
}; 
FastIonSolver[params0_, OptionsPattern[]] := Module[
  {
    params = Association[params0],
    paramKeys, numE, PrintFun, makeNotebook, eigenstateTruncationProbability,
    maxStatesInTable, Duplicator, host, ln, terms, termNames, termSimplifier,
    startTime, numH, numEH, simpleHam, endTime, loadTime, numHam,
    basis, eigensys, eigenEnergies, titleTemplate, title, parsedStates,
    stateLabels, simplerStateLabels, truncatedStates, truncationTime,
    energyDiagram, appToFname, statesTable, DefaultIfMissing, b, 
    diffTableData, diffTable, nb, exportFname, tinyexportFname, tinyExport,
    nbFname, uncertaintySentence, sigmaRMS, energyUncertaintyTemplate, onlyEnergies, originalParams, epiThings
  },
  (
    paramKeys = Keys[params];
    Which[
      MemberQ[paramKeys, nE],
      numE = params[nE],
      MemberQ[paramKeys, "nE"],
      numE = params["nE"],
      MemberQ[paramKeys, "numE"],
      numE = params["numE"],
      True,
      (
        Print["The number of electrons was not specified in params, exiting ..."];
        Return[Null];
      )
    ];
    NiceSort[params_] := (
      Association@SortBy[Normal@params, Position[paramSymbols, #[[1]]] &]);
    energyUncertaintyTemplate = StringTemplate["\nWhen the Hamiltonian was fitted to experimental data, the root mean square deviation between the data and the calculated energies was `energySigma` cm^-1.\nThe red points in the diagram above (if present) indicate the energies of the known experimental data."];
    sigmaRMS     = OptionValue["Energy Uncertainty in K"];
    uncertaintySentence = If[
      sigmaRMS === 0.,
      "",
      energyUncertaintyTemplate[<|"energySigma" -> sigmaRMS|>]
    ];
    PrintFun     = OptionValue["PrintFun"];
    makeNotebook = OptionValue["MakeNotebook"];
    onlyEnergies = OptionValue["EnergiesOnly"];
    eigenstateTruncationProbability = OptionValue["eigenstateTruncationProbability"];
    maxStatesInTable                = OptionValue["Max Eigenstates in Table"];
    Duplicator[aList_] := Flatten[{#, #} & /@ aList];
    
    host      = OptionValue["Host"];
    ln        = theLanthanides[[numE]];
    terms     = AllowedNKSLJTerms[Min[numE, 14 - numE]];
    termNames = First /@ terms;
    
    (* For labeling the states, the degeneracy in some of the terms is elided *)
    PrintFun["> Calculating simpler term labels ..."];
    termSimplifier = Table[termN -> If[StringLength[termN] == 3,
      StringTake[termN, {1, 2}],
      termN
      ],
      {termN, termNames}
    ];
    originalParams = NiceSort[params]; 
    If[numE>7,
      (
        PrintFun["> Conjugating the parameters accounting for the hole-particle equivalence ..."];
        params = HoleElectronConjugation[params];
        params[t2Switch] = (numE-2)/(numE-7);
      ),
      params[t2Switch] = 1;
    ];
    
    (* Import the symbolic Hamiltonian *)
    PrintFun["> Loading the symbolic Hamiltonian for this configuration ..."];
    startTime = Now;
    numH = 14 - numE;
    numEH = Min[numE, numH];
    {symGroup, symmetrySimplifier} = OptionValue["SymmetrySimplifier"];
    (* If the necessary symbolicHamiltonian is defined load it, if not, then calculate it *)
    simpleHam = If[
      ValueQ[symbolicHamiltonians[{symGroup, numEH}]],
      symbolicHamiltonians[{symGroup, numEH}],
      SimplerEffectiveHamiltonian[numE,
        symmetrySimplifier,
        "PrependToFilename" -> symGroup <> "-",
        "Overwrite" -> False]
    ];
    endTime  = Now; 
    loadTime = QuantityMagnitude[endTime - startTime, "Seconds"];
    PrintFun[">> Loading the symbolic Hamiltonian took ", loadTime, " seconds."];
    
    (* Enforce the override to the spin-spin contribution to the magnetic interactions *)
    If[Not@MemberQ[paramKeys, \[Sigma]SS],
      (
        PrintFun["> params does not implicitly specify if the spin-spin contribution should be included. Including it by default ..."];
        params[\[Sigma]SS] = 1;
      )
    ];
    
    (*Everything that is not given is set to zero*)
    params = ParamPad[params, "PrintFun" -> PrintFun];
    PrintFun[params];
    numHam = ReplaceInSparseArray[simpleHam, params];
    If[Not[OptionValue["Sparse"]],
        numHam = Normal[numHam]
    ];
    PrintFun["> Calculating the SLJ basis ..."];
    basis = BasisLSJMJ[numE];
    
    (* Eigensolver *)
    If[onlyEnergies,
      (
        PrintFun["> Diagonalizing the numerical Hamiltonian to get just the energies ..."];
        startTime     = Now;
        eigenEnergies = Eigenvalues[numHam];
        endTime       = Now;
        diagonalTime  = QuantityMagnitude[endTime - startTime, "Seconds"];
        PrintFun[">> Diagonalization took ", diagonalTime, " seconds."];
        eigenEnergies = Sort[eigenEnergies];
        eigenEnergies = eigenEnergies - eigenEnergies[[1]];
        If[MemberQ[Keys[params], \[Epsilon]],
              eigenEnergies += params[\[Epsilon]]
        ];
        (* Energies are doubly degenerate in the case of odd number of electrons, keep only one as long as the option \"Remove Kramers\" is set to True *)
        If[And[OddQ[numE], OptionValue["Remove Kramers"]], 
            (
            PrintFun["> Since there's an odd number of electrons energies come in pairs, taking just one for each pair ..."];
            eigenEnergies = eigenEnergies[[;; ;; 2]];
            )
        ];
        Return[eigenEnergies];
      ),
      (
        PrintFun["> Diagonalizing the numerical Hamiltonian ..."];
        startTime    = Now;
        eigensys     = Eigensystem[numHam];
        endTime      = Now;
        diagonalTime = QuantityMagnitude[endTime - startTime, "Seconds"];
        PrintFun[">> Diagonalization took ", diagonalTime, " seconds."];
        eigensys     = Chop[eigensys];
        eigensys     = Transpose[eigensys];  
      )
    ];

    (* Shift the baseline energy *)
    eigensys = ShiftedLevels[eigensys];
    (* Sort according to energy *)
    eigensys = SortBy[eigensys, First];
    energyShift = 0;
    energyShift = If[
      MemberQ[Keys[params], \[Epsilon]],
      Chop[Round[params[\[Epsilon]],0.1]],
      0.
    ];
    If[\[Epsilon]=!=0.,
      eigensys = ({#[[1]]+energyShift, #[[2]] } & /@ eigensys)
    ];
    (* Grab just the energies *)
    eigenEnergies = Chop[First /@ eigensys];
    
    
    (* Energies are doubly degenerate in the case of odd number of electrons, keep only one as long as the option \"Remove Kramers\" is set to True *)
    If[And[OddQ[numE], OptionValue["Remove Kramers"]], 
        (
        PrintFun["> Since there's an odd number of electrons energies come in pairs, taking just one for each pair ..."];
        eigenEnergies = eigenEnergies[[;; ;; 2]];
        )
    ];
    
    titleTemplate = StringTemplate[
      "Energy Level Diagram of \!\(\*SuperscriptBox[\(`ion`\), \(\(3\)\(+\)\)]\)"];
    title = titleTemplate[<|"ion" -> ln|>];
    parsedStates = ParseStates[eigensys, basis];
    If[And[OddQ[numE],OptionValue["Remove Kramers"]], 
      parsedStates = parsedStates[[;; ;; 2]]
    ];
    
    stateLabels = #[[-1]] & /@ parsedStates;
    simplerStateLabels = ((#[[2]] /. termSimplifier) <> ToString[#[[3]], InputForm]) & /@ parsedStates;
    
    PrintFun[">> Truncating eigenvectors to given probability ..."];
    startTime = Now;
    truncatedStates = ParseStatesByProbabilitySum[eigensys, basis, 
      eigenstateTruncationProbability, 
      0.01];
    endTime = Now;
    truncationTime = QuantityMagnitude[endTime - startTime, "Seconds"];
    PrintFun[">>> Truncation took ", truncationTime, " seconds."];
    
    If[makeNotebook,
    (
      hamImage = HamTeX[numE, "T2" -> (params[T2]  === 0.)];
      epiThings = Which[
        host == "LaF3",
        (
          carnalKey = "appendix:" <> ln <> ":RawTable";
          expData = {#[[2]], #[[1]], #[[3]]} & /@ Carnall[carnalKey];
          expData = Select[First /@ expData, NumericQ];
          {Red, PointSize[0.002], Point[{#, 0.5}]} & /@ expData
        ),
        True,
        {}
      ];
      PrintFun["> Putting together results in a notebook ..."];
      energyDiagram = Framed[
            EnergyLevelDiagram[eigensys, "Title" -> title, 
            "Explorer" -> OptionValue["Explorer"],
            "Background" -> White,
            "Epilog"->epiThings],
        Background -> White,
        FrameMargins -> 50];
      appToFname = OptionValue["Append to Filename"];
    
      statesTable = {Round[#[[1]]], #[[2]]} & /@ truncatedStates[[;;Min[Length[eigensys], maxStatesInTable]]];
      statesTableHeader = {"Energy/\!\(\*SuperscriptBox[\(cm\), \(-1\)]\)", "\[Psi]"};
      statesTable = Prepend[statesTable, statesTableHeader];
      statesTable = Grid[statesTable,
            Frame      -> All,
            Spacings   -> {2, 2}, 
            FrameStyle -> Blue,
            Dividers   -> {{False, True, False}, {True, True}}
            ];
      DefaultIfMissing[expr_]:= If[FreeQ[expr, Missing[]], expr,"NA"];
      diffTableData = Transpose[{simplerStateLabels, Round[#,1] & /@ eigenEnergies}];
      diffTable     = TableForm[diffTableData, 
        TableHeadings -> {None, {"coarse label", 
        "E/\!\(\*SuperscriptBox[\(cm\), \(-1\)]\)"}}];
      PrintFun[{"numE",numE}];
      (* now place all of this in a new notebook *)
      exportParams = Select[originalParams, And[# =!= 0, # =!= 0.] &];
      infoBlurb = Cell[TextData[{
 "The entire calculated electronic structure is contained in the file that \
may be loaded in the line below.\nIt loads an association with the following \
keys:\n\n",
 StyleBox["\[OpenCurlyDoubleQuote]energies/K\[CloseCurlyDoubleQuote]",
  FontWeight->"Bold"],
 " : A list with all the calculated energies, explicitly including Kramers\
\[CloseCurlyQuote] degeneracy.\n",
 StyleBox["\[OpenCurlyDoubleQuote]lnSymbol\[CloseCurlyDoubleQuote]",
  FontWeight->"Bold"],
 ": The symbol representing the current lanthanide.\n",
 StyleBox["\[OpenCurlyDoubleQuote]coarseStateLabels\[CloseCurlyDoubleQuote]",
  FontWeight->"Bold"],
 ": A list with the largest LSJ contribution to the eigenvectors.\n",
 StyleBox["\[OpenCurlyDoubleQuote]eigenSystem\[CloseCurlyDoubleQuote]:",
  FontWeight->"Bold"],
 " A list of lists, the first element being the energy and the second the \
eigenvector in the given basis.\n",
 StyleBox["\[OpenCurlyDoubleQuote]truncatedStates\[CloseCurlyDoubleQuote]",
  FontWeight->"Bold"],
 ": A truncated version of the eigenvectors.\n",
 StyleBox["\[OpenCurlyDoubleQuote]params\[CloseCurlyDoubleQuote]",
  FontWeight->"Bold"],
 ": The fitted parameters of the semi-empirical Hamiltonian.\n",
 StyleBox["\[OpenCurlyDoubleQuote]fitResult\[CloseCurlyDoubleQuote]",
  FontWeight->"Bold"],
 ": An association containing the full result of the solver that produced the \
parameters used here.\n"
}], "Text", TextAlignment->Center];
      nb = CreateDocument[
      {
        TextCell[
          Style[
            DisplayForm[RowBox[{SuperscriptBox[host <> ":" <> ln, "3+"], "(", SuperscriptBox["f", numE], ")"}]]
          ],
          "Title",
          TextAlignment -> Center
        ],
        TextCell[Show[hamImage, ImageSize -> 600],
          "Output",
          TextAlignment -> Center
        ],
        TextCell[exportParams,
          "Output",
          TextAlignment -> Center
        ],
        TextCell["Calculated Energy Diagram",
          "Section",
          TextAlignment -> Center
        ],
        TextCell[energyDiagram,
          TextAlignment -> Center
        ],
        TextCell["Multiplet Assignments & Energy Levels",
          "Section",
          TextAlignment -> Center
        ],
        TextCell[
          "The following table shows the multiplet assignments and energies of the eigenstates.\nThe energies are rounded to the nearest cm^-1 .\nThe coarse label for each state corresponds to the LSJ term with the largest contribution to the corresponding eigenvector."<> uncertaintySentence,
          "Text",
          TextAlignment -> Center
        ],
        TextCell[diffTable,
            "Output",
            TextAlignment -> Center
        ],
        TextCell["Truncated Eigenstates",
            "Section",
            TextAlignment -> Center
        ],
        TextCell[
            "These are some of the eigenstates, truncated each to add up to at least a total probability of " <> ToString[eigenstateTruncationProbability] <> ".", 
            "Text",
            TextAlignment -> Center
        ],
        TextCell[statesTable,
            "Output",
            TextAlignment -> Center
        ]
      },
      WindowSelected -> True,
      WindowTitle -> ln <> " in " <> host <> appToFname,
      WindowSize -> {1600, 800}];
      If[OptionValue["SaveData"],
      (
        exportFname = FileNameJoin[{workDir,OptionValue["OutputDirectory"], ln <> " in " <> host <> appToFname <> ".mx"}];
        SelectionMove[nb, After, Notebook];
        NotebookWrite[nb, Cell["Reload Data", "Section", TextAlignment -> Center]];
        NotebookWrite[nb, infoBlurb];
        NotebookWrite[nb,
          Cell[(
            "lnStructure = Import[FileNameJoin[{NotebookDirectory[],\"" <> StringSplit[exportFname,"/"][[-1]] <> "\"}]];"
            ),"Input"
          ]
        ];
        NotebookWrite[nb,
          Cell[(
            "Manipulate[First[MinimalBy[lnStructure[\"truncatedStates\"], Abs[First[#] - energy] &]], {energy,0}]"
            ),"Input"]
        ];
        exportAssoc = <|
          "energies/K" -> eigenEnergies,
          "lnSymbol" -> ln,
          "coarseStateLabels" -> simplerStateLabels,
          "eigenSystem" -> eigensys,
          "basis" -> basis,
          "truncatedStates" -> truncatedStates,
          "params" -> params
        |>; 
        If[OptionValue["AddToExport"] =!= <||>,
          exportAssoc = Join[exportAssoc, OptionValue["AddToExport"]];
        ]; 
        Export[exportFname, exportAssoc];
        tinyexportFname = FileNameJoin[
          {workDir, OptionValue["OutputDirectory"], ln <> " in " <> host  <> appToFname <> " - tiny.m"}
        ];
        tinyExport = <|"ln"->ln,
                    "eigenEnergies"-> eigenEnergies,
                    "simplerStateLabels" -> simplerStateLabels,
                    "params" -> params|>;
        If[OptionValue["AddToExport"] =!= <||>,
          tinyExport = Join[tinyExport, OptionValue["AddToExport"]];
        ];
        Export[tinyexportFname, tinyExport];
      )
      ];
      NotebookWrite[nb,
      Cell["\n\n\nThis notebook was generated on " <> DateString[] <> ".", "Text", TextAlignment -> Right]
      ];
      (* Move the cursor to the top of the notebook *)
      SelectionMove[nb, Before, Notebook];
      If[OptionValue["NotebookSave"],
        (
          nbFname = FileNameJoin[{workDir,
                OptionValue["OutputDirectory"],
                ln <> " in " <> host <> appToFname <> ".nb"}];
          PrintFun[">> Saving notebook to ", nbFname, " ..."];
          NotebookSave[nb, nbFname];
        )
      ];
    )
    ];
    Return[{eigenEnergies, ln, simplerStateLabels, eigensys, basis, truncatedStates}];
  )
];

MagneticDipoleTransitions::usage = "MagneticDipoleTransitions[params, fname] calculates the magnetic dipole transition rates and oscillator strengths for the lanthanide ion describe by params. This is done using the data saved in fname which corresponds to output files from FastIonSolverLaF3. This function saves a tabular file, a raw data file, and a CSV file. The tabular file contains the following columns: 
    \[Psi]i:simple, (* main contribution to the wavefuction |i>*)
    \[Psi]f:simple, (* main contribution to the wavefuction |j>*)
    \[Psi]i:idx,    (* index of the wavefuction |i>*)
    \[Psi]f:idx,    (* index of the wavefuction |j>*)
    Ei/K,           (* energy of the initial state in K *)
    Ef/K,           (* energy of the final state in K *)
    \[Lambda]/nm,   (* transition wavelength in nm *)
    \[CapitalDelta]\[Lambda]/nm, (* uncertainty in the transition wavelength in nm *)
    \[Tau]/s,       (* radiative lifetime in s *)
    AMD/s^-1        (* magnetic dipole transition rate in s^-1 *)

The raw data file contains the following keys: 
    - AMD, (* Magnetic dipole transition rates in 1/s *)
    - fMD, (* Oscillator strengths from ground to excited states *)
    - Radiative lifetimes, (* Radiative lifetimes in s *)
    - Transition Energies / K, (* Transition energies in K *)
    - Transition Wavelengths in nm. (* Transition wavelengths in nm *)

The CSV file contains the same information as the tabular file. 

The function also creates a notebook with a Manipulate that allows the user to select a wavelength interval and a lifetime power of ten. The results notebook is saved in the examples directory.

The function takes the following options: 
    - \"Make Notebook\" -> True or False. If True, the results notebook is created. Default is True.
    - \"Print Function\" -> PrintTemporary or Print. The function used to print the progress of the calculation. Default is PrintTemporary.
    - \"Host\" -> \"CrystalHost\". The host material. Default is LaF3.
    - \"Wavelength Range\" -> {50,2000}. The range of wavelengths in nm for the Manipulate object in the created notebook. Default is {50,2000}.

The function returns an association containing the following keys: \"AMD/s^-1\", \"fMD\", \"radLifetimes/s\", \"transitionEnergies/K\", \"transitionWavelengths/nm\".";
Options[MagneticDipoleTransitions] = {
        "Make Notebook" -> True, 
        "Close Notebook" -> True,
        "Print Function" -> PrintTemporary,
        "Host" -> "LaF3",
        "Energy Uncertainty in K" -> 1.,
        "Wavelength Range" -> {50, 2000}};
MagneticDipoleTransitions[params_Association, fname_String, OptionsPattern[]]:= Module[
    {numE, host, PrintFun, \[Lambda]Range,
    header, ln, energy\[CapitalDelta], b},
  (
    paramKeys    = Keys[params];
    Which[
        MemberQ[paramKeys, nE],
            numE = params[nE],
        MemberQ[paramKeys, "nE"],
            numE = params["nE"],
        MemberQ[paramKeys, "numE"],
            numE = params["numE"],
        True,
        (
        Print["The number of electrons was not specified in params, exiting ..."];
        Return[Null];
        )
    ];
    Print[{"numE",numE}];
    host           = OptionValue["Host"];
    \[Lambda]Range = OptionValue["Wavelength Range"];
    PrintFun       = OptionValue["Print Function"];
    {\[Lambda]min, \[Lambda]max} = OptionValue["Wavelength Range"];
    
    header    = {"\[Psi]i:simple","\[Psi]f:simple","\[Psi]i:idx","\[Psi]f:idx","Ei/K","Ef/K","\[Lambda]/nm","\[CapitalDelta]\[Lambda]/nm","\[Tau]*n^3/s","AMD/n^3/s^-1"};
    ln        = {"Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb"}[[numE]];
    
    PrintFun["Loading the data from ", fname, " ..."];
    lnStructure   = Import[fname];
    eigenEnergies = lnStructure["energies/K"];
    simplerStateLabels = lnStructure["coarseStateLabels"];
    eigensys           = lnStructure["eigenSystem"];
    fittedParams       = lnStructure["params"];
    
    (* Some of the above are not needed here *)
    Clear[lnStructure];
    Clear[truncatedStates];
    Clear[basis];
    If[OddQ[numE],
        eigenEnergies      = eigenEnergies[[;;;;2]];
        simplerStateLabels = simplerStateLabels[[;;;;2]];
        eigensys           = eigensys[[;;;;2]];
    ];
    eigenEnergies = eigenEnergies - eigenEnergies[[1]];

    magIon = <||>;
    (* PrintFun["Calculating the magnetic dipole line strength array..."];
    magIon["Line Strength"] = MagDipLineStrength[eigensys, numE, "Reload MagOp" -> False, "Units" -> "SI"]; *)
    
    PrintFun["Calculating the M1 spontaneous transition rates ..."];
    magIon["AMD/n^3/s^-1"] = MagDipoleRates[eigensys, numE, "Units"->"SI","Lifetime"->False];
    magIon["AMD/n^3/s^-1"] = magIon["AMD/n^3/s^-1"]/.{0.->Indeterminate};

    PrintFun["Calculating the oscillator strength for transitions from the ground state ..."];
    magIon["fMD/n"] = GroundMagDipoleOscillatorStrength[eigensys, numE];

    PrintFun["Calculating the natural radiative lifetimes ..."];
    magIon["radLifetimes*n^3/s"] = 1./magIon["AMD/n^3/s^-1"];

    PrintFun["Calculating the transition energies in K ..."];
    transitionEnergies = Outer[Subtract, First/@eigensys, First/@eigensys];
    magIon["transitionEnergies/K"] = ReplaceDiagonal[transitionEnergies, Indeterminate];

    PrintFun["Calculating the transition wavelengths in nm ..."];
    magIon["transitionWavelengths/nm"] = 10^7/magIon["transitionEnergies/K"];
    
    energy\[CapitalDelta] = OptionValue["Energy Uncertainty in K"];
    PrintFun["Estimating the uncertainties in \[Lambda]/nm assuming a " <> ToString[energy\[CapitalDelta]] <> " K uncertainty in energies."];
    
    \[Lambda]uncertainty = energy\[CapitalDelta] * Sqrt[2.] * 10^-7 * magIon["transitionWavelengths/nm"]^2;

    PrintFun["Formatting a tabular output file ..."];
    startTime = Now;
    numEigenvecs    = Length[eigensys];
    roundedEnergies = Round[eigenEnergies, 1.];
    simpleFromTo    = Outer[{#1, #2}&, simplerStateLabels,  simplerStateLabels];
    fromTo          = Outer[{#1, #2}&, Range[numEigenvecs], Range[numEigenvecs]];
    energyPairs     = Outer[{#1, #2}&, roundedEnergies,     roundedEnergies];
    allTransitions  = {simpleFromTo,
        fromTo,
        energyPairs,
        magIon["transitionWavelengths/nm"],
        \[Lambda]uncertainty,
        magIon["AMD/n^3/s^-1"],
        magIon["radLifetimes*n^3/s"]
    };
    allTransitions = (Flatten/@Transpose[Flatten[#,1]&/@allTransitions]);
    allTransitions = Select[allTransitions, #[[3]]!=#[[4]]&];
    allTransitions = Select[allTransitions, #[[10]]>0&];
    allTransitions = Transpose[allTransitions];
    endTime = Now;
    PrintFun[">> Formatting the tabular output took ", QuantityMagnitude[endTime-startTime,"Seconds"]," seconds."];

    (*round things up*)
    PrintFun["Rounding wavelengths according to estimated uncertainties ..."];
    starTime = Now;
    {roundedWaves, roundedDeltas} = Transpose[MapThread[RoundValueWithUncertainty,{allTransitions[[7]],allTransitions[[8]]}]];
    allTransitions[[7]]          = roundedWaves;
    allTransitions[[8]]          = roundedDeltas;
    endTime = Now;
    PrintFun[">> Rounding the wavelengths took ", QuantityMagnitude[endTime-startTime,"Seconds"]," seconds."];

    PrintFun["Rounding lifetimes and transition rates to three significant figures ..."];
    startTime = Now;
    allTransitions[[9]]          = RoundToSignificantFigures[#,3]&/@(allTransitions[[9]]);
    allTransitions[[10]]         = RoundToSignificantFigures[#,3]&/@(allTransitions[[10]]);
    finalTable                   = Transpose[allTransitions];
    finalTable                   = Prepend[finalTable,header];
    endTime = Now;
    PrintFun[">> Rounding the lifetimes and transition rates took ", QuantityMagnitude[endTime-startTime,"Seconds"]," seconds."];

    (* tabular output *)
    basename     = ln <> " in " <> host <> " - ql - " <> "M1 - tabular.zip";
    exportFname  = FileNameJoin[{DirectoryName[fname], basename}];
    PrintFun["Exporting tabular data to " <> exportFname <> " ..."];
    starTime = Now;
    exportKey    = StringReplace[basename,".zip"->".m"];
    Export[exportFname, <|exportKey -> finalTable|>];
    endTime = Now;
    PrintFun[">> Exporting the tabular data took ", QuantityMagnitude[endTime-starTime,"Seconds"]," seconds."];

    (* raw data output *)
    basename       = ln <> " in " <> host <> " - ql - " <> "M1 - raw.zip";
    rawexportFname = FileNameJoin[{DirectoryName[fname], basename}];
    starTime = Now;
    PrintFun["Exporting raw data as an association to "<>exportFname<>" ..."];
    rawexportKey   = StringReplace[basename, ".zip"->".m"];
    Export[rawexportFname, <|rawexportKey -> magIon|>];
    endTime = Now;
    PrintFun[">> Exporting the raw data took ", QuantityMagnitude[endTime-starTime,"Seconds"]," seconds."];

    (* csv output *)
    PrintFun["Formatting and exporting a CSV output..."];
    startTime = Now;
    csvOut = Table[
        StringJoin[Riffle[ToString[#,CForm]&/@finalTable[[i]],","]],
    {i,1,Length[finalTable]}
    ];
    csvOut      = StringJoin[Riffle[csvOut, "\n"]];
    endTime = Now;
    PrintFun[">> Formatting the CSV output took ", QuantityMagnitude[endTime-startTime,"Seconds"]," seconds."];
    basename    = ln <> " in " <> host <> " - ql - " <> "M1.csv";
    exportFname = FileNameJoin[{DirectoryName[fname], basename}];
    PrintFun["Exporting csv data to "<>exportFname<>" ..."];
    starTime = Now;
    Export[exportFname, csvOut, "Text"];
    endTime = Now;
    PrintFun[">> Exporting the csv data took ", QuantityMagnitude[endTime-starTime,"Seconds"]," seconds."];

    If[OptionValue["Make Notebook"],
    (
        PrintFun["Creating a notebook with a Manipulate to select a wavelength interval and a lifetime power of ten ..."];
        startTime      = Now;
        finalTable     = Rest[finalTable];
        finalTable     = SortBy[finalTable,#[[7]]&];
        opticalTable   = Select[finalTable,\[Lambda]min<=#[[7]]<=\[Lambda]max&];
        pows           = Sort[DeleteDuplicates[(MantissaExponent[#[[9]]][[2]]-1)&/@opticalTable]];
        
        man            = Manipulate[
        (
            {\[Lambda]min,\[Lambda]max} = \[Lambda]int;
            table = Select[opticalTable,And[(\[Lambda]min<=#[[7]]<=\[Lambda]max),
                        (MantissaExponent[#[[9]]][[2]]-1)==log10\[Tau]]&];
            tab   = TableForm[table,TableHeadings->{None,header}];
            Column[{{"\[Lambda]min="<>ToString[\[Lambda]min]<>" nm","\[Lambda]max="<>ToString[\[Lambda]max]<>" nm",log10\[Tau]},tab}]
        ),
          {{\[Lambda]int,\[Lambda]Range,"\[Lambda] interval"},
              \[Lambda]Range[[1]],
              \[Lambda]Range[[2]],
              50,
              ControlType->IntervalSlider
          },
          {{log10\[Tau],pows[[-1]]},
              pows
          },
          TrackedSymbols  :> {\[Lambda]int,log10\[Tau]},
          SaveDefinitions -> True
        ];

        nb = CreateDocument[{
            TextCell[Style[DisplayForm[RowBox[{"Magnetic Dipole Transitions", "\n", SuperscriptBox[host<>":"<>ln,"3+"],"(",SuperscriptBox["f",numE],")"}]]],"Title",TextAlignment->Center],
            (* TextCell["Magnetic Dipole Transition Lifetimes","Section",TextAlignment->Center], *)
            TextCell[man,"Output",TextAlignment->Center]
          },
          WindowSelected -> True,
          WindowTitle    -> "M1 - "<>ln<>" in "<>host,
          WindowSize     -> {1600,800}
        ];
        reloadDataTemplate = StringTemplate@If[EvenQ[numE],
          "lnStructure = Import[FileNameJoin[{NotebookDirectory[],\"`ln` in `host` - ql.mx\"}]];",
          "lnStructure = Import[FileNameJoin[{NotebookDirectory[],\"`ln` in `host` - ql.mx\"}]];
eigenEnergies = lnStructure[\"energies/K\"][[;;;;2]];
coarseStateLabels = lnStructure[\"coarseStateLabels\"][[;;;;2]];
eigensys = lnStructure[\"eigenSystem\"][[;;;;2]];
truncatedStates = lnStructure[\"truncatedStates\"][[;;;;2]];
Clear[lnStructure];" 
        ];
        SelectionMove[nb, After, Notebook];
        NotebookWrite[nb, Cell["Reload Data", "Section", TextAlignment -> Center]];
        If[OddQ[numE],
          NotebookWrite[nb, Cell["NOTE: the transition data elides Kramers' degeneracy, only one state from each pair is included.","Text"]];
        ];
        magInfoBlurb=Cell[TextData[{
 "The entire calculated data for magnetic dipole transitions may be loaded in \
the line below.\nIt loads an association with the following keys:\n\n",
 StyleBox["\[OpenCurlyDoubleQuote]AMD/n^3/s^-1\[CloseCurlyDoubleQuote]",
  FontWeight->"Bold"],
 " : A 2D array where the [[i,j]] element gives the transition rate (in 1/s) \
(divided by the cube of the refractive index) from initial state i to final \
state j.\n",
 StyleBox["\[OpenCurlyDoubleQuote]fMD/n\[CloseCurlyDoubleQuote]",
  FontWeight->"Bold"],
 ": A list with the [[i]] element giving the magnetic dipole oscillator \
strength (divided by the refractive index) between the ground state and state \
i.\n",
 StyleBox["\[OpenCurlyDoubleQuote]radLifetimes*n^3/s\[CloseCurlyDoubleQuote]",
  
  FontWeight->"Bold"],
 ": A 2D array where the [[i,j]] element gives the radiative lifetime (in s) \
(multiplied by the cube of the refractive index) between state i and state j.\
\n",
 StyleBox["\[OpenCurlyDoubleQuote]transitionEnergies/K\[CloseCurlyDoubleQuote]\
:",
  FontWeight->"Bold"],
 " A 2D array where the [[i,j]] element gives the energy difference (in Kayser) between \
state i and state j.\n",
 StyleBox["\[OpenCurlyDoubleQuote]transitionWavelengths/nm\
\[CloseCurlyDoubleQuote]",
  FontWeight->"Bold"],
 ": A 2D array where the [[i,j]] element gives the vacuum transition wavelength (in nm) \
between state i and state j.\n",
"\n",
"In addition, a comma-separated version of this data is saved in this same directory.\n"
}], "Text", TextAlignment->Center];
        infoBlurb = Cell[TextData[{
 "The entire calculated electronic structure is contained in the file that \
may be loaded in the line below.\nIt loads an association with the following \
keys:\n\n",
 StyleBox["\[OpenCurlyDoubleQuote]energies/K\[CloseCurlyDoubleQuote]",
  FontWeight->"Bold"],
 " : A list with all the calculated energies (in Kayser), explicitly including Kramers\
\[CloseCurlyQuote] degeneracy.\n",
 StyleBox["\[OpenCurlyDoubleQuote]lnSymbol\[CloseCurlyDoubleQuote]",
  FontWeight->"Bold"],
 ": The symbol representing the current lanthanide.\n",
 StyleBox["\[OpenCurlyDoubleQuote]coarseStateLabels\[CloseCurlyDoubleQuote]",
  FontWeight->"Bold"],
 ": A list with the largest LSJ contribution to the eigenvectors.\n",
 StyleBox["\[OpenCurlyDoubleQuote]eigenSystem\[CloseCurlyDoubleQuote]:",
  FontWeight->"Bold"],
 " A list of lists, the first element being the energy and the second the \
eigenvector in the given basis.\n",
 StyleBox["\[OpenCurlyDoubleQuote]truncatedStates\[CloseCurlyDoubleQuote]",
  FontWeight->"Bold"],
 ": A truncated version of the eigenvectors.\n",
 StyleBox["\[OpenCurlyDoubleQuote]params\[CloseCurlyDoubleQuote]",
  FontWeight->"Bold"],
 ": The fitted parameters of the semi-empirical Hamiltonian.\n",
 StyleBox["\[OpenCurlyDoubleQuote]fitResult\[CloseCurlyDoubleQuote]",
  FontWeight->"Bold"],
 ": An association containing the full result of the solver that produced the \
parameters used here."
}], "Text", TextAlignment->Center];
        NotebookWrite[nb, infoBlurb];
        NotebookWrite[nb, Cell[(
              reloadDataTemplate[<|"host"->host,"ln"->ln|>]
              ),"Input"]];
        NotebookWrite[nb, magInfoBlurb];
        NotebookWrite[nb, Cell[(
              "magTransitions = Import[FileNameJoin[{NotebookDirectory[],\"" <> StringSplit[rawexportFname,"/"][[-1]] <> "\"}],\""<>rawexportKey<>"\"];"
              ),"Input"]];
        NotebookWrite[nb, Cell[TextData[{"\n\n\n",
        StyleBox["This notebook was generated on " <> DateString[] <> ".",
        FontSlant->"Italic"]}], "Text", TextAlignment->Right]];
        SelectionMove[nb, Before, Notebook];
        nbFname = FileNameJoin[{DirectoryName[fname], "M1 - "<>ln<>" in "<>host<>".nb"}];
        endTime = Now;
        PrintFun[">> Creating the notebook took ", QuantityMagnitude[endTime-startTime,"Seconds"]," seconds."];
        PrintFun[">> Saving notebook to ", nbFname, " ..."];
        startTime = Now;
        NotebookSave[nb, nbFname];
        endTime = Now;
        PrintFun[">> Saving the notebook took ", QuantityMagnitude[endTime-startTime,"Seconds"]," seconds."];
        If[OptionValue["Close Notebook"],
            NotebookClose[nb];
        ];
    )
    ];

    Return[magIon];
  )
];

Options[FitLaF3] = {
  "fitOrder" -> "Pr Nd Dy Ce Sm Ho Er Tm Yb Tb Eu Gd Pm",
  "notify" -> False,
  "truncationEnergy" -> Infinity
  };  
FitLaF3[OptionsPattern[]]:=(
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
        "pathfinder-constraints-without-spinspin-and-no-truncation"
    ];
    filePrefix = "variant-" <> filePrefix;
    \[Sigma]exp = 1.0;
    trunEnergy = OptionValue["truncationEnergy"];
    Print[trunEnergy];
    fitOrder = OptionValue["fitOrder"];
    fitOrder = Position[theLanthanides, #][[1, 1]] & /@ StringSplit[fitOrder, " "];
    doThese = fitOrder;
    
    linearParams       = Join[casimirSymbols, cfSymbols, TSymbols, pseudoMagneticSymbols, marvinSymbols];
    quadraticParams    = {};
    interpolatedParams = Join[slaterSymbols,{\[Zeta]}];
    
    varModels = <||>;
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
            independentVars /. LoadLaF3Parameters[ln, "Vintage" -> "Carnall", "With Uncertainties" -> False]];
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
        If[
          Not[truncationEnergy == Infinity],
        (
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
        ),
        (
          eigenEnergies = thisSolution["energies"];
          fullRMSsquared = thisSolution["bestRMS"];
        )
        ];
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
        varModels[ln] = varModel; 
        If[OptionValue["notify"],
          Run["~/Scripts/pushover \"finished fitting " <> ln <> "\""];
        ]
    ),
    {numE0, doThese}];
    Return[{varModels, truncatedSols}];
);

LoadLiYF4Parameters::usage="Data from Cheng, Jun, Wen Jun, Chen Yonghu, Yin Min, and Duan Changkui. Crystal-Field Analyses for Trivalent Lanthanide Ions in LiYF4. Journal of Rare Earths 34, no. 10 (2016): 1048-52.";
LoadLiYF4Parameters[ln_] := (
  paramsChengLiYF4 = <|"Ce" -> <|\[Zeta] -> 630., B02 -> 354., B04 -> -1043., 
    B44 -> -1249., B06 -> -65., B46 -> -1069.|>, 
  "Pr" -> <|F0 -> 0., F2 -> 68955., F4 -> 50505., 
    F6 -> 33098., \[Zeta] -> 748., \[Alpha] -> 
      23.3, \[Beta] -> -644., \[Gamma] -> 1413., M0 -> 1.88, P2 -> 244.,
      B02 -> 512., B04 -> -1127., B44 -> -1239., B06 -> -85., 
    B46 -> -1205.|>, 
  "Nd" -> <|F0 -> 0., F2 -> 72952., F4 -> 52681., 
    F6 -> 35476., \[Zeta] -> 877., \[Alpha] -> 
      21., \[Beta] -> -579., \[Gamma] -> 1446., T2 -> 210., T3 -> 41., 
    T4 -> 74., T6 -> -293., T7 -> 321., T8 -> 205., M0 -> 1.85, 
    P2 -> 304., B02 -> 391., B04 -> -1031., B44 -> -1271., B06 -> -28.,
      B46 -> -1046.|>, 
  "Sm" -> <|F0 -> 0., F2 -> 79515., F4 -> 56766., 
    F6 -> 40078., \[Zeta] -> 1168., \[Alpha] -> 
      20.5, \[Beta] -> -616., \[Gamma] -> 1565., T2 -> 282., T3 -> 26., 
    T4 -> 71., T6 -> -257., T7 -> 314., T8 -> 328., M0 -> 2.38, 
    P2 -> 336., B02 -> 370., B04 -> -757., B44 -> -941., B06 -> -67., 
    B46 -> -895.|>, 
  "Eu" -> <|F0 -> 0., F2 -> 82573., F4 -> 59646., 
    F6 -> 43203., \[Zeta] -> 1329., \[Alpha] -> 
      21.6, \[Beta] -> -482., \[Gamma] -> 1140., T2 -> 370., T3 -> 40., 
    T4 -> 40., T6 -> -300., T7 -> 380., T8 -> 370., M0 -> 2.41, 
    P2 -> 332., B02 -> 339., B04 -> -733., B44 -> -1067., B06 -> -36., 
    B46 -> -764.|>, 
  "Tb" -> <|F0 -> 0., F2 -> 90972., F4 -> 64499., 
    F6 -> 45759., \[Zeta] -> 1702., \[Alpha] -> 
      17.6, \[Beta] -> -581., \[Gamma] -> 1792., T2 -> 330., T3 -> 40., 
    T4 -> 45., T6 -> -365., T7 -> 320., T8 -> 349., M0 -> 2.7, 
    P2 -> 482., B02 -> 413., B04 -> -867., B44 -> -1114., B06 -> -41., 
    B46 -> -736.|>, 
  "Dy" -> <|F0 -> 0, F2 -> 90421., F4 -> 63928., 
    F6 -> 46657., \[Zeta] -> 1895., \[Alpha] -> 
      17.9, \[Beta] -> -628., \[Gamma] -> 1790., T2 -> 326., T3 -> 23., 
    T4 -> 83., T6 -> -294., T7 -> 403., T8 -> 340., M0 -> 4.46, 
    P2 -> 610., B02 -> 360., B04 -> -737., B44 -> -943., B06 -> -35., 
    B46 -> -700.|>, 
  "Ho" -> <|F0 -> 0., F2 -> 93512., F4 -> 66084., 
    F6 -> 49765., \[Zeta] -> 2126., \[Alpha] -> 
      17.2, \[Beta] -> -596., \[Gamma] -> 1839., T2 -> 365., T3 -> 37., 
    T4 -> 95., T6 -> -274., T7 -> 331., T8 -> 343., M0 -> 3.92, 
    P2 -> 582., B02 -> 386., B04 -> -629., B44 -> -841., B06 -> -33., 
    B46 -> -687.|>, 
  "Er" -> <|F0 -> 0., F2 -> 97326., F4 -> 67987., 
    F6 -> 53651., \[Zeta] -> 2377., \[Alpha] -> 
      18.1, \[Beta] -> -599., \[Gamma] -> 1870., T2 -> 380., T3 -> 41., 
    T4 -> 69., T6 -> -356., T7 -> 239., T8 -> 390., M0 -> 4.41, 
    P2 -> 795., B02 -> 325., B04 -> -749., B44 -> -1014., B06 -> -19., 
    B46 -> -635.|>, 
  "Tm" -> <|F0 -> 0., T2 -> 0., F2 -> 101938., F4 -> 71553., 
    F6 -> 51359., \[Zeta] -> 2632., \[Alpha] -> 
      17.3, \[Beta] -> -665., \[Gamma] -> 1936., M0 -> 4.93, P2 -> 730.,
      B02 -> 339., B04 -> -627., B44 -> -913., B06 -> -39., 
    B46 -> -584.|>, 
  "Yb" -> <|\[Zeta] -> 2916., B02 -> 446., B04 -> -560., B44 -> -843., 
    B06 -> -23., B46 -> -512.|>
    |>;
  Return[paramsChengLiYF4[ln]];
);

Options[FitLiYF4] = {
  "fitOrder" -> "Er Nd Eu Ho Sm Pr Tm Yb Ce Tb Dy Gd Pm",
  "notify" -> False};
FitLiYF4[OptionsPattern[]]:=(
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
    filePrefix = "variant-" <> filePrefix;
    (* This according to Cheng probably following Gorller-Walrand *)
    magSimplifier = {
        M2 -> 56/100 M0,
        M4 -> 38/100 M0,
        P4 -> 3/4 P2,
        P6 -> 1/2 P2};
    fitOrder = "Er Nd Eu Ho Sm Pr Tm Yb Ce Tb Dy Gd Pm";
    fitOrder = Position[theLanthanides, #][[1, 1]] & /@ StringSplit[fitOrder, " "];
    doThese = fitOrder;
    
    linearParams       = Join[casimirSymbols, cfSymbols, TSymbols, pseudoMagneticSymbols, marvinSymbols];
    quadraticParams    = {};
    interpolatedParams = Join[slaterSymbols,{\[Zeta]}];
    
    varModels = <||>;
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
        
        If[OptionValue["notify"],
            Run["~/Scripts/pushover \"finished fitting " <> ln <> "\""];
        ];
        varModels[ln] = varModel;
    ),
    {numE0, doThese}];
    Return[{varModels, truncatedSols}];
    );


FitLaF3RestrictedOrthogonal[]:=(
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
        "restricted-orthogonal-pathfinder-constraints-with-spinspin-and-no-truncation",
        "restricted-orthogonal-pathfinder-constraints-without-spinspin-and-no-truncation"
    ];
    \[Sigma]exp = 1.0;
    trunEnergy  = Infinity;
    fitOrder = "Pr Nd Dy Ce Sm Ho Er Tm Yb Tb Eu Gd Pm";
    fitOrder = "Pr Nd Dy Ce Sm Ho Er Tm Yb Tb Eu Gd";
    fitOrder = Position[theLanthanides, #][[1, 1]] & /@ StringSplit[fitOrder, " "];
    doThese = fitOrder;
    
    linearParams       = Join[{\[Alpha]p,\[Beta]p,\[Gamma]p}, cfSymbols, TSymbols, pseudoMagneticSymbols, marvinSymbols];
    quadraticParams    = {};
    interpolatedParams = Join[{E1p,E2p,E3p},{\[Zeta]}];
    
    varModels = <||>;
    If[Not[restart],
    varModel = AssociationThread[{B02, B04, B06, B22, B24, B26, B44, B46, B66, E1p,
        E2p, E3p, M0, P2, T2p, T3, T4, T6, T7, 
        T8, \[Alpha]p, \[Beta]p, \[Gamma]p, \[Zeta]}, 
        ConstantArray[<|"data" -> {}, "function" -> Null|>, 24]];
    truncatedSols = <||>;
    ];
    
    Do[
    (
        numE = numE0;
        ln   = theLanthanides[[numE]];
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
        
        mostyOvalues = Association@Table[(param -> varModel[param]["function"][numE]), {param, Keys[varModel]}];
        nonOvalues   = FromMostlyOrthogonalToNonOrthogonal[mostyOvalues, numE];
        
        Do[
          (
            constraints[fixedCon] = nonOvalues[fixedCon] 
          ), 
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
            independentVars /. LoadLaF3Parameters[ln, "Vintage" -> "Carnall", "With Uncertainties" -> False]];
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
        mostlyOparams = FromNonOrthogonalToMostlyOrthogonal[bestParams, numE];
        Do[
        (
            If[Not[MemberQ[Keys[varModel], pVar]],
            Continue[];
            ];
            varModel[pVar]["data"] = 
            DeleteDuplicates@
            Append[varModel[pVar]["data"], {numE, mostlyOparams[pVar]}];
            numData = Length[varModel[pVar]["data"]];
            varModel[pVar]["data"] = Select[varModel[pVar]["data"], FreeQ[#, Missing] &];
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
        {pVar, Keys[varModel]}
        ];
        varModels[ln] = varModel; 
        Run["~/Scripts/pushover \"finished fitting " <> ln <> "\""];
    ),
    {numE0, doThese}];
    Return[{varModels, truncatedSols}];
);

FitLaF3MostlyOrthogonal::usage="FitLaF3MostlyOrthogonal[] fits the LaF3 data using the mostly-orthogonal operator basis. 
The function admits two options: \"FitOrder\" and \"MaxIterations\". The first one is a string with the order of the lanthanides to be fitted, and the second one is the maximum number of iterations for the fitting routine. The default values are \"Pr Nd Dy Ce Sm Ho Er Tm Yb Tb Eu Gd Pm\" and 1000, respectively.
The results of fitting every ion are saved to files in the directory \"mostly-orthogonal-pathfinder-constraints-with-spinspin-and-no-truncation\" or \"mostly-orthogonal-pathfinder-constraints-without-spinspin-and-no-truncation\". The results are saved in the form of an association with the following keys:
- simplifier
- excludeDataIndices
- startValues
- freeIonSymbols
- truncationEnergy
- numE
- expData
- problemVars
- maxIterations
- hamDim
- constraints
- allVars
- freeBies
- compiledIntermediateFname
- presentDataIndices
- degreesOfFreedom
- solWithUncertainty
- truncatedDim
- fittedLevels
- actualSteps
- paramSols
- rmsHistory
- Appendix
- timeTaken/s
- bestParams
- bestParamsWithConstraints
- bestRMSwithTruncation
- energiesWithTruncation
- excludedDataIndices
- bestRMSwithFullDiagonalization
- energiesWithFullDiagonalization
The function returns a list with two elements. The first element being an association with containing the regression models that were built to the several model parameters during fitting, and the second element being an association with the results (as associations with the keys listed above) of the fitting for each lanthanide ion.";
Options[FitLaF3MostlyOrthogonal] = {
  "FitOrder" -> "Pr Nd Dy Ce Sm Ho Er Tm Yb Tb Eu Gd Pm",
  "MaxIterations" -> 1000};
FitLaF3MostlyOrthogonal[OptionsPattern[]]:=Module[
  {
    withSpinSin,
    restart,
    saveEigenvectors,
    filePrefix,
    \[Sigma]exp,
    trunEnergy,
    fitOrder,
    doThese,
    linearParams,
    quadraticParams,
    interpolatedParams,
    varModels,
    varModel,
    truncatedSols,
    numE0,
    numE,
    ln,
    pmSol,
    pVar,
    bestP,
    constrainedMsandPs,
    pVars,
    constraints,
    fixedConstraints,
    carnalKey,
    expData,
    presentDataIndices,
    excludedDataIndices,
    independentVars,
    startValues,
    thisSolution,
    bestParams,
    params,
    eigenEnergies,
    fullDiffs,
    numFreeVars,
    numData,
    fullRMSsquared,
    order,
    b
},
(
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
        "mostly-orthogonal-pathfinder-constraints-with-spinspin-and-no-truncation",
        "mostly-orthogonal-pathfinder-constraints-without-spinspin-and-no-truncation"
    ];
    \[Sigma]exp = 1.0;
    trunEnergy  = Infinity;
    fitOrder = OptionValue["FitOrder"];
    fitOrder = Position[theLanthanides, #][[1, 1]] & /@ StringSplit[fitOrder, " "];
    doThese  = fitOrder;
    
    linearParams       = Join[{\[Alpha]p, \[Beta]p, \[Gamma]p}, cfSymbols, TSymbols, pseudoMagneticSymbols, marvinSymbols];
    quadraticParams    = {};
    interpolatedParams = Join[{E1p, E2p, E3p}, {\[Zeta]}];
    
    varModels = <||>;
    If[Not[restart],
    varModel = AssociationThread[{
        B02, B04, B06, B22, B24, B26, B44, B46, B66, 
        E1p, E2p, E3p, 
        M0, P2, 
        T2p, T3, T4, T6, T7, T8,
        \[Alpha]p, \[Beta]p, \[Gamma]p,
        \[Zeta]}, 
        ConstantArray[<|"data" -> {}, "function" -> Null|>, 24]];
    truncatedSols = <||>;
    ];
    
    Do[
    (
        numE = numE0;
        ln   = theLanthanides[[numE]];
        PrintAndDefine[Style["Working on "<>ln<>" ...", Red]];
        (* remember to modify this one, too!*)
        If[ln == "Pm",
            (
                pmSol = <||>;
                pmSol["hamDim"] = Binomial[14, numE];
                pmSol["numE"] = numE;
                pmSol["problemVars"] = variedSymbolsMO["Nd"];
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
                {pVar, variedSymbolsMO["Nd"]}];
                pmSol["bestParamsWithConstraints"][\[Epsilon]] = 0;
                constrainedMsandPs = ("MagneticSimplifier" /. 
                    Options[MostlyOrthogonalFit]) /. 
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
        pVars = variedSymbolsMO[ln];
        constraints = Association@caseConstraintsMO[ln];
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
            independentVars /. FromNonOrthogonalToMostlyOrthogonal[LoadLaF3Parameters[ln, "Vintage" -> "Carnall", "With Uncertainties" -> False], numE]];
        thisSolution = MostlyOrthogonalFit[numE,
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
                Options[MostlyOrthogonalFit]) /. (\[Sigma]SS -> _) :> (\[Sigma]SS ->
                    withSpinSin)),
            "MaxIterations" -> OptionValue["MaxIterations"],
            "SaveToLog" -> False,
            "MaxHistory" -> 1000,
            "FilePrefix" -> filePrefix <> "/" <> ln <> "-" <> ToString[trunEnergy],
            "TruncationEnergy" -> trunEnergy,
            "PrintFun" -> PrintAndDefine];
        thisSolution["bestParams"] = Append[thisSolution["bestParams"], (nE -> numE)];
        thisSolution["bestParamsWithConstraints"][nE] = numE;
        If[(2 <= numE <= 12),
            (constrainedMsandPs = ("MagneticSimplifier" /. 
                Options[MostlyOrthogonalFit]) /. thisSolution["bestParamsWithConstraints"];
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
        numFreeVars = Length[thisSolution["problemVars"]];
        numData = Length[thisSolution["presentDataIndices"]];
        Which[
          OddQ[numE],
          (
            Print["ATTENTION: Odd number of electrons, disregarding Kramers' degeneracy for calculation of sigma ..."];
            fullRMSsquared = Sqrt[Total[fullDiffs[[;;;;2]]^2]/(numData/2 - numFreeVars - 1)];
            thisSolution["bestRMS"] = Sqrt[
              thisSolution["bestRMS"]^2 * 
              thisSolution["degreesOfFreedom"] /
              (numData - 2*numFreeVars - 2)];
            thisSolution["degreesOfFreedom"] = numData/2 - numFreeVars - 1;
          ),
          True,
          (
            fullRMSsquared = Sqrt[Total[fullDiffs^2]/(thisSolution["degreesOfFreedom"])]
          )  
        ];
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
        varModels[ln] = varModel; 
        Run["~/Scripts/pushover \"finished fitting " <> ln <> "\""];
    ),
    {numE0, doThese}];
    Return[{varModels, truncatedSols}];
)
];

FromNonOrthogonalToMostlyOrthogonal::usage="FromNonOrthogonalToMostlyOrthogonal[nonOparams, nE] takes  an association of parameters corresponding to the non-orthogonal version of the Hamiltonian and returns an association for the parameters in the mostly-orthogonal representation (E0p, E1p, E2p, \[Alpha]p, \[Beta]p, \[Gamma]p). 
Mostly orthogonal is meant here in the sense that the Pk and Mk parameters are not orthogonal within themselves although they are to the other parameters.
The function takes the following arguments: 
    - nonOparams: an association with keys equal to non-orthogonal parameter symbols and associated values. 
    - nE: the number of electrons in the corresponding configuration.
The function returns and association for the parameters that can be reasonably computed given the input.    
";
FromNonOrthogonalToMostlyOrthogonal[nonOparams_, nE_] := Module[
  {paramChanger, mostlyOparams, nnE = nE},
  (
    paramChanger = <|
      E0p -> F0 - (4 F2)/195 - (2 F4)/143 - (100 F6)/5577 - (24 \[Alpha])/13 - \[Beta]/13 - (6 \[Gamma])/65, 
      E1p -> (14 F2)/405 + (7 F4)/297 + (350 F6)/11583 + (4 \[Alpha])/5 + \[Beta]/30 + \[Gamma]/25, 
      E2p -> F2/2025 - F4/3267 + (175 F6)/1656369, 
      E3p -> F2/135 + (2 F4)/1089 - (175 F6)/42471 - T2/(35 Sqrt[2]) + (nnE T2)/(70 Sqrt[2]) - (2 \[Alpha])/5, 
      \[Alpha]p -> (4 \[Alpha])/5, 
      \[Beta]p  -> -4 \[Alpha] - \[Beta]/6, 
      \[Gamma]p -> (8 \[Alpha])/5 + \[Beta]/15 + (2 \[Gamma])/25,
      B02 -> B02,
      B04 -> B04,
      B06 -> B06,
      B12 -> B12,
      B14 -> B14,
      B16 -> B16,
      B22 -> B22,
      B24 -> B24,
      B26 -> B26,
      B34 -> B34,
      B36 -> B36,
      B44 -> B44,
      B46 -> B46,
      B56 -> B56,
      B66 -> B66,
      S12 -> S12,
      S14 -> S14,
      S16 -> S16,
      S22 -> S22,
      S24 -> S24,
      S26 -> S26,
      S34 -> S34,
      S36 -> S36,
      S44 -> S44,
      S46 -> S46,
      S56 -> S56,
      S66 -> S66,
      \[Zeta] -> \[Zeta],
      T2p -> T2,
      T3 -> T3,
      T4 -> T4,
      T6 -> T6,
      T7 -> T7,
      T8 -> T8,
      M0 -> M0,
      M2 -> M2,
      M4 -> M4,
      P2 -> P2,
      P4 -> P4,
      P6 -> P6
      |>;
      paramChanger = Select[paramChanger, ContainsAll[Keys[nonOparams], Variables[#]] &];
    mostlyOparams = Expand /@ (paramChanger /. nonOparams);
    Return[mostlyOparams];
)
];

FromMostlyOrthogonalToNonOrthogonal::usage="FromMostlyOrthogonalToNonOrthogonal[mostlyOparams_, nE] takes  an association of parameters corresponding to the mostly-orthogonal version of the Hamiltonian and returns an association for the parameters in the non-orthogonal representation (F0, F2, F4, F6, \[Alpha], \[Beta], \[Gamma], ...). Mostly orthogonal is meant here in the sense that the Pk and Mk parameters are not orthogonal within themselves although they are to the other parameters.
The function takes the following arguments: 
    - mostlyOparams: an association with keys equal to mostly-orthogonal parameter symbols and associated values. 
    - nE: the number of electrons in the corresponding configuration.
The function returns and association for the parameters that can be reasonably computed given the input. 
";
FromMostlyOrthogonalToNonOrthogonal[mostlyOparams_, nE_]:=Module[
  {paramChanger, nonOparams, nnE = nE},
  (
    paramChanger = <|
        F0 -> 1/91 (54 E1p+91 E0p+78 \[Gamma]p),
        F2 -> 15/392 (140 E1p+20020 E2p+1540 E3p+770 \[Alpha]p-70 \[Gamma]p+22 Sqrt[2] T2p-11 Sqrt[2] nnE T2p),
        F4 -> 99/490 (70 E1p-9100 E2p+280 E3p+140 \[Alpha]p-35 \[Gamma]p+4 Sqrt[2] T2p-2 Sqrt[2] nnE T2p),
        F6 -> (5577 (20 E1p+700 E2p-140 E3p-70 \[Alpha]p-10 \[Gamma]p-2 Sqrt[2] T2p + Sqrt[2] nnE T2p))/7000,
        \[Zeta]  -> \[Zeta],
        \[Alpha] -> (5 \[Alpha]p)/4,
        \[Beta]  -> -6  (5 \[Alpha]p + \[Beta]p),
        \[Gamma] -> 5/2 (2 \[Beta]p + 5 \[Gamma]p),
        T2 -> T2p,
        T3 -> T3,
        T4 -> T4,
        T6 -> T6,
        T7 -> T7,
        T8 -> T8,
        M0 -> M0,
        M2 -> M2,
        M4 -> M4,
        P2 -> P2,
        P4 -> P4,
        P6 -> P6,
        B02 -> B02,
        B04 -> B04,
        B06 -> B06,
        B12 -> B12,
        B14 -> B14,
        B16 -> B16,
        B22 -> B22,
        B24 -> B24,
        B26 -> B26,
        B34 -> B34,
        B36 -> B36,
        B44 -> B44,
        B46 -> B46,
        B56 -> B56,
        B66 -> B66,
        S12 -> S12,
        S14 -> S14,
        S16 -> S16,
        S22 -> S22,
        S24 -> S24,
        S26 -> S26,
        S34 -> S34,
        S36 -> S36,
        S44 -> S44,
        S46 -> S46,
        S56 -> S56,
        S66 -> S66
      |>;
      paramChanger = Select[paramChanger,ContainsAll[Keys[mostlyOparams],Variables[#]]&];
    nonOparams = Expand/@(paramChanger/.mostlyOparams);
    Return[nonOparams];
  )
];

Options[SimulatedAnnealingSolver] = {
  "corral" -> 0.01,
  "maxIterations" -> 200,
  "superCycles" -> 1,
  "maxTimeInMin" -> 5
};
SimulatedAnnealingSolver[numE_,OptionsPattern[]]:=(
  (* corral controls how much the parameters are constrained based on the original solution *)
  corral        = OptionValue["corral"]; 
  maxIterations =  OptionValue["maxIterations"];
  (* when using the simulated annealing method this number of attempts are tried out*)
  superCycles   = OptionValue["superCycles"];
  (* each of those attempts is time constrained *)
  maxTimeInMin  = OptionValue["maxTimeInMin"]; 
  (* such that that part of the solver takes about superCycles * maxTimeInMin minutes *)
  ln = theLanthanides[[numE]];
  
  Print["Loading previously found solution ..."];
  truncatedSols = Import["/Users/juan/ZiaLab/Codebase/qlanth/LaF3-final-fit-mostly-orthogonal.m"]["truncatedSols"];
  {compiledHam, intermediateBasis} = Import[truncatedSols[numE]["compiledIntermediateFname"]];
  sol = truncatedSols[numE]["bestParamsWithConstraints"];
  problemVars = Append[truncatedSols[numE]["problemVars"], \[Epsilon]];
  expData     = First/@truncatedSols[numE]["expData"];
  presentDataIndices  = truncatedSols[numE]["presentDataIndices"];
  excludedDataIndices = truncatedSols[numE]["excludedDataIndices"];
  relevantDataIndices = Complement[presentDataIndices,excludedDataIndices];
  degsOfFree          = If[OddQ[numE],
        Length[relevantDataIndices]/2 - Length[problemVars],
        Length[relevantDataIndices]   - Length[problemVars]
  ];
  Which[
    ln == "Er",
    CostFun[
      B02num_?NumericQ,B04num_?NumericQ,B06num_?NumericQ,B22num_?NumericQ,
      B24num_?NumericQ,B26num_?NumericQ,B44num_?NumericQ,B46num_?NumericQ,B66num_?NumericQ,
      E1pnum_?NumericQ,E2pnum_?NumericQ,E3pnum_?NumericQ,
      M0num_?NumericQ,P2num_?NumericQ,
      T2pnum_?NumericQ,T3num_?NumericQ,T4num_?NumericQ,T6num_?NumericQ,T7num_?NumericQ,T8num_?NumericQ,
      \[Alpha]pnum_?NumericQ,\[Beta]pnum_?NumericQ,\[Gamma]pnum_?NumericQ,
      \[Zeta]num_?NumericQ,\[Epsilon]num_?NumericQ]:=(
    eigenValues = Sort@Eigenvalues@compiledHam[B02num,B04num,B06num,B22num,B24num,B26num,B44num,B46num,B66num,E1pnum,E2pnum,E3pnum,M0num,P2num,T2pnum,T3num,T4num,T6num,T7num,T8num,\[Alpha]pnum,\[Beta]pnum,\[Gamma]pnum,\[Zeta]num];
      eigenValues = (eigenValues-eigenValues[[1]]+\[Epsilon]num)[[;;Length[expData]]];
      diffs       = (eigenValues-expData)[[relevantDataIndices]];
      theCost     = Total[diffs^2];
      Return[theCost]
    );,
    ln == "Nd",
    CostFun[
      B02num_?NumericQ,B04num_?NumericQ,B06num_?NumericQ,B22num_?NumericQ,
      B24num_?NumericQ,B26num_?NumericQ,B44num_?NumericQ,B46num_?NumericQ,B66num_?NumericQ,
      E1pnum_?NumericQ,E2pnum_?NumericQ,E3pnum_?NumericQ,
      M0num_?NumericQ,P2num_?NumericQ,
      T2pnum_?NumericQ,T3num_?NumericQ,T4num_?NumericQ,T6num_?NumericQ,T7num_?NumericQ,T8num_?NumericQ,
      \[Alpha]pnum_?NumericQ,\[Beta]pnum_?NumericQ,\[Gamma]pnum_?NumericQ,
      \[Zeta]num_?NumericQ,\[Epsilon]num_?NumericQ]:=(
    eigenValues = Sort@Eigenvalues@compiledHam[B02num,B04num,B06num,B22num,B24num,B26num,B44num,B46num,B66num,E1pnum,E2pnum,E3pnum,M0num,P2num,T2pnum,T3num,T4num,T6num,T7num,T8num,\[Alpha]pnum,\[Beta]pnum,\[Gamma]pnum,\[Zeta]num];
      eigenValues = (eigenValues-eigenValues[[1]]+\[Epsilon]num)[[;;Length[expData]]];
      diffs       = (eigenValues-expData)[[relevantDataIndices]];
      theCost     = Total[diffs^2];
      Return[theCost]
    ); 
  ];

  Options[ProgressNotebook] = {"Basic" -> True};
  ProgressNotebook[
    OptionsPattern[]] := (nb = 
      Which[OptionValue["Basic"], 
      CreateDocument[({Dynamic[
          TextCell[
            If[Length[paramSols] > 0, 
            TableForm[
              Prepend[Transpose[{stringPartialVars, 
                paramSols[[-1]]}], {"RMS", rmsHistory[[-1]]}]], ""], 
            "Output"], 
          TrackedSymbols :> {paramSols, stringPartialVars},
          UpdateInterval -> 0.5
          ]}), WindowSize -> {600, 1000}, WindowSelected -> True, 
        TextAlignment -> Center, WindowTitle -> "Solver Progress"], 
      True, CreateDocument[({"", Dynamic[
          Framed[progressMessage], UpdateInterval -> 0.5], 
          Dynamic[GraphicsColumn[{ListPlot[rmsHistory, 
              PlotMarkers -> "OpenMarkers", Frame -> True, 
              FrameLabel -> {"Iteration", "RMS"}, ImageSize -> 800, 
              AspectRatio -> 1/3, FrameStyle -> Directive[Thick, 15], 
              PlotLabel -> 
              If[Length[rmsHistory] != 0, rmsHistory[[-1]], ""]], 
            ListPlot[(#/#[[1]]) & /@ Transpose[paramSols], 
              Joined -> True, PlotRange -> {All, {-5, 5}}, 
              Frame -> True, ImageSize -> 800, AspectRatio -> 1, 
              FrameStyle -> Directive[Thick, 15], 
              FrameLabel -> {"Iteration", "Params"}]}], 
          TrackedSymbols :> {rmsHistory, paramSols},
          UpdateInterval -> 0.5], 
          Dynamic[TextCell[
            If[Length[paramSols] > 0, 
            TableForm[Transpose[{stringPartialVars, paramSols[[-1]]}]],
              ""], "Output"], 
          TrackedSymbols :> {paramSols, stringPartialVars},
          UpdateInterval -> 0.5
          ]}), WindowSize -> {600, 1000}, WindowSelected -> True, 
        TextAlignment -> Center, WindowTitle -> "Solver Progress"]];
    Return[nb];
    );
   
  AddToList[list_, element_, maxSize_, addOnlyNew_ : True] := 
    Module[{tempList = If[
        addOnlyNew, If[Length[list] == 0,
        {element}, 
        If[element != list[[-1]], Append[list, element], list]], 
        Append[list, element]]}, 
    If[Length[tempList] > maxSize, 
      Drop[tempList, Length[tempList] - maxSize], tempList]];

  runningInteractive = True;
  stringPartialVars = ToString /@ problemVars;
  paramSols = {};
  rmsHistory = {};
  openNotebooks = 
    If[runningInteractive, ("WindowTitle" /. NotebookInformation[#]) & /@
      Notebooks[], {}];
  If[Not[MemberQ[openNotebooks, "Solver Progress"]], 
    ProgressNotebook["Basic" -> False]];

  bestSigma = Infinity;
  bestSol = {};
  atlas = {};

  maxHistory = Infinity;
  steps = 0;
  paramSols = {};
  rmsHistory = {};
  maxTimeInSec = maxTimeInMin*60;
  seeds = {};
  Do[(
    corrals = Table[
      (
      centerValue = sol[var];
      varName = ToExpression[(ToString[var] <> "opt")];
      If[centerValue < 0,
        (((1 + corral)*centerValue) <= 
          varName <= ((1 - corral) centerValue)),
        (((1 + corral)*centerValue) >= 
          varName >= ((1 - corral) centerValue))
        ]
      ), {var, problemVars}];
    corrals[[-1]] = (-30. <= \[Epsilon]opt <= 30.);
    
    startValues = Table[(
      centerValue = sol[var];
      If[centerValue > 0,
        RandomReal[{centerValue*(1 - corral), centerValue*(1 + corral)}],
        RandomReal[{centerValue*(1 + corral), 
          centerValue*(1 - corral)}]
        ]), {var, problemVars}];
    startValues[[-1]] = 0.;
    seed = RandomInteger[{1, 10000}];
    seeds = Append[seeds, seed];
    
    Print["Running FindMinimum with the simulated annealing method ..."];
    minSol = 
    TimeConstrained[
      NMinimize[
      Join[{CostFun[B02opt, B04opt, B06opt, B22opt, B24opt, B26opt, 
          B44opt, B46opt, B66opt, E1popt, E2popt, E3popt, M0opt, P2opt, 
          T2popt, T3opt, T4opt, T6opt, T7opt, 
          T8opt, \[Alpha]popt, \[Beta]popt, \[Gamma]popt, \[Zeta]opt, \
  \[Epsilon]opt]}, corrals],
      {B02opt, B04opt, B06opt, B22opt, B24opt, B26opt, B44opt, B46opt, 
        B66opt, E1popt, E2popt, E3popt, M0opt, P2opt, T2popt, T3opt, 
        T4opt, T6opt, T7opt, 
        T8opt, \[Alpha]popt, \[Beta]popt, \[Gamma]popt, \[Zeta]opt, \
  \[Epsilon]opt},
      Method -> {"SimulatedAnnealing", "RandomSeed" -> seed},
      MaxIterations -> maxIterations,
      StepMonitor :>
        (
        steps += 1;
        currentSqSum = 1/2*theCost;
        currentRMS = Sqrt[currentSqSum/degsOfFree];
        currentParams = {B02opt, B04opt, B06opt, B22opt, B24opt, 
          B26opt, B44opt, B46opt, B66opt, E1popt, E2popt, E3popt, 
          M0opt, P2opt, T2popt, T3opt, T4opt, T6opt, T7opt, 
          T8opt, \[Alpha]popt, \[Beta]popt, \[Gamma]popt, \[Zeta]opt, \
  \[Epsilon]opt};
        AbortProtect[(
          
          paramSols = 
            AddToList[paramSols, currentParams, maxHistory, True];
          
          rmsHistory = 
            AddToList[rmsHistory, currentRMS, maxHistory, True];
          atlas = Append[atlas, {currentRMS, paramSols}];
          )];
        If[currentRMS < bestSigma,
          AbortProtect[(
            bestSol = currentParams;
            bestSigma = currentRMS
            )]
          ]
        )
      ],
      maxTimeInSec,
      {bestSigma, AssociationThread[problemVars, bestSol]}
      ];
    )
  , {supercycle, 1, superCycles}];

  startValues = Values@minSol[[2]];
  Print["Running FindMinimum with automatic method, polishing the best solution found through simulated annealing ..."];
  paramSols = {};
  rmsHistory = {};  
  refinedSol  = TimeConstrained[
  FindMinimum[
    CostFun[B02opt,B04opt,B06opt,B22opt,B24opt,B26opt,B44opt,B46opt,B66opt,
            E1popt,E2popt,E3popt,
            M0opt,P2opt,
            T2popt,T3opt,T4opt,T6opt,T7opt,T8opt,
            \[Alpha]popt,\[Beta]popt,\[Gamma]popt,
            \[Zeta]opt,
            \[Epsilon]opt],
    Transpose[{{B02opt,B04opt,B06opt,B22opt,B24opt,B26opt,B44opt,B46opt,B66opt,
                E1popt,E2popt,E3popt,
                M0opt,P2opt,
                T2popt,T3opt,T4opt,T6opt,T7opt,T8opt,
                \[Alpha]popt,\[Beta]popt,\[Gamma]popt,
                \[Zeta]opt,
                \[Epsilon]opt},
                startValues}],
    Method        -> Automatic,
    MaxIterations -> maxIterations,
    StepMonitor :>
  (
    steps+=1;
    currentSqSum  = 1/2*theCost;
    currentRMS    = Sqrt[currentSqSum/degsOfFree];
    currentParams = {B02opt,B04opt,B06opt,B22opt,B24opt,B26opt,B44opt,B46opt,B66opt,
                     E1popt,E2popt,E3popt,
                     M0opt,P2opt,
                     T2popt,T3opt,T4opt,T6opt,T7opt,T8opt,
                     \[Alpha]popt,\[Beta]popt,\[Gamma]popt,
                     \[Zeta]opt,
                     \[Epsilon]opt};
    AbortProtect[(
      paramSols     = AddToList[paramSols,currentParams,maxHistory,True];
      rmsHistory    = AddToList[rmsHistory,currentRMS,maxHistory,True];
    )];
    If[
      currentRMS < bestSigma,
      AbortProtect[
        (
        bestSol=currentParams;
        bestSigma = currentRMS
        )
      ]
    ]
  )
  ],
  maxTimeInSec,
  {bestSigma,AssociationThread[problemVars,bestSol]}
  ];

  Return[{
    "bestSigma" -> bestSigma,
    "bestSolution" -> AssociationThread[problemVars,bestSol],
    "Simulated Annealing Atlas" -> atlas
    }];
);