Needs["qlanth`"];
Needs["misc`"];
Needs["qplotter`"];
Needs["qonstants`"]
LoadCarnall[];

workDir = DirectoryName[$InputFileName];

FastIonSolverLaF3::usage = "This function solves the energy levels of the given trivalent lanthanide in LaF3. The values for the Hamiltonian are simply taken from the values quoted by Carnall. It uses precomputed symbolic matrices for the Hamiltonian so it's faster than the previous alternatives.

The function returns a list with nine elements 
{rmsDifference, carnallEnergies, eigenEnergies, ln, carnallAssignments, simplerStateLabels, eigensys, basis, truncatedStates}.

Where:
1. rmsDifference is the root mean squared difference between the calculated values and those quoted by Carnall
2. carnallEnergies are the quoted calculated energies from Carnall;
3. eigenEnergies are the calculated energies (in the case of an odd number of electrons the Kramers degeneracy may have been removed from this list according to the option \"Remove Kramers\");
4. ln is simply a string labelling the corresponding lanthanide;
5. carnallAssignments is a list of strings providing the multiplet assignments that Carnall assumed;
6. simplerStateLabels is a list of strings providing the multiplet assignments that this function assumes;
7. eigensys is a list of tuples where the first element is the energy corresponding to the eigenvector given as the second element (in the case of an odd number of electrons the Kramers degeneracy may have been removed from this list according to the option \"Remove Kramers\");
8. basis is a list that specifies the basis in which the Hamiltonian was constructed and diagonalized, equal to BasisLSJMJ[numE];
9. Same as eigensys but the eigenvectors have been truncated so that the truncated version adds up to at least a total probability of eigenstateTruncationProbability.
";
Options[FastIonSolverLaF3] = {
  "MakeNotebook" -> True,
  "NotebookSave" -> True,
  "HTMLSave" -> False,
  "eigenstateTruncationProbability" -> 0.9,
  "Include spin-spin" -> True,
  "Max Eigenstates in Table" -> 100,
  "Sparse" -> True,
  "PrintFun" -> Print,
  "SaveData" -> True,
  "paramFiddle" -> {},
  "Append to Filename" -> "",
  "Remove Kramers" -> True,
  "OutputDirectory" -> "calcs",
  "Explorer" -> False
}; 
FastIonSolverLaF3[numE_, OptionsPattern[]] := Module[
  {makeNotebook, eigenstateTruncationProbability, host,
  ln, terms, termNames, carnallEnergies, eigenEnergies, simplerStateLabels,
  eigensys, basis, assignmentMatches, stateLabels, carnallAssignments},
  (
    PrintFun = OptionValue["PrintFun"];
    makeNotebook = OptionValue["MakeNotebook"];
    eigenstateTruncationProbability = OptionValue["eigenstateTruncationProbability"];
    maxStatesInTable = OptionValue["Max Eigenstates in Table"];
    Duplicator[aList_] := Flatten[{#, #} & /@ aList];
    host = "LaF3";
    paramFiddle = OptionValue["paramFiddle"];
    ln = theLanthanides[[numE]];
    terms = AllowedNKSLJTerms[Min[numE, 14 - numE]];
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
    params = LoadParameters[ln, "Free Ion" -> False];
    If[numE>7,
      (
        PrintFun["> Conjugating the parameters accounting for the hole-particle equivalence ..."];
        params = HoleElectronConjugation[params];
        params[t2Switch] = 0;
      ),
      params[t2Switch] = 1;
    ];

    Do[params[key] = paramFiddle[key], 
      {key, Keys[paramFiddle]}
    ];

    (* Import the symbolic Hamiltonian *)
    PrintFun["> Loading the symbolic Hamiltonian for this configuration ..."];
    startTime = Now;
    numH = 14 - numE;
    numEH = Min[numE, numH];
    C2vsimplifier = {B12 -> 0, B14 -> 0, B16 -> 0, B34 -> 0, B36 -> 0, 
      B56 -> 0,
      S12 -> 0, S14 -> 0, S16 -> 0, S22 -> 0, S24 -> 0, S26 -> 0, 
      S34 -> 0, S36 -> 0,
      S44 -> 0, S46 -> 0, S56 -> 0, S66 -> 0, T11p -> 0, T11 -> 0, 
      T12 -> 0, T14 -> 0, T15 -> 0,
      T16 -> 0, T18 -> 0, T17 -> 0, T19 -> 0};
    simpleHam = If[
      ValueQ[symbolicHamiltonians[numEH]],
      symbolicHamiltonians[numEH],
      SimplerSymbolicHamMatrix[numE, C2vsimplifier, "PrependToFilename" -> "C2v-", "Overwrite" -> False]
    ];
    endTime  = Now; 
    loadTime = QuantityMagnitude[endTime - startTime, "Seconds"];
    PrintFun[">> Loading the symbolic Hamiltonian took ", loadTime, " seconds."];

    (*Enforce the override to the spin-spin contribution to the magnetic interactions*)
    params[\[Sigma]SS] = If[OptionValue["Include spin-spin"], 1, 0];

    (*Everything that is not given is set to zero*)
    params = ParamPad[params, "Print" -> False];
    PrintFun[params];
    (* numHam = simpleHam /. params; *)
    numHam = ReplaceInSparseArray[simpleHam, params];
    If[Not[OptionValue["Sparse"]],
        numHam = Normal[numHam]
    ];
    PrintFun["> Calculating the SLJ basis ..."];
    basis = BasisLSJMJ[numE];

    (* Eigensolver *)
    PrintFun["> Diagonalizing the numerical Hamiltonian ..."];
    startTime = Now;
    eigensys  = Eigensystem[numHam];
    endTime   = Now;
    diagonalTime = QuantityMagnitude[endTime - startTime, "Seconds"];
    PrintFun[">> Diagonalization took ", diagonalTime, " seconds."];
    eigensys = Chop[eigensys];
    eigensys = Transpose[eigensys];

    (*Shift the baseline energy*)
    eigensys = ShiftedLevels[eigensys];
    (*Sort according to energy*)
    eigensys = SortBy[eigensys, First];
    (*Grab just the energies*)
    eigenEnergies = First /@ eigensys;

    (*Energies are doubly degenerate in the case of odd number of electrons, keep only one*)
    If[And[OddQ[numE], OptionValue["Remove Kramers"]], 
        (
        PrintFun["> Since there's an odd number of electrons energies come in pairs, taking just one for each pair ..."];
        eigenEnergies = eigenEnergies[[;; ;; 2]];
        )
    ];

    (*Compare against the data quoted by Bill Carnall*)
    PrintFun["> Comparing against the data from Carnall ..."];
    mainKey            = StringTemplate["appendix:`Ln`:Association"][<|"Ln" -> ln|>];
    lnData             = Carnall[mainKey];
    carnalKeys         = lnData // Keys;
    repetitions        = Length[lnData[#]["Calc (1/cm)"]] & /@ carnalKeys;
    carnallAssignments = First /@ Carnall["appendix:" <> ln <> ":RawTable"];
    carnalKey          = StringTemplate["appendix:`Ln`:Calculated"][<|"Ln" -> ln|>];
    carnallEnergies    = Carnall[carnalKey];
    If[And[OddQ[numE], Not[OptionValue["Remove Kramers"]]],
    (
      PrintFun[">> The number of eigenstates and the number of quoted states don't match, removing the last state ..."];
      carnallAssignments = Duplicator[carnallAssignments];
      carnallEnergies    = Duplicator[carnallEnergies];
    )
    ];

    (* For the difference take as many energies as quoted by Bill*)
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
      assignmentMatches = 
      If[StringContainsQ[#[[1]], #[[2]]], "\[Checkmark]", "X"] & /@ 
        Transpose[{carnallAssignments, simplerStateLabels[[;; Length[carnallAssignments]]]}];
      assignmentMatches = {{"\[Checkmark]", 
        Count[assignmentMatches, "\[Checkmark]"]}, {"X", 
        Count[assignmentMatches, "X"]}};
      labelComparison = (If[StringContainsQ[#[[1]], #[[2]]], "\[Checkmark]", "X"] & /@ 
        Transpose[{carnallAssignments, 
        simplerStateLabels[[;; Length[carnallAssignments]]]}]);
      labelComparison = 
      PadRight[labelComparison, Length[simplerStateLabels], "-"];

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
      diffs = Pick[diffs,notBad];
      diffHistogram = Histogram[diffs,
        Frame -> True,
        ImageSize -> 800, 
        AspectRatio -> 1/3, FrameStyle -> Directive[16], 
        FrameLabel -> {"(qlanth-carnall)/Ky", "Freq"}
      ];
      rmsDifference = Sqrt[Total[diffs^2/Length[diffs]]];
      labelTempate = StringTemplate["\!\(\*SuperscriptBox[\(`ln`\), \(\(3\)\(+\)\)]\)"];
      diffData   = diffs;
      diffLabels = simplerStateLabels[[;;Length[notBad]]];
      diffLabels = Pick[diffLabels, notBad];
      diffPlot   = Framed[
        ListLabelPlot[diffData,
        diffLabels,
        Frame -> True,
        PlotRange -> All,
        ImageSize -> 1200,
        AspectRatio -> 1/3,
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
        TextCell["Energy Diagram",
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
        TextCell[diffHistogram, TextAlignment -> Center],
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
        exportFname = FileNameJoin[{workDir,OptionValue["OutputDirectory"], ln <> " in " <> "LaF3" <> appToFname <> ".m"}];
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
      If[OptionValue["HTMLSave"],
      (
        htmlFname = FileNameJoin[{workDir,OptionValue["OutputDirectory"], "html", ln <> " in " <> "LaF3" <> appToFname <> ".html"}];
        PrintFun[">> Saving html version to ", htmlFname, " ..."];
        Export[htmlFname, nb];
      )
      ];
    )
    ];

    Return[{rmsDifference, carnallEnergies, eigenEnergies, ln, carnallAssignments, simplerStateLabels, eigensys, basis, truncatedStates}];
  )
];

MagneticDipoleTransitions::usage = "MagneticDipoleTransitions[numE] calculates the magnetic dipole transitions for the lanthanide ion numE in LaF3. The output is a tabular file, a raw data file, and a CSV file. The tabular file contains the following columns: 
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
    - \"Host\" -> \"LaF3\". The host material. Default is LaF3.
    - \"Wavelength Range\" -> {50,2000}. The range of wavelengths in nm for the Manipulate object in the created notebook. Default is {50,2000}.

The function returns an association containing the following keys: Line Strength, AMD, fMD, Radiative lifetimes, Transition Energies / K, Transition Wavelengths in nm.";
Options[MagneticDipoleTransitions] = {
        "Make Notebook" -> True, 
        "Close Notebook" -> True,
        "Print Function" -> PrintTemporary,
        "Host" -> "LaF3",
        "Wavelength Range" -> {50,2000}};
MagneticDipoleTransitions[numE_Integer, OptionsPattern[]]:= (
  host           = OptionValue["Host"];
  \[Lambda]Range = OptionValue["Wavelength Range"];
  PrintFun       = OptionValue["Print Function"];
  {\[Lambda]min, \[Lambda]max} = OptionValue["Wavelength Range"];

  header    = {"\[Psi]i:simple","\[Psi]f:simple","\[Psi]i:idx","\[Psi]f:idx","Ei/K","Ef/K","\[Lambda]/nm","\[CapitalDelta]\[Lambda]/nm","\[Tau]/s","AMD/s^-1"};
  ln        = {"Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb"}[[numE]];
  {rmsDifference,carnallEnergies,eigenEnergies,ln,
  carnallAssignments,simplerStateLabels,eigensys,basis,truncatedStates} = Import["./examples/"<>ln<>" in LaF3 - example.m"];

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
  PrintFun["Calculating the magnetic dipole line strength array..."];
  magIon["Line Strength"] = magIon;MagDipLineStrength[eigensys, numE, "Reload MagOp" -> False, "Units" -> "SI"];

  PrintFun["Calculating the M1 spontaneous transition rates ..."];
  magIon["AMD"] = MagDipoleRates[eigensys, numE, "Units"->"SI","Lifetime"->False];
  magIon["AMD"] = magIon["AMD"]/.{0.->Indeterminate};

  PrintFun["Calculating the oscillator strength for transition from the ground state ..."];
  magIon["fMD"]=GroundStateOscillatorStrength[eigensys, numE];

  PrintFun["Calculating the natural radiative lifetims ..."];
  magIon["Radiative lifetimes"]= 1/magIon["AMD"];

  PrintFun["Calculating the transition energies in K ..."];
  transitionEnergies=Outer[Subtract,First/@eigensys,First/@eigensys];
  magIon["Transition Energies / K"]=ReplaceDiagonal[transitionEnergies,Indeterminate];

  PrintFun["Calculating the transition wavelengths in nm ..."];
  magIon["Transition Wavelengths in nm"] =10^7/magIon["Transition Energies / K"];

  PrintFun["Estimating the uncertainties in \[Lambda]/nm assuming a 1 K uncertainty in energies."];
  (*Assuming an uncertainty of 1 K in both energies used to calculate the wavelength*)
  \[Lambda]uncertainty=Sqrt[2]*magIon["Transition Wavelengths in nm"]^2*10^-7;

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
  basename     = ln <> " in " <> host <> " - example - " <> "MD1 - tabular.zip";
  exportFname  = FileNameJoin[{"./examples",basename}];
  PrintFun["Exporting tabular data to "<>exportFname<>" ..."];
  exportKey    = StringReplace[basename,".zip"->".m"];
  Export[exportFname, <|exportKey->finalTable|>];

  (* raw data output *)
  basename    = ln <> " in " <> host <> " - example - " <> "MD1 - raw.zip";
  rawexportFname = FileNameJoin[{"./examples",basename}];
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
  basename    = ln <> " in " <> host <> " - example - " <> "MD1.csv";
  exportFname = FileNameJoin[{"./examples", basename}];
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
      WindowTitle    -> "MD1 - "<>ln<>" in "<>host,
      WindowSize     -> {1600,800}
      ];
      SelectionMove[nb, After, Notebook];
      NotebookWrite[nb, Cell["Reload Data", "Section", TextAlignment -> Center]];
      NotebookWrite[nb, Cell[(
            "magTransitions = Import[FileNameJoin[{NotebookDirectory[],\"" <> StringSplit[rawexportFname,"/"][[-1]] <> "\"}],\""<>rawexportKey<>"\"];"
            ),"Input"]];
      SelectionMove[nb, Before, Notebook];
      nbFname = FileNameJoin[{workDir,"examples","MD1 - "<>ln<>" in "<>"LaF3"<>".nb"}];
      PrintFun[">> Saving notebook to ",nbFname," ..."];
      NotebookSave[nb, nbFname];
      If[OptionValue["Close Notebook"],
          NotebookClose[nb];
      ];
  )
  ];

  Return[magIon];
)