Needs["qlanth`"];
Needs["misc`"];
Needs["qplotter`"];
Needs["qonstants`"];
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
  eigensys, basis, assignmentMatches, stateLabels, carnallAssignments},
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
    params = LoadParameters[ln, "Free Ion" -> False];
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
      T11p -> 0, T11 -> 0, T12 -> 0, T14 -> 0, T15 -> 0, T16 -> 0, T18 -> 0, T17 -> 0, T19 -> 0};
    (* If the necessary symbolicHamiltonian is defined load it, if not, then calculate it *)
    simpleHam = If[
      ValueQ[symbolicHamiltonians[numEH]],
      symbolicHamiltonians[numEH],
      SimplerSymbolicHamMatrix[numE,
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
    params = ParamPad[params, "Print" -> False];
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
          ListLabelPlot[
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
        exportFname = FileNameJoin[{workDir,OptionValue["OutputDirectory"], ln <> " in " <> "LaF3" <> appToFname <> ".mx"}];
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
  "Print Function" -> PrintTemporary
  };
MagneticDipoleTransitionsLaF3Carnall[numE_Integer, OptionsPattern[]]:= (
  host           = OptionValue["Host"];
  \[Lambda]Range = OptionValue["Wavelength Range"];
  PrintFun       = OptionValue["Print Function"];
  {\[Lambda]min, \[Lambda]max} = OptionValue["Wavelength Range"];
  
  header    = {"\[Psi]i:simple","\[Psi]f:simple","\[Psi]i:idx","\[Psi]f:idx","Ei/K","Ef/K","\[Lambda]/nm","\[CapitalDelta]\[Lambda]/nm","\[Tau]/s","AMD/s^-1"};
  ln        = {"Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb"}[[numE]];
  {rmsDifference,carnallEnergies,eigenEnergies,ln,
  carnallAssignments,simplerStateLabels,eigensys,basis,truncatedStates} = Import["./examples/"<>ln<>" in LaF3 - example.mx"];
  
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
  magIon["fMD"] = GroundStateOscillatorStrength[eigensys, numE];
  
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
  basename    = ln <> " in " <> host <> " - example - " <> "M1.csv";
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
      nbFname = FileNameJoin[{workDir,"examples","M1 - "<>ln<>" in "<>"LaF3"<>".nb"}];
      PrintFun[">> Saving notebook to ",nbFname," ..."];
      NotebookSave[nb, nbFname];
      If[OptionValue["Close Notebook"],
          NotebookClose[nb];
      ];
  )
  ];

  Return[magIon];
  )

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
  "Remove Kramers"     -> True,
  "Append to Filename" -> "",
  "OutputDirectory"    -> "calcs",
  "Max Eigenstates in Table" -> 100,
  "eigenstateTruncationProbability" -> 0.9,
  "SymmetrySimplifier" -> {"C2v",
    {
      B12 -> 0,  B14 -> 0, B16 -> 0, B34 -> 0, B36 -> 0, B56 -> 0,
      S12 -> 0,  S14 -> 0, S16 -> 0, S22 -> 0, S24 -> 0, S26 -> 0, S34 -> 0, S36 -> 0, S44 -> 0,  S46 -> 0, S56 -> 0, S66 -> 0,
      T11p -> 0, T11 -> 0, T12 -> 0, T14 -> 0, T15 -> 0, T16 -> 0, T18 -> 0, T17 -> 0, T19 -> 0
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
    energyDiagram, appToFname, statesTable, DefaultIfMissing, 
    diffTableData, diffTable, nb, exportFname, tinyexportFname, tinyExport,
    nbFname
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
    PrintFun     = OptionValue["PrintFun"];
    makeNotebook = OptionValue["MakeNotebook"];
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
    
    If[numE>7,
      (
        PrintFun["> Conjugating the parameters accounting for the hole-particle equivalence ..."];
        params = HoleElectronConjugation[params];
        params[t2Switch] = 0;
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
      SimplerSymbolicHamMatrix[numE,
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
    params = ParamPad[params, "Print" -> False];
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
      PrintFun["> Putting together results in a notebook ..."];
      energyDiagram = Framed[
            EnergyLevelDiagram[eigensys, "Title" -> title, 
            "Explorer" -> OptionValue["Explorer"],
            "Background" -> White],
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
      diffTableData = Transpose[{simplerStateLabels, eigenEnergies}];
      diffTable     = TableForm[diffTableData, 
        TableHeadings -> {None, {"coarse label", 
        "E/\!\(\*SuperscriptBox[\(cm\), \(-1\)]\)"}}];

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
        TextCell[diffTable,
            "Output",
            TextAlignment -> Center
        ],
        TextCell["Truncated Eigenstates",
            "Section",
            TextAlignment -> Center
        ],
        TextCell[
            "These are some of the resultant eigenstates which add up to at least a total probability of " <> ToString[eigenstateTruncationProbability] <> ".", 
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
        NotebookWrite[nb,
          Cell[(
            "{eigenEnergies, ln, simplerStateLabels, eigensys, basis, truncatedStates} = Import[FileNameJoin[{NotebookDirectory[],\"" <> StringSplit[exportFname,"/"][[-1]] <> "\"}]];"
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
          {eigenEnergies, ln, simplerStateLabels, eigensys, basis, truncatedStates}
        ];
        tinyexportFname = FileNameJoin[
          {workDir, OptionValue["OutputDirectory"], ln <> " in " <> host  <> appToFname <> " - tiny.m"}
        ];
        tinyExport = <|"ln"->ln,
                    "eigenEnergies"-> eigenEnergies,
                    "simplerStateLabels" -> simplerStateLabels|>;
        Export[tinyexportFname, tinyExport];
      )
      ];
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
    - Line Strength, (* Line strength array *)
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

The function returns an association containing the following keys: Line Strength, AMD, fMD, Radiative lifetimes, Transition Energies / K, Transition Wavelengths in nm.";
Options[MagneticDipoleTransitions] = {
        "Make Notebook" -> True, 
        "Close Notebook" -> True,
        "Print Function" -> PrintTemporary,
        "Host" -> "LaF3",
        "Energy Uncertainty in K" -> 1.,
        "Wavelength Range" -> {50, 2000}};
MagneticDipoleTransitions[params_Association, fname_String, OptionsPattern[]]:= Module[
    {numE, host, PrintFun, \[Lambda]Range,
    header, ln, energy\[CapitalDelta]},
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
    host           = OptionValue["Host"];
    \[Lambda]Range = OptionValue["Wavelength Range"];
    PrintFun       = OptionValue["Print Function"];
    {\[Lambda]min, \[Lambda]max} = OptionValue["Wavelength Range"];
    
    header    = {"\[Psi]i:simple","\[Psi]f:simple","\[Psi]i:idx","\[Psi]f:idx","Ei/K","Ef/K","\[Lambda]/nm","\[CapitalDelta]\[Lambda]/nm","\[Tau]/s","AMD/s^-1"};
    ln        = {"Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb"}[[numE]];
    {eigenEnergies, ln, simplerStateLabels, eigensys, basis, truncatedStates} = Import[fname];
    
    (* Some of the above are not needed here *)
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
    magIon["AMD"] = MagDipoleRates[eigensys, numE, "Units"->"SI","Lifetime"->False];
    magIon["AMD"] = magIon["AMD"]/.{0.->Indeterminate};

    PrintFun["Calculating the oscillator strength for transitions from the ground state ..."];
    magIon["fMD"] = GroundStateOscillatorStrength[eigensys, numE];

    PrintFun["Calculating the natural radiative lifetimes ..."];
    magIon["Radiative lifetimes"] = 1./magIon["AMD"];

    PrintFun["Calculating the transition energies in K ..."];
    transitionEnergies = Outer[Subtract, First/@eigensys, First/@eigensys];
    magIon["Transition Energies / K"] = ReplaceDiagonal[transitionEnergies, Indeterminate];

    PrintFun["Calculating the transition wavelengths in nm ..."];
    magIon["Transition Wavelengths in nm"] = 10^7/magIon["Transition Energies / K"];
    
    energy\[CapitalDelta] = OptionValue["Energy Uncertainty in K"];
    PrintFun["Estimating the uncertainties in \[Lambda]/nm assuming a " <> ToString[energy\[CapitalDelta]] <> " K uncertainty in energies."];
    
    \[Lambda]uncertainty = energy\[CapitalDelta] * Sqrt[2.] * 10^-7 * magIon["Transition Wavelengths in nm"]^2;

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
        magIon["Transition Wavelengths in nm"],
        \[Lambda]uncertainty,
        magIon["AMD"],
        magIon["Radiative lifetimes"]
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
        SelectionMove[nb, After, Notebook];
        NotebookWrite[nb, Cell["Reload Data", "Section", TextAlignment -> Center]];
        NotebookWrite[nb, Cell[(
              "magTransitions = Import[FileNameJoin[{NotebookDirectory[],\"" <> StringSplit[rawexportFname,"/"][[-1]] <> "\"}],\""<>rawexportKey<>"\"];"
              ),"Input"]];
        SelectionMove[nb, Before, Notebook];
        nbFname = FileNameJoin[{DirectoryName[fname], "M1 - "<>ln<>" in "<>"LaF3"<>".nb"}];
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