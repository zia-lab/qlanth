BeginPackage["misc`"];

ExportToH5;
FlattenBasis;
RecoverBasis;
FlowMatching;
SuperIdentity;
RobustMissingQ;
ReplaceDiagonal;

GreedyMatching;
HelperNotebook;
StochasticMatching;
ExtractSymbolNames;
GetModificationDate;
TextBasedProgressBar;
ToPythonSparseFunction;

FirstOrderPerturbation;
SecondOrderPerturbation;
RoundValueWithUncertainty;

ToPythonSymPyExpression;
RoundToSignificantFigures;
RobustMissingQ;

BlockMatrixMultiply;
BlockAndIndex;
TruncateBlockArray;
BlockArrayDimensionsArray;
ArrayBlocker;
BlockTranspose;
RemoveTrailingDigits;

Begin["`Private`"];

  RemoveTrailingDigits[s_String] := StringReplace[s, RegularExpression["\\d+$"] -> ""];

  BlockTranspose[anArray_]:=(
    Map[Transpose, Transpose[anArray], {2}]
  );

  BlockMatrixMultiply::usage="BlockMatrixMultiply[A,B] gives the matrix multiplication of A and B, with A and B having a compatible block structure that allows for matrix multiplication into a congruent block structure.";
  BlockMatrixMultiply[Amat_,Bmat_]:=Module[{rowIdx,colIdx,sumIdx},
  (
    Table[
      Sum[Amat[[rowIdx,sumIdx]].Bmat[[sumIdx,colIdx]],{sumIdx,1,Dimensions[Amat][[2]]}],
    {rowIdx,1,Dimensions[Amat][[1]]},
    {colIdx,1,Dimensions[Bmat][[2]]}
    ]
  )
  ];

  BlockAndIndex::usage="BlockAndIndex[blockSizes, index] takes a list of bin widths and index. The function return in which block the index would be, were the bins to be layed out from left to right. The function also returns the position within the bin in which it is accomodated. The function returns these two numbers as a list of two elements {blockIndex, blockSubIndex}";
  BlockAndIndex[blockSizes_List, index_Integer]:=Module[{accumulatedBlockSize,blockIndex, blockSubIndex},
  (
    accumulatedBlockSize = Accumulate[blockSizes];
    If[accumulatedBlockSize[[-1]]-index<0,
      Print["Index out of bounds"];
      Abort[]
    ];
    blockIndex    = Flatten[Position[accumulatedBlockSize-index,n_ /; n>=0]][[1]];
    blockSubIndex = Mod[index-accumulatedBlockSize[[blockIndex]],blockSizes[[blockIndex]],1];
    Return[{blockIndex,blockSubIndex}]
  )
  ];

  TruncateBlockArray::usage="TruncateBlockArray[blockArray, truncationIndices, blockWidths] takes a an array of blocks and selects the columns and rows corresponding to truncationIndices. The indices being given in what would be the ArrayFlatten[blockArray] version of the array. They blocks in the given array may be SparseArray. This is equivalent to FlattenArray[blockArray][truncationIndices, truncationIndices] but may be more efficient blockArray is sparse.";
  TruncateBlockArray[blockArray_,truncationIndices_,blockWidths_]:=Module[
  {truncatedArray,blockCol,blockRow,blockSubCol,blockSubRow},(
  truncatedArray = Table[
    {blockCol,blockSubCol} = BlockAndIndex[blockWidths,fullColIndex];
    {blockRow,blockSubRow} = BlockAndIndex[blockWidths,fullRowIndex];
    blockArray[[blockRow,blockCol]][[blockSubRow,blockSubCol]],
  {fullColIndex,truncationIndices},
  {fullRowIndex,truncationIndices}
  ];
  Return[truncatedArray]
  )
  ];

  BlockArrayDimensionsArray::usage="BlockArrayDimensionsArray[blockArray] returns the array of block sizes in a given blocked array.";
  BlockArrayDimensionsArray[blockArray_]:=(
    Map[Dimensions,blockArray,{2}]
  );

  ArrayBlocker::usage="ArrayBlocker[anArray, blockSizes] takes a flat 2d array and a congruent 2D array of block sizes, and with them it returns the original array with the block structure imposed by blockSizes. The resulting array satisfies ArrayFlatten[blockedArray] == anArray, and also Map[Dimensions, blockedArray,{2}] == blockSizes.";
  ArrayBlocker[anArray_,blockSizes_]:=Module[{rowStart,colStart,colEnd,numBlocks,blockedArray,blockSize,rowEnd,aBlock,idxRow,idxCol},(
    rowStart  = 1;
    colStart  = 1;
    colEnd    = 1;
    numBlocks = Length[blockSizes];
    blockedArray = Table[(
      blockSize = blockSizes[[idxRow,idxCol]];
      rowEnd    = rowStart+blockSize[[1]]-1;
      colEnd    = colStart+blockSize[[2]]-1;
      aBlock    = anArray[[rowStart;;rowEnd,colStart;;colEnd]];
      colStart  = colEnd+1;
      If[idxCol==numBlocks,
        rowStart=rowEnd+1;
        colStart=1;
      ];
      aBlock
    ),
    {idxRow,1,numBlocks},
    {idxCol,1,numBlocks}
    ];
    Return[blockedArray]
  )
  ];

  ReplaceDiagonal::usage = 
    "ReplaceDiagonal[matrix, repValue] replaces all the diagonal of the given array to the given value. The array is assumed to be square and the replacement value is assumed to be a number. The returned value is the array with the diagonal replaced. This function is useful for setting the diagonal of an array to a given value. The original array is not modified. The given array may be sparse.";
  ReplaceDiagonal[matrix_, repValue_] := 
    ReplacePart[matrix, 
      Table[{i, i} -> repValue, {i, 1, Length[matrix]}]];

  Options[RoundValueWithUncertainty] = {"SetPrecision" -> False};
  RoundValueWithUncertainty::usage = "RoundValueWithUncertainty[x,dx] given a number x together with an uncertainty dx this function rounds x to the first significant figure of dx and also rounds dx to have a single significant figure. 
  The returned value is a list with the form {roundedX, roundedDx}.
  The option \"SetPrecision\" can be used to control whether the Mathematica precision of x and dx is also set accordingly to these rules, otherwise the rounded numbers still have the original precision of the input values. 
  If the position of the first significant figure of x is after the position of the first significant figure of dx, the function returns {0,dx} with dx rounded to one significant figure.";
  RoundValueWithUncertainty[x_, dx_, OptionsPattern[]] := Module[
    {xExpo, dxExpo, sigFigs, roundedX, roundedDx, returning},
    (
      xExpo   = RealDigits[x][[2]];
      dxExpo  = RealDigits[dx][[2]];
      sigFigs = (xExpo - dxExpo) + 1;
      {roundedX, roundedDx} = If[sigFigs <= 0,
        {0., N@RoundToSignificantFigures[dx, 1]},
        N[
        {
          RoundToSignificantFigures[x, xExpo - dxExpo + 1],
          RoundToSignificantFigures[dx, 1]}
        ]
        ];
      returning = If[
        OptionValue["SetPrecision"],
        {SetPrecision[roundedX, Max[1, sigFigs]],
        SetPrecision[roundedDx, 1]},
        {roundedX, roundedDx}
        ];
      Return[returning]
      )
    ];

  RoundToSignificantFigures::usage = 
    "RoundToSignificantFigures[x, sigFigs] rounds x so that it only has \
  sigFigs significant figures.";
  RoundToSignificantFigures[x_, sigFigs_] := 
    Sign[x]*N[FromDigits[RealDigits[x, 10, sigFigs]]];

  RobustMissingQ[expr_] := (FreeQ[expr, _Missing] === False);

  TextBasedProgressBar[progress_, totalIterations_, prefix_:""] := Module[
      {progMessage},
      progMessage = ToString[progress] <> "/" <> ToString[totalIterations];
      If[progress < totalIterations,
          WriteString["stdout", StringJoin[prefix, progMessage, "\r"]],
          WriteString["stdout", StringJoin[prefix, progMessage, "\n"]]
      ];
  ];

  FirstOrderPerturbation::usage="Given the eigenVectors of a matrix A (which doesn't need to be given) together with a corresponding perturbation matrix perMatrix, this function calculates the first derivative of the eigenvalues with respect to the scale factor of the perturbation matrix. In the sense that the eigenvalues of the matrix A + \[Beta] perMatrix are to first order equal to \[Lambda] + \[Delta]_i \[Beta], where the \[Delta]_i are the returned values. This assuming that the eigenvalues are non-degenerate.";
  FirstOrderPerturbation[eigenVectors_, 
    perMatrix_] := (Diagonal[
    eigenVectors . perMatrix . Transpose[eigenVectors]])

  SecondOrderPerturbation::usage="Given the eigenValues and eigenVectors of a matrix A (which doesn't need to be given) together with a corresponding perturbation matrix perMatrix, this function calculates the second derivative of the eigenvalues with respect to the scale factor of the perturbation matrix. In the sense that the eigenvalues of the matrix A + \[Beta] perMatrix are to second order equal to \[Lambda] + \[Delta]_i \[Beta] + \[Delta]_i^{(2)}/2 \[Beta]^2, where the \[Delta]_i^{(2)} are the returned values. The eigenvalues and eigenvectors are assumed to be given in the same order, i.e. the ith eigenvalue corresponds to the ith eigenvector. This assuming that the eigenvalues are non-degenerate.";
  SecondOrderPerturbation[eigenValues_, eigenVectors_, perMatrix_] := (
    dim = Length[perMatrix];
    eigenBras = Conjugate[eigenVectors];
    eigenKets = eigenVectors;
    matV = Abs[eigenBras . perMatrix . Transpose[eigenKets]]^2;
    OneOver[x_, y_] := If[x == y, 0, 1/(x - y)];
    eigenDiffs = Outer[OneOver, eigenValues, eigenValues, 1];
    pProduct = Transpose[eigenDiffs]*matV;
    Return[2*(Total /@ Transpose[pProduct])];
    )

  SuperIdentity::usage="SuperIdentity[args] returns the arguments passed to it. This is useful for defining a function that does nothing, but that can be used in a composition.";
  SuperIdentity[args___] := {args};

  FlattenBasis::usage="FlattenBasis[basis] takes a basis in the standard representation and separates out the strings that describe the LS part of the labels and the additional numbers that define the values of J MJ and MI. It returns a list with two elements {flatbasisLS, flatbasisNums}. This is useful for saving the basis to an h5 file where the strings and numbers need to be separated.";
  FlattenBasis[basis_] := Module[{flatbasis, flatbasisLS, flatbasisNums},
    (
      flatbasis = Flatten[basis];
      flatbasisLS = flatbasis[[1 ;; ;; 4]];
      flatbasisNums = Select[flatbasis, Not[StringQ[#]] &];
      Return[{flatbasisLS, flatbasisNums}]
      )
    ];

  RecoverBasis::usage="RecoverBasis[{flatBasisLS, flatbasisNums}] takes the output of FlattenBasis and returns the original basis. The input is a list with two elements {flatbasisLS, flatbasisNums}.";
  RecoverBasis[{flatbasisLS_, flatbasisNums_}] := Module[{recBasis},
    (
    recBasis = {{{#[[1]], #[[2]]}, #[[3]]}, #[[4]]} & /@ (Flatten /@ 
        Transpose[{flatbasisLS, 
          Partition[Round[2*#]/2 & /@ flatbasisNums, 3]}]);
    Return[recBasis];
    )
    ]

  ExtractSymbolNames[expr_Hold] := Module[
    {strSymbols},
    strSymbols = ToString[expr, InputForm];
    StringCases[strSymbols, RegularExpression["\\w+"]][[2 ;;]]
    ]

  ExportToH5::usage = 
    "ExportToH5[fname, Hold[{symbol1, symbol2, ...}]] takes an .h5 filename and a held list of symbols and export to the .h5 file the values of the symbols with keys equal the symbol names. The values of the symbols cannot be arbitrary, for instance a list with mixes numbers and string will fail, but an Association with mixed values exports ok. Do give it a try. 
    If the file is already present in disk, this function will overwrite it by default. If the value of a given symbol contains symbolic numbers, e.g. \[Pi], these will be converted to floats in the exported file.";
  Options[ExportToH5] = {"Overwrite" -> True};
  ExportToH5[fname_String, symbols_Hold, OptionsPattern[]] := (
    If[And[FileExistsQ[fname], OptionValue["Overwrite"]],
    (
      Print["File already exists, overwriting ..."];
      DeleteFile[fname];
      )
    ];
    symbolNames = ExtractSymbolNames[symbols];
    Do[(Print[symbolName];
      Export[fname, ToExpression[symbolName], {"Datasets", symbolName}, 
      OverwriteTarget -> "Append"]
      ), {symbolName, symbolNames}]
    )

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
    ]

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
    ]

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
    ]

  HelperNotebook::usage="HelperNotebook[nbName] creates a separate notebook and returns a function that can be used to print to the bottom of it. The name of the notebook, nbName, is optional and defaults to OUT.";
  HelperNotebook[nbName_:"OUT"] := 
  Module[{screenDims, screenWidth, screenHeight, nbWidth, leftMargin, 
    PrintToOutputNb}, (
    screenDims = 
      SystemInformation["Devices", "ScreenInformation"][[1, 2, 2]];
    screenWidth = screenDims[[1, 2]];
    screenHeight = screenDims[[2, 2]];
    nbWidth = Round[screenWidth/3];
    leftMargin = screenWidth - nbWidth;
    outputNb = CreateDocument[{}, WindowTitle -> nbName,
      WindowMargins -> {{leftMargin, Automatic}, {Automatic, 
          Automatic}}, WindowSize -> {nbWidth, screenHeight}];
    PrintToOutputNb[text_] := 
      (
          SelectionMove[outputNb, After, Notebook];
          NotebookWrite[outputNb, Cell[BoxData[ToBoxes[text]], "Output"]];
      );
    Return[PrintToOutputNb]
    )
    ]

  GetModificationDate::usage="GetModificationDate[fname] returns the modification date of the given file.";
  GetModificationDate[theFileName_] := FileDate[theFileName, "Modification"];

  (*Helper function to convert Mathematica expressions to standard form*)
  StandardFormExpression[expr0_] := Module[{expr=expr0}, ToString[expr, InputForm]];

  (*Helper function to translate to Python/SymPy expressions*)
  ToPythonSymPyExpression::usage="ToPythonSymPyExpression[expr] converts a Mathematica expression to a SymPy expression. This is a little iffy and might break if the expression includes Mathematica functions that haven't been given a SymPy equivalent.";
  ToPythonSymPyExpression[expr0_] := Module[{standardForm, expr=expr0},
    standardForm = StandardFormExpression[expr];
    StringReplace[standardForm, {
      "Power[" -> "Pow(",
      "Sqrt[" -> "sqrt(",
      "[" -> "(",
      "]" -> ")",
      "\\" -> "",
      "I" -> "1j",
      (*Remove special Mathematica backslashes*)
      "/" -> "/" (*Ensure division is represented with a slash*)}]];

  ToPythonSparseFunction[sparseArray_SparseArray, funName_] := 
    Module[{data, rowPointers, columnIndices, dimensions, pyCode, vars, 
      varList, dataPyList, 
      colIndicesPyList},(*Extract unique symbolic variables from the \
  SparseArray*)
    vars = Union[Cases[Normal[sparseArray], _Symbol, Infinity]];
    varList = StringRiffle[ToString /@ vars, ", "];
    (*varList=ToPythonSymPyExpression/@varList;*)
    (*Convert data to SymPy compatible strings*)
    dataPyList = 
      StringRiffle[
      ToPythonSymPyExpression /@ Normal[sparseArray["NonzeroValues"]], 
      ", "];
    colIndicesPyList = 
      StringRiffle[
      ToPythonSymPyExpression /@ (Flatten[
          Normal[sparseArray["ColumnIndices"]] - 1]), ", "];
    (*Extract sparse array properties*)
    rowPointers = Normal[sparseArray["RowPointers"]];
    dimensions = Dimensions[sparseArray];
    (*Create Python code string*)pyCode = StringJoin[
      "#!/usr/bin/env python3\n\n",
      "from scipy.sparse import csr_matrix\n",
      "from sympy import *\n",
      "import numpy as np\n",
      "\n",
      "sqrt = np.sqrt\n",
      "\n",
      "def ", funName, "(",
      varList,
      "):\n",
      "    data = np.array([", dataPyList, "])\n",
      "    indices = np.array([",
      colIndicesPyList,
      "])\n",
      "    indptr = np.array([", 
      StringRiffle[ToString /@ rowPointers, ", "], "])\n",
      "    shape = (", StringRiffle[ToString /@ dimensions, ", "],
      ")\n",
      "    return csr_matrix((data, indices, indptr), shape=shape)"];
    pyCode
    ];

End[];
EndPackage[];