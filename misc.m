BeginPackage["misc`"];

ExportToH5;
FlattenBasis;
RecoverBasis;
FlowMatching;
GreedyMatching;
HelperNotebook;
StochasticMatching;
ExtractSymbolNames;
GetModificationDate;
ToPythonSymPyExpression;

Begin["`Private`"];

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

FlowMatching::usage="FlowMatching[aList, bList] returns a list of pairs of elements from aList and bList that are closest to each other, this is returned in a list together with a mapping of indices from the aList to those in bList to which they were matched. The option \"alistLabels\" can be used to specify labels for the elements in aList. The option \"blistLabels\" can be used to specify labels for the elements in bList. If these options are used, the function returns a list with three elements the pairs of matched elements, the pairs of corresponding matched labels, and the mapping of indices. This is basically a wrapper around Mathematica's FindMinimumCostFlow function. By default the option \"noMatched\" is zero, and this means that all elements of aList must be matched to elements of bList. If this is not the case, the option \"noMatched\" can be used to specify how many elements of aList can be left unmatched. By default the cost function is Abs[#1-#2]&, but this can be changed with the option \"CostFun\", this function needs to take two arguments.";
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
   cFun = OptionValue["CostFun"];
   nA      = Length[aValues];
   nB      = Length[bValues];
   (*Build up the edges costs and capacities*)
   (*From source to the nodes representing the values of the first \
list*)
   edgesSourceToA = ("source" \[DirectedEdge] {"A", #}) & /@ Range[1, nA];
   capacitySourceToA = ConstantArray[1, nA];
   costSourceToA = ConstantArray[0, nA];
   
   (*From all the elements of A to all the elements of B*)
   midLayer = Table[{{"A", i} \[DirectedEdge] ({"B", j}), 1, cFun[aValues[[i]], bValues[[j]]]}, {i, 1, nA}, {j, 1, nB}];
   midLayer = Flatten[midLayer, 1];
   {midLayerEdges, midCapacities, midCosts} = Transpose[midLayer];
   
   (*From the elements of B to the sink*)
   edgesBtoSink = ({"B", #} \[DirectedEdge] "sink") & /@ Range[1, nB];
   capacityBtoSink = ConstantArray[1, nB];
   costBtoSink = ConstantArray[0, nB];
   
   (*Put it all together*)
   allCapacities = Join[capacitySourceToA, midCapacities, capacityBtoSink];
   allCosts = Join[costSourceToA, midCosts, costBtoSink];
   allEdges = Join[edgesSourceToA, midLayerEdges, edgesBtoSink];
   graph = Graph[allEdges, EdgeCapacity -> allCapacities, 
     EdgeCost -> allCosts];
   
   (*Solve it*)
   flow = FindMinimumCostFlow[graph, "source", "sink", nA - OptionValue["notMatched"], "OptimumFlowData"];
   (*Collect the pairs of matched indices*)
   pairedIndices = Select[flow["EdgeList"], And[Not[#[[1]] === "source"], Not[#[[2]] === "sink"]] &];
   pairedIndices = {#[[1, 2]], #[[2, 2]]} & /@ pairedIndices;
   (*Collect the pairs of matched values*)
   bestValues = {aValues[[#[[1]]]], bValues[[#[[2]]]]} & /@ pairedIndices;
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