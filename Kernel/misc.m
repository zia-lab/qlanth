BeginPackage["misc`"];
(* Needs["MaTeX`"]; *)

ArrayBlocker;
BlockAndIndex;
BlockArrayDimensionsArray;
BlockMatrixMultiply;
BlockTranspose;

EllipsoidBoundingBox;
ExportToH5;
FirstOrderPerturbation;
FlattenBasis;

GetModificationDate;
HamTeX;
HelperNotebook;

RecoverBasis;
RemoveTrailingDigits;
ReplaceDiagonal;
RobustMissingQ;

RoundToSignificantFigures;
RoundValueWithUncertainty;
SecondOrderPerturbation;
SuperIdentity;

TextBasedProgressBar;
ToPythonSparseFunction;
ToPythonSymPyExpression;
TruncateBlockArray;

Begin["`Private`"];

  ExtractSymbolNames;

  RemoveTrailingDigits[s_String] := StringReplace[s, RegularExpression["\\d+$"] -> ""];

  BlockTranspose::usage="BlcockTranspose[anArray] takes a 2D array with a congruent block structure and returns the transposed array with the same block structure.";
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

  BlockAndIndex::usage="BlockAndIndex[{w1, w2, ...}, idx] takes a list of block lengths wi and an index idx. The function returns in which block idx would be in a a list defined by {Range[1,w1], Range[1+w1,w1+w2], ...}. The function also returns the position within the bin in which the given index would be found in. The function returns these two numbers as a list of two elements {blockIndex, blockSubIndex}.";
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

  TruncateBlockArray::usage="TruncateBlockArray[blockArray, truncationIndices, blockWidths] takes an array of blocks and selects the columns and rows corresponding to truncationIndices. The indices being given in what would be the ArrayFlatten[blockArray] version of the array. The blocks in the given array may be SparseArray. This is equivalent to FlattenArray[blockArray][truncationIndices, truncationIndices] but may be more efficient if blockArray is sparse.";
  TruncateBlockArray[blockArray_,truncationIndices_,blockWidths_]:=Module[
  {truncatedArray,blockCol,blockRow,blockSubCol,blockSubRow},(
  truncatedArray = Table[
    {blockCol,blockSubCol} = BlockAndIndex[blockWidths,fullColIndex];
    {blockRow,blockSubRow} = BlockAndIndex[blockWidths,fullRowIndex];
    blockArray[[blockRow,blockCol]][[blockSubRow,blockSubCol]],
  {fullRowIndex,truncationIndices},
  {fullColIndex,truncationIndices}
  ];
  Return[truncatedArray]
  )
  ];

  BlockArrayDimensionsArray::usage="BlockArrayDimensionsArray[blockArray] returns the array of block sizes in a given blocked array.";
  BlockArrayDimensionsArray[blockArray_]:=(
    Map[Dimensions,blockArray,{2}]
  );

  ArrayBlocker::usage="ArrayBlocker[anArray, blockSizes] takes a flat 2d array and a congruent 2D array of block sizes, and with them it returns the original array with the block structure imposed by blockSizes. The resulting array satisfies ArrayFlatten[blockedArray] == anArray, and also Map[Dimensions, blockedArray,{2}] == blockSizes.";
  ArrayBlocker[anArray_, blockSizes_] := Module[{rowStart,colStart,colEnd,numBlocks,blockedArray,blockSize,rowEnd,aBlock,idxRow,idxCol},(
    rowStart  = 1;
    colStart  = 1;
    colEnd    = 1;
    case      = Length[Dimensions[blockSizes]];
    Which[
      case == 3,
      (
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
      {idxCol,1,numBlocks}];
      ),
      case == 1,
      (
        expandedSizes = Table[
          {blockSizes[[idxRow]], blockSizes[[idxCol]]},
          {idxRow,1,Length[blockSizes]},
          {idxCol,1,Length[blockSizes]}
          ]; 
        numBlocks = Length[expandedSizes];
        blockedArray = Table[(
          blockSize = expandedSizes[[idxRow,idxCol]];
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
      {idxCol,1,numBlocks}];
      )
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

  FirstOrderPerturbation::usage="Given the eigenVectors of a matrix A (which doesn't need to be given) together with a corresponding perturbation matrix perMatrix, this function calculates the first derivative of the eigenvalues with respect to the scale factor of the perturbation matrix. In the sense that the eigenvalues of the matrix A + \[Beta] perMatrix are to first order equal to \[Lambda]_i + \[Delta]_i \[Beta], where the \[Delta]_i are the returned values. This assuming that the eigenvalues are non-degenerate.";
  FirstOrderPerturbation[eigenVectors_, 
    perMatrix_] := (Chop@Diagonal[
    Conjugate@eigenVectors . perMatrix . Transpose[eigenVectors]])

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

  ExtractSymbolNames[expr_Hold] := Module[
    {strSymbols},
    strSymbols = ToString[expr, InputForm];
    StringCases[strSymbols, RegularExpression["\\w+"]][[2 ;;]]
    ]

  ExportToH5::usage = 
    "ExportToH5[fname, Hold[{symbol1, symbol2, ...}]] takes an .h5 filename and a held list of symbols and export to the .h5 file the values of the symbols with keys equal the symbol names. The values of the symbols cannot be arbitrary, for instance a list which mixes numbers and strings will fail, but an Association with mixed values exports ok. Do give it a try. 
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
      "Sqrt[" -> "sq(",
      "[" -> "(",
      "]" -> ")",
      "\\" -> "",
      "I" -> "1j",
      "^" -> "**",
      (*Remove special Mathematica backslashes*)
      "/" -> "/" (*Ensure division is represented with a slash*)}]];

  ToPythonSparseFunction::usage="ToPythonSparseFunction[array, funName] takes a 2D array (whose elements are linear combinations of a set of variables) and returns a string which defines a Python function that can be used to evaluate the given array in Python. In Python the name of the function is equal to funName. The given array may be a list of lists or a 2D SparseArray.";
  ToPythonSparseFunction[arrayInput_, funName_] := 
    Module[{
      sparseArray = arrayInput,
      rowPointers,
      dimensions,
      pyCode,
      vars, 
      varList,
      dataPyList, 
      colIndicesPyList
      },
    (
      If[Head[sparseArray] === List,
        sparseArray = SparseArray[sparseArray]
        ];
      (* Extract unique symbolic variables from the SparseArray *)
      vars       = Union[Cases[Normal[sparseArray], _Symbol, Infinity]];
      varList    = StringRiffle[ToString /@ vars, ", "];
      (* Convert data to SymPy compatible strings *)
      dataPyList = StringRiffle[ToPythonSymPyExpression /@ Normal[sparseArray["NonzeroValues"]], ",\n        "];
      colIndicesPyList = StringRiffle[
        ToPythonSymPyExpression /@ (Flatten[Normal[sparseArray["ColumnIndices"]] - 1]), ", "];
      (* Extract sparse array properties *)
      rowPointers = Normal[sparseArray["RowPointers"]];
      dimensions = Dimensions[sparseArray];
      (*Create Python code string*)
      pyCode = StringJoin[
        "#!/usr/bin/env python3\n\n",
        "from scipy.sparse import csr_matrix\n",
        "import numpy as np\n",
        "\n",
        "sq = np.sqrt\n",
        "\n",
        "def ", funName, "(",
        varList,
        "):\n",
        "    data = np.array([\n        ", dataPyList, "\n        ])\n",
        "    indices = np.array([",
        colIndicesPyList,
        "])\n",
        "    indptr = np.array([", 
        StringRiffle[ToString /@ rowPointers, ", "], "])\n",
        "    shape = (", StringRiffle[ToString /@ dimensions, ", "],
        ")\n",
        "    return csr_matrix((data, indices, indptr), shape=shape)\n"];
      Return[pyCode];
    )
    ];

  Options[HamTeX] = {"T2" -> False};
  HamTeX::usage="HamTeX[numE] returns the LaTeX code for the semi-empirical Hamiltonian for f^numE. The option \"T2\" can be used to specify whether the T2 term should be included in the Hamiltonian for the f^12 configuration. The default is False and the option is ignored if the number of electrons is not 12.";
  HamTeX[nE_, OptionsPattern[]] := (
    tex = Which[
      MemberQ[{1, 13}, nE],
        "\\hat{H}&=\\zeta \\sum_{i=1}^{n}\\left(\\hat{s}_i \\cdot \
    \\hat{l}_i\\right) \
    +\\sum_{i=1}^n\\sum_{k=2,4,6}\\sum_{q=-k}^{k}B_q^{(k)}\\mathcal{C}(i)_\
    q^{(k)} + \\epsilon",
      nE == 2,
        "\\hat{H}&=\\sum_{k=2,4,6}F^{(k)}\\hat{f}_k 
    +\\alpha \\hat{L}^2
    +\\beta \\,\\mathcal{C}\\left(\\mathcal{G}(2)\\right)
    +\\gamma \\,\\mathcal{C}\\left(\\mathcal{SO}(7)\\right)\\\\
    &\\quad + \\zeta \\sum_{i=1}^{n}\\left(\\hat{s}_i \\cdot \
    \\hat{l}_i\\right) 
    +\\sum_{k=0,2,4}M^{(k)}\\hat{m}_k
    +\\sum_{k=2,4,6}P^{(k)}\\hat{p}_k \\\\
    &\\quad\\quad\\quad\\quad\\quad+\\sum_{i=1}^n\\sum_{k=2,4,6}\\sum_{q=-\
    k}^{k}B_q^{(k)}\\mathcal{C}(i)_q^{(k)} + \\epsilon",
      And[nE == 12, OptionValue["T2"]],
        "\\hat{H}&=\\sum_{k=2,4,6}F^{(k)}\\hat{f}_k 
        +T^{(2)}\\hat{t}_2
        +\\alpha \\hat{L}^2
        +\\beta \\,\\mathcal{C}\\left(\\mathcal{G}(2)\\right)
        +\\gamma \\,\\mathcal{C}\\left(\\mathcal{SO}(7)\\right)\\\\
        &\\quad\\quad\\quad + \\zeta \\sum_{i=1}^{n}\\left(\\hat{s}_i \\cdot \
        \\hat{l}_i\\right) 
        +\\sum_{k=0,2,4}M^{(k)}\\hat{m}_k
        +\\sum_{k=2,4,6}P^{(k)}\\hat{p}_k \\\\
        &\\quad\\quad\\quad\\quad\\quad\\quad\\quad+\\sum_{i=1}^n\\sum_{k=2,4,\
        6}\\sum_{q=-k}^{k}B_q^{(k)}\\mathcal{C}(i)_q^{(k)} + \\epsilon",
      And[nE == 12, Not@OptionValue["T2"]],
        "\\hat{H}&=\\sum_{k=2,4,6}F^{(k)}\\hat{f}_k 
        +\\alpha \\hat{L}^2
        +\\beta \\,\\mathcal{C}\\left(\\mathcal{G}(2)\\right)
        +\\gamma \\,\\mathcal{C}\\left(\\mathcal{SO}(7)\\right)\\\\
        &\\quad + \\zeta \\sum_{i=1}^{n}\\left(\\hat{s}_i \\cdot \
        \\hat{l}_i\\right) 
        +\\sum_{k=0,2,4}M^{(k)}\\hat{m}_k
        +\\sum_{k=2,4,6}P^{(k)}\\hat{p}_k \\\\
        &\\quad\\quad\\quad+\\sum_{i=1}^n\\sum_{k=2,4,6}\\sum_{q=-k}^{k}B_q^{(\
    k)}\\mathcal{C}(i)_q^{(k)} + \\epsilon",
      True,
        "\\hat{H}&=\\sum_{k=2,4,6}F^{(k)}\\hat{f}_k 
        +\\sum_{k=2,3,4,6,7,8}T^{(k)}\\hat{t}_k
        +\\alpha \\hat{L}^2
        +\\beta \\,\\mathcal{C}\\left(\\mathcal{G}(2)\\right)
        +\\gamma \\,\\mathcal{C}\\left(\\mathcal{SO}(7)\\right)\\\\
        &\\quad\\quad\\quad\\quad + \\zeta \\sum_{i=1}^{n}\\left(\\hat{s}_i \
        \\cdot \\hat{l}_i\\right) 
        +\\sum_{k=0,2,4}M^{(k)}\\hat{m}_k
        +\\sum_{k=2,4,6}P^{(k)}\\hat{p}_k \\\\
        &\\quad\\quad\\quad\\quad\\quad\\quad\\quad\\quad+\\sum_{i=1}^n\\sum_{\
        k=2,4,6}\\sum_{q=-k}^{k}B_q^{(k)}\\mathcal{C}(i)_q^{(k)} + \\epsilon"
      ];
    Return[StringJoin[{"\\begin{aligned}\n", tex, "\n\\end{aligned}"}]]
    )

  EllipsoidBoundingBox::usage = "EllipsoidBoundingBox[A,\[Kappa]] gives the coordinate intervals that contain the ellipsoid determined by r^T.A.r==\[Kappa]^2. The matrix A must be square NxN, symmetric, and positive definite. The function returns a list with N pairs of numbers, each pair being of the form {-x_i, x_i}.";
  EllipsoidBoundingBox[Amat_,\[Kappa]_]:=Module[
    {invAmat, stretchFactors, boundingPlanes, quad},
      (
      invAmat = Inverse[Amat];
      stretchFactors = Sqrt[1/Diagonal[invAmat]];
      boundingPlanes = DiagonalMatrix[stretchFactors].invAmat;
      (* The solution is proportional to \[Kappa] *)
      boundingPlanes = \[Kappa] * boundingPlanes;
      boundingPlanes = Max /@ Transpose[boundingPlanes];
      Return[{-#,#}& /@ boundingPlanes]
    )
  ];

End[];
EndPackage[];