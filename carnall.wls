#!/usr/bin/env wolframscript

(* This script parses the data in Carnall.xls to produce Carnall.m*)

ions={"Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm"};
templates=StringTemplate/@StringSplit["appendix:`ion`:Association appendix:`ion`:Calculated appendix:`ion`:RawTable appendix:`ion`:Headings"," "];

(* How many unique eigenvalues, after removing Kramer's degeneracy *)
fullSizes = AssociationThread[ions, 
    {91, 182, 1001, 1001, 3003, 1716, 3003, 1001, 1001, 182, 91}]

moduleDir = DirectoryName[$InputFileName];
carnall=Import[FileNameJoin[{moduleDir,"data","Carnall.xls"}]][[2]];
carnallErr=Import[FileNameJoin[{moduleDir,"data","Carnall.xls"}]][[3]];

elementNames=carnall[[1]][[2;;]];
carnall=carnall[[2;;]];
carnallErr=carnallErr[[2;;]];
carnall=Transpose[carnall];
carnallErr=Transpose[carnallErr];
paramNames=ToExpression/@carnall[[1]][[1;;]];
carnall=carnall[[2;;]];
carnallErr=carnallErr[[2;;]];
carnallData=Table[(
data=carnall[[i]];
data=(#[[1]]->#[[2]])&/@Select[Transpose[{paramNames,data}],#[[2]]!=""&];
elementNames[[i]]->data),{i,1,13}];
carnallData=Association[carnallData];
carnallNotes=Table[(
    data=carnallErr[[i]];
    elementName=elementNames[[i]];
    dataFun=(
        #[[1]] -> If[#[[2]]=="[]",
        "Not allowed to vary in fitting.",
        If[#[[2]]=="[R]",
            "Ratio constrained by: " <> <|"Eu"->"F4/F2=0.713; F6/F2=0.512",
                "Gd"->"F4/F2=0.710]",
                "Tb"->"F4/F2=0.707"|>[elementName],
            If[#[[2]]=="i",
                "Interpolated",
                #[[2]]
            ]
        ]
        ]) &;
    data = dataFun /@ Select[Transpose[{paramNames,data}],#[[2]]!=""&];
    elementName->data
    ),
{i,1,13}
];
carnallNotes=Association[carnallNotes];

annotatedData=Table[If[NumberQ[#[[1]]],Tooltip[#[[1]],#[[2]]],""] & /@ Transpose[{paramNames/.carnallData[element],
           paramNames/.carnallNotes[element]
    }],
{element,elementNames}
];
annotatedData=Transpose[annotatedData];

Carnall=<|"data"->carnallData,
    "annotations"->carnallNotes,
    "paramSymbols"->paramNames,
    "elementNames"->elementNames,
    "rawData"->carnall,
    "rawAnnotations"->carnallErr,
    "includedTableIons"->ions,
    "annnotatedData"->annotatedData
|>;

Do[(
    carnallData=Import[FileNameJoin[{moduleDir,"data","Carnall.xls"}]][[i]];
    headers=carnallData[[1]];
    calcIndex=Position[headers,"Calc (1/cm)"][[1,1]];
    headers=headers[[2;;]];
    carnallLabels=carnallData[[1]];
    carnallData=carnallData[[2;;]];
    carnallTerms=DeleteDuplicates[First/@carnallData];
    parsedData=Table[(
        rows=Select[carnallData,#[[1]]==term&];
        rows=#[[2;;]]&/@rows;
        rows=Transpose[rows];
        rows=Transpose[{headers,rows}];
        rows=Association[(#[[1]]->#[[2]])&/@rows];
        term->rows
        ),
    {term,carnallTerms}
    ];
    carnallAssoc=Association[parsedData];
    carnallCalcEnergies=#[[calcIndex]]&/@carnallData;
    carnallCalcEnergies=If[NumberQ[#],#,Missing[]]&/@carnallCalcEnergies;
    ion=ions[[i-3]];
    carnallCalcEnergies=PadRight[carnallCalcEnergies, fullSizes[ion], Missing[]];
    keys=#[<|"ion"->ion|>]&/@templates;
    Carnall[keys[[1]]]=carnallAssoc;
    Carnall[keys[[2]]]=carnallCalcEnergies;
    Carnall[keys[[3]]]=carnallData;
    Carnall[keys[[4]]]=headers;
    ),
{i,4,14}
];

goodions = Select[ions,#!="Pm"&];
expData  = Select[Transpose[Carnall["appendix:"<>#<>":RawTable"]][[1+Position[Carnall["appendix:"<>#<>":Headings"],"Exp (1/cm)"][[1,1]]]],NumberQ]&/@goodions;
Carnall["All Experimental Data"]=AssociationThread[goodions,expData];

exportFname = FileNameJoin[{moduleDir,"data","Carnall.m"}];
Print["Exporting to "<>exportFname];
Export[exportFname,Carnall];
