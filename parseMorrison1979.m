#!/usr/bin/env wolfram

scriptDir = DirectoryName[$InputFileName];
AppendTo[$Path, scriptDir];
SetDirectory[scriptDir];

cfVars = {B02, B04, B06, B22, B24, B26, B44, B46, B66};
data = Import["./data/Morrison - 1979.xls"][[4]];
numDataRows = Position[Transpose[data][[1]], ""][[1, 1]] - 1;
data = data[[;; numDataRows]];
header = data[[1]];
parsingHeader = header[[4 ;; -2]];
parsingHeader = (If[MemberQ[StringSplit[#, "-"], "Im"], I, 1]*
      ToExpression[
       StringJoin[
        Permute[Characters[StringSplit[#, "-"][[1]]], {1, 3, 
          2}]]]) & /@ parsingHeader;
clyde = Association[
   Table[(data[[idx]][[1]] -> 
      AssociationThread[
       cfVars -> 
        Coefficient[data[[idx]][[4 ;; -2]] . parsingHeader, 
         cfVars]]), {idx, 2, numDataRows}]
   ];
Export["./data/morrison1979.m", clyde]
