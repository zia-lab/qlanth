
BeginPackage["qplotter`"];

GetColor;
IndexMappingPlot;
ListLabelPlot::usage="ListLabelPlot[data, labels] takes a list of numbers with corresponding labels. The data is grouped according to the labels and a ListPlot is created with them so that each group has a different color and their corresponding label is shown in the horizontal axis.";

Begin["`Private`"];

Options[IndexMappingPlot] = Options[Graphics];
IndexMappingPlot::usage = 
  "IndexMappingPlot[pairs] take a list of pairs of integers and creates a visual representation of how they are paired. The first indices being depicted in the bottom and the second indices being depicted on top.";
IndexMappingPlot[pairs_, opts : OptionsPattern[]] := Module[{width, height}, (
   width = Max[First /@ pairs];
   height = width/3;
   Return[
    Graphics[{{Tooltip[Point[{#[[1]], 0}],#[[1]]]}, Tooltip[Point[{#[[2]], height}],#[[2]]], 
        Line[{{#[[1]], 0}, {#[[2]], height}}]} & /@ pairs, opts, ImageSize -> 800]]
   )
  ]

TickCompressor[fTicks_] := 
 Module[{avgTicks, prevTickLabel, groupCounter, groupTally, idx, 
  tickPosition, tickLabel, avgPosition, groupLabel}, (avgTicks = {};
  prevTickLabel = fTicks[[1, 2]];
  groupCounter = 0;
  groupTally = 0;
  idx = 1;
  Do[({tickPosition, tickLabel} = tick;
    If[
     tickLabel === prevTickLabel,
     (groupCounter += 1;
      groupTally += tickPosition;
      groupLabel = tickLabel;),
     (
      avgPosition = groupTally/groupCounter;
      avgTicks = Append[avgTicks, {avgPosition, groupLabel}];
      groupCounter = 1;
      groupTally = tickPosition;
      groupLabel = tickLabel;
      )
     ];
    If[idx != Length[fTicks],
     prevTickLabel = tickLabel;
     idx += 1;]
    ), {tick, fTicks}];
  If[Or[Not[prevTickLabel === tickLabel], groupCounter > 1],
   (
    avgPosition = groupTally/groupCounter;
    avgTicks = Append[avgTicks, {avgPosition, groupLabel}];
    )
   ];
  Return[avgTicks];)]

GetColor[s_Style] := s /. Style[_, c_] :> c
GetColor[_] := Black

Options[ListLabelPlot] = Append[Options[ListPlot], "TickCompression"->True];
ListLabelPlot[data_, labels_, opts : OptionsPattern[]] := Module[
  {uniqueLabels, pallete, groupedByTerm, groupedKeys, scatterGroups, 
   groupedColors, frameTicks, compTicks, bottomTicks, topTicks},
  (
   uniqueLabels  = DeleteDuplicates[labels];
   pallete = Table[ColorData["Rainbow", i], {i, 0, 1, 
      1/(Length[uniqueLabels] - 1)}];
   uniqueLabels  = (#[[1]] -> #[[2]]) & /@ 
     Transpose[{RandomSample[uniqueLabels], pallete}];
   uniqueLabels  = Association[uniqueLabels];
   groupedByTerm = GroupBy[Transpose[{labels, Range[Length[data]], data}], First];
   groupedKeys   = Keys[groupedByTerm];
   scatterGroups = Transpose[Transpose[#][[2 ;; 3]]] & /@ Values[groupedByTerm];
   groupedColors = uniqueLabels[#] & /@ groupedKeys;
   frameTicks    = {Transpose[{Range[Length[data]], 
    Style[Rotate[#, Pi/2], uniqueLabels[#]] & /@ labels}], 
     Automatic};
    If[OptionValue["TickCompression"], (
        compTicks = TickCompressor[frameTicks[[1]]];
        bottomTicks = 
            MapIndexed[
            If[EvenQ[First[#2]], {#1[[1]], 
                Tooltip[Style["\[SmallCircle]", GetColor[#1[[2]]]],#1[[2]]]
                }, #1] &, compTicks];
        topTicks = 
            MapIndexed[
            If[OddQ[First[#2]], {#1[[1]], 
                Tooltip[Style["\[SmallCircle]", GetColor[#1[[2]]]],#1[[2]]]
                }, #1] &, compTicks];
        frameTicks = {{Automatic, Automatic}, {bottomTicks, topTicks}};)
   ];
   ListPlot[scatterGroups,
    opts,
    Frame->True,
    PlotStyle -> groupedColors,
    FrameTicks -> frameTicks]
   )
  ]

End[];

EndPackage[];