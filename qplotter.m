
BeginPackage["qplotter`"];

GetColor;
IndexMappingPlot;
ListLabelPlot;
AutoGraphicsGrid;
SpectrumPlot;
WaveToRGB;

Begin["`Private`"];

  AutoGraphicsGrid::usage="AutoGraphicsGrid[graphsList] takes a list of graphics and creates a GraphicsGrid with them. The number of columns and rows is chosen automatically so that the grid has a squarish shape.";
  Options[AutoGraphicsGrid] = Options[GraphicsGrid];
  AutoGraphicsGrid[graphsList_, opts : OptionsPattern[]] := 
    (
      numGraphs = Length[graphsList];
      width = Floor[Sqrt[numGraphs]];
      height = Ceiling[numGraphs/width];
      groupedGraphs = Partition[graphsList, width, width, 1, Null];
      GraphicsGrid[groupedGraphs, opts]
    )

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

  ListLabelPlot::usage="ListLabelPlot[data, labels] takes a list of numbers with corresponding labels. The data is grouped according to the labels and a ListPlot is created with them so that each group has a different color and their corresponding label is shown in the horizontal axis.";
  Options[ListLabelPlot] = Join[Options[ListPlot], {"TickCompression"->True,
  "LabelLevels"->1}];
  ListLabelPlot[data_, labels_, opts : OptionsPattern[]] := Module[
    {uniqueLabels, pallete, groupedByTerm, groupedKeys, scatterGroups, 
    groupedColors, frameTicks, compTicks, bottomTicks, topTicks},
    (
    uniqueLabels  = DeleteDuplicates[labels];
    pallete = Table[ColorData["Rainbow", i], {i, 0, 1, 
        1/(Length[uniqueLabels] - 1)}];
    uniqueLabels  = (#[[1]] -> #[[2]]) & /@ Transpose[{RandomSample[uniqueLabels], pallete}];
    uniqueLabels  = Association[uniqueLabels];
    groupedByTerm = GroupBy[Transpose[{labels, Range[Length[data]], data}], First];
    groupedKeys   = Keys[groupedByTerm];
    scatterGroups = Transpose[Transpose[#][[2 ;; 3]]] & /@ Values[groupedByTerm];
    groupedColors = uniqueLabels[#] & /@ groupedKeys;
    frameTicks    = {Transpose[{Range[Length[data]], 
      Style[Rotate[#, 90 Degree], uniqueLabels[#]] & /@ labels}], 
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
      Frame      -> True,
      AxesStyle  -> {Directive[Black, Dotted], Automatic},
      PlotStyle  -> groupedColors,
      FrameTicks -> frameTicks]
    )
    ]

  WaveToRGB::usage="WaveToRGB[wave, gamma] takes a wavelength in nm and returns the corresponding RGB color. The gamma parameter is optional and defaults to 0.8. The wavelength wave is assumed to be in nm. If the wavelength is below 380 the color will be the same as for 380 nm. If the wavelength is above 750 the color will be the same as for 750 nm. The function returns an RGBColor object. REF: https://www.noah.org/wiki/wave_to_rgb_in_Python. ";
  WaveToRGB[wave_, gamma_ : 0.8] := (
    wavelength = (wave);
    Which[
      wavelength < 380,
        wavelength = 380,
      wavelength > 750,
        wavelength = 750
    ];
    Which[380 <= wavelength <= 440,
    (
      attenuation = 0.3 + 0.7*(wavelength - 380)/(440 - 380);
      R = ((-(wavelength - 440)/(440 - 380))*attenuation)^gamma;
      G = 0.0;
      B = (1.0*attenuation)^gamma;
      ),
    440 <= wavelength <= 490,
    (
      R = 0.0;
      G = ((wavelength - 440)/(490 - 440))^gamma;
      B = 1.0;
      ),
    490 <= wavelength <= 510,
    (
      R = 0.0;
      G = 1.0;
      B = (-(wavelength - 510)/(510 - 490))^gamma;
      ),
    510 <= wavelength <= 580,
    (
      R = ((wavelength - 510)/(580 - 510))^gamma;
      G = 1.0;
      B = 0.0;
      ),
    580 <= wavelength <= 645,
    (
      R = 1.0;
      G = (-(wavelength - 645)/(645 - 580))^gamma;
      B = 0.0;
      ),
    645 <= wavelength <= 750,
    (
      attenuation = 0.3 + 0.7*(750 - wavelength)/(750 - 645);
      R = (1.0*attenuation)^gamma;
      G = 0.0;
      B = 0.0;
      ),
    True,
    (
      R = 0;
      G = 0;
      B = 0;
      )];
    Return[RGBColor[R, G, B]]
    )

  FuzzyRectangle::usage = "FuzzyRectangle[xCenter, width, ymin, height, color] creates a polygon with a fuzzy edge. The polygon is centered at xCenter and has a full horizontal width of width. The bottom of the polygon is at ymin and the height is height. The color of the polygon is color. The left edge and the right edge of the resulting polygon will be transparent and the middle will be colored. The polygon is returned as a list of polygons.";
  FuzzyRectangle[xCenter_, width_, ymin_, height_, color_, intensity_:1] := Module[
    {intenseColor, nocolor, ymax, polys},
  (
    nocolor = Directive[Opacity[0], color];
    ymax = ymin + height;
    intenseColor = Directive[Opacity[intensity], color];
    polys = {
      Polygon[{
        {xCenter - width/2, ymin},
        {xCenter, ymin},
        {xCenter, ymax},
        {xCenter - width/2, ymax}},
        VertexColors -> {
          nocolor,
          intenseColor,
          intenseColor,
          nocolor,
          nocolor}],
      Polygon[{
        {xCenter, ymin},
        {xCenter + width/2, ymin},
        {xCenter + width/2, ymax},
        {xCenter, ymax}},
        VertexColors -> {
          intenseColor,
          nocolor,
          nocolor,
          intenseColor,
          intenseColor}]
      };
    Return[polys]
    );
  ]

  Options[SpectrumPlot] = Options[Graphics];
  Options[SpectrumPlot] = Join[Options[SpectrumPlot], {"Intensities" -> {},"Tooltips" -> True, "Comments" -> {}, "SpectrumFunction" -> WaveToRGB}];
  SpectrumPlot::usage="SpectrumPlot[lines, widthToHeightAspect, lineWidth] takes a list of spectral lines and creates a visual representation of them. The lines are represented as fuzzy rectangles with a width of lineWidth and a height that is determined by the overall condition that the width to height ratio of the resulting graph is widthToHeightAspect. The color of the lines is determined by the wavelength of the line. The function assumes that the lines are given in nm.
  If the lineWidth parameter is a single number, then every line shares that width. If the lineWidth parameter is a list of numbers, then each line has a different width. The function returns a Graphics object. The function also accepts any options that Graphics accepts. The background of the plot is black by default. The plot range is set to the minimum and maximum wavelength of the given lines.
  Besides the options for Graphics the function also admits the option Intensities. This option is a list of numbers that determines the intensity of each line. If the Intensities option is not given, then the lines are drawn with full intensity. If the Intensities option is given, then the lines are drawn with the given intensity. The intensity is a number between 0 and 1.
  The function also admits the option \"Tooltips\". If this option is set to True, then the lines will have a tooltip that shows the wavelength of the line. If this option is set to False, then the lines will not have a tooltip. The default value for this option is True.
  If \"Tooltips\" is set to True and the option \"Comments\" is a non-empty list, then the tooltip will append the wavelength and the values in the comments list for the tooltips.
  The function also admits the option \"SpectrumFunction\". This option is a function that takes a wavelength and returns a color. The default value for this option is WaveToRGB.
  ";
  SpectrumPlot[lines_, widthToHeightAspect_, lineWidth_, opts : OptionsPattern[]] := Module[
    {minWave, maxWave, height, fuzzyLines},
  (
    colorFun = OptionValue["SpectrumFunction"];
    {minWave, maxWave} = MinMax[lines];
    height     = (maxWave - minWave)/widthToHeightAspect;
    fuzzyLines = Which[
      NumberQ[lineWidth] && Length[OptionValue["Intensities"]] == 0,
        FuzzyRectangle[#, lineWidth, 0, height, colorFun[#]] &/@ lines,
      Not[NumberQ[lineWidth]] && Length[OptionValue["Intensities"]] == 0,
        MapThread[FuzzyRectangle[#1, #2, 0, height, colorFun[#1]] &, {lines, lineWidth}],
      NumberQ[lineWidth] && Length[OptionValue["Intensities"]] > 0,
        MapThread[FuzzyRectangle[#1, lineWidth, 0, height, colorFun[#1], #2] &, {lines, OptionValue["Intensities"]}],
      Not[NumberQ[lineWidth]] && Length[OptionValue["Intensities"]] > 0,
        MapThread[FuzzyRectangle[#1, #2, 0, height, colorFun[#1], #3] &, {lines, lineWidth, OptionValue["Intensities"]}]
    ];
    comments = Which[
      Length[OptionValue["Comments"]] > 0, 
        MapThread[StringJoin[ToString[#1]<>" nm","\n",ToString[#2]]&,
        {lines, OptionValue["Comments"]}],
      Length[OptionValue["Comments"]] == 0, 
        ToString[#] <>" nm" & /@ lines,
      True,
      {}
    ];
    If[OptionValue["Tooltips"],
      fuzzyLines = MapThread[Tooltip[#1, #2] &, {fuzzyLines, comments}];
    ];
    graphicsOpts = FilterRules[{opts}, Options[Graphics]];
    Graphics[fuzzyLines, 
        graphicsOpts,
        Background -> Black,
        PlotRange -> {{minWave, maxWave}, {0, height}}]
    )
  ];

End[];

EndPackage[];