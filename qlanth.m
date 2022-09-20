(* -----------------------------------------------------------------
This  code was authored by Chris Dodson and turned into a package by
David Lizarazo in 2022.
+------------------------------------------------------------------+
|                                                                  |
|                      __                     __     __            |
|            ____ _   / /  ____ _   ____     / /_   / /_           |
|           / __ `/  / /  / __ `/  / __ \   / __/  / __ \          |
|          / /_/ /  / /  / /_/ /  / / / /  / /_   / / / /          |
|          \__, /  /_/   \__,_/  /_/ /_/   \__/  /_/ /_/           |
|            /_/                                                   |
|                                                                  |
|                                                                  |
+------------------------------------------------------------------+
It  implements  the  solution  of  several  Hamiltonians that are of
interest  in  understanding  the  electronic structure of lanthanide
ions in crystals.

REFERENCES:
1. Velkov, Multi-Electron Coefficients of Fractional Parentage for the p, d, and f Shells.
2. Nielson and Koster, Spectroscopic Coefficients for the Pn, Dn, and Fn Configurations.
3. Wybourne and Smentek, Optical Spectroscopy of Lanthanides.
4. 
----------------------------------------------------------------- *)

BeginPackage["qlanth`"]

(* Physical Constants*)
bohrRadius = 5.29177210903 10^-9;
ee = 1.602176634 * 10^-9;

specAlphabet = "SPDFGHIKLMNOQRTUV";
moduleDir = DirectoryName[$InputFileName];

theLanthanides = {"Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb"};

printL::usage="printL[L] give the string representation of a given angular momentum."
findSL::usage = "findSL[LS] gives the spin and orbital angular momentum that corresponds to the provided string LS."
findSL\[Gamma]::usage = "findSL\[Gamma][LS] gives the spin and orbital angular momentum that corresponds to the provided string LS. The provided lists contains terms that belong to one of the f^n configurations. The last element of this list provides a digit that differentiates between different possibly degenerate terms."
printSLJ::usage = "Given a list with three elements {S,L,J} this function returns a symbol where the spin multiplicity is presented as a superscript, the orbital angular momentum as its corresponding spectroscopic letter, and J as a subscript. Function does not check to see if the given J is compatible with the given S and L."
printSLJM::usage = "Given a list with four elements {S,L,J,MJ} this function returns a symbol where the spin multiplicity is presented as a superscript, the orbital angular momentum as its corresponding spectroscopic letter, and {J, MJ} as a subscript. No attempt is made to guarantee that the given input is consistent."

fnTerms::usage = "This list contains the labels of f^n configurations. Each element of the list has four elements {LS, seniority, RacahW, RacahU}."
NKLSterms::usage = "NKLSterms[n, L, S] give all the terms that correspond to the given number of electrons n, with total orbital angular momentum L, and total spin angular momentum S."
findNKLSterm::usage = "findNKLSterm[SL] for a given string SL that only represent the LS part of a spectroscopic term, this function returns all the terms that are compatible with it. This is only for f^n configurations. The provided terms might belong to more than one configuration."
AA::usage = "Given a list with the parts corresponding to the quasi-rational representation of a number, this function parses the result into a regular number."
CFPterms::usage = "CFPterms[n] gives all the daughter and parent terms, together with the corresponding coefficients of fractional parentage, that correspond to the the f^n configuration. 
CFPterms[n, SL] gives all the daughter and parent terms, together with the corresponding coefficients of fractional parentage, that are compatible with the given string SL in the f^n configuration 
CFPterms[n, L, S] gives all the daughter and parent terms, together with the corresponding coefficients of fractional parentage, that correspond to the given total orbital angular momentum L and total spin S n the f^n configuration. L being an integer, and S being integer or half-integer.
In all cases the output is in the shape of a list with enclosed lists having the format {daughter_term, {parent_term_1, CFP_1}, {parent_term_2, CFP_2}, ...}.
Only the one-body coefficients for f-electrons are provided.
In all cases it must be that 1 <= n <= 7.
Data for these was parsed from Velkov, \"Multi-electron coefficients of fractional parentage for the p, d, and f shells\".
"

AllowedSLterms::usage = "AllowedSLterms[n] returns a list with the allowed terms in the f^n configuration, the terms are given as lists in the format {S,L}. This list may have redundancies which are compatible with the degeneracies that might correspond to the given case."
AllowedNKSLterms::usage = "AllowedNKSLterms[n] returns a list with the allowed terms in the f^n configuration, the terms are given as strings in spectroscopic notation. An integer as the last is used to distinguish cases with degenaracy."
maxJ::usage = "maxJ[n] gives the maximum J=S+L that corresponds to the configuration f^n."
minJ::usage = "minJ[n] gives the minimum J=S+L that corresponds to the configuration f^n."
AllowedSLJterms::usage = "AllowedSLJterms[n] returns a list with the allowed {S,L,J} terms in the f^n configuration, the terms are given as lists in the format {S,L,J}. This list may have redundancies which are compatible with the degeneracies that might correspond to the given case."

AllowedNKSLJterms::usage = "AllowedNKSLJterms[n] returns a list with the allowed {SL,J} terms in the f^n configuration, the terms are given as lists in the format {SL,J} where SL is a string in spectroscopic notation."
AllowedNKSLforJterms::usage = "AllowedNKSLforJterms[n, J] gives the terms that correspond to the given total angular momentum J in the f^n configuration. The result is a list whose elements are lists of length 2, the first element being the SL term in spectroscopic notation, and the second element being J."
AllowedMforJ::usage = "AllowedMforJ[j] is shorthand for Range[-j,j,1]."
AllowedNKSLJMIMforJIterms::usage = "AllowedNKSLJMIMforJIterms[n, J, I] returns the states that belong to the f^n configuration with a total angular momentun L+S equal to J and with a nuclear angular momentum I. The returned list has elements of the form {{{SL (string in spectroscopic notation), J}, MJ}, MI}."
AllowedJ::usage = "AllowedJ[n] returns the total angular momenta J that appear in the f^n configuration."

AllowedSLJMterms::usage = "AllowedSLJMterms[n] returns a list with all the states that correspond to the configuration f^n. A list is returned whose elements are lists of the form {S, L, J, MJ}."
AllowedNKSLJMforJMterms::usage = "AllowedNKSLJMforJMterms[n, J, MJ] returns a list with all the terms that contain states of the f^n configuration that have a total angular momentum J, and a projection along the z-axis MJ. The returned list has elements of the form {SL (string in spectroscopic notation), J, MJ}."
AllowedNKSLJMforJterms::usage = "AllowedNKSLJMforJterms[n, J] returns a list with all the states that have a total angular momentum J. The returned list has elements of the form {{SL (string in spectroscopic notation), J}, MJ}."

ShiftedLevels::usage = "
ShiftedLevels[originalLevels] takes a list of levels of the form
{{energy_1, coeff_vector_1, basis}, 
{energy_2, coeff_vector_2, basis},
...}} 
and returns the same input except that now to every energy the minimum of all of them has been subtracted.";
LoadGuillotParameters::usage = "";

CFP::usage = "CFP[{n, NKSL}] provides a list whose first element echoes NKSL and whose other elements are lists with two elements the first one being the symbol of a parent term and the second being the corresponding coefficient of fractional parentage. n must satisfy 1 <= n <=7"

(* ############################################################################################## *)
(* ############################################ Data ############################################ *)

Print["Loading coefficients of fractional parentage ..."]
CFPfname = FileNameJoin[{moduleDir,"data","CFPs.m"}]
CFP = Import[CFPfname]
CFP::usage = "";

Print["Loading reduced matrix elements for unit tensor operators ..."]
ReducedUktablefname = FileNameJoin[{moduleDir,"data","ReducedUktables.m"}]
ReducedUktable = Import[ReducedUktablefname]
ReducedUktable::usage="";

Print["Loading matrix elements for the electrostatic interaction ..."]
ElectrostaticMatrixTablefname = FileNameJoin[{moduleDir,"data","ElectrostaticMatrixTable.m"}]
ElectrostaticMatrixTable = Import[ElectrostaticMatrixTablefname]
ElectrostaticMatrix::usage="";

(* ############################################ Data ############################################ *)
(* ############################################################################################## *)

(* ############################################################################################## *)
(* ################################### Racah Algebra Goodies #################################### *)

ReducedUk::usage="ReducedUk[n, l, SL, SpLp,k] gives the reduced matrix element of the spherical tensor U^(k).";
ReducedUk[n_,l_,SL_,SpLp_,k_]:=Module[
  {Uk,S,L,Sp,Lp,Sb,Lb,parentSL,cfpSL,cfpSpLp},
  {S,L}=findSL[SL];
  {Sp,Lp}=findSL[SpLp];
  cfpSL=CFP[{n,SL}];
  cfpSpLp=CFP[{n,SpLp}];
  Uk=0;
  If[S==Sp,
    For[nn=2,nn<=Length[cfpSL],nn++,
      parentSL=Part[cfpSL,nn,1];
      {Sb,Lb}=findSL[parentSL];
      fun=If[
        #[[1]]!=parentSL,
        0,
        ((-1)^(Lb+L+l+k)
        * #[[2]]
        * cfpSL[[nn,2]]
        * SixJSymbol[{L,Lp,k},{l,l,Lb}])]&;
      summands=Map[fun,cfpSpLp[[2;;]]];
      Uk=Uk+Total[summands];
    ];
    n * Sqrt[(2L+1) (2Lp+1)] * Uk,
  0]
]

Ck::usage="Diagonal reduced matrix element <l||C^(k)\[VerticalSeparator]\[VerticalSeparator]l> where the Subscript[C, q]^(k) are reduced spherical harmonics.";
Ck[l_,k_]:=(-1)^l *(2 l+1) * ThreeJSymbol[{l,0},{l,0},{k,0}]

fk::usage="Slater integral. See Wybourne (1965) eqn. 2-93.";
fk[n_,l_,NKSL_,NKSLp_,k_]:=Module[
  {S,L,PossibleStates,Sb,Lb},
  {S,L}=findSL[NKSL];
  PossibleStates=StringCases[AllowedNKSLterms[n],ToString[2 S+1]~~__];
  PossibleStates=DeleteCases[PossibleStates,{}];
  PossibleStates=Flatten[PossibleStates];
  PossibleStates=
  If[{S,L} == findSL[NKSLp],
    (1/2 ( Abs[Ck[l,k]])^2 * (1/(2L+1) Sum[(ReducedUktable[{n,l,x,NKSL,k}]*ReducedUktable[{n,l,x,NKSLp,k}]),{x,PossibleStates}]
     - KroneckerDelta[NKSL,NKSLp] (n (4 l+2-n))/((2 l+1) (4 l+1)))),
  0]
]

fK::usage="Non-reduced Slater integral f^k = Subscript[f, k] * Subscript[D, k]";
fK[n_,l_,NKSL_,NKSLp_,k_]:= Dk[k] * fk[n,l,NKSL,NKSLp,k]

Dk::usage="Ratio between the reduced and non-reduced Slater direct (Subscript[F, k] and F^k) and exchange(Subscript[G, k]and G^k) integrals. Subscript[D, k]:= (Subscript[G, k](ff))/(G^k (ff))=(Subscript[F, k](ff))/(F^k (ff)). k must be even. See table 6-3 in Cowan (1981), and also section 2-7 of Wybourne (1965).";
Dk[k_]:={1,225,1089,184041/25}[[k/2+1]]

ElectrostaticMatrix::usage="See Wybourne (1965) section 2-9 'Electrostatic Interactions in f^n Configurations'"
ElectrostaticMatrix[n_,NKSL_,NKSLp_]:=Module[{f0,f2,f4,f6,e0,e1,e2,e3,eevals,EMatrix,l},
  l=3;
  Ek={E0,E1,E2,E3};
  f0=fK[n,l,NKSL,NKSLp,0];
  f2=fK[n,l,NKSL,NKSLp,2];
  f4=fK[n,l,NKSL,NKSLp,4];
  f6=fK[n,l,NKSL,NKSLp,6];
  e0=f0;
  e1=9/7 f0+1/42 f2+1/77 f4+1/462 f6;
  e2=143/42 f2-130/77 f4+35/462 f6;
  e3=11/42 f2+4/77 f4-7/462 f6;
  eevals={e0,e1,e2,e3};
  EMatrix=e0 E0+e1 E1+e2 E2+e3 E3
]

(* ################################### Racah Algebra Goodies #################################### *)
(* ############################################################################################## *)


(* ############################################################################################## *)
(* ########################################## Printers ########################################## *)

printL[L_] := 
 If[StringQ[L], L, StringTake[specAlphabet, {L + 1}]]

findSL[SL_] := 
 If[StringQ[SL], {(ToExpression[StringTake[SL, 1]] - 1)/2, 
   Part[StringPosition[specAlphabet, StringTake[SL, {2}]], 1, 
     1] - 1}, SL]

findSL\[Gamma][SL_] := 
 If[StringQ[SL], 
  If[StringLength[SL] == 2, 
    {
    (ToExpression[StringTake[SL, 1]] - 1)/2, 
    Part[StringPosition[specAlphabet, StringTake[SL, {2}]], 1, 1] - 1,
    0}, 
    {
    (ToExpression[StringTake[SL, 1]] - 1)/2, 
    Part[StringPosition[specAlphabet, StringTake[SL, {2}]], 1, 1] - 1,
    ToExpression[StringTake[SL, {3}]]}], SL]

printSLJ[SLJ_] := 
  RowBox[{SuperscriptBox[" ", 2 SLJ[[1]] + 1], 
     SubscriptBox[printL[SLJ[[2]]], SLJ[[3]]]}] // DisplayForm;

printSLJM[SLJM_] := 
  RowBox[{SuperscriptBox[" ", 2 SLJM[[1]] + 1], 
     SubscriptBox[printL[SLJM[[2]]], {SLJM[[3]], SLJM[[4]]}]}] // 
   DisplayForm;

(* ########################################## Printers ########################################## *)
(* ############################################################################################## *)

(* ############################################################################################## *)
(* ############################ Coefficients of Fractional Parentage ############################ *)

LabelsTextData = 
  FileNameJoin[{moduleDir, "data", "NielsonKosterLabels_f6_f7.txt"}];
fNtextLabels = Import[LabelsTextData];

NielsonKosterLabels = Partition[StringSplit[fNtextLabels], 3];
TermLabels = Map[Part[#, {1}] &, NielsonKosterLabels];
Seniority = Map[Part[# , {2}] &, NielsonKosterLabels];
RacahW = 
  Map[StringTake[
     Flatten[StringCases[Part[# , {3}], 
       "(" ~~ DigitCharacter ~~ DigitCharacter ~~ DigitCharacter ~~ 
        ")"]], {2, 4}] &, NielsonKosterLabels];
RacahU = 
  Map[StringTake[
     Flatten[StringCases[Part[# , {3}], 
       "(" ~~ DigitCharacter ~~ DigitCharacter ~~ ")"]], {2, 3}] &, 
   NielsonKosterLabels];


fNterms = Join[TermLabels, Seniority, RacahW, RacahU, 2];
fNterms = Sort[fNterms];
Clear[fNtextLabels, TermLabels, Seniority, RacahU, RacahW]

NKLSterms[n_, L_, S_] := Module[{NKterms, SL},
  SL = StringJoin[ToString[2 S + 1], printL[L]];
  NKterms = {{}};
  Map[If[! StringFreeQ[First[#], SL], 
     If[ToExpression[Part[#, 2]] <= n, 
      NKterms = Join[NKterms, {#}, 1]]] &, fNterms];
  NKterms = DeleteCases[NKterms, {}];
  If[Length[NKterms] == 1, 
   NKterms[[1]][[1]] = StringTake[Part[NKterms, 1, 1], 2]];
  NKterms]

findNKLSterm[SL_] := Module[{NKterms, n},
  n = 7;
  NKterms = {{}};
  Map[If[! StringFreeQ[First[#], SL], 
     If[ToExpression[Part[#, 2]] <= n, 
        NKterms = Join[NKterms, {#}, 1]]] &,
     fNterms];
  NKterms = DeleteCases[NKterms, {}];
  NKterms]

LabelsTextData = 
  FileNameJoin[{moduleDir, "data", "NielsonKosterLabels_f6_f7.txt"}];
fNtextLabels = Import[LabelsTextData];

(* ############################################################################################## *)
(* ############################ Coefficients of Fractional Parentage ############################ *)


CFPtextData = Import[FileNameJoin[{moduleDir, "data", "B1F_ALL.TXT"}]];
fConfigs = StringSplit[CFPtextData, 
  "[ONE PARTICLE FRACTIONAL PARENTAGE COEFFICIENTS "];
CFPtable = {};

AA[quasinum_] :=
 (radical = Product[Prime[nn - 1] ^ Part[quasinum, nn], {nn, 2, Length[quasinum]}];
 (Part[quasinum, 1] * Sqrt[radical]) /. {A -> 10, B -> 11, C -> 12, D -> 13})

(* This table parses the text file with the one-body coefficients of fractional parentage *)

CFPtable = Table[
 (fNx = StringReplace[Part[fConfigs, nn], "-" -> " -"];
 daughterLabelSpots = 
  StringPosition[fNx, 
   Shortest[
    StartOfLine ~~ DigitCharacter ~~ LetterCharacter ~~ ___ ~~ "["], 
   Overlaps -> False];
 daughterLabels = Map[StringDrop[#, -1] &, StringTake[fNx, daughterLabelSpots]];
 daughterLabelLines = 
  StringPosition[fNx, 
   Shortest[
    StartOfLine ~~ DigitCharacter ~~ LetterCharacter ~~ __ ~~ 
     EndOfLine], Overlaps -> False];
 startDaughters = Map[Last, daughterLabelLines + 2];
 stopDaughters = 
  Delete[Append[Map[First, daughterLabelLines - 2], 
    StringLength[fNx]], 1];
 daughterLines = 
  Join[Partition[startDaughters, 1], Partition[stopDaughters, 1], 2];
 testing = 
  Map[StringSplit, 
   StringSplit[StringTake[fNx, daughterLines], EndOfLine]];
 testing2 = Map[DeleteCases[#, {}] &, testing];
 CFPs = {};
 For[mm = 1, mm <= Length[testing2], mm++,
  tt = Part[testing2, mm];
  pLabels = Map[Extract[#, 1] &, tt];
  pValues = Map[AA, ToExpression[Map[Drop[#, 2] &, tt]]];
  CFPs = 
   Append[CFPs, Join[Partition[pLabels, 1], Partition[pValues, 1], 2]];
  ];
 CFPconfig = Join[Partition[daughterLabels, 1], CFPs, 2];
CFPconfig)
 , {nn, 2, 7}
 ];

CFPtable = Join[{{{"2F", {"1S", 1}}}}, CFPtable]; (*for f^1*)

CFPterms[n_] := Part[CFPtable, n]
CFPterms[n_, SL_] := Module[{NKterms, CFPconfig},
  NKterms = {{}};
  CFPconfig = Part[CFPtable, n];
  Map[If[StringFreeQ[First[#], SL], Null , 
     NKterms = Join[NKterms, {#}, 1]] &, CFPconfig];
  NKterms = DeleteCases[NKterms, {}]]
CFPterms[n_, L_, S_] := Module[{NKterms, SL, CFPconfig},
  SL = StringJoin[ToString[2 S + 1], printL[L]];
  NKterms = {{}};
  CFPconfig = Part[CFPtable, n];
  Map[If[StringFreeQ[First[#], SL], Null , 
     NKterms = Join[NKterms, {#}, 1]] &, CFPconfig];
  NKterms = DeleteCases[NKterms, {}]]

(* ############################ Coefficients of Fractional Parentage ############################ *)
(* ############################################################################################## *)

(* ############################################################################################## *)
(* ############################### Allowed SL, SLJ, and SLJM Terms ############################## *)


AllowedSLterms[n_] := Map[findSL[First[#]] &, CFPterms[n]]

AllowedNKSLterms[n_] := Map[First, CFPterms[n]]

maxJ[n_] := Max[Map[Total, AllowedSLterms[n]]]

minJ[n_] := Min[Map[Abs[Part[#, 1] - Part[#, 2]] &, AllowedSLterms[n]]]

AllowedSLJterms[n_] := Module[{allowedSL, allowedSLJ},
    allowedSL = AllowedSLterms[n];
    allowedSLJ = {}; 
    For[nn = 1, nn <= Length[allowedSL],
      termSL = allowedSL[[ nn]];
      termsSLJ = 
    Table[{termSL[[1]], termSL[[2]], J}, {J, 
      Abs[termSL[[1]] - termSL[[2]]], Total[termSL]}];
      allowedSLJ = Join[allowedSLJ, termsSLJ];
      nn++];
    SortBy[allowedSLJ, Last]]

AllowedNKSLJterms[n_] := Module[{allowedSL, allowedNKSL, allowedSLJ},
    allowedNKSL = AllowedNKSLterms[n];
    allowedSL = AllowedSLterms[n];
    allowedSLJ = {}; 
    For[nn = 1, nn <= Length[allowedSL],
      termSL = allowedSL[[nn]];
      termNKSL = allowedNKSL[[nn]];
      termsSLJ = Table[{termNKSL, J}, 
          {J, Abs[termSL[[1]] - termSL[[2]]], Total[termSL]}];
        allowedSLJ = Join[allowedSLJ, termsSLJ];
        nn++];
    SortBy[allowedSLJ, Last]]

AllowedNKSLforJterms[n_, J_] := 
 Module[{allowedSL, allowedNKSL, allowedSLJ},
    allowedNKSL = AllowedNKSLterms[n];
    allowedSL = AllowedSLterms[n];
    allowedSLJ = {}; 
    For[nn = 1, nn <= Length[allowedSL],
      termSL = allowedSL[[ nn]];
      termNKSL = allowedNKSL[[ nn]];
      termsSLJ = 
    If[Abs[termSL[[1]] - termSL[[2]]] <= J <= 
      Total[termSL], {{termNKSL, J}}, {}];
      allowedSLJ = Join[allowedSLJ, termsSLJ];
      nn++];
    allowedSLJ]

AllowedSLJMterms[n_] := 
 Module[{allowedSLJ, allowedSLJM, termSLJ, termsSLJM},
    allowedSLJ = AllowedSLJterms[n];
    allowedSLJM = {}; 
    For[nn = 1, nn <= Length[allowedSLJ], nn++,
      termSLJ = allowedSLJ[[ nn]];
      termsSLJM = 
    Table[{termSLJ[[1]], termSLJ[[2]], termSLJ[[3]], 
      M}, {M, - termSLJ[[3]], termSLJ[[3]]}];
      allowedSLJM = Join[allowedSLJM, termsSLJM];];
    SortBy[allowedSLJM, Last]]

AllowedNKSLJMforJMterms[n_, J_, MJ_] :=  
 Module[{allowedSL, allowedNKSL, allowedSLJM},
    allowedNKSL = AllowedNKSLterms[n];
  allowedSL = AllowedSLterms[n];
  allowedSLJM = {}; 
    For[nn = 1, nn <= Length[allowedSL],
      termSL = allowedSL[[ nn]];
      termNKSL = allowedNKSL[[ nn]];
      termsSLJ = 
    If[Abs[termSL[[1]] - termSL[[2]]] <= J <= Total[termSL] && 
      Abs[MJ] <= J, {{termNKSL, J, MJ}}, {}];
      allowedSLJM = Join[allowedSLJM, termsSLJ];
      nn++];
    allowedSLJM]

AllowedNKSLJMforJterms[n_, J_] :=  
 Module[{allowedSL, allowedNKSL, allowedSLJM, Mvals},
    allowedNKSL = AllowedNKSLterms[n];
  allowedSL = AllowedSLterms[n];
  allowedSLJM = {}; 
  Mvals = Range[-J, J];
    For[nn = 1, nn <= Length[allowedSL],
      termSL = allowedSL[[ nn]];
      termNKSL = allowedNKSL[[ nn]];
      termsSLJ = 
    If[Abs[termSL[[1]] - termSL[[2]]] <= J <= 
      Total[termSL], {{termNKSL, J}}, {}];
      allowedSLJM = Join[allowedSLJM, termsSLJ];
      nn++];
    Tuples[{allowedSLJM, Mvals}]]

AllowedMforJ[j_] := Range[-j, j, 1];

AllowedNKSLJMIMforJIterms[n_, J_, \[CapitalIota]_] := 
  Tuples[{AllowedNKSLJMforJterms[n, J], Table[MI, {MI, -\[CapitalIota], \[CapitalIota]}]}];

AllowedJ[n_] := Table[J, {J, minJ[n], maxJ[n]}];

(* ############################### Allowed SL, SLJ, and SLJM Terms ############################## *)
(* ############################################################################################## *)

(* ############################################################################################## *)
(* ###############################              SOLVERS            ############################## *)

ShiftedLevels[originalLevels_] := Module[{GroundEnergy},
  GroundEnergy = First[First[Sort[originalLevels]]];
  Map[{#[[1]] - GroundEnergy, #[[2]], #[[3]]} &, originalLevels]]

LoadGuillotParameters::usage = "LoadGuillotParameters[Ln] loads into the current session the parameters for \!\(\*SuperscriptBox[\(Ln\), \(\(3\)\(+\)\)]\). The function returns the following list of quantities : {nf, E0, E1 , E2, E3 , \[Zeta], \[Alpha], \[Beta], \[Gamma], T2, T3, T4, T6, T7, T8, M0, M2, M4, P0, P2, P4, P6, B02, B06}";
LoadGuillotParameters::argerr = "Not a lanthanide."

LoadGuillotParameters[Ln_]:=Module[{},
  Which[
    Ln=="Ce" || Ln==1,
      (
      nf = 1; Fsup2 = 0; Fsup4 = 0; Fsup6 = 0;
      \[Zeta] = 647.3; \[Alpha] = 0; \[Beta] = 0; \[Gamma] = 0;
      T2 = 0; T3 = 0; T4 = 0; T6 = 0; T7 = 0; T8 = 0;
      M0 = 0; P2 = 0;
      eOrbitalRad = 1/4.8;
      Get[FileNameJoin[{moduleDir, "params", Ln, "CeLaF3_CFParameters.txt"}]];
      ),

    Ln=="Pr" || Ln==2,
      (
      nf = 2; Fsup2 = 68878; Fsup4 = 50347; Fsup6=32901;
      \[Zeta] = 751.7; \[Alpha] =16.23; \[Beta] = -566.6; \[Gamma] = 1371;
      T2 = 0; T3 = 0; T4 = 0; T6 = 0; T7 = 0; T8 = 0;
      M0 = 2.08; P2 = -88.6;
      eOrbitalRad = 1/5.5;
      ),

    Ln=="Nd" || Ln==3,
      (
      nf = 3; Fsup2 = 73018; Fsup4 = 52789; Fsup6 = 35757;
      \[Zeta] = 885.3; \[Alpha] = 21.34; \[Beta] = -593.0; \[Gamma] = 1445;
      T2 = 298; T3 = 35; T4 = 59; T6 = -285; T7 = 332; T8 = 305;
      M0 = 2.11; P2 = 192;
      eOrbitalRad = 1/6.2;
      ),

    Ln=="Pm" || Ln==4,
      (
      nf=4; Fsup2 = 76400;Fsup4 = 54900; Fsup6 = 37700;
      \[Zeta] = 1025; \[Alpha] = 20.50; \[Beta] = -560; \[Gamma] = 1475;
      T2 = 300; T3 = 35; T4 = 58; T6 = -310; T7 = 350; T8 = 320;
      M0 = 2.4; P2 = 275;
      eOrbitalRad = 1/6.9;
      ),

    Ln=="Sm" || Ln==5,
      (
      nf = 5; Fsup2 = 79805; Fsup4 = 57175; Fsup6 = 40250; 
      \[Zeta] = 1176; \[Alpha] = 20.16; \[Beta] = -566.9; \[Gamma] = 1500;
      T2 = 300; T3 = 36; T4 = 56; T6 = -347; T7 = 373; T8 = 348;
      M0 = 2.60; P2 = 357; 
      eOrbitalRad = 1/7.5;
      ),

    Ln=="Eu" || Ln==6,
      (
      nf = 6; Fsup2 = 83125; Fsup4 = 59268; Fsup6 = 42560;
      \[Zeta] = 1338; \[Alpha] = 20.16; \[Beta] = -566.9; \[Gamma] = 1500;
      T2 = 300; T3 = 40; T4 = 60; T6 = -300; T7 = 370; T8 = 320;
      M0 = 2.1; P2 = 360;
      eOrbitalRad = 1/8.4;
      ),

    Ln=="Gd" || Ln==7,
      (
      nf = 7; Fsup2 = 85669; Fsup4 = 60825; Fsup6 = 44776;
      \[Zeta] = 1508; \[Alpha] = 18.92; \[Beta] = -600; \[Gamma] = 1575;
      T2 = 300; T3 = 42; T4 = 62; T6 = -295; T7 = 350; T8 = 310;
      M0 = 3.22; P2 = 676;
      eOrbitalRad = 1/9.2;
      ),

    Ln=="Tb" || Ln==8,
      (
      nf = 8; Fsup2 = 88995; Fsup4 = 62919; Fsup6 = 47252;
      \[Zeta] = 1707; \[Alpha] = 18.40; \[Beta] = -590.9; \[Gamma] = 1650;
      T2 = 320; T3 = 40; T4 = 50; T6 = -395; T7 = 303; T8 = 317;
      M0 = 2.39; P2 = 373;
      eOrbitalRad = 1/10.0;
      ),

    Ln=="Dy" || Ln==9,
      (
      nf = 9; Fsup2 = 91903; Fsup4 = 64372; Fsup6 = 49386;
      \[Zeta] = 1913; \[Alpha] = 18.02; \[Beta] = -633.4; \[Gamma] = 1790;
      T2 = 329; T3 = 36; T4 = 127; T6 = -314; T7 = 404; T8 = 315;
      M0 = 3.39; P2 = 719;
      eOrbitalRad = 1/10.9;
      Get[FileNameJoin[{moduleDir, "params", "FreeIon_CFParameters.txt"}]];
      ),

    Ln=="Ho" || Ln==10,
      (
      nf = 10; Fsup2 = 94564; Fsup4 = 66397; Fsup6 = 52022;
      \[Zeta] = 2145; \[Alpha] = 17.15; \[Beta] = -607.9; \[Gamma] = 1800;
      T2 = 400; T3 = 37; T4 = 107; T6 = -264; T7 = 316; T8 = 336;
      M0 = 2.54; P2 = 605;
      eOrbitalRad = 1/11.8;
      Get[FileNameJoin[{moduleDir, "params", "FreeIon_CFParameters.txt"}]];
      ),

    Ln=="Er" || Ln==11,
      (
      nf = 11; Fsup2 = 97483; Fsup4 = 67904; Fsup6 = 54010;
      \[Zeta] = 2376; \[Alpha] = 17.79; \[Beta] = -582.1; \[Gamma] = 1800;
      T2 = 400; T3 = 43; T4 = 73; T6 = -271; T7 = 308; T8 = 299;
      M0 = 3.86; P2 = 594; 
      eOrbitalRad = 1/12.7;
      Get[FileNameJoin[{moduleDir, "params", Ln, "ErYVO4_Parameters.txt"}]];
      ),

    Ln=="Tm" || Ln==12,
      (
      nf = 12; E0 = 17443; E1 = 7018; E2 = 33.38; E3 = 671.5;
      \[Zeta] = 2619.7; \[Alpha] = 17.06; \[Beta] = -658; \[Gamma] = 0;
      T2 = 0; T3 = 0; T4 = 0; T6 = 0; T7 = 0; T8 = 0;
      M0 = 0; P2 = 0;
      eOrbitalRad = 1/13.6;
      Get[FileNameJoin[{moduleDir, "params", Ln, "TmYAG_GuillotNoelParameters.txt"}]];
      ),

    Ln=="Yb" || Ln==13,
      (
      nf = 13; Fsup2 = 0; Fsup4 = 0; Fsup6 = 0;
      luFnames = FileNames[All, FileNameJoin[{moduleDir, "params", "Yb"}]];
      shortNames = FileNameSplit[#][[-1]] & /@ luFnames;
      luFname = DialogInput[Prepend[MapThread[Button[#1,DialogReturn[#2]]&,{shortNames,luFnames}],TextCell[Style["These are the files available for Yb, choose one:",12]]]];
      Get[luFname];
      ),
    
    True,
    (Message[LoadGuillotParameters::argerr];
    Abort[];)
  ];

  Fsup0 = 0; T11 = 0; T12 = 0; T14 = 0; T15 = 0; T16 = 0; T17 = 0; T18 = 0; T19 = 0;
  M2 = 0.56 * M0; M4 = 0.31 * M0; 
  P0 = 0; P4 = 0.5 * P2; P6 = 0.1 * P2;
  gs = 2.002319304386;
  gI = 0.987;
  eV2Icm = 1/8065.5439;
  \[Beta]BohrMag = 5.7883818012 * 10^-5 / eV2Icm;
  \[Beta]n = 3.1524512550 * 10^-8 / eV2Icm;
  {nf, E0, E1 , E2, E3 , \[Zeta], \[Alpha], \[Beta], \[Gamma],
  T2, T3, T4, T6, T7, T8, M0, M2, M4, P0, P2, P4, P6, B02, B06}
]


(*nf is number of electrons, and IiN is nuclear spin val*)
SolveStates::usage = 
"SolveStates[nf, IiN] solves the energy values and states for an atom with n f-electrons with a nucleus of spin IiN. This function assumes that all the parameters that determine the Hamiltonian for a specific ion in a specific host have already been loaded into the session; these can be loaded by using the function LoadGuillotParameters.
This function also requires files for pre-computed energy matrix tables that provide the symbols EnergyMatrixTable[_, _, _, _, _].
To account for configurations f^n with n > 7, particle-hole dualities are enforced for \[Zeta] and T_i.
The unit for the returned energies is cm^-1.
-----------------------
References:
1. Sign inversion for \[Zeta]: Wybourne, Spectroscopic Properties of Rare Earths. 
2. Sign inversion for {T2, T3, T4, T6, T7, T8}: Hansen and Judd, Matrix Elements of Scalar Three Electron Operators for the Atomic f Shell.";

EnergyMatrixFileName[n_, IiN_] := (
  FileNameJoin[{moduleDir, "hams", 
    "f" <> ToString[n] <> "_I_" <> ToString[2*IiN + 1] <> "_Zhang_EnergyMatrixTables"
    }]
  )

SolveStates::notNum = "Necessary parameters have not been loaded into the session."

SolveStates[nf_, IiN_] := Module[
   {n, ii, jj, JMvals},
   (*#####################################*)
   (*hole-particle equivalence enforcement*)
   n = nf;
   allLoaded = And@@(NumberQ /@ {\[Zeta], T2, T3, T4, T6, T7, T8});
   If[Not[allLoaded],
    Message[SolveStates::notNum];
    Abort[];
    ];
   If[nf > 7,
    (n = 14 - nf;(*hole-particle equivalence*)
     {\[Zeta], T2, T3, T4, T6, T7, T8} = -{\[Zeta], T2, T3, T4, T6, T7, T8};
     )];
   (*hole-particle equivalence enforcement*)
   (*#####################################*)

   (* Load symbolic expressions for energy sub-matrices. *)
   Get[EnergyMatrixFileName[n, IiN]];
   (* Get[StringJoin[NotebookDirectory[], 
     StringJoin["f", ToString[n], "_I_", ToString[2*IiN + 1], 
      "_Zhang_EnergyMatrixTables"]]]; *)
   
   (*Patch together the entire matrix representation in block-diagonal form.*)
   EnergyMatrix = 
    ConstantArray[0, {Length[AllowedJ[n]], Length[AllowedJ[n]]}];
   Do[EnergyMatrix[[jj, ii]] = 
      EnergyMatrixTable[n, AllowedJ[n][[ii]], AllowedJ[n][[jj]], IiN, 
       IiN];,
    {ii, 1, Length[AllowedJ[n]]},
    {jj, 1, Length[AllowedJ[n]]}];
   EnergyMatrix = ArrayFlatten[EnergyMatrix];
   
   (* Print["The energy matrix has dimensions:", 
    Dimensions[EnergyMatrix]]; *)
   
   (*Solve for eigenvalues and eigenvectors.*)
   EnergyLevels = {};
   {EigenvalueJM, EigenvectorJM} = Eigensystem[EnergyMatrix];
   EigenvalueJM = Re[EigenvalueJM]; (*There might be a very small imaginary part.*)
   (* Print[{Dimensions@EigenvalueJM, Dimensions@EigenvectorJM}]; *)
   
   (*Parse the results for the eigenvectors in terms of the ordered basis being used.*)
   EnergyStates = {};
   Do[EnergyStates = Join[EnergyStates,
      EnergyStatesTable[n, AllowedJ[n][[nn]], IiN]]
    , {nn, 1, Length[AllowedJ[n]]}
    ];
   Do[EnergyLevels = Join[EnergyLevels,
       {{EigenvalueJM[[nn]], EigenvectorJM[[nn]], EnergyStates}}];
    , {nn, 1, Length[EigenvalueJM]}];
   If[nf > 7,
    ({\[Zeta], T2, T3, T4, T6, T7, T8} = -{\[Zeta], T2, T3, T4, T6, T7, T8};
     )];
   EnergyLevels];

(* ###############################              SOLVERS            ############################## *)
(* ############################################################################################## *)

(* ############################################################################################## *)
(* ###############################               MISC              ############################## *)


EnergyLevelDiagram::usage = "EnergyLevelDiagram[states] takes states and produces a visualization of its energy spectrum.
The resultant visualization can be navigated by clicking and dragging to zoom in on a region, or by clicking and dragging horizontally while pressing Ctrl. Double-click to reset the view."

Options[EnergyLevelDiagram]={"Title"->"", "ImageSize"->1000, "AspectRatio" -> 1/8};

EnergyLevelDiagram[states_,OptionsPattern[]]:=(
energies = First/@states;
ExploreGraphics@ListPlot[Tooltip[{{#,0},{#,1}}, {Quantity[#/8065.54429,"eV"], Quantity[#,1/"Centimeters"]}] &/@ energies,
  Joined -> True,
  PlotStyle -> Black,
  AspectRatio -> OptionValue["AspectRatio"],
  ImageSize -> OptionValue["ImageSize"],
  Frame -> True,
  PlotRange -> {All,{0,1}},
  FrameTicks -> {{None,None},{Automatic,Automatic}},
  FrameStyle -> Directive[15,Dashed,Thin],
  PlotLabel -> Style[OptionValue["Title"],15,Bold],
  FrameLabel -> {"\!\(\*FractionBox[\(E\), SuperscriptBox[\(cm\), \(-1\)]]\)"}]
)

ExploreGraphics::usage = 
  "Pass a Graphics object to explore it. Zoom by clicking and dragging a rectangle. Pan by clicking and dragging while pressing Ctrl. Click twice to reset view. 
  Based on ZeitPolizei @ https://mathematica.stackexchange.com/questions/7142/how-to-manipulate-2d-plots";

OptAxesRedraw::usage = 
  "Option for ExploreGraphics to specify redrawing of axes. Default False.";
Options[ExploreGraphics] = {OptAxesRedraw -> False};

ExploreGraphics[graph_Graphics, opts : OptionsPattern[]] := With[
   {gr = First[graph],
    opt = 
     DeleteCases[Options[graph], 
      PlotRange -> PlotRange | AspectRatio | AxesOrigin -> _],
    plr = PlotRange /. AbsoluteOptions[graph, PlotRange],
    ar = AspectRatio /. AbsoluteOptions[graph, AspectRatio],
    ao = AbsoluteOptions[AxesOrigin],
    rectangle = {Dashing[Small], 
       Line[{#1, {First[#2], Last[#1]}, #2, {First[#1], 
          Last[#2]}, #1}]} &,
    optAxesRedraw = OptionValue[OptAxesRedraw]},
   DynamicModule[{dragging = False, first, second, rx1, rx2, ry1, ry2,
      range = plr}, {{rx1, rx2}, {ry1, ry2}} = plr;
    Panel@
     EventHandler[
      Dynamic@Graphics[
        If[dragging, {gr, rectangle[first, second]}, gr],
        PlotRange -> Dynamic@range,
        AspectRatio -> ar,
        AxesOrigin -> If[optAxesRedraw,
          Dynamic@Mean[range\[Transpose]], ao],
        Sequence @@ opt],
      {{"MouseDown", 1} :> (
         first = MousePosition["Graphics"]
         ),
       {"MouseDragged", 1} :> (
         dragging = True;
         second = MousePosition["Graphics"]
         ),
       "MouseClicked" :> (
         If[CurrentValue@"MouseClickCount" == 2,
           range = plr];
         ),
       {"MouseUp", 1} :> If[dragging,
         dragging = False;
         
         range = {{rx1, rx2}, {ry1, ry2}} = 
           Transpose@{first, second};
         range[[2]] = {0, 1}],
       {"MouseDown", 2} :> (
         first = {sx1, sy1} = MousePosition["Graphics"]
         ),
       {"MouseDragged", 2} :> (
         second = {sx2, sy2} = MousePosition["Graphics"];
         rx1 = rx1 - (sx2 - sx1);
         rx2 = rx2 - (sx2 - sx1);
         ry1 = ry1 - (sy2 - sy1);
         ry2 = ry2 - (sy2 - sy1);
         range = {{rx1, rx2}, {ry1, ry2}};
         range[[2]] = {0, 1};
         )}]]];

(* ###############################               MISC              ############################## *)
(* ############################################################################################## *)

EndPackage[]
