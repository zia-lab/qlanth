(* -----------------------------------------------------------------
This  code was authored by Chris Dodson and turned into a package by
David Lizarazo in 2022-2023.
+------------------------------------------------------------------+
|                                                                  |
|                      __                     __     __            |
|            ____ _   / /  ____ _   ____     / /_   / /_           |
|           / __ `/  / /  / __ `/  / __ \   / __/  / __ \          |
|          / /_/ /  / /  / /_/ /  / / / /  / /_   / / / /          |
|          \__, /  /_/   \__, _/  /_/ /_/   \__/  /_/ /_/           |
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

On[Assert];

Off[SixJSymbol::tri];

(* Physical Constants*)
bohrRadius = 5.29177210903 * 10^-9;
ee = 1.602176634 * 10^-19;

specAlphabet = "SPDFGHIKLMNOQRTUV";
moduleDir = DirectoryName[$InputFileName];

SimplifyFun = Expand;

theLanthanides = {"Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb"};
theActinides = {"Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"};
paramAtlas = "
E0:
E1:
E2:
E3:

\[Zeta]: spin-orbit strength parameter.

F0:
F2:
F4:
F6:

M0: 0th Marvin integral
M2: 2nd Marvin integral
M4: 4th Marvin integral

T2:
T3:
T4:
T6:
T7:
T8:

T11:
T12:
T14:
T15:
T16:
T18:
T17:
T19:

P0: 0th parameter for the two-body electrostatically correlated spin-orbit interaction.
P2: 2th parameter for the two-body electrostatically correlated spin-orbit interaction.
P4: 4th parameter for the two-body electrostatically correlated spin-orbit interaction.
P6: 6th parameter for the two-body electrostatically correlated spin-orbit interaction.

gs: electronic gyromagnetic ratio
gI: nuclear spin gyromagnetic ratio

\[Alpha]: (Trees parameter).
\[Beta]:
\[Gamma]:

B02: crystal field parameter B_0^2
B04: crystal field parameter B_0^4
B06: crystal field parameter B_0^6
B12: crystal field parameter B_1^2
B14: crystal field parameter B_1^4
B16: crystal field parameter B_1^6
B22: crystal field parameter B_2^2
B24: crystal field parameter B_2^4
B26: crystal field parameter B_2^6
B34: crystal field parameter B_3^4
B36: crystal field parameter B_3^6
B44: crystal field parameter B_4^4
B46: crystal field parameter B_4^6
B56: crystal field parameter B_5^6
B66: crystal field parameter B_6^6

S12: crystal field parameter S_1^2
S14: crystal field parameter S_1^4
S16: crystal field parameter S_1^6
S22: crystal field parameter S_2^2
S24: crystal field parameter S_2^4
S26: crystal field parameter S_2^6
S34: crystal field parameter S_3^4
S36: crystal field parameter S_3^6
S44: crystal field parameter S_4^4
S46: crystal field parameter S_4^6
S56: crystal field parameter S_5^6
S66: crystal field parameter S_6^6

\[Epsilon]: ground level baseline shift
";
paramSymbols = StringSplit[paramAtlas, "\n"];
paramSymbols = Select[paramSymbols, # != ""& ];
paramSymbols = ToExpression[StringSplit[#, ":"][[1]]] & /@ paramSymbols;
Protect /@ paramSymbols;

Options[ParamPad] = {"Print" -> True}

ParamPad::usage = 
  "ParamPad[params] takes an Association params whose keys are a subset of paramSymbols. The function returns a new Association where all the keys not present in the given association, will also be present in the returned Association, with their values set to zero.
  The function additionally takes an option \"Print\" that if set to True, will print the symbols that were not present in the given Association, default is True.";
ParamPad[params_, OptionsPattern[]] := (
  notPresentSymbols = Complement[paramSymbols, Keys[params]];
  If[OptionValue["Print"], Print["Symbols not in given params: ", notPresentSymbols]];
  newParams = Transpose[{paramSymbols, ConstantArray[0, Length[paramSymbols]]}];
  newParams = (#[[1]] -> #[[2]]) & /@ newParams;
  newParams = Association[newParams];
  newParams = Join[newParams, params];
  Return[newParams];
  )

cfSymbols = {B02, B04, B06, B12, B14, B16, B22, B24, B26, B34, B36, 
   B44, B46, B56, B66, S12, S14, S16, S22, S24, S26, S34, S36, S44, 
   S46, S56, S66};

printL::usage = "printL[L] give the string representation of a given angular momentum.";
findSL::usage = "findSL[LS] gives the spin and orbital angular momentum that corresponds to the provided string LS.";
findSL\[Gamma]::usage = "findSL\[Gamma][LS] gives the spin and orbital angular momentum that corresponds to the provided string LS. The provided lists contains terms that belong to one of the f^n configurations. The last element of this list provides a digit that differentiates between different possibly degenerate terms.";
printSLJ::usage = "Given a list with three elements {S, L, J} this function returns a symbol where the spin multiplicity is presented as a superscript, the orbital angular momentum as its corresponding spectroscopic letter, and J as a subscript. Function does not check to see if the given J is compatible with the given S and L.";
printSLJM::usage = "Given a list with four elements {S, L, J, MJ} this function returns a symbol where the spin multiplicity is presented as a superscript, the orbital angular momentum as its corresponding spectroscopic letter, and {J, MJ} as a subscript. No attempt is made to guarantee that the given input is consistent.";

fnTerms::usage = "This list contains the labels of f^n configurations. Each element of the list has four elements {LS, seniority, RacahW, RacahU}.";
NKLSterms::usage = "NKLSterms[n, L, S] give all the terms that correspond to the given number of electrons n, with total orbital angular momentum L, and total spin angular momentum S.";
findNKLSterm::usage = "findNKLSterm[SL] for a given string SL that only represent the LS part of a spectroscopic term, this function returns all the terms that are compatible with it. This is only for f^n configurations. The provided terms might belong to more than one configuration.";
AA::usage = "Given a list with the parts corresponding to the quasi-rational representation of a number, this function parses the result into a regular number.";
CFPterms::usage = "CFPterms[n] gives all the daughter and parent terms, together with the corresponding coefficients of fractional parentage, that correspond to the the f^n configuration.

CFPterms[n, SL] gives all the daughter and parent terms, together with the corresponding coefficients of fractional parentage, that are compatible with the given string SL in the f^n configuration.

CFPterms[n, L, S] gives all the daughter and parent terms, together with the corresponding coefficients of fractional parentage, that correspond to the given total orbital angular momentum L and total spin S n the f^n configuration. L being an integer, and S being integer or half-integer.

In all cases the output is in the shape of a list with enclosed lists having the format {daughter_term, {parent_term_1, CFP_1}, {parent_term_2, CFP_2}, ...}.
Only the one-body coefficients for f-electrons are provided.
In all cases it must be that 1 <= n <= 7.

Data for these was parsed from Velkov, \"Multi-electron coefficients of fractional parentage for the p, d, and f shells\".";

AllowedSLterms::usage = "AllowedSLterms[n] returns a list with the allowed terms in the f^n configuration, the terms are given as lists in the format {S, L}. This list may have redundancies which are compatible with the degeneracies that might correspond to the given case.";
AllowedNKSLterms::usage = "AllowedNKSLterms[n] returns a list with the allowed terms in the f^n configuration, the terms are given as strings in spectroscopic notation. The integers in the last positions are used to distinguish cases with degeneracy.";
maxJ::usage = "maxJ[n] gives the maximum J = S+L that corresponds to the configuration f^n.";
minJ::usage = "minJ[n] gives the minimum J = S+L that corresponds to the configuration f^n.";
AllowedSLJterms::usage = "AllowedSLJterms[n] returns a list with the allowed {S, L, J} terms in the f^n configuration, the terms are given as lists in the format {S, L, J}. This list may have repeated elements which account for possible degeneracies of the related term.";

AllowedNKSLJterms::usage = "AllowedNKSLJterms[n] returns a list with the allowed {SL, J} terms in the f^n configuration, the terms are given as lists in the format {SL, J} where SL is a string in spectroscopic notation.";
AllowedNKSLforJterms::usage = "AllowedNKSLforJterms[n, J] gives the terms that correspond to the given total angular momentum J in the f^n configuration. The result is a list whose elements are lists of length 2, the first element being the SL term in spectroscopic notation, and the second element being J.";
AllowedMforJ::usage = "AllowedMforJ[j] is shorthand for Range[-j, j, 1].";
AllowedNKSLJMIMforJIterms::usage = "AllowedNKSLJMIMforJIterms[n, J, I] returns the states that belong to the f^n configuration with a total angular momentun L+S equal to J and with a nuclear angular momentum I. The returned list has elements of the form {{{SL (string in spectroscopic notation), J}, MJ}, MI}.";
AllowedJ::usage = "AllowedJ[n] returns the total angular momenta J that appear in the f^n configuration.";

AllowedSLJMterms::usage = "AllowedSLJMterms[n] returns a list with all the states that correspond to the configuration f^n. A list is returned whose elements are lists of the form {S, L, J, MJ}.";
AllowedNKSLJMforJMterms::usage = "AllowedNKSLJMforJMterms[n, J, MJ] returns a list with all the terms that contain states of the f^n configuration that have a total angular momentum J, and a projection along the z-axis MJ. The returned list has elements of the form {SL (string in spectroscopic notation), J, MJ}.";
AllowedNKSLJMforJterms::usage = "AllowedNKSLJMforJterms[n, J] returns a list with all the states that have a total angular momentum J. The returned list has elements of the form {{SL (string in spectroscopic notation), J}, MJ}.";

ShiftedLevels::usage = "
ShiftedLevels[originalLevels] takes a list of levels of the form
{{energy_1, coeff_vector_1}, 
{energy_2, coeff_vector_2},
...}} 
and returns the same input except that now to every energy the minimum of all of them has been subtracted.";
LoadGuillotParameters::usage = "";

CFP::usage = "CFP[{n, NKSL}] provides a list whose first element echoes NKSL and whose other elements are lists with two elements the first one being the symbol of a parent term and the second being the corresponding coefficient of fractional parentage. n must satisfy 1 <= n <= 7";

chenDeltasRaw = 
  Import[FileNameJoin[{moduleDir, "data", "the-chen-deltas.csv"}]];
chenDeltasRaw = chenDeltasRaw[[2 ;;]];
chenDeltas = <||>;
Off[Power::infy];
Do[
  ({right, wrong} = {chenDeltasRaw[[row]][[4 ;;]], 
     chenDeltasRaw[[row + 1]][[4 ;;]]};
   key = chenDeltasRaw[[row]][[1 ;; 3]];
   repRule = (#[[1]] -> #[[2]]*#[[1]]) & /@ 
     Transpose[{{M0, M2, M4, P2, P4, P6}, right/wrong}];
   chenDeltas[key] = <|"right" -> right, "wrong" -> wrong, 
     "repRule" -> repRule|>;
   chenDeltas[{key[[1]], key[[3]], key[[2]]}] = <|"right" -> right, 
     "wrong" -> wrong, "repRule" -> repRule|>;
   ),
  {row, 1, Length[chenDeltasRaw], 2}];
On[Power::infy]

(* ############################################################################################## *)
(* ################################### Racah Algebra Goodies #################################### *)

ReducedUk::usage = "ReducedUk[n, l, SL, SpLp, k] gives the reduced matrix element of the spherical tensor U^(k).";
ReducedUk[n_, l_, SL_, SpLp_, k_]:= Module[
  {Uk, S, L, Sp, Lp, Sb, Lb, parentSL, cfpSL, cfpSpLp, Ukval},
  {S, L}   = findSL[SL];
  {Sp, Lp} = findSL[SpLp];
  cfpSL   = CFP[{n, SL}];
  cfpSpLp = CFP[{n, SpLp}];
  Uk      = 0;
  If[S==Sp,
    For[nn = 2,
        nn <= Length[cfpSL],
        nn ++,
        (
        parentSL = cfpSL[[nn, 1]];
        {Sb, Lb}  = findSL[parentSL];
        fun = If[
                #[[1]]!= parentSL,
                0,
                ((-1)^(Lb + L + l + k)
                * #[[2]]
                * cfpSL[[nn, 2]]
                * SixJay[{L, Lp, k}, {l, l, Lb}])
                ]&;
        summands = Map [fun, cfpSpLp[[2;;]]];
        Uk      += Total[summands];
        )
    ];
    Ukval = n * Sqrt[(2L+1) * (2Lp+1)] * Uk;
  ,
  Ukval = 0;
  ];
  Return[Ukval];
]

Ck::usage = "Diagonal reduced matrix element <l||C^(k)\[VerticalSeparator]\[VerticalSeparator]l> where the Subscript[C, q]^(k) are reduced spherical harmonics.";
Ck[l_, k_] := (-1)^l *(2 l+1) * ThreeJay[{l, 0}, {l, 0}, {k, 0}]

fk::usage = "Slater integral. See Wybourne (1965) eqn. 2-93.";
fk[n_, l_, NKSL_, NKSLp_, k_]:= Module[
  {S, L, PossibleStates, Sb, Lb, SL},
  {S, L}         = findSL[NKSL];
  PossibleStates = StringCases[AllowedNKSLterms[n], ToString[2 S+1]~~__];
  PossibleStates = DeleteCases[PossibleStates, {}];
  PossibleStates = Flatten[PossibleStates];
  PossibleStates = If[{S, L}==findSL[NKSLp],
      (1/2
      * ( Abs[Ck[l, k]])^2 
      * ( 1/(2L+1) 
        * Sum[( ReducedUkTable[{n, l, SL, NKSL , k}]
              * ReducedUkTable[{n, l, SL, NKSLp, k}]),
              {SL, PossibleStates}]
        - (KroneckerDelta[NKSL, NKSLp] 
          * (n (4 l+2-n)) 
          / (2 l+1)
          / (4 l+1))
        )
      ),
      0
    ]
]

fK::usage = "Non-reduced Slater integral f^k = Subscript[f, k] * Subscript[D, k]";
fK[n_, l_, NKSL_, NKSLp_ ,k_]:= (Dk[k] * fk[n, l, NKSL, NKSLp, k])

Dk::usage = "Ratio between the reduced and non-reduced Slater direct (Subscript[F, k] and F^k) and exchange(Subscript[G, k]and G^k) integrals. Subscript[D, k]:= (Subscript[G, k](ff))/(G^k (ff)) = (Subscript[F, k](ff))/(F^k (ff)). k must be even. See table 6-3 in Cowan (1981), and also section 2-7 of Wybourne (1965).";
Dk[k_] := {1, 225, 1089, 184041/25}[[k/2+1]]

FtoE::usage = "FtoE[F0, F2, F4, F6] calculates the corresponding {E0, E1, E2, E3} values.
See eqn. 2-80 in Wybourne. Note that in that equation the subscripted Slater integrals are used but since this function assumes the the input values are superscripted Slater integrals, it is necessary to convert them using Dk.";
FtoE[F0_, F2_, F4_, F6_] := (Module[
  {E0, E1, E2, E3},
  E0 = (F0 - 10/Dk[2]*F2 - 33/Dk[4]*F4 - 286/Dk[6]*F6);
  E1 = (70/Dk[2]*F2 + 231/Dk[4]*F4 + 2002/Dk[6]*F6)/9;
  E2 = (1/Dk[2]*F2 - 3/Dk[4]*F4 + 7/Dk[6]*F6)/9;
  E3 = (5/Dk[2]*F2 + 6/Dk[4]*F4 - 91/Dk[6]*F6)/3;
  Return[{E0, E1, E2, E3}];
]
);

EtoF::usage = "EtoF[E0, E1, E2, E3] calculates the corresponding {0, F2, F4, F6} values.";
EtoF[E0_, E1_, E2_, E3_] := (Module[
  {F0, F2, F4, F6},
  F0 = 1/7 (7 E0 + 9 E1);
  F2 = 75/14 (E1 + 143 E2 + 11 E3);
  F4 = 99/7 (E1 - 130 E2 + 4 E3);
  F6 = 5577/350 (E1 + 35 E2 - 7 E3);
  Return[{F0, F2, F4, F6}];
]
);

ElectrostaticMatrix::usage = "See Wybourne (1965) section 2-9 'Electrostatic Interactions in f^n Configurations'";
ElectrostaticMatrix[n_, NKSL_, NKSLp_]:= Module[
  {f0, f2, f4, f6, e0, e1, e2, e3, EMatrix, l},
  l = 3;
  Ek = {E0, E1, E2, E3};
  f0 = fK[n, l, NKSL, NKSLp, 0];
  f2 = fK[n, l, NKSL, NKSLp, 2];
  f4 = fK[n, l, NKSL, NKSLp, 4];
  f6 = fK[n, l, NKSL, NKSLp, 6];
  e0 = f0;
  e1 = 9/7*f0 + f2/42 + f4/77 + f6/462;
  e2 = 143/42*f2 - 130/77*f4 + 35/462*f6;
  e3 = 11/42*f2 + 4/77*f4 - 7/462*f6;
  EMatrix = e0*E0 + e1*E1 + e2*E2 + e3*E3;
  Return[EMatrix];
]

(* ################################### Racah Algebra Goodies #################################### *)
(* ############################################################################################## *)

(* ############################################################################################## *)
(* ####################################### Table Generation Functions ########################### *)

GenerateCFP::usage = "Generates the tables for the coefficients of fractional parentage. Results are exported to the file CFP.m";
GenerateCFP[export_:True]:= (
  CFP = Table[
    {n, NKSL} -> First[CFPterms[n, NKSL]],
    {n, 1, 7},
    {NKSL, AllowedNKSLterms[n]}];
  CFP = Association[CFP];
  If[export,
    Export[FileNameJoin[{moduleDir, "data", "CFPs.m"}], CFP];
  ];
  Return[CFP];
)

GenerateCFPAssoc::usage = "GenerateCFPAssoc[export] converts the coefficients of fractional parentage into an association in which zero values are explicit. If export is True, the association is exported to the file CFPAssoc.m.";
GenerateCFPAssoc[export_:True]:= (
  CFPAssoc = Association[];
  Do[
    (daughterTerms = AllowedNKSLterms[\[CapitalNu]];
    parentTerms    = AllowedNKSLterms[\[CapitalNu] - 1];
    Do[
      (
      cfps = CFP[{\[CapitalNu], daughter}];
      cfps = cfps[[2 ;;]];
      parents = First /@ cfps;
      Do[
        (
        key = {\[CapitalNu], daughter, parent};
        cfp = If[
          MemberQ[parents, parent],
          (
            idx = Position[parents, parent][[1, 1]];
            cfps[[idx]][[2]]
            ),
          (
            0
            )
          ];
        CFPAssoc[key] = cfp;
        ),
        {parent, parentTerms}
        ]
      ),
      {daughter, daughterTerms}
      ]
    ),
    {\[CapitalNu], 2, 7}
    ];
  If[export,
    (
    CFPAssocfname = FileNameJoin[{moduleDir, "data", "CFPAssoc.m"}];
    Export[CFPAssocfname, CFPAssoc];
    )
  ];
  Return[CFPAssoc];
)

GenerateElectrostaticMatrixTable[export_:True]:= (
  ElectrostaticMatrixTable = Table[
    {n, SL, SpLp} -> SimplifyFun[ElectrostaticMatrix[n, SL, SpLp]], 
    {n, 1, 7},
    {SL, AllowedNKSLterms[n]}, 
    {SpLp, AllowedNKSLterms[n]}
  ];
  ElectrostaticMatrixTable = Association[Flatten[ElectrostaticMatrixTable]];
  If[export,
    Export[FileNameJoin[{moduleDir, "data", "ElectrostaticMatrixTable.m"}], ElectrostaticMatrixTable];
  ];
  Return[ElectrostaticMatrixTable];
)
GenerateElectrostaticMatrixTable::usage = "GenerateElectrostaticMatrixTable[] can be used to generate the table for the electrostatic interaction. Results are exported to the file ElectrostaticMatrixTable.m.";

GenerateReducedUkTable::usage = "GenerateReducedUkTable[] can be used to generate the table of reduced matrix elements for the unit tensor operators Uk. Results are exported to the file ReducedUkTable.m.";
GenerateReducedUkTable[export_:True]:= (
  Off[SixJSymbol::tri];
  ReducedUkTable = Table[
    {n, 3, SL, SpLp, k} -> SimplifyFun[ReducedUk[n, 3, SL, SpLp, k]],
    {n, 1, 7},
    {SL, AllowedNKSLterms[n]}, 
    {SpLp, AllowedNKSLterms[n]}, 
    {k, {0, 2, 4, 6}}
  ];
  ReducedUkTable = Association[Flatten[ReducedUkTable]];
  ReducedUkTableFname = FileNameJoin[{moduleDir, "data", "ReducedUkTable.m"}];
  If[export,
    Export[ReducedUkTableFname, ReducedUkTable];
  ];
  On[SixJSymbol::tri];
  Return[ReducedUkTable];
)

SixJay::usage = "SixJay[{j1, j2, j3}, {j4, j5, j6}] provides the value for SixJSymbol[{j1, j2, j3}, {j4, j5, j6}] with the addition that it remembers values it has already computed.";
SixJay[{j1_, j2_, j3_}, {j4_, j5_, j6_}]:= (
  SixJay[{j1, j2, j3}, {j4, j5, j6}] = SixJSymbol[{j1, j2, j3}, {j4, j5, j6}]
);

ThreeJay::usage = "ThreeJay[{j1, m1}, {j2, m2}, {j3, m3}] gives the value of the Wigner 3j-symbol and memorizes previously computed values.";
ThreeJay[{j1_, m1_}, {j2_, m2_}, {j3_, m3_}] := (
   ThreeJay[{j1, m1}, {j2, m2}, {j3, m3}] = ThreeJSymbol[{j1, m1}, {j2, m2}, {j3, m3}]
);

ReducedVk1::usage = "ReducedVk1";
ReducedVk1[n_, SL_, SpLp_, k_]:= Module[
  {Vk1, S, L, Sp, Lp, Sb, Lb, parentSL, s, l, cfpSL, cfpSpLp, prefact, sign},
  {s, l}   = {1/2, 3};
  {S, L}   = findSL[SL];
  {Sp, Lp} = findSL[SpLp];
  cfpSL   = CFP[{n, SL}];
  cfpSpLp = CFP[{n, SpLp}];
  Vk1     = 0;
  For[nn = 2,
    nn<= Length[cfpSL],
    nn++,
    parentSL = cfpSL[[n, 1]];
    {Sb, Lb}  = findSL[parentSL];
    Vk1fun   = (If[Part[#, 1] != parentSL,
                0,
                (#[[2]] * cfpSL[[nn, 2]]
                * SixJay[{S, Sp, 1}, {s, s, Sb}]
                * SixJay[{L, Lp, k}, {l, l, Lb}]
                * (-1)^(Sb + Lb + S + L + l + k - s)
                )]&
               );
    Vk1adds  = Map[Vk1fun, cfpSpLp[[2;;]]];
    Vk1     += Total[Vk1adds];
  ];
  prefact = n Sqrt[s(s+1)(2s+1)(2S+1)(2L+1)(2Sp+1)(2Lp+1)];
  Return[prefact*Vk1];
]

GenerateVk1Table::usage = "GenerateVk1Table[nmax, export calculates values for Vk1 and returns an Association where the keys are lists of the form {n, SL, SpLp, 1}. If export is set to True, then the result is exported to the data subfolder for the folder in which this package is in."
GenerateVk1Table[nmax_:7, export_:False]:= (
  Off[SixJSymbol::tri];
  numValues = Total[Length[AllowedNKSLterms[#]]*Length[AllowedNKSLterms[#]]&/@Range[1, nmax]];
  Print["Calculating " <> ToString[numValues] <> " values for Vk1."];

  counter = 1;
  Monitor[
    ReducedVk1Table = Table[
      (counter = counter+1;
      {n, SL, SpLp, 1} -> SimplifyFun[ReducedVk1[n, SL, SpLp, 1]]
      ),
      {n, 1, nmax},
      {SL, AllowedNKSLterms[n]}, 
      {SpLp, AllowedNKSLterms[n]}
    ],
    Row[{ProgressIndicator[counter, {0, numValues}], " ", counter}]
  ];

  ReducedVk1Table = Association[ReducedVk1Table];

  exportFname = FileNameJoin[{moduleDir, "data", "ReducedVk1Table.m"}];
  If[export,
  (Print["Exporting to file "<>ToString[exportFname]];
  Export[exportFname, ReducedVk1Table];
  )
  ];
  Return[ReducedVk1Table];
)

(* ######################################################################### *)
(* ############################# Spin Orbit ################################ *)

SpinOrbit::usage = "SpinOrbit[n, SL, SpLp, J] returns the matrix element ζ <SL, J|L⋅S|SpLp, J> within an f^n configuration. These are given as a function of ζ. This function requires that the Association ReducedVk1Table be defined.";
SpinOrbit[n_, SL_, SpLp_, J_]:= Module[
  {S, L, Sp, Lp, l, sign, prefact},
  l       = 3;
  {S, L}   = findSL[SL];
  {Sp, Lp} = findSL[SpLp];
  prefact = Sqrt[l*(l+1)*(2l+1)] * SixJay[{L, Lp, 1}, {Sp, S, J}];
  sign    = (-1)^(J + L + Sp);
  Return[sign * prefact * \[Zeta] * ReducedVk1Table[{n, SL, SpLp, 1}]];
]

GenerateSpinOrbitTable::usage = "GenerateSpinOrbitTable[nmax, export] computes the matrix values for the spin-orbit interaction for f^n configurations up to n = nmax. The function returns an Association whose keys are lists of the form {n, SL, SpLp, J}. If export is set to True, then the result is exported to the data subfolder for the folder in which this package is in. It requires ReducedVk1Table to be defined.";
GenerateSpinOrbitTable[nmax_:7, export_:False]:= (
  SpinOrbitTable = Table[
    {n, SL, SpLp, J} -> SpinOrbit[n, SL, SpLp, J], 
    {n, 1, nmax},
    {J, minJ[n], maxJ[n]},
    {SL, Map[First, AllowedNKSLforJterms[n, J]]},
    {SpLp, Map[First, AllowedNKSLforJterms[n, J]]}
  ];
  SpinOrbitTable = Association[SpinOrbitTable];

  exportFname = FileNameJoin[{moduleDir, "data", "SpinOrbitTable.m"}];
  If[export,
  (Print["Exporting to file "<>ToString[exportFname]];
  Export[exportFname, SpinOrbitTable];
  )
  ];
  Return[SpinOrbitTable];
)

GenerateThreeBodyTable::usage = "GenerateThreeBodyTable[nmax, export] computes the matrix values for the three-body interaction for f^n configurations up to n = nmax. The function returns an Association whose keys are lists of the form {n, SL, SpLp}. If export is set to True, then the result is exported to the data subfolder for the folder in which this package is in.";
GenerateThreeBodyTable[nmax_ : 7, export_ : False] := (
  ti3fname = FileNameJoin[{moduleDir, "data", "ti3.m"}];
  ti3 = Import[ti3fname];
  
  ti3PsiPsipStates[SL_, SpLp_, tnp_, k_] := 
   Module[{Statesn1, tkvals, PositionofState, kposition, ii, Snt, 
     Sntt, Lnt, Lntt, jj, S, L, Sp, Lp, tnk},
    Statesn1 = 
     Reverse[Sort[
       Join[Thread[{AllowedNKSLterms[3], AllowedNKSLterms[3]}], 
        Permutations[AllowedNKSLterms[3], {2}]]]];
    For[ii = 1, ii <= Length[Statesn1], ii++, {
      {Snt, Lnt} = findSL[Statesn1[[ii, 1]]];
      {Sntt, Lntt} = findSL[Statesn1[[ii, 2]]];
      If[Snt == Sntt && Lnt == Lntt, 
       0, {Statesn1 = Delete[Statesn1, ii], ii = ii - 1}]}];
    tkvals = {2, 3, 4, 6, 7, 8, 11, 12, 14, 15, 16, 17, 18, 19};
    kposition = Position[tkvals, k][[1, 1]];
    {S, L} = findSL[SL];
    {Sp, Lp} = findSL[SpLp];
    tnk = 0;
    If[S == Sp && L == Lp,
     PositionofState = 0;
     jj = 1;
     While[PositionofState == 0 && jj <= Length[Statesn1], 
      If[Statesn1[[jj, 1]] == SL && Statesn1[[jj, 2]] == SpLp, 
       PositionofState = jj, PositionofState = 0]; jj++];
     If[PositionofState == 0, tnk = tnk + 0, 
      tnk = tnp[[kposition, PositionofState]]], 0];
     tnk
    ];
  
  ti[n_, SL_, SpLp_, k_] := 
   Module[{nn, S, L, Sp, Lp, cfpSL, cfpSpLp, parentSL, parentSpLp, 
     tnk},
    {S, L} = findSL[SL];
    {Sp, Lp} = findSL[SpLp];
    tnk = 0;
    If[S == Sp && L == Lp,
     (cfpSL = CFP[{n, SL}];
      cfpSpLp = CFP[{n, SpLp}];
      For[nn = 2, nn <= Length[cfpSL], nn++,
       parentSL = Part[cfpSL, nn, 1];
       tnk = (tnk
          + Sum[
           (cfpSL[[nn, 2]] 
             *cfpSpLp[[mm, 2]] 
             *tktable[{n - 1, parentSL, Part[cfpSpLp, mm, 1], k}]
            ),
           {mm, 2, Length[cfpSpLp]}]);
       ]
      )
     ];
    Return[n/(n - 3) tnk ];
    ];

  (*Calculate the matrix elements of t^i*)

  tktable = Association[];
  kvalues = {2, 3, 4, 6, 7, 8, 11, 12, 14, 15, 16, 17, 18, 19};
  (*Initialize n=1,2*)
  Do[
   (
    Do[
     (
      tktable[{n, SL, SpLp, k}] = Which[
         n <= 2,
         0,
         n == 3,
         SimplifyFun[ti3PsiPsipStates[SL, SpLp, ti3, k]],
         True,
          SimplifyFun[ti[n,  SL, SpLp, k]]
         ];
      ),
     {SL,     AllowedNKSLterms[n]}, 
     {SpLp, AllowedNKSLterms[n]}, 
     {k, kvalues}
     ];
    PrintTemporary[
     StringJoin["\[ScriptF]", ToString[n], " configuration complete"]];
    ),
    {n, 1, nmax}
   ];

  (*Now use those matrix elements to determine their sum as weighted by their corresponding strenghts Ti*)
  Tparams = {T2, T3, T4, T6, T7, T8, T11, T12, T14, T15, T16, T17, 
    T18, T19};
  ThreeBodyTable = Association[];
  Do[
   Do[
    ThreeBodyTable[{n, SL, SpLp}] = Sum[(
        tktable[{n, SL, SpLp, kvalues[[kk]]}]
         *Tparams[[kk]]),
       {kk, 1, 14}];, 
    {SL, AllowedNKSLterms[n]}, 
    {SpLp, AllowedNKSLterms[n]}
    ];
   PrintTemporary[
    StringJoin["\[ScriptF]", ToString[n], " matrix complete"]];,
    {n, 1, nmax}];
  If[export,
   (
    threeBodyTablefname = 
     FileNameJoin[{moduleDir, "data", "ThreeBodyTable.m"}];
    Export[threeBodyTablefname, ThreeBodyTable];
    )
   ];
  Return[ThreeBodyTable];
  )

(* ############################# Spin Orbit ################################ *)
(* ######################################################################### *)

(* ######################################################################### *)
(* ################ Magnetically-Correlated Corrections #################### *)


T222PsiPsipStates::usage="T222PsiPsipStates[SL, SpLp] returns the reduced matrix element of the scalar component of the double tensor T22 for the given SL terms SL, SpLp.
Data used here for m0, m2, m4 is from Table I of Judd, BR, HM Crosswhite, and Hannah Crosswhite. Intra-Atomic Magnetic Interactions for f Electrons. Physical Review 169, no. 1 (1968): 130.
";
T222PsiPsipStates[SL_, SpLp_] := 
  Module[{statePosition, PsiPsipStates, m0, m2, m4, Tkk2m},
   PsiPsipStates = {
    {"3P", "3P"}, {"3P", "3F"}, 
    {"3F", "3P"}, {"3F", "3F"}, 
    {"3F", "3H"}, {"3H", "3F"},
    {"3H", "3H"}
    };
   m0 = {-12, 
      3 8/Sqrt[3], 3 8/Sqrt[3], 
      -(4/3) Sqrt[14], 
      2 8/3 Sqrt[11/2], 2 8/3 Sqrt[11/2], 
      4/3 Sqrt[143]
      };
   m2 = {-24, 
     8/Sqrt[3], 8/Sqrt[3], 
     8 4/3 Sqrt[14], 
     -(23/11) 8/3 Sqrt[11/2], -(23/11) 8/3 Sqrt[11/2],
     -(34/11) 4/3 Sqrt[143]};
   m4 = {-(300/11), 
      -(100/11) 8/Sqrt[3], -(100/11) 8/Sqrt[3], 
      -(200/11) 4/3 Sqrt[14], 
      -(325/121) 8/3 Sqrt[11/2], 
      -(325/121) 8/3 Sqrt[11/2], 
      -(1325/1573) 4/3 Sqrt[143]
      };
   statePosition = Position[PsiPsipStates, {SL, SpLp}];
   If[Length[statePosition] == 0,
    Return[0],
    statePosition = statePosition[[1, 1]]];
   Tkk2m = (m0[[statePosition]] * M0 + m2[[statePosition]] * M2 + 
      m4[[statePosition]] * M4);
   Return[Tkk2m];
   ];

T112PsiPsipStates::usage="T112PsiPsipStates[SL, SpLp] returns the reduced matrix element of the scalar component of the double tensor T11 for the given SL terms SL, SpLp.
Data used here for m0, m2, m4 is from Table II of Judd, BR, HM Crosswhite, and Hannah Crosswhite. Intra-Atomic Magnetic Interactions for f Electrons. Physical Review 169, no. 1 (1968): 130.
";
T112PsiPsipStates[SL_, SpLp_] := 
  Module[{statePosition, PsiPsipStates, m0, m2, m4, Tkk2m},
   PsiPsipStates = {
      {"1S", "3P"}, {"3P", "1S"},
      {"3P", "3P"}, {"3P", "1D"},
      {"1D", "3P"}, {"1D", "3F"},
      {"3F", "1D"}, {"3F", "3F"}, 
      {"3F", "1G"}, {"1G", "3F"},
      {"1G", "3H"}, {"3H", "1G"},
      {"3H", "3H"}, {"3H", "1I"},
      {"1I", "3H"}
      };
   m0 = {6, 6, 
        -36,
        -Sqrt[(2/15)] 23, -Sqrt[(2/15)] 23,
        Sqrt[2/5] 23, Sqrt[2/5] 23, 
        -2 Sqrt[14] 15, 
        -Sqrt[11] 6, -Sqrt[11] 6, 
        Sqrt[2/5] 39, Sqrt[2/5] 39, 
        -8 1/Sqrt[55] 132, 
        -5 Sqrt[26], -5 Sqrt[26]
        };
   m2 = {2, 2, 
        -72, 
        -Sqrt[(2/15)] 14, 
        -Sqrt[(2/15)] 14, 
        Sqrt[2/5] 6, Sqrt[2/5] 6,
        -2 Sqrt[14], 
        Sqrt[11] 64/33, Sqrt[11] 64/33, 
        -Sqrt[(2/5)] 728/33, -Sqrt[(2/5)] 728/33, 
        8/Sqrt[55] 23, 
        -(30/11) Sqrt[26], -(30/11) Sqrt[26]
        };
   m4 = {10/11, 10/11, 
      -(900/11), 
      -Sqrt[(2/15)] 115/11, -Sqrt[(2/15)] 115/11, 
      -(195/11) Sqrt[2/5], -(195/11) Sqrt[2/5], 
      2 Sqrt[14] 10/11, 
      -Sqrt[11] 1240/363, -Sqrt[11] 1240/363, 
      -(3175/363) Sqrt[2/5], -(3175/363) Sqrt[2/5], 
      8/Sqrt[55] 130/11, 
      -(375/1573) Sqrt[26], -(375/1573) Sqrt[26]
      };
   statePosition = Position[PsiPsipStates, {SL, SpLp}];
   If[Length[statePosition] == 0,
    Return[0],
    statePosition = statePosition[[1, 1]]];
   Tkk2m = (m0[[statePosition]] * M0 + m2[[statePosition]] * M2 + 
      m4[[statePosition]] * M4);
   Return[Tkk2m];
   ];

(* ################ Magnetically-Correlated Corrections #################### *)
(* ######################################################################### *)

(* ######################################################################### *)
(* ############################ Reduced t^11 ############################### *)

t112PsiPsipStates::usage="t112PsiPsipStates[SL, SpLp] returns the reduced matrix element of the double tensor operator t11 for the corresponding given terms {SL, SpLp}.

The values for p0, p2, and p4 are the coefficients obtained from Table VII of \"Judd, BR, HM Crosswhite, and Hannah Crosswhite. “Intra-Atomic Magnetic Interactions for f Electrons.” Physical Review 169, no. 1 (1968): 130.\"

"
t112PsiPsipStates[SL_, SpLp_]:= Module[
  {statePosition, PsiPsipStates, p0, p2, p4, p6, t112p},

  PsiPsipStates = {{"1S", "3P"}, {"3P", "1S"}, {"3P", "3P"},
                   {"3P", "1D"}, {"1D", "3P"}, {"1D", "3F"},
                   {"3F", "1D"}, {"3F", "3F"}, {"3F", "1G"},
                   {"1G", "3F"}, {"1G", "3H"}, {"3H", "1G"},
                   {"3H", "3H"}, {"3H", "1I"}, {"1I", "3H"}};
  p0 = {-2, -2,
        -1,
        Sqrt[15/2], Sqrt[15/2],
        -Sqrt[10], -Sqrt[10],
        -Sqrt[14],
        Sqrt[11], Sqrt[11],
        -Sqrt[10], -Sqrt[10],
        -Sqrt[55], 
        Sqrt[13/2], Sqrt[13/2]};
  p2 = {-105, -105,
      -45, 
      Sqrt[15/2] 32, Sqrt[15/2] 32,
      -Sqrt[10] 9/2, -Sqrt[10] 9/2,
      Sqrt[14] 10,
      -Sqrt[11] 20, -Sqrt[11] 20,
      Sqrt[10] 55/2, Sqrt[10] 55/2,
      Sqrt[55] 25, 
      0, 0};
  p4 = {-231, -231,
      -33, 
      -Sqrt[(15/2)] 33, -Sqrt[(15/2)] 33, 
      Sqrt[10] 66, Sqrt[10] 66,
      Sqrt[14] 33, 
      Sqrt[11] 32, Sqrt[11] 32,
      -Sqrt[10] 23, -Sqrt[10] 23, 
      Sqrt[55] 51, 
      -Sqrt[(13/2)] 21, -Sqrt[(13/2)] 21};
  p6 = {-429, -429,
      1287,
      Sqrt[15/2] (-286), Sqrt[15/2] (-286),
      Sqrt[10] (-(429/2)), Sqrt[10] (-(429/2)),
      Sqrt[14] 286,
      Sqrt[11] (-104), Sqrt[11] (-104),
      Sqrt[10] (-(65/2)), Sqrt[10] (-(65/2)),
      Sqrt[55] 13, 
      Sqrt[13/2] (-6), Sqrt[13/2] (-6)};

  statePosition = Position[PsiPsipStates, {SL, SpLp}];
  t112p = If[Length[statePosition]==0,
    0,
    (
      statePosition = statePosition[[1, 1]];
      ( p0[[statePosition]] P0
      + p2[[statePosition]] P2
      + p4[[statePosition]] P4
      + p6[[statePosition]] P6)
    )];
  Return[t112p];
]

(* ############################ Reduced t^11 ############################### *)
(* ######################################################################### *)

(* ######################################################################### *)
(* ############################ Reduced AiZi ############################### *)

AiZi2PsiPsipStates::usage="AiZi2PsiPsipStates[SL, SpLp] returns the sum of the reduced matrix elements corresponding to the operator T11, t11, and -a14 z13 for the corresponding given terms {SL, SpLp}.

The value for a13 is taken from Table VIII of \"Judd, BR, HM Crosswhite, and Hannah Crosswhite. “Intra-Atomic Magnetic Interactions for f Electrons.” Physical Review 169, no. 1 (1968): 130.\"

The values for the reduced matrix elements of z13 are obtained from Table IX of the same paper.
";
AiZi2PsiPsipStates[SL_, SpLp_] :=
Module[{statePosition, PsiPsipStates, a13, z13, aizi2p}, 
  PositionofState = 0;
  kk = 1;
  PsiPsipStates = {
    {"1S", "3P"}, {"3P", "1S"}, 
    {"3P", "3P"}, {"3P", "1D"}, 
    {"1D", "3P"}, {"1D", "3F"}, 
    {"3F", "1D"}, {"3F", "3F"}, 
    {"3F", "1G"}, {"1G", "3F"}, 
    {"1G", "3H"}, {"3H", "1G"}, 
    {"3H", "3H"}, {"3H", "1I"}, 
    {"1I", "3H"}};
  a13 = (-33 M0 + 3 M2 + 15/11 M4 
        - 6 P0 + 3/2 (35 P2 + 77 P4 + 143 P6));
  z13 = {2, 2, 
    1, 
    1/Sqrt[1080] (-90), 
    1/Sqrt[1080] (-90), 
    Sqrt[2/405] 45,
    Sqrt[2/405] 45,
    Sqrt[14],
    1/Sqrt[891] (-99),
    1/Sqrt[891] (-99), 
    990/Sqrt[98010],
    990/Sqrt[98010], 
    55/Sqrt[55],
    -2574/Sqrt[1019304],
    -2574/Sqrt[1019304]};
  statePosition = Position[PsiPsipStates, {SL, SpLp}];
  If[Length[statePosition] == 0,
   Return[0],
   statePosition = statePosition[[1, 1]]];
  
  aizi2p = (
      T112PsiPsipStates[SL, SpLp]
     + t112PsiPsipStates[SL, SpLp]
     - a13/6 z13[[statePosition]]);
  aizi2p = SimplifyFun[aizi2p];
  Return[aizi2p];
  ];

aizi::usage="aizi[n, SL, SpLp] calculates the reduced matrix elements of the (spin-other-orbit + ECSO) operator for the f^n configuration corresponding to the terms SL and SpLp.

It does this by using equation (4) of \"Judd, BR, HM Crosswhite, and Hannah Crosswhite. “Intra-Atomic Magnetic Interactions for f Electrons.” Physical Review 169, no. 1 (1968): 130.\"

The values for the coefficients of fractional parentage originate from the tables of Velkov.
";
aizi[n_, SL_, SpLp_]:= Module[
  {s, l, nn, S, L, Sp, Lp, cfpSL, cfpSpLp, parentSL, Sb, Lb, parentSpLp, aizival},
  {s, l} = {1/2, 3};
  {S, L} = findSL[SL];
  {Sp, Lp} = findSL[SpLp];
  cfpSL   = CFP[{n, SL}];
  cfpSpLp = CFP[{n, SpLp}];

  aizival = Table[
    (
      parentSL = cfpSL[[nn, 1]];
      {Sb, Lb} = findSL[parentSL];
      dval = Sum[(
        cfpSpLp[[mm, 2]] 
        * SixJay[{S, 1, Sp}, {findSL[cfpSpLp[[mm, 1]]][[1]], s, Sb}] 
        * SixJay[{L, 1, Lp}, {findSL[cfpSpLp[[mm, 1]]][[2]], l, Lb}] 
        * AiZiTable[{n-1, parentSL, cfpSpLp[[mm, 1]]}]
      ),
        {mm, 2, Length[cfpSpLp]}
        ];
      phase = (-1)^(Sb+Lb+s+l+Sp+Lp);
      dval *= phase * cfpSL[[nn, 2]];
      dval
    ),
    {nn, 2, Length[cfpSL]}
  ];
  aizival = Total[aizival];
  aizival *= n/(n-2) * Sqrt[(2 S + 1) (2 Sp + 1) (2 L+ 1) (2Lp + 1)];
  Return[aizival];
];

GenerateAiZiTable[nmax_, export_:False, progressIndicator_:False]:= (
  If[progressIndicator,
  (
    numItersai = Association[Table[n->Length[AllowedNKSLterms[n]]^2, {n, 1, nmax}]];
    counters   = Association[Table[n->0, {n, 1, nmax}]];
    totalIters = Total[Values[numItersai[[1;;nmax]]]];
    template1  = StringTemplate["Iteration `numiter` of `totaliter`"];
    template2  = StringTemplate["`remtime` min remaining"];template3 = StringTemplate["Iteration speed = `speed` ms/it"];
    template4  = StringTemplate["Time elapsed = `runtime` min"];
    PrintTemporary[
      Dynamic[
        Pane[
          Grid[{
                {Superscript["f", n]},
                {template1[<|"numiter"->numiter, "totaliter"->totalIters|>]},
                {template4[<|"runtime"->Round[QuantityMagnitude[UnitConvert[(Now-startTime), "min"]], 0.1]|>]},
                {template2[<|"remtime"->Round[QuantityMagnitude[UnitConvert[(Now-startTime)/(numiter)*(totalIters-numiter), "min"]], 0.1]|>]},
                {template3[<|"speed"->Round[QuantityMagnitude[Now-startTime, "ms"]/(numiter), 0.01]|>]}, {ProgressIndicator[Dynamic[numiter], {1, totalIters}]}
              },
              Frame->All], Full, Alignment->Center]]];
  )
  ];
  AiZiTable = <||>;
  numiter = 1;
  startTime = Now;
  Do[
    (
      numiter+= 1;
      AiZiTable[{n, SL, SpLp}] = Which[
        n==1,
        0,
        n==2,
        SimplifyFun[AiZi2PsiPsipStates[SL, SpLp]],
        True,
        SimplifyFun[aizi[n,  SL, SpLp]]
      ];
    ),
    {n,    1, nmax},
    {SL,   AllowedNKSLterms[n]}, 
    {SpLp, AllowedNKSLterms[n]}
  ];
  If[export,
    (fname = FileNameJoin[{moduleDir, "data", "ReducedAiZiTable.m"}];
    Export[fname, AiZiTable];
    )
  ];
  Return[AiZiTable];
);

(* ############################ Reduced AiZi ############################### *)
(* ######################################################################### *)

(* ######################################################################### *)
(* ############################# T11 and T22 ############################### *)

T11n::usage="T11n[n, SL, SpLp] calculate the reduced matrix element of the T11 operator for the f^n configuration corresponding to the terms SL and SpLp. It is essentially the same as T22n with a different value of t.

It does this by using equation (4) of \"Judd, BR, HM Crosswhite, and Hannah Crosswhite. “Intra-Atomic Magnetic Interactions for f Electrons.” Physical Review 169, no. 1 (1968): 130.\"
";
T11n[n_, SL_, SpLp_]:= Module[{s, l, t, nn, S, L, Sp, Lp, cfpSL, cfpSpLp, parentSL, Sb, Lb, parentSpLp, Tnkk, dval, phase}, 
  {s, l} = {1/2, 3};
  {S, L} = findSL[SL];
  {Sp, Lp} = findSL[SpLp];
  t = 1;
  cfpSL = CFP[{n, SL}];
  cfpSpLp = CFP[{n, SpLp}];
  
  Tnkk = Table[
    (parentSL = cfpSL[[nn, 1]];
     {Sb, Lb} = findSL[parentSL];
     dval = Sum[(
        cfpSpLp[[mm, 2]]
         * SixJay[{S, t, Sp},
          {findSL[cfpSpLp[[mm, 1]]][[1]], s, Sb}]
         * SixJay[{L, t, Lp},
          {findSL[cfpSpLp[[mm, 1]]][[2]], l, Lb}]
         * T11Table[{n - 1, parentSL, cfpSpLp[[mm, 1]]}]),
       {mm, 2, Length[cfpSpLp]}];
     phase = (-1)^(Sb + Lb + s + l + Sp + Lp);
     dval *= phase * cfpSL[[nn, 2]];
     dval),
    {nn, 2, Length[cfpSL]}];
  Tnkk = Total[Tnkk];
  Tnkk *= n/(n - 2) * Sqrt[(2 S + 1) (2 Sp + 1) (2 L + 1) (2 Lp + 1)];
  Return[Tnkk];
  ];

T22n::usage="T22n[n, SL, SpLp] calculates the reduced matrix element of the T22 operator for the f^n configuration corresponding to the terms SL and SpLp.

It does this by using equation (4) of \"Judd, BR, HM Crosswhite, and Hannah Crosswhite. “Intra-Atomic Magnetic Interactions for f Electrons.” Physical Review 169, no. 1 (1968): 130.\"
";
T22n[n_, SL_, SpLp_]:= Module[{s, l, t, nn, S, L, Sp, Lp, cfpSL, cfpSpLp, parentSL, Sb, Lb, parentSpLp, Tnkk, dval, phase}, 
  {s, l} = {1/2, 3};
  {S, L} = findSL[SL];
  {Sp, Lp} = findSL[SpLp];
  t = 2;
  cfpSL   = CFP[{n, SL}];
  cfpSpLp = CFP[{n, SpLp}];
  
  Tnkk = Table[
    (parentSL = cfpSL[[nn, 1]];
     {Sb, Lb} = findSL[parentSL];
     dval = Sum[(
        cfpSpLp[[mm, 2]]
         * SixJay[{S, t, Sp},
          {findSL[cfpSpLp[[mm, 1]]][[1]], s, Sb}]
         * SixJay[{L, t, Lp},
          {findSL[cfpSpLp[[mm, 1]]][[2]], l, Lb}]
         * T22Table[{n - 1, parentSL, cfpSpLp[[mm, 1]]}]),
       {mm, 2, Length[cfpSpLp]}];
     phase = (-1)^(Sb + Lb + s + l + Sp + Lp);
     dval *= phase * cfpSL[[nn, 2]];
     dval),
    {nn, 2, Length[cfpSL]}];
  Tnkk = Total[Tnkk];
  Tnkk *= n/(n - 2) * Sqrt[(2 S + 1) (2 Sp + 1) (2 L + 1) (2 Lp + 1)];
  Return[Tnkk];
  ];

GenerateT11Table::usage="GenerateT11Table[nmax] calculates the reduced matrix elements of the double tensor operator T11 up to the configuration f^nmax. If the optional value \"Export\" is set to True, the resultant Association is saved to the data folder. If the optional value \"Progress Indicator\" is set to True then the function displays a progress indicator as the matrix elements are calculated.";
Options[GenerateT11Table]={"Export"->False, "Progress Indicator"->True};
GenerateT11Table[nmax_, OptionsPattern[]]:= (
  If[OptionValue["Progress Indicator"], (
    numItersai = Association[Table[n->Length[AllowedNKSLterms[n]]^2, {n, 1, nmax}]];
    counters = Association[Table[n->0, {n, 1, nmax}]];
    totalIters = Total[Values[numItersai[[1;;nmax]]]];
    template1 = StringTemplate["Iteration `numiter` of `totaliter`"];
    template2 = StringTemplate["`remtime` min remaining"];template3 = StringTemplate["Iteration speed = `speed` ms/it"];
    template4 = StringTemplate["Time elapsed = `runtime` min"];
    PrintTemporary[
      Dynamic[
        Pane[
          Grid[{{Superscript["f", n]},
          {template1[<|"numiter"->numiter, "totaliter"->totalIters|>]},
          {template4[<|"runtime"->Round[QuantityMagnitude[UnitConvert[(Now-startTime), "min"]], 0.1]|>]},
          {template2[<|"remtime"->Round[QuantityMagnitude[UnitConvert[(Now-startTime)/(numiter)*(totalIters-numiter), "min"]], 0.1]|>]},
          {template3[<|"speed"->Round[QuantityMagnitude[Now-startTime, "ms"]/(numiter), 0.01]|>]},
          {ProgressIndicator[Dynamic[numiter], {1, totalIters}]}
          },
          Frame->All],
          Full,
          Alignment->Center]
        ]
      ];
  )
  ];
  T11Table = <||>;
  numiter = 1;
  startTime = Now;
  Do[
    (
    numiter+= 1;
    T11Table[{n, SL, SpLp}] = Which[
    n==1,
    0,
    n==2,
    SimplifyFun[T112PsiPsipStates[SL, SpLp]],
    True,
    SimplifyFun[T11n[n,  SL, SpLp]]
    ];
    ),
    {n, 1, nmax},
    {SL, AllowedNKSLterms[n]},
    {SpLp, AllowedNKSLterms[n]}
  ];
  If[OptionValue["Export"],
    (
    fname = FileNameJoin[{moduleDir, "data", "ReducedT11Table.m"}];
    Export[fname, T11Table];
    )
  ];
  Return[T11Table];
);

GenerateT22Table[nmax_, export_, progressIndicator_]:= (
  If[progressIndicator,
  (
    numItersai = Association[Table[n->Length[AllowedNKSLterms[n]]^2, {n, 1, nmax}]];
    counters = Association[Table[n->0, {n, 1, nmax}]];
    totalIters = Total[Values[numItersai[[1;;nmax]]]];
    template1 = StringTemplate["Iteration `numiter` of `totaliter`"];
    template2 = StringTemplate["`remtime` min remaining"];template3 = StringTemplate["Iteration speed = `speed` ms/it"];
    template4 = StringTemplate["Time elapsed = `runtime` min"];
    PrintTemporary[
      Dynamic[
        Pane[
          Grid[{{Superscript["f", n]},
                {template1[<|"numiter"->numiter, "totaliter"->totalIters|>]},
                {template4[<|"runtime"->Round[QuantityMagnitude[UnitConvert[(Now-startTime), "min"]], 0.1]|>]},
                {template2[<|"remtime"->Round[QuantityMagnitude[UnitConvert[(Now-startTime)/(numiter)*(totalIters-numiter), "min"]], 0.1]|>]},
                {template3[<|"speed"->Round[QuantityMagnitude[Now-startTime, "ms"]/(numiter), 0.01]|>]},
                {ProgressIndicator[Dynamic[numiter], {1, totalIters}]}},
                Frame->All],
                Full,
                Alignment->Center]
              ]
            ];
  )
  ];
  T22Table = <||>;
  startTime = Now;
  numiter = 1;
  Do[
    (numiter+= 1;
    T22Table[{n, SL, SpLp}] = Which[
    n==1,
    0,
    n==2,
    SimplifyFun[T222PsiPsipStates[SL, SpLp]],
    True,
    SimplifyFun[T22n[n,  SL, SpLp]]
    ];),
    {n, 1, nmax},
    {SL, AllowedNKSLterms[n]}, 
    {SpLp, AllowedNKSLterms[n]}
  ];
  If[export,
  (
    fname = FileNameJoin[{moduleDir, "data", "ReducedT22Table.m"}];
    Export[fname, T22Table];
  )
  ];
  Return[T22Table];
);

(* ############################# T11 and T22 ############################### *)
(* ######################################################################### *)

(* ######################################################################### *)
(* ############################# Spin-Spin ################################# *)

SpinSpin::usage="SpinSpin[n, SL, SpLp, J] returns the matrix element <|SL,J|spin-spin|SpLp,J|> for the spin-spin operator within the configuration f^n. This matrix element is independent of MJ. This is obtained by querying the relevant reduced matrix element by querying the Association T22Table and putting in the adequate phase and 6-j symbol.

This is calculated according to equation (3) in \"Judd, BR, HM Crosswhite, and Hannah Crosswhite. “Intra-Atomic Magnetic Interactions for f Electrons.” Physical Review 169, no. 1 (1968): 130.\"
\".
";
SpinSpin[n_, SL_, SpLp_, J_] := Module[
  {S, L, Sp, Lp, \[Alpha], val},
   \[Alpha] = 2;
   {S, L} = findSL[SL];
   {Sp, Lp} = findSL[SpLp];
   val = ((-1)^(Sp + L + J) 
        * SixJay[{Sp, Lp, J}, {L, S, \[Alpha]}]
        * T22Table[{n, SL, SpLp}]);
   Return[val]
    ];

GenerateSpinSpinTable::usage="GenerateSpinSpinTable[nmax] generates the matrix elements in the |LSJ> basis for the (spin-other-orbit + electrostatically-correlated-spin-orbit) operator. It returns an Association where the keys are of the form {n, SL, SpLp, J}. If the option \"Export\" is set to True then the resulting object is saved to the data folder. Since this is a scalar operator, there is no MJ dependence. This dependence only comes into play when the crystal field contribution is taken into account.";
Options[GenerateSpinSpinTable] = {"Export"->False};
GenerateSpinSpinTable[nmax_, OptionsPattern[]] :=
 (
  SpinSpinTable = <||>;
  PrintTemporary[Dynamic[n]];
  Do[
   SpinSpinTable[{n, SL, SpLp, J}] = (SpinSpin[n, SL, SpLp, J]);,
   {n, 1, nmax},
   {J, minJ[n], maxJ[n]},
   {SL,   First /@ AllowedNKSLforJterms[n, J]},
   {SpLp, First /@ AllowedNKSLforJterms[n, J]}
   ];
  If[OptionValue["Export"],
   (fname = FileNameJoin[{moduleDir, "data", "SpinSpinTable.m"}];
    Export[fname, SpinSpinTable];
    )
   ];
  Return[SpinSpinTable];
  );

(* ############################# Spin-Spin ################################# *)
(* ######################################################################### *)


(* ######################################################################### *)
(* ##### Spin-Other-Orbit and Electrostatically-Correlated-Spin-Orbit ###### *)

SOOandECSO::usage="SOOandECSO[n, SL, SpLp, J] returns the matrix element <|SL,J|spin-spin|SpLp,J|> for the combined effects of the spin-other-orbit interaction and the electrostatically-correlated-spin-orbit (which originates from configuration interaction effects) within the configuration f^n. This matrix element is independent of MJ. This is obtained by querying the relevant reduced matrix element by querying the Association AiZiTable and putting in the adequate phase and 6-j symbol. The AiZiTable puts together the reduced matrix elements from three operators.

This is calculated according to equation (3) in \"Judd, BR, HM Crosswhite, and Hannah Crosswhite. “Intra-Atomic Magnetic Interactions for f Electrons.” Physical Review 169, no. 1 (1968): 130.\".
";
SOOandECSO[n_, SL_, SpLp_, J_]:= Module[
  {S, Sp, L, Lp, \[Alpha], val},
  \[Alpha] = 1;
  {S, L}   = findSL[SL];
  {Sp, Lp} = findSL[SpLp];
  val = ((-1)^(Sp + L + J) 
      * SixJay[{Sp, Lp, J}, {L, S, \[Alpha]}] 
      * AiZiTable[{n, SL, SpLp}]);
  Return[val];
]

Prescaling = {P2 -> P2/225, P4 -> P4/1089, P6 -> 25 P6 / 184041};
GenerateSOOandECSOTable::usage="GenerateSOOandECSOTable[nmax] generates the matrix elements in the |LSJ> basis for the (spin-other-orbit + electrostatically-correlated-spin-orbit) operator. It returns an Association where the keys are of the form {n, SL, SpLp, J}. If the option \"Export\" is set to True then the resulting object is saved to the data folder. Since this is a scalar operator, there is no MJ dependence. This dependence only comes into play when the crystal field contribution is taken into account.";
Options[GenerateSOOandECSOTable] = {"Export"->False}
GenerateSOOandECSOTable[nmax_, OptionsPattern[]]:= (
  SOOandECSOTable = <||>;
  Do[
    SOOandECSOTable[{n, SL, SpLp, J}] = (SOOandECSO[n, SL, SpLp, J] /. Prescaling);,
  {n, 1, nmax},
  {J, minJ[n], maxJ[n]},
  {SL,   First /@ AllowedNKSLforJterms[n, J]},
  {SpLp, First /@ AllowedNKSLforJterms[n, J]}
  ];
  If[OptionValue["Export"],
  (
    fname = FileNameJoin[{moduleDir, "data", "SOOandECSOTable.m"}];
    Export[fname, SOOandECSOTable];
  )
  ];
  Return[SOOandECSOTable];
);

(* ##### Spin-Other-Orbit and Electrostatically-Correlated-Spin-Orbit ###### *)
(* ######################################################################### *)

(* ######################################################################### *)
(* ############################ Crystal Field ############################## *)

Cqk::usage = "Cqk[nf, q, k, NKSL, J, M, NKSLp, Jp, Mp].";
Cqk[nf_, q_, k_, NKSL_, J_, M_, NKSLp_, Jp_, Mp_] := Module[
  {S, Sp, L, Lp, l, val},
  l = 3;
  {S, L} = findSL[NKSL];
  {Sp, Lp} = findSL[NKSLp];
  f1 = ThreeJay[{J, -M}, {k, q}, {Jp, Mp}];
  val = If[f1==0,
    0,
    (f2 = SixJay[{L, J, S}, {Jp, Lp, k}] ;
     If[f2==0,
      0,
      (
       f3 = ReducedUkTable[{nf, l, NKSL, NKSLp, k}];
       If[f3==0,
        0,
        (
         ( (-1)^(J - M + S + Lp + J + k)
           * Sqrt[(2 J + 1)*(2 Jp + 1)]
           * f1
           * f2  
           * f3
           * Ck[l, k]
         )
         )]
       )]
     )
    ];
  Cqk[nf, q, k, NKSL, J, M, NKSLp, Jp, Mp] = val
  ]

Bqk::usage = "";
Bqk[q_, 2] := {B02/2, B12, B22}[[q + 1]];
Bqk[q_, 4] := {B04/2, B14, B24, B34, B44}[[q + 1]];
Bqk[q_, 6] := {B06/2, B16, B26, B36, B46, B56, B66}[[q + 1]];

Sqk::usage = "";
Sqk[q_, 2] := {Sm22, Sm12, S02,  S12,  S22}[[q + 3]];
Sqk[q_, 4] := {Sm44, Sm34, Sm24, Sm14, S04,  S14,  S24, S34, S44}[[q + 5]];
Sqk[q_, 6] := {Sm66, Sm56, Sm46, Sm36, Sm26, Sm16, S06, S16, S26, S36, S46, S56, S66}[[q + 7]];

CrystalField::usage = "CrystalField[n, NKSL, J, M, NKSLp, Jp, Mp] gives the general expression for the matrix element of the crystal field Hamiltonian parametrized with Bqk and Sqk coefficients as a sum over spherical harmonics Cqk.";
CrystalField[n_, NKSL_, J_, M_, NKSLp_, Jp_, Mp_] := (
  Sum[(Bqk[q, k] * (Cqk[n, q, k, NKSL, J, M, NKSLp, Jp, Mp] 
                 + (-1)^q * Cqk[n, -q, k, NKSL, J, M, NKSLp, Jp, Mp]) 
    +I*Sqk[q, k] * (Cqk[n, q, k, NKSL, J, M, NKSLp, Jp, Mp] 
                 - (-1)^q * Cqk[n, -q, k, NKSL, J, M, NKSLp, Jp, Mp])
      ),
   {k, {2, 4, 6}},
   {q, 0, k}
     ]
  )

TotalCFIters::usage = 
  "TotalIters[i, j] returns total number of function evaluations for calculating all the matrix elements for the \!\(\*SuperscriptBox[\(f\), \(i\)]\) to the \!\(\*SuperscriptBox[\(f\), \(j\)]\) configurations.";
TotalCFIters[i_, j_] := (
  numIters = {196, 8281, 132496, 1002001, 4008004, 9018009, 11778624};
  Return[Total[numIters[[i ;; j]]]];
  )

GenerateCrystalFieldTable::usage = "GenerateCrystalFieldTable[nmax, export, progressIndicator, exportCompressed] computes the matrix values for the crystal field interaction for f^n configurations up to n = nmax. The function returns an Association whose keys are lists of the form {n, NKSL, J, M, NKSLp, Jp, Mp}. If export is set to True, then the result is exported to the data subfolder for the folder in which this package is in. If progressIndicator is set to True then an interactive progress indicator is shown. If exportCompressed is set to true the exported values are compressed.";
GenerateCrystalFieldTable[nmax_:7, export_:False, progressIndicator_:True, exportCompressed_:True]:= (
  Off[ClebschGordan::phy];
  Off[ClebschGordan::tri];
  Off[SixJSymbol::phy];
  Off[SixJSymbol::tri];
  CrystalFieldTables = <||>;
  ExportFun = If[exportCompressed,
    ExportMZip,
    Export
    ];
  numiter = 1;
  template1 = StringTemplate["Iteration `numiter` of `totaliter`"];
  template2 = StringTemplate["`remtime` min remaining"];
  template3 = StringTemplate["Iteration speed = `speed` ms/it"];
  template4 = StringTemplate["Time elapsed = `runtime` min"];
  totalIter = TotalCFIters[1, nmax];
  freebies = 0;
  startTime = Now;
  If[progressIndicator,
  PrintTemporary[
    Dynamic[
      Pane[
      Grid[{{Superscript["f", n]},
        {template1[<|"numiter" -> numiter, "totaliter" -> totalIter|>]},
        {template4[<|
            "runtime" -> Round[QuantityMagnitude[UnitConvert[(Now - startTime), "min"]], 0.1]|>]},
        {template2[<|
            "remtime" -> Round[QuantityMagnitude[UnitConvert[(Now - startTime)/(numiter - freebies) * (totalIter - numiter), "min"]], 0.1]|>]},
        {template3[<|
            "speed"   -> Round[QuantityMagnitude[Now - startTime, "ms"]/(numiter-freebies), 0.01]|>]},
        {ProgressIndicator[Dynamic[numiter], {1, totalIter}]}
        },
        Frame -> All
        ],
        Full,
        Alignment -> Center
    ]
    ]
    ];
  ];
  Do[
    (
    exportFname = FileNameJoin[{moduleDir, "data", "CrystalFieldTable_f"<>ToString[n]<>".m"}];
    If[FileExistsQ[exportFname],
      CrystalFieldTable = Import[exportFname];
      CrystalFieldTables[n] = CrystalFieldTable;
      Print["File exists, skipping calculation and importing file ..."];
      numiter+= TotalCFIters[n, n];
      freebies+= TotalCFIters[n, n];
      Continue[];
    ];
    CrystalFieldTable = Table[
      (numiter+= 1;
      {n, NKSL, J, M, NKSLp, Jp, Mp} -> CrystalField[n, NKSL, J, M, NKSLp, Jp, Mp]
      ),
      {J, minJ[n], maxJ[n]},
      {M, AllowedMforJ[J]},
      {Jp, minJ[n], maxJ[n]},
      {NKSL , First /@ AllowedNKSLforJterms[n, J]},
      {NKSLp, First /@ AllowedNKSLforJterms[n, Jp]},
      {Mp, AllowedMforJ[Jp]}
      ];
    CrystalFieldTable = Association[CrystalFieldTable];

    If[export, (
      Print["Exporting to file "<>ToString[exportFname]];
      ExportFun[exportFname, CrystalFieldTable];
      )
      ];
    CrystalFieldTables[n] = CrystalFieldTable;
    ),
  {n, 1, nmax}
  ];
  Return[CrystalFieldTables];
)

(* ############################ Crystal Field ############################## *)
(* ######################################################################### *)

(* ######################################################################### *)
(* ########### Configuration-Interaction via Casimir Operators ############# *)

GG2U::usage = "GG2U is an Association whose keys are labels for the irreducible representations of group G2 and whose values are the eigenvalues (multiplied by \[Beta]/12) of the corresponding Casimir operator.
Reference: Wybourne, \"Spectroscopic Properties of Rare Earths\", table 2-6.";
GG2U = Association[{
    "00" -> 0,
    "10" -> \[Beta]/12*6 ,
    "11" -> \[Beta]/12*12,
    "20" -> \[Beta]/12*14,
    "21" -> \[Beta]/12*21,
    "22" -> \[Beta]/12*30,
    "30" -> \[Beta]/12*24,
    "31" -> \[Beta]/12*32,
    "40" -> \[Beta]/12*36}
   ];

GR7W::usage = "GR7W is an Association whose keys are labels for the irreducible representations of group R7 and whose values are the eigenvalues (multiplied by \[Gamma] / 5) of the corresponding Casimir operator.
Reference: Wybourne, \"Spectroscopic Properties of Rare Earths\", table 2-7.";
GR7W := Association[
   {"000" -> 0,
    "100" -> \[Gamma]/5*3,
    "110" -> \[Gamma]/5*5,
    "111" -> \[Gamma]/5*6,
    "200" -> \[Gamma]/5*7,
    "210" -> \[Gamma]/5*9,
    "211" -> \[Gamma]/5*10,
    "220" -> \[Gamma]/5*12,
    "221" -> \[Gamma]/5*13,
    "222" -> \[Gamma]/5*15}
   ];

GG2NKSL::usage = "GG2NKSL[SL, SpLp] returns the matrix element of the configuration interaction term corresponding to the Casimir operator G2, for the two given terms SL and SpLp.";
GG2NKSL[SL_, SpLp_] := (
 If[SL==SpLp,
 GG2U[findNKLSterm[SL][[1]][[4]]],
 0]
)

GR7NKSL::usage = "GR7NKSL[SL, SpLp] returns the matrix element of the configuration interaction term corresponding to the Casimir operator R7, for the two given terms SL and SpLp.";
GR7NKSL[SL_, SpLp_] := (
 If[SL==SpLp,
 GR7W[findNKLSterm[SL][[1]][[3]]],
 0]
)

TwoBodyNKSL::usage = "TwoBodyNKSL[SL, SpLp] returns the matrix element for configuration interaction as approximated by the Cassimir operators of the groups SO(3), G2, and R7. SL and SpLp are strings that represent terms under LS coupling.";
TwoBodyNKSL[SL_, SpLp_] := Module[
  {S, L},
  {S, L} = findSL[SL];
  val = (If[SL==SpLp,
    \[Alpha] * L * (L + 1) +
    GR7W[Part[First[findNKLSterm[SL]], 3]] + 
    GG2U[Part[First[findNKLSterm[SL]], 4]],
    0]);
  TwoBodyNKSL[S, L] = val;
  Return[val];
  ]

(* ########### Configuration-Interaction via Casimir Operators ############# *)
(* ######################################################################### *)

MagneticInteractions[{n_, SLJ_, SLJp_, J_}, chenDelta_, includeSpinSpin_] := (
  ss = If[includeSpinSpin,
    SpinSpinTable[{n, SLJ, SLJp, J}],
    0];
  sooandecso = SOOandECSOTable[{n, SLJ, SLJp, J}];
  total = ss + sooandecso;
  key = {n, SLJ, SLJp, J};
  If[chenDelta,
    If[
      MemberQ[Keys[chenDeltas], {n, SLJ, SLJp}],
    (
      {S, L} = findSL[SLJ];
      {Sp, Lp} = findSL[SLJp];
      phase = (-1)^(Sp + L + J);
      Msixjay = SixJay[{Sp, Lp, J},{L, S, 2}];
      Psixjay = SixJay[{Sp, Lp, J},{L, S, 1}];
      {M0v, M2v, M4v, P2v, P4v, P6v} = chenDeltas[{n, SLJ, SLJp}]["wrong"];
      total = phase*Msixjay(M0v*M0+M2v*M2+M4v*M4);
      total +=  phase*Psixjay(P2v*P2+P4v*P4+P6v*P6);
      total = total /. Prescaling;
    )
    ]
    ];
  Return[total];
)

Options[EnergyMatrix] = {"Sparse"->True, "ChenDeltas"->False, "Include Spin-Spin"->True};
EnergyMatrix::usage = "EnergyMatrix[n, J, J', I, I'] provides the matrix element <J, I|H|J', I'> within the f^n configuration, it does this by adding the following interactions: Coulomb, spin-orbit, spin-other-orbit, electrostatically-correlated-spin-orbit, spin-spin, three-body interactions, and crystal-field.";
EnergyMatrix[n_, J_, Jp_, Ii_, Ip_, CFTable_, OptionsPattern[]]:= (
  eMatrix = Table[
    (*Condition for a scalar matrix op*)
    subKron = (KroneckerDelta[NKSLJM[[4]], NKSLJMp[[4]]]
            * KroneckerDelta[J, Jp]
            * KroneckerDelta[NKSLJM[[3]], NKSLJMp[[3]]]);
    matValue = If[subKron==0,
      0,
      SimplifyFun[
        (ElectrostaticMatrixTable[{n, NKSLJM[[1]], NKSLJMp[[1]]}]
          + TwoBodyNKSL[NKSLJM[[1]], NKSLJMp[[1]]]
          + SpinOrbitTable[{n, NKSLJM[[1]], NKSLJMp[[1]], NKSLJM[[2]]}]
          + MagneticInteractions[{n, NKSLJM[[1]], NKSLJMp[[1]], NKSLJM[[2]]}, 
                  OptionValue["ChenDeltas"], OptionValue["Include Spin-Spin"]]
          + ThreeBodyTable[{n, NKSLJM[[1]], NKSLJMp[[1]]}]
        )
      ]
    ];
    matValue += CFTable[{n, NKSLJM[[1]], NKSLJM[[2]], NKSLJM[[3]], NKSLJMp[[1]], NKSLJMp[[2]], NKSLJMp[[3]]}];
    matValue
    ,
  {NKSLJMp, Partition[Flatten[AllowedNKSLJMIMforJIterms[n, Jp, Ip]], 4]},
  {NKSLJM , Partition[Flatten[AllowedNKSLJMIMforJIterms[n, J, Ii]], 4]}
  ];
If[OptionValue["Sparse"],
 eMatrix = SparseArray[eMatrix]];
Return[eMatrix]
)

EnergyStates[n_, J_, Ii_]:= AllowedNKSLJMIMforJIterms[n, J, Ii];

Options[TabulateEnergyMatrixTable] = {"Sparse"->True, "ChenDeltas"->False, "Include Spin-Spin"->True};
TabulateEnergyMatrixTable::usage = "TabulateEnergyMatrixTable[n, I] returns a list with three elements {EnergyMatrixTable, EnergyStatesTable, AllowedM}. EnergyMatrixTable is an Association with keys equal to lists of the form {n, J, Jp, Ii, Ii}. EnergyStatesTable is an Association with keys equal to lists of the form {n, J, Ii}. AllowedM is another Association with keys equal to lists of the form {n, J} and values equal to lists equal to the corresponding values of MJ. It's unnecessary (and it won't work in this implementation) to give n > 7 given the equivalency between electron and hole configurations.";
TabulateEnergyMatrixTable[n_, Ii_, CFTable_, OptionsPattern[]]:= (
  EnergyMatrixTable = <||>;
  EnergyStatesTable = <||>;
  AllowedM = <||>;
  totalIterations = Length[AllowedJ[n]]^2;
  template1 = StringTemplate["Iteration `numiter` of `totaliter`"];
  template2 = StringTemplate["`remtime` min remaining"];
  template4 = StringTemplate["Time elapsed = `runtime` min"];
  numiter   = 0;
  startTime = Now;
  temp = PrintTemporary[
    Dynamic[
      Grid[{
          {template1[<|"numiter"->numiter, "totaliter"->totalIterations|>]},
          {template2[<|"remtime"->Round[QuantityMagnitude[UnitConvert[(Now-startTime)/(Max[1,numiter])*(totalIterations-numiter), "min"]], 0.1]|>]},
          {template4[<|"runtime"->Round[QuantityMagnitude[UnitConvert[(Now-startTime), "min"]], 0.1]|>]},
          {ProgressIndicator[numiter,{1,totalIterations}]}
    }]
  ]
  ];
  Do[
    (EnergyMatrixTable[{n, J, Jp, Ii, Ii}] = EnergyMatrix[n, J, Jp, Ii, Ii, CFTable, "Sparse"->OptionValue["Sparse"], "ChenDeltas"->OptionValue["ChenDeltas"], "Include Spin-Spin"->OptionValue["Include Spin-Spin"]];
     EnergyStatesTable[{n, J, Ii}]  = EnergyStates[n, J, Ii];
     AllowedM[{n, J}]               = Table[M, {J, minJ[n], maxJ[n]}, {M, -J, J}];
     numiter += 1;
    ),
    {Jp, AllowedJ[n]},
    {J, AllowedJ[n]}
  ];
  NotebookDelete[temp];
  Return[{EnergyMatrixTable, EnergyStatesTable, AllowedM}];
)

Options[TabulateManyEnergyMatrixTables] = {"Overwrite"->False, "Sparse"->True, "ChenDeltas"->False, "Include Spin-Spin"->True};
TabulateManyEnergyMatrixTables::usage = "TabulateManyEnergyMatrixTables[{n1, n2, ...}, {I1, I2, ...}] calculates the tables of matrix elements for the requested f^n_i configurations with the given nuclear spin I_i. The function does not return the matrices themselves. It instead returns an Association whose keys are lists of the form {n, I} and whose values are the filenames where the output of TabulateEnergyMatrixTables was saved to. When these files are loaded with Get, the following three symbols are thus defined: EnergyMatrixTable, EnergyStatesTable, and AllowedM.
EnergyMatrixTable is an Association whose keys are of the form {n, J, Jp, Ii, Ii} and whose values are matrix elements.";
TabulateManyEnergyMatrixTables[ns_, Iis_, OptionsPattern[]]:= (
  overwrite = OptionValue["Overwrite"];
  fNames = <||>;
  Do[
    (
    CFdataFilename = FileNameJoin[{moduleDir, "data", "CrystalFieldTable_f"<>ToString[n]<>".zip"}];
    PrintTemporary["Importing CrystalFieldTable from ", CFdataFilename, " ..."];
    CrystalFieldTable = ImportMZip[CFdataFilename];
    Do[(
      PrintTemporary["#------- n = ", n, " | I = ", Ii, " -------#"];
      exportFname = EnergyMatrixFileName[n, Ii];
      fNames[{n, Ii}] = exportFname;
      If[FileExistsQ[exportFname] && Not[overwrite],
        Continue[]];
      {EnergyMatrixTable, EnergyStatesTable, AllowedM} = TabulateEnergyMatrixTable[n, Ii, CrystalFieldTable, "Sparse"->OptionValue["Sparse"], "ChenDeltas"->OptionValue["ChenDeltas"], "Include Spin-Spin"->OptionValue["Include Spin-Spin"]];
      If[FileExistsQ[exportFname]&&overwrite,
        DeleteFile[exportFname]];
      Save[exportFname, {EnergyMatrixTable, EnergyStatesTable, AllowedM}];
      ),
    {Ii, Iis}
    ];
    ClearAll[CrystalFieldTable];
    ),
  {n, ns}
  ];
Return[fNames];
)

HamMatrixAssembly::usage="HamMatrixAssembly[n, IiN] returns the Hamiltonian matrix for the f^n_i configuration with nuclear spin I_i. The matrix is returned as a SparseArray."

HamMatrixAssembly[nf_, IiN_] := Module[
  {n, ii, jj, JMvals},
  (*#####################################*)
  (*hole-particle equivalence enforcement*)
  n = nf;
  numFittedLevels = Total[If[#, 1, 0] & /@ picker];
  allVars = {E0, E1, E2, E3, \[Zeta], F0, F2, F4, F6, M0, M2, M4, T2, 
    T3, T4, T6, T7, T8, P0, P2, P4, P6, gs, 
    gI, \[Alpha], \[Beta], \[Gamma], B02, B04, B06, B12, B14, B16, 
    B22, B24, B26, B34, B36, B44, B46, B56, B66, S12, S14, S16, S22, 
    S24, S26, S34, S36, S44, S46, S56, S66, T11, T12, T14, T15, T16, 
    T17, T18, T19};
  params0 = AssociationThread[allVars, allVars];

  (*hole-particle equivalence*)
  If[nf > 7, (n = 14 - nf;
    params = HoleElectronConjugation[params0];),
    params = params0;];

  (*Load symbolic expressions for energy sub-matrices.*)
  Get[EnergyMatrixFileName[n, IiN]];
  (*Patch together the entire matrix representation in block-
  diagonal form.*)
  EnergyMatrix = 
   ConstantArray[0, {Length[AllowedJ[n]], Length[AllowedJ[n]]}];
  Do[EnergyMatrix[[jj, ii]] = 
     EnergyMatrixTable[{n, AllowedJ[n][[ii]], AllowedJ[n][[jj]], IiN, 
       IiN}];,
       {ii, 1, Length[AllowedJ[n]]},
       {jj, 1, Length[AllowedJ[n]]}
   ];
  EnergyMatrix = ArrayFlatten[EnergyMatrix];
  EnergyMatrix = ReplaceInSparseArray[EnergyMatrix, params];
  Return[EnergyMatrix];
  ]

OperatorsMatrix::usage="OperatorsMatrix[n, Ii, ops] returns the matrix representation of the sum of the operators in the list ops for the f^n configuration with nuclear spin I_i.
ops may be any subset of: {\"Electrostatic\", \"Two-body CI\", \"Spin-orbit\", \"Spin-other-orbit and electrostatically-correlated-spin-orbit\", \"Spin-spin\", \"Three-body CI\"}";
OperatorsMatrix[n_,Ii_,ops_]:=
(
operatorCalls=<|
    "Electrostatic"->"ElectrostaticMatrixTable[{"<>ToString[n]<>",NKSLJM[[1]],NKSLJMp[[1]]}]",
    "Two-body CI"->"TwoBodyNKSL[NKSLJM[[1]],NKSLJMp[[1]]]",
    "Spin-orbit"->"SpinOrbitTable[{"<>ToString[n]<>",NKSLJM[[1]],NKSLJMp[[1]],NKSLJM[[2]]}]",
    "Spin-other-orbit and electrostatically-correlated-spin-orbit"->"SOOandECSOTable[{"<>ToString[n]<>",NKSLJM[[1]],NKSLJMp[[1]],NKSLJM[[2]]}]",
    "Spin-spin"->"SpinSpinTable[{"<>ToString[n]<>",NKSLJM[[1]],NKSLJMp[[1]],NKSLJM[[2]]}]",
    "Three-body CI"->"ThreeBodyTable[{"<>ToString[n]<>",NKSLJM[[1]],NKSLJMp[[1]]}]"
    |>;
opMatrixTemplate=If[MemberQ[ops,"Crystal-field"],
If[Length[ops]>1,
  StringTemplate["opMatrix = Table[
  (
  (*Condition for a scalar matrix op*)
  subKron=(
    KroneckerDelta[NKSLJM[[4]],NKSLJMp[[4]]]
  * KroneckerDelta[`J`,`Jp`]
  * KroneckerDelta[NKSLJM[[3]],NKSLJMp[[3]]]);
  matValue = If[subKron==0, 0, SimplifyFun[(`iteratorSum`)]];
  matValue += CrystalFieldTable[{"<>ToString[n]<>", NKSLJM[[1]], NKSLJM[[2]], NKSLJM[[3]], NKSLJMp[[1]], NKSLJMp[[2]], NKSLJMp[[3]]}];
      matValue
  ),
  {NKSLJMp, Partition[Flatten[AllowedNKSLJMIMforJIterms["<>ToString[n]<>",`Jp`,`Ip`]],4]},
  {NKSLJM,  Partition[Flatten[AllowedNKSLJMIMforJIterms["<>ToString[n]<>",`J`,`Ii`]],4]}];"],
  StringTemplate["opMatrix = Table[
  (
  CrystalFieldTable[{"<>ToString[n]<>", NKSLJM[[1]], NKSLJM[[2]], NKSLJM[[3]], NKSLJMp[[1]], NKSLJMp[[2]], NKSLJMp[[3]]}]
  ),
  {NKSLJMp, Partition[Flatten[AllowedNKSLJMIMforJIterms["<>ToString[n]<>",`Jp`,`Ip`]],4]},
  {NKSLJM,  Partition[Flatten[AllowedNKSLJMIMforJIterms["<>ToString[n]<>",`J`,`Ii`]],4]}];"]
],
  StringTemplate["opMatrix = Table[
  (
  (*Condition for a scalar matrix op*)
  subKron=(
    KroneckerDelta[NKSLJM[[4]],NKSLJMp[[4]]]
  * KroneckerDelta[`J`,`Jp`]
  * KroneckerDelta[NKSLJM[[3]],NKSLJMp[[3]]]);
  matValue = If[subKron==0, 0, SimplifyFun[(`iteratorSum`)]];
  matValue
  ),
  {NKSLJMp, Partition[Flatten[AllowedNKSLJMIMforJIterms["<>ToString[n]<>",`Jp`,`Ip`]],4]},
  {NKSLJM,  Partition[Flatten[AllowedNKSLJMIMforJIterms["<>ToString[n]<>",`J`,`Ii`]],4]}];"]
];

If[MemberQ[ops,"Crystal-field"],
  (
  CFdataFilename=FileNameJoin[{moduleDir,"data","CrystalFieldTable_f"<>ToString[n]<>".zip"}];PrintTemporary["Importing CrystalFieldTable from ",CFdataFilename," ..."];
  CrystalFieldTable=ImportMZip[CFdataFilename];
  operators=Delete[ops,Position[ops,"Crystal-field"][[1,1]]];
  ),
  operators=ops;
  ];

OperatorMatrixAssoc=<||>;
If[Length[operators]!=0,
  (
  iteratorSum=StringJoin[Riffle[operatorCalls[#]&/@operators," + "]];
  OperatorMatrix[J_,Jp_,Ip_,CFTable_:0]:=(
  tableExpression=opMatrixTemplate[<|"iteratorSum"->iteratorSum,
    "J"->ToString[J,InputForm],
    "Jp"->ToString[Jp,InputForm],
    "Ii"->ToString[Ii,InputForm],
    "Ip"->ToString[Ip,InputForm]|>];
  ToExpression[tableExpression];
  Return[opMatrix]
  );
  ),
  OperatorMatrix[J_,Jp_,Ip_,CFTable_:0]:=(
    tableExpression=opMatrixTemplate[<|"J"->ToString[J,InputForm],
      "Jp"->ToString[Jp,InputForm],
      "Ii"->ToString[Ii,InputForm],
      "Ip"->ToString[Ip,InputForm]|>];
    ToExpression[tableExpression];
    Return[opMatrix]
  );
];

Do[(
  OperatorMatrixAssoc[{n,J,Jp,Ii,Ii}]=OperatorMatrix[J,Jp,Ii];
  ),
{Jp,AllowedJ[n]},
{J,AllowedJ[n]}];

operatorMatrix=ConstantArray[0,{Length[AllowedJ[n]],Length[AllowedJ[n]]}];
Do[operatorMatrix[[jj,ii]]=OperatorMatrixAssoc[{n,AllowedJ[n][[ii]],AllowedJ[n][[jj]],Ii,Ii}];
  ,{ii,1,Length[AllowedJ[n]]}
  ,{jj,1,Length[AllowedJ[n]]}];
operatorMatrix=ArrayFlatten[operatorMatrix];
Return[operatorMatrix]
)

(* ####################################### Table Generation Functions ########################### *)
(* ############################################################################################## *)

(* ############################################################################################## *)
(* ########################################## Printers ########################################## *)

printL[L_] := 
 If[StringQ[L], L, StringTake[specAlphabet, {L + 1}]]

findSL[SL_]:= (
  findSL[SL] = If[StringQ[SL],
    {(ToExpression[StringTake[SL, 1]]-1)/2,
     StringPosition[specAlphabet, StringTake[SL, {2}]][[1, 1]]-1},
    SL]
    )

findSL\[Gamma][SL_] := 
 If[StringQ[SL], 
  If[StringLength[SL]==2, 
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

RSJnoIlabelPercent::usage = "Needs documentation and perhaps some options.";
RSJnoIlabelPercent[\[CapitalPsi]_]:= Module[
  {magFun, allRS, label},
  magFun = DisplayForm[
            SubscriptBox[
              StringTake[#[[1, 1]][[1]], 2],
              {#[[1, 1, 2]], #[[1]][[-1]], #[[-1]]}
              ]
              ]&;
  magRS[a_, b_]:= {Abs[a]^2, Map[magFun, b]};
  RSfunc = {Total[#][[1]], Flatten[#[[-1, -1]]]} &;
  allRS = Map[RSfunc, 
    GatherBy[
      Thread[magRS[\[CapitalPsi][[2]], \[CapitalPsi][[3]]]], 
      #[[2]] &]];
  allRS = SortBy[allRS, -#[[1]] &];
  label = If[Part[allRS, 1, 1]>= 0.5,
    DisplayForm[RowBox[
                  {
                  allRS[[1, 2]],
                  StringJoin["(", ToString[Round[100 * allRS[[1, 1]]]], "%)"]
                  }
                  ]
                ],
    DisplayForm[RowBox[
                  {
                  allRS[[1, 2]],
                  StringJoin["(", ToString[Round[100 * allRS[[1, 1]]]], "%) + "],
                  allRS[[2, 2]],
                  StringJoin["(", ToString[Round[100 * allRS[[2, 1]]]], "%) + "],
                  allRS[[3, 2]],
                  StringJoin["(", ToString[Round[100 * allRS[[3, 1]]]], "%) + "],
                  allRS[[4, 2]],
                  StringJoin["(", ToString[Round[100 * allRS[[4, 1]]]], "%) + "],
                  allRS[[5, 2]],
                  StringJoin["(", ToString[Round[100 * allRS[[5, 1]]]], "%)"]
                  }
                ]
              ]
          ];
  Return[label];
];

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

NKLSterms[n_, L_, S_] := Module[
  {NKterms, SL},
  SL = StringJoin[ToString[2 S + 1], printL[L]];
  NKterms = {{}};
  Map[If[! StringFreeQ[First[#], SL], 
     If[ToExpression[Part[#, 2]] <= n, 
      NKterms = Join[NKterms, {#}, 1]]] &, fNterms];
  NKterms = DeleteCases[NKterms, {}];
  If[Length[NKterms]==1, 
   NKterms[[1]][[1]] = StringTake[Part[NKterms, 1, 1], 2]];
  NKterms]

findNKLSterm[SL_] := Module[
  {NKterms, n},
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
CFPterms[n_, SL_] := Module[
  {NKterms, CFPconfig},
  NKterms = {{}};
  CFPconfig = Part[CFPtable, n];
  Map[If[StringFreeQ[First[#], SL], Null , 
     NKterms = Join[NKterms, {#}, 1]] &, CFPconfig];
  NKterms = DeleteCases[NKterms, {}]]
CFPterms[n_, L_, S_] := Module[
  {NKterms, SL, CFPconfig},
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


AllowedSLterms[n_] := Map[findSL[First[#]] &, CFPterms[Min[n,14-n]]]

AllowedNKSLterms[n_] := Map[First, CFPterms[Min[n,14-n]]]

maxJ[n_] := Max[Map[Total, AllowedSLterms[Min[n,14-n]]]]

minJ[n_] := Min[Map[Abs[Part[#, 1] - Part[#, 2]] &, AllowedSLterms[Min[n,14-n]]]]

AllowedSLJterms[n_] := Module[
  {allowedSL, allowedSLJ},
    allowedSL = AllowedSLterms[n];
    allowedSLJ = {}; 
    For[nn = 1, nn <= Length[allowedSL],
      termSL = allowedSL[[nn]];
      termsSLJ = 
    Table[{termSL[[1]], termSL[[2]], J}, {J, 
      Abs[termSL[[1]] - termSL[[2]]], Total[termSL]}];
      allowedSLJ = Join[allowedSLJ, termsSLJ];
      nn++];
    SortBy[allowedSLJ, Last]]

AllowedNKSLJterms[n_] := Module[
  {allowedSL, allowedNKSL, allowedSLJ},
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

AllowedNKSLforJterms[n_, J_] := Module[
    {allowedSL, allowedNKSL, allowedSLJ, nn, termSL, termNKSL, termsSLJ},
    allowedNKSL = AllowedNKSLterms[n];
    allowedSL = AllowedSLterms[n];
    allowedSLJ = {}; 
    For[nn = 1,
      nn <= Length[allowedSL],
      termSL = allowedSL[[nn]];
      termNKSL = allowedNKSL[[nn]];
      termsSLJ = If[Abs[termSL[[1]] - termSL[[2]]] <= J <= Total[termSL], 
          {{termNKSL, J}},
          {}
          ];
      allowedSLJ = Join[allowedSLJ, termsSLJ];
      nn++
      ];
    Return[allowedSLJ]
  ];

AllowedSLJMterms[n_] := Module[
    {allowedSLJ, allowedSLJM, termSLJ, termsSLJM, nn},
    allowedSLJ = AllowedSLJterms[n];
    allowedSLJM = {}; 
    For[nn = 1,
      nn <= Length[allowedSLJ],
      nn++,
      termSLJ = allowedSLJ[[nn]];
      termsSLJM = 
    Table[{termSLJ[[1]], termSLJ[[2]], termSLJ[[3]], 
      M}, {M, - termSLJ[[3]], termSLJ[[3]]}];
      allowedSLJM = Join[allowedSLJM, termsSLJM];
      ];
    Return[SortBy[allowedSLJM, Last]];
    ]

AllowedNKSLJMforJMterms[n_, J_, MJ_] := Module[
  {allowedSL, allowedNKSL, allowedSLJM},
  allowedNKSL = AllowedNKSLterms[n];
  allowedSL   = AllowedSLterms[n];
  allowedSLJM = {}; 
  For[nn = 1,
    nn <= Length[allowedSL],
    termSL = allowedSL[[nn]];
    termNKSL = allowedNKSL[[nn]];
    termsSLJ = If[(Abs[termSL[[1]] - termSL[[2]]] 
                   <= J 
                   <= Total[termSL] 
                   && (Abs[MJ] <= J)
                  ),
                  {{termNKSL, J, MJ}},
                  {}];
    allowedSLJM = Join[allowedSLJM, termsSLJ];
    nn++];
    Return[allowedSLJM];
    ]

AllowedNKSLJMforJterms[n_, J_] := Module[
  {allowedSL, allowedNKSL, allowedSLJM, Mvals},
  allowedNKSL = AllowedNKSLterms[n];
  allowedSL = AllowedSLterms[n];
  allowedSLJM = {}; 
  Mvals = Range[-J, J];
  For[nn = 1,
    nn <= Length[allowedSL],
    termSL = allowedSL[[nn]];
    termNKSL = allowedNKSL[[nn]];
    termsSLJ = If[Abs[termSL[[1]] - termSL[[2]]] <= J <= 
              Total[termSL], {{termNKSL, J}}, {}
              ];
    allowedSLJM = Join[allowedSLJM, termsSLJ];
    nn++];
  Return[Tuples[{allowedSLJM, Mvals}]];
  ]

AllowedMforJ[j_] := Range[-j, j, 1];

AllowedNKSLJMIMforJIterms[n_, J_, \[CapitalIota]_] := 
  Tuples[{AllowedNKSLJMforJterms[n, J], Table[MI, {MI, -\[CapitalIota], \[CapitalIota]}]}];

AllowedJ[n_] := Table[J, {J, minJ[n], maxJ[n]}];

(* ############################### Allowed SL, SLJ, and SLJM Terms ############################## *)
(* ############################################################################################## *)

(* ############################################################################################## *)
(* ###############################              SOLVERS            ############################## *)

ShiftedLevels[originalLevels_] := Module[
  {groundEnergy, shifted},
  groundEnergy = Sort[originalLevels][[1,1]];
  shifted      = Map[{#[[1]] - groundEnergy, #[[2]]} &, originalLevels];
  Return[shifted];
  ]

LoadGuillotParameters::usage = "LoadGuillotParameters[Ln] loads into the current session the parameters for \!\(\*SuperscriptBox[\(Ln\), \(\(3\)\(+\)\)]\). In addition of defining values inside of the session, the function also returns a list with the following quantities : {nf, E0, E1 , E2, E3 , \[Zeta], \[Alpha], \[Beta], \[Gamma], T2, T3, T4, T6, T7, T8, M0, M2, M4, P0, P2, P4, P6, B02, B06}";
LoadGuillotParameters::argerr = "Not a lanthanide."

LoadGuillotParameters[Ln_]:= Module[
  {paramAssoc},
  Which[
    (* The hard-coded values used here are from Table I, "Spectra of lanthanides", Carnall et al, 1989.*)
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
      nf = 2; Fsup2 = 68878; Fsup4 = 50347; Fsup6 = 32901;
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
      nf = 4; Fsup2 = 76400;Fsup4 = 54900; Fsup6 = 37700;
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
      luFname = DialogInput[Prepend[MapThread[Button[#1, DialogReturn[#2]]&, {shortNames, luFnames}], TextCell[Style["These are the files available for Yb, choose one:", 12]]]];
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
  Return[{nf, E0, E1 , E2, E3 , \[Zeta], \[Alpha], \[Beta], \[Gamma],
  T2, T3, T4, T6, T7, T8, M0, M2, M4, P0, P2, P4, P6, B02, B06}];
];

Options[LoadParameters] = {
    "Source"->"Carnall",
    "Free Ion"->False,
    "gs"->2.002319304386,
    "gI"->0.987};
LoadParameters[Ln_, OptionsPattern[]]:= Module[
  {source, params},
  (
    source = OptionValue["Source"];
    params = Which[source=="Carnall",
            (Association[Carnall["data"][Ln]])
            ];
    (*If a free ion then all the parameters from the crystal field are set to zero*)
    If[OptionValue["Free Ion"],
      Do[params[cfSymbol] = 0, {cfSymbol, cfSymbols}]
    ];
    params[F0] = 0;
    params[M2] = 0.56 * params[M0]; (* See Carnall 1989, Table I, caption, probably fixed based on HF values*)
    params[M4] = 0.31 * params[M0]; (* See Carnall 1989, Table I, caption, probably fixed based on HF values*)
    params[P0] = 0;
    params[P4] = 0.5 * params[P2];  (* See Carnall 1989, Table I, caption, probably fixed based on HF values*)
    params[P6] = 0.1 * params[P2];  (* See Carnall 1989, Table I, caption, probably fixed based on HF values*)
    params[gs] = OptionValue["gs"];
    params[gI] = OptionValue["gI"];
    params[nf] = Round[params[nf]];
    {params[E0], params[E1], params[E2], params[E3]} = FtoE[params[F0], params[F2], params[F4], params[F6]];
    params[E0] = 0;
    Return[params];
  )
];


EnergyMatrixFileName::usage = "EnergyMatrixFileName[nf, IiN] gives the filename for the energy matrix table for an atom with nf electrons and a nucleus of spin IiN.";
EnergyMatrixFileName[n_, IiN_] := (
  fname = FileNameJoin[{moduleDir,
      "hams",
      StringJoin[{"f", ToString[n], "_I_", ToString[2*IiN + 1], "_EnergyMatrixTable.m"}]}];
  Return[fname];
  );

HoleElectronConjugation::usage = "HoleElectronConjugation[params] takes the parameters that define a configuration and converts them so that they may be interpreted under its conjugate configuration.";

HoleElectronConjugation[params_] := (
  newParams = params;
  flipSignsOf = {\[Zeta], T2, T3, T4, T6, T7, T8};
  flipSignsOf = Join[flipSignsOf, cfSymbols];
  Do[
    (newParams[flipper] = -newParams[flipper]), 
  {flipper, flipSignsOf}];
  Return[newParams];
  )


SolveStates::usage = "SolveStates[nf, IiN, params] solves the energy values and states for an atom with n f-electrons with a nucleus of spin IiN. params is an Association with the parameters of the specific ion under study.
This function requires files for pre-computed energy matrix tables that provide the symbols EnergyMatrixTable[_, _, _, _, _].
The optional parameter \"maxEigenvalues\" (default: \"All\") specifies the number of eigenvalues to be returned. If maxE is \"All\" then all eigenvalues are returned. If maxE is positive then the k largest (in absolute value) eigenvalues are returned. If maxE is negative then the k smallest (in absolute value) eigenvalues are returned.
To account for configurations f^n with n > 7, particle-hole dualities are enforced for \[Zeta] and T_i.
The unit for the returned energies is cm^-1.

Parameters
----------
nf (int) : Number of f-electrons.
IiN (int) : Nuclear spin.
params (Association) : Parameters of the ion under study.

Returns
-------
{eigenstates, basis} (list): eigenstates is a list wher each element is a list with two elements, the first element being the energy eigenvalue and the second being a list that represents the eigenvector in the computational basis. basis is a list of lists that represent the computational basis. The elements of the basis are lists of the form {{{SL, J}, mJ}, I}, where SL is given a string.

Options
-------
\"Return Symbolic Matrix\" (bool) : If True then the function returns instead a list with the three elements {levels, basis, symbolicMatrix}.
\"maxEigenvalues\" (int) : Number of eigenvalues to be returned. If \"All\" then all eigenvalues are returned. If positive then the k largest (in absolute value) eigenvalues are returned. If negative then the k smallest (in absolute value) eigenvalues are returned.
-----------------------
References:
1. Sign inversion for \[Zeta]: Wybourne, Spectroscopic Properties of Rare Earths. 
2. Sign inversion for {T2, T3, T4, T6, T7, T8}: Hansen and Judd, Matrix Elements of Scalar Three Electron Operators for the Atomic f Shell.";

Options[SolveStates] = {"Return Symbolic Matrix" -> False,
                        "maxEigenvalues" -> "All"};
SolveStates[nf_, IiN_, params0_, OptionsPattern[]]:= Module[
  {n, ii, jj, JMvals},
  maxEigen = OptionValue["maxEigenvalues"];
  (*#####################################*)
  (*hole-particle equivalence enforcement*)
  n = nf;
  If[nf>7, 
  (
    n = 14-nf; (*hole-particle equivalence*)
    params = HoleElectronConjugation[params0];
  ),
    params = params0;
  ];
  (*hole-particle equivalence enforcement*)
  (*#####################################*)
  (*Load symbolic expressions for energy sub-matrices.*)
  Get[EnergyMatrixFileName[n, IiN]];
  (*Patch together the entire matrix representation in block-diagonal form.*)
  EnergyMatrix = ConstantArray[0, {Length[AllowedJ[n]], Length[AllowedJ[n]]}];
  Do[EnergyMatrix[[jj, ii]] = EnergyMatrixTable[{n, AllowedJ[n][[ii]], AllowedJ[n][[jj]], IiN, IiN}];,
    {ii, 1, Length[AllowedJ[n]]},
    {jj, 1, Length[AllowedJ[n]]}
    ];
  EnergyMatrix = ArrayFlatten[EnergyMatrix];
  symbolicMatrix = EnergyMatrix;
  EnergyMatrix = ReplaceInSparseArray[EnergyMatrix, params];
  problemSize = Dimensions[EnergyMatrix][[1]];
  If[maxEigen!="All",
  (If[Abs[maxEigen]>problemSize,
    maxEigen="All"]
    )
    ];
  PrintTemporary["The energy matrix has dimensions:", Dimensions[EnergyMatrix]];
  (*Solve for eigenvalues and eigenvectors.*)
  {EigenvalueJM, EigenvectorJM} = If[maxEigen=="All",
                  Eigensystem[EnergyMatrix],
                  Eigensystem[EnergyMatrix, maxEigen]
                  ];
  EigenvalueJM = Re[EigenvalueJM];
  (*There might be a very small imaginary part.*)
  (*Parse the results for the eigenvectors in terms of the ordered basis being used.*)
  basis = {};
  Do[basis = Join[basis, EnergyStatesTable[{n, AllowedJ[n][[nn]], IiN}]],
    {nn, 1, Length[AllowedJ[n]]}
    ];
  levels = {};
  Do[levels = Join[levels, {{EigenvalueJM[[nn]], EigenvectorJM[[nn]]}}];,
    {nn, 1, Length[EigenvalueJM]}
    ];
  If[OptionValue["Return Symbolic Matrix"],
  Return[{levels, basis, symbolicMatrix}]
  ];
  Return[{levels, basis}];
];

Options[FitToHam] = {
   "LogSolution" -> False,
   "TweakSigns" -> True,
   "ProgressWindow" -> True};
FitToHam::usage = "
FitToHam can be use to fit the parameters of a given model Hamiltonian to a set of observed energy levels. The starting values are taken from Carnall's paper on LaF3, and if the option \"TweakSigns\" is set to True, then the signs of the parameters for the crystal field parameters are randomly switched. If the option \"LogSolution\" is set to True, then each of the found solutions is written to disk. If the option \"ProgressWindow\" is set to True, then a progress window is opened that shows the progress of the average root mean squared error together with a view of the changes in the fitted parameters, in this plot at most 200 values for all the parameters are shown, and all the values are normalized using the first element of the last 200 values.

Parameters
----------
numElectons (int):
ionSymbol (str): the atomic symbol for the corresponding ion.
hamMatrix (List[List] or SparseArray): representing the model Hamiltonian with all its corresponding parameters.
expData (List): energy values for the different levels in increasing order.
simplifier (List): with replacement rules that are applied to modelMatrix. These replacement rules can be used to exclude certain parameters from the fitting, or to impose certain symmetry constraints.
picker (List): a list of booleans that determine wheter the corresponding energy in expData will be used or not for fitting.
numSuperCycles (int): the solver proceeds in stages, this value determines how many times this set of stages will be repeated.
solverIterations (List): a list of integers that determines how many iterations will be run for each of the steps in the solver sequence. Must have as many values as varConcerto does.

Returns
-------
solutions (List): a list of lists where the first element echoes back the starting values for the fit parameters and the second element equals a list where the first element is the minimum rms error and the second is an Association corresponding to the fitted parameters. 
";
FitToHam[numElectrons_,
  ionSymbol_,
  hamMatrix_,
  expData_,
  simplifier_,
  picker_,
  numSuperCycles_,
  numSols_,
  hostName_:"",
  solverIterations_ : {100, 200, 100, 100, 400, 10},
  OptionsPattern[]] := (
  Off[FindMinimum::lstol];
  Off[FindMinimum::cvmit];
  numFittedLevels = Total[If[#, 1, 0] & /@ picker];
  startParams = Association[Carnall["data"][ionSymbol]];
  (*For improved efficiency the given matrix needs to be converted to a CompiledFunction.*)
  fitMatrix = ReplaceInSparseArray[hamMatrix, simplifier];
  reps = 
   Transpose[{{B02, B04, B06, B22, B24, B26, B44, B46, B66, E1, E2, 
      E3, M0, P2, T2, T3, T4, T6, T7, 
      T8, \[Alpha], \[Beta], \[Gamma], \[Zeta]},
     {B02v, B04v, B06v, B22v, B24v, B26v, B44v, B46v, B66v, E1v, E2v, 
      E3v, M0v, 
      P2v, -T2v, -T3v, -T4v, -T6v, -T7v, -T8v, \[Alpha]v, \[Beta]v, \[Gamma]v, -\[Zeta]v}}];
  reps = (#[[1]] -> #[[2]]) & /@ reps;
  funMatrix = Normal[fitMatrix];
  fMatrix = 
   Compile[{B02v, B04v, B06v, B22v, B24v, B26v, B44v, B46v, B66v, E1v,
      E2v, E3v, M0v, P2v, T2v, T3v, T4v, T6v, T7v, 
     T8v, \[Alpha]v, \[Beta]v, \[Gamma]v, \[Zeta]v}, 
    Evaluate[funMatrix /. reps]];
  
  {numIter0, numIter1, numIter2, numIter3, numIter4, numIter5} = solverIterations;
  (*This association is used to define the different stages in the solver.*)
  varConcerto = {
    0 -> {numIter0, {E1, E2, 
       E3, \[Alpha], \[Beta], \[Gamma], \[Zeta]}},
    1 -> {numIter1, modelVars},
    2 -> {numIter2, {B02, B04, B06, B22, B24, B26, B44, B46, B66}}, 
    3 -> {numIter3, {M0, P2, T2, T3, T4, T6, T7, T8}},
    4 -> {numIter4, modelVars},
    5 -> {numIter5, \[Epsilon]}};
  varConcerto = Association[varConcerto];
  
  (* This template *)
  solverTemplate = StringTemplate["
ansatz=`chosenVarListR`;
vars=Transpose[{`chosenVarList`,ansatz}];
partialSol=FindMinimum[
SumOSquares[`B02`,`B04`,`B06`,`B22`,`B24`,`B26`,`B44`,`B46`,`B66`,`E1`\
,`E2`,`E3`,`M0`,`P2`,`T2`,`T3`,`T4`,`T6`,`T7`,`T8`,`\[Alpha]`,`\[Beta]\
`,`\[Gamma]`,`\[Zeta]`, `\[Epsilon]`],
vars,
MaxIterations->`maxIters`,
Method->\"QuasiNewton\",
StepMonitor:>(
rmsHistory=AddToList[rmsHistory,SumOSquares[`B02`,`B04`,`B06`,`B22`,`\
B24`,`B26`,`B44`,`B46`,`B66`,`E1`,`E2`,`E3`,`M0`,`P2`,`T2`,`T3`,`T4`,`\
T6`,`T7`,`T8`,`\[Alpha]`,`\[Beta]`,`\[Gamma]`,`\[Zeta]`,`\[Epsilon]`],maxHistory];
paramSols=AddToList[paramSols,{`B02`,`B04`,`B06`,`B22`,`B24`,`B26`,`\
B44`,`B46`,`B66`,`E1`,`E2`,`E3`,`M0`,`P2`,`T2`,`T3`,`T4`,`T6`,`T7`,`\
T8`,`\[Alpha]`,`\[Beta]`,`\[Gamma]`,`\[Zeta]`,`\[Epsilon]`},maxHistory];
)
];
`chosenVarListR`= Last/@partialSol[[2]];
"];
  
  AddToList[list_, element_, maxSize_] := Module[{
     tempList = Append[list, element]},
    If[Length[tempList] > maxSize,
     Drop[tempList,
      Length[tempList] - maxSize],
     tempList]
    ];

  progressTemplate = StringTemplate["\[Sigma]=`rms` | `elapsedTime`"];
  ProgressNotebook[] := (nb = CreateDocument[(
       Dynamic[
        GraphicsColumn[{
          ListPlot[rmsHistory,
           PlotMarkers -> "OpenMarkers",
           Frame -> True,
           FrameLabel -> {"Iteration", "RMS"},
           ImageSize -> 800,
           AspectRatio -> 1/3,
           FrameStyle -> Directive[Thick, 15], 
           PlotLabel -> If[Length[rmsHistory] != 0, progressTemplate[<|"rms" -> rmsHistory[[-1]], "elapsedTime" -> (Now - startTime)|>], ""]
            ], 
          ListPlot[(#/#[[1]]) & /@ Transpose[paramSols],
           Joined -> True,
           PlotRange -> {All, {-5, 5}},
           Frame -> True,
           ImageSize -> 800,
           AspectRatio -> 1,
           FrameStyle -> Directive[Thick, 15],
           FrameLabel -> {"Iteration", "Params"}]}], 
        TrackedSymbols :> {rmsHistory, paramSols}
        ]),
      WindowSize -> {590, 750},
      WindowSelected -> True,
      WindowTitle -> "Solver Progress"];
    Return[nb];
    );
  
  LogSol[init_, final_, solHistory_, prefix_:""] := (
    fname = prefix <> "-sols-" <> CreateUUID[] <> ".m";
    optimum = final[[1]];
    Print["Saving solution to: ", fname];
    exporter = 
     Association[{"start" -> init, 
       "bestRMS" -> optimum,
       "solHistory" -> solHistory,
       "prefix" -> prefix,
       "end" -> final[[2]]}];
    Export[fname, exporter]);
  
  solutions = {};
  Do[
   (
    ClearAll[SumOSquares];
    SumOSquares[
      B02v_?NumericQ, B04v_?NumericQ, B06v_?NumericQ,
      B22v_?NumericQ, B24v_?NumericQ, B26v_?NumericQ,
      B44v_?NumericQ, B46v_?NumericQ, B66v_?NumericQ,
      E1v_?NumericQ, E2v_?NumericQ, E3v_?NumericQ,
      M0v_?NumericQ, P2v_?NumericQ, T2v_?NumericQ,
      T3v_?NumericQ, T4v_?NumericQ, T6v_?NumericQ,
      T7v_?NumericQ, T8v_?NumericQ, \[Alpha]v_?NumericQ,
      \[Beta]v_?NumericQ, \[Gamma]v_?NumericQ,
      \[Zeta]v_?NumericQ, \[Epsilon]v_?NumericQ] := (
      mat = 
       fMatrix[B02v, B04v, B06v, B22v, B24v, B26v, B44v, B46v, B66v, 
        E1v, E2v, E3v, M0v, P2v, T2v, T3v, T4v, T6v, T7v, 
        T8v, \[Alpha]v, \[Beta]v, \[Gamma]v, \[Zeta]v];
      eigenvals = Sort[Eigenvalues[mat, Method -> "Banded"]];
      eigenvals = Chop[eigenvals - Min[eigenvals]];
      eigenvals = eigenvals[[;; ;; 2]];
      eigenvals = eigenvals[[;; 134]];
      eigenvals = Pick[eigenvals, picker];
      rms = Sqrt[Total[(eigenvals - expData + \[Epsilon]v)^2]/numFittedLevels];
      Return[rms];
      );
    
    numCycles = 4;
    modelVars = {B02, B04, B06, B22, B24, B26, B44, B46, B66, E1, E2, 
      E3, M0, P2, T2, T3, T4, T6, T7, 
      T8, \[Alpha], \[Beta], \[Gamma], \[Zeta], \[Epsilon]};
    maxHistory = 200;
    rmsHistory = {};
    paramSols = {};
    startTime = Now;
    If[OptionValue["ProgressWindow"],
     nb = ProgressNotebook[];
     ];
    maxChosenVars = Length[modelVars];
    minChosenVars = Ceiling[maxChosenVars/2];
    
    ansatz = (# -> startParams[#]) & /@ {B02, B04, B06, B22, B24, B26,
        B44, B46, B66, E1, E2, E3, M0, P2, T2, T3, T4, T6, T7, 
       T8, \[Alpha], \[Beta], \[Gamma], \[Zeta], F2, F4, F6};
    ansatz = Association[ansatz];
    
    (*Convert from Fs to Es*)
    ansatz[E0] = 0;
    ansatz[E1] = 
     14/405*ansatz[F2] + 7/297*ansatz[F4] + 350/11583*ansatz[F6];
    ansatz[E2] = -1/2025 ansatz[F2] - 1/3267 ansatz[F4] + 
      175/1656369*ansatz[F6];
    ansatz[E3] = 
     1/135 ansatz[F2] + 2/1089*ansatz[F4] - 175/42471*ansatz[F6];
    ansatz = 
     ansatz[#] & /@ {B02, B04, B06, B22, B24, B26, B44, B46, B66, E1, 
       E2, E3, M0, P2, T2, T3, T4, T6, T7, 
       T8, \[Alpha], \[Beta], \[Gamma], \[Zeta]};
    
    If[OptionValue["TweakSigns"],
     Do[
      ansatz[[i]] = 
       RandomChoice[{1, -1}]*ansatz[[i]], {i, {1, 2, 3, 4, 5, 6, 7, 8,
         9, 13, 14, 15, 16, 17, 18, 19, 20}}
      ]
     ];
    
    {B02r, B04r, B06r, B22r, B24r, B26r, B44r, B46r, B66r, E1r, E2r, 
      E3r, M0r, P2r, T2r, T3r, T4r, T6r, T7r, 
      T8r, \[Alpha]r, \[Beta]r, \[Gamma]r, \[Zeta]r} = ansatz[[;; 24]];
    \[Epsilon]r = 0.001;
    ansatz = Append[ansatz, \[Epsilon]r];
    startingPoint = 
     Association[(#[[1]] -> #[[2]]) & /@ 
       Transpose[{modelVars, ansatz}]];
    
    loopMessage = StringTemplate["Cycle `numCycle`/`allCycles`"];
    
    Do[
     Do[
      (
       {maxIters, chosenVars} = varConcerto[cycle];
       frozenVars = Complement[modelVars, chosenVars];
       frozenReps = (ToString[#] -> ToString[#] <> "r") & /@ 
         frozenVars;
       varyingReps = (ToString[#] -> ToString[#]) & /@ chosenVars;
       varyingRepsAnsatzBridge = (ToString[#] <> "r") & /@ 
         chosenVars;
       stringReps = Association[Join[frozenReps, varyingReps]];
       stringReps["chosenVarList"] = ToString[chosenVars];
       stringReps["chosenVarListR"] = 
        ToString[varyingRepsAnsatzBridge];
       stringReps["maxIters"] = 
        ToString[solverIterations[[cycle + 1]]];
       solverString = solverTemplate[stringReps];
       pt = PrintTemporary[
        loopMessage[<|"numCycle" -> (cycle+1), 
          "allCycles" -> (numCycles+1)|>]," : ","Using: ", chosenVars];
       ToExpression[solverString];
       NotebookDelete[pt];
       ),
      {cycle, 0, numCycles}
      ],
     {superCycle, 1, numSuperCycles}
     ];
    If[OptionValue["ProgressWindow"],
     NotebookClose[nb];
     ];
    If[OptionValue["LogSolution"],
     LogSol[startingPoint, partialSol, paramSols, ionSymbol<>"-"<>hostName];
     ];
    AppendTo[solutions, {startingPoint, partialSol}];
    ),
   {bigReps, 1, numSols}
   ];
  Return[solutions]
  )

(* ###############################              SOLVERS            ############################## *)
(* ############################################################################################## *)

(* ############################################################################################## *)
(* ############################        EIGEN-VECTOR ANALYSIS          ########################### *)

PrettySaunders::usage = "PrettySaunders[SL, J] produces a human-redeable symbol for the given Russel-Saunders term."
PrettySaundersSLJ[{{SL_, J_}, MJ_}] := (
  If[StringQ[SL],
   {S, L} = findSL[SL],
   {S, L} = SL
   ];
  Return[RowBox[
    {AdjustmentBox[Style[2*S + 1, Smaller], BoxBaselineShift -> -1, 
      BoxMargins -> 0], AdjustmentBox[printL[L], BoxMargins -> -0.2],
     AdjustmentBox[
      Style[InputForm[J], Small, FontTracking -> "Narrow"], 
      BoxBaselineShift -> 1, BoxMargins -> {{0.7, 0}, {0.4, 0.4}}]
     }] // DisplayForm])

PrettySaundersSLJmJ[{{SL_, J_}, mJ_}] := (If[
   StringQ[SL], 
   ({S, L} = findSL[SL];
    L = StringTake[SL, {2, -1}];
    ), 
   {S, L} = SL];
  Return[
   RowBox[{AdjustmentBox[Style[2*S + 1, Smaller], 
       BoxBaselineShift -> -1, BoxMargins -> 0], 
      AdjustmentBox[printL[L], BoxMargins -> -0.2], 
      AdjustmentBox[
       Style[{InputForm[J], mJ}, Small, FontTracking -> "Narrow"], 
       BoxBaselineShift -> 1, 
       BoxMargins -> {{0.7, 0}, {0.4, 0.4}}]}] // DisplayForm])

BasisVecInRusselSaunders::usage = "BasisVecInRusselSaunders[basisVec] takes a basis vector in the format {{{LSstring, Jval}, mJval}, nucSpin} and returns a human-readable symbol for the corresponding Russel-Saunders term. The nuclear spin is ignored."
BasisVecInRusselSaunders[basisVec_] := (
  {{{LSstring, Jval}, mJval}, nucSpin} = basisVec;
  Ket[PrettySaunders[LSstring, Jval], mJval]
  )

LSJMJTemplate = 
  StringTemplate[
   "\!\(\*TemplateBox[{\nRowBox[{\"`LS`\", \",\", \nRowBox[{\"J\", \
\"=\", \"`J`\"}], \",\", \nRowBox[{\"mJ\", \"=\", \"`mJ`\"}]}]},\n\
\"Ket\"]\)"];
BasisVecInLSJMJ::usage = "BasisVecInLSJMJ[basisVec] takes a basis vector in the format {{{LSstring, Jval}, mJval}, nucSpin} and returns a human-readable symbol for the corresponding LSJMJ term in the form |LS, J=..., mJ=...> The nuclear spin is ignored."
BasisVecInLSJMJ[basisVec_] := (
   {{{LSstring, Jval}, mJval}, nucSpin} = basisVec;
   LSJMJTemplate[<|
     "LS" -> LSstring,
     "J" -> ToString[Jval, InputForm], 
     "mJ" -> ToString[mJval, InputForm]|>]
   );

ParseStates::usage = "ParseStates[states, basis] takes a list of eigenstates in terms of their coefficients in the given basis and returns a list of the same states in terms of their energy, LSJMJ symbol, J, mJ, S, L, LSJ symbol, and LS symbol. The LS symbol returned corresponds to the term with the largest coefficient in the given basis.";
ParseStates[states_, basis_, OptionsPattern[]] := (
  parsedStates = Table[(
     {energy, eigenVec} = state;
     maxTermIndex = Ordering[Abs[eigenVec]][[-1]];
     {{{LSstring, Jval}, mJval}, nucSpin} = basis[[maxTermIndex]];
     LSJsymbol = Subscript[LSstring, {Jval, mJval}];
     LSJMJsymbol = LSstring <> ToString[Jval, InputForm];
     {S, L} = findSL[LSstring];
     {energy, LSstring, Jval, mJval, S, L, LSJsymbol, LSJMJsymbol}
     ),
    {state, fstates}];
  Return[parsedStates]
  )

ParseStatesByNumBasisVecs::usage = "ParseStatesByNumBasisVecs[states, basis, numBasisVecs] takes a list of eigenstates in terms of their coefficients in the given basis and returns a list of the same states in terms of their energy and the coefficients of the numBasisVecs most significant basis vectors.";
ParseStatesByNumBasisVecs[states_, basis_, numBasisVecs_, roundTo_ : 0.01] := (
  parsedStates = Table[(
     {energy, eigenVec} = state;
     energy = Chop[energy];
     probs = Round[Abs[eigenVec^2], roundTo];
     amplitudes = Round[eigenVec, roundTo];
     ordering = Ordering[probs];
     chosenIndices = ordering[[-numBasisVecs ;;]];
     majorComponents = basis[[chosenIndices]];
     majorProbabilities = amplitudes[[chosenIndices]];
     majorComponents = BasisVecInLSJMJ /@ majorComponents;
     majorRep = majorProbabilities . majorComponents;
     {energy, majorRep}
     ),
    {state, fstates}];
  Return[parsedStates]
  )

FindThresholdPosition::usage = "FindThresholdPosition[list, threshold] returns the position of the first element in list that is greater than threshold. If no such element exists, it returns the length of list. The elements of the given list must be in ascending order.";
FindThresholdPosition[list_, threshold_] := 
 Module[{position}, 
  position = Position[list, _?(# > threshold &), 1, 1];
  thrPos = If[Length[position] > 0,
    position[[1, 1]],
    Length[list]];
  If[thrPos == 0, Return[1], Return[thrPos+1]]
  ]

ParseStatesByProbabilitySum::usage = "ParseStatesByProbabilitySum[states, basis, probSum] takes a list of eigenstates in terms of their coefficients in the given basis and returns a list of the same states in terms of their energy and the coefficients of the basis vectors that sum to at least probSum.";
ParseStatesByProbabilitySum[states_, basis_, probSum_, roundTo_ : 0.01] :=(
    parsedStates = Table[({energy, eigenVec} = state;
      energy = Chop[energy];
      amplitudes = Round[eigenVec, roundTo];
      probs = Round[Abs[eigenVec^2], roundTo];
      ordering = Reverse[Ordering[probs]];
      orderedProbs = probs[[ordering]];
      accProb = Accumulate[orderedProbs];
      thresholdIndex = Min[Length[accProb], FindThresholdPosition[accProb, probSum]];
      chosenIndices = ordering[[;; thresholdIndex]];
      majorComponents = basis[[chosenIndices]];
      majorProbabilities = amplitudes[[chosenIndices]];
      notnullProbs = 
      Flatten[Position[majorProbabilities, _?(Abs[#] != 0 &)]];
      majorComponents = PrettySaundersSLJmJ[#[[1]]] & /@ majorComponents;
      majorProbabilities = majorProbabilities[[notnullProbs]];
      majorComponents = Ket /@ majorComponents[[notnullProbs]];
      majorRep = majorProbabilities . majorComponents;
    {energy, majorRep}), {state, states}];
 Return[parsedStates])

(* ############################        EIGEN-VECTOR ANALYSIS          ########################### *)
(* ############################################################################################## *)

(* ############################################################################################## *)
(* ###############################               MISC              ############################## *)

SymbToNum::usage = "SymbToNum[expr, numAssociation] takes an expression expr and returns what results after making the replacements defined in the given replacementAssociation. If replacementAssociation doesn't define values for expected keys, they are taken to be zero.";
SymbToNum[expr_, replacementAssociation_]:= (
  includedKeys = Keys[replacementAssociation];
  (*If a key is not defined, make its value zero.*)
  fullAssociation = Table[(
    If[MemberQ[includedKeys, key],
      ToExpression[key]->replacementAssociation[key],
      ToExpression[key]->0
    ]
  ),
  {key, paramSymbols}];
  Return[expr/.fullAssociation];
)

SimpleConjugate::usage = "SimpleConjugate[expr] takes an expression and applies a simplified version of the conjugate in that all it does is that it replaces the imaginary unit I with -I. It assumes that every other symbol is real so that it remains the same under complex conjugation. Among other expressions it is valid for any rational or polynomial expression with complex coefficients and real variables.";
SimpleConjugate[expr_] := expr /. Complex[a_, b_] :> a - I b;

ExportMZip::usage = 
  "ExportMZip[filename, object] exports the given object and compresses it. It first exports in .m format, it then compresses that file in to a zip file, and finally the .m file is deleted. The filename must be a full path and end with .m. This probably won't work on a PC.";
ExportMZip[filename_, object_] := (
   zipTemplate = StringTemplate["cd \"`sourceDir`\"; zip \"`dest`\" \"`source`\""];
   delTemplate = StringTemplate["rm \"`rmFname`\""];
   Export[filename, object];
   zipFilename = StringReplace[filename, ".m" -> ".zip"];
   splitName = FileNameSplit[zipFilename];
   zipFilename = splitName[[-1]];
   sourceDir = FileNameJoin[splitName[[1 ;; -2]]];
   zipCmd = 
    zipTemplate[<|"sourceDir" -> sourceDir, "dest" -> zipFilename, 
      "source" -> FileNameSplit[filename][[-1]]|>];
   delCmd = delTemplate[<|"rmFname" -> filename|>];
   Run[zipCmd];
   Run[delCmd];
   );

ImportMZip::usage = 
  "ImportMZip[filename] decompresses the provided filename, and imports the enclosed .m file that it is assumed to contain. After being imported the uncompressed file is deleted from disk. The provided filename must be a full path and end with .zip. This probably won't work on a PC.";
ImportMZip[filename_] := (
  splitName     = FileNameSplit[filename];
  sourceDir     = FileNameJoin[splitName[[1 ;; -2]]];
  delTemplate   = StringTemplate["rm \"`rmFname`\""];
  unzipTemplate = StringTemplate["cd \"`sourceDir`\"; unzip \"`source`\""];
  unzipCmd      = unzipTemplate[<|"sourceDir" -> sourceDir, 
                                  "source" -> FileNameSplit[filename][[-1]]|>];
  Run[unzipCmd];
  mFilename = StringReplace[filename, ".zip" -> ".m"];
  imported  = Import[mFilename];
  delCmd    = delTemplate[<|"rmFname" -> mFilename|>];
  Run[delCmd];
  Return[imported];
  );

EnergyLevelDiagram::usage = "EnergyLevelDiagram[states] takes states and produces a visualization of its energy spectrum.
The resultant visualization can be navigated by clicking and dragging to zoom in on a region, or by clicking and dragging horizontally while pressing Ctrl. Double-click to reset the view.";
Options[EnergyLevelDiagram] = {
  "Title"->"", 
  "ImageSize"->1000, 
  "AspectRatio" -> 1/8, 
  "Background"->"Automatic",
  "Epilog"->{}
  };
EnergyLevelDiagram[states_, OptionsPattern[]]:= (
  energies = First/@states;
  epi = OptionValue["Epilog"];
  ExploreGraphics@ListPlot[Tooltip[{{#, 0}, {#, 1}}, {Quantity[#/8065.54429, "eV"], Quantity[#, 1/"Centimeters"]}] &/@ energies,
    Joined       -> True,
    PlotStyle    -> Black,
    AspectRatio  -> OptionValue["AspectRatio"],
    ImageSize    -> OptionValue["ImageSize"],
    Frame        -> True,
    PlotRange    -> {All, {0, 1}},
    FrameTicks   -> {{None, None}, {Automatic, Automatic}},
    FrameStyle   -> Directive[15, Dashed, Thin],
    PlotLabel    -> Style[OptionValue["Title"], 15, Bold],
    Background   -> OptionValue["Background"],
    FrameLabel   -> {"\!\(\*FractionBox[\(E\), SuperscriptBox[\(cm\), \(-1\)]]\)"},
    Epilog       -> epi]
)

ExploreGraphics::usage = 
  "Pass a Graphics object to explore it. Zoom by clicking and dragging a rectangle. Pan by clicking and dragging while pressing Ctrl. Click twice to reset view. 
  Based on ZeitPolizei @ https://mathematica.stackexchange.com/questions/7142/how-to-manipulate-2d-plots";

OptAxesRedraw::usage = 
  "Option for ExploreGraphics to specify redrawing of axes. Default False.";
Options[ExploreGraphics] = {OptAxesRedraw -> False};

ExploreGraphics[graph_Graphics, opts : OptionsPattern[]] := With[
   {gr  = First[graph],
    opt = DeleteCases[Options[graph], 
          PlotRange -> PlotRange | AspectRatio | AxesOrigin -> _],
    plr = PlotRange /. AbsoluteOptions[graph, PlotRange],
    ar  = AspectRatio /. AbsoluteOptions[graph, AspectRatio],
    ao  = AbsoluteOptions[AxesOrigin],
    rectangle = {Dashing[Small], 
       Line[{#1, 
            {First[#2], Last[#1]},
            #2, 
            {First[#1], Last[#2]},
            #1}]} &,
    optAxesRedraw = OptionValue[OptAxesRedraw]},
   DynamicModule[
      {dragging=False, first, second, rx1, rx2, ry1, ry2,
      range = plr}, 
      {{rx1, rx2}, {ry1, ry2}} = plr;
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
         If[CurrentValue@"MouseClickCount"==2,
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

ReplaceInSparseArray::usage = "ReplaceInSparseArray[sparseArray, rules] takes a sparse array that may contain symbolic quantities and returns a sparse array in which the given replacement rules have been used.";
ReplaceInSparseArray[s_SparseArray, rule_] := (With[{
    elem = s["NonzeroValues"]/.rule,
    def  = s["Background"]/.rule
    },
    srep = SparseArray[Automatic, s["Dimensions"], def, {1, {s["RowPointers"], s["ColumnIndices"]}, elem}];
  ];
  Return[srep];
  );

Options[ParseTeXLikeSymbol] = {"Form" -> "List"};
ParseTeXLikeSymbol::usage = 
  "ParseTeXLikeSymbol[string] parses a string for a symbol given in \
latex notation and returns a corresponding mathematica symbol. The \
string may have expressions for several symbols, they need to be \
separated by single spaces. In addition the _ and ^ symbols used in \
LaTeX notation need to have arguments that are enclosed in \
parenthesis, for example \"x_2\" is invalid, instead \"x_{2}\" should \
have been given.";
ParseTeXLikeSymbol[bigString_, OptionsPattern[]] := (
   form = OptionValue["Form"];
   (*parse greek*)
   symbols = Table[(
      str = StringReplace[string, {"\\alpha" -> "\[Alpha]",
         "\\beta" -> "\[Beta]",
         "\\gamma" -> "\[Gamma]",
         "\\psi" -> "\[Psi]"}];
      symbol = Which[
        StringContainsQ[str, "_"] && Not[StringContainsQ[str, "^"]],
        (
         (*yes sub no sup*)
         mainSymbol = StringSplit[str, "_"][[1]];
         mainSymbol = ToExpression[mainSymbol];
         
         subPart = 
          StringCases[str, 
            RegularExpression@"\\{(.*?)\\}" -> "$1"][[1]];
         Subscript[mainSymbol, subPart]
         ),
        Not[StringContainsQ[str, "_"]] && StringContainsQ[str, "^"],
        (
         (*no sub yes sup*)
         mainSymbol = StringSplit[str, "^"][[1]];
         mainSymbol = ToExpression[mainSymbol];
         
         supPart = 
          StringCases[str, 
            RegularExpression@"\\{(.*?)\\}" -> "$1"][[1]];
         Superscript[mainSymbol, supPart]),
        StringContainsQ[str, "_"] && StringContainsQ[str, "^"],
        (
         (*yes sub yes sup*)
         mainSymbol = StringSplit[str, "_"][[1]];
         mainSymbol = ToExpression[mainSymbol];
         {subPart, supPart} = 
          StringCases[str, RegularExpression@"\\{(.*?)\\}" -> "$1"];
         Subsuperscript[mainSymbol, subPart, supPart]
         ),
        True,
        ((*no sup or sub*)
         str)
        ];
      symbol
      ),
     {string, StringSplit[bigString, " "]}];
   Which[
    form == "Row",
    Return[Row[symbols]],
    form == "List",
    Return[symbols]
    ]
   );
Options[LabeledGrid]={
    ItemSize->Automatic,
    Alignment->Center,
    Frame->All,
    "Separator"->",",
    "Pivot"->""
};
LabeledGrid::usage="LabeledGrid[data, rowHeaders, columnHeaders] provides a grid of given data interpreted as a matrix of values whose rows are labeled by rowHeaders and whose columns are labeled by columnHeaders. When hovering with the mouse over the grid elements, the row and column labels are displayed with the given separator between them.";
LabeledGrid[data_,rowHeaders_,columnHeaders_,OptionsPattern[]]:=Module[
    {gridList=data,rowHeads=rowHeaders,colHeads=columnHeaders},
    (
    separator=OptionValue["Separator"];
    pivot=OptionValue["Pivot"];
    gridList=Table[
            Tooltip[
              data[[rowIdx,colIdx]],
              DisplayForm[
                RowBox[{rowHeads[[rowIdx]],
                        separator,
                        colHeads[[colIdx]]}
                      ]
                      ]
                    ],
        {rowIdx,Dimensions[data][[1]]},
        {colIdx,Dimensions[data][[2]]}];
    gridList=Transpose[Prepend[gridList,colHeads]];
    rowHeads=Prepend[rowHeads,pivot];
    gridList=Prepend[gridList,rowHeads]//Transpose;
    Grid[gridList,
        Frame->OptionValue[Frame],
        Alignment->OptionValue[Alignment],
        Frame->OptionValue[Frame],
        ItemSize->OptionValue[ItemSize]
        ]
)
]

Options[HamiltonianForm]={
    "Separator"->"",
    "Pivot"->""
}
HamiltonianForm::usage="HamiltonianForm[hamMatrix, basisLabels] takes the matrix representation of a hamiltonian together with a set of symbols representing the ordered basis in which the operator is represented. With this it creates a displayed form that has adequately labeled row and columns together with informative values when hovering over the matrix elements using the mouse cursor.";
HamiltonianForm[hamMatrix_, basisLabels_, OptionsPattern[]]:=(
    braLabels=DisplayForm[RowBox[{"\[LeftAngleBracket]",#,"\[RightBracketingBar]"}]]&/@basisLabels;
    ketLabels=DisplayForm[RowBox[{"\[LeftBracketingBar]",#,"\[RightAngleBracket]"}]]&/@basisLabels;
    LabeledGrid[hamMatrix,braLabels,ketLabels,"Separator"->OptionValue["Separator"],"Pivot"->OptionValue["Pivot"]]
    )

(* ###############################               MISC              ############################## *)
(* ############################################################################################## *)

(* ############################################################################################## *)
(* ####################################### Ion-Solvers ########################################## *)


ErbiumSolver::usage = "ErbiumSolver[numSuperCycles, numSols, solverIterations] runs the solver for the Erbium problem. The solver is run for numSuperCycles super cycles, each of which consists of numSols solutions.

Inspect the code to modify a number of other parameters that can also be tweaked.

The solver returns a list of numSols solutions, each of which is a list of the form {startingPoint, solution}. Where startingPoint is the starting point for the solution, and solution is a list whose first element is the minimum value for the RMS errror and whose second element is an list of replacement rules for the Hamiltonian parameters.

This solver is rather ad-hoc, but it works well enough for the Erbium problem. But is left here since it will surely be helpful in other cases.

The matrices that are needed for the solver are imported from the file ./data/Erbium-Solver-Toolbits.m. This file contains the ansatz data from Carnall for Er3+ in LaF3. It also contains a template used to call FindMinimum with different chose varialbes. And most importantly, it contains a function that allows for efficient creation of the matrix that is diagonalized in the solver.

The starting values for the solution parameters are taken from Carnall's data for LaF3:Er3+. The different solutions returned by this function differ in how the signs of {B02, B04, B06, B22, B24, B26, B44, B46, B66, M0, P2, T2, T3, T4, T6, T7, T8} are randomly changed from the Carnall data on lanthanum fluoride.

The solver uses the following steps: 1) First only E1, E2, E3, α, β, γ, ζ are allowed to vary. 2) Then all parameters are allowed to vary. 3) Then only the crystal field parameters are allowed to vary. 4) Then only {M0,P2,T2,T3,T4,T6,T7,T8} are allowed to vary. 5) Finally all parameters are allowed to vary again. This sequency is repeated numSuperCycles times.

The solver uses the QuasiNewton method option in FindMinimum.

The solverIterations argument is a list of the number of iterations for each of the 5 solver steps. The default is solverIterations = {100, 200, 100, 100, 400}. 

The experimental data used here is from Gruber et al. - 1993 - Energy levels and correlation crystal-field effects in Er3+ doped garnets.

The solver creates a progress notebook that can be used to monitor the progress of the solver. The progress notebook can be closed at any time, the solver will continue to run in the background.

The following is the string used for the solverTemplate template:

  ansatz=`chosenVarListR`;
  vars=Transpose[{
  `chosenVarList`,
  ansatz}];
  partialSol=FindMinimum[
  SumOSquares[`B02`,`B04`,`B06`,`B22`,`B24`,`B26`,`B44`,`B46`,`B66`,
  `E1`,`E2`,`E3`,`M0`,`P2`,`T2`,`T3`,`T4`,`T6`,`T7`,`T8`,`\[Alpha]`,`\
  \[Beta]`,`\[Gamma]`,`\[Zeta]`],
  vars,
  MaxIterations->maxIters,
  Method->\" QuasiNewton \",
  StepMonitor:>(
  rmsHistory=AddToList[rmsHistory,SumOSquares[`B02`,`B04`,`B06`,`B22`,`\
  B24`,`B26`,`B44`,`B46`,`B66`,`E1`,`E2`,`E3`,`M0`,`P2`,`T2`,`T3`,`T4`,`\
  T6`,`T7`,`T8`,`\[Alpha]`,`\[Beta]`,`\[Gamma]`,`\[Zeta]`],maxHistory];
  paramSols=AddToList[paramSols,{`B02`,`B04`,`B06`,`B22`,`B24`,`B26`,`\
  B44`,`B46`,`B66`,`E1`,`E2`,`E3`,`M0`,`P2`,`T2`,`T3`,`T4`,`T6`,`T7`,`\
  T8`,`\[Alpha]`,`\[Beta]`,`\[Gamma]`,`\[Zeta]`},maxHistory];
  )
  ];
  `chosenVarListR`= Last/@partialSol[[2]];

";

ErbiumSolver[numSuperCycles_, numSols_, solverIterations_:{100, 200, 100, 100, 400}]:=(

{numIter0, numIter1, numIter2, numIter3, numIter4} = solverIterations;
AddToList[list_, element_, maxSize_] := 
Module[{tempList = Append[list, element]}, 
  If[Length[tempList] > maxSize, 
  Drop[tempList, Length[tempList] - maxSize], tempList]];

ProgressNotebook[] := (
  nb = CreateDocument[
    (
    Dynamic[GraphicsColumn[{
        ListPlot[rmsHistory,
        PlotMarkers -> "OpenMarkers",
        Frame -> True,
        FrameLabel -> {"Iteration", "RMS"},
        ImageSize -> 800,
        AspectRatio -> 1/3,
        FrameStyle -> Directive[Thick, 15],
        
        PlotLabel -> 
          If[Length[rmsHistory] != 0, rmsHistory[[-1]], ""]],
        ListPlot[(#/#[[1]]) & /@ Transpose[paramSols],
        Joined -> True,
        PlotRange -> {All, {-5, 5}},
        Frame -> True,
        ImageSize -> 800,
        AspectRatio -> 1,
        FrameStyle -> Directive[Thick, 15],
        FrameLabel -> {"Iteration", "Params"}]
        }], TrackedSymbols :> {rmsHistory, paramSols}]
    ),
    WindowSize -> {590, 750},
    WindowSelected -> True,
    WindowTitle -> "Solver Progress"];
  Return[nb];
  );
solutions = {};

Do[(
  (
   trunk = Import["./data/Er-Solver-Toolbits.m"];
   carnallAnsatz = trunk["carnallAnsatz"];
   fMatrix = trunk["fMatrix"];
   expData = trunk["expData"];
   picker = trunk["picker"];
   solverTemplate = trunk["solverTemplate"];
   LogSol[init_, final_, prefix_:""] := (
     fname = prefix <> "sols-" <> CreateUUID[] <> ".m";
     optimum = final[[1]];
     Print["Saving to", fname];
     exporter = 
      Association[{"start" -> init, 
                   "bestRMS" -> optimum, 
                   "end" -> final[[2]]}];
     Export[fname, exporter]
     );
   ClearAll[SumOSquares];
   SumOSquares[B02v_?NumericQ, B04v_?NumericQ, B06v_?NumericQ, 
     B22v_?NumericQ, B24v_?NumericQ, B26v_?NumericQ, B44v_?NumericQ, 
     B46v_?NumericQ, B66v_?NumericQ, E1v_?NumericQ, E2v_?NumericQ, 
     E3v_?NumericQ, M0v_?NumericQ, P2v_?NumericQ, T2v_?NumericQ, 
     T3v_?NumericQ, T4v_?NumericQ, T6v_?NumericQ, T7v_?NumericQ, 
     T8v_?NumericQ, \[Alpha]v_?NumericQ, \[Beta]v_?
      NumericQ, \[Gamma]v_?NumericQ, \[Zeta]v_?NumericQ] := (
     mat = 
      fMatrix[B02v, B04v, B06v, B22v, B24v, B26v, B44v, B46v, B66v, 
       E1v, E2v, E3v, M0v, P2v, T2v, T3v, T4v, T6v, T7v, 
       T8v, \[Alpha]v, \[Beta]v, \[Gamma]v, \[Zeta]v];
     eigenvals = Sort[Eigenvalues[mat, Method -> "Banded"]];
     eigenvals = Chop[eigenvals - Min[eigenvals]];
     eigenvals = eigenvals[[;; ;; 2]];
     eigenvals = eigenvals[[;; 134]];
     eigenvals = Pick[eigenvals, picker];
     rms = Sqrt[Total[(eigenvals - expData)^2]/117.];
     Return[rms];
     );
   );
  numCycles = 4;
  modelVars = {B02, B04, B06, B22, B24, B26, B44, B46, B66, E1, E2, 
    E3, M0, P2, T2, T3, T4, T6, T7, 
    T8, \[Alpha], \[Beta], \[Gamma], \[Zeta]};
  maxHistory = 200;
  rmsHistory = {};
  paramSols = {};
  
  nb = ProgressNotebook[];
  maxChosenVars = Length[modelVars];
  minChosenVars = Ceiling[maxChosenVars/2];
  tweakSigns = True;
  
  (*build initial ansatz*)
  ansatz = (# -> carnallAnsatz[#]) & /@ {B02, B04, B06, B22, B24, B26,
      B44, B46, B66, E1, E2, E3, M0, P2, T2, T3, T4, T6, T7, 
     T8, \[Alpha], \[Beta], \[Gamma], \[Zeta], F2, F4, F6};
  ansatz = Association[ansatz];
  ansatz[E0] = 0;
  ansatz[E1] = 
   14*ansatz[F2]/405 + 7*ansatz[F4]/297 + 350*ansatz[F6]/11583;
  ansatz[E2] = -ansatz[F2]/2025 - ansatz[F4]/3267 + 
    175*ansatz[F6]/1656369;
  ansatz[E3] = 
   ansatz[F2]/135 + 2*ansatz[F4]/1089 - 175*ansatz[F6]/42471;
  ansatz = 
   ansatz[#] & /@ {B02, B04, B06, B22, B24, B26, B44, B46, B66, E1, 
     E2, E3, M0, P2, T2, T3, T4, T6, T7, 
     T8, \[Alpha], \[Beta], \[Gamma], \[Zeta]};
  If[tweakSigns,
   Do[ansatz[[i]] = 
     RandomChoice[{1, -1}]*ansatz[[i]], {i, {1, 2, 3, 4, 5, 6, 7, 8, 
      9, 13, 14, 15, 16, 17, 18, 19, 20}}]
   ];
  {B02r, B04r, B06r, B22r, B24r, B26r, B44r, B46r, B66r, E1r, E2r, 
    E3r, M0r, P2r, T2r, T3r, T4r, T6r, T7r, 
    T8r, \[Alpha]r, \[Beta]r, \[Gamma]r, \[Zeta]r} = ansatz[[;; 24]];
  startingPoint = 
   Association[(#[[1]] -> #[[2]]) & /@ Transpose[{modelVars, ansatz}]];

  varConcerto = {
    0 -> {numIter0, {E1, E2, E3, \[Alpha], \[Beta], \[Gamma], \[Zeta]}},
    1 -> {numIter1, modelVars},
    2 -> {numIter2, {B02, B04, B06, B22, B24, B26, B44, B46, B66}},
    3 -> {numIter3, {M0, P2, T2, T3, T4, T6, T7, T8}},
    4 -> {numIter4, modelVars}};
  varConcerto = Association[varConcerto];
  loopMessage = StringTemplate["Cycle `numCycle`/`allCycles`"];
  Off[FindMinimum::lstol];
  Off[FindMinimum::cvmit];
  
  Do[
   Do[(
     {maxIters, chosenVars} = varConcerto[cycle];
     frozenVars = Complement[modelVars, chosenVars];
     frozenReps = (ToString[#] -> ToString[#] <> "r") & /@ 
       frozenVars;
     varyingReps = (ToString[#] -> ToString[#]) & /@ chosenVars;
     varyingRepsAnsatzBridge = (ToString[#] <> "r") & /@ chosenVars;
     stringReps = Association[Join[frozenReps, varyingReps]];
     stringReps["chosenVarList"] = ToString[chosenVars];
     stringReps["chosenVarListR"] = 
      ToString[varyingRepsAnsatzBridge];
     solverString = solverTemplate[stringReps];
     Print[
      loopMessage[<|"numCycle" -> cycle, 
        "allCycles" -> numCycles|>]];
     Print["Using: ", chosenVars];
     ToExpression[solverString];
     ),
    {cycle, 0, numCycles}],
   {superCycle, 1, numSuperCycles}
   ];
  NotebookClose[nb];
  Print[partialSol];
  AppendTo[solutions,{startingPoint, partialSol}];
  ),
 {bigReps, 1, numSols}
 ];
 Return[solutions]
);

(* ####################################### Ion-Solvers ########################################## *)
(* ############################################################################################## *)


(* ############################################################################################## *)
(* ############################################ Data ############################################ *)

PrintTemporary["Loading data from Carnall for LaF3 ..."];
carnallFname = FileNameJoin[{moduleDir, "data", "Carnall.m"}];
Carnall = Import[carnallFname];
Carnall::usage = "Association of data from Carnall et al (1989) with the following keys: {data, annotations, paramSymbols, elementNames, rawData, rawAnnotations, annnotatedData, appendix:Pr:Association, appendix:Pr:Calculated, appendix:Pr:RawTable, appendix:Headings}";

PrintTemporary["Loading table with coefficients of fractional parentage ..."];
CFPfname = FileNameJoin[{moduleDir, "data", "CFPs.m"}]
If[!FileExistsQ[CFPfname],
  (PrintTemporary[">> CFPs.m not found, generating ..."];
    CFP = GenerateCFP[];
  ),
  CFP = Import[CFPfname];
]
CFP::usage = "CFP[n, termSymbol] provides the coefficients of fractional parentage for the given term termSymbol and number of electrons n.";

PrintTemporary["Loading association with coefficients of fractional parentage ..."];
CFPAfname = FileNameJoin[{moduleDir, "data", "CFPAssoc.m"}]
If[!FileExistsQ[CFPAfname],
  (PrintTemporary[">> CFPAssoc.m not found, generating ..."];
    CFPAssoc = GenerateCFPAssoc[True];
  ),
  CFPAssoc = Import[CFPAfname];
]
CFPAssoc::usage = " CFPAssoc is an association where keys are of lists of the form {num_electrons, daugherTerm, parentTerm} and values are the corresponding coefficients of fractional parentage. The terms given in string-spectroscopic notation. If a certain daughter term does not have a parent term, the value is 0.";

PrintTemporary["Loading table of reduced matrix elements for unit tensor operators ..."]
ReducedUkTableFname = FileNameJoin[{moduleDir, "data", "ReducedUkTable.m"}]
If[!FileExistsQ[ReducedUkTableFname],
  (PrintTemporary[">> ReducedUkTable.m not found, generating ..."];
    ReducedUkTable = GenerateReducedUkTable[];
  ),
  ReducedUkTable = Import[ReducedUkTableFname];
]
ReducedUkTable::usage = "ReducedUkTable[{n, l = 3, SL, SpLp, k}] provides reduced matrix elements of the spherical tensor operator Uk. See Cowan (1981) section 11-9 \"Unit Tensor Operators\".";

PrintTemporary["Loading table of matrix elements for the electrostatic interaction ..."];
ElectrostaticMatrixTablefname = FileNameJoin[{moduleDir, "data", "ElectrostaticMatrixTable.m"}]
If[!FileExistsQ[ElectrostaticMatrixTablefname],
  (PrintTemporary[">> ElectrostaticMatrixTable.m not found, generating ..."];
    ElectrostaticMatrixTable = GenerateElectrostaticMatrixTable[];
  ),
  ElectrostaticMatrixTable = Import[ElectrostaticMatrixTablefname];
]

ElectrostaticMatrixTable::usage = "ElectrostaticMatrixTable[{n, SL, SpLp}] provides the calculated result of ElectrostaticMatrix[n, SL, SpLp]."; 

PrintTemporary["Loading table of matrix elements for Vk1 ..."];
ReducedVk1TableFname = FileNameJoin[{moduleDir, "data", "ReducedVk1Table.m"}];
If[!FileExistsQ[ReducedVk1TableFname],
  (PrintTemporary[">> ReducedVk1Table.m not found, generating ..."];
    ReducedVk1Table = GenerateVk1Table[7, True];
  ),
  ReducedVk1Table = Import[ReducedVk1TableFname];
]

PrintTemporary["Loading table of matrix elements for spin-orbit ..."];
SpinOrbitTableFname = FileNameJoin[{moduleDir, "data", "SpinOrbitTable.m"}];
If[!FileExistsQ[SpinOrbitTableFname],
  (PrintTemporary[">> SpinOrbitTable.m not found, generating ..."];
    SpinOrbitTable = GenerateSpinOrbitTable[7, True];
  ),
  SpinOrbitTable = Import[SpinOrbitTableFname];
]

PrintTemporary["Loading table of reduced AiZi matrix elements ..."];
AiZiTableFname = FileNameJoin[{moduleDir, "data", "ReducedAiZiTable.m"}];
If[!FileExistsQ[AiZiTableFname],
  (PrintTemporary[">> ReducedAiZiTable.m not found, generating ..."];
    AiZiTable = GenerateAiZiTable[7, True, True];
  ),
  AiZiTable = Import[AiZiTableFname];
];

PrintTemporary["Loading table of matrix elements for spin-other-orbit and electrostatically-correlated-spin-orbit ..."];
SOOandECSOTableFname = FileNameJoin[{moduleDir, "data", "SOOandECSOTable.m"}];
If[!FileExistsQ[SpinOrbitTableFname],
  (PrintTemporary[">> SOOandECSOTable.m not found, generating ..."];
    SOOandECSOTable = GenerateSOOandECSOTable[7, True];
  ),
  SOOandECSOTable = Import[SOOandECSOTableFname];
]

PrintTemporary["Loading table of reduced T22 matrix elements ..."];
T22TableFname = FileNameJoin[{moduleDir, "data", "ReducedT22Table.m"}];
If[!FileExistsQ[T22TableFname],
  (PrintTemporary[">> ReducedT22Table.m not found, generating ..."];
    T22Table = GenerateT22Table[7, True, True];
  ),
  T22Table = Import[T22TableFname];
]

PrintTemporary["Loading table of matrix elements for spin-spin ..."];
SpinSpinTableFname = FileNameJoin[{moduleDir, "data", "SpinSpinTable.m"}];
If[!FileExistsQ[SpinSpinTableFname],
  (PrintTemporary[">> SpinSpinTable.m not found, generating ..."];
    SpinSpinTable = GenerateSpinSpinTable[7, True];
  ),
  SpinSpinTable = Import[SpinSpinTableFname];
]

PrintTemporary["Loading table of matrix elements for three-body configuration-interaction effects ..."];
ThreeBodyFname = FileNameJoin[{moduleDir, "data", "ThreeBodyTable.m"}];
If[!FileExistsQ[ThreeBodyFname],
  (PrintTemporary[">> ThreeBody.m not found, generating ..."];
    ThreeBodyTable = GenerateThreeBodyTable[7, True];
  ),
  ThreeBodyTable = Import[ThreeBodyFname];
]

(* ############################################ Data ############################################ *)
(* ############################################################################################## *)



EndPackage[]
