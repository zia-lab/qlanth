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
ee = 1.602176634 * 10^-9;

specAlphabet = "SPDFGHIKLMNOQRTUV";
moduleDir = DirectoryName[$InputFileName];

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
";
paramSymbols = StringSplit[paramAtlas, "\n"];
paramSymbols = Select[paramSymbols, # != ""& ];
paramSymbols = ToExpression[StringSplit[#, ":"][[1]]] & /@ paramSymbols;
Protect /@ paramSymbols;

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
CFPterms[n, SL] gives all the daughter and parent terms, together with the corresponding coefficients of fractional parentage, that are compatible with the given string SL in the f^n configuration 
CFPterms[n, L, S] gives all the daughter and parent terms, together with the corresponding coefficients of fractional parentage, that correspond to the given total orbital angular momentum L and total spin S n the f^n configuration. L being an integer, and S being integer or half-integer.
In all cases the output is in the shape of a list with enclosed lists having the format {daughter_term, {parent_term_1, CFP_1}, {parent_term_2, CFP_2}, ...}.
Only the one-body coefficients for f-electrons are provided.
In all cases it must be that 1 <= n <= 7.
Data for these was parsed from Velkov, \"Multi-electron coefficients of fractional parentage for the p, d, and f shells\".";

AllowedSLterms::usage = "AllowedSLterms[n] returns a list with the allowed terms in the f^n configuration, the terms are given as lists in the format {S, L}. This list may have redundancies which are compatible with the degeneracies that might correspond to the given case.";
AllowedNKSLterms::usage = "AllowedNKSLterms[n] returns a list with the allowed terms in the f^n configuration, the terms are given as strings in spectroscopic notation. An integer as the last is used to distinguish cases with degenaracy.";
maxJ::usage = "maxJ[n] gives the maximum J = S+L that corresponds to the configuration f^n.";
minJ::usage = "minJ[n] gives the minimum J = S+L that corresponds to the configuration f^n.";
AllowedSLJterms::usage = "AllowedSLJterms[n] returns a list with the allowed {S, L, J} terms in the f^n configuration, the terms are given as lists in the format {S, L, J}. This list may have redundancies which are compatible with the degeneracies that might correspond to the given case.";

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
{{energy_1, coeff_vector_1, basis}, 
{energy_2, coeff_vector_2, basis},
...}} 
and returns the same input except that now to every energy the minimum of all of them has been subtracted.";
LoadGuillotParameters::usage = "";

CFP::usage = "CFP[{n, NKSL}] provides a list whose first element echoes NKSL and whose other elements are lists with two elements the first one being the symbol of a parent term and the second being the corresponding coefficient of fractional parentage. n must satisfy 1 <= n <= 7";

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
    {n, SL, SpLp} -> Simplify[ElectrostaticMatrix[n, SL, SpLp]], 
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
    {n, 3, SL, SpLp, k} -> Simplify[ReducedUk[n, 3, SL, SpLp, k]],
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
      {n, SL, SpLp, 1} -> Simplify[ReducedVk1[n, SL, SpLp, 1]]
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

(* ############################# Spin Orbit ################################ *)
(* ######################################################################### *)

(* ######################################################################### *)
(* ################ Magnetically-Correlated Corrections #################### *)

T222PsiPsipStates[SL_, SpLp_]:= Module[{jj, PositionofState, \[Alpha], PsiPsipStates, m0, m2, m4, m, Tkk2m},
  PositionofState = 0;
  jj = 1;
  \[Alpha]==2;
  PsiPsipStates = {{"3P", "3P"},
                   {"3P", "3F"}, {"3F", "3P"},
                   {"3F", "3F"}, {"3F", "3H"},
                   {"3H", "3F"}, {"3H", "3H"}};
  m0 = {-12, 3 8/Sqrt[3], 3 8/Sqrt[3], -(4/3) Sqrt[14], 2 8/3 Sqrt[11/2], 2 8/3 Sqrt[11/2], 4/3 Sqrt[143]};
  m2 = {-24, 8/Sqrt[3], 8/Sqrt[3], 8 4/3 Sqrt[14], -(23/11) 8/3 Sqrt[11/2], -(23/11) 8/3 Sqrt[11/2], -(34/11) 4/3 Sqrt[143]};
  m4 = {-(300/11), -(100/11) 8/Sqrt[3], -(100/11) 8/Sqrt[3], -(200/11) 4/3 Sqrt[14], -(325/121) 8/3 Sqrt[11/2], -(325/121) 8/3 Sqrt[11/2], -(1325/1573) 4/3 Sqrt[143]};
  m = {m0, m2, m4};
  While[jj<= Length[PsiPsipStates]&&PositionofState==0,
    PositionofState = If[PsiPsipStates[[jj, 1]]==SL && PsiPsipStates[[jj, 2]]==SpLp,
                      jj,
                      0];
    jj++
  ];
  Tkk2m = If[PositionofState==0, 
            0, 
            (m[[1, PositionofState]] * M0
            +m[[2, PositionofState]] * M2
            +m[[3, PositionofState]] * M4)];
  Return[Tkk2m];
];

T112PsiPsipStates[SL_, SpLp_]:= Module[
  {jj, PositionofState, \[Alpha], PsiPsipStates, m0, m2, m4, m, Tkk2m},
  PositionofState = 0;
  jj = 1;
  \[Alpha]==1;
  PsiPsipStates= {{"1S", "3P"},
                  {"3P", "1S"}, {"3P", "3P"},
                  {"3P", "1D"}, {"1D", "3P"},
                  {"1D", "3F"}, {"3F", "1D"},
                  {"3F", "3F"}, {"3F", "1G"},
                  {"1G", "3F"}, {"1G", "3H"},
                  {"3H", "1G"}, {"3H", "3H"},
                  {"3H", "1I"}, {"1I", "3H"}};
  m0 = {6, 6, -36, -Sqrt[(2/15)] 23, -Sqrt[(2/15)] 23, Sqrt[2/5] 23, Sqrt[2/5] 23, -2 Sqrt[14] 15, -Sqrt[11] 6, -Sqrt[11] 6, Sqrt[2/5] 39, Sqrt[2/5] 39, -8 1/Sqrt[55] 132, -5 Sqrt[26], -5 Sqrt[26]};
  m2 = {2, 2, -72, -Sqrt[(2/15)] 14, -Sqrt[(2/15)] 14, Sqrt[2/5] 6, Sqrt[2/5] 6, -2 Sqrt[14], Sqrt[11] 64/33, Sqrt[11] 64/33, -Sqrt[(2/5)] 728/33, -Sqrt[(2/5)] 728/33, 8/Sqrt[55] 23, -(30/11) Sqrt[26], -(30/11) Sqrt[26]};
  m4 = {10/11, 10/11, -(900/11), -Sqrt[(2/15)] 115/11, -Sqrt[(2/15)] 115/11, -(195/11) Sqrt[2/5], -(195/11) Sqrt[2/5], 2 Sqrt[14] 10/11, -Sqrt[11] 1240/363, -Sqrt[11] 1240/363, -(3175/363) Sqrt[2/5], -(3175/363) Sqrt[2/5], 8/Sqrt[55] 130/11, -(375/1573) Sqrt[26], -(375/1573) Sqrt[26]};

  m = {m0, m2, m4};
  While[jj<= Length[PsiPsipStates]&&PositionofState==0,
    PositionofState = If[PsiPsipStates[[jj, 1]]==SL && PsiPsipStates[[jj, 2]]==SpLp,
                      jj,
                      0];
    jj++
  ];
  Tkk2m = If[PositionofState==0, 
          0, 
          (m[[1, PositionofState]] M0
          +m[[2, PositionofState]] M2
          +m[[3, PositionofState]] M4)
          ];
  Return[Tkk2m];
];

(* ################ Magnetically-Correlated Corrections #################### *)
(* ######################################################################### *)

(* ######################################################################### *)
(* ############################ Reduced t^11 ############################### *)

t112PsiPsipStates[SL_, SpLp_]:= Module[
  {jj, PositionofState, \[Alpha], PsiPsipStates, p0, p2, p4, p6, p, t112p},
  PositionofState = 0;
  jj = 1;
  \[Alpha]==1;
  PsiPsipStates = {{"1S", "3P"}, {"3P", "1S"}, {"3P", "3P"},
                   {"3P", "1D"}, {"1D", "3P"}, {"1D", "3F"},
                   {"3F", "1D"}, {"3F", "3F"}, {"3F", "1G"},
                   {"1G", "3F"}, {"1G", "3H"}, {"3H", "1G"},
                   {"3H", "3H"}, {"3H", "1I"}, {"1I", "3H"}};
  p0 = {-2, -2, -1, Sqrt[15/2], Sqrt[15/2], -Sqrt[10], -Sqrt[10],
      -Sqrt[14], Sqrt[11], Sqrt[11], -Sqrt[10], -Sqrt[10],
      -Sqrt[55], Sqrt[13/2], Sqrt[13/2]};
  p2 = {-105, -105, -45, Sqrt[15/2] 32, Sqrt[15/2] 32, -Sqrt[10] 9/2,
      -Sqrt[10] 9/2, Sqrt[14] 10, -Sqrt[11] 20, -Sqrt[11] 20,
      Sqrt[10] 55/2, Sqrt[10] 55/2, Sqrt[55] 25, -Sqrt[(13/2)] 0,
      -Sqrt[(13/2)] 0};
  p4 = {-231, -231, -33, -Sqrt[(15/2)] 33, -Sqrt[(15/2)] 33, Sqrt[10] 66,
      Sqrt[10] 66, Sqrt[14] 33, Sqrt[11] 32, Sqrt[11] 32, -Sqrt[10] 23,
      -Sqrt[10] 23, Sqrt[55] 51, -Sqrt[(13/2)] 21, -Sqrt[(13/2)] 21};
  p6 = {-429, -429, 1287, Sqrt[15/2] (-286), Sqrt[15/2] (-286),
      Sqrt[10] (-(429/2)), Sqrt[10] (-(429/2)), Sqrt[14] 286,
      Sqrt[11] (-104), Sqrt[11] (-104), Sqrt[10] (-(65/2)),
      Sqrt[10] (-(65/2)), Sqrt[55] 13, Sqrt[13/2] (-6), Sqrt[13/2] (-6)};

  p = {p0, p2, p4, p6};
  While[jj<= Length[PsiPsipStates]&&PositionofState==0,
        PositionofState = If[PsiPsipStates[[jj, 1]]==SL && PsiPsipStates[[jj, 2]]==SpLp,
                          jj,
                          0];
        jj++
  ];
  t112p = If[PositionofState==0,
          0, 
          (p[[1, PositionofState]] P0
          +p[[2, PositionofState]] P2
          +p[[3, PositionofState]] P4
          +p[[4, PositionofState]] P6)];
  Return[t112p];
]

(* ############################ Reduced t^11 ############################### *)
(* ######################################################################### *)

(* ######################################################################### *)
(* ############################ Reduced AiZi ############################### *)

AiZi2PsiPsipStates[SL_, SpLp_]:= Module[
  {PositionofState, kk, PsiPsipStates, a13, z13, aizi2p},
  PositionofState = 0;
  kk = 1;
  PsiPsipStates = {{"1S", "3P"}, {"3P", "1S"}, {"3P", "3P"},
                 {"3P", "1D"}, {"1D", "3P"}, {"1D", "3F"},
                 {"3F", "1D"}, {"3F", "3F"}, {"3F", "1G"},
                 {"1G", "3F"}, {"1G", "3H"}, {"3H", "1G"},
                 {"3H", "3H"}, {"3H", "1I"}, {"1I", "3H"}};
  a13 = -33 M0+3 M2+15/11 M4-6 P0+3/2 (35 P2+77 P4+143 P6);
  z13 = {2, 2, 1, 1/Sqrt[1080] (-90), 1/Sqrt[1080] (-90), Sqrt[2/405] 45, Sqrt[2/405] 45, Sqrt[14] 1, 1/Sqrt[891] (-99 ), 1/Sqrt[891] (-99 ), 990/Sqrt[98010], 990/Sqrt[98010], 55/Sqrt[55], -2574/Sqrt[1019304], -2574/Sqrt[1019304]};
  While[kk<= Length[PsiPsipStates]&&PositionofState==0,
    PositionofState = If[PsiPsipStates[[kk, 1]]==SL && PsiPsipStates[[kk, 2]]==SpLp,
                      kk,
                      0];
    kk++
  ];
  aizi2p = If[PositionofState==0, 
            0,
            Simplify[(
              T112PsiPsipStates[SL, SpLp]
              +t112PsiPsipStates[SL, SpLp]
              -a13/6 z13[[PositionofState]]
            )
            ]
  ];
  Return[aizi2p];
];

aizi[n_, SL_, SpLp_]:= Module[
  {s, l, nn, S, L, Sp, Lp, cfpSL, cfpSpLp, parentSL, Sb, Lb, parentSpLp, aizival},
  {s, l} = {1/2, 3};
  {S, L} = findSL[SL];
  {Sp, Lp} = findSL[SpLp];
  cfpSL = CFP[{n, SL}];
  cfpSpLp = CFP[{n, SpLp}];
  For[(
    aizival = 0;nn = 2),
    nn<= Length[cfpSL],
    nn++,
    (
      parentSL = cfpSL[[nn, 1]];
      {Sb, Lb} = findSL[parentSL];
      dval = Sum[(
        cfpSpLp[[mm, 2]] 
        *SixJay[{S, 1, Sp}, {findSL[cfpSpLp[[mm, 1]]][[1]], s, Sb}] 
        *SixJay[{L, 1, Lp}, {findSL[cfpSpLp[[mm, 1]]][[2]], l, Lb}] 
        *AiZiTable[{n-1, parentSL, Part[cfpSpLp, mm, 1]}]
      ),
        {mm, 2, Length[cfpSpLp]}
        ];
      dval = (cfpSL[[nn, 2]] * (-1)^(Sb+Lb+s+l+Sp+Lp) * dval);
      aizival+= dval;
    )
  ];
  Return[n/(n-2) * Sqrt[(2 S+1) (2 Sp+1) (2 L+1) (2Lp+1)] * aizival];
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
        AiZi2PsiPsipStates[SL, SpLp],
        True,
        Simplify[aizi[n,  SL, SpLp]]
      ];
    ),
    {n, 1, nmax},
    {SL, AllowedNKSLterms[n]}, 
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

T11n[n_, SL_, SpLp_]:= Module[
  {s, l, nn, S, L, Sp, Lp, cfpSL, cfpSpLp, parentSL, Sb, Lb, parentSpLp, Tnkk, dval},
  s = 1/2;
  l = 3;
  {S, L} = findSL[SL];
  {Sp, Lp} = findSL[SpLp];
  cfpSL = CFP[{n, SL}];
  cfpSpLp = CFP[{n, SpLp}];
  For[(Tnkk=0; nn=2;),
  nn<= Length[cfpSL],
  nn++,
  parentSL = cfpSL[[nn, 1]];
  {Sb, Lb} = findSL[parentSL];
  dval = Sum[(
    cfpSpLp[[mm, 2]] 
    * SixJay[{S, 1, Sp}, {findSL[cfpSpLp[[mm, 1]]][[1]], s, Sb}]
    * SixJay[{L, 1, Lp}, {findSL[cfpSpLp[[mm, 1]]][[2]], l, Lb}] 
    * T11Table[{n-1, parentSL, Part[cfpSpLp, mm, 1]}]
    ),
    {mm, 2, Length[cfpSpLp]}];
  dval = (cfpSL[[nn, 2]]
        * (-1)^(Sb+Lb+s+l+Sp+Lp)
        * dval);
  Tnkk+= dval;
  ];
  Tnkk = n/(n-2) * Sqrt[(2 S+1) (2 Sp+1) (2 L+1) (2Lp+1)] *Tnkk;
  Return[Tnkk];
];

T22n[n_, SL_, SpLp_]:= Module[
  {s, l, nn, S, L, Sp, Lp, cfpSL, cfpSpLp, parentSL, Sb, Lb, parentSpLp, Tnkk, dval},
  s = 1/2;
  l = 3;
  {S, L} = findSL[SL];
  {Sp, Lp} = findSL[SpLp];
  cfpSL = CFP[{n, SL}];
  cfpSpLp = CFP[{n, SpLp}];
  For[(Tnkk = 0;nn = 2),
  nn<= Length[cfpSL],
  nn++,
  parentSL = Part[cfpSL, nn, 1];
  {Sb, Lb} = findSL[parentSL];
  dval = Sum[(cfpSpLp[[mm, 2]] 
            * SixJay[{S, 2, Sp}, {findSL[Part[cfpSpLp, mm, 1]][[1]], s, Sb}] 
            * SixJay[{L, 2, Lp}, {findSL[Part[cfpSpLp, mm, 1]][[2]], l, Lb}] 
            * T22Table[{n-1, parentSL, Part[cfpSpLp, mm, 1]}]
  ),
  {mm, 2, Length[cfpSpLp]}
  ];
  dval = (dval
      *cfpSL[[nn, 2]]
      *(-1)^(Sb+Lb+s+l+Sp+Lp));
  Tnkk+= dval;
  ];
  Tnkk = n/(n-2) * Sqrt[(2 S+1) (2 Sp+1) (2 L+1) (2Lp+1)] *Tnkk;
  Return[Tnkk];
];

GenerateT11Table[nmax_, export_, progressIndicator_]:= (
  If[progressIndicator, (
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
    Simplify[T112PsiPsipStates[SL, SpLp]],
    True,
    Simplify[T11n[n,  SL, SpLp]]
    ];
    ),
    {n, 1, nmax},
    {SL, AllowedNKSLterms[n]},
    {SpLp, AllowedNKSLterms[n]}
  ];
  If[export,
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
    Simplify[T222PsiPsipStates[SL, SpLp]],
    True,
    Simplify[T22n[n,  SL, SpLp]]
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

SpinSpin[n_, SL_, SpLp_, J_] := Module[
  {S, L, Sp, Lp, \[Alpha]},
   \[Alpha] = 2;
   {S, L} = findSL[SL];
   {Sp, Lp} = findSL[SpLp];
   (-1)^(Sp + L + J) SixJay[{Sp, Lp, J}, {L, S, \[Alpha]}]* 
    T22Table[{n, SL, SpLp}]];

GenerateSpinSpinTable[nmax_, export_] :=
 (
  SpinSpinTable = <||>;
  PrintTemporary[Dynamic[n]];
  Do[
   SpinSpinTable[{n, SL, SpLp, J}] = (SpinSpin[n, SL, SpLp, J]);,
   {n, 1, nmax},
   {J, minJ[n], maxJ[n]},
   {SL, Map[First, AllowedNKSLforJterms[n, J]]}, {SpLp, 
    Map[First, AllowedNKSLforJterms[n, J]]}
   ];
  If[export,
   (fname = FileNameJoin[{moduleDir, "data", "SpinSpinTable.m"}];
    Export[fname, SpinSpinTable];
    )
   ];
  Return[SpinSpinTable];
  );

(* ############################# Spin-Spin ################################# *)
(* ######################################################################### *)


(* ######################################################################### *)
(* ########################## Spin-Other-Orbit ############################# *)

SOOandECSO[n_, SL_, SpLp_, J_]:= Module[
  {S, Sp, L, Lp},
  {S, L}   = findSL[SL];
  {Sp, Lp} = findSL[SpLp];
  Return[(-1)^(Sp+L+J) * SixJay[{Sp, Lp, J}, {L, S, 1}] * AiZiTable[{n, SL, SpLp}]];
]

GenerateSOOandECSOTable[nmax_, export_:False]:= (
  SOOandECSOTable = <||>;
  Do[
    SOOandECSOTable[{n, SL, SpLp, J}] = (SOOandECSO[n, SL, SpLp, J]/.{P2->P2/225, P4->P4/1089, P6->25 P6/184041});,
  {n, 1, nmax},
  {J, minJ[n], maxJ[n]},
  {SL, Map[First, AllowedNKSLforJterms[n, J]]},
  {SpLp, Map[First, AllowedNKSLforJterms[n, J]]}
  ];
  If[export,
  (
    fname = FileNameJoin[{moduleDir, "data", "SOOandECSOTable.m"}];
    Export[fname, SOOandECSOTable];
  )
  ];
  Return[SOOandECSOTable];
);

(* ########################## Spin-Other-Orbit ############################# *)
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
Bqk[q_, 2] := {B02, B12, B22}[[q + 1]];
Bqk[q_, 4] := {B04, B14, B24, B34, B44}[[q + 1]];
Bqk[q_, 6] := {B06, B16, B26, B36, B46, B56, B66}[[q + 1]];

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

AlphaL::usage = "AlphaL[SL, SpLp] returns the matrix element for the Trees operator for the two given terms SL and SpLp.";
AlphaL[SL_, SpLp_] := Module[
  {S, L},
  (
    {S, L} = findSL[SL];
    If[L==findSL[SpLp][[2]], 
      \[Alpha] L (L + 1),
      0]
  )
  ]

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
  val = (If[L==findSL[SpLp][[2]],
    \[Alpha] * L * (L + 1), 
    0] + 
   If[SL==SpLp, 
    GR7W[Part[First[findNKLSterm[SL]], 3]] + 
    GG2U[Part[First[findNKLSterm[SL]], 4]],
    0]);
  TwoBodyNKSL[S, L] = val;
  Return[val];
  ]

(* ########### Configuration-Interaction via Casimir Operators ############# *)
(* ######################################################################### *)

Options[EnergyMatrix] = {"Sparse"->True};
EnergyMatrix::usage = "EnergyMatrix[n, J, J', I, I'] provides the matrix element <J, I|H|J', I'> within the f^n configuration, it does this by adding the following interactions: Coulomb, spin-orbit, spin-other-orbit, electrostatically-correlated-spin-orbit, spin-spin, and crystal-field.";
EnergyMatrix[n_, J_, Jp_, Ii_, Ip_, CFTable_, OptionsPattern[]]:= (
  eMatrix = Table[
  subKron = ( KroneckerDelta[NKSLJM[[4]], NKSLJMp[[4]]]
          * KroneckerDelta[J, Jp]
          * KroneckerDelta[NKSLJM[[3]], NKSLJMp[[3]]]);
  matValue = If[subKron==0,
    0,
    Simplify[
      (subKron
      * ( ElectrostaticMatrixTable[{n, NKSLJM[[1]], NKSLJMp[[1]]}]
        + TwoBodyNKSL[NKSLJM[[1]], NKSLJMp[[1]]]
        + SpinOrbitTable[{n, NKSLJM[[1]], NKSLJMp[[1]], NKSLJM[[2]]}]
        + SOOandECSOTable[{n, NKSLJM[[1]], NKSLJMp[[1]], NKSLJM[[2]]}]
        + SpinSpinTable[{n, NKSLJM[[1]], NKSLJMp[[1]], NKSLJM[[2]]}]
       )
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

Options[TabulateEnergyMatrixTable] = {"Sparse"->True}
TabulateEnergyMatrixTable::usage = "TabulateEnergyMatrixTable[n, I] returns a list with three elements {EnergyMatrixTable, EnergyStatesTable, AllowedM}. EnergyMatrixTable is an Association with keys equal to lists of the form {n, J, Jp, Ii, Ii}. EnergyStatesTable is an association with keys equal to lists of the form {n, J, Ii}. AllowedM is a list with keys equal to lists of the form {n, J} and values equal to lists equal to the corresponding values of MJ.";
TabulateEnergyMatrixTable[n_, Ii_, CFTable_, OptionsPattern[]]:= (
  EnergyMatrixTable = <||>;
  EnergyStatesTable = <||>;
  AllowedM = <||>;
  Do[
    (EnergyMatrixTable[{n, J, Jp, Ii, Ii}] = EnergyMatrix[n, J, Jp, Ii, Ii, CFTable, "Sparse"->OptionValue["Sparse"]];
     EnergyStatesTable[{n, J, Ii}]       = EnergyStates[n, J, Ii];
     AllowedM[{n, J}]                   = Table[M, {J, minJ[n], maxJ[n]}, {M, -J, J}];
    ),
    {Jp, AllowedJ[n]},
    {J, AllowedJ[n]}
  ];
  Return[{EnergyMatrixTable, EnergyStatesTable, AllowedM}];
)

Options[TabulateManyEnergyMatrixTables] = {"Overwrite"->False, "Sparse"->True};
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
      {EnergyMatrixTable, EnergyStatesTable, AllowedM} = TabulateEnergyMatrixTable[n, Ii, CrystalFieldTable, "Sparse"->OptionValue["Sparse"]];
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


AllowedSLterms[n_] := Map[findSL[First[#]] &, CFPterms[n]]

AllowedNKSLterms[n_] := Map[First, CFPterms[n]]

maxJ[n_] := Max[Map[Total, AllowedSLterms[n]]]

minJ[n_] := Min[Map[Abs[Part[#, 1] - Part[#, 2]] &, AllowedSLterms[n]]]

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
  shifted      = Map[{#[[1]] - groundEnergy, #[[2]], #[[3]]} &, originalLevels];
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

SolveStates::usage = 
"SolveStates[nf, IiN, params] solves the energy values and states for an atom with n f-electrons with a nucleus of spin IiN. params is an association with the parameters of the specific ion under study.
This function also requires files for pre-computed energy matrix tables that provide the symbols EnergyMatrixTable[_, _, _, _, _].
To account for configurations f^n with n > 7, particle-hole dualities are enforced for \[Zeta] and T_i.
The unit for the returned energies is cm^-1.
-----------------------
References:
1. Sign inversion for \[Zeta]: Wybourne, Spectroscopic Properties of Rare Earths. 
2. Sign inversion for {T2, T3, T4, T6, T7, T8}: Hansen and Judd, Matrix Elements of Scalar Three Electron Operators for the Atomic f Shell.";

Options[SolveStates] = {"Return Symbolic Matrix" -> False};
SolveStates[nf_, IiN_, params_, OptionsPattern[]]:= Module[
  {n, ii, jj, JMvals},
  (*#####################################*)
  (*hole-particle equivalence enforcement*)
  n = nf;
  (*If[nf>7, (n = 14-nf;(*hole-particle equivalence*)
  {\[Zeta], T2, T3, T4, T6, T7, T8} = -{\[Zeta], T2, T3, T4, T6, T7, T8};
  )];*)
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
  SymbolicMatrix = EnergyMatrix;
  EnergyMatrix = ReplaceInSparseArray[EnergyMatrix, params];
  Print["The energy matrix has dimensions:", Dimensions[EnergyMatrix]];
  (*Solve for eigenvalues and eigenvectors.*)
  EnergyLevels = {};
  {EigenvalueJM, EigenvectorJM} = Eigensystem[EnergyMatrix];
  EigenvalueJM = Re[EigenvalueJM];
  (*There might be a very small imaginary part.*)
  (*Print[{Dimensions@EigenvalueJM, Dimensions@EigenvectorJM}];*)
  (*Parse the results for the eigenvectors in terms of the ordered basis being used.*)
  EnergyStates = {};
  Do[EnergyStates = Join[EnergyStates, EnergyStatesTable[{n, AllowedJ[n][[nn]], IiN}]],
    {nn, 1, Length[AllowedJ[n]]}
    ];
  Do[EnergyLevels = Join[EnergyLevels, {{EigenvalueJM[[nn]], EigenvectorJM[[nn]], EnergyStates}}];,
    {nn, 1, Length[EigenvalueJM]}
    ];
  If[OptionValue["Return Symbolic Matrix"],
  Return[{SymbolicMatrix, EnergyLevels}]
  ];
  Return[EnergyLevels];
];

(* ###############################              SOLVERS            ############################## *)
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
  "ExportMZip[filename, object] exports the given object and compresses it. It first exports in .m format, it then compresses that file in to a zip file, and finally the .m file is deleted. The filename must be a full path, and end with .m. This probably won't work on a PC.";
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
  "ImportMZip[filename] decompresses the provided filename, and imports the enclosed .m file that it is assumed to contain. After being imported the uncompressed file is deleted from disk. The provided filename bust be a full path, and end with .zip. This probably won't work on a PC.";
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
(* ############################################ Data ############################################ *)

PrintTemporary["Loading data from Carnall ..."];
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

PrintTemporary["Loading table of reduce AiZi matrix elements ..."];
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

PrintTemporary["Loading table of reduce T22 matrix elements ..."];
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

(* ############################################ Data ############################################ *)
(* ############################################################################################## *)

EndPackage[]
