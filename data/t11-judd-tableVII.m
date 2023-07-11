(*
These are the reduced matrix elements of the t11 operator for f2.
1st column psi
2nd column psi'
3rd column P^2
4th column P^4
5th column P^6

Taken from Table VII in Judd, Crosswhite, and Crosswhite, 1968, "Intra-Atomic Magnetic Interactions for f Electrons"
*)
f2t11redmat = {
{"1S", "3P", -2        , -105        , -231          , -429          },
{"3P", "3P", -1        , -45         , -33           , 1287          },
{"3P", "1D", Sqrt[15/2], 16*Sqrt[30] , -33*Sqrt[15/2], -143*Sqrt[30] },
{"1D", "3F", -Sqrt[10] , -9*Sqrt[5/2], 66*Sqrt[10]   , -429*Sqrt[5/2]},
{"3F", "3F", -Sqrt[14] , 10*Sqrt[14] , 33*Sqrt[14]   , 286*Sqrt[14]  },
{"3F", "1G", Sqrt[11]  , -20*Sqrt[11], 32*Sqrt[11]   , -104*Sqrt[11] },
{"1G", "3H", -Sqrt[10] , 55*Sqrt[5/2], -23*Sqrt[10]  , -65*Sqrt[5/2] },
{"3H", "3H", -Sqrt[55] , 25*Sqrt[55] , 51*Sqrt[55]   , 13*Sqrt[55]   },
{"3H", "1I", Sqrt[13/2], 0           , -21*Sqrt[13/2], -3*Sqrt[26]   }
}