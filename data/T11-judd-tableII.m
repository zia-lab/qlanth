(*
These are the reduced matrix elements of the T11 operator for f2.
1st column psi
2nd column psi'
3rd column M^0
4th column M^2
5th column M^4

Taken from Table II in Judd, Crosswhite, and Crosswhite, 1968, "Intra-Atomic Magnetic Interactions for f Electrons"
*)
f2T11redmat = {
{"1S", "3P", 6             , 2                  , 10/11                },
{"3P", "3P", -36           , -72                , -900/11              },
{"3P", "1D", -23*Sqrt[2/15], -14*Sqrt[2/15]     , (-23*Sqrt[10/3])/11  },
{"1D", "3F", 23*Sqrt[2/5]  , 6*Sqrt[2/5]        , (-39*Sqrt[10])/11    },
{"3F", "3F", -30*Sqrt[14]  , -2*Sqrt[14]        , (20*Sqrt[14])/11     },
{"3F", "1G", -6*Sqrt[11]   , 64/(3*Sqrt[11])    , -1240/(33*Sqrt[11])  },
{"1G", "3H", 39*Sqrt[2/5]  , (-728*Sqrt[2/5])/33, (-635*Sqrt[10])/363  },
{"3H", "3H", -96*Sqrt[11/5], 184/Sqrt[55]       , (208*Sqrt[5/11])/11  },
{"3H", "1I", -5*Sqrt[26]   , (-30*Sqrt[26])/11  , (-375*Sqrt[2/13])/121}
}