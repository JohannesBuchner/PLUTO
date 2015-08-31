(* ****************************************************************
    pload.m
    Provide basic capabilities to read PLUTO binary files
	
	Last Modified:  Jan 23, 2013
	Author:         A. Mignone (mignone@ph.unito.it)
   **************************************************************** *)

readFloatMatrix[fileName_,nx_,ny_,nv_,type_] := Module[{file, obj},
  file = OpenRead[fileName, BinaryFormat -> True];
  If [type == 4, obj = BinaryReadList[file, "Real32", nx*ny*nv]];
  If [type == 8, obj = BinaryReadList[file, "Real64", nx*ny*nv]];
  Close[file];
  Partition[Partition[obj, nx],ny]
 ]
 
(* *******************************************************
         Read Grid information from grid.out
   ******************************************************* *)
 
file = OpenRead["grid.out"];
For[i=1, i < 10, i++, 
  {str=Read[file,String], 
   If[UnsameQ[StringTake[str,1],"#"], Break[]]  
  }]; 
   
nx    = Read[file,Number];
xgrid = Partition[ReadList[file,Number,nx*3],3];

ny    = Read[file,Number];
ygrid = Partition[ReadList[file,Number,ny*3],3];

Close[file]

(* *******************************************************
     Read .out file and get various information
   ******************************************************* *)
   
nout = Input["Enter output number"]
type = Input["Enter precision type (4 = real/8 = double)"]

If [type == 4, {file = OpenRead["flt.out"], suf = ".flt"} ]
If [type == 8, {file = OpenRead["dbl.out"], suf = ".dbl"} ]

For[i=0, i <= nout, i++, {str=Read[file,String]}]; 
Close[file]

strlist = StringSplit[str]
t       = strlist[[2]]
dt      = strlist[[3]]
mode    = strlist[[5]]

n    = Length[strlist]
nvar = n-6
Print["Sim. time: "<>ToString[t]]
Print["Grid size: "<>ToString[nx]<>" x "<>ToString[ny]]
Print["File mode: "<>mode]
Print["Variables: "<>ToString[nvar]]
For[i=7, i <= n, i++, 
   {name = strlist[[i]];
    Print[ToString[i-6]<>" -> "<>strlist[[i]]];
   }
]; 

varnames = strlist[[7;;n]]

(* *******************************************************
        Set file name(s) and read data
   ******************************************************* *)

If [mode === "single_file",
  If [nout < 10, fname = "data.000"<>ToString[nout]<>suf, 
    If [nout < 100, fname = "data.00"<>ToString[nout]<>suf,
      If [nout < 1000, fname = "data.0"<>ToString[nout]<>suf,
        If [nout < 1000, fname = "data."<>ToString[nout]<>suf];
      ];
    ];
  ];
  data = readFloatMatrix[fname,nx,ny,nvar,type];
];


If [mode === "multiple_files",
  data = {};
  For[n=1, n <= nvar, n++,
     {If [nout < 10, fname = ".000"<>ToString[nout]<>suf, 
        If [nout < 100, fname = ".00"<>ToString[nout]<>suf,
          If [nout < 1000, fname = ".0"<>ToString[nout]<>suf,
            If [nout < 1000, fname = "."<>ToString[nout]<>suf];
          ];
        ];
      ];
      fname = varnames[[n]]<>fname;
      Print[fname];
      data = Join[data,readFloatMatrix[fname,nx,ny,1,type]];
    };
  ];
]


(* ** Pload binary data ** *)
(*
xmin = 0.5*(xgrid[[ 1,2]] + xgrid[[ 1,3]])
xmax = 0.5*(xgrid[[nx,2]] + xgrid[[nx,3]])

ymin = 0.5*(ygrid[[ 1,2]] + ygrid[[ 1,3]])
ymax = 0.5*(ygrid[[ny,2]] + ygrid[[ny,3]])
*)
xmin = xgrid[[1,2]]
xmax = xgrid[[nx,3]]

ymin = ygrid[[1,2]]
ymax = ygrid[[ny,3]]

nv = Input["Enter variable to plot"]

ArrayPlot[data[[nv]], ColorFunction -> "Rainbow", 
          DataReversed -> True, 
          DataRange -> {{xmin,xmax},{ymin,ymax}}, 
          Frame -> True, FrameTicks -> Automatic]

(* the following should plot a vertical cut along y at i = 40 of rho 
ListPlot[data[[1, 1 ;; ny, 40]]]  
*)


(* the following should plot a horizontal cut along x at j = 20 of vx 
ListPlot[data[[2,20]]]
*)

