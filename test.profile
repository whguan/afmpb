Flat profile:

Each sample counts as 0.01 seconds.
  %   cumulative   self              self     total           
 time   seconds   seconds    calls   s/call   s/call  name    
 11.85      0.48     0.48  6264801     0.00     0.00  sLapGreenFunction
  9.63      0.87     0.39   749758     0.00     0.00  NSingularCoeff
  7.16      1.16     0.29  2170106     0.00     0.00  dYukGreenFunction
  6.67      1.43     0.27                             exp.L
  5.93      1.67     0.24  2088328     0.00     0.00  sYukGreenFunction
  4.94      1.87     0.20        1     0.20     0.61  SetSelfMatrix
  4.69      2.06     0.19   228486     0.00     0.00  lgndr
  4.20      2.23     0.17     2487     0.00     0.00  LapLocalToTarget
  3.95      2.39     0.16     2548     0.00     0.00  YukLocalToTarget
  3.70      2.54     0.15    35947     0.00     0.00  DirectEvaluation
  3.21      2.67     0.13  1856214     0.00     0.00  dLapGreenFunction
  2.72      2.78     0.11   115510     0.00     0.00  ribesl_
  2.22      2.87     0.09      520     0.00     0.00  LapLocalToLocal
  1.98      2.95     0.08    12742     0.00     0.00  LapExponentialToLocalPhase1
  1.98      3.03     0.08     3611     0.00     0.00  dYukSourceToMultipole
  1.98      3.11     0.08      507     0.00     0.00  YukLocalToLocal
  1.73      3.18     0.07     3815     0.00     0.00  dLapSourceToMultipole
  1.73      3.25     0.07    22320     0.00     0.00  rotz2x
  1.48      3.31     0.06    16492     0.00     0.00  rotz2y
  1.48      3.37     0.06    15396     0.00     0.00  MakeDList
  1.23      3.42     0.05    46955     0.00     0.00  YukMultipoleToExponentialPhase2
  1.23      3.47     0.05    15232     0.00     0.00  MakeUList
  1.23      3.52     0.05    13396     0.00     0.00  YukExponentialToLocalPhase1
  1.23      3.57     0.05     9003     0.00     0.00  YukExponentialToLocalPhase2
  1.23      3.62     0.05     1481     0.00     0.00  LapMultipoleToMultipole
  1.23      3.67     0.05     1358     0.00     0.00  YukMultipoleToMultipole
  1.23      3.72     0.05                             __intel_memset
  0.99      3.76     0.04    26421     0.00     0.00  LapMultipoleToExponentialPhase1
  0.99      3.80     0.04     8745     0.00     0.00  LapExponentialToLocalPhase2
  0.74      3.83     0.03    52173     0.00     0.00  LapMultipoleToExponentialPhase2
  0.74      3.86     0.03    23491     0.00     0.00  YukMultipoleToExponentialPhase1
  0.74      3.89     0.03      620     0.00     0.00  YukExponentialToLocal
  0.74      3.92     0.03                             exp
  0.49      3.94     0.02     5892     0.00     0.00  roty2z
  0.49      3.96     0.02     3904     0.00     0.00  sLapSourceToMultipole
  0.49      3.98     0.02     1273     0.00     0.00  CloseCoeff
  0.25      3.99     0.01   107136     0.00     0.00  in
  0.25      4.00     0.01   103554     0.00     0.00  dgamma_
  0.25      4.01     0.01    63796     0.00     0.00  UpdateList
  0.25      4.02     0.01     3370     0.00     0.00  sYukSourceToMultipole
  0.25      4.03     0.01      646     0.00     0.00  LapExponentialToLocal
  0.25      4.04     0.01      277     0.00     0.00  PartitionBox
  0.25      4.05     0.01                             _intel_fast_memset.J
  0.00      4.05     0.00    31025     0.00     0.00  PushStack
  0.00      4.05     0.00    18993     0.00     0.00  IfAdjacent
  0.00      4.05     0.00    17014     0.00     0.00  AggregateSweep
  0.00      4.05     0.00    11545     0.00     0.00  PopStack
  0.00      4.05     0.00     8515     0.00     0.00  LapMultipoleToExponential
  0.00      4.05     0.00     7369     0.00     0.00  YukMultipoleToExponential
  0.00      4.05     0.00     6090     0.00     0.00  ProcessList4
  0.00      4.05     0.00     5783     0.00     0.00  DisAggregateSweep
  0.00      4.05     0.00     1840     0.00     0.00  PopAll
  0.00      4.05     0.00     1134     0.00     0.00  BuildMergedList2
  0.00      4.05     0.00      756     0.00     0.00  BuildFinerList
  0.00      4.05     0.00       86     0.00     0.00  ProcessList13
  0.00      4.05     0.00       46     0.00     0.07  AdapFMMCompute
  0.00      4.05     0.00       36     0.00     0.00  lgndrgt1
  0.00      4.05     0.00       24     0.00     0.00  BuildList
  0.00      4.05     0.00       12     0.00     0.00  fom_
  0.00      4.05     0.00       12     0.00     0.00  fstrtn
  0.00      4.05     0.00       12     0.00     0.00  gmres_
  0.00      4.05     0.00       11     0.00     0.27  MatrixVectorMultiply
  0.00      4.05     0.00       10     0.00     0.00  givens_
  0.00      4.05     0.00        8     0.00     0.00  numthetafour
  0.00      4.05     0.00        5     0.00     0.00  ylcshftcoef
  0.00      4.05     0.00        5     0.00     0.00  ympshftcoef
  0.00      4.05     0.00        4     0.00     0.00  bnlcft
  0.00      4.05     0.00        4     0.00     0.00  numthetahalf
  0.00      4.05     0.00        4     0.00     0.00  yhfstrtn
  0.00      4.05     0.00        4     0.00     0.00  ymkexps
  0.00      4.05     0.00        4     0.00     0.00  ymkfexp
  0.00      4.05     0.00        4     0.00     0.00  yrlscini
  0.00      4.05     0.00        3     0.00     0.00  BuildGraph
  0.00      4.05     0.00        3     0.00     0.00  DestroyGraph
  0.00      4.05     0.00        3     0.00     0.00  FMMClean
  0.00      4.05     0.00        3     0.00     0.00  LapFMMClean
  0.00      4.05     0.00        3     0.00     0.00  LapFMMInit
  0.00      4.05     0.00        3     0.00     0.00  frmini
  0.00      4.05     0.00        3     0.00     0.00  lapvwts
  0.00      4.05     0.00        3     0.00     0.00  mkexps
  0.00      4.05     0.00        3     0.00     0.00  mkfexp
  0.00      4.05     0.00        3     0.00     0.00  rlscini
  0.00      4.05     0.00        3     0.00     0.00  rotgen
  0.00      4.05     0.00        2     0.00     0.00  wctime
  0.00      4.05     0.00        1     0.00     0.00  BuildDirectList13
  0.00      4.05     0.00        1     0.00     0.00  CleanDirectList13
  0.00      4.05     0.00        1     0.00     0.00  ComputeSolvationEnergy
  0.00      4.05     0.00        1     0.00     0.00  ParseCommandLine
  0.00      4.05     0.00        1     0.00     0.00  ProcessElementGeometry
  0.00      4.05     0.00        1     0.00     0.00  ProcessPQRFile
  0.00      4.05     0.00        1     0.00     0.00  ReadMeshFile
  0.00      4.05     0.00        1     0.00     0.00  ReadPQRFile
  0.00      4.05     0.00        1     0.00     0.00  RemoveIsolatedNodes
  0.00      4.05     0.00        1     0.00     0.07  SetGNDGN
  0.00      4.05     0.00        1     0.00     0.00  SetParameters
  0.00      4.05     0.00        1     0.00     0.07  SetRightHandSide
  0.00      4.05     0.00        1     0.00     3.62  SolvePoissonBoltzmann
  0.00      4.05     0.00        1     0.00     0.00  WritePotential
  0.00      4.05     0.00        1     0.00     0.00  YukFMMClean
  0.00      4.05     0.00        1     0.00     0.00  YukFMMInit
  0.00      4.05     0.00        1     0.00     0.00  implu_
  0.00      4.05     0.00        1     0.00     0.00  yhfrmini
  0.00      4.05     0.00        1     0.00     0.00  yhrotgen
  0.00      4.05     0.00        1     0.00     0.00  yukvwts

 %         the percentage of the total running time of the
time       program used by this function.

cumulative a running sum of the number of seconds accounted
 seconds   for by this function and those listed above it.

 self      the number of seconds accounted for by this
seconds    function alone.  This is the major sort for this
           listing.

calls      the number of times this function was invoked, if
           this function is profiled, else blank.
 
 self      the average number of milliseconds spent in this
ms/call    function per call, if this function is profiled,
	   else blank.

 total     the average number of milliseconds spent in this
ms/call    function and its descendents per call, if this 
	   function is profiled, else blank.

name       the name of the function.  This is the minor sort
           for this listing. The index shows the location of
	   the function in the gprof listing. If the index is
	   in parenthesis it shows where it would appear in
	   the gprof listing if it were to be printed.

		     Call graph (explanation follows)


granularity: each sample hit covers 2 byte(s) for 0.25% of 4.05 seconds

index % time    self  children    called     name
                                                 <spontaneous>
[1]     91.1    0.00    3.69                 main [1]
                0.00    3.62       1/1           SolvePoissonBoltzmann [2]
                0.00    0.07       1/1           SetGNDGN [31]
                0.00    0.00       2/2           wctime [90]
                0.00    0.00       1/1           ParseCommandLine [94]
                0.00    0.00       1/1           SetParameters [100]
                0.00    0.00       1/1           ProcessPQRFile [96]
                0.00    0.00       1/1           ReadPQRFile [98]
                0.00    0.00       1/1           ReadMeshFile [97]
                0.00    0.00       1/1           RemoveIsolatedNodes [99]
                0.00    0.00       1/1           ProcessElementGeometry [95]
                0.00    0.00       1/1           WritePotential [101]
                0.00    0.00       1/1           ComputeSolvationEnergy [93]
-----------------------------------------------
                0.00    3.62       1/1           main [1]
[2]     89.4    0.00    3.62       1         SolvePoissonBoltzmann [2]
                0.00    2.94      11/11          MatrixVectorMultiply [5]
                0.20    0.41       1/1           SetSelfMatrix [9]
                0.00    0.07       1/1           SetRightHandSide [32]
                0.00    0.00       1/3           BuildGraph [55]
                0.00    0.00       1/1           YukFMMInit [57]
                0.00    0.00      12/12          gmres_ [71]
                0.00    0.00       1/3           LapFMMInit [83]
                0.00    0.00       1/24          BuildList [68]
                0.00    0.00       1/1           BuildDirectList13 [91]
                0.00    0.00       1/1           CleanDirectList13 [92]
                0.00    0.00       1/1           YukFMMClean [102]
                0.00    0.00       1/3           LapFMMClean [82]
                0.00    0.00       1/3           DestroyGraph [80]
                0.00    0.00       1/3           FMMClean [81]
-----------------------------------------------
                                  46             AdapFMMCompute [3]
                0.00    0.07       1/46          SetGNDGN [31]
                0.00    0.07       1/46          SetRightHandSide [32]
                0.00    2.94      44/46          MatrixVectorMultiply [5]
[3]     75.8    0.00    3.07      46+46      AdapFMMCompute [3]
                0.15    2.92      46/46          DisAggregateSweep <cycle 1> [7]
                0.00    0.00      24/1481        LapMultipoleToMultipole [41]
                0.00    0.00      22/1358        YukMultipoleToMultipole [42]
                                  46             AdapFMMCompute [3]
-----------------------------------------------
[4]     75.8    0.15    2.92      46+64874   <cycle 1 as a whole> [4]
                0.15    1.14   35947             DirectEvaluation <cycle 1> [6]
                0.00    1.00    5783             DisAggregateSweep <cycle 1> [7]
                0.00    0.00    6090             ProcessList4 <cycle 1> [63]
                0.00    0.00      86             ProcessList13 <cycle 1> [66]
-----------------------------------------------
                0.00    2.94      11/11          SolvePoissonBoltzmann [2]
[5]     72.5    0.00    2.94      11         MatrixVectorMultiply [5]
                0.00    2.94      44/46          AdapFMMCompute [3]
-----------------------------------------------
                                 831             ProcessList13 <cycle 1> [66]
                               35116             ProcessList4 <cycle 1> [63]
[6]     31.9    0.15    1.14   35947         DirectEvaluation <cycle 1> [6]
                0.48    0.00 6264801/6264801     sLapGreenFunction [10]
                0.29    0.00 2170106/2170106     dYukGreenFunction [12]
                0.24    0.00 2088328/2088328     sYukGreenFunction [14]
                0.13    0.00 1856214/1856214     dLapGreenFunction [23]
                               17014             AggregateSweep <cycle 1> [8]
-----------------------------------------------
                                5737             AggregateSweep <cycle 1> [8]
                0.15    2.92      46/46          AdapFMMCompute [3]
[7]     24.7    0.00    1.00    5783         DisAggregateSweep <cycle 1> [7]
                0.16    0.06    2548/2548        YukLocalToTarget [15]
                0.03    0.18     620/620         YukExponentialToLocal [16]
                0.01    0.20     646/646         LapExponentialToLocal [17]
                0.17    0.03    2487/2487        LapLocalToTarget [18]
                0.09    0.00     520/520         LapLocalToLocal [27]
                0.08    0.00     507/507         YukLocalToLocal [30]
                                6090             ProcessList4 <cycle 1> [63]
-----------------------------------------------
                                5315             AggregateSweep <cycle 1> [8]
                               17014             DirectEvaluation <cycle 1> [6]
[8]     19.2    0.00    0.78   17014+5315    AggregateSweep <cycle 1> [8]
                0.08    0.08    3611/3611        dYukSourceToMultipole [20]
                0.00    0.13    7369/7369        YukMultipoleToExponential [21]
                0.00    0.13    8515/8515        LapMultipoleToExponential [24]
                0.07    0.04    3815/3815        dLapSourceToMultipole [26]
                0.01    0.08    3370/3370        sYukSourceToMultipole [28]
                0.02    0.04    3904/3904        sLapSourceToMultipole [36]
                0.05    0.00    1336/1358        YukMultipoleToMultipole [42]
                0.05    0.00    1457/1481        LapMultipoleToMultipole [41]
                                5737             DisAggregateSweep <cycle 1> [7]
                                5315             AggregateSweep <cycle 1> [8]
-----------------------------------------------
                0.20    0.41       1/1           SolvePoissonBoltzmann [2]
[9]     15.1    0.20    0.41       1         SetSelfMatrix [9]
                0.39    0.00  749758/749758      NSingularCoeff [11]
                0.02    0.00    1273/1273        CloseCoeff [50]
-----------------------------------------------
                0.48    0.00 6264801/6264801     DirectEvaluation <cycle 1> [6]
[10]    11.9    0.48    0.00 6264801         sLapGreenFunction [10]
-----------------------------------------------
                0.39    0.00  749758/749758      SetSelfMatrix [9]
[11]     9.6    0.39    0.00  749758         NSingularCoeff [11]
-----------------------------------------------
                0.29    0.00 2170106/2170106     DirectEvaluation <cycle 1> [6]
[12]     7.2    0.29    0.00 2170106         dYukGreenFunction [12]
-----------------------------------------------
                                                 <spontaneous>
[13]     6.7    0.27    0.00                 exp.L [13]
-----------------------------------------------
                0.24    0.00 2088328/2088328     DirectEvaluation <cycle 1> [6]
[14]     5.9    0.24    0.00 2088328         sYukGreenFunction [14]
-----------------------------------------------
                0.16    0.06    2548/2548        DisAggregateSweep <cycle 1> [7]
[15]     5.3    0.16    0.06    2548         YukLocalToTarget [15]
                0.00    0.03   27516/107136      in [22]
                0.02    0.00   26562/228486      lgndr [19]
-----------------------------------------------
                0.03    0.18     620/620         DisAggregateSweep <cycle 1> [7]
[16]     5.2    0.03    0.18     620         YukExponentialToLocal [16]
                0.05    0.00   13396/13396       YukExponentialToLocalPhase1 [39]
                0.05    0.00    9003/9003        YukExponentialToLocalPhase2 [40]
                0.03    0.00    7790/15396       MakeDList [35]
                0.03    0.00    7810/15232       MakeUList [38]
                0.01    0.00    2988/5892        roty2z [49]
                0.01    0.00    2910/22320       rotz2x [33]
                0.00    0.00     558/1134        BuildMergedList2 [53]
-----------------------------------------------
                0.01    0.20     646/646         DisAggregateSweep <cycle 1> [7]
[17]     5.1    0.01    0.20     646         LapExponentialToLocal [17]
                0.08    0.00   12742/12742       LapExponentialToLocalPhase1 [29]
                0.04    0.00    8745/8745        LapExponentialToLocalPhase2 [45]
                0.03    0.00    7606/15396       MakeDList [35]
                0.02    0.00    7422/15232       MakeUList [38]
                0.01    0.00    2904/5892        roty2z [49]
                0.01    0.00    2888/22320       rotz2x [33]
                0.00    0.01     576/1134        BuildMergedList2 [53]
-----------------------------------------------
                0.17    0.03    2487/2487        DisAggregateSweep <cycle 1> [7]
[18]     4.9    0.17    0.03    2487         LapLocalToTarget [18]
                0.03    0.00   33045/228486      lgndr [19]
-----------------------------------------------
                0.02    0.00   26562/228486      YukLocalToTarget [15]
                0.03    0.00   33045/228486      LapLocalToTarget [18]
                0.03    0.00   38025/228486      sYukSourceToMultipole [28]
                0.03    0.00   41439/228486      dYukSourceToMultipole [20]
                0.04    0.00   44624/228486      sLapSourceToMultipole [36]
                0.04    0.00   44791/228486      dLapSourceToMultipole [26]
[19]     4.7    0.19    0.00  228486         lgndr [19]
-----------------------------------------------
                0.08    0.08    3611/3611        AggregateSweep <cycle 1> [8]
[20]     4.1    0.08    0.08    3611         dYukSourceToMultipole [20]
                0.00    0.05   40996/107136      in [22]
                0.03    0.00   41439/228486      lgndr [19]
-----------------------------------------------
                0.00    0.13    7369/7369        AggregateSweep <cycle 1> [8]
[21]     3.3    0.00    0.13    7369         YukMultipoleToExponential [21]
                0.05    0.00   46955/46955       YukMultipoleToExponentialPhase2 [37]
                0.03    0.00   23491/23491       YukMultipoleToExponentialPhase1 [47]
                0.03    0.00    7990/16492       rotz2y [34]
                0.02    0.00    7864/22320       rotz2x [33]
-----------------------------------------------
                0.00    0.00       5/107136      ympshftcoef [59]
                0.00    0.00       5/107136      ylcshftcoef [58]
                0.00    0.03   27516/107136      YukLocalToTarget [15]
                0.00    0.04   38614/107136      sYukSourceToMultipole [28]
                0.00    0.05   40996/107136      dYukSourceToMultipole [20]
[22]     3.2    0.01    0.12  107136         in [22]
                0.11    0.01  115510/115510      ribesl_ [25]
-----------------------------------------------
                0.13    0.00 1856214/1856214     DirectEvaluation <cycle 1> [6]
[23]     3.2    0.13    0.00 1856214         dLapGreenFunction [23]
-----------------------------------------------
                0.00    0.13    8515/8515        AggregateSweep <cycle 1> [8]
[24]     3.2    0.00    0.13    8515         LapMultipoleToExponential [24]
                0.04    0.00   26421/26421       LapMultipoleToExponentialPhase1 [44]
                0.03    0.00    8502/16492       rotz2y [34]
                0.03    0.00   52173/52173       LapMultipoleToExponentialPhase2 [46]
                0.03    0.00    8658/22320       rotz2x [33]
-----------------------------------------------
                0.11    0.01  115510/115510      in [22]
[25]     3.0    0.11    0.01  115510         ribesl_ [25]
                0.01    0.00  103554/103554      dgamma_ [51]
-----------------------------------------------
                0.07    0.04    3815/3815        AggregateSweep <cycle 1> [8]
[26]     2.6    0.07    0.04    3815         dLapSourceToMultipole [26]
                0.04    0.00   44791/228486      lgndr [19]
-----------------------------------------------
                0.09    0.00     520/520         DisAggregateSweep <cycle 1> [7]
[27]     2.2    0.09    0.00     520         LapLocalToLocal [27]
-----------------------------------------------
                0.01    0.08    3370/3370        AggregateSweep <cycle 1> [8]
[28]     2.2    0.01    0.08    3370         sYukSourceToMultipole [28]
                0.00    0.04   38614/107136      in [22]
                0.03    0.00   38025/228486      lgndr [19]
-----------------------------------------------
                0.08    0.00   12742/12742       LapExponentialToLocal [17]
[29]     2.0    0.08    0.00   12742         LapExponentialToLocalPhase1 [29]
-----------------------------------------------
                0.08    0.00     507/507         DisAggregateSweep <cycle 1> [7]
[30]     2.0    0.08    0.00     507         YukLocalToLocal [30]
-----------------------------------------------
                0.00    0.07       1/1           main [1]
[31]     1.7    0.00    0.07       1         SetGNDGN [31]
                0.00    0.07       1/46          AdapFMMCompute [3]
                0.00    0.00       1/3           BuildGraph [55]
                0.00    0.00       1/3           LapFMMInit [83]
                0.00    0.00       1/24          BuildList [68]
                0.00    0.00       1/3           FMMClean [81]
                0.00    0.00       1/3           LapFMMClean [82]
                0.00    0.00       1/3           DestroyGraph [80]
-----------------------------------------------
                0.00    0.07       1/1           SolvePoissonBoltzmann [2]
[32]     1.7    0.00    0.07       1         SetRightHandSide [32]
                0.00    0.07       1/46          AdapFMMCompute [3]
                0.00    0.00       1/3           BuildGraph [55]
                0.00    0.00       1/3           LapFMMInit [83]
                0.00    0.00       1/24          BuildList [68]
                0.00    0.00       1/3           DestroyGraph [80]
                0.00    0.00       1/3           FMMClean [81]
                0.00    0.00       1/3           LapFMMClean [82]
-----------------------------------------------
                0.01    0.00    2888/22320       LapExponentialToLocal [17]
                0.01    0.00    2910/22320       YukExponentialToLocal [16]
                0.02    0.00    7864/22320       YukMultipoleToExponential [21]
                0.03    0.00    8658/22320       LapMultipoleToExponential [24]
[33]     1.7    0.07    0.00   22320         rotz2x [33]
-----------------------------------------------
                0.03    0.00    7990/16492       YukMultipoleToExponential [21]
                0.03    0.00    8502/16492       LapMultipoleToExponential [24]
[34]     1.5    0.06    0.00   16492         rotz2y [34]
-----------------------------------------------
                0.03    0.00    7606/15396       LapExponentialToLocal [17]
                0.03    0.00    7790/15396       YukExponentialToLocal [16]
[35]     1.5    0.06    0.00   15396         MakeDList [35]
-----------------------------------------------
                0.02    0.04    3904/3904        AggregateSweep <cycle 1> [8]
[36]     1.4    0.02    0.04    3904         sLapSourceToMultipole [36]
                0.04    0.00   44624/228486      lgndr [19]
-----------------------------------------------
                0.05    0.00   46955/46955       YukMultipoleToExponential [21]
[37]     1.2    0.05    0.00   46955         YukMultipoleToExponentialPhase2 [37]
-----------------------------------------------
                0.02    0.00    7422/15232       LapExponentialToLocal [17]
                0.03    0.00    7810/15232       YukExponentialToLocal [16]
[38]     1.2    0.05    0.00   15232         MakeUList [38]
-----------------------------------------------
                0.05    0.00   13396/13396       YukExponentialToLocal [16]
[39]     1.2    0.05    0.00   13396         YukExponentialToLocalPhase1 [39]
-----------------------------------------------
                0.05    0.00    9003/9003        YukExponentialToLocal [16]
[40]     1.2    0.05    0.00    9003         YukExponentialToLocalPhase2 [40]
-----------------------------------------------
                0.00    0.00      24/1481        AdapFMMCompute [3]
                0.05    0.00    1457/1481        AggregateSweep <cycle 1> [8]
[41]     1.2    0.05    0.00    1481         LapMultipoleToMultipole [41]
-----------------------------------------------
                0.00    0.00      22/1358        AdapFMMCompute [3]
                0.05    0.00    1336/1358        AggregateSweep <cycle 1> [8]
[42]     1.2    0.05    0.00    1358         YukMultipoleToMultipole [42]
-----------------------------------------------
                                                 <spontaneous>
[43]     1.2    0.05    0.00                 __intel_memset [43]
-----------------------------------------------
                0.04    0.00   26421/26421       LapMultipoleToExponential [24]
[44]     1.0    0.04    0.00   26421         LapMultipoleToExponentialPhase1 [44]
-----------------------------------------------
                0.04    0.00    8745/8745        LapExponentialToLocal [17]
[45]     1.0    0.04    0.00    8745         LapExponentialToLocalPhase2 [45]
-----------------------------------------------
                0.03    0.00   52173/52173       LapMultipoleToExponential [24]
[46]     0.7    0.03    0.00   52173         LapMultipoleToExponentialPhase2 [46]
-----------------------------------------------
                0.03    0.00   23491/23491       YukMultipoleToExponential [21]
[47]     0.7    0.03    0.00   23491         YukMultipoleToExponentialPhase1 [47]
-----------------------------------------------
                                                 <spontaneous>
[48]     0.7    0.03    0.00                 exp [48]
-----------------------------------------------
                0.01    0.00    2904/5892        LapExponentialToLocal [17]
                0.01    0.00    2988/5892        YukExponentialToLocal [16]
[49]     0.5    0.02    0.00    5892         roty2z [49]
-----------------------------------------------
                0.02    0.00    1273/1273        SetSelfMatrix [9]
[50]     0.5    0.02    0.00    1273         CloseCoeff [50]
-----------------------------------------------
                0.01    0.00  103554/103554      ribesl_ [25]
[51]     0.2    0.01    0.00  103554         dgamma_ [51]
-----------------------------------------------
                0.01    0.00   63796/63796       BuildMergedList2 [53]
[52]     0.2    0.01    0.00   63796         UpdateList [52]
-----------------------------------------------
                0.00    0.00     558/1134        YukExponentialToLocal [16]
                0.00    0.01     576/1134        LapExponentialToLocal [17]
[53]     0.2    0.00    0.01    1134         BuildMergedList2 [53]
                0.01    0.00   63796/63796       UpdateList [52]
-----------------------------------------------
                0.01    0.00     277/277         BuildGraph [55]
[54]     0.2    0.01    0.00     277         PartitionBox [54]
-----------------------------------------------
                0.00    0.00       1/3           SetGNDGN [31]
                0.00    0.00       1/3           SolvePoissonBoltzmann [2]
                0.00    0.00       1/3           SetRightHandSide [32]
[55]     0.2    0.00    0.01       3         BuildGraph [55]
                0.01    0.00     277/277         PartitionBox [54]
                0.00    0.00      21/24          BuildList [68]
-----------------------------------------------
                                                 <spontaneous>
[56]     0.2    0.01    0.00                 _intel_fast_memset.J [56]
-----------------------------------------------
                0.00    0.00       1/1           SolvePoissonBoltzmann [2]
[57]     0.0    0.00    0.00       1         YukFMMInit [57]
                0.00    0.00       5/5           ympshftcoef [59]
                0.00    0.00       5/5           ylcshftcoef [58]
                0.00    0.00       5/8           numthetafour [73]
                0.00    0.00       4/4           yrlscini [79]
                0.00    0.00       4/4           ymkfexp [78]
                0.00    0.00       4/4           ymkexps [77]
                0.00    0.00       1/1           yhrotgen [105]
                0.00    0.00       1/1           yhfrmini [104]
                0.00    0.00       1/4           numthetahalf [75]
                0.00    0.00       1/1           yukvwts [106]
-----------------------------------------------
                0.00    0.00       5/5           YukFMMInit [57]
[58]     0.0    0.00    0.00       5         ylcshftcoef [58]
                0.00    0.00       5/107136      in [22]
-----------------------------------------------
                0.00    0.00       5/5           YukFMMInit [57]
[59]     0.0    0.00    0.00       5         ympshftcoef [59]
                0.00    0.00       5/107136      in [22]
-----------------------------------------------
                0.00    0.00    8624/31025       BuildList [68]
                0.00    0.00   22401/31025       BuildFinerList [65]
[60]     0.0    0.00    0.00   31025         PushStack [60]
-----------------------------------------------
                0.00    0.00    7448/18993       BuildList [68]
                0.00    0.00   11545/18993       BuildFinerList [65]
[61]     0.0    0.00    0.00   18993         IfAdjacent [61]
-----------------------------------------------
                0.00    0.00   11545/11545       BuildFinerList [65]
[62]     0.0    0.00    0.00   11545         PopStack [62]
-----------------------------------------------
                                6090             DisAggregateSweep <cycle 1> [7]
[63]     0.0    0.00    0.00    6090         ProcessList4 <cycle 1> [63]
                               35116             DirectEvaluation <cycle 1> [6]
                                  86             ProcessList13 <cycle 1> [66]
-----------------------------------------------
                0.00    0.00     895/1840        BuildFinerList [65]
                0.00    0.00     945/1840        BuildList [68]
[64]     0.0    0.00    0.00    1840         PopAll [64]
-----------------------------------------------
                0.00    0.00     756/756         BuildList [68]
[65]     0.0    0.00    0.00     756         BuildFinerList [65]
                0.00    0.00   22401/31025       PushStack [60]
                0.00    0.00   11545/11545       PopStack [62]
                0.00    0.00   11545/18993       IfAdjacent [61]
                0.00    0.00     895/1840        PopAll [64]
-----------------------------------------------
                                  86             ProcessList4 <cycle 1> [63]
[66]     0.0    0.00    0.00      86         ProcessList13 <cycle 1> [66]
                                 831             DirectEvaluation <cycle 1> [6]
-----------------------------------------------
                0.00    0.00      36/36          yrlscini [79]
[67]     0.0    0.00    0.00      36         lgndrgt1 [67]
-----------------------------------------------
                                 895             BuildList [68]
                0.00    0.00       1/24          SetGNDGN [31]
                0.00    0.00       1/24          SolvePoissonBoltzmann [2]
                0.00    0.00       1/24          SetRightHandSide [32]
                0.00    0.00      21/24          BuildGraph [55]
[68]     0.0    0.00    0.00      24+895     BuildList [68]
                0.00    0.00    8624/31025       PushStack [60]
                0.00    0.00    7448/18993       IfAdjacent [61]
                0.00    0.00     945/1840        PopAll [64]
                0.00    0.00     756/756         BuildFinerList [65]
                                 895             BuildList [68]
-----------------------------------------------
                                 152             fom_ [69]
                0.00    0.00      12/12          gmres_ [71]
[69]     0.0    0.00    0.00      12+152     fom_ [69]
                                 152             fom_ [69]
-----------------------------------------------
                0.00    0.00       3/12          LapFMMInit [83]
                0.00    0.00       9/12          rotgen [89]
[70]     0.0    0.00    0.00      12         fstrtn [70]
-----------------------------------------------
                0.00    0.00      12/12          SolvePoissonBoltzmann [2]
[71]     0.0    0.00    0.00      12         gmres_ [71]
                0.00    0.00      12/12          fom_ [69]
                0.00    0.00      10/10          givens_ [72]
                0.00    0.00       1/1           implu_ [103]
-----------------------------------------------
                0.00    0.00      10/10          gmres_ [71]
[72]     0.0    0.00    0.00      10         givens_ [72]
-----------------------------------------------
                0.00    0.00       3/8           LapFMMInit [83]
                0.00    0.00       5/8           YukFMMInit [57]
[73]     0.0    0.00    0.00       8         numthetafour [73]
-----------------------------------------------
                0.00    0.00       1/4           yhrotgen [105]
                0.00    0.00       3/4           rotgen [89]
[74]     0.0    0.00    0.00       4         bnlcft [74]
-----------------------------------------------
                0.00    0.00       1/4           YukFMMInit [57]
                0.00    0.00       3/4           LapFMMInit [83]
[75]     0.0    0.00    0.00       4         numthetahalf [75]
-----------------------------------------------
                0.00    0.00       4/4           yhrotgen [105]
[76]     0.0    0.00    0.00       4         yhfstrtn [76]
-----------------------------------------------
                0.00    0.00       4/4           YukFMMInit [57]
[77]     0.0    0.00    0.00       4         ymkexps [77]
-----------------------------------------------
                0.00    0.00       4/4           YukFMMInit [57]
[78]     0.0    0.00    0.00       4         ymkfexp [78]
-----------------------------------------------
                0.00    0.00       4/4           YukFMMInit [57]
[79]     0.0    0.00    0.00       4         yrlscini [79]
                0.00    0.00      36/36          lgndrgt1 [67]
-----------------------------------------------
                0.00    0.00       1/3           SetGNDGN [31]
                0.00    0.00       1/3           SolvePoissonBoltzmann [2]
                0.00    0.00       1/3           SetRightHandSide [32]
[80]     0.0    0.00    0.00       3         DestroyGraph [80]
-----------------------------------------------
                0.00    0.00       1/3           SetGNDGN [31]
                0.00    0.00       1/3           SolvePoissonBoltzmann [2]
                0.00    0.00       1/3           SetRightHandSide [32]
[81]     0.0    0.00    0.00       3         FMMClean [81]
-----------------------------------------------
                0.00    0.00       1/3           SetGNDGN [31]
                0.00    0.00       1/3           SolvePoissonBoltzmann [2]
                0.00    0.00       1/3           SetRightHandSide [32]
[82]     0.0    0.00    0.00       3         LapFMMClean [82]
-----------------------------------------------
                0.00    0.00       1/3           SetGNDGN [31]
                0.00    0.00       1/3           SolvePoissonBoltzmann [2]
                0.00    0.00       1/3           SetRightHandSide [32]
[83]     0.0    0.00    0.00       3         LapFMMInit [83]
                0.00    0.00       3/12          fstrtn [70]
                0.00    0.00       3/3           rotgen [89]
                0.00    0.00       3/3           frmini [84]
                0.00    0.00       3/4           numthetahalf [75]
                0.00    0.00       3/3           lapvwts [85]
                0.00    0.00       3/3           rlscini [88]
                0.00    0.00       3/8           numthetafour [73]
                0.00    0.00       3/3           mkfexp [87]
                0.00    0.00       3/3           mkexps [86]
-----------------------------------------------
                0.00    0.00       3/3           LapFMMInit [83]
[84]     0.0    0.00    0.00       3         frmini [84]
-----------------------------------------------
                0.00    0.00       3/3           LapFMMInit [83]
[85]     0.0    0.00    0.00       3         lapvwts [85]
-----------------------------------------------
                0.00    0.00       3/3           LapFMMInit [83]
[86]     0.0    0.00    0.00       3         mkexps [86]
-----------------------------------------------
                0.00    0.00       3/3           LapFMMInit [83]
[87]     0.0    0.00    0.00       3         mkfexp [87]
-----------------------------------------------
                0.00    0.00       3/3           LapFMMInit [83]
[88]     0.0    0.00    0.00       3         rlscini [88]
-----------------------------------------------
                0.00    0.00       3/3           LapFMMInit [83]
[89]     0.0    0.00    0.00       3         rotgen [89]
                0.00    0.00       9/12          fstrtn [70]
                0.00    0.00       3/4           bnlcft [74]
-----------------------------------------------
                0.00    0.00       2/2           main [1]
[90]     0.0    0.00    0.00       2         wctime [90]
-----------------------------------------------
                0.00    0.00       1/1           SolvePoissonBoltzmann [2]
[91]     0.0    0.00    0.00       1         BuildDirectList13 [91]
-----------------------------------------------
                0.00    0.00       1/1           SolvePoissonBoltzmann [2]
[92]     0.0    0.00    0.00       1         CleanDirectList13 [92]
-----------------------------------------------
                0.00    0.00       1/1           main [1]
[93]     0.0    0.00    0.00       1         ComputeSolvationEnergy [93]
-----------------------------------------------
                0.00    0.00       1/1           main [1]
[94]     0.0    0.00    0.00       1         ParseCommandLine [94]
-----------------------------------------------
                0.00    0.00       1/1           main [1]
[95]     0.0    0.00    0.00       1         ProcessElementGeometry [95]
-----------------------------------------------
                0.00    0.00       1/1           main [1]
[96]     0.0    0.00    0.00       1         ProcessPQRFile [96]
-----------------------------------------------
                0.00    0.00       1/1           main [1]
[97]     0.0    0.00    0.00       1         ReadMeshFile [97]
-----------------------------------------------
                0.00    0.00       1/1           main [1]
[98]     0.0    0.00    0.00       1         ReadPQRFile [98]
-----------------------------------------------
                0.00    0.00       1/1           main [1]
[99]     0.0    0.00    0.00       1         RemoveIsolatedNodes [99]
-----------------------------------------------
                0.00    0.00       1/1           main [1]
[100]    0.0    0.00    0.00       1         SetParameters [100]
-----------------------------------------------
                0.00    0.00       1/1           main [1]
[101]    0.0    0.00    0.00       1         WritePotential [101]
-----------------------------------------------
                0.00    0.00       1/1           SolvePoissonBoltzmann [2]
[102]    0.0    0.00    0.00       1         YukFMMClean [102]
-----------------------------------------------
                0.00    0.00       1/1           gmres_ [71]
[103]    0.0    0.00    0.00       1         implu_ [103]
-----------------------------------------------
                0.00    0.00       1/1           YukFMMInit [57]
[104]    0.0    0.00    0.00       1         yhfrmini [104]
-----------------------------------------------
                0.00    0.00       1/1           YukFMMInit [57]
[105]    0.0    0.00    0.00       1         yhrotgen [105]
                0.00    0.00       4/4           yhfstrtn [76]
                0.00    0.00       1/4           bnlcft [74]
-----------------------------------------------
                0.00    0.00       1/1           YukFMMInit [57]
[106]    0.0    0.00    0.00       1         yukvwts [106]
-----------------------------------------------

 This table describes the call tree of the program, and was sorted by
 the total amount of time spent in each function and its children.

 Each entry in this table consists of several lines.  The line with the
 index number at the left hand margin lists the current function.
 The lines above it list the functions that called this function,
 and the lines below it list the functions this one called.
 This line lists:
     index	A unique number given to each element of the table.
		Index numbers are sorted numerically.
		The index number is printed next to every function name so
		it is easier to look up where the function in the table.

     % time	This is the percentage of the `total' time that was spent
		in this function and its children.  Note that due to
		different viewpoints, functions excluded by options, etc,
		these numbers will NOT add up to 100%.

     self	This is the total amount of time spent in this function.

     children	This is the total amount of time propagated into this
		function by its children.

     called	This is the number of times the function was called.
		If the function called itself recursively, the number
		only includes non-recursive calls, and is followed by
		a `+' and the number of recursive calls.

     name	The name of the current function.  The index number is
		printed after it.  If the function is a member of a
		cycle, the cycle number is printed between the
		function's name and the index number.


 For the function's parents, the fields have the following meanings:

     self	This is the amount of time that was propagated directly
		from the function into this parent.

     children	This is the amount of time that was propagated from
		the function's children into this parent.

     called	This is the number of times this parent called the
		function `/' the total number of times the function
		was called.  Recursive calls to the function are not
		included in the number after the `/'.

     name	This is the name of the parent.  The parent's index
		number is printed after it.  If the parent is a
		member of a cycle, the cycle number is printed between
		the name and the index number.

 If the parents of the function cannot be determined, the word
 `<spontaneous>' is printed in the `name' field, and all the other
 fields are blank.

 For the function's children, the fields have the following meanings:

     self	This is the amount of time that was propagated directly
		from the child into the function.

     children	This is the amount of time that was propagated from the
		child's children to the function.

     called	This is the number of times the function called
		this child `/' the total number of times the child
		was called.  Recursive calls by the child are not
		listed in the number after the `/'.

     name	This is the name of the child.  The child's index
		number is printed after it.  If the child is a
		member of a cycle, the cycle number is printed
		between the name and the index number.

 If there are any cycles (circles) in the call graph, there is an
 entry for the cycle-as-a-whole.  This entry shows who called the
 cycle (as parents) and the members of the cycle (as children.)
 The `+' recursive calls entry shows the number of function calls that
 were internal to the cycle, and the calls entry for each member shows,
 for that member, how many times it was called from other members of
 the cycle.


Index by function name

   [3] AdapFMMCompute         [66] ProcessList13          [69] fom_
   [8] AggregateSweep         [63] ProcessList4           [84] frmini
  [91] BuildDirectList13      [96] ProcessPQRFile         [70] fstrtn
  [65] BuildFinerList         [60] PushStack              [72] givens_
  [55] BuildGraph             [97] ReadMeshFile           [71] gmres_
  [68] BuildList              [98] ReadPQRFile           [103] implu_
  [53] BuildMergedList2       [99] RemoveIsolatedNodes    [22] in
  [92] CleanDirectList13      [31] SetGNDGN               [85] lapvwts
  [50] CloseCoeff            [100] SetParameters          [19] lgndr
  [93] ComputeSolvationEnergy [32] SetRightHandSide       [67] lgndrgt1
  [80] DestroyGraph            [9] SetSelfMatrix          [86] mkexps
   [6] DirectEvaluation        [2] SolvePoissonBoltzmann  [87] mkfexp
   [7] DisAggregateSweep      [52] UpdateList             [73] numthetafour
  [81] FMMClean              [101] WritePotential         [75] numthetahalf
  [61] IfAdjacent             [16] YukExponentialToLocal  [25] ribesl_
  [17] LapExponentialToLocal  [39] YukExponentialToLocalPhase1 [88] rlscini
  [29] LapExponentialToLocalPhase1 [40] YukExponentialToLocalPhase2 [89] rotgen
  [45] LapExponentialToLocalPhase2 [102] YukFMMClean      [49] roty2z
  [82] LapFMMClean            [57] YukFMMInit             [33] rotz2x
  [83] LapFMMInit             [30] YukLocalToLocal        [34] rotz2y
  [27] LapLocalToLocal        [15] YukLocalToTarget       [10] sLapGreenFunction
  [18] LapLocalToTarget       [21] YukMultipoleToExponential [36] sLapSourceToMultipole
  [24] LapMultipoleToExponential [47] YukMultipoleToExponentialPhase1 [14] sYukGreenFunction
  [44] LapMultipoleToExponentialPhase1 [37] YukMultipoleToExponentialPhase2 [28] sYukSourceToMultipole
  [46] LapMultipoleToExponentialPhase2 [42] YukMultipoleToMultipole [90] wctime
  [41] LapMultipoleToMultipole [43] __intel_memset       [104] yhfrmini
  [35] MakeDList              [56] _intel_fast_memset.J   [76] yhfstrtn
  [38] MakeUList              [74] bnlcft                [105] yhrotgen
   [5] MatrixVectorMultiply   [23] dLapGreenFunction      [58] ylcshftcoef
  [11] NSingularCoeff         [26] dLapSourceToMultipole  [77] ymkexps
  [94] ParseCommandLine       [12] dYukGreenFunction      [78] ymkfexp
  [54] PartitionBox           [20] dYukSourceToMultipole  [59] ympshftcoef
  [64] PopAll                 [51] dgamma_                [79] yrlscini
  [62] PopStack               [48] exp                   [106] yukvwts
  [95] ProcessElementGeometry [13] exp.L                   [4] <cycle 1>
