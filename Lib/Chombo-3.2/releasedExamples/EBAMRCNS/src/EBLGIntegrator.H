#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif


#if (CH_SPACEDIM == 2)

#define QRHO  0
#define QVELX 1
#define QVELY 2
#define QPRES 3
#define QENTR 4
#define QC 5
#define QINTERN 6
#define QCVTEMP 7
#define QNUM  8
#define QSLOPE 8

#define CRHO  0
#define CMOMX 1
#define CMOMY 2
#define CENG  3
#define CNUM  4
#define FNUM  4


#endif

#if (CH_SPACEDIM == 3)

#define QRHO  0
#define QVELX 1
#define QVELY 2
#define QVELZ 3
#define QPRES 4
#define QENTR 5
#define QC    6
#define QINTERN 7
#define QCVTEMP 8
#define QNUM  9
#define QSLOPE 9


#define CRHO  0
#define CMOMX 1
#define CMOMY 2
#define CMOMZ 3
#define CENG  4
#define CNUM  5
#define FNUM  5
#endif
