#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>

#include "EBArith.H"
#include "EBArithF_F.H"
#include "EBISLayout.H"
#include "BoxIterator.H"
#include "VoFIterator.H"
#include "FaceIterator.H"
#include "BaseFab.H"
#include "PolyGeom.H"
#include "BaseEBCellFactory.H"
#include "BaseIFFactory.H"
#include <iomanip>
#include "CH_Timer.H"
#include "Stencils.H"
#include "EBIndexSpace.H"
#include "EBLevelDataOps.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "EBLoadBalance.H"
#include "NamespaceHeader.H"

Real     EBArith::s_valMax = 0.0;
VolIndex EBArith::s_vofMax = VolIndex(IntVect::Zero, 0);
Real     EBArith::s_minVolFrac = 0.0;

/******/
int
EBArith::
orderScript(int icomp, int inorm, int ncomp)
{
  return icomp + inorm*ncomp;
}
void
EBArith::
getLeastSquaresGradStenAllVoFsRad(VoFStencil&          a_stencil,
                                  Real&                a_weight,
                                  const RealVect&      a_normal,
                                  const RealVect&      a_centroid,
                                  const VolIndex&      a_vof,
                                  const EBISBox&       a_ebisBox,
                                  const RealVect&      a_dx,
                                  const ProblemDomain& a_domain,
                                  int                  a_ivar,
                                  int                  a_rad)
{
#if   CH_SPACEDIM == 2
  //  int numNeighbors = 24;
  int minStenSize  = 3;
#elif CH_SPACEDIM == 3
  //  int numNeighbors = 26;
  int minStenSize  = 7;
#else
  THIS_IS_AN_ERROR_MESSAGE__THIS_WILL_ONLY_COMPILE_WHEN_CH_SPACEDIM_IS_2_OR_3;
#endif
  Box domainBox = a_domain.domainBox();
  //  Vector<IntVect> ivSten(numNeighbors);

  IntVect iv0 = a_vof.gridIndex();
  int radius = a_rad;
  Vector<VolIndex> volSten;
  getAllVoFsInMonotonePath(volSten, a_vof, a_ebisBox, radius);


  bool dropOrder = false;
  int stenSize = volSten.size();
  if (stenSize < minStenSize)
    {
      dropOrder = true;
    }

  if (!dropOrder)
    {
      RealVect x0;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          x0[idir] = a_dx[idir] * (0.5 + a_centroid[idir] + iv0[idir]);
        }

      Vector<RealVect> xp(stenSize);
      for (int isten = 0; isten < stenSize; isten++)
        {
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              xp[isten][idir] = a_dx[idir] * (0.5 + volSten[isten].gridIndex()[idir]);
            }
        }

      Vector<RealVect> invATransAdeltaX(stenSize,RealVect::Zero);
      bool detZero = false;
      EBArith::calculateWeightingMatrix(x0, xp, invATransAdeltaX, detZero);

      a_stencil.clear();
      a_weight = 0.0;

      //if (!detZero)
      // {
          for (int isten = 0; isten < stenSize; isten++)
            {
              Real dphidnWeight = 0.0;
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  dphidnWeight -= invATransAdeltaX[isten][idir] * a_normal[idir];
                }
              
              a_stencil.add(volSten[isten],dphidnWeight, a_ivar);
              a_weight -= dphidnWeight;
            }
          // }
    }
  else
    {
      // getLeastSquaresGradStenAllVoFs(a_stencil, a_weight, a_normal, a_centroid,
      //                                a_vof, a_ebisBox,a_dx, a_domain, a_ivar);
      a_stencil.clear();
      a_weight = 0.0;
    }
}
/******/
// void
// EBArith::
// compareError(Vector<Real>&                            a_orders,
//              const Vector< LevelData<EBCellFAB>* >&   a_errorFine,
//              const Vector< LevelData<EBCellFAB>* >&   a_errorCoar,
//              const Vector< DisjointBoxLayout >&       a_gridsFine,
//              const Vector< DisjointBoxLayout >&       a_gridsCoar,
//              const Vector< EBISLayout >&              a_ebislFine,
//              const Vector< EBISLayout >&              a_ebislCoar,
//              const Vector<int>&                       a_refRat,
//              const Box&                               a_coarseDom,
//              int                                      a_testverbosity,
//              fstream*                                 a_fout,
//              Vector<string> a_names)
// {
//   CH_TIME("EBArith::compareError");

//   const Vector<int>& refRat = a_refRat;
//   const int ncomp = a_errorFine[0]->nComp();
//   const int nnorm = 3;
//   Real* normsCoar = new Real[ncomp*nnorm];
//   Real* normsFine = new Real[ncomp*nnorm];
//   a_orders.resize(ncomp*nnorm, 0.0);
//   Real* orders    = &(a_orders[0]);
//   for (int icomp = 0; icomp < ncomp; icomp++)
//     {
//       orders[icomp] = 0.0;
//       for (int inorm = 0; inorm < nnorm; inorm++)
//         {
//           normsCoar[EBArith::orderScript(icomp, inorm, ncomp)] = 0;
//           normsFine[EBArith::orderScript(icomp, inorm, ncomp)] = 0;
//         }
//     }
//   int testverbosity = a_testverbosity;
//   if (testverbosity > 1)
//   {
//     if (a_fout == NULL)
//     {
//       pout() << "==============================================" << endl;
//     }
//     else
//     {
//       *a_fout << "==============================================" << endl;
//     }
//   }
//   for (int comp = 0; comp < ncomp; comp++)
//     {
//       if (testverbosity > 1)
//       {
//         if (a_fout == NULL)
//         {
//           if (a_names.size() > comp)
//             pout() << "Comparing error in variable  " << a_names[comp] << endl;
//           else
//             pout() << "Comparing error in variable  " << comp << endl;
//         }
//         else
//         {
//           if (a_names.size() > comp)
//             (*a_fout) << "Comparing error in variable  " << a_names[comp] << endl;
//           else
//             (*a_fout) << "Comparing error in variable  " << comp << endl;
//         }
//       }
//       if (testverbosity > 1)
//       {
//         if (a_fout == NULL)
//         {
//           pout() << "==============================================" << endl;
//         }
//         else
//         {
//           (*a_fout) << "==============================================" << endl;
//         }
//       }
//       for (int itype = 0; itype < 3; itype++)
//         {
//           EBNormType::NormMode normtype;
//           if (itype == 0)
//             {
//               normtype = EBNormType::OverBoth;
//               if (testverbosity > 1)
//               {
//                 if (a_fout == NULL)
//                 {
//                   pout() << endl << "Using all uncovered cells." << endl  ;
//                 }
//                 else
//                 {
//                   (*a_fout) << endl << "Using all uncovered cells." << endl  ;
//                 }
//               }
//             }
//           else if (itype == 1)
//             {
//               normtype = EBNormType::OverOnlyRegular;
//               if (testverbosity > 1)
//               {
//                 if (a_fout == NULL)
//                 {
//                   pout() << endl << "Using only regular cells." << endl ;
//                 }
//                 else
//                 {
//                   (*a_fout) << endl << "Using only regular cells." << endl ;
//                 }
//               }
//             }
//           else
//             {
//               normtype = EBNormType::OverOnlyIrregular;
//               if (testverbosity > 1)
//               {
//                 if (a_fout == NULL)
//                 {
//                   pout() << endl << "Using only irregular cells." << endl;
//                 }
//                 else
//                 {
//                   (*a_fout) << endl << "Using only irregular cells." << endl;
//                 }
//               }
//             }

//           for (int inorm = 0; inorm <= 2; inorm++)
//             {

//               if (inorm == 0)
//                 {
//                   if (testverbosity > 1)
//                   {
//                     if (a_fout == NULL)
//                     {
//                       pout() << endl << "Using max norm." << endl;
//                     }
//                     else
//                     {
//                       (*a_fout) << endl << "Using max norm." << endl;
//                     }
//                   }
//                 }
//               else
//                 {
//                   if (testverbosity > 1)
//                   {
//                     if (a_fout == NULL)
//                     {
//                       pout() << endl << "Using L-" << inorm << "norm." << endl;
//                     }
//                     else
//                     {
//                       (*a_fout) << endl << "Using L-" << inorm << "norm." << endl;
//                     }
//                   }
//                 }
//               Real coarnorm = EBArith::norm(a_errorCoar,
//                                             a_gridsCoar, a_ebislCoar,
//                                             refRat,
//                                             comp, inorm, normtype);
//               Real finenorm = EBArith::norm(a_errorFine,
//                                             a_gridsFine, a_ebislFine,
//                                             refRat,
//                                             comp, inorm, normtype);
//               if (testverbosity > 1)
//               {
//                 if (a_fout == NULL)
//                 {
//                   pout() << "Coarse Error Norm = " << coarnorm << endl;
//                 }
//                 else
//                 {
//                   (*a_fout) << "Coarse Error Norm = " << coarnorm << endl;
//                 }
//               }
//               if (testverbosity > 1)
//               {
//                 if (a_fout == NULL)
//                 {
//                   pout() << "Fine   Error Norm = " << finenorm << endl;
//                 }
//                 else
//                 {
//                   (*a_fout) << "Fine   Error Norm = " << finenorm << endl;
//                 }
//               }
//               if (itype == 0)
//                 {
//                   normsCoar[EBArith::orderScript(comp,inorm,ncomp)] = coarnorm;
//                   normsFine[EBArith::orderScript(comp,inorm,ncomp)] = finenorm;
//                 }
//               if ((Abs(finenorm) > 1.0e-10) && (Abs(coarnorm) > 1.0e-10))
//                 {
//                   Real order = log(Abs(coarnorm/finenorm))/log(2.0);
//                   //if (a_fout == NULL)
//                   //{
//                   //  pout() << "Order of scheme = " << order << endl;
//                   //}
//                   //else
//                   //{
//                   //  (*a_fout) << "Order of scheme = " << order << endl;
//                   //}
//                   if (itype == 0)
//                     {
//                       orders[EBArith::orderScript(comp,inorm,ncomp)] = order;
//                     }
//                 }
//             }
//         }
//       if (testverbosity > 1)
//       {
//         if (a_fout == NULL)
//         {
//           pout() << "==============================================" << endl ;;
//         }
//         else
//         {
//           (*a_fout) << "==============================================" << endl ;;
//         }
//       }
//     }

//   //output in latex format to be safe
//   int nfine = a_coarseDom.size(0);
//   nfine *= 2;
//   if (testverbosity > 0)
//     {
//       if (a_fout == NULL)
//       {
//         pout()    << setw(12)
//                   << setprecision(6)
//                   << setiosflags(ios::showpoint)
//                   << setiosflags(ios::scientific);
//       }
//       else
//       {
//         (*a_fout) << setw(12)
//                   << setprecision(6)
//                   << setiosflags(ios::showpoint)
//                   << setiosflags(ios::scientific);
//       }
//       for (int inorm = 0; inorm <= 2; inorm++)
//         {
//           if (a_fout == NULL)
//           {
//             pout() << "\\begin{table}[p]" << endl;
//             pout() << "\\begin{center}" << endl;
//             pout() << "\\begin{tabular}{|c|c|c|c|} \\hline" << endl;
//             pout() << "Variable & Coarse Error & Fine Error & Order\\\\" << endl;;
//             pout() << "\\hline \\hline " << endl;
//           }
//           else
//           {
//             (*a_fout) << "\\begin{table}[p]" << endl;
//             (*a_fout) << "\\begin{center}" << endl;
//             (*a_fout) << "\\begin{tabular}{|c|c|c|c|} \\hline" << endl;
//             (*a_fout) << "Variable & Coarse Error & Fine Error & Order\\\\" << endl;;
//             (*a_fout) << "\\hline \\hline " << endl;
//           }
//           for (int icomp = 0; icomp < ncomp; icomp++)
//             {
//               int iindex = EBArith::orderScript(icomp,inorm,ncomp);
//               if (a_fout == NULL)
//               {
//                 if (a_names.size() > icomp)
//                   pout() << a_names[icomp] << " &    \t ";
//                 else
//                   pout() << "var" << icomp << " &    \t ";

//                 pout()  << setw(12)
//                        << setprecision(6)
//                        << setiosflags(ios::showpoint)
//                        << setiosflags(ios::scientific)
//                        << normsCoar[iindex]  << " & "
//                        << setw(12)
//                        << setprecision(6)
//                        << setiosflags(ios::showpoint)
//                        << setiosflags(ios::scientific)
//                        << normsFine[iindex] << " & "
//                        << setw(12)
//                        << setprecision(6)
//                        << setiosflags(ios::showpoint)
//                        << setiosflags(ios::scientific)
//                        << orders[iindex];
//                 pout() << " \\\\ " << endl <<   "\\hline"  <<  endl;
//               }
//               else
//               {
//                 (*a_fout) << "var" << icomp << " &    \t "
//                           << setw(12)
//                           << setprecision(6)
//                           << setiosflags(ios::showpoint)
//                           << setiosflags(ios::scientific)
//                           << normsCoar[iindex]  << " & "
//                           << setw(12)
//                           << setprecision(6)
//                           << setiosflags(ios::showpoint)
//                           << setiosflags(ios::scientific)
//                           << normsFine[iindex] << " & "
//                           << setw(12)
//                           << setprecision(6)
//                           << setiosflags(ios::showpoint)
//                           << setiosflags(ios::scientific)
//                           << orders[iindex];
//                 (*a_fout) << " \\\\ " << endl <<   "\\hline"  <<  endl;
//               }
//             }
//           if (a_fout == NULL)
//           {
//             pout() << "\\end{tabular}" << endl;
//             pout() << "\\end{center}" << endl;
//             pout() << "\\caption{Solution error convergence rates using L-" << inorm << " norm. " << endl;
//             pout() << "$h_f = \\frac{1}{" << nfine << "}$ and $h_c = 2 h_f$, $D = " << SpaceDim << "$ }"  << endl;
//             pout() << "\\end{table}" << endl;
//             pout() << endl << endl;
//           }
//           else
//           {
//             (*a_fout) << "\\end{tabular}" << endl;
//             (*a_fout) << "\\end{center}" << endl;
//             (*a_fout) << "\\caption{Solution error convergence rates using L-" << inorm << " norm. " << endl;
//             (*a_fout) << "$h_f = \\frac{1}{" << nfine << "}$ and $h_c = 2 h_f$, $D = " << SpaceDim << "$ }"  << endl;
//             (*a_fout) << "\\end{table}" << endl;
//             (*a_fout) << endl << endl;
//           }
//         }
//     }
//   //latex output

//   delete normsFine;
//   delete normsCoar;
// }
void
EBArith::
compareError(Vector<Real>&                            a_orders,
             const Vector< LevelData<EBCellFAB>* >&   a_errorFine,
             const Vector< LevelData<EBCellFAB>* >&   a_errorCoar,
             const Vector< DisjointBoxLayout >&       a_gridsFine,
             const Vector< DisjointBoxLayout >&       a_gridsCoar,
             const Vector< EBISLayout >&              a_ebislFine,
             const Vector< EBISLayout >&              a_ebislCoar,
             const Vector<int>&                       a_refRat,
             const Box&                               a_coarseDom,
             int                                      a_testverbosity,
             fstream*                                 a_fout,
             Vector<string> a_names)
{
  CH_TIME("EBArith::compareError");

  const Vector<int>& refRat = a_refRat;
  const int ncomp = a_errorFine[0]->nComp();
  const int nnorm = 3;
  Real* normsCoar = new Real[ncomp*nnorm];
  Real* normsFine = new Real[ncomp*nnorm];
  a_orders.resize(ncomp*nnorm, 0.0);
  Real* orders    = &(a_orders[0]);
  for (int icomp = 0; icomp < ncomp; icomp++)
    {
      orders[icomp] = 0.0;
      for (int inorm = 0; inorm < nnorm; inorm++)
        {
          normsCoar[EBArith::orderScript(icomp, inorm, ncomp)] = 0;
          normsFine[EBArith::orderScript(icomp, inorm, ncomp)] = 0;
        }
    }
  int testverbosity = a_testverbosity;
  if (testverbosity > 1)
  {
    if (a_fout == NULL)
    {
      pout() << "==============================================" << endl;
    }
    else
    {
      *a_fout << "==============================================" << endl;
    }
  }
  for (int comp = 0; comp < ncomp; comp++)
    {
      if (testverbosity > 1)
      {
        if (a_fout == NULL)
        {
          if (a_names.size() > comp)
            pout() << "Comparing error in variable  " << a_names[comp] << endl;
          else
            pout() << "Comparing error in variable  " << comp << endl;
        }
        else
        {
          if (a_names.size() > comp)
            (*a_fout) << "Comparing error in variable  " << a_names[comp] << endl;
          else
            (*a_fout) << "Comparing error in variable  " << comp << endl;
        }
      }
      if (testverbosity > 1)
      {
        if (a_fout == NULL)
        {
          pout() << "==============================================" << endl;
        }
        else
        {
          (*a_fout) << "==============================================" << endl;
        }
      }
      for (int itype = 0; itype < 3; itype++)
        {
          EBNormType::NormMode normtype;
          if (itype == 0)
            {
              normtype = EBNormType::OverBoth;
              if (testverbosity > 1)
              {
                if (a_fout == NULL)
                {
                  pout() << endl << "Using all uncovered cells." << endl  ;
                }
                else
                {
                  (*a_fout) << endl << "Using all uncovered cells." << endl  ;
                }
              }
            }
          else if (itype == 1)
            {
              normtype = EBNormType::OverOnlyRegular;
              if (testverbosity > 1)
              {
                if (a_fout == NULL)
                {
                  pout() << endl << "Using only regular cells." << endl ;
                }
                else
                {
                  (*a_fout) << endl << "Using only regular cells." << endl ;
                }
              }
            }
          else
            {
              normtype = EBNormType::OverOnlyIrregular;
              if (testverbosity > 1)
              {
                if (a_fout == NULL)
                {
                  pout() << endl << "Using only irregular cells." << endl;
                }
                else
                {
                  (*a_fout) << endl << "Using only irregular cells." << endl;
                }
              }
            }

          for (int inorm = 0; inorm <= 2; inorm++)
            {

              if (inorm == 0)
                {
                  if (testverbosity > 1)
                  {
                    if (a_fout == NULL)
                    {
                      pout() << endl << "Using max norm." << endl;
                    }
                    else
                    {
                      (*a_fout) << endl << "Using max norm." << endl;
                    }
                  }
                }
              else
                {
                  if (testverbosity > 1)
                  {
                    if (a_fout == NULL)
                    {
                      pout() << endl << "Using L-" << inorm << "norm." << endl;
                    }
                    else
                    {
                      (*a_fout) << endl << "Using L-" << inorm << "norm." << endl;
                    }
                  }
                }
              Real coarnorm = EBArith::norm(a_errorCoar,
                                            a_gridsCoar, a_ebislCoar,
                                            refRat,
                                            comp, inorm, normtype);
              Real finenorm = EBArith::norm(a_errorFine,
                                            a_gridsFine, a_ebislFine,
                                            refRat,
                                            comp, inorm, normtype);
              if (testverbosity > 1)
              {
                if (a_fout == NULL)
                {
                  pout() << "Coarse Error Norm = " << coarnorm << endl;
                }
                else
                {
                  (*a_fout) << "Coarse Error Norm = " << coarnorm << endl;
                }
              }
              if (testverbosity > 1)
              {
                if (a_fout == NULL)
                {
                  pout() << "Fine   Error Norm = " << finenorm << endl;
                }
                else
                {
                  (*a_fout) << "Fine   Error Norm = " << finenorm << endl;
                }
              }
              if (itype == 0)
                {
                  normsCoar[EBArith::orderScript(comp,inorm,ncomp)] = coarnorm;
                  normsFine[EBArith::orderScript(comp,inorm,ncomp)] = finenorm;
                }
              if ((Abs(finenorm) > 1.0e-10) && (Abs(coarnorm) > 1.0e-10))
                {
                  Real order = log(Abs(coarnorm/finenorm))/log(2.0);
                  //if (a_fout == NULL)
                  //{
                  //  pout() << "Order of scheme = " << order << endl;
                  //}
                  //else
                  //{
                  //  (*a_fout) << "Order of scheme = " << order << endl;
                  //}
                  if (itype == 0)
                    {
                      orders[EBArith::orderScript(comp,inorm,ncomp)] = order;
                    }
                }
            }
        }
      if (testverbosity > 1)
      {
        if (a_fout == NULL)
        {
          pout() << "==============================================" << endl ;;
        }
        else
        {
          (*a_fout) << "==============================================" << endl ;;
        }
      }
    }

  //output in latex format to be safe
  int ncoar = a_coarseDom.size(0);
  int nmedi = 2*ncoar;
  int nfine = 2*nmedi;
  if (testverbosity > 0)
    {
      if (a_fout == NULL)
      {
        pout()    << setw(12)
                  << setprecision(6)
                  << setiosflags(ios::showpoint)
                  << setiosflags(ios::scientific);
      }
      else
      {
        (*a_fout) << setw(12)
                  << setprecision(6)
                  << setiosflags(ios::showpoint)
                  << setiosflags(ios::scientific);
      }
      for (int inorm = 0; inorm <= 2; inorm++)
        {
          if (a_fout == NULL)
          {
            pout() << "\\begin{table}[p]" << endl;
            pout() << "\\begin{center}" << endl;
            pout() << "\\begin{tabular}{|cccc|} \\hline" << endl;
            pout() << "Variable & $e_{4h \\rightarrow 2h}$ & Order & $e_{2h \\rightarrow h}$\\\\" << endl;;
            pout() << "\\hline " << endl;
          }
          else
          {
            (*a_fout) << "\\begin{table}[p]" << endl;
            (*a_fout) << "\\begin{center}" << endl;
            (*a_fout) << "\\begin{tabular}{|cccc|} \\hline" << endl;
           pout() << "Variable & $e_{\\frac{1}{" << ncoar << "}\\rightarrow\\frac{1}{" << nmedi << "}}$ & Order & $e_{\\frac{1}{" << nmedi << "}\\rightarrow\\frac{1}{" << nfine << "}}$\\\\" << endl;;
            (*a_fout) << "\\hline " << endl;
          }
          for (int icomp = 0; icomp < ncomp; icomp++)
            {
              int iindex = EBArith::orderScript(icomp,inorm,ncomp);
              if (a_fout == NULL)
              {
                if (a_names.size() > icomp)
                  pout() << setw(14) << a_names[icomp] << " &    \t ";
                else
                  pout() << "var" << icomp << " &    \t ";

                pout() << setw(12)
                       << setprecision(6)
                       << setiosflags(ios::showpoint)
                       << setiosflags(ios::scientific)
                       << normsCoar[iindex]  << " & "
                       << setw(12)
                       << setprecision(3)
                       << setiosflags(ios::showpoint)
                       << setiosflags(ios::scientific)
                       << orders[iindex] << " & "
                       << setw(12)
                       << setprecision(6)
                       << setiosflags(ios::showpoint)
                       << setiosflags(ios::scientific)
                       << normsFine[iindex];
                pout() << " \\\\ " << endl;
              }
              else
              {
                (*a_fout) << "var" << icomp << " &    \t "
                          << setw(12)
                          << setprecision(6)
                          << setiosflags(ios::showpoint)
                          << setiosflags(ios::scientific)
                          << normsCoar[iindex]  << " & "
                          << setw(12)
                          << setprecision(3)
                          << setiosflags(ios::showpoint)
                          << orders[iindex] << " & "
                          << setw(12)
                          << setprecision(6)
                          << setiosflags(ios::showpoint)
                          << setiosflags(ios::scientific)
                          << normsFine[iindex];
                (*a_fout) << " \\\\ " << endl;
              }
            }
          if (a_fout == NULL)
          {
            pout() << "\\hline " << endl;
            pout() << "\\end{tabular}" << endl;
            pout() << "\\end{center}" << endl;
            pout() << "\\caption{Solution error convergence rates using $L_" << inorm << "$-norm for $h=\\frac{1}{" << nfine << "}$.} " << endl;
            pout() << "\\end{table}" << endl;
            pout() << endl << endl;
          }
          else
          {
            (*a_fout) << "\\end{tabular}" << endl;
            (*a_fout) << "\\end{center}" << endl;
            (*a_fout) << "\\caption{Solution error convergence rates using $L_" << inorm << "$-norm for $h=\\frac{1}{" << nfine << "}$.} " << endl;
            (*a_fout) << "\\end{table}" << endl;
            (*a_fout) << endl << endl;
          }
        }
    }
  //latex output

  delete normsFine;
  delete normsCoar;
}

void
EBArith::
shrinkIVS(IntVectSet& a_ivs, const int& a_numShrink)
{
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (SideIterator  sit; sit.ok(); ++sit)
        {
          IntVect ivshift = IntVect::Zero;
          ivshift[idir] = sign(sit())*a_numShrink;
          IntVectSet shiftedIVS = a_ivs;
          shiftedIVS.shift(ivshift);
          a_ivs &= shiftedIVS;
        }
    }
}


void
EBArith::
timeInterpolate(LevelData<EBCellFAB>&       a_U,
                const LevelData<EBCellFAB>& a_UOld,
                const LevelData<EBCellFAB>& a_UNew,
                const DisjointBoxLayout&    a_grids,
                const Real&                 a_time,
                const Real&                 a_told,
                const Real&                 a_tnew)
{
  CH_assert(a_time >= a_told);
  CH_assert(a_time <= a_tnew);
  CH_assert(a_told < (a_tnew-1.0e-12));
  Real newFac = (a_time-a_told)/(a_tnew-a_told);
  Real oldFac = 1.0-newFac;
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      a_U[dit()].axby(a_UOld[dit()], a_UNew[dit()], oldFac, newFac);
    }
}
Real
EBArith::
getDiagWeight(  VoFStencil&     a_vofStencil,
                const VolIndex& a_vof,
                int             a_ivar)
{
  //has to be zero because adding to this later
  Real retval = 0;
  bool found = false;
  for (int ivof = 0; ivof  < a_vofStencil.size(); ivof++)
    {
      if ((a_vofStencil.vof(ivof) == a_vof) && (a_vofStencil.variable(ivof) == a_ivar))
        {
          found = true;
          //additive in case there are more than one entry with vof == a_vof
          retval += a_vofStencil.weight(ivof);
        }
    }
  if (!found)
    {
      //      MayDay::Warning("no diagonal weight, probably an empty cell");
      retval = 1;
    }
  return retval;
}

void
EBArith::
getMultiColors(Vector<IntVect>& a_colors)
{

#if CH_SPACEDIM==2
  a_colors.resize(4);
  a_colors[0] = IntVect::Zero;//(0,0)
  a_colors[1] = IntVect::Unit;//(1,1)
  a_colors[2] = IntVect::Zero + BASISV(1);//(0,1)
  a_colors[3] = IntVect::Zero + BASISV(0);//(1,0)
#elif CH_SPACEDIM==3
  a_colors.resize(8);
  a_colors[0] = IntVect::Zero;//(0,0,0)
  a_colors[1] = IntVect::Zero + BASISV(0) + BASISV(1);//(1,1,0)
  a_colors[2] = IntVect::Zero + BASISV(1) + BASISV(2);//(0,1,1)
  a_colors[3] = IntVect::Zero + BASISV(0) + BASISV(2);//(1,0,1)
  a_colors[4] = IntVect::Zero + BASISV(1);//(0,1,0)
  a_colors[5] = IntVect::Zero + BASISV(0);//(1,0,0)
  a_colors[6] = IntVect::Zero + BASISV(2);//(0,0,1)
  a_colors[7] = IntVect::Unit;//(1,1,1)
#endif
}
void
EBArith::
interpolateCFH(EBCellFAB&                  a_phi,
               const int  &                a_idir,
               const Side::LoHiSide&       a_hiorlo,
               const EBISBox&              a_ebisBox,
               const Real&                 a_dxfine,
               const Real&                 a_dxcoar,
               const IntVectSet&           a_interpIVS)
{
  Real halfdxcoar = a_dxcoar/2.0;
  Real halfdxfine = a_dxfine/2.0;
  Real xg = halfdxcoar -   halfdxfine;
  Real xc = halfdxcoar +   halfdxfine;
  Real xf = halfdxcoar + 3*halfdxfine;
  Real hf = a_dxfine;
  Real denom = xf*xc*hf;
  for (int ivar = 0; ivar < a_phi.nComp(); ivar++)
    {
      for (VoFIterator vofit(a_interpIVS, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& VoFGhost = vofit();

          IntVect ivGhost  = VoFGhost.gridIndex();
          IntVect ivClose =  ivGhost;
          IntVect ivFar   =  ivGhost;

          Vector<VolIndex> farVoFs;
          Vector<VolIndex> closeVoFs = a_ebisBox.getVoFs(VoFGhost,
                                                         a_idir,
                                                         flip(a_hiorlo),
                                                         1);
          bool hasClose = (closeVoFs.size() > 0);
          bool hasFar = false;
          Real phic = 0.0;
          Real phif = 0.0;
          if (hasClose)
            {
              const int& numClose = closeVoFs.size();
              for (int iVof=0;iVof<numClose;iVof++)
                {
                  const VolIndex& vofClose = closeVoFs[iVof];
                  phic += a_phi(vofClose,0);
                }
              phic /= Real(numClose);

              farVoFs = a_ebisBox.getVoFs(VoFGhost,
                                          a_idir,
                                          flip(a_hiorlo),
                                          2);
              hasFar   = (farVoFs.size()   > 0);
              if (hasFar)
                {
                  const int& numFar = farVoFs.size();
                  for (int iVof=0;iVof<numFar;iVof++)
                    {
                      const VolIndex& vofFar = farVoFs[iVof];
                      phif += a_phi(vofFar,0);
                    }
                  phif /= Real(numFar);
                }
            }

          Real phiGhost;
          if (hasClose && hasFar)
            {
              // quadratic interpolation  phi = ax^2 + bx + c
              Real A = (phif*xc - phic*xf)/denom;
              Real B = (phic*hf*xf - phif*xc*xc + phic*xf*xc)/denom;

              phiGhost = A*xg*xg + B*xg;
            }
          else if (hasClose)
            {
              //linear interpolation
              Real slope =  phic/xc;
              phiGhost   =  slope*xg;
            }
          else
            {
              phiGhost = 0.0; //nothing to interpolate from
            }
          a_phi(VoFGhost, ivar) = phiGhost;
        }
    }
}
void
EBArith::
getDir1Dir2(int& a_dir1, int& a_dir2, const int& a_idir)
{
#if CH_SPACEDIM==2
  a_dir1 = 0;
  a_dir2 = 1;
#elif CH_SPACEDIM==3
  if (a_idir == 0)
    {
      a_dir1 = 1;
      a_dir2 = 2;
    }
  else if (a_idir == 1)
    {
      a_dir1 = 0;
      a_dir2 = 2;
    }
  else if (a_idir == 2)
    {
      a_dir1 = 0;
      a_dir2 = 1;
    }
  else
    {
      MayDay::Error("bogus idir");
    }
#else
  bogus_spacedim();
#endif
}

int
EBArith::
getExtrapolationStencil(VoFStencil&     a_stencil,
                        const RealVect& a_dist,
                        const RealVect& a_dx,
                        const VolIndex& a_startVoF,
                        const EBISBox&  a_ebisBox,
                        int a_noExtrapThisDir,
                        IntVectSet*    a_cfivsPtr,
                        int ivar)
{
  int order = 2;
  a_stencil.clear();

  //zeroth order Taylor series
  a_stencil.add(a_startVoF,1.0,ivar);

  //do taylor series extrapolation stencil
  //get stencils for derivatives
  int derivorder;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (idir != a_noExtrapThisDir)
        {
          VoFStencil firstDSten;
          derivorder = EBArith::getFirstDerivStencil(firstDSten,
                                                     a_startVoF, a_ebisBox,
                                                     idir, a_dx[idir], a_cfivsPtr, ivar);
          order = Min(order, derivorder);

          firstDSten *= a_dist[idir];
          a_stencil += firstDSten;

          VoFStencil secondDSten;
          derivorder = EBArith::getSecondDerivStencil(secondDSten,
                                                      a_startVoF, a_ebisBox,
                                                      idir, a_dx[idir], a_cfivsPtr, ivar);
          order = Min(order, derivorder);
          secondDSten *= (0.5*a_dist[idir]*a_dist[idir]);
          a_stencil += secondDSten;
        }
#if CH_SPACEDIM==3
      int dir1, dir2;
      EBArith::getDir1Dir2(dir1, dir2, idir);
      if ((dir1 != a_noExtrapThisDir) && (dir2 != a_noExtrapThisDir))
        {
          VoFStencil mixedDSten;
          derivorder = EBArith::getMixedDerivStencil(mixedDSten,
                                                     a_startVoF, a_ebisBox,
                                                     dir1, dir2, a_dx[dir1], a_dx[dir2], a_cfivsPtr, ivar);
          order = Min(order, derivorder);
          mixedDSten *= (a_dist[dir1]*a_dist[dir2]);
          a_stencil += mixedDSten;
        }
#endif
    }

#if CH_SPACEDIM==2
  int dir1 = 0;
  int dir2 = 1;
  if ((dir1 != a_noExtrapThisDir) && (dir2 != a_noExtrapThisDir))
    {
      VoFStencil mixedDSten;
      derivorder = EBArith::getMixedDerivStencil(mixedDSten,
                                                 a_startVoF, a_ebisBox,
                                                 dir1, dir2, a_dx[dir1], a_dx[dir2], a_cfivsPtr, ivar);
      order = Min(order, derivorder);
      mixedDSten *= (a_dist[dir1]*a_dist[dir2]);
      a_stencil += mixedDSten;
    }
#endif

  return order;
}

int
EBArith::
getFirstOrderExtrapolationStencil(VoFStencil&     a_stencil,
                                  const RealVect& a_dist,
                                  const RealVect& a_dx,
                                  const VolIndex& a_startVoF,
                                  const EBISBox&  a_ebisBox,
                                  int a_noExtrapThisDir,
                                  IntVectSet*    a_cfivsPtr,
                                  int ivar)
{
  int order = 2;
  a_stencil.clear();

  //zeroth order Taylor series
  a_stencil.add(a_startVoF,1.0,ivar);

  //do taylor series extrapolation stencil
  //get stencils for derivatives
  int derivorder;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (idir != a_noExtrapThisDir)
        {
          VoFStencil firstDSten;
          derivorder = EBArith::getFirstDerivStencilWidthOne(firstDSten,
                                                     a_startVoF, a_ebisBox,
                                                     idir, a_dx[idir], a_cfivsPtr, ivar);
          order = Min(order, derivorder);

          firstDSten *= a_dist[idir];
          a_stencil += firstDSten;

        }
    }
  return order;
}

// int
// EBArith::
// get1stOrderExtrapolationStencil(VoFStencil&     a_stencil,
//                                 const RealVect& a_dist,
//                                 const RealVect& a_dx,
//                                 const VolIndex& a_startVoF,
//                                 const EBISBox&  a_ebisBox,
//                                 int a_noExtrapThisDir,
//                                 IntVectSet*    a_cfivsPtr,
//                                 int ivar)
// {
//   int order = 2;
//   a_stencil.clear();

//   //zeroth order Taylor series
//   a_stencil.add(a_startVoF,1.0,ivar);

//   //do taylor series extrapolation stencil
//   //get stencils for derivatives
//   int derivorder;
//   for (int idir = 0; idir < SpaceDim; idir++)
//     {
//       if (idir != a_noExtrapThisDir)
//         {
//           VoFStencil firstDSten;
//           derivorder = EBArith::getFirstDerivStencil(firstDSten,
//                                                      a_startVoF, a_ebisBox,
//                                                      idir, a_dx[idir], a_cfivsPtr, ivar);
//           order = Min(order, derivorder);

//           firstDSten *= a_dist[idir];
//           a_stencil += firstDSten;

//         }
//     }

//   return order;
// }

/****/
void
EBArith::
loHiCenter(Box&                 a_loBox,
           int&                 a_hasLo,
           Box&                 a_hiBox,
           int&                 a_hasHi,
           Box&                 a_centerBox,
           const ProblemDomain& a_eblg,
           const Box&           a_inBox,
           const int&           a_dir,
           IntVectSet*          a_cfivsPtr)
{
  bool checkCFIVS = false;
  if (a_cfivsPtr != NULL)
    {
      checkCFIVS = true;
    }

  CH_TIME("EBArith::loHiCenter");
  // Make a copy of the input box which can be modified
  Box inBox = a_inBox;

  inBox &= a_eblg;

  a_centerBox = inBox;
  a_centerBox.grow(a_dir, 1);
  a_centerBox &= a_eblg;
  a_centerBox.grow(a_dir,-1);

  // See if this chops off the high side of the input box
  Box tmp = inBox;
  tmp.shift(a_dir,1);
  tmp &= a_eblg;
  tmp.shift(a_dir,-1);

  // If so, set up the high, one-sided difference box, a_hiBox
  if (tmp != inBox)
    {
      a_hasHi = 1;
      a_hiBox = adjCellHi(tmp,a_dir);
    }
  else if (checkCFIVS)
    {
      Box sideBox = adjCellHi(a_centerBox, a_dir);
      IntVectSet ivsInter(sideBox);
      ivsInter &= (*a_cfivsPtr);
      if (ivsInter.isEmpty())
        {
          a_hasHi = 0;
        }
      else
        {
          a_hasHi = 1;
          a_hiBox = sideBox;
          a_hiBox.shift(a_dir, -1);
          a_centerBox.growDir(a_dir,Side::Hi,-1);
        }
    }
  else
    {
      a_hasHi = 0;
    }

  // See if this chops off the low side of the input box
  tmp = inBox;
  tmp.shift(a_dir,-1);
  tmp &= a_eblg;
  tmp.shift(a_dir,1);

  // If so, set up the low, one-sided difference box, a_loBox
  if (tmp != inBox)
    {
      a_hasLo = 1;
      a_loBox = adjCellLo(tmp,a_dir);
    }
  else if (checkCFIVS)
    {
      Box sideBox = adjCellLo(a_centerBox, a_dir);
      IntVectSet ivsInter(sideBox);
      ivsInter &= (*a_cfivsPtr);
      if (ivsInter.isEmpty())
        {
          a_hasLo = 0;
        }
      else
        {
          a_hasLo = 1;
          a_loBox = sideBox;
          a_loBox.shift(a_dir, 1);
          a_centerBox.growDir(a_dir,Side::Lo,-1);
        }
    }
  else
    {
      a_hasLo = 0;
    }
}
/***/
/// Can be used instead of loHiCenter when the center box isn't needed
/**/
void
EBArith::
loHi(Box&                 a_loBox,
     int&                 a_hasLo,
     Box&                 a_hiBox,
     int&                 a_hasHi,
     const ProblemDomain& a_eblg,
     const Box&           a_inBox,
     const int&           a_dir)
{
  CH_TIME("EBArith::loHi");
  // Make a copy of the input box which can be modified
  Box inBox = a_inBox;

  inBox &= a_eblg;

  // See if this chops off the high side of the input box
  Box tmp = inBox;
  tmp.shift(a_dir,1);
  tmp &= a_eblg;
  tmp.shift(a_dir,-1);

  // If so, set up the high, one-sided difference box, a_hiBox
  if (tmp != inBox)
    {
      a_hasHi = 1;
      a_hiBox = adjCellHi(tmp,a_dir);
    }
  else
    {
      a_hasHi = 0;
    }

  // See if this chops off the low side of the input box
  tmp = inBox;
  tmp.shift(a_dir,-1);
  tmp &= a_eblg;
  tmp.shift(a_dir,1);

  // If so, set up the low, one-sided difference box, a_loBox
  if (tmp != inBox)
    {
      a_hasLo = 1;
      a_loBox = adjCellLo(tmp,a_dir);
    }
  else
    {
      a_hasLo = 0;
    }
}

//version that does not fill ebislCoar
bool EBArith::
getCoarserLayouts(DisjointBoxLayout&       a_dblCoar,
                  ProblemDomain&           a_domainCoar,
                  const DisjointBoxLayout& a_dblFine,
                  const ProblemDomain&     a_domainFine,
                  int                      a_refToCoar,
                  int                      a_maxBoxSize,
                  bool&                    a_layoutChanged,
                  int                      a_testRef)
{
  CH_TIME("EBArith::getCoarserLayouts(no ebisl)");
  // check to see if domain is coarsenable by this amount
  //never want to coarsen down to 1x1
  int testRef = a_testRef*a_refToCoar;
  ProblemDomain testBox = coarsen(a_domainFine, testRef);
  testBox.refine(testRef);
  a_layoutChanged = true;
  if (testBox != a_domainFine)
    {
      //not coarsenable.
      return false;
    }

  a_domainCoar = coarsen(a_domainFine, a_refToCoar);

  int fac = 2;

  if ((a_dblFine.coarsenable(fac*a_refToCoar) && a_dblFine.isClosed() && (a_dblFine.size() > 0)))
    {
      //if we can, just coarsen the grids
      coarsen(a_dblCoar, a_dblFine, a_refToCoar);
      a_layoutChanged = false;
    }
  else if (a_dblFine.size() != 1)
    {//check to see if we are covering all uncovered IntVects with multiple boxes
      unsigned long long numPtsLeft = a_domainFine.domainBox().numPts();
      for (LayoutIterator lit = a_dblFine.layoutIterator(); lit.ok(); ++lit)
        {
          unsigned long long  ptsGrid = a_dblFine[lit()].numPts();
          numPtsLeft -= ptsGrid;
        }
      if (numPtsLeft == 0)
         {
          //if we are covering domain with more than one box, make one big split-up box
          Vector<Box> boxes;
          domainSplit(a_domainCoar, boxes, a_maxBoxSize);
          mortonOrdering(boxes);

          Vector<int> procs;
          LoadBalance(procs, boxes);
          a_dblCoar.define(boxes, procs, a_domainCoar);
        }
      else
        {
          //we are out of ideas
          return false;
        }
    }
  else
    {
      //we are out of ideas
      return false;
    }

  //if we got here, then we have coarser stuff
  return true;
}

void
EBArith::
getVoFsDir(bool& a_hasClose, VolIndex& a_closeVoF,
           bool& a_hasFar,   VolIndex& a_farVoF,
           const EBISBox& a_ebisBox,
           const VolIndex& a_vof,
           int a_idir, Side::LoHiSide a_sd,
           IntVectSet*    a_cfivsPtr)
{
  a_hasClose = false;
  a_hasFar   = false;
  bool checkCFIVS = false;
  if (a_cfivsPtr != NULL)
    {
      checkCFIVS = true;
    }
  //get faces on both sides to see in which direction we can do diffs
  Vector<FaceIndex> closeFaces = a_ebisBox.getFaces(a_vof, a_idir, a_sd);

  //boundary faces and multi-valued faces are to be one-sided away from
  a_hasClose = ((closeFaces.size() == 1) && (!closeFaces[0].isBoundary()));
  if (a_hasClose)
    {
      a_closeVoF = closeFaces[0].getVoF(a_sd);
      if (checkCFIVS && (*a_cfivsPtr).contains(a_closeVoF.gridIndex()))
        {
          a_hasClose = false;
        }
      Vector<FaceIndex> farFaces = a_ebisBox.getFaces(a_closeVoF, a_idir, a_sd);
      a_hasFar = ((farFaces.size() == 1) && (!farFaces[0].isBoundary()));
      if (a_hasFar)
        {
          a_farVoF =  farFaces[0].getVoF(a_sd);
          if (checkCFIVS && (*a_cfivsPtr).contains(a_farVoF.gridIndex()))
            {
              a_hasFar = false;
            }
        }
    }
}

///
/**
   Gets the stencil to take the first derivative of  cell centered data.
*/
int
EBArith::
getFirstDerivStencil(VoFStencil&      a_sten,
                     const VolIndex&  a_vof,
                     const EBISBox&   a_ebisBox,
                     const int&       a_idir,
                     const Real&      a_dx,
                     IntVectSet*    a_cfivsPtr,
                     int ivar)

{
  CH_assert(a_dx > 0.0);
  CH_assert(a_idir >= 0);
  CH_assert(a_idir <  SpaceDim);
  int order;
  bool hasLo, hasLower, hasHi, hasHigher;
  VolIndex loVoF, lowerVoF, hiVoF, higherVoF;
  EBArith::getVoFsDir(hasLo, loVoF, hasLower, lowerVoF,   a_ebisBox, a_vof, a_idir, Side::Lo, a_cfivsPtr);
  EBArith::getVoFsDir(hasHi, hiVoF, hasHigher, higherVoF, a_ebisBox, a_vof, a_idir, Side::Hi, a_cfivsPtr);

  //clear any residual stencil info
  a_sten.clear();

  if (hasHi && hasLo)
    {
      //if we have vofs on both sides, we do our friend the centered difference
      order = 2;
      a_sten.add(hiVoF,  1.0, ivar);
      a_sten.add(loVoF, -1.0, ivar);
      a_sten *= (1.0/(2.*a_dx));
    }
  else if (hasHi)
    {
      //if we only have vofs on the high side, we see if we can take a higher
      //order one sided diff.   If not, we take simple one-sided diff and drop order
      if (hasHigher)
        {
          a_sten.add(hiVoF,       4.0, ivar);
          a_sten.add(higherVoF,  -1.0, ivar);
          a_sten.add(a_vof,      -3.0, ivar);
          a_sten *= (1.0/(2.*a_dx));
          order = 2;
        }
      else
        {
          a_sten.add(hiVoF,       1.0, ivar);
          a_sten.add(a_vof,      -1.0, ivar);
          a_sten *= (1.0/(a_dx));
          order = 1;
        }
    }
  else if (hasLo)
    {
      //if we only have vofs on the low side, we see if we can take a higher
      //order one sided diff.   If not, we take simple one-sided diff and drop order
      if (hasLower)
        {
          a_sten.add(loVoF,      -4.0, ivar);
          a_sten.add(lowerVoF,    1.0, ivar);
          a_sten.add(a_vof,       3.0, ivar);
          a_sten *= (1.0/(2.*a_dx));
          order = 2;
        }
      else
        {
          a_sten.add(a_vof,       1.0, ivar);
          a_sten.add(loVoF,      -1.0, ivar);
          a_sten *= (1.0/(a_dx));
          order = 1;
        }
    }
  else
    {
      //no vofs on either side.   return cleared stencil and order=0
      order = 0;
    }

  return order;
}


///
/**
   Gets the stencil to take the first derivative of  cell centered data.
   no high order one sided stuff (to keep stencil close)
*/
int
EBArith::
getFirstDerivStencilWidthOne(VoFStencil&      a_sten,
                             const VolIndex&  a_vof,
                             const EBISBox&   a_ebisBox,
                             const int&       a_idir,
                             const Real&      a_dx,
                             IntVectSet*    a_cfivsPtr,
                             int ivar)

{
  CH_assert(a_dx > 0.0);
  CH_assert(a_idir >= 0);
  CH_assert(a_idir <  SpaceDim);
  int order;
  bool hasLo, hasLower, hasHi, hasHigher;
  VolIndex loVoF, lowerVoF, hiVoF, higherVoF;
  EBArith::getVoFsDir(hasLo, loVoF, hasLower, lowerVoF,   a_ebisBox, a_vof, a_idir, Side::Lo, a_cfivsPtr);
  EBArith::getVoFsDir(hasHi, hiVoF, hasHigher, higherVoF, a_ebisBox, a_vof, a_idir, Side::Hi, a_cfivsPtr);

  //clear any residual stencil info
  a_sten.clear();

  if (hasHi && hasLo)
    {
      //if we have vofs on both sides, we do our friend the centered difference
      order = 2;
      a_sten.add(hiVoF,  1.0, ivar);
      a_sten.add(loVoF, -1.0, ivar);
      a_sten *= (1.0/(2.*a_dx));
    }
  else if (hasHi)
    {
      a_sten.add(hiVoF,       1.0, ivar);
      a_sten.add(a_vof,      -1.0, ivar);
      a_sten *= (1.0/(a_dx));
      order = 1;
    }
  else if (hasLo)
    {
      a_sten.add(a_vof,       1.0, ivar);
      a_sten.add(loVoF,      -1.0, ivar);
      a_sten *= (1.0/(a_dx));
      order = 1;
    }
  else
    {
      //no vofs on either side.   return cleared stencil and order=0
      order = 0;
    }

  return order;
}


///
/**
   Gets the stencil to take the first derivative of  cell centered data.
*/
int
EBArith::
getSecondDerivStencil(VoFStencil&      a_sten,
                      const VolIndex&  a_vof,
                      const EBISBox&   a_ebisBox,
                      const int&       a_idir,
                      const Real&      a_dx,
                      IntVectSet*    a_cfivsPtr,
                      int ivar)
{
  CH_assert(a_dx > 0.0);
  CH_assert(a_idir >= 0);
  CH_assert(a_idir <  SpaceDim);
  bool hasLo, hasLower, hasHi, hasHigher;
  VolIndex loVoF, lowerVoF, hiVoF, higherVoF;
  EBArith::getVoFsDir(hasLo, loVoF, hasLower, lowerVoF,   a_ebisBox, a_vof, a_idir, Side::Lo, a_cfivsPtr);
  EBArith::getVoFsDir(hasHi, hiVoF, hasHigher, higherVoF, a_ebisBox, a_vof, a_idir, Side::Hi, a_cfivsPtr);
  int order;

  //clear any residual stencil info
  a_sten.clear();


  if (hasHi && hasLo)
    {
      //if we have vofs on both sides, we do our friend the centered difference
      order = 2;
      a_sten.add(hiVoF,  1.0, ivar);
      a_sten.add(loVoF,  1.0, ivar);
      a_sten.add(a_vof, -2.0, ivar);
      a_sten *= (1.0/(a_dx*a_dx));
    }
  else if (hasHi)
    {
      //no vof on the low side
      //if we can, shift stencil one to the high side
      if (hasHigher)
        {
          a_sten.add(higherVoF,   1.0, ivar);
          a_sten.add(a_vof,       1.0, ivar);
          a_sten.add(hiVoF,      -2.0, ivar);
          a_sten *= (1.0/(a_dx*a_dx));
          order = 1;
        }
      else
        {
          //need 3 points
          order = 0;
        }
    }
  else if (hasLo)
    {
      //no vof on the high side
      //if we can, shift stencil one to the low side
      if (hasLower)
        {
          a_sten.add(lowerVoF,    1.0, ivar);
          a_sten.add(a_vof,       1.0, ivar);
          a_sten.add(loVoF,      -2.0, ivar);
          a_sten *= (1.0/(a_dx*a_dx));
          order = 1;
        }
      else
        {
          //need 3 points
          order = 0;
        }
    }
  else
    {
      //no vofs on either side.   return cleared stencil and order=0
      order = 0;
    }

  return order;
}
/****/
///
/**
   Gets the stencil to take a mixed derivative of  cell centered data.
*/
int
EBArith::
getMixedDerivStencil(VoFStencil&      a_sten,
                     const VolIndex&  a_vof,
                     const EBISBox&   a_ebisBox,
                     const int&       a_dir1,
                     const int&       a_dir2,
                     const Real&      a_dx1,
                     const Real&      a_dx2,
                     IntVectSet*    a_cfivsPtr,
                     int ivar)
{
  CH_assert(a_dx1 > 0.0);
  CH_assert(a_dx2 > 0.0);
  CH_assert(a_dir1 >= 0);
  CH_assert(a_dir2 >= 0);
  CH_assert(a_dir1 <  SpaceDim);
  CH_assert(a_dir2 <  SpaceDim);
  bool checkCFIVS = false;
  if (a_cfivsPtr != NULL)
    {
      checkCFIVS = true;
    }
  int radius = 1;
  Vector<VolIndex> vofList;

  EBArith::getAllVoFsInMonotonePath(vofList, a_vof, a_ebisBox, radius);


  //clear any residual stencil info
  a_sten.clear();

  //see how many corners have all vofs available for a mixed stencil derivative
  //and average the available stencils
  int numSten = 0;
  IntVect iv = a_vof.gridIndex();
  IntVect ivOne, ivTwo, ivCor;
  VolIndex vofOne, vofTwo, vofCor;
  for (SideIterator sit1; sit1.ok(); ++sit1)
    {
      int sign1 = sign(sit1());
      ivOne = iv + sign1*BASISV(a_dir1);
      bool isOneHere = EBArith::isVoFHere(vofOne, vofList, ivOne);
      if (isOneHere)
        {
          for (SideIterator sit2; sit2.ok(); ++sit2)
            {
              int sign2 = sign(sit2());
              ivTwo = iv + sign2*BASISV(a_dir2);
              bool isTwoHere = EBArith::isVoFHere(vofTwo, vofList, ivTwo);
              if (isTwoHere)
                {
                  ivCor = iv + sign2*BASISV(a_dir2) +  sign1*BASISV(a_dir1);
                  bool isCorHere = EBArith::isVoFHere(vofCor, vofList, ivCor);

                  if (isCorHere && checkCFIVS)
                    {
                      if ( (*a_cfivsPtr).contains(vofCor.gridIndex()) ||
                          (*a_cfivsPtr).contains(vofOne.gridIndex()) ||
                          (*a_cfivsPtr).contains(vofTwo.gridIndex()))
                        {
                          isCorHere = false;
                        }
                    }
                  if (isCorHere)
                    {
                      //if we have all vofs, construct the stencil.
                      Real rsign = Real(sign1*sign2);

                      VoFStencil mixedStencil;
                      mixedStencil.add(a_vof,   1.0, ivar);
                      mixedStencil.add(vofCor,  1.0, ivar);
                      mixedStencil.add(vofOne, -1.0, ivar);
                      mixedStencil.add(vofTwo, -1.0, ivar);

                      mixedStencil *= rsign/(a_dx1*a_dx2);

                      a_sten += mixedStencil;
                      numSten++;
                    }
                }
            }
        }
    }
  int order = 0;
  if (numSten > 0)
    {
      a_sten *= 1.0/numSten;

      if (numSten == 4)
        {
          order = 2;
        }
      else
        {
          order = 1;
        }
    }

  return order;
}
/****/


void
EBArith::defineCFIVS(LayoutData<IntVectSet>&   a_cfivs,
                     const DisjointBoxLayout&  a_grids,
                     const ProblemDomain&      a_probDom)
{
  a_cfivs.define(a_grids);
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      Box grownBox = grow(a_grids.get(dit()), 1);
      grownBox &= a_probDom;
      a_cfivs[dit()] = IntVectSet(grownBox);
      for (LayoutIterator lit = a_grids.layoutIterator(); lit.ok(); ++lit)
        {
          a_cfivs[dit()] -= a_grids[lit()];
        }
    }
}

Real
EBArith::extrapFaceGradToOutflow(const FaceIndex&      a_bndryFace,
                                 const Side::LoHiSide& a_side,
                                 const int&            a_idir,
                                 const EBGraph&        a_ebGraph,
                                 const EBFaceFAB&      a_faceData,
                                 const int&            a_comp)
{
  Real extrapValue = -1.e99;
  Side::LoHiSide flipSide = flip(a_side);
  const VolIndex & closeVoF = a_bndryFace.getVoF(flipSide);
  Vector<FaceIndex> nearFaces = a_ebGraph.getFaces(closeVoF, a_idir, flipSide);
  bool hasNearFace = ((nearFaces.size() == 1)  && !nearFaces[0].isBoundary());
  bool hasFarFace  = false;
  bool hasFarFarFace  = false;
  Real nearVal = 0.0;
  Real farVal  = 0.0;
  Real farFarVal  = 0.0;
  if (hasNearFace)
    {
      nearVal = a_faceData(nearFaces[0], 0);
      const VolIndex & nextVoF = nearFaces[0].getVoF(flipSide);
      Vector<FaceIndex> farFaces = a_ebGraph.getFaces(nextVoF, a_idir, flipSide);
      hasFarFace = ((farFaces.size() == 1) && !farFaces[0].isBoundary());
      if (hasFarFace)
        {
          farVal = a_faceData(farFaces[0], 0);
          const VolIndex & nextNextVoF = farFaces[0].getVoF(flipSide);
          Vector<FaceIndex> farFarFaces = a_ebGraph.getFaces(nextNextVoF, a_idir, flipSide);
          hasFarFarFace = ((farFarFaces.size() == 1) && !farFarFaces[0].isBoundary());
          if (hasFarFarFace)
            {
              farFarVal = a_faceData(farFarFaces[0], 0);
            }
        }
    }

  // if (hasNearFace && hasFarFace && hasFarFarFace)
  //   {
  //     extrapValue = 3.*nearVal - 3.*farVal + farFarVal;
  //   }
  // else if (hasNearFace && hasFarFace)
  // // if (hasNearFace && hasFarFace)
  //   {
  //     extrapValue = 2.*nearVal - farVal;
  //   }
  // else if (hasNearFace)
  if (hasNearFace)
    {
      extrapValue = nearVal;//all that's needed to preserve constant gradient
    }
  else
    {
      extrapValue = 0.0; //for want of a better option.
    }
  return extrapValue;
}
Real
EBArith::extrapFaceValueToDomain(const FaceIndex&      a_bndryFace,
                                 const Side::LoHiSide& a_side,
                                 const int&            a_idir,
                                 const EBGraph&        a_ebGraph,
                                 const EBFaceFAB&      a_faceData,
                                 const int&            a_comp)
{
  Real extrapValue = -1.e99;
  Side::LoHiSide flipSide = flip(a_side);
  const VolIndex & closeVoF = a_bndryFace.getVoF(flipSide);
  Vector<FaceIndex> nearFaces = a_ebGraph.getFaces(closeVoF, a_idir, flipSide);
  bool hasNearFace = ((nearFaces.size() == 1)  && !nearFaces[0].isBoundary());
  bool hasFarFace  = false;
  bool hasFarFarFace  = false;
  Real nearVal = 0.0;
  Real farVal  = 0.0;
  Real farFarVal  = 0.0;
  if (hasNearFace)
    {
      nearVal = a_faceData(nearFaces[0], 0);
      const VolIndex & nextVoF = nearFaces[0].getVoF(flipSide);
      Vector<FaceIndex> farFaces = a_ebGraph.getFaces(nextVoF, a_idir, flipSide);
      hasFarFace = ((farFaces.size() == 1) && !farFaces[0].isBoundary());
      if (hasFarFace)
        {
          farVal = a_faceData(farFaces[0], 0);
          const VolIndex & nextNextVoF = farFaces[0].getVoF(flipSide);
          Vector<FaceIndex> farFarFaces = a_ebGraph.getFaces(nextNextVoF, a_idir, flipSide);
          hasFarFarFace = ((farFarFaces.size() == 1) && !farFarFaces[0].isBoundary());
          if (hasFarFarFace)
            {
              farFarVal = a_faceData(farFarFaces[0], 0);
            }
        }
    }

  if (hasNearFace && hasFarFace && hasFarFarFace)
    {
      extrapValue = 3.*nearVal - 3.*farVal + farFarVal;
    }
  else if (hasNearFace && hasFarFace)
    {
      extrapValue = 2.*nearVal - farVal;
    }
  else if (hasNearFace)
    {
      extrapValue = nearVal;//all that's needed to preserve constant gradient
    }
  else
    {
      extrapValue = 0.0; //for want of a better option.
    }
  return extrapValue;
}
Real
EBArith::extrapFaceVelToOutflow(const FaceIndex&      a_bndryFace,
                                const Side::LoHiSide& a_side,
                                const int&            a_idir,
                                const EBGraph&        a_ebGraph,
                                const EBFaceFAB&      a_faceData,
                                const int&            a_comp)
{
  Real extrapValue = -1.e99;
  Side::LoHiSide flipSide = flip(a_side);
  const VolIndex & closeVoF = a_bndryFace.getVoF(flipSide);
  Vector<FaceIndex> nearFaces = a_ebGraph.getFaces(closeVoF, a_idir, flipSide);
  bool hasNearFace = ((nearFaces.size() == 1)  && !nearFaces[0].isBoundary());
  bool hasFarFace  = false;
  bool hasFarFarFace  = false;
  Real nearVal = 0.0;
  Real farVal  = 0.0;
  Real farFarVal  = 0.0;
  if (hasNearFace)
    {
      nearVal = a_faceData(nearFaces[0], 0);
      const VolIndex & nextVoF = nearFaces[0].getVoF(flipSide);
      Vector<FaceIndex> farFaces = a_ebGraph.getFaces(nextVoF, a_idir, flipSide);
      hasFarFace = ((farFaces.size() == 1) && !farFaces[0].isBoundary());
      if (hasFarFace)
        {
          farVal = a_faceData(farFaces[0], 0);
          const VolIndex & nextNextVoF = farFaces[0].getVoF(flipSide);
          Vector<FaceIndex> farFarFaces = a_ebGraph.getFaces(nextNextVoF, a_idir, flipSide);
          hasFarFarFace = ((farFarFaces.size() == 1) && !farFarFaces[0].isBoundary());
          if (hasFarFarFace)
            {
              farFarVal = a_faceData(farFarFaces[0], 0);
            }
        }
    }

  // if (hasNearFace && hasFarFace && hasFarFarFace)
  //   {
  //     extrapValue = 3.*nearVal - 3.*farVal + farFarVal;
  //   }
  // else if (hasNearFace && hasFarFace)
  //   {
  //     extrapValue = 2.*nearVal - farVal;
  //   }
  if (hasNearFace && hasFarFace)
    {
      //quadratic extrapolation on v with the bc dv/dx=0 applied
      extrapValue = (4.*nearVal - farVal)/3.;
    }
  else if (hasNearFace)
    {
      extrapValue = nearVal;
    }
  else
    {
      extrapValue = 0.0; //for want of a better option.
    }
  return extrapValue;
}
Real
EBArith::interpolateVel(const EBFaceFAB& a_vel,
                        const EBISBox&   a_ebisBox,
                        const FaceIndex& a_face)
{

  //returns sum(weights*a_vels over face interpolation stencil)
  //to give velocity at face centroid
  const ProblemDomain& domain = a_ebisBox.getDomain();
  //like most things, this will break
  //if the coarse-fine interface intersects the
  //embedded boundary
  IntVectSet cfivs;
  FaceStencil sten = EBArith::getInterpStencil(a_face, cfivs, a_ebisBox, domain);
  Real velFace = 0.;
  for (int isten = 0; isten < sten.size(); isten++)
    {
      const Real&      weight = sten.weight(isten);
      const FaceIndex& stFace = sten.face(isten);
      Real velSten = a_vel(stFace, 0);
      velFace += weight*velSten;
    }
  return velFace;
}


Real
EBArith::deInterpolateVel(const Real&      a_centroidVel,
                          const EBFaceFAB& a_vel,
                          const EBISBox&   a_ebisBox,
                          const FaceIndex& a_face)
{

  const ProblemDomain& domain = a_ebisBox.getDomain();
  //like most things, this will break
  //if the coarse-fine interface intersects the
  //embedded boundary
  IntVectSet cfivs;
  FaceStencil sten = EBArith::getInterpStencil(a_face, cfivs, a_ebisBox, domain);
  Real velFace   = a_centroidVel;
  Real weightFace = -1;
  bool found = false;
  for (int isten = 0; isten < sten.size(); isten++)
    {
      const Real&      weight = sten.weight(isten);
      const FaceIndex& stFace = sten.face(isten);
      if (stFace == a_face)
        {
          found = true;
          weightFace = weight;
        }
      else
        {
          Real velSten = a_vel(stFace, 0);
          velFace -= weight*velSten;
        }
    }
  if (found && (weightFace != 0.0))
    {
      velFace /= weightFace;
    }
  else
    {
      //current face is not part of the stencil.  cannot think of a case where
      //this should happen.  return input vel as a default
      MayDay::Warning("wacky face interpolation stencil");
      velFace = a_centroidVel;
    }
  return velFace;
}

void
EBArith::
getLeastSquaresGradSten(VoFStencil&     a_stencil,
                        Real&           a_weight,
                        const VolIndex& a_vof,
                        const EBISBox&  a_ebisBox,
                        const RealVect& a_dx,
                        const ProblemDomain& a_domain,
                        int a_ivar)
{
  const RealVect normal   = a_ebisBox.normal(a_vof);
  const RealVect centroid = a_ebisBox.bndryCentroid(a_vof);
  getLeastSquaresGradSten(a_stencil, a_weight, normal, centroid,
                          a_vof, a_ebisBox, a_dx, a_domain, a_ivar);
  // getLeastSquaresGradStenAllVoFs(a_stencil, a_weight, normal, centroid,
  //                                a_vof, a_ebisBox,a_dx, a_domain, a_ivar);
  // getLeastSquaresGradStenAllQuad(a_stencil, a_weight, normal, centroid,
  //                                a_vof, a_ebisBox, a_dx, a_domain, a_ivar);
}
// This is the original way of doing the least squares stencil when eb normal points out of domain
// Only use stencil of all available vofs when symmetry is needed
void
EBArith::
getLeastSquaresGradSten(VoFStencil&     a_stencil,
                        Real&           a_weight,
                        const RealVect& a_normal  ,
                        const RealVect& a_centroid,
                        const VolIndex& a_vof,
                        const EBISBox&  a_ebisBox,
                        const RealVect& a_dx,
                        const ProblemDomain& a_domain,
                        int a_ivar)
{
  bool needSymStencil = true;
  for (int idir=0; idir<SpaceDim; idir++)
    {
      if (needSymStencil)
        {
          if (abs(a_normal[idir]) != 1. && a_normal[idir] != 0.)
            {
              needSymStencil = false;
            }
        }
    }

  if (needSymStencil)
    {
      getLeastSquaresGradStenAllQuad(a_stencil, a_weight, a_normal,
                                     a_centroid, a_vof, a_ebisBox,
                                     a_dx,a_domain, a_ivar, true);
      // getLeastSquaresGradStenAllVoFs(a_stencil, a_weight, a_normal,
      //                                a_centroid, a_vof, a_ebisBox,
      //                                a_dx,a_domain, a_ivar);
    }
  else
    {
      IntVect quadrant;
      const RealVect& normal = a_ebisBox.normal(a_vof);
      IntVect iv0 = a_vof.gridIndex();

      Box domain = a_domain.domainBox();
      for (int idir=0; idir<SpaceDim; idir++)
        {
          // if (iv0[idir]==domain.smallEnd(idir))
          //   {
          //     quadrant[idir]=1;
          //   }
          // else if (iv0[idir]==domain.bigEnd(idir))
          //   {
          //     quadrant[idir]=-1;
          //   }
          // else
          //   {
              if (normal[idir] < 0)
                {
                  quadrant[idir]=-1;
                }
              else
                {
                  quadrant[idir]=1;
                }
            // }
        }

      getLeastSquaresGradSten(a_stencil, a_weight, a_normal, a_centroid,
                              quadrant, a_vof, a_ebisBox, a_dx, a_domain,
                              a_ivar);

      if (a_stencil.size()==0)
        {
          getLeastSquaresGradStenAllQuad(a_stencil, a_weight, a_normal,
                                          a_centroid, a_vof, a_ebisBox,
                                          a_dx,a_domain, a_ivar);
           if (a_stencil.size()==0)
             {
               getLeastSquaresGradStenAllVoFs(a_stencil, a_weight, a_normal,
                                              a_centroid, a_vof, a_ebisBox,
                                              a_dx,a_domain, a_ivar);
             }
           else
             {
               pout()<<"leastSquares stencil dropping order "<<a_vof<<endl;
             }
        }
    }
}
// This is the actual solution to the least squares problem with the inversion of matrices and solution in flux form
void
EBArith::
getLeastSquaresGradSten(VoFStencil&     a_stencil,
                        Real&           a_weight,
                        const RealVect& a_normal,
                        const RealVect& a_centroid,
                        const IntVect&  a_quadrant,
                        const VolIndex& a_vof,
                        const EBISBox&  a_ebisBox,
                        const RealVect& a_dx,
                        const ProblemDomain& a_domain,
                        int a_ivar)
{
  IntVect iv0 = a_vof.gridIndex();

#if   CH_SPACEDIM == 2
  int stenSize = 3;
#elif CH_SPACEDIM == 3
  int stenSize = 7;
#else
  THIS_IS_AN_ERROR_MESSAGE__THIS_WILL_ONLY_COMPILE_WHEN_CH_SPACEDIM_IS_2_OR_3;
#endif
  Box domainBox = a_domain.domainBox();
  Vector<IntVect> ivSten(stenSize);

  ivSten[0] = iv0 + a_quadrant[0]*BASISV(0)                                                    ;
  ivSten[1] = iv0                           + a_quadrant[1]*BASISV(1)                          ;
  ivSten[2] = iv0 + a_quadrant[0]*BASISV(0) + a_quadrant[1]*BASISV(1)                          ;
#if CH_SPACEDIM == 3
  ivSten[3] = iv0                                                     + a_quadrant[2]*BASISV(2);
  ivSten[4] = iv0 + a_quadrant[0]*BASISV(0)                           + a_quadrant[2]*BASISV(2);
  ivSten[5] = iv0                           + a_quadrant[1]*BASISV(1) + a_quadrant[2]*BASISV(2);
  ivSten[6] = iv0 + a_quadrant[0]*BASISV(0) + a_quadrant[1]*BASISV(1) + a_quadrant[2]*BASISV(2);
#endif

  bool dropOrder = false;

  Vector<VolIndex> volSten(stenSize);
  for (int isten = 0; isten < stenSize; isten++)
    {
      //cp: it needs to be populated anyways
      if (a_ebisBox.getDomain().contains(ivSten[isten]))
        {
          if (a_ebisBox.getVoFs(ivSten[isten]).size() > 0)
            {
              volSten[isten] = a_ebisBox.getVoFs(ivSten[isten])[0];
            }
          else
            {
              dropOrder = true;
              volSten[isten] = VolIndex(IntVect(D_DECL(0,0,0)),0);
              // break;
            }
        }
      else
        {
          volSten[isten] = VolIndex(IntVect(D_DECL(0,0,0)),0);
          dropOrder = true;
          // break;
        }
    }

  std::vector<bool> monotonePath(stenSize);
  int nConn = stenSize;

  //restrictive because if one of the cells in the quadrant is not there then drop order
  for (int isten = 0; isten < stenSize; isten++)
    {
      monotonePath[isten] = EBArith::monotonePathVoFToCellVoF(volSten[isten],
                                                              a_vof,
                                                              ivSten[isten],
                                                              a_ebisBox);
      if (!monotonePath[isten])
        {
          dropOrder = true;
          nConn--;
        }
    }

  if (!dropOrder)
    {
      RealVect x0;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          x0[idir] = a_dx[idir] * (0.5 + a_centroid[idir] + iv0[idir]);
        }

      Vector<RealVect> xp(stenSize);
      for (int isten = 0; isten < stenSize; isten++)
        {
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              xp[isten][idir] = a_dx[idir] * (0.5 + ivSten[isten][idir]);
            }
        }

      Vector<RealVect> invATransAdeltaX(stenSize,RealVect::Zero);
      bool detZero = false;
      EBArith::calculateWeightingMatrix(x0, xp, invATransAdeltaX, detZero);

      a_stencil.clear();
      a_weight = 0.0;

      //if (!detZero)
        {
          for (int isten = 0; isten < stenSize; isten++)
            {
              Real dphidnWeight = 0.0;
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  dphidnWeight -= invATransAdeltaX[isten][idir] * a_normal[idir];
                }
              
              a_stencil.add(volSten[isten],dphidnWeight, a_ivar);
              a_weight -= dphidnWeight;
            }
        }
    }
  else
    {
      bool deadCell = true;
      if (nConn < 1)
        {
          deadCell = true;
        }
      else
        {
          //CP changed
          //collect all potentially usable points for the stencil
          RealVect x0;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              x0[idir] = a_dx[idir] * (0.5 + a_centroid[idir] + iv0[idir]);
            }

          Vector<RealVect> xp;
          Vector<int> volStenIdx;
          int ns = 0;
          for (int isten = 0; isten < stenSize; isten++)
            {
              if (monotonePath[isten])
                {
                  xp.push_back(RealVect::Zero);
                  ns++;
                  volStenIdx.push_back(isten);
                  for (int idir = 0; idir < SpaceDim; idir++)
                    {
                      xp[ns-1][idir] = a_dx[idir] * (0.5 + ivSten[isten][idir]);
                    }
                }
            }

          //determine which dimension to choose
          IntVect dimm = IntVect::Zero;
          IntVect mainDir = IntVect::Zero;
          //int diagDir[SpaceDim][2*SpaceDim - 1];
          int minDiag;
          int nDiagN = 2*(SpaceDim-1)-1; // the number of diagonal neighbors in one direction. Is this calculation correct? It is for 2D and 3D
          // mainDir and diagDir both contain neighbors' indices into monotonePath
          // main direction is its side neighbor, diagDir is the corner neighbors
          // how do I get diagDir? by looking at my 3D visualization box
          // by setting minDiag, we control, in the case of mainDir missing, how many corner neighbors must exist to make this dimension "valid". To always require the side neighbor, set this > 3 in 3D and > 1 in 2D

#if   CH_SPACEDIM == 2
          int diagDir[2][1] =
          {
            {
              2
            },
            {
              2
            }
          };
          mainDir   = IntVect(0,1);
          minDiag= 2; //minDiag = 2 this will always require the face neighbor
#elif CH_SPACEDIM == 3
          mainDir   = IntVect(0,1,3);
          int diagDir[3][3] =
          {
            {
              2,4,6
            },
            {
              2,5,6
            },
            {4,5,6
            }
          };
          minDiag= 2;
#else
          THIS_IS_AN_ERROR_MESSAGE__THIS_WILL_ONLY_COMPILE_WHEN_CH_SPACEDIM_IS_2_OR_3;
#endif
          for (int iDir = 0; iDir < SpaceDim; iDir++)
            {
              if (monotonePath[mainDir[iDir]])
                {
                  int diagCount = 0;
                  for (int jDiag = 0; jDiag< nDiagN; jDiag++)
                    {
                      if (monotonePath[diagDir[iDir][jDiag] ]) diagCount++;
                    }
                  if (diagCount >= 1) dimm[iDir]=1;
                }
              else
                // face neighbor covered, counting number of diagonal neighbors
                {
                  int diagCount = 0;
                  for (int jDiag = 0; jDiag< nDiagN; jDiag++)
                    {
                      if (monotonePath[diagDir[iDir][jDiag] ]) diagCount++;
                    }
                  if (diagCount >= minDiag) dimm[iDir]=1;
                }
            }
          if (dimm.sum()<1)
            {
              deadCell = true;
            }
          else
            {
              //calculate weights, corresponding to all stencil points
              Vector<RealVect> invATransAdeltaX(ns,RealVect::Zero);
              EBArith::calculateWeightingMatrixRed(x0,xp,dimm,invATransAdeltaX, deadCell);

              a_stencil.clear();
              a_weight = 0.0;

              for (int isten = 0; isten < xp.size(); isten++)
                {
                  Real dphidnWeight = 0.0;
                  for (int idir = 0; idir < SpaceDim; idir++)
                    {
                      dphidnWeight -= invATransAdeltaX[isten][idir] * a_normal[idir];
                    }

                  a_stencil.add(volSten[volStenIdx[isten]],dphidnWeight, a_ivar);
                  a_weight -= dphidnWeight;
                }
            }
        }
      if (deadCell)
        {
          a_stencil.clear();
          a_weight = 0.0;
        }
    }
}
void
EBArith::
getLeastSquaresGradStenAllVoFs(VoFStencil&          a_stencil,
                               Real&                a_weight,
                               const RealVect&      a_normal,
                               const RealVect&      a_centroid,
                               const VolIndex&      a_vof,
                               const EBISBox&       a_ebisBox,
                               const RealVect&      a_dx,
                               const ProblemDomain& a_domain,
                               int                  a_ivar)
{
  IntVect iv0 = a_vof.gridIndex();

#if   CH_SPACEDIM == 2
  int numNeighbors = 8;
  int minStenSize  = 3;
#elif CH_SPACEDIM == 3
  int numNeighbors = 26;
  int minStenSize  = 7;
#else
  THIS_IS_AN_ERROR_MESSAGE__THIS_WILL_ONLY_COMPILE_WHEN_CH_SPACEDIM_IS_2_OR_3;
#endif
  Box domainBox = a_domain.domainBox();
  Vector<IntVect> ivSten(numNeighbors);

  ivSten[0] = iv0 + BASISV(0)            ;
  ivSten[1] = iv0             + BASISV(1);
  ivSten[2] = iv0 + BASISV(0) + BASISV(1);
  ivSten[3] = iv0 - BASISV(0)            ;
  ivSten[4] = iv0 - BASISV(0) + BASISV(1);
  ivSten[5] = iv0             - BASISV(1);
  ivSten[6] = iv0 + BASISV(0) - BASISV(1);
  ivSten[7] = iv0 - BASISV(0) - BASISV(1);
#if CH_SPACEDIM == 3
  ivSten[ 8] = iv0                         + BASISV(2);
  ivSten[ 9] = iv0                         - BASISV(2);
  ivSten[10] = iv0 + BASISV(0)             + BASISV(2);
  ivSten[11] = iv0 + BASISV(0)             - BASISV(2);
  ivSten[12] = iv0 - BASISV(0)             + BASISV(2);
  ivSten[13] = iv0 - BASISV(0)             - BASISV(2);
  ivSten[14] = iv0             + BASISV(1) + BASISV(2);
  ivSten[15] = iv0             + BASISV(1) - BASISV(2);
  ivSten[16] = iv0             - BASISV(1) + BASISV(2);
  ivSten[17] = iv0             - BASISV(1) - BASISV(2);
  ivSten[18] = iv0 + BASISV(0) + BASISV(1) + BASISV(2);
  ivSten[19] = iv0 - BASISV(0) + BASISV(1) + BASISV(2);
  ivSten[20] = iv0 + BASISV(0) - BASISV(1) + BASISV(2);
  ivSten[21] = iv0 + BASISV(0) + BASISV(1) - BASISV(2);
  ivSten[22] = iv0 - BASISV(0) - BASISV(1) + BASISV(2);
  ivSten[23] = iv0 - BASISV(0) + BASISV(1) - BASISV(2);
  ivSten[24] = iv0 + BASISV(0) - BASISV(1) - BASISV(2);
  ivSten[25] = iv0 - BASISV(0) - BASISV(1) - BASISV(2);
#endif

  bool dropOrder = false;

  Vector<VolIndex> volSten;
  for (int inbr=0; inbr<numNeighbors; inbr++)
    {
      if (a_ebisBox.getDomain().contains(ivSten[inbr]) && a_ebisBox.getVoFs(ivSten[inbr]).size() > 0)
        {
          VolIndex vof;
          bool vofIsThere = EBArith::monotonePathVoFToCellVoF(vof,
                                                              a_vof,
                                                              ivSten[inbr],
                                                              a_ebisBox);
          if (vofIsThere) volSten.push_back(vof);
        }
    }

  int stenSize = volSten.size();
  if (stenSize < minStenSize)
    {
      dropOrder = true;
    }

  if (!dropOrder)
    {
      RealVect x0;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          x0[idir] = a_dx[idir] * (0.5 + a_centroid[idir] + iv0[idir]);
        }

      Vector<RealVect> xp(stenSize);
      for (int isten = 0; isten < stenSize; isten++)
        {
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              xp[isten][idir] = a_dx[idir] * (0.5 + volSten[isten].gridIndex()[idir]);
            }
        }

      Vector<RealVect> invATransAdeltaX(stenSize,RealVect::Zero);
      bool detZero = false;
      EBArith::calculateWeightingMatrix(x0, xp, invATransAdeltaX, detZero);

      a_stencil.clear();
      a_weight = 0.0;

      //if (!detZero)
      // {
          for (int isten = 0; isten < stenSize; isten++)
            {
              Real dphidnWeight = 0.0;
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  dphidnWeight -= invATransAdeltaX[isten][idir] * a_normal[idir];
                }
              
              a_stencil.add(volSten[isten],dphidnWeight, a_ivar);
              a_weight -= dphidnWeight;
            }
          // }
    }
  else
    {
      a_stencil.clear();
      a_weight = 0.0;
    }
}
void
EBArith::
getLeastSquaresGradStenAllQuad(VoFStencil&          a_stencil,
                               Real&                a_weight,
                               const RealVect&      a_normal,
                               const RealVect&      a_centroid,
                               const VolIndex&      a_vof,
                               const EBISBox&       a_ebisBox,
                               const RealVect&      a_dx,
                               const ProblemDomain& a_domain,
                               int                  a_ivar,
                               bool                 a_doSymmetric)
{
  // if we've gotten here, then we shouldn't have a stencil
  a_stencil.clear();
  a_weight = 0.;
  // number of valid stencils so far
  int numValidStencils = 0;

  // calculate the quadrant
  IntVect quadrant;
  IntVect iv0 = a_vof.gridIndex();

  for (int idir=0; idir<SpaceDim; idir++)
    {
      if (a_normal[idir] < 0)
        {
          quadrant[idir]=-1;
        }
      else
        {
          quadrant[idir]=1;
        }
    }

  // next stencil candidate
  VoFStencil curStencil;
  IntVect curQuadrant;
  Real curWeight = 0.;
  bool validCurQuadrant;

  // in symmetric case and near domain boundaries this quadrant will not
  //  have been tried, so we try it here
  getLeastSquaresGradSten(curStencil, curWeight, a_normal, a_centroid,
                          quadrant, a_vof, a_ebisBox, a_dx, a_domain,
                          a_ivar);

  if (curStencil.size() != 0)
    {
      numValidStencils += 1;
      a_stencil += curStencil;
      a_weight += curWeight;
    }

  // cycle through directions, switching the sign of quadrant's
  //  components looking for a better one
  for (int idir=0; idir<SpaceDim; idir++)
    {
      curQuadrant = IntVect::Unit;
      // first try flipping this direction
      curQuadrant[idir] *= -1;
      validCurQuadrant = curQuadrant != quadrant;
      if (!a_doSymmetric)
        {
          validCurQuadrant = validCurQuadrant && (curQuadrant != -quadrant);
        }
      if (validCurQuadrant)
        {
          curStencil.clear();
          curWeight = 0.;
          getLeastSquaresGradSten(curStencil, curWeight, a_normal, a_centroid,
                                  curQuadrant, a_vof, a_ebisBox, a_dx, a_domain,
                                  a_ivar);

          if (curStencil.size() != 0)
            {
              numValidStencils += 1;
              a_stencil += curStencil;
              a_weight += curWeight;
            }
        }
#if CH_SPACEDIM == 3
      // now try flipping a second direction
      int jdir = (idir+1)%SpaceDim;
      curQuadrant[jdir] *= -1;
      validCurQuadrant = curQuadrant != quadrant;
      if (!a_doSymmetric)
        {
          validCurQuadrant = validCurQuadrant && (curQuadrant != -quadrant);
        }
      if (validCurQuadrant)
        {
          curStencil.clear();
          curWeight = 0.;
          getLeastSquaresGradSten(curStencil, curWeight, a_normal, a_centroid,
                                  curQuadrant, a_vof, a_ebisBox, a_dx, a_domain,
                                  a_ivar);

          if (curStencil.size() != 0)
            {
              numValidStencils += 1;
              a_stencil += curStencil;
              a_weight += curWeight;
            }
        }
#endif
    }
  // lastly, try flipping all directions
  curQuadrant = IntVect::Unit;
  curQuadrant *= -1;
  validCurQuadrant = curQuadrant != quadrant;
  if (!a_doSymmetric)
    {
      validCurQuadrant = validCurQuadrant && (curQuadrant != -quadrant);
    }
  if (validCurQuadrant)
    {
      curStencil.clear();
      curWeight = 0.;
      getLeastSquaresGradSten(curStencil, curWeight, a_normal, a_centroid,
                              curQuadrant, a_vof, a_ebisBox, a_dx, a_domain,
                              a_ivar);

      if (curStencil.size() != 0)
        {
          numValidStencils += 1;
          a_stencil += curStencil;
          a_weight += curWeight;
        }
    }

  if (numValidStencils > 1)
    {
      a_stencil *= 1./numValidStencils;
      a_weight *= 1./numValidStencils;
    }
}

void
EBArith::calculateWeightingMatrix(RealVect           x0,
                                  Vector<RealVect>&  xp,
                                  Vector<RealVect>&  weightMatrix,
                                  bool&              detZero)
{
  int stenSize = xp.size();

  Vector<RealVect> deltaX = xp;
  for (int isten = 0; isten < stenSize; isten++)
    {
      deltaX[isten] -= x0;
    }

  Vector<RealVect> aTransA(SpaceDim,RealVect::Zero), invATransA(SpaceDim,RealVect::Zero);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (int jdir = 0; jdir < SpaceDim; jdir++)
        {
          for (int isten = 0; isten < stenSize; isten++)
            {
              aTransA[idir][jdir] = aTransA[idir][jdir]
                                  + deltaX[isten][idir]*deltaX[isten][jdir];
            }
        }
    }

  Real det;
#if   CH_SPACEDIM == 2
  det = aTransA[0][0] * aTransA[1][1] - aTransA[0][1] * aTransA[1][0];
  if (det < 1.e-15 && det > -1.e-15)
    {
      detZero = true;
    }
  else
    {
      invATransA[0][0] =  aTransA[1][1] / det;
      invATransA[0][1] = -aTransA[0][1] / det;
      invATransA[1][0] = -aTransA[1][0] / det;
      invATransA[1][1] =  aTransA[0][0] / det;
    }
#elif CH_SPACEDIM == 3
  det = aTransA[0][0] * ( aTransA[1][1] * aTransA[2][2]
                        - aTransA[1][2] * aTransA[2][1])
      + aTransA[0][1] * ( aTransA[1][2] * aTransA[2][0]
                        - aTransA[1][0] * aTransA[2][2])
      + aTransA[0][2] * ( aTransA[1][0] * aTransA[2][1]
                        - aTransA[1][1] * aTransA[2][0]);

  if (det < 1.e-15 && det > -1.e-15)
    {
      detZero = true;
    }
  else
    {
      invATransA[0][0] = ( aTransA[1][1] * aTransA[2][2]
                           - aTransA[1][2] * aTransA[2][1]) / det;
      invATransA[0][1] = ( aTransA[1][2] * aTransA[2][0]
                           - aTransA[1][0] * aTransA[2][2]) / det;
      invATransA[0][2] = ( aTransA[1][0] * aTransA[2][1]
                           - aTransA[1][1] * aTransA[2][0]) / det;
      invATransA[1][0] = ( aTransA[2][1] * aTransA[0][2]
                           - aTransA[2][2] * aTransA[0][1]) / det;
      invATransA[1][1] = ( aTransA[2][2] * aTransA[0][0]
                           - aTransA[2][0] * aTransA[0][2]) / det;
      invATransA[1][2] = ( aTransA[2][0] * aTransA[0][1]
                           - aTransA[2][1] * aTransA[0][0]) / det;
      invATransA[2][0] = ( aTransA[0][1] * aTransA[1][2]
                           - aTransA[0][2] * aTransA[1][1]) / det;
      invATransA[2][1] = ( aTransA[0][2] * aTransA[1][0]
                           - aTransA[0][0] * aTransA[1][2]) / det;
      invATransA[2][2] = ( aTransA[0][0] * aTransA[1][1]
                           - aTransA[0][1] * aTransA[1][0]) / det;
    }
#else
  THIS_IS_AN_ERROR_MESSAGE__THIS_WILL_ONLY_COMPILE_WHEN_CH_SPACEDIM_IS_2_OR_3;
#endif

  //if (!detZero)
    {
      weightMatrix.resize(stenSize,RealVect::Zero);
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          for (int isten = 0; isten < stenSize; isten++)
            {
              for (int jdir = 0; jdir < SpaceDim; jdir++)
                {
                  weightMatrix[isten][idir] += invATransA[idir][jdir] * deltaX[isten][jdir];
                }
            }
        }
    }
}


void
EBArith::calculateWeightingMatrixRed(RealVect           x00,
                                     Vector<RealVect>&  xpp,
                                     IntVect            dimm,
                                     Vector<RealVect>&  weightMatrix,
                                     bool&              deadRed)
//CP: do the same thing for a reduced system, where some neighbors in the normal leastSquare stencil are covered
//some dimensions might also have vanished. these need to be recorded outside
//dimm[idir]==1, idir is valid, dimm[idir]==0, idir is missing

{

  int stenSize = xpp.size();

  //now cast the problem to reduced dimension: make x0, xp
  int nr = 0;
  RealVect x0=RealVect::Zero;
  Vector<RealVect> xp(stenSize,RealVect::Zero);
  IntVect dirN;
  for (int idir = 0; idir< SpaceDim; idir++)
    {
      if (dimm[idir]>0)
        {
          nr++;
          x0[nr-1] = x00[idir];
          dirN[nr-1] = idir;
          for (int isten = 0; isten < stenSize; isten++)
            {
              xp[isten][nr-1]=xpp[isten][idir];
            }
        }
    }

  Vector<RealVect> deltaX = xp;
  for (int isten = 0; isten < stenSize; isten++)
    {
      deltaX[isten] -= x0;
    }

  Vector<RealVect> aTransA(SpaceDim,RealVect::Zero), invATransA(SpaceDim,RealVect::Zero);
  //CP: using the fact that nr <= SpaceDim

  for (int idir = 0; idir < nr; idir++)
    {
      for (int jdir = 0; jdir < nr; jdir++)
        {
          for (int isten = 0; isten < stenSize; isten++)
            {
              aTransA[idir][jdir] = aTransA[idir][jdir]
                                  + deltaX[isten][idir]*deltaX[isten][jdir];
            }
        }
    }

  Real det;
  if (nr == 1)
    {
      for (int isten = 0; isten< stenSize; isten++)
        {
          // this is worth more consideration when there is only one dimension
          // should all cells be included?
          invATransA[0][0] =  invATransA[0][0] + deltaX[isten][0]*deltaX[isten][0];
        }
      invATransA[0][0]=1/ invATransA[0][0];
      if ((invATransA[0][0]) == 0.0) deadRed = true;
    }
  else if (nr == 2)
    {
      det = aTransA[0][0] * aTransA[1][1] - aTransA[0][1] * aTransA[1][0];

      invATransA[0][0] =  aTransA[1][1] / det;
      invATransA[0][1] = -aTransA[0][1] / det;
      invATransA[1][0] = -aTransA[1][0] / det;
      invATransA[1][1] =  aTransA[0][0] / det;
      if ((det) == 0.0) deadRed = true;
    }
  else if (nr == 3)
    {
      det = aTransA[0][0] * ( aTransA[1][1] * aTransA[2][2]
                              - aTransA[1][2] * aTransA[2][1])
        + aTransA[0][1] * ( aTransA[1][2] * aTransA[2][0]
                            - aTransA[1][0] * aTransA[2][2])
        + aTransA[0][2] * ( aTransA[1][0] * aTransA[2][1]
                            - aTransA[1][1] * aTransA[2][0]);
      if (det< 0)
        {
          // this matrix has no inverse. one dimension has no solution
          // what to do?
          CH_assert(det >= 0);
        }
      else if (det == 0)
        {
          deadRed = true;
        }
      else
        {
          invATransA[0][0] = ( aTransA[1][1] * aTransA[2][2]
                               - aTransA[1][2] * aTransA[2][1]) / det;
          invATransA[0][1] = ( aTransA[1][2] * aTransA[2][0]
                               - aTransA[1][0] * aTransA[2][2]) / det;
          invATransA[0][2] = ( aTransA[1][0] * aTransA[2][1]
                               - aTransA[1][1] * aTransA[2][0]) / det;
          invATransA[1][0] = ( aTransA[2][1] * aTransA[0][2]
                               - aTransA[2][2] * aTransA[0][1]) / det;
          invATransA[1][1] = ( aTransA[2][2] * aTransA[0][0]
                               - aTransA[2][0] * aTransA[0][2]) / det;
          invATransA[1][2] = ( aTransA[2][0] * aTransA[0][1]
                               - aTransA[2][1] * aTransA[0][0]) / det;
          invATransA[2][0] = ( aTransA[0][1] * aTransA[1][2]
                               - aTransA[0][2] * aTransA[1][1]) / det;
          invATransA[2][1] = ( aTransA[0][2] * aTransA[1][0]
                               - aTransA[0][0] * aTransA[1][2]) / det;
          invATransA[2][2] = ( aTransA[0][0] * aTransA[1][1]
                               - aTransA[0][1] * aTransA[1][0]) / det;
        }
    }
  else
    {
      CH_assert(nr<=3 && nr>0);
    }

  weightMatrix.resize(stenSize,RealVect::Zero);
  for (int idir = 0; idir < nr; idir++)
    {
      for (int isten = 0; isten < stenSize; isten++)
        {
          for (int jdir = 0; jdir < nr; jdir++)
            {
              weightMatrix[isten][dirN[idir]] += invATransA[idir][jdir] * deltaX[isten][jdir];
              // direction offset: set to original dimension
            }
        }
    }
}


void
EBArith::dataRayCast(bool&               a_dropOrder,
                     Vector<VoFStencil>& a_pointStencils,
                     Vector<Real>&       a_distanceAlongLine,
                     const RealVect&     a_normal,
                     const RealVect&     a_bndryCentroid,
                     const VolIndex&     a_vof,
                     const EBISBox&      a_ebisBox,
                     const RealVect&     a_dx,
                     const IntVectSet&   a_cfivs,
                     int a_ivar,
                     int a_numPoints)
{
  a_dropOrder = false;

  CH_TIME("EBArith::johanStencil");
  IntVect iv = a_vof.gridIndex();

  //line starts at xb = location in physical space
  //of boundary centroid
  RealVect xb;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      xb[idir] = a_dx[idir]*(Real(iv[idir]) + 0.5 + a_bndryCentroid[idir]);
    }

  // Find InterpDirection and hiLo
  Real nMax = 0.0;
  int nMaxDir = 0;

  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      if (Abs(a_normal[idir]) > nMax)
        {
          nMax = Abs(a_normal[idir]);
          nMaxDir = idir;
        }
    }
  //sometimes normals can be zero
  if (Abs(nMax) < 1.0e-15)
    {
      a_dropOrder = true;
      return;
    }
  int hiLo = 1;
  if (a_normal[nMaxDir] < 0.0)
    {
      hiLo = -1;
    }

  // Find Intersection of planes and ray
  //and the cells in which they live
  Vector<RealVect> intersectLoc(a_numPoints);
  Vector<IntVect>   intersectIV(a_numPoints);
  //equation of line
  //y = y_b + (ny/nx)*(x-xb)
  //we know x because that is the location of the
  //next cell center over in the nMaxDir direction
  //hiLo in the stencil depends on where the intersection point is in
  //relation to the cell center of the intersecting cell
  a_distanceAlongLine.resize(a_numPoints);
  a_pointStencils.resize(a_numPoints);
  const Box region = a_ebisBox.getRegion();

  for (int iinter = 0; iinter < a_numPoints; iinter++)
    {
      intersectIV[iinter] = iv + (iinter+1)*hiLo*BASISV(nMaxDir);
      // check whether intersectIV occurs outside of iv[idir!=nMaxDir]
      for (int idir=0; idir<SpaceDim; idir++)
        {
          if (idir != nMaxDir)
            {
              // what direction are we looking in?
              int isign = (a_normal[idir]<0.0) ? -1 : 1;
              // how far do we go in idir as a result of moving
              //  (iinter+1) cells in nMaxDir?
              Real xDist = Abs((Real(iinter+1)-hiLo*a_bndryCentroid[nMaxDir])
                               *(a_normal[idir]/a_normal[nMaxDir]));
              // how far from the boundary centroid to the idir cell-edge
              Real xEdgeDist = Abs(Real(isign)*0.5 - a_bndryCentroid[idir]);
              if (xDist > xEdgeDist)
                { // we're outside of iv[idir]. calculate by how many cells
                  //  and adjust intersectIV accordingly
                  intersectIV[iinter][idir] += isign*int(1+floor(xDist-xEdgeDist));
                }
            }
        }

      if (!region.contains(intersectIV[iinter]))
        {
          a_dropOrder = true;
          return;
        }
      if (a_ebisBox.numVoFs(intersectIV[iinter]) != 1)
        {
          a_dropOrder = true;
          return;
        }
      VolIndex centerVoF(intersectIV[iinter], 0);

      //small cells as centers can be bad
      if (a_ebisBox.volFrac(centerVoF) < s_minVolFrac)
        {
          a_dropOrder = true;
          return;
        }

      Real xMaxDir  = a_dx[nMaxDir]*(Real(intersectIV[iinter][nMaxDir]) + 0.5);
      intersectLoc[iinter][nMaxDir] = xMaxDir;
      a_distanceAlongLine[iinter] = (xMaxDir-xb[nMaxDir])*(xMaxDir-xb[nMaxDir]);
      RealVect extrapDist = RealVect::Zero; //the initialization is important
      //get location of intersection and the distance along the line
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if (idir != nMaxDir)
            {
              Real normalRat = a_normal[idir]/a_normal[nMaxDir];
              Real distDir =  normalRat*(xMaxDir-xb[nMaxDir]);
              Real spaceLoc  =  xb[idir] + distDir;
              intersectLoc[iinter][idir] = spaceLoc;
              a_distanceAlongLine[iinter] += distDir*distDir;

              //the nMaxDir value is set with RealVect::Zero initialization
              //no extrapolation in the nmax dir
              Real ccLoc = a_dx[idir]*(Real(intersectIV[iinter][idir]) + 0.5);
              extrapDist[idir] = spaceLoc - ccLoc;
            }
        }

      a_distanceAlongLine[iinter] = sqrt(a_distanceAlongLine[iinter]);

      int order = getExtrapolationStencil(a_pointStencils[iinter], extrapDist,
                                          a_dx, centerVoF, a_ebisBox, nMaxDir,
                                          NULL, a_ivar);

      //the returned order is the order of the derivs taken.  can tolerate 1 or 2
      if (order == 0)
        {
          a_dropOrder = true;
          return;
        }
    }//end loop over intersection points

  return;
}
void
EBArith::johanStencil(bool&               a_dropOrder,
                      Vector<VoFStencil>& a_pointStencils,
                      Vector<Real>&       a_distanceAlongLine,
                      const VolIndex&     a_vof,
                      const EBISBox&      a_ebisBox,
                      const RealVect&     a_dx,
                      const IntVectSet&   a_cfivs,
                      int a_ivar)
{
  //CH_TIME("EBArith::johanStencil");
  IntVect iv = a_vof.gridIndex();

  // Do Johansen-Colella cast-a-ray algorithm
  RealVect bndryCentroid = a_ebisBox.bndryCentroid(a_vof);
  RealVect normal     = a_ebisBox.normal(a_vof);
  //splitting up stuff this way to facillitate multifluid
  //which can have multiple normals and boundary centroids per cell.
  johanStencil(a_dropOrder, a_pointStencils, a_distanceAlongLine,
               normal, bndryCentroid,
               a_vof, a_ebisBox, a_dx, a_cfivs, a_ivar);
}
void
EBArith::johanStencil(bool&               a_dropOrder,
                      Vector<VoFStencil>& a_pointStencils,
                      Vector<Real>&       a_distanceAlongLine,
                      const RealVect&     a_normal,
                      const RealVect&     a_bndryCentroid,
                      const VolIndex&     a_vof,
                      const EBISBox&      a_ebisBox,
                      const RealVect&     a_dx,
                      const IntVectSet&   a_cfivs,
                      int a_ivar)
{
  int numExtrapPoints = 2;
  EBArith::dataRayCast(a_dropOrder, a_pointStencils, a_distanceAlongLine, a_normal, a_bndryCentroid,
              a_vof, a_ebisBox, a_dx, a_cfivs, a_ivar, numExtrapPoints);
}
/******/
RealVect
EBArith::
getFaceLocation(const FaceIndex& a_face,
                const RealVect&  a_dx,
                const RealVect&  a_probLo)
{
  const IntVect& iv = a_face.gridIndex(Side::Hi);
  RealVect loc = a_probLo;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (idir == a_face.direction())
        {
          loc[idir] += a_dx[idir]*Real(iv[idir]);
        }
      else
        {
          loc[idir] += a_dx[idir]*(Real(iv[idir]) + 0.5);
        }
    }
  return loc;
}
RealVect
EBArith::
getVofLocation(const VolIndex& a_vof,
               const RealVect& a_dx,
               const RealVect& a_probLo)
{
  const IntVect& iv = a_vof.gridIndex();
  RealVect loc = getIVLocation(iv,a_dx,a_probLo);
  return loc;
}
RealVect
EBArith::
getIVLocation(const IntVect&  a_iv,
              const RealVect& a_dx,
              const RealVect& a_probLo)
{
  RealVect loc = a_probLo;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      loc[idir] += a_dx[idir]*(Real(a_iv[idir]) + 0.5);
    }
  return loc;
}
RealVect
EBArith::
getDomainNormal(int a_idir, Side::LoHiSide a_side)
{
  RealVect normal = BASISREALV(a_idir);
  if (a_side == Side::Hi)
    {
      normal *= -1.0;
    }
  return normal;
}
/***/
void
EBArith::
defineFluxInterpolant(LevelData<BaseIFFAB<Real> >& a_fluxInterpolant,
                     LayoutData<IntVectSet>     & a_irregSetsGrown,
                      const DisjointBoxLayout    & a_dbl,
                      const EBISLayout           & a_ebisl,
                      const ProblemDomain        & a_domain,
                      const int                  & a_ncomp,
                      const int                  & a_faceDir)
{
  a_irregSetsGrown.define(a_dbl);

  for (DataIterator dit = a_dbl.dataIterator(); dit.ok(); ++dit)
    {
      const Box&     grid = a_dbl.get(dit());
      Box grownBox2 = grow(grid, 2);
      grownBox2 &= a_domain;
      const EBISBox& ebisBox = a_ebisl[dit()];
      if (!ebisBox.isAllCovered())
        {
          IntVectSet& grownIVS = a_irregSetsGrown[dit()];

          //need to do this in case there is an irregular cells
          //just outside the box.
          IntVectSet irregivs = ebisBox.getIrregIVS(grownBox2);
          grownIVS = irregivs;
          Box grownBox1 = grid;
          //grow in all directions != faceDir
          for (int jdir = 0; jdir < SpaceDim; jdir++)
            {
              if (jdir != a_faceDir)
                {
                  grownIVS.grow(jdir, 1);
                  grownBox1.grow(jdir, 1);
                }
            }
          //restrict to domain
          grownBox1 &= a_domain;
          grownIVS &= grownBox1;
        }
    }

  //define flux interpolant stuff
  BaseIFFactory<Real> faceFactory(a_ebisl, a_irregSetsGrown, a_faceDir);
  IntVect ivghost = IntVect::Zero;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (a_faceDir != idir)
        {
          ivghost += BASISV(idir);
        }
    }
  a_fluxInterpolant.define(a_dbl, a_ncomp, ivghost, faceFactory);
}

/****/
void
EBArith::
interpolateFluxToCentroids(LevelData<EBFluxFAB>&       a_centroidFlux,
                           const LevelData<EBFluxFAB>& a_faceCentFlux,
                           const DisjointBoxLayout&    a_grids,
                           const EBISLayout&           a_ebisl,
                           const ProblemDomain&        a_domain)
{
  CH_TIME("EBArith::interpolateFluxToCentroids(level)");
  int ibox = 0;
  int idummy = 0;
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      EBFluxFAB&       centroidFluxFAB = a_centroidFlux[dit()];
      const EBFluxFAB& faceCentFluxFAB = a_faceCentFlux[dit()];
      const EBISBox&   ebisBox = a_ebisl[dit()];
      const Box& box = a_grids.get(dit());
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          EBFaceFAB&       centroidFaceFAB = centroidFluxFAB[idir];
          const EBFaceFAB& faceCentFaceFAB = faceCentFluxFAB[idir];
          interpolateFluxToCentroids(centroidFaceFAB, faceCentFaceFAB,  box, ebisBox, a_domain, idir);
          idummy++;
        }
      ibox++;
    }
}
/****/
void
EBArith::
interpolateFluxToCentroids(EBFaceFAB&                  a_centroidFlux,
                           const EBFaceFAB&            a_faceCentFlux,
                           const Box&                  a_box,
                           const EBISBox&              a_ebisBox,
                           const ProblemDomain&        a_domain,
                           const int&                  a_idir)
{
  CH_TIME("EBArith::interpolateFluxToCentroids");
  FaceStop::WhichFaces stopCrit = FaceStop::SurroundingWithBoundary;
  IntVectSet cfivs;

  //first copy all the fluxes over then interpolate on irregular cells
  a_centroidFlux.copy(a_faceCentFlux);
  IntVectSet irregIVS = a_ebisBox.getIrregIVS(a_box);

  for (FaceIterator faceit(irregIVS, a_ebisBox.getEBGraph(), a_idir, stopCrit);
      faceit.ok(); ++faceit)
    {
      const FaceIndex& face = faceit();
      FaceStencil sten = getInterpStencil(face, cfivs, a_ebisBox, a_domain);
      for (int ivar = 0; ivar < a_centroidFlux.nComp(); ivar++)
        {
          Real newFlux = 0.0;
          for (int isten = 0; isten < sten.size(); isten++)
            {
              const FaceIndex& stenFace = sten.face(isten);
              Real weight     = sten.weight(isten);
              Real interpFlux = a_faceCentFlux(stenFace, ivar);
              newFlux += weight*interpFlux;
            }
          a_centroidFlux(face, ivar) = newFlux;
        }
    }
}
/****/
void EBArith::
irregNorm(Real& a_ebIrregNorm,
          const BaseIVFAB<Real>& a_ebiError,
          const IntVectSet& a_ivsIrreg,
          const EBISBox& a_ebisBox,
          const int& a_comp,
          const int& a_normtype)
{
  CH_assert(a_normtype >= 0);
  CH_assert(a_normtype <= 2);
  a_ebIrregNorm = 0.0;
  if (a_normtype == 0)
    {
      for (VoFIterator vofit(a_ivsIrreg, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real bdarea = a_ebisBox.bndryArea(vof);
          if (bdarea > 0)
            {
              Real valVoF = a_ebiError(vof, a_comp);
              a_ebIrregNorm = Max(Abs(valVoF), a_ebIrregNorm);
            }
        }
    }
  else
    {
      //integral norm
      Real areaTot = 0.0;
      Real normTot = 0.0;
      for (VoFIterator vofit(a_ivsIrreg, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real valVoF = a_ebiError(vof, a_comp);
          Real bdarea = a_ebisBox.bndryArea(vof);
          areaTot += bdarea;
          if (a_normtype == 1)
            normTot += Abs(valVoF)*bdarea;
          else
            normTot += valVoF*valVoF*bdarea;
        }
      if (areaTot > 1.0e-8)
        normTot /= areaTot;
      if (a_normtype == 2)
        normTot = sqrt(normTot);
      a_ebIrregNorm = normTot;
    }

}
IntVect ebcoarsen (const IntVect& b,
                   int        refinement_ratio)
{
  return coarsen(b, refinement_ratio);
}

Box ebrefine (const Box& b,
              int        refinement_ratio)
{
  return refine(b, refinement_ratio);
}
ProblemDomain ebrefine (const ProblemDomain& b,
                        int        refinement_ratio)
{
  return refine(b, refinement_ratio);
}

Box ebcoarsen (const Box& b,
               int        refinement_ratio)
{
  return coarsen(b, refinement_ratio);
}

ProblemDomain ebcoarsen (const ProblemDomain& b,
                         int        refinement_ratio)
{
  return coarsen(b, refinement_ratio);
}

void ebrefine(DisjointBoxLayout& output,
              const DisjointBoxLayout& input,
              int refinement)
{
  return refine(output, input, refinement);
}

void ebcoarsen(DisjointBoxLayout& output,
               const DisjointBoxLayout& input,
               int refinement)
{
  return coarsen(output, input, refinement);
}
/*****************************/
void
EBArith::
getKVolRZ(Real&           a_kvol,
          Real&           a_cellVol,
          const EBISBox&  a_ebisBox,
          const Real&     a_dx,
          const VolIndex& a_vof)
{
  int rdir = 0;
  int rindex   =  a_vof.gridIndex()[rdir];
  Real cellRad = a_dx*(Real(rindex) + 0.5);
  a_cellVol =cellRad*a_dx*a_dx;

  //kvol = (1/2*cellvol)*int(r^2 nr dl)
  Real openContrib = 0.;
  for (SideIterator sit; sit.ok(); ++sit)
    {
      Vector<FaceIndex> rfaces = a_ebisBox.getFaces(a_vof, rdir, sit());
      int isign = sign(sit());
      for (int iface = 0; iface < rfaces.size(); iface++)
        {
          const FaceIndex& face = rfaces[iface];
          Real faceRad = cellRad + 0.5*isign*a_dx;
          Real faceRad2 = a_dx*Real(face.gridIndex(Side::Hi)[rdir]);
          CH_assert(Abs(faceRad - faceRad2) < 1.0e-9);
          Real areaFrac = a_ebisBox.areaFrac(face);
          openContrib += isign*areaFrac*faceRad*faceRad;
        }
    }
  //add in contribution from irregular face using midpoint rule
  Real bndryArea  = a_ebisBox.bndryArea(a_vof);
  Real rnormal    = a_ebisBox.normal(a_vof)[rdir];
  Real rbndryCent = a_ebisBox.bndryCentroid(a_vof)[rdir];
  Real radCent = cellRad + a_dx*rbndryCent;
  Real bdryContrib= -bndryArea*rnormal*radCent*radCent;

  Real factor = a_dx/(2.*a_cellVol);

  openContrib *=  factor;
  bdryContrib *=  factor;
  a_kvol =  openContrib + bdryContrib;
  //the two is from r^2/2
  //the dx comes from dl
}
/****/

/*****************************/
void
EBArith::
getCompVolRZ(Real&           a_compVol,
             const EBISBox&  a_ebisBox,
             const Real&     a_dx,
             const VolIndex& a_vof,
             bool a_verbose)
{
  int rdir = 0;
  int rindex   =  a_vof.gridIndex()[rdir];
  if (a_verbose)
    {
      pout()     << setprecision(6)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific);
    }
  Real cellRad = a_dx*(Real(rindex) + 0.5);
  Real openContrib = 0.;
  if (a_verbose)
    {
      pout() << endl << "vof = " << a_vof.gridIndex() << endl;
    }
  for (SideIterator sit; sit.ok(); ++sit)
    {
      Vector<FaceIndex> rfaces = a_ebisBox.getFaces(a_vof, rdir, sit());
      int isign = sign(sit());
      if (a_verbose)
        {
          pout()  << "side = " << isign;
        }
      for (int iface = 0; iface < rfaces.size(); iface++)
        {
          const FaceIndex& face = rfaces[iface];
          Real faceRad = cellRad + 0.5*isign*a_dx;
          Real areaFrac = a_ebisBox.areaFrac(face);
          openContrib += isign*areaFrac*faceRad*faceRad;
          if (a_verbose)
            {

              pout() << ", faceRad  = " << faceRad;
              pout() << ", areaFrac = " << areaFrac;
            }
        }
      if (a_verbose)
        {
          pout() << endl;
        }
    }
  //add in contribution from irregular face using midpoint rule
  Real bndryArea  = a_ebisBox.bndryArea(a_vof);
  Real rnormal    = a_ebisBox.normal(a_vof)[rdir];
  Real rbndryCent = a_ebisBox.bndryCentroid(a_vof)[rdir];
  Real radCent    = a_dx*(rindex + 0.5 + rbndryCent);


  Real bdryContrib= -bndryArea*rnormal*radCent*radCent;

  Real factor = a_dx/2.;

  openContrib *=  factor;
  bdryContrib *=  factor;
  a_compVol =  openContrib + bdryContrib;

  if (a_verbose)
    {
      pout() << "bndryArea = " << bndryArea  ;
      pout() << ", rnormal = " << rnormal    ;
      pout() << ", rbndryCent = " << rbndryCent ;
      pout() << ", radCent = " << radCent    ;
      pout() << ", openContrib = " << openContrib  ;
      pout() << ", bdryContrib = " << bdryContrib  ;
      pout() << ", compVol = " << a_compVol  ;
      pout() << endl;
    }
}
/****/

/*****************************/
void
EBArith::
getKVolRZNoDx(Real&           a_kvol,
              Real&           a_cellVol,
              const EBISBox&  a_ebisBox,
              const VolIndex& a_vof)
{
  int rdir = 0;
  int rindex   =  a_vof.gridIndex()[rdir];
  Real cellRad = (Real(rindex) + 0.5);

  a_cellVol =cellRad;

  //kvol = (1/2*cellvol)*int(r^2 nr dl)
  Real openContrib = 0.;
  for (SideIterator sit; sit.ok(); ++sit)
    {
      Vector<FaceIndex> rfaces = a_ebisBox.getFaces(a_vof, rdir, sit());
      int isign = sign(sit());
      for (int iface = 0; iface < rfaces.size(); iface++)
        {
          const FaceIndex& face = rfaces[iface];
          Real faceRad = cellRad + 0.5*Real(isign);
          Real areaFrac = a_ebisBox.areaFrac(face);
          openContrib += isign*areaFrac*faceRad*faceRad;
        }
    }
  //add in contribution from irregular face using midpoint rule
  Real bndryArea  = a_ebisBox.bndryArea(a_vof);
  Real rnormal    = a_ebisBox.normal(a_vof)[rdir];
  Real rbndryCent = a_ebisBox.bndryCentroid(a_vof)[rdir];
  Real radCent = cellRad + rbndryCent;
  Real bdryContrib= -bndryArea*rnormal*radCent*radCent;

  Real factor = 1.0/(2.*a_cellVol);

  openContrib *=  factor;
  bdryContrib *=  factor;
  a_kvol =  openContrib + bdryContrib;
}
/****/
Real
EBArith::norm(Real&                a_volume,
              const EBCellFAB&     a_dataOne,
              const Box&           a_grid,
              const EBISBox&       a_ebisBox,
              const int&           a_comp,
              const int&           a_pval,
              EBNormType::NormMode a_mode)
{
  Real normval, sum;
  RealVect dx = RealVect::Unit;
  //nothing excluded (no finer levels as far as i know)
  IntVectSet ivsExclude;
  volWeightedSum(sum, a_volume, a_dataOne,
                 a_grid, a_ebisBox, dx, ivsExclude, a_comp, a_pval, a_mode);

  if (a_pval == 0)
    {
      normval = sum;
    }
  else if (a_pval == 1)
    {
      if (a_volume > 0.0)
        {
          normval = Abs(sum)/a_volume;
        }
      else
        {
          normval = 0.0;
        }
    }
  else
    {
      Real denom = a_pval;
      Real exponent = 1.0/denom;
      if (a_volume > 0.0)
        {
          sum /= a_volume;
          normval = pow(sum, exponent);
        }
      else
        {
          normval = 0.0;
        }
    }

  return normval;
}

Real
EBArith::norm(const EBCellFAB&     a_dataOne,
              const Box&           a_grid,
              const EBISBox&       a_ebisBox,
              const int&           a_comp,
              const int&           a_pval,
              EBNormType::NormMode a_mode)
{
  Real volume, normval, sum;
  RealVect dx = RealVect::Unit;
  //nothing excluded (no finer levels as far as i know)
  IntVectSet ivsExclude;
  volWeightedSum(sum, volume, a_dataOne,
                 a_grid, a_ebisBox, dx, ivsExclude, a_comp, a_pval, a_mode);

  if (a_pval == 0)
    {
      normval = sum;
    }
  else if (a_pval == 1)
    {
      if (volume > 0.0)
        {
          normval = Abs(sum)/volume;
        }
      else
        {
          normval = 0.0;
        }
    }
  else
    {
      Real denom = a_pval;
      Real exponent = 1.0/denom;
      if (volume > 0.0)
        {
          sum /= volume;
          normval = pow(sum, exponent);
        }
      else
        {
          normval = 0.0;
        }
    }

  return normval;
}

bool
EBArith::getAdjacentFace(FaceIndex& a_adjacentFace,
                         const FaceIndex& a_face,
                         const EBISBox& a_ebisBox,
                         const ProblemDomain& a_domainBox,
                         const int& a_idir,
                         const Side::LoHiSide& a_side)
{
  bool uniqueFace = true;
  int faceDir = a_face.direction();
  VolIndex loFaceVoF = a_face.getVoF(Side::Lo);
  VolIndex hiFaceVoF = a_face.getVoF(Side::Hi);
  //figure out if we have the connectivity to find a face
  //with which to interpolate.  If so,
  //find the other face
  if (!a_face.isBoundary())
    {
      Vector<FaceIndex> loTanFaces = a_ebisBox.getFaces(loFaceVoF, a_idir, a_side);
      Vector<FaceIndex> hiTanFaces = a_ebisBox.getFaces(hiFaceVoF, a_idir, a_side);
      if ((loTanFaces.size() != 1) || (hiTanFaces.size() != 1))
        {
          uniqueFace = false;
        }
      else if ((loTanFaces[0].isBoundary()) && (hiTanFaces[0].isBoundary()))
        {
          uniqueFace = false;
        }
      else
        {
          const FaceIndex& loTanFace = loTanFaces[0];
          const FaceIndex& hiTanFace = hiTanFaces[0];
          VolIndex loOtherVoF = loTanFace.getVoF(a_side);
          VolIndex hiOtherVoF = hiTanFace.getVoF(a_side);
          if (!(a_ebisBox.isConnected(loOtherVoF, hiOtherVoF)))
            {
              uniqueFace = false;
            }
          if (uniqueFace)
            {
              Vector<FaceIndex> otherFaces = a_ebisBox.getFaces(loOtherVoF, faceDir, Side::Hi);
              if (otherFaces.size() != 1)
                {
                  uniqueFace = false;
                }
              else
                {
                  a_adjacentFace = otherFaces[0];
                }
            }
        }
    } //end if !face.isBoundary
  else
    {
      //boundary face.
      IntVect loVoFiv = loFaceVoF.gridIndex();
      IntVect hiVoFiv = hiFaceVoF.gridIndex();
      IntVect loDomiv = a_domainBox.domainBox().smallEnd();
      IntVect hiDomiv = a_domainBox.domainBox().bigEnd();
      Vector<FaceIndex> otherFaces;

      if (hiVoFiv[faceDir] == loDomiv[faceDir])
        {
          Vector<FaceIndex> hiTanFaces = a_ebisBox.getFaces(hiFaceVoF, a_idir, a_side);
          if (hiTanFaces.size() != 1)
            {
              uniqueFace = false;
            }
          else if (hiTanFaces[0].isBoundary())
            {
              //in 3d,  can be a boundary in two directions
              uniqueFace = false;
            }
          else
            {
              VolIndex hiOtherVoF = hiTanFaces[0].getVoF(a_side);
              Vector<FaceIndex> otherFaces = a_ebisBox.getFaces(hiOtherVoF, faceDir, Side::Lo);
              if (otherFaces.size() != 1)
                {
                  uniqueFace = false;
                }
              else
                {
                  a_adjacentFace = otherFaces[0];
                }
            }

        }
      else if (loVoFiv[faceDir] == hiDomiv[faceDir])
        {
          Vector<FaceIndex> loTanFaces = a_ebisBox.getFaces(loFaceVoF, a_idir, a_side);
          if (loTanFaces.size() != 1)
            {
              uniqueFace = false;
            }
          else if (loTanFaces[0].isBoundary())
            {
              //in 3d,  can be a boundary in two directions
              uniqueFace = false;
            }
          else
            {
              VolIndex loOtherVoF = loTanFaces[0].getVoF(a_side);
              Vector<FaceIndex> otherFaces = a_ebisBox.getFaces(loOtherVoF, faceDir, Side::Hi);
              if (otherFaces.size() != 1)
                {
                  uniqueFace = false;
                }
              else
                {
                  a_adjacentFace = otherFaces[0];
                }
            }

        }
      else
        {
          MayDay::Error("boundary face confusion");
        }

    }
  return uniqueFace;
}

/**************/
//If the face is covered, one of the vofs
//should have been covered which means
//i should not have gotten a face
//so I won't check that.
//If a vof is outside the domain, then
//we have a boundary face an those don't
//get done here.
// check conditions for using a nontrivial O(h**2) flux stencil.
// The conditions are:
// (1) Is this face irregular ?
// (2) If so, are the components of the boundary normal in the directions
//      tangent to the face of the same sign ?
// (3) If (1) or (2) are true, are the faces that will be used to
//      computed the centered flues regular, and are the vofs making up
//      those faces connected appropriately to the vofs making up
//      face at which the flux is being evaluated ?
// If any of these conditions fail, then we use the centered difference
// flux. This corresponds to dropping order for some of the
// irregular faces, but is still second order accurate for regular faces
// mumble set of codimension two mumble mumble
/**************/
void
EBArith::computeInterpStencil(FaceStencil& a_thisStencil,
                              const FaceIndex& a_thisFace,
                              const EBISBox& a_ebisBox,
                              const ProblemDomain& a_domain,
                              const int& a_dir)
{
  a_thisStencil.clear();
  const IntVect& ivlo = a_thisFace.gridIndex(Side::Lo);
  const IntVect& ivhi = a_thisFace.gridIndex(Side::Hi);

  //boundary faces are done separately
  if (!a_domain.contains(ivlo) || !a_domain.contains(ivhi))
    return;

  IntVect kperp = IntVect::Zero;
  for (int d = 0; d < SpaceDim-1; ++d)
    {
      kperp[d] = (a_dir + d + 1) % SpaceDim;
    }

  bool MakeSimpleStencil = true;
  Real areaFrac = a_ebisBox.areaFrac(a_thisFace);
  if (areaFrac < 1.0)
    {
      // Condition (1) is met.
      VolIndex vofLo = a_thisFace.getVoF(Side::Lo);
      VolIndex vofHi = a_thisFace.getVoF(Side::Hi);

      RealVect rvnormLo = a_ebisBox.normal(vofLo);
      RealVect rvnormHi = a_ebisBox.normal(vofHi);
      //this checks that the normals are monotonic
      //in the non-facedir direction
      bool montonNorm = true;
      bool facesOk = true;
      for (int d = 0; d < SpaceDim-1; ++d)
        {
          montonNorm = montonNorm
            && ( rvnormLo[d]*rvnormHi[d] >= 0.0 );
        }
      if (montonNorm)
        {
          // Condition (2) is met.
          //find perpendicular directions and signs
          RealVect rvnorm = RealVect::Zero;
          IntVect  ivsign =  IntVect::Zero;
          Tuple<IntVect,SpaceDim-1> ivPerp;
          for (int d = 0; d < SpaceDim-1; ++d)
            {
              rvnorm[d] = 0.5 * (rvnormLo[d] + rvnormHi[d]);
              if (rvnorm[d] > 0.0)
                ivsign[d] = 1;
              else
                ivsign[d] = -1;

              ivPerp[d] = IntVect::Zero;
              ivPerp[d].setVal(kperp[d],ivsign[d]);
              //by convention
              const IntVect& ivface = vofHi.gridIndex();
              ivPerp[d] += ivface;
              //this stuff  checks to see
              //if there are vofs in the perpendicular direction
              //that are regular.   if so, we can do the johansen
              //interpolation between faces thing without fear.
              //If there is no regular face nearby, this gets
              //godawful complicated and we just drop to the
              //simple stencil.
              if (ivsign[d] != 0)
                {
                  const IntVect& perpIV = ivPerp[d];
                  IntVect othrIV = perpIV - BASISV(a_dir);
                  if (!a_ebisBox.getRegion().contains(perpIV) ||
                     !a_ebisBox.getRegion().contains(othrIV))
                    {
                      facesOk = false;
                    }
                  else if (a_ebisBox.isCovered(perpIV) ||
                          a_ebisBox.isCovered(othrIV))
                    {
                      facesOk = false;
                    }
                  else
                    {
                      facesOk = facesOk && (a_ebisBox.isRegular(perpIV) ||
                                            a_ebisBox.isRegular(othrIV));
                    }
                }
            }
          if (facesOk)
            {
              // Get vofs for faces at ivPerp to check connectivity.
              Tuple<FaceIndex,SpaceDim-1> perpFaces;
              for (int d = 0; d < SpaceDim-1; ++d)
                {
                  if (ivsign[d] != 0)
                    {
                      perpFaces[d] =
                        a_ebisBox.getAllFaces(ivPerp[d],a_dir,Side::Lo)[0];
                      facesOk = facesOk &&
                        a_ebisBox.isConnected(perpFaces[d].getVoF(Side::Lo),vofLo);
                      facesOk = facesOk &&
                        a_ebisBox.isConnected(perpFaces[d].getVoF(Side::Hi),vofHi);
                    }
                }
              // Condition (3) is met - construct the stencil.
              // Compute the center of area for the face.
              RealVect faceCent = a_ebisBox.centroid(a_thisFace);
              Real thisCoef = 1.0;
              for (int d = 0; d < SpaceDim-1; ++d)
                {
                  thisCoef -= Abs(faceCent[kperp[d]]);
                }
              a_thisStencil.add(a_thisFace, thisCoef);

              for (int d = 0; d < SpaceDim-1; ++d)
                {
                  Real thisXC = Abs(faceCent[kperp[d]]);
                  if (ivsign[d] != 0)
                    {
                      a_thisStencil.add(perpFaces[d], thisXC);
                    }
                } //end loop over directions
              MakeSimpleStencil = false;
            } // end if (facesOK) (all conditons met
        } //end if (perp norms are monotonic)
    } //end isFaceIrregular (if (areaFrac <1))

  if (MakeSimpleStencil)
    {
      a_thisStencil.add(a_thisFace, 1.0);
    }
}
/**************/
//If the face is covered, one of the vofs
//should have been covered which means
//i should not have gotten a face
//so I won't check that.
//If a vof is outside the domain, then
//we have a boundary face an those don't
//get done here.
// check conditions for using a nontrivial O(h**2) flux stencil.
// The conditions are:
// (1) Is this face irregular ?
// (2) If so, are the components of the boundary normal in the directions
//      tangent to the face of the same sign ?
// (3) If (1) or (2) are true, are the faces that will be used to
//      computed the centered flues regular, and are the vofs making up
//      those faces connected appropriately to the vofs making up
//      face at which the flux is being evaluated ?
// If any of these conditions fail, then we use the centered difference
// flux. This corresponds to dropping order for some of the
// irregular faces, but is still second order accurate for regular faces
// mumble set of codimension two mumble mumble
/**************/
void
EBArith::computeGradFluxStencil(VoFStencil& a_thisStencil,
                                const FaceIndex& a_thisFace,
                                const EBISBox& a_ebisBox,
                                const ProblemDomain& a_domain,
                                const int& a_dir)
{

  a_thisStencil.clear();
  const IntVect& ivlo = a_thisFace.gridIndex(Side::Lo);
  const IntVect& ivhi = a_thisFace.gridIndex(Side::Hi);

  //boundary faces are done separately
  if (!a_domain.contains(ivlo) || !a_domain.contains(ivhi))
    return;

  IntVect kperp = IntVect::Zero;
  for (int d = 0; d < SpaceDim-1; ++d)
    {
      kperp[d] = (a_dir + d + 1) % SpaceDim;
    }

  bool MakeSimpleStencil = true;
  Real areaFrac = a_ebisBox.areaFrac(a_thisFace);
  if (areaFrac < 1.0)
    {
      // Condition (1) is met.
      VolIndex vofLo = a_thisFace.getVoF(Side::Lo);
      VolIndex vofHi = a_thisFace.getVoF(Side::Hi);

      RealVect rvnormLo = a_ebisBox.normal(vofLo);
      RealVect rvnormHi = a_ebisBox.normal(vofHi);
      //this cheecks that the normals are monotonic
      //in the non-facedir direction
      bool montonNorm = true;
      bool facesOk = true;
      for (int d = 0; d < SpaceDim-1; ++d)
        {
          montonNorm = montonNorm
            && ( rvnormLo[d]*rvnormHi[d] >= 0.0 );
        }
      if (montonNorm)
        {
          // Condition (2) is met.
          //find perpendicular directions and signs
          RealVect rvnorm = RealVect::Zero;
          IntVect  ivsign =  IntVect::Zero;
          Tuple<IntVect,SpaceDim-1> ivPerp;
          for (int d = 0; d < SpaceDim-1; ++d)
            {
              rvnorm[d] = 0.5 * (rvnormLo[d] + rvnormHi[d]);
              if (rvnorm[d] > 0.0)
                ivsign[d] = 1;
              else
                ivsign[d] = -1;

              ivPerp[d] = IntVect::Zero;
              ivPerp[d].setVal(kperp[d],ivsign[d]);
              //by convention
              const IntVect& ivface = vofHi.gridIndex();
              ivPerp[d] += ivface;
              //this stuff  checks to see
              //if there are vofs in the perpendicular direction
              //that are regular.   if so, we can do the johansen
              //interpolation between faces thing without fear.
              //If there is no regular face nearby, this gets
              //godawful complicated and we just drop to the
              //simple stencil.
              if (ivsign[d] != 0)
                {
                  const IntVect& perpIV = ivPerp[d];
                  IntVect othrIV = perpIV - BASISV(a_dir);
                  if (!a_ebisBox.getRegion().contains(perpIV) ||
                     !a_ebisBox.getRegion().contains(othrIV))
                    {
                      facesOk = false;
                    }
                  else if (a_ebisBox.isCovered(perpIV) ||
                          a_ebisBox.isCovered(othrIV))
                    {
                      facesOk = false;
                    }
                  else
                    {
                      facesOk = facesOk && (a_ebisBox.isRegular(perpIV) ||
                                            a_ebisBox.isRegular(othrIV));
                    }
                }
            }
          if (facesOk)
            {
              // Get vofs for faces at ivPerp to check connectivity.
              Tuple<FaceIndex,SpaceDim-1> perpFaces;
              for (int d = 0; d < SpaceDim-1; ++d)
                {
                  if (ivsign[d] != 0)
                    {
                      perpFaces[d] =
                        a_ebisBox.getAllFaces(ivPerp[d],a_dir,Side::Lo)[0];
                      facesOk = facesOk &&
                        a_ebisBox.isConnected(perpFaces[d].getVoF(Side::Lo),vofLo);
                      facesOk = facesOk &&
                        a_ebisBox.isConnected(perpFaces[d].getVoF(Side::Hi),vofHi);
                    }
                }
              // Condition (3) is met - construct the stencil.
              // Compute the center of area for the face.
              RealVect faceCent = a_ebisBox.centroid(a_thisFace);
              Real thisCoef = 1.0;
              for (int d = 0; d < SpaceDim-1; ++d)
                {
                  thisCoef -= Abs(faceCent[kperp[d]]);
                }
              a_thisStencil.add(a_thisFace.getVoF(Side::Lo), -thisCoef);
              a_thisStencil.add(a_thisFace.getVoF(Side::Hi), thisCoef);
              for (int d = 0; d < SpaceDim-1; ++d)
                {
                  if (ivsign[d] != 0)
                    {
                      Real thisXC = Abs(faceCent[kperp[d]]);

                      for (SideIterator jsit; jsit.ok(); ++jsit)
                        {
                          Real coef = Real(sign(jsit()))*thisXC;
                          const VolIndex& vof = perpFaces[d].getVoF(jsit());
                          a_thisStencil.add(vof, coef);
                        }
                    }
                } //end loop over directions
              MakeSimpleStencil = false;
            } // end if (facesOK) (all conditons met
        } //end if (perp norms are monotonic)
    } //end isFaceIrregular (if (areaFrac <1))
  if (MakeSimpleStencil)
    {
      // We don't need / have the more elaborate stencil
      // - instead use simple centered differences.
      for (SideIterator jsit; jsit.ok(); ++jsit)
        {
          Real weight = Real(sign(jsit()));
          VolIndex vof = a_thisFace.getVoF(jsit());
          a_thisStencil.add(vof, weight);
        }
    } // End stencil construction for simple face.
}

void
EBArith::volWeightedSum(Real&                a_norm,
                        Real&                a_volume,
                        const EBCellFAB&     a_src,
                        const Box&           a_region,
                        const EBISBox&       a_ebisBox,
                        const RealVect&      a_dx,
                        const IntVectSet&    a_ivsExclude,
                        const int&           a_comp,
                        const int&           a_pval,
                        EBNormType::NormMode a_mode)
{
  CH_TIME("EBArith::volWeightedSum");
  Real cellVol = 1.0;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      cellVol *= a_dx[idir];
    }

  CH_assert(a_pval >= -2);
  CH_assert(a_comp >= 0);
  CH_assert(a_comp < a_src.nComp());
  CH_assert(a_src.getRegion().contains(a_region));
  CH_assert(a_ebisBox.getRegion().contains(a_region));
  int idoirr = 0;
  int idoreg = 0;
  if (a_mode == EBNormType::OverBoth)
    {
      idoirr = 1;
      idoreg = 1;
    }
  else if (a_mode == EBNormType::OverOnlyRegular)
    {
      idoirr = 0;
      idoreg = 1;
    }
  else if (a_mode == EBNormType::OverOnlyIrregular)
    {
      idoirr = 1;
      idoreg = 0;
    }

  a_norm = 0.;
  a_volume = 0.;

  IntVectSet ivsRegion(a_region);
  ivsRegion -= a_ivsExclude;
  for (VoFIterator vofit(ivsRegion, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      const IntVect& iv = vof.gridIndex();
      if (((idoirr == 1) && a_ebisBox.isIrregular(iv)) ||
         ((idoreg == 1) && a_ebisBox.isRegular(iv)))
        {
          Real val     = Abs(a_src(vof,a_comp));
          if (Abs(val) >s_valMax)
            {
              s_valMax = Abs(val);
              s_vofMax = vof;
            }
          Real volfrac =  a_ebisBox.volFrac(vof);

          a_volume += cellVol*volfrac;
          if (a_pval == 0) // max a_norm
            {
              if (volfrac > 0)
                {
                  a_norm = Max(a_norm, val);
                }
            }
          else if (a_pval == 1) //L1 norm
            {
              a_norm += Abs(cellVol*volfrac*val);
            }
          else if (a_pval == -1) //just doing integral of a_src  where a_src
                                //is already multiplied by volume fraction
            {
              Real noAbsVal = a_src(vof,a_comp);
              a_norm += cellVol*noAbsVal;
            }
          else if (a_pval == -2) //just doing integral of a_src  where a_src
                                //is NOT already multiplied by volume fraction
            {
              Real noAbsVal = a_src(vof,a_comp);
              a_norm += volfrac*cellVol*noAbsVal;
            }
          else // a_pval > 1
            {
              Real rpval = a_pval;
              Real integrand= pow(val, rpval);
              a_norm += cellVol*volfrac*integrand;
            }
        }
    }
}


void
EBArith::computeSum(Real&                a_norm,
                    Real&                a_volume,
                    const EBCellFAB&     a_src,
                    const Box&           a_region,
                    const EBISBox&       a_ebisBox,
                    const RealVect&      a_dx,
                    const IntVectSet&    a_ivsExclude,
                    const int&           a_comp,
                    EBNormType::NormMode a_mode)
{
  Real cellVol = 1.0;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      cellVol *= a_dx[idir];
    }

  CH_assert(a_comp >= 0);
  CH_assert(a_comp < a_src.nComp());
  CH_assert(a_src.getRegion().contains(a_region));
  CH_assert(a_ebisBox.getRegion().contains(a_region));
  int idoirr = 0;
  int idoreg = 0;
  if (a_mode == EBNormType::OverBoth)
    {
      idoirr = 1;
      idoreg = 1;
    }
  else if (a_mode == EBNormType::OverOnlyRegular)
    {
      idoirr = 0;
      idoreg = 1;
    }
  else if (a_mode == EBNormType::OverOnlyIrregular)
    {
      idoirr = 1;
      idoreg = 0;
    }

  a_norm = 0.;
  a_volume = 0.;

  IntVectSet ivsRegion(a_region);
  ivsRegion -= a_ivsExclude;
  for (VoFIterator vofit(ivsRegion, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      const IntVect& iv = vof.gridIndex();
      if (((idoirr == 1) && a_ebisBox.isIrregular(iv)) ||
         ((idoreg == 1) && a_ebisBox.isRegular(iv)))
        {
          Real val     = a_src(vof,a_comp);

          Real volfrac =  a_ebisBox.volFrac(vof);

          a_volume += cellVol*volfrac;
          a_norm += cellVol*volfrac*val;

        }
    }
}



void
EBArith::computeUnweightedSum(Real&                a_norm,
                              Real&                a_volume,
                              const EBCellFAB&     a_src,
                              const Box&           a_region,
                              const EBISBox&       a_ebisBox,
                              const RealVect&      a_dx,
                              const IntVectSet&    a_ivsExclude,
                              const int&           a_comp,
                              EBNormType::NormMode a_mode)
{
  Real cellVol = 1.0;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      cellVol *= a_dx[idir];
    }

  CH_assert(a_comp >= 0);
  CH_assert(a_comp < a_src.nComp());
  CH_assert(a_src.getRegion().contains(a_region));
  CH_assert(a_ebisBox.getRegion().contains(a_region));
  int idoirr = 0;
  int idoreg = 0;
  if (a_mode == EBNormType::OverBoth)
    {
      idoirr = 1;
      idoreg = 1;
    }
  else if (a_mode == EBNormType::OverOnlyRegular)
    {
      idoirr = 0;
      idoreg = 1;
    }
  else if (a_mode == EBNormType::OverOnlyIrregular)
    {
      idoirr = 1;
      idoreg = 0;
    }

  a_norm = 0.;
  a_volume = 0.;

  IntVectSet ivsRegion(a_region);
  ivsRegion -= a_ivsExclude;
  for (VoFIterator vofit(ivsRegion, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      const IntVect& iv = vof.gridIndex();
      if (((idoirr == 1) && a_ebisBox.isIrregular(iv)) ||
         ((idoreg == 1) && a_ebisBox.isRegular(iv)))
        {
          Real val     = a_src(vof,a_comp);

          Real volfrac =  a_ebisBox.volFrac(vof);

          a_volume += cellVol*volfrac;
          a_norm += val;

        }
    }
}


Real
EBArith::dotProduct(const EBCellFAB& a_dataOne,
                    const EBCellFAB& a_dataTwo,
                    const Box& a_region,
                    const EBISBox& a_ebisBox,
                    const int& a_comp)
{
  //check everything in sight
  CH_assert(a_comp >= 0);
  CH_assert(a_comp < a_dataOne.nComp());
  CH_assert(a_comp < a_dataTwo.nComp());
  CH_assert(a_dataOne.getRegion().contains(a_region));
  CH_assert(a_dataTwo.getRegion().contains(a_region));
  CH_assert(a_ebisBox.getRegion().contains(a_region));

  Real dotgrid = 0.0;
  IntVectSet ivsregion(a_region);
  for (VoFIterator vofit(ivsregion, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      const Real& valone = a_dataOne(vof,a_comp);
      const Real& valtwo = a_dataTwo(vof,a_comp);
      dotgrid += valone*valtwo;
    }
  return dotgrid;
}
/***************/
/***************/
Real
EBArith::dotProduct(const BoxLayoutData<EBCellFAB >& a_dataOne,
                    const BoxLayoutData<EBCellFAB >& a_dataTwo,
                    const BoxLayout& a_dblIn,
                    const EBISLayout& a_ebisl,
                    const int& a_comp)
{
  BaseEBCellFactory<Real> testFactory(a_ebisl);
  Real rhodot = 0.0;
  //calculate the single-processor dot product
  //lots of assertions  have to hold
  //(same number of vars, grid contained in all the fabs etc)
  DataIterator dit = a_dataOne.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      rhodot += dotProduct(a_dataOne[dit()], a_dataTwo[dit()],
                           a_dblIn.get(dit()), a_ebisl[dit()], a_comp);
    }

  // now for the multi-processor fandango
  //gather all the rhodots onto a vector and add them up
  int baseProc = 0;
  Vector<Real> dotVec;
  gather(dotVec, rhodot, baseProc);

  Real rhodotTot = 0.0;
  if (procID() == baseProc)
    {
      CH_assert(dotVec.size() == numProc());
      for (int ivec = 0; ivec < dotVec.size(); ivec++)
        {
          rhodotTot += dotVec[ivec];
        }
    }
  //broadcast the sum to all processors.
  broadcast(rhodotTot, baseProc);

  //return the total
  return rhodotTot;
}
/***************/
/*****************************/
void
EBArith::
getInterpStencil2D(FaceStencil&         a_sten,
                   const FaceIndex&     a_face,
                   const IntVectSet&    a_coarseFineIVS,
                   const EBISBox&       a_ebisBox,
                   const ProblemDomain& a_domainBox)
{
  CH_assert(SpaceDim == 2);
  a_sten.clear();
  RealVect faceCentroid = a_ebisBox.centroid(a_face);
  int faceDir = a_face.direction();
  int tanDir = 1 - faceDir;
  Side::LoHiSide tanFaceSide;
  if (faceCentroid[tanDir] > 0.)
    {
      tanFaceSide = Side::Hi;
    }
  else if (faceCentroid[tanDir] < 0.)
    {
      tanFaceSide = Side::Lo;
    }
  else
    {
      tanFaceSide = Side::Lo;
      // MayDay::Error("EBArith: getInterpStencil2D faceCentroid[tanDir] = 0");
    }

  FaceIndex otherFace;
  bool uniqueFace = EBArith::getAdjacentFace(otherFace, a_face,
                                             a_ebisBox, a_domainBox,
                                             tanDir, tanFaceSide);

  //drop order of interpolation to zero if the other face is not unique
  //or if this face or the otherface is on the coarse fine interface
  bool dropOrder = false;
  dropOrder = !uniqueFace;
  if (!dropOrder)
    {
      for (SideIterator sit; sit.ok(); ++sit)
        {
          if (a_coarseFineIVS.contains(   a_face.gridIndex(sit())) ||
             a_coarseFineIVS.contains(otherFace.gridIndex(sit())))
            {
              dropOrder = true;
            }
        }
    }

  if (dropOrder)
    {
      a_sten.add(a_face, 1.0);
    }
  else
    {
      CH_assert(a_face.isDefined());
      CH_assert(otherFace.isDefined());
      Real dist = Abs(faceCentroid[tanDir]);
      Real thisWeight = 1.0 - dist;
      Real otherWeight =  dist;
      a_sten.add(a_face, thisWeight);
      a_sten.add(otherFace, otherWeight);
    }
}
/*****************************/
void
EBArith::
getInterpStencil3D(FaceStencil&         a_sten,
                   const FaceIndex&     a_face,
                   const IntVectSet&    a_coarseFineIVS,
                   const EBISBox&       a_ebisBox,
                   const ProblemDomain& a_domainBox)
{
  CH_assert(SpaceDim == 3);
  a_sten.clear();

  int faceDir = a_face.direction();
  RealVect faceCentroid = a_ebisBox.centroid(a_face);
  //first check to see if this is a genuinely 3d face
  //centroid > 0 in all directions
  bool really3D = true;
  int izerod = 7;
  Real tol = PolyGeom::getTolerance();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (idir != faceDir)
        {
          if (Abs(faceCentroid[idir]) < tol)
            {
              really3D= false;
              izerod = idir;
            }
        }
    }
  if (!really3D)
    {
      //use 2D stencil because we have a centroid that is zero
      //in the given direction
      int tanDir = 3 - izerod - faceDir;
      Side::LoHiSide tanFaceSide;
      if (faceCentroid[tanDir] > 0)
        tanFaceSide = Side::Hi;
      else
        tanFaceSide = Side::Lo;

      bool dropOrder = false;
      FaceIndex otherFace;
      bool uniqueFace = EBArith::getAdjacentFace(otherFace, a_face,
                                                 a_ebisBox, a_domainBox,
                                                 tanDir, tanFaceSide);
      dropOrder = !uniqueFace;

      if ( !dropOrder )
        {
          for (SideIterator sit; sit.ok(); ++sit)
            {
              if (a_coarseFineIVS.contains(   a_face.gridIndex(sit())) ||
                 a_coarseFineIVS.contains(otherFace.gridIndex(sit())))
                {
                  dropOrder = true;
                }
            }
        }

      if (dropOrder)
        {
          a_sten.add(a_face, 1.0);
        }
      else
        {
          CH_assert(a_face.isDefined());
          CH_assert(otherFace.isDefined());
          Real dist = Abs(faceCentroid[tanDir]);
          Real thisWeight = 1 - dist;
          Real otherWeight =  dist;
          a_sten.add(a_face, thisWeight);
          a_sten.add(otherFace, otherWeight);
        }
    }
  else
    {
      int tanDirs[2];
      Side::LoHiSide tanFaceSides[2];
      {
        int itan = 0;
        for (int idir = 0; idir < 3; idir++)
          {
            if (idir != faceDir)
              {
                tanDirs[itan] = idir;
                if (faceCentroid[tanDirs[itan]] > 0)
                  tanFaceSides[itan] = Side::Hi;
                else
                  tanFaceSides[itan] = Side::Lo;
                itan++;
              }
          }
      }

      bool dropOrder = false;
      FaceIndex face01, face10, face11;
      int ixDir = tanDirs[0];
      int iyDir = tanDirs[1];
      Side::LoHiSide xSide = tanFaceSides[0];
      Side::LoHiSide ySide = tanFaceSides[1];
      //figure out whether to drop order.
      //if not,  get the faces involved
      VolIndex vofLo = a_face.getVoF(Side::Lo);
      VolIndex vofHi = a_face.getVoF(Side::Hi);

      //first get face10 which is the in
      //the  ixdir direction from the input face
      if (!dropOrder)
        {
          bool uniqueFace10
            = EBArith::getAdjacentFace(face10, a_face, a_ebisBox, a_domainBox, ixDir, xSide);
          dropOrder = !uniqueFace10;
          //now get face11 which is the in
          //the  ixdir,iydir direction from the input face (diagonal)
          if (uniqueFace10)
            {
              bool uniqueFace11
                = EBArith::getAdjacentFace(face11, face10, a_ebisBox, a_domainBox, iyDir, ySide);
              dropOrder = !uniqueFace11;
            }
        }
      //the  ixdir,iydir direction from the input face (diagonal)
      //now get face01 which is the in
      //the  iydir direction from the input face
      if (!dropOrder)
        {
          //first get face01 which is the in
          //the  iydir direction from the input face
          bool uniqueFace01
            = EBArith::getAdjacentFace(face01, a_face, a_ebisBox, a_domainBox, iyDir, ySide);
          dropOrder = !uniqueFace01;
          //now get face11 which is the in
          //the  ixdir,iydir direction from the input face (diagonal)
          //compute temp face and see if it is the same as
          //the one computed from the other direction.
          //if not , drop order
          FaceIndex face11temp;
          if (uniqueFace01)
            {
              bool uniqueFace11
                = EBArith::getAdjacentFace(face11temp, face01, a_ebisBox, a_domainBox, ixDir, xSide);
              dropOrder = !uniqueFace11;
              if ((!dropOrder) && !(face11temp == face11) )
                {
                  dropOrder = true;
                }
            }

          //finally if any of the stencil faces are in the coarse-fine interface, drop order
          if (!dropOrder)
            {
              for (SideIterator sit; sit.ok(); ++sit)
                {
                  if ((a_coarseFineIVS.contains(a_face.gridIndex(sit()))) ||
                     (a_coarseFineIVS.contains(face01.gridIndex(sit()))) ||
                     (a_coarseFineIVS.contains(face10.gridIndex(sit()))) ||
                     (a_coarseFineIVS.contains(face11.gridIndex(sit()))))
                    {
                      dropOrder = true;
                    }
                }
            }
        }
      ///////////////////////
      //construct the stencils
      ///////////////////////
      if (dropOrder)
        {
          a_sten.add(a_face, 1.0);
        }
      else
        {
          FaceIndex face00 = a_face;
          Real xbar = Abs(faceCentroid[ixDir]);
          Real ybar = Abs(faceCentroid[iyDir]);
          Real f00coef = 1.0 - xbar - ybar + xbar*ybar;
          Real f10coef = xbar - xbar*ybar;
          Real f01coef = ybar - xbar*ybar;
          Real f11coef = xbar*ybar;

          for (SideIterator sit; sit.ok(); ++sit)
            {
              if (!face00.isBoundary())
                CH_assert(face00.cellIndex(sit()) >= 0);
              if (!face01.isBoundary())
                CH_assert(face01.cellIndex(sit()) >= 0);
              if (!face10.isBoundary())
                CH_assert(face10.cellIndex(sit()) >= 0);
              if (!face11.isBoundary())
                CH_assert(face11.cellIndex(sit()) >= 0);
            }
          CH_assert(face00.isDefined());
          CH_assert(face01.isDefined());
          CH_assert(face10.isDefined());
          CH_assert(face11.isDefined());

          a_sten.add(face00, f00coef);
          a_sten.add(face01, f01coef);
          a_sten.add(face10, f10coef);
          a_sten.add(face11, f11coef);
        }
    }
}
/*****************************/
void
EBArith::
computeCoveredFaces(Vector<VolIndex>&     a_coveredFace,
                    IntVectSet&           a_coveredSets,
                    IntVectSet&           a_irregIVS,
                    const int&            a_idir,
                    const Side::LoHiSide& a_sd,
                    const EBISBox&        a_ebisBox,
                    const Box&            a_region)
{
  //first compute the sets where where the covered faces exist
  //start with all irregular cells and subtract off  cells
  //whose vofs all have faces in the given direction
  a_irregIVS =  a_ebisBox.getIrregIVS(a_region);
  a_coveredSets = a_irregIVS;
  a_coveredFace.resize(0);
  for (IVSIterator ivsit(a_irregIVS); ivsit.ok(); ++ivsit)
    {
      const IntVect& iv = ivsit();
      Vector<VolIndex> vofs = a_ebisBox.getVoFs(iv);
      bool allVoFsHaveFaces = true;
      for (int ivof = 0; ivof < vofs.size(); ivof++)
        {
          const VolIndex& vof = vofs[ivof];
          Vector<FaceIndex> faces = a_ebisBox.getFaces(vof, a_idir, a_sd);
          if (faces.size() == 0)
            {
              allVoFsHaveFaces = false;
              a_coveredFace.push_back(vof);
            }
        }
      if (allVoFsHaveFaces)
        a_coveredSets -= iv;
    }
}
/*****************************/
FaceStencil
EBArith::
getInterpStencil(const FaceIndex&     a_face,
                 const IntVectSet&    a_cfivs,
                 const EBISBox&       a_ebisBox,
                 const ProblemDomain& a_domainBox)
{
  FaceStencil sten;
#if CH_SPACEDIM==2
  getInterpStencil2D(sten, a_face, a_cfivs, a_ebisBox, a_domainBox);

#elif CH_SPACEDIM==3
  getInterpStencil3D(sten, a_face, a_cfivs, a_ebisBox, a_domainBox);
#else
  bogus_ch_spacedim_macro();
#endif
  //EBArith::computeInterpStencil(sten, a_face, a_ebisBox, a_domainBox, a_face.direction());
  return sten;
}
/***************/
void
EBArith::meanOverHierarchy(Real&                                   a_mean,
                           const Vector<LevelData<EBCellFAB> *>&   a_divu,
                           const Vector<DisjointBoxLayout>&        a_grids,
                           const Vector<EBISLayout>&               a_ebisl,
                           const Vector<int>&                      a_refRat,
                           const ProblemDomain&                    a_domainCoarsest,
                           const RealVect&                         a_dxCoarsest,
                           const int& a_comp,
                           int a_pval)
{
  Real sumTotal = 0;
  Real volTotal = 0;
  RealVect dxLev = a_dxCoarsest;
  CH_assert((a_pval == -1) || (a_pval == -2));
  int pval= a_pval;
  for (int ilev = 0; ilev < a_divu.size(); ilev++)
    {
      Real sumLevel, volLevel;
      IntVectSet ivsExclude;
      if (ilev < a_divu.size()-1)
        {
          for (LayoutIterator lit = a_grids[ilev+1].layoutIterator(); lit.ok(); ++lit)
            {
              ivsExclude |= a_grids[ilev+1].get(lit());
            }
        }

      volWeightedSum(sumLevel, volLevel, *a_divu[ilev], a_grids[ilev], a_ebisl[ilev],
                     dxLev, ivsExclude, a_comp, pval, EBNormType::OverBoth);

      sumTotal += sumLevel;
      volTotal += volLevel;
      dxLev /= a_refRat[ilev];
    }
  if (volTotal > 0)
    {
      a_mean = sumTotal/volTotal;
    }
  else
    {
      a_mean = 0;
    }
}

void
EBArith::volWeightedSum(Real&                            a_norm,
                        Real&                            a_volume,
                        const BoxLayoutData<EBCellFAB >& a_src,
                        const BoxLayout&                 a_region,
                        const EBISLayout&                a_ebisl,
                        const RealVect&                  a_dx,
                        const IntVectSet&                a_ivsExclude,
                        const int&                       a_comp,
                        const int&                       a_pval,
                        EBNormType::NormMode             a_mode)
{
  CH_assert(a_comp >= 0);
  CH_assert(a_pval >= -2);
  CH_assert(a_comp < a_src.nComp());

  Real normlocal = 0.;
  Real vollocal =  0.;
  DataIterator dit = a_src.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      Real volgrid = 0.;
      Real normgrid = 0.;
      const Box& layoutRegion =  a_region.get(dit());
      const Box& srcRegion    =  a_src[dit()].getRegion();
      const EBCellFAB& srcfab = a_src[dit()];
      if (!srcRegion.contains(layoutRegion))
        {
          MayDay::Error("Incompatible layouts in EBArith::volWeightedSum");
        }

      volWeightedSum(normgrid, volgrid,
                     srcfab,
                     a_region.get(dit()), a_ebisl[dit()], a_dx, a_ivsExclude,
                     a_comp, a_pval, a_mode);

      if (volgrid > 0.)
        {
          vollocal += volgrid;
          if (a_pval == 0)
            {
              normlocal = Max(normlocal, normgrid);
            }
          else
            {
              normlocal += normgrid;
            }
        }
    }
  Real volTot = 0.0;
  Real normTot = 0.0;

  int baseProc = 0;
  Vector<Real> normVec;
  Vector<Real> volVec;
  gather(normVec, normlocal, baseProc);
  gather( volVec,  vollocal, baseProc);

  if (procID() == baseProc)
    {
      CH_assert(volVec.size() == numProc());
      CH_assert(normVec.size() == numProc());
      for (int ivec = 0; ivec < numProc(); ivec++)
        {
          if (volVec[ivec] > 0.)
            {
              volTot += volVec[ivec];
              if (a_pval == 0)
                {
                  normTot = Max(normTot, normVec[ivec]);
                }
              else
                {
                  normTot += normVec[ivec];
                }
            }
        }
    }

  //broadcast the sum to all processors.
  broadcast(normTot, baseProc);
  broadcast(volTot, baseProc);

  //return the totals
  a_norm   = normTot;
  a_volume = volTot;
}
/***************/
void
EBArith::sumBndryArea(Real&               a_area,
                      const BoxLayout&    a_region,
                      const EBISLayout&   a_ebisl)
{
  a_area = 0.0;
  DataIterator dit = a_region.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      const EBISBox& ebisbox = a_ebisl[dit()];
      const Box&  gridRegion =  a_region.get(dit());
      IntVectSet    ivsBndry = ebisbox.boundaryIVS(gridRegion);
      EBGraph        ebGraph = ebisbox.getEBGraph();
      VoFIterator vofIt(ivsBndry, ebGraph);

      for (vofIt.reset(); vofIt.ok(); ++vofIt)
        {
          VolIndex vof = vofIt();
          int nFace = ebisbox.numFacePhase(vof);
          for (int iface=0; iface<nFace; iface++)
            {
              a_area += ebisbox.bndryArea(vof, iface);
            }
        }
    }
}
/***************/
Real
EBArith::norm(const Vector< LevelData<EBCellFAB>* >& a_src,
              const Vector< DisjointBoxLayout >&     a_grids,
              const Vector< EBISLayout >&            a_ebisl,
              const Vector<int>&                     a_refRatio,
              const int& a_comp,
              const int& a_pval,
              EBNormType::NormMode a_mode)
{
  CH_TIME("EBArith::norm(vector)");
  CH_assert(a_comp >= 0);
  CH_assert(a_pval >= -2);

  Real sum = 0;
  Real volume = 0;
  RealVect dxLev = RealVect::Unit;
  int maxlev = a_src.size() - 1;
  for (int ilev  = 0; ilev <= maxlev; ilev++)
    {
      //don't count stuff covered by finer levels
      IntVectSet ivsExclude;
      if (ilev < maxlev)
        {
          //put next finer grids into an IVS
          for (LayoutIterator lit = a_grids[ilev+1].layoutIterator(); lit.ok(); ++lit)
            {
              ivsExclude |= a_grids[ilev+1].get(lit());
            }
          //coarsen back down to this level.
          //ivs will now hold the image of the next finer level
          ivsExclude.coarsen(a_refRatio[ilev]);
        }

      Real sumLev, volLev;
      volWeightedSum(sumLev, volLev, *(a_src[ilev]),  a_grids[ilev],
                     a_ebisl[ilev], dxLev, ivsExclude,
                     a_comp, a_pval, a_mode);

      sum += sumLev;
      volume += volLev;
      dxLev /= a_refRatio[ilev];
    }

  Real normval;
  if (a_pval == 0)
    {
      normval = sum;
    }
  else if ( (a_pval == 1) || (a_pval == -1) || (a_pval == -2))
    {
      if (volume > 0.0)
        {
          normval = sum/volume;
        }
      else
        {
          normval = 0.0;
        }
    }
  else
    {
      Real denom = a_pval;
      Real exponent = 1.0/denom;
      if (volume > 0.0)
        {
          sum /= volume;
          normval = pow(sum, exponent);
        }
      else
        {
          normval = 0.0;
        }
    }
  return normval;
}
/***************/
Real
EBArith::norm(const BoxLayoutData<EBCellFAB >& a_src,
              const BoxLayout& a_region,
              const EBISLayout& a_ebisl,
              const int& a_comp,
              const int& a_pval,
              EBNormType::NormMode a_mode)
{
  CH_TIME("EBArith::norm");
  CH_assert(a_comp >= 0);
  CH_assert(a_pval >= 0);
  CH_assert(a_comp < a_src.nComp());

  Real sum, volume;
  RealVect dx = RealVect::Unit;
  //nothing excluded (no finer levels as far as i know)
  IntVectSet ivsExclude;
  volWeightedSum(sum, volume, a_src,  a_region,
                 a_ebisl, dx, ivsExclude, a_comp, a_pval, a_mode);
  Real normval;
  if (a_pval == 0)
    {
      normval = sum;
    }
  else if (a_pval == 1)
    {
      if (volume > 0.0)
        {
          normval = Abs(sum)/volume;
        }
      else
        {
          normval = 0.0;
        }
    }
  else
    {
      Real denom = a_pval;
      Real exponent = 1.0/denom;
      if (volume > 0.0)
        {
          sum /= volume;
          normval = pow(sum, exponent);
        }
      else
        {
          normval = 0.0;
        }
    }
  return normval;
}

/***************/
/***************/
Real
EBArith::norm(Real& a_volume, const BoxLayoutData<EBCellFAB >& a_src,
              const BoxLayout& a_region,
              const EBISLayout& a_ebisl,
              const int& a_comp,
              const int& a_pval,
              EBNormType::NormMode a_mode)
{
  CH_assert(a_comp >= 0);
  CH_assert(a_pval >= 0);
  CH_assert(a_comp < a_src.nComp());

  Real sum;
  RealVect dx = RealVect::Unit;
  //nothing excluded (no finer levels as far as i know)
  IntVectSet ivsExclude;
  volWeightedSum(sum, a_volume, a_src,  a_region,
                 a_ebisl, dx, ivsExclude, a_comp, a_pval, a_mode);
  Real normval;
  if (a_pval == 0)
    {
      normval = sum;
    }
  else if (a_pval == 1)
    {
      if (a_volume > 0.0)
        {
          normval = Abs(sum)/a_volume;
        }
      else
        {
          normval = 0.0;
        }
    }
  else
    {
      Real denom = a_pval;
      Real exponent = 1.0/denom;
      if (a_volume > 0.0)
        {
          sum /= a_volume;
          normval = pow(sum, exponent);
        }
      else
        {
          normval = 0.0;
        }
    }
  return normval;
}
/***************/
bool
EBArith::isVoFHere(VolIndex& a_vof2, int& a_whichVoF,
                   const Vector<VolIndex>& a_vofsStencil,
                   const IntVect& a_cell2)
{
  bool found = false;
  for (int isten = 0; isten < a_vofsStencil.size(); isten++)
    {
      if (a_vofsStencil[isten].gridIndex() == a_cell2)
        {
          if (found == true)
            {
              //vof not unique.  return false;
              return false;
            }
          else
            {
              found = true;
              a_vof2 = a_vofsStencil[isten];
              a_whichVoF = isten;
            }
        }
    }

  return found;
}

/***************/
bool
EBArith::isVoFHere(VolIndex& a_vof2,
                   const Vector<VolIndex>& a_vofsStencil,
                   const IntVect& a_cell2)
{
  int whichVoF;
  bool found = isVoFHere(a_vof2, whichVoF, a_vofsStencil, a_cell2);

  return found;
}
/***************/
void
EBArith::monotonePathVoFToCellMultiVoFs(Vector<VolIndex>& a_vofs,
                                        const Vector<VolIndex>& a_vofsInPath,
                                        const IntVect& a_cell2,
                                        const EBISBox& a_ebisBox)
{
  Vector<VolIndex> vofsInCell = a_ebisBox.getVoFs(a_cell2);
  for (int ivof=0; ivof < vofsInCell.size(); ivof++)
    {
      for (int isten = 0; isten < a_vofsInPath.size(); isten++)
        {
          if (a_vofsInPath[isten].gridIndex() == a_cell2 &&
             vofsInCell[ivof].cellIndex() == a_vofsInPath[isten].cellIndex())
//              a_ebisBox.volFrac(vofsInCell[ivof]) == a_ebisBox.volFrac(a_vofsInPath[isten]) &&
//              a_ebisBox.normal(vofsInCell[ivof])  == a_ebisBox.normal(a_vofsInPath[isten]))
            {
              a_vofs.push_back(vofsInCell[ivof]);
            }
        }
    }
}
/***************/
bool
EBArith::monotonePathVoFToCellVoF(VolIndex& a_vof2,
                                  const VolIndex& a_vof1,
                                  const IntVect& a_cell2,
                                  const EBISBox& a_ebisBox)
{
  IntVect diffVect = a_cell2 - a_vof1.gridIndex();
  int imaxdiff = 0;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (Abs(diffVect[idir]) > imaxdiff)
        imaxdiff = Abs(diffVect[idir]);
    }

  IntVect timesMoved = IntVect::Zero;
  IntVect pathSign   = IntVect::Zero;
  Vector<VolIndex> vofsStencil;
  getAllVoFsInMonotonePath(vofsStencil, timesMoved,
                           pathSign, a_vof1, a_ebisBox,
                           imaxdiff);

  bool found = isVoFHere(a_vof2,  vofsStencil, a_cell2);

  return found;
}
/*******/
void
EBArith::
getAllVoFsInMonotonePath(Vector<VolIndex>& a_vofList,
                         const VolIndex&   a_vof,
                         const EBISBox&    a_ebisBox,
                         const int&        a_redistRad)
{
  Vector<VolIndex> vofsStencil;
  IntVect timesMoved = IntVect::Zero;
  IntVect pathSign   = IntVect::Zero;
  EBArith::getAllVoFsInMonotonePath(vofsStencil, timesMoved,
                                    pathSign, a_vof, a_ebisBox,
                                    a_redistRad);
  a_vofList = vofsStencil;
}
/*******/
void
EBArith::
getAllVoFsInMonotonePath(Vector<VolIndex>& a_vofList,
                         const IntVect&    a_timesMoved,
                         const IntVect&    a_pathSign,
                         const VolIndex&   a_vof,
                         const EBISBox&    a_ebisBox,
                         const int&        a_redistRad)
{
  const ProblemDomain& domain = a_ebisBox.getDomain();
  if (domain.contains(a_vof.gridIndex()))
    {
      //check to see if we have already added it.
      //if not, add it
      bool found = false;
      for (int ivof = 0; ivof < a_vofList.size(); ivof++)
        {
          if (a_vofList[ivof] == a_vof)
            {
              found = true;
            }
        }
      if (!found)
        {
          a_vofList.push_back(a_vof);
        }
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          //only move redist radius times in a direction
          if (a_timesMoved[idir] < a_redistRad)
            {
              IntVect newTimesMoved = a_timesMoved;
              newTimesMoved[idir]++;
              //pathSign preserves monotonicity
              //we set pathsign to -1 once you start in the low
              //direction, 1 in the high direction.  It starts at zero

              //if we started in the low direction or have not started
              //add vofs on low side
              if ((a_pathSign[idir] == -1) || (a_pathSign[idir] == 0))
                {
                  IntVect newSign = a_pathSign;
                  newSign[idir] = -1;
                  Vector<FaceIndex> facesLo =
                    a_ebisBox.getFaces(a_vof, idir, Side::Lo);
                  for (int iface = 0; iface < facesLo.size(); iface++)
                    {
                      VolIndex newVoF = facesLo[iface].getVoF(Side::Lo);
                      getAllVoFsInMonotonePath(a_vofList, newTimesMoved, newSign,
                                               newVoF, a_ebisBox, a_redistRad);
                    }
                }
              //if we started in the high direction or have not started
              //add vofs on high side
              if ((a_pathSign[idir] == 1) || (a_pathSign[idir] == 0))
                {
                  IntVect newSign = a_pathSign;
                  newSign[idir] = 1;
                  Vector<FaceIndex> facesHi =
                    a_ebisBox.getFaces(a_vof, idir, Side::Hi);
                  for (int iface = 0; iface < facesHi.size(); iface++)
                    {
                      VolIndex newVoF = facesHi[iface].getVoF(Side::Hi);
                      getAllVoFsInMonotonePath(a_vofList, newTimesMoved, newSign,
                                               newVoF, a_ebisBox, a_redistRad);
                    }
                }
            } //end if (we are less than redist radius away)
        }//end loop over directions
    }
}
#include "NamespaceFooter.H"
