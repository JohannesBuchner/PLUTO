#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "Divergence.H"
#include "DivergenceF_F.H"
#include "CellToEdge.H"
#include "NamespaceHeader.H"

// ---------------------------------------------------------
void levelDivergenceMAC(LevelData<FArrayBox>& a_div,
                        const LevelData<FluxBox>& a_uEdge,
                        const Real a_dx)
{

  // silly way to do this until i figure out a better
  // way to make this dimensionally-independent
  CH_assert (a_uEdge.nComp() >= a_div.nComp());


  DataIterator dit = a_div.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {
    a_div[dit()].setVal(0.0);

    const FluxBox& thisFluxBox = a_uEdge[dit()];
    Box cellBox(thisFluxBox.box());
    // just to be sure we don't accidentally trash memory...
    cellBox &= a_div[dit()].box();

    // now loop over coordinate directions and add to divergence
    for (int dir=0; dir<SpaceDim; dir++)
    {
      const FArrayBox& uEdgeDir = thisFluxBox[dir];

      FORT_DIVERGENCE(CHF_CONST_FRA(uEdgeDir),
                      CHF_FRA(a_div[dit()]),
                      CHF_BOX(cellBox),
                      CHF_CONST_REAL(a_dx),
                      CHF_INT(dir));
    }

  }

}

void simpleDivergenceMAC( FArrayBox& a_div, const FluxBox& a_uEdge,
                          const Real a_dx)
{

  a_div.setVal(0.0);
  const Box& cellBox = a_uEdge.box();

  // now loop over coordinate directions and increment divergence
  for (int dir=0; dir<SpaceDim; dir++)
  {
    const FArrayBox& uEdgeDir = a_uEdge[dir];

    FORT_DIVERGENCE(CHF_CONST_FRA(uEdgeDir),
                    CHF_FRA(a_div),
                    CHF_BOX(cellBox),
                    CHF_CONST_REAL(a_dx),
                    CHF_INT(dir));
  }
}



#include "NamespaceFooter.H"
