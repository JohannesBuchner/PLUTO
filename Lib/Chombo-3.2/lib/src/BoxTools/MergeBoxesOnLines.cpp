#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// Merge boxes on lines class
// this is used for line solves with AMR
// ---------
//  class MergeBoxesOnLines
//
///////////////////////////////////////////////////////////////////////////////

// Include files:

#include "Box.H"
#include "MayDay.H"
#include "parstream.H"
#include "MergeBoxesOnLines.H"
#include "NamespaceHeader.H"
void MergeBoxesOnLines::
mergeBoxes(Vector<Vector<Box> >& a_vvbox,
           const int&   a_dir)
{
  for (int ilev = 0; ilev < a_vvbox.size(); ilev++)
    {
      mergeBoxes(a_vvbox[ilev],a_dir);
    }
}
void MergeBoxesOnLines::
mergeBoxes(Vector<Box>& a_vboxin,
           const int&   a_dir)
{

  //loop through all the boxes and chop them based on their neighbors
  bool keepChopping = (a_vboxin.size()>1);
  while (keepChopping)
    {
      keepChopping = false;
      for (SideIterator sit; sit.ok(); ++sit)
        {
          for (int dir = 0; dir < SpaceDim; dir++)
            {
              if (dir != a_dir)
                {
                  Vector<Box> vboxtemp;
                  for (int ibox = 0; ibox < a_vboxin.size(); ibox++)
                    {
                      Box& boxThis = a_vboxin[ibox];
                      Box growThis = boxThis;
                      growThis.grow(1);
                      growThis.grow(a_dir,-1);
                      for (int jbox = 0; jbox < a_vboxin.size(); jbox++)
                        {
                          if (ibox!=jbox)
                            {
                              const Box& boxNeigh = a_vboxin[jbox];

                              //check if this is a chop candidate (only chop based on a_dir neighbor boxes)
                              if (!boxNeigh.intersects(growThis))
                                {
                                  //get an adjacent box to the neighbor in the tangent direction dir
                                  Box sideNeigh = adjCellBox(boxNeigh,dir,sit(),1);
                                  //grow this in the line direction (a_dir)
                                  sideNeigh.grow(a_dir,1);
                                  if (sideNeigh.intersects(boxThis))
                                    {
                                      int chop_pnt = sideNeigh.bigEnd(dir);
                                      if (chop_pnt>boxThis.smallEnd(dir))
                                        {
                                          //chop this box based on the neighbor (results in two boxes)
                                          //pout() << "Chopping ibox = " << ibox << std::endl;
                                          Box newBox = boxThis.chop(dir,chop_pnt);
                                          if (sit()==Side::Lo)
                                            {
                                              boxThis.growHi(dir,1);
                                              newBox.growLo(dir,-1);
                                            }
                                          //add the new box to vboxtemp
                                          if (!newBox.isEmpty())
                                            {
                                              keepChopping = true;
                                              vboxtemp.push_back(newBox);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                      //add this box to the new list
                      vboxtemp.push_back(boxThis);
                    }
                  a_vboxin = vboxtemp;
                }
            }
        }
    }

  //loop through all the boxes and merge them based on the neighbors
  bool keepMerging = (a_vboxin.size()>1);
  while (keepMerging)
    {
      keepMerging = false;
      Vector<Box> vboxtemp;
      for (int ibox = 0; ibox < a_vboxin.size(); ibox++)
        {
          Box& boxThis = a_vboxin[ibox];
          Box newBoxThis = boxThis;
          for (int jbox = 0; jbox < a_vboxin.size(); jbox++)
            {
              if (ibox<jbox)
                {
                  Box& boxNeigh = a_vboxin[jbox];
                  if (!boxNeigh.isEmpty())
                    {
                      //get an adjacent box to the neighbor in the tangent direction dir
                      Box growNeigh = boxNeigh;
                      //grow this in the line direction (a_dir)
                      growNeigh.grow(a_dir,1);
                      if (growNeigh.intersects(boxThis))
                        {
                          //pout() << "Merging ibox = " << ibox << std::endl;
                          keepMerging = true;
                          //merge this box with the neighbor, keeping the neighbor
                          newBoxThis = minBox(newBoxThis,boxNeigh);
                          boxNeigh = Box();//makes the neighbor empty
                        }
                    }
                }
            }
          //add this box to the new list only if it's not empty
          if (!newBoxThis.isEmpty())
            {
              vboxtemp.push_back(newBoxThis);
            }
        }
      a_vboxin = vboxtemp;
    }

  //lexigraphically sort the box vector to keep some locality
  a_vboxin.sort();

}

MergeBoxesOnLines::MergeBoxesOnLines()
{
}
MergeBoxesOnLines::~MergeBoxesOnLines()
{
}

#include "NamespaceFooter.H"
