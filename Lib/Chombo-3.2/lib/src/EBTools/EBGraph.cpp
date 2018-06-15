#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  ANAG, LBNL, DTG

#include "EBGraph.H"
#include "DebugOut.H"
#include "EBDebugOut.H"
#include "BoxIterator.H"
#include "VoFIterator.H"
#include "EBArith.H"
#include "FaceIterator.H"
#include "NamespaceHeader.H"

bool    EBGraphImplem::s_verbose = false;
IntVect EBGraphImplem::s_ivDebug = IntVect(D_DECL(11, 5, 3));
Box     EBGraphImplem::s_doDebug = Box(IntVect(D_DECL(0, 0, 0)), IntVect(D_DECL(31, 31, 15)));

/*******************************/
Vector<FaceIndex> EBGraph::getMultiValuedFaces(const int&  a_idir,
                                               const Box&  a_box) const
{
  return
    m_implem->getMultiValuedFaces(a_idir, a_box, *this);
}


/*******************************/
Vector<FaceIndex> EBGraphImplem::getMultiValuedFaces(const int&     a_idir,
                                                     const Box&     a_box,
                                                     const EBGraph& a_ebgraph) const
{
  Vector<FaceIndex> multiValuedFaces;
  Box ghostRegion = a_box;
  ghostRegion.grow(a_idir, 1);
  const IntVectSet ivsMulti = a_ebgraph.getMultiCells(ghostRegion);
  //Use faceiterator to only stop at faces once
  FaceIterator faceit(ivsMulti, a_ebgraph, a_idir, FaceStop::SurroundingWithBoundary);
  for (faceit.reset(); faceit.ok(); ++faceit)
    {
      const IntVect& ivHi = faceit().gridIndex(Side::Hi);
      const IntVect& ivLo = faceit().gridIndex(Side::Lo);
      if (a_box.contains(ivLo) || a_box.contains(ivHi))
        {
          bool isMulti = false;
          if (a_ebgraph.getDomain().contains(ivHi))
            {
              const Vector<FaceIndex> faces = getAllFaces(ivHi, a_idir, Side::Lo);
              isMulti = (faces.size() > 1);
            }
          else if (a_ebgraph.getDomain().contains(ivLo))
            {
              const Vector<FaceIndex> faces = getAllFaces(ivLo, a_idir, Side::Hi);
              isMulti = (faces.size() > 1);
             }
          else
            {
              pout() << "EBGraph::getMultiValuedFaces --  Domain does not contain either ivHi or ivLo" << endl;
              MayDay::Error("EBGraph::getMultiValuedFaces --  domain does not contain either ivHi or ivLo");
            }
          if (isMulti)
            {
              multiValuedFaces.push_back(faceit());
            }
        }
    }
  return multiValuedFaces;
}

/*******************************/
void EBGraph::getRegNextToMultiValued(IntVectSet&    a_vofsToChange,
                                      const EBGraph& a_ghostGraph) const
{
  m_implem->getRegNextToMultiValued(a_vofsToChange, a_ghostGraph);
}

/*******************************/
void EBGraphImplem::getRegNextToMultiValued(IntVectSet&    a_vofsToChange,
                                            const EBGraph& a_ghostGraph) const
{
  Box ghostRegion = a_ghostGraph.getRegion();
  Box region = getRegion();
  CH_assert(ghostRegion.contains(region));

  //loop through multiply-valued vofs in the grown region
  //if any of the vofs next to the multiply-valued vofs
  //are regular, collect them so they can be changed to
  //irregular vofs with unit vol fracs and area fracs.
  a_vofsToChange = IntVectSet();
  IntVectSet multiIVS = a_ghostGraph.getMultiCells(ghostRegion);
  for (VoFIterator vofit(multiIVS, a_ghostGraph);
      vofit.ok(); ++vofit)
    {
      const VolIndex& multiVoF = vofit();
      const IntVect& multiIV = multiVoF.gridIndex();
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          for (SideIterator sit; sit.ok(); ++sit)
            {
              IntVect otherIV = multiIV + sign(sit())*BASISV(idir);
              if ((region.contains(otherIV) && isRegular(otherIV)))
                {
                  a_vofsToChange |= otherIV;
                }
            }
        }
    }
}

/*******************************/
void EBGraph::addFullIrregularVoFs(const IntVectSet& a_vofsToChange,
                                   const EBGraph&    a_ghostGraph)
{
  m_implem->addFullIrregularVoFs(a_vofsToChange, a_ghostGraph);
}

/*******************************/
void EBGraph::addEmptyIrregularVoFs(const IntVectSet& a_vofsToChange)
{
  m_implem->addEmptyIrregularVoFs(a_vofsToChange);
}

/*******************************/
void EBGraphImplem::addEmptyIrregularVoFs(const IntVectSet& a_vofsToChange)
{
  if (!a_vofsToChange.isEmpty())
    {
      CH_assert(isDefined());
      CH_assert(isDomainSet());
      //this is all supposed to be called for covered vofs
      //to be changed to empty irregular vofs
      CH_assert(!isAllRegular());
      if (isAllCovered())
        {
          m_tag = HasIrregular;
          if (m_irregIVS != NULL) delete m_irregIVS;
          if (m_multiIVS != NULL) delete m_multiIVS;
          m_multiIVS = new IntVectSet(DenseIntVectSet(m_region, false));
          m_irregIVS = new IntVectSet(DenseIntVectSet(m_region, false));
          m_graph.define(m_region, 1);
          for (BoxIterator bit(m_region); bit.ok(); ++bit)
            {
              m_graph(bit(), 0).defineAsCovered();
            }
        }

      //  //now for changing vofs
      for (IVSIterator ivsit(a_vofsToChange); ivsit.ok(); ++ivsit)
        {
          const IntVect&  iv = ivsit();
          //needs to be a covered cell to start with
          CH_assert(isCovered(iv));
          VolIndex vof(iv, 0);

          GraphNode coveredNode;
          coveredNode.defineAsCovered();

          //create node and its arcs (no arcs)
          //the coarse-fine info is created later.
          //this operation must be done before all that
          //the finer ones can be set trivially but not the
          //coarse ones.
          GraphNodeImplem node;
          node.m_finerNodes = coveredNode.refine(vof);
          node.m_coarserNode = 0;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              for (SideIterator sit; sit.ok(); ++sit)
                {

                  int gNodeIndex = IrregNode::index(idir, sit());
                  Vector<int>& nodeArcsDir = node.m_arc[gNodeIndex];
                  nodeArcsDir.resize(0);
                }
            }

          //finally add node into graph.
          //again coarse and fine info have to be added on later.
          m_graph(iv, 0).addIrregularNode(node);
          (*m_irregIVS) |= iv;
          if (m_graph(iv, 0).size() > 1)
            {
              MayDay::Error("that vof was already irregular");
            }
        }
    }
}

/*******************************/
void EBGraphImplem::addFullIrregularVoFs(const IntVectSet& a_vofsToChange,
                                         const EBGraph&    a_ghostGraph)
{
  if (!a_vofsToChange.isEmpty())
    {
      CH_assert(isDefined());
      CH_assert(isDomainSet());
      //this is all supposed to be called for regular vofs
      //to be changed to full irregular vofs
      CH_assert(!isAllCovered());
      if (isAllRegular())
        {
          m_tag = HasIrregular;
          if (m_irregIVS != NULL) delete m_irregIVS;
          if (m_multiIVS != NULL) delete m_multiIVS;
          m_multiIVS = new IntVectSet(DenseIntVectSet(m_region, false));
          m_irregIVS = new IntVectSet(DenseIntVectSet(m_region, false));
          m_graph.define(m_region, 1);
          for (BoxIterator bit(m_region); bit.ok(); ++bit)
            {
              m_graph(bit(), 0).defineAsRegular();
            }
        }

      //  //now for changing vofs
      for (IVSIterator ivsit(a_vofsToChange); ivsit.ok(); ++ivsit)
        {
          const IntVect&  iv = ivsit();
          //needs to be a regular cell to start with
          CH_assert(isRegular(iv));
          VolIndex vof(iv, 0);

          GraphNode regularNode;
          regularNode.defineAsRegular();

          //create node and its arcs
          //the coarse-fine info is created later.
          //this operation must be done before all that
          //the finer ones can be set trivially but not the
          //coarse ones.
          GraphNodeImplem node;
          node.m_finerNodes = regularNode.refine(vof);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              for (SideIterator sit; sit.ok(); ++sit)
                {

                  int gNodeIndex = IrregNode::index(idir, sit());
                  Vector<int>& nodeArcsDir = node.m_arc[gNodeIndex];
                  nodeArcsDir.resize(0);
                  //find which vof node is connected to in each direction.
                  //cannot use isConnected here because it will always return
                  //true since one vof is still regular
                  IntVect  otherIV = iv + sign(sit())*BASISV(idir);
                  if (m_domain.contains(otherIV))
                    {
                      Vector<VolIndex> otherVoFs = a_ghostGraph.getVoFs(otherIV);
                      bool found = false;
                      for (int iother = 0; iother < otherVoFs.size(); iother++)
                        {
                          const VolIndex& otherVoF = otherVoFs[iother];
                          Vector<FaceIndex> otherFaces = a_ghostGraph.getFaces(otherVoF, idir, flip(sit()));
                          //there are rare cases where the number of other faces is greater than 1
                          for (int iface = 0; iface < otherFaces.size(); iface++)
                            {
                              nodeArcsDir.push_back(otherVoF.cellIndex());
                              found = true;
                            }
                        }
                      if (!found)
                        {
                          pout() << iv << " was specified as regular  but is not connected to anything in " << otherIV << endl;
                          MayDay::Error("former regular vof not connected to anything");
                        }
                    }
                  else
                    {
                      //boundary arc
                      nodeArcsDir.resize(1, -1);
                    }

                }
            }
          //finally add node into graph.
          //again coarse and fine info have to be added on later.
          m_graph(iv, 0).addIrregularNode(node);
          (*m_irregIVS) |= iv;
          if (m_graph(iv, 0).size() > 1)
            {
              MayDay::Error("that vof was already irregular");
            }
        }
    }
}

/*******************************/
IntVectSet EBGraphImplem::getIrregCells(const Box& a_subbox) const
{
  static IntVectSet emptySet;
  if (m_irregIVS == NULL) return emptySet;
  return (*m_irregIVS) & a_subbox;
}

/*******************************/
IntVectSet EBGraphImplem::getMultiCells(const Box& a_subbox) const
{
  static IntVectSet emptySet;
  if (m_multiIVS == NULL) return emptySet;
  return (*m_multiIVS) & a_subbox;
}

/*******************************/
void EBGraphImplem::setToAllRegular()
{
  CH_assert(isDefined());
  CH_assert(isDomainSet());
  m_tag = AllRegular;
  if (m_irregIVS != NULL) delete m_irregIVS;
  if (m_multiIVS != NULL) delete m_multiIVS;
  m_irregIVS = NULL;
  m_multiIVS = NULL;
}

/*******************************/
void EBGraphImplem::setToAllCovered()
{
  CH_assert(isDefined());
  CH_assert(isDomainSet());
  m_tag = AllCovered;
  if (m_irregIVS != NULL) delete m_irregIVS;
  if (m_multiIVS != NULL) delete m_multiIVS;
  m_irregIVS = NULL;
  m_multiIVS = NULL;
}

/*******************************/
void EBGraphImplem::checkGraph(const BaseFab<int>&      a_regIrregCovered,
                               const Vector<IrregNode>& a_irregGraph,
                               const Box&               a_validRegion,
                               const ProblemDomain&     a_domain)
{
  ///account for regular and covered cells
  for (BoxIterator bit(a_validRegion); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      //-1 is covered, 0 is irregular, and 1 is regular
      int valCen = a_regIrregCovered(iv, 0);
      //the only valid values are 0, 1 and -1
      if (!((valCen == 1) || (valCen == -1) || (valCen == 0)))
        {
          pout() << "EBGraph::checkgraph error 1" << endl;
          MayDay::Error("invalid flag in graph regirregcover3ed");
        }
      //check to see that there is never a regular cell next
      //to a covered cell
      if ((valCen == 1) || (valCen == -1))
        {
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              for (SideIterator sit; sit.ok(); ++sit)
                {
                  IntVect ivshift = BASISV(idir);
                  ivshift[idir] *= sign(sit());
                  IntVect ivneigh = iv + ivshift;
                  if (a_domain.contains(ivneigh))
                    {
                      int valNeigh = a_regIrregCovered(ivneigh, 0);
                      //only way this can happen is if one value is -1 and one is 1
                      //which implies regular next to covered.
                      bool regNextToCovered = (valCen*valNeigh < 0);
                      if (regNextToCovered)
                        {
                          pout() << "EBGraph::checkgraph error 2" << endl;
                          pout() << "graph inconsistent between cells " << iv
                                 << "and " << ivneigh <<  endl;
                          MayDay::Error("regular cell put next to covered cell");
                        }
                    }
                }
            }
        }
    }

  //now check to see that the graph is internally consistent.
  //and that the real data associated with it is reasonable.
  Real eps = 1.0e-6;
  for (int ivec = 0; ivec < a_irregGraph.size(); ivec++)
    {
      const IrregNode& inputNode = a_irregGraph[ivec];
      if ((inputNode.m_volFrac < -eps) || (inputNode.m_volFrac > 1.+eps))
        {
          pout() << "EBGraph::checkgraph error 3 at cell " << inputNode.m_cell << endl;
          MayDay::Error("invalid volume fraction");
        }
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if ((inputNode.m_volCentroid[idir]   < -0.5-eps) || (inputNode.m_volCentroid[idir]   > 0.5+eps))
            {
              pout() << "EBGraph::checkgraph error 4 at cell " << inputNode.m_cell << endl;
              MayDay::Error("invalid volume centroid");
            }
          if ((inputNode.m_bndryCentroid[idir] < -0.5-eps) || (inputNode.m_bndryCentroid[idir] > 0.5+eps))
            {
              pout() << "EBGraph::checkgraph error 5 at cell " << inputNode.m_cell << endl;
              MayDay::Error("invalid boundary centroid");
            }
        }
      //now check arcs and whether stuff associated with that is internally consistent
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          for (SideIterator sit; sit.ok(); ++sit)
            {
              int vecIndex = IrregNode::index(idir, sit());
              const Vector<int>&      arcs  = inputNode.m_arc[vecIndex];
              const Vector<Real>&     areas = inputNode.m_areaFrac[vecIndex];
              const Vector<RealVect>& cents = inputNode.m_faceCentroid[vecIndex];
              if ((arcs.size() != areas.size() ) || (arcs.size() != cents.size()))
                {
                  pout() << "EBGraph::checkgraph error 6 at cell " << inputNode.m_cell
                         << "faces for direction " << idir
                         << "inconsistent with lengths of areas or centroids" << endl;
                  MayDay::Error("inconsistent face information");
                }

              for (int ivec = 0; ivec < arcs.size(); ivec++)
                {
                  const int     &  arc  =  arcs[ivec];
                  const Real    &  area = areas[ivec];
                  const RealVect&  cent = cents[ivec];
                  if ((area < -eps) || (area > 1.+eps))
                    {
                      pout() << "EBGraph::checkgraph error 7 at cell " << inputNode.m_cell << endl;
                      MayDay::Error("invalid area fraction");
                    }
                  for (int jdir = 0; jdir < SpaceDim; jdir++)
                    {
                      if ((jdir != idir) && (((cent[jdir]   < -0.5-eps)) || ((cent[jdir]   > 0.5+eps))))
                        {
                          pout() << "EBGraph::checkgraph error 8 at cell " << inputNode.m_cell << endl;
                          MayDay::Error("invalid face centroid");
                        }
                    }

                  //now see if arc value is consistent with what the basefab thinks is there
                  IntVect ivshift = BASISV(idir);
                  ivshift[idir] *= sign(sit());
                  IntVect ivneigh = inputNode.m_cell + ivshift;
//                   if (inputNode.m_cell[0]==38 && inputNode.m_cell[1]==52 && inputNode.m_cell[2]==40
//                      && ivneigh[0]==38 && ivneigh[1]==52 && ivneigh[2]==41)
//                     {
//                       pout() << "this is the cell" << endl;
//                       pout() << "a_domain " << a_domain << endl;
//                       pout() << "inputNode.m_cell " << inputNode.m_cell << endl;
//                       pout() << "ivneigh " << ivneigh << endl;
//                     }
                  if (a_domain.contains(ivneigh))
                    {
                      int valNeigh = a_regIrregCovered(ivneigh, 0);
                      int otherNodeInd = arc;
//                       if (inputNode.m_cell[0]==38 && inputNode.m_cell[1]==52 && inputNode.m_cell[2]==40
//                          && ivneigh[0]==38 && ivneigh[1]==52 && ivneigh[2]==41)
//                         {
//                           pout() << " valNeigh " << valNeigh << endl;
//                           pout() << " otherNodeInd " << otherNodeInd << endl;
//                         }
                      if (otherNodeInd < -2)
                        {
                          pout() << "EBGraph::checkgraph error 8.5 at cell " << inputNode.m_cell << endl;
                          pout() << "arc has invalid value " << endl;
                          MayDay::Error("arc list whacked ");
                        }

                      if (valNeigh == -1)
                        {
                          pout() << "EBGraph::checkgraph error 9 at cell " << inputNode.m_cell << endl;
                          pout() << "arc reaches into covered cell " << endl;
                          MayDay::Error("arc list inconsistent with basefab ");
                        }

                      if (otherNodeInd == -1)
                        {
                          pout() << "EBGraph::checkgraph error 10 at cell " << inputNode.m_cell << endl;
                          pout() << "arc does not reach  outside domain but is set to -1 " << endl;
                          MayDay::Error("arc list inconsistent with domain ");
                        }
                      else if (otherNodeInd == -2)
                        {
                          if ((a_validRegion.contains(ivneigh)) && (valNeigh != 1))
                            {
                              pout() << "EBGraph::checkgraph error 11 at cell " << inputNode.m_cell << endl;
                              pout() << "arc points to a regular cell but the cell is not regular" << endl;
                              MayDay::Error("arc list inconsistent with basefab ");
                            }
                        }
                      else if (otherNodeInd >= 0)
                        {
//                           if ((a_validRegion.contains(ivneigh)) && (valNeigh != 0) && (area != 1.))//it's ok to have a regular cell over there
                          if ((a_validRegion.contains(ivneigh)) && (valNeigh != 0) && ((area < 1.-eps) || (area > 1.+eps)))//it's ok to have a regular cell over there
//                           if ((a_validRegion.contains(ivneigh)) && (valNeigh != 0) && (area < 1.-eps))//it's ok to have a regular cell over there
                            {
                              pout() << "EBGraph::checkgraph error 12 at cell " << inputNode.m_cell << endl;
                              pout() << "arc points to an irregular cell but the cell is not irregular" << endl;
                              MayDay::Error("arc list inconsistent with basefab ");
                            }
                        }

                    }
                  else
                    {
                      //check to see if arc knows that it is reaching outside the domain
                      if (arc != -1)
                        {
                          pout() << "EBGraph::checkgraph error 13 at cell " << inputNode.m_cell << endl;
                          pout() << "arc reaches outside domain but is not set to -1 " << endl;
                          MayDay::Error("arc list inconsistent with domain ");
                        }
                    }
                }
            }
        }
    }
}

/*******************************/
void EBGraphImplem::buildGraph(const BaseFab<int>&      a_regIrregCovered,
                               const Vector<IrregNode>& a_irregGraph,
                               const Box&               a_validRegion,
                               const ProblemDomain&     a_domain)
{
  define(a_validRegion);
  setDomain(a_domain);

  //remember that this forces a dense representation.
  m_tag = HasIrregular;
  if (m_irregIVS != NULL) delete m_irregIVS;
  if (m_multiIVS != NULL) delete m_multiIVS;
  m_multiIVS = new IntVectSet(DenseIntVectSet(m_region, false));
  m_irregIVS = new IntVectSet(DenseIntVectSet(m_region, false));
  m_graph.define(m_region, 1);

  //set regular and covered cells
  for (BoxIterator bit(m_region); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();
      if (a_regIrregCovered(iv, 0) == 1) //regular cell
        {
          m_graph(iv, 0).defineAsRegular();
        }
      else if (a_regIrregCovered(iv, 0) == -1) //covered cell
        {
          m_graph(iv, 0).defineAsCovered();
        }
      else if (a_regIrregCovered(iv, 0) != 0)
        {
          MayDay::Error("invalid flag");
        }
    }

  //now for irregular cells
  //add the vofs
  for (int ivecIrreg = 0; ivecIrreg < a_irregGraph.size(); ivecIrreg++)
    {
      const IrregNode& inputNode = a_irregGraph[ivecIrreg];
      const IntVect& iv =inputNode.m_cell;

      GraphNodeImplem newImplem;
      newImplem.m_nodeInd = ivecIrreg;

      m_graph(iv, 0).addIrregularNode(newImplem, inputNode.m_cellIndex);
      (*m_irregIVS) |= iv;
      if (m_graph(iv, 0).size() > 1)
        {
          (*m_multiIVS) |= iv;
        }
    }

  //add the faces
  for (int ivecIrreg = 0; ivecIrreg < a_irregGraph.size(); ivecIrreg++)
    {
      const IrregNode& inputNode = a_irregGraph[ivecIrreg];
      const IntVect& iv =inputNode.m_cell;
      Vector<GraphNodeImplem>& vecNodes =
        *(m_graph(iv, 0).m_cellList);
      //pick out which node we are talking about
      //by maching its nodeInd with ivecIrreg
      bool foundNode = false;
      GraphNodeImplem* nodePtr = NULL;
      for (int ivecGraph = 0; ivecGraph < vecNodes.size(); ivecGraph++)
        {
          if (vecNodes[ivecGraph].m_nodeInd == ivecIrreg)
            {
              foundNode = true;
              nodePtr = &(vecNodes[ivecGraph]);
            }
        }
      if (!foundNode)
        {
          MayDay::Error("EBGraph: internal error in construction");
        }
      //now add the arcs in the input to the node
      GraphNodeImplem& node = *nodePtr;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          for (SideIterator sit; sit.ok(); ++sit)
            {
              int irregIndex = IrregNode::index(idir, sit());
              int gNodeIndex = IrregNode::index(idir, sit());
              const Vector<int>& irregArcs = inputNode.m_arc[irregIndex];
              Vector<int>& nodeArcs = node.m_arc[gNodeIndex];
              nodeArcs.resize(irregArcs.size());
              for (int iarc = 0; iarc < irregArcs.size(); iarc++)
                {
                  int otherNodeInd = irregArcs[iarc];
                  if (otherNodeInd == -1)
                    {
                      //if otherNodeInd == -1, boundary arc.
                      //just make the arc in our node = -1
                      //to signify the same
                      nodeArcs[iarc] = -1;
                    }
                  else if (otherNodeInd == -2)
                    {
                      //this means that the vof is connected
                      //to a regular vof,
                      //which always have a cell index of 0
                      nodeArcs[iarc] =  0;
                    }
                  else
                    {
                      nodeArcs[iarc] =  otherNodeInd;
                    }
                }
            }
        }

    }
}

/*******************************/
const Box& EBGraphImplem::getRegion() const
{
  return m_region;
}

/*******************************/
const ProblemDomain& EBGraphImplem::getDomain() const
{
  return m_domain;
}

/*******************************/
EBGraphImplem::EBGraphImplem(const Box& a_box)
  :m_irregIVS(NULL),  m_multiIVS(NULL), m_isMaskBuilt(false)
{
  define(a_box);
}

/*******************************/
void EBGraphImplem::define(const Box& a_region)
{
  CH_assert(!a_region.isEmpty());

  m_tag = AllRegular;
  m_region = a_region;
  if (m_irregIVS != NULL) delete m_irregIVS;
  if (m_multiIVS != NULL) delete m_multiIVS;
  m_irregIVS = NULL;
  m_multiIVS = NULL;
  m_mask.clear();
  m_isMaskBuilt = false;
  m_isDefined= true;
}

/*******************************/
void EBGraphImplem::setDomain(const ProblemDomain& a_domain)
{
  m_isDomainSet = true;
  m_domain = a_domain;
}

/*******************************/
EBGraphImplem::EBGraphImplem()
  :m_irregIVS(NULL),  m_multiIVS(NULL)
{
  m_isDefined = false;
  m_isDomainSet = false;
}

/*******************************/
EBGraphImplem::~EBGraphImplem()
{
  if (m_irregIVS != NULL) delete m_irregIVS;
  if (m_multiIVS != NULL) delete m_multiIVS;
}

/*******************************/
bool EBGraphImplem::isDefined() const
{
  return m_isDefined;
}

/*******************************/
bool EBGraphImplem::isDomainSet() const
{
  return m_isDomainSet;
}

/*******************************/
const BaseFab<int>& EBGraphImplem::getMask(int& a_regIrregCovered) const
{
  if (this->isAllRegular())
    {
      a_regIrregCovered = 1;
    }
  else if (this->isAllCovered())
    {
      a_regIrregCovered = -1;
    }
  else
    {
      a_regIrregCovered = 0;
      if (!m_isMaskBuilt)
        {
          Box maskBox = m_region & m_domain;
          m_mask.define(maskBox, 1);
          fillIntMask(m_mask);
          m_isMaskBuilt = true;
        }
    }
  return m_mask;
}

/*******************************/
Vector<VolIndex> EBGraphImplem::getVoFs(const IntVect& a_iv) const
{
  CH_assert(isDefined());
  CH_assert(isDomainSet());
  Vector<VolIndex> retvec;
  if (m_tag == AllRegular)
    {
      retvec.push_back(VolIndex(a_iv, 0));
    }
  else if (m_tag == AllCovered)
    {
      //return an empty vector
    }
  else if (m_tag == HasIrregular)
    {
      CH_assert(m_region.contains(a_iv));
      CH_assert(m_domain.contains(a_iv));
      const GraphNode& node = m_graph(a_iv, 0);
      retvec = node.getVoFs(a_iv);
    }
  return retvec;
}

/*******************************/
bool EBGraphImplem::isRegular(const IntVect& a_iv) const
{
  CH_assert(isDefined());
  CH_assert(isDomainSet());
  bool retval;
  if (m_tag == AllRegular)
    {
      retval = true;
    }
  else if (m_tag == AllCovered)
    {
      retval = false;
    }
  else if (m_tag == HasIrregular)
    {
      //CH_assert(m_region.contains(a_iv)); //picked up my m_graph already
      //CH_assert(m_domain.contains(a_iv));
      const GraphNode& node = m_graph(a_iv, 0);
      retval = node.isRegular();
    }
  else
    {
      retval = false;
      MayDay::Error("EBGraphImplem::isRegular:Bogus Tag");
    }
  return retval;
}

/*******************************/
Vector<FaceIndex> EBGraphImplem::getAllFaces(const IntVect&        a_iv,
                                             const int&            a_idir,
                                             const Side::LoHiSide& a_sd) const
{
  Vector<FaceIndex> retval(0);
  Vector<VolIndex> vofs = getVoFs(a_iv);
  for (int ivof= 0; ivof < vofs.size(); ivof++)
    {
      retval.append(getFaces(vofs[ivof], a_idir, a_sd));
    }
  return retval;
}

/*******************************/
bool EBGraphImplem::isIrregular(const IntVect& a_iv) const
{
  CH_assert(isDefined());
  CH_assert(isDomainSet());
  bool retval;
  if (m_tag == AllRegular)
    {
      retval = false;
    }
  else if (m_tag == AllCovered)
    {
      retval = false;
    }
  else if (m_tag == HasIrregular)
    {
      CH_assert(m_region.contains(a_iv));
      CH_assert(m_domain.contains(a_iv));
      const GraphNode& node = m_graph(a_iv, 0);
      retval = node.isIrregular();
    }
  else
    {
      retval = false;
      MayDay::Error("EBGraphImplem::isIrregular:Bogus Tag");
    }
  return retval;
}

/*******************************/
bool EBGraphImplem::isAllCovered() const
{
  return (m_tag == AllCovered);
}

/*******************************/
bool EBGraphImplem::isAllRegular() const
{
  return m_tag == AllRegular;
}

/*******************************/
int EBGraphImplem::size(const Box&      a_region,
                        const Interval& a_comps) const
{
  CH_assert(isDefined());
  CH_assert(isDomainSet());
  //regular irregular covered code
  int linearSize = sizeof(int);
  if (!isRegular(a_region) && !isCovered(a_region))
    {
      for (BoxIterator bit(a_region); bit.ok(); ++bit)
        {
          const GraphNode& node = m_graph(bit(), 0);
          int nodeSize = node.linearSize();
          linearSize += nodeSize;
        }
    }
  return linearSize;
}

/*******************************/
void EBGraphImplem::linearOut(void*           a_buf,
                              const Box&      a_region,
                              const Interval& a_comps) const
{
  CH_assert(isDefined());
  CH_assert(isDomainSet());
  int secretCode;
  if (isCovered(a_region))
    {
      secretCode = 0;
    }
  else if (isRegular(a_region))
    {
      secretCode = 1;
    }
  else
    {
      secretCode = 2;
    }
  //  int linearSize = sizeof(int);

  int* intbuf = (int*) a_buf;
  *intbuf = secretCode;
  intbuf++;

  if (!isRegular(a_region) && !isCovered(a_region))
    {
      unsigned char* buffer = (unsigned char*) intbuf;
      for (BoxIterator bit(a_region); bit.ok(); ++bit)
        {
          const GraphNode& node = m_graph(bit(), 0);
          int nodeSize = node.linearSize();
          node.linearOut(buffer);
          buffer += nodeSize;
          //          linearSize += nodeSize;
        }
    }
}

/*******************************/
void EBGraphImplem::linearIn(void*           a_buf,
                             const Box&      a_region,
                             const Interval& a_comps)
{
  CH_assert(isDefined());
  CH_assert(isDomainSet());
  CH_assert(isDefined());
  CH_assert(isDomainSet());
  int* intbuf = (int*) a_buf;
  int secretCode = *intbuf;
  intbuf++;
  //  int linearSize = sizeof(int);

  if (secretCode == 0)
    {
      //all covered input
      EBGraphImplem ebgraphSrc(a_region);
      ebgraphSrc.setDomain(m_domain);
      ebgraphSrc.setToAllCovered();
      copy(a_region, a_comps, a_region, ebgraphSrc, a_comps);
    }
  else if (secretCode == 1)
    {
      //all regular input
      EBGraphImplem ebgraphSrc(a_region);
      ebgraphSrc.setDomain(m_domain);
      ebgraphSrc.setToAllRegular();
      copy(a_region, a_comps, a_region, ebgraphSrc, a_comps);
    }
  else
    {
      CH_assert(secretCode==2);
      unsigned char* buffer = (unsigned char*) intbuf;

      if (isAllRegular() || isAllCovered())
        {
          m_tag = HasIrregular;
          if (m_irregIVS != NULL) delete m_irregIVS;
          if (m_multiIVS != NULL) delete m_multiIVS;
          m_multiIVS = new IntVectSet(DenseIntVectSet(m_region, false));
          m_irregIVS = new IntVectSet(DenseIntVectSet(m_region, false));
          m_graph.define(m_region, 1);
        }
      for (BoxIterator bit(a_region); bit.ok(); ++bit)
        {
          GraphNode& node = m_graph(bit(), 0);
          node.linearIn(buffer);
          if (node.isIrregular()>0) (*m_irregIVS)|=bit();
          if (node.size()>1) (*m_multiIVS)|=bit();
          int nodeSize = node.linearSize();
          buffer += nodeSize;
          //          linearSize += nodeSize;
        }
    }
}

/*******************************/
bool EBGraph::hasIrregular() const
{
  return m_implem->hasIrregular();
}

/*******************************/
bool EBGraphImplem::hasIrregular() const
{
  return m_tag == HasIrregular;
}

/*******************************/
bool EBGraphImplem::isCovered(const IntVect& a_iv) const
{
  CH_assert(isDefined());
  CH_assert(isDomainSet());
  bool retval;
  if (m_tag == AllRegular)
    {
      retval = false;
    }
  else if (m_tag == AllCovered)
    {
      retval = true;
    }
  else if (m_tag == HasIrregular)
    {
      //CH_assert(m_region.contains(a_iv)); this check picked up by m_graph
      //CH_assert(m_domain.contains(a_iv));
      const GraphNode& node = m_graph(a_iv, 0);
      retval = node.isCovered();
    }
  else
    {
      retval = false;
      MayDay::Error("EBGraphImplem::isIrregular:Bogus Tag");
    }
  return retval;
}

/*******************************/
bool EBGraphImplem::isCovered(const Box& a_box) const
{
  CH_assert(isDefined());
  CH_assert(isDomainSet());
  bool retval;
  if (m_tag == AllRegular)
    {
      retval = false;
    }
  else if (m_tag == AllCovered)
    {
      retval = true;
    }
  else if (m_tag == HasIrregular)
    {
      CH_assert(m_region.contains(a_box));
      CH_assert(m_domain.contains(a_box));
      retval = true;
      BoxIterator bit(a_box);
      for (bit.reset(); bit.ok() && retval; ++bit)
        {
          if (!isCovered(bit()))
            {
              retval = false;
            }
        }
    }
  else
    {
      retval = false;
      MayDay::Error("EBGraphImplem::isCovered:Bogus Tag");
    }
  return retval;
}

/*******************************/
bool EBGraphImplem::isRegular(const Box& a_box) const
{
  CH_assert(isDefined());
  CH_assert(isDomainSet());
  bool retval;
  if (m_tag == AllRegular)
    {
      retval = true;
    }
  else if (m_tag == AllCovered)
    {
      retval = false;
    }
  else if (m_tag == HasIrregular)
    {
      CH_assert(m_region.contains(a_box));
      CH_assert(m_domain.contains(a_box));
      BoxIterator bit(a_box);
      retval = true;
      for (bit.reset(); bit.ok() && retval; ++bit)
        {
          if (!isRegular(bit()))
            {
              retval = false;
            }
        }
    }
  else
    {
      retval = false;
      MayDay::Error("EBGraphImplem::isRegular:Bogus Tag");
    }
  return retval;
}

/*******************************/
Vector<FaceIndex> EBGraphImplem::getFaces(const VolIndex&       a_vof,
                                          const int&            a_idir,
                                          const Side::LoHiSide& a_sd) const
{
  CH_assert(isDefined());
  CH_assert(isDomainSet());
  CH_assert(m_region.contains(a_vof.gridIndex()));
  CH_assert(m_domain.contains(a_vof.gridIndex()));
  CH_assert((a_idir >= 0) && (a_idir < SpaceDim));
  CH_assert((a_sd == Side::Lo) || (a_sd == Side::Hi));

  Vector<FaceIndex> retvec;
  if (m_tag == AllRegular)
    {
      IntVect otherIV = a_vof.gridIndex()
        + sign(a_sd)*BASISV(a_idir);
      int otherCellInd = 0;
      if (!m_domain.contains(otherIV))
        {
          otherCellInd = -1;
        }
      VolIndex otherVoF(otherIV, otherCellInd);
      retvec.push_back(FaceIndex(a_vof, otherVoF, a_idir));
    }
  else if (m_tag == AllCovered)
    {
      //return empty vector
    }
  else if (m_tag == HasIrregular)
    {
      const IntVect& iv = a_vof.gridIndex();
      const GraphNode& node = m_graph(iv, 0);
      retvec = node.getFaces(a_vof, a_idir, a_sd, m_domain);
    }

  return retvec;
}

void EBGraph::fillMask(BaseFab<char>& a_mask) const
{
  m_implem->fillMask(a_mask);
}

void EBGraph::fillCellTypeMask(BaseFab<char>& a_mask) const
{
  m_implem->fillCellTypeMask(a_mask);
}

void EBGraph::fillIntMask(BaseFab<int>& a_mask) const
{
  m_implem->fillIntMask(a_mask);
}

void EBGraphImplem::fillMask(BaseFab<char>& a_mask) const
{
  Box b = a_mask.box() & m_graph.box();
  if (b.isEmpty()) return;
  BoxIterator bit(b);
  for (; bit.ok(); ++bit)
    {
      const GraphNode& g = m_graph(bit(), 0);
      CH_assert(g.size() < 128);
      a_mask(bit(), 0) = (char)(g.size());
    }
}

void EBGraphImplem::fillIntMask(BaseFab<int>& a_mask) const
{
  if (m_tag == AllRegular )
    {
      a_mask.setVal(1);
    }
  else if (m_tag == AllCovered )
    {
      a_mask.setVal(-1);
    }
  else
    {
      Box b = a_mask.box();
      for (BoxIterator bit(b); bit.ok(); ++bit)
        {
          if (isRegular(bit()))
            {
              a_mask(bit(), 0) =  1;
            }
          else if (isCovered(bit()))
            {
              a_mask(bit(), 0) = -1;
            }
          else if (numVoFs(bit()) > 1)
            {
              //multivalued need to be treated as covered as far as the mask is concerned
              a_mask(bit(), 0) = -1;
            }
          else
            {
              a_mask(bit(), 0) = 0;
            }
        }
    }
}

/*******************************/
void EBGraphImplem::fillCellTypeMask(BaseFab<char>& a_mask) const
{
  // This is added for VisIt support.
  //  0       outside priblem domain
  //  1       covered
  //  2       regular
  //  3       irregular
  if (m_tag == AllCovered )
    {
      a_mask.setVal(1);
    }
  else
  if (m_tag == AllRegular )
    {
      a_mask.setVal(2);
    }
  else
    {
      Box b = a_mask.box() & m_graph.box();
      if (b.isEmpty()) return;
      for (BoxIterator bit(b); bit.ok(); ++bit)
        {
          if (isCovered(bit()))
            {
              a_mask(bit(), 0) = 1;
            }
          else
          if (isRegular(bit()))
            {
              a_mask(bit(), 0) = 2;
            }
          else
          if (isIrregular(bit()))
            {
              a_mask(bit(), 0) = 3;
            }
          else
            {
              a_mask(bit(), 0) = 0;
            }
        }
    }
}

Vector<FaceIndex> EBGraph:: getIrregFaces(const Box& a_box,
                                          int        a_dir) const
{
  return m_implem->getIrregFaces(a_box, a_dir);
}

Vector<FaceIndex>  EBGraphImplem::getIrregFaces(const Box& a_box,
                                                int        a_dir) const
{
  Vector<FaceIndex> faces;
  if (m_tag == AllRegular || m_tag == AllCovered)
    {
      // do nothing
    }
  else
    {
      Box b = a_box & m_domain;
      IntVectSet ivs = this->getIrregCells(b);
      for (IVSIterator it(ivs); it.ok(); ++it)
        {
          const GraphNode& node = m_graph(it(), 0);
          IntVect hi = it() + BASISV(a_dir);
          CH_assert(node.isIrregular());
          if (b.contains(hi) && m_graph(hi, 0).isIrregular())
          {
            faces.append(node.getFaces(it(), a_dir, Side::Hi, m_domain));
          }
        }
    }
  return faces;
}

/*******************************/
FaceIndex EBGraphImplem::coarsen(const FaceIndex& a_fineFace) const
{
  VolIndex loVoFCoar = coarsen(a_fineFace.getVoF(Side::Lo));
  VolIndex hiVoFCoar = coarsen(a_fineFace.getVoF(Side::Hi));
  return FaceIndex(loVoFCoar, hiVoFCoar, a_fineFace.direction());
}

/*******************************/
Vector<FaceIndex> EBGraphImplem::refine(const FaceIndex&     a_coarFace,
                                        const EBGraphImplem& a_fineGraph) const
{
#ifndef NDEBUG
  // This code was originally:
  //     CH_assert( CH_XD::refine(m_domain, 2) == a_fineGraph.m_domain );
  // This failed to compile with gcc 4.1.x because this is a "const" method
  ProblemDomain fineDomain = m_domain;
  fineDomain.refine(2);
  CH_assert( fineDomain == a_fineGraph.m_domain );
#endif
  Vector<FaceIndex> retval;
  IntVect ivLoCoar = a_coarFace.gridIndex(Side::Lo);
  IntVect ivHiCoar = a_coarFace.gridIndex(Side::Hi);
  int direction = a_coarFace.direction();
  if (m_region.contains(ivLoCoar) && m_region.contains(ivHiCoar))
    {
      //interior face
      Vector<VolIndex> loVoFFine = refine(a_coarFace.getVoF(Side::Lo));
      Vector<VolIndex> hiVoFFine = refine(a_coarFace.getVoF(Side::Hi));
      int idir = a_coarFace.direction();
      for (int ilo = 0; ilo < loVoFFine.size(); ilo++)
        {
          for (int ihi = 0; ihi < hiVoFFine.size(); ihi++)
            {
              if (a_fineGraph.isConnected(loVoFFine[ilo], hiVoFFine[ihi]))
                {
                  FaceIndex newFace(loVoFFine[ilo], hiVoFFine[ihi], idir);
                  retval.push_back(newFace);
                }
            }
        }
    }
  else if (m_region.contains(ivLoCoar))
    {
      //boundary face on the hi side of the domain
      Vector<VolIndex> loVoFsFine = refine(a_coarFace.getVoF(Side::Lo));
      Box fineRegion = ebrefine(m_region, 2);
      for (int ivof = 0; ivof < loVoFsFine.size(); ivof++)
        {
          VolIndex& loVoFFine = loVoFsFine[ivof];
          IntVect ivHiFine = loVoFFine.gridIndex() + BASISV(direction);
          if (!fineRegion.contains(ivHiFine))
            {
              Vector<FaceIndex> fineFaces = a_fineGraph.getFaces(loVoFFine, direction, Side::Hi);
              retval.append(fineFaces);
            }
        }
    }
  else if (m_region.contains(ivHiCoar))
    {
      //boundary face on the low side of the domain
      Vector<VolIndex> hiVoFsFine = refine(a_coarFace.getVoF(Side::Hi));
      Box fineRegion = ebrefine(m_region, 2);
      for (int ivof = 0; ivof < hiVoFsFine.size(); ivof++)
        {
          VolIndex& hiVoFFine = hiVoFsFine[ivof];
          IntVect ivLoFine = hiVoFFine.gridIndex() - BASISV(direction);
          if (!fineRegion.contains(ivLoFine) )
            {
              Vector<FaceIndex> fineFaces = a_fineGraph.getFaces(hiVoFFine, direction, Side::Lo);
              retval.append(fineFaces);
            }
        }
    }
  else
    {
      MayDay::Error("neither vof of the face is inside the domain");
    }
  return retval;
}

/*******************************/
bool EBGraphImplem::isConnected(const VolIndex& a_vof1,
                                const VolIndex& a_vof2) const
{
  CH_TIME("EBGraphImplem::isConnected");

  CH_assert(isDefined());
  const IntVect& iv1 = a_vof1.gridIndex();
  const IntVect& iv2 = a_vof2.gridIndex();

  CH_assert(m_region.contains(iv1));
  CH_assert(m_region.contains(iv2));
  VolIndex vofLo, vofHi;
  //use the lexical greaterThan because it will
  //always give the correct answer since only one component
  //is different.
  if (iv1.lexGT(iv2))
    {
      vofLo = a_vof2;
      vofHi = a_vof1;
    }
  else
    {
      vofLo = a_vof1;
      vofHi = a_vof2;
    }
  int direction;
  bool dirfound;
  const IntVect& ivLo = vofLo.gridIndex();
  const IntVect& ivHi = vofHi.gridIndex();

  dirfound = false;
  for (int idir = 0; ((idir<SpaceDim) && !dirfound); idir++)
    {
      if ((ivHi - ivLo) == BASISV(idir))
        {
          direction = idir;
          dirfound = true;
        }
    }
  //if not neigboring intvects, no way it can be connected
  if (!dirfound) return false;

  Vector<FaceIndex> faces = getFaces(vofLo, direction, Side::Hi);
  bool voffound = false;
  for (int iface = 0; iface < faces.size(); iface++)
    {
      const FaceIndex& face = faces[iface];
      if (face.getVoF(Side::Hi) == vofHi)
        {
          voffound = true;
        }
    }
  return voffound;
}

// // a) assumes neither 1 or 2 is covered
// // b) assumes both are within this graphImplem
// // c) assumes are within 1 cell of each other.
// bool
// EBGraphImplem::isNeighbourConnected(const VolIndex& a_vof1,
//                                     const VolIndex& a_vof2) const
// {
//   const IntVectSet& ivs = *m_irregIVS;
//   if (!ivs.contains(a_vof1.gridIndex()) || !ivs.contains(a_vof2.gridIndex()))
//     return true ;// returns true.  regular cell connected to all neighbours.

//   IntVect delta = a_vof1.gridIndex() - a_vof2.gridIndex();
//   int sum = D_TERM(abs(delta[0]), + abs(delta[1]), +abs(delta[2]));
//   if (sum > 1) return false;

//   D_TERM(int dir=0;,
//          if (delta[0]==0) dir=1;,
//          if (delta[1]==0) dir=2;)

//   int index1 = a_vof1.cellIndex();
//   int index2 = a_vof2.cellIndex();
//   const GraphNodeImplem& node = m_graph(a_vof1.gridIndex(), 0)->m_cellList[a_vof1.cellIndex()];
//   const Vector<int>& arcs = node.m_arc[0];
//   if (delta[dir] > 0)
//     arcs = node.m_arc[dir];
//   else
//     arcs = node.m_arc[CH_SPACEDIM + dir];

//   for (int i=0; i<arcs.size(); ++i)
//     {
//       if (arcs[i] ==  a_vof2.cellIndex()) return true;
//     }
//   return false;
// }

/*******************************/
VolIndex EBGraphImplem::coarsen(const VolIndex& a_fineVoF) const
{
  VolIndex retval;
  IntVect ivfine = a_fineVoF.gridIndex();
  IntVect ivcoar = ebcoarsen(ivfine, 2);
  if (!m_domain.contains(ivfine))
    {
      //boundary vof
      retval = VolIndex(ivcoar, -1);
    }
  else if ((m_tag == AllRegular ) || (m_tag == AllCovered ))
    {
      retval = VolIndex(ivcoar, 0);
    }
  else
    {
      CH_assert(m_tag == HasIrregular);
      const IntVect& iv = a_fineVoF.gridIndex();
      const GraphNode& node = m_graph(iv, 0);
      retval = node.coarsen(a_fineVoF);
    }
  return retval;
}

/*******************************/
Vector<VolIndex> EBGraphImplem::refine(const VolIndex& a_coarVoF) const
{
  CH_assert(isDefined());
  CH_assert(isDomainSet());

  Vector<VolIndex> retval(0);
  IntVect ivCoar = a_coarVoF.gridIndex();
  CH_assert(m_domain.contains(ivCoar));

  if (m_tag == AllRegular)
    {
      if (s_verbose)
        {
          pout() << "tag is all regular" << endl;
        }
      const IntVect& iv = a_coarVoF.gridIndex();
      Box refbox = ebrefine(Box(iv,iv),2);
      BoxIterator bit(refbox);
      for (bit.reset(); bit.ok(); ++bit)
        {
          retval.push_back(VolIndex(bit(), 0));
        }
    }
  else if (m_tag == AllCovered)
    {
      //return empty vector
      if (s_verbose)
        {
          pout() << "tag is all covered" << endl;
        }
    }
  else
    {
      CH_assert(m_tag == HasIrregular);
      const IntVect& iv = a_coarVoF.gridIndex();
      const GraphNode& node = m_graph(iv, 0);
      if (s_verbose)
        {
          pout() << "tag is has irregular" << endl;
          if (node.isIrregular())
            {
              pout() << " node has " << node.m_cellList->size() << " elements" << endl;
              const Vector<GraphNodeImplem>& nodes = *node.m_cellList;
              for (int inode = 0; inode < nodes.size(); inode++)
                {
                  const GraphNodeImplem node = nodes[inode];
                  pout() << "finernodes size  = " <<  node.m_finerNodes.size() << endl;
                }
            }
        }
      retval = node.refine(a_coarVoF);
    }
  return retval;
}

/*******************************/
void EBGraphImplem::copy(const Box&           a_regionFrom,
                         const Interval&      a_intDst,
                         const Box&           a_regionTo,
                         const EBGraphImplem& a_source,
                         const Interval&      a_intSrc)
{
  CH_TIME("EBGraphImplem::copy");
  CH_assert(isDefined());
  CH_assert(isDomainSet());
  if (isRegular(a_regionTo) && a_source.isRegular(a_regionFrom))
    {
      return;
    }
  else if (isCovered(a_regionTo) && a_source.isCovered(a_regionFrom))
    {
      return;
    }
  else if (a_source.isCovered(a_regionFrom) && a_regionTo.contains(m_region))
    {
      setToAllCovered();
      return;
    }
  else if (a_source.isRegular(a_regionFrom) && a_regionTo.contains(m_region))
    {
      setToAllRegular();
      return;
    }
  else if (isAllRegular() && a_source.isAllCovered())
    {
      //define the basefab as all regular and set the region to
      //covered in the intersection
      m_tag = HasIrregular;
      if (m_irregIVS != NULL) delete m_irregIVS;
      if (m_multiIVS != NULL) delete m_multiIVS;
      m_multiIVS = new IntVectSet(DenseIntVectSet(m_region, false));
      m_irregIVS = new IntVectSet(DenseIntVectSet(m_region, false));
      m_graph.define(m_region, 1);
      GraphNode regularNode;
      regularNode.defineAsRegular();
      m_graph.setVal(regularNode);
      Box interBox = m_region & a_regionTo;
      for (BoxIterator bit(interBox); bit.ok(); ++bit)
        {
          m_graph(bit(), 0).defineAsCovered();
        }
    }
  else if (isAllCovered() && a_source.isAllRegular())
    {
      //define the basefab as all covered and set the region to
      //regular in the intersection
      m_tag = HasIrregular;
      if (m_irregIVS != NULL) delete m_irregIVS;
      if (m_multiIVS != NULL) delete m_multiIVS;
      m_multiIVS = new IntVectSet(DenseIntVectSet(m_region, false));
      m_irregIVS = new IntVectSet(DenseIntVectSet(m_region, false));
      m_graph.define(m_region, 1);
      GraphNode  coveredNode;
      coveredNode.defineAsCovered();
      m_graph.setVal(coveredNode);
      Box interBox = m_region & a_regionTo;
      for (BoxIterator bit(interBox); bit.ok(); ++bit)
        {
          m_graph(bit(), 0).defineAsRegular();
        }
    }
  else
    {
      //one or both has irregular cells.
      //because i am sick of combinatorics,
      //use basefab copy to transfer the data.

      //see if we need to generate a source fab
      BaseFab<GraphNode>* srcFabPtr;
      bool needToDelete;
      if (a_source.hasIrregular())
        {
          srcFabPtr = (BaseFab<GraphNode>*)&a_source.m_graph;
          needToDelete = false;
        }
      else
        {
          needToDelete = true;
          srcFabPtr = new BaseFab<GraphNode>(a_regionFrom, 1);
          GraphNode srcVal;
          if (a_source.isAllRegular())
            {
              srcVal.defineAsRegular();
            }
          else
            {
              //this really has to be true
              CH_assert(a_source.isAllCovered());
              srcVal.defineAsCovered();
            }
          srcFabPtr->setVal(srcVal);
        }
      //define our graph if i need to. leave alone otherwise
      if (isAllRegular() || isAllCovered())
        {
          if (m_irregIVS != NULL) delete m_irregIVS;
          if (m_multiIVS != NULL) delete m_multiIVS;
          m_multiIVS = new IntVectSet(DenseIntVectSet(m_region, false));
          m_irregIVS = new IntVectSet(DenseIntVectSet(m_region, false));
          m_graph.define(m_region, 1);
        }

      //copy the data
      m_tag = HasIrregular;
      m_graph.copy(a_regionFrom, a_intSrc,
                   a_regionTo, *srcFabPtr, a_intDst);

      //if we needed to new the basefab, clean it up
      if (needToDelete)
        {
          delete srcFabPtr;
        }
    }
  //  now fix up the IntVectSets to match the information
  if (a_source.hasIrregular())
    {
      IntVectSet ivsInterIrreg = (*a_source.m_irregIVS) & a_regionTo;
      if (!ivsInterIrreg.isEmpty())
        {
          if (m_irregIVS == NULL)
          {
            m_irregIVS = new IntVectSet(DenseIntVectSet(m_region, false));
            m_multiIVS = new IntVectSet(DenseIntVectSet(m_region, false));
          }
          for (IVSIterator it(ivsInterIrreg); it.ok(); ++it)
            {
              IntVect iv = it();
              (*m_irregIVS) |= iv;
              if (numVoFs(iv) > 1) // this will be correct since we already
                // did a m_graph copy operation
                {
                  (*m_multiIVS) |= iv;
                }
            }
        }
    }
}

/*******************************/
void EBGraphImplem::coarsenVoFs(const EBGraphImplem& a_fineGraph,
                                const Box&           a_coarRegion)
{
  CH_TIME("EBGraphImplem::coarsenVoFs");

  //this also defines the boxes
  m_region = a_coarRegion;
  m_domain = ebcoarsen(a_fineGraph.getDomain(), 2);
  m_isDomainSet = true;
  m_isDefined = true;
  if (a_fineGraph.isCovered(ebrefine(m_region, 2)))
    {
      m_tag = AllCovered;
    }
  else if (a_fineGraph.isRegular(ebrefine(m_region, 2)))
    {
      m_tag = AllRegular;
    }
  else
    {
      Box fineRegion = a_coarRegion;
      fineRegion.refine(2);
      int numFineVoFs = a_fineGraph.numVoFs(fineRegion);
      int numCoarVoFs = 0;
      m_tag = HasIrregular;
      m_region = a_coarRegion;
      m_domain = ebcoarsen(a_fineGraph.m_domain, 2);
      if (m_irregIVS != NULL) delete m_irregIVS;
      if (m_multiIVS != NULL) delete m_multiIVS;
      m_multiIVS = new IntVectSet(DenseIntVectSet(m_region, false));
      m_irregIVS = new IntVectSet(DenseIntVectSet(m_region, false));
      m_graph.define(m_region, 1);
      for (BoxIterator bit(m_region); bit.ok(); ++bit)
        {
          CH_TIME("EBGraphImplem::coarsenVoFs_BoxIterator");
          bool blab = ((bit() == s_ivDebug) && (m_domain.domainBox() == s_doDebug));
          blab = false;
          Box fineBox = ebrefine(Box(bit(), bit()), 2);
          if (a_fineGraph.isRegular(fineBox))
            {
              m_graph(bit(), 0).defineAsRegular();
              if (blab)
                pout() << s_ivDebug << " is regular" << endl;
            }
          else if (a_fineGraph.isCovered(fineBox))
            {
              m_graph(bit(), 0).defineAsCovered();
              if (blab)
                pout() << s_ivDebug << " is covered" << endl;
            }
          else
            {
              //get sets of all connected vofs in the box
              //the number of sets will be the number of
              //vofs in the coarse cell
              if (blab)
                pout() << s_ivDebug << " is covered" << endl;

              Vector<Vector<VolIndex> > fineVoFSets
                = a_fineGraph.getVoFSets(fineBox);

              for (int iset = 0; iset < fineVoFSets.size(); iset++)
                {
                  numCoarVoFs++;
                  GraphNodeImplem newImplem;
                  newImplem.m_finerNodes = fineVoFSets[iset];
                  m_graph(bit(), 0).addIrregularNode(newImplem);
                  (*m_irregIVS) |= bit();
                  if (m_graph(bit(), 0).size() > 1)
                    {
                      (*m_multiIVS) |= bit();
                    }
                }
            }
        }
      if (numCoarVoFs > numFineVoFs) MayDay::Error("Coarsening generated more VoFs");
    }
}

/*******************************/
Vector<Vector<VolIndex> >  EBGraphImplem::getVoFSets(const Box& a_region) const
{
  Vector<Vector<VolIndex> > retval;
  //gather all  vofs
  Vector<VolIndex> allVoFs;
  for (BoxIterator bit(a_region); bit.ok(); ++bit)
    {
      allVoFs.append(getVoFs(bit()));
    }

  std::vector<bool> beenAdded(allVoFs.size(), false);
  for (int ivof = 0; ivof < allVoFs.size(); ivof++)
    {
      if (!beenAdded[ivof])
        {
          const VolIndex& thisVoF = allVoFs[ivof];
          //we have a vof to start with.
          //now find all vofs connected to this inside
          //this cell----not necessarily directly---
          //hence the inner loop of kvof.
          Vector<VolIndex> thisVoFSet(1, thisVoF);
          beenAdded[ivof] = true;
          bool doneAdding = false;
          while (!doneAdding)
            {
              int oldSetSize = thisVoFSet.size();
              //i can start at ivof+1 here because 0 to ivof
              //has always been added.
              for (int jvof = ivof+1; jvof < allVoFs.size(); jvof++)
                {
                  const VolIndex& testVoF = allVoFs[jvof];
                  //this might be a tad confusing to people because
                  //the length of the vector can be changing inside
                  //the loop.  fortran would puke horribly here.
                  for (int kvof = 0; kvof < thisVoFSet.size(); kvof++)
                    {
                      //because of the nature of push_back,
                      //the vector's pointer can change inside so beware.
                      //need to hold the vof by value here
                      VolIndex innerVoF = thisVoFSet[kvof];
                      if (!beenAdded[jvof] && isConnected(innerVoF, testVoF))
                        {
                          thisVoFSet.push_back(testVoF);
                          beenAdded[jvof] = true;
                        }
                    }
                }
              doneAdding = (thisVoFSet.size() == oldSetSize);
            } //end while !doneadding
          //add the finished bunch to the list of  vofs.
          retval.push_back(thisVoFSet);
        } //end this bunch of connected vofs
    }

  //can't hurt to check in case some wacky case was missed.
  for (int ivof = 0; ivof < allVoFs.size(); ivof++)
    {
      if (!beenAdded[ivof])
        {
          MayDay::Error("We seem to have missed a vof in getVoFSets");
        }
    }

  return retval;
}

/*******************************/
Vector<int> EBGraphImplem::coarsenFaces(const VolIndex&       a_coarVoF,
                                        const EBGraphImplem&  a_coarGhostGraph,
                                        const EBGraphImplem&  a_fineGraph,
                                        const int&            a_idir,
                                        const Side::LoHiSide& a_sd)
{
  CH_TIME("EBGraphImplem::coarsenFaces_1");

  Vector<int> retval;

  IntVect coarIV = a_coarVoF.gridIndex();
  IntVect otherIV= coarIV + sign(a_sd)*BASISV(a_idir);
  Vector<VolIndex> theseFineVoFs = a_coarGhostGraph.refine(a_coarVoF);
  if (m_domain.contains(otherIV))
    {
      CH_TIME("EBGraphImplem::coarsenFaces_1_contains_otherIV");

      //interior faces.
      //just get all possible vofs to connect to and
      //check connectivity
      Vector<VolIndex> otherCoarVoFs = a_coarGhostGraph.getVoFs(otherIV);
      for (int iotherCoar = 0; iotherCoar < otherCoarVoFs.size(); iotherCoar++)
        {
          const VolIndex& otherCoarVoF = otherCoarVoFs[iotherCoar];
          Vector<VolIndex> otherFineVoFs = a_coarGhostGraph.refine(otherCoarVoF);
          bool addThisFace = false;
          for (int iotherFine = 0; iotherFine < otherFineVoFs.size(); iotherFine++)
            {
              const VolIndex& otherFineVoF = otherFineVoFs[iotherFine];
              for (int ithisFine = 0; ithisFine < theseFineVoFs.size(); ithisFine++)
                {
                  const VolIndex& thisFineVoF = theseFineVoFs[ithisFine];
                  int thesetwocon = -1;
                  if (a_fineGraph.isConnected(thisFineVoF, otherFineVoF))
                    {
                      addThisFace = true;
                      thesetwocon = 1;
                    }
                }
            }
          if (addThisFace)
            {
              retval.push_back(iotherCoar);
            }
        }
    }
  else
    {
      CH_TIME("EBGraphImplem::coarsenFaces_1_no_otherIV");

      //boundary faces.
      //if there are any boundary faces on the fine level,
      //make one here too
      bool hasBoundaryFaceFine = false;
      for (int ithis = 0; ithis < theseFineVoFs.size() && !hasBoundaryFaceFine ; ithis++)
        {
          const VolIndex& thisVoF = theseFineVoFs[ithis];
          Vector<FaceIndex> fineFaces =
            a_fineGraph.getFaces(thisVoF,a_idir, a_sd);
          for (int iface = 0; iface < fineFaces.size() && !hasBoundaryFaceFine; iface++)
            {
              if (fineFaces[iface].isBoundary())
                {
                  hasBoundaryFaceFine = true;
                }
            }
        }
      if (hasBoundaryFaceFine)
        {
          //remember that -1 is the secret code for boundary face arcs
          retval.push_back(-1);
        }
    }

  return retval;
}

/*******************************/
void EBGraphImplem::coarsenFaces(const EBGraphImplem& a_coarGhostGraph,
                                 const EBGraphImplem& a_fineGraph)
{
  CH_TIME("EBGraphImplem::coarsenFaces_2");

  if (hasIrregular())
    {
      CH_assert(a_coarGhostGraph.getDomain() == m_domain);
      for (BoxIterator bit(m_region); bit.ok(); ++bit)
        {
          if (isIrregular(bit()))
            {
              Vector<VolIndex> vofsCoar = getVoFs(bit());
              Vector<GraphNodeImplem>& nodes =
                *(m_graph(bit(), 0).m_cellList);
              for (int ivof = 0; ivof < vofsCoar.size(); ivof++)
                {
                  const VolIndex& vofCoar= vofsCoar[ivof];
                  GraphNodeImplem& node = nodes[vofCoar.cellIndex()];
                  for (int idir = 0; idir < SpaceDim; idir++)
                    {
                      for (SideIterator sit; sit.ok(); ++sit)
                        {
                          Vector<int> coarArcs =
                            coarsenFaces(vofsCoar[ivof],
                                         a_coarGhostGraph,
                                         a_fineGraph,
                                         idir, sit());

                          int nodeind = IrregNode::index(idir, sit());
                          node.m_arc[nodeind] = coarArcs;
                        }
                    }
                }
            }
        }
    }
}

/*******************************/
void EBGraphImplem::fixFineToCoarse(EBGraphImplem& a_fineGraph) const
{
  if (hasIrregular())
    {
      for (BoxIterator bit(m_region); bit.ok(); ++bit)
        {
          if (isIrregular(bit()))
            {
              const IntVect& ivCoar = bit();

              const Vector<GraphNodeImplem>&
                vofsCoar = *(m_graph(ivCoar, 0).m_cellList);

              int numVofsCoar = vofsCoar.size();

              for (int icoar = 0; icoar < numVofsCoar; icoar++)
                {
                  VolIndex vofCoar(ivCoar, icoar);
                  Vector<VolIndex> vofsFine = refine(vofCoar);

                  for (int ifine = 0; ifine < vofsFine.size(); ifine++)
                    {
                      const IntVect& ivFine = vofsFine[ifine].gridIndex();

                      if (a_fineGraph.isIrregular(ivFine))
                        {
                          int cellIndexFine = vofsFine[ifine].cellIndex();

                          Vector<GraphNodeImplem>&
                            nodesFine = *(a_fineGraph.m_graph(ivFine, 0).m_cellList);

                          nodesFine[cellIndexFine].m_coarserNode = icoar;
                        }
                      else if ((numVofsCoar > 1) && (a_fineGraph.isRegular(ivFine)))
                        {
                          GraphNode& nodeFine = a_fineGraph.m_graph(ivFine, 0);

                          nodeFine.m_cellList = new(nodeFine.alloc()) Vector<GraphNodeImplem>(1);

                          (*(nodeFine.m_cellList))[0].m_isRegular   = true;
                          (*(nodeFine.m_cellList))[0].m_coarserNode = icoar;
                        }
                    }
                }
            }
        }
    }
}

/*******************************/
long long EBGraphImplem::numVoFs(const IntVect& a_iv) const
{
//   Vector<VolIndex> vofs = getVoFs(a_iv);
//   return vofs.size();
  CH_assert(isDefined());
  CH_assert(isDomainSet());
  long long count = 0;
  if (m_tag == AllRegular)
    {
      return 1;
    }
  else if (m_tag == AllCovered)
    {
      return 0;
    }
  else
    {
      CH_assert(m_region.contains(a_iv));
      CH_assert(m_domain.contains(a_iv));
      count = m_graph(a_iv, 0).size();
    }
  return count;
}

/*******************************/
EBGraph::EBGraph(const Box& a_box, int a_comps)
  : m_implem( RefCountedPtr<EBGraphImplem>( new EBGraphImplem(a_box) ) )
{
}

/*******************************/
EBGraph::EBGraph(const Box& a_box)
  :  m_implem( RefCountedPtr<EBGraphImplem>( new EBGraphImplem(a_box) ) )
{
}

/*******************************/
EBGraph::EBGraph(const EBGraph& a_ebiin)
  :  m_implem( a_ebiin.m_implem )
{
}

/*******************************/
EBGraph::EBGraph()
  : m_implem( RefCountedPtr<EBGraphImplem>( new EBGraphImplem() ) )
{
}

/*******************************/
EBGraph::~EBGraph()
{
}

/*******************************/
long long EBGraph::numVoFs(const Box& a_box) const
{
  return m_implem->numVoFs(a_box);
}

/*******************************/
long long EBGraphImplem::numVoFs(const Box& a_box) const
{
  CH_assert(isDefined());
  CH_assert(isDomainSet());
  long long count = 0;
  if (m_tag == AllRegular)
    {
      return a_box.numPts();
    }
  else if (m_tag == AllCovered)
    {
      return 0;
    }
  else
    {
      Box locRegion = a_box & m_region;
      count = 0;
      for (BoxIterator bit(locRegion); bit.ok(); ++bit)
        {
          count += numVoFs(bit());
        }
    }
  return count;
}

/*******************************/
/* This is going to be expensive.  Let's see if pointer comparison is enough.
bool
EBGraphImplem::operator<(const EBGraphImplem & that) const
{
// Data members:
  Box m_region;
  ProblemDomain m_domain;
  TAG m_tag;
  BaseFab<GraphNode> m_graph;
  bool m_isDefined;
  bool m_isDomainSet;
  IntVectSet* m_irregIVS;
  IntVectSet* m_multiIVS;
}
*/

/*******************************/
void EBGraph::define(const Box& a_box)
{
  m_implem->define(a_box);
}

/*******************************/
int EBGraph::size(const Box& R, const Interval& comps) const
{
  return  (m_implem->size(R, comps));
}

/*******************************/
void EBGraph::linearOut(void* buf, const Box& R, const Interval& comps) const
{
  m_implem->linearOut(buf, R, comps);
}

/*******************************/
void EBGraph::linearIn(void* buf, const Box& R, const Interval& comps)
{
  m_implem->linearIn(buf, R, comps);
}

/*******************************/
template <> void dataTypes(Vector<hid_t>& a_types,
                           const EBGraph& a_dummySpecializationArg)
{
#ifdef CH_USE_HDF5
  a_types.resize(1); a_types[0] = H5T_NATIVE_INT;
#endif
}

template <> void dataSize(const EBGraph&  a_item,
                          Vector<int>&    a_sizes,
                          const Box&      a_box,
                          const Interval& a_comps)
{
  a_sizes[0] = a_item.size(a_box, a_comps) / sizeof(int);
}

template <> void write(const EBGraph&  a_item,
                       Vector<void*>&  a_allocatedBuffers,
                       const Box&      a_box,
                       const Interval& a_comps)
{
  a_item.linearOut(a_allocatedBuffers[0], a_box, a_comps);
}

template <> void read(EBGraph&        a_item,
                      Vector<void*>&  a_allocatedBuffers,
                      const Box&      a_box,
                      const Interval& a_comps)
{
  a_item.linearIn(a_allocatedBuffers[0], a_box, a_comps);
}

template <> const char* name(const EBGraph& a_dummySpecializationArg)
{
  static const char* name = "EBGraph";
  return name;
}

#include "NamespaceFooter.H"
