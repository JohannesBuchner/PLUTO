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

#include "GraphNode.H"
#include "BoxIterator.H"
#include "parstream.H"
#include "EBArith.H"
#include "NamespaceHeader.H"

BArena GraphNode::s_Arena("GraphNode Vectors");

/*******************************/
GraphNodeImplem::GraphNodeImplem(const GraphNodeImplem& a_impin)
{

  for (int iarc = 0; iarc < 2*SpaceDim; iarc++)
    {
      m_arc[iarc] = a_impin.m_arc[iarc];
    }

  m_isRegular   = a_impin.m_isRegular;
  m_isValid     = a_impin.m_isValid;
  m_coarserNode = a_impin.m_coarserNode;
  m_nodeInd     = a_impin.m_nodeInd;
  m_finerNodes  = a_impin.m_finerNodes;
}

/*******************************/
GraphNodeImplem& GraphNodeImplem::operator=(const GraphNodeImplem& a_impin)
{
  if (&a_impin != this)
    {
      for (int iarc = 0; iarc < 2*SpaceDim; iarc++)
        {
          m_arc[iarc] = a_impin.m_arc[iarc];
        }

      m_isRegular   = a_impin.m_isRegular;
      m_isValid     = a_impin.m_isValid;
      m_coarserNode = a_impin.m_coarserNode;
      m_nodeInd     = a_impin.m_nodeInd;
      m_finerNodes  = a_impin.m_finerNodes;
    }

  return *this;
}

/*******************************/
int GraphNode::size() const
{
  int retval;
  if (isRegular())
    {
      retval =  1;
    }
  else if (isCovered())
    {
      retval =  0;
    }
  else
    {
      retval =  m_cellList->size();
    }
  return retval;
}

void* GraphNode::alloc()
{
  void* rtn = s_Arena.alloc(sizeof(Vector<GraphNodeImplem>));
#ifdef CH_USE_MEMORY_TRACKING
  s_Arena.bytes += sizeof(Vector<GraphNodeImplem>);
  CH_assert(s_Arena.bytes > 0);
  if (s_Arena.bytes > s_Arena.peak)
    {
      s_Arena.peak = s_Arena.bytes;
    }
#endif
  return rtn;
}

/*******************************/
void GraphNode::addIrregularNode(const GraphNodeImplem& a_nodein, int cellIndex)
{
  if (!hasValidCellList())
    {
      m_cellList = new(alloc()) Vector<GraphNodeImplem>();
    }

  if (m_cellList->size() < cellIndex+1)
    {
      m_cellList->resize(cellIndex+1);
    }

  (*m_cellList)[cellIndex] = a_nodein;
}

void GraphNode::addIrregularNode(const GraphNodeImplem& a_nodein)
{
  if (!hasValidCellList())
    {
      m_cellList = new(alloc()) Vector<GraphNodeImplem>();
    }

  m_cellList->push_back(a_nodein);
}

/*******************************/
Vector<VolIndex> GraphNode::getVoFs(const IntVect& a_iv) const
{
  Vector<VolIndex> retvec;
  if (isCovered())
    {
      //return empty vector
    }
  else if (isRegular())
    {
      retvec.push_back(VolIndex(a_iv, 0));
    }
  else
    {
      const Vector<GraphNodeImplem>& vofVec = *m_cellList;
      for (int ivec = 0; ivec < vofVec.size(); ivec++)
        {
          VolIndex vof(a_iv, ivec);
          retvec.push_back(vof);
        }
    }
  return retvec;
}

Vector<FaceIndex> GraphNode::getFaces(const IntVect&        a_this,
                                      const int&            a_idir,
                                      const Side::LoHiSide& a_sd,
                                      const ProblemDomain&  a_domain) const
{
  static Vector<FaceIndex> emptyVec;
  static Vector<FaceIndex> regularVec(1);

  IntVect otherIV = a_this +sign(a_sd)*BASISV(a_idir);
  if (isRegular())
    {
      VolIndex vof(a_this, 0);
      FaceIndex& face = regularVec[0];
      // if node is regular, the other iv must be single valued
      int otherCellIndex = 0;
      if (!a_domain.contains(otherIV))
        {
          otherCellIndex = -1;
        }
      VolIndex otherVoF(otherIV, otherCellIndex);
      face.define(vof, otherVoF, a_idir);
      return regularVec;
    }
  else if (isCovered())
    {
      //return empty vector
      return emptyVec;
    }
  else
    {
      const Vector<GraphNodeImplem>& nodeVec = *m_cellList;
      if (nodeVec.size()==1)
        {
          const GraphNodeImplem& node =  nodeVec[0];
          const Vector<int>& arcs = node.m_arc[node.index(a_idir, a_sd)];
          if (arcs.size()==0)
            {
              return emptyVec;
            }
          VolIndex vof(a_this, 0);
          if (arcs.size()==1)
            {
              FaceIndex& face = regularVec[0];
              VolIndex otherVoF(otherIV, arcs[0]);
              face.define(vof, otherVoF, a_idir);
              return regularVec;
            }
          Vector<FaceIndex> faces;
          for (int a=0; a<arcs.size(); a++)
            {
              VolIndex otherVoF(otherIV, arcs[a]);
              faces.push_back(FaceIndex(vof, otherVoF, a_idir));
            }
          return faces;
        }
      Vector<FaceIndex> faces;
      for (int v = 0; v<nodeVec.size(); ++v)
        {
          const GraphNodeImplem& node =  nodeVec[v];
          const Vector<int>& arcs = node.m_arc[node.index(a_idir, a_sd)];
          VolIndex vof(a_this, v);
          for (int a=0; a<arcs.size(); a++)
            {
              VolIndex otherVoF(otherIV, arcs[a]);
              faces.push_back(FaceIndex(vof, otherVoF, a_idir));
            }
        }
      return faces;
    }
  return emptyVec;  //make compiler happy
}

/*******************************/
Vector<FaceIndex> GraphNode::getFaces(const VolIndex&       a_vof,
                                      const int&            a_idir,
                                      const Side::LoHiSide& a_sd,
                                      const ProblemDomain&  a_domain) const
{
  static Vector<FaceIndex> emptyVec;
  static Vector<FaceIndex> regularVec(1);

  IntVect otherIV = a_vof.gridIndex() +sign(a_sd)*BASISV(a_idir);

  if (isRegular())
    {
      FaceIndex& face = regularVec[0];
      // if node is regular, the other iv must be single valued
      int otherCellIndex = 0;
      if (!a_domain.contains(otherIV))
        {
          otherCellIndex = -1;
        }
      VolIndex otherVoF(otherIV, otherCellIndex);
      face.define(a_vof, otherVoF, a_idir);
      return regularVec;
    }
  else if (isCovered())
    {
      //return empty vector
      return emptyVec;
    }
  else
    {
      const Vector<GraphNodeImplem>& nodeVec = *m_cellList;
      const GraphNodeImplem& node =  nodeVec[a_vof.cellIndex()];
      const Vector<int>& arcs = node.m_arc[node.index(a_idir, a_sd)];
      if (arcs.size()==0)
        {
          return emptyVec;
        }
      if (arcs.size()==1)
        {
          FaceIndex& face = regularVec[0];
          VolIndex otherVoF(otherIV, arcs[0]);
          face.define(a_vof, otherVoF, a_idir);
          return regularVec;
        }

      Vector<FaceIndex> retvec;
      //cell index of the list is the same as the
      //index into the vector. if the input cell
      //index is too big (or < 0), we can tell by Vector
      //going out of bounds

      for (int ivec = 0; ivec < arcs.size(); ivec++)
        {
          VolIndex otherVoF(otherIV, arcs[ivec]);
          retvec.push_back(FaceIndex(a_vof, otherVoF, a_idir));
        }
      return retvec;
    }
  return emptyVec;
}

/*******************************/
Vector<VolIndex> GraphNode::refine(const VolIndex& a_coarVoF) const
{
  Vector<VolIndex> retvec;
  if (isCovered())
    {
      //return empty vector
    }
  else if (isRegular())
    {
      const IntVect& iv = a_coarVoF.gridIndex();
      Box refbox = ebrefine(Box(iv,iv),2);
      BoxIterator bit(refbox);
      for (bit.reset(); bit.ok(); ++bit)
        {
          retvec.push_back(VolIndex(bit(), 0));
        }
    }
  else
    {
      //irregular node.
      //cell index of the list is the same as the
      //index into the vector. if the input cell
      //index is too big (or < 0), we can tell by Vector
      //going out of bounds
      const Vector<GraphNodeImplem>& nodeVec = *m_cellList;
      const GraphNodeImplem& node =  nodeVec[a_coarVoF.cellIndex()];
      retvec = node.m_finerNodes;
    }
  return retvec;
}

/*******************************/
const GraphNode& GraphNode::operator=(const GraphNode& a_nodein)
{
  if (this != &a_nodein)
    {
      clear();
      setDefaults();
      //if node is regular or covered, just copy the pointer
      //otherwise, append the list of nodes
      if ((a_nodein.isRegularWithSingleValuedParent()) || (a_nodein.isCovered()))
        {
          m_cellList = a_nodein.m_cellList;
        }
      else
        {
          m_cellList = new(alloc()) Vector<GraphNodeImplem>();
          (*m_cellList) = (*a_nodein.m_cellList);
        }
    }
  return *this;
}

/*******************************/
GraphNode::GraphNode(const GraphNode& a_nodein)
{
  //if node is regular or covered, just copy the pointer
  //otherwise, append the list of nodes
  if ((a_nodein.isRegularWithSingleValuedParent()) || (a_nodein.isCovered()))
    {
      m_cellList = a_nodein.m_cellList;
    }
  else
    {
      m_cellList = new(alloc()) Vector<GraphNodeImplem>();
      (*m_cellList) = *(a_nodein.m_cellList);
    }
}

/*******************************/
VolIndex GraphNode::coarsen(const VolIndex& a_fineVoF) const
{
  IntVect ivCoar = ebcoarsen(a_fineVoF.gridIndex(), 2);
  int cellIndexCoar = 0;

  if (isRegularWithSingleValuedParent() || isCovered())
    {
      // Already set correctly
    }
  else if (isRegularWithMultiValuedParent())
    {
      cellIndexCoar = (*m_cellList)[0].m_coarserNode;
    }
  else
    {
      const Vector<GraphNodeImplem>& nodes = *m_cellList;
      int inode = a_fineVoF.cellIndex();

      cellIndexCoar = nodes[inode].m_coarserNode;
    }

  return VolIndex(ivCoar, cellIndexCoar);
}

/*******************************/
int GraphNode::linearSize() const
{
  int retval;

  if (isRegularWithSingleValuedParent() || isCovered())
    {
      // regular/irregular covered
      retval = sizeof(int);
    }
  else
    {
      // regular/irregular covered then
      // number of vofs
      retval = 2*sizeof(int);
      // node data
      const Vector<GraphNodeImplem>& nodes = *m_cellList;
      for (int inode = 0; inode < nodes.size(); inode++)
        {
          retval +=  nodes[inode].linearSize();
        }
    }

  return retval;
}

/*******************************/
void GraphNode::linearOut(void*  a_buf) const
{
  int secretCode;

  if (isRegularWithSingleValuedParent())
    {
      //secret code for regular with single-valued parent.
      //trying to not have two separate secret codes in one class,
      //this matches the m_cellList val for regular
      secretCode = 1;
    }
  else if (isCovered())
    {
      //secret code for covered
      //trying to not have two separate secret codes in one class,
      //this matches the m_cellList val for covered
      secretCode =  0;
    }
  else
    {
      CH_assert(hasValidCellList());
      //secret code for irregular or regular with multi-valued parent.
      secretCode =  2;
    }

  int* intbuf = (int *) a_buf;

  //regular/irregular covered
  *intbuf = secretCode;
  intbuf++;

  if (hasValidCellList())
    {
      int nvofs = m_cellList->size();
      //number of vofs
      *intbuf =  nvofs;
      intbuf++;

      //using intbuf for the two ints we just extracted
      unsigned char* buffer=(unsigned char*) intbuf;

      //now put in the actual nodes
      const Vector<GraphNodeImplem>& nodes = *m_cellList;
      for (int inode = 0; inode < nodes.size(); inode++)
        {
          nodes[inode].linearOut(buffer);

          int nodeSize = nodes[inode].linearSize();
          buffer += nodeSize;
        }
    }
}

/*******************************/
void GraphNode::linearIn(void* a_buf)
{
  int* intbuf = (int *) a_buf;

  int secretCode = *intbuf;
  intbuf++;
  if (secretCode == 1)
    {
      //secret code for regular with single-valued parent.
      //trying to not have two separate secret codes in one class,
      //this matches the m_cellList val for regular
      defineAsRegular();
    }
  else if (secretCode == 0)
    {
      //secret code for covered
      //trying to not have two separate secret codes in one class,
      //this matches the m_cellList val for covered
      defineAsCovered();
    }
  else
    {
      //secret code for irregular or regular with multi-valued parent.
      CH_assert(secretCode == 2);

      //regular/irregular covered
      //number of vofs
      int nvofs = *intbuf;
      intbuf++;

      //using intbuf for the two ints we just extracted
      unsigned char* buffer = (unsigned char*) intbuf;

      //now pull out the actual nodes
      for (int inode = 0; inode < nvofs; inode++)
        {
          GraphNodeImplem newNode;
          newNode.linearIn(buffer);

          addIrregularNode(newNode);

          int nodeSize = newNode.linearSize();
          buffer += nodeSize;
        }
    }
}

/*******************************/
int GraphNodeImplem::linearSize() const
{
  int linSize = 0;

  //isRegular flag
  linSize += sizeof(int);

  //isValid flag
  linSize += sizeof(int);

  //arc sizes
  for (int iarc = 0; iarc < 2*SpaceDim; iarc++)
    {
      int thisArcSize = m_arc[iarc].size();
      //space for each int in each vector +
      //the size of the vector
      linSize += sizeof(int)*(thisArcSize + 1);
    }

  //coarser node
  linSize += sizeof(int);

  //finer nodess size
  linSize += sizeof(int);

  //finer nodes
  for (int inode = 0; inode < m_finerNodes.size(); inode++)
    {
      linSize += m_finerNodes[inode].linearSize();
    }
  return linSize;
}

/*******************************/
void GraphNodeImplem::linearOut(void*  a_buf) const
{
  int* intbuf = (int*) a_buf;
  int linSize = 0;

  // isRegular flag
  *intbuf = m_isRegular;
  linSize += sizeof(int);
  intbuf++;

  // isValid flag
  *intbuf = m_isValid;
  linSize += sizeof(int);
  intbuf++;

  //space for each int in each vector + size of the vector
  for (int iarc = 0; iarc < 2*SpaceDim; iarc++)
    {
      const Vector<int>& thisArc = m_arc[iarc];
      *intbuf = thisArc.size();
      intbuf++;
      linSize += sizeof(int);
      for (int ivec = 0; ivec < thisArc.size(); ivec++)
        {
          *intbuf = thisArc[ivec];

          intbuf++;
          linSize += sizeof(int);
        }
    }

  //coarser node
  *intbuf = m_coarserNode;
  linSize += sizeof(int);
  intbuf++;

  //finer nodess size
  *intbuf = m_finerNodes.size();
  linSize += sizeof(int);
  intbuf++;

  //finer nodes
  char* charbuf = (char*) intbuf;
  for (int inode = 0; inode < m_finerNodes.size(); inode++)
    {
      m_finerNodes[inode].linearOut(charbuf);
      int thisSize = m_finerNodes[inode].linearSize();
      linSize += thisSize;
      charbuf += thisSize;
    }
}

/*******************************/
void GraphNodeImplem::linearIn(void* a_buf)
{
  int* intbuf = (int*) a_buf;
  int linSize = 0;

  // isRegular flag
  m_isRegular = *intbuf;
  linSize += sizeof(int);
  intbuf++;

  // isValid flag
  m_isValid = *intbuf;
  linSize += sizeof(int);
  intbuf++;

  //i am only outputting the arcs.
  //space for each int in each vector + size of the vector
  for (int iarc = 0; iarc < 2*SpaceDim; iarc++)
    {
      Vector<int>& thisArc = m_arc[iarc];
      int thisArcSize = *intbuf;
      intbuf++;
      linSize += sizeof(int);
      thisArc.resize(thisArcSize);
      for (int ivec = 0; ivec < thisArc.size(); ivec++)
        {

          thisArc[ivec] = *intbuf;
          intbuf++;
          linSize += sizeof(int);
        }
    }

  //coarser node
  m_coarserNode = *intbuf;
  linSize += sizeof(int);
  intbuf++;

  //finer nodess size
  int fineNodeSize = *intbuf;
  linSize += sizeof(int);
  intbuf++;

  //finer nodes
  m_finerNodes.resize(fineNodeSize);
  char* charbuf = (char*) intbuf;
  for (int inode = 0; inode < m_finerNodes.size(); inode++)
    {
      m_finerNodes[inode].linearIn(charbuf);
      int thisSize = m_finerNodes[inode].linearSize();
      linSize += thisSize;
      charbuf += thisSize;
    }
}

#include "NamespaceFooter.H"
