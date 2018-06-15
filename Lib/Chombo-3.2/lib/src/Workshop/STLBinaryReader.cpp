#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "STLBinaryReader.H"
#include <stdint.h>
#include "STLUtil.H"
#include "KDTree.H"
#include "CH_Timer.H"

#include "NamespaceHeader.H"

using namespace STLUtil;

/*
 * Reads binary STL files and generates a mesh
 * see http://www.cplusplus.com/doc/tutorial/files/
 */

/// Constructor - read from standard input
STLBinaryReader::STLBinaryReader()
{
  m_header = NULL;

  ReadData(cin,0);
}

/// Constructor - read from file name
STLBinaryReader::STLBinaryReader(const string& a_filename)
{
  m_header = NULL;

  ifstream curFile;
  curFile.open(a_filename.c_str(),ios::in|ios::binary);
  if (!curFile.good() || !curFile.is_open())
  {
    MayDay::Abort("STLBinaryReader - unable to open file");
  }

  ReadData(curFile,0);

  curFile.close();
}

/// Destructor
STLBinaryReader::~STLBinaryReader()
{
  delete m_header;
}

/// Return pointer to header string
string* STLBinaryReader::GetHeader() const
{
  return m_header;
}


/// Return number of elements
void STLBinaryReader::GetNtri(int& a_ntri) const
{
  a_ntri = m_ntri;
}

/// Return whether number of elements from header matches file
void STLBinaryReader::GetNtriMatch(bool& a_ntriMatch) const
{
  a_ntriMatch = m_ntriMatch;
}

/// Return pointer to the mesh
RefCountedPtr<STLMesh> STLBinaryReader::GetMesh() const
{
  return m_stlmesh;
}

void STLBinaryReader::ReadData(istream&   a_file,
                               const int offset)
{

  CH_TIMERS("STLBinaryReader::ReadData");
  CH_TIMER("read in mesh",tmesh);
  CH_TIMER("set new data",tdat);
  CH_TIMER("check existing vertices",tvert);
  CH_TIMER("check existing edges",tedge);

  CH_START(tmesh);

  RefCountedPtr<STLMesh> temp(new STLMesh());
  m_stlmesh = temp;
  //m_header = new string(80);
  m_stlmesh->tol = 1.0e-10; // maybe do something fancier later, for now just constant

  char* memblock;
  a_file.seekg(0,ios::end); // go to end of file
  ifstream::pos_type fsize = a_file.tellg(); // get total size of file
  a_file.seekg(0,ios::beg); // go to beginning of file
  size_t isize = fsize;
  memblock = (char *)malloc(isize);
  a_file.read(memblock, fsize); // read it into memory

  int chunksize;
  int readsize=0;

  // read header
  //char hstr[80] = " ";
  //memcpy( hstr , memblock , chunksize );
  chunksize = sizeof(uint8_t)*80;
  m_header = new string( memblock , chunksize );
  //m_header->assign( memblock , chunksize ); // copy the header string
  readsize += chunksize; // increment progress in file

  uint32_t ntri;
  chunksize = sizeof(uint32_t);
  memcpy( &ntri , memblock+readsize , chunksize ); // copy first 4 bytes
  m_ntri = (int) ntri;
  readsize += chunksize; // increment progress in file

  CH_STOP(tmesh);

  // create KDtree for nodes and edgemap for map to speed up searching for connectivity info
  int KDError;
  KDTree * nodetree = (KDTree *) KDCreate( SpaceDim , &KDError );
  if (KDError!=0)
    pout() << "KDCreate returned an error" << endl;
  Vector<int>* data = new Vector<int>(3*m_ntri);
  KDSetGlobalData( nodetree , (void *) data );
  Real rvarray[SpaceDim];
  Real closestarray[SpaceDim];
  int treeInsertIndex = 0;

  // create map for mesh edges (note, only 1st 2 indices of intvect matter here
  if (SpaceDim<2)
    pout() << "Error, STLBinaryReader does not support SpaceDim<2" << endl;
#ifndef STL_UNORDERED_MAP
  typedef map<IntVect, int, IVCompareSWO> MeshEdgeMap;
#else
  typedef unordered_map<IntVect, int, IVHash> MeshEdgeMap;
#endif
  MeshEdgeMap meshedgemap;

  // read triangles
  float tridat[12];
  RealVect normal;
  Vector<RealVect> verts(3);
  chunksize = sizeof(float)*12 + sizeof(uint16_t);
  int itri=0;
  pout() << "# of triangles: " << m_ntri << endl;
  while (readsize<fsize && itri<m_ntri) // force only m_ntri triangles... what we want?
  {
    CH_START(tdat);

    // read a triangle
    memcpy( tridat , memblock+readsize , sizeof(float)*12 );
    readsize += chunksize; // increment progress in file

    // cast to reals
    normal[0] = (Real) tridat[0];
    normal[1] = (Real) tridat[1];
    normal[2] = (Real) tridat[2];
    verts[0][0] = (Real) tridat[3];
    verts[0][1] = (Real) tridat[4];
    verts[0][2] = (Real) tridat[5];
    verts[1][0] = (Real) tridat[6];
    verts[1][1] = (Real) tridat[7];
    verts[1][2] = (Real) tridat[8];
    verts[2][0] = (Real) tridat[9];
    verts[2][1] = (Real) tridat[10];
    verts[2][2] = (Real) tridat[11];

    // check for degenerate triangles
    if ((verts[0]-verts[1]).vectorLength() < m_stlmesh->tol || \
        (verts[1]-verts[2]).vectorLength() < m_stlmesh->tol || \
        (verts[2]-verts[0]).vectorLength() < m_stlmesh->tol)
    {
      // pout() << "STLBinaryReader: Building mesh: Warning, encountered degenerate triangle\n";
      // pout() << " itri = " << itri << " corners = "; PRV(verts[0]); PRV(verts[1]); PRV(verts[2]); pout() << "\n";
    }

    // initialize triangle
    m_stlmesh->triangles.corners.resize(itri+1);
    m_stlmesh->triangles.corners[itri].resize(3);
    m_stlmesh->triangles.normal.resize(itri+1);

    // insert normal
    m_stlmesh->triangles.normal[itri] = normal;

    for (int i = 0; i < SpaceDim; i++)
      m_stlmesh->triangles.corners[itri][i]=-1; // initialize to -1
        
    CH_STOP(tdat);
    CH_START(tvert);

    // see if vertices exist already
    /*
    bool condition;
    for (int ivertg = 0; ivertg < m_stlmesh->vertices.vertex.size(); ivertg++)
    {
      for (int ivertl = 0; ivertl < SpaceDim; ivertl++)
      {
        // if all dimensions match within tol, then it's the same point
        condition = true;
        for (int j = 0; j < SpaceDim; j++)
          condition = condition && Abs(verts[ivertl][j]-m_stlmesh->vertices.vertex[ivertg][j])<m_stlmesh->tol;
        if (condition)
          m_stlmesh->triangles.corners[itri][ivertl] = ivertg;
      }
    }*/
    for (int ivertl = 0; ivertl < SpaceDim; ivertl++)
    {
      // set up for KDTree search
      for (int j=0; j<SpaceDim; j++)
      {
        rvarray[j] = verts[ivertl][j];
        closestarray[j] = INFINITY;
      }
      int* ivertg_ptr;

      // do the search
      if (treeInsertIndex>0)
      {
        KDError = KDNearestNeighbor( nodetree , rvarray , closestarray , (void**) &ivertg_ptr , NULL , 0 );
        if (KDError!=0)
          pout() << "KDNearestNeighbor returned an error (STLBinaryReader)" << endl;

        // set the triangle's corner to the found node if it is the same as this one
        bool condition = true;
        for (int j = 0; j < SpaceDim; j++)
          condition = condition && Abs(verts[ivertl][j]-m_stlmesh->vertices.vertex[*ivertg_ptr][j])<m_stlmesh->tol;
        if (condition)
          m_stlmesh->triangles.corners[itri][ivertl] = *ivertg_ptr;
      }
    }

    CH_STOP(tvert);
    CH_START(tdat);

    // if vertices don't exist, add them
    for (int ivertl = 0; ivertl < SpaceDim; ivertl++)
    {
      if (m_stlmesh->triangles.corners[itri][ivertl]==-1)
      {
        m_stlmesh->vertices.vertex.push_back(verts[ivertl]);
        m_stlmesh->triangles.corners[itri][ivertl] = m_stlmesh->vertices.vertex.size()-1;
        // and initialize a place in connect.vertexToTriangle
        m_stlmesh->connect.vertexToTriangle.resize( m_stlmesh->vertices.vertex.size() );

        // add to nodetree
        for (int j=0; j<SpaceDim; j++)
          rvarray[j] = verts[ivertl][j];
        (*data)[treeInsertIndex] = m_stlmesh->vertices.vertex.size()-1;
        KDError = KDInsert( nodetree , rvarray , &(*data)[treeInsertIndex] );
        if (KDError!=0)
          pout() << "KDInsert returned an error (STLBinaryReader)" << endl;
        treeInsertIndex++;
        
      }
    }

    CH_STOP(tdat);
    CH_START(tedge);

    // see if edges exist already, if not add them
    /*
    Vector<int> tmpedge(2);
    Vector<int> tmpedg2(2);
    for (int iedgel = 0; iedgel < SpaceDim; iedgel ++)
    {
      // create edges
      tmpedge[0] = m_stlmesh->triangles.corners[itri][ (iedgel)   % 3 ];
      tmpedge[1] = m_stlmesh->triangles.corners[itri][ (iedgel+1) % 3 ];
      tmpedg2[0]=tmpedge[1]; tmpedg2[1]=tmpedge[0]; // swap nodes
      bool foundedge = false; // flag to see if we need to add a new edge

      for (int iedgeg = 0; iedgeg < m_stlmesh->edges.edge.size(); iedgeg++)
      {
        // if all indices are the same, it's the same edge
        if ( m_stlmesh->edges.edge[iedgeg].stdVector() == tmpedge.stdVector() || \
             m_stlmesh->edges.edge[iedgeg].stdVector() == tmpedg2.stdVector() )
        {
          foundedge=true;
          if (m_stlmesh->connect.edgeToTriangle[iedgeg][0]==-1)
            m_stlmesh->connect.edgeToTriangle[iedgeg][0] = itri;
          else if (m_stlmesh->connect.edgeToTriangle[iedgeg][1]==-1)
            m_stlmesh->connect.edgeToTriangle[iedgeg][1] = itri;
          else
          {
            //MayDay::Abort("STLBinaryReader: Building mesh: edge has more than two triangles connected to it");
            pout() << "STLBinaryReader: Building mesh: edge has more than two triangles connected\n";
            printf(" iedgel=%i, iedgeg=%i, itri=%i, connected tri1=%i, connected tri2=%i\n",iedgel,iedgeg,itri,m_stlmesh->connect.edgeToTriangle[iedgeg][0],m_stlmesh->connect.edgeToTriangle[iedgeg][1]);
            //pout() << " old left="; PRV(m_stlmesh->vertices.vertex[ tmpedge[0] ]);
            //pout() << " right="; PRV(m_stlmesh->vertices.vertex[ tmpedge[1] ]); pout() << "\n";
            //pout() << " current left="; PRV(verts[ (iedgel) % 3 ]);
            //pout() << " right="; PRV(verts[ (iedgel+1) % 3 ]); pout() << "\n";
            pout() << " this triangle: "; PRV(verts[0]); PRV(verts[1]); PRV(verts[2]); pout() << "\n";
            int itri1 = m_stlmesh->connect.edgeToTriangle[iedgeg][0];
            pout() << " old triangle1: "; PRV(m_stlmesh->vertices.vertex[ m_stlmesh->triangles.corners[itri1][0] ]);
              PRV(m_stlmesh->vertices.vertex[ m_stlmesh->triangles.corners[itri1][1] ]);
              PRV(m_stlmesh->vertices.vertex[ m_stlmesh->triangles.corners[itri1][2] ]); pout() << "\n";
            itri1 = m_stlmesh->connect.edgeToTriangle[iedgeg][1];
            pout() << " old triangle2: "; PRV(m_stlmesh->vertices.vertex[ m_stlmesh->triangles.corners[itri1][0] ]);
              PRV(m_stlmesh->vertices.vertex[ m_stlmesh->triangles.corners[itri1][1] ]);
              PRV(m_stlmesh->vertices.vertex[ m_stlmesh->triangles.corners[itri1][2] ]); pout() << "\n";
            //m_stlmesh->PrintMesh(); pout() << "\n";
          }
        }
      }

      // add edge to list and update connectivity
      if (!foundedge)
      {
        m_stlmesh->edges.edge.push_back(tmpedge);
        Vector<int> tmpe2t(2);
        tmpe2t[0]=itri; tmpe2t[1]=-1;
        m_stlmesh->connect.edgeToTriangle.push_back(tmpe2t);
      }
    }
    */
    
    // see if edges exist already
    IntVect tmpedge = IntVect::Zero;
    IntVect tmpedge2 = IntVect::Zero;
    for (int iedgel = 0; iedgel < SpaceDim; iedgel++)
    {
      tmpedge[0] = m_stlmesh->triangles.corners[itri][ (iedgel)   % 3 ];
      tmpedge[1] = m_stlmesh->triangles.corners[itri][ (iedgel+1) % 3 ];
      tmpedge2[0]=tmpedge[1]; tmpedge2[1]=tmpedge[0]; // swap nodes

      // do search
      MeshEdgeMap::iterator memit = meshedgemap.find( tmpedge );
      MeshEdgeMap::iterator memit2 = meshedgemap.find( tmpedge2 );
      // set memit to the found edge
      if (memit2 != meshedgemap.end())
      {
        memit = memit2; 
        tmpedge = tmpedge2;
      }

      // found the edge in edgemap
      if (memit != meshedgemap.end())
      {
        int iedgeg = memit->second;
        if (m_stlmesh->connect.edgeToTriangle[iedgeg][0]==-1)
          m_stlmesh->connect.edgeToTriangle[iedgeg][0] = itri;
        else if (m_stlmesh->connect.edgeToTriangle[iedgeg][1]==-1)
          m_stlmesh->connect.edgeToTriangle[iedgeg][1] = itri;
        else
        {
          // pout() << "STLBinaryReader: Building mesh: edge has more than two triangles connected" << endl;
        }
      }
      else // didn't find the edge
      {
        Vector<int> tmpedgevec(2); 
        tmpedgevec[0] = tmpedge[0];
        tmpedgevec[1] = tmpedge[1];
        m_stlmesh->edges.edge.push_back(tmpedgevec);
        Vector<int> tmpe2t(2);
        tmpe2t[0]=itri; tmpe2t[1]=-1;
        m_stlmesh->connect.edgeToTriangle.push_back(tmpe2t);
        // add to meshedgemap
        pair<MeshEdgeMap::iterator, bool> rval = meshedgemap.insert(make_pair(tmpedge,m_stlmesh->edges.edge.size()-1));
        if (rval.second==false)
        {
          // pout() << "STLBinaryReader: Building mesh: we thought the edge wasn't in the map, but it is" << endl;
        }
      }
    }

    CH_STOP(tedge);
    CH_START(tdat);

    // add vertex to triangle connectivity
    m_stlmesh->connect.vertexToTriangle[ m_stlmesh->triangles.corners[itri][0] ].push_back(itri);
    m_stlmesh->connect.vertexToTriangle[ m_stlmesh->triangles.corners[itri][1] ].push_back(itri);
    m_stlmesh->connect.vertexToTriangle[ m_stlmesh->triangles.corners[itri][2] ].push_back(itri);

    CH_STOP(tdat);

    itri++;
  }
  pout() << "# of vertices: " << m_stlmesh->vertices.vertex.size() << endl;

  // get rid of the KDTree
  delete data;
  KDError = KDFree( nodetree );
  if (KDError!=0)
    pout() << "KDFree returned an error (STLBinaryReader)" << endl;
  KDTreeFinalize();

  //memblock -= fsize; // go back to beginning of memory block
  free(memblock);
  m_ntriMatch = m_stlmesh->triangles.corners.size()==m_ntri;
  //m_stlmesh->PrintMesh();

}

#include "NamespaceFooter.H"
