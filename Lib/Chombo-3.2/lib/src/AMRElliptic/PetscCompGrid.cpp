#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifdef CH_USE_PETSC

#include "PetscCompGrid.H"
#include "IntVectSet.H"

#include "NamespaceHeader.H"
 
std::ostream& operator<< (std::ostream& os, GID_type a_type)
{
  switch (a_type)
    {
    case GHOST: os << "ghost";                       break;
    case FINE_COVERED: os << "covered (with fine)" ; break;
    case DIRI_BC: os << "BC (Diri)";                 break;
    case NEUM_BC: os << "BC (Neum)";                 break;
    case ANY_BC: os << "BC";                         break;
    case UNKNOWN: os << "unknown";                   break;
    default: os << "index " << (int)a_type;          break;
    }
  return os;
}

//
// My BC function
//
CompBC::~CompBC()
{
  PetscFree(m_Rcoefs);
}

CompBC::CompBC(int a_order, IntVect a_nGhosts) : m_Rcoefs(0), m_isDefined(false)
{
  define(a_order,a_nGhosts);
}

void 
CompBC::define(int a_order, IntVect a_nGhosts)
{  
  if (m_Rcoefs) PetscFree(m_Rcoefs);
  if (a_nGhosts[0]!=1 && a_nGhosts[0]!=2)
    MayDay::Error("Unsupported number of ghosts in CompBC");

  m_nGhosts = a_nGhosts; m_nSources = a_order+1;
  PetscMalloc(m_nSources*m_nGhosts[0]*sizeof(PetscReal),&m_Rcoefs);

  m_isDefined = false; // needs to be set
}

PetscReal CompBC::getCoef(int a_iSrc, int a_iGhost)
{
  if (!m_isDefined) createCoefs();
  return m_Rcoefs[a_iGhost*m_nSources + a_iSrc];
}
//
void 
ConstDiriBC::createCoefs()
{  
  m_isDefined=true;

  if (m_nSources==1) 
    {
      m_Rcoefs[0] = -1.;
      if (m_nGhosts[0]==2) m_Rcoefs[1] = -3.;   
    }
  else if (m_nSources==2) 
    { // s = 6, 18
      m_Rcoefs[0] = -5./2; m_Rcoefs[1] = 1./2;
      if (m_nGhosts[0]==2) 
        {
          m_Rcoefs[2] = -21./2; m_Rcoefs[3] = 5./2;
        }
    }
  else if (m_nSources==3) 
    { // s = 12, 48
      m_Rcoefs[0] = -13./3; m_Rcoefs[1] = 5./3; m_Rcoefs[2] = -1./3;
      if (m_nGhosts[0]==2) 
        {
          m_Rcoefs[3] = -70./3; m_Rcoefs[4] = 32./3;  m_Rcoefs[5] = -7./3; 
        }
    }
  else if (m_nSources==4) 
    { // s = 60, 300
      m_Rcoefs[0] = -77./12; m_Rcoefs[1] = 43./12; m_Rcoefs[2] = -17./12; m_Rcoefs[3] = 3./12;
      if (m_nGhosts[0]==2) 
        {
          m_Rcoefs[4] = -505./12; m_Rcoefs[5] = 335./12;  m_Rcoefs[6] = -145./12; m_Rcoefs[7] = 27./12; 
        }
    }
  else
    MayDay::Error("Unsupported degree in ConstDiriBC");
}

void 
ConstDiriBC::operator()( FArrayBox&           a_state,
                         const Box&           a_valid,
                         const ProblemDomain& a_domain,
                         Real                 a_dx,
                         bool                 a_homogeneous)
{
  const Box& domainBox = a_domain.domainBox();
  
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (!a_domain.isPeriodic(idir))
        {
          for (SideIterator sit; sit.ok(); ++sit)
            {
              Side::LoHiSide side = sit();              
              if (a_valid.sideEnd(side)[idir] == domainBox.sideEnd(side)[idir])
                {
                  // Dirichlet BC
                  int isign = sign(side);
                  Box toRegion = adjCellBox(a_valid, idir, side, 1);
                  toRegion &= a_state.box();
                  for (BoxIterator bit(toRegion); bit.ok(); ++bit)
                    {
                      IntVect ivTo = bit();                     
                      IntVect ivClose = ivTo - isign*BASISV(idir);
                      for (int ighost=0;ighost<m_nGhosts[0];ighost++,ivTo += isign*BASISV(idir))
                        {
                          for (int icomp = 0; icomp < a_state.nComp(); icomp++) a_state(ivTo, icomp) = 0.0;
                          IntVect ivFrom = ivClose;
                          for (int i=0;i<m_nSources;i++,ivFrom -= isign*BASISV(idir))
                            {
                              for (int icomp = 0; icomp < a_state.nComp(); icomp++)
                                {
                                  a_state(ivTo, icomp) += m_Rcoefs[ighost*m_nSources + i]*a_state(ivFrom, icomp);
                                }
                            }
                        }
                    }
                } // if ends match
            } // end loop over sides
        } // if not periodic in this direction
    } // end loop over directions    
}

//
// PetscCompGrid
//
void
PetscCompGrid::clean()
{
  if (m_mat) 
    {
      MatDestroy(&m_mat);
      m_mat = 0; // this can get called multiple times
      if (m_domains.size() > 1)
        {
          const Box& stencilBox = m_FCStencils.box();
          for (BoxIterator bit(stencilBox); bit.ok(); ++bit)
            {
              IntVect offset = bit();
              delete m_FCStencils(offset, 0);
            }   
        }
    }
}

//
PetscCompGrid::~PetscCompGrid()
{
  clean();
}

// Petsc composite grid solver - a matrix with solve methods that builds itself with a Chombo operator and a hierarchy of grids.
void
PetscCompGrid::define( const ProblemDomain &a_cdomain,
                       Vector<DisjointBoxLayout> &a_grids, 
                       Vector<int> &a_refRatios, 
                       BCHolder a_bc,
                       const RealVect &a_cdx,
                       int a_numLevels/* =-1 */, int a_ibase/* =0 */)
{
  CH_TIME("PetscCompGrid::define");
  int maxiLev;
  if (a_numLevels<=0) maxiLev = a_grids.size() - 1;
  else maxiLev = a_numLevels - 1;
  const int numLevs = maxiLev - a_ibase + 1;
  m_bc = a_bc;

  PetscCompGrid::clean(); // this is virtual so lets not step on derived classes data
  
  if (m_verbose>5) 
    {
      pout() << "PetscCompGrid::define: a_numLevels=" << a_numLevels << ", numLevs=" << 
	numLevs << ", a_ibase=" << a_ibase << endl; 
    }
  if (numLevs>1)
    {
      int degree = 3, nref = a_refRatios[0]; // assume same coarsening for all directions and levels
      IntVect interpUnit = IntVect::Unit;
      Box stencilBox( -m_CFStencilRad*interpUnit,
                      m_CFStencilRad*interpUnit );
      m_FCStencils.define(stencilBox,1);
      for (BoxIterator bit(stencilBox); bit.ok(); ++bit)
        {
          IntVect offset = bit();
          m_FCStencils(offset, 0) =
            new FourthOrderInterpStencil(offset, nref, degree);
        }
    }

  // copy in data
  m_domains.resize(numLevs);
  m_grids.resize(numLevs);
  m_refRatios.resize(numLevs);
  m_GIDs.resize(numLevs);
  m_crsSupportGIDs.resize(numLevs);
  m_fineCoverGIDs.resize(numLevs);
  m_dxs.resize(numLevs);

  RealVect dx = a_cdx;
  ProblemDomain dom = a_cdomain;
  for (int ilev=0; ilev<a_ibase ; ilev++) 
    {
      dom.refine(a_refRatios[ilev]);
      dx /= a_refRatios[ilev];
    }
  for (int ilev=a_ibase, ii=0; ilev < numLevs ; ilev++, ii++)
    {
      m_domains[ii] = dom;
      m_grids[ii] = a_grids[ilev];
      m_dxs[ii] = dx;
      if (ilev != numLevs-1) 
        {
          m_refRatios[ii] = a_refRatios[ilev];
          dom.refine(a_refRatios[ilev]);
          dx /= a_refRatios[ilev];
        }
      if (m_verbose>5) 
	{
	  pout() << "PetscCompGrid::define: level=" << ilev << ", " << 
	    m_grids[ii].size() << " patches" << endl; 
	}
    }

  // iterate over real cells, get gids, make ops, start numbering from coarsest
  PetscInt my0 = 0;
  for (int ilev=0;ilev<numLevs;ilev++)
    { 
      const DisjointBoxLayout& dbl = m_grids[ilev];
      // 3 is a hack for treb (?); refRat*(rad+1) happens at corners (not clear def of "radius")
      IntVect nProcessGhost = (ilev==0) ? /*getGhostVect()*/ 3*IntVect::Unit : 
	m_refRatios[ilev-1]*(m_CFStencilRad+1)*IntVect::Unit;
      m_GIDs[ilev] = RefCountedPtr<LevelData<BaseFab<PetscInt> > >
        (new LevelData<BaseFab<PetscInt> >(dbl,1,nProcessGhost));
      // set gids == GHOST
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
        {
          BaseFab<PetscInt>& gidfab = (*m_GIDs[ilev])[dit];
          gidfab.setVal(GHOST); 
        }

      // zero covered if not finest
      if (ilev!=numLevs-1) 
        {
          // zero out fine cover
          DisjointBoxLayout dblCoarsenedFine;
          Copier copier;
          coarsen(dblCoarsenedFine, m_grids[ilev+1], m_refRatios[ilev]); // coarsens entire grid
          copier.define(dblCoarsenedFine, dbl, IntVect::Zero);
          for (CopyIterator it(copier, CopyIterator::LOCAL); it.ok(); ++it)
            {
              const MotionItem& item = it();
              (*m_GIDs[ilev])[item.toIndex].setVal(FINE_COVERED, item.toRegion, 0);           
            }
          for (CopyIterator it(copier, CopyIterator::TO); it.ok(); ++it)
            {
              const MotionItem& item = it();
              (*m_GIDs[ilev])[item.toIndex].setVal(FINE_COVERED, item.toRegion, 0);
            }
        }
      
      // add extra covered
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
        {
          BaseFab<PetscInt>& gidfab = (*m_GIDs[ilev])[dit]; // make explicit that we modify this
	  addExtraCovered(FINE_COVERED,ilev,dit(),gidfab);
	}
      // count
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
        {
          Box region = dbl[dit]; // no ghosts
          BaseFab<PetscInt>& gidfab = (*m_GIDs[ilev])[dit];
          for (BoxIterator bit(region); bit.ok(); ++bit)
            {
              const IntVect& iv = bit(); 
              if (gidfab(iv,0) == GHOST) my0++; // not covered, so real
            }
        }
    }

  // create glabal matrix, not used until close
#ifdef CH_MPI
  MPI_Comm wcomm = Chombo_MPI::comm;
#else
  MPI_Comm wcomm = PETSC_COMM_SELF;
#endif
  MatCreate(wcomm,&m_mat);
  MatSetSizes(m_mat,my0*m_dof,my0*m_dof,PETSC_DECIDE,PETSC_DECIDE);
  MatSetBlockSize(m_mat,m_dof);
  MatSetType(m_mat,MATAIJ);
  MatSetFromOptions(m_mat);
#ifdef CH_MPI
  PetscInt result, petscdata = my0;
  MPI_Datatype mtype;
  PetscDataTypeToMPIDataType(PETSC_INT,&mtype);
  MPI_Scan( &petscdata, &result, 1, mtype, MPI_SUM, Chombo_MPI::comm );
  m_gid0 = result - my0;
#else
  m_gid0 = 0;
#endif

  // set GIDs
  my0 = m_gid0;
  for (int ilev=0;ilev<numLevs;ilev++)
    {
      const DisjointBoxLayout& dbl = m_grids[ilev];
      // count and set ids
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
        {
          Box region = dbl[dit]; // no ghosts
          BaseFab<PetscInt>& gidfab = (*m_GIDs[ilev])[dit];
          for (BoxIterator bit(region); bit.ok(); ++bit)
            {
              const IntVect& iv = bit(); 
              if (gidfab(iv,0) == GHOST) gidfab(iv,0) = my0++; // not covered, so real
              else CH_assert(gidfab(iv,0) == FINE_COVERED);   // fine covered
            }
        }
      m_GIDs[ilev]->exchange();
    }

  // construct covering gids maps
  for (int ilev=1;ilev<numLevs;ilev++)
    {
      Interval inter = m_GIDs[ilev]->interval();
      DisjointBoxLayout crsenedFineDBL,refinedCrsDBL;
      LevelData<BaseFab<PetscInt> > *pl;
      int refRatio = m_refRatios[ilev-1];

      // form m_crsSupportGIDs[ilev-1]:
      const DisjointBoxLayout& fdbl = m_grids[ilev];
      coarsen(crsenedFineDBL,fdbl,refRatio);
      // rad+1 needed to get stencils coarse cell coarse cover
      const IntVect nCrsSuppCells = (m_CFStencilRad+1)*IntVect::Unit;
      pl = new LevelData<BaseFab<PetscInt> >(crsenedFineDBL,1,nCrsSuppCells);
      m_crsSupportGIDs[ilev-1] = RefCountedPtr<LevelData<BaseFab<PetscInt> > >(pl);
      // set gids == UNKNOWN
      for (DataIterator dit = fdbl.dataIterator(); dit.ok(); ++dit)
        {
          BaseFab<PetscInt>& gidfab = (*pl)[dit];
          gidfab.setVal(UNKNOWN); // enum this.  BCs and covered by this 
        }
      Copier rcopier(m_grids[ilev-1],crsenedFineDBL,nCrsSuppCells);
      m_GIDs[ilev-1]->copyTo(inter,*m_crsSupportGIDs[ilev-1],inter,rcopier);

      // from m_fineCoverGIDs[ilev]: fine support of stencils for coarse grid ilev-1
      const DisjointBoxLayout& cdbl = m_grids[ilev-1];
      refine(refinedCrsDBL,cdbl,refRatio);
      //
      if (getGhostVect()[0]>refRatio)
	{
	  MayDay::Error("PetscCompGrid::define getGhostVect>refRatio - not setup for this");
	}
      // ghost cells for fine support governed by CF radius (12=nr*nr(rad+1), could use 1 with nesting radius==4 maybe
      const IntVect nFineProcGhosts = refRatio*refRatio*(m_CFStencilRad+1)*IntVect::Unit;
      pl = new LevelData<BaseFab<PetscInt> >(refinedCrsDBL,1,nFineProcGhosts); 
      m_fineCoverGIDs[ilev] = RefCountedPtr<LevelData<BaseFab<PetscInt> > >(pl);
      // set gids == UNKNOWN
      for (DataIterator dit = cdbl.dataIterator(); dit.ok(); ++dit)
        {
          BaseFab<PetscInt>& gidfab = (*pl)[dit];
          gidfab.setVal(UNKNOWN); // enum this. fine grid that is not there (true data here)
        }
      Copier pcopier(fdbl,refinedCrsDBL,nFineProcGhosts);
      m_GIDs[ilev]->copyTo(inter,*m_fineCoverGIDs[ilev],inter,pcopier);
      if (m_verbose>3) 
	{
          pout() << "PetscCompGrid::define: m_fineCoverGIDs["<< ilev-1 <<"] boxes: " << endl; 
	  for (DataIterator dit = cdbl.dataIterator(); dit.ok(); ++dit)
	    {
	      BaseFab<PetscInt>& gidfab = (*pl)[dit];
              if (m_verbose>5)
                {
                  pout() << "PetscCompGrid::define: refined box="
                         << gidfab.box() 
                         << ", coarse box = " << m_grids[ilev-1][dit] << endl;
                }
	    }
	}
    }
}

// main method of this class -- make a matrix
#undef __FUNCT__
#define __FUNCT__ "createMatrix"
PetscErrorCode
PetscCompGrid::createMatrix()
{
  CH_TIME("PetscCompGrid::createMatrix");
  PetscErrorCode ierr,ilev,idx;
  PetscInt nloc,max_size;
  int nGrids=m_domains.size();
  IntVect nghost = getGhostVect();
  Vector<StencilTensor> stenVect;
  PetscInt *d_nnz, *o_nnz;
  PetscFunctionBeginUser;

  // set preallocation
#if defined(PETSC_USE_LOG)
  PetscLogEventBegin(m_event0,0,0,0,0);
#endif
  ierr = MatGetLocalSize(m_mat,&nloc,PETSC_NULL);CHKERRQ(ierr);
  nloc /= m_dof;
  /* count nnz */
  PetscMalloc((nloc+1)*sizeof(PetscInt)*m_dof, &d_nnz);
  PetscMalloc((nloc+1)*sizeof(PetscInt)*m_dof, &o_nnz);
  stenVect.resize(nloc); 
  // add data
  for (ilev=0,idx=0;ilev<nGrids;ilev++)
    {
      const DisjointBoxLayout& dbl = m_grids[ilev];
      Box testbox,dombox(m_domains[ilev].domainBox());
      for (int dir=0; dir<CH_SPACEDIM; ++dir) 
        {
          if (m_domains[ilev].isPeriodic(dir)) {
            dombox.grow(dir,nghost[dir]); // nullify
          }
        }
      testbox = dombox;
      testbox.grow(-nghost); // non-boundary cells
      max_size=0;
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
        { 
          BaseFab<PetscInt>& gidfab = (*m_GIDs[ilev])[dit];
          Box region = dbl[dit];
          for (BoxIterator bit(region); bit.ok(); ++bit)
            {
              const IntVect& iv = bit();
              if ( gidfab(iv,0) >= 0)
                { 
                  // doit 
                  StencilTensor &sten2 = stenVect[idx];
                  createOpStencil(iv,ilev,dit(),sten2);
                  if (!testbox.contains(iv))
                    {
                      applyBCs(iv,ilev,dit(),dombox,sten2);
                    }
                  if (ilev != 0)
                    {
		      InterpToCoarse(iv,ilev,dit(),sten2);
                    }
                  if (ilev != nGrids-1)
                    {
                      InterpToFine(iv,ilev,dit(),sten2);
                    }
                  // collect prealloc sizes
                  o_nnz[idx*m_dof] = d_nnz[idx*m_dof] = sten2.size();
                  if (o_nnz[idx*m_dof]>max_size) { max_size = o_nnz[idx*m_dof];}
                  if (d_nnz[idx*m_dof]>nloc) d_nnz[idx*m_dof] = nloc;
                  d_nnz[idx*m_dof] *= m_dof;
                  o_nnz[idx*m_dof] *= m_dof;
                  for (int i=1; i<m_dof; ++i) 
                    {
                      d_nnz[idx*m_dof+i] = d_nnz[idx*m_dof];
                      o_nnz[idx*m_dof+i] = o_nnz[idx*m_dof];
                    }
                  idx++;
                } // real
            } // cell
        } // patch
      if (m_verbose>0) 
        {
          if (m_verbose>1) 
            {
#ifdef CH_MPI
              MPI_Comm wcomm = Chombo_MPI::comm;
              PetscInt n = max_size;
              ierr = MPI_Allreduce(&n, &max_size, 1, MPIU_INT, MPI_MAX, wcomm);CHKERRQ(ierr);
#endif        
            }
          if (m_verbose>3)
            {
              pout() << "\t PetscCompGrid::createMatrix level " << ilev+1 << 
                "/" << nGrids << ". domain " << m_domains[ilev] 
                     << ". max. stencil size: " 
                     << max_size << endl;
            }
        }
    } // level
  CH_assert(idx==nloc);
  
  ierr = MatSeqAIJSetPreallocation(m_mat,0,d_nnz);CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(m_mat,0,d_nnz,0,o_nnz);CHKERRQ(ierr);

  // assemble
  for (ilev=0,idx=0;ilev<nGrids;ilev++)
    {
      const DisjointBoxLayout& dbl = m_grids[ilev];
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
        {
          BaseFab<PetscInt>& gidfab = (*m_GIDs[ilev])[dit];
          Box region = dbl[dit];
          for (BoxIterator bit(region); bit.ok(); ++bit)
            {
              const IntVect& iv = bit();
              if ( gidfab(iv,0) >= 0)
                {
                  ierr = AddStencilToMat(iv,ilev,dit(),stenVect[idx],m_mat);CHKERRQ(ierr);
                  idx++;
                }
            }
        }
    }
  CH_assert(idx==nloc);
  PetscFree(d_nnz);  PetscFree(o_nnz);

  // total assesmbly
  ierr = MatAssemblyBegin(m_mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(m_mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  
  if (m_writeMatlab)
    {
      PetscViewer viewer;
      char suffix[30];
      sprintf(suffix, "A%d.m",nGrids);
      ierr = PetscViewerASCIIOpen(((PetscObject) m_mat)->comm, suffix, &viewer);CHKERRQ(ierr);
      ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
      ierr = MatView(m_mat,viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&viewer);
    }
#if defined(PETSC_USE_LOG)
  PetscLogEventEnd(m_event0,0,0,0,0);
#endif
  PetscFunctionReturn(0);
}

static bool s_check_row_sum = false;

void 
PetscCompGrid::applyBCs( IntVect a_iv, int a_ilev, const DataIndex &a_dummy, Box a_dombox, 
                         StencilTensor &a_sten )
{
  CH_TIME("PetscCompGrid::applyBCs");
  int nSources = -1;
  Vector<Vector<StenScalarNode> > new_vals;
  
  // count degrees, add to 'corners'
  Vector<IntVectSet> corners(CH_SPACEDIM);
  StencilTensor::const_iterator end2 = a_sten.end(); 
  for (StencilTensor::const_iterator it = a_sten.begin(); it != end2; ++it) 
    {
      const IntVect &jiv = it->first.iv();
      if (!a_dombox.contains(jiv) && it->first.level()==a_ilev) 
        {
          int degree = CH_SPACEDIM-1;
          for (int dir=0;dir<CH_SPACEDIM;++dir)
            {
              if (jiv[dir] >= a_dombox.smallEnd(dir) && jiv[dir] <= a_dombox.bigEnd(dir)) degree--;
            }
          CH_assert(degree>=0);
          corners[degree] |= jiv;
        }
    }
  // move ghosts starting at with high degree corners and cascade to lower degree
  for (int ideg=CH_SPACEDIM-1;ideg>=0;ideg--)
    {
      for (IVSIterator ivit(corners[ideg]); ivit.ok(); ++ivit)
        {
          const IntVect &jiv = ivit(); // get rid of this bc ghost
          // get layer of ghost
          Box gbox(IntVect::Zero,IntVect::Zero); gbox.shift(jiv);
          Box inter = a_dombox & gbox;    CH_assert(inter.numPts()==0);
          int igid = -1; // find which layer of ghost I am
          do{
            igid++;
            gbox.grow(1);
            inter = a_dombox & gbox;
          }while (inter.numPts()==0);
          for (int dir=0;dir<CH_SPACEDIM;++dir)
            {
              if (jiv[dir] < a_dombox.smallEnd(dir) || jiv[dir] > a_dombox.bigEnd(dir)) 
		{
		  // have a BC, get coefs on demand
		  if (nSources == -1) { 
		    IntVect ghostVect = getGhostVect();
		    RefCountedPtr<BCFunction> bc = m_bc.getBCFunction();
		    ConstDiriBC *mybc = dynamic_cast<ConstDiriBC*>(&(*bc));
		    if (mybc)
		      {
			CH_assert(mybc->nGhosts()==ghostVect);
			nSources = mybc->nSources();
			new_vals.resize(ghostVect[0]);
			for (int i = 0 ; i < ghostVect[0]; i++)
			  {
			    new_vals[i].resize(nSources); 
			    for (int j = 0 ; j < nSources; j++ )
			      {
				new_vals[i][j].second.setValue(mybc->getCoef(j,i));
				new_vals[i][j].first.setLevel(a_ilev);
			      }
			  }
		      }
		    else
		      {
			CH_assert(ghostVect[0]==1); 
			Real denom=2.;   new_vals.resize(1); new_vals[0].resize(2);
			new_vals[0][0].second.setValue(-5./denom); new_vals[0][1].second.setValue(1./denom);
			new_vals[0][0].first.setLevel(a_ilev);  new_vals[0][1].first.setLevel(a_ilev);
		      }		    
		    // for debugging -- neumann
		    if (s_check_row_sum)
		      {
			new_vals.resize(ghostVect[0]);
			for (int idx=0;idx<ghostVect[0];idx++)
			  {
			    new_vals[idx].resize(1);
			    new_vals[idx][0].second.setValue(1.);
			  }
			nSources = 1;
		      }
		  } // init coefs
		  
                  IndexML kill(jiv,a_ilev);
                  int isign = 1;
                  if (jiv[dir] < a_dombox.smallEnd(dir)) isign = -1;
                  IntVect biv = jiv;
                  biv.shift(dir,-isign*(igid+1));                 
                  for (int j=0; j<nSources; j++, biv.shift(dir,-isign))
                    {
                      new_vals[igid][j].first.setIV(biv);
                      if (ideg>0 && !a_dombox.contains(biv)) // this could be a new stencil value for high order BCs
                        {
                          corners[ideg-1] |= biv;
                        }
                      else CH_assert(a_dombox.contains(biv));
                    }
                  StencilProject(kill,new_vals[igid],a_sten);
                  int nrm = a_sten.erase(kill); CH_assert(nrm==1);
                  break;
                }
            }
        } // ghosts
    } // degree
  
  // debug
  if (0)
    {
      double summ=0.;
      StencilTensor::const_iterator end = a_sten.end(); 
      for (StencilTensor::const_iterator it = a_sten.begin(); it != end; ++it) {
        for (int i=0; i<m_dof; ++i) 
          {
            for (int j=0; j<m_dof; ++j) 
              {
                summ += it->second.value(i,j);
              }
          }
      }
      if (Abs(summ)>1.e-10)
        {
          for (StencilTensor::const_iterator it = a_sten.begin(); it != end; ++it) {
            pout() << it->first << ": \n";
            for (int i=0; i<m_dof; ++i) 
              {
                pout() << "\t\t";
                for (int j=0; j<m_dof; ++j) 
                  {
                    pout() << it->second.value(i,j) << ", ";
                  }
                pout() << endl;
              }
          }
          pout() << "\t ERROR summ (if Neumann): " << summ << endl; 
        }
    }
}

// BCs have been removed (should not care).
//
void 
PetscCompGrid::InterpToCoarse(IntVect a_iv,int a_ilev,const DataIndex &a_di, StencilTensor &a_sten)
{
  CH_TIME("PetscCompGrid::InterpToCoarse");
  int nref=m_refRatios[a_ilev-1];
  int nfine = (int)pow((Real)nref,SpaceDim);
  BaseFab<PetscInt>& supgidfab = (*m_crsSupportGIDs[a_ilev-1])[a_di];
  BaseFab<PetscInt>& gidfab = (*m_GIDs[a_ilev])[a_di];
  ProblemDomain cdom = m_domains[a_ilev-1];
  std::list<StencilTensor::iterator> killList;

  // remove ghosts on domain
  StencilTensor::iterator end = a_sten.end(); 
  for (StencilTensor::iterator it = a_sten.begin(); it != end; ++it) 
    {
      const IndexML &mlivJ = it->first;
      if (mlivJ.level()==a_ilev && gidfab(mlivJ.iv(),0) == GHOST)
        {
          const IntVect &jiv = mlivJ.iv(); // kill
          IntVect cjiv = jiv; cjiv.coarsen(nref);
          const IntVect fivBase = cjiv*nref;
          const IntVect fividx = jiv - fivBase;
          IntVect bndryIndex = getCFStencil(cdom,cjiv);
          const FourthOrderInterpStencil *interp = m_FCStencils(bndryIndex,0);
          const FArrayBox &coarseToFineFab = interp->m_coarseToFineFab;
          const Vector<int> &coarseBaseIndices = interp->m_coarseBaseIndices;
          Vector<StenScalarNode> new_vals(interp->m_stencilSize);
          for (PetscInt cidx=0,kk=0;cidx < interp->m_stencilSize;cidx++)
            {
              IntVect civ = cjiv;
              for (int dir=0;dir<CH_SPACEDIM;++dir,kk++) 
                {
                  civ[dir] += coarseBaseIndices[kk];              
                }
              Real val = coarseToFineFab(fividx,cidx);
              NodeDefine(new_vals[cidx],civ,a_ilev-1,val);
            }
          StencilProject(mlivJ,new_vals,a_sten);
          killList.push_back (it);
        }
    }

  // get rid of coarse covered with simple averaging
  end = a_sten.end(); 
  for (StencilTensor::iterator it = a_sten.begin(); it != end; ++it) 
    {
      const IndexML &mlivJ = it->first;
      IntVect civ = mlivJ.iv();
      if (mlivJ.level()==a_ilev-1 && supgidfab(civ,0)<0)
        {
          CH_assert(supgidfab(civ,0) == FINE_COVERED);    
          Real fineInterp = 1./(Real)nfine;
          const IntVect fivBase = civ*nref; // fine box low end to interp to
          Box fbox(IntVect::Zero,(nref-1)*IntVect::Unit); fbox.shift(fivBase);
          Vector<StenScalarNode> new_vals(nfine);
          BoxIterator bit(fbox);
          for (int i=0; bit.ok(); ++bit,i++)
            {
              NodeDefine(new_vals[i],bit(),a_ilev,fineInterp);
            }
          StencilProject(mlivJ,new_vals,a_sten);
          killList.push_back (it);
        }
    }
  // erase killed
  for (std::list<StencilTensor::iterator>::iterator lit=killList.begin(); 
       lit != killList.end(); 
       ++lit)
    {
      a_sten.erase(*lit);
    }
}

// F: talk to fine grid covers -- refluxing
void 
PetscCompGrid::InterpToFine(IntVect a_iv, int a_ilev,const DataIndex &a_di, StencilTensor &a_sten )
{
  CH_TIME("PetscCompGrid::InterpToFine");
  int nref=m_refRatios[a_ilev];
  BaseFab<PetscInt>& gidfab = (*m_GIDs[a_ilev])[a_di];
  BaseFab<PetscInt>& covergidfab = (*m_fineCoverGIDs[a_ilev+1])[a_di];
  Box coverfabbox = covergidfab.box();
  std::list<StencilTensor::iterator> killList;
  
  // get rid of coarse covered
  StencilTensor::iterator end = a_sten.end(); 
  for (StencilTensor::iterator it = a_sten.begin(); it != end; ++it) 
    {
      const IndexML &mlivJ = it->first;
      const IntVect &jiv = mlivJ.iv();
      if (mlivJ.level()==a_ilev && gidfab(jiv,0)<0)     
        {
          // simple averaging
          const IntVect fivBase = jiv*nref; // fine box low end to interp to
          Box fbox(IntVect::Zero,(nref-1)*IntVect::Unit); fbox.shift(fivBase);
          int nfine = (int)pow((Real)nref,SpaceDim);
          Vector<StenScalarNode> new_vals(nfine);
          Real fineInterp = 1./(Real)nfine;
          BoxIterator bit(fbox);
          for (int i=0; bit.ok(); ++bit,i++)
            {
              NodeDefine(new_vals[i],bit(),a_ilev+1,fineInterp);
            }
          StencilProject(mlivJ,new_vals,a_sten);
          killList.push_back (it);
        }
    }
  // erase killed - would like to have StencilProject do this
  for (std::list<StencilTensor::iterator>::iterator lit=killList.begin(); 
       lit != killList.end(); 
       ++lit)
    {
      a_sten.erase(*lit);
    }
}

// add row of matrix (stencil)
#undef __FUNCT__
#define __FUNCT__ "AddStencilToMat"
PetscErrorCode 
PetscCompGrid::AddStencilToMat(IntVect a_iv, int a_ilev,const DataIndex &a_di, StencilTensor &a_sten, Mat a_mat)
{
  PetscErrorCode ierr;
  PetscScalar vals[4096];
  PetscInt cols[4096],iLevel=a_ilev,vidx[STENCIL_MAX_DOF],gid,ncols=m_dof*a_sten.size();
  BaseFab<PetscInt>& this_gidfabJ = (*m_GIDs[iLevel])[a_di];
  double summ=0.,abssum=0.;
  PetscFunctionBeginUser;

  if (a_sten.size()*m_dof > 4096) MayDay::Error("PetscCompGrid::AddStencilToMat buffer (4096) too small");
  // get cols & vals
  gid = this_gidfabJ(a_iv,0)*m_dof;
  for (int ni=0;ni<m_dof;ni++,gid++) vidx[ni] = gid;
  StencilTensor::const_iterator end = a_sten.end(); 
  int ci=0,jj=0;
  for (StencilTensor::const_iterator it = a_sten.begin(); it != end; ++it,jj++) 
    {
      const IndexML &ivJ = it->first;
      int jLevel = ivJ.level();
      BaseFab<PetscInt>& gidfabJ = (jLevel==iLevel) ? this_gidfabJ : 
        (jLevel==iLevel+1) ? (*m_fineCoverGIDs[iLevel+1])[a_di] : 
        (*m_crsSupportGIDs[iLevel-1])[a_di];
                  
      if (!gidfabJ.box().contains(ivJ.iv())){
      	pout() << "ERROR row  " << IndexML(a_iv,a_ilev) << ", this box:" << 
      	  m_grids[a_ilev][a_di] << ", fine box (failure): " << gidfabJ.box() << 
      	  " failed to find "<< ivJ << endl;
      	jj=0;
      	for (StencilTensor::const_iterator it2 = a_sten.begin(); it2 != end; ++it2) pout()<<++jj<<") j="<<it2->first<<endl;
      	MayDay::Error("PetscCompGrid::AddStencilToMat failed to find cell");
      }

      PetscInt gidj = gidfabJ(ivJ.iv(),0)*m_dof;
      if (gidj<0) 
        {
          pout() << "\tFAILED TO FIND iv=" << ivJ << "\t gidj type:" << (GID_type)gidj << "\tiLevel=" << iLevel << "\tiv=" << a_iv << endl;
          MayDay::Error("PetscCompGrid::AddStencilToMat failed to find gid");
        }
      const Real *vv = it->second.getVals();
      for (int nj=0;nj<m_dof;nj++,gidj++,ci++) cols[ci] = gidj;   // columns
      for (int ni=0;ni<m_dof;ni++) 
        {
          for (int nj=0;nj<m_dof;nj++) 
            {
              double tt = vv[ni*m_dof + nj];
              vals[ni*ncols + jj*m_dof + nj] = tt;
              summ += tt;
              abssum += Abs(tt);
            }
        }
    }
  
  ierr = MatSetValues(a_mat,m_dof,vidx,ncols,cols,vals,INSERT_VALUES);CHKERRQ(ierr);

  // debug
  if (s_check_row_sum && Abs(summ)/abssum>1.e-5)
    {
      pout() << " error gid=" << vidx[0]+1 << " row iv =" << a_iv << endl;
      for (int idx=0;idx<ncols;idx++)
        {
          if (vals[idx]!=0.0)
            {
              PetscInt gidj = cols[idx];
              pout() << "\t" << gidj+1 << ", v = " << vals[idx] << endl;
            }
        }
      SETERRQ4(PETSC_COMM_WORLD,1,"sum=%e on level %d, iv=[%d,%d]",summ,iLevel,a_iv[0],a_iv[1]);
    }

  PetscFunctionReturn(0);
}

//
IntVect PetscCompGrid::getCFStencil(const ProblemDomain &a_cdom, const IntVect a_ivc)
{
  const Box& coarseDomainBox = a_cdom.domainBox();
  const IntVect& coarseDomainLo = coarseDomainBox.smallEnd();
  const IntVect& coarseDomainHi = coarseDomainBox.bigEnd();
  IntVect dist = IntVect::Zero;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (!a_cdom.isPeriodic(idir))
        {
          int offLo = coarseDomainLo[idir] - a_ivc[idir] - 1;
          int offHi = coarseDomainHi[idir] - a_ivc[idir] + 1;
          if (offLo < 0 && offHi >0) // condition means ivc in coarseDomain
            {
              if ((offLo >= -m_CFStencilRad) &&
                  (offHi <= m_CFStencilRad))
                { // both of these:  very narrow domain, you are in trouble
                  MayDay::Error("FourthOrderFineInterp::define bad boxes");
                }
              if (offLo >= -m_CFStencilRad) // -1 or -2
                dist[idir] = offLo;
              if (offHi <= m_CFStencilRad) // 1 or 2
                dist[idir] = offHi;
              // Otherwise, dist[idir] = 0.
            }
        }
    }
  return dist;
}
#undef __FUNCT__
#define __FUNCT__ "putChomboInPetsc"
PetscErrorCode 
PetscCompGrid::putChomboInPetsc( const Vector<LevelData<FArrayBox> * > &a_rhs, Vec a_b )const
{
  CH_TIME("PetscCompGrid::putChomboInPetsc");
  PetscInt gid,nGrids=a_rhs.size(),*idxj,idx;
  PetscErrorCode ierr;
  PetscScalar *vals;
  PetscFunctionBeginUser;
  // iterate over real cells, get gids, make ops
  for (int ilev=nGrids-1;ilev>=0;ilev--)
    {
      const int nc=a_rhs[ilev]->nComp();
      const DisjointBoxLayout &dbl = a_rhs[ilev]->getBoxes();
      for (DataIterator dit = a_rhs[ilev]->dataIterator(); dit.ok(); ++dit)
        {
          Box region = dbl[dit];
          const BaseFab<PetscInt>& gidfab = (*m_GIDs[ilev])[dit];
          FArrayBox& phiFAB = (*a_rhs[ilev])[dit];
          idx = nc*region.numPts();
	  PetscMalloc(idx*sizeof(PetscInt), &idxj);
	  PetscMalloc(idx*sizeof(PetscScalar),&vals);
          idx = 0;
          for (BoxIterator bit(region); bit.ok(); ++bit)
            {
              const IntVect& iv = bit();
              if ( (gid=gidfab(iv,0)) >= 0)
                {
                  for (PetscInt i=nc*gid,n=0;n<nc;n++,i++)
                    {
                      idxj[idx] = i;
                      vals[idx] = phiFAB(iv,n);
                      idx++;
                    }
                }
            }
          ierr = VecSetValues(a_b,idx,idxj,vals,INSERT_VALUES);CHKERRQ(ierr);
          PetscFree(vals);
          PetscFree(idxj);
        }
    }
  ierr = VecAssemblyBegin(a_b);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(a_b);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "putPetscInChombo"
PetscErrorCode 
PetscCompGrid::putPetscInChombo( Vec a_x, Vector<LevelData<FArrayBox> * > &a_phi )const
{
  CH_TIME("PetscCompGrid::putPetscInChombo");
  PetscInt gid,nGrids=a_phi.size();
  PetscErrorCode ierr;
  const int nc=a_phi[0]->nComp(), my0eq=nc*m_gid0;
  PetscScalar *avec;
  PetscFunctionBeginUser;
  ierr = VecGetArray(a_x,&avec);CHKERRQ(ierr);
  // iterate over real cells, get gids, make ops
  for (int ilev=nGrids-1;ilev>=0;ilev--)
    {
      const DisjointBoxLayout &dbl = a_phi[ilev]->getBoxes();
      for (DataIterator dit = a_phi[ilev]->dataIterator(); dit.ok(); ++dit)
        {
          Box region = dbl[dit];
          const BaseFab<PetscInt>& gidfab = (*m_GIDs[ilev])[dit];
          FArrayBox& phiFAB = (*a_phi[ilev])[dit];
          for (BoxIterator bit(region); bit.ok(); ++bit)
            {
              const IntVect& iv = bit();
              if ( (gid=gidfab(iv,0)) >= 0)
                {
                  for (int i=nc*gid-my0eq,n=0;n<nc;n++,i++) phiFAB(iv,n) = avec[i];
                }
              else
                {
                  for (int n=0;n<nc;n++) phiFAB(iv,n) = 0.;
                }
            }
        }
      //a_phi[ilev]->exchange();
    }
  ierr = VecRestoreArray(a_x,&avec);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#include "NamespaceFooter.H"

#endif
