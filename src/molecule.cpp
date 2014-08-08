#include "molecule.h"
#include "readutil.h"

#include <fstream>

#define MAX_POLAR_DEV (1e-2)//(2e-4)//(1e-4)  //tolerance for iterative solution of F and H (S. Liu)
#define MAX_POLAR_DEV_SQR (MAX_POLAR_DEV*MAX_POLAR_DEV) 
#define MAXMONOPOLE 5
#ifdef __ACCURATE__
#define MAX_POL_ROUNDS 1000  //number of iterations (S. Liu)
#else
#define MAX_POL_ROUNDS 2
#endif

// for generating sphere points
const int NMIN = 1000;
const int NMAX = 72000;//100000;
const int NP_PER_OSCILLATION = 5;//10;

// Static Variables
double CMolecule::m_sdiel;
double CMolecule::KAPPA;
double CMolecule::SPHERETOL = 0.0;

double CMolecule::m_interRcutoff = DBL_MAX, CMolecule::m_interactRcutoff = DBL_MAX;
double CMolecule::INTRA_RCUTOFF = 30; //10000;
double CMolecule::MAXDEV = 0.0;
int CMolecule::NUM_POINTS_I,CMolecule::NUM_POINTS_X ;
bool CMolecule::m_bGrad = false;
int CMolecule::N_MOL = 0;
double CMolecule::m_total;
int CMolecule::m_unit = 1;

////////////////////// static functions ////////////////////////////
void 
CMolecule::initConstants(REAL kappa, REAL sdiel)
{
  N_MOL = 0; 
  m_sdiel = sdiel;
  KAPPA = kappa;

  CSolExpCenter::initConstants(sdiel, kappa);
  CXFormA::initConstants();

  /*
  NUM_POINTS_I = (int)(2  * _PSQ_ * NP_PER_OSCILLATION); 
  if(NUM_POINTS_I <  NMIN) NUM_POINTS_I = NMIN;
  if(NUM_POINTS_I >  NMAX) NUM_POINTS_I = NMAX;
  */
#ifdef __ACCURATE__
  NUM_POINTS_I = 250000;  //number of surface points for calculation of Imat (S. Liu)
#else
  NUM_POINTS_I = 72000;
#endif

  NUM_POINTS_X = (int)(0.5* _PSQ_ * NP_PER_OSCILLATION);
  if(NUM_POINTS_X <  NMIN) NUM_POINTS_X = NMIN;
  if(NUM_POINTS_X >  NMAX) NUM_POINTS_X = NMAX;

#ifdef __DIFF__
  NUM_POINTS_X = 1;
#endif
#ifdef __ACCURATE__
  NUM_POINTS_X = 10000;
#endif


  //#endif

  cout <<"SPHERE TOLERANCE (for extracting charges) = "<<SPHERETOL<<endl;
  cout <<" NUM_POINTS_I "<<NUM_POINTS_I<<" NUM_POINTS_X "<<NUM_POINTS_X<<endl;
  cout <<" MIN SEPDIST for using CXFormN "<<MINSEPDIST<<endl;
}

void 
CMolecule::initMutualConstants(const vector<CMolecule*> & mols, REAL interRCutoff, 
			       REAL interactRCutoff, bool bGrad)
{
  m_bGrad = bGrad;

  setInterRcutoff(interRCutoff);
  setInteractRcutoff(interactRCutoff);
  
  int nmol = mols.size(); 
  if(bGrad) 
    {  
      for(int i=0; i < nmol; i++)
	for(int ki=0; ki < mols[i]->getNKS(); ki++)
	  {
	    mols[i]->getKS(ki).initGradient(nmol); 
	  }
    }
  

}

// reset static variables concerning molecules in the system 
void 
CMolecule::resetMolSystem()
{

  deleteConstants();

  N_MOL = 0; 
  m_unit = 1;
  m_bInfinite = false;
  m_interRcutoff = 10;
  //m_mols.clear();

}

void 
CMolecule::deleteConstants()
{}

// prepare a list of spheres per molecule for gradient polarization
void 
CMolecule::generateInterMolPolList(vector<CMolecule*> & mols)
{
  int nmol = mols.size();
//  cout<<"nmol:"<<nmol<<endl;
  for(int i=0; i < nmol; i++) mols[i]->clearInterMolPolList();

  for(int i=0; i < nmol; i++)
    for(int ki=0; ki<mols[i]->getNKS(); ki++)
      {
	vector<bool> bSphereAdded(nmol, false);
	const vector<CFullSphereID> plist = mols[i]->getpKS(ki)->getInterPolList(); 
	for(int k=0; k<plist.size(); k++)
	  {
	    int j  = plist[k].mid();
	    if( bSphereAdded[j] ) continue;
	    else
	      {
		mols[j]->addInterMolPolList( CFullSphereID(i,ki) );
		bSphereAdded[j] = true;
	      }
	  }
	/*
	plist = mols[i]->getpKS(ki)->getInterOverlapPolList(); 
	for(int k=0; k<plist.size(); k++)
	  {
	    int j  = plist[k].mid();
	    if( bSphereAdded[j] ) continue;
	    else
	      {
		mols[j]->addInterMolPolList( CFullSphereID(i,ki) );
		bSphereAdded[j] = true;
	      }
	  }
	*/

      }

}

// prepare a list of spheres per molecule for gradient polarization
void 
CMolecule::generateIntraMolPolList(vector<CMolecule*> & mols )
{
  int nmol = mols.size();
  for(int j=0; j < nmol; j++) mols[j]->clearIntraMolPolLists();

  for(int j=0; j < nmol; j++)
    for(int kj=0; kj<mols[j]->getNKS(); kj++)

      // add spheres that have inter-pol neigbours
      if( !mols[j]->getKS(kj).IsEmptyInterPolList() ) 
	mols[j]->addIntraMolPolList( kj );
  // add spheres that have no inter-pol neighbors, but have inter-interact neighbors
      else if (! mols[j]->getKS(kj).IsNoInteractionList())
        mols[j]->addIntraMolInteractOnlyList( kj );

  }

bool
CMolecule::checkGradSpheres(const vector<CMolecule*> & mols, int j, int m, int km)
{

  if(m==j)
    {
      vector<int> intraMolPolList = mols[j]->getIntraMolPolList();
      for (int h=0; h < intraMolPolList.size(); h++) // intra
	{
	  int kj  = intraMolPolList[h];
	  if(kj == km) return true;
	}
    }
  else
    {
      vector<CFullSphereID> interMolPolList = mols[j]->getInterMolPolList();
      for (int h=0; h < interMolPolList.size(); h++) // inter
	{
	  int i   = interMolPolList[h].mid();
	  int ki  = interMolPolList[h].kid();
	  if(m==i && km == ki) return true; 
	}
    }
  
  return false;
  
}

void 
CMolecule::writeMolsPQR(const char * fname, const vector<CMolecule*> & mols)
{

  ofstream fout(fname);
  char buf[50];
  for(int i = 0; i < mols.size(); i++) 
    {
      char aname0[5] = " O  ";
      char aname1[5] = " N  ";
      int moltype = mols[i]->getMolType(); 
      if (moltype!=0 && moltype !=1) 
	{ cout <<"cannot handle > 2 moltype yet"<<endl; return; }
      
      CPnt rcen = mols[i]->getRCen();

      for(int ki = 0; ki < mols[i]->getNKS(); ki++)
	{

	  CPnt pos = rcen + mols[i]->getKS(ki).getCenRot();
	  fout << "ATOM  ";
	  
	  // atom index
	  fout.width(5);
	  fout.setf(ios::right, ios::adjustfield);
	  fout << ki << ' ';
	  
	  fout.setf(ios::left, ios::adjustfield);
	 
	  if(moltype == 0 ) 
	    fout << aname0 << ' ';
	  else if (moltype == 1 )
	    fout << aname1 << ' ';

	  fout << "CAV" << ' ' << " ";
	  
	  fout.setf(ios::right, ios::adjustfield);
	  fout.width(4);
	  fout << i << "    ";
	  
	  sprintf(buf, "%8.3f%8.3f%8.3f%6.2f%6.2f", pos.x(), pos.y(), 
		  pos.z(), 0.0, mols[i]->getKS(ki).getRad());
	  fout << buf << endl;;
	  
	}// endki
      
    }// endi
  fout.close();
}

void
CMolecule::saveConfig(const char * fname, const vector<CMolecule*> & mols)
{

  ofstream fout(fname);
  
  for(int i = 0; i < mols.size(); i++) 
    {
      CPnt rcen = mols[i]->getRCen();
      CQuat rot = mols[i]->getOrient(); 
      int moltype = mols[i]->getMolType();
      CPnt rI = rot.imag();
      fout <<rcen.x()<<" "<<rcen.y()<<" "<<rcen.z()<<" "<<rot.real()<<" "
	   <<rI.x()<<" "<<rI.y()<<" "<<rI.z()<<" "<<moltype<<endl;
      
      
    }
  
  fout.close();
}




////////////////////// member functions ////////////////////////////

// create a molecule from scratch with just sphere centers and radii
// use for generating imat
CMolecule::CMolecule(CPnt rcen, const vector<CPnt> &cens, const vector<double> &radii, 
		     const vector<double> &chg, const vector<CPnt> &cpos, double idiel,
		     const vector<vector<CPnt> > &SPxes,  const  vector<int> &nSPx, 
		     const vector<vector<int> >&neighs,
		     const vector< vector<int> > &intraPolLists_near, 
		     const vector<CMolCell> &molcell )
  : m_rot(false), m_p(N_POLES), m_idiel(idiel), m_bKappa(false), m_bAggregateM(false),
    m_molcells(molcell) 
{

  cout <<" =============== Molecule "<<N_MOL<<"================" <<endl;

  m_rot.reset(CQuat(), m_p);
  
  // set geometry
  m_rcen = rcen;
  vector<CPnt> cens_r = cens; 
  resetPos(cens_r, rcen);

  m_chg  = chg;
  m_cpos = cpos;
  resetPos(m_cpos, rcen);
  
  m_maxR = computeMaxRadius(cens_r, radii);  
  
  vector<int> clabel;
  assignCharges(cens_r, radii, m_cpos, clabel);

  //create solvent-exposed expansion centers
  m_nks = cens_r.size();
  m_k.resize(m_nks);

  const int N = NUM_POINTS_I;
  vector<double> rad2(m_nks);

  for(int ki = 0; ki < getNKS(); ki++)
    rad2[ki] = radii[ki]*radii[ki] *1.0;

   // setup cells
  for(int ci=0; ci<m_molcells.size(); ci++) m_cellCens.push_back(m_molcells[ci].getCen());

  int ki, tid;
#pragma omp parallel shared(rad2) private(ki, tid)
 {   
#pragma omp for 
   for(ki = 0; ki < getNKS(); ki++)
     {
       
       //      printf("================Solvent Center %d =================\n", ki);
#ifdef __OMP
       tid = omp_get_thread_num();
       //       printf("Thread %d ki = %d\n", tid, ki);
#endif
       
       vector<CPnt> posAssigned, allPosKi(m_cpos.size());
       vector<double> chgAssigned;
       vector<CPnt> SPE, SPB,SPdum, NPdum;;
       //, SPx, 
       //vector<int> neigh;
       //int nSPx;
       //findNeighbors(ki, cens_r, radii, neigh);	
       getSpherePoints(cens_r, radii, rad2,SPdum, NPdum, ki, SPE, SPB, neighs[ki]);
       cout <<"SPE "<<SPE.size()<<" SPB "<<SPB.size()<<endl;
       //getXFormSpherePoints(cens_r, radii, rad2, SPdum, NPdum, neigh, ki, SPx, nSPx, NUM_POINTS_X);
       
       // extract charge and positions       
	extractCharges(ki, cens_r[ki], m_cpos, m_chg, clabel, posAssigned, chgAssigned, allPosKi);
       
	//double start = read_timer();
	m_k[ki] = new CSolExpCenter(ki, cens_r[ki], radii[ki], SPE, SPB, 
				   SPxes[ki], nSPx[ki], neighs[ki], intraPolLists_near[ki], m_idiel, m_bKappa, m_chg,
				    allPosKi, chgAssigned, posAssigned);


     }
   
 }//endparallel

 //update static variables
 m_id = N_MOL;
 N_MOL++;
   
 cout <<"Molecule constructed."<<endl;
 cout <<" ==========================================" <<endl;
}


/////////////////////////////////////////////////////////////
// for self-polarization
//
CMolecule::CMolecule(CPnt rcen, const vector<CPnt> &cens, const vector<double> &radii, 
		     const vector<double> &chg, const vector<CPnt> &cpos, double idiel, 
		     vector<REAL*> &iMats, REAL intraRcutoff, 
		     const vector<vector<CPnt> > &SPxes,  const  vector<int> &nSPx, 
		     const vector<vector<int> >&neighs, 
		     const vector< vector<int> > &intraPolLists_near, 
		     const vector<CMolCell> &molcell )
  : m_rot(false), m_p(N_POLES), m_idiel(idiel), m_bKappa(false), m_bAggregateM(false),
    m_molcells(molcell) 
{

  cout <<" =============== Molecule "<<N_MOL<<"================" <<endl;

  assert(iMats.size() == cens.size());

  m_rot.reset(CQuat(), m_p);
  
  // set geometry
  m_rcen = rcen;
  vector<CPnt> cens_r = cens; 
  resetPos(cens_r, rcen);
  m_maxR = computeMaxRadius(cens_r, radii);

  m_chg  = chg;
  m_cpos = cpos;
  resetPos(m_cpos, rcen);
   
  vector<int> clabel;
  assignCharges(cens_r, radii, m_cpos, clabel);

  //create solvent-exposed expansion centers
  m_nks = cens_r.size();
  m_k.resize(m_nks);

  // setup cells
  for(int ci=0; ci<m_molcells.size(); ci++) m_cellCens.push_back(m_molcells[ci].getCen());

  int ki, tid;
#pragma omp parallel private(ki, tid)
 {   
#pragma omp for 
   for(ki = 0; ki < getNKS(); ki++)
     {
       
       //       printf("================Solvent Center %d =================\n", ki);
#ifdef __OMP
       tid = omp_get_thread_num();
       //       printf("Thread %d ki = %d\n", tid, ki);
#endif
       
       vector<CPnt> posAssigned, allPosKi(m_cpos.size());
       vector<double> chgAssigned;
       
       // extract charge and positions       
	extractCharges(ki, cens_r[ki], m_cpos, m_chg, clabel, posAssigned, chgAssigned, allPosKi);
       
	//double start = read_timer();
	m_k[ki] = new CSolExpCenter(ki, cens_r[ki], radii[ki], 
				    SPxes[ki], nSPx[ki], neighs[ki], intraPolLists_near[ki],
				    m_idiel, m_bKappa, m_chg,
				    allPosKi, chgAssigned, posAssigned, iMats[ki]);
     }
   
 }//endparallel

 //update static variables
 m_id = N_MOL;
 N_MOL++;
   
 cout <<"Molecule "<<N_MOL-1<<" constructed."<<endl;
 //cout <<" ==========================================" <<endl;
}

/////////////////////////////////////////////////////////////
// using pointers to imat
// read in F, H, SPx

CMolecule::CMolecule(int moltype, CPnt rcen, const vector<CPnt> &cens, const vector<double> &radii, 
		     const vector<double> &chg, const vector<CPnt> &cpos, double idiel, 
		     const vector<REAL*> &iMats, REAL intraRcutoff,
		     const vector<vector<CPnt> > &SPxes, const vector<int> &nSPx,
		     const vector<vector<int> >&neighs,
		     const vector< vector<int> > &intraPolLists_near, 
		     const vector<CMulExpan*> &Fself, const vector<CMulExpan*> &Hself,
		     const vector<CLocalExpan*> & LFs_intraSelf, const vector<CLocalExpan*> &LHs_intraSelf,
		     const vector<CLocalExpan*> & LFs_intraSelf_far, const vector<CLocalExpan*> &LHs_intraSelf_far,
 		     const vector<vector<REAL> > &qSolvedFself, 
		     const vector<vector<REAL> > &qSolvedHself, 
		     const vector<double> &totalFself, const vector<double> &totalHself, 
		     const vector<CMolCell> &molcell)
  : m_rot(false), m_p(N_POLES), m_idiel(idiel), m_bKappa(false), m_bAggregateM(false), m_moltype(moltype),
    m_molcells(molcell) 
{

  //  cout <<" =============== Molecule "<<N_MOL<<"================" <<endl;

  assert(iMats.size() == cens.size());

  m_rot.reset(CQuat(), m_p);
  
  // set geometry
  m_rcen = rcen;
  vector<CPnt> cens_r = cens; 
  resetPos(cens_r, rcen);
  m_maxR = computeMaxRadius(cens_r, radii);
  //  cout <<" mol "<<N_MOL<<" "<<rcen<<" "<<cens_r[0]<<" "<<radii[0]<<endl;

  m_chg  = chg;
  m_cpos = cpos;
  resetPos(m_cpos, rcen);
   
  vector<int> clabel;
  assignCharges(cens_r, radii, m_cpos, clabel);

  //create solvent-exposed expansion centers
  m_nks = cens_r.size();
  m_k.resize(m_nks);

  // setup cells
  for(int ci=0; ci<m_molcells.size(); ci++) m_cellCens.push_back(m_molcells[ci].getCen());
  
  int ki, tid;
#pragma omp parallel private(ki, tid)
  {   
#pragma omp for 
   for(ki = 0; ki < getNKS(); ki++)
     {       
       vector<CPnt> posAssigned, allPosKi(m_cpos.size());
       vector<double> chgAssigned;

       // extract charge and positions       
	extractCharges(ki, cens_r[ki], m_cpos, m_chg, clabel, posAssigned, chgAssigned, allPosKi);
       
	m_k[ki] = new CSolExpCenter(ki, cens_r[ki], radii[ki],
				    SPxes[ki],nSPx[ki], neighs[ki], intraPolLists_near[ki],
				    m_idiel, m_bKappa, m_chg,
				    allPosKi, chgAssigned, posAssigned, 
				    iMats[ki], Fself[ki], Hself[ki],
				    LFs_intraSelf[ki], LHs_intraSelf[ki],
				    LFs_intraSelf_far[ki], LHs_intraSelf_far[ki],
				    &(qSolvedFself[ki][0]), &(qSolvedHself[ki][0]),
				    &(totalFself[ki]), &(totalHself[ki]),
				    &m_orient, &m_rot, m_bGrad);
     }
   
 }//endparallel

 //update static variables
 m_id = N_MOL;
 N_MOL++;

}


//this function will generate lists for intra-molecular polarization
//if the distance between sphere ki and kj is less than intraRcutoff,
//they will be put into intraPolLists_near, where the local expansion
// will be performed numerically (see Eq. 27a, 28b in 2010 JCTC paper).
//Otherwise they will be  put into intraPolLists_far, where the local
//expansion will be performed analytically (S. Liu)
void
CMolecule::generateMolTypeIntraPolLists(const vector<CPnt> &cens, const vector<double> &radii,
					REAL intraRcutoff,
					vector< vector<int> > &intraPolLists_near,
					vector< vector<int> > &intraPolLists_far)
{

  const int nks = cens.size();

  intraPolLists_near.clear(); intraPolLists_near.resize(nks);
  intraPolLists_far.clear();  intraPolLists_far.resize(nks);
 
  for(int kj=0;kj<nks; kj++)
   {
     CPnt cj = cens[kj];
     REAL rj = radii[kj];

     for(int ki=0; ki<kj; ki++)
       {
	 CPnt ci = cens[ki];
	 REAL ri = radii[ki];
	 CPnt P = ci-cj; // absolute, no pbc distance for intra
	 REAL rho = P.norm();
	 bool bNear = (  rho < (ri + rj + intraRcutoff) );

	 if(bNear)
	   {
	     intraPolLists_near[ki].push_back(kj);            
	     intraPolLists_near[kj].push_back(ki);
	   } 
	 else
	   {
	     intraPolLists_far[ki].push_back(kj);            
	     intraPolLists_far[kj].push_back(ki);
	   } 
            
       }//end-ki                                                                                                    
   }//end-kj          

}

// generate a temporary molecule just to compute quantities for each moltype
// note: rcen, cens, cpos are in labframe
void
CMolecule::computeMolTypeValues(CPnt rcen, const vector<CPnt> &cens, const vector<double> &radii, 
				const vector<double> &chg, const vector<CPnt> &cpos, double idiel, 
				REAL intraRcutoff,
				const vector<vector<CPnt> > &SPxes, const vector<int> &nSPx,
				const vector<vector<int> >&neighs,
				const vector<CMulExpan*> &Fself, const vector<CMulExpan*> &Hself,
				const vector<vector<REAL> > &qSolvedFself, 
				const vector<vector<REAL> > &qSolvedHself, 
				const vector<double> &totalFself, const vector<double> &totalHself, 
				const vector< vector<int> > &intraPolLists_near, 
				const vector< vector<int> > &intraPolLists_far, 
				vector<CLocalExpan*> &LFs_intraSelf, vector<CLocalExpan*> &LHs_intraSelf,
				vector<CLocalExpan*> &LFs_intraSelf_far, vector<CLocalExpan*> &LHs_intraSelf_far)
{

  // set geometry
  vector<CPnt> cens_r = cens;
  vector<CPnt> cpos_r = cpos; 
  resetPos(cens_r, rcen);
  resetPos(cpos_r, rcen);
   
  vector<int> clabel;
  assignCharges(cens_r, radii, cpos_r, clabel);

  //create solvent-exposed expansion centers
  int nks = cens_r.size();

  // create expcenters
  vector<CSolExpCenter*> KS(nks);

#pragma omp parallel for 
  for(int ki = 0; ki < nks; ki++)
     {

       vector<CPnt> posAssigned, allPosKi(cpos_r.size());
       vector<double> chgAssigned;
       
       // extract charge and positions       
	extractCharges(ki, cens_r[ki], cpos_r, chg, clabel, posAssigned, chgAssigned, allPosKi);
		
	KS[ki] = new CSolExpCenter(ki, cens_r[ki], radii[ki],
				   SPxes[ki], nSPx[ki], neighs[ki], intraPolLists_near[ki],
				   Fself[ki], Hself[ki],
				   &(qSolvedFself[ki][0]), &(qSolvedHself[ki][0]),
				   &(totalFself[ki]), &(totalHself[ki]));
     }
  
  // initialise expansions' scale and order
  //
  LFs_intraSelf_far.resize(nks);
  LHs_intraSelf_far.resize(nks);
  LFs_intraSelf.resize(nks);
  LHs_intraSelf.resize(nks);
  
  for(int ki=0; ki<nks; ki++)
    {
      double lscale = KS[ki]->getLScale();
      LFs_intraSelf[ki] = new CLocalExpan(CRange(0, N_POLES), lscale);
      LHs_intraSelf[ki] = new CLocalExpan(CRange(0, N_POLES), lscale);
      LFs_intraSelf_far[ki] = new CLocalExpan(CRange(0, N_POLES), lscale);
      LHs_intraSelf_far[ki] = new CLocalExpan(CRange(0, N_POLES), lscale);
      
    }
  
  //
  // perform xforms to get local expansions from self values
  //
  //'xform' means re-expansion (S. Liu)
#pragma omp parallel for
  for(int ki=0; ki<nks; ki++)
    {
     CLocalExpan tLF, tLH;

     CSolExpCenter* pKi = KS[ki];
     CPnt cen_ki = pKi->getCen();
     REAL ri = pKi->getRad();

     // first compute local expansions due to far field
     vector<int> plist = intraPolLists_far[ki];

     for(int k=0; k<plist.size(); k++)
       {
	 // source 
	 const int kj = plist[k];
	 const CSolExpCenter *pKj = KS[kj];
	 const CPnt cen_kj = pKj->getCen();
	 const REAL rj = pKj->getRad();
	 
	 // generate xform                                                        
	 const CPnt P = cen_ki - cen_kj; // use absolute, not pbc distance for intra      
	 const REAL rho = P.norm();
	 const REAL sepdist = rho-ri-rj;

	 assert(sepdist  > 5);
	 // assume that far field will always be using CXFormA	 
 	 int pF = N_POLES; //CXFormBase::computeOrder(rho, pKj->getTQFself(), rj);
	 int pH = N_POLES; //CXFormBase::computeOrder(rho, pKj->getTQHself(), rj);
	 
	 CXFormAIntra XA(*pKj, *pKi);
	 XA.reset(P, pH, pF);
	 
	 XA.xformF(pKj->getFself(), tLF, true);
	 XA.xformH(pKj->getHself(), tLH, true);
	 
	 *(LFs_intraSelf_far[ki]) += tLF;
	 *(LHs_intraSelf_far[ki]) += tLH;

       }
     
     // then add contributions due to near field
     *(LFs_intraSelf[ki]) = *(LFs_intraSelf_far[ki]);
     *(LHs_intraSelf[ki]) = *(LHs_intraSelf_far[ki]);
     
     plist = intraPolLists_near[ki];
     
     for(int k=0; k<plist.size(); k++)
       {
	 // source                                                                                                      
	 const int kj = plist[k];
	 const CSolExpCenter *pKj = KS[kj];
	 const CPnt cen_kj = pKj->getCen();
	 const REAL rj = pKj->getRad();
	 
	 // generate xform                                                                                              
	 const CPnt P = cen_ki - cen_kj; // use absolute, not pbc distance for intra      
	 const REAL rho = P.norm();
	 const REAL sepdist = rho-ri-rj;
	 
	 if(sepdist > 0)
	   {
	     if( useXFormN(sepdist)  )
	       {
		 /*
		 CXFormNIntra::xformFH(*pKj, *pKi, tLF, tLH, P,
				       pKj->getQSolvedFself(), pKj->getQSolvedHself(),
				       pKj->getTQFself(), pKj->getTQHself());
		 */
		 CXFormNIntra::xformFH(*pKj, *pKi, tLF, tLH, P,
				       &(qSolvedFself[kj][0]), &(qSolvedHself[kj][0]),
				       totalFself[kj], totalHself[kj]);
		 
	       }
	     else
	       {
		 int pF = N_POLES;// CXFormBase::computeOrder(rho, pKj->getTQFself(), rj);
		 int pH = N_POLES;// CXFormBase::computeOrder(rho, pKj->getTQHself(), rj);
		 
		 CXFormAIntra XA(*pKj, *pKi);
		 XA.reset(P, pH, pF);
		 
		 XA.xformF(pKj->getFself(), tLF, true);
		 XA.xformH(pKj->getHself(), tLH, true);
	       }
	   }
	 else
	   {
	     /*
	     CXFormNIntra::xformFH(*pKj, *pKi, tLF, tLH, P,
				   pKj->getQSolvedFself(), pKj->getQSolvedHself(),
				   pKj->getTQFself(), pKj->getTQHself());
	     */
	     CXFormNIntra::xformFH(*pKj, *pKi, tLF, tLH, P,   
				   &(qSolvedFself[kj][0]), &(qSolvedHself[kj][0]),
                                   totalFself[kj], totalHself[kj]);
	     
	   }
	 
	 *(LFs_intraSelf[ki]) += tLF;
	 *(LHs_intraSelf[ki]) += tLH;

       }
     //     cout <<ki<<" LHs_intraSelf\n"<<*(LHs_intraSelf[ki])<<endl;
     //cout <<ki<<" LFs_intraSelf\n"<<*(LFs_intraSelf[ki])<<endl;
     

   }// end-ki


 // for(int ki=0; ki<nks; ki++)         delete KS[ki]; // why seg fault? later
  return;
}

/////////////////////////////////////////////////////////////
// for queries - only need H
CMolecule::CMolecule(CPnt rcen, const vector<CPnt> &cens, const vector<double> &radii, 
		     const vector<CMulExpan*> Hself, const vector<CMolCell> &molcell )
  : m_rot(false), m_p(N_POLES), m_idiel(1.0), m_bKappa(false), m_bAggregateM(false),
    m_molcells(molcell) 
{

  // set geometry
  m_rcen = rcen;
  vector<CPnt> cens_r = cens; 
  resetPos(cens_r, rcen);
  m_maxR = computeMaxRadius(cens_r, radii);

  //create solvent-exposed expansion centers
  m_nks = cens_r.size();
  m_k.resize(m_nks);

#pragma omp parallel
 {   
#pragma omp for 
   for(int ki = 0; ki < getNKS(); ki++)
     {
#ifdef __OMP
       int tid = omp_get_thread_num();
#endif
       m_k[ki] = new CSolExpCenter(ki, cens_r[ki], radii[ki], Hself[ki]);
     }
   
 }//endparallel
 
 //update static variables
 m_id = N_MOL;
 //m_mols.push_back(this);
 N_MOL++;
   
 cout <<"Molecule "<<N_MOL-1<<" constructed."<<endl;
 cout <<" ==========================================" <<endl;
}

CMolecule::~CMolecule()
{
  // delete centers
  for(int ki=0; ki< m_nks; ki++)
    delete m_k[ki]; 
}

// generate xform points for all spheres
void 
CMolecule::generateMolSPX(const vector<CPnt> &scens,const vector<double> &srad, 
			  vector<vector<CPnt> > &SPxes, vector<int> &nSPx, 
			  vector<vector<int> >&neighs)
{
  int ncen = scens.size();
  SPxes.resize(ncen);
  nSPx.resize(ncen);
  neighs.resize(ncen);

  vector<double> rad2(ncen);
  for(int ki = 0; ki < ncen; ki++) rad2[ki] = srad[ki]*srad[ki];

#pragma omp for 
   for(int ki = 0; ki < ncen; ki++)
     {
       neighs[ki].clear();
       SPxes[ki].clear();

       vector<CPnt>  SPEdum, SPBdum, SPx, SPdum, NPdum;
       int SPExSize;
       findNeighbors(ki, scens, srad, neighs[ki]);	
       getXFormSpherePoints(scens, srad, rad2, SPdum, NPdum, neighs[ki], ki, 
			    SPxes[ki], nSPx[ki], NUM_POINTS_X);
     }
}

void 
CMolecule::generateMolExposedChargesFH(const vector<CMulExpan*> Fself, const vector<CMulExpan*> Hself, 
				       const vector<CPnt> &scens,const vector<double> &srad, 
				       const vector<vector<CPnt> > &SPEx, const vector<int> & nSPx, 
				       vector<vector<double> > &qSolvedF,vector<vector<double> > &qSolvedH, 
				       vector<double> &totalF, vector<double> &totalH)
{
  int ncen = scens.size();
  qSolvedF.resize(ncen); qSolvedH.resize(ncen);
  totalF.resize(ncen); totalH.resize(ncen);

  double constH[N_POLES];
 
  for(int ki=0; ki<ncen; ki++) 
    {
      vector<double> bessel_I(N_POLES+1, 1.0);
      CMExpan::BESSEL(bessel_I, KAPPA * srad[ki]);
      for(int n=0; n<N_POLES; n++) constH[n] = CSolExpCenter::CONST2[n] / bessel_I[n]; 

      qSolvedF[ki].clear(); qSolvedH[ki].clear(); 
      CSolExpCenter::computeExposedSurfaceChargesFH_(SPEx[ki], *Fself[ki], *Hself[ki], constH, nSPx[ki], 
						     qSolvedF[ki], qSolvedH[ki], totalF[ki], totalH[ki], N_POLES);

      //      printf("%d %8.3f %8.3f %5.2f %5.2f\n" , ki, qSolvedF[ki][0], qSolvedH[ki][0], totalF[ki], totalH[ki]);

    }
}

void 
CMolecule::generateMolCells(const vector<CPnt> & scens_,const vector<double> &srad, CPnt rcen, vector<CMolCell> &molcells)
{
  
  // set scens wrt to rcen
  vector<CPnt> scens = scens_;// scen is a working copy
  for(int ki=0; ki<scens.size(); ki++) scens[ki] -= rcen; 

  // quick way to generate 8 cells
  const double maxR = CMolecule::computeMaxRadius(scens, srad);
  cout <<"Max R computed in generateMolCells = "<<maxR<<endl;
  const int nx = 2;
  const int ny = 2; 
  const int nz = 2;
  const double dx = 2*maxR/double(nx);  // cell length
  const double dy = 2*maxR/double(ny);  
  const double dz = 2*maxR/double(nz);
  double maxd = (dx > dy ? dx : dy); maxd = (maxd > dz ? maxd : dz);
  const double rc = 0.5 * maxd * sqrt(3.0); // most encompassing cell radius
  cout  << "rc = "<<rc<<endl; 
  const double xs = - 0.5*nx * dx;  
  const double ys = - 0.5*ny * dy;  
  const double zs = - 0.5*nz * dz;

  molcells.clear();
  for(int ix=0; ix<2; ix++)
    for(int iy=0; iy<2; iy++) 
      for(int iz=0; iz<2; iz++)
	{
	  double cx = xs + (ix+0.5)*dx;
	  double cy = ys + (iy+0.5)*dy;
	  double cz = zs + (iz+0.5)*dz;
	  CPnt cen = CPnt(cx, cy, cz);

	  molcells.push_back( CMolCell( cen, rc) );
	}

  for(int ki=0; ki<scens.size(); ki++) 
    {
      int ix = int( (scens[ki].x() - xs) / dx );
      int iy = int( (scens[ki].y() - ys) / dy );
      int iz = int( (scens[ki].z() - zs) / dz );
      int key = CMolCell::getKey( ix, iy, iz);
      assert(key < molcells.size() ) ;
      molcells[key].addSphere( ki );

    }


  // debug
  for(int c=0; c<8; c++)
    {
      cout <<" --- cell "<<c<<" "<<rcen + molcells[c].getCen()<<" "<<molcells[c].getSphereListSize()<<endl;
      vector<int> spherelist = molcells[c].getSphereList(); 
      // for(int h=0; h < spherelist.size(); h++ ) 
      //	cout <<" "<<spherelist[h]<<endl;
    }


}

// generate surface points for sphere i (wrt to center i)
// assigned to 'exposed' or 'buried' based on neighboring spheres
void 
CMolecule::findNeighbors(int ki, const vector<CPnt> &cens,const vector<double> &radii, 
			  vector<int> &neigh)
{
  neigh.clear();

  for (int kj = 0; kj < cens.size(); kj++)
    {
      if( ki!=kj )
	{
	  double sumRad = radii[ki]+radii[kj];
	  if( (cens[ki]-cens[kj]).normsq() < sumRad * sumRad ) 
	    {
	      neigh.push_back(kj);
	    }
	}
    }
}

// 
void 
CMolecule::getSpherePoints(const vector<CPnt> &cens,const vector<double> &radii, 
			   const vector<double> &rad2,  const vector<CPnt> &SP, const vector<CPnt> &NP, 
			   int ki, vector<CPnt> &SPE,
			    vector<CPnt> &SPB, vector<int> neigh)
{
  SPE.clear();
  SPB.clear();

  // generate sphere points
  vector<REAL> th, ph;
  spherePts(NUM_POINTS_I, th, ph);
  
  for (int n = 0; n < th.size(); n++)
    {
      CPnt q = SphToCart(CSpPnt(1, th[n], ph[n]));
      
      CPnt p = radii[ki]*q;

      // Do not choose any points which are located very close to the 
      // Z-axis (singularity issues
 
     //   if (q.x()*q.x() + q.y()*q.y() < 1e-10) continue;
      
      // Check whether point is inside / outside neighboring exposed spheres
      int k = 0;
      for (k = 0; k < neigh.size(); k++)
	{
	  int kj = neigh[k];	  
	  double distsq = (p + cens[ki] - cens[kj]).normsq();
	  
	  if ( distsq < rad2[kj])
	    {
	      SPB.push_back(p);
	      break;
	    }
	}
      
      // also check if inside boundary
      if (k == neigh.size())  
	{
	    SPE.push_back(p);
	}

    }
  
  return;

}

// itay's method
// generate surface points for sphere i (wrt to center i)
// assigned to 'exposed+interface' (SPE) or 'buried' (SPB) 
// based on neighboring spheres
void 
CMolecule::getXFormSpherePoints(const vector<CPnt> &cens, const vector<double> &radii, 
				const vector<double> &rad2,  const vector<CPnt> &SP, const vector<CPnt> &NP, 
				const vector<int> &neigh, int ki,
				vector<CPnt> &SPx, int &nSPx, int Nin) 
{
  
  SPx.clear();
  vector<CPnt> SPEx;

  int N;
  const int minimumPoints = 200;
  const double cutoff = 1.5; 
  const double maxrad = 10; 
  const double gradient = Nin / ((maxrad-cutoff)*(maxrad-cutoff));
  if (radii[ki] < cutoff ) N = minimumPoints;
  else {
    N = int ( minimumPoints + gradient * (radii[ki]-cutoff)*(radii[ki]-cutoff) );
    if (N > Nin) N = Nin;
  }
 
  // generate sphere points
  vector<REAL> th, ph;
  spherePts(N, th, ph);
  

  for (int n = 0; n < th.size(); n++)
    {
      CPnt q = SphToCart(CSpPnt(1, th[n], ph[n]));
      CPnt p = radii[ki]*q;

      // Do not choose points close to Z-axis (singularity issues)
      // if (q.x()*q.x() + q.y()*q.y() < 1e-10) continue;
      
      // Check whether point is inside / outside neighboring spheres
      bool bInside = false;
      for (int k = 0; k < neigh.size(); k++)
	{
	  int j = neigh[k];
	  double distsq = (p + cens[ki] - cens[j]).normsq();
	  
	  if( distsq < rad2[j] )
	    {
	      bInside = true;
	      break;
	    }	  
	}
      
      if (!bInside)   SPEx.push_back(p); 
	  
    }
  

  nSPx = N;
  SPx.insert(SPx.begin(), SPEx.begin(), SPEx.end() );

  return;

}

// assign charges to respective centers
void 
CMolecule::assignCharges(const vector<CPnt> &cens, const vector<double> &radii, 
			 const vector<CPnt> &cpos, vector<int> &clabel)
{


  clabel.resize(cpos.size());
  
  for(int i=0; i<cpos.size(); i++)
    {
      clabel[i] = -1; 

      for(int ki=0; ki < cens.size(); ki++)
	{
	  double dist = (cpos[i]-cens[ki]).norm() + SPHERETOL; 
	  if( dist <= radii[ki]) 
	    {
	      clabel[i] = ki;
	      break;
	    }
	}

      if(clabel[i] == -1 ) 
	{  cout <<"Error! charge "<<i<<" at "<<cpos[i]<<" not assigned to any center!"<<endl; exit(1); }

    }


  return;
}


// use rotated centers (getCenRot)
// just generate pol-lists and update bInterXForm, 
// do not generate inter-xforms, intermap - xforms will be computed on the fly
//this function will generate inter-molecular polarization list (S. Liu)
bool
CMolecule::generateInterXFormsForPolarize_LowMemory(vector<CMolecule*> & mols)
{
  const int nmol = mols.size();

  //  cout <<"Generating Interxform using Low Memory Mode"<<endl;
  if(nmol == 1)   return true;

  if (m_bInfinite)
    {
      cout<<"cannot handle infinite yet for multicenters : "<<endl;
      exit(1);
    }
 
  int i_max = (m_bInfinite ? m_unit : nmol);

  // clear away old pol lists if they exist
  for(int i=0; i<nmol; i++)
    {
      mols[i]->setInterXForm(false);
      for(int ki=0;ki<mols[i]->getNKS(); ki++) 
	  mols[i]->m_k[ki]->clearInterPolList();
    }

  
  // Assign transforms

  for(int j=1; j<nmol; j++)
    {
      // make a quick list of j's neighbors
      vector<int> neigh;
      for(int i=0; i<j; i++)
	if(getSepDist(mols[i], mols[j]) < m_interactRcutoff) 
	  neigh.push_back( i );

      if(neigh.size()==0) continue;

      // now go through all kj to list each sphere's neighbors
      int NKj = mols[j]->getNKS();
      CPnt rcen_j = mols[j]->getRCen();
      
      for(int kj=0; kj<NKj; kj++)
	{
	  CSolExpCenter* pKj = mols[j]->getpKS(kj);
	  CPnt mcen_j = pKj->getCenRot();
	  REAL rj = pKj->getRad();

	  for(int ii=0; ii<neigh.size(); ii++)
	    {
	      int i = neigh[ii];
	      int NKi = mols[i]->getNKS();
	      CPnt rcen_i = mols[i]->getRCen();
	      
	      for(int ki=0; ki<NKi; ki++)
		{
		  CSolExpCenter* pKi = mols[i]->getpKS(ki);
		  CPnt mcen_i = pKi->getCenRot();
		  REAL ri = pKi->getRad();

		  CPnt P = CSystem::pbcPos((rcen_i + mcen_i) - (rcen_j + mcen_j));
		  REAL sepdist = P.norm() -ri - rj;
		  
		  //cout <<ki<<" "<<kj<<" "<<rho<<" "<<ri<<" "<<rj<<" "<<" "<<rho-(ri+rj)<<endl;
		  // assert(rho > (ri+rj) );
		  //if(minDist > rho-(ri+rj)) {minDist = rho-(ri+rj); minKi=ki; minKj=kj;}

		  if(sepdist <=0)
		    {
		      cout <<"sepdist <=0 !"<<endl;
		      cout <<i<<" "<<ki<<" "<<j<<" "<<kj<<endl;
		      cout <<rcen_i<<" "<<rcen_j<<" "<<mcen_i<<" "<<mcen_j<<endl;
		      cout <<P.norm()<<" "<<ri<<" "<<rj<<" "<<" "<<sepdist<<endl;
		      return false;
		    }
		  //		  assert(sepdist > 0);

		  if(  sepdist < m_interactRcutoff )
		    {
		      // update interact list (reexpandDTA but don't polarize)
		      if(sepdist > m_interRcutoff)
			{
			  pKi->addInteractList(CFullSphereID(j,kj));
			  pKj->addInteractList(CFullSphereID(i,ki));
			}
		      // update polarization lists
		      else 
			//if(sepdist > 0) 
			{
			  pKi->addInterPolList(CFullSphereID(j,kj));
			  pKj->addInterPolList(CFullSphereID(i,ki));
			}
		      /*		      
		      else {
			pKi->addInterOverlapPolList(CFullSphereID(j,kj));
			pKj->addInterOverlapPolList(CFullSphereID(i,ki));
		      }
		      */

		      if(!mols[i]->getbInterXForm()) mols[i]->setInterXForm(true);
		      if(!mols[j]->getbInterXForm()) mols[j]->setInterXForm(true);
			
		    }//if-<minsep
		 
		  //else 		      cout << "no mut xforms"<<endl;
		 
		}//ki
	    }//i

	}//kj      
    }//j

  // cout <<"InterR Cutoff = "<<m_interRcutoff<<" "<<" Number of interXForms "<<nx<<endl;

  generateInterMolPolList(mols);
  generateIntraMolPolList(mols);

  return true;
}

//this function will do the first round of mutual polarization, using F and H values 
//from self-polarization to calculate new F and H multipole expansion coefficient (S. Liu)
REAL
CMolecule::recomputeFromSelfVal_LowMemory(const vector<CMolecule*> & mols, int i, int ki, bool bUpdateFarField)
{
  //generate re-expansions from self-polarization value (S. Liu)
  reexpandFromSelfVal_LowMemory(mols, i, ki,bUpdateFarField);
 
  REAL dev = m_k[ki]->solveSurfaceCharges(m_k[ki]->getFself(), m_k[ki]->getHself(), true);
#if __DEBUGDIE__
  bool bH   = m_k[ki]->getH().isBlownup() ; if(bH) cout <<"died in selfrecompute :H"<<endl;
  bool bLHS = m_k[ki]->getLHS().isBlownup(); if(bLHS) cout <<"died in selfrecompute :LHS"<<endl;
  bool bLFS = m_k[ki]->getLFS().isBlownup(); if(bLFS) cout <<"died in selfrecompute :LFS"<<endl;
  
  if( bH || bLHS || bLFS )
    {
      cout      <<m_id<<" "<<ki<<endl;
      cout <<"H\n"<<m_k[ki]->getH()<<endl;
      cout <<"LS\n"<<m_k[ki]->getLS()<<endl;      cout <<"LFN\n"<<m_k[ki]->getLFS()<<endl;
      cout <<"LHintra\n"<<m_k[ki]->getLH_intraSelf()<<endl; 
      CMolecule::writeMolsPQR("died.recomputeFromSelfVal.pqr", mols);
      CMolecule::saveConfig("died.recomputeFromSelfVal.config", mols);
      exit(1);
    }
#endif
  return dev;
}

REAL
CMolecule::recompute_LowMemory(const vector<CMolecule*> & mols, int i, int ki, bool bUpdateFarField)
{
  reexpand_LowMemory(mols, i, ki,bUpdateFarField);

  REAL dev = m_k[ki]->solveSurfaceCharges();
#if __DEBUGDIE__
  if(m_k[ki]->getH().isBlownup())
    {
      cout <<"died in recompute"<<i<<" "<<ki<<endl;
      cout <<"H\n"<<m_k[ki]->getH()<<endl;
      cout <<"LS\n"<<m_k[ki]->getLS()<<endl; 
      cout <<"LFN\n"<<m_k[ki]->getLFS()<<endl;
      cout <<"LHintra\n"<<m_k[ki]->getLH_intraSelf()<<endl;

      CMolecule::writeMolsPQR("died.recompute.pqr", mols);  
      CMolecule::saveConfig("died.recompute.config", mols);  
      exit(1);
    }
#endif
  return dev;
}

REAL
CMolecule::recompute_LowMemory(int ki, bool bUpdateFarField)
{
  
  reexpand_LowMemory(ki,bUpdateFarField);

  REAL dev;
  
  dev = m_k[ki]->solveSurfaceCharges();
 
  return dev;
}

// reexpand intra and inter centers
void
CMolecule::reexpand_LowMemory(const vector<CMolecule*> & mols, int i, int ki, bool bUpdateFarField)
{
 
  // sets LHS and LFS
  reexpandLSFromList_LowMemory(mols, i, ki);

  reexpandIntra_Near_LowMemory(ki, m_k[ki]->getLS());
  
  
  // add LHS_Far and LFS_Far
  if(bUpdateFarField)
    {
       cout <<" CMolecule::reexpand_LowMemory : Not updated"<<endl; exit(1);
      //reexpandIntra_Far_LowMemory(ki); 
      // m_k[ki]->getLHS() += m_k[ki]->getLHS_Far();  
      // m_k[ki]->getLFS() += m_k[ki]->getLFS_Far(); 
 
    }
  else
    {
      m_k[ki]->getLHS() += m_k[ki]->getLH_intraSelf_far();
      m_k[ki]->getLFS() += m_k[ki]->getLF_intraSelf_far(); 
    }
}

// reexpand intra and inter centers
void
CMolecule::reexpandFromSelfVal_LowMemory(const vector<CMolecule*> & mols, int i, int ki, bool bUpdateFarField)
{
 
  // sets LHS and LFS
  reexpandLSFromSelfVal_LowMemory(mols, i, ki);
  
  m_k[ki]->getLHS() = m_k[ki]->getLS();       
  m_k[ki]->getLHS() +=  m_k[ki]->getLH_intraSelf();

  m_k[ki]->getLFS() = m_k[ki]->getLF_intraSelf();  
}

// reexpand intra centers
void
CMolecule::reexpand_LowMemory(int ki, bool bUpdateFarField)
{
 
  reexpandIntra_Near_LowMemory(ki, CLocalExpan( CRange(),m_k[ki]->getLScale()) );
  
  // sets LHS_Far and LFS_Far
  /*
  if(bUpdateFarField)
    {
      reexpandIntra_Far_LowMemory(ki);
    }
  */

}

// reexpand m_LS (from SELECTED external centers)
// memory saving version - generate xform on the fly
void
CMolecule::reexpandLSFromList_LowMemory(const vector<CMolecule*> & mols, int i, int ki)
{

  CSolExpCenter *pKi = getpKS(ki);
  CPnt cen_ki = getRCen() + pKi->getCenRot();
  REAL ri = pKi->getRad();

  CLocalExpan LSrot( CRange(pKi->getOrder()), pKi->getLScale() );
  CLocalExpan tL;
  const vector<CFullSphereID> plist = pKi->getInterPolList(); 
  for(int k=0; k<plist.size(); k++)
    {

      // source
      int j  = plist[k].mid();
      int kj = plist[k].kid();
      const CSolExpCenter *pKj = mols[j]->getpKS(kj);

      assert (! pKj->IsEmptyInterPolList());
      CPnt cen_kj = mols[j]->getRCen() + pKj->getCenRot();
      REAL rj = pKj->getRad();

      // generate xform -only use CXFormA 
      CPnt P = CSystem::pbcPos(cen_ki - cen_kj);
      
      int pj2i = CXFormBase::computeOrder(P.norm(), pKj->getTQH(), pKj->getRad());

      CXFormA xform(*pKj, *pKi); //CXFormA is for analytical re-expansion (S. Liu)
      xform.reset(P,pj2i); 
      xform.xformH( pKj->getHrot(), tL, true);
      LSrot += tL;
    }
    
  /*
  plist = pKi->getInterOverlapPolList(); 
  for(int k=0; k<plist.size(); k++)
    {
      // source
      int j  = plist[k].mid();
      int kj = plist[k].kid();
      CXFormBase::xformHrotMono(mols[j]->getKS(kj), *pKi, tL);
      LSrot += tL;
    }
  */
  //cout<<"in reexpand: m_LS\n"<<pK->getLS()<<endl;
  m_rot.rotateWithXi(LSrot, pKi->getLS(), 1, m_p, false);

}

//generate LS from Hselfrot
void
CMolecule::reexpandLSFromSelfVal_LowMemory(const vector<CMolecule*> & mols, int i, int ki)
{

  const CSolExpCenter *pKi = getpKS(ki);
  const CPnt cen_ki = getRCen() + pKi->getCenRot();
  const REAL ri = pKi->getRad();
  CLocalExpan LSrot( CRange(pKi->getOrder()), pKi->getLScale() );
  CLocalExpan tL;

  const vector<CFullSphereID> &plist = pKi->getInterPolList(); 
  for(int k=0; k<plist.size(); k++)
    {
      // source
      const int j  = plist[k].mid();
      const int kj = plist[k].kid();
      const CSolExpCenter *pKj = mols[j]->getpKS(kj);
      const CPnt cen_kj = mols[j]->getRCen() + pKj->getCenRot();
      const REAL rj = pKj->getRad();

      // generate xform -only use CXFormA 
      const CPnt P = CSystem::pbcPos(cen_ki - cen_kj);

      // debug 
      const double sepdist = P.norm()-ri-rj;
      assert( sepdist > 0 && sepdist < m_interRcutoff);
      
      const int pj2i = CXFormBase::computeOrder(P.norm(), pKj->getTQHself(), rj);

      CXFormA xform(*pKj, *pKi);
      xform.reset(P,pj2i); 
      xform.xformH( pKj->getHselfRot(), tL, true);
#if __DEBUGDIE__ 
     if(tL.isBlownup())
	{
	  cout <<"Blown:"<<endl;
	  cout <<"rceni: "<<getRCen()<<" rcenj"<<mols[j]->getRCen()<<endl;
	  cout <<"orienti:"<<pKi->getOrient()<<endl;
	  cout <<"orientj:"<<pKj->getOrient()<<endl;
 
	  cout <<"centers:"<< cen_ki<<" "<<cen_kj<<"P= "<<P<<" sepdist:"<<sepdist<<" pj2i="<<pj2i<<endl;
	  cout <<"tL blown in reexpandLSFromSelfVal_LowMemory: "<<j<<" "<<kj<<"->"<<i<<" "<<ki<<endl;
	  cout <<"tL\n"<<tL<<endl;
	  cout <<"HselfRot_kj\n"<<pKj->getHselfRot()<<endl;
	  CMulExpan tM;
	  pKj->getRotCoeff().rotateWithXi(pKj->getHself(), tM, 1, N_POLES, true);
	  cout <<"Rotate Hself_kj\n"<<tM<<endl;
	  cout <<"Hself_kj\n"<<pKj->getHself()<<endl;
	}
#endif

      LSrot += tL;
    }
    
  /*
  plist = pKi->getInterOverlapPolList(); 
  for(int k=0; k<plist.size(); k++)
    {
      // source
      int j  = plist[k].mid();
      int kj = plist[k].kid();
      CXFormBase::xformHrotMono(mols[j]->getKS(kj), *pKi, tL);
      LSrot += tL;
    }
  */
 
  m_rot.rotateWithXi(LSrot, getKS(ki).getLS(), 1, pKi->getOrder(), false);

}

//debug only
//generate LS from Hselfrot
void
CMolecule::reexpandLSFromSelfVal_debug(const vector<CMolecule*> & mols, 
				       int i, int ki, int j, int kj)
{
  const CSolExpCenter *pKi = getpKS(ki);
  const CPnt cen_ki = getRCen() + pKi->getCenRot();
  const REAL ri = pKi->getRad();
  CLocalExpan LSrot( CRange(pKi->getOrder()), pKi->getLScale() );
  CLocalExpan tL;

  const CSolExpCenter *pKj = mols[j]->getpKS(kj);
  const CPnt cen_kj = mols[j]->getRCen() + pKj->getCenRot();
  const REAL rj = pKj->getRad();
  
  // generate xform -only use CXFormA 
  const CPnt P = CSystem::pbcPos(cen_ki - cen_kj);
  
  // debug 
  const double sepdist = P.norm()-ri-rj;
  const int pj2i = CXFormBase::computeOrder(P.norm(), pKj->getTQHself(), rj);
  /*    
  cout <<"rceni: "<<getRCen()<<" rcenj"<<mols[j]->getRCen()<<endl;
      cout <<"centers:"<< cen_ki<<" "<<cen_kj<<"P= "<<P<<" sepdist:"<<sepdist<<" pj2i="<<pj2i<<endl;
      cout <<"orienti:"<<pKi->getOrient()<<endl;
      cout <<"orientj:"<<pKj->getOrient()<<endl;
  */
      assert( sepdist > 0);
      //      assert(sepdist < m_interRcutoff);

  CXFormA xform(*pKj, *pKi);
  xform.reset(P,pj2i); 
  xform.xformH( pKj->getHselfRot(), tL, true);
  if(tL.isBlownup())
    {
      cout <<"rceni: "<<getRCen()<<" rcenj"<<mols[j]->getRCen()<<endl;
      cout <<"centers:"<< cen_ki<<" "<<cen_kj<<"P= "<<P<<" sepdist:"<<sepdist<<" pj2i="<<pj2i<<endl;
      cout <<"orienti:"<<pKi->getOrient()<<endl;
      cout <<"orientj:"<<pKj->getOrient()<<endl;
      cout <<"tL blown in reexpandLSFromSelfVal_LowMemory: "<<j<<" "<<kj<<"->"<<i<<" "<<ki<<endl;
      cout <<"tL\n"<<tL<<endl;
      cout <<"HselfRot_kj\n"<<pKj->getHselfRot()<<endl;
      CMulExpan tM;
      pKj->getRotCoeff().rotateWithXi(pKj->getHself(), tM, 1, N_POLES, true);
      cout <<"Rotate Hself_kj\n"<<tM<<endl;
      cout <<"Hself_kj\n"<<pKj->getHself()<<endl;
    }
  
  LSrot += tL;
  
  
  m_rot.rotateWithXi(LSrot, getKS(ki).getLS(), 1, pKi->getOrder(), false);
  
}


// reexpand LFS and LHS of near spheres using latest values
void
CMolecule::reexpandIntra_Near_LowMemory(int ki, const CLocalExpan &LH0)
{
  CLocalExpan tLF, tLH;

  CSolExpCenter* pKi = m_k[ki]; 
  CPnt cen_ki = pKi->getCen();
  REAL ri = pKi->getRad();
  
  // get initial values for LH and LF
  pKi->getLFS().reset(  CRange(pKi->getOrder() ) );
  pKi->getLHS().reset(  CRange(pKi->getOrder() ) );
  pKi->getLHS()+= LH0;
  
  const  vector<int> plist = pKi->getIntraPolList_near(); 
    
  for(int k=0; k<plist.size(); k++)
    {
      // source
      int kj = plist[k];
      const CSolExpCenter *pKj = m_k[kj];
      CPnt cen_kj = pKj->getCen();
      REAL rj = pKj->getRad();
      
      // generate xform
      CPnt P = cen_ki - cen_kj; // use absolute not pbc distance for intra
      REAL rho = P.norm();
      REAL sepdist = rho-ri-rj;
  
      if(sepdist > 0) 
	{	  
	  if( useXFormN(sepdist) )
	    {
          //note that here is using CXFormNIntra, which will do the re-expansion
          //numerically, rather than analytically (S. Liu)
	      CXFormNIntra::xformFH(*m_k[kj], *m_k[ki], tLF, tLH, P); 
	      pKi->getLFS() += tLF;
	      pKi->getLHS() += tLH;	  
	    }
	  else
	    {
	      int pF = CXFormBase::computeOrder(rho, pKj->getTQF(), rj);
	      int pH = CXFormBase::computeOrder(rho, pKj->getTQH(), rj);
	      CXFormAIntra XA(*m_k[kj], *m_k[ki]); 
	      XA.reset(P, pH, pF); 
	      
	      XA.xformF(m_k[kj]->getF(), tLF, true);
	      pKi->getLFS() += tLF;
	      XA.xformH(m_k[kj]->getH(), tLH, true);
	      pKi->getLHS() += tLH;
	    }
	}
      else
	{
	  CXFormNIntra::xformFH(*m_k[kj], *m_k[ki], tLF, tLH, P); 
	  pKi->getLFS() += tLF;
	  pKi->getLHS() += tLH;
	}   
    } // end kj
  
  //  cout <<"Rexpanded LH\n"<<m_k[ki]->getLH()<<endl;
  return;

}
/*
//  reexpansion from nearby other spheres in same cavity  
void 
CMolecule::reexpandIntra_Far_LowMemory(int ki )
{

  //  tA = 0.0; 
  CLocalExpan tL;

  CSolExpCenter* pKi = m_k[ki]; 

  CPnt cen_ki = pKi->getCen();
  REAL ri = pKi->getRad();

  // get initial values for LH and LF
  pKi->getLHS_Far().reset(CRange(m_k[ki]->getOrder() ) );
  pKi->getLFS_Far().reset(CRange(m_k[ki]->getOrder() ) );

  const vector<int> plist = pKi->getIntraPolList_Far(); 

  for(int k=0; k<plist.size(); k++)
    {
      // source
      int kj = plist[k];
      const CSolExpCenter *pKj = m_k[kj];
      CPnt cen_kj = pKj->getCen();
      REAL rj = pKj->getRad();

      // generate xform
      CPnt P = cen_ki - cen_kj;
      REAL rho = P.norm();

      int pF = CXFormBase::computeOrder(rho, pKj->getTQF(), rj);
      int pH = CXFormBase::computeOrder(rho, pKj->getTQH(), rj);
      //double start = read_timer();
      CXFormAIntra xform(*m_k[kj], *m_k[ki]);
      xform.reset(P, pH, pF); 
      //      double end = read_timer();time += (end-start);

      xform.xformH(m_k[kj]->getH(), tL, true);
      pKi->getLHS_Far() += tL;
   
      xform.xformF(m_k[kj]->getF(), tL, true);
      pKi->getLFS_Far() += tL;
	
      //      double end = read_timer();tA += (end-start);
    }
  
  //  cout <<"Rexpanded LH\n"<<m_k[ki]->getLH()<<endl;
  return;

}
*/

// Precompute the sum of the products of dT(i,j)*A(i) for all molecules ( i = 0 -> N_MOL)
void
CMolecule::prepareDTA_LowMemory(const vector<CMolecule*> & mols, int j, vector<CGradExpan> &tG)
{

  CGradExpan tG2;
  
  vector<CFullSphereID> interMolPolList = mols[j]->getInterMolPolList();
  vector<int> intraMolPolList = mols[j]->getIntraMolPolList();
  int j_intersize = interMolPolList.size();
  int j_intrasize = intraMolPolList.size(); 
  tG.resize( j_intersize + j_intrasize);
  
  int k;
  // first do external spheres to mol_j
  for (k = 0; k <j_intersize; k++)
    {
      int i   = interMolPolList[k].mid();
      int ki  = interMolPolList[k].kid();

      CSolExpCenter* pKi = mols[i]->getpKS(ki);
      CPnt cen_ki = mols[i]->getRCen() + pKi->getCenRot();
      
      tG[k].reset(    pKi->getOrder()  );
      tG[k].setScale( pKi->getLScale() );
      const vector<CFullSphereID> plist = pKi->getInterPolList(); 
      for(int h=0; h<plist.size(); h++)
	{
	  int m  = plist[h].mid();
	  
	  if (j != m)  continue; // only consider when m = j
	  
	  // generate xform
	  int km = plist[h].kid();
	  
	  //assert( checkGradSpheres(mols, j, m, km));

	  const CSolExpCenter *pKm = mols[m]->getpKS(km);
	  CPnt cen_km = mols[m]->getRCen() + pKm->getCenRot();
	  CPnt P = cen_ki - cen_km;
	  int pm2i = CXFormBase::computeOrder(P.norm(), pKm->getTQH(), pKm->getRad());
	  CXFormA XA(*pKm, *pKi, true);
	  XA.reset(P,pm2i); 
	  
	  XA.xformH(mols[m]->getKS(km).getHrot(), tG2, true);
	  tG[k] += tG2;
	  
	}// end h
      /*  
      plist = pKi->getInterOverlapPolList(); 
      for(int h=0; h<plist.size(); h++)
	{
	  int m  = plist[h].mid();

	  // only consider when m = j
	  if (m != j) continue; 
	  
	  int km = plist[h].kid();
	  CXFormBase::xformHrotMonoG(mols[m]->getKS(km), *pKi, tG2); 
	  tG[k] -= tG2;
	  
	}// end h-overlap
      */      
    }//end k-inter

  //==========================================================
  // then consider spheres in mol_j; i.e. i == j => del_(j,kj)T(i,ki <- m,km) = del_(i,ki)T(i,ki <- m,km)
  for (int hh = 0; hh <j_intrasize; hh++, k++)
    {
      int ki  = intraMolPolList[hh];

      CSolExpCenter* pKi = mols[j]->getpKS(ki);
      CPnt cen_ki = mols[j]->getRCen() + pKi->getCenRot();
      
      tG[k].reset(    pKi->getOrder()  );
      tG[k].setScale( pKi->getLScale() );
      const vector<CFullSphereID> plist = pKi->getInterPolList(); 
      for(int h=0; h<plist.size(); h++)
	{
	  int m  = plist[h].mid();
	  int km = plist[h].kid();
	  assert( checkGradSpheres(mols, j, m, km));
	  
	  // generate xform
	  const CSolExpCenter *pKm = mols[m]->getpKS(km);
	  CPnt cen_km = mols[m]->getRCen() + pKm->getCenRot();
	  CPnt P = cen_ki - cen_km;
	  int pm2i = CXFormBase::computeOrder(P.norm(), pKm->getTQH(), pKm->getRad());
	  CXFormA XA(*pKm, *pKi, true);
	  XA.reset(P,pm2i); 
	  XA.xformH(mols[m]->getKS(km).getHrot(), tG2, true);
	  tG[k] += tG2;
	  
	}// end-h

      const vector<CFullSphereID> irlist = pKi->getInteractList();
      for(int h=0; h<irlist.size(); h++)
	{
	  int m  = irlist[h].mid();
	  int km = irlist[h].kid();
	  // generate xform                                                                                                 
	  const CSolExpCenter *pKm = mols[m]->getpKS(km);
	  CPnt cen_km = mols[m]->getRCen() + pKm->getCenRot();
	  CPnt P = cen_ki - cen_km;
	  int pm2i = CXFormBase::computeOrder(P.norm(), pKm->getTQH(), pKm->getRad());
	  CXFormA XA(*pKm, *pKi, true);
	  XA.reset(P,pm2i);
	  XA.xformH(mols[m]->getKS(km).getHrot(), tG2, true); // use self val                                           
	  tG[k] += tG2;

	}// end-h                                                                                                           



      /*
	plist = pKi->getInterOverlapPolList(); 
	for(int h=0; h<plist.size(); h++)
	  {
	    int m  = plist[h].mid();
	    int km = plist[h].kid();
	    
	    CXFormBase::xformHrotMonoG(mols[m]->getKS(km), *pKi, tG2); 
	    tG[k] += tG2;
	    
	  }// end-h-overlap
      */
    } // end hh-intra
  return;
}

// Precompute the sum of the products of dT(i,j)*A(i) for all molecules ( i = 0 -> N_MOL)
// spheres in mol_j; i.e. i == j => del_(j,kj)T(i,ki <- m,km) = del_(i,ki)T(i,ki <- m,km)
// 1) only for i=j case
// 2) using selfH values
// 3) include interactlist because we need it for force
void
CMolecule::prepareDTA_iself_LowMemory(const vector<CMolecule*> & mols, int i, int ki)
{
  
  CGradExpan tG2, tG;

  CSolExpCenter* pKi = mols[i]->getpKS(ki);
  CPnt cen_ki = mols[i]->getRCen() + pKi->getCenRot();
  tG.reset(    pKi->getOrder()  );
  tG.setScale( pKi->getLScale() );
  
  // first reexpand H from ki's interpollist
  const vector<CFullSphereID> plist = pKi->getInterPolList();
  for(int h=0; h<plist.size(); h++)
    {
      int m  = plist[h].mid();
      int km = plist[h].kid();
      
      // generate xform                                                                                           
      const CSolExpCenter *pKm = mols[m]->getpKS(km);
      CPnt cen_km = mols[m]->getRCen() + pKm->getCenRot();
      CPnt P = cen_ki - cen_km;
      int pm2i = CXFormBase::computeOrder(P.norm(), pKm->getTQHself(), pKm->getRad());
      CXFormA XA(*pKm, *pKi, true);
      XA.reset(P,pm2i);
      XA.xformH(mols[m]->getKS(km).getHselfRot(), tG2, true); // use self val
      tG += tG2;
      
    }// end-h   
  
  // next reexpand H from ki's interactlist
  const vector<CFullSphereID> irlist = pKi->getInteractList();
  for(int h=0; h<irlist.size(); h++)
    {
      int m  = irlist[h].mid();
      int km = irlist[h].kid();
      // generate xform      
      const CSolExpCenter *pKm = mols[m]->getpKS(km);
      CPnt cen_km = mols[m]->getRCen() + pKm->getCenRot();
      CPnt P = cen_ki - cen_km;
      int pm2i = CXFormBase::computeOrder(P.norm(), pKm->getTQHself(), pKm->getRad());
      CXFormA XA(*pKm, *pKi, true);
      XA.reset(P,pm2i);
      XA.xformH(mols[m]->getKS(km).getHselfRot(), tG2, true); // use self val
      tG += tG2;
      
    }// end-h   
  
  pKi->getgLHN() = tG;

  return;
}

// debug
// Precompute the sum of the products of dT(i,j)*A(i) for all molecules ( i = 0 -> N_MOL)
void
CMolecule::prepareDTAintra_LowMemory(const vector<CMolecule*> & mols, int j, vector<CGradExpan> &tG)
{
  
  CGradExpan tG2;
  vector<int> intraMolPolList = mols[j]->getIntraMolPolList();
  int j_intrasize = intraMolPolList.size();
  tG.resize( j_intrasize);

  // then consider spheres in mol_j; i.e. i == j => del_(j,kj)T(i,ki <- m,km) = del_(i,ki)T(i,ki <- m,km)               
  for (int k = 0; k <j_intrasize; k++)
    {
      int ki  = intraMolPolList[k];

      CSolExpCenter* pKi = mols[j]->getpKS(ki);
      CPnt cen_ki = mols[j]->getRCen() + pKi->getCenRot();

      tG[k].reset(    pKi->getOrder()  );
      tG[k].setScale( pKi->getLScale() );
      
      // first reexpand H from ki's interpollist
      const vector<CFullSphereID> plist = pKi->getInterPolList();
      for(int h=0; h<plist.size(); h++)
          {
            int m  = plist[h].mid();
            int km = plist[h].kid();
            assert( checkGradSpheres(mols, j, m, km));

            // generate xform                                                                                           
            const CSolExpCenter *pKm = mols[m]->getpKS(km);
            CPnt cen_km = mols[m]->getRCen() + pKm->getCenRot();
            CPnt P = cen_ki - cen_km;
            int pm2i = CXFormBase::computeOrder(P.norm(), pKm->getTQHself(), pKm->getRad());
            CXFormA XA(*pKm, *pKi, true);
            XA.reset(P,pm2i);
	    XA.xformH(mols[m]->getKS(km).getHselfRot(), tG2, true); // use self val
            tG[k] += tG2;

          }// end-h   
                                                                     
    } // end k-intra                                          
  return;
}


//this function will do the re-expansion of F and H's gradients (S. Liu)
void 
CMolecule::reexpandGrad_LowMemory(const vector<CMolecule*> & mols, 
				  const CGradExpan &tG_DTA, CGradExpan &tG, 
				  int i, int ki, int j, bool bUpdateFarField, 
				  const vector<CGradExpan> &tGHrots,  
				  map<CFullSphereID, int, classCompareCFullSphereID> &gHmap)
{
  CSolExpCenter* pKi = mols[i]->getpKS(ki);
  CPnt cen_ki = mols[i]->getRCen() + pKi->getCenRot();
  CGradExpan tGrot = tG_DTA;

  vector<CFullSphereID> plist;
  CGradExpan tG2;

  // if i=j, ki's external spheres to polarize grad = ki's interpollist
  // otherwise, we consider ki's interpolist, but only include spheres from mol j
  // i.e. we ignore spheres in ki's interpolist that do not belong to j
  if(i==j) plist = pKi->getInterPolList(); 
  else 
    {
      const vector<CFullSphereID> interPolList = pKi->getInterPolList(); 
      for(int h=0; h < interPolList.size(); h++) 
	if( interPolList[h].mid() == j ) 
	  plist.push_back( interPolList[h] );
    }
  for(int h=0; h<plist.size(); h++)
    {
      int m  = plist[h].mid();
      int km = plist[h].kid();

      // generate xform
      const CSolExpCenter *pKm = mols[m]->getpKS(km);
      CPnt cen_km = mols[m]->getRCen() + pKm->getCenRot();
      CPnt P =CSystem::pbcPos(cen_ki - cen_km);
      int pm2i = CXFormBase::computeOrder(P.norm(), pKm->getTGH(), pKm->getRad());
      CXFormA XA(*pKm, *pKi, false);
      XA.reset(P,pm2i); 
      XA.xformH(tGHrots[ gHmap[CFullSphereID(m,km)] ], tG2, true);

      tGrot += tG2;
    }
  /*
  //==========================================================
  // now do overlap list
  if(i==j) plist = pKi->getInterOverlapPolList(); 
  else 
    {
      plist.clear();
      vector<CFullSphereID> interPolList = pKi->getInterOverlapPolList(); 
      for(int h=0; h < interPolList.size(); h++) 
	if( interPolList[h].mid() == j ) 
	  plist.push_back( interPolList[h] );
    }
  for(int h=0; h<plist.size(); h++)
    {
      int m  = plist[h].mid();
      int km = plist[h].kid();

      CXFormBase::xformGHrotMono(mols[m]->getKS(km), *pKi, tG2);
      tGrot += tG2;
    }
  
  //==========================================================
  */
  // if(bUpdateFarField) { do far field updates ... tG += ...}

  // rotate final tG to sphere ki's original frame
   m_rot.rotateWithXi(tGrot, tG, 1, m_p, false);

  return;
}

void 
CMolecule::reexpandIntraGrad_LowMemory(const vector<CMolecule*> & mols, 
				   CGradExpan &tGF, CGradExpan &tGH, 
				       int i, int ki, int j,  bool bUpdateFarField, 
				       const vector<CGradExpan> &tGFs,
				       const vector<CGradExpan> &tGHs,
				       map<CFullSphereID, int, classCompareCFullSphereID> &gHmap) 
{

  CSolExpCenter* pKi = m_k[ki];
  CPnt cen_ki = pKi->getCen();
  REAL ri = pKi->getRad();
   
  REAL scale = pKi->getLScale();
  int p = m_k[ki]->getOrder();
  tGF.setScale( scale ); tGF.reset( p );
  tGH.setScale( scale ); tGH.reset( p );

  CGradExpan tG2, tLGF, tLGH;
  vector<int> plist;
  if(i==j) plist = mols[j]->getIntraMolPolList();
  else 
    {
      vector<CFullSphereID> interMolPolList = mols[j]->getInterMolPolList();
      for(int h=0; h < interMolPolList.size(); h++) 
	if( interMolPolList[h].mid() == i) 
	  plist.push_back( interMolPolList[h].kid() );
    }
  
  for(int h=0; h<plist.size(); h++)
    {
      int km = plist[h];
      if (km == ki) continue;
      assert ( checkGradSpheres(mols, j, i, km) );

      const CSolExpCenter *pKm = m_k[km];
      
      // generate xform
      CPnt cen_km = pKm->getCen();
      REAL rm = pKm->getRad();
      CPnt P = cen_ki - cen_km;
      REAL rho = P.norm();
      REAL sepdist = rho-ri-rm;
      
      if(sepdist > 0) 
	{
	  
	  if( useXFormN(sepdist)  )
	    {
	      CXFormNIntra::xformGFH(*m_k[km], *m_k[ki], tLGF, tLGH, P); 
	      tGF += tLGF;
	      tGH += tLGH;
	    }
	  else
	    {
	      int pF = CXFormBase::computeOrder(rho, pKm->getTGF(), rm);
	      int pH = CXFormBase::computeOrder(rho, pKm->getTGH(), rm);
	      CXFormAIntra XA(*m_k[km], *m_k[ki]); 
	      XA.reset(P, pH, pF); 
	      int key = gHmap[CFullSphereID(m_id,km)];
	      XA.xformF(tGFs[key], tG2, true);
	      tGF += tG2;
	      XA.xformH(tGHs[key], tG2, true);
	      tGH += tG2;	      
	    }
	}

      else // overlapping neighbors - use cxformN
	{
	  CXFormNIntra::xformGFH(*m_k[km], *m_k[ki], tLGF, tLGH, P); 
	  tGF += tLGF;
	  tGH += tLGH;
	}   
	  
    }



  // if(bUpdateFarField) { do far field updates ... tG += ...}

  return;
}

//this function will calculate gradients of F and H iteratively (S. Liu)
double
CMolecule::recomputeGrad_LowMemory(const vector<CMolecule*> & mols, const CGradExpan &tG_DTA, 
				   int i, int ki, int j, bool bUpdateFarField,  vector<CGradExpan> &tGFs,
				    vector<CGradExpan> &tGHs, vector<CGradExpan> &tGHrots, 
				   map<CFullSphereID, int, classCompareCFullSphereID> &gHmap)
{
  CGradExpan LGF, LGH, LGH2;
  REAL dev;

  reexpandGrad_LowMemory(mols, tG_DTA, LGH, i, ki, j, bUpdateFarField, tGHrots, gHmap);
  
  if(j==i) m_k[ki]->getgLHN() = LGH;

  
  reexpandIntraGrad_LowMemory(mols, LGF, LGH2, i, ki, j, bUpdateFarField, tGFs, tGHs, gHmap);
  
  LGH += LGH2;
  
  int key = gHmap[CFullSphereID(m_id,ki)];
  assert(key < tGFs.size()); 
  assert(key < tGHs.size()); 
  assert(key < tGHrots.size()); 
  dev   = m_k[ki]->solveSurfaceGradient(j, tGFs[key], tGHs[key], tGHrots[key], LGF, LGH);
 
#if __DEBUGDIE__
  CGradExpan gh = tGHs[key];
  if( fabs(gh[0][0]) >= MAXMONOPOLE || fabs(gh[1][0]) >= MAXMONOPOLE  || fabs(gh[2][0]) >= MAXMONOPOLE 
      || isnan(gh[0][0]) || isnan(gh[1][0]) || isnan(gh[2][0]))
    
    {
      cout <<"died in recomputegrad"<<endl;
      CMolecule::writeMolsPQR("died.recomputeGrad.pqr", mols);
      CMolecule::saveConfig("died.recomputeGrad.config", mols);
      exit(1);
    }
#endif
  return dev; 
  
}

//this function will extract local expansions for interaction energy calculation 
//purpose (S. Liu)
REAL 
CMolecule::interactCenters_LowMemory(CMolecule* moli, CMolecule* molj)
{
  
  REAL pot = 0.0;
  int dummyP = 1;

  int i = moli->getID();
  int j = molj->getID();
  CPnt rcen_i = moli->getRCen();

  for (int ki=0; ki < moli->getNKS(); ki++)
    {
      const CSolExpCenter *pKi = moli->getpKS(ki);
      int nQki = pKi->getSPExSize();
      CPnt mcen_i = pKi->getCen();
	  
      CLocalExpan Lki( CRange(0,1), pKi->getLScale() ); 
	  
      for (int kj=0; kj < molj->getNKS(); kj++)
	{
	  const CSolExpCenter *pKj = molj->getpKS(kj);
	  CLocalExpan tL;
	  
	  CPnt P = (rcen_i + mcen_i) - (molj->getRCen() + pKj->getCen() );
	  REAL rho = P.norm();
	  int pi2j = CXFormBase::computeOrder(rho, pKi->getTQH(), pKi->getRad());
	  int pj2i = CXFormBase::computeOrder(rho, pKj->getTQH(), pKj->getRad());
	  int minp = (pi2j < pj2i ? pi2j : pj2i);
	  bool bJ2I; 

	  CXFormA XA(*pKj,*pKi);
	  XA.reset(P, minp); 
	  bJ2I = minp < pi2j;
	  	  
//	  if( bJ2I ) // add to ki's local expansion if it's cheaper (transform pKj) 
//	    {
             // for debug purpose, S. Liu
            // cout << "bFor is true"<<endl;
	      XA.xformH( pKj->getH(), tL, true);
	      Lki += tL;
//	    }
	  
//	  else // otherwise, transform pKi 
//	    {
            //  cout << "bFor is false"<<endl;
//	      XA.xformH( pKi->getH(), tL, false);
       
//	      pot += inprod(pKj->getH(), tL);
//	    }
	  
	}// end-kj
	  
      pot += inprod(pKi->getH(), Lki);
      
    }//end-ki
    
  return pot;
}

// Iteratively solve for MP of all molecules until converged
//this function will perform self-polarization (S. Liu)
void
CMolecule::polarize_self(bool bPot, int farFieldFreq)
{

  cout<<"self polarize with farfield freq = "<< farFieldFreq<<endl;

  double  dev = 0.0;   
  int ct;
  const int MAX_POL_ROUNDS_SELF = 200;

  m_total = m_nks;
  REAL itot = 1.0/m_total;  
  cout <<"m_total = "<<m_total<<endl;

  if(m_nks == 1) 
    {
      dev = recompute_LowMemory(0, true);
      ct = 1;
      while (dev*itot > MAX_POLAR_DEV && ct < MAX_POL_ROUNDS_SELF) 
	{
	  dev = recompute_LowMemory(0, true);
	  ct++;
	}
    }
  
  else 
    {      
      double start = read_timer();
      
      int NKis = getNKS();
      
      int ki, tid;
      REAL maxdev = 0.0, maxdev_new;

#pragma omp parallel for reduction(+:dev) shared(NKis)
      for(ki=0; ki < NKis; ki++) 
	{
	  REAL dev_ki = 0; 

	  dev_ki = recompute_LowMemory(ki, true); // update far field for first round
	  dev += (getKS(ki).getRad()*getKS(ki).getRad()) * dev_ki;
	  if(maxdev < dev_ki) maxdev = dev_ki;
	}

      ct = 1; // no. of molecule iterations
      printf("After 1st round : ct %d maxdev %e total dev: %e avg %e\n", ct, maxdev, dev, itot*dev);
            
      while ( maxdev > MAX_POLAR_DEV_SQR && ct < MAX_POL_ROUNDS_SELF) 
	
	{	  
	  dev = 0.0; 
	  maxdev_new = 0.0;
	  int count = 0;

#pragma omp parallel for reduction(+:dev) shared(NKis)
	  for(ki=0; ki < NKis; ki++)    
	    {
	      REAL dev_ki;
	      dev_ki = getKS(ki).getDev() ; 
	      if( dev_ki > 0.1*maxdev || ct % 5 == 0)
		{
		  dev_ki = recompute_LowMemory(ki, true);
		  count ++;
		}
	      
	      if(maxdev_new < dev_ki ) maxdev_new = dev_ki;
	      dev += (getKS(ki).getRad()*getKS(ki).getRad()) * dev_ki;

	    } // endki

	  printf("%d) count %d maxdev: %e dev: %e avg: %e\n",ct, count, maxdev_new, dev, dev*itot);
	  maxdev = maxdev_new;
	  
	  ct++;
	  
	  if (ct == MAX_POL_ROUNDS_SELF)
	    {
	      cout << "Polarization does not converge!!! dev=" 
		   << dev << " avg " << itot*dev<<" maxdev "<<maxdev <<" ct "<< ct << endl;
	      // exit(0);
	    }
	  
	  double vm, rss;
	  
	}//endwhile
      
      
      double end = read_timer();
      cout <<"Time taken to polarize [s]: "<<end-start<<
	", ct: "<< ct<< ", per cycle: "<<(end-start)/ double(ct)<<endl;  
      printf("Potential converged : ct %d total dev: %e avg %e\n", ct, dev, itot*dev);
      cout <<"avg moldev = "<<maxdev<<endl;      
      
    } // endif_nscenter


  if (bPot)
    return;


}

// Iteratively solve for MP of all molecules until converged
// only update farfield infrequently 
//this function will do mutual polarization (S. Liu)
void
CMolecule::polarize_mutual(vector<CMolecule*> & mols, bool bPot, int farFieldFreq)
{

//   cout <<"Polarize Mutual (Parallel)"<<endl;

  int nmol = mols.size();
  const int i_max = (m_bInfinite ? m_unit : nmol);
//  cout <<"i_max="<<i_max<<endl; 
  //////////// ONE MOLECULE ////////////
  if(i_max == 1 )     return;
  
  //////////// MULTIPLE  MOLECULES ////////////
  double  dev = 0.0;
  int ct, mtotal=0;
  REAL itot;
  vector<double> maxdevs(i_max, 0.0);
  
  double start = read_timer();

#pragma omp parallel
 {

   // if centers will be polarized: update Hrot using self polarized values
#pragma omp for schedule(dynamic)
 for(int i = 0; i < i_max; i++)
   {
     if( ! mols[i]->getbInterXForm() ) continue;

     mols[i]->rotateRotCoeff();
     for(int ki=0; ki < mols[i]->getNKS(); ki++)
       {
	 if( !mols[i]->getKS(ki).IsNoInteractionList() ) 
	   mols[i]->getKS(ki).rotateHself();
       }
   }
 
 // first round; use self values
#pragma omp for reduction(+:dev) reduction(+:mtotal)  schedule(dynamic)
 for (int i = 0; i < i_max; i++)
   {
     REAL dmol_i = 0.0;
     
     if( mols[i]->getbInterXForm() )
       for(int ki=0; ki < mols[i]->getNKS(); ki++){
	 
	 if( mols[i]->getKS(ki).IsEmptyInterPolList() ) continue;
	 
	 REAL dev_ki = mols[i]->recomputeFromSelfVal_LowMemory(mols, i, ki, false); 

	 dmol_i += dev_ki;
	 if(maxdevs[i] < dev_ki) maxdevs[i]  = dev_ki;
	 mtotal++;
	 
	 
       }// end if-bxform, ki loop      
     
     dev += dmol_i;
     
   } //end i-loop
 
 }//endparallel
 
 MAXDEV = 0.0;
 for(int i=0; i<i_max; i++) if (maxdevs[i] > MAXDEV) MAXDEV = maxdevs[i]; // record max deviation
 
 ct = 1; // no. of molecule iterations
 double e1 = read_timer();
 itot = (mtotal > 0 ? 1.0/mtotal : 0);
 while ( dev*itot > MAX_POLAR_DEV_SQR && ct < MAX_POL_ROUNDS) 
   {             
     dev = 0.0;
     maxdevs.assign(i_max, 0.0);
     int count=0;
     
#pragma omp parallel for reduction(+:dev) schedule(dynamic)
     for (int i = 0; i < i_max; i++)
       {
	 REAL dmol_i = 0.0;
	 
	 if( mols[i]->getbInterXForm() )
	   for(int ki=0; ki < mols[i]->getNKS(); ki++) {

	     if( mols[i]->getKS(ki).IsEmptyInterPolList() ) continue;
	     
	     REAL dev_ki = mols[i]->getKS(ki).getDev();

	     if(  dev_ki > 0.1*MAXDEV)
	       {
		 dev_ki = mols[i]->recompute_LowMemory(mols, i, ki, false);
                 count++;

	       }// endif-dev
	      
	      if(maxdevs[i] < dev_ki) maxdevs[i]  = dev_ki;
	      dmol_i += dev_ki; 

	    }//end ki loop

	  dev += dmol_i;	      
	  //printf("molecule %d : dmol_i: %e\n", i, dmol_i);
	  
	} // end-i

      MAXDEV=0.0;
      for(int i=0; i<i_max; i++) if (maxdevs[i] > MAXDEV) MAXDEV = maxdevs[i];

      //Incrementing  indices	
      ct ++ ;
    } // end-while
  
 
  double end = read_timer();
 if (bPot) return;
 
  //--------------------------------------------------------------------
  // Solving for gradient (dH_j)
 double startg = read_timer();
 
#pragma omp parallel for schedule(dynamic)
  for (int j = 0; j < i_max; j++) 
    {

     //for debug (S. Liu)
//       cout << "molecule index"<<j<<endl;
      double startj = read_timer();

      vector<CFullSphereID> interMolPolList = mols[j]->getInterMolPolList();
      int j_intersize = interMolPolList.size();

      // if molecule j has no inter polarization neighbors, 
      // compute DTA for j's spheres that have interaction neighbor
      if( j_intersize == 0)
	{
	  if(mols[j]->getbInterXForm())
	    for(int kj=0; kj<mols[j]->getNKS(); kj++)
	      {
		if(! mols[j]->getKS(kj).IsNoInteractionList()) 
		  prepareDTA_iself_LowMemory(mols, j, kj);
	      }	    
	}
     
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // start polarizing grad if there are ext. spheres nearby
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      else
	{
      //for debug (S. Liu)
//      cout<< "there are ext"<<endl;
	  
	  vector<int> intraMolPolList = mols[j]->getIntraMolPolList();
	  int j_intrasize = intraMolPolList.size();
	  int j_size = j_intersize + j_intrasize;
	  vector<REAL> d(j_size);
	  itot =1.0/double(j_size);
	  
	  //      prepare DTA[i,ki] = sum_m,km( del_j(T_ki,km) * H_km )
#ifdef __NOGRADPOL__
	  vector<CGradExpan> tGe;
	  prepareDTAintra_LowMemory(mols, j, tGe); // debug!
#else
	  vector<CGradExpan> tG;
	  prepareDTA_LowMemory(mols, j, tG);
#endif
	  
	  
	  // polarize inter and intra(j) spheres wrt molecule j
	  
	  
	  //===============================================
	  vector<CGradExpan> tGHs( j_size ), tGFs( j_size ), tGHrots( j_size );  
	  map<CFullSphereID, int, classCompareCFullSphereID> gHmap;
	  
      //for debug (S. Liu)
//      cout<< "polarize inter and intra"<<endl;
	  int k = 0;
	  // reset gradients that will be used
	  for (k=0; k < j_intersize; k++) // inter
	    {
	      const int i   = interMolPolList[k].mid();
	      const int ki  = interMolPolList[k].kid();
	      const int p = mols[i]->getKS(ki).getOrder();
	      const double mscale = mols[i]->getKS(ki).getMScale();
	      
	      tGFs[k].reset( p );
	      tGFs[k].setScale( mscale );
	      tGHs[k].reset( p );
	      tGHs[k].setScale( mscale );
	      tGHrots[k].reset( p );
	      tGHrots[k].setScale( mscale );
	      
	      gHmap[CFullSphereID(i,ki)] = k;
	      
	    }
      //for debug (S. Liu)
//      cout<< "point 1"<<endl;
	  for (int h=0; h < j_intrasize; h++,k++) // intra
	    {
	      const int kj  = intraMolPolList[h];
	      const int p = mols[j]->getKS(kj).getOrder();
	      const double mscale = mols[j]->getKS(kj).getMScale();
	      tGFs[k].reset( p );
	      tGFs[k].setScale( mscale );
	      tGHs[k].reset( p );
	      tGHs[k].setScale( mscale );
	      tGHrots[k].reset( p );
	      tGHrots[k].setScale( mscale );
	      gHmap[CFullSphereID(j,kj)] = k;
	    }
      //for debug (S. Liu)
//      cout<< "point 2"<<endl;
	  
	  REAL dev = 0.0;
	  
#ifndef __NOGRADPOL__
	  for (k=0; k < j_intersize; k++) // inter
	    {
	      int i   = interMolPolList[k].mid();
	      int ki  = interMolPolList[k].kid();
	      d[k] = mols[i]->recomputeGrad_LowMemory(mols, tG[k], i, ki, j, false, tGFs, tGHs, tGHrots, gHmap);
	      dev += d[k];
	    }
#endif
	  
	  for (int h=0; h < j_intrasize; h++, k++) // intra
	    {
	      
	      int kj  = intraMolPolList[h];
	      
#ifdef __NOGRADPOL__
	      mols[j]->getKS(kj).getgLHN() = tGe[h];// debug bypass grad polarization
#else
	      d[k] = mols[j]->recomputeGrad_LowMemory(mols, tG[k], j, kj, j, false, tGFs, tGHs, tGHrots, gHmap);
	      dev += d[k];
#endif
	    }
	  
	  //      printf("%d) dev: %e avg: %e\n",ct,  dev, dev*itot);
	  
#ifndef __NOGRADPOL__      
	  ct = 1;
	  k = 0;
	  int i   = interMolPolList[k].mid();
	  int ki  = interMolPolList[k].kid();
      //for debug (S. Liu)
//      cout<< "point 3"<<endl;
	  
	  
	  while (dev*itot > MAX_POLAR_DEV_SQR && ct < MAX_POL_ROUNDS)
	    {
      //for debug (S. Liu)
//      cout<< "ct ="<<ct<<endl;
	      dev -= d[k];
	      d[k] = mols[i]->recomputeGrad_LowMemory(mols, tG[k], i, ki, j, false,tGFs, tGHs, tGHrots, gHmap);
	      dev += d[k];
	      
	      //Incrementing indices
	      k++; 
	      if(k==j_size) 
		{
		  k = 0;
		  //      printf("%d) dev: %e avg: %e\n",ct,  dev, dev*itot);
		  ct++;
		  
		}
	      
	      if(k < j_intersize) // read from intermollist
		{
		  i  = interMolPolList[k].mid();
		  ki  = interMolPolList[k].kid();
		}
	      else // read from intramollist
		{
		  i  = j;
		  ki  = intraMolPolList[k-j_intersize];
		}
	      
	    }//end while
	  // save a copy of gH for mol j's spheres (to be used in force/torque)
	  for (int h=0; h < j_intrasize; h++, k++) // intra
	    {
	      int kj  = intraMolPolList[h];
	      int key = gHmap[CFullSphereID(j,kj)]; 
	      assert(key < tGHs.size()); 
	      mols[j]->getKS(kj).getGH() = tGHs[key];

	    }
#endif
	  
	  // compute gLHN for mol j's spheres with interactlist but not polarized 
	  vector<int> intraMolInteractionList = mols[j]->getIntraMolInteractOnlyList();
	  for(int h=0; h < intraMolInteractionList.size(); h++)      
	    {
	      int kj = intraMolInteractionList[h]; 
	      prepareDTA_iself_LowMemory(mols, j, kj);
	    }
	  
	}// endif-inter=0
      
    }//end-delj
  

  /* //~~~~~~~~~~~~~~ expensive step!!! ignore if not used for torque ~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~compute sphere's gradient wrt own molecule on grid points ~~~~
  for (int i = 0; i < i_max; i++)
    for(int ki=0; ki < mols[i]->getNKS(); ki++)
      mols[i]->getpKS(ki)->computeExposedSurfaceGradientH(i);
  */
  
#if __DEBUGDIE__
  // debug
  double max = 100;
  for (int i=0; i<mols.size(); i++)
    for (int ki=0; ki<mols[i]->getNKS(); ki++)
    {
      CSolExpCenter KS=mols[i]->getKS(ki);

      CGradExpan gLHN = KS.getgLHN();
      bool bLS = false, bgLHN=false;
      
      
     if(!mols[i]->getKS(ki).IsEmptyInterPolList() )
       {
	 CLocalExpan LS = KS.getLS();
	 if(LS[0] > max || isnan(LS[0]) ) {  cout <<"LS > max"<<endl; bLS = true;}
       }
     

      if( (mols[i]->getInterMolPolList()).size() > 0)
	if(fabs(gLHN[0][0]) > max || fabs(gLHN[1][0]) > max || fabs(gLHN[2][0]) > max
	   || isnan(gLHN[0][0]) || isnan(gLHN[1][0]) || isnan(gLHN[2][0]) ) 
	  {cout <<"gLHN > max"<<endl; bgLHN = true;}
     
	
	if(bLS || bgLHN) 
	{
	  CMolecule::writeMolsPQR("died.L.pqr", mols);
	  CMolecule::saveConfig("died.L.config", mols);
	}
    }
#endif  

  double endg = read_timer();
  return;
}






// Generate a set of N points uniformly distributed on a unit sphere
void 
CMolecule::spherePts(int N, vector<REAL> & th, vector<REAL> & ph)
{
  th.resize(N,0.0);
  ph.resize(N,0.0);

  for (int i = 1; i <= N; i++)      
    {
      REAL h = -1.0 + 2.0*(i-1.0)/(N-1.0);
      th[i-1] = acos(h);

      if (i == 1 or i == N)
	ph[i-1] = 0;
      else 
	ph[i-1] = (ph[i-2] + 3.6/sqrt(N*(1.0 - h*h)));
	
      while(ph[i-1] > 2*M_PI)
	ph[i-1] -= 2*M_PI;

      while(ph[i-1] < -2*M_PI)
	ph[i-1] += 2*M_PI;
    }
}

// MOLECULE-RELATED FUNCTIONS


// extract charge and positions       
void 
CMolecule::extractCharges(int ki, CPnt cenKi, const vector<CPnt> &cpos, const vector<double> &chg, 
			  const vector<int> & clabel, 
			  vector<CPnt> &posAssigned, vector<double> & chgAssigned)
{	
  posAssigned.clear();
  chgAssigned.clear();

  for(int i=0; i<cpos.size();i++)
    if(clabel[i]==ki) 
      {
	posAssigned.push_back( cpos[i] - cenKi);
	chgAssigned.push_back( chg[i] );
      }
}

void 
CMolecule::extractCharges(int ki, CPnt cenKi, const vector<CPnt> &cpos, const vector<double> &chg,
			  const vector<int> & clabel, 
			  vector<CPnt> &posAssigned, vector<double> & chgAssigned, vector<CPnt> &allPosKi)
{

  extractCharges(ki, cenKi, cpos, chg, clabel, posAssigned, chgAssigned);

  for(int i=0; i<cpos.size(); i++) allPosKi[i] = cpos[i] - cenKi;

}




void 
CMolecule::writeMolExpansions(char* runname) const 
{


  for(int ki=0; ki < m_nks; ki++)
    {
      char fname1[MAX_CHARNUM]; 
      char fname2[MAX_CHARNUM]; 
      int n;

      n = sprintf(fname1, "%s.%d.%d.F.exp", runname, m_id,ki); 
      assert(n<=MAX_CHARNUM);
      writeExpansion(m_k[ki]->getF(), 1000, fname1);

      n = sprintf(fname2, "%s.%d.%d.H.exp", runname, m_id, ki);
      assert(n<=MAX_CHARNUM);
      writeExpansion(m_k[ki]->getH(), 1000,fname2);
    }

  return;
}
