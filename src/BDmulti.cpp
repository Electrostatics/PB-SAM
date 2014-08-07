#include  <cstdlib>
#include <iostream>
#include <fstream>
#include <cfloat>
#include <stdio.h>
#include <time.h>


#include "BDmulti.h"
#include "readutil.h"

#define WRITEFREQ 500
#define DISTCUTOFF_TIME 10 // angstroms
REAL CBDmulti::MINDIST, CBDmulti::MAXDIST;
int CBDmulti::m_nMolType;
vector< vector<CMolContact> > CBDmulti::MOLCONTACTLIST;

using namespace std;

// static functions
void 
CBDmulti::initConstants(int nMolType)
{
  
  if(DISTCUTOFF_TIME > CMolecule::m_interactRcutoff )
    {
      MAXDIST = DISTCUTOFF_TIME; 
      MINDIST = CMolecule::m_interactRcutoff;
    } 
  else
    {
      MAXDIST = CMolecule::m_interactRcutoff;  
      MINDIST = DISTCUTOFF_TIME;
    } 
  m_nMolType = nMolType;
  MOLCONTACTLIST.resize(1); // for 2 species, we just need one def for protein 1 
}

void
CBDmulti::addMolContact(int mol1type, const vector<CMolContact> &molcontactlist)
{
  MOLCONTACTLIST[0] = molcontactlist ;
}

// member function


// construct from array placement with random orientation
CBDmulti::CBDmulti(int np1, int np2, 
		   const vector<char*> &molfnames1, const vector<char*> &molfnames2, 
		   REAL idiel) : m_np1(np1), m_np2(np2)
{

  initMolTypes(molfnames1, molfnames2, idiel);

  double vm, rss;
  process_mem_usage(vm, rss);  cout << "VM: " << vm << "; RSS: " << rss << endl;	   

  cout <<"BD initialization complete for "<<m_mols.size()<<" molecules." <<endl;

}

// read in and precompute values for each species
void
CBDmulti::initMolTypes(const vector<char*> &molfnames1,
		      const vector<char*> &molfnames2, REAL idiel)
{
  char* pqrfile1  = molfnames1[0];
  char* imatpath1 = molfnames1[1];
  char* exppath1  = molfnames1[2];
  char* expname1  = molfnames1[3];

  char* pqrfile2  = molfnames2[0];
  char* imatpath2 = molfnames2[1];
  char* exppath2  = molfnames2[2];
  char* expname2  = molfnames2[3];

  // read in charge + centers                                                                         
  vector<CPnt> scen1, POS1, scen2, POS2;
  vector<double> srad1, CHG1,  srad2, CHG2;
  readpqr(pqrfile1, POS1, CHG1 ,srad1, scen1);
  readpqr(pqrfile2, POS2, CHG2 ,srad2, scen2);

  const int ncen1 = scen1.size();
  const int ncen2 = scen2.size();

  m_initialrcen1 = CPnt(0,0,0), m_initialrcen2 = CPnt(0,0,0);
  for(int i=0; i<ncen1; i++) m_initialrcen1 += scen1[i]; m_initialrcen1 /= scen1.size();
  for(int i=0; i<ncen2; i++) m_initialrcen2 += scen2[i]; m_initialrcen2 /= scen2.size();
  cout <<"initialrcen1 : " <<m_initialrcen1<<endl;
  cout <<"initialrcen2 : " <<m_initialrcen2<<endl;
  

  //////////////// read in solved expansions //////////////                                           
  readMats(m_iMats1, N_POLES, imatpath1, ncen1);
  readMats(m_iMats2, N_POLES, imatpath2, ncen2);
  REAL intraRcutoff_dum;
  readExpansionsP(m_iF1, m_iH1, N_POLES, exppath1, expname1, ncen1, intraRcutoff_dum);
  readExpansionsP(m_iF2, m_iH2, N_POLES, exppath2, expname2, ncen2, intraRcutoff_dum);

  /////////////// generate values for each moltype ///////////////                                   
  const double intraRcutoff = 30;  //<<<<<<<<                                                        
  vector< vector<int> > intraPolLists_far1, intraPolLists_far2;

  // protein 1  
  //'generateMolSPX' will generate surface points (S. Liu)
  CMolecule::generateMolSPX(scen1, srad1, m_SPxes1, m_nSPx1, m_neighs1 );
  //'generateMolExposedChargesFH' will calculate surface charge on each sphere
  // represented by multipole expansions of F and H (S. Liu)
  CMolecule::generateMolExposedChargesFH(m_iF1, m_iH1, scen1, srad1, m_SPxes1, m_nSPx1,
                                         m_qSolvedF1, m_qSolvedH1, m_totalF1, m_totalH1);
  CMolecule::generateMolCells(scen1, srad1, m_initialrcen1, m_molcells1 );

  //'generateMolTypeIntraPolLists' will generate lists for intra-molecular polarization (S. Liu)
  CMolecule::generateMolTypeIntraPolLists(scen1, srad1, intraRcutoff,
                                          m_intraPolLists_near1, intraPolLists_far1);
  //see "molecule.cpp" for interpretations of 'computeMolTypeValues'  (S. Liu)
  CMolecule::computeMolTypeValues(m_initialrcen1, scen1, srad1, CHG1, POS1, idiel,
                                  intraRcutoff,
                                  m_SPxes1, m_nSPx1, m_neighs1,
                                  m_iF1, m_iH1, m_qSolvedF1, m_qSolvedH1, m_totalF1, m_totalH1,
                                  m_intraPolLists_near1,
                                  intraPolLists_far1,
                                  m_LFs_intraSelf1, m_LHs_intraSelf1,
                                  m_LFs_intraSelf_far1, m_LHs_intraSelf_far1);

  // protein 2                                                                                       
  CMolecule::generateMolSPX(scen2, srad2, m_SPxes2, m_nSPx2, m_neighs2 );
  CMolecule::generateMolExposedChargesFH(m_iF2, m_iH2, scen2, srad2, m_SPxes2, m_nSPx2,
                                         m_qSolvedF2, m_qSolvedH2, m_totalF2, m_totalH2);
  CMolecule::generateMolCells(scen2, srad2, m_initialrcen2, m_molcells2 );

  CMolecule::generateMolTypeIntraPolLists(scen2, srad2, intraRcutoff,
                                          m_intraPolLists_near2, intraPolLists_far2);

  CMolecule::computeMolTypeValues(m_initialrcen2, scen2, srad2, CHG2, POS2, idiel,
                                  intraRcutoff,
                                  m_SPxes2, m_nSPx2, m_neighs2,
                                  m_iF2, m_iH2, m_qSolvedF2, m_qSolvedH2, m_totalF2, m_totalH2,
                                  m_intraPolLists_near2,
                                  intraPolLists_far2,
                                  m_LFs_intraSelf2, m_LHs_intraSelf2,
                                  m_LFs_intraSelf_far2, m_LHs_intraSelf_far2);
    
  // create molecules
  m_mols.clear();
  
  for(int i=0; i < m_np1; i++)
    {
      m_mols.push_back( new CMolecule(0, m_initialrcen1, scen1, srad1, CHG1, POS1, idiel,
				      m_iMats1, intraRcutoff,
				      m_SPxes1, m_nSPx1, m_neighs1, m_intraPolLists_near1,
				      m_iF1, m_iH1, m_LFs_intraSelf1, m_LHs_intraSelf1,
				      m_LFs_intraSelf_far1, m_LHs_intraSelf_far1,
				      m_qSolvedF1, m_qSolvedH1, m_totalF1, m_totalH1, m_molcells1));
    }

  for(int i=0; i < m_np2; i++)
    {
      m_mols.push_back( new CMolecule(1, m_initialrcen2, scen2, srad2, CHG2, POS2, idiel,
				      m_iMats2, intraRcutoff,
				      m_SPxes2, m_nSPx2, m_neighs2, m_intraPolLists_near2,
				      m_iF2, m_iH2, m_LFs_intraSelf2, m_LHs_intraSelf2,
				      m_LFs_intraSelf_far2, m_LHs_intraSelf_far2,
				      m_qSolvedF2, m_qSolvedH2, m_totalF2, m_totalH2, m_molcells2));
    }

  assert(m_mols.size() == m_np1+m_np2);
  if(m_np1 > 0) cout << "index 0 to " <<m_np1-1<<" are moltype 0"<<endl;
  if(m_np2 > 0) cout << "index "<<m_np1<<" to "<<m_np1+m_np2-1<<" are moltype 1"<<endl; 
  // =========================================================
  //"interRCutoff", which gives the value of "m_interRCutoff" in class "CMolecule",
  //as you will see in function CMolecule::generateInterXFormsForPolarize_LowMemory,
  // if the distance between sphere i in one molecule and sphere j is less than 
  // 'interRCutoff', we will not polarize the gradients of their expansions, which 
  // will be used for force calculation (S. Liu)
  const REAL interRCutoff = 10; // <<<<<<<<<<<<<<<<<<<<<<<<
  //interactRCutoff, which gives the valule of "m_interactRCutoff" in class CMolecule
  // is the cutoff radius for doing mutual polarization, i.e., if the distance 
  //between two molecules is greater than interacRCutoff, this pair won't be mutually 
  //polarized (S. Liu)
  const REAL interactRCutoff = 100; // <<<<<<<<<<<<<<<<<<<<<<<<
  CMolecule::initMutualConstants(m_mols, interRCutoff, interactRCutoff, true );
  cout << "interRcutoff = "<<interRCutoff<<endl;
  cout << "interactRcutoff = "<<interactRCutoff<<endl;
  // =========================================================  
  
  // initialize diffusion parameters
  m_Dtr.resize(m_nMolType);
  m_Dr.resize(m_nMolType);
  
  m_Dtr[0] = 0.015; m_Dtr[1] = 0.015; // 1BRS values, translational diffusion coefficient
  m_Dr[0] = 4.5e-5; m_Dr[1] = 4e-5;//rotational diffusion coefficient (S. Liu)
  /* 
  // debug
  CPnt cen0 = CPnt(-118.136, 4.67571, 93.6392); //i
  CQuat Q0 = CQuat(-0.679896,CPnt(-0.574849, 0.0429259, -0.453263));
  CPnt cen1 = CPnt(-94.3519, 32.3896, 114.463); //j
  CQuat Q1 = CQuat(-0.11112, CPnt(-0.77886, 0.565594, -0.24725));
  m_mols[0]->setPos(cen0);
  m_mols[0]->setOrient(Q0);
  m_mols[1]->setPos(cen1);
  m_mols[1]->setOrient(Q1);

  m_mols[0]->rotateRotCoeff();
  m_mols[1]->rotateRotCoeff();
  
  for(int ki=0; ki<m_mols[0]->getNKS(); ki++)
    for(int kj=0; kj<m_mols[1]->getNKS(); kj++)
      {
	//	cout <<ki<<" "<<kj<<endl;
	m_mols[1]->getKS(kj).rotateExpansions();
	m_mols[0]->reexpandLSFromSelfVal_debug(m_mols, 0, ki, 1, kj);
      }
include
  cout << "j>i done"<<endl;

  for(int ki=0; ki<m_mols[0]->getNKS(); ki++)
    for(int kj=0; kj<m_mols[1]->getNKS(); kj++)
      {
        m_mols[0]->getKS(ki).rotateExpansions();
        m_mols[1]->reexpandLSFromSelfVal_debug(m_mols, 1, kj, 0, ki);
      }
  cout <<"i>j done"<<endl;
  

  m_mols[1]->getKS(7).rotateExpansions();                                                                   
  m_mols[0]->reexpandLSFromSelfVal_debug(m_mols, 0, 7, 1, 7);                                              
  */
  return;

 }


// initialize positions, orientation and moltypes on a lattice
void 
CBDmulti::resetLattice()
{
  // generate all positions and orientation in the array, regardless of moltype
  const int number_of_monomers = m_np1 + m_np2;
  const int number_per_side = int(ceil( pow(number_of_monomers, 1.0/3.0)));
  vector<CPnt> rcens;
  vector<CQuat> rots;

  // determine the size of the biggest protein 
  double maxR = DBL_MIN;
  for (int i=0; i<m_mols.size(); i++)
    { 
      double maxr_i = m_mols[i]->getMaxR();
      if(maxr_i > maxR) maxR = maxr_i; 
    }
  REAL fact = CSystem::BOXLENGTH / double(number_per_side); 
  cout <<"maxR for all moltype "<<maxR<<endl;
  cout <<"c2c distance between molecules "<<fact<<endl; 
  assert (fact > 2*maxR); 

  const REAL offset = 0.5*CSystem::BOXLENGTH;  
  int n = 0; 
  for(int i=0; i<number_per_side; i++)
    {
      REAL px = -offset + (0.5+i)*fact;
      for(int j=0; j<number_per_side ; j++)
	{
	  REAL py = -offset + (0.5+j)*fact;
	  for(int k=0; k<number_per_side; k++)
	    {
	      REAL pz = -offset + (0.5+k)*fact;
	      rcens.push_back( CPnt(px,py,pz) );
	      rots.push_back (  CQuat::chooseRandom() );
	    }
	}
    }

  // assign molecules randomly to the generated positions/orientation
  vector<bool> bOccupied(number_of_monomers, false);
  for(int i=0; i<number_of_monomers; i++)
    {
      bool bSpotFound = false;
      while ( !bSpotFound )
	{
	  int k = (int) floor(drand48()*  number_of_monomers );
	  if( bOccupied[k] == false) 
	    {
	      m_mols[i]->setPos( rcens[k] );
	      m_mols[i]->setOrient( rots[k] ); 
	      bSpotFound = true; 
	      bOccupied[k] = true;
	    }
	}
    }

  // check that nothing collides (use direct method instead of cell)
  for(int i=0; i < number_of_monomers; i++)
    if( m_mols[i]->isCollided(m_mols) )  
      {
	cout <<"Error : initial config clashes"<<i<<endl;
	exit(1);
      }

  return;
}


void
CBDmulti::restartConfig(char* configfname)
{
  vector<CPnt> rcens;
  vector<CQuat> rots; 
  vector<int> moltypes, moltypeid, nps;
  readConfigWithMolType(configfname, rcens, rots, moltypeid, moltypes, nps); 

  // for now, two species
  assert (nps.size() == 2);
  assert (nps.size() == moltypes.size());

  // check that the number of particles agree
  assert (moltypes[0] == 0 && moltypes[1] == 1);
  assert (nps[0] == m_np1 && nps[1] == m_np2);
  const int number_of_monomers = m_np1 + m_np2;
  assert (rcens.size() == number_of_monomers && rots.size() == number_of_monomers);

  // set position and orientation  
  for(int i=0; i<number_of_monomers; i++)
    {
      m_mols[i]->setPos( rcens[i] );
      m_mols[i]->setOrient( rots[i] );
    }
  

  // check that nothing collides (use direct method instead of cell)                                             
  for(int i=0; i < number_of_monomers; i++)
    if( m_mols[i]->isCollided(m_mols) )
      {
        cout <<"Error : initial config clashes in restart "<<i<<endl;
        exit(1);
      }




  return;

}


// bd run for 1 trajectory
double 
CBDmulti::runFirstDockTime(double maxtime, char* runname)
{
  
  char cfgname[200];
  sprintf(cfgname, "restart_%s.config", runname);

  int nmol = m_mols.size();

  double totaltime = 0;
  double t = 0;
  int n = 0; 
  bool bDocked = false;  

  // reset molecule positions and orientations
  // resetLattice();

  CMolecule::writeMolsPQR("initial.pqr", m_mols);                                                               
  double start = read_timer();
  double avgtime = 0.0;
  
  while(!bDocked && t < maxtime )
    //while(!bDocked && n < 1)
    {
      REAL minDist = computeMinSepDist();

      bool bWrite = (n % WRITEFREQ == 0); 
      if(bWrite) cout <<n<<") Timings ";

      double t1 = read_timer();
      REAL dt = propagate(minDist, bWrite); //this is the function that did force calculation 
      //and moved molecules (S. Liu)
      t += dt;     

      for(int i=0; i<nmol; i++)
	{
	  if(m_mols[i]->getMolType() != 0 ) continue;

	  for(int j=0; j<nmol; j++)
	    {
	      if(m_mols[j]->getMolType() != 1 ) continue;

	      if (IsDocked(m_mols[i], m_mols[j]))
		{
		  bDocked = true;
		  break;
		}
	      if(bDocked) break;
	    }
	}

      double t2 = read_timer(); 
      avgtime += t2-t1;

      if (bWrite)
	{      
	  cout <<n<<") Timings: dock "<<t2-t1
	       <<" Time for this "<<WRITEFREQ<<" steps: "<<avgtime/(int) WRITEFREQ<<"s"<< endl;
	  avgtime = 0.0;
	  cout <<n<< ") mindist "<<minDist <<" dt "<<dt<<"[ps]; simulated "<<t*PS_TO_NS<<" [ns]"<<endl; 
	  CMolecule::saveConfig(cfgname, m_mols);
	  
	}
      
      n++;
      
    }//end while
  
  //  CMolecule::writeMolsPQR("finished.pqr", m_mols);

  double end = read_timer();
  cout <<"Docked?  " <<bDocked<<" out of "<<n<<" steps; comp time for this traj: "<<end-start<<endl;

  return t * PS_TO_NS;  
}

double
CBDmulti::propagate(double minDist, bool bWriteTime)
{

  // compute minimum distance to decide timestep and whether to compute forces
  bool bComputeForce = (minDist < CMolecule::m_interactRcutoff ) ;  
  double dt = compute_dt(minDist);
 
  // compute forces
  const int nmol = m_mols.size();
  vector<CPnt> force(nmol), torque(nmol);
  vector<CPnt> dR(nmol), dO(nmol); 
#ifndef __DIFF__

  double t1 = read_timer();

  if(bComputeForce) {
  // this is the function that computes forces and torques (S. Liu)
  bool bBlown = CMolecule::computeForces(m_mols, force, torque); 
  if(!bBlown) 
    {
      cout <<"died somewherein computeforces"<<endl; 
      CMolecule::saveConfig("died.config", m_mols);
      CMolecule::writeMolsPQR("died.pqr", m_mols);                                                                     
      exit(1);
    }
  }
#else
  force.clear(); force.resize(nmol);
  torque.clear(); torque.resize(nmol);
#endif

  double t2 = read_timer();      
  // move molecules
  for(int i=0; i<nmol; i++) 
    {
      int m = m_mols[i]->getMolType();
      dR[i] = (m_Dtr[m]*dt*FACT2)*force[i];
      dO[i] = (m_Dr[m]*dt*IKbT*COUL_K)*torque[i];
      
    }// end-i

  makeMove(dR, dO, dt); //this is the function that moves molecules (S. Liu)
  double t3 = read_timer();

  if(bWriteTime) cout <<"  computeForces "<<t2-t1<<" move "<<t3-t2<<endl;

  return dt;      
}

//this is the function that moves molecules (S. Liu)
void
CBDmulti::makeMove(const vector<CPnt> & dR, const vector<CPnt> & dO, REAL dt) 
{
  
  const int maxc = 500; // for now, stay put if it doesnt move after maxc trials
  
  for(int i=0; i<m_mols.size(); i++)
    {
      int c;
      int m = m_mols[i]->getMolType(); 

       // translate
      CPnt dR_; 
      bool bTranslated = false; c=0;
      while(!bTranslated && c < maxc)
	{
	  c++;

      //the displacement is the sum of contribution from electrostatic forces and 
      //random displacement, CBDmulti::getRandVec will generate a random displacement
      //whose distribution is Gaussian with standard deviation 2*dt*m_Dtr[m] (S. Liu)
	  dR_ = dR[i] + CBDmulti::getRandVec(sqrt(2*dt*m_Dtr[m]));
	  m_mols[i]->translate(dR_);
	
#ifdef __CELL__
	  bool bCollided =  m_mols[i]->isCollided_cell(m_mols); 
#else
	  bool bCollided =  m_mols[i]->isCollided(m_mols); 
#endif
	  if(bCollided)	
	      m_mols[i]->untranslate(); //if one molecule collides with another, the move 
          //will be given up and a new move will be generated (S. Liu)
	  else 
	    bTranslated = true;
	}

      //      if (c == maxc) { cout << "stuck in translate " << endl; exit(1);}
            
      
      // rotate 
	      
      CQuat Q;   //CQuat is the class for quarterion (S. Liu)
      CPnt dO_;
      bool bRotated = false; c = 0;
      while(!bRotated && c < maxc)
	{
	  c++;
	  dO_ = dO[i] + CBDmulti::getRandVec(sqrt(2*dt*m_Dr[m]));
	  Q = CQuat(dO_, dO_.norm());
	  
#ifdef __CELL__
	  if(! m_mols[i]->willRotCollide_cell(m_mols, Q) ) 
#else
	    if(! m_mols[i]->willRotCollide(m_mols, Q) )
#endif
	    {
	      m_mols[i]->rotate(Q);	
	      bRotated = true;
	    }
	}
        
      //      if (c == maxc) 	{ cout << "stuck in rotate " <<i<< endl; exit(1); }
       
      
    }

    return;
  
}

inline REAL 
CBDmulti::compute_dt(double d) const
{

  if (d - DISTCUTOFF_TIME > 0)
    return 2.0 + (d - DISTCUTOFF_TIME)/15;
  else
    return 2.0;
}


// assume DISTCUTOFF_TIME = CMolecule::m_interRcutoff to simplify
REAL 
CBDmulti::computeMinSepDist()     //function for calculation of minimum separation between molecuels in a configuration (S. Liu)
{ 
  double minsepdist = DBL_MAX;
  double minM2Msepdist = DBL_MAX;
  const int nmol = m_mols.size();
  bool bTrackingS2S = false;
  
  for(int i=0; i<nmol; i++)
      for(int j=i+1; j<nmol; j++)
	{
	  double dM2M = CMolecule::getSepDist(m_mols[i], m_mols[j]); 

	  // first check if they are close enough that we care about
	  if(dM2M < MAXDIST )
	    {
	      if(!bTrackingS2S) bTrackingS2S = true;
	      // evaluate it until it reaches the small limit that we care about
	      for(int ki=0; ki<m_mols[i]->getNKS(); ki++)
		{

		  for(int kj=0; kj<m_mols[j]->getNKS(); kj++)
		    {
		      double d = CMolecule::getSepDist(m_mols[i], ki, m_mols[j], kj);
		      if( d < minsepdist )
			{
			  if( d < MINDIST ) { return MINDIST; }
			  else minsepdist = d;
			}
		    }// end kj 

		}// end k
	    }
	  // otherwise keep track of mol-to-mol distance just in case
	  else if (!bTrackingS2S)
	    {
	      if(dM2M < minM2Msepdist) minM2Msepdist  = dM2M;
	    }
	  
	}//end (i,j)
    
  
  if(bTrackingS2S) return  minsepdist;
  else  return (minM2Msepdist > 0 ?  minM2Msepdist : 0); 
  
}

/*
// assume DISTCUTOFF_TIME = CMolecule::m_interRcutoff to simplify
REAL 
CBDmulti::computeMinSepDist()     
{ 
  double minsepdist = DBL_MAX;
  double minM2Msepdist = DBL_MAX;
  const int nmol = m_mols.size();
  bool bTrackingS2S = false;
  //bool bMinReached = false;
  
  //#pragma omp parallel shared(bMinReached, bTrackingS2S, minsepdist, minM2Msepdist)
  {
    for(int i=0; i<nmol; i++)
      for(int j=i+1; j<nmol; j++)
	{
	  double dM2M = CMolecule::getSepDist(m_mols[i], m_mols[j]); 
	  // first check if they are close enough that we care about
	  if(dM2M < MAXDIST )
	    {
	      if(!bTrackingS2S) bTrackingS2S = true;
	      // evaluate it until it reaches the small limit that we care about
	      for(int ki=0; ki<m_mols[i]->getNKS(); ki++)
		{
		  //#pragma omp for
		  for(int kj=0; kj<m_mols[j]->getNKS(); kj++)
		    {
		      double d = CMolecule::getSepDist(m_mols[i], ki, m_mols[j], kj);
		      if( d < minsepdist )
			{
			  if( d < MINDIST ) { return MINDIST; }
			  else minsepdist = d;
			}
		    }// end kj (parallel for)

		}// end k
	    }
	  // otherwise keep track of mol-to-mol distance just in case
	  else if (!bTrackingS2S)
	    {
	      //#pragma omp master 
	      {
		if(dM2M < minM2Msepdist) minM2Msepdist  = dM2M;
	      }
	    }
	  
	}//end (i,j)
    
  }// end parallel
  
  if(bTrackingS2S) return  minsepdist;
  else  return (minM2Msepdist > 0 ?  minM2Msepdist : 0); 
  
}
*/

// general, for proteins to bind multiple partners                                        
bool
CBDmulti::IsDocked(CMolecule* mol1, CMolecule* mol2) const
{
  int m1 = mol1->getMolType();

  int nDockDef = MOLCONTACTLIST[m1].size();
  for(int h=0; h<nDockDef; h++)
    {
      CMolContact mcon = MOLCONTACTLIST[m1][h];

      if(mcon.getMol2Type()!= mol2->getMolType() )  continue;
      //      assert(mcon.getMol2Type()== mol2->getMolType() ) ;                          

      int ncon = 0;
      vector<CContact> clist = mcon.getContactList();
      for(int k=0; k<clist.size() && ncon < mcon.getNContact(); k++)
        {
          if( CMolecule::getSepDist(mol1, clist[k].getID1(), mol2, clist[k].getID2())
              <= clist[k].getDist() ) ncon++;

          if( ncon == mcon.getNContact()) return true;
        }


    } // end all definitions                                                              

  return false;
}
