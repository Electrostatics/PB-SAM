#include  <cstdlib>
#include <iostream>
#include <fstream>
#include <cfloat>
#include <stdio.h>
#include <time.h>

#include "BDnam.h"
#include "readutil.h"
#include "contact.h"


using namespace std;

// Static / General Variables
const char CBDnam::STATUS_STRING[4][10] = {"ESCAPED", "DOCKED", "RUNNING", "STUCK"};
vector< vector<CMolContact> > CBDnam::MOLCONTACTLIST;
int CBDnam::m_nMolType;

/******************************************************************/
/******************************************************************/
/**
 * Initialize constants for a nam run
 ******************************************************************/
void
CBDnam::initConstants(int nMolType)
{
  m_nMolType = nMolType;
  MOLCONTACTLIST.resize(1); // for NAM, we just need one def for protein 1
}

/******************************************************************/
/******************************************************************/
/**
 * BDnam add molecule to the contact list
 ******************************************************************/
void
CBDnam::addMolContact(int mol1type, const vector<CMolContact> &molcontactlist)
{
  MOLCONTACTLIST[0] = molcontactlist ;
}

/******************************************************************/
/******************************************************************/
/**
 * BDnam constructor
 ******************************************************************/
CBDnam::CBDnam(vector<char*> molfnames1, vector<char*> molfnames2, REAL idiel)
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
	
  CPnt initialrcen1(0,0,0), initialrcen2(0,0,0);
  for(int i=0; i<ncen1; i++) initialrcen1 += scen1[i]; initialrcen1 /= scen1.size();
  for(int i=0; i<ncen2; i++) initialrcen2 += scen2[i]; initialrcen2 /= scen2.size();
  cout <<"initialrcen1 : " <<initialrcen1<<endl;
  cout <<"initialrcen2 : " <<initialrcen2<<endl;
  
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
  CMolecule::generateMolSPX(scen1, srad1, m_SPxes1, m_nSPx1, m_neighs1 );
  CMolecule::generateMolExposedChargesFH(m_iF1, m_iH1, scen1, srad1, m_SPxes1, m_nSPx1,
                                         m_qSolvedF1, m_qSolvedH1, m_totalF1, m_totalH1);
  CMolecule::generateMolCells(scen1, srad1, initialrcen1, m_molcells1 );
  
  CMolecule::generateMolTypeIntraPolLists(scen1, srad1, intraRcutoff,
                                          m_intraPolLists_near1, intraPolLists_far1);
  
  CMolecule::computeMolTypeValues(initialrcen1, scen1, srad1, CHG1, POS1, idiel,
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
  CMolecule::generateMolCells(scen2, srad2, initialrcen2, m_molcells2 );
	
  CMolecule::generateMolTypeIntraPolLists(scen2, srad2, intraRcutoff,
                                          m_intraPolLists_near2, intraPolLists_far2);
  
  CMolecule::computeMolTypeValues(initialrcen2, scen2, srad2, CHG2, POS2, idiel,
                                  intraRcutoff,
                                  m_SPxes2, m_nSPx2, m_neighs2,
                                  m_iF2, m_iH2, m_qSolvedF2, m_qSolvedH2, m_totalF2, m_totalH2,
                                  m_intraPolLists_near2,
                                  intraPolLists_far2,
                                  m_LFs_intraSelf2, m_LHs_intraSelf2,
                                  m_LFs_intraSelf_far2, m_LHs_intraSelf_far2);
	
  // create 2 molecules
  m_mols.clear();
	
  m_mols.push_back( new CMolecule(0, initialrcen1, scen1, srad1, CHG1, POS1, idiel,
																	m_iMats1, intraRcutoff,
																	m_SPxes1, m_nSPx1, m_neighs1, m_intraPolLists_near1,
																	m_iF1, m_iH1, m_LFs_intraSelf1, m_LHs_intraSelf1,
																	m_LFs_intraSelf_far1, m_LHs_intraSelf_far1,
																	m_qSolvedF1, m_qSolvedH1, m_totalF1, m_totalH1, m_molcells1));
  
  m_mols.push_back( new CMolecule(1, initialrcen2, scen2, srad2, CHG2, POS2, idiel,
																	m_iMats2, intraRcutoff,
																	m_SPxes2, m_nSPx2, m_neighs2, m_intraPolLists_near2,
																	m_iF2, m_iH2, m_LFs_intraSelf2, m_LHs_intraSelf2,
																	m_LFs_intraSelf_far2, m_LHs_intraSelf_far2,
																	m_qSolvedF2, m_qSolvedH2, m_totalF2, m_totalH2, m_molcells2));

  REAL interRCutoff = 10;  // <<<<<<<<<<<<<<<<<<<<<<<<
  REAL interactRCutoff = 100;  // <<<<<<<<<<<<<<<<<<<<<<<<
  CMolecule::initMutualConstants(m_mols,interRCutoff, interactRCutoff, true );
  cout << "interRcutoff = "<<interRCutoff<<endl;
  cout << "interactRcutoff = "<<interactRCutoff<<endl;
	
  // initialize diffusion parameters
  m_Dtr = 0.03; m_Dr1 = 4.5e-5;  m_Dr2 = 4e-5;// 1brs
	
  double vm, rss;
  process_mem_usage(vm, rss);  cout << "Memory used to construct molecules: VM: " << vm << "; RSS: " << rss << endl;
  cout <<"BD construction complete"<<endl;
}


/******************************************************************/
/******************************************************************/
/**
 * BDnam run function: generally what happens here:
 1. two molecules position and orientation are set
 2. one is in the center, the other at some random
 distance from it in the simulation box
 both have random orientations
 3. Then the timestep is set, and for the time that
 they are not STUCK, DOCKED or ESCAPED, a BD simulation
 advances
 ******************************************************************/
CBDnam::STATUS
CBDnam::run(int &scount, char* runname)
{
  char fname[MAX_CHARNUM];
  sprintf(fname, "run_%s.log", runname);
  ofstream fout(fname);
	
  CQuat Q;
  // MOLECULE 1 (ORIGIN)
  m_mols[0]->setPos(CPnt());
  Q = CQuat::chooseRandom();
  m_mols[0]->setOrient(Q);
	
  // MOLECULE 2
  m_mols[1]->setPos(b_DIST*randOrient()); // <<ACTUAL CODE
  Q = CQuat::chooseRandom(); //<ACTUAL CODE
  m_mols[1]->setOrient(Q);
  
#ifdef __CELL__
  if( m_mols[1]->isCollided_cell(m_mols) )
	{cout <<"initial config collided_cell"<<endl;}
#else
  if( m_mols[1]->isCollided(m_mols) )
	{cout <<"initial config collided"<<endl;}
#endif
	
  vector<CPnt> force(2), torque(2);
  vector<REAL> pot;
  REAL fact = FACT2;
  STATUS status = RUNNING;
  scount = 0;
  CPnt dR2, dO1, dO2;
  double totaltime = 0;
  double start, end;

  while (status == RUNNING) // Removed upper limit && scount < 300000)
	{
		REAL dt = compute_dt();
		REAL dist = CMolecule::getC2CDist(m_mols[0], m_mols[1]);
		
#ifndef __DIFF__
		if (dist < f_DIST)
		{
			start = read_timer();
			
			bool bForce = 	CMolecule::computeForces(m_mols, force, torque);
			if(!bForce )
	    {cout <<"failed computeforce in step "<<scount<<endl; exit(1); }
			end = read_timer();
			totaltime += (end-start);
			
			dR2 = (m_Dtr*dt*fact)*force[1];
			dO1 = (m_Dr1*dt*IKbT*COUL_K)*torque[0];
			dO2 = (m_Dr2*dt*IKbT*COUL_K)*torque[1];
		}
		else
		{
			dR2.zero();
			dO1.zero();
			dO2.zero();
		}
		
#else
		dR2.zero();
		dO1.zero();
		dO2.zero();
#endif
	
		if (scount % 1000 == 0)
		{
			fout << scount << ") "
			<< m_mols[0]->getRCen() << m_mols[1]->getRCen()
			<< "-> " << CMolecule::getC2CDist(m_mols[0], m_mols[1]) << endl;
			fout << force[1] << " " << torque[1] <<endl;
			fout<<"dt: "<<dt<<" dR2: "<<dR2<<endl;
		}
		
		status = makeMove(dR2, dO1, dO2, dt);
		scount++;
	}
  
  totaltime /= double(scount);
  fout << CBDnam::STATUS_STRING[status] << endl;
  fout.close();
  return status;
}

/******************************************************************/
/******************************************************************/
/**
 * Function to move the second molecule
 Inputs: dR2 - differential distance step change
 between the two molecules
 d01 - differential rotation of the first molecue
 d02 - diff rotation of second
 dt - diff timestep
 ******************************************************************/
CBDnam::STATUS
CBDnam::makeMove(const CPnt & dR2, const CPnt & dO1,
								 const CPnt & dO2, REAL dt)
{
  REAL bdist = CMolecule::getC2CDist(m_mols[0], m_mols[1]);
	
  int c = 0;
  const int maxc = 500;
  bool bTranslated = false;
  while(!bTranslated && c < maxc)
	{
		c++;
		CPnt dR2_ = dR2 + CBDnam::getRandVec(sqrt(2*dt*m_Dtr));
		m_mols[1]->translate(dR2_);
		
		
#ifdef __CELL__
		if( m_mols[1]->isCollided_cell(m_mols) )  	m_mols[1]->untranslate();
#else
		if( m_mols[1]->isCollided(m_mols) )  	m_mols[1]->untranslate();
#endif
		else
			bTranslated = true;
	}
  if (c == maxc) { cout << "stuck in translate " << endl; return STUCK; }
	
  if (escaped(q_DIST)) return ESCAPED;
  
  // now rotate molecules
  CQuat Q1, Q2;
  CPnt dO1_, dO2_;
  
  dO1_ = dO1 + CBDnam::getRandVec(sqrt(2*dt*m_Dr1));
  dO2_ = dO2 + CBDnam::getRandVec(sqrt(2*dt*m_Dr2));
  Q1 = CQuat(dO1_, dO1_.norm());
  Q2 = CQuat(dO2_, dO2_.norm());
	
  REAL adist = CMolecule::getC2CDist(m_mols[0], m_mols[1]);
  if (adist > f_DIST) // too far away, store the rotation
	{
		m_rot1 = Q1*m_rot1;
		m_rot2 = Q2*m_rot2;
	}
  else  if (bdist > f_DIST) // just got within fdist, now rotate stored rotation
	{
		m_mols[0]->rotate(m_rot1);
		m_mols[1]->rotate(m_rot2);
		
		m_rot1.identity();
		m_rot2.identity();
	}
  else
	{
		// rotate mol1
		bool bRotated = false; c = 0;
		while(!bRotated && c < maxc)
		{
			c++;
			dO1_ = dO1 + CBDnam::getRandVec(sqrt(2*dt*m_Dr1));
			Q1 = CQuat(dO1_, dO1_.norm());
			
#ifdef __CELL__
			if (! m_mols[0]->willRotCollide_cell(m_mols, Q1) )
#else
				if (! m_mols[0]->willRotCollide(m_mols, Q1) )
#endif
				{
					m_mols[0]->rotate(Q1);
					bRotated = true;
				}
		}
		if (c == maxc)
		{
			cout << "stuck in rotate1 " << endl;
			m_rot1 = Q1*m_rot1;
			return STUCK;
		}

		// rotate mol2
		bRotated = false; c = 0;
		while(!bRotated && c < maxc)
		{
			c++;
			CPnt dO2_ = dO2 + CBDnam::getRandVec(sqrt(2*dt*m_Dr2));
			Q2 = CQuat(dO2_, dO2_.norm());
#ifdef __CELL__
			if (! m_mols[1]->willRotCollide_cell(m_mols, Q2) )
#else
				if (! m_mols[1]->willRotCollide(m_mols, Q2) )
#endif
	      {
					m_mols[1]->rotate(Q2);
					bRotated = true;
				}
		}
		if (c == maxc)
		{
			cout << "stuck in rotate2 " << endl;
			m_rot2 = Q2*m_rot2;
			return STUCK;
		}
		
	}
  if (IsDocked(m_mols[0], m_mols[1]))
    return DOCKED;
	
  else
    return RUNNING;
}

/******************************************************************/
/******************************************************************/
/**
 * This function checks to see whether two molecules
 have docked general, for proteins to bind multiple partners
 ******************************************************************/
bool
CBDnam::IsDocked(CMolecule* mol1, CMolecule* mol2) const
{
  int m1 = mol1->getMolType();
	
  int nDockDef = MOLCONTACTLIST[m1].size();
  for(int h=0; h<nDockDef; h++)
	{
		CMolContact mcon = MOLCONTACTLIST[m1][h];
		if(mcon.getMol2Type()!= mol2->getMolType() )  continue;
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

/******************************************************************/
/******************************************************************/
/**
 * Depending on the distance between the two molecules,
 we will modify our timestep as they get closer, forces
 are larger and we need a smaller timestep
 ******************************************************************/
inline REAL
CBDnam::compute_dt() const
{
  REAL d = CMolecule::getSepDist(m_mols[0], m_mols[1]);
	
  if (d-5 > 0)
    return 2.0 + (d-5)/15;
  else
    return 2.0;
}

