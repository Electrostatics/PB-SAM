#include "molecule.h"
bool CMolecule::m_bInfinite = false;
double CMolContactRigid::angleTolerance;
double CMolContactRigid::maxC2CDist2 = DBL_MIN;

REAL
CMolecule::getC2CDist2(const CMolecule *mol1, const CMolecule *mol2)

{
  assert(mol1->getID() != mol2->getID());
  return (CSystem::pbcPos(mol1->getRCen() - mol2->getRCen())).normsq();
}
REAL
CMolecule::getC2CDist(const CMolecule *mol1, const CMolecule *mol2)

{
  assert(mol1->getID() != mol2->getID());
  return (CSystem::pbcPos(mol1->getRCen() - mol2->getRCen())).norm();
}

REAL
CMolecule::getS2SDist(const CMolecule *mol1, int ki, const CMolecule *mol2, int kj)
{
  assert(mol1->getID() != mol2->getID());
  return (CSystem::pbcPos( mol1->getRCen() + mol1->getKS(ki).getCenRot()
													- mol2->getRCen() - mol2->getKS(kj).getCenRot() ).norm());
}

REAL
CMolecule::getSepDist(const CMolecule *mol1,const CMolecule *mol2)
{
  return getC2CDist(mol1, mol2) - mol1->getMaxR() - mol2->getMaxR();
}

REAL
CMolecule::getSepDist(const CMolecule *mol1, int ki, const CMolecule *mol2, int kj)
{
  return getS2SDist(mol1, ki, mol2, kj)
	- mol1->getKS(ki).getRad()
	- mol2->getKS(kj).getRad();
}
bool
CMolecule::willRotCollide(const vector<CMolecule*> &mols, CQuat dQ) const
{
  const int nks = getNKS();
  const int nmol = mols.size();
  const CPnt rcen_i = getRCen();
  const CQuat Q =  dQ * getOrient(); // Order of multiplying matters!
	
  vector<CPnt> testCenRot(nks);
	
  for(int ki=0; ki<nks; ki++)
    testCenRot[ki] = Q * getKS(ki).getCen();
  
  //#pragma omp parallel
  {
		
    for(int j=0; j<nmol; j++)
		{
			if(j==m_id || getSepDist(this, mols[j]) > 0) continue;;
			
			CPnt rcen_j = mols[j]->getRCen();
			
			for(int ki=0; ki<nks; ki++)
			{
				double ri = getKS(ki).getRad();
				CPnt cki = rcen_i + testCenRot[ki];
				
				//#pragma omp for
				for(int kj=0; kj<mols[j]->getNKS(); kj++)
	      {
					double allowedOverlap = 0.0;
					if( ri  + mols[j]->getKS(kj).getRad()
						 - CSystem::pbcPos(cki - rcen_j - mols[j]->getKS(kj).getCenRot()).norm()
						 > allowedOverlap )
						return true;
	      }// end kj (parallel)
				
			}//end ki
			
		}// end j
		
  }//end parallel
  return false;
}

#ifdef __CELL__
bool
CMolecule::willRotCollide_cell(const vector<CMolecule*> &mols, CQuat dQ) const
{
  
  const int nks = getNKS();
  const int nmol = mols.size();
  const CPnt rcen_i = getRCen();
  const CQuat Q =  dQ * getOrient(); // Order of multiplying matters!
  
  const vector<CMolCell> &molcelli = getMolCells();
  const vector<CPnt> &cellCensi = getCellCens();
	
	
  vector<CPnt> testCenRot(nks);
  vector<CPnt> testCellCenRot(molcelli.size());
	
  for(int ci=0; ci<molcelli.size(); ci++)    testCellCenRot[ci] = Q * molcelli[ci].getCen();
  for(int ki=0; ki<nks; ki++)    testCenRot[ki] = Q * getKS(ki).getCen();
  
  //#pragma omp parallel
  {
		
    for(int j=0; j<nmol; j++)
		{
			if(j==m_id || getSepDist(this, mols[j]) > 0) continue;;
			
			const vector<CMolCell> &molcellj = mols[j]->getMolCells();
			const vector<CPnt> &cellCensj = mols[j]->getCellCens();
			
			const CPnt rcen_j = mols[j]->getRCen();
			
			for(int ci=0; ci<cellCensi.size(); ci++)
			{
				CPnt p_ci = rcen_i + testCellCenRot[ci];
				double rad_ci = molcelli[ci].getRad();
				
				for(int cj=0; cj<cellCensj.size(); cj++)
	      {
					CPnt p_cj = rcen_j + cellCensj[cj];
					double rad_sum = rad_ci + molcellj[cj].getRad();
					
					if( CSystem::pbcPos(p_ci - p_cj).normsq() > rad_sum*rad_sum )
						continue;
					
					// if cells overlap, then check each continuent sphere
					vector<int> sphereListI = molcelli[ci].getSphereList();
					vector<int> sphereListJ = molcellj[cj].getSphereList();
					
					for(int k=0; k<molcelli[ci].getSphereListSize(); k++)
					{
						int ki = sphereListI[k];
						double ri = getKS(ki).getRad();
						CPnt cki = rcen_i + testCenRot[ki];
						
						for(int h=0; h<molcellj[cj].getSphereListSize(); h++)
						{
							int kj =  sphereListJ[h];
							double allowedOverlap = 0.0;
							
							if( ri  + mols[j]->getKS(kj).getRad()
								 - CSystem::pbcPos(cki - rcen_j - mols[j]->getKS(kj).getCenRot()).norm()
								 > allowedOverlap )
								return true;
							
						}//end kj
					}// end ki
					
	      }//end cj
			}//end ci
			
		}// end j
    
  }//end parallel
	
  return false;
}
#endif

bool
CMolecule::isCollided(const vector<CMolecule*> &mols) const
{
  const int nks = getNKS();
  const int nmol = mols.size();
  const CPnt rcen_i = getRCen();
	
	for(int j=0; j<nmol; j++)
	{
		if(j==m_id || getSepDist(this, mols[j]) > 0) continue;
		
		for(int ki=0; ki<nks; ki++)
		{
			
			//#pragma omp parallel for
			for(int kj=0; kj<mols[j]->getNKS(); kj++)
	    {
	      double allowedOverlap = 0.0;
	      if( getSepDist(this, ki, mols[j], kj) < - allowedOverlap) {
					cout <<"collided in isCollided: i="<<m_id<<" ki="<<ki<<" j="<<j<<" kj="<<kj<<endl;
					cout <<"pos_ki: "<<getRCen() + getKS(ki).getCenRot()<<" pos_kj: "<<mols[j]->getRCen() + mols[j]->getKS(kj).getCenRot()<<endl;
					return true;
	      }
	    }
		}
	}// end j
	return false;
}

#ifdef __CELL__
// don't use : has bug
// for one, cell positions are not correct when restarting
bool
CMolecule::isCollided_cell(const vector<CMolecule*> &mols) const
{
  const int nks = getNKS();
  const int nmol = mols.size();
  const CPnt rcen_i = getRCen();
  
  const vector<CMolCell> &molcelli = getMolCells();
  const vector<CPnt> &cellCensi = getCellCens();
  
  for(int j=0; j<nmol; j++)
	{
		if(j==m_id || getSepDist(this, mols[j]) > 0) {continue;}
		
		const vector<CMolCell> &molcellj = mols[j]->getMolCells();
		const vector<CPnt> &cellCensj = mols[j]->getCellCens();
		
		const CPnt rcen_j = mols[j]->getRCen();
		
		for(int ci=0; ci<cellCensi.size(); ci++)
		{
			CPnt p_ci = rcen_i + cellCensi[ci];
			double rad_ci = molcelli[ci].getRad();
			
			for(int cj=0; cj<cellCensj.size(); cj++)
	    {
	      CPnt p_cj = rcen_j + cellCensj[cj];
	      double rad_sum = rad_ci + molcellj[cj].getRad();
	      
	      if( CSystem::pbcPos(p_ci - p_cj).normsq() > rad_sum*rad_sum ) continue;
				
	      // if cells overlap, then check each continuent sphere
	      const vector<int> &sphereListI = molcelli[ci].getSphereList();
	      const vector<int> &sphereListJ = molcellj[cj].getSphereList();
				
	      for(int k=0; k<molcelli[ci].getSphereListSize() ; k++)
				{
					int ki = sphereListI[k];
					
					double ri = getKS(ki).getRad();
					CPnt cki = rcen_i + getKS(ki).getCenRot();
					
					for(int h=0; h<sphereListJ.size(); h++)
					{
						int kj = 	sphereListJ[h];
						double allowedOverlap = 0;
						
						if( ri  + mols[j]->getKS(kj).getRad()
							 - CSystem::pbcPos(cki - rcen_j - mols[j]->getKS(kj).getCenRot()).norm()
							 > allowedOverlap )
							
							return true;
						
					}//end kj
					
				}// end ki
	      
	    }//end cj
		}//end ci
		
	}// end j
	
	return false;
}
#endif

// static function
// finds the closest surface point SP[minID] to P
// if dot product = negative -> antiparallel = inside

bool
CMolecule::IsInsideSurface(CPnt P, const vector<CPnt> &SP, const vector<CPnt> &NP)
{
	
  if(SP.size()==0) return false;
	
  double  minsq = DBL_MAX;
  int minID = 0;
  for (int i = 0; i < SP.size(); i++)
	{
		double distsq = (SP[i]- P).normsq();
		if (distsq < minsq)
		{
			minsq = distsq;
			minID = i;
		}
	}
  return ( dot( SP[minID]- P,  NP[minID] )  > 0.0 );
}

bool
CMolecule::OutsideAllBoundaries(const vector<CMolecule*> & mols, CPnt P)
{
  bool bOut = true;
  for(int i=0; i<mols.size() && bOut ; i++)
	{
		CPnt rcen = mols[i]->getRCen();
		if( CSystem::pbcPos(P - rcen).norm() < mols[i]->getMaxR() )
			for(int ki=0; ki<mols[i]->getNKS(); ki++)
			{
				CSolExpCenter *pKi = mols[i]->getpKS(ki);
				if( CSystem::pbcPos(P - rcen - pKi->getCen()).norm() < pKi->getRad() )
	      {
					bOut = false;
					break;
	      }
			}
	}
  return bOut;
}

//member function
bool
CMolecule::IsInsideSurface(CPnt P) const
{
  double  minsq = DBL_MAX;
  int minID = 0;
  for (int i = 0; i < m_SP.size(); i++)
	{
		double distsq = (m_SP[i]- P).normsq();
		if (distsq < minsq)
		{
			minsq = distsq;
			minID = i;
		}
	}
  
  return ( dot( m_SP[minID]- P,  m_NP[minID] )  > 0.0 );
	
}
// calculate maxR from dielectric spheres
double
CMolecule::computeMaxRadius(const vector<CPnt> &cens, const vector<REAL> &radii)
{
  double maxR = 0;
  for(int ki=0; ki<cens.size(); ki++)
	{
		REAL dist = cens[ki].norm() + radii[ki] ;
		if(dist > maxR) maxR = dist;
	}
  return maxR;
}

// calculate maxR from furthestmost charge points
void
CMolecule::computeMaxRadius(const vector<CPnt> &cpos)
{
  m_maxR = 0;
  for(int k=0; k<cpos.size(); k++)
	{
		REAL dist = cpos[k].norm() + SPHERETOL ;
		if(dist > m_maxR) m_maxR = dist;
	}
}

/////////////////////////////////////////////////////////////////////
// full xform is not available
// - have to generate xform for non-pol neigbors
/////////////////////////////////////////////////////////////////////
//this function will compute total interaction energy (S. Liu)
REAL
CMolecule::computeTotalIntEnergy(const vector<CMolecule*> & mols)
{
  int number_of_monomers = mols.size();
  int i_max = (m_bInfinite ? m_unit : number_of_monomers);
	
  REAL pot = 0.0;
  int dummyP = 1;
	
  cout <<"COMPUTING TOTAL INTERACTION ENERGY ..."<<endl;
	
#pragma omp parallel for reduction(+:pot)
  for (int i = 0; i < i_max; i++)
    for (int j = i+1; j < i_max; j++)
    {
      REAL c2cdist = (mols[i]->getRCen() - mols[j]->getRCen()).norm();
      bool bAggI = mols[i]->isAggregateM();
      bool bAggJ = mols[j]->isAggregateM();
			
      cout <<"mols i: "<<i<<" j: "<<j<<"-> ";
			cout <<"interactCenters"<<endl;
			double deltaE = interactCenters_LowMemory(mols[i], mols[j]);
			cout <<"this pair interaction energy is:"<<deltaE<<endl;
			pot += interactCenters_LowMemory(mols[i], mols[j]);
    }
  
  return pot;
  
}


// moli = as a molecule
// molj = as centers
// creates xforms on the fly
REAL
CMolecule::interactMolWithCenters(CMolecule* moli, CMolecule* molj)
{
	
  CXFormBase* xform;
  REAL pot = 0.0;
  int dummyP = 1;
	
  int i = moli->getID();
  int j = molj->getID();
  CPnt rcen_i = moli->getRCen();
  //REAL ri = pKi->getRad();
	
  CLocalExpan Li( CRange(0,1), moli->getMaxR() );
	
  for (int kj=0; kj < molj->getNKS(); kj++)
	{
		const CSolExpCenter *pKj = molj->getpKS(kj);
		REAL rj = pKj->getRad();
		CLocalExpan tL;
		
		CPnt P = CSystem::pbcPos( rcen_i - (molj->getRCen() + pKj->getCen() )  );
		REAL rho = P.norm();
		int pi2j = CXFormBase::computeOrder(rho, moli->getMolTQ(), moli->getMaxR());
		int pj2i = CXFormBase::computeOrder(rho, pKj->getTQH(), pKj->getRad());
		int minp = (pi2j < pj2i ? pi2j : pj2i);
		//int minN = pKj->getSPExSize();
		bool bJ2I;
		CXFormA XA(pKj->getH(), moli->getMolM());
		XA.reset(P, minp);
		bJ2I = minp < pi2j;
		
		if( bJ2I ) // add to moli's local expansion if it's cheaper (transform pKj)
		{
			XA.xformH( pKj->getH(), tL, true);
			Li += tL;
		}
		
		else // otherwise, transform moli
		{
			XA.xformH( moli->getMolM(), tL, false);
			pot += inprod(pKj->getH(), tL);
		}
	}// end-kj
  
  pot += inprod(moli->getMolM(), Li);
  return pot;
}

REAL
CMolecule::interactMols(CMolecule* moli, CMolecule* molj)
{
  REAL pot;
  CLocalExpan tL;
	
  CPnt P = CSystem::pbcPos(moli->getRCen() - molj->getRCen());
  REAL rho = P.norm();
  int pi2j = CXFormBase::computeOrder(rho, moli->getMolTQ(), moli->getMaxR());
  int pj2i = CXFormBase::computeOrder(rho, molj->getMolTQ(), molj->getMaxR());
  int minp = (pi2j < pj2i ? pi2j : pj2i);
	
  CXFormA xform( molj->getMolM(), moli->getMolM() );
  xform.reset(P, minp);
  if(minp < pi2j)
	{
		xform.xformH( molj->getMolM(), tL, true); // if minp < pi2j: transform j2i (i.e. fwd)
		pot = inprod(moli->getMolM(), tL);
	}
  else
	{
		xform.xformH( moli->getMolM(), tL, false); // if minp > pi2j: transform i2j (i.e. reverse)
		pot = inprod(molj->getMolM(), tL);
	}
	
  return pot;
}

REAL
CMolecule::computePotInSpace(const vector<CMolecule*> & mols, CPnt P)
{
  int i_max = (m_bInfinite ? m_unit : N_MOL);
  REAL pot = 0.0;
	
	
#pragma omp parallel for reduction(+:pot) shared (P)
  for (int i = 0; i < i_max; i++)
	{
		CLocalExpan tL;
		CPnt   rcen = mols[i]->getRCen();
		CSpPnt s    = CartToSph( CSystem::pbcPos(P-rcen) );
		REAL   maxR = mols[i]->getMaxR();
		if( s.rho() > maxR && mols[i]->isAggregateM())
		{
			int p = CXFormBase::computeOrder(s.rho(), mols[i]->getMolTQ(), maxR);
			tL = CLExpan(1.0, s, true, p, maxR);
			pot += inprod(mols[i]->getMolM(), tL);
			
		}
		else
		{
			for (int ki=0; ki < mols[i]->getNKS(); ki++)
	    {
	      const CSolExpCenter *pKi = mols[i]->getpKS(ki);
	      CSpPnt s_ki = CartToSph(CSystem::pbcPos(P - rcen - pKi->getCen()) );
	      int p = CXFormBase::computeOrder(s_ki.rho(), pKi->getTQH(), pKi->getRad());
				
	      tL = CLExpan(1.0, s_ki, true, p, pKi->getLScale() );
	      pot += inprod(pKi->getH(), tL);
	    }
			
		}
	}
  
  return pot;
}

const REAL
CMolecule::computePot() const
{
  REAL pot = 0.0;
	
#pragma omp parallel for reduction(+:pot)
  for(int ki=0; ki< getNKS(); ki++)
	{
		const CExpCenter* pK = getpKS(ki);
		pot += pK->computePot();
	}
	
  return pot;
}

//this function will compute force and torque on a specific molecule (S. Liu)
void
CMolecule::computeMol_Force_and_Torque(CPnt &force, CPnt &torque) const
{
	
  force = CPnt();
  torque = CPnt();
	
  if(!getbInterXForm()) return;
	
  for(int ki=0; ki< getNKS(); ki++)
	{
		if( m_k[ki]->IsNoInteractionList()) continue;
		
		CPnt f_ki = m_k[ki]->computeForceOn_0();
		force  += f_ki;
		
		torque += cross (m_k[ki]->getCen(), f_ki); // (1)
		//torque += m_k[ki]->computeTorqueOn_0(m_id); // (2) ->HACK torque dominated by (1)
		
	}
	
  // convert to labframe
  force  = conj(m_orient) * force;
  torque = conj(m_orient) * torque;
  
  return;
}

//this function will compute forces and torques for all molecules (S. Liu)
bool
CMolecule::computeForces(vector<CMolecule*> & mols,
												 vector<CPnt> & force, vector<CPnt> & torque)
{
  if(mols.size() == 1 ) return true;
	
  bool bInterXFS =  CMolecule::generateInterXFormsForPolarize_LowMemory(mols);
  if(! bInterXFS ) {
    cout <<"error in generateInterXFormsForPolarize_LowMemory"<<endl;
    return false;
  }
  CMolecule::polarize_mutual(mols,!m_bGrad, 1000);
  int i_max = (m_bInfinite ? m_unit : N_MOL);
	
  force.clear(); force.resize(i_max);
  torque.clear(); torque.resize(i_max);
  
  for (int i = 0; i < i_max; i++)
	{
		mols[i]->computeMol_Force_and_Torque(force[i],torque[i]);
#if __DEBUGDIE__
		double f = force[i].norm();
		if(fabs(f) > 5 || isnan(f) )
		{
			cout <<"died in computeforces; mol "<<i<<endl;
			CMolecule::writeMolsPQR("died.compforces.pqr", mols);
			CMolecule::saveConfig("died.compforces.config", mols);
			return false; //exit(1);
		}
#endif
		
	}
  
  return true;
}

void
CMolecule::computeMolMultipole()
{
  m_molM.reset(CRange(0,N_POLES));
  m_molM.setScale(1.0/m_maxR);
	
  m_molTQ = 0.0;
  // currently using doing M2M numerically
  for(int ki=0; ki < getNKS(); ki++)
	{
		CSolExpCenter *pKi = m_k[ki];
		CPnt cen = pKi->getCen();
		vector<CPnt> SP = pKi->getSPx();
		vector<REAL> Q  = pKi->getQSolvedH();
		
		for(int k=0; k<pKi->getSPExSize(); k++)
			m_molM += CMExpan(Q[k], CartToSph(SP[k]+cen), m_bKappa, N_POLES, m_maxR);
		
		m_molTQ += pKi->getTQH();
	}
  return;
}


// to worth the effort, only aggregate if molecule has more then 2 spheres,
// and well-separated from at least ONE other molecule
void
CMolecule::aggregateMolMultipoles_Conditional(vector<CMolecule*> & mols)
{
	
  for(int i=0; i<mols.size(); i++)
	{
		mols[i]->setAggregateM(false);
		
		//      if(mols[i]->getNKS() > 2)
		for(int j=0; j<mols.size(); j++)
			if( (mols[i]->getRCen() - mols[j]->getRCen()).norm() > mols[i]->getMaxR() + mols[j]->getMaxR() )
	    {
	      mols[i]->setAggregateM(true);
	      mols[i]->computeMolMultipole();
	      break;
	    }
	}
  
}

// aggregate all
void
CMolecule::aggregateMolMultipoles(vector<CMolecule*> & mols)
{
	
#pragma omp parallel for
  for(int i=0; i<mols.size(); i++)
	{
		mols[i]->setAggregateM(true);
		mols[i]->computeMolMultipole();
	}
  
}


// Compute the maximal distance between a point and a set of points
REAL
CMolecule::maxDist(const CPnt & pnt, const vector<CPnt> & pts)
{
  REAL max = 0.0;
  for (int i = 0; i < pts.size(); i++)
	{
		double distsq = (pts[i]-pnt).normsq();
		if (distsq > max)
			max = distsq;
	}
	
  return sqrt(max);
}

void
CMolecule::rotate(const CQuat & dQ)
{
	
  m_orient *= dQ;
  
  for(int ki=0; ki<getNKS(); ki++)
		m_k[ki]->rotateCenters();
	
#ifdef __CELL__
  for(int ci=0; ci<m_molcells.size(); ci++)
    m_cellCens[ci] =   m_orient * m_molcells[ci].getCen();
#endif
	
}

void
CMolecule::rotateRotCoeff()
{
  m_rot.reset(m_orient, m_p);
}

void
CMolecule::printMolConfig(const CMolecule *mol, char* fname)
{
  ofstream fout(fname);
  CPnt rcen = mol->getRCen();
  for(int ki=0; ki<mol->getNKS(); ki++)
	{
		CPnt p = CSystem::pbcPos(rcen + mol->getKS(ki).getCenRot());
		REAL rad = mol->getKS(ki).getRad();
		fout <<p.x()<<" "<<p.y()<<" "<<p.z()<<" "<<rad<<endl;
	}
}

// contacts using rigid body criterion
void
CMolecule::readMolContactRigid(const char* fname, int &mol1type, vector<CMolContactRigid> &molcontactlist, double addDist)
{
  cout <<"Reading molcontact file "<<fname<<endl;
	
  molcontactlist.clear();
  ifstream fin(fname);
  if (!fin.is_open())
	{
		cout << "Could not open file " << fname << endl;
		exit(0);
	}
	
  int mol2type, ndef;
  fin >> mol1type;
  fin >> ndef;
	
  double rideal,rideal2;
  double  o1Ax, o1Ay, o1Az, o2Ax, o2Ay, o2Az;// with P
  double  o1Bx, o1By, o1Bz, o2Bx, o2By, o2Bz;// z-orient
  
  while (!fin.eof())
	{
		fin >> mol2type; // what the partner molecule should be
		if(fin.eof()) break;
		fin >>rideal
	  >>o1Ax>>o1Ay>>o1Az
	  >>o2Ax>>o2Ay>>o2Az
	  >>o1Bx>>o1By>>o1Bz
	  >>o2Bx>>o2By>>o2Bz;
		rideal2 = (rideal+addDist)*(rideal+addDist);
		
		if(  rideal2 > CMolContactRigid::maxC2CDist2   ) CMolContactRigid::maxC2CDist2 = rideal2;
		molcontactlist.push_back(CMolContactRigid(mol2type, rideal2,
																							CPnt(o1Ax,o1Ay,o1Az), CPnt(o1Bx,o1By,o1Bz),
																							CPnt(o2Ax,o2Ay,o2Az), CPnt(o2Bx,o2By,o2Bz) ) );
	}
	
  cout <<"mol definitions read "<<molcontactlist.size()<<endl;
  fin.close();
}

void
CMolecule::initMolContactRigid(vector<CMolecule*> &mols, int mol1type, vector<CMolContactRigid> &molcontactlist)
{
  for(int i=0; i<mols.size(); i++)
	{
		if(mol1type != mols[i]->getMolType() )
		{
			cout <<"MOLTYPE DIFFERENT!"<<mol1type <<" "<<mols[i]->getMolType()<<endl;
			exit(1);//continue;
		}
		
		mols[i]->setMolContactList( molcontactlist );
	}
}

bool
CMolecule::updateDockStatus(vector<CMolecule*> &mols, bool bDebug, int n)
{
	
  bool bAnyDocked = false;
	
  for(int i=0; i<mols.size(); i++)
    mols[i]->resetDockStatus();
  
  for(int i=0; i<mols.size(); i++)
	{
		for(int j=i+1; j<mols.size(); j++)
		{
			bool bThisDocked = checkDocked(mols[i], mols[j], bDebug)  ;
			if(!bAnyDocked && bThisDocked) bAnyDocked = true;
			
		} //end-j
		
	}//end-i
  return bAnyDocked;
}

bool
CMolecule::checkDocked(CMolecule* moli,  CMolecule* molj, bool bDebug )
{
	
  CPnt P = CSystem::pbcPos(molj->getRCen() - moli->getRCen());
  double c2cdist2 = P.normsq();
	
  if(c2cdist2 >  CMolContactRigid::maxC2CDist2 ) return false;
	
  P.normalize();
	
  for(int ii=0; ii<moli->getNDockSides(); ii++)
  {
		CMolContactRigid molcon = moli->getMolContactRigid(ii);
		//-------------------------
		// needs retool for multiple mol types
		int jj = moli->getReciprocalSide(ii);
		//-------------------------
		
		if(c2cdist2 > molcon.getRideal2() )
		{
			continue;
		}
		
		CQuat Q1 = moli->getOrient();
		CQuat Q2 = molj->getOrient();
		
		if( dot(Q1*molcon.getVec1A(), P)  < CMolContactRigid::angleTolerance ) {
			continue;
		}
		
		if( dot(Q2*molcon.getVec2A(), P)  < CMolContactRigid::angleTolerance ) {
			continue;
		}
		
		if( dot(Q1*molcon.getVec1B(), Q2*molcon.getVec2B() )  < CMolContactRigid::angleTolerance )
		{
			continue;
		}
		
		moli->setBDocked(ii, true);
		molj->setBDocked(jj, true);
		
		moli->addDockedNeigh(molj->getID());
		molj->addDockedNeigh(moli->getID());
		
		if(bDebug)
			cout <<"DOCKED: "<< moli->getID() <<" "<<molj->getID()<<" "<<" mode: "<<ii<<endl;
		
		return true;
		
	}//end-ii
  
  return false;
  
}

bool
CMolecule::checkDocked_debug(CMolecule* moli,  CMolecule* molj, bool bDebug )
{
	
  CPnt P = CSystem::pbcPos(molj->getRCen() - moli->getRCen());
  double c2cdist2 = P.normsq();
	
  if(c2cdist2 >  CMolContactRigid::maxC2CDist2 ) return false;
	
  P.normalize();
	
  for(int ii=0; ii<moli->getNDockSides(); ii++)
  {
		CMolContactRigid molcon = moli->getMolContactRigid(ii);
		
		//if( molcon.getMol2Type() != molj->getMolType()) continue;
		//if( moli->getBDocked(ii) ) continue;
		
		//-------------------------
		// needs retool for multiple mol types
		int jj = moli->getReciprocalSide(ii);
		//if( molj->getBDocked(jj) ) continue;
		//-------------------------
		
		bool bc2cdist = c2cdist2 > molcon.getRideal2();
		
		CQuat Q1 = moli->getOrient();
		CQuat Q2 = molj->getOrient();
		bool b1AP = dot(Q1*molcon.getVec1A(), P) < CMolContactRigid::angleTolerance;
		bool b2AP = dot(Q2*molcon.getVec2A(), P) < CMolContactRigid::angleTolerance;
		bool b1B2B = dot(Q1*molcon.getVec1B(), Q2*molcon.getVec2B() )  < CMolContactRigid::angleTolerance;
		
		if(bc2cdist)
			cout << "1AP:"<<b1AP<<" 2AP:"<<b2AP<<" 1B2B:" <<b1B2B<<endl;
		
		if( !bc2cdist || !b1AP || !b2AP || !b1B2B)
			continue;
		
		moli->setBDocked(ii, true);
		molj->setBDocked(jj, true);
		
		moli->addDockedNeigh(molj->getID());
		molj->addDockedNeigh(moli->getID());
		
		if(bDebug)
			cout <<"DOCKED: "<< moli->getID() <<" "<<molj->getID()<<" "<<" mode: "<<ii<<endl;
		
		return true;
	}//end-ii
  return false;
}

void
CMolecule::updateDockStatistics(const vector<CMolecule*> &mols, char* runname, int n, double t,
																vector<ofstream*> &douts, int nSpeciesPerBin)

{
  const int nmol = mols.size();
  const int nbin = nmol;
	
  char fname[MAX_CHARNUM];
  vector<bool> bCounted(nmol, false);
  vector<int> histogram(nbin,0);
  bool bError = false;
	
  for(int i=0; i<nmol; i++)
	{
		int size = 0;
		updateSpecies(mols, i, bCounted, size);
		if(size == 0 ) continue; // already counted
		
		assert(size <= nbin);
		
		histogram[ size-1 ]++;
		
	}
	
  //------------------------------------------
  // prints out molecules that have docked
  vector<CMolecule*> dockedMols;
  for(int i=0; i<nmol; i++)
	{
		if(mols[i]->getDockedNeigh().size() > 0)
		{
			dockedMols.push_back( mols[i] );
		}
	}
	
  if(dockedMols.size() > 0)
	{
		sprintf(fname, "config.%s.%d.pqr", runname,n);
		writeMolsPQR(fname, dockedMols);
	}// endif
	
  //------------------------------------------
	
}
void
CMolecule::updateSpecies(const vector<CMolecule*> &mols, int j, vector<bool> &bCounted, int & size)
{
  if(bCounted[j]) return;
	
  bCounted[j] = true;
  size++;
  
  vector<int> dockedNeigh = mols[j]->getDockedNeigh();
  for(int n=0; n<dockedNeigh.size(); n++)
	{
		int m = dockedNeigh[n];
		// assert(m < mols.size());
		updateSpecies(mols, m, bCounted, size);
	}
  
  return;
}

//// member functions ////
// reset at every timestep (if we don't have any retaining docking forces)
void
CMolecule::resetDockStatus()
{
  m_dockedNeigh.clear();
  m_bDocked.assign(getNDockSides(), false);
}

