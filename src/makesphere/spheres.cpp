#include <vector>
#include <cfloat>
#include <fstream>
#include "spheres.h"
#define NTRIALS 9
#define MAX_NTRIALS 40
#define MCTRIALS 2000
#define NUM_POINTS 25001
#define NMAX 50000

//!  The getInterfaceAtoms function
/*! The function that uses the molecule contact information
to determine which atoms lay at the contact interface. Outputs 
the positions ( in labframe ) of atoms in contact, ignoring 
all hydrogens. 
		\param pos1 a vector of coordinates for all atoms in the molecule
		\param pos2 a vector of coordinates for atoms in the second molecule
		\param rad1 a vector of atom radii for all atoms in the molecule
		\param rad2 a vector of atom radii for atoms in the second molecule
		\param conpos1 a vector of atoms that are in contact with the 2nd mol */
void getInterfaceAtoms( const vector <CPnt> &pos1,const vector <CPnt> &pos2,
											 const vector <REAL> &rad1,const vector <REAL> &rad2,
											 vector <CPnt> & conpos1, double cutoff )
{	
	// find contact atoms for Protein 1 that contacts at 
	// least one heavy atom in protein 2
	conpos1.clear(  );
	for ( int i=0; i < pos1.size( ); i++)
	{
		if (  rad1[i] < 1.5 ) continue; // ignore hydrogen atoms
		
		for ( int j=0; j < pos2.size( ); j++)
		{
			if (  rad2[j] < 1.5 )  continue; // ignore hydrogen atoms
			
			if ( (pos1[i]  - pos2[j] ).norm() <= cutoff)
			{
				conpos1.push_back( pos1[i] );
				break;
			}
		}
	}	
	return;
}

//!  The classifyVertPoints function
/*! The function that uses the molecule contact information
to determine which MSMS points are at the molecule interface
and which are not.
		\param interfaceAtomPos a vector of coordinates for molecule atoms at interface
		\param vertPts a vector of all vertex points
		\param IPind a vector of indices of all interface points
		\param atomVertCutoff a cutoff for all vertex points */
void classifyVertPoints(const vector<CPnt> & interfaceAtomPos,const vector<CPnt> & vertPts,
												vector<int> & IPind, double atomVertCutoff)
{
  IPind.clear();
  vector<CPnt>::const_iterator it2;
  for (int i=0; i< vertPts.size(); i++)
	{
		it2 = interfaceAtomPos.begin();
		for (; it2 != interfaceAtomPos.end(); it2++)
		{
			if ( (vertPts[i] - *it2).norm() <= atomVertCutoff )
	    {
	      IPind.push_back( i );
	      break;
	    }
		}
	}
  cout <<"No. of Interface Vert Points: "<<IPind.size()<< " No. of NonInterface Vert Points: "<<vertPts.size() - IPind.size()<<endl;
}

//!  The findCenters function
/*! The function that chooses centers for charges in the molecule, 
for a system with contacts.
		\param P a vector of cartesian coordinates of all atoms in the molecule
		\param arad a vector of floating points of all atom radii in the molecule
		\param SP a vector of coordinates of all MSMS surface points in the mol 
		\param NP a vector of all MSMS surface points in the molecule 
		\param cens a vector of centers for CG spheres
		\param R a vector of CG sphere radii
		\param tolSP the tolerance of how much the sphere can protrude from SES  */
void findCenters(const vector<CPnt> & P, const vector<double> & arad,
								 const vector<CPnt> & SP, const vector<CPnt> & NP,
								 vector<CPnt> & cens, vector<double> & R, double tolSP)
{
  cens.clear();
  int np = P.size();
  // initialize ind
  vector<int> ind(np);
  for (int i = 0; i < np; i++) ind[i] = i;
	
  int N, N1;
  vector<int> ind1, ind2;
  ind2.reserve(np);
  CPnt cn;
  double r;
  int j =-1;
  int ct =0;
  // Find centers to encompass charges
  while(ind.size() !=0 && ct < np)
	{
		j++;
		cens.resize(j+1);
		R.resize(j+1);
		int Nmax = 0;
		cout<<"Finding center no. "<<j<<":"<<endl;
		// Try at least 10 times to find the best choice for next center
		int m=-1;
		while (m < NTRIALS || ( m>= NTRIALS && Nmax==0 && m < MAX_NTRIALS) )
		{
			m++;
			middle(P, arad, SP, NP, cn, r, N, tolSP, ind, ind2);
			cout<<"trial "<<m<<": center = "<<cn<<" bound points :"<<ind2.size()<<" "<<N<<endl;
			if (N > Nmax)
	    {
	      cens[j] = cn;
	      R[j] = sqrt(r);
	      ind1 = ind2;
	      Nmax = N;
	    }
		}
		if(Nmax == 0) // pick a center from ind
		{
			ind1.clear(); ind1.push_back( ind[(int) floor(drand48()* ind.size()  )]);
			Nmax = 1;
		}
		
		if( Nmax == 1)
		{
			cens[j] = P[ind1[0]]; // set sphere center as position of bounded atom
			R[j] = arad[ind1[0]]; // set sphere radius as vdw of bounded atom
		}
		
		// Remove the indices assigned to new center from list of unassigned points
		subtract(ind1, ind);
		cout<<"  Chosen center : "<<cens[j]<<" radius : "<<R[j]<<endl;
		cout<<"  "<<j<<") points bounded : "<<ind1.size()<<"  "<<"points remaining : "<<ind.size()<<endl;
		ct++;
	}//end while
	
  cout<<"Total Centers = "<<cens.size()<<endl;
  if(ind.size() > 0)
	{
		cout <<"Remaining Centers not assigned: "<<endl;
		for(int i=0; i<ind.size(); i++)
		{
			CPnt p = P[ind[i]];
			cout <<p.x()<<" "<<p.y()<<" "<<p.z()<<" "<<arad[ind[i]]<<endl;
		}
		
	}
  return;
}

//!  The findCenters_interface function
/*! The function that chooses centers for charges in the molecule, 
for a system with contacts.
		\param P a vector of cartesian coordinates of all atoms in the molecule
		\param arad a vector of floating points of all atom radii in the molecule
		\param SP a vector of coordinates of all MSMS surface points in the mol 
		\param NP a vector of all MSMS surface points in the molecule 
		\param IPind a vector of indices for each MSMS point at the interface
		\param cens a vector of centers for CG spheres
		\param R a vector of CG sphere radii
		\param tolSP the tolerance of how much the sphere can protrude from SES
		\param tolIP the interface tolerance only for contact situation  */
void findCenters_interface(const vector<CPnt> & P, const vector<double> & arad,
													 const vector<CPnt> & SP, const vector<CPnt> & NP,
													 const vector<int> & IPind,
													 vector<CPnt> & cens, vector<double> & R,
													 double tolSP, double tolIP)
{
  cens.clear();
  int np = P.size();
  // initialize ind
  vector<int> ind(np);
  for (int i = 0; i < np; i++) ind[i] = i;
	
  int N, N1;
  vector<int> ind1, ind2;
  ind2.reserve(np);
  CPnt cn;
  double r;
  int j =-1;
  int ct =0;
  // Find centers to encompass charges
  while(ind.size() !=0 && ct < np)
	{
		j++;
		cens.resize(j+1);
		R.resize(j+1);
		int Nmax = 0;
		cout<<"Finding center no. "<<j<<":"<<endl;
		// Try at least 10 times to find the best choice for next center
		int m=-1;
		while (m < NTRIALS || ( m>= NTRIALS && Nmax==0 && m < MAX_NTRIALS) )
		{
			m++;
			middle_interface(P, arad, SP, NP, IPind, cn, r, N, tolSP, tolIP, ind, ind2);
			cout<<"trial "<<m<<": center = "<<cn<<" bound points :"<<ind2.size()<<" "<<N<<endl;
			if (N > Nmax)
	    {
	      cens[j] = cn;
	      R[j] = sqrt(r);
	      ind1 = ind2;
	      Nmax = N;
	    }
		}
		
		if(Nmax == 0) // pick a center from ind
		{
			ind1.clear(); ind1.push_back( ind[(int) floor(drand48()* ind.size()  )]);
			Nmax = 1;
		}
		
		if( Nmax == 1)
		{
			cens[j] = P[ind1[0]]; // set sphere center as position of bounded atom
			R[j] = arad[ind1[0]]; // set sphere radius as vdw of bounded atom
		}
		// Remove the indices assigned to new center from list of unassigned points
		subtract(ind1, ind);
		cout<<"  Chosen center : "<<cens[j]<<" radius : "<<R[j]<<endl;
		cout<<"  "<<j<<") points bounded : "<<ind1.size()<<"  "<<"points remaining : "<<ind.size()<<endl;
		ct++;
	}//end while
	
  cout<<"Total Centers = "<<cens.size()<<endl;
  if(ind.size() > 0)
	{
		cout <<"Remaining Centers not assigned: "<<endl;
		for(int i=0; i<ind.size(); i++)
		{
			CPnt p = P[ind[i]];
			cout <<p.x()<<" "<<p.y()<<" "<<p.z()<<" "<<arad[ind[i]]<<endl;
		}
		
	}
  return;
}


//!  The middle function
/*! The function that finds the center and radius of the ball that 
contains the most charge points and is bounded within surface S with tolerance. 
Use an MC-SA approach. tight tol for interface, looser otherwise
 Tries to maximize points listed in indmax, bounded points will be returned in ind2
		\param P a vector of cartesian coordinates of all atoms in the molecule
		\param arad a vector of floating points of all atom radii in the molecule
		\param SP a vector of coordinates of all MSMS surface points in the mol 
		\param NP a vector of all MSMS surface points in the molecule 
		\param cens a CPnt object for center location
		\param R a floating point for center radius
		\param N an integer for the CG sphere index
		\param tolSP the tolerance of how much the sphere can protrude from SES
		\param indmax a vector of indices for all unbounded atoms
		\param ind2 a vector of indices bounded by created CG sphere*/
void middle(const vector<CPnt> & P, const vector<double> & arad, const vector<CPnt> & SP, const vector<CPnt> & NP,
						CPnt & cen, double & r, int & N, double tolSP,
						const vector<int> & indmax, vector<int> & ind2)
{
  float beta = 2.0;
  int sz = indmax.size();
  const REAL sphereTol =  1.0;
	
  // initialize cen, N, r
  cen = P[ indmax[(int) floor(drand48()* sz  )]];// pick a center from indmax
  N = 0;
  r = 0;
  // monte carlo search of new center
  for (int n = 0; n < 1200; n++)
	{
		CPnt tcen; // trial center
		REAL tmin; // trial radius^2
		int tN = 0; // trial N
		
		// find the current minimum distance between current cen and surface to scale rand vector
		REAL cmin = DBL_MAX;
		for (int i = 0; i < SP.size(); i++)
		{
			double distsq = (SP[i]-cen).normsq();
			if (distsq < cmin)   cmin = distsq;
		}
		REAL scale = sqrt(cmin);//( sqrt(cmin) < 1.0? sqrt(cmin) : 1.0);
		tcen = cen + ( scale * normRand() * randOrient());
    
		// Find the minimal distance from trial center to surface
		int minID = 0;
		tmin = DBL_MAX;
		for (int i = 0; i < SP.size(); i++)
		{
			double distsq = (SP[i]-tcen).normsq();
			if (distsq < tmin)
	    {
	      tmin = distsq;
	      minID = i;
	    }
		}
		// reject trial if tcen is outside boundary
		if( dot( SP[minID]-tcen,  NP[minID] )  < 0.0 ) continue;
		// extend radius by tolSP
		tmin = sqrt(tmin) + tolSP;
		tmin *= tmin;
		// Count the number of indmax points that fall within this ball
		for (int i = 0; i < sz; i++)
		{
			double dist = (P[indmax[i]]-tcen).norm() + arad[indmax[i]]; 
			if (dist < sqrt(tmin))
				tN++;
		}
		// Apply MC criterion to the change in number of bounded charges
		if (expf(beta*(tN - N)) > drand48())
		{
			cen = tcen;
			N = tN;
			r = tmin;
		}
		
		if ((n+1) % 100 == 0)
			beta *= 1.1;
	}//end n-loop
	
  // Create a list of all points contained inside the bounding sphere
  ind2.resize(0);
  for (int i = 0; i < sz; i++)
	{
		double dist = (P[indmax[i]]-cen).norm() + arad[indmax[i]];//sphereTol;
		//double d = arad[indmax[i]] > sphereTol ? sphereTol : arad[indmax[i]];
		//      double dist = (P[indmax[i]]-cen).norm() + d;
		if (  dist  < sqrt(r))
			ind2.push_back(indmax[i]);
	}
}


//!  The middle_interface function
/*! The function that finds the center and radius of the ball that 
contains the most charge points and is bounded within surface S with tolerance. 
Use an MC-SA approach. tight tol for interface, looser otherwise
 Tries to maximize points listed in indmax, bounded points will be returned in ind2
		\param P a vector of cartesian coordinates of all atoms in the molecule
		\param arad a vector of floating points of all atom radii in the molecule
		\param SP a vector of coordinates of all MSMS surface points in the mol 
		\param NP a vector of all MSMS surface points in the molecule 
		\param IPind a vector of indices for each MSMS point at the interface
		\param cens a CPnt object for center location
		\param R a floating point for center radius
		\param N an integer for the CG sphere index
		\param tolSP the tolerance of how much the sphere can protrude from SES
		\param tolIP the interface tolerance only for contact situation
		\param indmax a vector of indices for all unbounded atoms
		\param ind2 a vector of indices bounded by created CG sphere*/
void middle_interface(const vector<CPnt> & P, const vector<double> & arad,
											const vector<CPnt> & SP, const vector<CPnt> & NP,
											const vector<int> & IPind,
											CPnt & cen, double & r, int & N, double tolSP, double tolIP,
											const vector<int> & indmax, vector<int> & ind2)
{
  float beta = 2.0;
  int sz = indmax.size();
	
  const REAL sphereTol =  1.0;
	
  // initialize cen, N, r
  cen = P[ indmax[(int) floor(drand48()* sz  )]];// pick a center from indmax
  N = 0;
  r = 0;
	
  // monte carlo search of new center
  for (int n = 0; n < 1200; n++)
	{
		CPnt tcen; // trial center
		REAL tmin; // trial radius^2
		int tN = 0; // trial N
		
		// find the current minimum distance between current cen and surface to scale rand vector
		REAL cmin = DBL_MAX;
		for (int i = 0; i < SP.size(); i++)
		{
			double distsq = (SP[i]-cen).normsq();
			if (distsq < cmin)   cmin = distsq;
		}
		
		REAL scale = sqrt(cmin);//( sqrt(cmin) < 1.0? sqrt(cmin) : 1.0);
		tcen = cen + ( scale * normRand() * randOrient());
    
		// Find the minimal distance from trial center to surface for SP points
		int minID = 0;
		REAL tmin_s = DBL_MAX;
		
		for (int i = 0; i < SP.size(); i++)
		{
			double distsq = (SP[i]-tcen).normsq();
			if (distsq < tmin_s)
	    {
	      tmin_s = distsq;
	      minID = i;
	    }
		}
		
		// reject trial if tcen is outside boundary
		if( dot( SP[minID]-tcen,  NP[minID] )  < 0.0 ) continue;
		
		// extend radius by tolSP
		tmin_s = sqrt(tmin_s) + tolSP;
		tmin_s *= tmin_s;
		
		// Find the minimal distance from trial center to surface for SP points
		REAL tmin_i = DBL_MAX;
		for (int i = 0; i < IPind.size(); i++)
		{
			int k=IPind[i];
			double distsq = (SP[k]-tcen).normsq();
			if (distsq < tmin_i)	      tmin_i = distsq;
		}
		// extend radius by tolSP
		tmin_i = sqrt(tmin_i) + tolIP;
		tmin_i *= tmin_i;
		
		tmin = (tmin_i < tmin_s ? tmin_i : tmin_s);
		
		// Count the number of indmax points that fall within this ball
		for (int i = 0; i < sz; i++)
		{
			double dist = (P[indmax[i]]-tcen).norm() + arad[indmax[i]]; //sphereTol;
			if (dist < sqrt(tmin))
				tN++;
		}		
		// Apply MC criterion to the change in number of bounded charges
		if (expf(beta*(tN - N)) > drand48())
		{
			cen = tcen;
			N = tN;
			r = tmin;
		}
		
		if ((n+1) % 100 == 0)
			beta *= 1.1;
	}//end n-loop
	
  // Create a list of all points contained inside the bounding sphere
  ind2.resize(0);
  for (int i = 0; i < sz; i++)
	{
		double dist = (P[indmax[i]]-cen).norm() + arad[indmax[i]];//sphereTol;
		//double d = arad[indmax[i]] > sphereTol ? sphereTol : arad[indmax[i]];
		//      double dist = (P[indmax[i]]-cen).norm() + d;
		if (  dist  < sqrt(r))
			ind2.push_back(indmax[i]);
	}
}

void middle_pt(const vector<CPnt> & P, const vector<CPnt> & SP, const vector<CPnt> & NP,
							 CPnt & cen, double & r, int & N, double tolSP,
							 const vector<int> & indmax, vector<int> & ind2)
{
	
  cout <<"calling middle pt"<<endl;
  float beta = 2.0;
  int sz = indmax.size();
	
  // pick a center from indmax
	cen = P[ indmax[(int) floor(drand48()* sz  )]];
	
  // monte carlo search of new center
  N = 0;
  r = 0;
	
  // find radius
  REAL minSP;
  minSP = DBL_MAX;
  for (int i = 0; i < SP.size(); i++)
	{
		double distsq = (SP[i]-cen).normsq();
		if (distsq < minSP) minSP = distsq;
	}
	
  minSP = (sqrt(minSP) + tolSP);
  r = minSP *minSP;
  
  // Create a list of all points contained inside the bounding sphere
  ind2.resize(0);
  for (int i = 0; i < indmax.size(); i++)
	{
		double distsq = (P[indmax[i]]-cen).normsq();
		if (distsq < r)
			ind2.push_back(indmax[i]);
	}
  
  N = ind2.size();
  
}


//!  The subtract function
/*! The function that removes all objects in A from the B list
		\param A a vector of ints to remove
		\param B a vector of ints to remove from  */
void subtract(const vector<int> & A, vector<int> & B)
{
  if (A.size() == 0 || B.size() == 0)
    return;
	
  vector<int>::const_iterator ait = A.begin();
  vector<int>::iterator bit = B.begin();
  int ib = 0;
  int ia = 0;
  //  while (bit != B.end())
  while (bit != B.end() && ait != A.end())
	{
		if (*ait < *bit)
			ait++;
		else if (*ait > *bit)
			bit++;
		else
			bit = B.erase(bit);
	}
}

bool IsInsideSurface(CPnt P, const vector<CPnt> & SP, const vector<CPnt> & NP, double & minsq)
{
  minsq = DBL_MAX;
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

//!  The IsInsideSurface function
/*! The function that determines whether a sphere is 
exposed 
		\param P a cartesian coordinate object of the CG sphere point
		\param SP a vector of coordinates of all MSMS surface points in the mol 
		\param NP a vector of all MSMS surface points in the molecule 
		\return bExposed a boolean of whether or not the sphere is inside 
				the surface, is true for yes, false for no.  */
bool IsInsideSurface(CPnt P, const vector<CPnt> & SP, const vector<CPnt> & NP)
{
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

//!  The spherePts function
/*! The function that generates a set of N points uniformly 
distributed on a unit sphere.
		\param N an integer of number of points on surf of the sphere
		\param th a vector of theta angles for all sphere surface points
		\param ph a vector of phi angles for all sphere surface points */
void
spherePts(int N, vector<REAL> & th, vector<REAL> & ph)
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


// generate surface points for sphere i (wrt to center i)
// assigned to 'exposed' or 'buried' based on neighboring spheres
void
getSpherePoints(int ki, const vector<int> &neigh, const vector<CPnt> &cens,
								const vector<double> &radii, const vector<double> &rad2,
								const vector<CPnt> & SP,  const vector<CPnt> & NP,
								vector<CPnt> &SPE, vector<CPnt> &SPB)
{
  SPE.clear();
  SPB.clear();
  // generate sphere points
  vector<REAL> th, ph;
  int N = (int)ceil(radii[ki]*radii[ki]*NUM_POINTS);
  spherePts(N, th, ph);
  
  double minsq; // dummy
  for (int n = 0; n < th.size(); n++)
	{
		CPnt q = SphToCart(CSpPnt(1, th[n], ph[n]));
		CPnt p = radii[ki]*q;
		
		// Do not choose points located very close to Z-axis (singularity issues)
		//   if (q.x()*q.x() + q.y()*q.y() < 1e-10) continue;
		// First check whether point is outside neighboring spheres
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
		
		// Then check that point is outside boundary
		if (k == neigh.size() )
		{
			if( !IsInsideSurface( p + cens[ki], SP, NP, minsq) )
				SPE.push_back(p);
			else SPB.push_back(p);
		}
	}
  return;
}

//!  The findNeighbors function
/*! The function that determines the neighboring CG 
spheres for sphere ki
		\param ki an integer index of the sphere of interest
		\param cens a vector of coordinates for all CG spheres
		\param radii a vector of radii for all CG spheres
		\param neigh a vector of neighbor sphere indices for ki  */
void
findNeighbors(int ki, const vector<CPnt> &cens,
							const vector<double> &radii, vector<int> &neigh)
{
  neigh.clear();
	
  for (int kj = 0; kj < cens.size(); kj++)
	{
		if( ki!=kj )
		{
			double sumRad = radii[ki]+radii[kj];
			if( (cens[ki]-cens[kj]).normsq() < sumRad * sumRad )
				neigh.push_back(kj);
		}
	}
	
}

//!  The IsSphereExposed function
/*! The function that determines whether a sphere is 
exposed 
		\param ki an integer index of the sphere of interest
		\param neigh a vector of neighbor sphere indices for ki 
		\param cens a vector of coordinates for all CG spheres
		\param radii a vector of radii for all CG spheres
		\param rad2 a vector of the radii^2 of each CG sphere
		\param SP a vector of coordinates of all MSMS surface points in the mol 
		\param NP a vector of all MSMS surface points in the molecule 
		\return bExposed a boolean of whether or not the sphere is exposed
				is true for yes, false for no.  */
bool IsSphereExposed(int ki, const vector<int> &neigh, const vector<CPnt> &cens,
										 const vector<double> &radii, const vector<double> &rad2,
										 const vector<CPnt> & SP,  const vector<CPnt> & NP)
{
  // generate sphere points
  vector<REAL> th, ph;
  int N = (int)ceil(radii[ki]*radii[ki]*NUM_POINTS);
  if(N > NMAX) N = NMAX;
  spherePts(N, th, ph);
  bool bExposed = false;
  double minsq; // dummy
  for (int n = 0; n < th.size(); n++)
	{
		CPnt q = SphToCart(CSpPnt(1, th[n], ph[n]));
		CPnt p = radii[ki]*q;
		
		// First check whether point is outside neighboring spheres
		int k = 0;
		for (k = 0; k < neigh.size(); k++)
		{
			int kj = neigh[k];
			double distsq = (p + cens[ki] - cens[kj]).normsq();
			if ( distsq < rad2[kj]) break;
		}
		
		// If not bound by any sphere, then check that point is also outside boundary
		if (k == neigh.size() )
			if (!IsInsideSurface( p+ cens[ki], SP, NP) )
	    {
	      bExposed = true;
	      break;
	    }
	}
  return bExposed;
}

// check if any sphere point is exposed
bool IsSphereExposedPercent(int ki, const vector<int> &neigh, const vector<CPnt> &cens,
														const vector<double> &radii, const vector<double> &rad2,
														const vector<CPnt> & SP,  const vector<CPnt> & NP)
{
  const double expRatio = 0.05;
  // generate sphere points
  vector<REAL> th, ph;
  int N = (int)ceil(radii[ki]*radii[ki]*NUM_POINTS);
  if(N > NMAX) N = NMAX;
	
  spherePts(N, th, ph);
  bool bExposed = false;
  int nExposed = 0;
  double minsq; // dummy
  for (int n = 0; n < th.size(); n++)
	{
		CPnt q = SphToCart(CSpPnt(1, th[n], ph[n]));
		CPnt p = radii[ki]*q;
		
		// First check whether point is outside neighboring spheres
		int k = 0;
		for (k = 0; k < neigh.size(); k++)
		{
			int kj = neigh[k];
			double distsq = (p + cens[ki] - cens[kj]).normsq();
			
			if ( distsq < rad2[kj]) break;
		}
		
		// If not bound by any sphere, then check that point is also outside boundary
		if (k == neigh.size() )
			if (!IsInsideSurface( p+ cens[ki], SP, NP) )
	    {
	      nExposed ++;
	    }
	}
  return ( nExposed / double(N) > expRatio );
}

