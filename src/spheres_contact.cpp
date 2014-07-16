#include <vector>
#include <cfloat>
#include <fstream>
#include "spheres.h"

#define NTRIALS 9						//!< The number of trials to attempt for sphere generation
#define MAX_NTRIALS 40			//!< An upper limit to number of trials to attempt for sphere generation

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
void classifyVertPoints( const vector<CPnt> & interfaceAtomPos,
												const vector<CPnt> & vertPts, 
												vector<int> & IPind, double atomVertCutoff )
{
	IPind.clear(  );
	
	vector<CPnt>::const_iterator it2;
	for ( int i=0; i< vertPts.size( ); i++)
	{
		it2 = interfaceAtomPos.begin(  );
		for ( ; it2 != interfaceAtomPos.end( ); it2++)
		{
			if (  (vertPts[i] - *it2 ).norm() <= atomVertCutoff )
			{ 
				IPind.push_back(  i  );
				break;
			}	  
		}
	}
	
	cout <<"No. of Interface Vert Points: "<<IPind.size(  )<< 
	" No. of NonInterface Vert Points: "<<vertPts.size(  ) - IPind.size()<<endl;
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
void findCenters_interface( const vector<CPnt> & P, const vector<double> & arad,
													 const vector<CPnt> & SP, const vector<CPnt> & NP, 
													 const vector<int> & IPind, 
													 vector<CPnt> & cens, vector<double> & R, 
													 double tolSP, double tolIP )
{
	cens.clear(  ); 	
	int np = P.size(  );		// number of atoms in the system
	
	// initialize ind with indices for each of the atoms in the system
	vector<int> ind( np );
	for ( int i = 0; i < np; i++ ) ind[i] = i;
	
	int N, N1;
	vector<int> ind1, ind2;
	ind2.reserve( np );
	CPnt cn;
	double r;
	int j =-1;				// Index of CG centers
	int ct =0;				// count of 
	
	// Find centers to encompass charges, while the size of
	// ind is not zero and the count is less than number of atoms in sys
	while( ind.size( ) !=0 && ct < np)		
	{
		j++;
		cens.resize( j+1 );
		R.resize( j+1 );
		int Nmax = 0;
		cout<<"Finding center no. "<<j<<":"<<endl;
		
		// Try at least 10 times to find the best choice for next center
		// And continue if 
		int m=-1;
		while ( m < NTRIALS || ( m>= NTRIALS && Nmax==0 && m < MAX_NTRIALS ) )
		{
			m++;
			middle_interface( P, arad, SP, NP, IPind, cn, r, N, tolSP, tolIP, ind, ind2 );
			cout<<"trial "<<m<<": center = "<<cn<<" bound points :"<<ind2.size(  )<<" "<<N<<endl; 
			if ( N > Nmax )
			{
				cens[j] = cn;
				R[j] = sqrt( r );
				ind1 = ind2;
				Nmax = N;
			}
		}
		
		if( Nmax == 0 ) // If no spheres were bound pick a CG sphere from ind
		{
			ind1.clear(  ); ind1.push_back( ind[(int) floor(drand48()* ind.size()  )]); 
			Nmax = 1; 
		}
		
		if(  Nmax == 1 ) 
		{
			cens[j] = P[ind1[0]]; // set sphere center as position of bounded atom
			R[j] = arad[ind1[0]]; // set sphere radius as vdw of bounded atom
		}
		
		// Remove the indices assigned to new center from list of unassigned points
		subtract( ind1, ind );
		
		cout<<"  Chosen center : "<<cens[j]<<" radius : "<<R[j]<<endl;
		cout<<"  "<<j<<" ) points bounded : "<<ind1.size( )<<"  "
		<<"points remaining : "<<ind.size()<<endl;
		
		ct++;
		
	}//end while
	
	cout<<"Total Centers = "<<cens.size(  )<<endl;
	
	if( ind.size( ) > 0) // Print out unassigned centers
	{ 
		cout <<"Remaining Centers not assigned: "<<endl;
		for( int i=0; i<ind.size( ); i++)
		{
			CPnt p = P[ind[i]];
			cout <<p.x(  )<<" "<<p.y()<<" "<<p.z()<<" "<<arad[ind[i]]<<endl;
		}
	}  
	return;
} // end findCenters_interface

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
void middle_interface( const vector<CPnt> & P, const vector<double> & arad, 
											const vector<CPnt> & SP, const vector<CPnt> & NP,
											const vector<int> & IPind, 
											CPnt & cen, double & r, int & N, double tolSP, double tolIP,
											const vector<int> & indmax, vector<int> & ind2 )
{
	float beta = 2.0;
	int sz = indmax.size(  );	
	const REAL sphereTol =  1.0;
	
	// initialize cen, N, r 
	cen = P[ indmax[( int ) floor(drand48()* sz  )]];// pick a center from indmax
	N = 0;
	r = 0;
	
	// monte carlo search of new center
	for ( int n = 0; n < 1200; n++ )
	{
		CPnt tcen; // trial center 
		REAL tmin; // trial radius^2
		int tN = 0; // trial N
		
		// find the current minimum distance between 
		// current cen and surface to scale rand vector
		REAL cmin = DBL_MAX;
		for ( int i = 0; i < SP.size( ); i++)
		{
			double distsq = ( SP[i]-cen ).normsq();
			if ( distsq < cmin )   cmin = distsq;
		}
		
		REAL scale = sqrt( cmin );
		tcen = cen + (  scale * normRand( ) * randOrient());
		
		// Find the minimal distance from trial center to surface for SP points
		int minID = 0;
		REAL tmin_s = DBL_MAX;
		
		for ( int i = 0; i < SP.size( ); i++)
		{
			double distsq = ( SP[i]-tcen ).normsq();
			if ( distsq < tmin_s )
			{
				tmin_s = distsq;
				minID = i;
			}
		}
		
		// reject trial if tcen is outside boundary
		if(  dot( SP[minID]-tcen,  NP[minID]  )  < 0.0 ) continue; 
		
		// extend radius by tolSP
		tmin_s = sqrt( tmin_s ) + tolSP;
		tmin_s *= tmin_s;
		
		// Find the minimal distance from trial center to surface for interface points
		REAL tmin_i = DBL_MAX;
		for ( int i = 0; i < IPind.size( ); i++)
		{
			int k=IPind[i];
			double distsq = ( SP[k]-tcen ).normsq();
			if ( distsq < tmin_i )	      tmin_i = distsq;
		}
		// extend radius by tolIP
		tmin_i = sqrt( tmin_i ) + tolIP;
		tmin_i *= tmin_i;
		
		tmin = ( tmin_i < tmin_s ? tmin_i : tmin_s );
		
		// Count the number of indmax points that fall within this ball
		for ( int i = 0; i < sz; i++ )
		{
			double dist = ( P[indmax[i]]-tcen ).norm() + arad[indmax[i]]; //sphereTol;
			if ( dist < sqrt(tmin ))
				tN++;
		}
		
		// Apply MC criterion to the change in number of bounded charges
		if ( expf(beta*(tN - N )) > drand48())
		{
			cen = tcen;
			N = tN;
			r = tmin;
		}
		
		if ( (n+1 ) % 100 == 0)
			beta *= 1.1;
		
	}//end n-loop
	
	// Create a list of all points contained inside the bounding sphere
	ind2.resize( 0 );  
	for ( int i = 0; i < sz; i++ )
	{
		double dist = ( P[indmax[i]]-cen ).norm() + arad[indmax[i]];//sphereTol;
		//double d = arad[indmax[i]] > sphereTol ? sphereTol : arad[indmax[i]];
		//      double dist = ( P[indmax[i]]-cen ).norm() + d;
		if (   dist  < sqrt(r ))
			ind2.push_back( indmax[i] );
	}
} // end middle_interface
