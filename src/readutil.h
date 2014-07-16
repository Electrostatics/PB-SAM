#ifndef READUTIL_H
#define READUTIL_H

#include <vector>
#include <cstring>
#include "util.h"
#include "expansion.h"
 
#define GSIZE 401										//!< The length of one side of the DelPhi grid
#define GSIZE_SQ ( GSIZE*GSIZE )		//!< The square of GSIZE
#define GSIZE_CB ( GSIZE_SQ*GSIZE )	//!< GSIZE to the third power

const int MAX_CHARNUM = 300;				//!< maximum # of characters for filename ( including path )

class CMolecule;

//!  GridPnt structure 
/*! A  structure containing information
		on a grid point froma DelPhi grid  */
struct GridPnt
{
	GridPnt(  ) {}
	GridPnt( int i_, int j_, int k_ ) : i(i_), j(j_), k(k_) {}
	bool operator==( const GridPnt & g )
	{ return ( i == g.i && j == g.j && k == g.k ); }
	
	int i, j, k;
};	// end struct GridPnt


/*#########################################################*/
// various functions to read / write files

void readgrid(const char * fname, float * V, float * scale, 
							CPnt & pGridRef)
void readSurface( const char * fname, vector<CPnt> & SP, vector<CPnt> & NP  );

void readpqrAtoms( const char * fname, vector<CPnt> & pnt, vector<REAL> & ch, 
									vector<REAL> & R )
void readpqrCenters(const char * fname, vector<CPnt> & scen, vector<REAL> & srad, bool bNew)
void readpqr( const char * fname, vector<CPnt> & pnt, vector<REAL> & ch, 
						 vector<REAL> & rad, vector<CPnt> & cen );

void printM( double * A, int m, int n );
void writeMat_Unformatted( double *mat, int p, char* fname );
void readMat_Unformatted( double *mat, int maxp, char* fname );
void readMat_Formatted( double *mat, int maxp, char* fname );
void readMats( vector<REAL*> &iMats, int p, char* imatpath, int ncen );
void writeMat_Formatted( double *mat, int p, char* fname );

void writeExpansion( const CExpan &E, REAL rcut, char* fname );
void readExpansion( vector<double> &V, int maxp, char* fname, REAL & scale, REAL &rcut );
void readExpansionsP( vector<CMulExpan*> &iF, vector<CMulExpan*>  &iH , 
										 int p, char* exppath, char* runname, int ncen, REAL & intraRcutoff );

void generateConfigPQR( const char* configfname, const char* cenfname, 
											 const char* outfname );
void readConfig( const char* fname, vector<CPnt> &rcens, vector<CQuat> &rots );											 
											 
// Unused

/*

void readCenters( const char* fname, vector<CPnt> &cens, vector<double> &radii );

void readmdp( const char * fname, vector<CPnt> & pnt, vector<REAL> & ch, 
						 vector<REAL> & rad, vector<CPnt> & cen );
void writemdp( const char * ofname, const vector<CPnt> & scen, double * srad,
							const vector<CPnt> & pnt, const vector<double> ch );
							
void readOpenDX( const char * fname, vector<double> &V, int &sizex,int &sizey,int &sizez, 
								double &scalex,double &scaley,double &scalez, CPnt &pGridRef );

void readCompFile( char* fname, vector<CPnt> &pts, vector<REAL> &pot, int potcolumn );


void readConfigWithMolType( const char* fname, vector<CPnt> &rcens, 
													 vector<CQuat> &rots, vector<int> &moltypeid,													 
													 vector<int> &moltypes, vector<int> &nps );

void readExpansions( vector<CMulExpan> &iF, vector<CMulExpan>  &iH , 
										int p, char* exppath, char* runname, int ncen, REAL & intraRcutoff );
*/
											 
#endif
