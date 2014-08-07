#ifndef READUTIL_H
#define READUTIL_H

#include <vector>
#include <cstring>
#include "util.h"
#include "expansion.h"

#define GSIZE 401//301
#define GSIZE_SQ (GSIZE*GSIZE)
#define GSIZE_CB (GSIZE_SQ*GSIZE)

const int MAX_CHARNUM = 300; // maximum # of characters for filename (including path)

class CMolecule;

// A point on a DelPhi grid
struct GridPnt
{
  GridPnt() {}
  GridPnt(int i_, int j_, int k_) : i(i_), j(j_), k(k_) {}
  bool operator==(const GridPnt & g)
  { return (i == g.i && j == g.j && k == g.k); }

  int i, j, k;
};

// various functions to read / write files
void  readSurface(const char * fname, vector<CPnt> & SP, vector<CPnt> & NP );
void readmdp(const char * fname, vector<CPnt> & pnt, vector<REAL> & ch, 
	     vector<REAL> & rad, vector<CPnt> & cen);
void writemdp(const char * ofname, const vector<CPnt> & scen, double * srad,
	      const vector<CPnt> & pnt, const vector<double> ch);
void readpqr(const char * fname, vector<CPnt> & pnt, vector<REAL> & ch, 
	     vector<REAL> & rad, vector<CPnt> & cen);
void printM(double * A, int m, int n);
void writeMat_Unformatted(double *mat, int p, char* fname);
void readMat_Unformatted(double *mat, int maxp, char* fname);
void readMat_Formatted(double *mat, int maxp, char* fname);
void readMats(vector<REAL*> &iMats, int p, char* imatpath, int ncen);
void writeMat_Formatted(double *mat, int p, char* fname);
void readCenters(const char* fname, vector<CPnt> &cens, vector<double> &radii);
void readOpenDX(const char * fname, vector<double> &V, int &sizex,int &sizey,int &sizez, 
		double &scalex,double &scaley,double &scalez, CPnt &pGridRef);

void writeExpansion(const CExpan &E, REAL rcut, char* fname);
void readExpansion(vector<double> &V, int maxp, char* fname, REAL & scale, REAL &rcut);
void readExpansions(vector<CMulExpan> &iF, vector<CMulExpan>  &iH , 
		    /*vector<CLocalExpan>  &iLF, vector<CLocalExpan> &iLH,*/ 
		    int p, char* exppath, char* runname, int ncen, REAL & intraRcutoff);
void readExpansionsP(vector<CMulExpan*> &iF, vector<CMulExpan*>  &iH , 
		    int p, char* exppath, char* runname, int ncen, REAL & intraRcutoff);
void readCompFile(char* fname, vector<CPnt> &pts, vector<REAL> &pot, int potcolumn);
void readConfig(const char* fname, vector<CPnt> &rcens, vector<CQuat> &rots);

void readConfigWithMolType(const char* fname, vector<CPnt> &rcens, vector<CQuat> &rots, vector<int> &moltypeid,
                           vector<int> &moltypes, vector<int> &nps);

void generateConfigPQR(const char* configfname, const char* cenfname, const char* outfname);
#endif
