#include <vector>
#include "util.h"

void findCenters(const vector<CPnt> & P, const vector<double> & arad,
								 const vector<CPnt> & SP, const vector<CPnt> & NP,
								 vector<CPnt> & cens, vector<double> & R, double tolSP);

void middle(const vector<CPnt> & P, const vector<double> & arad, 
						const vector<CPnt> & SP, const vector<CPnt> & NP, 
						CPnt & cen, double & r, int & N, double tolSP,
						const vector<int> & indmax, vector<int> & ind2);

void subtract(const vector<int> & A, vector<int> & B);

bool IsInsideSurface(CPnt P, const vector<CPnt> & SP, const vector<CPnt> & NP);

void findNeighbors(int ki, const vector<CPnt> &cens, 
									 const vector<double> &radii, vector<int> &neigh);

bool IsSphereExposed(int ki, const vector<int> &neigh, const vector<CPnt> &cens, 
										 const vector<double> &radii, const vector<double> &rad2, 
										 const vector<CPnt> & SP,  const vector<CPnt> & NP);

/*#########################################################*/
/*#########################################################*/
// Functions for makesphere with contacts
/*#########################################################*/
/*#########################################################*/

void getInterfaceAtoms(const vector <CPnt> &pos1,const vector <CPnt> &pos2,
											 const vector <REAL> &rad1,const vector <REAL> &rad2,
											 vector <CPnt> & conpos1, double cutoff);
void classifyVertPoints(const vector<CPnt> & interfaceAtomPos,const vector<CPnt> & vertPts, 
												vector<int> & IPind, double atomVertCutoff);
void findCenters_interface(const vector<CPnt> & P, const vector<double> & arad,
													 const vector<CPnt> & SP, const vector<CPnt> & NP, 
													 const vector<int> & IPind, 
													 vector<CPnt> & cens, vector<double> & R, 
													 double tolSP, double tolIP);
void middle_interface(const vector<CPnt> & P, const vector<double> & arad, 
											const vector<CPnt> & SP, const vector<CPnt> & NP,
											const vector<int> & IPind, 
											CPnt & cen, double & r, int & N, double tolSP, double tolIP,
											const vector<int> & indmax, vector<int> & ind2);

/*
void middle_buried(const vector<CPnt> & P, const vector<CPnt> & SP, CPnt & cen, 
										double & r, int & N, const vector<int> & indmax, vector<int> & ind2);
void middle_surface(const vector<CPnt> & SP, const vector<CPnt> & NP,
											CPnt & cen, double & r, int & N, CPnt scen, double tolSP,  
											const vector<int> & indmax,  vector<int> & ind2);
void middle_surface(const vector<CPnt> & SP, const vector<CPnt> & NP, const vector<int> & IPind,
											CPnt & cen, double & r, int & N, CPnt scen, double tolSP, double tolIP, 
											const vector<int> & indmax,  vector<int> & ind2);
void findSolventCenters(const vector<CPnt> & SP,  const vector<int> & IPind, 
													vector<CPnt> & cens, vector<double> & R, const char* outfname); 
void middle_pt(const vector<CPnt> & P, const vector<CPnt> & SP, 
							 const vector<CPnt> & NP, 
							 CPnt & cen, double & r, int & N, double tolSP,
							 const vector<int> & indmax, vector<int> & ind2);
void getSpherePoints(int ki, const vector<int> &neigh, const vector<CPnt> &cens,
										 const vector<double> &radii, const vector<double> &rad2, 
										 const vector<CPnt> & SP,  const vector<CPnt> & NP, 
										 vector<CPnt> &SPE, vector<CPnt> &SPB);
bool IsInsideSurface(CPnt P, const vector<CPnt> & SP, 
										 const vector<CPnt> & NP, double & minsq);
bool IsSphereExposedPercent(int ki, const vector<int> &neigh, 
														const vector<CPnt> &cens, 
														const vector<double> &radii, const vector<double> &rad2, 
														const vector<CPnt> & SP,  const vector<CPnt> & NP);
 */

