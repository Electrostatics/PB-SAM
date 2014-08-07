#include  <cstdlib>
#include <iostream>
#include <fstream>
#include <cfloat>
#include <stdio.h>
#include <string>
#include <ctime>
#include "molecule.h"
#include "system.h"
#include "readutil.h"

#define FACT2 5.604586e2 //(COUL_K*IKbT)

const double AVOGADRO_NUM = 6.0221415e23;
const double ELECT_CHG = 1.60217646e-19;
const double PERMITTIVITY_VAC = 8.854187817e-12;
const double KB = 1.3806503e-23;
const double ANGSTROM = 1e-10;
const double LITRE = 1e-3;
using namespace std;                       

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// param 1 = flag (spol)
// param 2 = number of threads
// param 3 = salt conc (M)
// param 4 = pqr filename
// param 5 = imat path
// param 6 = runname
// param 7 = interior dielectric 

int main_selfPolarize(int argc, char ** argv)

{
  cout << "running self-polarization ..." << endl;
  cout.precision(10);
  int num_threads = atoi(argv[2]);
#ifdef __OMP
  omp_set_num_threads(num_threads);
#endif

  double start = read_timer();

  //////////////// system parameters //////////////
  const double salt_conc = atof(argv[3]);
  const double idiel = ( argc == 7 ? 4.0 :  atof(argv[7]) );
  const double sdiel = 78;
  const double temp = 298.15;
  const double kappa = ANGSTROM * sqrt(  (2 * salt_conc * AVOGADRO_NUM / LITRE * ELECT_CHG * ELECT_CHG) 
		       / (sdiel* PERMITTIVITY_VAC * KB * temp ) );
  cout <<"EPS_SOLVENT:"<<sdiel<<" EPS_IN: "<<idiel<<" KAPPA: "<<kappa<<endl;

  CSystem::initConstants(kappa, sdiel, temp);
  
  //////////////// read in charge + centers //////////////
  vector<CPnt> scen, POS;
  vector<double> srad, CHG;
  const char * configfile  = argv[4];
  readpqr(configfile, POS, CHG ,srad, scen); // as PQR 

  //////////////// // define center as center of spheres ////////////////
  const int ncen = scen.size();// no of unique spheres
  CPnt initialrcen(0,0,0);
  for(int i=0; i<ncen; i++) initialrcen += scen[i]; initialrcen /= scen.size(); 
  cout <<"rcen : " <<initialrcen<<" ; no. of centers:  "<<ncen<<endl;

  //////////////// read in interaction matrices //////////////
  vector<REAL*> iMats;
  char* imatpath = argv[5];
  readMats(iMats, N_POLES, imatpath, ncen);
  
  
  /////////////// generate a copy of spx ///////////////
  vector<vector<CPnt> > SPxes;
  vector<int> nSPx; 
  vector<vector<int> > neighs;
  CMolecule::generateMolSPX(scen, srad, SPxes, nSPx, neighs );
  
  //////////////// create molecule  //////////////
  const double intraRcutoff = 30;  //set it big so that all spheres are considered near

  vector< vector<int> > intraPolLists_far, intraPolLists_near;
  CMolecule::generateMolTypeIntraPolLists(scen, srad, intraRcutoff,
                                          intraPolLists_near, intraPolLists_far);

  assert (intraPolLists_near.size() == ncen ); // all spheres are 'near' for selfpol

  CMolecule mol(initialrcen, scen, srad, CHG, POS, idiel, iMats, CMolecule::INTRA_RCUTOFF, 
		SPxes, nSPx, neighs, intraPolLists_near);
  double vm, rss;
  process_mem_usage(vm, rss);  cout << "After creating molecules: VM: " << vm << "; RSS: " << rss << endl;
    
  //////////////// self polarization  //////////////

  bool bPot = true;
  mol.polarize_self(bPot);

  double end = read_timer();

  cout <<"Total Solve Time [s] = "<<end-start<<endl;
  
  mol.writeMolExpansions(argv[6]);
    
  ////////////////  clean up memory    //////////////
  for(int i=0; i<ncen; i++) delete [] iMats[i];
  CSystem::deleteConstants();

  return 0;
}

// param 1 = flag (mat)
// param 2 = number of threads
// param 3 = salt conc (M)
// param 4 = pqr filename

int main_makeImat(int argc, char ** argv)
{

  cout <<"Building Surface Integral Matrices ..."<<endl;
  int num_threads = atoi(argv[2]);
#ifdef __OMP
  omp_set_num_threads(num_threads);
#endif

  double start = read_timer();

  // system parameters
  const double salt_conc = atof(argv[3]);
  const double idiel = 4.0;
  const double sdiel = 78;
  const double temp = 298.15;
  const double kappa = ANGSTROM * sqrt(  (2 * salt_conc * AVOGADRO_NUM / LITRE * ELECT_CHG * ELECT_CHG) 
		       / (sdiel* PERMITTIVITY_VAC * KB * temp ) );
  cout <<"EPS_SOLVENT:"<<sdiel<<" EPS_IN: "<<idiel<<" KAPPA: "<<kappa<<endl;

  CSystem::initConstants(kappa, sdiel, temp);
  
  // read in charge + centers 
  vector<CPnt> scen, POS;
  vector<double> srad, CHG;
  const char * configfile  = argv[4];
  readpqr(configfile, POS, CHG ,srad, scen); // as PQR 
  
  // define center as center of spheres
  CPnt rcen(0,0,0); 
  for(int i=0; i<scen.size(); i++) rcen += scen[i]; rcen /= scen.size(); 
  cout <<"rcen = "<<rcen<<endl;

  /////////////// generate a copy of moltype values ///////////////
  vector<vector<CPnt> > SPxes;
  vector<int> nSPx; 
  vector<vector<int> > neighs;
  CMolecule::generateMolSPX(scen, srad, SPxes, nSPx, neighs );

  const double intraRcutoff = 30;  //<<<<<<<<        
  vector< vector<int> > intraPolLists_far, intraPolLists_near;
  CMolecule::generateMolTypeIntraPolLists(scen, srad, intraRcutoff,
                                          intraPolLists_near, intraPolLists_far);

 // create molecule
  CMolecule mol(rcen, scen, srad, CHG, POS, idiel, SPxes, nSPx, neighs, intraPolLists_near);
  double end = read_timer();
  
  cout <<"Total Time to compute Imat [s] = "<<end-start<<endl;
  
  return 0;
}
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char ** argv)
{
  if (strncmp(argv[1], "spol", 4) == 0)
    return main_selfPolarize(argc, argv);  
  else if (strncmp(argv[1], "mat", 3) == 0)
    return main_makeImat(argc, argv);  
  
    return 0;
}
