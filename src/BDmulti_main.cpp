#include "BDmulti.h"

/******************************************************************/
/******************************************************************/
/**
 *  This function runs the core of BDmulti, inputs:
 \param no of thread
 \param salt conc
 \param no. of trajectories
 \param distance for docking contact
 \param number of monomers ( type 0 )
 \param number of monomers ( type 1 )
 \param boxlength
 \param contactfilename
 \param maxtime
 \param runname
 \param system ( 0 = 2s, 1= 107, 2=354, 3=712 )
 \param restartfilename <<<limit to 12
 \param Mol0 pqr file
 \param Mol0 Imat path
 \param Mol0 expansion path
 \param Mol0 expname
 \param Mol1 pqr file
 \param Mol1 Imat path
 \param Mol1 expansion path
 \param Mol1 expname
 ******************************************************************/
int bd_multi( int argc, char ** argv )
{
  int num_threads = atoi(argv[1]);
  const double salt_conc = atof(argv[2]);
  const int ntraj = atoi(argv[3]);
  const double consepdist = atof(argv[4]);
  const int np1 = atoi(argv[5]);
  const int np2 = atoi(argv[6]);
  const double boxlength = atof(argv[7]);
  const char* contactfile = argv[8];
  const double maxtime = atof(argv[9]);
  char* runname = argv[10];
  const int systemIndex  = atoi(argv[11]);
	
  for(int i=0; i<argc; i++) cout <<argv[i]<<" "; cout<<"variables read"<<endl;
  cout.precision(5);
  long int idum=(unsigned)time(NULL);
  srand48(idum);  //this will initialize random number generator for BD simulation (S. Liu)
  cout <<"Random Seed = "<<idum<<endl;
#ifdef __OMP
  omp_set_num_threads(num_threads);
#endif
  // system parameters
  const double idiel = 4.0;
  const double sdiel = 78.0;
  const double temperature = 298.15;
  const double kappa = ANGSTROM * sqrt(  (2 * salt_conc * AVOGADRO_NUM / LITRE
																					* ELECT_CHG * ELECT_CHG)
																			 / (sdiel* PERMITTIVITY_VAC * KB * temperature ) );
	
  cout <<"EPS_SOLVENT:"<<sdiel<<" EPS_IN: "<<idiel<<" KAPPA: "<<kappa<<endl;
  cout <<" Np1 ="<<np1<<" Np2 = "<<np2<<" Boxlength ="<<boxlength<<endl;
  cout <<" -> Species 0 Paritial Concentration: "<<np1/pow(boxlength,3.0)<<" [A-3]"<<endl;
  cout <<" -> Species 1 Paritial Concentration: "<<np2/pow(boxlength,3.0)<<" [A-3]"<<endl;
	
  //initialize constants of the system, "true" means here we are using periodic boundary
  //condition (S. Liu)
  CSystem::initConstants(kappa, sdiel, temperature, true, boxlength);
  const int nmoltype=2;
  CBDmulti::initConstants(nmoltype);
  // read in contacts
  CContact::setSepDist(consepdist);
  vector<CMolContact> molcontactlist;
	
  int moltype;
  CMolContact::readMolContact(contactfile, moltype, molcontactlist);
  assert (moltype == 0); // since only two species, we always check species 0's contacts
  CBDmulti::addMolContact(moltype, molcontactlist);
  
  // create molecules
  vector<char*> molfnames1, molfnames2;
	molfnames1.resize(4);
	molfnames1[0] = argv[12];  //pqr file
	molfnames1[1] = argv[13]; //Imat path
	molfnames1[2] = argv[14]; //expansion path
	molfnames1[3] = argv[15]; //expname
	molfnames2.resize(4);
	molfnames2[0] = argv[16];  //pqr file
	molfnames2[1] = argv[17]; //Imat path
	molfnames2[2] = argv[18]; //expansion path
	molfnames2[3] = argv[19]; //expname
	
  // create BD system
  bool bRestart;
  char* configfname;
  //"bd" is a pointer that points to an object of class CBDmulti (S. Liu)
  CBDmulti* bd;
  
  if(argc==13)
	{
		//restart from an existing configuration file (S. Liu)
		bRestart = true;
		configfname = argv[12];
		cout <<"Initializing from restart file "<<configfname<<endl;
	}
  else
	{
		//restart from an initial configuration where all molecules are placed on a lattice (S. Liu)
		bRestart = false;
		cout <<"Initializing from lattice "<<endl;
	}
  //create "bd" with a constructor (S. Liu)
  bd = new CBDmulti(np1, np2, molfnames1,molfnames2, idiel);
  
  // start BD
  char outfname[MAX_CHARNUM];
  sprintf(outfname, "mfpt.%s.dat", runname);
  ofstream fout(outfname);
  double start = read_timer();
  int cdock = 0;
  for (int r=0; r<ntraj; r++)
	{
		double starti = read_timer();
		
		if(bRestart)
			bd->restartConfig(configfname);
		else
			bd->resetLattice(); //place molecules on a lattice (S. Liu)
    //"runFirstDockTime" will perform the BD simulation until one trajectory ends (S. Liu)
    double docktime = bd->runFirstDockTime( maxtime, runname );
    double endi = read_timer();
    if(docktime < maxtime ) cdock++;
    fout <<r<<" "<<docktime<<"[ns]; computation time ="<<endi-starti<< endl;
	}
  
  cout << "Total Docked "<<cdock<<" out of "<<ntraj<<" trajectories"<<endl;
  delete bd;
  double end = read_timer();
  
  cout <<"Program runtime = "<<end-start<<" s"<< endl;
  cout <<"Program completed successfully"<<endl;
  return 0;
} // end BD multi

/******************************************************************/
/******************************************************************/
/**
 * Main for BD multi
 ******************************************************************/
int main( int argc, char ** argv )
{
	return bd_multi( argc, argv );
}

