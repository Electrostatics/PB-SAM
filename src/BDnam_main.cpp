#include "BDnam.h"


/******************************************************************/
/******************************************************************/
/**
 *  	BDNAM portion: initalizes many simulation parameters
 Calls on BDnam.cpp to gather input from external files
 and run BD
 \param no of thread
 \param salt conc
 \param no. of trajectories
 \param distance for docking contact
 \param output filename
 \param runname
 \param contact file location
 \param Mol1 pqr file
 \param Mol1 Imat pat
 \param Mol1 expansion path
 \param Mol1 expname
 \param Mol2 pqr file
 \param Mol2 Imat path
 \param Mol2 expansion path
 \param Mol2 expname
 ******************************************************************/
int bdNam(int argc, char ** argv)
{
	
  int num_threads = atoi(argv[1]);
  const double salt_conc = atof(argv[2]);
  int nTraj = atoi( argv[3]);
  double consepdist = atof(argv[4]);
  char* outfname = argv[5];
  char* runname = argv[6];
  char* contactfile = argv[7];
	
  for(int i=0; i<argc; i++) cout <<argv[i]<<" "; cout<<endl;

  cout.precision(5);
  long int idum=(unsigned)time(NULL);
  srand48(idum);
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
	
  CSystem::initConstants(kappa, sdiel, temperature); // non-periodic BC
  
  const int nmoltype=2;
  CBDnam::initConstants(nmoltype);
	
  // read in contacts
  CContact::setSepDist(consepdist);
  vector<CMolContact> molcontactlist;
  int moltype;
  CMolContact::readMolContact(contactfile, moltype, molcontactlist);
  assert (moltype == 0); // for nam;
  CBDnam::addMolContact(moltype, molcontactlist);
  vector<char*> molfnames1, molfnames2;
	
	molfnames1.resize(4);
	molfnames1[0] = argv[8];  //pqr file
	molfnames1[1] = argv[9]; //Imat path
	molfnames1[2] = argv[10]; //expansion path
	molfnames1[3] = argv[11]; //expname
	molfnames2.resize(4);
	molfnames2[0] = argv[12];  //pqr file
	molfnames2[1] = argv[13]; //Imat path
	molfnames2[2] = argv[14]; //expansion path
	molfnames2[3] = argv[15]; //expname

  CBDnam bd(molfnames1, molfnames2, idiel);
  ofstream fout(outfname);
  double start = read_timer();
  int cdock = 0;
  int scount;
  bool bFirstDock = false;
	
  for (int i = 0; i < nTraj; i++)
	{
		double starti = read_timer();
		CBDnam::STATUS status = bd.run( scount, runname );
		
		if (status == CBDnam::DOCKED)
		{
			cdock++;
			if(!bFirstDock)
	    {
	      char fname[MAX_CHARNUM];
	      sprintf(fname, "%s.mol0.dock.out", runname);
	      CMolecule::printMolConfig(bd.getMol(0), fname);
	      sprintf(fname, "%s.mol1.dock.out", runname);
	      CMolecule::printMolConfig(bd.getMol(1), fname);
	      bFirstDock = true;
	    }
			
		}
		
		double endi = read_timer();
		fout <<i<<" "<< CBDnam::STATUS_STRING[status] << " " << (int)status
	  << " ct " << scount <<" traj time ="<<endi-starti<<" s "<< endl;
	}
  fout <<" "<<cdock << " " << nTraj<<" delta = "<< (cdock / double(nTraj)) << endl;

  double end = read_timer();
  printf("Program run time = %.4lf seconds\n", end-start );
	CSystem::deleteConstants();
  return 0;
}

/******************************************************************/
/******************************************************************/
/**  *  	BDNAM main
 ******************************************************************/
int main( int argc, char ** argv ) {
	
	if ( argc!= 16 )
	{
		cout << "Incorrect usage: ./nam [# threads] [salt conc] ";
		cout << " [# traj] [Dock dist] [ofname] [runname] [contact file]";
		cout << " [Mol1 pqr] [Mol1 Imat] [Mol1 Exp] [Mol1 Expname] ";
		cout << " [Mol2 pqr] [Mol2 Imat] [Mol2 Exp] [Mol2 Expname] "<< endl;
		exit(0);
	}
	return bdNam( argc, argv );
	return 0;
}
