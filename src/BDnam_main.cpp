#include "BDnam.h"


// param 1 : no of thread
// param 2 : salt conc
// param 3 : no. of trajectories
// param 4 : distance for docking contact
// param 5 : outfname
// param 6 : runname
// param 7 : contact file location
// param 8 : system (0 = 2s, 1=1brs_s5i1)

int bdNam(int argc, char ** argv)
{

  int num_threads = atoi(argv[1]);
  const double salt_conc = atof(argv[2]);
  int nTraj = atoi( argv[3]); 
  double consepdist = atof(argv[4]);
  char* outfname = argv[5];
  char* runname = argv[6];
  char* contactfile = argv[7];  
  int systemIndex = atoi(argv[8]);

  for(int i=0; i<argc; i++) cout <<argv[i]<<" "; cout<<endl;
//  char* hostname;
//  gethostname(hostname, 50);
//  cout << "Hostname: "<<hostname<<endl;

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
    molfnames1[0] = argv[9];  //pqr file
    molfnames1[1] = argv[10]; //Imat path
    molfnames1[2] = argv[11]; //expansion path
    molfnames1[3] = argv[12]; //expname

/*
    molfnames1[0] = "/home/shuleliu/mytrial/Example_1BRS/config/1BRS_chainA_p1.5d3_s5i1.pqr"; 
    molfnames1[1] = "/home/shuleliu/mytrial/Example_1BRS/DATA_s5i1/makeIMat/IMAT_P30_250k/chainA";
    molfnames1[2] = "/home/shuleliu/mytrial/Example_1BRS/DATA_s5i1/DATA_SELFPOL/IMAT_250k_based/SOLEXP_p30/salt_i0.05/chainA";
    molfnames1[3] = "1BRS_A_s5i1_i0.05_p30.0"; //expname
*/
    molfnames2.resize(4);
    molfnames2[0] = argv[13];  //pqr file
    molfnames2[1] = argv[14]; //Imat path
    molfnames2[2] = argv[15]; //expansion path
    molfnames2[3] = argv[16]; //expname

/*
    molfnames2[0] = "/home/shuleliu/mytrial/Example_1BRS/config/1BRS_chainD_p1.5d3_s5i1.pqr"; 
    molfnames2[1] = "/home/shuleliu/mytrial/Example_1BRS/DATA_s5i1/makeIMat/IMAT_P30_250k/chainD";
    molfnames2[2] = "/home/shuleliu/mytrial/Example_1BRS/DATA_s5i1/DATA_SELFPOL/IMAT_250k_based/SOLEXP_p30/salt_i0.05/chainD";
    molfnames2[3] = "1BRS_D_s5i1_i0.05_p30.0"; //expname
    break;
 */

      
  
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

int main(int argc, char ** argv)
{

  return bdNam(argc, argv);

  return 0;
}
