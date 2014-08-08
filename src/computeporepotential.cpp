#include "pore.h"

int computeporepotential(int argc, char** argv)
{
  const double salt_conc = atof(argv[1]);
  const double rpore = atof(argv[2]);
  const int nbin = atoi(argv[3]);
	
  const int np1 = 1;
  const int np2 = 1;
  const double boxlength = 114;
	
  // system parameters
  const double idiel = 2.1;
  const double sdiel = 80.0;
  const double temperature = 353;
  const double kappa = ANGSTROM * sqrt(  (2 * salt_conc * AVOGADRO_NUM / LITRE
																					* ELECT_CHG * ELECT_CHG)
																			 / (sdiel* PERMITTIVITY_VAC * KB * temperature ) );
	
  CSystem::initConstants(kappa, sdiel, temperature, false, boxlength);   //no PBC for potential calculation for the moment
  const int nmoltype=2;
  CComputePore::initConstants(nmoltype);
	
	
  // create molecules
  vector<char*> molfnames1, molfnames2;
	
	molfnames1.resize(4);
	molfnames1[0] = "/home/shuleliu/mytrial/Example_1BRS/config/cylinder6A12ring_onepore/cylinder6A12ring_onepore_p2d3_s5_polarize.pqr";
	molfnames1[1] = "/home/shuleliu/mytrial/Example_1BRS/Imat_nafion/cylinder6A12ring_onepore_p2d3_s5";
	molfnames1[2] = "/home/shuleliu/mytrial/Example_1BRS/SELFPOLARIZE_nafion/cylinder6A12ring_onepore_p2d3_s5";
	molfnames1[3] = "cylinder6A12ring_onepore_p2d3_s57.5_p30.0"; //expname
	
	molfnames2.resize(4);
	molfnames2[0] = "/home/shuleliu/mytrial/Example_1BRS/config/unitcharge/unitcharge_withcenter.pqr";
	molfnames2[1] = "/home/shuleliu/mytrial/Example_1BRS/Imat_nafion/unitcharge";
	molfnames2[2] = "/home/shuleliu/mytrial/Example_1BRS/SELFPOLARIZE_nafion/unitcharge";
	molfnames2[3] = "unitcharge7.5_p30.0"; //expname
	
	CComputePore* unitcharge;
	unitcharge = new CComputePore(np1,np2,molfnames1,molfnames2,idiel);
	
	const char* pfname = "unitpotential.dat";
	ofstream fout(pfname);
	fout.precision(8);
	
	double deltax = rpore/nbin;
	for(int i=0;i<nbin;i++)
	{
    unitcharge->setconfig(i, deltax);
    double unitpot=unitcharge->calculatePot(unitcharge->m_mols);
    fout<< deltax*i << unitpot<<endl;
	}
	
	fout.close();
	return 0;
}
int main(int argc, char ** argv)
{
	return computeporepotential(argc, argv);
}
