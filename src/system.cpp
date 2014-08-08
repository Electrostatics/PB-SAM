#include "system.h"
#include "molecule.h"

double CSystem::KAPPA;
double CSystem::m_sdiel;
double CSystem::m_temp;
double CSystem::BOXLENGTH;
bool CSystem::bPBC=false;

/******************************************************************/
/******************************************************************/
/**
 * initConstants initializes constants of the system from user input:
 ******************************************************************/
void
CSystem::initConstants(double kappa, double sdiel, double temp,
											 bool bPBC_, double boxlength)
{
  KAPPA = kappa;
  m_sdiel = sdiel;
  m_temp = temp;
  CMolecule::initConstants(kappa,sdiel);
  bPBC = bPBC_;
  BOXLENGTH = (bPBC ? boxlength : DBL_MAX);
}

/******************************************************************/
/******************************************************************/
/**
 *  deletes any static array on heap
 ******************************************************************/
void
CSystem::deleteConstants()
{
  CMolecule::deleteConstants();
  return;
}

/******************************************************************/
/******************************************************************/
/**
 * Function to wrap a cartesian coordinate into its position in the
 first box.
 ******************************************************************/
CPnt
CSystem::pbcPos(CPnt p)
{
  if(!bPBC) return p;
  return CPnt(p.x()-BOXLENGTH*round(p.x()/BOXLENGTH),
							p.y()-BOXLENGTH*round(p.y()/BOXLENGTH),
							p.z()-BOXLENGTH*round(p.z()/BOXLENGTH));
	/* 	return CPnt( p.x( )-floor( ( p.x() + BOXLENGTH/2.0 ) / BOXLENGTH ) * BOXLENGTH,
	 p.y( )-floor( ( p.y() + BOXLENGTH/2.0 ) / BOXLENGTH ) * BOXLENGTH,
	 p.z( )-floor( ( p.z() + BOXLENGTH/2.0 ) / BOXLENGTH ) * BOXLENGTH ); */
}
