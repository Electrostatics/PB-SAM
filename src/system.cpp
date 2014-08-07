#include "system.h"
#include "molecule.h"

// Static members of the CSystem class, must include an external definition
double CSystem::KAPPA;
double CSystem::m_sdiel;
double CSystem::m_temp;
double CSystem::BOXLENGTH;
bool CSystem::bPBC=false;

/******************************************************************/
/******************************************************************//**
* initConstants initializes constants of the system from user input:
******************************************************************/
void 
CSystem::initConstants( double kappa, double sdiel, double temp,
		       bool bPBC_, double boxlength )
{ 
	// Just copying inputs from BDmain and saving as system
	// constants
	KAPPA = kappa; 
	m_sdiel = sdiel; 
	m_temp = temp;	
	// Setting similar constants within CMolecule
	CMolecule::initConstants( kappa,sdiel );

	// Periodic boundary conditions, if it is there ( true )
	// then set the boxlength equal to boxlength,
	// else set BOXLENGTH equal to DBL_MAX
	// which stands for the max debye length
	bPBC = bPBC_;
	BOXLENGTH = boxlength;
	//BOXLENGTH = ( bPBC ? boxlength : DBL_MAX );
}	 // end InitConstants


/******************************************************************/
/******************************************************************//**
*  deletes any static array on heap
******************************************************************/
//void
//CSystem::deleteConstants(  )
//{
//	CMolecule::deleteConstants(  );
//	return;
//}	// end deleteConstants

/******************************************************************/
/******************************************************************//**
* Function to wrap a cartesian coordinate into its position in the
first box.
******************************************************************/
CPnt
CSystem::pbcPos( CPnt p )
{
	if( !bPBC ) return p; 
	return CPnt( p.x( )-floor( ( p.x() + BOXLENGTH/2.0 ) / BOXLENGTH ) * BOXLENGTH,
			 p.y( )-floor( ( p.y() + BOXLENGTH/2.0 ) / BOXLENGTH ) * BOXLENGTH,
			 p.z( )-floor( ( p.z() + BOXLENGTH/2.0 ) / BOXLENGTH ) * BOXLENGTH );
	//return CPnt( p.x( )-BOXLENGTH*round(p.x()/BOXLENGTH),
	//		p.y(  )-BOXLENGTH*round(p.y()/BOXLENGTH),
	//		p.z(  )-BOXLENGTH*round(p.z()/BOXLENGTH));
}

