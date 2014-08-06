#ifndef _SYSTEM_H_
#define _SYSTEM_H_

#include <cfloat>
#include "util.h"

//!  The CSystem class
/*!	The class that contains all the constants for a system  */
class CSystem
{
public:
	static double KAPPA;			//!< A floating point of the inverse debye length
	static double m_sdiel;		//!< A floating point for the solvent dielectric
	static double m_temp;			//!< A floating point of the system temperature
	static double BOXLENGTH;	//!< A floating point of the box length for PBCs
	static bool bPBC;					//!< A boolean indicating whether or not periodic boundary conds. are implemented

	//!  The CSystem initConstants function
	/*!	The function to initialize system constants  
	\param kappa a floating point of the inverse debye length
	\param sdiel a floating point for the solvent dielectric
	\param temp a floating point of the system temperature 
	\param bPBC a boolean indicating whether or not periodic boundary conds. are implemented
	\param boxlength a floating point of the box length for PBCs  */	
	static void initConstants( double kappa, double sdiel, double temp,
														bool bPBC=false, double boxlength=DBL_MAX );

	//!  The CSystem deleteConstants function
	/*!	The function to clear system constants  */	
	//static void deleteConstants(  );

	//!  The CSystem setBoxLength function
	/*!	The function to set the boxlength as that given from in [Angstroms]  
			\param bl a floating point of the boxlength to set */ 
	static void setBoxLength( double bl ) {BOXLENGTH = bl; }

	//!  The CSystem pbcPos function
	/*!	The function to move a point to its position in original
			image of a periodic system
			\param p an XYZ cartesian object to move
			\return an XYZ object of position in first periodic box*/		
	static CPnt pbcPos( CPnt p );
};

#endif
