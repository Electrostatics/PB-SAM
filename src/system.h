#ifndef _SYSTEM_H_
#define _SYSTEM_H_

#include "util.h"
#include <cfloat>

class CSystem
{
 public:
  static void initConstants(double kappa, double sdiel, double temp,
			    bool bPBC=false, double boxlength=DBL_MAX);
  static void deleteConstants();
  static void setBoxLength(double bl) {BOXLENGTH = bl; }
  static CPnt pbcPos(CPnt p);

  static double KAPPA, m_sdiel, m_temp, BOXLENGTH;
  static bool bPBC;
};

#endif
