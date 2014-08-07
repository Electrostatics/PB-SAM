#ifndef _EXPCENTER_H_
#define _EXPCENTER_H_
#include "expansion.h"
#include "rotcoeff.h"
#include "hash.h"

#define _PSQH_ N_POLES*(N_POLES+1)/2
#define _PSQH2_ _PSQH_ * _PSQH_
#define _PQUADH_ _PSQH_ * (_PSQH_+1) / 2
#define _PSQ_ N_POLES*N_POLES
#define _PQUAD_ _PSQ_*_PSQ_

////////////////////////////////////////////////
// CExpCenter
// Base Class for Expansion centers
////////////////////////////////////////////////
class CExpCenter
{

 public:
  static void initConstants(double sdiel, double kappa);
  
  CExpCenter();
  CExpCenter(int ki, CPnt cen, double rad);

  CExpCenter(int ki, CPnt cen, double rad, double idiel, bool bKappa,
	     const vector<double> &chgAssigned, const vector<CPnt> &posAssigned,
	     CQuat* orient=NULL, CRotCoeff* rot=NULL);

  void setOrder(const int p, const CRotCoeff &rot);
  void incOrder(const CRotCoeff & rot);
  virtual const CPnt computeTorqueOn_0(int i) const { /* later, for fixed charge */};
  const REAL computePot() const;

  const double getRad() const { return m_rad;}
  const CPnt getCen() const { return m_cen;}
  const CPnt getCenRot() const { return m_cenRot;}
  const CQuat getOrient() const {return *m_orient;}
  const int getOrder() const {return m_p;}
  const double getLScale() const {return m_lscale;}
  const double getMScale() const {return m_mscale;}
  const bool getbKappa() const {return m_bKappa;}

  const CMulExpan * const getpH() const  {return &m_H;}
  CMulExpan * const getpH() {return &m_H;}
  CMulExpan & getH() {return m_H;}
  const CMulExpan getH() const {return m_H;}
  CMulExpan & getHrot() {return m_Hrot;}
  const CMulExpan getHrot() const {return m_Hrot;}
  const double getHrotMono() const {return m_Hrot[0];}

  //  CLocalExpan & getL() {return m_L;}
  // const CLocalExpan & getL() const {return m_L;}
  CLocalExpan & getLS() {return m_LS;}
  const CLocalExpan & getLS() const {return m_LS;}
  const CMulExpan getrE() const {return m_rE;}

  CGradExpan & getgLHN() {return m_gLHN;}
  const CGradExpan getgLHN() const {return m_gLHN;}

  CExpCenter & operator=(const CExpCenter &M);

  // debug
  const CRotCoeff getRotCoeff() const {return *m_rot;} 

 protected:
  void setFixedCharges();
  //  void incRotate(const CRotCoeff &rot);


  static double m_sdiel, m_kappa;

  CPnt m_cen, m_cenRot;
  double m_rad, m_lscale, m_mscale;
  double m_idiel;
  bool m_bKappa;
  vector<double> m_ch;
  vector<CPnt> m_pos;

  //  CLocalExpan m_L; // only due to external source
  CLocalExpan m_LS; // only due to external source (selected)
  CMulExpan m_E, m_rE, m_H, m_Hrot; 
  //vector<CGradExpan> m_gH, m_gF, m_gHrot;
  CGradExpan m_gH /*, m_gF, m_gHrot*/;
  CGradExpan m_gLHN;

  CQuat* m_orient;
  CRotCoeff* m_rot;


  int m_p;
  int m_id;

};


////////////////////////////////////////////////
// CSolExpCenter
// Expansion centers with partially exposed surface
////////////////////////////////////////////////

class CSolExpCenter : public CExpCenter
{

 public:
  CSolExpCenter(int ki, CPnt cen, double rad, const vector<CPnt> &SPE, const vector<CPnt> &SPB,
		const vector<CPnt> &SPx, const int &nSPx, const vector<int> &neigh,
		const vector<int> &intraPolList_near,
		double idiel, bool bKappa, const vector<double> &allchg, const vector<CPnt> &allpos, 
		const vector<double> &chgAssigned, const vector<CPnt> &posAssigned);

  CSolExpCenter(int ki, CPnt cen, double rad, 
		const vector<CPnt> &SPx, const int &nSPx, const vector<int> &neigh,
		const vector<int> &intraPolList_near,
		double idiel, bool bKappa, const vector<double> &allchg, const vector<CPnt> &allpos, 
		const vector<double> &chgAssigned, const vector<CPnt> &posAssigned,REAL* iMat);

  CSolExpCenter(int ki, CPnt cen, double rad,
		const vector<CPnt> &SPx, const int &nSPx, const vector<int> &neigh,
		const vector<int> &intraPolList_near,
		double idiel, bool bKappa, const vector<double> &allchg, const vector<CPnt> &allpos,
		const vector<double> &chgAssigned, const vector<CPnt> &posAssigned, REAL* iMat,
		const CMulExpan *Fself, const CMulExpan *Hself,
		const CLocalExpan* LF_intraSelf, const CLocalExpan* LH_intraSelf,
		const CLocalExpan* LF_intraSelf_far, const CLocalExpan* LH_intraSelf_far,
		const double *qSolvedFself, const double  *qSolvedHself,
		const double *totalFself, const double  *totalHself,
		CQuat* orient, CRotCoeff* rot, bool bGrad);

  CSolExpCenter(int ki, CPnt cen, double rad,
		const vector<CPnt> &SPx, const int &nSPx, const vector<int> &neigh,
		const vector<int> &intraPolList_near,
		const CMulExpan *Fself, const CMulExpan *Hself,
		const double *qSolvedFself, const double  *qSolvedHself,
		const double *totalFself, const double  *totalHself);
    
    
  // for queries
  CSolExpCenter(int ki, CPnt cen, double rad, const CMulExpan *Hself, const vector<CPnt> &SPx=vector<CPnt>(0), 
		const int &nSPx=0, const vector<int> & intraPolList_near=vector<int>(0));

 ~CSolExpCenter()
   {   if(m_bOwnMat) delete [] m_IMat; }

  // public functions
  static void initConstants(double sdiel, double kappa);
  static void build_IMat(const vector<CPnt> &SPE, const vector<CPnt> &SPB, REAL* IMat);
  static void computeExposedSurfaceChargesFH_(const vector<CPnt> &SPx, 
					      const CMulExpan & F, const CMulExpan & H,
					      double *constH, double nSPx, 
					      vector<double> &qSolvedF, vector<double> &qSolvedH, 
					      double &totalF, double &totalH, int p);

  double solveSurfaceCharges();
  double solveSurfaceCharges(const CMulExpan &Finit, const CMulExpan &Hinit, bool bRotateH);
  double solveSurfaceGradient(int j, CGradExpan &gF, CGradExpan &gH, CGradExpan &gHrot, 
			      const CGradExpan & LGF_j, const CGradExpan & LGH_j);
  void initGradient(int nmol);
  //void resetGradient(int p);
  //void setGradientRange(int p);
  void resetExpansions(const CMulExpan &iF, const CMulExpan &iH);
  void setToSelfPolValues(bool bResetCharges);
  const vector<REAL> & getExposedSurfaceChargesH() const {return m_qSolvedH; } 
  const vector<REAL> & getExposedSurfaceChargesF() const {return m_qSolvedF; }
  //void calTQ_forPole();

  void rotateCenters();
  void rotateHself();

#ifdef __NOPOL_UNCHANGED__
  void rotateCurrentH();
  CGradExpan rotateGH() const; 
#endif

  REAL getDev() const {return m_dev;}
  CMulExpan & getF() {return m_F;}
  const CMulExpan getF() const {return m_F;}

  CGradExpan & getGH() {return m_gH;}
  const CGradExpan getGH() const {return m_gH;}

  const CPnt getGHMono() const {return CPnt(m_gH[0][0],m_gH[1][0],m_gH[2][0]);}

  CLocalExpan & getLHS() {return m_LHS;}
  const CLocalExpan getLHS() const {return m_LHS;}
  CLocalExpan & getLFS() {return m_LFS;}
  const CLocalExpan getLFS() const {return m_LFS;}

  CLocalExpan & getLFS_Far() {return m_LFS_Far;}
  const CLocalExpan getLFS_Far() const {return m_LFS_Far;}
  CLocalExpan & getLHS_Far() {return m_LHS_Far;}
  const CLocalExpan getLHS_Far() const {return m_LHS_Far;}


  const vector<CPnt> & getSPx() const {return m_SPx;}

  const vector<double> & getQSolvedF() const {return m_qSolvedF;}
  const vector<double> & getQSolvedH() const {return m_qSolvedH;}
  const vector<CPnt> & getGSolvedF() const {return m_gSolvedF;}
  const vector<CPnt> & getGSolvedH() const {return m_gSolvedH;}
  const int getSPExSize() const {return m_SPExSize;}

  const REAL& getTQH() const {return m_totalH;}
  const REAL& getTQF() const {return m_totalF;}
  const REAL& getTGH() const {return m_maxgH;}
  const REAL& getTGF() const {return m_maxgF;}

  const CMulExpan &getFself() const {return *m_Fself;}
  const CMulExpan &getHself() const {return *m_Hself;}
  const CMulExpan &getHselfRot() const {return m_HselfRot;}
  const CLocalExpan getLF_intraSelf() const {return *m_LF_intraSelf;}
  const CLocalExpan getLH_intraSelf() const {return *m_LH_intraSelf;}
  const CLocalExpan getLF_intraSelf_far() const {return *m_LF_intraSelf_far;}
  const CLocalExpan getLH_intraSelf_far() const {return *m_LH_intraSelf_far;}
 

  const REAL* getQSolvedFself() const {return m_qSolvedFself;}
  const REAL* getQSolvedHself() const {return m_qSolvedHself;}
  const REAL getTQFself() const {return *m_totalFself;}
  const REAL getTQHself() const {return *m_totalHself;}

  CSolExpCenter & operator=(const CSolExpCenter &M);

  const vector<int> & getIntraPolList_near() {return m_intraPolList_near;} 

  void clearInterPolList() {
    m_interactList.clear(); m_interPolList.clear(); /*m_interOverlapPolList.clear();*/ }
  void addInterPolList(CFullSphereID c) {m_interPolList.push_back(c);}
  void addInteractList(CFullSphereID c) {m_interactList.push_back(c);}
  //  void addInterOverlapPolList(CFullSphereID c) {m_interOverlapPolList.push_back(c);}
  const vector<CFullSphereID> & getInterPolList() const {return m_interPolList;} 
  const vector<CFullSphereID> & getInteractList() const {return m_interactList;} 
  //  const vector<CFullSphereID> & getInterOverlapPolList() const {return m_interOverlapPolList;} 

#ifdef __NOPOL_UNCHANGED__
  void setbInterPolListChanged(bool bChanged) {m_bInterPolListChanged = bChanged;}
  bool IsInterPolListChanged() const {return m_bInterPolListChanged;}
#endif
  bool IsEmptyInterPolList() const { return m_interPolList.size() == 0 
				       /*&& m_interOverlapPolList.size() == 0*/;}
  bool IsNoInteractionList() const { return m_interPolList.size() == 0 &&
				       m_interactList.size() == 0;  } 
  
  bool isOnInterPolList(int j, int kj) const;

  virtual const CPnt computeTorqueOn_0(int i) const;
  const CPnt computeForceOn_0() const;

  static double CONST2[N_POLES];

 private:
// private functions 
  static int & IDK(int l, int s)    {return idk[id[l]+s];}
  static void computeIntegralE(double * Yrr,double * Yri,double * Yir,double * Yii, 
			       const vector<CPnt> & SPE, const vector<CPnt> & SPB);
  static void computeIntegralB(double * Yrr,double * Yri,double * Yir,double * Yii, 
			       const vector<CPnt> & SPE, const vector<CPnt> & SPB );

  virtual void setFixedCharges(const vector<double> &allchg, const vector<CPnt> &allpos);
  void initMyConstants();
  void initSurfaceCharges();
  void computeExposedSurfaceChargesH();
  void computeExposedSurfaceChargesF();
  void computeExposedSurfaceChargesFH();
  void computeExposedSurfaceGradientFH(const CGradExpan &gF, const CGradExpan &gH);
  void build_IMatE();

  void cal_DVFixed();
  void getSurfaceChargesH(const vector<CPnt> &SPx, vector<double> &qS);
  void getSurfaceChargesF(const vector<CPnt> &SPx, vector<double> &qS);
  void revertF(double * F, int p, int D) const ;
  void revertH(double * H, int p, int D) const ;

  double solveFH_N(const vector<double> &LH, const vector<double>  &LF, int pm, int pn);
  double solveFH(vector<double> &F, vector<double>  &H,
		 const vector<double> &LF, const vector<double>  &LH, int pm, int pn, int D) const;
  //  double solveFH(const vector<double> &LH, const vector<double>  &LF, int pm, int pn);
  void compute_FHbase(const double* LH, const double *LF, double* fbase, double *hbase, int pm) const;
  void compute_FHbase3(const double* LH, const double *LF, double* fbase, double *hbase, int pm) const;
  void compute_fx(double *x, const double *fbase, const double *F, const double *H, int p, int D) const;
  void compute_hx(double *x, const double *hbase, const double *F, const double *H, int p, int D) const;

  static double computeDev(const double *M1, const double *M2, int p, int D);
  static double computeDev2(const double *M1, const double *M2, int p, int D);

  // utilities
  static void printY(double *Y) ;
  static void printYFull(double *Y) ;
  static void printMat(double *mat) ;
  static void printMat(const vector<double> &mat) ;
  void testY(double *Y1, double *Y2) const;
  void testYFull(double *Y1, double *Y2) const;
  void checkSolveDifferenceDirichlet();
  void checkSolveDifferenceVonNeumann();
  //void checkSolveDifferenceBuriedCharges();
  //void printCharges();

  // private variables
  static double CONST3[N_POLES];
  static int id[N_POLES], idk[_PSQH_];
  static double & getYY(double * YY, int n, int m, int l, int s) {return YY[ IDK(l,s) + id[n] + m];}

  vector<CPnt> m_SPE, m_SPB;
  vector<int> m_neigh;//, m_intraPolList_Near,m_intraPolList_Far;
  const vector<int> & m_intraPolList_near;
  vector<CFullSphereID> m_interPolList, m_interactList; //m_interOverlapPolList;
  vector<double> m_qSolvedF, m_qSolvedH;
  vector<CPnt> m_gSolvedF, m_gSolvedH;
  REAL* m_IMat;
  double m_Dfix[_PSQ_], m_Vfix[_PSQ_];
  REAL m_totalH, m_totalF, m_maxgH, m_maxgF, m_dev;

  REAL m_dA;

  double m_constH1[N_POLES], m_constF1[N_POLES], m_constH2[N_POLES], m_constF2[N_POLES], 
    m_constLH1[N_POLES], m_constLF1[N_POLES],m_constLH2[N_POLES], m_constLF2[N_POLES];
  double m_constInvH[N_POLES], m_constChgH[N_POLES];

  CLocalExpan m_Lfix, /*m_LH, m_LF,*/ m_LHS, m_LFS, m_LFS_Far, m_LHS_Far;
  CMulExpan m_Efix, m_F, m_Ev, m_HselfRot; 
  bool m_bRead, m_bOwnMat, m_bGrad;

  const CMulExpan *m_Fself, *m_Hself;
  const CLocalExpan *m_LF_intraSelf, *m_LH_intraSelf, *m_LF_intraSelf_far, *m_LH_intraSelf_far;
  const double *m_qSolvedFself, *m_qSolvedHself;
  const double *m_totalFself, *m_totalHself;
  const vector<CPnt> &m_SPx;
  const int m_SPExSize; 
  const int &m_nSPx;
#ifdef __NOPOL_UNCHANGED__
  bool m_bInterPolListChanged;
#endif

};

inline void
CExpCenter::setOrder(const int p, const CRotCoeff &rot)
{
  assert(p <= N_POLES);
  cout <<"setorder is used" <<endl;
  if (m_p < p)  
    while ( m_p < p)    
      incOrder(rot);
  
  else if (m_p > p)
    {
      // (LATER) m_rT[ki].setRange(p);
      m_H.setRange(p); //CHECK LATER
      m_p = p;
    }
  
  /* (LATER)
     if(ki == 0)
     for(int k=0; k<m_ncenter; k++)
     m_pG[k].setRange(p);
  */
}
/*
inline void
CExpCenter::incOrder(const CRotCoeff &rot)
{
  m_p++;

  assert(m_p <= rot.getOrder());
  incRotate(rot);
  CMulExpan & getF() {return m_F;}
  const CMulExpan getF() const {return m_F;}

  const CMulExpan &getFself() const {return *m_Fself;}
  const CMulExpan &getHself() const {return *m_Hself;}

  int getGFsize() const {return m_gF.size();}

  CGradExpan & getGF(int j) {return m_gF[j];}
  const CGradExpan getGF(int j) const {return m_gF[j];}
  CGradExpan & getGH(int j) {return m_gH[j];}
  const CGradExpan getGH(int j) const {return m_gH[j];}
  CGradExpan & getGHrot(int j) {return m_gHrot[j];}
  const CGradExpan getGHrot(int j) const {return m_gHrot[j];}


}
*/
#endif

