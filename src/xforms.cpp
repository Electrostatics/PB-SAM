#include <cmath>
#include "xforms.h"
#include "lautil.h"
#include "molecule.h"

CXFormA::CXFormA(const CMulExpan & M1, const CMulExpan & M2, bool bGrad)  
  : m_pM1(&M1), m_pM2(&M2) , m_pH(1) ,  m_transH(bGrad), m_rot(bGrad)
{

  REAL scale1 = 1.0/(m_pM1->getScale()); //center kj (in)
  REAL scale2 = 1.0/(m_pM2->getScale()) ; //center ki (out)
  m_transH.initScale( scale1 ,  scale2 );  

}

CXFormA::CXFormA(const CSolExpCenter & C1, const CSolExpCenter & C2, bool bGrad)  
  : m_pM1(C1.getpH()), m_pM2(C2.getpH()), m_pH(1) , m_transH(bGrad), m_rot(bGrad)
{

  REAL scale1 = C1.getLScale(); //center kj (in)
  REAL scale2 = C2.getLScale(); //center ki (out)
  m_transH.initScale( scale1 ,  scale2 );  
}

void
CXFormA::initConstants()
{
  CRotCoeff::initConstants();
  CTransCoeff::initConstants();
}

// Initialize to a MtoL transform by the vector P.
void 
CXFormA::reset(const CPnt & P, int p)
{
  assert(p > 0 && p <= N_POLES);
  
  m_pH = p;
  CSpPnt S = CartToSph(P);
  
  m_rot.reset(S.theta(), S.phi(), M_PI, p);// itay's (correct)
  
  m_transH.reset(S.rho(), CRExpan::KAPPA, p);

  /*
  // not used - pole is precalculate using computeorder in generatexform  
  // calculate relative error of using p-1 vs. p poles  
  m_rot.rotate(*m_pM1, m_tM1, 1, p, true);
  m_rot.rotate(*m_pM2, m_tM2, 1, p, true);

  m_trans.translate(m_tM1, m_tL, m_p, false);
  m_base = inprod(m_tM2, m_tL);	  
  m_resid1 = (p > 1 ? m_trans.computeError(m_tM1, m_tM2, p-1) : 0.0);
  m_resid2 = m_trans.computeError(m_tM1, m_tM2, p);
  compRelError();
  */
  
  // prepare matrix to convert gradient from spherical to cartesian coordinates   
  REAL sint, cost, sinp ,cosp;
  m_rot.getParams(sint, cost, sinp, cosp);
  REAL ir = 1.0/S.rho();
  
  if (m_rot.isSingular())
    {
      m_R[0] = CPnt(0.0, ir*cost, 0.0);
      m_R[1] = CPnt(0.0, 0.0, ir);
      m_R[2] = CPnt(cost, 0.0, 0.0);
    }
  else 
    {
  
      REAL irst = ir/sint, irct = ir*cost;
            
      m_R[0] = CPnt(sint*cosp, irct*cosp, -irst*sinp);
      m_R[1] = CPnt(sint*sinp, irct*sinp, irst*cosp);
      m_R[2] = CPnt(cost, -ir*sint, 0.0);
    }
  
}


// Compute the derivatives of the transformed MP coeffs
// ( d_T(i,j).A(j) )  
void 
CXFormA::xformH(const CExpan & Min, CGradExpan & Gout, bool bFor)
{
  assert(Min.getRange().p2() >= m_pH);
  CMulExpan tM1;
  CLocalExpan tL;

  // Derivative with respect to rho
  m_rot.rotate(Min, tM1, true);
  m_transH.dTranslate(tM1, tL, !bFor); 
  m_rot.rotate(tL, Gout[dRHO], false);
   
  // Derivative with respect to theta and phi (first part)
  m_transH.translate(tM1, tL, !bFor);
  m_rot.dRotateT(tL, Gout[dTHETA], false);
  m_rot.dRotateP(tL, Gout[dPHI], false);
    
  // Derivative with respect to theta (second part)
  m_rot.dRotateT(Min, tM1, true);
  m_transH.translate(tM1, tL, !bFor);
  m_rot.rotate(tL, tM1, false);
  Gout[dTHETA] += tM1;
  
  // Derivative with respect to phi (second part)
  m_rot.dRotateP(Min, tM1, true);
  m_transH.translate(tM1, tL, !bFor);
  m_rot.rotate(tL, tM1, false);
  Gout[dPHI] += tM1;

  sphToCart(Gout);
}



// Transform the MP coeffs
void 
CXFormA::xformH(const CMulExpan & Min, CLocalExpan & Mout, bool bFor)
{

  //  cout <<"called CXFormA:: xformH "<<endl;
  CMulExpan tM1;
  CLocalExpan tL;
  
  //  cout<<Min.getRange().p2()<<" "<<m_p<<endl;
  assert(Min.getRange().p2() >= m_pH);
  
  //cout<<"before rotate "<<Min.getScale()<<Min<<endl;

  m_rot.rotate(Min, tM1, 1, m_pH, true);
  //  if(tM1.isBlownup()) cout <<"IsBlown in xformH tM1\n"<<tM1<<endl;
  //  cout<<"before translate "<<m_tM1.getScale()<<m_tM1<<endl;
 
  m_transH.translate(tM1, tL, !bFor); 
  //  if(tL.isBlownup()) cout <<"IsBlown in xformH tL\n"<<tL<<endl;
  //cout<<"after translate "<<m_tL.getScale()<<m_tL<<endl;

  m_rot.rotate(tL, Mout,  1, m_pH, false);
  //  if(Mout.isBlownup()) cout <<"IsBlown in xformH Mout\n"<<Mout<<endl;
  //cout<<"after rotate "<<Mout.getScale()<<Mout<<endl;

}

// Transform the MP coeff triplet (for T.d_A)
// Why use TriExpan and not CGradExpan? 
void 
CXFormA::xformH(const CTriExpan & Gin, CTriExpan & Gout, bool bFor)
{
  CLocalExpan tL[3];
  CMulExpan tM[3];
  tM[0] = CMulExpan(Gin[0].getVector(), Gin.getRange(), Gin.getScale()); 
  tM[1] = CMulExpan(Gin[1].getVector(), Gin.getRange(), Gin.getScale()); 
  tM[2] = CMulExpan(Gin[2].getVector(), Gin.getRange(), Gin.getScale()); 

  xformH(tM[0], tL[0], bFor);
  xformH(tM[1], tL[1], bFor);
  xformH(tM[2], tL[2], bFor);

  Gout[0] = CExpan(tL[0].getVector(), tL[0].getRange(), tL[0].getScale()); 
  Gout[1] = CExpan(tL[1].getVector(), tL[0].getRange(), tL[0].getScale()); 
  Gout[2] = CExpan(tL[2].getVector(), tL[0].getRange(), tL[0].getScale()); 

}



/*
void
CXFormA::sphToCart(CPnt & p)
{
  p = CPnt(dot(m_R[0],p), dot(m_R[1],p), dot(m_R[2],p));
}
*/

void
CXFormA::sphToCart(CGradExpan & G)
{
  for (int n = 0; n < G.getRange().p2(); n++)
    {
      // m = 0 case
      CPnt ar(G[dRHO](n,0), G[dTHETA](n,0), G[dPHI](n,0));
      CPnt br(dot(m_R[0],ar), dot(m_R[1],ar), dot(m_R[2],ar));
      
      G[dX](n,0) = br.x();
      G[dY](n,0) = br.y();
      G[dZ](n,0) = br.z();
      
      for (int m = 1; m <= n; m++)
	{
	  CPnt ar(G[dRHO](n,2*m-1), G[dTHETA](n,2*m-1), 
		  G[dPHI](n,2*m-1));
	  CPnt ai(G[dRHO](n,2*m), G[dTHETA](n,2*m), 
		  G[dPHI](n,2*m));
	  
	  CPnt br(dot(m_R[0],ar), dot(m_R[1],ar), dot(m_R[2],ar));
	  CPnt bi(dot(m_R[0],ai), dot(m_R[1],ai), dot(m_R[2],ai));
	  
	  G[dX](n,2*m-1) = br.x(); G[dX](n,2*m) = bi.x();
	  G[dY](n,2*m-1) = br.y(); G[dY](n,2*m) = bi.y();
	  G[dZ](n,2*m-1) = br.z(); G[dZ](n,2*m) = bi.z();

	}
    }
}

/////////////////////////////////////////////////////////////////

CXFormAIntra::CXFormAIntra(const CSolExpCenter & C1, const CSolExpCenter & C2)  
  : CXFormA(C1, C2), m_transF(false)
{

  REAL scale1 = C1.getLScale(); //center kj (in)
  REAL scale2 = C2.getLScale();; //center ki (out)
  m_transF.initScale( scale1 ,  scale2 );  
}

void
CXFormAIntra::reset(const CPnt & P, int pH, int pF)
{

  assert(pH > 0 && pH <= N_POLES);
  assert(pF > 0 && pF <= N_POLES);
  

  int maxp = ( pH > pF ? pH : pF);
  m_pH = pH;
  m_pF = pF;
  CSpPnt S = CartToSph(P);
  
  m_rot.reset(S.theta(), S.phi(), M_PI, maxp);// itay's (correct)
  
  m_transH.reset(S.rho(), CRExpan::KAPPA, pH);
  m_transF.reset(S.rho(), 0.0, pF);

}

// Transform the MP coeff triplet (for T.d_A)
// Why use TriExpan and not CGradExpan? 
void 
CXFormAIntra::xformF(const CTriExpan & Gin, CTriExpan & Gout, bool bFor)
{
  CLocalExpan tL[3];
  CMulExpan tM[3];
  tM[0] = CMulExpan(Gin[0].getVector(), Gin.getRange(), Gin.getScale()); 
  tM[1] = CMulExpan(Gin[1].getVector(), Gin.getRange(), Gin.getScale()); 
  tM[2] = CMulExpan(Gin[2].getVector(), Gin.getRange(), Gin.getScale()); 

  xformF(tM[0], tL[0], bFor);
  xformF(tM[1], tL[1], bFor);
  xformF(tM[2], tL[2], bFor);

  Gout[0] = CExpan(tL[0].getVector(), tL[0].getRange(), tL[0].getScale()); 
  Gout[1] = CExpan(tL[1].getVector(), tL[0].getRange(), tL[0].getScale()); 
  Gout[2] = CExpan(tL[2].getVector(), tL[0].getRange(), tL[0].getScale()); 

}
/*
void
CXFormAIntra::incOrderH()
{
  assert(m_pH < N_POLES);

  //  m_p++;//???
  m_pH++;
  m_transH.incOrder();
  if (m_rot.getOrder() < m_pH) m_rot.incOrder();

}
*/

/* not working yet
void
CXFormAIntra::decOrderH()
{
  assert (m_p > 1);
    
  m_p--;

  m_transH.decOrder();
  m_rot.decOrder();

}
*/

/*
void CXFormAIntra::setOrderF(int p)
{
  assert(p >= 1 && p <= N_POLES);
  if (m_pF < p)  
    while (m_pF < p)
      incOrderF();
  else if (m_pF > p)
    {
      cout <<"does not do dec yet "<<endl; exit(1);
    }
 
}
*/
/*
void CXFormAIntra::setOrderH(int p)
{
  assert(p >= 1 && p <= N_POLES);
  if (m_pH < p)  
    while (m_pH < p)
      incOrderH();
  else if (m_pH > p)
    {
      cout <<"does not do dec yet "<<endl; exit(1);
    }
 
}
*/
/*
void
CXFormAIntra::incOrderF()
{
  assert(m_pF < N_POLES);

  //  m_p++;//???
  m_pF++;
  m_transF.incOrder();
  if (m_rot.getOrder() < m_pF) m_rot.incOrder();

}
*/

// Transform the MP coeffs with kappa = 0
void 
CXFormAIntra::xformF(const CMulExpan & Min, CLocalExpan & Mout, bool bFor)
{
  //  cout <<"in xformf"<<endl;
  CMulExpan tM1;
  CLocalExpan tL;
  //  if(Min.getRange().p2() < m_pF) cout <<"AxformF "<<Min.getRange().p2()<<" "<< m_pF<<endl;
  assert(Min.getRange().p2() >= m_pF);
  m_rot.rotate(Min, tM1, 1, m_pF, true);
  m_transF.translate(tM1, tL, !bFor); 
  m_rot.rotate(tL, Mout, 1, m_pF, false);

}

/////////////////////////////////////////////////////////////////
// CXFormN : Numerical Transforms
/////////////////////////////////////////////////////////////////


// Constructors
CXFormN::CXFormN(const CSolExpCenter &C1, const CSolExpCenter &C2)
  : m_reset(false), m_qSolvedH1(C1.getQSolvedH()),  m_qSolvedH2(C2.getQSolvedH()), 
    m_gSolvedH1(C1.getGSolvedH()), m_gSolvedH2(C2.getGSolvedH())

{

  m_pM1       = C1.getpH();
  m_SPx1      = &(C1.getSPx()[0]);
  m_SPExSize1 = C1.getSPExSize();

  m_scale1    = C1.getLScale();
  m_tQH1       = &(C1.getTQH()); 
  m_tGH1       = &(C1.getTGH()); 

  m_pM2       = C2.getpH();
  m_SPx2      = &(C2.getSPx()[0]);
  m_SPExSize2 = C2.getSPExSize();


  m_scale2   = C2.getLScale();
  m_tQH2       = &(C2.getTQH()); 
  m_tGH2       = &(C2.getTGH()); 

  m_d1.resize(m_SPExSize2 );
  m_d2.resize(m_SPExSize1 );
  
 }
/*
// Constructors
CXFormN::CXFormN(const CSolExpCenter &C1, const CMolecule &M2, REAL scale)
  : m_reset(false)
{

  m_SPx1      = &(C1.getSPx()[0]);
  m_SPExSize1 = C1.getSPExSize();
  m_qSolvedH1  = &( C1.getQSolvedH()[0] );
  m_scale1    = C1.getLScale();
  m_tQH1       = &(C1.getTQH()); 

  m_SPx2      = NULL;
  m_SPExSize2 = 0;
  m_qSolvedH2  = NULL;
  m_scale2  = scale;
  m_tQH2       = NULL; 

  m_d2.resize(m_SPExSize1 );

 }
*/
// Functions
void
CXFormN::reset(const CPnt & P, int p)
{
  //  cout <<"reseting cxformN "<<P<<endl;
  // recalculate distance between surface charge points and target centers
  m_P = P;
  m_resetI = false;
  m_resetJ = false;

  return;
}

// get qsolvedExposed directly from CSolExpcenter
// Note : Min is not used - we get qsolvedH directly
void 
CXFormN::xformH(const CMulExpan & Min, CLocalExpan & Lout, bool bFor)
{
  int p;
  REAL rho = m_P.norm();

  // Transform J(1) to I(2)
  if(bFor)
    {  
      if(!m_resetJ )
	{
	  for(int i=0; i<m_SPExSize1; i++) m_d2[i] = m_SPx1[i] - m_P ; 
	  m_resetJ = true;
	}

      p = computeOrder(rho, *m_tQH1, m_scale1);

      Lout = CLocalExpan(m_qSolvedH1, m_d2, p, true, m_scale2);
    }

  // Transform I(2) to J(1)
  else
    {
      if(!m_resetI )
	{
	  for(int i=0; i<m_SPExSize2; i++) m_d1[i] = m_SPx2[i] + m_P ; 
	  m_resetI = true;
	}

      p = computeOrder(rho, *m_tQH2, m_scale2);


      Lout = CLocalExpan(m_qSolvedH2, m_d1, p, true, m_scale1);
    }

  return;
}

// Note : Gin is not used - we get gsolvedH directly
void 
CXFormN::xformH(const CTriExpan & Gin, CTriExpan & Gout, bool bFor)
{
  int p;
  REAL rho = m_P.norm();
  
  vector<CPnt> gSolved;
  // Transform J(1) to I(2)
  if(bFor)
    {  
      if(!m_resetJ )
	{
	  for(int i=0; i<m_SPExSize1; i++) m_d2[i] = m_SPx1[i] - m_P ; 
	  m_resetJ = true;
	}

      p = computeOrder(rho, *m_tGH1, m_scale1);
      //cout <<"fwd p = "<<p<<endl;
      //gSolved.assign(m_gSolvedH1, m_gSolvedH1+m_SPExSize1);
      Gout = CGradExpan(&(m_gSolvedH1[0]), m_d2, p, true, m_scale2, false);
    }

  // Transform I(2) to J(1)
  else
    {
      if(!m_resetI )
	{
	  for(int i=0; i<m_SPExSize2; i++) m_d1[i] = m_SPx2[i] + m_P ; 
	  m_resetI = true;
	}

      p = computeOrder(rho, *m_tGH2, m_scale2);
      //cout <<"tp p = "<<p<<endl;
      //gSolved.assign(m_gSolvedH2, m_gSolvedH2+m_SPExSize2);
      Gout = CGradExpan(&(m_gSolvedH2[0]), m_d1, p, true, m_scale1, false);
    }
  
  return;
}


int 
CXFormBase::computeOrder(REAL rho, REAL sourceTQ, REAL sourceScale)
{
  //return N_POLES;
  
  if(sourceTQ <1e-10) return 1;

  REAL c = rho/sourceScale- 1;
  /*
  if(c < 0 ) 
    {
      cout <<"Error! rho ("<<rho<<") cannot be smaller than bounding radius ("<<sourceScale<<")!"<<endl;
      exit(1);
    }
  */
  if( c <= 1 ) return N_POLES;

  int p = int (log(sourceTQ / (MAX_ERROR * sourceScale * (c-1)) ) / log(c) - 1 + 0.5); // 0.5 for round up
  //if(p < N_POLES) cout <<"in computeorder: rho "<<rho<<" sourceTQ "<<sourceTQ<<" sourceScale "<<sourceScale<<" c "<<c<<" p="<<p<<endl;

  if(p < 1 ) p = 1;
  if(p > N_POLES ) p = N_POLES;

  return p;

  
}

////////////////////
CXFormNIntra::CXFormNIntra(const CSolExpCenter &C1, const CSolExpCenter &C2)
  : CXFormN(C1, C2),   m_qSolvedF1( C1.getQSolvedF()), m_qSolvedF2(C2.getQSolvedF()),
    m_gSolvedF1(C1.getGSolvedF()),m_gSolvedF2(C2.getGSolvedF())
{
  m_tQF1       = &(C1.getTQF()); 
  m_tGF1       = &(C1.getTGF()); 

  m_tQF2       = &(C2.getTQF()); 
  m_tGF2       = &(C2.getTGF()); 
 }
 
//Min not used
void 
CXFormNIntra::xformF(const CMulExpan & Min, CLocalExpan & Lout, bool bFor) 
{
  int p;
  REAL rho = m_P.norm();
  vector<double> qSolved;

  if(bFor)
    {

      if(!m_resetJ )
	{
	  for(int i=0; i<m_SPExSize1; i++) m_d2[i] = m_SPx1[i] - m_P ; 
	  m_resetJ = true;
	}

      p = computeOrder(rho, *m_tQF1, m_scale1);
      Lout = CLocalExpan(m_qSolvedF1, m_d2, p, false, m_scale2);
    }

  // Transform I(2) to J(1)
  else
    {
     if(!m_resetI )
	{
	  for(int i=0; i<m_SPExSize2; i++) m_d1[i] = m_SPx2[i] + m_P ; 
	  m_resetI = true;
	}
      p = computeOrder(rho, *m_tQF2, m_scale2);
      Lout = CLocalExpan(m_qSolvedF2, m_d1, p, false, m_scale1);
    }

  return;
}

//Min not used
void 
CXFormNIntra::xformF(const CTriExpan & Gin, CTriExpan & Gout, bool bFor) 
{
  int p;
  REAL rho = m_P.norm();
  

  // Transform J(1) to I(2)
  if(bFor)
    {
      if(!m_resetJ )
	{
	  for(int i=0; i<m_SPExSize1; i++) m_d2[i] = m_SPx1[i] - m_P ; 
	  m_resetJ = true;
	}

      p = computeOrder(rho, *m_tGF1, m_scale1);
      Gout = CGradExpan(&(m_gSolvedF1[0]), m_d2, p, false, m_scale2, false);
    }

  // Transform I(2) to J(1)
  else
    {
     if(!m_resetI )
	{
	  for(int i=0; i<m_SPExSize2; i++) m_d1[i] = m_SPx2[i] + m_P ; 
	  m_resetI = true;
	}
      p = computeOrder(rho, *m_tGF2, m_scale2);
      Gout = CGradExpan(&(m_gSolvedF2[0]), m_d1, p, false, m_scale1, false);
    }
  
  return;
}

void 
CXFormNIntra::xformFH(const CSolExpCenter &source, const CSolExpCenter &target, 
		     CLocalExpan & LFout, CLocalExpan & LHout, CPnt P) 
{
  xformFH(source, target, LFout, LHout, P, 
	  &(source.getQSolvedF()[0]), &(source.getQSolvedH()[0]), 
	  source.getTQF(),source.getTQH());
} 

void 
CXFormNIntra::xformFH(const CSolExpCenter &source, const CSolExpCenter &target, 
		      CLocalExpan & LFout, CLocalExpan & LHout, CPnt P, 
		      const REAL *qSolvedF, const REAL* qSolvedH, 
		      REAL tQF, REAL tQH) 
{
  const vector<CPnt> &SPx       = source.getSPx();
  int SPExSize                  = source.getSPExSize();
  REAL scaleS                   = source.getLScale();
  REAL scaleT                   = target.getLScale();

  REAL rho = P.norm();
  vector<CPnt> d(SPExSize);
  if(SPExSize != SPx.size() )    { cout<<" spx "<<SPx.size()<<" " <<SPExSize<<endl; exit(1); }

  for(int i=0; i<SPExSize; i++) d[i] = SPx[i] - P ; 

  int p; 
  p = computeOrder(rho, tQF, scaleS);
  LFout = CLocalExpan(qSolvedF, d, p, false, scaleT);
  p = computeOrder(rho, tQH, scaleS);
  LHout = CLocalExpan(qSolvedH, d, p, true,  scaleT);


  return;
}

void 
CXFormNIntra::xformGFH(const CSolExpCenter &source, const CSolExpCenter &target, 
		     CTriExpan & LGFout, CTriExpan & LGHout, CPnt P) 
{
  const CPnt *gSolvedF = &(source.getGSolvedF()[0]);
  const CPnt *gSolvedH = &(source.getGSolvedH()[0]);
  REAL tGF             = source.getTGF();
  REAL tGH             = source.getTGH(); 

  const vector<CPnt> &SPx       = source.getSPx();
  int SPExSize                  = source.getSPExSize();
  REAL scaleS                   = source.getLScale();
  REAL scaleT                   = target.getLScale();
  REAL rho = P.norm();
  vector<CPnt> d(SPExSize);
  // if(SPExSize != SPx.size() )    { cout<<" spx "<<SPx.size()<<" " <<SPExSize<<endl; exit(1); }

  for(int i=0; i<SPExSize; i++) d[i] = SPx[i] - P ; 

  int p; 
  p = computeOrder(rho, tGF, scaleS);
  LGFout = CGradExpan(gSolvedF, d, p, false, scaleT, false); 

  p = computeOrder(rho, tGH, scaleS);
  LGHout = CGradExpan(gSolvedH, d, p, true,  scaleT, false);


  return;
}
/*
/////////////////////////
// reexpand monopole only (works for inter only because of pbc!) 

void 
CXFormBase::xformHrotMono(const CSolExpCenter &source, const CSolExpCenter &target, 
		       CLocalExpan &Mout)
{
  Mout = CLExpan( source.getHrotMono(), 
		  CartToSph( CSystem::pbcPos(target.getCenRot()-source.getCenRot()) ), 
		  target.getbKappa(),
		  1,
		  target.getLScale());

}

// convert Hrot_monopole to local gradient
void 
CXFormBase::xformHrotMonoG(const CSolExpCenter &source, const CSolExpCenter &target, 
			   CTriExpan &Mout)
{
  double Q  = source.getHrotMono();
  xformHrotMonoG(source,target,Mout,Q);
}


// works for inter only (due to pbc)
void 
CXFormBase::xformHrotMonoG(const CSolExpCenter &source, const CSolExpCenter &target, 
			   CTriExpan &Mout, double Q)
{
  CSpPnt sp = CartToSph( CSystem::pbcPos(target.getCenRot()-source.getCenRot()) );
  bool bKappa =  target.getbKappa();
  double scale = target.getLScale();

  double r = sp.rho();
  double fact = -(CMolecule::KAPPA * r + 1) / r;
  double fx, fy, fz;
  if( fabs(sp.theta()) < 1e-5) 
    { fx = 0.0, fy = 1.0/r; fz = 0.0;} 
  else
    {
      REAL sint = sin(sp.theta());
      REAL cost = cos(sp.theta());
      REAL cosp = cos(sp.phi());
      REAL sinp = sin(sp.phi()); 
      fx = sint*cosp; fy = cost * cosp / r; fz = - sint * sinp / r;
    }

  Mout[0] = CLExpan( Q * fact * fx, sp, bKappa,1,scale); 
  Mout[1] = CLExpan( Q * fact * fy, sp, bKappa,1,scale); 
  Mout[2] = CLExpan( Q * fact * fz, sp, bKappa,1,scale); 
}

// convert multipole gradient to local gradient
void 
CXFormBase::xformGHrotMono(const CSolExpCenter &source, const CSolExpCenter &target, 
			   CTriExpan &Mout)
{
  CPnt gH = source.getGHMono();
  xformGHrotMono(source, target, Mout, gH);
}
// convert multipole gradient to local gradient
void 
CXFormBase::xformGHrotMono(const CSolExpCenter &source, const CSolExpCenter &target, 
			   CTriExpan &Mout, CPnt gH)
{
  CSpPnt sp = CartToSph( CSystem::pbcPos(target.getCenRot()-source.getCenRot()) );
  bool bKappa =  target.getbKappa();
  double scale = target.getLScale();

  Mout[0] = CLExpan( gH.x(), sp, bKappa,1,scale); 
  Mout[1] = CLExpan( gH.y(), sp, bKappa,1,scale); 
  Mout[2] = CLExpan( gH.z(), sp, bKappa,1,scale); 
}
*/
