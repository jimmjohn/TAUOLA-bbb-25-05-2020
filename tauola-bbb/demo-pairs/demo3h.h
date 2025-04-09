#include <complex>
using std::complex;

float gfun(float s);

void curr3h(complex<float> *hadcur, const float *p1, const float *p2, const float *p3, 
            complex<float> f1, complex<float> f2, complex<float> f3, complex<float> f4, complex<float> f5);
// This function constructs unnormalized hadronic current from given form factors and momentas
// Form factors f1-f5 should already be multiplied by correct (Clebsh-Gordon) coefficients

void curr3pi(const float *p1, const float *p2, const float *p3, complex<float> *hadcur);


/*      FUNCTION FORM1(MNUM,QQ,S1,SDWA)

C
      IF     (MNUM.EQ.0) THEN
C ------------  3 pi hadronic state (a1)
C       FORMRO = FPIKM(SQRT(S1),AMPI,AMPI)
C       FORMRO = F3PI(1,QQ,S1,SDWA)
C       FORMA1 = FA1A1P(QQ)
C       FORM1 = FORMA1*FORMRO
      IF (IVER.EQ.0) THEN
           GAMAX=GAMA1*GFUN(QQ)/GFUN(AMA1**2)
           FORM1=AMA1**2*WIGNER(QQ,AMA1,GAMAX)*FPIKM(SQRT(S1),AMPI,AMPI)
c        FORM1 = F3PI(1,QQ,S1,SDWA)
      ELSE
        FORM1 = F3PI_RCHT(1,QQ,S1,SDWA)
      ENDIF

      ELSEIF (MNUM.EQ.4) THEN
C ------------ pi0 pi0 K-  (K*-pi0)
      XM2   = 1.402
      GAM2  = 0.174
      FORM1 = BWIGM(S1,AMKST,GAMKST,AMK,AMPI)
      FORM1 = WIGFOR(QQ,XM2,GAM2)*FORM1
c       FORMKS = BWIGM(S1,AMKST,GAMKST,AMPI,AMK)
c       FORMK1 = FK1AB(QQ,3)
c       FORM1 = FORMK1*FORMKS

      ELSEIF (MNUM.EQ.5) THEN
C ------------ K- pi- pi+ (rho0 K-)
      XM2   = 1.402
      GAM2  = 0.174
      FORM1 = WIGFOR(QQ,XM2,GAM2)*FPIKM(SQRT(S1),AMPI,AMPI)
c       FORMK1 = FK1AB(QQ,4)
c       FORMRO = FPIKM(SQRT(S1),AMPI,AMPI)
c       FORM1 = FORMK1*FORMRO
      ENDIF
end
*/
//      FUNCTION FORM5(MNUM,QQ,S1,S2)
//C ------------ K- pi- pi+
//        ELPHA=-0.2
//        FORM5=BWIGM(QQ,AMKST,GAMKST,AMPI,AMK)/(1+ELPHA)
//     $       *(       FPIKM(SQRT(S1),AMPI,AMPI)
//     $         +ELPHA*BWIGM(S2,AMKST,GAMKST,AMPI,AMK))
//      END
