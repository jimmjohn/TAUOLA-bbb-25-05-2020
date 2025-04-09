#include <math.h>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <complex>
#include "MEutils.h"
#include "demo3h.h"
using namespace std;


float gfun(float s){
    float mpi=0.13957;
    float mro=0.7759;
    float mpiz=0.134976;

    float m29=mpiz*mpiz*9;
    float l=(mro+mpi)*(mro+mpi);
    float GFUN;

    if(s<l){
        GFUN=4.1*(s-m29)*(s-m29)*(s-m29)*(1.-3.3*(s-m29)+5.8*(s-m29)*(s-m29));
    }
    else{
        GFUN=s*(1.623+10.38/s-9.32/(s*s)+0.65/(s*s*s));
    }
return GFUN;
}


void curr3h(complex<float> *hadcur, const float *p1, const float *p2, const float *p3, 
            complex<float> f1, complex<float> f2, complex<float> f3, complex<float> f4, complex<float> f5){
// This function constructs unnormalized hadronic current from given form factors and momentas
// Form factors f1-f5 should already be multiplied by correct (Clebsh-Gordon) coefficients
    float v1[4], v2[4], v3[4], v4[4], v5[4], Q[4], Q2;
    int i;
    for (i=0;i<4;i++){
        Q[i]=p1[i]+p2[i]+p3[i];
    }
      Q2   =Q[3]*Q[3]-Q[2]*Q[2]-Q[1]*Q[1]-Q[0]*Q[0];

    for (i=0;i<4;i++){
      v1[i]= p2[i]-p3[i]-Q[i]*(Q[3]*(p2[3]-p3[3])-Q[0]*(p2[0]-p3[0])-Q[1]*(p2[1]-p3[1])-Q[2]*(p2[2]-p3[2]))/Q2; 
      v2[i]= p3[i]-p1[i]-Q[i]*(Q[3]*(p3[3]-p1[3])-Q[0]*(p3[0]-p1[0])-Q[1]*(p3[1]-p1[1])-Q[2]*(p3[2]-p1[2]))/Q2;
      v3[i]= p1[i]-p2[i]-Q[i]*(Q[3]*(p1[3]-p2[3])-Q[0]*(p1[0]-p2[0])-Q[1]*(p1[1]-p2[1])-Q[2]*(p1[2]-p2[2]))/Q2;
      v4[i]= Q[i]; 
    }

    prod5(p1, p2, p3, v5);

//   F5 usually requires multiplication by -i/4.0/PI^2/FPI^2
//   it should be done outside this function !!!!!!!!!!!!!!!
//   same goes for F4 and multiplication by -i
    for (i=0;i<4;i++){
        hadcur[i]=v1[i]*f1+v2[i]*f2+v3[i]*f3+v4[i]*f4+v5[i]*f5;
    }
return;
}


void curr3pi(const float *p1, const float *p2, const float *p3, complex<float> *hadcur){
    float s1, s2, s3, Q2, Q[4];
    float norm=0.975/0.0933; // cos_cabibo/Fpi
    complex<float> f1, f2, f3, f4, f5;
    for (int i=0;i<4;i++){
        Q[i]=p1[i]+p2[i]+p3[i];
    }
    Q2   =Q[3]*Q[3]-Q[2]*Q[2]-Q[1]*Q[1]-Q[0]*Q[0];
    s1   =(p3[3]+p2[3])*(p3[3]+p2[3])-(p3[2]+p2[2])*(p3[2]+p2[2])-(p3[1]+p2[1])*(p3[1]+p2[1])-(p3[0]+p2[0])*(p3[0]+p2[0]);
    s2   =(p3[3]+p1[3])*(p3[3]+p1[3])-(p3[2]+p1[2])*(p3[2]+p1[2])-(p3[1]+p1[1])*(p3[1]+p1[1])-(p3[0]+p1[0])*(p3[0]+p1[0]);
    s3   =(p1[3]+p2[3])*(p1[3]+p2[3])-(p1[2]+p2[2])*(p1[2]+p2[2])-(p1[1]+p2[1])*(p1[1]+p2[1])-(p1[0]+p2[0])*(p1[0]+p2[0]);

//BaBar model
    float mpi=0.13957;
    float ma1=1.251;
    float gama1=0.599;
    float gamx;
    gamx=gama1*gfun(Q2)/gfun(ma1*ma1);
           
    f1 = float(2.0*sqrt(2.0)/3.0)*(ma1*ma1*wigner(Q2,ma1,gamx)*fpikm(sqrt(s1),mpi,mpi)); 	//CMPLX(COEF(1,MNUM))*FORM1(XMAA**2,XMRO1**2,XMRO2**2)
    f2 = -float(2.0*sqrt(2.0)/3.0)*(ma1*ma1*wigner(Q2,ma1,gamx)*fpikm(sqrt(s2),mpi,mpi)); 	//CMPLX(COEF(2,MNUM))*FORM2(XMAA**2,XMRO2**2,XMRO1**2)
    f3 = complex<float>(0.0,0.0);		//CMPLX(COEF(3,MNUM))*FORM3(XMAA**2,XMRO3**2,XMRO1**2)
    f4 = complex<float>(0.0,0.0);  		//(-1.0*UROJ)*CMPLX(COEF(4,MNUM))*FORM4(XMAA**2,XMRO1**2,XMRO2**2,XMRO3**2)
    f5 = complex<float>(0.0,0.0);  		//(-1.0)*UROJ/4.0/PI**2/FPI**2*CMPLX(COEF(5,MNUM))*FORM5(XMAA**2,XMRO1**2,XMRO2**2)

    f1=f1*norm;
    f2=f2*norm;

    curr3h(hadcur, p1, p2, p3, f1, f2, f3, f4, f5);
return;
}


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
