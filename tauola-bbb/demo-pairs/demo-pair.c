#include <math.h>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <complex>
#include "MEutils.h"
#include "demo-pair.h"
using namespace std;

void eepi(const float *pt, const float *pn, const float *p1, const float *p2, const float *p3, float &amplit, float *hv){
    
    double Mtau=1.777;
    double Mnu=0.0; // 0.01
    double me=0.000511;
    double Gf=0.0000116637;
    double fpi=0.0922;
    double vud=0.9742;

    
float p23[4]; // momentum of e- e+ pair
      p23[3]=p3[3]+p2[3];
      p23[2]=p3[2]+p2[2];
      p23[1]=p3[1]+p2[1];
      p23[0]=p3[0]+p2[0];
float m2ee=p23[3]*p23[3]-p23[2]*p23[2]-p23[1]*p23[1]-p23[0]*p23[0];

float col[4];
      col[3]=p1[3]*p23[3];
      col[2]=p1[2]*p23[2];
      col[1]=p1[1]*p23[1];
      col[0]=p1[0]*p23[0];
float pk=col[3]-col[2]-col[1]-col[0];

float ptk[4];
      ptk[3]=pt[3]*p23[3];
      ptk[2]=pt[2]*p23[2];
      ptk[1]=pt[1]*p23[1];
      ptk[0]=pt[0]*p23[0];
float ptp23=ptk[3]-ptk[2]-ptk[1]-ptk[0]; 

float kpn[4];
      kpn[3]=pn[3]*p23[3];
      kpn[2]=pn[2]*p23[2];
      kpn[1]=pn[1]*p23[1];
      kpn[0]=pn[0]*p23[0];
float pnk=kpn[3]-kpn[2]-kpn[1]-kpn[0];

float ptq[4];
      ptq[3]=pn[3]*pt[3];
      ptq[2]=pn[2]*pt[2];
      ptq[1]=pn[1]*pt[1];
      ptq[0]=pn[0]*pt[0];
float ptpn=ptq[3]-ptq[2]-ptq[1]-ptq[0];


float pp[4];
      pp[3]=p2[3]*p3[3];
      pp[2]=p2[2]*p3[2];
      pp[1]=p2[1]*p3[1];
      pp[0]=p2[0]*p3[0];
float p2p3=pp[3]-pp[2]-pp[1]-pp[0];

float tmp[4];
      tmp[3]=p2[3]*pt[3];
      tmp[2]=p2[2]*pt[2];
      tmp[1]=p2[1]*pt[1];
      tmp[0]=p2[0]*pt[0];
float p2pt=tmp[3]-tmp[2]-tmp[1]-tmp[0];
      tmp[3]=p3[3]*pt[3];
      tmp[2]=p3[2]*pt[2];
      tmp[1]=p3[1]*pt[1];
      tmp[0]=p3[0]*pt[0];
float p3pt=tmp[3]-tmp[2]-tmp[1]-tmp[0];
      tmp[3]=p2[3]*pn[3];
      tmp[2]=p2[2]*pn[2];
      tmp[1]=p2[1]*pn[1];
      tmp[0]=p2[0]*pn[0];
float p2pn=tmp[3]-tmp[2]-tmp[1]-tmp[0];
      tmp[3]=p3[3]*pn[3];
      tmp[2]=p3[2]*pn[2];
      tmp[1]=p3[1]*pn[1];
      tmp[0]=p3[0]*pn[0];
float p3pn=tmp[3]-tmp[2]-tmp[1]-tmp[0];
      tmp[3]=p2[3]*p1[3];
      tmp[2]=p2[2]*p1[2];
      tmp[1]=p2[1]*p1[1];
      tmp[0]=p2[0]*p1[0];
float p2p1=tmp[3]-tmp[2]-tmp[1]-tmp[0];
      tmp[3]=p3[3]*p1[3];
      tmp[2]=p3[2]*p1[2];
      tmp[1]=p3[1]*p1[1];
      tmp[0]=p3[0]*p1[0];
float p3p1=tmp[3]-tmp[2]-tmp[1]-tmp[0];
      tmp[3]=p1[3]*pn[3];
      tmp[2]=p1[2]*pn[2];
      tmp[1]=p1[1]*pn[1];
      tmp[0]=p1[0]*pn[0];
float p1pn=tmp[3]-tmp[2]-tmp[1]-tmp[0];
      tmp[3]=p1[3]*pt[3];
      tmp[2]=p1[2]*pt[2];
      tmp[1]=p1[1]*pt[1];
      tmp[0]=p1[0]*pt[0];
float p1pt=tmp[3]-tmp[2]-tmp[1]-tmp[0];

float ptsq=pt[3]*pt[3]-pt[2]*pt[2]-pt[1]*pt[1]-pt[0]*pt[0];
float p1sq=p1[3]*p1[3]-p1[2]*p1[2]-p1[1]*p1[1]-p1[0]*p1[0];


      
      // eq. 31 of https://arxiv.org/pdf/1306.1732.pdf
float d1,d2;   //denominator parts
      d1=m2ee-2.*ptp23;
      d2=m2ee+2.*pk;
float c1,c2,c3,c4,c5,c6,c7,c8,c9,a;      // 9 contractions needed

      a=(p2p3+me*me); // part of l_munu, used in all contractions
      
      c1= p2pt*p3pn+p3pt*p2pn-p2p3*ptpn-ptpn*a;
      c1=-(c1*2.+a*ptpn)*m2ee/d1/d1;
      
      c2=p2p1*p3pn+p2pn*p3p1-p1pn*a;
      c2=c2*4.*ptp23/d1/d2;
      
      c3=p2pt*p3pn+p2pn*p3pt-ptpn*a;
      c3=c3*4.*ptp23/d1/d1;
      
      c4=2.*p2p3-a;
      c4=-c4*2.*pnk*ptp23/d1/d1;
      
      c5=p2p1*p3pt+p2pt*p3p1-p1pt*a;
      c5=-c5*4.*pnk/d1/d2;
      
      c6=2.*p2pt*p3pt-ptsq*a;
      c6=-c6*4.*pnk/d1/d1;
      
      c7=p2p1*p3pt+p2pt*p3p1-p1pt*a;
      c7=c7*8.*ptpn/d1/d2;
      
      c8=2.*p2p1*p3p1-p1sq*a;
      c8=c8*4.*ptpn/d2/d2;
      
      c9=2.*p2pt*p3pt-ptsq*a;
      c9=c9*4.*ptpn/d1/d1;
      
float ce;
      ce=16.*Gf*Gf*Mtau*Mtau*fpi*fpi*vud*vud*(16.*M_PI*M_PI/137./137.);
      
      amplit=(c1+c2+c3+c4+c5+c6+c7+c8+c9)/m2ee/m2ee;
      amplit =amplit*ce;

hv[0]=1.0;
hv[1]=0.0;
hv[2]=0.0;
hv[3]=0.0;

}

void mumupi(const float *pt, const float *pn, const float *p1, const float *p2, const float *p3, float &amplit, float *hv){
    
    double Mtau=1.777;
    double Mnu=0.0; // 0.01
    double me=0.1056584; 
    double Gf=0.0000116637;
    double fpi=0.0922;
    double vud=0.9742;

    
float p23[4]; // momentum of pair
      p23[3]=p3[3]+p2[3];
      p23[2]=p3[2]+p2[2];
      p23[1]=p3[1]+p2[1];
      p23[0]=p3[0]+p2[0];
float m2ee=p23[3]*p23[3]-p23[2]*p23[2]-p23[1]*p23[1]-p23[0]*p23[0];

float col[4];
      col[3]=p1[3]*p23[3];
      col[2]=p1[2]*p23[2];
      col[1]=p1[1]*p23[1];
      col[0]=p1[0]*p23[0];
float pk=col[3]-col[2]-col[1]-col[0];

float ptk[4];
      ptk[3]=pt[3]*p23[3];
      ptk[2]=pt[2]*p23[2];
      ptk[1]=pt[1]*p23[1];
      ptk[0]=pt[0]*p23[0];
float ptp23=ptk[3]-ptk[2]-ptk[1]-ptk[0]; 

float kpn[4];
      kpn[3]=pn[3]*p23[3];
      kpn[2]=pn[2]*p23[2];
      kpn[1]=pn[1]*p23[1];
      kpn[0]=pn[0]*p23[0];
float pnk=kpn[3]-kpn[2]-kpn[1]-kpn[0];

float ptq[4];
      ptq[3]=pn[3]*pt[3];
      ptq[2]=pn[2]*pt[2];
      ptq[1]=pn[1]*pt[1];
      ptq[0]=pn[0]*pt[0];
float ptpn=ptq[3]-ptq[2]-ptq[1]-ptq[0];


float pp[4];
      pp[3]=p2[3]*p3[3];
      pp[2]=p2[2]*p3[2];
      pp[1]=p2[1]*p3[1];
      pp[0]=p2[0]*p3[0];
float p2p3=pp[3]-pp[2]-pp[1]-pp[0];

float tmp[4];
      tmp[3]=p2[3]*pt[3];
      tmp[2]=p2[2]*pt[2];
      tmp[1]=p2[1]*pt[1];
      tmp[0]=p2[0]*pt[0];
float p2pt=tmp[3]-tmp[2]-tmp[1]-tmp[0];
      tmp[3]=p3[3]*pt[3];
      tmp[2]=p3[2]*pt[2];
      tmp[1]=p3[1]*pt[1];
      tmp[0]=p3[0]*pt[0];
float p3pt=tmp[3]-tmp[2]-tmp[1]-tmp[0];
      tmp[3]=p2[3]*pn[3];
      tmp[2]=p2[2]*pn[2];
      tmp[1]=p2[1]*pn[1];
      tmp[0]=p2[0]*pn[0];
float p2pn=tmp[3]-tmp[2]-tmp[1]-tmp[0];
      tmp[3]=p3[3]*pn[3];
      tmp[2]=p3[2]*pn[2];
      tmp[1]=p3[1]*pn[1];
      tmp[0]=p3[0]*pn[0];
float p3pn=tmp[3]-tmp[2]-tmp[1]-tmp[0];
      tmp[3]=p2[3]*p1[3];
      tmp[2]=p2[2]*p1[2];
      tmp[1]=p2[1]*p1[1];
      tmp[0]=p2[0]*p1[0];
float p2p1=tmp[3]-tmp[2]-tmp[1]-tmp[0];
      tmp[3]=p3[3]*p1[3];
      tmp[2]=p3[2]*p1[2];
      tmp[1]=p3[1]*p1[1];
      tmp[0]=p3[0]*p1[0];
float p3p1=tmp[3]-tmp[2]-tmp[1]-tmp[0];
      tmp[3]=p1[3]*pn[3];
      tmp[2]=p1[2]*pn[2];
      tmp[1]=p1[1]*pn[1];
      tmp[0]=p1[0]*pn[0];
float p1pn=tmp[3]-tmp[2]-tmp[1]-tmp[0];
      tmp[3]=p1[3]*pt[3];
      tmp[2]=p1[2]*pt[2];
      tmp[1]=p1[1]*pt[1];
      tmp[0]=p1[0]*pt[0];
float p1pt=tmp[3]-tmp[2]-tmp[1]-tmp[0];

float ptsq=pt[3]*pt[3]-pt[2]*pt[2]-pt[1]*pt[1]-pt[0]*pt[0];
float p1sq=p1[3]*p1[3]-p1[2]*p1[2]-p1[1]*p1[1]-p1[0]*p1[0];


      
      // eq. 31 of https://arxiv.org/pdf/1306.1732.pdf
float d1,d2;   //denominator parts
      d1=m2ee-2.*ptp23;
      d2=m2ee+2.*pk;
float c1,c2,c3,c4,c5,c6,c7,c8,c9,a;      // 9 contractions needed

      a=(p2p3+me*me); // part of l_munu, used in all contractions
      
      c1= p2pt*p3pn+p3pt*p2pn-p2p3*ptpn-ptpn*a;
      c1=-(c1*2.+a*ptpn)*m2ee/d1/d1;
      
      c2=p2p1*p3pn+p2pn*p3p1-p1pn*a;
      c2=c2*4.*ptp23/d1/d2;
      
      c3=p2pt*p3pn+p2pn*p3pt-ptpn*a;
      c3=c3*4.*ptp23/d1/d1;
      
      c4=2.*p2p3-a;
      c4=-c4*2.*pnk*ptp23/d1/d1;
      
      c5=p2p1*p3pt+p2pt*p3p1-p1pt*a;
      c5=-c5*4.*pnk/d1/d2;
      
      c6=2.*p2pt*p3pt-ptsq*a;
      c6=-c6*4.*pnk/d1/d1;
      
      c7=p2p1*p3pt+p2pt*p3p1-p1pt*a;
      c7=c7*8.*ptpn/d1/d2;
      
      c8=2.*p2p1*p3p1-p1sq*a;
      c8=c8*4.*ptpn/d2/d2;
      
      c9=2.*p2pt*p3pt-ptsq*a;
      c9=c9*4.*ptpn/d1/d1;
      
float ce;
      ce=16.*Gf*Gf*Mtau*Mtau*fpi*fpi*vud*vud*(16.*M_PI*M_PI/137./137.);
      
      amplit=(c1+c2+c3+c4+c5+c6+c7+c8+c9)/m2ee/m2ee;
      amplit =amplit*ce;

hv[0]=1.0;
hv[1]=0.0;
hv[2]=0.0;
hv[3]=0.0;

}

void eemu(const float *pt, const float *pn, const float *p1, const float *p2, const float *p3, const float *p4, float &amplit, float *hv){
//we assume order: neutrino, muon, electron, positron
    double Gf=1.16637E-5;
    double me=0.000511;
    double mmu=0.10566;
    double Mtau=1.777;
    
double p34[4];
      p34[3]=p3[3]+p4[3];
      p34[2]=p3[2]+p4[2];
      p34[1]=p3[1]+p4[1];
      p34[0]=p3[0]+p4[0];
      
//double p234[4];
      float p234[4];
  
      p234[3]=p2[3]+p3[3]+p4[3];
      p234[2]=p2[2]+p3[2]+p4[2];
      p234[1]=p2[1]+p3[1]+p4[1];
      p234[0]=p2[0]+p3[0]+p4[0];
      
double pl[4];
      pl[3]=p2[3]+p1[3]+pn[3];
      pl[2]=p2[2]+p1[2]+pn[2];
      pl[1]=p2[1]+p1[1]+pn[1];
      pl[0]=p2[0]+p1[0]+pn[0];     
      
double pl2=pl[3]*pl[3]-pl[2]*pl[2]-pl[1]*pl[1]-pl[0]*pl[0];

double col[4];
      col[3]=p2[3]*p34[3];
      col[2]=p2[2]*p34[2];
      col[1]=p2[1]*p34[1];
      col[0]=p2[0]*p34[0];
      
double p2sq=p2[3]*p2[3]-p2[2]*p2[2]-p2[1]*p2[1]-p2[0]*p2[0];

double pk=col[3]-col[2]-col[1]-col[0];
double m2ee=p34[3]*p34[3]-p34[2]*p34[2]-p34[1]*p34[1]-p34[0]*p34[0];
double col2=col[3]*col[3]-col[2]*col[2]-col[1]*col[1]-col[0]*col[0];

double ptk[4];
      ptk[3]=pt[3]*p34[3];
      ptk[2]=pt[2]*p34[2];
      ptk[1]=pt[1]*p34[1];
      ptk[0]=pt[0]*p34[0];
      
double ptp34=ptk[3]*ptk[3]-ptk[2]*ptk[2]-ptk[1]*ptk[1]-ptk[0]*ptk[0];
ptp34=ptk[3]-ptk[2]-ptk[1]-ptk[0];
      
double PP[4];
      PP[3]=p2[3]/(pk+m2ee/2.)-pt[3]/(ptp34-m2ee/2.);
      PP[2]=p2[2]/(pk+m2ee/2.)-pt[2]/(ptp34-m2ee/2.);
      PP[1]=p2[1]/(pk+m2ee/2.)-pt[1]/(ptp34-m2ee/2.);
      PP[0]=p2[0]/(pk+m2ee/2.)-pt[0]/(ptp34-m2ee/2.);
      
double PPsq;
      PPsq=PP[3]*PP[3]-PP[2]*PP[2]-PP[1]*PP[1]-PP[0]*PP[0];
      
double p3PP;
      p3PP=p3[3]*PP[3]-p3[2]*PP[2]-p3[1]*PP[1]-p3[0]*PP[0];
      
double p4PP;
      p4PP=p4[3]*PP[3]-p4[2]*PP[2]-p4[1]*PP[1]-p4[0]*PP[0];

      

//amplit = 4*p2sq/(2*pk+m2ee)/(2*pk+m2ee)/m2ee/m2ee; // matrix element squared
//amplit=1/m2ee/m2ee;
      
double ampser,ce,ampser2;
      ce=4.*(16.*M_PI*M_PI/137./137.);
      ampser=(4.*p3PP*p4PP-m2ee*PPsq)/(2.*m2ee*m2ee)*ce;

double ak0=0.001*1.777;
hv[3]=1.0;
hv[0]=0.0;
hv[1]=0.0;
hv[2]=0.0;
//amplit = (float)nunul(p234,pn,p1,ak0,hv);
amplit = (float)nunul(p2,pn,p1,ak0,hv)*ampser;

}


double nunul(const float *pl,const float *pnu,const float *pnl,float ak0,float *hv){
 
double p234[4];
      p234[3]=pl[3]+pnu[3]+pnl[3];
      p234[2]=pl[2]+pnu[2]+pnl[2];
      p234[1]=pl[1]+pnu[1]+pnl[1];
      p234[0]=pl[0]+pnu[0]+pnl[0];
      
    double Mtau=1.777;//sqrt(p234[3]*p234[3]-p234[2]*p234[2]-p234[1]*p234[1]-p234[0]*p234[0]);
//1.777;
    double Mnu=0.0; // 0.01
    double Gf=1.16637E-5;
    double alpha=137.03604;
    double Msq=Mtau*Mtau;
    double pr[4];
    double rxnl[3],rxnu[3],rxl[3];
    pr[3]=Mtau;
    for(int i=0;i<3;i++){
      pr[0]=0.0;
      pr[1]=0.0;
      pr[2]=0.0;
      pr[i]=Mtau;
      rxnl[i]=pr[3]*pnl[3]-pr[2]*pnl[2]-pr[1]*pnl[1]-pr[0]*pnl[0];
      rxnu[i]=pr[3]*pnu[3]-pr[2]*pnu[2]-pr[1]*pnu[1]-pr[0]*pnu[0];
      rxl[i]=pr[3]*pl[3]-pr[2]*pl[2]-pr[1]*pl[1]-pr[0]*pl[0];        
    }
//     QUASI TWO-BODY VARIABLES
    double u0,u3, w3, w0,up, um, wp, wm, yu, yw, eps2, eps, y, al;
    u0=pl[3]/Mtau;
    u3=sqrt(pl[0]*pl[0]+pl[1]*pl[1]+pl[2]*pl[2])/Mtau;
    w3=u3;
    w0=(pnu[3]+pnl[3])/Mtau;  //does this require ee pair energy?
    up=u0+u3;
    um=u0-u3;
    wp=w0+w3;
    wm=w0-w3;
    yu=log(up/um)/2.;
    yw=log(wp/wm)/2.;
    eps2=u0*u0-u3*u3;
    eps=sqrt(eps2);
    y=w0*w0-w3*w3;
    al=ak0/Mtau;
    
//     FORMFACTORS
    double f0,fp,fm,f3;
    f0=2.*u0/u3*(dilog(1.0-(um*wm/(up*wp))) - dilog(1.0-wm/wp)
       + dilog(1.0-um/up) - 2.*yu + 2.*log(up)*(yw+yu) )
       + 1.0/y*(2.*u3*yu + (1.0-eps2-2*y)*log(eps)) + 2.0 - 4.*(u0/u3*yu-1.0)*log(2.*al);
    fp= yu/(2*u3)*(1.+(1.-eps2)/y)+log(eps)/y;
    fm= yu/(2*u3)*(1.-(1.-eps2)/y)-log(eps)/y;
    f3= eps2*(fp+fm)/2.;

//     SCALAR PRODUCTS OF FOUR-MOMENTA      
    double lxnl,lxnu,nuxnl,mxnl,mxnu,mxl;
    lxnu=pl[3]*pnu[3]-pl[2]*pnu[2]-pl[1]*pnu[1]-pl[0]*pnu[0];
    lxnl=pl[3]*pnl[3]-pl[2]*pnl[2]-pl[1]*pnl[1]-pl[0]*pnl[0];
    nuxnl=pnu[3]*pnl[3]-pnu[2]*pnl[2]-pnu[1]*pnl[1]-pnu[0]*pnl[0];
    
    mxnu=Mtau*pnu[3];
    mxnl=Mtau*pnl[3];
    mxl =Mtau*pl[3];
//     DECAY DIFFERENTIAL WIDTH WITHOUT AND WITH POLARIZATION 
    double c3,am3,xm3;
    c3=1./(2.*alpha*M_PI)*64.*Gf*Gf;
//    if(itdkrc==0) c3=0.0;
    xm3=-(f0*lxnu*mxnl+fp*eps2*mxnu*mxnl+fm*lxnu*lxnl+f3*Msq*nuxnl);
    am3=0.0; // xm3*c3;
//    am3=xm3*c3;
//     V-A  AND  V+A COUPLINGS, BUT IN THE BORN PART ONLY
    double brak,gv,ga,born,xm3pol[3],am3pol[3],bornp[3];
    gv=1.;
    ga=-1.;
    brak= (gv+ga)*(gv+ga)*mxl*nuxnl
         +(gv-ga)*(gv-ga)*mxnl*lxnu
         -(gv*gv-ga*ga)*Mtau*Mnu*lxnl;
    born= 32.*(Gf*Gf/2.0)*brak;
    for(int i=0;i<3;i++){
        xm3pol[i]= -(f0*lxnu*rxnl[i]+fp*eps2*mxnu*rxnl[i] 
                     +fm*lxnu*(lxnl+(rxnl[i]*mxl-mxnl*rxl[i])/Msq)
                     +f3*(Msq*nuxnl+mxnu*rxnl[i]-rxnu[i]*mxnl));
        am3pol[i]=0.0; // xm3pol[i]*c3;
//     V-A  AND  V+A COUPLINGS, BUT IN THE BORN PART ONLY
        bornp[i]=born+((gv+ga)*(gv+ga)*Mtau*nuxnl*pl[i]
                     -(gv+ga)*(gv+ga)*Mtau*lxnu*pnu[i]
                     +(gv*gv-ga*ga)*Mnu*mxnl*pl[i]
                     -(gv*gv-ga*ga)*Mnu*mxl*pnu[i])*32.*(Gf*Gf/2.0);
        hv[i]=(bornp[i]+am3pol[i])/(born+am3)-1.0;
    }
    
    double tbh;
    tbh=born+am3;
    if(isnan(tbh)) {
        cout<<f0<<yw<<wm<<wp<<" "<<w0<<" "<<w3<<endl;
        exit(-1);
    }
    if(tbh/born < 0.1){
        cout<<"ERROR in nunul (matrix element)"<<tbh/born<<endl;
        tbh=0.0;
        return tbh;
    }
//    if(tbh!=0.0) cout<<"amplit= "<< tbh << endl;
    return tbh;
}

void demo_mu(const float *pt, const float *pn, const float *pnu_mu, const float *pmu, float &amplit, float *hv){
    double ak0=0.001*1.777;
    hv[3]=1.0;
    amplit=nunul(pmu,pn,pnu_mu,ak0,hv);
}
