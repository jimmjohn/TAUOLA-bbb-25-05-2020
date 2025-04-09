#include <complex>
using std::complex;

void eepi(const float *pt, const float *pn, const float *pe1, const float *pe2, const float *ppi, float &amplit, float *hv);

void mumupi(const float *pt, const float *pn, const float *pe1, const float *pe2, const float *ppi, float &amplit, float *hv);

void eemu(const float *pt, const float *pn, const float *pe1, const float *pe2, const float *pmu, const float *pnu_mu, float &amplit, float *hv);

void eee(const float *pt, const float *pn, const float *pe1, const float *pe2, const float *pmu, const float *pnu_mu, float &amplit, float *hv);

double nunul(const float *pl,const float *pnu,const float *pnl,float ak0,float *hv);

void demo_mu(const float *pt, const float *pn, const float *pmu, const float *pnu_mu, float &amplit, float *hv);
