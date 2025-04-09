#include <cstdlib>
#include <cstdio>
#include <complex>
#include "../tauola-c/ChannelForTauolaInterface.h"
#include "../tauola-c/ChannelForTauola.h"
#include "pipi0.h" //demo 2pi current+ME
#include "demo3h.h" // demo 3pi babar
#include "demo-pair.h" // pair emision for leptonic and pion decay
using namespace std;
using Tauolapp::ChannelForTauola;

void lfv_dam2pi(const float *pt, const float *pn, const float *pim1, const float *pim2, float &amplit, float *hv);

// Examples for channel (re-)initialization
void RedefExample() {
    vector<int> products(9);

    
    //**************************************************************************
    // demo pi- pi0 current, it's a copy fortran current using pipi0_curr      *
    //**************************************************************************
/*    
    products[0] = -1; //pi-
    products[1] =  2; //pi0 
    ChannelForTauola *pipi0_test1 = new ChannelForTauola(9.0, products, "  demo pi- pi0  " , pipi0_curr );
    Tauolapp::RegisterChannel( 86, pipi0_test1 );

    
    //**************************************************************************
    // demo pi- pi- pi+ current, it's a copy fortran current used by BaBar     *
    //**************************************************************************
    products[0] = -1; //pi-
    products[1] = -1; //pi-
    products[2] =  1; //pi+
    ChannelForTauola *demo3pi = new ChannelForTauola(9.0, products, "  demo pi- pi- pi+  " , curr3pi );    
    Tauolapp::RegisterChannel( 76, demo3pi );
*/

    //**************************************************************************
    // demo pi- e- e+ ME, based on paper arxiv.org/abs/1306.1732               *
    //**************************************************************************
    products[2] = -11; //e-
    products[1] =  11; //e+
    products[0] = -1; //pi-
    ChannelForTauola *piee = new ChannelForTauola(10.0, products, "  demo ee pi  " , eepi );
    Tauolapp::RegisterChannel( 85, piee );
// Optimizing presampler is important here
    SetPresampler3( piee, 0.0, 0.8, 0.16, 0.1, 0.7759, 0.1479, 0.001, 0.005);
    
    //**************************************************************************
    // demo pi- mu- mu+ ME, based on paper arxiv.org/abs/1306.1732             *
    //**************************************************************************
    products[2] = -13; //mu-
    products[1] =  13; //mu+
    products[0] = -1; //pi-
    ChannelForTauola *pimumu = new ChannelForTauola(10.0, products, "  demo mumu pi  " , mumupi );
    Tauolapp::RegisterChannel( 84, pimumu );
    SetPresampler3( pimumu, 0.0, 0.8, 0.7, 0.4, 0.7759, 0.1479, 0.27, 0.2);
    
    
    //**************************************************************************
    // demo mu nu_mu ee nu_tau ME, based on paper by S.Antropov                *
    //**************************************************************************
//Check order on matrix element, and in phase space calculation
    products[3] = -11; //e+
    products[2] = 11; //e-
//order is important!
    products[1] = -13; //mu-  as second due to phase space generating mass of three (muee)
    products[0] =  14; //nu_mu 
    ChannelForTauola *muee = new ChannelForTauola(10.0, products, "  demo ee mu  " , eemu );
    Tauolapp::RegisterChannel( 30, muee );
    SetPresampler4( muee, 0.0, 0.8, 1.0, 0.9, 0.15, 0.1);
    
    
    //**************************************************************************
    // demo e nu_e ee nu_tau   ME, based on paper by S.Antropov                *
    //**************************************************************************
//Check order on matrix element, and in phase space calculation
    products[3] = -11;//e-
    products[2] =  11; //e+
//order is important!
    products[1] = -11; //e-  as second due to phase space generating mass of three (eee)
    products[0] =  12; //nu_e
    ChannelForTauola *eee = new ChannelForTauola(10.0, products, "  demo ee e  " , eemu );
    Tauolapp::RegisterChannel( 29, eee ); //Functional form of eee is the same as for muee
    SetPresampler4( eee, 0.0, 0.8, 1.0, 0.9, 0.03, 0.07);  //it is hard to balance generation efficiency and number of overweighted events

}


// Set pointer for user channel redefinition
// This function is called from taumain.f
extern "C" void tauolaredef_() {

    // uncomment to run examples for re-initialization from C++:
    Tauolapp::SetUserRedefinitions(RedefExample);      // basic examples 
    //Tauolapp::SetUserRedefinitions(TestCommunication); //more technical examples
}
