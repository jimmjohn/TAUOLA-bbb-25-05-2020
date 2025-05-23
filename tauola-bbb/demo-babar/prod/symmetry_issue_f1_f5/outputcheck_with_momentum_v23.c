// 27.03.2014 JZ

#include <fstream>
#include <iostream>
#include <TVector3.h>
#include <TH1D.h>
#include <TMath.h>
#include <sstream>
#include <string>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector>

using namespace std;

//int main(int argc, char **argv)

void outputcheck()
 {
  fstream in[6];
  fstream out;

  in[0].open("../../tauola_test23_taum.output", ios::in);   //tau -
  in[1].open("../../tauola_test23_taul.output", ios::in);   //tau +



  TH1D *h1[6];
  TH1D *resonance[6];
  TH2D *h1_resonance[6];
  for(int i=0; i<6; i++) {
    h1[i] = new TH1D(Form("h1_%d", i), Form("h1_%d", i), 100, -1., 1. );
    resonance[i] = new TH1D(Form("resonance_%d", i), Form("resonance_%d", i), 400, 0., 2. );
    h1_resonance[i] = new TH2D(Form("h1_resonance_%d", i), Form("h1_resonance_%d", i), 400, 0., 2., 100, -1., 1. );
  }

for(int ij=0; ij<6; ij++) {

  if(ij==2 || ij==3 || ij==5){continue;} //skip tau+ and tau- for K K K

  if (in[ij].good() == true) cout<<"file open"<<endl;
  else cout<<"file is missing or unreadable"<<endl;

  string dane;
  int event_count=0;

  int id[12], parent[12], daughter[12], status[12], np1[12], np2[12];
  double px[12], py[12], pz[12], energy[12], mass[12];
  std::string name[12];

  TVector3 p1P, p2P, p3P, pnuP, pWP, p1N, p2N, p3N, pnuN, pWN , pTauP, pTauN;

  int print=0;

  while (in[ij].fail() == false)
  {
   getline(in[ij], dane);
   if (dane == "                            Event listing (standard)")
      {
        for (int i=1; i<=5; i++) getline(in[ij], dane);
//        for (int i=1; i<=18; i++) {out.write(& dane[0], dane.length()); out<<"\n"; getline(in, dane);}

        for(int i=0;i<7;i++) {
          getline(in[ij], dane);
          std::istringstream pi(dane);
          pi >> id[i] >> name[i] >> parent[i] >> daughter[i] >> status[i] >> np1[i] >> np2[i] >> px[i] >> py[i] >> pz[i] >> energy[i] >> mass[i];
          out << id[i] << " " << name[i] << " " << parent[i] << " " << daughter[i] << " " << status[i] << " " << np1[i] << " " << np2[i] << " "
          << px[i] << " " << py[i] << " " << pz[i] << " " << energy[i] << " " << mass[i] << "\n";
          if(daughter[i]==22){i--;}
        }

//         if(print<10) {
//           cout<<id[0]<<" "<<name[0]<<" "<<parent[0]<<" "<<daughter[0]<<" "<<status[0]<<" "<<np1[0]<<" "<<np2[0]<<" "
//           <<px[0]<<" "<<py[0]<<" "<<pz[0]<<" "<<energy[0]<<" "<<mass[0]<<endl;
//           cout<<id[1]<<" "<<name[1]<<" "<<parent[1]<<" "<<daughter[1]<<" "<<status[1]<<" "<<np1[1]<<" "<<np2[1]<<" "
//           <<px[1]<<" "<<py[1]<<" "<<pz[1]<<" "<<energy[1]<<" "<<mass[1]<<endl;
//           cout<<id[2]<<" "<<name[2]<<" "<<parent[2]<<" "<<daughter[2]<<" "<<status[2]<<" "<<np1[2]<<" "<<np2[2]<<" "
//           <<px[2]<<" "<<py[2]<<" "<<pz[2]<<" "<<energy[2]<<" "<<mass[2]<<endl;
//           cout<<id[3]<<" "<<name[3]<<" "<<parent[3]<<" "<<daughter[3]<<" "<<status[3]<<" "<<np1[3]<<" "<<np2[3]<<" "
//           <<px[3]<<" "<<py[3]<<" "<<pz[3]<<" "<<energy[3]<<" "<<mass[3]<<endl;
//           cout<<id[4]<<" "<<name[4]<<" "<<parent[4]<<" "<<daughter[4]<<" "<<status[4]<<" "<<np1[4]<<" "<<np2[4]<<" "
//           <<px[4]<<" "<<py[4]<<" "<<pz[4]<<" "<<energy[4]<<" "<<mass[4]<<endl;
//           cout<<id[5]<<" "<<name[5]<<" "<<parent[5]<<" "<<daughter[5]<<" "<<status[5]<<" "<<np1[5]<<" "<<np2[5]<<" "
//           <<px[5]<<" "<<py[5]<<" "<<pz[5]<<" "<<energy[5]<<" "<<mass[5]<<endl;
//           cout<<id[6]<<" "<<name[6]<<" "<<parent[6]<<" "<<daughter[6]<<" "<<status[6]<<" "<<np1[6]<<" "<<np2[6]<<" "
//           <<px[6]<<" "<<py[6]<<" "<<pz[6]<<" "<<energy[6]<<" "<<mass[6]<<endl;
//           cout<<id[7]<<" "<<name[7]<<" "<<parent[7]<<" "<<daughter[7]<<" "<<status[7]<<" "<<np1[7]<<" "<<np2[7]<<" "
//           <<px[7]<<" "<<py[7]<<" "<<pz[7]<<" "<<energy[7]<<" "<<mass[7]<<endl;
//           cout<<id[8]<<" "<<name[8]<<" "<<parent[8]<<" "<<daughter[8]<<" "<<status[8]<<" "<<np1[8]<<" "<<np2[8]<<" "
//           <<px[8]<<" "<<py[8]<<" "<<pz[8]<<" "<<energy[8]<<" "<<mass[8]<<endl;
//           cout<<id[9]<<" "<<name[9]<<" "<<parent[9]<<" "<<daughter[9]<<" "<<status[9]<<" "<<np1[9]<<" "<<np2[9]<<" "
//           <<px[9]<<" "<<py[9]<<" "<<pz[9]<<" "<<energy[9]<<" "<<mass[9]<<endl;

//           cout<<"------------------------------------------------------------------------------------------------------"<<endl;

//           print++;
//         }
//  //       What I read is correct - checked


        pTauP.SetXYZ(px[0], py[0], pz[0]);
        pTauN.SetXYZ(px[1], py[1], pz[1]);
        pnuP.SetXYZ(px[2], py[2], pz[2]);
        pWP.SetXYZ(px[3], py[3], pz[3]);
        p1P.SetXYZ(px[4], py[4], pz[4]);
        p2P.SetXYZ(px[5], py[5], pz[5]);
        p3P.SetXYZ(px[6], py[6], pz[6]);


        TLorentzVector vTauP(pTauP, energy[0]);
        TLorentzVector vTauN(pTauN, energy[1]);
        TLorentzVector vnuP(pnuP, energy[2]);
        TLorentzVector vWP(pWP, energy[3]);
        TLorentzVector v1P(p1P, energy[4]);
        TLorentzVector v2P(p2P, energy[5]);
        TLorentzVector v3P(p3P, energy[6]);


        TVector3 tauPboost_direction(-pTauP); // boost direction - opposite of hadronic system
        tauPboost_direction = tauPboost_direction.Unit();
        TVector3 tauNboost_direction(-pTauN); // boost direction - opposite of hadronic system
        tauNboost_direction = tauNboost_direction.Unit();

        TVector3 tauBoostP;
        if(ij==0){
          tauBoostP = tauNboost_direction * (pTauP.Mag() / energy[0]);
        } else if(ij==1){
          tauBoostP = tauNboost_direction * (pTauN.Mag() / energy[1]);
        } else {
          cout<<"Error in the boost direction"<<endl;
        }

        // Apply the boost to tau rest frame
        vTauP.Boost(tauBoostP);
        vTauN.Boost(tauBoostP);
        v1P.Boost(tauBoostP);
        v2P.Boost(tauBoostP);
        v3P.Boost(tauBoostP);
        vnuP.Boost(tauBoostP);
        vWP.Boost(tauBoostP);

         if(print<2 && ij==0){
           cout<<"Point1"<<endl;
           if(ij==0){
            cout<<"Momenta="<<vTauN.Px()<<" "<<vTauN.Py()<<" "<<vTauN.Pz()<<endl;
           } else if(ij==1){
            cout<<"Momenta="<<vTauN.Px()<<" "<<vTauN.Py()<<" "<<vTauN.Pz()<<endl;
           }
           cout<<"Momenta_vWP="<<vWP.Px()<<" "<<vWP.Py()<<" "<<vWP.Pz()<<endl;
           cout<<"Momenta_vnuP="<<vnuP.Px()<<" "<<vnuP.Py()<<" "<<vnuP.Pz()<<endl;
           print++;
         }


        double M2P = sqrt((pow((energy[5]+energy[6]),2)-pow((px[5]+px[6]),2)-pow((py[5]+py[6]),2)-pow((pz[5]+pz[6]),2)));

        // double modP1 = abs(pow(0.77590,2)-(pow(energy[2]+energy[4],2) - pow(px[2]+px[4],2) - pow(py[2]+py[4],2) - pow(pz[2]+pz[4],2)));
        // double modP2 = abs(pow(0.77590,2)-(pow(energy[3]+energy[4],2) - pow(px[3]+px[4],2) - pow(py[3]+py[4],2) - pow(pz[3]+pz[4],2)));
        // double modP3 = abs(pow(0.77590,2)-(pow(energy[2]+energy[3],2) - pow(px[2]+px[3],2) - pow(py[2]+py[3],2) - pow(pz[2]+pz[3],2)));

        double mPiSystemP = sqrt((pow(energy[3],2)-pow(px[3],2)-pow(py[3],2)-pow(pz[3],2)));

        resonance[ij+1]->Fill(mPiSystemP);

        // Define the boost direction (opposite direction of this vector)
        TVector3 boost_directionP = vnuP.Vect(); // boost direction - opposite of hadronic system
        boost_directionP = boost_directionP.Unit();

        double pnuenergyP = energy[2];
        double WenergyP = energy[3];

        // Compute the boost velocity vector v = p / E
        TVector3 boostP = boost_directionP * (vWP.P() / vWP.E());

        // Apply the boost
        v1P.Boost(boostP);
        v2P.Boost(boostP);
        v3P.Boost(boostP);
        //vnuP.Boost(boostP);
        vWP.Boost(boostP);



        if(print<2 && ij==0){
          cout<<"Point2"<<endl;
          cout<<"Momenta="<<vWP.Px()<<" "<<vWP.Py()<<" "<<vWP.Pz()<<endl;
        }


        // double modP1 = abs(sqrt(v1.Px()*v1.Px() + v1.Py()*v1.Py() + v1.Pz()*v1.Pz()));
        // double modP2 = abs(sqrt(v2.Px()*v2.Px() + v2.Py()*v2.Py() + v2.Pz()*v2.Pz()));
        // double modP3 = abs(sqrt(v3.Px()*v3.Px() + v3.Py()*v3.Py() + v3.Pz()*v3.Pz()));
        // double modPnu = sqrt(vnu.Px()*vnu.Px() + vnu.Py()*vnu.Py() + vnu.Pz()*vnu.Pz());
        // double modPW = sqrt(vW.Px()*vW.Px() + vW.Py()*vW.Py() + vW.Pz()*vW.Pz());


        //Modification suggestb by Prof. Was
        double modP1P = abs(pow(0.77590,2)-(pow(v1P.Energy()+v3P.Energy(),2) - pow(v1P.Px()+v3P.Px(),2) - pow(v1P.Py()+v3P.Py(),2) - pow(v1P.Pz()+v3P.Pz(),2)));
        double modP2P = abs(pow(0.77590,2)-(pow(v2P.Energy()+v3P.Energy(),2) - pow(v2P.Px()+v3P.Px(),2) - pow(v2P.Py()+v3P.Py(),2) - pow(v2P.Pz()+v3P.Pz(),2)));
        double modP3P = abs(pow(0.77590,2)-(pow(v1P.Energy()+v2P.Energy(),2) - pow(v1P.Px()+v2P.Px(),2) - pow(v1P.Py()+v2P.Py(),2) - pow(v1P.Pz()+v2P.Pz(),2)));
        double modPnuP = sqrt(vnuP.Px()*vnuP.Px() + vnuP.Py()*vnuP.Py() + vnuP.Pz()*vnuP.Pz());
        double modPWP = sqrt(vWP.Px()*vWP.Px() + vWP.Py()*vWP.Py() + vWP.Pz()*vWP.Pz());


        TVector3 pi1P(v1P.Px(), v1P.Py(), v1P.Pz());
        TVector3 pi2P(v2P.Px(), v2P.Py(), v2P.Pz());
        TVector3 pi3P(v3P.Px(), v3P.Py(), v3P.Pz());


        if(ij<4) {
         if(modP1P<modP2P) {
            TVector3 n_perpendicularP;
            n_perpendicularP = pi1P.Cross(pi2P);
            double betaP = n_perpendicularP.Angle(vnuP.Vect());
            h1[ij]->Fill(cos(betaP));
            h1_resonance[ij+1]->Fill(mPiSystemP, cos(betaP));
         } else {
            TVector3 n_perpendicularP;
            n_perpendicularP = pi2P.Cross(pi1P);
            double betaP = n_perpendicularP.Angle(vnuP.Vect());
            h1[ij]->Fill(cos(betaP));
            h1_resonance[ij+1]->Fill(mPiSystemP, cos(betaP));
         }
        } else {

        }


        //while (dane != "") {out.write(& dane[0], dane.length()); out<<"\n"; getline(in, dane);}
        out<<"\n";
        event_count++;
      }
  }

  in[ij].close();

  TVector3 t1, t2, t3;
  t1.SetXYZ(0,0,1);
  t2.SetXYZ(0,1,0);
  //cout<<"Angle=" <<t1.Angle(t2)<<endl;


}




  TCanvas *beta_dist1 = new TCanvas("beta_distribution1", "beta_distribution1", 800, 600);
  beta_dist1->cd();
  h1[0]->SetLineColor(kBlack);
  h1[0]->SetTitle("Beta distribution");
  h1[0]->GetXaxis()->SetTitle("cos(#beta)");
  h1[0]->GetYaxis()->SetTitle("Events");
  h1[0]->GetXaxis()->CenterTitle();
  h1[0]->GetYaxis()->CenterTitle();
  h1[0]->GetXaxis()->SetRangeUser(-1., 1.);
  h1[0]->GetYaxis()->SetRangeUser(0., 1400.);
  h1[0]->SetStats(0);
  h1[0]->SetLineWidth(2);
  h1[0]->SetMarkerStyle(1);
  h1[0]->Draw();
  h1[1]->SetLineColor(kRed);
  h1[1]->SetLineWidth(2);
  h1[1]->SetLineStyle(2);
  h1[1]->SetMarkerStyle(7);
  h1[1]->Draw("same");
  TLatex *latex1 = new TLatex();
  latex1->SetTextSize(0.05);
  latex1->SetTextColor(kBlack);
  latex1->SetTextAlign(12);
  //latex1->DrawLatex(-0.4, 200, "F1=CONJG(F1) && F2=CONJG(F2)");
  TLegend *legend1 = new TLegend(0.3, 0.35, 0.6, 0.5);  // Adjust position as needed
  legend1->AddEntry(h1[0], "#pi^{-} #pi^{-} #pi^{+}", "l");  // "l" for line, "p" for point, etc.
  legend1->AddEntry(h1[1], "#pi^{+} #pi^{+} #pi^{-}", "l");
  legend1->Draw();
  beta_dist1->SaveAs("beta_distribution_pi_pi_pi_2.png");


  TCanvas *beta_dist2 = new TCanvas("beta_distribution2", "beta_distribution2", 800, 600);
  beta_dist2->cd();
  h1[2]->SetLineColor(kOrange);
  h1[2]->SetTitle("Beta distribution");
  h1[2]->GetXaxis()->SetTitle("cos(#beta)");
  h1[2]->GetYaxis()->SetTitle("Events");
  h1[2]->GetXaxis()->CenterTitle();
  h1[2]->GetYaxis()->CenterTitle();
  h1[2]->GetXaxis()->SetRangeUser(-1., 1.);
  h1[2]->GetYaxis()->SetRangeUser(0., 1500.);
  h1[2]->SetStats(0);
  h1[2]->SetLineWidth(2);
  h1[2]->SetMarkerStyle(1);
  h1[2]->Draw();
  h1[3]->SetLineColor(kMagenta);
  h1[3]->SetLineWidth(2);
  h1[3]->SetLineStyle(2);
  h1[3]->SetMarkerStyle(7);
  h1[3]->Draw("same");
  TLatex *latex2 = new TLatex();
  latex2->SetTextSize(0.05);
  latex2->SetTextColor(kBlack);
  latex2->SetTextAlign(12);
 // latex2->DrawLatex(-0.4, 200, "F1=CONJG(F1) && F2=CONJG(F2)");
  TLegend *legend2 = new TLegend(0.3, 0.35, 0.6, 0.5);  // Adjust position as needed
  legend2->AddEntry(h1[2], "#pi^{0} #pi^{0} #pi^{-}", "l");
  legend2->AddEntry(h1[3], "#pi^{0} #pi^{0} #pi^{+}", "l");
  legend2->Draw();
  beta_dist2->SaveAs("beta_distribution_pi0_pi0_pi_2.png");

  TCanvas *beta_dist3 = new TCanvas("beta_distribution3", "beta_distribution3", 800, 600);
  beta_dist3->cd();
  h1[4]->SetLineColor(kGreen);
  h1[4]->SetTitle("Beta distribution");
  h1[4]->GetXaxis()->SetTitle("cos(#beta)");
  h1[4]->GetYaxis()->SetTitle("Events");
  h1[4]->GetXaxis()->CenterTitle();
  h1[4]->GetYaxis()->CenterTitle();
  h1[4]->GetXaxis()->SetRangeUser(-1., 1.);
  h1[4]->GetYaxis()->SetRangeUser(0., 1200.);
  h1[4]->SetStats(0);
  h1[4]->SetLineWidth(2);
  h1[4]->SetMarkerStyle(1);
  h1[4]->Draw();
  h1[5]->SetLineColor(kBlue);
  h1[5]->SetLineWidth(2);
  h1[5]->SetLineStyle(2);
  h1[5]->SetMarkerStyle(7);
  h1[5]->Draw("same");
  TLatex *latex3 = new TLatex();
  latex3->SetTextSize(0.05);
  latex3->SetTextColor(kBlack);
  latex3->SetTextAlign(12);
  //latex3->DrawLatex(-0.4, 200, "F1=CONJG(F1) && F2=CONJG(F2)");
  TLegend *legend3 = new TLegend(0.3, 0.35, 0.6, 0.5);  // Adjust position as needed
  legend3->AddEntry(h1[4], "K^{-} #pi^{-} #pi^{+}", "l");
  legend3->AddEntry(h1[5], "K^{+} #pi^{+} #pi^{-}", "l");
  legend3->Draw();
  beta_dist3->SaveAs("beta_distribution_K_pi_pi_2.png");


 TCanvas *resonance_dist = new TCanvas("resonance_distribution", "resonance_distribution", 800, 600);
  resonance[0]->SetLineColor(kBlack);
  resonance[0]->SetTitle("Resonance distribution");
  resonance[0]->GetXaxis()->SetTitle("M2");
  resonance[0]->GetXaxis()->CenterTitle();
  resonance[0]->GetYaxis()->SetTitle("Events");
  resonance[0]->GetYaxis()->CenterTitle();
  resonance[0]->GetXaxis()->SetRangeUser(0., 3.);
  resonance[0]->GetYaxis()->SetRangeUser(0., 5000.);
  resonance[0]->SetStats(0);
  resonance[0]->SetLineWidth(2);
  resonance[0]->SetMarkerStyle(1);
  resonance[0]->Draw();
  resonance[1]->SetLineColor(kRed);
  resonance[1]->SetLineWidth(2);
  resonance[1]->SetLineStyle(2);
  resonance[1]->SetMarkerStyle(7);
  resonance[1]->Draw("same");
  resonance[2]->SetLineColor(kGreen);
  resonance[2]->SetLineWidth(2);
  resonance[2]->SetMarkerStyle(1);
  resonance[2]->Draw("same");
  resonance[3]->SetLineColor(kBlue);
  resonance[3]->SetLineWidth(2);
  resonance[3]->SetLineStyle(2);
  resonance[3]->SetMarkerStyle(7);
  resonance[3]->Draw("same");
  resonance[4]->SetLineColor(kMagenta);
  resonance[4]->SetLineWidth(2);
  resonance[4]->SetMarkerStyle(1);
  resonance[4]->Draw("same");
  resonance[5]->SetLineColor(kCyan);
  resonance[5]->SetLineWidth(2);
  resonance[5]->SetLineStyle(2);
  resonance[5]->SetMarkerStyle(7);
  resonance[5]->Draw("same");

  TCanvas *h1_reso1 = new TCanvas("y", "y", 800, 600);
  h1_reso1->cd();
  h1_resonance[0]->Draw("COLZ");

  TCanvas *h1_reso2 = new TCanvas("x", "x", 800, 600);
  h1_reso2->cd();
  h1_resonance[1]->Draw("COLZ");


  TLegend *legendR = new TLegend(0.15, 0.5, 0.45, 0.8);  // Adjust position as needed
  legendR->AddEntry(resonance[0], "#pi^{-} #pi^{-} #pi^{+}", "l");  // "l" for line, "p" for point, etc.
  legendR->AddEntry(resonance[1], "#pi^{+} #pi^{+} #pi^{-}", "l");
  legendR->AddEntry(resonance[2], "#pi^{0} #pi^{0} #pi^{-}", "l");
  legendR->AddEntry(resonance[3], "#pi^{0} #pi^{0} #pi^{+}", "l");
  legendR->AddEntry(resonance[4], "K^{-} #pi^{-} #pi^{+}", "l");
  legendR->AddEntry(resonance[5], "K^{+} #pi^{+} #pi^{-}", "l");
  legendR->Draw();

  resonance_dist->SaveAs("resonance_distribution.png");


 }
