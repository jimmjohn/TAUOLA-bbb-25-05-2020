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

//commented lines alows, to make another file wih list og generated channel only
// if zero appears as number of efents in channel, that mean it's not initialized properly
/*
  fstream chan;
  chan.open("channels.txt", ios::out);
 */
  in[0].open("../../tauola_pm_pm_pl.output", ios::in);   //tau -
  in[1].open("../../tauola_pl_pl_pm.output", ios::in);   //tau +
  in[2].open("../../tauola_km_km_kl.output", ios::in);   //tau -
  in[3].open("../../tauola_kl_kl_km.output", ios::in);   //tau +
  in[4].open("../../tauola_km_pm_pl.output", ios::in);   //tau -
  in[5].open("../../tauola_kl_pl_pm.output", ios::in);   //tau +

  TH1D *h1[6];
  TH1D *resonance[6];
  for(int i=0; i<6; i++) {
    h1[i] = new TH1D(Form("h1_%d", i), Form("h1_%d", i), 100, -1., 1. );
    resonance[i] = new TH1D(Form("resonance_%d", i), Form("resonance_%d", i), 400, 0., 2. );
  }

for(int ij=0; ij<6; ij++) {

  if (in[ij].good() == true) cout<<"file open"<<endl;
  else cout<<"file is missing or unreadable"<<endl;

  string dane;
  int event_count=0;

  int id[5], parent[5], daughter[5], status[5], np1[5], np2[5];
  double px[5], py[5], pz[5], energy[5], mass[5];
  std::string name[5];

  TVector3 p1, p2, p3, pnu, pW;

  while (in[ij].fail() == false)
  {
   getline(in[ij], dane);
   if (dane == "                            Event listing (standard)")
      {
        for (int i=1; i<=7; i++) getline(in[ij], dane);
//        for (int i=1; i<=18; i++) {out.write(& dane[0], dane.length()); out<<"\n"; getline(in, dane);}

        for(int i=0;i<5;i++) {
          getline(in[ij], dane);
          std::istringstream pi(dane);
          pi >> id[i] >> name[i] >> parent[i] >> daughter[i] >> status[i] >> np1[i] >> np2[i] >> px[i] >> py[i] >> pz[i] >> energy[i] >> mass[i];
          out << id[i] << " " << name[i] << " " << parent[i] << " " << daughter[i] << " " << status[i] << " " << np1[i] << " " << np2[i] << " "
          << px[i] << " " << py[i] << " " << pz[i] << " " << energy[i] << " " << mass[i] << "\n";
        }

        pnu.SetXYZ(px[0], py[0], pz[0]);
        pW.SetXYZ(px[1], py[1], pz[1]);
        p1.SetXYZ(px[2], py[2], pz[2]);
        p2.SetXYZ(px[3], py[3], pz[3]);
        p3.SetXYZ(px[4], py[4], pz[4]);

        TLorentzVector v1(p1, energy[2]);
        TLorentzVector v2(p2, energy[3]);
        TLorentzVector v3(p3, energy[4]);
        TLorentzVector vnu(pnu, energy[0]);
        TLorentzVector vW(pW, energy[1]);

        double M2 = sqrt((pow((energy[3]+energy[4]),2)-pow((px[3]+px[4]),2)-pow((py[3]+py[4]),2)-pow((pz[3]+pz[4]),2)));

        // double modP1 = abs(pow(0.77590,2)-(pow(energy[2]+energy[4],2) - pow(px[2]+px[4],2) - pow(py[2]+py[4],2) - pow(pz[2]+pz[4],2)));
        // double modP2 = abs(pow(0.77590,2)-(pow(energy[3]+energy[4],2) - pow(px[3]+px[4],2) - pow(py[3]+py[4],2) - pow(pz[3]+pz[4],2)));
        // double modP3 = abs(pow(0.77590,2)-(pow(energy[2]+energy[3],2) - pow(px[2]+px[3],2) - pow(py[2]+py[3],2) - pow(pz[2]+pz[3],2)));

        double mPiSystem = sqrt((pow(energy[1],2)-pow(px[1],2)-pow(py[1],2)-pow(pz[1],2)));
          resonance[ij]->Fill(mPiSystem);


        // Define the boost direction (opposite direction of this vector)
        TVector3 boost_direction(pnu); // boost direction - opposite of hadronic system
        boost_direction = boost_direction.Unit();

        double pnuenergy = energy[0];
        double Wenergy = energy[1];

        // Compute the boost velocity vector v = p / E
        TVector3 boost = boost_direction * (pW.Mag() / Wenergy);

        // Apply the boost
        v1.Boost(boost);
        v2.Boost(boost);
        v3.Boost(boost);
        vnu.Boost(boost);
        vW.Boost(boost);

        cout<<"Momenta="<<vW.Px()<<" "<<vW.Py()<<" "<<vW.Pz()<<endl;

        // double modP1 = abs(sqrt(v1.Px()*v1.Px() + v1.Py()*v1.Py() + v1.Pz()*v1.Pz()));
        // double modP2 = abs(sqrt(v2.Px()*v2.Px() + v2.Py()*v2.Py() + v2.Pz()*v2.Pz()));
        // double modP3 = abs(sqrt(v3.Px()*v3.Px() + v3.Py()*v3.Py() + v3.Pz()*v3.Pz()));
        // double modPnu = sqrt(vnu.Px()*vnu.Px() + vnu.Py()*vnu.Py() + vnu.Pz()*vnu.Pz());
        // double modPW = sqrt(vW.Px()*vW.Px() + vW.Py()*vW.Py() + vW.Pz()*vW.Pz());


        //Modification suggestb by Prof. Was
        double modP1 = abs(pow(0.77590,2)-(pow(v1.Energy()+v3.Energy(),2) - pow(v1.Px()+v3.Px(),2) - pow(v1.Py()+v3.Py(),2) - pow(v1.Pz()+v3.Pz(),2)));
        double modP2 = abs(pow(0.77590,2)-(pow(v2.Energy()+v3.Energy(),2) - pow(v2.Px()+v3.Px(),2) - pow(v2.Py()+v3.Py(),2) - pow(v2.Pz()+v3.Pz(),2)));
        double modP3 = abs(pow(0.77590,2)-(pow(v1.Energy()+v2.Energy(),2) - pow(v1.Px()+v2.Px(),2) - pow(v1.Py()+v2.Py(),2) - pow(v1.Pz()+v2.Pz(),2)));
        double modPnu = sqrt(vnu.Px()*vnu.Px() + vnu.Py()*vnu.Py() + vnu.Pz()*vnu.Pz());
        double modPW = sqrt(vW.Px()*vW.Px() + vW.Py()*vW.Py() + vW.Pz()*vW.Pz());


        TVector3 pi1(v1.Px(), v1.Py(), v1.Pz());
        TVector3 pi2(v2.Px(), v2.Py(), v2.Pz());
        TVector3 pi3(v3.Px(), v3.Py(), v3.Pz());


        if(ij<4) {
         if(modP1<modP2) {
           TVector3 n_perpendicular = pi1.Cross(pi2);
           double beta = n_perpendicular.Angle(pnu);
           h1[ij]->Fill(cos(beta));
         } else {
           TVector3 n_perpendicular = pi2.Cross(pi1);
           double beta = n_perpendicular.Angle(pnu);
           h1[ij]->Fill(cos(beta));
         }
        } else {
          TVector3 n_perpendicular = pi1.Cross(pi2);
           double beta = n_perpendicular.Angle(pnu);
           h1[ij]->Fill(cos(beta));
        }


        //while (dane != "") {out.write(& dane[0], dane.length()); out<<"\n"; getline(in, dane);}
        out<<"\n";
        event_count++;
      }
  }

  in[ij].close();


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
  h1[0]->GetYaxis()->SetRangeUser(0., 1500.);
  h1[0]->SetStats(0);
  h1[0]->SetLineWidth(2);
  h1[0]->SetMarkerStyle(1);
  h1[0]->Draw();
  h1[1]->SetLineColor(kRed);
  h1[1]->SetLineWidth(2);
  h1[1]->SetLineStyle(2);
  h1[1]->SetMarkerStyle(7);
  h1[1]->Draw("same");
  TLegend *legend1 = new TLegend(0.3, 0.35, 0.6, 0.5);  // Adjust position as needed
  legend1->AddEntry(h1[0], "#pi^{-} #pi^{-} #pi^{+}", "l");  // "l" for line, "p" for point, etc.
  legend1->AddEntry(h1[1], "#pi^{+} #pi^{+} #pi^{-}", "l");
  legend1->Draw();
  beta_dist1->SaveAs("beta_distribution_pi_pi_pi_2.png");


  TCanvas *beta_dist2 = new TCanvas("beta_distribution2", "beta_distribution2", 800, 600);
  beta_dist2->cd();
  h1[2]->SetLineColor(kGreen);
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
  h1[3]->SetLineColor(kBlue);
  h1[3]->SetLineWidth(2);
  h1[3]->SetLineStyle(2);
  h1[3]->SetMarkerStyle(7);
  h1[3]->Draw("same");
  TLegend *legend2 = new TLegend(0.3, 0.35, 0.6, 0.5);  // Adjust position as needed
  legend2->AddEntry(h1[2], "K^{-} K^{-} K^{+}", "l");
  legend2->AddEntry(h1[3], "K^{+} K^{+} K^{-}", "l");
  legend2->Draw();
  beta_dist2->SaveAs("beta_distribution_K_K_K_2.png");

  TCanvas *beta_dist3 = new TCanvas("beta_distribution3", "beta_distribution3", 800, 600);
  beta_dist3->cd();
  h1[4]->SetLineColor(kGreen);
  h1[4]->SetTitle("Beta distribution");
  h1[4]->GetXaxis()->SetTitle("cos(#beta)");
  h1[4]->GetYaxis()->SetTitle("Events");
  h1[4]->GetXaxis()->CenterTitle();
  h1[4]->GetYaxis()->CenterTitle();
  h1[4]->GetXaxis()->SetRangeUser(-1., 1.);
  h1[4]->GetYaxis()->SetRangeUser(0., 1500.);
  h1[4]->SetStats(0);
  h1[4]->SetLineWidth(2);
  h1[4]->SetMarkerStyle(1);
  h1[4]->Draw();
  h1[5]->SetLineColor(kBlue);
  h1[5]->SetLineWidth(2);
  h1[5]->SetLineStyle(2);
  h1[5]->SetMarkerStyle(7);
  h1[5]->Draw("same");
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

  TLegend *legendR = new TLegend(0.15, 0.5, 0.45, 0.8);  // Adjust position as needed
  legendR->AddEntry(resonance[0], "#pi^{-} #pi^{-} #pi^{+}", "l");  // "l" for line, "p" for point, etc.
  legendR->AddEntry(resonance[1], "#pi^{+} #pi^{+} #pi^{-}", "l");
  legendR->AddEntry(resonance[2], "K^{-} K^{-} K^{+}", "l");
  legendR->AddEntry(resonance[3], "K^{+} K^{+} K^{-}", "l");
  legendR->AddEntry(resonance[4], "K^{-} #pi^{-} #pi^{+}", "l");
  legendR->AddEntry(resonance[5], "K^{+} #pi^{+} #pi^{-}", "l");
  legendR->Draw();

  resonance_dist->SaveAs("resonance_distribution.png");


 }
