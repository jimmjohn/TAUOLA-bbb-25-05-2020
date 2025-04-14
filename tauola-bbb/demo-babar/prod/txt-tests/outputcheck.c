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
  fstream in;
  fstream out;

//commented lines alows, to make another file wih list og generated channel only
// if zero appears as number of efents in channel, that mean it's not initialized properly
/*
  fstream chan;
  chan.open("channels.txt", ios::out);
 */
  in.open("../../tauola.output", ios::in);
  out.open("out.txt",	 ios::out);
  if (in.good() == true) cout<<"file open"<<endl;
  else cout<<"file is missing or unreadable"<<endl;
  if (out.good() == false) cout<<"cannot create file for writing"<<endl;
  string dane;
  int event_count=0;

  TVector3 p1, p2, p3, p4, pnu, pW;

  TH1D *h1 = new TH1D("beta", "beta", 100, -1.1, 1.1 );

  TH1D *p1x = new TH1D("p1x", "p1x", 40, -0.8, -0.6 );
  TH1D *p1y = new TH1D("p1y", "p1y", 40, -0.8, -0.6 );
  TH1D *p1z = new TH1D("p1z", "p1z", 40, -0.8, -0.6 );


  TH1D *p2x = new TH1D("p2x", "p2x", 40, -0.8, -0.6 );
  TH1D *p2y = new TH1D("p2y", "p2y", 40, -0.8, -0.6 );
  TH1D *p2z = new TH1D("p2z", "p2z", 40, -0.8, -0.6 );

  TH1D *p3x = new TH1D("p3x", "p3x", 40, -0.8, -0.6 );
  TH1D *p3y = new TH1D("p3y", "p3y", 40, -0.8, -0.6 );
  TH1D *p3z = new TH1D("p3z", "p3z", 40, -0.8, -0.6 );


  TH1D *p4x = new TH1D("p4x", "p4x", 40, -0.8, -0.6 );
  TH1D *p4y = new TH1D("p4y", "p4y", 40, -0.8, -0.6 );
  TH1D *p4z = new TH1D("p4z", "p4z", 40, -0.8, -0.6 );

  TH1D *resonance = new TH1D("resonance", "resonance", 500, 0.3, 1 );


  int id[6], parent[6], daughter[6], status[6], np1[6], np2[6];
  double px[6], py[6], pz[6], energy[6], mass[6];
  std::string name[6];


  while (in.fail() == false)
  {
   getline(in, dane);
   if (dane == "                            Event listing (standard)")
      {
        for (int i=1; i<=7; i++) getline(in, dane);
//        for (int i=1; i<=18; i++) {out.write(& dane[0], dane.length()); out<<"\n"; getline(in, dane);}

        for(int i=0;i<6;i++) {
          getline(in, dane);
          std::istringstream pi(dane);
          pi >> id[i] >> name[i] >> parent[i] >> daughter[i] >> status[i] >> np1[i] >> np2[i] >> px[i] >> py[i] >> pz[i] >> energy[i] >> mass[i];
          out << id[i] << " " << name[i] << " " << parent[i] << " " << daughter[i] << " " << status[i] << " " << np1[i] << " " << np2[i] << " "
          << px[i] << " " << py[i] << " " << pz[i] << " " << energy[i] << " " << mass[i] << "\n";
        }

        //Calculate the cos(beta)

        pnu.SetXYZ(px[0], py[0], pz[0]);
        pW.SetXYZ(px[1], py[1], pz[1]);
        p1.SetXYZ(px[2], py[2], pz[2]);
        p2.SetXYZ(px[3], py[3], pz[3]);
        p3.SetXYZ(px[4], py[4], pz[4]);
        p4.SetXYZ(px[5], py[5], pz[5]);

      // Fill histograms with momenta
      p1x->Fill(p1.Px());
      p1y->Fill(p1.Py());
      p1z->Fill(p1.Pz());
      p2x->Fill(p2.Px());
      p2y->Fill(p2.Py());
      p2z->Fill(p2.Pz());
      p3x->Fill(p3.Px());
      p3y->Fill(p3.Py());
      p3z->Fill(p3.Pz());
      p4x->Fill(p4.Px());
      p4y->Fill(p4.Py());
      p4z->Fill(p4.Pz());

        // // Fill histograms with momenta
        // p1x->Fill(px[2]);
        // p1y->Fill(py[2]);
        // p1z->Fill(pz[2]);
        // p2x->Fill(px[3]);
        // p2y->Fill(py[3]);
        // p2z->Fill(pz[3]);
        // p3x->Fill(px[4]);
        // p3y->Fill(py[4]);
        // p3z->Fill(pz[4]);


        //double mPiSystem = sqrt((pow(energy[1],2)-pow(px[1],2)-pow(py[1],2)-pow(pz[1],2)));
        double Mpi3 = sqrt((pow((energy[3]+energy[4]+energy[5]),2)-pow((px[3]+px[4]+px[5]),2)-pow((py[3]+py[4]+py[5]),2)-pow((pz[3]+pz[4]+pz[5]),2)));
        //resonance->Fill(mPiSystem);
        resonance->Fill(Mpi3);



        // Define the boost direction (opposite direction of this vector)
        TVector3 boost_direction(pnu); // Example boost direction
        boost_direction = boost_direction.Unit(); // Opposite direction


        double pnuenergy = energy[0];
        double Wenergy = energy[1];

        // Compute the boost velocity vector v = p / E
        TVector3 boost = boost_direction * (pW.Mag() / Wenergy);

        // Convert initial vectors to TLorentzVector (assuming mass = 1)
        TLorentzVector v1(p1, energy[2]);
        TLorentzVector v2(p2, energy[3]);
        TLorentzVector v3(p3, energy[4]);
        TLorentzVector v4(p4, energy[1]);
        TLorentzVector v5(pW, energy[1]);
        TLorentzVector v6(pW, energy[1]);


        TLorentzVector k1(p1, energy[2]);
        TLorentzVector k2(p2, energy[3]);
        TLorentzVector k3(p3, energy[4]);
        TLorentzVector k4(p4, energy[4]);


        // Apply the boost
        v1.Boost(boost);
        v2.Boost(boost);
        v3.Boost(boost);

        v4.Boost(boost);

        cout<<"Momenta="<<v4.Px()<<" "<<v4.Py()<<" "<<v4.Pz()<<endl;


        //Jim checking
        TVector3 pi1(v1.Px(), v1.Py(), v1.Pz());
        TVector3 pi2(v2.Px(), v2.Py(), v2.Pz());

        //jim comment from here for ananya beta

        TVector3 pi3(v3.Px(), v3.Py(), v3.Pz());

        TVector3 n_perpendicular = pi2.Cross(pi3);
        TVector3 n_Lab = -boost_direction;

        double beta = n_perpendicular.Angle(n_Lab);





        TVector3 zaxis(0,0,1);
        TVector3 boost_direction2(zaxis); // Example boost direction
        boost_direction2 = boost_direction2.Unit(); // Opposite direction
        TVector3 boostLab = boost_direction2 * 0.86602540378;



        v5.Boost(boostLab);
        k1.Boost(boostLab);
        k2.Boost(boostLab);
        k3.Boost(boostLab);
        cout<<"Momenta after boost="<<v5.Px()<<" "<<v5.Py()<<" "<<v5.Pz()<<" "<<v5.Energy()<<endl;



        TVector3 WLab(v5.Px(), v5.Py(), v5.Pz());


        TVector3 oppWLab = -WLab;
        TVector3 boost_direction3 = oppWLab.Unit();;
        TVector3 boostWrest = boost_direction3 * (WLab.Mag() / v5.Energy());

        cout<<"WLab.Mag="<<WLab.Mag()/ v5.Energy()<<endl;


        k1.Boost(boostWrest);
        k2.Boost(boostWrest);
        k3.Boost(boostWrest);
        v5.Boost(boostWrest);

        cout<<"Momenta after boost to rest frame="<<v5.Px()<<" "<<v5.Py()<<" "<<v5.Pz()<<" "<<v5.Energy()<<endl;
        cout<<"Momenta after boost to rest frame="<<k1.Px()<<" "<<k1.Py()<<" "<<k1.Pz()<<" "<<k1.Energy()<<endl;
        cout<<"Momenta after boost to rest frame="<<k2.Px()<<" "<<k2.Py()<<" "<<k2.Pz()<<" "<<k2.Energy()<<endl;
        cout<<"Momenta after boost to rest frame="<<k3.Px()<<" "<<k3.Py()<<" "<<k3.Pz()<<" "<<k3.Energy()<<endl;


        TVector3 pi1rest(k1.Px(), k1.Py(), k1.Pz());
        TVector3 pi2rest(k2.Px(), k2.Py(), k2.Pz());
        TVector3 pi3rest(k3.Px(), k3.Py(), k3.Pz());


        // Compute cross product
        TVector3 cross_product = pi1rest.Cross(pi2rest);

        //double beta = cross_product.Angle(oppWLab);//Ananya beta

        //TVector3 cross_product2 = p2.Cross(p3);

        // Compute angle beta with the z-axis

        //double Mag1= cross_product.Mag();
        //double Mag2= p3.Mag();

        //double value = cross_product.Dot(p3) / (Mag1 * Mag2);

        //double beta = value;
        //double beta = cross_product2.Angle(cross_product);
        //double beta = cross_product2.Angle(TVector3(0, 0, -1));
        //double beta = p1.Angle(-p2);

        //cout<<beta<<":"<<cos(beta)<<endl;
        // if(cos(beta) > 0.5) {
        //   h1->Fill(cos(beta));
        // } else {
        //   h1->Fill(cos(beta));
        // }
        h1->Fill(cos(beta));

        //double scalar_triple_product = (p1.Cross(p2)).Dot(p3);
        //std::cout << "Scalar Triple Product: " << scalar_triple_product << std::endl;









        //while (dane != "") {out.write(& dane[0], dane.length()); out<<"\n"; getline(in, dane);}
        out<<"\n";
        event_count++;

        // for (int i=1; i<6; i++)
        // {
        //  getline(in, dane);
        //  if (dane == " ***************************************************************************")
        //     {
        //      getline(in, dane);
        //      out.write(& dane[0], dane.length()); out<<"\n"; getline(in, dane);
        //      out.write(& dane[0], dane.length()); out<<"\n"; getline(in, dane);
        //      out<<"\n";
        //     }
        // }
      }
//     if (dane == " *NCHAN    NOEVTS  PART.WIDTH     ERROR       ROUTINE    DECAY MODE        *")
//      {
//       for(int i=1; i<(event_count + 1) ;i++)
//       {
//        getline(in, dane);
//         if(i==event_count)
//           {
//            out.write(& dane[0], dane.length());
// //           chan.write(& dane[0], dane.length());
// //           chan<<"\n";
//            out<<"\n \n";
//            out<<"________________________________________next channel____________________________________________________________";
// 	   out<<"\n";
// 	  }
//       }
//      }
  }
  in.close();
  out.close();


  // TVector3 v1(5, 6, 7);
  // TVector3 v2(1, 2, 3);
  // TVector3 v3(-6, -8, -10);


  // TVector3 cross_product = v1.Cross(v2);
  // cout<<"cross product: "<<cross_product.X()<<" "<<cross_product.Y()<<" "<<cross_product.Z()<<endl;
  // cout<<"scalr product: "<<cross_product.Dot(v3)<<endl;


  TCanvas *c1 = new TCanvas("c1", "c1", 1000, 600);
  c1->Divide(3, 4);
  c1->cd(1);
  gPad->SetLogy();
  p1x->Draw();
  c1->cd(2);
  gPad->SetLogy();
  p1y->Draw();
  c1->cd(3);
  gPad->SetLogy();
  p1z->Draw();
  c1->cd(4);
  gPad->SetLogy();
  p2x->Draw();
  c1->cd(5);
  gPad->SetLogy();
  p2y->Draw();
  c1->cd(6);
  gPad->SetLogy();
  p2z->Draw();
  c1->cd(7);
  gPad->SetLogy();
  p3x->Draw();
  c1->cd(8);
  gPad->SetLogy();
  p3y->Draw();
  c1->cd(9);
  gPad->SetLogy();
  p3z->Draw();
  c1->cd(10);
  gPad->SetLogy();
  p4x->Draw();
  c1->cd(11);
  gPad->SetLogy();
  p4y->Draw();
  c1->cd(12);
  gPad->SetLogy();
  p4z->Draw();
  c1->Update();
  c1->SaveAs("momentum_distribution.png");
  TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
  h1->Draw();
  c2->SaveAs("beta_distribution.png");


  TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
  resonance->Draw();

  TLatex *latex = new TLatex();
  latex->SetTextSize(0.05);
  latex->SetTextColor(kRed);
  latex->DrawLatex(0.5, 200, Form("PREZ=0.0"));
  c3->SaveAs("resonance_distribution.png");



  //return 0;
 }
