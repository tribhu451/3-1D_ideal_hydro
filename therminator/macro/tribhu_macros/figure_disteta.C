/********************************************************************************
 *                                                                              *
 *             THERMINATOR 2: THERMal heavy-IoN generATOR 2                     *
 *                                                                              *
 * Version:                                                                     *
 *      Release, 2.0.3, 1 February 2011                                         *
 *                                                                              *
 * Authors:                                                                     *
 *      Mikolaj Chojnacki   (Mikolaj.Chojnacki@ifj.edu.pl)                      *
 *      Adam Kisiel         (kisiel@if.pw.edu.pl)                               *
 *      Wojciech Broniowski (Wojciech.Broniowski@ifj.edu.pl)                    *
 *      Wojciech Florkowski (Wojciech.Florkowski@ifj.edu.pl)                    *
 *                                                                              *
 * Project homepage:                                                            *
 *      http://therminator2.ifj.edu.pl/                                         *
 *                                                                              *
 * For the detailed description of the program and further references           *
 * to the description of the model please refer to                              *
 * http://arxiv.org/abs/1102.0273                                               *
 *                                                                              *
 * This code can be freely used and redistributed. However if you decide to     *
 * make modifications to the code, please, inform the authors.                  *
 * Any publication of results obtained using this code must include the         *
 * reference to arXiv:1102.0273 and the published version of it, when           *
 * available.                                                                   *
 *                                                                              *
 ********************************************************************************/

#include <TH1D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TStyle.h>
#include "../events2chain.C"
#include "../model2legend.C"
#include "../hist2xml.C"

#define  _FIGURE_NAME_ "fig_disteta"
#define  _N_HISTOGRAMS_ 5


void Base_hist(TH1D* h0)
{

        h0->GetXaxis()->SetRangeUser(-6.,6.);
	h0->GetXaxis()->SetLabelSize(0.045);
	h0->GetXaxis()->SetTitle("#eta");
	h0->GetXaxis()->SetTitleSize(0.04);
	h0->GetXaxis()->CenterTitle();
	h0->GetXaxis()->SetTitleOffset(1.2);
	h0->GetYaxis()->SetRangeUser(0,10000);
	h0->GetYaxis()->SetLabelSize(0.045);
	h0->GetYaxis()->SetTitle("#frac{dN_{ch}}{d#eta}");
	h0->GetYaxis()->CenterTitle();
	h0->GetYaxis()->SetTitleSize(0.04);
	h0->GetYaxis()->SetNdivisions(510);


}

void Hist_prpts(TH1D* h, int aMStyle,Color_t aMColor,const double aMSize, Color_t aLColor)
{
	
	h->SetMarkerStyle(aMStyle);
	h->SetMarkerColor(aMColor);
	h->SetLineColor(aLColor);
	h->SetMarkerSize(aMSize);


}

void graph_prpts(TGraph* g, Color_t aMColor, int aMs)
{
  g->SetMarkerColor(aMColor);
  g->SetMarkerStyle(aMs);
     	
}


void figure_disteta(TString aEventDir = "../events/", Int_t aEventFiles = 1)
{

     gStyle->SetOptStat(0);
     gStyle->SetErrorX(0);
     
// ##########################################################################
// # READ ROOT FILES
// ##########################################################################
  static ParticleCoor Particle;
  Int_t   Events;
  TChain* Chain = events2chain(aEventDir, aEventFiles, &Particle, &Events);


    // TCanvas
    TCanvas* c0 = new TCanvas("c0","",800,600);
    c0->cd();
    c0->SetMargin(0.14,0.002,0.11,0.002); 


  Int_t   XBins  = 60;
  Float_t XMin   = -6.0;
  Float_t XMax   =  6.0;
  Float_t dX     = (XMax - XMin) / XBins;

  // Base Histogram
  TH1D* h0 = new TH1D("h0", "",XBins,XMin,XMax);
     Base_hist(h0);

  // HISTOGRAM
  TH1D* h1 = new TH1D("h1", "",XBins,XMin,XMax);


  Float_t P, Eta;
  Int_t   pid;



   TGraphAsymmErrors* gr1 = new TGraphAsymmErrors("AuAu_200GeV_Centrality_0-3.txt","%lg %lg %lg %lg");



  for(Int_t i=0; i<Chain->GetEntries(); i++)
    {
      Chain->GetEntry(i);
      P   = TMath::Sqrt(Particle.px*Particle.px + Particle.py*Particle.py + Particle.pz*Particle.pz);
      if(P == Particle.pz)
	continue;
      pid = Particle.pid;
      Eta = 0.5 * TMath::Log((P+Particle.pz) / (P-Particle.pz));
      if((pid == 211) || (pid == -211) || (pid == -321) || (pid == 321) || (pid == -2212) || (pid == 2212)|| 
         (pid==-10321)||(pid==10321)||(pid==-3324)||(pid==3324)||(pid==-323)||(pid==323)||(pid==-3114)||
            (pid==3114)||(pid==-3224)||(pid==3224)||(pid==-10323)||(pid==10323)||
	     (pid==-3112)||(pid==3112)||(pid==-67719)||(pid==67719)||(pid==-3334)||(pid==3334)||(pid==13214)||
              (pid==13114)||(pid==-46653)||(pid==46653)||(pid==-3312)||(pid==3312)||
	         (pid==-3222)||(pid==3222)||(pid==-9000211)||(pid==9000211))
         {
	  h1->Fill(Eta,1.);
	}
    }
 h1->Scale(1.0 / (Events * dX));
int tMaxBin = h1->GetMaximumBin();   h0->SetMaximum((h1->GetBinContent(tMaxBin) + h1->GetBinError(tMaxBin)) * 1.15);


Hist_prpts(h1, 20,kRed,0.01, kBlack);
h1->SetFillColor(kRed-2);

  gr1->SetMarkerColor(4);
   gr1->SetMarkerStyle(21);

h0->Draw();
h1->Draw("e3same");
gr1->Draw("psame");


}


