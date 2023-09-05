#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include <fstream>
using namespace std;


void ascii3 ()

{
	
	// int n = h4->GetNbinsX();
// for (int i=1; i<=n; i++)
// {
	// printf ("%g %g \n",
		// (h4->GetBinLowEdge(i) + h4->GetBinWidth(i)/2)*1000, //*1000, //
		// h4->GetBinContent(i)//*0.01 //*1250 // 
// );
// }
	const char * HN = "2";
	TH1F::SetDefaultSumw2(kTRUE);
	TH2F::SetDefaultSumw2(kTRUE);
	TFile *f_g4 = 
	new TFile("C:/Geant4/workwithg4/TestEm10_b/Release/salice_c_100_100.root");//PbF2plexi  PbF2matrix PbF2matrix_wr  PbF2multilayer_wr PbF2matrix_wr_e
		
			TH1F *h1 = (TH1F*)f_g4->Get(HN);
		
		TFile *f_g42 = 
	 new TFile("C:/Geant4/workwithg4/TestEm10_b/Release/salice_c_100_100.root"); //PbF2plexi_rad PMMAplexi PbF2plexi_wr PbF2multilayer_rad PbF2matrix_rad PbF2matrix_rad_e
		
	TH1F *h2 = (TH1F*)f_g42->Get(HN);
	
	TFile *f_g43 = 
		new TFile("C:/Geant4/workwithg4/TestEm10_b/Release/salice_c_100_100.root"); // PbF2multilayer_wr
		
	TH1F *h3= (TH1F*)f_g43->Get(HN);	
	
		TFile *f_g44 = 
		new TFile("C:/Geant4/workwithg4/TestEm10_b/Release/salice_c_100_100.root"); // PbF2multilayer_wr
		
	TH1F *h4= (TH1F*)f_g44->Get(HN);
	
	//h8->Sumw2(kTRUE);
	//TH1F *hg2 = (TH1F*)f_g4->Get("h6");
	//hg2->GetXaxis()->SetLimits(0.07, 4);
	//Scale(hg2->GetXaxis(), 1000);
	// double g4el = h11->Integral();
	// double g4el1 = h8->Integral();

	TH1F *hg1 = (TH1F*)h1->Clone("hg1");
	TH1F *hg2 = (TH1F*)h2->Clone("hg2");
	TH1F *hg3 = (TH1F*)h3->Clone("hg3");
	TH1F *hg4 = (TH1F*)h4->Clone("hg4");

    hg1->Sumw2(kTRUE);

  int b = 5; //20   5
	 hg2->Rebin(b);
	 hg1->Rebin(b);
	  hg3->Rebin(b);
	  hg4->Rebin(b);

	  // hg1->Rebin(5);
	 
   // for (int i=0;i<44446;i++) {
      // hg2->Fill(hg2t->GetRandom());
   // }
		double Nsc = 3e4;//1e5; ��������������� �� ����� ����������� ���������� runBeamOn => ������� �� ��������� ����� ������� �� 1 �
	hg1->Scale(1. / Nsc);
	hg2->Scale(1. / Nsc);
	hg3->Scale(1. / Nsc);
	hg4->Scale(1. / Nsc);
	
	// cout << g4el << endl;
	// cout << g4el1 << endl;
	
  // h15->Scale(1. / g4el);

//norm to width of bins
	// for (int i = 1; i <= hg1->GetNbinsX(); i++) {
	// 	hg1->SetBinContent(i, float(hg1->GetBinContent(i)) / hg1->GetBinWidth(i)*(hg1->GetBinLowEdge(i) + hg1->GetBinWidth(i)/2)*0.00175);
	// 	hg1->SetBinError(i, float(hg1->GetBinError(i)) / hg1->GetBinWidth(i));
	// }
	
	
	// for (int i = 1; i <= hg2->GetNbinsX(); i++) {
	// 	hg2->SetBinContent(i, float(hg2->GetBinContent(i)) / hg2->GetBinWidth(i)*(hg2->GetBinLowEdge(i) + hg2->GetBinWidth(i)/2));
	// 	hg2->SetBinError(i, float(hg2->GetBinError(i)) / hg2->GetBinWidth(i));
	// } to make spectral dW
	for (int i = 1; i <= hg1->GetNbinsX(); i++) {
		hg1->SetBinContent(i, float(hg1->GetBinContent(i)) / hg1->GetBinWidth(i));
		hg1->SetBinError(i, float(hg1->GetBinError(i)) / hg1->GetBinWidth(i));
	}
	
	
	for (int i = 1; i <= hg2->GetNbinsX(); i++) {
		hg2->SetBinContent(i, float(hg2->GetBinContent(i)) / hg2->GetBinWidth(i));
		hg2->SetBinError(i, float(hg2->GetBinError(i)) / hg2->GetBinWidth(i));
	}
		
	for (int i = 1; i <= hg3->GetNbinsX(); i++) {
		hg3->SetBinContent(i, float(hg3->GetBinContent(i)) / hg3->GetBinWidth(i));
		hg3->SetBinError(i, float(hg3->GetBinError(i)) / hg3->GetBinWidth(i));
	}
	
	for (int i = 1; i <= hg4->GetNbinsX(); i++) {
		hg4->SetBinContent(i, float(hg4->GetBinContent(i)) / hg4->GetBinWidth(i));
		hg4->SetBinError(i, float(hg4->GetBinError(i)) / hg4->GetBinWidth(i));
	}
	
	//hg2->Scale(1. / g4el);
	
		// hg1->Add(h_w_back,+1); 
	// hg2->Add(h_tt_back,+1);
	
	
	// hg2->GetXaxis()->SetRangeUser(5500., 8500);
	// hg2->SetStats(0);
   // hg2->SetTitle("Number of Cherenkov photons in SiPM #23  ");//// Number of Cherenkov photons in all SiPMs
   // hg2->GetXaxis()->SetTitle("Number of photons");//#Theta [rad]   XTR energy [keV]
   // hg2->GetYaxis()->SetTitle("Number of Events");//
	
   
   
   hg1->GetXaxis()->SetRangeUser(0, 33);//(0.09, 0.511); (0.41, 0.605)//diapazon po x
   //hg1->GetYaxis()->SetRangeUser(0., 1.);//diapazon po y
   hg1->SetStats(0);//delete stats if 0
   hg1->SetTitle(0);//no title if 0
   //hg1->SetTitle("Cherenkov photons registred by the forth layer of SiPMs");//  sum of Cherenkov photon energies in all SiPMs //Cherenkov photons in all SiPMs
	//axis label
	 hg1->GetXaxis()->SetTitle("#hbar#omega [MeV]");//#lambda [um] global time ns  Electron energy [MeV] // #lambda [um]// #Theta [rad]   XTR energy [keV]  
   hg1->GetYaxis()->SetTitle("#frac{dN}{d#hbar#omega} x10^{-2}") ;//#frac{dN}{d#lambda} #frac{N_{ph}}{N_{gamma}} #frac{N_ph}{N_gamma}Counts #lambda #frac{dN}{d#hbar#omega} [#frac{photon}{eV}]// probab. to deposit energy in event [#frac{1}{MeV}]
	//colors etc
    hg1->SetLineColor(kBlack);
    hg1->SetMarkerStyle(0);
	hg1->SetMarkerSize(0.8);
	hg1->SetLineWidth(2);
	hg1->SetMarkerColor(kBlack);	

	
	hg2->SetLineColor(kRed);
    hg2->SetMarkerStyle(0);
	hg2->SetMarkerSize(0.8);
	hg2->SetLineWidth(2);
	hg2->SetMarkerColor(kBlue);
	
    hg3->SetLineColor(kMagenta);
    hg3->SetMarkerStyle(0);
	hg3->SetLineStyle(0);
	hg3->SetMarkerSize(0.8);
	hg3->SetLineWidth(2);
	hg3->SetMarkerColor(kMagenta);
	
	hg4->SetLineColor(kMagenta);
    hg4->SetMarkerStyle(0);
	hg4->SetMarkerSize(0.8);
	hg4->SetLineWidth(2);
	hg4->SetMarkerColor(kMagenta);
	//size of labels and numbers
	hg1->GetXaxis()->SetTitleSize(0.08);
	hg1->GetYaxis()->SetTitleSize(0.08);
	hg1->GetXaxis()->SetLabelSize(0.08);
	hg1->GetYaxis()->SetLabelSize(0.08);
	hg1->GetXaxis()->SetTitleOffset(1.1);//perecrytie of numbers with axis label
	hg1->GetXaxis()->SetLabelOffset(0.02);//+-

	hg1->GetYaxis()->SetTitleOffset(0.75);//+-
	hg1->GetYaxis()->SetTitleFont(42);
	hg1->GetXaxis()->SetTitleFont(42);
	hg1->GetYaxis()->SetLabelFont(42);
	hg1->GetXaxis()->SetLabelFont(42);

	//drawing
   TCanvas *c_cf6 = new TCanvas("c6", "c6", 1200, 800);//size of  hist
   //borders 
   c_cf6->SetLeftMargin(0.23);
   c_cf6->SetRightMargin(0.08);
   c_cf6->SetTopMargin(0.05);
c_cf6->SetBottomMargin(0.2);
//drawing of the hist
   hg1->Draw("E1 hist");
   hg2->Draw("E1 hist Same");//там же, где и hg1 нарисует
//    hg3->Draw("E1 hist Same");
   //hg4->Draw("E1 hist Same");
     	

   //отрисовка легенды
   TLegend *leg5 = new TLegend(0.6, 0.55, 0.82, 0.65, NULL, "brNDC");// координаты или размеры легенды
   TLegendEntry  * entry5;
   entry5 = leg5->AddEntry(hg1, "Geant4 800 MeV Ge(1,1,1)", "pl");//тип отображения легенды во вторых кавычках, p-точка,l-линия
//    entry5 = leg5->AddEntry(hg2, "e^{+} 2.6 GeV", "pl"); //Setup 3 PMMA
//    entry5 = leg5->AddEntry(hg3, "e^{+} 3.0 GeV", "pl");
  // entry5 = leg5->AddEntry(hg4, "e^{-}", "pl");
   leg5->SetTextSize(0.07);
   leg5->SetTextFont(42);
   leg5->SetFillColor(0);
   leg5->SetFillStyle(0);
   leg5->SetBorderSize(0);
   leg5->Draw();
   
	int n = hg1->GetNbinsX();
	for (int i=1; i<=n; i++)
	{
		printf ("%g %g \n",
			(hg1->GetBinLowEdge(i) + hg1->GetBinWidth(i)/2), //*1000, //
			hg1->GetBinContent(i)//*0.01 //*1250 // 
	);
	}

   // TLatex* ltx1 = new TLatex(); 
	// ltx1->SetTextSize(0.08);
    // ltx1->SetTextAngle(0.);
    // ltx1->DrawLatex(500.0,.01,Form("your text"));
	
	
	
	// TCanvas *c_cf5 = new TCanvas("c5", "c5", 1200, 800);
   // hg2->Draw("E1 hist");
   
     	

   
   // TLegend *leg5 = new TLegend(0.6, 0.55, 0.82, 0.65, NULL, "brNDC");
   // TLegendEntry  * entry5;
   // entry5 = leg5->AddEntry(hg1, "Geant4", "pl");
   // leg5->SetTextSize(0.04);
   // leg5->SetTextFont(42);
   // leg5->SetFillColor(0);
   // leg5->SetFillStyle(0);
   // leg5->SetBorderSize(0);
   // leg5->Draw();
   
   // TLatex* ltx1 = new TLatex(); 
	// ltx1->SetTextSize(0.08);
    // ltx1->SetTextAngle(0.);
    // ltx1->DrawLatex(5500.0,5.04,Form("cut value 1mm"));
	
	// hg2->Print("base");
	// hg2->Print("range");
	// hg2->Print("all");

}

//����������� ��� ��� ������� �� ������� ���� ��������(2):
//.x C:/Users/sanch/Documents/Geant4/workwithg4/PbF2calorim/CalorimPbF2_b/RelWithDebInfo/ascii2.c
//.x C:/Users/sanch/Documents/workwithg4/BCompton_b/RelWithDebInfo/asciiBC.c; > C:\Users\sanch\Documents\DAT\20gevespectEl_250bins_1304.dat
//������ ����� ������ � ���� ���������� ����
//.x C:/Geant4/workwithg4/TestEm10_b/Release/ascii3.c; > C:/Geant4/workwithg4/TestEm10_b/Release/graph.txt
//salice_si_110_100
//salice_ge_111_5