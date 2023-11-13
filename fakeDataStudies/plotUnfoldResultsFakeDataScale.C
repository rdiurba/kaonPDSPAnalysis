
#include "TH1.h"
#include "TGraph.h"
#include "TH2.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TVectorD.h"
#include "TTimeStamp.h"
#include <fstream>
#include "TMinuit.h"
#include "TString.h"
#include <vector>
#include <string.h>
#include "TLatex.h"
#include "TPaveStats.h"
#include "TDatime.h"
#include "TColor.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TPrincipal.h"
#include "TDecompChol.h"
#include <iomanip>
double pearsonBasic_chi2(std::vector<double> data, std::vector<double> mc, std::vector<double> err){
double chi2=0;
for(int j=0; j<data.size(); j++){
//double denom=err.at(j)*err.at(j);
double denom=mc.at(j);
double num=(data.at(j)-mc.at(j))*(data.at(j)-mc.at(j));
if (denom>0) chi2=chi2+num/denom;
}
return chi2;
}



double neymanBasic_chi2(std::vector<double> data, std::vector<double> mc, std::vector<double> err){
double chi2=0;
for(int j=0; j<data.size(); j++){
//double denom=err.at(j)*err.at(j);
double denom=data.at(j);
double num=(data.at(j)-mc.at(j))*(data.at(j)-mc.at(j));
if (denom>0) chi2=chi2+num/denom;
}
return chi2;
}


double pearson_chi2(std::vector<double> data, std::vector<double> mc, std::vector<double> err){
double chi2=0;
for(int j=0; j<data.size(); j++){
double denom=err.at(j)*err.at(j);
//double denom=mc.at(j);
double num=(data.at(j)-mc.at(j))*(data.at(j)-mc.at(j));
if (denom>0) chi2=chi2+num/denom;
}
return chi2;
}

double neyman_chi2(std::vector<double> data, std::vector<double> mc, std::vector<double> err){
double chi2=0;
for(int j=0; j<data.size(); j++){
//double denom=err.at(j)*err.at(j);
double denom=mc.at(j);

double num=(data.at(j)-mc.at(j))*(data.at(j)-mc.at(j));
if (denom>0) chi2=chi2+num/denom;
}
return chi2;
}




double loglike_chi2(std::vector<double> data, std::vector<double> mc, std::vector<double> err_data, std::vector<double> err_MC){
double chi2=0;
for(int j=0; j<data.size(); j++){
double d=data.at(j);
double e=mc.at(j);
double error=TMath::Power((err_MC.at(j)*err_MC.at(j)+err_data.at(j)+err_data.at(j)),0.5);
chi2=chi2+((d-e)*(d-e))/(error*error);
}
double norm_fac=23/2;
double lngamma=-norm_fac*TMath::Log(2*3.14)-norm_fac*TMath::Log(0.1)-0.5*chi2;


return -lngamma;

}




double poisson_chi2(std::vector<double> data, std::vector<double> mc, std::vector<double> err){

double lngamma=0;
for(int j=0; j<data.size(); j++){
double d=data.at(j);
double e=mc.at(j);




if (d>0 && e>0){
lngamma=lngamma+(e-d)+d*TMath::Log(d/e);
}
else lngamma=lngamma+(e-d);
}
return 2.f*lngamma;
}
void plotUnfoldResultsFakeData(double coeff)
{
//std::cout<<std::setprecision(2);
gROOT->LoadMacro("protoDUNEStyle.C");
gROOT->SetStyle("protoDUNEStyle");
gROOT->ForceStyle();
gStyle->SetTitleX(0.35);
gStyle->SetOptFit(111);
//gStyle->SetPadRightMargin(0.15);
//gStyle->SetPadLeftMargin(0.15);
TCanvas c1=TCanvas();
int index=coeff*10-1;
//TFile f("/dune/data/users/rdiurba/kaonana_run5770.root");
TFile f("kaonUnfoldStandardGenRecoPlots.root");
//TFile fThrow("kaonSystShiftOutput_1kThrows.root");
TFile f14("kaonUnfoldDatadrivenG4RW_13.root");
TFile f12("kaonUnfoldDatadrivenG4RW_11.root");
TFile f8("kaonUnfoldDatadrivenG4RW_7.root");
TFile f6(Form("kaonUnfoldStandardGenSampleScale_%1.1f.root",coeff));
//TFile fStat("kaonUnfoldWithStatShifts.root");
//TFile fStat("kaonUnfoldWithStatShifts.root");
TFile fMCStat("kaonUnfoldWithStatShiftsStandardGenTraining.root");
//TFile fResponse("kaonStatResponseShiftOutput_1kThrows_fix_v2.root");
TFile fMC("kaonUnfoldStandardGenTrainingRecoPlots.root");
TFile fMC2("kaonUnfoldStandardGenTrainingRecoPlots_FlipIndep.root");
//TFile fMC("kaonUnfoldStandardGenTrainingRecoPlots.root");
TFile fTruth("/dune/data/users/rdiurba/rootDump/kaon_cross_section_out.root");
TGraph* g4Truth=(TGraph*)fTruth.Get("inel_KE");
TGraph* g46=(TGraph*)fTruth.Get("inel_KE");
TGraph* g48=(TGraph*)fTruth.Get("inel_KE");
TGraph* g412=(TGraph*)fTruth.Get("inel_KE");
TGraph* g414=(TGraph*)fTruth.Get("inel_KE");
g46->SetLineColor(kBlue); g4Truth->SetLineColor(kRed);
g412->SetLineColor(kBlue); g414->SetLineColor(kCyan);
TH1D* recoXSec=(TH1D*)fMC.Get("crossSectionUnfoldData");
TH1D* recoXSec6=(TH1D*)f6.Get("crossSectionUnfoldTest");
TH1D* recoXSec8=(TH1D*)f8.Get("crossSectionUnfoldTest");
TH1D* recoXSec12=(TH1D*)f12.Get("crossSectionUnfoldTest");
TH1D* recoXSec14=(TH1D*)f14.Get("crossSectionUnfoldTest");
recoXSec6->SetLineColor(kBlue); recoXSec->SetLineColor(kRed);
recoXSec12->SetLineColor(kBlue); recoXSec14->SetLineColor(kCyan);
recoXSec6->SetMarkerColor(kBlue); recoXSec->SetMarkerColor(kRed);
recoXSec12->SetMarkerColor(kBlue); recoXSec14->SetMarkerColor(kCyan);
std::cout<<"RW G4"<<std::endl;
for (int i=0;i<g4Truth->GetN();i++){ 
g46->GetY()[i] *= coeff ;
//std::cout<<g46->GetY()[i]<<std::endl;
}
std::cout<<"Plotting"<<std::endl;
recoXSec->GetXaxis()->SetTitle("Unfolded Kinetic Energy [MeV]");
recoXSec->SetTitle("Fake Data Cross Section");
recoXSec->GetYaxis()->SetTitle("#sigma [mbarn]");
recoXSec->GetYaxis()->CenterTitle();
recoXSec->GetXaxis()->CenterTitle();
recoXSec->GetYaxis()->SetRangeUser(100,1000);
g46->SetLineStyle(3);
recoXSec6->SetLineStyle(3);
recoXSec->Draw("E0 P");
recoXSec6->SetMarkerStyle(kFullSquare);
recoXSec6->Draw("E0 P SAME");
/*
recoXSec8->Draw("E0 P SAME");
recoXSec12->Draw("E0 P SAME");
recoXSec14->Draw("E0 P SAME");
*/
g4Truth->Draw("SAME");
g46->Draw("SAME L");
//g48->Draw("SAME");
//g412->Draw("SAME");
//g414->Draw("SAME");
  TLegend *lErr = new TLegend(0.3,0.7,0.7,0.9);
  lErr->AddEntry(recoXSec6,Form("Fake Data: Weight Int. Points by %1.2f",coeff),"lp");
//  lErr->AddEntry(g48,"Fake Data: #sigma_{G4}*0.8","lp");
  lErr->AddEntry(recoXSec,"Fake Data: #sigma_{G4}*1","lp");
//  lErr->AddEntry(g412,"Fake Data: #sigma_{G4}*1.2","lp");
//  lErr->AddEntry(g414,"Fake Data: #sigma_{G4}*1.4","lp");
lErr->Draw("SAME");
c1.Print(Form("xSecUnfoldedFakeDataSampleScaleStandardGen_%1.2f.png",coeff));






TH1D* recoXSecFlip=(TH1D*)fMC2.Get("crossSectionUnfoldData");
recoXSecFlip->SetLineColor(kRed);
recoXSecFlip->SetMarkerColor(kRed);
recoXSecFlip->GetXaxis()->SetTitle("Unfolded Kinetic Energy [MeV]");
recoXSecFlip->SetTitle("6 GeV/c Sample: K^{+}-Ar Tot. Inel.");
recoXSecFlip->GetYaxis()->SetTitle("#sigma [mbarn]");



recoXSecFlip->GetXaxis()->SetTitle("Unfolded Kinetic Energy [MeV]");
recoXSecFlip->SetTitle("6 GeV/c Sample: K^{+}-Ar Tot. Inel.");
recoXSecFlip->GetYaxis()->SetTitle("#sigma [mbarn]");
//TH1D* copyRecoXSec=(TH1D*)recoXSec->Clone("copy_XSec");
recoXSecFlip->GetYaxis()->SetRangeUser(300,600);
//XSec->GetXaxis()->SetRangeUser(4500,6000);
recoXSec->GetXaxis()->CenterTitle();
recoXSec->GetYaxis()->CenterTitle();
g4Truth->SetLineColor(kGreen);
recoXSecFlip->Draw("E0 P");
TH1D* copyRecoXSecFlip=(TH1D*)recoXSecFlip->Clone("copyRecoXSecFlip");
//trueXSec->Draw("E0 P SAME");
g4Truth->Draw("SAME");
//copyRecoXSecFlip->Draw("E0 P SAME");
//dataXSec->Draw("E0 P SAME");
//lMC->Draw("SAME");
TLatex tL;
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c1.Print("finalPlots/xSecKaonMCAll6GeVStandardGenFlip.png");







}

