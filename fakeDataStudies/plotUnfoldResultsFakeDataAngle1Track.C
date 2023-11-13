
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
TFile f16(Form("kaonUnfoldFakeDataScaleForward_4.000_maxTracks1.root"));
TFile f6(Form("kaonUnfoldFakeDataScaleForward_2.000_maxTracks1.root"));
TFile f8(Form("kaonUnfoldFakeDataScatterExtend_0.00.root"));
TFile f12(Form("kaonUnfoldFakeDataScatterExtend_2.00.root"));
TFile f14(Form("kaonUnfoldFakeDataFlatAngle_0.975_1.root"));
TFile f24(Form("kaonUnfoldFakeDataFlatAngle_0.975_100.root"));
TFile fMCStat("kaonUnfoldWithStatShiftsStandardGenTraining.root");

TFile fMC("kaonUnfoldStandardGenTrainingRecoPlots.root");
TFile fMC2("kaonUnfoldStandardGenTrainingRecoPlots_FlipIndep.root");
//TFile fMC("kaonUnfoldStandardGenTrainingRecoPlots.root");
TFile fTruth("/dune/data/users/rdiurba/rootDump/kaon_cross_section_out.root");
TGraph* g4Truth=(TGraph*)fTruth.Get("inel_KE");
TGraph* g46=(TGraph*)fTruth.Get("inel_KE");
//g46->SetLineColor(kBlue); 
g4Truth->SetLineColor(kRed);

//g412->SetLineColor(kBlue); g414->SetLineColor(kCyan);
TH1D* recoXSec=(TH1D*)fMC.Get("crossSectionUnfoldData");
TH1D* recoXSec6=(TH1D*)f6.Get("crossSectionUnfoldTest");
TH1D* recoXSec8=(TH1D*)f8.Get("crossSectionUnfoldTest");
TH1D* recoXSec12=(TH1D*)f12.Get("crossSectionUnfoldTest");
TH1D* recoXSec14=(TH1D*)f14.Get("crossSectionUnfoldTest");
TH1D* recoXSec16=(TH1D*)f16.Get("crossSectionUnfoldTest");
TH1D* recoXSec24=(TH1D*)f24.Get("crossSectionUnfoldTest");
TH1D* histLeadingAngle=(TH1D*)f14.Get("histLeadingAngle");
TH1D* weightLeadingAngle=(TH1D*)f14.Get("weightLeadingAngle");


recoXSec6->SetLineColor(kBlack); recoXSec->SetLineColor(kRed); recoXSec16->SetLineColor(kMagenta);
recoXSec12->SetLineColor(kBlue); recoXSec14->SetLineColor(kMagenta); recoXSec8->SetLineColor(kGreen);
recoXSec6->SetMarkerColor(kBlack); recoXSec->SetMarkerColor(kRed); recoXSec8->SetMarkerColor(kGreen);
recoXSec12->SetMarkerColor(kBlue); recoXSec14->SetMarkerColor(kMagenta); recoXSec16->SetMarkerColor(kMagenta);
recoXSec24->SetMarkerColor(kBlack); 
std::cout<<"RW G4"<<std::endl;
/*
for (int i=0;i<g4Truth->GetN();i++){ 
g46->GetY()[i] *= coeff ;
//std::cout<<g46->GetY()[i]<<std::endl;
}*/
std::cout<<"Plotting"<<std::endl;
recoXSec->GetXaxis()->SetTitle("Unfolded Kinetic Energy [MeV]");
std::string chi2=recoXSec->GetTitle();
recoXSec->SetTitle("Fake Data with Endpoint/Angle Reweights");
recoXSec->GetYaxis()->SetTitle("#sigma [mbarn]");
recoXSec->GetYaxis()->CenterTitle();
recoXSec->GetXaxis()->CenterTitle();
recoXSec->GetYaxis()->SetRangeUser(300,600);
g46->SetLineStyle(3);
recoXSec6->SetLineStyle(3);
recoXSec24->SetLineStyle(3);
recoXSec8->SetLineStyle(4);
recoXSec12->SetLineStyle(5);
recoXSec14->SetLineStyle(6);
recoXSec16->SetLineStyle(6);
recoXSec->Draw("HIST");
recoXSec6->SetMarkerStyle(kFullSquare);
recoXSec6->Draw("HIST SAME");
recoXSec8->SetMarkerStyle(kStar);
recoXSec8->Draw("HIST  SAME");
recoXSec12->SetMarkerStyle(kCross);
recoXSec12->Draw("HIST SAME");
recoXSec16->SetMarkerStyle(kPlus);
recoXSec16->Draw("HIST  SAME");
recoXSec14->SetMarkerStyle(kPlus);
//recoXSec14->Draw("HIST  SAME");
recoXSec24->SetMarkerStyle(kFullSquare);
//recoXSec24->Draw("HIST SAME");
/*
recoXSec8->Draw("E0 P SAME");
recoXSec12->Draw("E0 P SAME");
recoXSec14->Draw("E0 P SAME");
*/
//g4Truth->Draw("SAME");
//g46->Draw("SAME L");
//g48->Draw("SAME");
//g412->Draw("SAME");
//g414->Draw("SAME");
  TLegend *lErr = new TLegend(0.15,0.67,0.85,0.87);
  lErr->AddEntry(recoXSec,Form("Nominal Fake Data: #chi^{2}=%s",chi2.c_str()),"l");
  lErr->AddEntry(recoXSec6,Form("Sec. Kaon with cos(#theta)>0.995 for all 1 K+0 tracks evts. reweight by 2: #chi^{2}=%s",recoXSec6->GetTitle()),"l");
  lErr->AddEntry(recoXSec16,Form("Sec. Kaon with cos(#theta)>0.995 for all 1 K+0 tracks evts. reweight by 4: #chi^{2}=%s",recoXSec16->GetTitle()),"l");
  //lErr->AddEntry(recoXSec14,Form("Flat Dist. for Evt. with One Kaon 0.975<cos(#theta)<1: #chi^{2}=%s",recoXSec14->GetTitle()),"l");
  //lErr->AddEntry(recoXSec24,Form("Flat Dist. for Evt. with Leading Kaon 0.975<cos(#theta)<1: #chi^{2}=%s",recoXSec24->GetTitle()),"l");
  lErr->AddEntry(recoXSec8,Form("Reweight of Extended Trks by 0: #chi^{2}=%s",recoXSec8->GetTitle()),"l");
  lErr->AddEntry(recoXSec12,Form("Reweight of Extended Trks by 2: #chi^{2}=%s",recoXSec12->GetTitle()),"l");

  TLegend *lErr2 = new TLegend(0.2,0.7,0.8,0.9);
  lErr2->AddEntry(recoXSec,Form("Nominal Fake Data: #chi^{2}=%s",chi2.c_str()),"l");
  //lErr2->AddEntry(recoXSec6,Form("Sec. Kaon with cos(#theta)>0.995 reweight times 2:  #chi^{2}=%s",recoXSec6->GetTitle()),"l");
  //lErr2->AddEntry(recoXSec16,Form("Sec. Kaon with cos(#theta)>0.995 reweight times 0:  #chi^{2}=%s",recoXSec16->GetTitle()),"l");
  lErr2->AddEntry(recoXSec14,Form("Flat Dist. for Evt. with One Kaon 0.975<cos(#theta)<1: #chi^{2}=%s",recoXSec14->GetTitle()),"l");
  lErr2->AddEntry(recoXSec24,Form("Flat Dist. for Evt. with Leading Kaon 0.975<cos(#theta)<1: #chi^{2}=%s",recoXSec24->GetTitle()),"l");
  lErr2->AddEntry(recoXSec8,Form("Reweight of Extended Trks by 0: #chi^{2}=%s",recoXSec8->GetTitle()),"l");
  lErr2->AddEntry(recoXSec12,Form("Reweight of Extended Trks by 2: #chi^{2}=%s",recoXSec12->GetTitle()),"l");


//  lErr->AddEntry(g48,"Fake Data: #sigma_{G4}*0.8","lp");

//  lErr->AddEntry(g412,"Fake Data: #sigma_{G4}*1.2","lp");
//  lErr->AddEntry(g414,"Fake Data: #sigma_{G4}*1.4","lp");
lErr->Draw("SAME");
c1.Print(Form("xSecUnfoldedFakeDataStandardGenScaleAngle1Tracks.png"));

recoXSec->Draw("HIST");

recoXSec8->Draw("HIST  SAME");
recoXSec12->Draw("HIST  SAME");
recoXSec14->Draw("HIST  SAME");
recoXSec24->Draw("HIST SAME");
lErr2->Draw("SAME");
c1.Print(Form("xSecUnfoldedFakeDataStandardGenFlatAngle.png"));


gPad->SetLeftMargin(0.15);
histLeadingAngle->SetTitle("Fake Data Leading Angle");
histLeadingAngle->GetXaxis()->SetTitle("Cos(#theta) with Parent");
histLeadingAngle->GetYaxis()->SetTitle("Number of Events");
histLeadingAngle->GetXaxis()->CenterTitle();
histLeadingAngle->GetYaxis()->CenterTitle();
histLeadingAngle->SetLineColor(kRed);
histLeadingAngle->SetMarkerColor(kRed);
histLeadingAngle->Draw("E0 P");
c1.Print(Form("fakeDataHistLeadingAngle.png"));

gPad->SetLeftMargin(0.15);
weightLeadingAngle->SetTitle("Fake Data Leading Angle");
weightLeadingAngle->GetXaxis()->SetTitle("Cos(#theta) with Parent");
weightLeadingAngle->GetYaxis()->SetTitle("Weight to Flatten Distribution");
weightLeadingAngle->GetXaxis()->CenterTitle();
weightLeadingAngle->GetYaxis()->CenterTitle();
weightLeadingAngle->SetLineColor(kRed);
weightLeadingAngle->SetMarkerColor(kRed);
weightLeadingAngle->Draw("HIST");
c1.Print(Form("fakeDataWeightLeadingAngle.png"));

















/*



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


*/




}

