	
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
void plotUnfoldResultsStandardGen()
{

std::vector<double> ke={2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0};
std::vector<double> hA={470.921, 472.224, 469.255, 471.893, 468.27, 469.748, 470.461, 469.203, 468.056, 470.117, 470.208, 471.076, 467.991, 468.944, 470.487, 470.61, 470.454, 471.044, 469.709, 470.325, 470.078, 469.3, 469.424, 472.243, 471.07, 471.459, 470.052, 469.378, 470.538, 469.378, 470.422, 471.563, 468.354, 471.239, 468.186, 470.493, 469.722, 469.722, 466.753, 467.868, 470.292, 468.011, 470.156, 467.012, 470.389, 470.519, 468.322, 469.858, 471.809, 465.943, 469.307, 467.862, 469.838, 470.007, 467.674, 471.154, 468.801, 467.933, 470.707, 468.795, 471.952};
std::vector<double> hN={468.711, 468.069, 467.064, 466.124, 469.988, 471.167, 468.069, 466.974, 469.676, 468.4, 465.619, 468.704, 470.526, 465.224, 468.445, 466.578, 464.42, 469.391, 469.832, 468.186, 464.536, 467.421, 468.225, 467.933, 465.211, 467.907, 468.588, 466.974, 469.761, 469.028, 466.423, 467.304, 466.487, 464.115, 470.072, 467.518, 469.631, 466.922, 466.753, 467.7, 468.075, 467.181, 469.877, 470.104, 466.889, 469.67, 468.938, 465.729, 467.661, 466.144, 467.084, 467.635, 469.657, 468.652, 470.214, 466.896, 468.011, 471.005, 466.124, 466.786, 466.954};
TH1D* hA_hist=new TH1D("hA_hist","hA_hist", ke.size(),1950,8050);
TH1D* hN_hist=new TH1D("hN_hist","hN_hist", ke.size(),1950,8050);
for(int i=0; i<ke.size(); i++){
hA_hist->SetBinContent(i+1,hA.at(i));
hN_hist->SetBinContent(i+1,hN.at(i));

}

//std::cout<<std::setprecision(2);
gROOT->LoadMacro("protoDUNEStyle.C");
gROOT->SetStyle("protoDUNEStyle");
gROOT->ForceStyle();
gStyle->SetTitleX(0.35);
gStyle->SetOptFit(111);
//gStyle->SetPadRightMargin(0.15);
//gStyle->SetPadLeftMargin(0.15);
TCanvas c1=TCanvas();

TFile fEff("effHistogram.root");
TH1D* effHist=(TH1D*)fEff.Get("effHist");
TH1D* effHist_close=(TH1D*)fEff.Get("effHist_close");
effHist->GetXaxis()->CenterTitle();
effHist->GetYaxis()->CenterTitle();
effHist->GetYaxis()->SetRangeUser(0.2,1.4);
effHist_close->GetXaxis()->CenterTitle();
effHist_close->GetYaxis()->CenterTitle();
effHist_close->GetXaxis()->SetTitle("True Cos(#theta) for sole K^{+} secondary and beam K^{+}");
effHist->GetXaxis()->SetTitle("True Cos(#theta) for sole K^{+} secondary and beam K^{+}");
effHist_close->GetYaxis()->SetTitle("Fraction with Vertex Correctly Selected");
effHist->GetYaxis()->SetTitle("Fraction with Vertex Correctly Selected");
effHist_close->GetYaxis()->SetRangeUser(0.2,1.4);



TLatex tL;
effHist->SetLineColor(kRed);
effHist_close->SetLineColor(kRed);

effHist->SetMarkerColor(kRed);
effHist_close->SetMarkerColor(kRed);


effHist->SetTitle("6 GeV/c Sample: Evts. with 1 K^{+/-}");

effHist_close->SetTitle("6 GeV/c Sample: Evts. with 1 K^{+/-}");

effHist->Draw("E0 P ");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c1.Print("finalPlots/efficiencyPlot6GeV.png");


effHist_close->Draw("E0 P");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c1.Print("finalPlots/efficiencyClosePlot6GeV.png");

TH1D* effHistOnePart_close=(TH1D*)fEff.Get("effHistAnyPart_close");
effHistOnePart_close->SetTitle("6 GeV/c Sample: Evts. with 1 Track");
effHistOnePart_close->GetXaxis()->CenterTitle();
effHistOnePart_close->SetMarkerColor(kRed);
effHistOnePart_close->SetLineColor(kRed);
effHistOnePart_close->GetYaxis()->CenterTitle();
effHistOnePart_close->GetXaxis()->CenterTitle();
effHistOnePart_close->GetYaxis()->CenterTitle();
effHistOnePart_close->GetXaxis()->SetTitle("True Cos(#theta) for sole secondary and beam K^{+}");
effHistOnePart_close->GetYaxis()->SetTitle("Fraction with Vertex Correctly Selected");
//effHist->GetXaxis()->SetTitle("Cos(#theta) for sole K^{+} secondary and beam K^{+}");
effHistOnePart_close->GetYaxis()->SetRangeUser(0.0,1.4);
effHistOnePart_close->Draw("E0 P");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c1.Print("finalPlots/efficiencyCloseOnePartPlot6GeV.png");

TH1D* totalHistOnePart=(TH1D*)fEff.Get("total_histAnyPart");
totalHistOnePart->SetTitle("6 GeV/c Sample: Evts. With 1 Track");
totalHistOnePart->GetXaxis()->CenterTitle();
totalHistOnePart->SetMarkerColor(kRed);
totalHistOnePart->SetLineColor(kRed);
totalHistOnePart->GetYaxis()->CenterTitle();
totalHistOnePart->GetXaxis()->CenterTitle();
totalHistOnePart->GetYaxis()->CenterTitle();
totalHistOnePart->GetXaxis()->SetTitle("True Cos(#theta) for sole secondary and beam K^{+}");
totalHistOnePart->GetYaxis()->SetTitle("Number of Events");
//effHist->GetXaxis()->SetTitle("Cos(#theta) for sole K^{+} secondary and beam K^{+}");
totalHistOnePart->GetYaxis()->SetRangeUser(0.0,totalHistOnePart->GetMaximum()*1.4);
totalHistOnePart->Draw("E0 P");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c1.Print("finalPlots/totalHistOnePartPlot6GeV.png");



TH1D* totalHist=(TH1D*)fEff.Get("total_hist");
totalHist->SetTitle("6 GeV/c Sample: Evts. with 1 K^{+/-}");
totalHist->GetXaxis()->CenterTitle();
totalHist->SetMarkerColor(kRed);
totalHist->SetLineColor(kRed);
totalHist->GetYaxis()->CenterTitle();
totalHist->GetXaxis()->CenterTitle();
totalHist->GetYaxis()->CenterTitle();
totalHist->GetXaxis()->SetTitle("True Cos(#theta) for sole secondary and beam K^{+}");
totalHist->GetYaxis()->SetTitle("Number of Events");
//effHist->GetXaxis()->SetTitle("Cos(#theta) for sole K^{+} secondary and beam K^{+}");
totalHist->GetYaxis()->SetRangeUser(0.0,totalHist->GetMaximum()*1.4);
totalHist->Draw("E0 P");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c1.Print("finalPlots/totalHistPlot6GeV.png");


TH1D* dEOnePart_close=(TH1D*)fEff.Get("keDeltaEAnyPart_close");
dEOnePart_close->SetTitle("6 GeV/c Sample: Evts. With 1 Track");
dEOnePart_close->GetXaxis()->CenterTitle();
dEOnePart_close->SetMarkerColor(kRed);
dEOnePart_close->SetLineColor(kRed);
dEOnePart_close->GetYaxis()->CenterTitle();
dEOnePart_close->GetXaxis()->CenterTitle();
dEOnePart_close->GetYaxis()->CenterTitle();
dEOnePart_close->GetXaxis()->SetTitle("#DeltaKE between Secondary Kaon and Parent with true cos(#theta)>0.995 [MeV]");
dEOnePart_close->GetYaxis()->SetTitle("Number of Events");
//effHist->GetXaxis()->SetTitle("Cos(#theta) for sole K^{+} secondary and beam K^{+}");
dEOnePart_close->GetYaxis()->SetRangeUser(0.0,dEOnePart_close->GetMaximum()*1.4);
dEOnePart_close->Draw("E0 P");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c1.Print("finalPlots/dEOnePartClosePlot6GeV.png");





TH1D* effHistOnePart=(TH1D*)fEff.Get("effHistAnyPart");
effHistOnePart->GetYaxis()->SetTitle("Fraction with Vertex Correctly Selected");
effHistOnePart->SetTitle("6 GeV/c Sample: Evts. with 1 Track");
effHistOnePart->GetXaxis()->CenterTitle();
effHistOnePart->SetMarkerColor(kRed);
effHistOnePart->SetLineColor(kRed);
effHistOnePart->GetYaxis()->CenterTitle();
effHistOnePart->GetXaxis()->CenterTitle();
effHistOnePart->GetYaxis()->CenterTitle();
effHistOnePart->GetXaxis()->SetTitle("True Cos(#theta) for sole secondary and beam K^{+}");
//effHist->GetXaxis()->SetTitle("Cos(#theta) for sole K^{+} secondary and beam K^{+}");
effHistOnePart->GetYaxis()->SetRangeUser(0.0,1.4);
effHistOnePart->Draw("E0 P");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c1.Print("finalPlots/efficiencyOnePartPlot6GeV.png");







//TFile f("/dune/data/users/rdiurba/kaonana_run5770.root");
TFile f("kaonUnfoldStandardGenRecoPlots.root");
//TFile fThrow("kaonSystShiftOutput_1kThrows.root");
TFile fThrow("kaonSystShiftOutputStandardGeni_all_1000.root");
//TFile fStat("kaonUnfoldWithStatShifts.root");
//TFile fStat("kaonUnfoldWithStatShifts.root");
//TFile fMCStat("kaonUnfoldWithStatShiftsStandardGenTraining.root");
TFile fResponseMCMC("kaonTrainingRespShiftOutputStandardGen_all_1000.root");
TFile fResponseMC("kaonRespShiftOutputStandardGen_all_1000.root");
TFile fMC("kaonUnfoldStandardGenTrainingRecoPlots.root");
//TFile fMC("kaonUnfoldStandardGenTrainingRecoPlots.root");
TFile fTruth("/dune/data/users/rdiurba/rootDump/kaon_cross_section_out.root");
TGraph* g4Truth=(TGraph*)fTruth.Get("inel_KE");
TH1D* fakes=(TH1D*)f.Get("fakes");
TH1D* fakesIncident=(TH1D*)f.Get("fakesIncident");
TH1D* missInc=(TH1D*)f.Get("h1d_miss_incidentE");
TH1D* missInt=(TH1D*)f.Get("h1d_miss_interactingE");
TH2D* covMatSyst=(TH2D*)fThrow.Get("covMat");
TH2D* covMatMCRespStat=(TH2D*)fResponseMC.Get("covMat");
TH2D* covMatMCRespStatMC=(TH2D*)fResponseMCMC.Get("covMat");
TH2D* covMatStat=(TH2D*)f.Get("covMatStat");
TH2D* covMatMCStat=(TH2D*)fMC.Get("covMatStat");
//TH2D* covMatStatResponse=(TH2D*)fResponse.Get("covMatResponse");
//TH2D* covMatMCStatResponse=(TH2D*)fResponse.Get("covMatMCResponse");
TH1D* h1d_DataIntSyst=(TH1D*)fThrow.Get("h1d_DataInt");
TH1D* h1d_DataIncSyst=(TH1D*)fThrow.Get("h1d_DataInc");


TH1D* h1d_DataIntResp=(TH1D*)fResponseMC.Get("h1d_DataInt");
TH1D* h1d_DataIncResp=(TH1D*)fResponseMC.Get("h1d_DataInc");


TH1D* h1d_RecoIntResp=(TH1D*)fResponseMCMC.Get("h1d_DataInt");
TH1D* h1d_RecoIncResp=(TH1D*)fResponseMCMC.Get("h1d_DataInc");




TH1D* h1d_DataIntStat=(TH1D*)f.Get("h1d_DataIntStat");
TH1D* h1d_DataIncStat=(TH1D*)f.Get("h1d_DataIncStat");
TH1D* h1d_MCIntStat=(TH1D*)fMC.Get("h1d_DataIntStat");
TH1D* h1d_MCIncStat=(TH1D*)fMC.Get("h1d_DataIncStat");

//TH1D* h1d_DataIntStatResponse=(TH1D*)fResponse.Get("h1d_DataIntStatResponse");
//TH1D* h1d_DataIncStatResponse=(TH1D*)fResponse.Get("h1d_DataIncStatResponse");


//TH1D* h1d_MCIntStatResponse=(TH1D*)fResponse.Get("h1d_MCIntStatResponse");
//TH1D* h1d_MCIncStatResponse=(TH1D*)fResponse.Get("h1d_MCIncStatResponse");

TH1D* h1d_SameMCIntStat=(TH1D*)fMC.Get("h1d_MCIntStat");
TH1D* h1d_SameMCIncStat=(TH1D*)fMC.Get("h1d_MCIncStat");
TH1D* recoUnfoldSameIntE=(TH1D*)fMC.Get("fullCorrRecoHist");
TH1D* recoUnfoldSameIncE=(TH1D*)fMC.Get("fullCorrIncHist");

TH1D* recoSameXSec=(TH1D*)fMC.Get("crossSectionUnfoldMC");
TH2D* covMatSameMCStat=(TH2D*)fMC.Get("covMatMCStat");
TH1D* recoXSec=(TH1D*)fMC.Get("crossSectionUnfoldData");
TH1D* dataXSec=(TH1D*)f.Get("crossSectionUnfoldData");
TH1D* recoIntOnlyXSec=(TH1D*)f.Get("crossSectionHist");
TH1D* dataIntOnlyXSec=(TH1D*)f.Get("crossSectionData");

TH1D* trueXSec=(TH1D*)f.Get("crossSectionTrue");



TH2D* response=(TH2D*)f.Get("response");
TH2D* responseIncident=(TH2D*)f.Get("responseIncident");
TH1D* trueIntE=(TH1D*)f.Get("h1d_true_interactingE");
TH1D* trueIncE=(TH1D*)f.Get("h1d_true_incidentE");
TH1D* recoIncE=(TH1D*)fMC.Get("h1d_reco_incidentE_data");
TH1D* dataIncE=(TH1D*)f.Get("h1d_reco_incidentE_data");
TH1D* recoIntE=(TH1D*)fMC.Get("h1d_reco_interactingE_data");
TH1D* dataIntE=(TH1D*)f.Get("h1d_reco_interactingE_data");
TH1D* recoUnfoldIntE=(TH1D*)fMC.Get("fullCorrData");
TH1D* dataUnfoldIntE=(TH1D*)f.Get("fullCorrData");



TH1D* totUncertainty=(TH1D*)dataUnfoldIntE->Clone("totUncertainty");
TH1D* systUncertainty=(TH1D*)dataUnfoldIntE->Clone("systUncertainty");
TH1D* statResUncertainty=(TH1D*)dataUnfoldIntE->Clone("statResUncertainty");
TH1D* statHistUncertainty=(TH1D*)dataUnfoldIntE->Clone("statHistUncertainty");

systUncertainty->SetLineColor(kRed); systUncertainty->SetLineStyle(2);
statResUncertainty->SetLineColor(kGreen); statResUncertainty->SetLineStyle(3);
statHistUncertainty->SetLineColor(kBlue); statHistUncertainty->SetLineStyle(4);
TH1D* recoUnfoldIncE=(TH1D*)fMC.Get("fullCorrIncData");
TH1D* dataUnfoldIncE=(TH1D*)f.Get("fullCorrIncData");

for(int i=0; i<dataUnfoldIncE->GetNbinsX(); ++i){

dataXSec->SetBinError(i+1, TMath::Sqrt(covMatMCRespStat->GetBinContent(i+1,i+1)+covMatSyst->GetBinContent(i+1,i+1)+covMatStat->GetBinContent(i+1,i+1)));

totUncertainty->SetBinContent(i+1,(dataXSec->GetBinError(i+1)/dataXSec->GetBinContent(i+1)));
systUncertainty->SetBinContent(i+1,TMath::Sqrt(covMatSyst->GetBinContent(i+1,i+1))/dataXSec->GetBinContent(i+1));
statResUncertainty->SetBinContent(i+1,TMath::Sqrt(covMatMCRespStat->GetBinContent(i+1,i+1))/dataXSec->GetBinContent(i+1));
statHistUncertainty->SetBinContent(i+1,TMath::Sqrt(covMatStat->GetBinContent(i+1,i+1))/dataXSec->GetBinContent(i+1));


std::cout<<dataXSec->GetBinContent(i+1)<<"& "<<dataXSec->GetBinError(i+1)<<" & "<<TMath::Sqrt(covMatSyst->GetBinContent(i+1,i+1))<<" & "<<TMath::Sqrt(covMatStat->GetBinContent(i+1,i+1))<<"&"<<TMath::Sqrt(covMatMCRespStat->GetBinContent(i+1,i+1))<<std::endl;//" & "<<TMath::Sqrt(covMatStatResponse->GetBinContent(i+1,i+1))<<std::endl;


recoXSec->SetBinError(i+1, TMath::Sqrt(covMatMCRespStatMC->GetBinContent(i+1,i+1)+covMatMCStat->GetBinContent(i+1,i+1)));

//std::cout<<covMatSyst->GetBinContent(i+1,i+1)<<","<<covMatStat->GetBinContent(i+1,i+1)<<","<<TMath::Sqrt(covMatSyst->GetBinContent(i+1,i+1)+covMatStat->GetBinContent(i+1,i+1))<<std::endl;

dataUnfoldIntE->SetBinError(i+1, TMath::Sqrt(TMath::Power(h1d_DataIntResp->GetBinError(i+1),2)+TMath::Power(h1d_DataIntStat->GetBinError(i+1),2)+TMath::Power(h1d_DataIntSyst->GetBinError(i+1),2)));
dataUnfoldIncE->SetBinError(i+1, TMath::Sqrt(TMath::Power(h1d_DataIncResp->GetBinError(i+1),2)+TMath::Power(h1d_DataIncStat->GetBinError(i+1),2)+TMath::Power(h1d_DataIncSyst->GetBinError(i+1),2)));

//std::cout<<covMatMCStatResponse->GetBinContent(i+1,i+1)<<","<<covMatStatResponse->GetBinContent(i+1,i+1)<<std::endl;

recoUnfoldIntE->SetBinError(i+1, TMath::Sqrt(TMath::Power(h1d_RecoIntResp->GetBinError(i+1),2)+TMath::Power(h1d_MCIntStat->GetBinError(i+1),2)));
recoUnfoldIncE->SetBinError(i+1, TMath::Sqrt(TMath::Power(h1d_RecoIncResp->GetBinError(i+1),2)+TMath::Power(h1d_MCIncStat->GetBinError(i+1),2)));










}


dataUnfoldIntE->SetLineColor(kBlack);
dataUnfoldIncE->SetLineColor(kBlack);
dataIntE->SetLineColor(kBlack);
dataIncE->SetLineColor(kBlack);
dataUnfoldIntE->SetMarkerColor(kBlack);
dataUnfoldIncE->SetMarkerColor(kBlack);
dataIntE->SetMarkerColor(kBlack);
dataIncE->SetMarkerColor(kBlack);
dataXSec->SetLineColor(kBlack);
dataXSec->SetMarkerColor(kBlack);
dataIntOnlyXSec->SetLineColor(kBlack);
dataIntOnlyXSec->SetMarkerColor(kBlack);

recoUnfoldIntE->SetLineColor(kRed);
recoUnfoldIncE->SetLineColor(kRed);
recoUnfoldSameIntE->SetLineColor(kRed);
recoUnfoldSameIncE->SetLineColor(kRed);
recoUnfoldSameIntE->SetMarkerColor(kRed);
recoUnfoldSameIncE->SetMarkerColor(kRed);
recoIntE->SetLineColor(kRed);
recoIncE->SetLineColor(kRed);
recoUnfoldIntE->SetMarkerColor(kRed);
recoUnfoldIncE->SetMarkerColor(kRed);
recoIntE->SetMarkerColor(kRed);
recoIncE->SetMarkerColor(kRed);
recoXSec->SetLineColor(kRed);
recoXSec->SetMarkerColor(kRed);
recoSameXSec->SetLineColor(kRed);
recoSameXSec->SetMarkerColor(kRed);
recoIntOnlyXSec->SetLineColor(kRed);
recoIntOnlyXSec->SetMarkerColor(kRed);


recoIntE->SetMarkerStyle(4);
recoIncE->SetMarkerStyle(4);
recoXSec->SetMarkerStyle(4);
recoUnfoldIntE->SetMarkerStyle(4);
recoUnfoldIncE->SetMarkerStyle(4);
recoUnfoldSameIntE->SetMarkerStyle(4);
recoUnfoldSameIncE->SetMarkerStyle(4);
recoSameXSec->SetMarkerStyle(4);
trueIntE->SetLineColor(kBlue);
trueIntE->SetMarkerColor(kBlue);
trueIncE->SetLineColor(kBlue);
trueIncE->SetMarkerColor(kBlue);
trueXSec->SetMarkerColor(kBlue);
trueXSec->SetLineColor(kBlue);

g4Truth->SetLineColor(kGreen);
g4Truth->SetMarkerColor(kGreen);

gStyle->SetOptFit(0);
  TLegend *lErr = new TLegend(0.3,0.7,0.5,0.9);
  lErr->SetTextFont(133);
  lErr->SetTextSize(25);
  lErr->AddEntry(totUncertainty,"Total Uncertainty","l");
  lErr->AddEntry(systUncertainty,"Systematic Uncertainty","l");
  lErr->AddEntry(statResUncertainty,"Res. Mat. Stat. Uncertainty","l");
  lErr->AddEntry(statHistUncertainty,"Input Hist. Stat. Uncertainty","l");



  TLegend *l = new TLegend(0.4,0.7,0.65,0.9);
  l->SetTextFont(133);
  l->SetTextSize(25);
  l->AddEntry(dataIntE,"Data (stat.+syst.)","lp");
  l->AddEntry(recoIntE,"Simulation (stat.)","lp");


  TLegend *ldata = new TLegend(0.45,0.7,0.65,0.9);
  ldata->SetTextFont(133);
  ldata->SetTextSize(25);
  ldata->AddEntry(dataIntE,"Data (stat.+syst.)","l");



  TLegend *lt = new TLegend(0.45,0.7,0.65,0.9);
  lt->SetTextFont(133);
  lt->SetTextSize(25);
  lt->AddEntry(dataIntE,"Data (stat.)","lp");
  lt->AddEntry(recoIntE,"Simulation","lp");

  TLegend *lMC= new TLegend(0.25,0.7,0.45,0.9);
  lMC->SetTextFont(133);
  lMC->SetTextSize(25);
  lMC->AddEntry(recoIntE,"Simulation (stat.)","lp");

  TLegend *lTrue= new TLegend(0.45,0.7,0.65,0.9);
  lTrue->SetTextFont(133);
  lTrue->SetTextSize(25);
  lTrue->AddEntry(trueIntE,"True Sim. (stat.)","l");



  TLegend *lTrue2= new TLegend(0.45,0.7,0.65,0.9);
  lTrue2->SetTextFont(133);
  lTrue2->SetTextSize(25);
  lTrue2->AddEntry(trueIntE,"True Sim. (stat.)","l");







recoIntE->GetXaxis()->SetTitle("Reconstructed Interacting Kinetic Energy [MeV]");
recoIntE->GetYaxis()->SetTitle("Number of Slices");
recoIntE->SetTitle("6 GeV/c Sample: Reco. Interacting Points");

recoIntE->Scale(dataIntE->GetEntries()/recoIntE->GetEntries());
//dataIntE->Scale(1.f/dataIntE->GetEntries());
recoIntE->GetYaxis()->SetRangeUser(0,recoIntE->GetMaximum()*2.0);
recoIntE->GetYaxis()->CenterTitle();
recoIntE->GetXaxis()->CenterTitle();
recoIntE->Draw("HIST");
dataIntE->Draw("E0 P SAME");
lt->Draw("SAME");
    //TLatex tL;
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c1.Print("finalPlots/interactingEnergyReco6GeVStandardGen.png");


recoIncE->GetXaxis()->SetTitle("Reconstructed Inc. Kinetic Energy [MeV]");
recoIncE->GetYaxis()->SetTitle("Number of Slices");
recoIncE->SetTitle("6 GeV/c Sample: Reco. Incident Points");

recoIncE->Scale(dataIncE->GetEntries()/recoIncE->GetEntries());
//dataIncE->Scale(1.f/dataIncE->GetEntries());

recoIncE->GetYaxis()->SetRangeUser(0,recoIncE->GetMaximum()*2.0);
recoIncE->GetYaxis()->CenterTitle();
recoIncE->GetXaxis()->CenterTitle();
recoIncE->Draw("HIST");
dataIncE->Draw("E0 P SAME");
lt->Draw("SAME");

    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c1.Print("finalPlots/incidentEnergyReco6GeVStandardGen.png");


  l->AddEntry(trueIntE,"True Sim. (stat.)","l");
 // lMC->AddEntry(trueIntE,"True Sim. (stat.)","l");

trueIntE->GetXaxis()->SetTitle("Unfolded Interacting Kinetic Energy [MeV]");
trueIntE->GetYaxis()->SetTitle("Number of Slices");
TH1D* missIntCopy=(TH1D*)missInt->Clone("missIntCopy");
TH1D* missIncCopy=(TH1D*)missInc->Clone("missIncCopy");
missIntCopy->Scale(1.f/trueIntE->Integral());
missIncCopy->Scale(1.f/trueIncE->Integral());
trueIntE->Scale(dataUnfoldIntE->Integral()/trueIntE->Integral());

recoUnfoldIntE->Scale(dataUnfoldIntE->Integral()/recoUnfoldIntE->Integral());
//dataUnfoldIntE->Scale(1.f/dataUnfoldIntE->Integral());


trueIntE->SetTitle("6 GeV/c Sample: Int. Point Energy");
trueIntE->GetYaxis()->SetRangeUser(0,trueIntE->GetMaximum()*2.0);
trueIntE->GetYaxis()->CenterTitle();
trueIntE->GetXaxis()->CenterTitle();
trueIntE->Draw("HIST");
recoUnfoldIntE->Draw("E1 P SAME");
dataUnfoldIntE->Draw("E0 P SAME");
l->Draw("SAME");

    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c1.Print("finalPlots/interactingEnergyUnfolded6GeVStandardGen.png");

lTrue2->AddEntry(missIntCopy,"Unselected (stat.)","l");
trueIntE->GetXaxis()->SetTitle("True Kinetic Energy [MeV]");
trueIntE->Draw("E0 P");
missIntCopy->Draw("E0 P SAME");
//recoUnfoldIntE->Draw("E0 P SAME");
//dataUnfoldIntE->Draw("E0 P SAME");
lTrue2->Draw("SAME");

    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c1.Print("finalPlots/interactingEnergyMCOnly6GeVStandardGen.png");



trueIncE->GetXaxis()->SetTitle("Unfolded Inc. Kinetic Energy [MeV]");
trueIncE->GetYaxis()->SetTitle("Number of Slices");
trueIncE->Scale(dataUnfoldIncE->Integral()/trueIncE->Integral());
recoUnfoldIncE->Scale(dataUnfoldIncE->Integral()/recoUnfoldIncE->Integral());
//dataUnfoldIncE->Scale(1.f/dataUnfoldIncE->Integral());
trueIncE->SetTitle("6 GeV/c Sample: Inc. Slices");
trueIncE->GetYaxis()->SetRangeUser(0,trueIncE->GetMaximum()*2.0);
trueIncE->GetYaxis()->CenterTitle();
trueIncE->GetXaxis()->CenterTitle();
trueIncE->Draw("HIST");
recoUnfoldIncE->Draw("E1 P SAME");
dataUnfoldIncE->Draw("E0 P SAME");
l->Draw("SAME");

    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c1.Print("finalPlots/incidentEnergyUnfolded6GeVStandardGen.png");
trueIncE->GetXaxis()->SetTitle("True Kinetic Energy [MeV]");
trueIncE->Draw("E0 P");
missIncCopy->Draw("E0 P SAME");
lTrue2->Draw("SAME");

    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c1.Print("finalPlots/incidentEnergyMCOnly6GeVStandardGen.png");






g4Truth->SetLineStyle(1);




  lTrue->AddEntry(g4Truth,"Geant4 v4.10.6 Bertini","l");
  l->AddEntry(g4Truth,"Geant4 v4.10.6 Bertini","l");
  lMC->AddEntry(g4Truth,"Geant4 v4.10.6 Bertini","l");
  ldata->AddEntry(g4Truth,"Geant4 v4.10.6 Bertini","l");
  hA_hist->SetLineColor(kRed);
hN_hist->SetLineColor(kBlue);
 hA_hist->SetLineStyle(2);
  hN_hist->SetLineStyle(3); 
 ldata->AddEntry(hA_hist,"GENIE v3.2 hA2018","l");
  ldata->AddEntry(hN_hist,"GENIE v3.2 hN2018","l");
recoXSec->GetXaxis()->SetTitle("Unfolded Kinetic Energy [MeV]");
recoXSec->SetTitle("6 GeV/c Sample: K^{+}-Ar Tot. Inel.");
recoXSec->GetYaxis()->SetTitle("#sigma [mbarn]");
TH1D* copyRecoXSec=(TH1D*)recoXSec->Clone("copy_XSec");
copyRecoXSec->GetYaxis()->SetRangeUser(300,800);
//trueXSec->GetXaxis()->SetRangeUser(4500,6000);
copyRecoXSec->GetXaxis()->CenterTitle();
copyRecoXSec->GetYaxis()->CenterTitle();
copyRecoXSec->Draw("E0 P");
g4Truth->Draw("SAME");
trueXSec->Draw("HIST SAME");
recoXSec->Draw("E0 P SAME");
dataXSec->Draw("E0 P SAME");
dataXSec->Fit("pol0","NO");
recoXSec->Fit("pol0","NO");
l->Draw("SAME");

    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c1.Print("finalPlots/xSecKaonUnfoldedAll6GeVStandardGen.png");




recoXSec->GetXaxis()->SetTitle("Unfolded Kinetic Energy [MeV]");
recoXSec->SetTitle("6 GeV/c Sample: K^{+}-Ar Tot. Inel.");
recoXSec->GetYaxis()->SetTitle("#sigma [mbarn]");
//TH1D* copyRecoXSec=(TH1D*)recoXSec->Clone("copy_XSec");
recoXSec->GetYaxis()->SetRangeUser(300,600);
//XSec->GetXaxis()->SetRangeUser(4500,6000);
recoXSec->GetXaxis()->CenterTitle();
recoXSec->GetYaxis()->CenterTitle();


recoXSec->Draw("E0 P");
//trueXSec->Draw("E0 P SAME");
g4Truth->Draw("SAME");
copyRecoXSec->Draw("E0 P SAME");
//dataXSec->Draw("E0 P SAME");
lMC->Draw("SAME");

    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c1.Print("finalPlots/xSecKaonMCAll6GeVStandardGen.png");




trueXSec->GetXaxis()->SetTitle("True Kinetic Energy [MeV]");
trueXSec->SetTitle("6 GeV/c Sample: K^{+}-Ar Tot. Inel.");
trueXSec->GetYaxis()->SetTitle("#sigma [mbarn]");
TH1D* trueRecoXSec=(TH1D*)trueXSec->Clone("copy_TrueXSec");
trueXSec->GetYaxis()->SetRangeUser(300,600);
//trueXSec->GetXaxis()->SetRangeUser(4500,6000);
trueXSec->GetXaxis()->CenterTitle();
trueXSec->GetYaxis()->CenterTitle();


trueXSec->Draw("E0 P");
//trueXSec->Draw("E0 P SAME");
g4Truth->Draw("SAME");
trueRecoXSec->Draw("E0 P SAME");
//dataXSec->Draw("E0 P SAME");
lTrue->Draw("SAME");

    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c1.Print("finalPlots/xSecKaonTrueAll6GeVStandardGen.png");


dataXSec->GetXaxis()->SetTitle("Unfolded Kinetic Energy [MeV]");
dataXSec->SetTitle("6 GeV/c Sample: K^{+}-Ar Tot. Inel.");
dataXSec->GetYaxis()->SetTitle("#sigma [mbarn]");
TH1D* dataRecoXSec=(TH1D*)dataXSec->Clone("data_XSec");
//trueXSec->GetXaxis()->SetRangeUser(4500,6000);
dataRecoXSec->GetYaxis()->SetRangeUser(200,650);
dataRecoXSec->GetXaxis()->CenterTitle();
dataRecoXSec->GetYaxis()->CenterTitle();
dataRecoXSec->Draw("E0 P");
g4Truth->Draw("SAME");
hA_hist->Draw("SAME");
hN_hist->Draw("SAME");
dataXSec->Draw("E0 P SAME");
//dataXSec->Draw("E0 P SAME");
ldata->Draw("SAME");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ProtoDUNE-SP}");
//l->Draw("SAME");
c1.Print("finalPlots/xSecKaonUnfoldedDataOnly6GeVStandardGen.png");


recoIntOnlyXSec->GetXaxis()->SetTitle("Unfolded Interacting Kinetic Energy [MeV]");
recoIntOnlyXSec->SetTitle("6 GeV/c Sample: K^{+}-Ar Tot. Inel.");
recoIntOnlyXSec->GetYaxis()->SetTitle("#sigma [mbarn]");
TH1D* copyRecoIntOnlyXSec=(TH1D*)recoIntOnlyXSec->Clone("copy_IntOnlyXSec");
//trueXSec->GetXaxis()->SetRangeUser(4500,6000);
copyRecoIntOnlyXSec->GetXaxis()->CenterTitle();
copyRecoIntOnlyXSec->GetYaxis()->CenterTitle();
copyRecoIntOnlyXSec->Draw("E0 P");
g4Truth->Draw("SAME");
trueXSec->Draw("E0 P SAME");
recoIntOnlyXSec->Draw("E0 P SAME");
dataIntOnlyXSec->Draw("E0 P SAME");
l->Draw("SAME");
c1.Print("finalPlots/xSecKaonUnfoldedIntOnly6GeVStandardGen.png");




  TLegend *lSame = new TLegend(0.15,0.7,0.85,0.9);
  lSame->SetTextFont(133);
  lSame->SetTextSize(25);
  lSame->AddEntry(trueIntE,"True Sim.","l");
  lSame->AddEntry(recoUnfoldSameIncE,"Sim. from Response Mat. (raw hist. stat.)","lp");

trueIncE->GetXaxis()->SetTitle("Unfolded Inc. Kinetic Energy [MeV]");
trueIncE->GetYaxis()->SetTitle("Number of Slices");
trueIncE->Scale(dataUnfoldIncE->Integral()/trueIncE->Integral());
recoUnfoldSameIncE->Scale(dataUnfoldIncE->Integral()/recoUnfoldSameIncE->Integral());
//dataUnfoldIncE->Scale(1.f/dataUnfoldIncE->Integral());
trueIncE->SetTitle("6 GeV/c Sample: Inc. Slices");
trueIncE->GetYaxis()->SetRangeUser(0,trueIncE->GetMaximum()*2.0);
trueIncE->GetYaxis()->CenterTitle();
trueIncE->GetXaxis()->CenterTitle();
trueIncE->Draw("HIST");
recoUnfoldSameIncE->Draw("E0 P SAME");
lSame->Draw("SAME");

    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c1.Print("finalPlots/incidentEnergySameUnfolded6GeVStandardGen.png");


trueIntE->GetXaxis()->SetTitle("Unfolded Int. Kinetic Energy [MeV]");
trueIntE->GetYaxis()->SetTitle("Number of Slices");
trueIntE->Scale(dataUnfoldIntE->Integral()/trueIntE->Integral());
recoUnfoldSameIntE->Scale(dataUnfoldIntE->Integral()/recoUnfoldSameIntE->Integral());
//dataUnfoldIncE->Scale(1.f/dataUnfoldIncE->Integral());
trueIntE->SetTitle("6 GeV/c Sample: Int. Slices");
trueIntE->GetYaxis()->SetRangeUser(0,trueIntE->GetMaximum()*2.0);
trueIntE->GetYaxis()->CenterTitle();
trueIntE->GetXaxis()->CenterTitle();
trueIntE->Draw("HIST");
recoUnfoldSameIntE->Draw("E0 P SAME");
lSame->Draw("SAME");

    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c1.Print("finalPlots/interactingEnergySameUnfolded6GeVStandardGen.png");




trueXSec->GetXaxis()->SetTitle("Unfolded Kinetic Energy [MeV]");
trueXSec->SetTitle("6 GeV/c Sample: K^{+}-Ar Tot. Inel.");
trueXSec->GetYaxis()->SetTitle("#sigma [mbarn]");
//trueXSec->GetXaxis()->SetRangeUser(4500,6000);
trueXSec->GetYaxis()->SetRangeUser(200,650);
trueXSec->GetXaxis()->CenterTitle();
trueXSec->GetYaxis()->CenterTitle();
trueXSec->Draw("HIST");
recoSameXSec->Draw("E0 P SAME");
lSame->Draw("SAME");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ProtoDUNE-SP}");
//l->Draw("SAME");
c1.Print("finalPlots/xSecKaonUnfoldedMCOnlySame6GeVStandardGen.png");

TH2D* responseCopy=new TH2D("responseCopy","resonpseCopy",5,0,5,4,0,4);
TH2D* responseIncCopy=new TH2D("responseIncCopy","responseIncCopy",5,0,5,4,0,4);
//int totBins=5;
//double edges[7]={3000,5000,5200,5400,5600,5900,6800};
////double edges[6]={4000,5000,5240,5440,5700,6650};
double edges[5]={4480,5080,5340,5610,6170};



for(int i=0;i<response->GetNbinsX();i++){
double rowValue=0;
double rowValueInc=0;
rowValue=missInt->GetBinContent(i+1);
rowValueInc=missInc->GetBinContent(i+1);

for(int j=0;j<response->GetNbinsX();j++){
rowValue=rowValue+response->GetBinContent(j+1,i+1);
rowValueInc=rowValueInc+responseIncident->GetBinContent(j+1,i+1);


}
responseCopy->SetBinContent(1,i+1,missInt->GetBinContent(i+1)/rowValue);
responseIncCopy->SetBinContent(1,i+1,missInc->GetBinContent(i+1)/rowValueInc);
for(int j=0;j<response->GetNbinsX();j++){
if(response->GetBinContent(j+1,i+1)==0){
response->SetBinContent(j+1,i+1,0.00001);

}
if(responseIncident->GetBinContent(j+1,i+1)==0){
responseIncident->SetBinContent(j+1,i+1,0.00001);


}

response->SetBinContent(j+1,i+1,response->GetBinContent(j+1,i+1)/rowValue);
responseCopy->SetBinContent(j+2,i+1,response->GetBinContent(j+1,i+1));

responseIncident->SetBinContent(j+1,i+1,responseIncident->GetBinContent(j+1,i+1)/rowValueInc);
responseIncCopy->SetBinContent(j+2,i+1,responseIncident->GetBinContent(j+1,i+1));

}



 responseCopy->GetXaxis()->SetBinLabel(1,"Unselected");
//responseCopy->GetYaxis()->SetBinLabel(i+1,"Unselected");

responseIncCopy->GetXaxis()->SetBinLabel(1,"Unselected");
//responseIncCopy->GetYaxis()->SetBinLabel(i+1,"Unselected");




 responseCopy->GetXaxis()->SetBinLabel(i+2,Form("%1.0f-%1.0f",edges[i],edges[i+1])); 
responseCopy->GetYaxis()->SetBinLabel(i+1,Form("%1.0f-%1.0f",edges[i],edges[i+1]));


responseIncCopy->GetXaxis()->SetBinLabel(i+2,Form("%1.0f-%1.0f",edges[i],edges[i+1]));
responseIncCopy->GetYaxis()->SetBinLabel(i+1,Form("%1.0f-%1.0f",edges[i],edges[i+1]));









}







gStyle->SetPaintTextFormat("1.3f");


responseCopy->SetMarkerSize(1.5);
responseIncCopy->SetMarkerSize(1.5);


gStyle->SetPalette(kBird);
responseCopy->GetXaxis()->SetTitle("Reconstructed Kinetic Energy [MeV]");
responseCopy->GetYaxis()->SetTitle("True Kinetic Energy [MeV]");
responseCopy->GetXaxis()->CenterTitle();
responseCopy->GetYaxis()->CenterTitle();
responseCopy->GetYaxis()->LabelsOption("v");
responseIncCopy->GetYaxis()->LabelsOption("v");
gPad->SetLeftMargin(0.18);
gPad->SetRightMargin(0.15);
responseCopy->GetYaxis()->SetTitleOffset(1.9);
responseCopy->GetZaxis()->SetRangeUser(-1.0,1.0);
responseCopy->SetTitle("6 GeV/c Sample: Int. Resp. Matrix");
responseCopy->Draw("COLZ TEXT0");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c1.Print("finalPlots/responseMatrix6GeVStandardGen.png");




gStyle->SetPaintTextFormat("1.3f");



responseIncCopy->GetXaxis()->SetTitle("Reconstructed Kinetic Energy [MeV]");
responseIncCopy->GetYaxis()->SetTitle("True Kinetic Energy [MeV]");
responseIncCopy->GetXaxis()->CenterTitle();
responseIncCopy->GetYaxis()->CenterTitle();
responseIncCopy->GetZaxis()->SetRangeUser(-1.0,1.0);
gPad->SetLeftMargin(0.18);
gPad->SetRightMargin(0.15);
responseIncCopy->GetYaxis()->SetTitleOffset(1.9);

responseIncCopy->SetTitle("6 GeV/c Sample: Inc. Resp. Matrix");
responseIncCopy->Draw("COLZ TEXT0");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c1.Print("finalPlots/responseIncidentMatrix6GeVStandardGen.png");
gPad->SetLeftMargin(0);

trueIncE->GetXaxis()->SetTitle("True Incident Kinetic Energy [MeV]");
trueIncE->SetTitle("MC True Incident");
trueIncE->GetYaxis()->SetTitle("Number of Entries");
trueIncE->Draw("HIST");
//missIncCopy->Draw("HIST SAME");
c1.Print("finalPlots/incEnergyTrue6GeVStandardGen.png");

TCanvas c2=TCanvas();
totUncertainty->GetXaxis()->SetTitle("Unfolded Kinetic Energy [MeV]");
totUncertainty->GetYaxis()->SetTitle("Relative Uncertainty");
totUncertainty->SetTitle("6 GeV/c Sample");
totUncertainty->GetXaxis()->CenterTitle(); totUncertainty->GetYaxis()->CenterTitle();
gPad->SetLeftMargin(0.15);
totUncertainty->GetYaxis()->SetTitleOffset(1.4);
totUncertainty->SetTitle("6 GeV/c Sample: K^{+}-Ar Tot. Inel.");
totUncertainty->GetYaxis()->SetRangeUser(0,0.3);
totUncertainty->Draw("HIST");
systUncertainty->Draw("HIST SAME");
statHistUncertainty->Draw("HIST SAME");
statResUncertainty->Draw("HIST SAME");
lErr->Draw("SAME");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c2.Print("finalPlots/relativeUncertaintyData6GeVStandardGen.png");




  int totBins=dataXSec->GetNbinsX();
   TMatrixD erecoData(totBins,totBins);
   TMatrixD erecoMC(totBins, totBins);
   TMatrixD erecoDataStat(totBins, totBins);
   TMatrixD erecoSameMC(totBins,totBins);
TH2D* covMatDraw=(TH2D*)covMatStat->Clone("CovMatAll");
TH2D* corrMC=(TH2D*)covMatStat->Clone("CorrMatMC");
TH2D* corrData=(TH2D*)covMatDraw->Clone("CorrData");
TH2D* corrDataStat=(TH2D*)covMatStat->Clone("CorrDatatStat");

   for (int i = 0; i < totBins; ++i) {
    for (int j = 0; j < totBins; ++j) {
    covMatDraw->SetBinContent(i+1,j+1, covMatStat->GetBinContent(i+1,j+1)+covMatSyst->GetBinContent(i+1,j+1)+covMatMCRespStat->GetBinContent(i+1,j+1));
    covMatMCStat->SetBinContent(i+1,j+1,covMatMCStat->GetBinContent(i+1,j+1)+covMatMCRespStatMC->GetBinContent(1+i,j+1));
    covMatStat->SetBinContent(1+i,j+1,covMatStat->GetBinContent(i+1,j+1)/*+covMatMCRespStat->GetBinContent(1+i,j+1)*/);
}}
   for (int i = 0; i < totBins; ++i) {
    std::cout<<TMath::Sqrt(covMatDraw->GetBinContent(i+1,i+1))<<","<<TMath::Sqrt(covMatStat->GetBinContent(i+1,i+1))<<","<<TMath::Sqrt(covMatSyst->GetBinContent(i+1,i+1))<<","<<TMath::Sqrt(covMatMCRespStat->GetBinContent(i+1,i+1))<<std::endl;
    for (int j = 0; j < totBins; ++j) {
  //  covMatMCStat->SetBinContent(1+i,j+i,covMatMCStat->GetBinContent(i+1,j+1)+covMatMCRespStatMC->GetBinContent(i+1,j+1));
    
    erecoMC[i][j] = covMatMCStat->GetBinContent(i+1,j+1);//+covMatMCRespStatMC->GetBinContent(i+1,j+1);
    erecoData[i][j] = covMatStat->GetBinContent(i+1, j+1)+covMatSyst->GetBinContent(i+1, j+1)+covMatMCRespStat->GetBinContent(i+1,j+1);
    erecoDataStat[i][j]=covMatStat->GetBinContent(i+1,j+1); 
    erecoSameMC[i][j]=covMatSameMCStat->GetBinContent(i+1,j+1);
    std::cout<<TMath::Sqrt(erecoData[i][j])<<std::endl;
    corrMC->SetBinContent(i+1,j+1,covMatMCStat->GetBinContent(i+1,j+1)/TMath::Sqrt(covMatMCStat->GetBinContent(i+1,i+1)*covMatMCStat->GetBinContent(j+1,j+1)));
   corrDataStat->SetBinContent(i+1,j+1,covMatStat->GetBinContent(i+1,j+1)/TMath::Sqrt(covMatStat->GetBinContent(i+1,i+1)*covMatStat->GetBinContent(j+1,j+1)));
  corrData->SetBinContent(i+1,j+1,covMatDraw->GetBinContent(i+1,j+1)/TMath::Sqrt(covMatDraw->GetBinContent(i+1,i+1)*covMatDraw->GetBinContent(j+1,j+1)));
 }}





gStyle->SetPalette(kBird);
covMatDraw->GetXaxis()->SetTitle("Unfolded Kinetic Energy [MeV]");
covMatDraw->GetYaxis()->SetTitle("Unfolded Kinetic Energy [MeV]");
covMatDraw->GetXaxis()->CenterTitle();
covMatDraw->GetYaxis()->CenterTitle();
gPad->SetLeftMargin(0.15);
gPad->SetRightMargin(0.15);
covMatDraw->GetYaxis()->SetTitleOffset(1.3);

covMatDraw->SetTitle("6 GeV/c Sample: Covariance Matrix");
covMatDraw->Draw("COLZ");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ProtoDUNE-SP}");
//covMatDraw->Draw("COLZ");
c2.Print("finalPlots/covarianceData6GeVStandardGen.png");
tL.SetNDC();
tL.DrawLatex(0.2,0.94,"#bf{DUNE:ProtoDUNE-SP}");
covMatDraw->SetMarkerSize(1.5);
covMatStat->SetMarkerSize(1.5);
covMatSyst->SetMarkerSize(1.5);
covMatMCRespStat->SetMarkerSize(1.5);
covMatDraw->GetZaxis()->SetRangeUser(-covMatDraw->GetMaximum(),covMatDraw->GetMaximum());
covMatStat->GetZaxis()->SetRangeUser(-covMatStat->GetMaximum(),covMatStat->GetMaximum());
covMatSyst->GetZaxis()->SetRangeUser(-covMatSyst->GetMaximum(),covMatSyst->GetMaximum());
covMatMCRespStat->GetZaxis()->SetRangeUser(-covMatMCRespStat->GetMaximum(),covMatMCRespStat->GetMaximum());

covMatMCStat->GetZaxis()->SetRangeUser(-covMatMCStat->GetMaximum(),covMatMCStat->GetMaximum());

covMatDraw->Draw("COLZ TEXT0.f");
tL.SetNDC();
tL.DrawLatex(0.2,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c2.Print("finalPlots/covarianceDataText6GeVStandardGen.png");

covMatStat->GetXaxis()->SetTitle("Unfolded Kinetic Energy [MeV]");
covMatStat->GetYaxis()->SetTitle("Unfolded Kinetic Energy [MeV]");
covMatStat->SetTitle("6 GeV/c Sample: Stat. Data Cov.");
covMatStat->GetXaxis()->CenterTitle();
covMatStat->GetYaxis()->CenterTitle();
gPad->SetLeftMargin(0.15);
gPad->SetRightMargin(0.15);
covMatStat->GetYaxis()->SetTitleOffset(1.3);
tL.DrawLatex(0.2,0.94,"#bf{DUNE:ProtoDUNE-SP}");
covMatStat->Draw("COLZ TEXT0.f");
tL.SetNDC();
tL.DrawLatex(0.2,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c2.Print("finalPlots/covarianceDataStatText6GeVStandardGen.png");


covMatMCRespStat->GetXaxis()->SetTitle("Unfolded Kinetic Energy [MeV]");
covMatMCRespStat->GetYaxis()->SetTitle("Unfolded Kinetic Energy [MeV]");
covMatMCRespStat->SetTitle("6 GeV/c Sample: Stat. Resp. Cov.");
covMatMCRespStat->GetXaxis()->CenterTitle();
covMatMCRespStat->GetYaxis()->CenterTitle();
gPad->SetLeftMargin(0.15);
gPad->SetRightMargin(0.15);
covMatMCRespStat->GetYaxis()->SetTitleOffset(1.3);
tL.DrawLatex(0.2,0.94,"#bf{DUNE:ProtoDUNE-SP}");
covMatMCRespStat->Draw("COLZ TEXT0.f");
tL.SetNDC();
tL.DrawLatex(0.2,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c2.Print("finalPlots/covarianceRespStatText6GeVStandardGen.png");


covMatSyst->GetXaxis()->SetTitle("Unfolded Kinetic Energy [MeV]");
covMatSyst->GetYaxis()->SetTitle("Unfolded Kinetic Energy [MeV]");
covMatSyst->SetTitle("6 GeV/c Sample: Syst. Cov.");
covMatSyst->GetXaxis()->CenterTitle();
covMatSyst->GetYaxis()->CenterTitle();
gPad->SetLeftMargin(0.15);
gPad->SetRightMargin(0.15);
covMatSyst->GetYaxis()->SetTitleOffset(1.3);
tL.DrawLatex(0.2,0.94,"#bf{DUNE:ProtoDUNE-SP}");
covMatSyst->Draw("COLZ TEXT0.f");
tL.SetNDC();
tL.DrawLatex(0.2,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c2.Print("finalPlots/covarianceDataSystText6GeVStandardGen.png");



covMatMCStat->GetXaxis()->SetTitle("Unfolded Kinetic Energy [MeV]");
covMatMCStat->GetYaxis()->SetTitle("Unfolded Kinetic Energy [MeV]");
covMatMCStat->GetXaxis()->CenterTitle();
covMatMCStat->GetYaxis()->CenterTitle();
gPad->SetLeftMargin(0.15);
covMatMCStat->GetYaxis()->SetTitleOffset(1.3);

covMatMCStat->SetTitle("6 GeV/c Sample: Sim. Stat. Cov. Mat.");
covMatMCStat->Draw("COLZ TEXT0.f");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ProtoDUNE-SP}");
//covMatMCStat->Draw("COLZ TEXT0.f");
c2.Print("finalPlots/covarianceMCStatOnly6GeVStandardGen.png");




corrMC->GetZaxis()->SetRangeUser(-1.0,1.0);
corrDataStat->GetZaxis()->SetRangeUser(-1.0,1.0);
corrData->GetZaxis()->SetRangeUser(-1.0,1.0);


corrMC->SetMarkerSize(1.5);
corrDataStat->SetMarkerSize(1.5);
corrData->SetMarkerSize(1.5);
/*
corrMC->GetZaxis()->SetLabelFont(133);
   corrMC->GetZaxis()->SetLabelSize(50);
   corrMC->GetZaxis()->SetTitleSize(25);
   //correlation->GetZaxis()->SetLabelColor(kGray);
   corrMC->GetZaxis()->SetTitleFont(133);

corrDataStat->GetZaxis()->SetLabelFont(133);
   corrDataStat->GetZaxis()->SetLabelSize(50);
   corrDataStat->GetZaxis()->SetTitleSize(25);
   //correlation->GetZaxis()->SetLabelColor(kGray);
   corrDataStat->GetZaxis()->SetTitleFont(133);

corrData->GetZaxis()->SetLabelFont(133);
   corrData->GetZaxis()->SetLabelSize(50);
   corrData->GetZaxis()->SetTitleSize(25);
   //correlation->GetZaxis()->SetLabelColor(kGray);
   corrData->GetZaxis()->SetTitleFont(133);
*/

corrMC->GetXaxis()->SetTitle("Unfolded Kinetic Energy [MeV]");
corrMC->GetYaxis()->SetTitle("Unfolded Kinetic Energy [MeV]");
corrMC->GetXaxis()->CenterTitle();
corrMC->GetYaxis()->CenterTitle();
corrMC->GetYaxis()->SetTitleOffset(1.3);
gPad->SetLeftMargin(0.15);
gPad->SetRightMargin(0.12);
gStyle->SetPaintTextFormat("1.3f");
corrMC->SetTitle("6 GeV/c Sample: Sim. Stat. Corr.");
corrMC->Draw("COLZ TEXT0.f");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ProtoDUNE-SP}");
//corrMC->Draw("COLZ TEXT0.f");
c2.Print("finalPlots/correlationMCStatOnly6GeVStandardGen.png");

corrDataStat->GetXaxis()->SetTitle("Unfolded Kinetic Energy [MeV]");
corrDataStat->GetYaxis()->SetTitle("Unfolded Kinetic Energy [MeV]");
corrDataStat->GetXaxis()->CenterTitle();
corrDataStat->GetYaxis()->CenterTitle();
gPad->SetLeftMargin(0.15);
gPad->SetRightMargin(0.12);
corrDataStat->GetYaxis()->SetTitleOffset(1.3);

corrDataStat->SetTitle("6 GeV/c Sample: Stat. Corr. Matrix");
corrDataStat->Draw("COLZ TEXT0.f");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ProtoDUNE-SP}");
//corrDataStat->Draw("COLZ TEXT0.f");
c2.Print("finalPlots/correlationDataStatOnly6GeVStandardGen.png");

corrData->GetXaxis()->SetTitle("Unfolded Kinetic Energy [MeV]");
corrData->GetYaxis()->SetTitle("Unfolded Kinetic Energy [MeV]");
corrData->GetXaxis()->CenterTitle();
corrData->GetYaxis()->CenterTitle();
gPad->SetLeftMargin(0.15);
gPad->SetRightMargin(0.12);
corrData->GetYaxis()->SetTitleOffset(1.3);

corrData->SetTitle("6 GeV/c Sample: Corr. Matrix");
corrData->Draw("COLZ TEXT0.f");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ProtoDUNE-SP}");
//corrData->Draw("COLZ TEXT0.f");
c2.Print("finalPlots/correlationData6GeVStandardGen.png");



 
TDecompChol mat_decomp(erecoData);

bool didItWork=mat_decomp.Decompose();
//std::cout<<didItWork<<std::endl;
auto new_mat=mat_decomp.Invert();

TDecompChol mat_decompStat(erecoDataStat);
auto new_matStat=mat_decompStat.Invert();


TDecompChol mat_decompMC(erecoMC);
bool didItWorkMC=mat_decompMC.Decompose();
auto new_matMC=mat_decompMC.Invert();

TDecompChol mat_decompSameStat(erecoSameMC);
auto new_matSameMC=mat_decompSameStat.Invert();

double chi2=0;
double chi2MC=0;
double chi2SameMC=0;
double chi2Stat=0;
double chi2GENIE=0;
double chi2Truth=0;
    for(int bin=0; bin<dataXSec->GetNbinsX(); ++bin){
     chi2Truth=chi2Truth+TMath::Power(trueXSec->GetBinContent(bin+1)-g4Truth->Eval(dataXSec->GetBinCenter(bin+1)),2)/TMath::Power(trueXSec->GetBinError(bin+1),2);
     for(int bin2=0; bin2<dataXSec->GetNbinsX(); ++bin2){

    //std::cout<<bin2+1<<std::endl;
    double diffBin1=dataXSec->GetBinContent(bin+1)-g4Truth->Eval(dataXSec->GetBinCenter(bin+1));
    double diffBin2=dataXSec->GetBinContent(bin2+1)-g4Truth->Eval(dataXSec->GetBinCenter(bin2+1));
   // double diffBin2=averageData->GetBinContent(bin2+1)-averageMC->GetBinContent(bin2+1);
    double bin_cont=diffBin1*(new_mat)(bin, bin2)*diffBin2;
    double bin_contStat=diffBin1*(new_matStat)(bin, bin2)*diffBin2;
    //std::cout<<g4Truth->Eval(dataXSec->GetBinCenter(bin+1))<<std::endl;

    chi2=chi2+bin_cont;
    chi2Stat=chi2Stat+bin_contStat;
    double diffBin1MC=recoXSec->GetBinContent(bin+1)-g4Truth->Eval(recoXSec->GetBinCenter(bin+1));
    double diffBin2MC=recoXSec->GetBinContent(bin2+1)-g4Truth->Eval(recoXSec->GetBinCenter(bin2+1));
    double bin_contMC=diffBin1MC*(new_matMC)(bin,bin2)*diffBin2MC;
    chi2MC=chi2MC+bin_contMC;

    double diffBin1GENIE=dataXSec->GetBinContent(bin+1)-hA_hist->Interpolate(dataXSec->GetBinCenter(bin+1));
    double diffBin2GENIE=dataXSec->GetBinContent(bin2+1)-hA_hist->Interpolate(dataXSec->GetBinCenter(bin2+1));
    double bin_contGENIE=diffBin1GENIE*(new_mat)(bin,bin2)*diffBin2GENIE;
    chi2GENIE=chi2GENIE+bin_contGENIE;


    double diffBin1SameMC=recoSameXSec->GetBinContent(bin+1)-g4Truth->Eval(recoSameXSec->GetBinCenter(bin+1));
    double diffBin2SameMC=recoSameXSec->GetBinContent(bin2+1)-g4Truth->Eval(recoSameXSec->GetBinCenter(bin2+1));
    double bin_contSameMC=diffBin1SameMC*(new_matSameMC)(bin,bin2)*diffBin2SameMC;
    chi2SameMC=chi2SameMC+bin_contSameMC;





}
}





std::cout<<"Chi2 (G4, GENIE, Data Stat. to G4,MC Stat to G4, dof): "<<chi2<<','<<chi2GENIE<<","<<chi2Stat<<","<<chi2MC<<","<<dataXSec->GetNbinsX()<<std::endl;
std::cout<<"Chi2 truth "<<chi2Truth<<std::endl;
std::cout<<"Chi2 Same MC+Resp: "<<chi2SameMC<<std::endl;
std::vector<double> data, mc, err;
for(int i=0; i<dataXSec->GetNbinsX(); ++i){
data.push_back(dataXSec->GetBinContent(i+1));
mc.push_back(g4Truth->Eval(dataXSec->GetBinCenter(i+1)));
err.push_back(dataXSec->GetBinError(i+1));









}
/*
TFile fEff("effHistogram.root");
TH1D* effHist=(TH1D*)fEff.Get("effHist");
TH1D* effHist_close=(TH1D*)fEff.Get("effHist_close");
effHist->GetXaxis()->CenterTitle();
effHist->GetYaxis()->CenterTitle();
effHist_close->GetXaxis()->CenterTitle();
effHist_close->GetYaxis()->CenterTitle();



effHist->Draw("e0 P");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c1.Print("efficiencyPlot6GeV.png");


effHist_close->Draw("e0 P");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c1.Print("efficiencyClosePlot6GeV.png");



std::cout<<pearson_chi2(data, mc, err)<<","<<mc.size()<<std::endl;


std::vector<double> data2, mc2, err2;
double dataAvg;
for(int i=2; i<8; ++i){
data2.push_back(dataXSec->GetBinContent(i+1));
mc2.push_back(g4Truth->Eval(dataXSec->GetBinCenter(i+1)));
err2.push_back(dataXSec->GetBinError(i+1));

dataAvg=dataAvg+dataXSec->GetBinContent(i+1)/6;
}



std::cout<<pearson_chi2(data2, mc2, err2)<<","<<mc2.size()<<std::endl;
std::cout<<dataAvg<<std::endl;

std::vector<double> data3, mc3, err3;

double par;
double start=0.8;
for(int j=0; j<50; j++){
double dataAvg=0;
data3.clear();
mc3.clear();
err3.clear();
for(int i=0; i<dataXSec->GetNbinsX(); ++i){
par=j*0.005+start;
data3.push_back(dataXSec->GetBinContent(i+1));
mc3.push_back(par*g4Truth->Eval(dataXSec->GetBinCenter(i+1)));
err3.push_back(dataXSec->GetBinError(i+1));

dataAvg=dataAvg+dataXSec->GetBinContent(i+1)/6;
}
std::cout<<pearson_chi2(data3, mc3, err3)<<','<<j<<","<<par<<std::endl;
}
*/
}

