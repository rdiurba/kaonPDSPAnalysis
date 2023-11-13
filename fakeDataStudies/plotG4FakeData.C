
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
void plotUnfoldResults()
{

gROOT->LoadMacro("protoDUNEStyle.C");
gROOT->SetStyle("protoDUNEStyle");
gROOT->ForceStyle();
gStyle->SetTitleX(0.35);
gStyle->SetOptFit(111);
//gStyle->SetPadRightMargin(0.15);
//gStyle->SetPadLeftMargin(0.15);
TCanvas c1=TCanvas();

//TFile f("/dune/data/users/rdiurba/kaonana_run5770.root");
TFile f("kaonUnfoldRecoPlots.root");
//TFile fThrow("kaonSystShiftOutput_1kThrows.root");
TFile fThrow("kaonSystShiftOutput_1kThrows_caloShift.root");
TFile fStat("kaonUnfoldWithStatShifts.root");
TFile fResponse("kaonStatResponseShiftOutput_1kThrows.root");
TFile fTruth("/dune/data/users/rdiurba/rootDump/kaon_cross_section_out.root");
TGraph* g4Truth=(TGraph*)fTruth.Get("inel_KE");
TH2D* covMatSyst=(TH2D*)fThrow.Get("covMat");
TH2D* covMatStat=(TH2D*)fStat.Get("covMatStat");
TH2D* covMatMCStat=(TH2D*)fStat.Get("covMatMCStat");
TH2D* covMatStatResponse=(TH2D*)fResponse.Get("covMatResponse");
TH2D* covMatMCStatResponse=(TH2D*)fResponse.Get("covMatMCResponse");
TH1D* h1d_DataIntSyst=(TH1D*)fThrow.Get("h1d_DataInt");
TH1D* h1d_DataIncSyst=(TH1D*)fThrow.Get("h1d_DataInc");

TH1D* h1d_DataIntStat=(TH1D*)fStat.Get("h1d_DataIntStat");
TH1D* h1d_DataIncStat=(TH1D*)fStat.Get("h1d_DataIncStat");


TH1D* h1d_DataIntStatResponse=(TH1D*)fResponse.Get("h1d_DataIntStatResponse");
TH1D* h1d_DataIncStatResponse=(TH1D*)fResponse.Get("h1d_DataIncStatResponse");


TH1D* h1d_MCIntStatResponse=(TH1D*)fResponse.Get("h1d_MCIntStatResponse");
TH1D* h1d_MCIncStatResponse=(TH1D*)fResponse.Get("h1d_MCIncStatResponse");


TH1D* recoXSec=(TH1D*)f.Get("crossSectionUnfoldMC");
TH1D* dataXSec=(TH1D*)f.Get("crossSectionUnfoldData");
TH1D* recoIntOnlyXSec=(TH1D*)f.Get("crossSectionHist");
TH1D* dataIntOnlyXSec=(TH1D*)f.Get("crossSectionData");

TH1D* trueXSec=(TH1D*)f.Get("crossSectionTrue");
TH2D* response=(TH2D*)f.Get("response");
TH2D* responseIncident=(TH2D*)f.Get("responseIncident");
TH1D* trueIntE=(TH1D*)f.Get("h1d_true_interactingE");
TH1D* trueIncE=(TH1D*)f.Get("h1d_true_incidentE");
TH1D* recoIncE=(TH1D*)f.Get("h1d_reco_incidentE");
TH1D* dataIncE=(TH1D*)f.Get("h1d_reco_incidentE_data");
TH1D* recoIntE=(TH1D*)f.Get("h1d_reco_interactingE");
TH1D* dataIntE=(TH1D*)f.Get("h1d_reco_interactingE_data");
TH1D* recoUnfoldIntE=(TH1D*)f.Get("fullCorrRecoHist");
TH1D* dataUnfoldIntE=(TH1D*)f.Get("fullCorrData");


TH1D* recoUnfoldIncE=(TH1D*)f.Get("fullCorrIncHist");
TH1D* dataUnfoldIncE=(TH1D*)f.Get("fullCorrIncData");

for(int i=0; i<dataUnfoldIncE->GetNbinsX(); ++i){

dataXSec->SetBinError(i+1, TMath::Sqrt(covMatStatResponse->GetBinContent(i+1,i+1)+covMatSyst->GetBinContent(i+1,i+1)+covMatStat->GetBinContent(i+1,i+1)));


recoXSec->SetBinError(i+1, TMath::Sqrt(covMatMCStatResponse->GetBinContent(i+1,i+1)+covMatMCStat->GetBinContent(i+1,i+1)));

std::cout<<covMatSyst->GetBinContent(i+1,i+1)<<","<<covMatStat->GetBinContent(i+1,i+1)<<","<<TMath::Sqrt(covMatSyst->GetBinContent(i+1,i+1)+covMatStat->GetBinContent(i+1,i+1))<<std::endl;

dataUnfoldIntE->SetBinError(i+1, TMath::Sqrt(TMath::Power(h1d_DataIntStat->GetBinError(i+1),2)+TMath::Power(h1d_DataIntSyst->GetBinError(i+1),2)));
dataUnfoldIncE->SetBinError(i+1, TMath::Sqrt(TMath::Power(h1d_DataIncStat->GetBinError(i+1),2)+TMath::Power(h1d_DataIncSyst->GetBinError(i+1),2)));

std::cout<<covMatMCStatResponse->GetBinContent(i+1,i+1)<<","<<covMatStatResponse->GetBinContent(i+1,i+1)<<std::endl;
//recoUnfoldIntE->SetBinError(i+1, TMath::Sqrt(TMath::Power(h1d_MCIntStatResponse->GetBinError(i+1),2)+TMath::Power(recoUnfoldIntE->GetBinError(i+1),2)));

//recoUnfoldIncE->SetBinError(i+1, TMath::Sqrt(TMath::Power(h1d_MCIncStatResponse->GetBinError(i+1),2)+TMath::Power(recoUnfoldIncE->GetBinError(i+1),2)));



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
recoIntE->SetLineColor(kRed);
recoIncE->SetLineColor(kRed);
recoUnfoldIntE->SetMarkerColor(kRed);
recoUnfoldIncE->SetMarkerColor(kRed);
recoIntE->SetMarkerColor(kRed);
recoIncE->SetMarkerColor(kRed);
recoXSec->SetLineColor(kRed);
recoXSec->SetMarkerColor(kRed);
recoIntOnlyXSec->SetLineColor(kRed);
recoIntOnlyXSec->SetMarkerColor(kRed);

trueIntE->SetLineColor(kBlue);
trueIntE->SetMarkerColor(kBlue);
trueIncE->SetLineColor(kBlue);
trueIncE->SetMarkerColor(kBlue);
trueXSec->SetMarkerColor(kBlue);
trueXSec->SetLineColor(kBlue);

g4Truth->SetLineColor(kGreen);
g4Truth->SetMarkerColor(kGreen);

  TLegend *l = new TLegend(0.45,0.7,0.65,0.9);
  l->SetTextFont(133);
  l->SetTextSize(25);
  l->AddEntry(dataIntE,"Data (stat.+syst.)","l");
  l->AddEntry(recoIntE,"Simulation (stat.)","l");

  TLegend *lt = new TLegend(0.45,0.7,0.65,0.9);
  lt->SetTextFont(133);
  lt->SetTextSize(25);
  lt->AddEntry(dataIntE,"Data (stat.)","l");
  lt->AddEntry(recoIntE,"Simulation (stat.)","l");

  TLegend *lMC= new TLegend(0.45,0.7,0.65,0.9);
  lMC->SetTextFont(133);
  lMC->SetTextSize(25);
  lMC->AddEntry(recoIntE,"Simulation (stat.)","l");



recoIntE->GetXaxis()->SetTitle("Reconstructed Interacting Kinetic Energy [MeV]");
recoIntE->GetYaxis()->SetTitle("Arbitrary Units");
recoIntE->SetTitle("Reco. Interacting Points");

recoIntE->Scale(1.f/recoIntE->GetEntries());
dataIntE->Scale(1.f/dataIntE->GetEntries());
recoIntE->GetYaxis()->SetRangeUser(0,recoIntE->GetMaximum()*2.0);
recoIntE->GetYaxis()->CenterTitle();
recoIntE->GetXaxis()->CenterTitle();
recoIntE->Draw("HIST");
dataIntE->Draw("E0 P SAME");
lt->Draw("SAME");
    TLatex tL;
    tL.SetNDC();
    tL.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c1.Print("interactingEnergyReco.png");


recoIncE->GetXaxis()->SetTitle("Reconstructed Inc. Kinetic Energy [MeV]");
recoIncE->GetYaxis()->SetTitle("Arbitrary Units");
recoIncE->SetTitle("Reco. Incident Points");

recoIncE->Scale(1.f/recoIncE->GetEntries());
dataIncE->Scale(1.f/dataIncE->GetEntries());

recoIncE->GetYaxis()->SetRangeUser(0,recoIncE->GetMaximum()*2.0);
recoIncE->GetYaxis()->CenterTitle();
recoIncE->GetXaxis()->CenterTitle();
recoIncE->Draw("HIST");
dataIncE->Draw("E0 P SAME");
lt->Draw("SAME");

    tL.SetNDC();
    tL.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c1.Print("incidentEnergyReco.png");


  l->AddEntry(trueIntE,"True Sim.","l");
  lMC->AddEntry(trueIntE,"True Sim.","l");

trueIntE->GetXaxis()->SetTitle("Unfolded Interacting Kinetic Energy [MeV]");
trueIntE->GetYaxis()->SetTitle("Arbitrary Units");
trueIntE->Scale(1.f/trueIntE->Integral());
recoUnfoldIntE->Scale(1.f/recoUnfoldIntE->Integral());
dataUnfoldIntE->Scale(1.f/dataUnfoldIntE->Integral());
std::cout<<trueIntE->Integral()<<','<<dataUnfoldIntE->Integral()<<std::endl;

trueIntE->SetTitle("Interacting Point Energy");
trueIntE->GetYaxis()->SetRangeUser(0,trueIntE->GetMaximum()*2.0);
trueIntE->GetYaxis()->CenterTitle();
trueIntE->GetXaxis()->CenterTitle();
trueIntE->Draw("HIST");
recoUnfoldIntE->Draw("E0 P SAME");
dataUnfoldIntE->Draw("E0 P SAME");
l->Draw("SAME");

    tL.SetNDC();
    tL.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c1.Print("interactingEnergyUnfolded.png");


trueIntE->Draw("HIST");
recoUnfoldIntE->Draw("E0 P SAME");
//dataUnfoldIntE->Draw("E0 P SAME");
lMC->Draw("SAME");

    tL.SetNDC();
    tL.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c1.Print("interactingEnergyMCOnly.png");



trueIncE->GetXaxis()->SetTitle("Unfolded Inc. Kinetic Energy [MeV]");
trueIncE->GetYaxis()->SetTitle("Arbitrary Units");
trueIncE->Scale(1.f/trueIncE->Integral());
recoUnfoldIncE->Scale(1.f/recoUnfoldIncE->Integral());
dataUnfoldIncE->Scale(1.f/dataUnfoldIncE->Integral());
trueIncE->SetTitle("Inc. Hit Energy");
trueIncE->GetYaxis()->SetRangeUser(0,trueIncE->GetMaximum()*2.0);
trueIncE->GetYaxis()->CenterTitle();
trueIncE->GetXaxis()->CenterTitle();
trueIncE->Draw("HIST");
recoUnfoldIncE->Draw("E0 P SAME");
dataUnfoldIncE->Draw("E0 P SAME");
l->Draw("SAME");

    tL.SetNDC();
    tL.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c1.Print("incidentEnergyUnfolded.png");

trueIncE->Draw("HIST");
recoUnfoldIncE->Draw("E0 P SAME");
lMC->Draw("SAME");

    tL.SetNDC();
    tL.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c1.Print("incidentEnergyUnfoldedMCOnly.png");



g4Truth->SetLineStyle(1);


  l->AddEntry(g4Truth,"G4 Input","l");

recoXSec->GetXaxis()->SetTitle("Unfolded Kinetic Energy [MeV]");
recoXSec->SetTitle("Measured Cross Section");
recoXSec->GetYaxis()->SetTitle("#sigma [mbarn]");
TH1D* copyRecoXSec=(TH1D*)recoXSec->Clone("copy_XSec");
copyRecoXSec->GetYaxis()->SetRangeUser(100,800);
//trueXSec->GetXaxis()->SetRangeUser(4500,6000);
copyRecoXSec->GetXaxis()->CenterTitle();
copyRecoXSec->GetYaxis()->CenterTitle();
copyRecoXSec->Draw("E0 P");
g4Truth->Draw("SAME");
trueXSec->Draw("E0 P SAME");
recoXSec->Draw("E0 P SAME");
dataXSec->Draw("E0 P SAME");
l->Draw("SAME");

    tL.SetNDC();
    tL.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c1.Print("xSecKaonUnfoldedAll.png");

recoXSec->GetXaxis()->SetTitle("Unfolded Kinetic Energy [MeV]");
recoXSec->SetTitle("Measured Cross Section");
recoXSec->GetYaxis()->SetTitle("#sigma [mbarn]");
//trueXSec->GetXaxis()->SetRangeUser(4500,6000);
copyRecoXSec->GetYaxis()->SetRangeUser(100,800);
copyRecoXSec->GetXaxis()->CenterTitle();
copyRecoXSec->GetYaxis()->CenterTitle();
copyRecoXSec->Draw("E0 P");
g4Truth->Draw("SAME");
recoXSec->Draw("E0 P SAME");
dataXSec->Draw("E0 P SAME");

    tL.SetNDC();
    tL.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
//l->Draw("SAME");
c1.Print("xSecKaonUnfoldedReco.png");


recoIntOnlyXSec->GetXaxis()->SetTitle("Unfolded Interacting Kinetic Energy [MeV]");
recoIntOnlyXSec->SetTitle("Measured Cross Section");
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
c1.Print("xSecKaonUnfoldedIntOnly.png");


response->GetXaxis()->SetTitle("Reconstructed Kinetic Energy [MeV]");
response->GetYaxis()->SetTitle("True Kinetic Energy [MeV]");
response->GetXaxis()->CenterTitle();
response->GetYaxis()->CenterTitle();
gPad->SetLeftMargin(0.15);
response->GetYaxis()->SetTitleOffset(1.4);

response->SetTitle("Smearing Interaction Point Matrix");
response->Draw("COLZ");
    tL.SetNDC();
    tL.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c1.Print("responseMatrix.png");


responseIncident->GetXaxis()->SetTitle("Reconstructed Kinetic Energy [MeV]");
responseIncident->GetYaxis()->SetTitle("True Kinetic Energy [MeV]");
responseIncident->GetXaxis()->CenterTitle();
responseIncident->GetYaxis()->CenterTitle();
gPad->SetLeftMargin(0.15);
responseIncident->GetYaxis()->SetTitleOffset(1.4);

responseIncident->SetTitle("Smearing Incident Hit Matrix");
responseIncident->Draw("COLZ");
    tL.SetNDC();
    tL.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c1.Print("responseIncidentMatrix.png");
gPad->SetLeftMargin(0);

trueIncE->GetXaxis()->SetTitle("True Incident Kinetic Energy [MeV]");
trueIncE->SetTitle("MC True Incident");
trueIncE->GetYaxis()->SetTitle("Number of Entries");
trueIncE->Draw("HIST");
c1.Print("incEnergyTrue.png");






covMatSyst->Add(covMatStat);
covMatSyst->Add(covMatStatResponse);


   TMatrixD* covMatrixTemp=new TMatrixD(dataXSec->GetNbinsX()+2,dataXSec->GetNbinsX()+2, covMatSyst->GetArray(),"");

   auto covMatrix=covMatrixTemp->GetSub(1,dataXSec->GetNbinsX(),1,dataXSec->GetNbinsX());
TDecompChol mat_decomp(covMatrix);

bool didItWork=mat_decomp.Decompose();
std::cout<<didItWork<<std::endl;
auto new_mat=mat_decomp.Invert();

double chi2;
    for(int bin=0; bin<dataXSec->GetNbinsX(); ++bin){
    for(int bin2=0; bin2<dataXSec->GetNbinsX(); ++bin2){

    std::cout<<bin2+1<<std::endl;
    double diffBin1=dataXSec->GetBinContent(bin+1)-g4Truth->Eval(dataXSec->GetBinCenter(bin+1));
    double diffBin2=dataXSec->GetBinContent(bin2+1)-g4Truth->Eval(dataXSec->GetBinCenter(bin2+1));
   // double diffBin2=averageData->GetBinContent(bin2+1)-averageMC->GetBinContent(bin2+1);
    double bin_cont=diffBin1*(new_mat)(bin, bin2)*diffBin2;
    std::cout<<g4Truth->Eval(dataXSec->GetBinCenter(bin+1))<<std::endl;

    chi2=chi2+bin_cont;
}
}

std::cout<<chi2<<','<<dataXSec->GetNbinsX()<<std::endl;
std::vector<double> data, mc, err;
for(int i=0; i<dataXSec->GetNbinsX(); ++i){
data.push_back(dataXSec->GetBinContent(i+1));
mc.push_back(g4Truth->Eval(dataXSec->GetBinCenter(i+1)));
err.push_back(dataXSec->GetBinError(i+1));


}



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
for(int i=2; i<8; ++i){
par=j*0.005+start;
data3.push_back(dataXSec->GetBinContent(i+1));
mc3.push_back(par*g4Truth->Eval(dataXSec->GetBinCenter(i+1)));
err3.push_back(dataXSec->GetBinError(i+1));

dataAvg=dataAvg+dataXSec->GetBinContent(i+1)/6;
}
std::cout<<pearson_chi2(data3, mc3, err3)<<','<<j<<","<<par<<std::endl;
}

}

