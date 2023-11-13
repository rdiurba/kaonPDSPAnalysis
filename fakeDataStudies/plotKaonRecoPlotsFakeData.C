
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
void plotKaonRecoPlots(int sel)
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
TFile f("kaonTreeSelDataOld.root");
TFile fMC("/dune/data/users/rdiurba/standard6GeVkaonTreeSel.root");
//TFile fMC("kaonana_mc_test.root");
//TFile fG4("/dune/data/users/rdiurba/kaonXSTruthG4.root","read")
TTree* tData=(TTree*)f.Get("ana");
TTree* tMC=(TTree*)fMC.Get("ana");

TH1D* beamPMC=new TH1D("beamPMC","beamPMC",40,5,7);
TH1D* beamPData=new TH1D("beamPData","beamPData",40,5,7);
TH1D* beamPDataLow=new TH1D("beamPDataLow","beamPDataLow",40,5,7);
TH1D* beamPDataHigh=new TH1D("beamPDataHigh","beamPDataHigh",40,5,7);
TH1D* beamRes=new TH1D("beamRes","beamRes",200,-5,5);
double meanShifts[12]={0,0,0,0,0,0,0,0,0,0,0,0};
double stdShifts[12]={0,0,0,0,0,0,0,0,0,0,0,0};
double sigmaShifts[12]={0,0,0,0,0,0,0,0,0,0,0,0};
tMC->Project("beamPMC","beam_inst_P",Form("selection_ID<=%d",sel));
tData->Project("beamPData","beam_inst_P",Form("selection_ID<=%d",sel));
tData->Project("beamPDataLow","beam_inst_P*0.988",Form("selection_ID<=%d",sel));
tData->Project("beamPDataHigh","beam_inst_P*1.012",Form("selection_ID<=%d",sel));
tMC->Project("beamPMC","beam_inst_P",Form("selection_ID<=%d",sel));
tData->Project("beamPData","beam_inst_P",Form("selection_ID<=%d",sel));
tData->Project("beamPDataLow","beam_inst_P*0.988",Form("selection_ID<=%d",sel));
tData->Project("beamPDataHigh","beam_inst_P*1.012",Form("selection_ID<=%d",sel));
double std=-5;
for (int i=0; i<12; ++i){
beamRes->Reset();
tMC->Project("beamRes",Form("(((1+%f*0.025)*beam_inst_P)-true_beam_startP)/true_beam_startP",std),Form("true_beam_PDG==321"));
std::cout<<std<<","<<beamRes->GetMean()<<","<<beamRes->GetStdDev()<<std::endl;
sigmaShifts[i]=std;
meanShifts[i]=beamRes->GetMean();
stdShifts[i]=beamRes->GetStdDev();
std=std+1;


}
TFile beamMap(Form("beamResMap_%d.root",sel),"recreate");
TGraph beamResMeans(12,sigmaShifts, meanShifts);
TGraph beamResStd(12, sigmaShifts, stdShifts);
beamResMeans.SetName("gMeans");
beamResStd.SetName("gWidths");
beamMap.cd();
beamResMeans.Write();
beamResStd.Write();
beamMap.Write();
beamMap.Close();


//tData->Project("beamPData","beam_inst_P",Form("abs(beam_inst_Y-reco_beam_startY)<1 && abs(beam_inst_Z-reco_beam_startZ+29)<4  && abs(beam_inst_X-reco_beam_startX+2.5)<0.5 && true_beam_PDG==321"));


//c1.SetLogy();
std::vector<double> dataBeam, mcBeam, errBeam;
for (int k=1; k<=beamPMC->GetNbinsX(); k++){
if(abs(beamPMC->GetBinCenter(k)-6)>1) continue;
dataBeam.push_back(beamPData->GetBinContent(k));
errBeam.push_back(beamPData->GetBinError(k));
mcBeam.push_back(beamPMC->GetBinContent(k));

} 
std::cout<<beamPData->Integral()<<","<<beamPMC->Integral()<<std::endl;
beamPData->Scale(1.f/beamPData->GetEntries());
beamPDataLow->Scale(1.f/beamPDataLow->GetEntries()); 
beamPDataHigh->Scale(1.f/beamPDataHigh->GetEntries());
std::cout<<beamPData->GetMaximumBin()<<std::endl;
std::cout<<beamPDataLow->GetMaximumBin()<<std::endl;
std::cout<<beamPDataHigh->GetMaximumBin()<<std::endl;

for(int i=0; i<beamPData->GetNbinsX(); i++){
double statErr=beamPData->GetBinError(i+1);
double systErr=0.5*(abs(beamPData->GetBinContent(i+1)-beamPDataHigh->GetBinContent(i+1))+abs(beamPData->GetBinContent(i+1)-beamPDataLow->GetBinContent(i+1)));
double totErr=TMath::Sqrt(statErr*statErr+systErr*systErr);
beamPData->SetBinError(i+1, totErr);


}


beamPMC->Scale(1.f/beamPMC->GetEntries());double tmp_chi2=pearsonBasic_chi2(dataBeam, mcBeam, errBeam);
beamPData->SetTitle(Form("#chi^{2}/d.o.f.=%f/%d",tmp_chi2,dataBeam.size()));
//beamPData->GetYaxis()->SetRangeUser(0.001,1);
beamPData->GetYaxis()->SetRangeUser(0,0.1);
beamPMC->SetLineColor(kRed);
beamPData->GetXaxis()->SetTitle("Momentum (GeV/c)");
beamPData->GetYaxis()->SetTitle("Arbitrary Units");
//beamPData->SetTitle("P_{beam} of Pion/Muon Candidates");
beamPData->SetTitle("P_{beam} of Candidate Kaons");
beamPData->GetXaxis()->CenterTitle();
beamPData->GetYaxis()->CenterTitle();
beamPData->Draw("e0 p");
    TLatex tL;
    tL.SetNDC();
    tL.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
beamPMC->Draw("HIST SAME ");
TLegend *leg = new TLegend(0.2,0.75,0.8,0.9);
leg->AddEntry(beamPData,Form("Data (#mu=%1.3f GeV/c, #sigma=%1.3f GeV/c)",beamPData->GetMean(),beamPData->GetStdDev()),"l");
leg->AddEntry(beamPMC,Form("Simulation (#mu=%1.3f GeV/c, #sigma=%1.3f GeV/c)",beamPMC->GetMean(),beamPMC->GetStdDev()),"l");
leg->Draw("SAME");
c1.Print(Form("beamMomGeneralFakeData_%d.png",sel));


TFile *fWrite=new TFile("plotBeamDistFakeData.root","recreate");
beamPData->Write();

beamPDataLow->Write();
beamPDataHigh->Write();

/*

int n_trklen=40; double trklen_min=0; double trklen_max=200;
TH1D* trkLen=new TH1D("trkLen","recoTrkLen",n_trklen,trklen_min,trklen_max);
TH1D* trkdEdx=new TH1D("trkdEdx","recodEdx",100,0,10);
TH1D* trkdEdxMC=new TH1D("trkdEdxMC","recodEdxMC",100,0,10);
TH1D* trkLenMC=new TH1D("trkLenMC","recoTrkLenMC",n_trklen,trklen_min,trklen_max);
TH1D* trkLenMCTrue=(TH1D*)f_rw.Get(Form("h1d_trklen_inelastic"));
TH1D* trkLenMCBkg=(TH1D*)f_rw.Get(Form("h1d_trklen_kaonbkg"));
TH1D* trkLenMCOther=(TH1D*)f_rw.Get(Form("h1d_trklen_other"));

TH1D* h1d_mpv_dEdx=(TH1D*)f_rw.Get(Form("hist_mpv_dEdx"));
TH1D* h1d_wire_dEdx=(TH1D*)f_rw.Get(Form("hist_wire_dEdx"));
TH1D* h1d_wire_pitch=(TH1D*)f_rw.Get(Form("hist_wire_pitch"));
TH1D* h1d_wire_dE=(TH1D*)f_rw.Get(Form("hist_wire_dE"));

TH1D* h1d_reco_incidentE=(TH1D*)f_rw.Get(Form("h1d_reco_incidentE"));
TH1D* h1d_reco_interactingE=(TH1D*)f_rw.Get(Form("h1d_reco_interactingE"));

TH1D* h1d_reco_cheat_incidentE=(TH1D*)f_rw.Get(Form("h1d_reco_cheat_incidentE"));
TH1D* h1d_reco_cheat_interactingE=(TH1D*)f_rw.Get(Form("h1d_reco_cheat_interactingE"));

TH1D* h1d_reco_wrongPar_incidentE=(TH1D*)f_rw.Get(Form("h1d_reco_wrongPar_incidentE"));
TH1D* h1d_reco_wrongPar_interactingE=(TH1D*)f_rw.Get(Form("h1d_reco_wrongPar_interactingE"));

TH1D* h1d_sideband_interactingE=(TH1D*)f_rw.Get(Form("h1d_reco_beyondFirstAPA_interactingE"));


TH1D* h1d_true_eff_interactingE=(TH1D*)f_rw.Get(Form("h1d_true_eff_interactingE"));
TH1D* h1d_pur_interactingE=(TH1D*)f_rw.Get("h1d_pur_interactingE");

TH1D* h1d_reco_notInel_incidentE=(TH1D*)h1d_reco_incidentE->Clone();
TH1D* h1d_reco_notInel_interactingE=(TH1D*)h1d_reco_interactingE->Clone();

h1d_reco_notInel_incidentE->Add(h1d_reco_cheat_incidentE,-1);
h1d_reco_notInel_incidentE->Add(h1d_reco_wrongPar_incidentE,-1);

h1d_reco_notInel_interactingE->Add(h1d_reco_cheat_interactingE,-1);
h1d_reco_notInel_interactingE->Add(h1d_reco_wrongPar_interactingE,-1);

h1d_reco_cheat_interactingE->SetFillColor(kRed);
h1d_reco_cheat_interactingE->SetLineColor(kRed);

h1d_reco_cheat_incidentE->SetFillColor(kRed);
h1d_reco_cheat_incidentE->SetLineColor(kRed);

h1d_reco_wrongPar_interactingE->SetFillColor(kGreen);
h1d_reco_wrongPar_interactingE->SetLineColor(kGreen);

h1d_reco_wrongPar_incidentE->SetFillColor(kGreen);
h1d_reco_wrongPar_incidentE->SetLineColor(kGreen);

h1d_reco_notInel_interactingE->SetFillColor(kBlue);
h1d_reco_notInel_interactingE->SetLineColor(kBlue);

h1d_reco_notInel_incidentE->SetFillColor(kBlue);
h1d_reco_notInel_incidentE->SetLineColor(kBlue);
std::cout<<h1d_sideband_interactingE->Integral()<<std::endl;
std::cout<<h1d_reco_incidentE->Integral()<<','<<h1d_reco_cheat_incidentE->Integral()<<','<<h1d_reco_wrongPar_incidentE->Integral()<<','<<h1d_reco_notInel_incidentE->Integral()<<std::endl;
std::cout<<h1d_reco_interactingE->Integral()<<','<<h1d_reco_cheat_interactingE->Integral()<<','<<h1d_reco_wrongPar_interactingE->Integral()<<','<<h1d_reco_notInel_interactingE->Integral()<<std::endl;
h1d_reco_interactingE->SetLineStyle(3);
  TLegend *l = new TLegend(0.2,0.7,0.6,0.9);
  l->AddEntry(h1d_reco_interactingE,"MC","l");
  l->AddEntry(h1d_reco_cheat_interactingE,"MC Kaon+Inelastic","f");
  l->AddEntry(h1d_reco_notInel_interactingE,"MC Kaon Non-Scattering Endpoint","f");
  l->AddEntry(h1d_reco_wrongPar_interactingE,"MC Non-Kaon Track","f");




THStack *hs = new THStack("hs","");
hs->Add(h1d_reco_cheat_interactingE);
hs->Add(h1d_reco_wrongPar_interactingE);
hs->Add(h1d_reco_notInel_interactingE);
h1d_reco_interactingE->GetYaxis()->SetRangeUser(0,1.5*h1d_reco_interactingE->GetMaximum());
h1d_reco_interactingE->GetXaxis()->SetTitle("Reconstructed KE [MeV]");
h1d_reco_interactingE->GetYaxis()->SetTitle("N_{int}");
h1d_reco_interactingE->GetXaxis()->CenterTitle();
h1d_reco_interactingE->GetYaxis()->CenterTitle();
h1d_reco_interactingE->SetTitle("Interacting Histogram");
h1d_reco_interactingE->Draw("HIST");
hs->Draw("HIST SAME");
TH1D* h1d_reco_clone=(TH1D*)h1d_reco_interactingE->Clone();
h1d_reco_clone->SetLineStyle(3);
h1d_reco_clone->Draw("HIST SAME");
l->Draw("SAME");
c1.Print("kaonRecoInteractingE.png");



THStack *hsIncident = new THStack("hsIncident","");
hsIncident->Add(h1d_reco_cheat_incidentE);
hsIncident->Add(h1d_reco_wrongPar_incidentE);
hsIncident->Add(h1d_reco_notInel_incidentE);
h1d_reco_incidentE->GetYaxis()->SetRangeUser(0,1.5*h1d_reco_incidentE->GetMaximum());
h1d_reco_incidentE->GetXaxis()->SetTitle("Reconstructed KE [MeV]");
h1d_reco_incidentE->GetYaxis()->SetTitle("N_{inc}");
h1d_reco_incidentE->GetXaxis()->CenterTitle();
h1d_reco_incidentE->GetYaxis()->CenterTitle();
h1d_reco_incidentE->SetTitle("Incident Histogram");
h1d_reco_incidentE->Draw("HIST");
hsIncident->Draw("HIST SAME");
h1d_reco_clone=(TH1D*)h1d_reco_incidentE->Clone();
h1d_reco_clone->SetLineStyle(3);
h1d_reco_clone->Draw("HIST SAME");
l->Draw("SAME");
c1.Print("kaonRecoIncidentE.png");


TFile fTrue("/dune/data/users/rdiurba/kaonXSTruthG4.root","read");
TTree* tTrue=(TTree*)fTrue.Get("RunInfo");
TH1D* hCrossSection= new TH1D("InelCross","Inel Cross Section",80,4000,8000);
tTrue->Project("InelCross","enerPrimGen","inelasticXS");
TH1D* h1d_reco_xsec=(TH1D*)f_data.Get(Form("h1d_reco_xsec"));
TH1D* h1d_recoMC_xsec=(TH1D*)f_rw.Get(Form("h1d_reco_xsec"));
//TH1D* h1d_recoMC_xsec=(TH1D*)f_fit.Get(Form("PostFitXSec/PostFitInelXSec"));
hCrossSection->SetTitle("Kaon Inelastic #sigma");
hCrossSection->GetYaxis()->SetTitle("#sigma [mbarn]");
hCrossSection->GetXaxis()->SetTitle("KE [MeV]");
hCrossSection->GetXaxis()->SetRangeUser(5000,5700);
hCrossSection->GetYaxis()->SetRangeUser(0,800);
hCrossSection->SetLineColor(kRed);


hCrossSection->SetLineColor(kRed);
hCrossSection->SetLineColor(kRed);

h1d_recoMC_xsec->SetLineColor(kBlack);
h1d_recoMC_xsec->SetMarkerColor(kBlack);


hCrossSection->GetXaxis()->CenterTitle();
hCrossSection->GetYaxis()->CenterTitle();
hCrossSection->Draw("HIST");
double react=0.868;
TH1D* hCrossSectionTrkLen= new TH1D("InelCrossTrkLen","Inel Cross SectionTrkLen",80,4000,8000);
tTrue->Project("InelCrossTrkLen","enerPrimGen",Form("inelasticXS*%f",react));
hCrossSectionTrkLen->SetLineColor(kBlue);
hCrossSectionTrkLen->SetMarkerColor(kBlue);
hCrossSectionTrkLen->SetLineStyle(2);
hCrossSectionTrkLen->Draw("SAME HIST");
h1d_reco_xsec->SetLineColor(kBlack);
h1d_reco_xsec->SetMarkerColor(kBlack);
h1d_reco_xsec->Draw("e0 p3 SAME");
h1d_recoMC_xsec->Draw("e0 p3 SAME");
c1.Print("kaonRecoInelXSecAll.png");



  TLegend *l2 = new TLegend(0.5,0.7,0.8,0.9);
  l2->AddEntry(hCrossSection,"Geant4 True","l");
  l2->AddEntry(hCrossSectionTrkLen,"MC-Data TrkLen Tune","l");
  l2->AddEntry(h1d_reco_xsec,"Data Raw Reco Cross Section","p");

hCrossSection->Draw("HIST");
hCrossSectionTrkLen->Draw("SAME HIST");
l2->Draw("SAME");
h1d_reco_xsec->Draw("e0 p3 SAME");
c1.Print("kaonRecoInelXSecDataOnly.png");

  TLegend *l3 = new TLegend(0.2,0.7,0.6,0.9);
  l3->AddEntry(hCrossSection,"Geant4 True","l");
  l3->AddEntry(h1d_recoMC_xsec,"MC Raw Reco Cross Section","p");


hCrossSection->Draw("HIST");
//hCrossSectionTrkLen->Draw("SAME HIST");
h1d_recoMC_xsec->Draw("e0 p3 SAME");
l3->Draw("SAME");
c1.Print("kaonRecoInelXSecMCOnly.png");



h1d_wire_dEdx->GetXaxis()->SetTitle("Wire Number");
h1d_wire_dEdx->GetYaxis()->SetTitle("dE/dx [MeV/cm]");
h1d_wire_dEdx->GetXaxis()->CenterTitle();
h1d_wire_dEdx->GetYaxis()->CenterTitle();
h1d_wire_dEdx->Draw("HIST");
c1.Print("kaon6GeVdEdxWire.png");

h1d_wire_pitch->GetXaxis()->SetTitle("Wire Number");
h1d_wire_pitch->GetYaxis()->SetTitle("Wire Pitch [cm]");
h1d_wire_pitch->GetXaxis()->CenterTitle();
h1d_wire_pitch->GetYaxis()->CenterTitle();
h1d_wire_pitch->Draw("HIST");
c1.Print("kaon6GeVTrkPitchWire.png");

h1d_wire_dE->SetTitle("Summed Energy Loss using Averages");
h1d_wire_dE->GetXaxis()->SetTitle("Reconstructed Wire Number");
h1d_wire_dE->GetYaxis()->SetTitle("Total E_{loss} [MeV]");
h1d_wire_dE->GetXaxis()->CenterTitle();
h1d_wire_dE->GetYaxis()->CenterTitle();
h1d_wire_dE->Draw("HIST");
c1.Print("kaon6GeVSumdEWire.png");

h1d_true_eff_interactingE->SetTitle("Efficiency");
h1d_true_eff_interactingE->GetXaxis()->SetTitle("True Interacting Kinetic Energy [MeV]");
h1d_true_eff_interactingE->GetYaxis()->SetTitle("Selection Efficiency");
h1d_true_eff_interactingE->GetXaxis()->CenterTitle();
h1d_true_eff_interactingE->GetYaxis()->CenterTitle();
h1d_true_eff_interactingE->GetYaxis()->SetRangeUser(0,1.0);
h1d_true_eff_interactingE->Draw("p e0");
c1.Print("kaonTrueSelEff.png");

h1d_pur_interactingE->SetTitle("Purity");
h1d_pur_interactingE->GetXaxis()->SetTitle("Reco Interacting Kinetic Energy [MeV]");
h1d_pur_interactingE->GetYaxis()->SetTitle("Selection Purity");
h1d_pur_interactingE->GetXaxis()->CenterTitle();
h1d_pur_interactingE->GetYaxis()->SetRangeUser(0,1.0);
h1d_pur_interactingE->GetYaxis()->CenterTitle();
h1d_pur_interactingE->Draw("p e0");
c1.Print("kaonRecoSelPur.png");

*/

}

