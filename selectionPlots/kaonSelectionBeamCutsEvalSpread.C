
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
#include "TTree.h"
#include "TLegend.h"
#include "THStack.h"
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
for(long unsigned int j=0; j<data.size(); j++){
double d=data.at(j);
double e=mc.at(j);




if (d>0 && e>0){
lngamma=lngamma+(e-d)+d*TMath::Log(d/e);
}
else lngamma=lngamma+(e-d);
}
return 2.f*lngamma;
}
void kaonSelectionBeamCutsEvalWithData( int selection_ID, bool data, int mom=6)
{

gROOT->LoadMacro("protoDUNEStyle.C");
gROOT->SetStyle("protoDUNEStyle");
gROOT->ForceStyle();
gStyle->SetTitleX(0.35);
gStyle->SetOptFit(111);
//gStyle->SetPadRightMargin(0.15);
//gStyle->SetPadLeftMargin(0.15);
TCanvas c1=TCanvas();
    c1.SetLeftMargin(0.12);
c1.SetBottomMargin(0.18);
std::string inelProcess="kaon+Inelastic";
//TFile f("/dune/data/users/rdiurba/kaonana_run5770.root");
std::string fileMC="/dune/data/users/rdiurba/standard6GeVkaonTreeSel.root";
std::string fileData="kaonTreeSelData.root";
if (mom==7)
{
fileMC="../7GeVKaonAnaMacros/7GeVkaonTreeSel.root";
fileData="../7GeVKaonAnaMacros/7GeVkaonTreeSelData.root";

}
TFile f(Form("%s",fileData.c_str()));//"kaonTreeSelData.root");
TFile fMC(Form("%s",fileMC.c_str()));///"/dune/data/users/rdiurba/standard6GeVkaonTreeSel.root");
TTree* t=(TTree*)f.Get("ana");
TTree* tMC=(TTree*)fMC.Get("ana");

int nBin=30;
double low=-4;
double high=2;






double minEnd, eff, pur;
std::vector<std::string> varArray={/*"reco_beam_calo_endZ","reco_beam_calo_startZ",*/"reco_beam_diffXY","reco_beam_startX-beam_inst_X","reco_beam_startY-beam_inst_Y","reco_beam_startZ-beam_inst_Z","beam_inst_XY","reco_beam_dCos"/*,"reco_beam_dirCos","reco_beam_endZ","reco_beam_len","reco_beam_alt_len","true_beam_endZ","beam_inst_P","reco_beam_calibrated_interactingEnergy","reco_beam_calibrated_dEdX_SCE","reco_beam_calibrated_dEdX_SCE*reco_beam_TrkPitch_SCE","beam_inst_XY"*/};
std::string var;
std::string varMC;
std::string varData;
for(long unsigned int i=0; i<varArray.size(); i++){






var=varArray.at(i);
varData=var;
varMC=var;
if (var.find("reco_beam_diffXY")!=std::string::npos){
low=-5;
high=5;
varMC="(reco_beam_diffXY-0.4308)/0.2288";
varData="(reco_beam_diffXY-4.4336)/0.4739";

}


if (var.find("reco_beam_calibrated_dEdX_SCE")!=std::string::npos){
low=0;
high=5;



}


if (var.find("reco_beam_calibrated_dEdX_SCE*reco_beam_TrkPitch_SCE")!=std::string::npos){
low=0;
high=3;



}


if (var.find("beam_inst_P")!=std::string::npos){
low=4;
high=8;



}

if(var.find("beam_inst_XY")!=std::string::npos){
low=0;
high=20;



}
if(var.find("reco_beam_calo_startZ")!=std::string::npos){
low=-4;
high=4;

}

if(var.find("reco_beam_calo_endZ")!=std::string::npos){
low=-4;
high=222;

}








if (var.find("reco_beam_startX-beam_inst_X")!=std::string::npos){
varMC="(reco_beam_startX-beam_inst_X+0.04133)/0.2566";
varData="(reco_beam_startX-beam_inst_X-3.9788)/0.2949";
low=-5;
high=5;


}
if (var.find("reco_beam_startY-beam_inst_Y")!=std::string::npos){
varMC="(reco_beam_startY-beam_inst_Y-0.3339)/0.2894";
varData="(reco_beam_startY-beam_inst_Y-1.8326)/1.0939";
low=-5;
high=5;
}


if (var.find("reco_beam_startZ-beam_inst_Z")!=std::string::npos){

varMC="(reco_beam_startZ-beam_inst_Z-29.6143)/0.4054";
varData="(reco_beam_startZ-beam_inst_Z-30.519)/1.0871";
high=5;
low=-5;

if(data){
high=35;
low=20;
}

}

if(mom==7){

if (var.find("reco_beam_diffXY")!=std::string::npos){
low=-5;
high=5;
varMC="(reco_beam_diffXY-0.39529)/0.20087";
varData="(reco_beam_diffXY-5.1616)/0.4379";

}





if (var.find("reco_beam_startX-beam_inst_X")!=std::string::npos){
varMC="(reco_beam_startX-beam_inst_X+0.04490)/0.2095";
varData="(reco_beam_startX-beam_inst_X-4.9189)/0.4086";
low=-5;
high=5;


}
if (var.find("reco_beam_startY-beam_inst_Y")!=std::string::npos){
varMC="(reco_beam_startY-beam_inst_Y-0.3300)/0.2512";
varData="(reco_beam_startY-beam_inst_Y-1.3076)/0.9548";
low=-5;
high=5;
}


if (var.find("reco_beam_startZ-beam_inst_Z")!=std::string::npos){

varMC="((reco_beam_startZ-beam_inst_Z)-29.6243)/0.4092";
varData="((reco_beam_startZ-beam_inst_Z)-31.1803)/1.2062";
high=5;
low=-5;


}

}

if (var.find("reco_beam_dCos")!=std::string::npos){
low=0.9;
high=1.0;

}

if (var.find("reco_beam_dirCos")!=std::string::npos){
low=0.9;
high=1.0;

}

if (var.find("reco_beam_endZ")!=std::string::npos){
nBin=46;
low=0.0;
high=230.0;
if(selection_ID>1){
nBin=138;
low=0.0;
high=690.0;
}

}


if (var.find("reco_beam_len")!=std::string::npos){
nBin=46;
low=0.0;
high=230.0;
if(selection_ID>1){
nBin=158;
low=0.0;
high=790.0;
}

}


if (var.find("true_beam_endZ")!=std::string::npos){
nBin=138;
low=0.0;
high=690.0;
if(selection_ID>1){
nBin=138;
low=0.0;
high=690.0;
}

}

std::cout<<nBin<<std::endl;
double edges[nBin+1];
for (int i=0; i<nBin+1; i++){
edges[i]=float(low)+float(float(high-low)*float(float(i)/(float(nBin))));
}


TH1D *h1d_reco_nocut_cheat_endZ=new TH1D("reco_nocut_cheat_endZ","reco_nocut_cheat_endZ",nBin,edges);
TH1D *h1d_reco_sel_endZ=new TH1D("reco_endZ","reco_endZ",nBin,edges);
TH1D *h1d_reco_sel_endZ_data=new TH1D("reco_endZ_data","reco_endZ_data",nBin,edges);



TH1D *h1d_reco_sel_endZ_p=new TH1D("reco_endZ_p","reco_endZ_p",nBin,edges);
TH1D *h1d_reco_sel_endZ_mu=new TH1D("reco_endZ_mu","reco_endZ_mu",nBin,edges);
TH1D *h1d_reco_sel_endZ_pi=new TH1D("reco_endZ_pi","reco_endZ_pi",nBin,edges);
TH1D *h1d_reco_sel_endZ_e=new TH1D("reco_endZ_e","reco_endZ_e",nBin,edges);
TH1D *h1d_reco_sel_endZ_other=new TH1D("reco_endZ_other","reco_endZ_other",nBin,edges);
TH1D *h1d_reco_sel_endZ_noreco=new TH1D("reco_endZ_noreco","reco_endZ_noreco",nBin,edges);
TH1D *h1d_reco_sel_cheat_endZ=new TH1D("reco_cheat_endZ","reco_cheat_endZ",nBin,edges);
TH1D *h1d_reco_sel_cheat_outOfVol_endZ=new TH1D("reco_cheat_outOfVol_endZ","reco_cheat_outOfVol_endZ",nBin,edges);
TH1D *h1d_reco_sel_cheat_endZ_broken=new TH1D("reco_cheat_endZ_broken","reco_cheat_endZ_broken",nBin,edges);
TH1D *h1d_reco_sel_cheat_other_endZ=new TH1D("reco_cheat_k_other_endZ","reco_cheat_k_other_endZ",nBin,edges);
TH1D *h1d_reco_sel_cheat_sec_endZ=new TH1D("reco_cheat_sec_endZ","reco_cheat_sec_endZ",nBin,edges);
TH1D *h1d_reco_pur_endZ=new TH1D("reco_pur_endZ","reco_pur_endZ",nBin,edges);
TH1D *h1d_reco_eff_endZ=new TH1D("reco_eff_endZ","reco_eff_endZ",nBin,edges);
TH1D *h1d_reco_pureff_endZ=new TH1D("reco_pureff_endZ","reco_pureff_endZ",23,low,115);
TH1D *h1d_reco_totEvents_endZ=new TH1D("reco_totEvents_endZ","reco_totEvents_endZ",23,low,115);

h1d_reco_sel_endZ->Reset();
h1d_reco_sel_cheat_outOfVol_endZ->Reset();
h1d_reco_sel_endZ_p->Reset();
h1d_reco_sel_endZ_mu->Reset();
h1d_reco_sel_endZ_e->Reset();
h1d_reco_sel_endZ_other->Reset();
h1d_reco_sel_endZ_noreco->Reset();
h1d_reco_sel_cheat_endZ->Reset();
h1d_reco_nocut_cheat_endZ->Reset();


if(var.find("reco_beam_calibrated_dEdX_SCE*reco_beam_TrkPitch_SCE")!=std::string::npos){
t->Project("reco_endZ_data",Form("reco_beam_dEdX_SCE*reco_beam_TrkPitch_SCE"),Form("selection_ID<=%d",selection_ID));
}
else if(var.find("reco_beam_calibrated_dEdX_SCE")!=std::string::npos){
t->Project("reco_endZ_data",Form("reco_beam_dEdX_SCE"),Form("selection_ID<=%d",selection_ID));
}

else t->Project("reco_endZ_data",Form("%s",varData.c_str()),Form("selection_ID<=%d",selection_ID));




if(selection_ID==5){

tMC->Project("reco_cheat_endZ",Form("%s",varMC.c_str()),Form("true_beam_endProcess==\"%s\" && (true_beam_endZ>30 && true_beam_endZ<222.1056) && reco_beam_true_byE_ID==true_beam_ID && selection_ID<%d",inelProcess.c_str(), selection_ID));

tMC->Project("reco_cheat_endZ_broken",Form("%s",varMC.c_str()),Form("true_beam_endProcess==\"%s\" && (true_beam_endZ>30 && true_beam_endZ<222.1056) && reco_beam_true_byE_ID==true_beam_ID && selection_ID<%d && reco_beam_daughter_sameID==1",inelProcess.c_str(), selection_ID));

tMC->Project("reco_cheat_outOfVol_endZ",Form("%s",varMC.c_str()),Form(" (true_beam_endZ<30 || true_beam_endZ>222.1056)  && reco_beam_true_byE_ID==true_beam_ID  && selection_ID<%d", selection_ID));
tMC->Project("reco_cheat_k_other_endZ",Form("%s",varMC.c_str()),Form(" (true_beam_endZ>30 && true_beam_endZ<222.1056) && true_beam_endProcess!=\"%s\" && reco_beam_true_byE_ID==true_beam_ID  && selection_ID<%d",inelProcess.c_str(), selection_ID));

tMC->Project("reco_cheat_sec_endZ",Form("%s",varMC.c_str()),Form("(reco_beam_true_byE_PDG)==321 && reco_beam_true_byE_ID!=true_beam_ID && selection_ID<%d ", selection_ID));
tMC->Project("reco_endZ_p",Form("%s",varMC.c_str()),Form("abs(reco_beam_true_byE_PDG)==2212 && selection_ID<%d  ", selection_ID));
tMC->Project("reco_endZ_mu",Form("%s",varMC.c_str()),Form("abs(reco_beam_true_byE_PDG)==13 && selection_ID<%d", selection_ID));
tMC->Project("reco_endZ_pi",Form("%s",varMC.c_str()),Form("abs(reco_beam_true_byE_PDG)==211 && selection_ID<%d  ", selection_ID));
tMC->Project("reco_endZ_e",Form("%s",varMC.c_str()),Form("(abs(reco_beam_true_byE_PDG)==11 || reco_beam_true_byE_PDG==22) && selection_ID<%d   ", selection_ID));
tMC->Project("reco_endZ_other",Form("%s",varMC.c_str()),Form("(reco_beam_true_byE_PDG>3000 || reco_beam_true_byE_PDG==-321)  && selection_ID<%d ", selection_ID));
tMC->Project("reco_endZ_noreco",Form("%s",varMC.c_str()),Form("selection_ID==%d ", selection_ID));
tMC->Project("reco_endZ",Form("%s",varMC.c_str()),Form("selection_ID<=%d",selection_ID));



}
else{
tMC->Project("reco_cheat_endZ",Form("%s",varMC.c_str()),Form("true_beam_endProcess==\"%s\" && (true_beam_endZ>30 && true_beam_endZ<222.1056) && reco_beam_true_byE_ID==true_beam_ID && selection_ID<=%d",inelProcess.c_str(), selection_ID));
tMC->Project("reco_cheat_endZ_broken",Form("%s",varMC.c_str()),Form("true_beam_endProcess==\"%s\" && reco_beam_true_byE_ID==true_beam_ID && selection_ID<=%d && reco_beam_daughter_sameID==1",inelProcess.c_str(), selection_ID));
tMC->Project("reco_cheat_outOfVol_endZ",Form("%s",varMC.c_str()),Form(" (true_beam_endZ<30 || true_beam_endZ>222.1056) && reco_beam_true_byE_ID==true_beam_ID  && selection_ID<=%d", selection_ID));
tMC->Project("reco_cheat_k_other_endZ",Form("%s",varMC.c_str()),Form(" (true_beam_endZ>30 && true_beam_endZ<222.1056) && true_beam_endProcess!=\"%s\" && reco_beam_true_byE_ID==true_beam_ID  && selection_ID<=%d",inelProcess.c_str(), selection_ID));
tMC->Project("reco_cheat_sec_endZ",Form("%s",varMC.c_str()),Form("(reco_beam_true_byE_PDG)==321 && reco_beam_true_byE_ID!=true_beam_ID && selection_ID<=%d ", selection_ID));
tMC->Project("reco_endZ_p",Form("%s",varMC.c_str()),Form("abs(reco_beam_true_byE_PDG)==2212 && selection_ID<=%d  ", selection_ID));
tMC->Project("reco_endZ_mu",Form("%s",varMC.c_str()),Form("abs(reco_beam_true_byE_PDG)==13 && selection_ID<=%d", selection_ID));
tMC->Project("reco_endZ_pi",Form("%s",varMC.c_str()),Form("abs(reco_beam_true_byE_PDG)==211 && selection_ID<=%d  ", selection_ID));
tMC->Project("reco_endZ_e",Form("%s",varMC.c_str()),Form("(abs(reco_beam_true_byE_PDG)==11 || reco_beam_true_byE_PDG==22) && selection_ID<=%d   ", selection_ID));
tMC->Project("reco_endZ_other",Form("%s",varMC.c_str()),Form("(reco_beam_true_byE_PDG>3000 || reco_beam_true_byE_PDG==-321 || (reco_beam_true_byE_PDG)==-999)  && selection_ID<=%d ", selection_ID));
tMC->Project("reco_endZ_noreco",Form("%s",varMC.c_str()),Form("selection_ID==%d ", selection_ID));
tMC->Project("reco_endZ",Form("%s",varMC.c_str()),Form("selection_ID<=%d",selection_ID));
}



//h1d_reco_sel_endZ->SetLineColor(kWhite);
//h1d_reco_sel_endZ->SetLineColor(0);

THStack* h1d_reco_stack_endZ=new THStack("hs","");
h1d_reco_sel_cheat_endZ->SetFillColor(kGreen+3);
h1d_reco_sel_cheat_outOfVol_endZ->SetFillColor(kViolet-2);
h1d_reco_sel_endZ_p->SetFillColor(kRed-10);
h1d_reco_sel_endZ_mu->SetFillColor(kYellow);
h1d_reco_sel_endZ_pi->SetFillColor(kGray);
h1d_reco_sel_endZ_e->SetFillColor(kOrange);
h1d_reco_sel_endZ_other->SetFillColor(kBlue-10);
h1d_reco_sel_cheat_other_endZ->SetFillColor(kGreen);
h1d_reco_sel_cheat_sec_endZ->SetFillColor(kCyan);


//c1.SetGrayscale();


h1d_reco_sel_cheat_endZ->SetLineColor(kGreen+3);
h1d_reco_sel_cheat_outOfVol_endZ->SetLineColor(kViolet-2);
h1d_reco_sel_endZ_p->SetLineColor(kRed-10);
h1d_reco_sel_endZ_mu->SetLineColor(kYellow);
h1d_reco_sel_endZ_pi->SetLineColor(kGray);
h1d_reco_sel_endZ_e->SetLineColor(kOrange);
h1d_reco_sel_endZ_other->SetLineColor(kBlue-10);
h1d_reco_sel_cheat_other_endZ->SetLineColor(kGreen);
h1d_reco_sel_cheat_sec_endZ->SetLineColor(kCyan);
h1d_reco_sel_endZ_noreco->SetFillStyle(4000);
h1d_reco_sel_endZ_noreco->SetLineColor(0);

h1d_reco_stack_endZ->Add(h1d_reco_sel_cheat_endZ);
h1d_reco_stack_endZ->Add(h1d_reco_sel_cheat_outOfVol_endZ);
h1d_reco_stack_endZ->Add(h1d_reco_sel_cheat_other_endZ);
h1d_reco_stack_endZ->Add(h1d_reco_sel_cheat_sec_endZ);
h1d_reco_stack_endZ->Add(h1d_reco_sel_endZ_p);
h1d_reco_stack_endZ->Add(h1d_reco_sel_endZ_mu);
h1d_reco_stack_endZ->Add(h1d_reco_sel_endZ_pi);
h1d_reco_stack_endZ->Add(h1d_reco_sel_endZ_e);
h1d_reco_stack_endZ->Add(h1d_reco_sel_endZ_other);
if (selection_ID==5) h1d_reco_stack_endZ->Add(h1d_reco_sel_endZ_noreco);


//else h1d_reco_stack_endZ->Add(h1d_reco_sel_endZ);
//h1d_reco_stack_endZ->SetLineColor(0);
//h1d_reco_stack_endZ->SetLineStyle(0);
std::cout<<h1d_reco_sel_endZ->GetEntries()<<','<<h1d_reco_sel_cheat_endZ->GetEntries()<<','<<h1d_reco_sel_cheat_other_endZ->GetEntries()<<','<<h1d_reco_sel_cheat_sec_endZ->GetEntries()<<','<<h1d_reco_sel_endZ_p->GetEntries()<<','<<h1d_reco_sel_endZ_pi->GetEntries()<<','<<h1d_reco_sel_endZ_mu->GetEntries()<<','<<h1d_reco_sel_endZ_e->GetEntries()<<','<<h1d_reco_sel_endZ_other->GetEntries()<<std::endl;


//h1d_reco_stack_endZ->GetXaxis()->SetTitle("Reco Endpoint Z [cm]");
//h1d_reco_stack_endZ->GetYaxis()->SetTitle("Total Selected");
std::string selection_str;
if (selection_ID==1){
selection_str="All Cuts";

}
else if (selection_ID==2){
selection_str="Var. Exist, Tracking, and Beam Quality Cuts";

}
else if (selection_ID==3){
selection_str="Calo Var. Exists and Beam Quality Cuts";
}
else if (selection_ID==4){
selection_str="Calo Var. Exists Cut";

}
else selection_str="No Cuts";


  TLegend *l = new TLegend(0.2,0.6,0.80,0.85);


h1d_reco_stack_endZ->SetTitle(Form("%d GeV/c Sample: No Selection Cuts",mom)); 
  if (selection_ID==4) h1d_reco_stack_endZ->SetTitle(Form("%d GeV/c Sample: Has Reco. Info",mom));
  if (selection_ID==3) h1d_reco_stack_endZ->SetTitle(Form("%d GeV/c Sample: Trk. in Fid. Vol.",mom));
  if (selection_ID==2) h1d_reco_stack_endZ->SetTitle(Form("%d GeV/c Sample: K^{+} Inc.+Int. Cand.",mom));
  if (selection_ID==1) h1d_reco_stack_endZ->SetTitle(Form("%d GeV/c Sample: K^{+} Int. Cand.",mom));



h1d_reco_sel_endZ->SetFillStyle(4000);
h1d_reco_sel_endZ->SetLineColor(0);

 h1d_reco_sel_endZ->Scale(h1d_reco_sel_endZ_data->GetEntries()/h1d_reco_sel_endZ->GetEntries());
  h1d_reco_sel_cheat_endZ->Scale(h1d_reco_sel_endZ_data->GetEntries()/h1d_reco_sel_endZ->GetEntries());
  h1d_reco_sel_cheat_outOfVol_endZ->Scale(h1d_reco_sel_endZ_data->GetEntries()/h1d_reco_sel_endZ->GetEntries());
  h1d_reco_sel_cheat_endZ_broken->Scale(h1d_reco_sel_endZ_data->GetEntries()/h1d_reco_sel_endZ->GetEntries());
  h1d_reco_sel_cheat_other_endZ->Scale(h1d_reco_sel_endZ_data->GetEntries()/h1d_reco_sel_endZ->GetEntries());
  h1d_reco_sel_cheat_sec_endZ->Scale(h1d_reco_sel_endZ_data->GetEntries()/h1d_reco_sel_endZ->GetEntries());
  h1d_reco_sel_endZ_p->Scale(h1d_reco_sel_endZ_data->GetEntries()/h1d_reco_sel_endZ->GetEntries());
  h1d_reco_sel_endZ_pi->Scale(h1d_reco_sel_endZ_data->GetEntries()/h1d_reco_sel_endZ->GetEntries());
  h1d_reco_sel_endZ_mu->Scale(h1d_reco_sel_endZ_data->GetEntries()/h1d_reco_sel_endZ->GetEntries());
  h1d_reco_sel_endZ_e->Scale(h1d_reco_sel_endZ_data->GetEntries()/h1d_reco_sel_endZ->GetEntries());
  h1d_reco_sel_endZ_other->Scale(h1d_reco_sel_endZ_data->GetEntries()/h1d_reco_sel_endZ->GetEntries());
  h1d_reco_sel_endZ_noreco->Scale(h1d_reco_sel_endZ_data->GetEntries()/h1d_reco_sel_endZ->GetEntries());



if(var.find("reco_beam_calibrated_dEdX_SCE")!=std::string::npos){

 double norm=h1d_reco_sel_endZ_data->Integral()/h1d_reco_sel_endZ->Integral();
  l->AddEntry(h1d_reco_sel_endZ_data,Form("Total Data Ent.: %3.0f",h1d_reco_sel_endZ_data->Integral()),"p");
  l->AddEntry(h1d_reco_sel_endZ,Form("Total Sim. Ent.: %3.1f",h1d_reco_sel_endZ->Integral()*norm),"l");
  l->AddEntry(h1d_reco_sel_cheat_endZ,Form("K^{+} inel.: %3.1f",h1d_reco_sel_cheat_endZ->Integral()*norm),"f");
  l->AddEntry(h1d_reco_sel_cheat_outOfVol_endZ,Form("K^{+} inel. out of Fid. Vol.: %3.1f",h1d_reco_sel_cheat_outOfVol_endZ->Integral()*norm),"f");
  //l->AddEntry(h1d_reco_sel_cheat_endZ_broken,Form("kaonBrokenTrk: %3.0f",h1d_reco_sel_cheat_endZ_broken->Integral()),"f");
  l->AddEntry(h1d_reco_sel_cheat_other_endZ,Form("K^{+} decay: %3.1f",h1d_reco_sel_cheat_other_endZ->Integral()*norm),"f");
  l->AddEntry(h1d_reco_sel_cheat_sec_endZ,Form("Sec. K^{+}: %3.1f",h1d_reco_sel_cheat_sec_endZ->Integral()*norm),"f");
  l->AddEntry(h1d_reco_sel_endZ_p,Form("Sec. p^{+/-}: %3.1f",h1d_reco_sel_endZ_p->Integral()*norm),"f");
  l->AddEntry(h1d_reco_sel_endZ_pi,Form("Sec. #pi^{+/-}: %3.1f",h1d_reco_sel_endZ_pi->Integral()*norm),"f");
  l->AddEntry(h1d_reco_sel_endZ_e,Form("Sec. e^{+/-}/#gamma: %3.1f",h1d_reco_sel_endZ_e->Integral()*norm),"f");
  l->AddEntry(h1d_reco_sel_endZ_mu,Form("Sec. or cosmic #mu^{+/-}: %3.1f",h1d_reco_sel_endZ_mu->Integral()*norm),"f");
  l->AddEntry(h1d_reco_sel_endZ_other,Form("Sec. other: %3.1f",h1d_reco_sel_endZ_other->Integral()*norm),"f");
  if (selection_ID==5)   l->AddEntry(h1d_reco_sel_endZ_noreco,Form("No Reco. Trk. or Calo.: %3.1f",h1d_reco_sel_endZ_noreco->Integral()),"f");
std::cout<<"GetEntries"<<h1d_reco_sel_endZ_data->GetIntegral()<<std::endl;
/*
 double norm=h1d_reco_sel_endZ_data->Integral(0,nBin+1)/h1d_reco_sel_endZ->Integral(0,nBin+1);
  l->AddEntry(h1d_reco_sel_endZ_data,Form("Total Data Ent.: %3.0f",h1d_reco_sel_endZ_data->Integral(0, nBin+1)),"p");
  l->AddEntry(h1d_reco_sel_endZ_data,Form("Data Overflow: %3.0f",h1d_reco_sel_endZ_data->Integral(0, nBin+1)-h1d_reco_sel_endZ_data->Integral(1,nBin)),"p");
  l->AddEntry(h1d_reco_sel_endZ,Form("Total Sim. Ent.: %3.1f",h1d_reco_sel_endZ->Integral(0,nBin+1)*norm),"l");
  l->AddEntry(h1d_reco_sel_cheat_endZ,Form("K^{+} inel.: %3.1f",h1d_reco_sel_cheat_endZ->Integral(1,nBin)*norm),"f");
 // l->AddEntry(h1d_reco_sel_cheat_endZ_broken,Form("kaonBrokenTrk: %3.0f",h1d_reco_sel_cheat_endZ_broken->Integral(1,nBin)),"f");
  l->AddEntry(h1d_reco_sel_cheat_other_endZ,Form("K^{+} decay: %3.1f",h1d_reco_sel_cheat_other_endZ->Integral(1,nBin)*norm),"f");
  l->AddEntry(h1d_reco_sel_cheat_sec_endZ,Form("Sec. Beam K^{+}: %3.1f",h1d_reco_sel_cheat_sec_endZ->Integral(1,nBin)*norm),"f");
  l->AddEntry(h1d_reco_sel_endZ_p,Form("Sec. p^{+/-}: %3.1f",h1d_reco_sel_endZ_p->Integral(1,nBin)*norm),"f");
  l->AddEntry(h1d_reco_sel_endZ_pi,Form("Sec. #pi^{+/-}: %3.1f",h1d_reco_sel_endZ_pi->Integral(1,nBin)*norm),"f");
  l->AddEntry(h1d_reco_sel_endZ_e,Form("Sec. e^{+/-}/#gamma: %3.1f",h1d_reco_sel_endZ_e->Integral(1,nBin)*norm),"f");
  l->AddEntry(h1d_reco_sel_endZ_mu,Form("Sec. or cosmic #mu^{+/-}: %3.1f",h1d_reco_sel_endZ_mu->Integral(1,nBin)*norm),"f");
  l->AddEntry(h1d_reco_sel_endZ_other,Form("Sec. other: %3.1f",h1d_reco_sel_endZ_other->Integral(1,nBin)*norm),"f");
//  l->AddEntry(h1d_reco_sel_endZ_data,Form("K^{+} inel. Overflow: %3.1f",h1d_reco_sel_cheat_endZ->GetEntries()-h1d_reco_sel_cheat_endZ->Integral(1,nBin)),"f");
  TH1D* h1d_reco_copy=(TH1D*)h1d_reco_sel_endZ_data->Clone("copy");
  h1d_reco_copy->SetLineColor(kWhite);
  h1d_reco_copy->SetFillColor(kWhite);
  h1d_reco_copy->SetMarkerColor(kWhite);
  l->AddEntry(h1d_reco_copy,Form("K^{+} inel. Overflow: %3.1f",norm*(h1d_reco_sel_cheat_endZ->Integral(0,nBin+1)-h1d_reco_sel_cheat_endZ->Integral(1,nBin))),"l");

  l->AddEntry(  h1d_reco_copy,Form("Bkg. Overflow: %3.1f",norm*(h1d_reco_sel_endZ->Integral(0,nBin+1)-h1d_reco_sel_endZ->Integral(1,nBin)-(h1d_reco_sel_cheat_endZ->Integral(0,nBin+1)-h1d_reco_sel_cheat_endZ->Integral(1,nBin)))),"l");
  if (selection_ID==4)   l->AddEntry(h1d_reco_sel_endZ_noreco,Form("No Reco. Trk.: %3.1f",h1d_reco_sel_endZ_noreco->GetEntries()),"f");
std::cout<<"GetEntries"<<h1d_reco_sel_endZ_data->GetIntegral()<<std::endl;
*/

} 

else{
   double norm=h1d_reco_sel_endZ_data->GetEntries()/h1d_reco_sel_endZ->GetEntries();

  l->AddEntry(h1d_reco_sel_endZ_data,Form("Total Data Evt.: %3.0f",h1d_reco_sel_endZ_data->GetEntries()),"p");
    if (selection_ID==5)   l->AddEntry(h1d_reco_sel_endZ,Form("Data No Reco. Trk.: %3.0f",h1d_reco_sel_endZ_data->GetBinContent(0)),"l");
  l->AddEntry(h1d_reco_sel_endZ,Form("Total Sim. Evt.: %3.1f",h1d_reco_sel_endZ->GetEntries()*norm),"l");
  l->AddEntry(h1d_reco_sel_cheat_endZ,Form("K^{+} inel.: %3.1f",h1d_reco_sel_cheat_endZ->GetEntries()*norm),"f");
  l->AddEntry(h1d_reco_sel_cheat_outOfVol_endZ,Form("K^{+} inel. out of Fid. Vol.: %3.1f",h1d_reco_sel_cheat_outOfVol_endZ->GetEntries()*norm),"f");
 // l->AddEntry(h1d_reco_sel_cheat_endZ_broken,Form("kaonBrokenTrk: %3.0f",h1d_reco_sel_cheat_endZ_broken->GetEntries()),"f");
  l->AddEntry(h1d_reco_sel_cheat_other_endZ,Form("K^{+} decay: %3.1f",h1d_reco_sel_cheat_other_endZ->GetEntries()*norm),"f");
  l->AddEntry(h1d_reco_sel_cheat_sec_endZ,Form("Sec. Beam K^{+}: %3.1f",h1d_reco_sel_cheat_sec_endZ->GetEntries()*norm),"f");
  l->AddEntry(h1d_reco_sel_endZ_p,Form("Sec. p^{+/-}: %3.1f",h1d_reco_sel_endZ_p->GetEntries()*norm),"f");
  l->AddEntry(h1d_reco_sel_endZ_pi,Form("Sec. #pi^{+/-}: %3.1f",h1d_reco_sel_endZ_pi->GetEntries()*norm),"f");
  l->AddEntry(h1d_reco_sel_endZ_e,Form("Sec. e^{+/-}/#gamma: %3.1f",h1d_reco_sel_endZ_e->GetEntries()*norm),"f");
  l->AddEntry(h1d_reco_sel_endZ_mu,Form("Sec. or cosmic #mu^{+/-}: %3.1f",h1d_reco_sel_endZ_mu->GetEntries()*norm),"f");
  l->AddEntry(h1d_reco_sel_endZ_other,Form("Sec. other: %3.1f",h1d_reco_sel_endZ_other->GetEntries()*norm),"f");
  if (selection_ID==5)   l->AddEntry(h1d_reco_sel_endZ_noreco,Form("No Reco. Trk. or Calo.: %3.1f",norm*h1d_reco_sel_endZ_noreco->GetEntries()),"f");
std::cout<<"GetEntries"<<h1d_reco_sel_endZ_data->GetIntegral()<<std::endl;

/*
 double norm=h1d_reco_sel_endZ_data->Integral(0,nBin+1)/h1d_reco_sel_endZ->Integral(0,nBin+1);
  l->AddEntry(h1d_reco_sel_endZ_data,Form("Total Data Ent.: %3.0f",h1d_reco_sel_endZ_data->Integral(0, nBin+1)),"p");
  l->AddEntry(h1d_reco_sel_endZ_data,Form("Data Overflow: %3.0f",h1d_reco_sel_endZ_data->Integral(0, nBin+1)-h1d_reco_sel_endZ_data->Integral(1,nBin)),"p");
  l->AddEntry(h1d_reco_sel_endZ,Form("Total Sim. Ent.: %3.1f",h1d_reco_sel_endZ->Integral(0,nBin+1)*norm),"l");
  l->AddEntry(h1d_reco_sel_cheat_endZ,Form("K^{+} inel.: %3.1f",h1d_reco_sel_cheat_endZ->Integral(1,nBin)*norm),"f");
 // l->AddEntry(h1d_reco_sel_cheat_endZ_broken,Form("kaonBrokenTrk: %3.0f",h1d_reco_sel_cheat_endZ_broken->Integral(1,nBin)),"f");
  l->AddEntry(h1d_reco_sel_cheat_other_endZ,Form("K^{+} decay: %3.1f",h1d_reco_sel_cheat_other_endZ->Integral(1,nBin)*norm),"f");
  l->AddEntry(h1d_reco_sel_cheat_sec_endZ,Form("Sec. Beam K^{+}: %3.1f",h1d_reco_sel_cheat_sec_endZ->Integral(1,nBin)*norm),"f");
  l->AddEntry(h1d_reco_sel_endZ_p,Form("Sec. p^{+/-}: %3.1f",h1d_reco_sel_endZ_p->Integral(1,nBin)*norm),"f");
  l->AddEntry(h1d_reco_sel_endZ_pi,Form("Sec. #pi^{+/-}: %3.1f",h1d_reco_sel_endZ_pi->Integral(1,nBin)*norm),"f");
  l->AddEntry(h1d_reco_sel_endZ_e,Form("Sec. e^{+/-}/#gamma: %3.1f",h1d_reco_sel_endZ_e->Integral(1,nBin)*norm),"f");
  l->AddEntry(h1d_reco_sel_endZ_mu,Form("Sec. or cosmic #mu^{+/-}: %3.1f",h1d_reco_sel_endZ_mu->Integral(1,nBin)*norm),"f");
  l->AddEntry(h1d_reco_sel_endZ_other,Form("Sec. other: %3.1f",h1d_reco_sel_endZ_other->Integral(1,nBin)*norm),"f");
//  l->AddEntry(h1d_reco_sel_endZ_data,Form("K^{+} inel. Overflow: %3.1f",h1d_reco_sel_cheat_endZ->GetEntries()-h1d_reco_sel_cheat_endZ->Integral(1,nBin)),"f");
  TH1D* h1d_reco_copy=(TH1D*)h1d_reco_sel_endZ_data->Clone("copy");
  h1d_reco_copy->SetLineColor(kWhite);
  h1d_reco_copy->SetFillColor(kWhite);
  h1d_reco_copy->SetMarkerColor(kWhite);
  l->AddEntry(h1d_reco_copy,Form("K^{+} inel. Overflow: %3.1f",norm*(h1d_reco_sel_cheat_endZ->Integral(0,nBin+1)-h1d_reco_sel_cheat_endZ->Integral(1,nBin))),"l");

  l->AddEntry(  h1d_reco_copy,Form("Bkg. Overflow: %3.1f",norm*(h1d_reco_sel_endZ->Integral(0,nBin+1)-h1d_reco_sel_endZ->Integral(1,nBin)-(h1d_reco_sel_cheat_endZ->Integral(0,nBin+1)-h1d_reco_sel_cheat_endZ->Integral(1,nBin)))),"l");
  if (selection_ID==4)   l->AddEntry(h1d_reco_sel_endZ_noreco,Form("No Reco. Trk.: %3.1f",h1d_reco_sel_endZ_noreco->GetEntries()),"f");
std::cout<<"GetEntries"<<h1d_reco_sel_endZ_data->GetIntegral()<<std::endl;

*/


}




//  h1d_reco_sel_endZ_data->Scale(1.f/h1d_reco_sel_endZ_data->GetEntries());

std::cout<<"GetEntries"<<h1d_reco_sel_endZ_data->Integral()<<std::endl;
std::cout<<"Get Entries "<<h1d_reco_sel_endZ->Integral()<<std::endl;
/*
  l->AddEntry(h1d_reco_sel_endZ,Form("Total: %3.0f",h1d_reco_sel_endZ->Integral()),"l");
  l->AddEntry(h1d_reco_sel_cheat_endZ,Form("kaonInel: %3.0f",h1d_reco_sel_cheat_endZ->Integral()),"f");
  l->AddEntry(h1d_reco_sel_cheat_endZ_broken,Form("kaonBrkTrk: %3.0f",h1d_reco_sel_cheat_endZ_broken->Integral()),"f");
  l->AddEntry(h1d_reco_sel_cheat_other_endZ,Form("kaonOtherEndpoint: %3.0f",h1d_reco_sel_cheat_other_endZ->Integral()),"f");
  l->AddEntry(h1d_reco_sel_cheat_sec_endZ,Form("trkIsSecKaon: %3.0f",h1d_reco_sel_cheat_sec_endZ->Integral()),"f");
  l->AddEntry(h1d_reco_sel_endZ_p,Form("p: %3.0f",h1d_reco_sel_endZ_p->Integral()),"f");
  l->AddEntry(h1d_reco_sel_endZ_pi,Form("#pi^{+/-}: %3.0f",h1d_reco_sel_endZ_pi->Integral()),"f");
  l->AddEntry(h1d_reco_sel_endZ_e,Form("e/#gamma: %3.0f",h1d_reco_sel_endZ_e->Integral()),"f");
  l->AddEntry(h1d_reco_sel_endZ_mu,Form("#mu: %3.0f",h1d_reco_sel_endZ_mu->Integral()),"f");
  l->AddEntry(h1d_reco_sel_endZ_other,Form("Other: %3.0f",h1d_reco_sel_endZ_other->Integral()),"f");

*/


   

//h1d_reco_stack_endZ->SetTitle(Form("%s with %s",var.c_str(),selection_str.c_str()));
h1d_reco_stack_endZ->SetMinimum(0);
h1d_reco_stack_endZ->SetMaximum(h1d_reco_sel_endZ->GetMaximum()*2.2);
//h1d_reco_sel_endZ->SetLineColor(0);
//h1d_reco_sel_endZ->SetLineWidth(0);
l->SetFillStyle(0);
l->SetLineWidth(0);
l->SetLineColor(kBlack);
l->SetTextFont(133);
l->SetTextSize(15);
 l->SetNColumns(2);

h1d_reco_stack_endZ->Draw("HIST");
h1d_reco_sel_endZ->Draw("SAME HIST");
h1d_reco_sel_endZ_data->Draw("SAME p e0");
l->Draw("SAME");

if(var.find("reco_beam_calibrated_dEdX_SCE")!=std::string::npos){
for(int i=0; i<h1d_reco_sel_endZ->GetNbinsX();++i){
h1d_reco_sel_endZ->SetBinError(i+1,h1d_reco_sel_endZ->GetBinError(i+1)*10);
h1d_reco_sel_endZ_data->SetBinError(i+1,h1d_reco_sel_endZ_data->GetBinError(i+1)*10);
}
}



std::cout<<"MC"<<std::endl;
TF1* f1=new TF1("f1","gaus",low,high);
TF1* f2=new TF1("f2","gaus",low,high);
h1d_reco_sel_endZ->Fit("f1","RS N");
std::cout<<"DATA"<<std::endl;
h1d_reco_sel_endZ_data->Fit("f2","RS N");
if (var.find("reco_beam_endZ")!=std::string::npos){
h1d_reco_stack_endZ->GetYaxis()->SetTitle("Number of Events");
h1d_reco_stack_endZ->GetXaxis()->SetTitle("Reconstructed Z Endpoint [cm]");
h1d_reco_stack_endZ->GetYaxis()->CenterTitle();
h1d_reco_stack_endZ->GetXaxis()->CenterTitle();
}
if (var.find("true_beam_endZ")!=std::string::npos){
h1d_reco_stack_endZ->GetYaxis()->SetTitle("Number of Events");
h1d_reco_stack_endZ->GetXaxis()->SetTitle("True Z Endpoint [cm]");
h1d_reco_stack_endZ->GetYaxis()->CenterTitle();
h1d_reco_stack_endZ->GetXaxis()->CenterTitle();
}
if (var.find("reco_beam_len")!=std::string::npos){
h1d_reco_stack_endZ->GetYaxis()->SetTitle("Number of Events");
h1d_reco_stack_endZ->GetXaxis()->SetTitle("Reconstructed Track Length [cm]");
h1d_reco_stack_endZ->GetYaxis()->CenterTitle();
h1d_reco_stack_endZ->GetXaxis()->CenterTitle();
}
if (var.find("reco_beam_alt_len")!=std::string::npos){
h1d_reco_stack_endZ->GetYaxis()->SetTitle("Number of Events");
h1d_reco_stack_endZ->GetXaxis()->SetTitle("Reconstructed Calib. Track Length [cm]");
h1d_reco_stack_endZ->GetYaxis()->CenterTitle();
h1d_reco_stack_endZ->GetXaxis()->CenterTitle();
}
if (var.find("reco_beam_startX-beam_inst_X")!=std::string::npos){
h1d_reco_stack_endZ->GetYaxis()->SetTitle("Number of Events");
h1d_reco_stack_endZ->GetXaxis()->SetTitle("#frac{#DeltaX_{TPC-beamline}-#mu}{#sigma}");
h1d_reco_stack_endZ->GetYaxis()->CenterTitle();
h1d_reco_stack_endZ->GetXaxis()->CenterTitle();
}
if (var.find("reco_beam_startY-beam_inst_Y")!=std::string::npos){
h1d_reco_stack_endZ->GetYaxis()->SetTitle("Number of Events");
h1d_reco_stack_endZ->GetXaxis()->SetTitle("#frac{#DeltaY_{TPC-beamline}-#mu}{#sigma}");
h1d_reco_stack_endZ->GetYaxis()->CenterTitle();
h1d_reco_stack_endZ->GetXaxis()->CenterTitle();
}
if (var.find("reco_beam_startZ-beam_inst_Z")!=std::string::npos){
h1d_reco_stack_endZ->GetYaxis()->SetTitle("Number of Events");
//h1d_reco_stack_endZ->GetXaxis()->SetTitle("#DeltaZ_{TPC-beamline} [cm]");
h1d_reco_stack_endZ->GetXaxis()->SetTitle("#frac{#DeltaZ_{TPC-beamline}-#mu}{#sigma}");

h1d_reco_stack_endZ->GetYaxis()->CenterTitle();
h1d_reco_stack_endZ->GetXaxis()->CenterTitle();
}

if (var.find("reco_beam_diffXY")!=std::string::npos){
h1d_reco_stack_endZ->GetYaxis()->SetTitle("Number of Events");
h1d_reco_stack_endZ->GetXaxis()->SetTitle("#frac{#DeltaR_{XY,TPC-beamline}-#mu}{#sigma}");
h1d_reco_stack_endZ->GetYaxis()->CenterTitle();
h1d_reco_stack_endZ->GetXaxis()->CenterTitle();
}


if (var.find("reco_beam_dCos")!=std::string::npos){
h1d_reco_stack_endZ->GetYaxis()->SetTitle("Number of Events");
h1d_reco_stack_endZ->GetXaxis()->SetTitle("cos(#theta_{#hat{r_{b}}*#hat{r_{T}}})");
h1d_reco_stack_endZ->GetYaxis()->CenterTitle();
h1d_reco_stack_endZ->GetXaxis()->CenterTitle();
}
if (var.find("reco_beam_dirCos")!=std::string::npos){
h1d_reco_stack_endZ->GetYaxis()->SetTitle("Number of Events");
h1d_reco_stack_endZ->GetXaxis()->SetTitle("cos(#theta_{TPC})");
h1d_reco_stack_endZ->GetYaxis()->CenterTitle();
h1d_reco_stack_endZ->GetXaxis()->CenterTitle();
}
if (var.find("reco_beam_calo_startZ")!=std::string::npos){
h1d_reco_stack_endZ->GetYaxis()->SetTitle("Number of Events");
h1d_reco_stack_endZ->GetXaxis()->SetTitle("Z_{start,SCE corr.} [cm]");
h1d_reco_stack_endZ->GetYaxis()->CenterTitle();
h1d_reco_stack_endZ->GetXaxis()->CenterTitle();
}
if (var.find("reco_beam_calo_endZ")!=std::string::npos){
h1d_reco_stack_endZ->GetYaxis()->SetTitle("Number of Events");
h1d_reco_stack_endZ->GetXaxis()->SetTitle("Z_{end,SCE corr.} [cm]");
h1d_reco_stack_endZ->GetYaxis()->CenterTitle();
h1d_reco_stack_endZ->GetXaxis()->CenterTitle();
}
if (var.find("reco_beam_calibrated_dEdX_SCE")!=std::string::npos){
h1d_reco_stack_endZ->GetYaxis()->SetTitle("Number of Collection Plane Slices");
h1d_reco_stack_endZ->GetXaxis()->SetTitle("dE/dx [MeV/c]");
h1d_reco_stack_endZ->GetYaxis()->CenterTitle();
h1d_reco_stack_endZ->GetXaxis()->CenterTitle();
}
if (var.find("reco_beam_calibrated_dEdX_SCE*reco_beam_TrkPitch_SCE")!=std::string::npos){
h1d_reco_stack_endZ->GetYaxis()->SetTitle("Number of Events");
h1d_reco_stack_endZ->GetXaxis()->SetTitle("dE");
h1d_reco_stack_endZ->GetYaxis()->CenterTitle();
h1d_reco_stack_endZ->GetXaxis()->CenterTitle();
h1d_reco_stack_endZ->GetXaxis()->SetRangeUser(0.0,3.0);
f1->Draw("lsame");
f2->Draw("lsame");
}


 h1d_reco_stack_endZ->GetYaxis()->SetTitleOffset(1.2);
 h1d_reco_stack_endZ->GetXaxis()->SetTitleOffset(1.2);
    TLatex tL;
    tL.SetNDC();
    tL.DrawLatex(0.2,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c1.Modified();




std::cout<<f1->GetParameter(0)<<','<<f1->GetParameter(1)<<','<<f1->GetParameter(2)<<std::endl;
std::cout<<f2->GetParameter(0)<<','<<f2->GetParameter(1)<<','<<f2->GetParameter(2)<<std::endl;
std::cout<<f1->GetChisquare()<<','<<f1->GetNDF()<<std::endl;
std::cout<<f2->GetChisquare()<<','<<f2->GetNDF()<<std::endl;
c1.Print(Form("validationPlots/kaonEval%s_%d_%d_wDataStandardSpread%dGeV.png",var.c_str(),selection_ID,data,mom));
c1.Print(Form("validationPlots/kaonEval%s_%d_%d_wDataStandardSpread%dGeV.pdf",var.c_str(),selection_ID,data,mom));

}







}

