
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
void kaonSelectionBeamCutsEval( int selection_ID, bool data, int mom=6)
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

std::string fileMC="/dune/data/users/rdiurba/standard6GeVkaonTreeSel.root";
std::string fileData="/dune/app/users/rdiurba/kaonAnaMacros/kaonTreeSelData.root";
if (mom==7)
{
fileMC="/dune/app/users/rdiurba/7GeVKaonAnaMacros/7GeVkaonTreeSel.root";
fileData="/dune/app/users/rdiurba/7GeVKaonAnaMacros/7GeVkaonTreeSelData.root";

}
TFile f(Form("%s",fileData.c_str()));//"kaonTreeSelData.root");
TFile fMC(Form("%s",fileMC.c_str()));///"/dune/data/users/rdiurba/standard6GeVkaonTreeSel.root");
//TFile fMC("kaonana_mc_test.root");
//TFile fG4("/dune/data/users/rdiurba/kaonXSTruthG4.root","read")
TTree* t=(TTree*)f.Get("ana");
TTree* tMC=(TTree*)fMC.Get("ana");


int nBin=100;
double low=-4;
double high=2;






double minEnd, eff, pur;
std::vector<std::string> varArray={"true_beam_endZ-reco_beam_calo_endZ","true_beam_len-reco_beam_alt_len"/*,"reco_beam_diffXY","beam_inst_X-reco_beam_startX","beam_inst_Y-reco_beam_startY","beam_inst_Z-reco_beam_startZ","reco_beam_dCos","reco_beam_dirCos","reco_beam_endZ","reco_beam_len","reco_beam_alt_len","true_beam_endZ","beam_inst_P", "reco_beam_calibrated_dEdX_SCE","beam_inst_X","beam_inst_Y","beam_inst_XY"*/};
if (!data){
varArray.push_back("reco_beam_calibrated_interactingEnergy-true_beam_traj_interactingEnergy");
}


std::string var;

for(long unsigned int i=0; i<varArray.size(); i++){






var=varArray.at(i);







if(var.find("beam_inst_X")!=std::string::npos){
low=-45;
high=-15;


}
if(var.find("beam_inst_Y")!=std::string::npos){
low=415;
high=445;


}


if(var.find("beam_inst_XY")!=std::string::npos){
low=0;
high=20;



}



if(var.find("reco_beam_diffXY")!=std::string::npos){
low=0;
high=7;


}


if(var.find("reco_beam_calibrated_interactingEnergy-true_beam_traj_interactingEnergy")!=std::string::npos){
low=-500;
high=500;

}


if(var.find("reco_beam_calibrated_dEdX")!=std::string::npos){
low=0;
high=5;

}

if (var.find("beam_inst_P")!=std::string::npos){
low=4;
high=8;
}



if (var.find("beam_inst_X-reco_beam_startX")!=std::string::npos){
low=-5;
high=5;

if(data){
low=-6;
high=-2;
}

}
if (var.find("beam_inst_Y-reco_beam_startY")!=std::string::npos){
low=-2;
high=2;

if(data){
low=-6;
high=2;
}

}
if (var.find("beam_inst_Z-reco_beam_startZ")!=std::string::npos){
low=-40;
high=10;

if(data){
low=-35;
high=-20;
}

}

if (var.find("reco_beam_dCos")!=std::string::npos){
low=0.97;
high=1.0;

}

if (var.find("reco_beam_dirCos")!=std::string::npos){
low=0.97;
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

if(var.find("true_beam_len-reco_beam_alt_len")!=std::string::npos){
nBin=60;
low=-30;
high=30;

}
if(var.find("true_beam_endZ-reco_beam_calo_endZ")!=std::string::npos){
nBin=60;
low=-30;
high=30;

}

std::cout<<nBin<<std::endl;
double edges[nBin+1];
for (int i=0; i<nBin+1; i++){
edges[i]=float(low)+float(float(high-low)*float(float(i)/(float(nBin))));
}


std::cout<<nBin<<std::endl;

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

else t->Project("reco_endZ_data",Form("%s",var.c_str()),Form("selection_ID<=%d",selection_ID));



if(selection_ID==5){
tMC->Project("reco_cheat_endZ",Form("%s",var.c_str()),Form("true_beam_endProcess==\"%s\" && (true_beam_endZ>30 && true_beam_endZ<222.1056) && reco_beam_true_byE_ID==true_beam_ID && selection_ID<%d",inelProcess.c_str(), selection_ID));

//tMC->Project("reco_cheat_endZ_broken",Form("%s",var.c_str()),Form("true_beam_endProcess==\"%s\" && (true_beam_endZ>30 && true_beam_endZ<222.1056) && reco_beam_true_byE_ID==true_beam_ID && selection_ID<%d && reco_beam_daughter_sameID==1",inelProcess.c_str(), selection_ID));

tMC->Project("reco_cheat_outOfVol_endZ",Form("%s",var.c_str()),Form(" (true_beam_endZ<30 || true_beam_endZ>222.1056)  && reco_beam_true_byE_ID==true_beam_ID  && selection_ID<%d", selection_ID));
tMC->Project("reco_cheat_k_other_endZ",Form("%s",var.c_str()),Form(" (true_beam_endZ>30 && true_beam_endZ<222.1056) && true_beam_endProcess!=\"%s\" && reco_beam_true_byE_ID==true_beam_ID  && selection_ID<%d",inelProcess.c_str(), selection_ID));

tMC->Project("reco_cheat_sec_endZ",Form("%s",var.c_str()),Form("(reco_beam_true_byE_PDG)==321 && reco_beam_true_byE_ID!=true_beam_ID && selection_ID<%d ", selection_ID));
tMC->Project("reco_endZ_p",Form("%s",var.c_str()),Form("abs(reco_beam_true_byE_PDG)==2212 && selection_ID<%d  ", selection_ID));
tMC->Project("reco_endZ_mu",Form("%s",var.c_str()),Form("abs(reco_beam_true_byE_PDG)==13 && selection_ID<%d", selection_ID));
tMC->Project("reco_endZ_pi",Form("%s",var.c_str()),Form("abs(reco_beam_true_byE_PDG)==211 && selection_ID<%d  ", selection_ID));
tMC->Project("reco_endZ_e",Form("%s",var.c_str()),Form("(abs(reco_beam_true_byE_PDG)==11 || reco_beam_true_byE_PDG==22) && selection_ID<%d   ", selection_ID));
tMC->Project("reco_endZ_other",Form("%s",var.c_str()),Form("(reco_beam_true_byE_PDG>3000 || reco_beam_true_byE_PDG==-321)  && selection_ID<%d ", selection_ID));
tMC->Project("reco_endZ_noreco",Form("%s",var.c_str()),Form("selection_ID==%d ", selection_ID));
tMC->Project("reco_endZ",Form("%s",var.c_str()),Form("selection_ID<=%d",selection_ID));



}
else{
tMC->Project("reco_cheat_endZ",Form("%s",var.c_str()),Form("true_beam_endProcess==\"%s\" && (true_beam_endZ>30 && true_beam_endZ<222.1056) && reco_beam_true_byE_ID==true_beam_ID && selection_ID<=%d",inelProcess.c_str(), selection_ID));
//tMC->Project("reco_cheat_endZ_broken",Form("%s",var.c_str()),Form("true_beam_endProcess==\"%s\" && reco_beam_true_byE_ID==true_beam_ID && selection_ID<=%d && reco_beam_daughter_sameID==1",inelProcess.c_str(), selection_ID));
tMC->Project("reco_cheat_outOfVol_endZ",Form("%s",var.c_str()),Form(" (true_beam_endZ<30 || true_beam_endZ>222.1056) && reco_beam_true_byE_ID==true_beam_ID  && selection_ID<=%d", selection_ID));
tMC->Project("reco_cheat_k_other_endZ",Form("%s",var.c_str()),Form(" (true_beam_endZ>30 && true_beam_endZ<222.1056) && true_beam_endProcess!=\"%s\" && reco_beam_true_byE_ID==true_beam_ID  && selection_ID<=%d",inelProcess.c_str(), selection_ID));
tMC->Project("reco_cheat_sec_endZ",Form("%s",var.c_str()),Form("(reco_beam_true_byE_PDG)==321 && reco_beam_true_byE_ID!=true_beam_ID && selection_ID<=%d ", selection_ID));
tMC->Project("reco_endZ_p",Form("%s",var.c_str()),Form("abs(reco_beam_true_byE_PDG)==2212 && selection_ID<=%d  ", selection_ID));
tMC->Project("reco_endZ_mu",Form("%s",var.c_str()),Form("abs(reco_beam_true_byE_PDG)==13 && selection_ID<=%d", selection_ID));
tMC->Project("reco_endZ_pi",Form("%s",var.c_str()),Form("abs(reco_beam_true_byE_PDG)==211 && selection_ID<=%d  ", selection_ID));
tMC->Project("reco_endZ_e",Form("%s",var.c_str()),Form("(abs(reco_beam_true_byE_PDG)==11 || reco_beam_true_byE_PDG==22) && selection_ID<=%d   ", selection_ID));
tMC->Project("reco_endZ_other",Form("%s",var.c_str()),Form("(reco_beam_true_byE_PDG>3000 || reco_beam_true_byE_PDG==-321 || (reco_beam_true_byE_PDG)==-999)  && selection_ID<=%d ", selection_ID));
tMC->Project("reco_endZ_noreco",Form("%s",var.c_str()),Form("selection_ID==%d ", selection_ID));
tMC->Project("reco_endZ",Form("%s",var.c_str()),Form("selection_ID<=%d",selection_ID));
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


  TLegend *l = new TLegend(0.25,0.6,0.80,0.85);
 /* if (selection_ID==4){ h1d_reco_stack_endZ->SetTitle("No Selection Cuts"); std::cout<<"HEY"<<std::endl;}
  if (selection_ID==3) h1d_reco_stack_endZ->SetTitle("Reco Trk. Req.");
  if (selection_ID==2) h1d_reco_stack_endZ->SetTitle("Beam Cut+Reco Trk. Req.");
  if (selection_ID==1) h1d_reco_stack_endZ->SetTitle("All Cuts for Inel. Cand.");
*/
h1d_reco_stack_endZ->SetTitle(Form("%d GeV/c Sample: No Selection Cuts",mom)); 
  if (selection_ID==4) h1d_reco_stack_endZ->SetTitle(Form("%d GeV/c Sample: Has Reco. Info",mom));
  if (selection_ID==3) h1d_reco_stack_endZ->SetTitle(Form("%d GeV/c Sample: Trk. in Fid. Vol.",mom));
  if (selection_ID==2) h1d_reco_stack_endZ->SetTitle(Form("%d GeV/c Sample: K^{+} Inc.+Int. Cand.",mom));
  if (selection_ID==1) h1d_reco_stack_endZ->SetTitle(Form("%d GeV/c Sample: K^{+} Int. Cand.",mom));
  double norm=1.f;
if (!data && var.find("reco_beam_calibrated_interactingEnergy-true_beam_traj_interactingEnergy")==std::string::npos) {
  l->AddEntry(h1d_reco_sel_endZ,Form("Total Sim. Evt.: %3.0f",h1d_reco_sel_endZ->GetEntries()),"lp");
  l->AddEntry(h1d_reco_sel_cheat_endZ,Form("K^{+} inel.: %3.1f",h1d_reco_sel_cheat_endZ->GetEntries()*norm),"f");
  l->AddEntry(h1d_reco_sel_cheat_outOfVol_endZ,Form("K^{+} inel. out of Fid. Vol.: %3.1f",h1d_reco_sel_cheat_outOfVol_endZ->Integral()*norm),"f");
 // l->AddEntry(h1d_reco_sel_cheat_endZ_broken,Form("kaonBrokenTrk: %3.0f",h1d_reco_sel_cheat_endZ_broken->GetEntries()),"f");
  l->AddEntry(h1d_reco_sel_cheat_other_endZ,Form("K^{+} decay: %3.1f",h1d_reco_sel_cheat_other_endZ->GetEntries()*norm),"f");
  l->AddEntry(h1d_reco_sel_cheat_sec_endZ,Form("Sec. Beam K^{+}: %3.1f",h1d_reco_sel_cheat_sec_endZ->GetEntries()*norm),"f");
  l->AddEntry(h1d_reco_sel_endZ_p,Form("Sec. p^{+/-}: %3.1f",h1d_reco_sel_endZ_p->GetEntries()*norm),"f");
  l->AddEntry(h1d_reco_sel_endZ_pi,Form("Sec. #pi^{+/-}: %3.1f",h1d_reco_sel_endZ_pi->GetEntries()*norm),"f");
  l->AddEntry(h1d_reco_sel_endZ_e,Form("Sec. e^{+/-}/#gamma: %3.1f",h1d_reco_sel_endZ_e->GetEntries()*norm),"f");
  l->AddEntry(h1d_reco_sel_endZ_mu,Form("Sec. or cosmic #mu^{+/-}: %3.1f",h1d_reco_sel_endZ_mu->GetEntries()*norm),"f");
  l->AddEntry(h1d_reco_sel_endZ_other,Form("Sec. other: %3.1f",h1d_reco_sel_endZ_other->GetEntries()),"f");
  if (selection_ID==4)   l->AddEntry(h1d_reco_sel_endZ_noreco,Form("No Reco. Trk. or Calo.: %3.0f",h1d_reco_sel_endZ_noreco->GetEntries()),"f");

}
else if(!data){
l->AddEntry(h1d_reco_sel_endZ,Form("Total Sim. Evt.: %3.0f",h1d_reco_sel_endZ->GetEntries()),"lp");
  l->AddEntry(h1d_reco_sel_cheat_endZ,Form("K^{+} inel.: %3.1f",h1d_reco_sel_cheat_endZ->GetEntries()*norm),"f");
  l->AddEntry(h1d_reco_sel_cheat_outOfVol_endZ,Form("K^{+} inel. out of Fid. Vol.: %3.1f",h1d_reco_sel_cheat_outOfVol_endZ->GetEntries()*norm),"f");
 // l->AddEntry(h1d_reco_sel_cheat_endZ_broken,Form("kaonBrokenTrk: %3.0f",h1d_reco_sel_cheat_endZ_broken->GetEntries()),"f");
  l->AddEntry(h1d_reco_sel_cheat_other_endZ,Form("K^{+} decay: %3.1f",h1d_reco_sel_cheat_other_endZ->GetEntries()*norm),"f");
  l->AddEntry(h1d_reco_sel_cheat_sec_endZ,Form("Sec. Beam K^{+}: %3.1f",h1d_reco_sel_cheat_sec_endZ->GetEntries()*norm),"f");
  l->AddEntry(h1d_reco_sel_endZ_p,Form("Sec. p^{+/-}: %3.1f",h1d_reco_sel_endZ_p->GetEntries()*norm),"f");
  l->AddEntry(h1d_reco_sel_endZ_pi,Form("Sec. #pi^{+/-}: %3.1f",h1d_reco_sel_endZ_pi->GetEntries()*norm),"f");
  l->AddEntry(h1d_reco_sel_endZ_e,Form("Sec. e^{+/-}/#gamma: %3.1f",h1d_reco_sel_endZ_e->GetEntries()*norm),"f");
  l->AddEntry(h1d_reco_sel_endZ_mu,Form("Sec. or cosmic #mu^{+/-}: %3.1f",h1d_reco_sel_endZ_mu->GetEntries()*norm),"f");
  l->AddEntry(h1d_reco_sel_endZ_other,Form("Sec. other: %3.0f",h1d_reco_sel_endZ_other->GetEntries()),"f");



}
else    l->AddEntry(h1d_reco_sel_endZ,Form("Total Data Evt.: %3.0f",h1d_reco_sel_endZ->GetEntries()),"p");

//h1d_reco_stack_endZ->SetTitle(Form("%s with %s",var.c_str(),selection_str.c_str()));
h1d_reco_stack_endZ->SetMinimum(0);
h1d_reco_stack_endZ->SetMaximum(h1d_reco_sel_endZ->GetMaximum()*1.8);
//h1d_reco_sel_endZ->SetLineColor(0);
//h1d_reco_sel_endZ->SetLineWidth(0);
l->SetFillStyle(0);
l->SetLineWidth(0);
l->SetLineColor(kBlack);
l->SetTextFont(133);
l->SetTextSize(15);
 l->SetNColumns(2);

h1d_reco_stack_endZ->Draw("HIST");
h1d_reco_sel_endZ->Draw("SAME");
l->Draw("SAME");
if (data){ 

h1d_reco_sel_endZ->Draw("p e0");
h1d_reco_stack_endZ->Draw("p e0");
l->Draw("SAME");
}
TF1* f1=new TF1("f1","gaus",low,high);
h1d_reco_sel_endZ->Fit("f1","RS N");
if (var.find("reco_beam_endZ")!=std::string::npos){
h1d_reco_stack_endZ->GetYaxis()->SetTitle("Number of Events");
h1d_reco_stack_endZ->GetXaxis()->SetTitle("Reconstructed Z Endpoint [cm]");
h1d_reco_stack_endZ->GetYaxis()->CenterTitle();
h1d_reco_stack_endZ->GetXaxis()->CenterTitle();
}
if (var.find("true_beam_endZ")!=std::string::npos){
h1d_reco_stack_endZ->GetYaxis()->SetTitle("Number of Events");
h1d_reco_stack_endZ->GetXaxis()->SetTitle("True K^{+} Z Endpoint [cm]");
h1d_reco_stack_endZ->GetYaxis()->CenterTitle();
h1d_reco_stack_endZ->GetXaxis()->CenterTitle();
}
if(var.find("true_beam_endZ-reco_beam_calo_endZ")!=std::string::npos){
h1d_reco_stack_endZ->GetYaxis()->SetTitle("Number of Events");
h1d_reco_stack_endZ->GetXaxis()->SetTitle("Z_{end,true}-Z_{end,reco,SCE corr.} [cm]");
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
if (var.find("beam_inst_X-reco_beam_startX")!=std::string::npos){
h1d_reco_stack_endZ->GetYaxis()->SetTitle("Number of Events");
h1d_reco_stack_endZ->GetXaxis()->SetTitle("#Delta_{X,beamline-TPC} [cm]");
h1d_reco_stack_endZ->GetYaxis()->CenterTitle();
h1d_reco_stack_endZ->GetXaxis()->CenterTitle();
}
if (var.find("beam_inst_Y-reco_beam_startY")!=std::string::npos){
h1d_reco_stack_endZ->GetYaxis()->SetTitle("Number of Events");
h1d_reco_stack_endZ->GetXaxis()->SetTitle("#Delta_{Y,beamline-TPC} [cm]");
h1d_reco_stack_endZ->GetYaxis()->CenterTitle();
h1d_reco_stack_endZ->GetXaxis()->CenterTitle();
}
if (var.find("beam_inst_Z-reco_beam_startZ")!=std::string::npos){
h1d_reco_stack_endZ->GetYaxis()->SetTitle("Number of Events");
h1d_reco_stack_endZ->GetXaxis()->SetTitle("#Delta_{Z,beamline-TPC} [cm]");
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
if (var.find("reco_beam_calibrated_interactingEnergy-true_beam_traj_interactingEnergy")!=std::string::npos){
h1d_reco_stack_endZ->GetYaxis()->SetTitle("Number of Events");
h1d_reco_stack_endZ->GetXaxis()->SetTitle("#DeltaKE_{reco-true} [MeV]");
h1d_reco_stack_endZ->GetYaxis()->CenterTitle();
h1d_reco_stack_endZ->GetXaxis()->CenterTitle();
}
if (var.find("true_beam_len-reco_beam_alt_len")!=std::string::npos){
h1d_reco_stack_endZ->GetYaxis()->SetTitle("Number of Events");
h1d_reco_stack_endZ->GetXaxis()->SetTitle("Length_{true}-Length_{reco,SCE corr.} [cm]");
h1d_reco_stack_endZ->GetYaxis()->CenterTitle();
h1d_reco_stack_endZ->GetXaxis()->CenterTitle();
}


 h1d_reco_stack_endZ->GetYaxis()->SetTitleOffset(1.2);
 h1d_reco_stack_endZ->GetXaxis()->SetTitleOffset(1.2);
    TLatex tL;
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ProtoDUNE-SP}");
c1.Modified();
std::cout<<f1->GetParameter(0)<<','<<f1->GetParameter(1)<<','<<f1->GetParameter(2)<<std::endl;
c1.Print(Form("validationPlots/kaonEval%s_%d_%dStandard%dGeV.png",var.c_str(),selection_ID,data,mom));

}






}

