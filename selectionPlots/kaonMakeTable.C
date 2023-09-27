
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




if (d>0 & e>0){
lngamma=lngamma+(e-d)+d*TMath::Log(d/e);
}
else lngamma=lngamma+(e-d);
}
return 2.f*lngamma;
}
void kaonMakeTable()
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
TFile fMC("/dune/data/users/rdiurba/standard6GeVkaonTreeSel.root");
TFile f("kaonTreeSelData.root");

TFile fMC7("../7GeVKaonAnaMacros/7GeVkaonTreeSel.root");
TFile f7("../7GeVKaonAnaMacros/7GeVkaonTreeSelData.root");;
//TFile fMC("kaonana_mc_test.root");
//TFile fG4("/dune/data/users/rdiurba/kaonXSTruthG4.root","read")
TTree* tData=(TTree*)f.Get("ana");
TTree* tMC=(TTree*)fMC.Get("ana");
TTree* tData7=(TTree*)f7.Get("ana");
TTree* tMC7=(TTree*)fMC7.Get("ana");







double allTrueKaonsInVol=tMC->GetEntries(Form("true_beam_endZ>30 && true_beam_endZ<222.1056 && true_beam_endProcess==\"%s\"",inelProcess.c_str()));
double all7GeVTrueKaonsInVol=tMC7->GetEntries(Form("true_beam_endZ>30 && true_beam_endZ<222.1056 && true_beam_endProcess==\"%s\"",inelProcess.c_str()));
std::vector<int> goodKaons;
std::vector<int>  outOfVolKaons;
std::vector<int>  decayKaons;
std::vector<int> kaonsRemaining;
std::vector<int>  secKaons;
std::vector<int>  secP;
std::vector<int>  secMu;
std::vector<int>  secEM;
std::vector<int>  secPi;
std::vector<int>  secOther;
std::vector<int>  all;



std::cout<<"6 GeV flow table"<<std::endl;
for(int selection_ID=1; selection_ID<6; selection_ID++){
all.push_back(tMC->GetEntries(Form("selection_ID<=%d",selection_ID)));
kaonsRemaining.push_back(tMC->GetEntries(Form("true_beam_endZ>30 && true_beam_endZ<222.1056 && true_beam_endProcess==\"%s\" &&  selection_ID<=%d",inelProcess.c_str(), selection_ID)));
goodKaons.push_back(tMC->GetEntries(Form("true_beam_endZ>30 && true_beam_endZ<222.1056 && true_beam_endProcess==\"%s\" && reco_beam_true_byE_ID==true_beam_ID && selection_ID<=%d",inelProcess.c_str(), selection_ID)));
outOfVolKaons.push_back(tMC->GetEntries(Form("(true_beam_endZ<30 || true_beam_endZ>222.1056) && reco_beam_true_byE_ID==true_beam_ID && selection_ID<=%d", selection_ID)));
decayKaons.push_back(tMC->GetEntries(Form("true_beam_endZ>30 && true_beam_endZ<222.1056 && true_beam_endProcess!=\"%s\" && reco_beam_true_byE_ID==true_beam_ID && selection_ID<=%d ",inelProcess.c_str(), selection_ID)));
secKaons.push_back(tMC->GetEntries(Form("(reco_beam_true_byE_PDG)==321 && reco_beam_true_byE_ID!=true_beam_ID && selection_ID<=%d ", selection_ID)));
secP.push_back(tMC->GetEntries(Form("abs(reco_beam_true_byE_PDG)==2212 && selection_ID<=%d  ", selection_ID)));
secMu.push_back(tMC->GetEntries(Form("abs(reco_beam_true_byE_PDG)==13 && selection_ID<=%d  ", selection_ID)));
secEM.push_back(tMC->GetEntries(Form("(abs(reco_beam_true_byE_PDG)==11 || reco_beam_true_byE_PDG==22) && selection_ID<=%d  ", selection_ID)));

secPi.push_back(tMC->GetEntries(Form("abs(reco_beam_true_byE_PDG)==211 && selection_ID<=%d  ", selection_ID)));
secOther.push_back(tMC->GetEntries(Form("(reco_beam_true_byE_PDG>3000 || reco_beam_true_byE_PDG==-321 || (reco_beam_true_byE_PDG)==-999)  && selection_ID<=%d ", selection_ID)));
}





std::vector<int>  goodKaons7;
std::vector<int>  kaonsRemaining7;
std::vector<int>  outOfVolKaons7;
std::vector<int>  decayKaons7;
std::vector<int>  secKaons7;
std::vector<int>  secP7;
std::vector<int>  secMu7;
std::vector<int>  secEM7;
std::vector<int>  secPi7;
std::vector<int>  secOther7;
std::vector<int>  all7;

std::cout<<"7 GeV flow table"<<std::endl;
for(int selection_ID=1; selection_ID<6; selection_ID++){
all7.push_back(tMC7->GetEntries(Form("selection_ID<=%d",selection_ID)));
goodKaons7.push_back(tMC7->GetEntries(Form("true_beam_endZ>30 && true_beam_endZ<222.1056 && true_beam_endProcess==\"%s\" && reco_beam_true_byE_ID==true_beam_ID && selection_ID<=%d",inelProcess.c_str(), selection_ID)));
kaonsRemaining7.push_back(tMC7->GetEntries(Form("true_beam_endZ>30 && true_beam_endZ<222.1056 && true_beam_endProcess==\"%s\" &&  selection_ID<=%d",inelProcess.c_str(), selection_ID)));
outOfVolKaons7.push_back(tMC7->GetEntries(Form("(true_beam_endZ<30 || true_beam_endZ>222.1056) && reco_beam_true_byE_ID==true_beam_ID && selection_ID<=%d", selection_ID)));
decayKaons7.push_back(tMC7->GetEntries(Form("true_beam_endZ>30 && true_beam_endZ<222.1056 && true_beam_endProcess!=\"%s\" && reco_beam_true_byE_ID==true_beam_ID && selection_ID<=%d",inelProcess.c_str(), selection_ID)));
secKaons7.push_back(tMC7->GetEntries(Form("(reco_beam_true_byE_PDG)==321 && reco_beam_true_byE_ID!=true_beam_ID && selection_ID<=%d ", selection_ID)));
secP7.push_back(tMC7->GetEntries(Form("abs(reco_beam_true_byE_PDG)==2212 && selection_ID<=%d  ", selection_ID)));
secMu7.push_back(tMC7->GetEntries(Form("abs(reco_beam_true_byE_PDG)==13 && selection_ID<=%d  ", selection_ID)));
secEM7.push_back(tMC7->GetEntries(Form("(abs(reco_beam_true_byE_PDG)==11 || reco_beam_true_byE_PDG==22) && selection_ID<=%d  ", selection_ID)));

secPi7.push_back(tMC7->GetEntries(Form("abs(reco_beam_true_byE_PDG)==211 && selection_ID<=%d  ", selection_ID)));
secOther7.push_back(tMC7->GetEntries(Form("(reco_beam_true_byE_PDG>3000 || reco_beam_true_byE_PDG==-321 || (reco_beam_true_byE_PDG)==-999)  && selection_ID<=%d ", selection_ID)));
}
std::cout<<"General Statistics"<<std::endl;
std::cout<<"6 GeV/$c$ sim."<<" & "<< all.at(4) <<" & "<<all.at(3) <<" & " <<all.at(2) <<" & " <<all.at(1) <<" & "<<all.at(0)<<std::endl;

std::cout<<"6 GeV/$c$ data"<<" & "<<tData->GetEntries("selection_ID<6")<<" & "<<tData->GetEntries("selection_ID<5")<<" & "<<tData->GetEntries("selection_ID<4")<<" & "<<tData->GetEntries("selection_ID<3")<<" & "<<tData->GetEntries("selection_ID<2")<<std::endl;


std::cout<<"7 GeV/$c$ sim."<<" & "<<all7.at(4) <<" & "<<all7.at(3) <<" & " <<all7.at(2) <<" & " <<all7.at(1) <<" & "<<all7.at(0)<<std::endl;

std::cout<<"7 GeV/$c$ data"<<" & "<<tData7->GetEntries("selection_ID<6")<<" & "<<tData7->GetEntries("selection_ID<5")<<" & "<<tData7->GetEntries("selection_ID<4")<<" & "<<tData7->GetEntries("selection_ID<3")<<" & "<<tData7->GetEntries("selection_ID<2")<<std::endl;



std::cout<<"Purity and Efficiency"<<std::endl;

std::cout<<"6 GeV/$c$ sample event efficiency"<<" & "<<100.0*float(kaonsRemaining.at(4))/allTrueKaonsInVol<<"\% & "<<100.0*float(kaonsRemaining.at(3))/allTrueKaonsInVol<<"\% & "<<100.0*float(kaonsRemaining.at(2))/allTrueKaonsInVol<<"\% & "<<100.0*float(kaonsRemaining.at(1))/allTrueKaonsInVol<<"\% & "<<100.0*float(kaonsRemaining.at(0))/allTrueKaonsInVol<<"\%"<<std::endl,

std::cout<<"6 GeV/$c$ sample event purity"<<" & "<<100.0*float(goodKaons.at(4))/float(all.at(4))<<"\% & "<<
100.0*float(goodKaons.at(3))/float(all.at(3))<<"\% & "<<
100.0*float(goodKaons.at(2))/float(all.at(2))<<"\% & "<<
100.0*float(goodKaons.at(1))/float(all.at(1))<<"\% & "<<
100.0*float(goodKaons.at(0))/float(all.at(0))<<"\% "<<std::endl;

std::cout<<"7 GeV/$c$ sample event efficiency"<<" & "<<100.0*float(kaonsRemaining7.at(4))/all7GeVTrueKaonsInVol<<"\% & "<<100.0*float(kaonsRemaining7.at(3))/all7GeVTrueKaonsInVol<<"\% & "<<100.0*float(kaonsRemaining7.at(2))/all7GeVTrueKaonsInVol<<"\% & "<<100.0*float(kaonsRemaining7.at(1))/all7GeVTrueKaonsInVol<<"\% & "<<100.0*float(kaonsRemaining7.at(0))/all7GeVTrueKaonsInVol<<"\%"<<std::endl,

std::cout<<"7 GeV/$c$ sample event purity "<<" & "<<100.0*float(goodKaons7.at(4))/float(all7.at(4))<<"\% & "<<
100.0*float(goodKaons7.at(3))/float(all7.at(3))<<"\% & "<<
100.0*float(goodKaons7.at(2))/float(all7.at(2))<<"\% & "<<
100.0*float(goodKaons7.at(1))/float(all7.at(1))<<"\% & "<<
100.0*float(goodKaons7.at(0))/float(all7.at(0))<<"\% "<<std::endl;

std::cout<<"Flow chart for 6 GeV"<<std::endl;
for(int selection_ID=1; selection_ID<5; selection_ID++){


std::cout<<selection_ID<<" & "<<all.at(4-selection_ID)<<" & "<<goodKaons.at(4-selection_ID)<<" & "<<outOfVolKaons.at(4-selection_ID)<<" & "<<decayKaons.at(4-selection_ID)<<" & "<<secKaons.at(4-selection_ID)<<" & "<<secP.at(4-selection_ID)<<" & "<<secPi.at(4-selection_ID)<<" & "<<secEM.at(4-selection_ID)<<" & "<<secMu.at(4-selection_ID)<<" & "<<secOther.at(4-selection_ID)<<std::endl;
}

std::cout<<"Flow chart for 7 GeV"<<std::endl;
for(int selection_ID=1; selection_ID<5; selection_ID++){


std::cout<<selection_ID<<" & "<<all7.at(4-selection_ID)<<" & "<<goodKaons7.at(4-selection_ID)<<" & "<<outOfVolKaons7.at(4-selection_ID)<<" & "<<decayKaons7.at(4-selection_ID)<<" & "<<secKaons7.at(4-selection_ID)<<" & "<<secP7.at(4-selection_ID)<<" & "<<secPi7.at(4-selection_ID)<<" & "<<secEM7.at(4-selection_ID)<<" & "<<secMu7.at(4-selection_ID)<<" & "<<secOther7.at(4-selection_ID)<<std::endl;
}





}

