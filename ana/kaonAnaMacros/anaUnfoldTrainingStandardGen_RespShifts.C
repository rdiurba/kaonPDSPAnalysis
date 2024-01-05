#define anaUnfoldTrainingStandardGen_RespShifts_cxx
#include "anaUnfoldTrainingStandardGen_RespShifts.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"


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
#include "TRandom.h"

double beta(double gamma){
  double value=TMath::Sqrt(1-(1.0/(gamma*gamma)));
  return value;
}

double gamma(double KE,double mass){
  double value=(double(KE)/mass)+1;
  return value;
}
auto densityEffect = [](long double beta, double gamma){

   double lar_C = 5.215, lar_x0 = 0.201, lar_x1 = 3, lar_a = 0.196, lar_k = 3;
   long double x = log10(beta * gamma);
   
   if( x >= lar_x1 ) return 2*log(10)*x - lar_C;

   else if ( lar_x0 <= x && x < lar_x1) return 2*log(10)*x - lar_C + lar_a * pow(( lar_x1 - x ) , lar_k );

   else return (long double) 0; //if x < lar_x0

};

  double BetheBloch(double energy, double mass) {
   double K,rho,Z,A, charge, me, I, gamma,  /*momentum ,*/wmax, pitch;
   long double beta;
    K = 0.307;
    rho = 1.4;
    charge = 1;
    Z = 18;
    A = 39.948;
    I = pow(10,-6)*10.5*18; //MeV
    me = 0.51; //MeV me*c^2
    pitch = 1;
    
    //momentum = sqrt( pow(energy,2) - pow(massicle,2));
    //beta = momentum/sqrt(pow(massicle,2) + pow(momentum,2));
    //gamma =  1/sqrt(1 - pow(beta,2));
    
    gamma = (energy + mass) / mass;
    beta = sqrt( 1 - 1/pow(gamma,2));

    wmax = 2*me*pow(beta,2)*pow(gamma,2)/(1+2*gamma*me/mass + pow(me,2)/pow(mass,2));
    
    
    double dEdX;
    //multiply by rho to have dEdX MeV/cm in LAr

    dEdX = pitch*(rho*K*Z*pow(charge,2))/(A*pow(beta,2))*(0.5*log(2*me*pow(gamma,2)*pow(beta,2)*wmax/pow(I,2)) - pow(beta,2) - densityEffect( beta, gamma )/2 );

   return dEdX;
  };


void anaUnfoldTrainingStandardGen_RespShifts::Loop(int nUniverses, std::string weightDecision)
{
//   In a ROOT session, you can do:
//      root> .L anaUnfoldTrainingStandardGen_RespShifts.C
//      root> anaUnfoldTrainingStandardGen_RespShifts t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = 0.6667*fChain->GetEntriesFast();
double eBin = 200.;
//double minE=2000;
//double maxE=8000;
int fSliceCut=1460;
//double range=maxE-minE;
//int totBins=(range)/eBin;
int totBins=6;
totBins=4;
//totBins=3;
//double edges[7]={3000,5000,5200,5400,5600,5900,6800};
//MC binning 
//double edges[6]={4700,5000,5200,5400,5620,6020};
//Data binning
double edges[5]={4480,5080,5340,5610,6170};
// Larger center bins
//double edges[4]={4000,5150,5450,5900};
TH1D* h1d_true_thinslice_incidentE=new TH1D("h1d_true_incidentE","h1d_true_incidentE",totBins, edges);
TH1D* h1d_true_thinslice_interactingE=new TH1D("h1d_true_interactingE","h1d_true_interactingE",totBins, edges);

TH1D* h1d_reco_incidentE=new TH1D("h1d_reco_incidentE","h1d_reco_incidentE",totBins, edges);
TH1D* h1d_reco_interactingE=new TH1D("h1d_reco_interactingE","h1d_reco_interactingE",totBins, edges);

TH1D* h1d_reco_incidentE_data=new TH1D("h1d_reco_incidentE_data","h1d_reco_incidentE_data",totBins, edges);
TH1D* h1d_reco_interactingE_data=new TH1D("h1d_reco_interactingE_data","h1d_reco_interactingE_data",totBins,edges);


TH1D* h1d_reco_interactingE_passBQT=new TH1D("h1d_reco_interactingE_passBQT","h1d_reco_interactingE_passBQT",totBins, edges);
TH1D* h1d_reco_interactingE_passBQT_data=new TH1D("h1d_reco_interactingE_passBQT_data","h1d_reco_interactingE_passBQT_data",totBins, edges);
TH1D* h1d_reco_cheat_interactingE=new TH1D("h1d_reco_cheat_interactingE","h1d_reco_cheat_interactingE",totBins, edges);


std::vector<double> xSecMultisim, multisimdEdXShift, multisimBeamShift, multisimDataInc, multisimDataInt;
std::vector<double>* xsec_unfoldDataVec=0x0;
std::vector<double>* xsec_unfoldMCVec=0x0;
std::vector<TH1D*> throwHists;
std::vector<TH1D*> throwHistsMC;
std::vector<std::vector<double>> xsec_unfoldDataVecVec;
std::vector<std::vector<double>> xsec_unfoldMCVecVec;
std::vector<TH1D> xsec_unfoldDataVecHist;
std::vector<TH1D> xsec_unfoldMCVecHist;
TFile *fWrite=new TFile(Form("kaonTrainingRespShiftOutputStandardGen_%s_%d.root",weightDecision.c_str(),nUniverses),"recreate");
TTree tree("systShift","systShiftTree");
tree.Branch("xsec_unfoldDataVec",&xsec_unfoldDataVec);
tree.Branch("xsec_unfoldMCVec",&xsec_unfoldMCVec);

TProfile* h1d_MCXSec=new TProfile("h1d_MCXSec","h1d_MCXSec",totBins, edges, "s");
TProfile* h1d_DataXSec=new TProfile("h1d_DataXSec","h1d_DataXSec",totBins, edges, "s");

TProfile* h1d_MCInt=new TProfile("h1d_MCInt","h1d_MCInt",totBins, edges, "s");
TProfile* h1d_DataInt=new TProfile("h1d_DataInt","h1d_DataInt",totBins, edges, "s");

TProfile* h1d_MCInc=new TProfile("h1d_MCInc","h1d_MCInc",totBins, edges, "s");
TProfile* h1d_DataInc=new TProfile("h1d_DataInc","h1d_DataInc",totBins, edges, "s");

TFile *fCV=new TFile(Form("kaonUnfoldStandardGenRecoPlots.root"),"r");
TH2D* responseCV=(TH2D*)fCV->Get("response");
TH1D* fakesCV=(TH1D*)fCV->Get("fake");
TH1D* missesCV=(TH1D*)fCV->Get("h1d_miss_interactingE");

  TFile shift_file("beamResMap.root");
  TGraph* fSystBeamShiftMeans = (TGraph*)shift_file.Get("gMeans");
  //fSystBeamShiftMeans->SetDirectory(0);
  TGraph* fSystBeamShiftWidths = (TGraph*)shift_file.Get("gWidths");
   auto t=new TTimeStamp();
   auto time=t->AsDouble();
    std::cout<<time<<std::endl;
gRandom->SetSeed(int(time-time+5));
TH1D* sigmas_misses=new TH1D("sigmas_misses","sigmas_misses",totBins, edges);
TH1D* sigmas_fakes=new TH1D("fakes","sigmas_fakes",totBins, edges);
TH2D* sigmas=new TH2D("sigmas","sigmas",totBins, edges,totBins, edges);
for(int universe=0; universe<nUniverses; universe++){
   sigmas->Reset();
   sigmas_fakes->Reset();
   sigmas_misses->Reset();
  for(int i=0; i<totBins; i++){
  sigmas_fakes->SetBinContent(i+1,gRandom->Gaus(0,1));
  sigmas_misses->SetBinContent(i+1,gRandom->Gaus(0,1));
  for(int j=0; j<totBins; j++){
  sigmas->SetBinContent(i+1,j+1,gRandom->Gaus(0,1));

}

}




   h1d_reco_incidentE->Reset();
   h1d_reco_interactingE->Reset();
   h1d_reco_incidentE_data->Reset();
   h1d_reco_interactingE_data->Reset();
   xsec_unfoldDataVec->clear();
   xsec_unfoldMCVec->clear();

   RooUnfoldResponse response(h1d_reco_interactingE, h1d_true_thinslice_interactingE);
   RooUnfoldResponse responseIncident(h1d_reco_incidentE, h1d_true_thinslice_incidentE);
std::cout<<"Running"<<std::endl;
response.UseOverflow();
responseIncident.UseOverflow();
std::cout<<response.UseOverflowStatus()<<std::endl;
   Long64_t nbytes = 0, nb = 0;
   nentries=fChain->GetEntries()*0.66667;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {





      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
  double beam_weight=1.0;


double initialKE=beam_inst_KE;
double interactingKE=beam_inst_KE;
double interactingKERW=(beam_inst_KE*1.f);

 double tot_weight=1.f;
   if (true_beam_traj_incidentEnergies->size()){
   // h1d_true_thinslice_incidentE->Fill(5000,true_beam_traj_incidentEnergies->size());
   size_t finalSliceIndex=0;  
   int finalSlice = (true_beam_traj_slice_index)->at(true_beam_traj_incidentEnergies->size()-1); 
   for (size_t j = 0; j < true_beam_traj_incidentEnergies->size(); j++) {
        int slice = (true_beam_traj_slice_index)->at(j);
        //std::cout<<true_beam_traj_slice_z->at(j)<<std::endl;
        if (true_beam_traj_slice_z->at(j)>220.0 || true_beam_traj_slice_z->at(j)<30  )  continue;
        //h1d_true_thinslice_incidentE->Fill(true_beam_traj_incidentEnergies->at(j));
        if (selection_ID>2) responseIncident.Miss(true_beam_traj_incidentEnergies->at(j), tot_weight);
        if(selection_ID<3 && (reco_beam_calibrated_dEdX_SCE->size()-1)<j) responseIncident.Miss(true_beam_traj_incidentEnergies->at(j), tot_weight);
        finalSliceIndex=j;
      }

}
if(selection_ID<3){



if (reco_beam_calibrated_dEdX_SCE->size()>0 && reco_beam_calo_wire->size()>0 ){



double currentdEdx=BetheBloch(beam_inst_KE, 493.67);
std::vector<double> copyCalibrated=*reco_beam_calibrated_dEdX_SCE;
for(long unsigned int index=0; index<reco_beam_calibrated_dEdX_SCE->size(); ++index){

//if (true_beam_slices->size()) std::cout<<"Slices: "<<true_beam_slices->size()<<","<<reco_beam_calibrated_dEdX_SCE->size()<<std::endl;

if (reco_beam_calibrated_dEdX_SCE->at(index)>20|| reco_beam_calibrated_dEdX_SCE->at(index)<0){
interactingKE=interactingKE-currentdEdx*reco_beam_TrkPitch_SCE->at(index);
interactingKERW=interactingKERW-currentdEdx*reco_beam_TrkPitch_SCE->at(index);
//if(index<reco_beam_calibrated_dEdX_SCE->size()-1) reco_beam_calibrated_incidentEnergies->push_back(interactingKE);
}
else{
 
currentdEdx=reco_beam_calibrated_dEdX_SCE->at(index);

interactingKE=interactingKE-currentdEdx*reco_beam_TrkPitch_SCE->at(index);
interactingKERW=interactingKERW-currentdEdx*reco_beam_TrkPitch_SCE->at(index);
//if(index<reco_beam_calibrated_dEdX_SCE->size()-1) reco_beam_calibrated_incidentEnergies->push_back(interactingKE);

}
if(index<reco_beam_calibrated_dEdX_SCE->size()-1){
if (reco_beam_calo_Z->at(index)>220.0 ||  reco_beam_calo_Z->at(index)<30.0) continue;
if (reco_beam_calo_Z->at(index)<220.0) h1d_reco_incidentE->Fill(interactingKE);
if(true_beam_traj_incidentEnergies->size()){
if (index<(true_beam_slices->size()-1) && true_beam_PDG==321 && true_beam_ID==reco_beam_true_byE_ID && reco_beam_calo_Z->at(index)<220.0) responseIncident.Fill(interactingKERW,true_beam_traj_incidentEnergies->at(index), tot_weight);

}
if(!true_beam_traj_incidentEnergies->size()) responseIncident.Fake(interactingKERW, tot_weight);
if ((index>=(true_beam_traj_slice_index->size()-1) || true_beam_PDG!=321 || true_beam_ID!=reco_beam_true_byE_ID) && reco_beam_calo_Z->at(index)<220.0) responseIncident.Fake(interactingKERW, tot_weight);
}



}

}

  }
   int recoNum=0;
   int truthNum=0;
  double resp_weight=1.f;

    recoNum=h1d_reco_interactingE_data->FindBin(interactingKERW)-1;
  truthNum=h1d_reco_interactingE_data->FindBin(true_beam_traj_interactingEnergy)-1;
  
//std::cout<<recoNum<<","<<truthNum<<std::endl;
   if (selection_ID==1){ 
 // std::cout<<"GOOD"<<std::endl;
  int resp_index=recoNum*totBins+truthNum;
  resp_weight=1-sigmas->GetBinContent(recoNum+1,truthNum+1)*(1/TMath::Sqrt(responseCV->GetBinContent(recoNum+1,truthNum+1)));
  if (weightDecision=="none") resp_weight=1;
  if(recoNum>totBins || truthNum>totBins) resp_weight=1;
  if(recoNum<0 || truthNum<0) resp_weight=1;
//  std::cout<<resp_weight<<","<<","<<recoNum<<","<<truthNum<<","<<responseCV->GetBinContent(recoNum+1,truthNum+1)<<std::endl;


  h1d_reco_interactingE->Fill(reco_beam_calibrated_interactingEnergy);
  if(true_beam_endZ>30 && true_beam_endZ<220.0 &&  true_beam_endProcess->find("kaon+Inelastic")!=std::string::npos && reco_beam_true_byE_ID==true_beam_ID /*&& reco_beam_daughter_sameID==0*/){ response.Fill(interactingKERW, true_beam_traj_interactingEnergy, tot_weight);

  }
  else{
//std::cout<<"MISS"<<std::endl;  
resp_weight=1-sigmas_fakes->GetBinContent(recoNum+1)*(1/TMath::Sqrt(fakesCV->GetBinContent(recoNum+1)));
if (weightDecision=="none") resp_weight=1;
 //std::cout<<"MISS"<<std::endl; 
 if(recoNum>totBins || truthNum>totBins) resp_weight=1;
  if(recoNum<0 || truthNum<0) resp_weight=1;
    response.Fake(interactingKERW, resp_weight*tot_weight);


   } 
   }


  else{  
resp_weight=1-sigmas_misses->GetBinContent(truthNum+1)*(1/TMath::Sqrt(missesCV->GetBinContent(truthNum+1)));
if (weightDecision=="none") resp_weight=1;
  //std::cout<<"FAKE"<<std::endl;
  if(recoNum>totBins || truthNum>totBins) resp_weight=1;
  if(recoNum<0 || truthNum<0) resp_weight=1;
  if (true_beam_endZ>30 && true_beam_endZ<220.0 && true_beam_endProcess->find("kaon+Inelastic")!=std::string::npos) response.Miss(true_beam_traj_interactingEnergy, resp_weight*tot_weight);


   }



   }

std::cout<<"Unfolding MC"<<std::endl;
std::cout<<h1d_reco_interactingE->GetEntries()<<std::endl;
//TH1D* h1d_recoCopy_interactingE=(TH1D*)h1d_reco_interactingE->Clone("h1d_reco_copy_interactingE");
//TH1D* h1d_recoCopy_incidentE=(TH1D*)h1d_reco_incidentE->Clone("h1d_reco_copy_incidentE");
//RooUnfoldBayes unfoldCorr(&response, h1d_recoCopy_interactingE,4);
//unfoldCorr.IncludeSystematics();
//TH1D* hRecoFullCorr=(TH1D*) unfoldCorr.Hreco(/*RooUnfold::kCovToy*/);
//RooUnfoldBayes unfoldInc(&responseIncident, h1d_recoCopy_incidentE,4);
//unfoldInc.IncludeSystematics();
//RooUnfold::ErrorTreatment withError=3;
//TH1D* hRecoInc=(TH1D*)unfoldInc.Hreco(/*RooUnfold::kCovToy*/);
//hRecoInc->SetName("fullCorrIncHist");

//RooUnfoldBayes unfoldCheat(&response, h1d_reco_cheat_interactingE,4);
//TH1D* hRecoFullCheat=(TH1D*) unfoldCheat.Hreco();
//TMatrixD* hMatrix=(TH1D*)unfoldCorr.Hreco();
//hRecoFullCorr->SetName("fullCorrRecoHist");
TCanvas c1=TCanvas();
//h1d_recoCopy_interactingE->Draw("HIST");

TFile f("/dune/data/users/rdiurba/standard6GeVkaonTreeSel.root");
TTree* t=(TTree*)f.Get("ana");
double reco_beam_calibrated_interactingEnergy_data=0;
double beam_inst_KE_data=0;
std::vector<double>* reco_beam_calibrated_dEdX_SCE_data=0x0;
std::vector<double>* reco_beam_TrkPitch_SCE_data=0x0;
std::vector<double>* reco_beam_calo_Z_data=0x0;
std::vector<double>* reco_beam_calo_wire_data=0x0;
int selection_ID_data=0;
t->SetBranchAddress("reco_beam_calibrated_interactingEnergy",&reco_beam_calibrated_interactingEnergy_data);
t->SetBranchAddress("reco_beam_dEdX_SCE",&reco_beam_calibrated_dEdX_SCE_data);
t->SetBranchAddress("reco_beam_TrkPitch_SCE",&reco_beam_TrkPitch_SCE_data);
t->SetBranchAddress("selection_ID",&selection_ID_data);
t->SetBranchAddress("beam_inst_KE",&beam_inst_KE_data);
t->SetBranchAddress("reco_beam_calo_Z",&reco_beam_calo_Z_data);
t->SetBranchAddress("reco_beam_calo_wire",&reco_beam_calo_wire_data);
   Long64_t nentries_data=t->GetEntries();

   for (Long64_t jentry=nentries_data*0.6667+1; jentry<nentries_data;jentry++) {
   t->GetEntry(jentry);

  if (selection_ID_data==1) h1d_reco_interactingE_data->Fill(reco_beam_calibrated_interactingEnergy_data);


if(selection_ID_data<3){


if (reco_beam_calibrated_dEdX_SCE_data->size()>0 && reco_beam_calo_wire_data->size()>0 ){
double initialKE=beam_inst_KE_data;
double interactingKE=beam_inst_KE_data;

//double currentdEdx=2.3;
double currentdEdx=BetheBloch(beam_inst_KE_data, 493.67);
std::vector<double> copyCalibrated_data=*reco_beam_calibrated_dEdX_SCE_data;
for(long unsigned int index=0; index<reco_beam_calibrated_dEdX_SCE_data->size(); ++index){
if (reco_beam_calibrated_dEdX_SCE_data->at(index)>20 || reco_beam_calibrated_dEdX_SCE_data->at(index)<0){
interactingKE=interactingKE-currentdEdx*reco_beam_TrkPitch_SCE_data->at(index);
}
else{
 
currentdEdx=reco_beam_calibrated_dEdX_SCE_data->at(index);

interactingKE=interactingKE-currentdEdx*reco_beam_TrkPitch_SCE_data->at(index);
}
if(index<reco_beam_calibrated_dEdX_SCE_data->size()-1){
if (reco_beam_calo_Z_data->at(index)<220.0 && reco_beam_calo_Z_data->at(index)>30.0) h1d_reco_incidentE_data->Fill(interactingKE);
}

}
}
}
}





std::cout<<"Unfolding Data"<<std::endl;
TH1D* h1d_recoCopy_interactingE_data=(TH1D*)h1d_reco_interactingE_data->Clone("h1d_reco_copy_interactingE_data");
TH1D* h1d_recoCopy_incidentE_data=(TH1D*)h1d_reco_incidentE_data->Clone("h1d_reco_copy_incidentE_data");
auto* R=response.Hresponse();
auto* RInc=responseIncident.Hresponse();
RInc->SetName("responseIncident");
R->SetName("response");
R->SetStats(0);
RooUnfoldBayes unfoldData(&response, h1d_recoCopy_interactingE_data,4);
RooUnfoldBayes unfoldIncData(&responseIncident, h1d_recoCopy_incidentE_data,4);
TH1D* hRecoFullData=(TH1D*) unfoldData.Hreco(/*RooUnfold::kCovToy*/);
TH1D* hRecoIncData=(TH1D*) unfoldIncData.Hreco(/*RooUnfold::kCovToy*/);
/*

       TH1D * xsec_unfoldMC = (TH1D*)hRecoFullCorr->Clone("crossSectionUnfoldMC");

     xsec_unfoldMC->Divide(hRecoInc);
      for (int i = 1; i <= xsec_unfoldMC->GetNbinsX(); ++i) {
        xsec_unfoldMC->SetBinContent(i, -1.*log(1. - xsec_unfoldMC->GetBinContent(i)));


        //std::cout<<h1d_true_thinslice_incidentE_scaled->GetBinContent(i)<<','<<hRecoFullData->GetBinContent(i)<<std::endl;
        //if (hRecoFullData->GetBinContent(i)>h1d_true_thinslice_incidentE_scaled->GetBinContent(i)) xsec_data->SetBinContent(i,0);

      }
      xsec_unfoldMC->Scale(1.E27/ (0.4979 * 1.4 * 6.022E23 / 39.948 ));
*/

       TH1D * xsec_unfoldData = (TH1D*)hRecoFullData->Clone("crossSectionUnfoldData");

     xsec_unfoldData->Divide(hRecoIncData);
      for (int i = 1; i <= xsec_unfoldData->GetNbinsX(); ++i) {
        xsec_unfoldData->SetBinContent(i, -1.*log(1. - xsec_unfoldData->GetBinContent(i)));


        //std::cout<<h1d_true_thinslice_incidentE_scaled->GetBinContent(i)<<','<<hRecoFullData->GetBinContent(i)<<std::endl;
        //if (hRecoFullData->GetBinContent(i)>h1d_true_thinslice_incidentE_scaled->GetBinContent(i)) xsec_data->SetBinContent(i,0);
   
      }
      xsec_unfoldData->Scale(1.E27/ (0.4979 * 1.4 * 6.022E23 / 39.948 ));
   std::vector<TH1D*> throws;
   std::vector<TH1D*> throwsMC;
  TH1D* temp_hist=new TH1D(Form("hist_%d",int(universe)),Form("hist_%d",int(universe)),totBins, edges); 
    TH1D* temp_histMC=new TH1D(Form("histMC_%d",int(universe)),Form("histMC_%d",int(universe)),totBins, edges); 

std::vector<double> tmpVec;
std::vector<double> tmpVecMC;
      for (int i = 1; i <= xsec_unfoldData->GetNbinsX(); ++i) {
        xsec_unfoldDataVec->push_back(xsec_unfoldData->GetBinContent(i));
        tmpVec.push_back(xsec_unfoldData->GetBinContent(i));
        h1d_DataXSec->Fill(xsec_unfoldData->GetBinCenter(i),xsec_unfoldData->GetBinContent(i));
        h1d_DataInc->Fill(xsec_unfoldData->GetBinCenter(i),hRecoIncData->GetBinContent(i));
        h1d_DataInt->Fill(xsec_unfoldData->GetBinCenter(i),hRecoFullData->GetBinContent(i));
        temp_hist->Fill(xsec_unfoldData->GetBinCenter(i),xsec_unfoldData->GetBinContent(i));
        
      }
     /* for (int i = 1; i <= xsec_unfoldMC->GetNbinsX(); ++i) {
        xsec_unfoldMCVec->push_back(xsec_unfoldMC->GetBinContent(i));
tmpVecMC.push_back(xsec_unfoldData->GetBinContent(i));
        h1d_MCXSec->Fill(xsec_unfoldMC->GetBinCenter(i),xsec_unfoldMC->GetBinContent(i));
        temp_histMC->Fill(xsec_unfoldMC->GetBinCenter(i),xsec_unfoldMC->GetBinContent(i));
        h1d_MCInc->Fill(xsec_unfoldMC->GetBinCenter(i),hRecoInc->GetBinContent(i));
        h1d_MCInt->Fill(xsec_unfoldMC->GetBinCenter(i),hRecoFullCorr->GetBinContent(i));


      }*/
 throwHists.push_back(temp_hist);
 //throwHistsMC.push_back(temp_histMC);  






 

std::cout<<temp_hist->GetBinContent(4)<<std::endl;

 xSecMultisim.push_back(xsec_unfoldData->Interpolate(5000));

 multisimDataInt.push_back(hRecoFullData->Interpolate(5000));
 multisimDataInc.push_back(hRecoIncData->Interpolate(5000));

TH1D* tempHist=(TH1D*)xsec_unfoldData->Clone(Form("toy_%d",universe));
//TH1D* tempHistMC=(TH1D*)xsec_unfoldMC->Clone(Form("toyMC_%d",universe));
xsec_unfoldDataVecHist.push_back(*tempHist);
//xsec_unfoldMCVecHist.push_back(*tempHistMC);


xsec_unfoldDataVecVec.push_back(tmpVec);
xsec_unfoldMCVecVec.push_back(tmpVecMC);
fWrite->cd();
tree.Fill();


}
std::cout<<time<<std::endl;
TH2D* covMat=new TH2D("covMat","covMat",totBins, edges,totBins, edges);


TH2D* covMatMC=new TH2D("covMatMC","covMatMC",totBins, edges,totBins, edges);


for(long unsigned int toy=0; toy<xSecMultisim.size(); toy++){
    auto dataVec=xsec_unfoldDataVecVec.at(toy);
    auto mcVec=xsec_unfoldMCVecVec.at(toy);
    TH1D tmpHist=xsec_unfoldDataVecHist.at(toy);
   // TH1D tmpHistMC=xsec_unfoldMCVecHist.at(toy);

   // std::cout<<tmpHist->GetBinContent(5)<<std::endl;
  for (Int_t bin=0; bin<totBins; ++bin){
    for (Int_t bin2=0; bin2<totBins; ++bin2){
    int binNum=(bin)*totBins+(bin2);
    //double diffBin1MC=tmpHistMC.GetBinContent(bin+1)-h1d_MCXSec->GetBinContent(bin+1);
    //double diffBin2MC=tmpHistMC.GetBinContent(bin2+1)-h1d_MCXSec->GetBinContent(bin2+1);
    //if(bin>40) std::cout<<diffBin1<<std::endl;
    //covMatMC->SetBinContent(bin+1,bin2+1,covMatMC->GetBinContent(bin+1,bin2+1)+((diffBin1MC*diffBin2MC)/(xSecMultisim.size()-1)));




    double diffBin1=tmpHist.GetBinContent(bin+1)-h1d_DataXSec->GetBinContent(bin+1);
    double diffBin2=tmpHist.GetBinContent(bin2+1)-h1d_DataXSec->GetBinContent(bin2+1);
   if(bin==4 && bin2==5) std::cout<<diffBin1<<','<<diffBin2<<','<<dataVec.at(bin)<<","<<dataVec.at(bin2)<<','<<h1d_DataXSec->GetBinContent(bin+1)<<','<<h1d_DataXSec->GetBinContent(bin2+1)<<","<<bin<<','<<bin2<<std::endl;
    covMat->SetBinContent(bin+1,bin2+1,covMat->GetBinContent(bin+1,bin2+1)+((diffBin1*diffBin2)/(xSecMultisim.size()-1)));
  }
  }
}




fWrite->cd();
tree.Write();
//h1d_MCXSec->Write();
h1d_DataXSec->Write();

//h1d_MCInc->Write();
h1d_DataInc->Write();

//h1d_MCInt->Write();
h1d_DataInt->Write();

covMat->Write();
//covMatMC->Write();
h1d_reco_incidentE->Write();
h1d_reco_incidentE_data->Write();
h1d_reco_interactingE->Write();
h1d_reco_interactingE_data->Write();

fWrite->Write();
fWrite->Close();
    std::cout<<time<<std::endl;
}
