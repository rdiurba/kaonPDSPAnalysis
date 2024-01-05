#define anaUnfoldFakeDataNewFidVol_cxx
#include "anaUnfoldFakeDataNewFidVol.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "TRandom.h"
#include "TDecompChol.h"

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



void anaUnfoldFakeDataNewFidVol::Loop(double  vol, double vol2)
{
//   In a ROOT session, you can do:
//      root> .L anaUnfoldFakeDataNewFidVol.C
//      root> anaUnfoldFakeDataNewFidVol t
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

   Long64_t nentries = fChain->GetEntriesFast();
double g4rw_index=1.f;
double eBin = 200.;
//double minE=2000;
//double maxE=8000;
int fSliceCut=1460;
//double range=maxE-minE;
//int totBins=(range)/eBin;
int totBins=4;
//double edges[7]={3000,5000,5200,5400,5600,5900,6800};
//double edges[6]={4000,5000,5240,5440,5700,6650};
//double edges[6]={4480,4980,5220,5410,5630,6200};
double edges[5]={4480,5080,5340,5610,6170};

TFile fTruth("/dune/data/users/rdiurba/rootDump/kaon_cross_section_out.root");
TGraph* g4Truth=(TGraph*)fTruth.Get("inel_KE");
//RooUnfoldResponse::UseOverflow();

TH1D* h1d_true_thinslice_incidentE=new TH1D("h1d_true_incidentE","h1d_true_incidentE",totBins, edges);
TH1D* h1d_true_thinslice_interactingE=new TH1D("h1d_true_interactingE","h1d_true_interactingE",totBins, edges);

TH1D* h1d_reco_incidentE=new TH1D("h1d_reco_incidentE","h1d_reco_incidentE",totBins, edges);
TH1D* h1d_reco_interactingE=new TH1D("h1d_reco_interactingE","h1d_reco_interactingE",totBins, edges);

TH1D* h1d_reco_interactingE_passBQT=new TH1D("h1d_reco_interactingE_passBQT","h1d_reco_interactingE_passBQT",totBins, edges);
TH1D* h1d_reco_interactingE_passBQT_data=new TH1D("h1d_reco_interactingE_passBQT_data","h1d_reco_interactingE_passBQT_data",totBins, edges);
TH1D* h1d_reco_cheat_interactingE=new TH1D("h1d_reco_cheat_interactingE","h1d_reco_cheat_interactingE",totBins, edges);
RooUnfoldResponse response(h1d_reco_interactingE, h1d_true_thinslice_interactingE);
RooUnfoldResponse responseIncident(h1d_reco_incidentE, h1d_true_thinslice_incidentE);

response.UseOverflow();
responseIncident.UseOverflow();
std::cout<<response.UseOverflowStatus()<<std::endl;
   Long64_t nbytes = 0, nb = 0;
   nentries=nentries*0.66667;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
           Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   double ediv_var=0.4287;

   double nobeam_var=0.953;

  /* Get Ediv */
  double ediv_weight=1;
  if (selection_ID==2){
  if (reco_beam_endZ > 220.0 && reco_beam_endZ<234.0) {
    ediv_weight=ediv_var;
  }
  else {
    ediv_weight = (1-ediv_var*0.55)/(1-0.55);
  }

  
  }

double nobeam_weight=1;
 if(true_beam_traj_incidentEnergies->size()>0 && true_beam_endZ<30){
  if(selection_ID==5) nobeam_weight=nobeam_var;
  else nobeam_weight=(1-0.45*nobeam_var)/(1-0.45);

}
double  weight=1.f;//nobeam_weight*ediv_weight;  

if(true_beam_endZ>30 && true_beam_endZ<220.0 && true_beam_endProcess->find("kaon+Inelastic")!=std::string::npos)  h1d_true_thinslice_interactingE->Fill(true_beam_traj_interactingEnergy);
   
   if (true_beam_traj_incidentEnergies->size()){
   // h1d_true_thinslice_incidentE->Fill(5000,true_beam_traj_incidentEnergies->size());
   size_t finalSliceIndex=0;  
   int finalSlice = (true_beam_traj_slice_index)->at(true_beam_traj_incidentEnergies->size()-1); 
   for (size_t j = 0; j < true_beam_traj_incidentEnergies->size(); j++) {
        //int slice = (true_beam_traj_slice_index)->at(j);
        //if (slice>fSliceCut) continue; 
        //std::cout<<true_beam_traj_slice_z->at(j)<<std::endl;
        if (true_beam_traj_slice_z->at(j)<30 || true_beam_traj_slice_z->at(j)>220.0 )  continue;
        h1d_true_thinslice_incidentE->Fill(true_beam_traj_incidentEnergies->at(j));
        if (selection_ID>2){ responseIncident.Miss(true_beam_traj_incidentEnergies->at(j),weight); 
}
        if(selection_ID<3 && (reco_beam_calibrated_dEdX_SCE->size()-1)<j){ responseIncident.Miss(true_beam_traj_incidentEnergies->at(j),weight);}
        finalSliceIndex=j;
      }

}
if(selection_ID<3){

h1d_reco_interactingE_passBQT->Fill(reco_beam_calibrated_interactingEnergy);

if (reco_beam_calibrated_dEdX_SCE->size()>0 && reco_beam_calo_wire->size()>0 ){
double initialKE=beam_inst_KE;
double interactingKE=beam_inst_KE;

//double currentdEdx=2.3;
double currentdEdx=BetheBloch(beam_inst_KE, 493.67);
std::vector<double> copyCalibrated=*reco_beam_calibrated_dEdX_SCE;
for(long unsigned int index=0; index<reco_beam_calibrated_dEdX_SCE->size(); ++index){

//if (true_beam_slices->size()) std::cout<<"Slices: "<<true_beam_slices->size()<<","<<reco_beam_calibrated_dEdX_SCE->size()<<std::endl;

if (reco_beam_calibrated_dEdX_SCE->at(index)>20 || reco_beam_calibrated_dEdX_SCE->at(index)<0){
interactingKE=interactingKE-currentdEdx*reco_beam_TrkPitch_SCE->at(index);
//if(index<reco_beam_calibrated_dEdX_SCE->size()-1) reco_beam_calibrated_incidentEnergies->push_back(interactingKE);
}
else{
 
currentdEdx=reco_beam_calibrated_dEdX_SCE->at(index);

interactingKE=interactingKE-currentdEdx*reco_beam_TrkPitch_SCE->at(index);

//if(index<reco_beam_calibrated_dEdX_SCE->size()-1) reco_beam_calibrated_incidentEnergies->push_back(interactingKE);

}
if(index<reco_beam_calibrated_dEdX_SCE->size()-1){
if (reco_beam_calo_Z->at(index)<30.0 || reco_beam_calo_Z->at(index)>220.0) continue;
if (reco_beam_calo_Z->at(index)>30.0 && reco_beam_calo_Z->at(index)<220.0) h1d_reco_incidentE->Fill(interactingKE);
if (true_beam_traj_incidentEnergies->size()){
if (index<(true_beam_slices->size()-1) && true_beam_PDG==321 && true_beam_ID==reco_beam_true_byE_ID && reco_beam_calo_Z->at(index)<220.0 && reco_beam_calo_Z->at(index)>30.0) responseIncident.Fill(interactingKE,true_beam_traj_incidentEnergies->at(index),weight);

}
if(!true_beam_traj_incidentEnergies->size()) responseIncident.Fake(interactingKE,weight);
if ((index>=(true_beam_traj_slice_index->size()-1) || true_beam_PDG!=321 || true_beam_ID!=reco_beam_true_byE_ID) && reco_beam_calo_Z->at(index)>30.0 && reco_beam_calo_Z->at(index)<220.0) responseIncident.Fake(interactingKE,weight);
}



}

}

  }
  if (selection_ID==1){ 


  h1d_reco_interactingE->Fill(reco_beam_calibrated_interactingEnergy);
  if(true_beam_endZ>30 && true_beam_endZ<220.0 && true_beam_endProcess->find("kaon+Inelastic")!=std::string::npos && reco_beam_true_byE_ID==true_beam_ID /* && reco_beam_daughter_sameID==0*/){ response.Fill(reco_beam_calibrated_interactingEnergy, true_beam_traj_interactingEnergy,weight);

  h1d_reco_cheat_interactingE->Fill(reco_beam_calibrated_interactingEnergy);
  }
  else{
    response.Fake(reco_beam_calibrated_interactingEnergy,weight);


   } 
   }


  else{
  if (true_beam_endZ>30  && true_beam_endZ<220.0 && true_beam_endProcess->find("kaon+Inelastic")!=std::string::npos){ response.Miss(true_beam_traj_interactingEnergy,weight);

}

   }



   }
TFile *fout = new TFile(Form("kaonUnfoldStandardGenFidVol_%1.1f_%1.1f.root",vol,vol2),"RECREATE");

std::cout<<h1d_reco_interactingE->GetEntries()<<std::endl;
      for (int i = 1; i <= h1d_reco_interactingE->GetNbinsX(); ++i) {
      h1d_reco_interactingE->SetBinError(i,TMath::Sqrt(h1d_reco_interactingE->GetBinContent(i)*(1-h1d_reco_interactingE->GetBinContent(i)/h1d_reco_incidentE->GetBinContent(i))));
      }


TH1D* h1d_recoCopy_interactingE=(TH1D*)h1d_reco_interactingE->Clone("h1d_reco_copy_interactingE");
TH1D* h1d_recoCopy_incidentE=(TH1D*)h1d_reco_incidentE->Clone("h1d_reco_copy_incidentE");
RooUnfoldBayes unfoldCorr(&response, h1d_recoCopy_interactingE,4);
unfoldCorr.Overflow();
unfoldCorr.IncludeSystematics();
unfoldCorr.SetNToys(1000);

TH1D* hRecoFullCorr=(TH1D*) unfoldCorr.Hreco(RooUnfold::kCovToy);
RooUnfoldBayes unfoldInc(&responseIncident, h1d_recoCopy_incidentE,4);
unfoldInc.Overflow();
unfoldInc.IncludeSystematics();
unfoldInc.SetNToys(1000);
//RooUnfold::ErrorTreatment withError=3;
TH1D* hRecoInc=(TH1D*)unfoldInc.Hreco(RooUnfold::kCovToy);
hRecoInc->SetName("fullCorrIncHist");

//RooUnfoldBayes unfoldCheat(&response, h1d_reco_cheat_interactingE,4);
//TH1D* hRecoFullCheat=(TH1D*) unfoldCheat.Hreco();
//TMatrixD* hMatrix=(TH1D*)unfoldCorr.Hreco();
hRecoFullCorr->SetName("fullCorrRecoHist");
TCanvas c1=TCanvas();
h1d_recoCopy_interactingE->Draw("HIST");
c1.Print("testRecoCopyInteractingE.png");



h1d_reco_interactingE->Write();
h1d_reco_cheat_interactingE->Write();
h1d_true_thinslice_incidentE->Write();
h1d_true_thinslice_interactingE->Write();

hRecoFullCorr->Write();
hRecoInc->Write();
       TH1D * xsec_hist = (TH1D*)hRecoFullCorr->Clone("crossSectionHist");

      xsec_hist->Divide(h1d_true_thinslice_incidentE);
      for (int i = 1; i <= xsec_hist->GetNbinsX(); ++i) {
        xsec_hist->SetBinContent(i, -1.*log(1. - xsec_hist->GetBinContent(i)));
//        if (hRecoFullCorr->GetBinContent(i)<500) xsec_hist->SetBinContent(i,0);
      }
      xsec_hist->Scale(1.E27/ (0.4979 * 1.4 * 6.022E23 / 39.948 ));



xsec_hist->GetXaxis()->SetTitle("True Interacting KE (MeV)");
xsec_hist->GetYaxis()->SetTitle("#sigma (mbarn)");
xsec_hist->Write();


       TH1D * xsec_true = (TH1D*)h1d_true_thinslice_interactingE->Clone("crossSectionTrue");

      xsec_true->Divide(h1d_true_thinslice_incidentE);
      for (int i = 1; i <= xsec_hist->GetNbinsX(); ++i) {
        xsec_true->SetBinContent(i, -1.*log(1. - xsec_true->GetBinContent(i)));
       if (h1d_true_thinslice_interactingE->GetBinContent(i)>h1d_true_thinslice_incidentE->GetBinContent(i)) xsec_true->SetBinContent(i,0);
      }
      xsec_true->Scale(1.E27/ (0.4979 * 1.4 * 6.022E23 / 39.948 ));
      for (int i = 1; i <= xsec_hist->GetNbinsX(); ++i){
        

           double recoErr=hRecoFullCorr->GetBinError(i)/hRecoFullCorr->GetBinContent(i);
           double trueIncErr=h1d_true_thinslice_incidentE->GetBinError(i)/h1d_true_thinslice_incidentE->GetBinContent(i);
           double trueErr=h1d_true_thinslice_interactingE->GetBinError(i)/h1d_true_thinslice_interactingE->GetBinContent(i);
           double tempReco1=TMath::Sqrt((recoErr*recoErr)+(trueIncErr*trueIncErr));
	   if(tempReco1>0){
           xsec_hist->SetBinError(i, tempReco1*xsec_hist->GetBinContent(i));
           }
           double tempTruth=TMath::Sqrt((trueErr*trueErr)+(trueIncErr*trueIncErr));
        if(tempTruth>0){
	xsec_true->SetBinError(i,tempTruth*xsec_true->GetBinContent(i));
        }

    }

xsec_true->GetXaxis()->SetTitle("True Interacting KE (MeV)");
xsec_true->GetYaxis()->SetTitle("#sigma (mbarn)");
xsec_true->Write();

TFile f("/dune/data/users/rdiurba/standard6GeVkaonTreeSel.root");
TTree* t=(TTree*)f.Get("ana");
double reco_beam_calibrated_interactingEnergy_data=0;
double beam_inst_KE_data=0;
std::vector<double>* reco_beam_calibrated_dEdX_SCE_data=0x0;
std::vector<double>* reco_beam_TrkPitch_SCE_data=0x0;
std::vector<double>* reco_beam_calo_Z_data=0x0;
std::vector<double>* reco_beam_calo_wire_data=0x0;
Int_t true_beam_ID;
Int_t reco_beam_true_byE_ID;
double reco_beam_endZ_data;
std::vector<std::vector<double>>* g4rw_full_grid_kplus_weights=0x0;
int selection_ID_data=0;
double true_beam_interactingEnergy=0;
t->SetBranchAddress("reco_beam_calibrated_interactingEnergy",&reco_beam_calibrated_interactingEnergy_data);
t->SetBranchAddress("reco_beam_calibrated_dEdX_SCE",&reco_beam_calibrated_dEdX_SCE_data);
t->SetBranchAddress("reco_beam_TrkPitch_SCE",&reco_beam_TrkPitch_SCE_data);
t->SetBranchAddress("selection_ID",&selection_ID_data);
t->SetBranchAddress("beam_inst_KE",&beam_inst_KE_data);
t->SetBranchAddress("reco_beam_calo_Z",&reco_beam_calo_Z_data);
t->SetBranchAddress("reco_beam_endZ",&reco_beam_endZ_data);
t->SetBranchAddress("reco_beam_calo_wire",&reco_beam_calo_wire_data);
t->SetBranchAddress("g4rw_full_grid_kplus_weights",&g4rw_full_grid_kplus_weights);
t->SetBranchAddress("reco_beam_true_byE_ID",&reco_beam_true_byE_ID);
t->SetBranchAddress("true_beam_ID",&true_beam_ID);
t->SetBranchAddress("true_beam_interactingEnergy",&true_beam_interactingEnergy);
   Long64_t nentries_data=t->GetEntries();
   //Long64_t nentries_start=0;
   Long64_t nentries_start=nentries_data*0.6667+1;
TH1D* h1d_reco_interactingE_data=new TH1D("h1d_reco_interactingE_data","h1d_reco_interactingE_data",totBins, edges);
TH1D* h1d_reco_incidentE_data=new TH1D("h1d_reco_incidentE_data","h1d_reco_incidentE_data",totBins, edges);
   for (Long64_t jentry=nentries_start; jentry<nentries_data;jentry++) {

   t->GetEntry(jentry);
double weight=1.f;//coeff;
 if (selection_ID_data==1  && reco_beam_endZ_data<220.0+0.239*vol2 && reco_beam_endZ_data>(30.0+1.751*vol)){ h1d_reco_interactingE_data->Fill(reco_beam_calibrated_interactingEnergy_data,weight);
}
if(selection_ID_data<3){
h1d_reco_interactingE_passBQT_data->Fill(reco_beam_calibrated_interactingEnergy_data);

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
if (reco_beam_calo_Z_data->at(index)<220.0+0.239*vol2  && reco_beam_calo_Z_data->at(index)>(30.0+1.751*vol)) h1d_reco_incidentE_data->Fill(interactingKE);
}

}
}
}
}






TH1D* h1d_recoCopy_interactingE_data=(TH1D*)h1d_reco_interactingE_data->Clone("h1d_reco_copy_interactingE_test");
TH1D* h1d_recoCopy_incidentE_data=(TH1D*)h1d_reco_incidentE_data->Clone("h1d_reco_copy_incidentE_test");

for(int k=0; k<h1d_recoCopy_interactingE_data->GetNbinsX(); k++){
//h1d_recoCopy_interactingE_data->SetBinError(k+1,0);
//h1d_recoCopy_incidentE_data->SetBinError(k+1,0);

}
auto* R=response.Hresponse();
auto* RInc=responseIncident.Hresponse();
RInc->SetName("responseIncident");
R->SetStats(0);

RooUnfoldBayes unfoldData(&response, h1d_recoCopy_interactingE_data,4);
RooUnfoldBayes unfoldIncData(&responseIncident, h1d_recoCopy_incidentE_data,4);
//unfoldData.SetVerbose(3);
unfoldData.SetNToys(1000);
unfoldIncData.SetNToys(1000);
unfoldData.IncludeSystematics();
unfoldIncData.IncludeSystematics();
//unfoldData.SetVerbose(2);
unfoldData.Overflow();
unfoldIncData.Overflow();
//unfoldData.RemoveZeros();
//unfoldIncData.RemoveZeros();
TH1D* hRecoFullData=(TH1D*) unfoldData.Hreco(RooUnfold::kCovToy);
TH1D* hRecoIncData=(TH1D*) unfoldIncData.Hreco(RooUnfold::kCovToy);
//unfoldData.Overflow();
//unfoldIncData.Overflow();

auto vmeasured=unfoldData.Vmeasured();
//std::cout<<h1d_recoCopy_interactingE_data->GetBinContent(10)<<','<<vmeasured.Print()<<std::endl;
h1d_recoCopy_interactingE_data->Print("all");
//auto emeas=unfoldData.GetCov();
//emeas.Print();
auto ereco=unfoldData.ErecoV(RooUnfold::kCovariance);
ereco.Print();
auto erecoM=unfoldData.Ereco(RooUnfold::kCovariance);
erecoM.Print();
hRecoFullData->SetName("fullCorrData");
hRecoIncData->SetName("fullCorrIncData");

       TH1D * xsec_data = (TH1D*)hRecoFullData->Clone("crossSectionTest");
      TH1D* h1d_true_thinslice_incidentE_scaled=(TH1D*)h1d_true_thinslice_incidentE->Clone("h1d_true_thinslice_incidentE_scaled");
      h1d_true_thinslice_incidentE_scaled->Scale(h1d_reco_interactingE_passBQT_data->Integral()/h1d_reco_interactingE_passBQT->Integral());
      std::cout<<hRecoFullData->Integral()/hRecoFullCorr->Integral()<<std::endl;
      xsec_data->Divide(h1d_true_thinslice_incidentE_scaled);
      for (int i = 1; i <= xsec_hist->GetNbinsX(); ++i) {
        xsec_data->SetBinContent(i, -1.*log(1. - xsec_data->GetBinContent(i)));
        //std::cout<<h1d_true_thinslice_incidentE_scaled->GetBinContent(i)<<','<<hRecoFullData->GetBinContent(i)<<std::endl;
        //if (hRecoFullData->GetBinContent(i)>h1d_true_thinslice_incidentE_scaled->GetBinContent(i)) xsec_data->SetBinContent(i,0);
      }
      xsec_data->Scale(1.E27/ (0.4979 * 1.4 * 6.022E23 / 39.948 ));


         for (int i = 1; i <= xsec_data->GetNbinsX(); ++i){


           double recoErr=hRecoFullData->GetBinError(i)/hRecoFullData->GetBinContent(i);
           double trueIncErr=h1d_true_thinslice_incidentE->GetBinError(i)/h1d_true_thinslice_incidentE->GetBinContent(i);
           double trueErr=h1d_true_thinslice_interactingE->GetBinError(i)/h1d_true_thinslice_interactingE->GetBinContent(i);
           double tempReco1=TMath::Sqrt((recoErr*recoErr)+(trueIncErr*trueIncErr));
	   if(tempReco1>0){
           xsec_data->SetBinError(i, tempReco1*xsec_data->GetBinContent(i));
           }
/*
           double diffErr=(TMath::Sqrt(hRecoFullData->GetBinError(i)*hRecoFullData->GetBinError(i)+TMath::Sqrt(h1d_true_thinslice_incidentE_scaled->GetBinContent(i))*TMath::Sqrt(h1d_true_thinslice_incidentE_scaled->GetBinContent(i))))/(h1d_true_thinslice_incidentE_scaled->GetBinContent(i)-hRecoFullData->GetBinContent(i));
          
           double numErr=TMath::Sqrt(h1d_true_thinslice_incidentE_scaled->GetBinContent(i))/h1d_true_thinslice_incidentE_scaled->GetBinContent(i);

           double inLogErr=TMath::Sqrt((diffErr*diffErr)+(numErr*numErr));

           double postLogErr=inLogErr;


           double percErr=postLogErr/(TMath::Log(h1d_true_thinslice_incidentE_scaled->GetBinContent(i)/(h1d_true_thinslice_incidentE_scaled->GetBinContent(i)-hRecoFullData->GetBinContent(i))));
           std::cout<<percErr<<std::endl;

	   if(percErr>0){
           xsec_data->SetBinError(i, percErr*xsec_data->GetBinContent(i));
           }
*/
       }


       TH1D * xsec_rawMC = (TH1D*)h1d_reco_interactingE->Clone("crossSectionRawMC");

      xsec_rawMC->Divide(h1d_reco_incidentE);
      for (int i = 1; i <= xsec_hist->GetNbinsX(); ++i) {
        xsec_rawMC->SetBinContent(i, -1.*log(1. - xsec_rawMC->GetBinContent(i)));
        //std::cout<<h1d_true_thinslice_incidentE_scaled->GetBinContent(i)<<','<<hRecoFullData->GetBinContent(i)<<std::endl;
        //if (hRecoFullData->GetBinContent(i)>h1d_true_thinslice_incidentE_scaled->GetBinContent(i)) xsec_data->SetBinContent(i,0);
      }
      xsec_rawMC->Scale(1.E27/ (0.4979 * 1.4 * 6.022E23 / 39.948 ));
      for (int i = 1; i <= xsec_hist->GetNbinsX(); ++i){
        

           double recoErr=h1d_reco_interactingE->GetBinError(i)/h1d_reco_interactingE->GetBinContent(i);
           double recoIncErr=h1d_reco_incidentE->GetBinError(i)/h1d_reco_incidentE->GetBinContent(i);
           //double trueIncErr=h1d_true_thinslice_incidentE->GetBinError(i)/h1d_true_thinslice_incidentE->GetBinContent(i);
           //double trueErr=h1d_true_thinslice_interactingE->GetBinError(i)/h1d_true_thinslice_interactingE->GetBinContent(i);
           double tempReco1=TMath::Sqrt((recoErr*recoErr)+(recoIncErr*recoIncErr));
        if(tempReco1>0){
	xsec_rawMC->SetBinError(i,tempReco1*xsec_rawMC->GetBinContent(i));
        }

    }

       TH1D * xsec_unfoldData = (TH1D*)hRecoFullData->Clone("crossSectionUnfoldTest");

     xsec_unfoldData->Divide(hRecoIncData);
      for (int i = 1; i <= xsec_hist->GetNbinsX(); ++i) {
        xsec_unfoldData->SetBinContent(i, -1.*log(1. - xsec_unfoldData->GetBinContent(i)));
        //std::cout<<h1d_true_thinslice_incidentE_scaled->GetBinContent(i)<<','<<hRecoFullData->GetBinContent(i)<<std::endl;
        //if (hRecoFullData->GetBinContent(i)>h1d_true_thinslice_incidentE_scaled->GetBinContent(i)) xsec_data->SetBinContent(i,0);
      }
      xsec_unfoldData->Scale(1.E27/ (0.4979 * 1.4 * 6.022E23 / 39.948 ));

    float rho = 1396; //kg/m^3
    float molar_mass = 39.95; //g/mol
    float g_per_kg = 1000; 
    float avogadro = 6.022e+23; //number/mol
    float number_density = rho*g_per_kg/molar_mass*avogadro;
    //float slab_width = 0.0047;//in m
    float slab_width = 0.004792;//in m
      for (int i = 1; i <= xsec_hist->GetNbinsX(); ++i) {
         float term1=hRecoFullData->GetBinError(i)/hRecoFullData->GetBinContent(i);
         float term2=hRecoIncData->GetBinError(i)/hRecoIncData->GetBinContent(i);
	  float totalError = (xsec_unfoldData->GetBinContent(i)) * (pow( ( (term1*term1) + (term2*term2) ),0.5)) * (1/number_density) * (1/slab_width)  *(1e26);
        std::cout<<totalError<<std::endl;
      // if(hRecoFullData->GetBinContent(i)>0) xsec_unfoldData->SetBinError(i, totalError);
        //std::cout<<h1d_true_thinslice_incidentE_scaled->GetBinContent(i)<<','<<hRecoFullData->GetBinContent(i)<<std::endl;
        //if (hRecoFullData->GetBinContent(i)>h1d_true_thinslice_incidentE_scaled->GetBinContent(i)) xsec_data->SetBinContent(i,0);
      }


       TH1D * xsec_unfoldMC = (TH1D*)hRecoFullCorr->Clone("crossSectionUnfoldMC");

     xsec_unfoldMC->Divide(hRecoInc);
      for (int i = 1; i <= xsec_hist->GetNbinsX(); ++i) {
        xsec_unfoldMC->SetBinContent(i, -1.*log(1. - xsec_unfoldMC->GetBinContent(i)));
        //std::cout<<h1d_true_thinslice_incidentE_scaled->GetBinContent(i)<<','<<hRecoFullData->GetBinContent(i)<<std::endl;
        //if (hRecoFullData->GetBinContent(i)>h1d_true_thinslice_incidentE_scaled->GetBinContent(i)) xsec_data->SetBinContent(i,0);

      }
      xsec_unfoldMC->Scale(1.E27/ (0.4979* 1.4 * 6.022E23 / 39.948 ));

//std::cout<<unfoldData.GetSmoothing()<<std::endl;
//std::cout<<unfoldData.GetRegParm()<<std::endl;
//std::cout<<unfoldData.GetIterations()<<std::endl;
//std::cout<<unfoldData.NToys()<<std::endl;
//xsec_unfoldData->Fit("pol0","O");
//xsec_unfoldMC->Fit("pol0","O");



TH2D* covMat=new TH2D("covMatStat","covMatStat",totBins, edges , totBins, edges );


TH2D* covMatMC=new TH2D("covMatMCStat","covMatMCStat",totBins, edges , totBins, edges );
   auto ttime=new TTimeStamp();
   auto time=ttime->AsDouble();


TProfile* h1d_MCXSec=new TProfile("h1d_MCXSecStat","h1d_MCXSecStat",totBins, edges , "s");
TProfile* h1d_DataXSec=new TProfile("h1d_DataXSecStat","h1d_DataXSecStat",totBins, edges , "s");


TProfile* h1d_MCInt=new TProfile("h1d_MCIntStat","h1d_MCIntStat",totBins, edges , "s");
TProfile* h1d_DataInt=new TProfile("h1d_DataIntStat","h1d_DataIntStat",totBins, edges , "s");

TProfile* h1d_MCInc=new TProfile("h1d_MCIncStat","h1d_MCIncStat",totBins, edges , "s");
TProfile* h1d_DataInc=new TProfile("h1d_DataIncStat","h1d_DataIncStat",totBins, edges , "s");

std::vector<std::vector<double>> xsec_unfoldDataVecVec;
std::vector<std::vector<double>> xsec_unfoldMCVecVec;
std::vector<TH1D> xsec_unfoldDataVecHist;
std::vector<TH1D> xsec_unfoldMCVecHist;
   gRandom->SetSeed(int(time+time*2)+int(time));
int totToys=1000;
for (Int_t k=0; k<totToys; k++){
//    gRandom->SetSeed(int(time+time*2));
    //double sigma=gRandom->Gaus(0,1);
    //double sigma2=gRandom->Gaus(0,1);
    //double sigma3=gRandom->Gaus(0,1);
    //double sigma4=gRandom->Gaus(0,1);
TH1D* h1d_recoCopyFluc_interactingE=(TH1D*)h1d_reco_interactingE->Clone("h1d_reco_copyFluc_interactingE");
TH1D* h1d_recoCopyFluc_incidentE=(TH1D*)h1d_reco_incidentE->Clone("h1d_reco_copyFluc_incidentE");
TH1D* h1d_recoCopyFluc_interactingE_data=(TH1D*)h1d_reco_interactingE_data->Clone("h1d_reco_copyFluc_interactingE_data");
TH1D* h1d_recoCopyFluc_incidentE_data=(TH1D*)h1d_reco_incidentE_data->Clone("h1d_reco_copyFluc_incidentE_data");
      for (int i = 1; i <= xsec_hist->GetNbinsX(); ++i) {

double sigma=gRandom->Gaus(0,1);
double sigma2=gRandom->Gaus(0,1);
double sigma3=gRandom->Gaus(0,1);
double sigma4=gRandom->Gaus(0,1);


double mcUncert=TMath::Sqrt(h1d_recoCopyFluc_interactingE->GetBinContent(i)-(h1d_recoCopyFluc_interactingE->GetBinContent(i)/h1d_recoCopyFluc_incidentE->GetBinContent(i)));
double dataUncert=TMath::Sqrt(h1d_recoCopyFluc_interactingE_data->GetBinContent(i)-(h1d_recoCopyFluc_interactingE_data->GetBinContent(i)/h1d_recoCopyFluc_incidentE_data->GetBinContent(i)));
//std::cout<<mcUncert<<','<<dataUncert<<std::endl;
//std::cout<<mcUncert*sigma<<std::endl;
//std::cout<<sigma<<','<<sigma2<<','<<sigma3<<','<<sigma4<<std::endl;
h1d_recoCopyFluc_interactingE->SetBinContent(i,mcUncert*sigma+h1d_recoCopyFluc_interactingE->GetBinContent(i));
h1d_recoCopyFluc_interactingE_data->SetBinContent(i,dataUncert*sigma2+h1d_recoCopyFluc_interactingE_data->GetBinContent(i));



h1d_recoCopyFluc_incidentE->SetBinContent(i,TMath::Sqrt(h1d_recoCopyFluc_incidentE->GetBinContent(i))*sigma3+h1d_recoCopyFluc_incidentE->GetBinContent(i));
h1d_recoCopyFluc_incidentE_data->SetBinContent(i,TMath::Sqrt(h1d_recoCopyFluc_incidentE_data->GetBinContent(i))*sigma4+h1d_recoCopyFluc_incidentE_data->GetBinContent(i));
	

	}
RooUnfoldBayes unfoldMCIntCopy(&response, h1d_recoCopyFluc_interactingE,4);
RooUnfoldBayes unfoldMCIncCopy(&responseIncident, h1d_recoCopyFluc_incidentE,4);
RooUnfoldBayes unfoldDataIntCopy(&response, h1d_recoCopyFluc_interactingE_data,4);
RooUnfoldBayes unfoldDataIncCopy(&responseIncident, h1d_recoCopyFluc_incidentE_data,4);
TH1D* recoMCIncCopy=(TH1D*) unfoldMCIncCopy.Hreco();
TH1D* recoMCIntCopy=(TH1D*) unfoldMCIntCopy.Hreco();
TH1D* recoDataIncCopy=(TH1D*) unfoldDataIncCopy.Hreco();
TH1D* recoDataIntCopy=(TH1D*) unfoldDataIntCopy.Hreco();
TH1D * xsec_unfoldMCFluc = (TH1D*)recoMCIntCopy->Clone("crossSectionUnfoldFluc");
TH1D * xsec_unfoldDataFluc = (TH1D*)recoDataIntCopy->Clone("crossSectionUnfoldDataFluc");
//std::cout<<recoMCIncCopy->GetBinContent(1)<<','<<recoMCIntCopy->GetBinContent(1)<<','<<recoDataIncCopy->GetBinContent(1)<<','<<recoDataIntCopy->GetBinContent(1)<<std::endl;
     xsec_unfoldMCFluc->Divide(recoMCIncCopy);
     xsec_unfoldDataFluc->Divide(recoDataIncCopy);
 // std::cout<<xsec_unfoldMCFluc->GetBinContent(1)<<','<<xsec_unfoldDataFluc->GetBinContent(1)<<std::endl;    
//std::cout<<recoMCIncCopy->GetBinContent(1)<<','<<recoMCIntCopy->GetBinContent(1)<<','<<recoDataIncCopy->GetBinContent(1)<<','<<recoDataIntCopy->GetBinContent(1)<<std::endl;

     for (int i = 1; i <= xsec_unfoldData->GetNbinsX(); ++i) {
       

        xsec_unfoldDataFluc->SetBinContent(i, -1.*log(1. - xsec_unfoldDataFluc->GetBinContent(i)));
        xsec_unfoldMCFluc->SetBinContent(i, -1.*log(1. - xsec_unfoldMCFluc->GetBinContent(i)));
   
      }
      xsec_unfoldDataFluc->Scale(1.E27/ (0.4979 * 1.4 * 6.022E23 / 39.948 ));
      xsec_unfoldMCFluc->Scale(1.E27/ (0.4979 * 1.4 * 6.022E23 / 39.948 ));
 	std::vector<double> tmpVec;
std::vector<double> tmpVecMC;

      for (int i = 1; i <=xsec_unfoldData->GetNbinsX(); ++i) {
        h1d_MCXSec->Fill(xsec_unfoldMCFluc->GetBinCenter(i),xsec_unfoldMCFluc->GetBinContent(i));
        h1d_DataXSec->Fill(xsec_unfoldDataFluc->GetBinCenter(i),xsec_unfoldDataFluc->GetBinContent(i));
        tmpVec.push_back(xsec_unfoldDataFluc->GetBinContent(i));
        tmpVecMC.push_back(xsec_unfoldMCFluc->GetBinContent(i));
        h1d_MCInc->Fill(xsec_unfoldMC->GetBinCenter(i),recoMCIncCopy->GetBinContent(i));
        h1d_MCInt->Fill(xsec_unfoldMC->GetBinCenter(i),recoMCIntCopy->GetBinContent(i));

        h1d_DataInc->Fill(xsec_unfoldMC->GetBinCenter(i),recoDataIncCopy->GetBinContent(i));
        h1d_DataInt->Fill(xsec_unfoldMC->GetBinCenter(i),recoDataIntCopy->GetBinContent(i));


      }
TH1D* tempHist=(TH1D*)xsec_unfoldDataFluc->Clone(Form("toy_%d",k));
TH1D* tempHistMC=(TH1D*)xsec_unfoldMCFluc->Clone(Form("toyMC_%d",k));
xsec_unfoldDataVecHist.push_back(*tempHist);
xsec_unfoldMCVecHist.push_back(*tempHistMC);
xsec_unfoldDataVecVec.push_back(tmpVec);
xsec_unfoldMCVecVec.push_back(tmpVecMC);
std::cout<<"On Toy:"<<k<<std::endl;
//    delete unfoldDataIncCopy;
//    delete unfoldDataIntCopy;
//    delete unfoldMCIncCopy;
//    delete unfoldMCIntCopy;
//std::cout<<tmpVec.at(4)<<','<<tmpVec.at(5)<<std::endl;
}

for(long unsigned int toy=0; toy<xsec_unfoldDataVecVec.size(); toy++){
    auto dataVec=xsec_unfoldDataVecVec.at(toy);
    auto mcVec=xsec_unfoldMCVecVec.at(toy);
    TH1D tmpHist=xsec_unfoldDataVecHist.at(toy);
    TH1D tmpHistMC=xsec_unfoldMCVecHist.at(toy);
 //   std::cout<<dataVec.at(4)<<','<<dataVec.at(5)<<std::endl;
    //std::cout<<tmpHist->GetBinContent(5)<<std::endl;
  for (Int_t bin=0; bin<totBins; ++bin){
    for (Int_t bin2=0; bin2<totBins; ++bin2){
    int binNum=(bin)*totBins+(bin2);
    double diffBin1MC=tmpHistMC.GetBinContent(bin+1)-h1d_MCXSec->GetBinContent(bin+1);
    double diffBin2MC=tmpHistMC.GetBinContent(bin2+1)-h1d_MCXSec->GetBinContent(bin2+1);
    //if(bin>40) std::cout<<diffBin1<<std::endl;
    covMatMC->SetBinContent(bin+1,bin2+1,covMatMC->GetBinContent(bin+1,bin2+1)+((diffBin1MC*diffBin2MC)/(xsec_unfoldDataVecVec.size()-1)));
    double diffBin1=tmpHist.GetBinContent(bin+1)-h1d_DataXSec->GetBinContent(bin+1);
    double diffBin2=tmpHist.GetBinContent(bin2+1)-h1d_DataXSec->GetBinContent(bin2+1);
//if (bin==2 && bin2==3)  std::cout<<tmpHist.GetBinContent(bin+1)<<","<<tmpHist.GetBinContent(bin2+1)<<","<<h1d_DataXSec->GetBinContent(bin+1)<<","<<diffBin1<<","<<diffBin2<<std::endl; 
//  if(bin==4 && bin2==5) std::cout<<diffBin1<<','<<diffBin2<<','<<dataVec.at(bin)<<","<<dataVec.at(bin2)<<','<<h1d_DataXSec->GetBinContent(bin+1)<<','<<h1d_DataXSec->GetBinContent(bin2+1)<<","<<bin<<','<<bin2<<std::endl;
    covMat->SetBinContent(bin+1,bin2+1,covMat->GetBinContent(bin+1,bin2+1)+((diffBin1*diffBin2)/(xsec_unfoldDataVecVec.size()-1)));
  }
  }
}
   TMatrixD erecoData(totBins,totBins);
   TMatrixD erecoMC(totBins,totBins);


   for (int i = 0; i < totBins; ++i) {
    for (int j = 0; j < totBins; ++j) {
    if(i==j) std::cout<<covMatMC->GetBinContent(i+1,j+1)<<std::endl;
    erecoMC[i][j] = covMatMC->GetBinContent(i+1, j+1);
    erecoData[i][j] = covMat->GetBinContent(i+1, j+1);
  }}
 TDecompChol invMC(erecoMC);
 TDecompChol invData(erecoData);

bool test=invMC.Decompose();
bool test2=invData.Decompose();
std::cout<<test<<","<<test2<<std::endl;
auto new_matMC=invMC.Invert();
auto new_matData=invData.Invert();


double chi2=0;
double chi2Inc=0;
double chi2XSec=0;
for(int bin=0; bin<h1d_DataXSec->GetNbinsX(); ++bin){
    for(int bin2=0; bin2<h1d_DataXSec->GetNbinsX(); ++bin2){
    double weight=1.f;//coeff;
    double weight2=1.f;//coeff;
    double diffBin1=xsec_unfoldData->GetBinContent(bin+1)-weight*g4Truth->Eval(xsec_unfoldData->GetBinCenter(bin+1));
    double diffBin2=xsec_unfoldData->GetBinContent(bin2+1)-weight2*g4Truth->Eval(xsec_unfoldData->GetBinCenter(bin2+1));
//    std::cout<<diffBin1<<","<<diffBin2<<","<<erecoData[bin][bin2]<<std::endl;

    double bin_cont=diffBin1*(new_matData)(bin,bin2)*diffBin2;
    chi2=chi2+bin_cont;
    double diffIncBin1=xsec_unfoldMC->GetBinContent(bin+1)-g4Truth->Eval(xsec_unfoldData->GetBinCenter(bin+1));
    double diffIncBin2=xsec_unfoldMC->GetBinContent(bin2+1)-g4Truth->Eval(xsec_unfoldData->GetBinCenter(bin2+1));
    double binInc_cont=diffIncBin1*(new_matMC)(bin,bin2)*diffIncBin2;
   // std::cout<<""<<std::endl;
   //std::cout<<diffIncBin2<<","<<diffIncBin2<<","<<new_matMC(bin,bin2)<<std::endl;


    chi2Inc=chi2Inc+binInc_cont;


}
}
std::cout<<chi2<<","<<chi2Inc<<std::endl;



h1d_DataXSec->SetTitle(Form("%1.3f",chi2));
xsec_unfoldData->SetTitle(Form("%1.3f",chi2));


xsec_unfoldData->Fit("pol0","NO");

fout->cd();
xsec_rawMC->Write();
hRecoFullData->Write();
h1d_reco_incidentE->Write();
h1d_reco_incidentE_data->Write();
xsec_data->Write();
hRecoFullData->Write();
hRecoIncData->Write();
h1d_reco_interactingE_data->Write();
xsec_unfoldData->Write();
xsec_unfoldMC->Write();
covMatMC->Write();
covMat->Write();
h1d_MCXSec->Write();
h1d_DataXSec->Write();

h1d_MCInc->Write();
h1d_DataInc->Write();


h1d_MCInt->Write();
h1d_DataInt->Write();


R->Write();
RInc->Write();


}
