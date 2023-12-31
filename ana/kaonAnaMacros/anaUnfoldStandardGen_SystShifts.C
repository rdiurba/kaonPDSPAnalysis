#define anaUnfoldStandardGen_SystShifts_cxx
#include "anaUnfoldStandardGen_SystShifts.h"
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


void anaUnfoldStandardGen_SystShifts::Loop(int nUniverses, std::string weightDecision, int oneShift)
{
//   In a ROOT session, you can do:
//      root> .L anaUnfoldStandardGen_SystShifts.C
//      root> anaUnfoldStandardGen_SystShifts t
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
//double edges[7]={3000,5000,5200,5400,5600,5900,6800};
//double edges[6]={4000,5000,5240,5440,5700,6650};
//double edges[6]={4480,4980,5220,5410,5630,6200};
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
double startShift=-1; double endShift=-1; double  dEdxshift=-1; double beamshift=-1; double  ediv_var=-1; double  nobeam_var=-1; double long_var=-1; double short_var=-1; double  g4rw_index=-1;
double  beamScraper_var=-1;

std::vector<double> xSecMultisim, multisimdEdXShift, multisimBeamShift, multisimDataInc, multisimDataInt;
std::vector<double>* xsec_unfoldDataVec=0x0;
std::vector<double>* xsec_unfoldMCVec=0x0;
std::vector<TH1D*> throwHists;
std::vector<TH1D*> throwHistsMC;
std::vector<std::vector<double>> xsec_unfoldDataVecVec;
std::vector<std::vector<double>> xsec_unfoldMCVecVec;
std::vector<TH1D> xsec_unfoldDataVecHist;
std::vector<TH1D> xsec_unfoldMCVecHist;
TFile *fWrite=new TFile(Form("kaonSystShiftOutputStandardGeni_%s_%d_%d.root",weightDecision.c_str(),nUniverses,oneShift),"recreate");
TTree tree("systShift","systShiftTree");
tree.Branch("xsec_unfoldDataVec",&xsec_unfoldDataVec);
tree.Branch("xsec_unfoldMCVec",&xsec_unfoldMCVec);

tree.Branch("dEdxshift",&dEdxshift,"dEdxshift/D");
tree.Branch("beamshift",&beamshift,"beamshift/D");
tree.Branch("ediv_var",&ediv_var,"ediv_var/D");
tree.Branch("nobeam_var",&nobeam_var,"nobeam_var/D");
tree.Branch("long_var",&long_var,"long_var/D");
tree.Branch("short_var",&short_var,"short_var/D");
tree.Branch("g4rw_index",&g4rw_index,"g4rw_index/D");
tree.Branch("beamScraper_var",&beamScraper_var,"beamScraper_var/D");


TProfile* h1d_MCXSec=new TProfile("h1d_MCXSec","h1d_MCXSec",totBins, edges, "s");
TProfile* h1d_DataXSec=new TProfile("h1d_DataXSec","h1d_DataXSec",totBins, edges, "s");

TProfile* h1d_MCInt=new TProfile("h1d_MCInt","h1d_MCInt",totBins, edges, "s");
TProfile* h1d_DataInt=new TProfile("h1d_DataInt","h1d_DataInt",totBins, edges, "s");

TProfile* h1d_MCInc=new TProfile("h1d_MCInc","h1d_MCInc",totBins, edges, "s");
TProfile* h1d_DataInc=new TProfile("h1d_DataInc","h1d_DataInc",totBins, edges, "s");


  TFile shift_file("beamResMap.root");
  TGraph* fSystBeamShiftMeans = (TGraph*)shift_file.Get("gMeans");
  //fSystBeamShiftMeans->SetDirectory(0);
  TGraph* fSystBeamShiftWidths = (TGraph*)shift_file.Get("gWidths");
   auto t=new TTimeStamp();
   auto time=t->AsDouble();
    std::cout<<time<<std::endl;
gRandom->SetSeed(int(time+time+1));

  double relative_nobeam_var=TMath::Sqrt(fChain->GetEntries("selection_ID==4"))/fChain->GetEntries("selection_ID==4");
  double relative_ediv_var=TMath::Sqrt(fChain->GetEntries("selection_ID==3"))/fChain->GetEntries("selection_ID==3");
  double relative_brkTrk_var=TMath::Sqrt(fChain->GetEntries("selection_ID<3 && reco_beam_daughter_sameID==1"))/fChain->GetEntries("selection_ID<3 && reco_beam_daughter_sameID==1");
if (weightDecision=="none") nUniverses=1;


if (abs(oneShift)==1) nUniverses=1;
for(int universe=0; universe<nUniverses; universe++){
   std::cout<<"Universe: "<<universe<<std::endl;
   dEdxshift=gRandom->Gaus(1,0.03);
   beamshift=gRandom->Gaus(1,0.012);
    startShift=gRandom->Gaus(0,1.751);
    endShift=gRandom->Gaus(0,0.239);
 
   std::cout<<"Generating variables"<<std::endl;
    ediv_var=gRandom->Gaus(1.0,1.00);
    if (ediv_var<0) ediv_var=0;
   //nobeam_var=gRandom->Gaus(0.953,0.1);
   nobeam_var=gRandom->Gaus(1.0,0.06);
   if (weightDecision=="none"){
   nobeam_var=1.0;//0.976;
   ediv_var=1.0;
   }
   //double nobeam_var=gRandom->Gaus(1.0,0.5);
   if (nobeam_var<0) nobeam_var=0;
   long_var=gRandom->Gaus(1,0.03);
   if (long_var<0) long_var=0;
   short_var=gRandom->Gaus(1,0.2);
   if (short_var<0) short_var=0;
   double extTrk_var=gRandom->Gaus(1,1);
   if (extTrk_var<0) extTrk_var=0;
   double brkTrk_var=gRandom->Gaus(1,1);
   if (brkTrk_var<0) brkTrk_var=0;

   beamScraper_var=1+abs(gRandom->Gaus(0,2.15));
   //if (beamScraper_var<1) beamScraper_var=1+(1-beamScraper_var);
   if (beamScraper_var<0) beamScraper_var=0;
   if (weightDecision=="none" || weightDecision=="beamScraper"){
   ediv_var=1.0;
   nobeam_var=1.0;

}   

   double bqt_var=gRandom->Gaus(0.92,0.1);
   //double bqt_var=gRandom->Gaus(1.0,0.23);
   if (bqt_var>2) bqt_var=1.99;
   if (bqt_var<0) bqt_var=0.01;


    g4rw_index=gRandom->Gaus(1,0.2);
   //double g4rw_var=0.8;
   double k0_var=gRandom->Gaus(1,0.2);
  
 /* Get Beam Weight */  
 double x_val = gRandom->Gaus(0,1);

 if (abs(oneShift)==1){
  startShift=1.751*oneShift;
  endShift=0.239*oneShift;
  dEdxshift=1+0.03*oneShift;
  beamshift=1+0.012*oneShift;

  short_var=1+0.2*oneShift;
  long_var=1+0.03*oneShift;
  beamScraper_var=1+2.15*abs(oneShift);
  g4rw_index=1+0.2*oneShift;
  k0_var=1+0.2*oneShift;
  brkTrk_var=1+1*oneShift;
  extTrk_var=1+1*oneShift;


   nobeam_var=1.0;
   ediv_var=1.0;

  if (weightDecision=="ediv" || weightDecision=="all") ediv_var=1.0+oneShift;
  if (weightDecision=="nobeam" || weightDecision=="all") nobeam_var=1.0+0.06*oneShift;
}


std::cout<<"Throws generated"<<std::endl;
   




   h1d_reco_incidentE->Reset();
   h1d_reco_interactingE->Reset();
   h1d_reco_incidentE_data->Reset();
   h1d_reco_interactingE_data->Reset();
   xsec_unfoldDataVec->clear();
   xsec_unfoldMCVec->clear();

   RooUnfoldResponse response(h1d_reco_interactingE, h1d_true_thinslice_interactingE);
   RooUnfoldResponse responseIncident(h1d_reco_incidentE, h1d_true_thinslice_incidentE);

response.UseOverflow();
responseIncident.UseOverflow();
std::cout<<response.UseOverflowStatus()<<std::endl;
   Long64_t nbytes = 0, nb = 0;
   nentries=fChain->GetEntries()*0.66667;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //std::cout<<"Running threw stuff"<<std::endl;    




      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
  double beam_weight=1.0;
  

  /* Get Ediv */
  double ediv_weight=1;
  if (selection_ID==2){
  if (reco_beam_endZ > 220.0 && reco_beam_endZ<234.0) {
    ediv_weight=ediv_var;
  }
  else {
    ediv_weight = (1-ediv_var*0.54025)/(1-0.54025);
  }

  
  }

  
double nobeam_weight=1;
  if(selection_ID>3) nobeam_weight=nobeam_var;
  else nobeam_weight=(1-0.55323*nobeam_var)/(1-0.55323);
double matched_weight=1;
if (true_beam_traj_incidentEnergies->size()>0 && selection_ID<4){
if (true_beam_endZ<30){
if (reco_beam_true_byE_ID==true_beam_ID) matched_weight=short_var;
else matched_weight=(1-0.04337*short_var)/(1-0.04337);
// 0.043373 was 0.374
}
else{
if (reco_beam_true_byE_ID==true_beam_ID) matched_weight=long_var;
else matched_weight=(1-0.8118*long_var)/(1-0.8118);
// 0.8118 was 0.783

}
}

if (matched_weight<0) matched_weight=0.00;
if (nobeam_weight<0) nobeam_weight=0.000;
double beamScraper_weight=1.f;
if(true_beam_endZ>0 && selection_ID<3 ) {

if ((beam_inst_KE-true_beam_traj_Eff)>200) beamScraper_weight=beamScraper_var;
else beamScraper_weight=(1-0.03959*beamScraper_var)/(1-0.03959);

}

double g4rw_weight=1.0;   
if(true_beam_endZ>-0.5 && true_beam_endZ<700){
std::vector<double> coeffs=g4rw_full_grid_kplus_coeffs->at(0);
g4rw_weight=coeffs.at(0)+coeffs.at(1)*g4rw_index+coeffs.at(2)*g4rw_index*g4rw_index+coeffs.at(3)*g4rw_index*g4rw_index*g4rw_index+coeffs.at(4)*g4rw_index*g4rw_index*g4rw_index*g4rw_index+coeffs.at(5)*g4rw_index*g4rw_index*g4rw_index*g4rw_index*g4rw_index+coeffs.at(6)*TMath::Power(g4rw_index,6)+coeffs.at(7)*TMath::Power(g4rw_index,7)+coeffs.at(8)*TMath::Power(g4rw_index,8)+coeffs.at(9)*TMath::Power(g4rw_index,9);
//if (g4rw_weight>2) g4rw_weight=2;
if (g4rw_weight<0.001) g4rw_weight=0.00;
//std::cout<<"Weight: "<<g4rw_weight<<std::endl;
coeffs.clear();
}
double k0_weight=1.f;
int numPart=0;
int numKaon=0;
       int  nKaonPlus=0; 
       int nKaon0=0;
       int nKaonMinus=0;
      for(long unsigned int i=0; i<true_beam_daughter_PDG->size(); i++){


        if (true_beam_daughter_PDG->at(i)==321) nKaonPlus++;
        if (true_beam_daughter_PDG->at(i)==310 || true_beam_daughter_PDG->at(i)==130) nKaon0++;
        if (true_beam_daughter_PDG->at(i)==-321) nKaonMinus++;
	} 

if ( true_beam_endZ<700 && true_beam_endZ>-0.5){
if (nKaon0==1 && nKaonMinus==0 && nKaonPlus==0) k0_weight=k0_var;
else k0_weight=(1-k0_var*3170.000/6578.000)/(1-3170.000/6578.000);
}
if (k0_weight<0) k0_weight=0.0000;


double extTrk_weight=1.f;
double brkTrk_weight=1.f;
if (true_beam_endZ>30 && true_beam_endZ<220 && reco_beam_true_byE_ID==true_beam_ID && selection_ID<3){
if (reco_beam_calo_endZ-true_beam_endZ>10) extTrk_weight=extTrk_var;
else extTrk_weight=(1-extTrk_var*0.14063)/(1-0.14063);


if (reco_beam_calo_endZ-true_beam_endZ<-10) brkTrk_weight=brkTrk_var;
else brkTrk_weight=(1-brkTrk_var*0.02982)/(1-0.02982);

}


if (brkTrk_weight<0) brkTrk_weight=0.0;
if (extTrk_weight<0) extTrk_weight=0.0;


   
	  double tot_weight=brkTrk_weight*extTrk_weight*beamScraper_weight*k0_weight*g4rw_weight*ediv_weight*nobeam_weight;
if (weightDecision=="none"){
tot_weight=ediv_weight*nobeam_weight;
beamshift=1.f;
dEdxshift=1.f;
startShift=0;
endShift=0;
}

else if (weightDecision=="ediv"){
tot_weight=ediv_weight;
beamshift=1.f;
dEdxshift=1.f;
startShift=0;
endShift=0;
}
else if(weightDecision=="calo"){
startShift=0;
endShift=0;
tot_weight=ediv_weight*nobeam_weight;
}
else if(weightDecision=="calodEdx"){
startShift=0;
endShift=0;
beamshift=1.f;
tot_weight=ediv_weight*nobeam_weight;
}
else if(weightDecision=="caloBeam"){
startShift=0;
endShift=0;
dEdxshift=1.f;
tot_weight=ediv_weight*nobeam_weight;
}
else if (weightDecision=="nobeam"){
tot_weight=nobeam_weight;
beamshift=1.f;
dEdxshift=1.f;
startShift=0;
endShift=0;
}
else if (weightDecision=="matchbeam"){
tot_weight=ediv_weight*nobeam_weight*matched_weight;
beamshift=1.f;
dEdxshift=1.f;
startShift=0;
endShift=0;
}
else if (weightDecision=="beamScraper"){
//std::cout>>beamScraper_weight
tot_weight=ediv_weight*nobeam_weight*beamScraper_weight;
beamshift=1.f;
dEdxshift=1.f;
startShift=0;
endShift=0;

}
else if (weightDecision=="g4rw"){
tot_weight=ediv_weight*nobeam_weight*g4rw_weight;
beamshift=1.f;
dEdxshift=1.f;
endShift=0;
startShift=0;

}
else if (weightDecision=="k0rw"){

tot_weight=ediv_weight*nobeam_weight*k0_weight;
beamshift=1.f;
dEdxshift=1.f;
endShift=0;
startShift=0;


}
else if(weightDecision=="g4rwTest"){
tot_weight=g4rw_weight;
beamshift=1.f;
dEdxshift=1.f;
startShift=0;
endShift=0;

}
else if (weightDecision=="brk"){
startShift=0;
endShift=0;
beamshift=1.f;
dEdxshift=1.f;
tot_weight=ediv_weight*nobeam_weight*brkTrk_weight;

}
else if(weightDecision=="ext"){
startShift=0;
endShift=0;
beamshift=1.f;
dEdxshift=1.f;
tot_weight=ediv_weight*nobeam_weight*extTrk_weight;
}
else if(weightDecision=="sceFront"){
endShift=0;
beamshift=1.f;
dEdxshift=1.f;
tot_weight=ediv_weight*nobeam_weight;



}
else if (weightDecision=="sceBack"){

startShift=0;
beamshift=1.f;
dEdxshift=1.f;
tot_weight=ediv_weight*nobeam_weight;




}


else weightDecision="all";  
//double tot_weight=g4rw_primary_grid_weights->at(g4rw_index); 
//std::cout<<"Generated weights"<<std::endl;   
double initialKE=beam_inst_KE;
double interactingKE=beam_inst_KE;
double interactingKERW=(beam_inst_KE*beamshift);

 
   if (true_beam_traj_incidentEnergies->size()){
   // h1d_true_thinslice_incidentE->Fill(5000,true_beam_traj_incidentEnergies->size());
   size_t finalSliceIndex=0;  
   int finalSlice = (true_beam_traj_slice_index)->at(true_beam_traj_incidentEnergies->size()-1); 
   for (size_t j = 0; j < true_beam_traj_incidentEnergies->size(); j++) {
        int slice = (true_beam_traj_slice_index)->at(j);
        //if (slice>fSliceCut) continue;
        //std::cout<<true_beam_traj_slice_z->at(j)<<std::endl;
        if (true_beam_traj_slice_z->at(j)>220.0 || true_beam_traj_slice_z->at(j)<30)  continue;
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
interactingKERW=interactingKERW-currentdEdx*reco_beam_TrkPitch_SCE->at(index)*dEdxshift;
//if(index<reco_beam_calibrated_dEdX_SCE->size()-1) reco_beam_calibrated_incidentEnergies->push_back(interactingKE);
}
else{
 
currentdEdx=reco_beam_calibrated_dEdX_SCE->at(index);

interactingKE=interactingKE-currentdEdx*reco_beam_TrkPitch_SCE->at(index);
interactingKERW=interactingKERW-currentdEdx*reco_beam_TrkPitch_SCE->at(index)*dEdxshift;
//if(index<reco_beam_calibrated_dEdX_SCE->size()-1) reco_beam_calibrated_incidentEnergies->push_back(interactingKE);

}
if(index<reco_beam_calibrated_dEdX_SCE->size()-1){
if (reco_beam_calo_Z->at(index)<30.0+startShift || reco_beam_calo_Z->at(index)>220.00+endShift) continue;
if (reco_beam_calo_Z->at(index)>30.0 && reco_beam_calo_Z->at(index)<220.0) h1d_reco_incidentE->Fill(interactingKE);
if(true_beam_traj_incidentEnergies->size()){
if (index<(true_beam_slices->size()-1) && true_beam_PDG==321 && true_beam_ID==reco_beam_true_byE_ID && reco_beam_calo_Z->at(index)<220.0) responseIncident.Fill(interactingKERW,true_beam_traj_incidentEnergies->at(index), tot_weight);

}
if(!true_beam_traj_incidentEnergies->size()) responseIncident.Fake(interactingKERW, tot_weight);
if ((index>=(true_beam_traj_slice_index->size()-1) || true_beam_PDG!=321 || true_beam_ID!=reco_beam_true_byE_ID) && reco_beam_calo_Z->at(index)<220.0) responseIncident.Fake(interactingKERW, tot_weight);
}



}

}

  }
  if (selection_ID==1 && reco_beam_endZ>30.0+startShift && reco_beam_endZ<220.00+endShift){ 


  h1d_reco_interactingE->Fill(reco_beam_calibrated_interactingEnergy);
  if( true_beam_endZ>30 && true_beam_endZ<220.0 && true_beam_endProcess->find("kaon+Inelastic")!=std::string::npos && reco_beam_true_byE_ID==true_beam_ID /*&& reco_beam_daughter_sameID==0*/){ response.Fill(interactingKERW, true_beam_traj_interactingEnergy, tot_weight);

  }
  else{
    response.Fake(interactingKERW, tot_weight);


   } 
   }


  else{
  if (true_beam_endZ>30 && true_beam_endZ<220.0 && true_beam_endProcess->find("kaon+Inelastic")!=std::string::npos) response.Miss(true_beam_traj_interactingEnergy, tot_weight);


   }



   }

std::cout<<"Unfolding MC"<<std::endl;


std::cout<<","<<h1d_reco_interactingE->GetEntries()<<std::endl;
TH1D* h1d_recoCopy_interactingE=(TH1D*)h1d_reco_interactingE->Clone("h1d_reco_copy_interactingE");
TH1D* h1d_recoCopy_incidentE=(TH1D*)h1d_reco_incidentE->Clone("h1d_reco_copy_incidentE");
RooUnfoldBayes unfoldCorr(&response, h1d_recoCopy_interactingE,4);
//unfoldCorr.IncludeSystematics();
TH1D* hRecoFullCorr=(TH1D*) unfoldCorr.Hreco(/*RooUnfold::kCovToy*/);
RooUnfoldBayes unfoldInc(&responseIncident, h1d_recoCopy_incidentE,4);
//unfoldInc.IncludeSystematics();
//RooUnfold::ErrorTreatment withError=3;
TH1D* hRecoInc=(TH1D*)unfoldInc.Hreco(/*RooUnfold::kCovToy*/);
hRecoInc->SetName("fullCorrIncHist");

//RooUnfoldBayes unfoldCheat(&response, h1d_reco_cheat_interactingE,4);
//TH1D* hRecoFullCheat=(TH1D*) unfoldCheat.Hreco();
//TMatrixD* hMatrix=(TH1D*)unfoldCorr.Hreco();
hRecoFullCorr->SetName("fullCorrRecoHist");
TCanvas c1=TCanvas();
//h1d_recoCopy_interactingE->Draw("HIST");

TFile f("kaonTreeSelData.root");
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

   for (Long64_t jentry=0; jentry<nentries_data;jentry++) {
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
std::cout<<"Response Stat:"<<R->Integral()<<","<<RInc->Integral()<<std::endl;
RooUnfoldBayes unfoldData(&response, h1d_recoCopy_interactingE_data,4);
RooUnfoldBayes unfoldIncData(&responseIncident, h1d_recoCopy_incidentE_data,4);
unfoldData.Overflow();
unfoldIncData.Overflow();
TH1D* hRecoFullData=(TH1D*) unfoldData.Hreco(/*RooUnfold::kCovToy*/);
TH1D* hRecoIncData=(TH1D*) unfoldIncData.Hreco(/*RooUnfold::kCovToy*/);


       TH1D * xsec_unfoldMC = (TH1D*)hRecoFullCorr->Clone("crossSectionUnfoldMC");

     xsec_unfoldMC->Divide(hRecoInc);
      for (int i = 1; i <= xsec_unfoldMC->GetNbinsX(); ++i) {
        xsec_unfoldMC->SetBinContent(i, -1.*log(1. - xsec_unfoldMC->GetBinContent(i)));


        //std::cout<<h1d_true_thinslice_incidentE_scaled->GetBinContent(i)<<','<<hRecoFullData->GetBinContent(i)<<std::endl;
        //if (hRecoFullData->GetBinContent(i)>h1d_true_thinslice_incidentE_scaled->GetBinContent(i)) xsec_data->SetBinContent(i,0);

      }
      xsec_unfoldMC->Scale(1.E27/ (0.4979 * 1.4 * 6.022E23 / 39.948 ));


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
        std::cout<<xsec_unfoldData->GetBinContent(i)<<std::endl;
         tmpVec.push_back(xsec_unfoldData->GetBinContent(i));
        h1d_DataXSec->Fill(xsec_unfoldData->GetBinCenter(i),xsec_unfoldData->GetBinContent(i));
        h1d_DataInc->Fill(xsec_unfoldData->GetBinCenter(i),hRecoIncData->GetBinContent(i));
        h1d_DataInt->Fill(xsec_unfoldData->GetBinCenter(i),hRecoFullData->GetBinContent(i));
        temp_hist->Fill(xsec_unfoldData->GetBinCenter(i),xsec_unfoldData->GetBinContent(i));
        
      }
      for (int i = 1; i <= xsec_unfoldMC->GetNbinsX(); ++i) {
        xsec_unfoldMCVec->push_back(xsec_unfoldMC->GetBinContent(i));
tmpVecMC.push_back(xsec_unfoldData->GetBinContent(i));
        h1d_MCXSec->Fill(xsec_unfoldMC->GetBinCenter(i),xsec_unfoldMC->GetBinContent(i));
        temp_histMC->Fill(xsec_unfoldMC->GetBinCenter(i),xsec_unfoldMC->GetBinContent(i));
        h1d_MCInc->Fill(xsec_unfoldMC->GetBinCenter(i),hRecoInc->GetBinContent(i));
        h1d_MCInt->Fill(xsec_unfoldMC->GetBinCenter(i),hRecoFullCorr->GetBinContent(i));


      }
 throwHists.push_back(temp_hist);
 throwHistsMC.push_back(temp_histMC);  






 

std::cout<<temp_hist->GetBinContent(4)<<std::endl;

 xSecMultisim.push_back(xsec_unfoldData->Interpolate(5000));
 multisimdEdXShift.push_back(dEdxshift);
 multisimBeamShift.push_back(x_val);
 multisimDataInt.push_back(hRecoFullData->Interpolate(5000));
 multisimDataInc.push_back(hRecoIncData->Interpolate(5000));

TH1D* tempHist=(TH1D*)xsec_unfoldData->Clone(Form("toy_%d",universe));
TH1D* tempHistMC=(TH1D*)xsec_unfoldMC->Clone(Form("toyMC_%d",universe));
xsec_unfoldDataVecHist.push_back(*tempHist);
xsec_unfoldMCVecHist.push_back(*tempHistMC);


xsec_unfoldDataVecVec.push_back(tmpVec);
xsec_unfoldMCVecVec.push_back(tmpVecMC);
fWrite->cd();
xsec_unfoldData->Write();
tree.Fill();


}
std::cout<<time<<std::endl;
TH2D* covMat=new TH2D("covMat","covMat",totBins, edges,totBins, edges);


TH2D* covMatMC=new TH2D("covMatMC","covMatMC",totBins, edges,totBins, edges);


for(long unsigned int toy=0; toy<xSecMultisim.size(); toy++){
    auto dataVec=xsec_unfoldDataVecVec.at(toy);
    auto mcVec=xsec_unfoldMCVecVec.at(toy);
    TH1D tmpHist=xsec_unfoldDataVecHist.at(toy);
    TH1D tmpHistMC=xsec_unfoldMCVecHist.at(toy);

   // std::cout<<tmpHist->GetBinContent(5)<<std::endl;
  for (Int_t bin=0; bin<totBins; ++bin){
    for (Int_t bin2=0; bin2<totBins; ++bin2){
    int binNum=(bin)*totBins+(bin2);
    double diffBin1MC=tmpHistMC.GetBinContent(bin+1)-h1d_MCXSec->GetBinContent(bin+1);
    double diffBin2MC=tmpHistMC.GetBinContent(bin2+1)-h1d_MCXSec->GetBinContent(bin2+1);
    //if(bin>40) std::cout<<diffBin1<<std::endl;
    covMatMC->SetBinContent(bin+1,bin2+1,covMatMC->GetBinContent(bin+1,bin2+1)+((diffBin1MC*diffBin2MC)/(xSecMultisim.size()-1)));




    double diffBin1=tmpHist.GetBinContent(bin+1)-h1d_DataXSec->GetBinContent(bin+1);
    double diffBin2=tmpHist.GetBinContent(bin2+1)-h1d_DataXSec->GetBinContent(bin2+1);
   if(bin==4 && bin2==5) std::cout<<diffBin1<<','<<diffBin2<<','<<dataVec.at(bin)<<","<<dataVec.at(bin2)<<','<<h1d_DataXSec->GetBinContent(bin+1)<<','<<h1d_DataXSec->GetBinContent(bin2+1)<<","<<bin<<','<<bin2<<std::endl;
    covMat->SetBinContent(bin+1,bin2+1,covMat->GetBinContent(bin+1,bin2+1)+((diffBin1*diffBin2)/(xSecMultisim.size()-1)));
  }
  }
}




fWrite->cd();
tree.Write();
h1d_MCXSec->Write();
h1d_DataXSec->Write();

h1d_MCInc->Write();
h1d_DataInc->Write();

h1d_MCInt->Write();
h1d_DataInt->Write();

covMat->Write();
covMatMC->Write();
h1d_reco_incidentE->Write();
h1d_reco_incidentE_data->Write();
h1d_reco_interactingE->Write();
h1d_reco_interactingE_data->Write();

fWrite->Write();
fWrite->Close();
    std::cout<<time<<std::endl;
}
