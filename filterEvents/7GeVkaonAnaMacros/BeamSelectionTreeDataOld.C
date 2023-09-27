#define BeamSelectionTreeDataOld_cxx
#include "BeamSelectionTreeDataOld.h"
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include "TMath.h"
#include "TProfile.h"
#include <iostream>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <vector>
#include <numeric>

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

void BeamSelectionTreeDataOld::Loop(int PDG)
{
//gSystem->Load("RooUnfold/libRooUnfold");
//   In a ROOT session, you can do:
//      root> .L BeamSelectionTreeDataOld.C
//      root> BeamSelectionTreeDataOld t
//      root> t.GetEntry(12); // Fill t DataOld members with entry number 12
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
// METHOD2: replace 
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
  if (fChain == 0) return;





   Long64_t nentries = fChain->GetEntriesFast();
//    nentries=100;
std::string fileString="kaon";
if(PDG==211) fileString="pion";
if(PDG==2212) fileString="proton";
TFile *fout = new TFile(Form("%sTreeSelDataOld.root",fileString.c_str()),"RECREATE");
TTree tree("ana","anaTree");
    int selection_ID=0;
double beam_inst_KE=0;
double reco_beam_dCos=-10, reco_beam_dirCos=-10;
double reco_beam_cosXZ=-10, reco_beam_cosYZ=-10;
double reco_beam_calibrated_interactingEnergy=-999;
std::vector<double>* reco_beam_calibrated_incidentEnergies=0x0;
double reco_beam_diffXY=-999;

tree.Branch("event",&event,"event/I");
tree.Branch("run",&run,"run/I");
tree.Branch("reco_beam_diffXY",&reco_beam_diffXY,"reco_beam_diffXY/D");
tree.Branch("reco_beam_true_byE_ID",&reco_beam_true_byE_ID,"reco_beam_true_byE_ID/I");
tree.Branch("reco_beam_true_by_PDG",&reco_beam_true_byE_PDG,"reco_beam_true_byE_PDG/I");
tree.Branch("true_beam_ID",&true_beam_ID,"true_beam_ID/I");
tree.Branch("true_beam_PDG",&true_beam_PDG,"true_beam_PDG/I");


tree.Branch("selection_ID",&selection_ID,"selection_ID/I");

tree.Branch("true_beam_interactingEnergy",&true_beam_interactingEnergy);



tree.Branch("reco_beam_startZ",&reco_beam_startZ, "reco_beam_startZ/D");
tree.Branch("true_beam_startZ",&true_beam_startZ, "true_beam_startZ/D");
tree.Branch("reco_beam_startY",&reco_beam_startY, "reco_beam_startY/D");
tree.Branch("true_beam_startY",&true_beam_startY, "true_beam_startY/D");
tree.Branch("reco_beam_startX",&reco_beam_startX, "reco_beam_startX/D");
tree.Branch("true_beam_startX",&true_beam_startX, "true_beam_startX/D");


tree.Branch("reco_beam_endY",&reco_beam_endY, "reco_beam_endY/D");
tree.Branch("true_beam_endY",&true_beam_endY, "true_beam_endY/D");
tree.Branch("reco_beam_endX",&reco_beam_endX, "reco_beam_endX/D");
tree.Branch("true_beam_endX",&true_beam_endX, "true_beam_endX/D");

tree.Branch("reco_beam_trackDirZ",&reco_beam_trackDirZ, "reco_beam_trackDirZ/D");
tree.Branch("reco_beam_trackDirY",&reco_beam_trackDirY, "reco_beam_trackDirY/D");
tree.Branch("reco_beam_trackDirX",&reco_beam_trackDirX, "reco_beam_trackDirX/D");

tree.Branch("reco_beam_trackEndDirZ",&reco_beam_trackEndDirZ, "reco_beam_trackEndDirZ/D");
tree.Branch("reco_beam_trackEndDirY",&reco_beam_trackEndDirY, "reco_beam_trackEndDirY/D");
tree.Branch("reco_beam_trackEndDirX",&reco_beam_trackEndDirX, "reco_beam_trackEndDirX/D");


tree.Branch("reco_beam_endZ",&reco_beam_endZ, "reco_beam_endZ/D");
tree.Branch("true_beam_endZ",&true_beam_endZ, "true_beam_endZ/D");
tree.Branch("reco_beam_len",&reco_beam_len, "reco_beam_len/D");
tree.Branch("reco_beam_alt_len",&reco_beam_alt_len, "reco_beam_alt_len/D");
tree.Branch("reco_beam_true_byE_endProcess",&reco_beam_true_byE_endProcess);
tree.Branch("true_beam_endProcess",&true_beam_endProcess);
tree.Branch("reco_beam_calibrated_dEdX_SCE",&reco_beam_calibrated_dEdX_SCE);
tree.Branch("reco_beam_dEdX_SCE",&reco_beam_dEdX_SCE);
tree.Branch("true_beam_indicentEnergies",&true_beam_incidentEnergies);
tree.Branch("reco_beam_TrkPitch_SCE",&reco_beam_TrkPitch_SCE);
tree.Branch("reco_beam_calo_wire",&reco_beam_calo_wire);
tree.Branch("reco_beam_calo_wire_z",&reco_beam_calo_wire_z);
tree.Branch("true_beam_slices",&true_beam_slices);
tree.Branch("true_beam_slices_dE",&true_beam_slices_deltaE);
tree.Branch("beam_inst_KE",&beam_inst_KE,"beam_inst_KE/D");
tree.Branch("beam_inst_P",&beam_inst_P,"beam_inst_P/D");

tree.Branch("beam_inst_X",&beam_inst_X,"beam_inst_X/D");
tree.Branch("beam_inst_Y",&beam_inst_Y,"beam_inst_Y/D");
tree.Branch("beam_inst_Z",&beam_inst_Z,"beam_inst_Z/D");
tree.Branch("beam_inst_dirX",&beam_inst_dirX,"beam_inst_dirX/D");
tree.Branch("beam_inst_dirY",&beam_inst_dirY,"beam_inst_dirY/D");
tree.Branch("beam_inst_dirZ",&beam_inst_dirZ,"beam_inst_dirZ/D");


tree.Branch("reco_beam_dCos",&reco_beam_dCos,"reco_beam_dCos/D");
tree.Branch("reco_beam_dirCos",&reco_beam_dirCos,"reco_beam_dirCos/D");

tree.Branch("reco_beam_cosXZ",&reco_beam_cosXZ,"reco_beam_cosXZ/D");
tree.Branch("reco_beam_cosYZ",&reco_beam_cosYZ,"reco_beam_cosYZ/D");

tree.Branch("reco_beam_type",&reco_beam_type,"reco_beam_type/I");
tree.Branch("reco_daughter_allTrack_len",&reco_daughter_allTrack_len);
tree.Branch("reco_daughter_allTrack_ID",&reco_daughter_allTrack_ID);
tree.Branch("reco_daughter_PFP_true_byHits_PDG",&reco_daughter_PFP_true_byHits_PDG);
tree.Branch("reco_daughter_PFP_true_byHits_ID",&reco_daughter_PFP_true_byHits_ID);
int reco_beam_nDaughters=-1; int reco_daughter_PFP_true_nDaughters=-1;
tree.Branch("reco_beam_nDaughters",&reco_beam_nDaughters);
tree.Branch("reco_daughter_PFP_true_nDaughters",&reco_daughter_PFP_true_nDaughters);

tree.Branch("reco_beam_calibrated_interactingEnergy",&reco_beam_calibrated_interactingEnergy);
tree.Branch("reco_beam_calibrated_incidentEnergies",&reco_beam_calibrated_incidentEnergies);
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if(jentry%1000==0) std::cout<<jentry<<"/"<<nentries<<std::endl;
      // if (Cut(ientry) < 0) continue;
     double trklen=reco_beam_alt_len; 
    double dZ=reco_beam_endZ-reco_beam_startZ;
    double dX=reco_beam_endX-reco_beam_startX;
    double dY=reco_beam_endY-reco_beam_startY;
    trklen=TMath::Sqrt(dX*dX+dY*dY+dZ*dZ);
reco_beam_nDaughters=reco_daughter_allTrack_ID->size();
reco_daughter_PFP_true_nDaughters=reco_daughter_PFP_true_byHits_ID->size();
    // trklen=reco_beam_len;
    beam_inst_KE=1000.f*(TMath::Sqrt(beam_inst_P*beam_inst_P+0.493*0.493)-0.493);
//double totalEvent=0;
std::string inelProcess="kaon+Inelastic";
   if(PDG==211){
     if(beam_inst_c0!=0 || beam_inst_c1!=0) continue;
beam_inst_KE=1000.f*(TMath::Sqrt(beam_inst_P*beam_inst_P+0.1396*0.1396)-0.1396);
   }
   else if(PDG==2212){
     if(beam_inst_c0==0 || beam_inst_c1==0) continue;
beam_inst_KE=1000.f*(TMath::Sqrt(beam_inst_P*beam_inst_P+0.9383*0.9383)-0.9383);
   } 
   else{
     if(beam_inst_c0!=1 || beam_inst_c1!=0) continue;


   }

reco_beam_diffXY=TMath::Sqrt(TMath::Power(reco_beam_startX-beam_inst_X,2)+TMath::Power(reco_beam_startY-beam_inst_Y,2));
    reco_beam_dCos=beam_inst_dirX*reco_beam_trackDirX+beam_inst_dirY*reco_beam_trackDirY+beam_inst_dirZ*reco_beam_trackDirZ;
    reco_beam_dirCos=reco_beam_trackEndDirX*reco_beam_trackDirX+reco_beam_trackEndDirY*reco_beam_trackDirY+reco_beam_trackEndDirZ*reco_beam_trackDirZ;
       reco_beam_cosXZ=(reco_beam_endZ-reco_beam_startZ)/(reco_beam_endX-reco_beam_startX);
       reco_beam_cosYZ=(reco_beam_endZ-reco_beam_startZ)/(reco_beam_endY-reco_beam_startY);

   //if(reco_beam_alt_len<0 || reco_beam_calo_wire->size()<1) selection_ID=4;
   if (beam_inst_P<4 || beam_inst_P>8) continue;
   if(reco_beam_len<-800 || reco_beam_calo_wire->size()<1) selection_ID=4;
   //if(reco_beam_len<-800 || reco_beam_calo_wire->size()<1) continue;
   else if (!(abs(beam_inst_Y-reco_beam_startY+1.91)<3.27  && abs(beam_inst_Z-reco_beam_startZ+30.7)<3.45  && abs(beam_inst_X-reco_beam_startX+4)<1 && reco_beam_dCos>0.9623 && abs(reco_beam_diffXY-4.433)<1.515)) selection_ID=3;
   else if (reco_beam_endZ>220.0) selection_ID=2;
   else if (reco_beam_dirCos<.99) selection_ID=1;
   else{ selection_ID=1;

     }

double initialKE=-999;
double interactingKE=-999;


if (reco_beam_dEdX_SCE->size()>0 && reco_beam_calo_wire->size()>0 ){
initialKE=beam_inst_KE;
interactingKE=beam_inst_KE;


double currentdEdx=BetheBloch(beam_inst_KE, 493.67);
std::vector<double> copyCalibrated=*reco_beam_dEdX_SCE;
for(long unsigned int index=0; index<reco_beam_dEdX_SCE->size(); ++index){

if (reco_beam_dEdX_SCE->at(index)>5 || reco_beam_dEdX_SCE->at(index)<0){
interactingKE=interactingKE-currentdEdx*reco_beam_TrkPitch_SCE->at(index);
//if(index<reco_beam_calibrated_dEdX->size()-1) reco_beam_calibrated_incidentEnergies->push_back(interactingKE);
}
else{
 
currentdEdx=reco_beam_dEdX_SCE->at(index);

interactingKE=interactingKE-currentdEdx*reco_beam_TrkPitch_SCE->at(index);
//if(index<reco_beam_calibrated_dEdX->size()-1) reco_beam_calibrated_incidentEnergies->push_back(interactingKE);

}
}
}

reco_beam_calibrated_interactingEnergy=interactingKE;


   tree.Fill();

}
   tree.Write();
}
