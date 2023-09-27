#define BeamSelectionTreeStandardGenOff_cxx
#include "BeamSelectionTreeStandardGenOff.h"
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include "TMath.h"
#include "TProfile.h"
#include "TVector3.h"
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



void BeamSelectionTreeStandardGenOff::Loop(int PDG)
{
//gSystem->Load("RooUnfold/libRooUnfold");
//   In a ROOT session, you can do:
//      root> .L BeamSelectionTreeStandardGenOff.C
//      root> BeamSelectionTreeStandardGenOff t
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
// METHOD2: replace 
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
  if (fChain == 0) return;





   Long64_t nentries = fChain->GetEntriesFast();
//    nentries=100;
std::string fileString="kaon";
if(PDG==211) fileString="pion";
if(PDG==2212) fileString="proton";
TFile *fout = new TFile(Form("/dune/data/users/rdiurba/offStandard6GeV%sTreeSel.root",fileString.c_str()),"RECREATE");
TTree tree("ana","anaTree");
    int selection_ID=0; int sample_ID=0;
double beam_inst_KE=0;
double reco_beam_dCos=-10, reco_beam_dirCos=-10, reco_beam_zEndPointCos=-10;
double reco_beam_cosXZ=-10, reco_beam_cosYZ=-10;
double reco_beam_calibrated_interactingEnergy=-999;
std::vector<double>* reco_beam_calibrated_incidentEnergies=0x0;
std::vector<int>* reco_daughter_true_byHits_ID=0x0;
std::vector<std::string>* reco_daughter_PFP_true_byHits_process=0x0;
std::vector<double>* true_beam_traj_incidentEnergies=0x0;
std::vector<double>* true_beam_traj_slice_z=0x0;
std::vector<int>* true_beam_traj_slice_index=0x0;
double true_beam_traj_interactingEnergy=-999;
double true_beam_mass=493.677;
int reco_beam_daughter_sameID=0;
double true_beam_len=-999;
double reco_beam_diffXY=-999;
tree.Branch("event",&event);
tree.Branch("subrun",&subrun);
tree.Branch("run",&run);
tree.Branch("true_beam_len",&true_beam_len,"true_beam_len/D");
tree.Branch("reco_beam_true_byE_ID",&reco_beam_true_byE_ID,"reco_beam_true_byE_ID/I");
tree.Branch("reco_beam_true_byE_PDG",&reco_beam_true_byE_PDG,"reco_beam_true_byE_PDG/I");
tree.Branch("reco_beam_true_byHits_ID",&reco_beam_true_byHits_ID,"reco_beam_true_byHits_ID/I");
tree.Branch("true_beam_ID",&true_beam_ID,"true_beam_ID/I");
tree.Branch("true_beam_PDG",&true_beam_PDG,"true_beam_PDG/I");
tree.Branch("true_beam_endP",&true_beam_endP,"true_beam_endP/D");
tree.Branch("true_beam_startP",&true_beam_startP,"true_beam_startP/D");
tree.Branch("true_beam_mass",&true_beam_mass,"true_beam_mass/D");
tree.Branch("selection_ID",&selection_ID,"selection_ID/I");
tree.Branch("sample_ID",&sample_ID,"sample_ID/I");

tree.Branch("reco_beam_daughter_sameID",&reco_beam_daughter_sameID,"reco_beam_daughter_sameID/I");
tree.Branch("true_beam_interactingEnergy",&true_beam_interactingEnergy);
tree.Branch("reco_beam_interactingEnergy",&reco_beam_interactingEnergy);


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


tree.Branch("true_beam_startDirZ",&true_beam_startDirZ, "true_beam_startDirZ/D");
tree.Branch("true_beam_startDirY",&true_beam_startDirY, "true_beam_startDirY/D");
tree.Branch("true_beam_startDirX",&true_beam_startDirX, "true_beam_startDirX/D");






tree.Branch("reco_beam_endZ",&reco_beam_endZ, "reco_beam_endZ/D");
tree.Branch("true_beam_endZ",&true_beam_endZ, "true_beam_endZ/D");
tree.Branch("reco_beam_len",&reco_beam_len, "reco_beam_len/D");
tree.Branch("reco_beam_alt_len",&reco_beam_alt_len, "reco_beam_alt_len/D");
tree.Branch("reco_beam_true_byE_endProcess",&reco_beam_true_byE_endProcess);
tree.Branch("true_beam_endProcess",&true_beam_endProcess);
tree.Branch("reco_beam_dEdX_NoSCE",&reco_beam_dEdX_NoSCE);
tree.Branch("reco_beam_calibrated_dEdX_NoSCE",&reco_beam_calibrated_dEdX_NoSCE);
tree.Branch("reco_beam_dEdX_SCE",&reco_beam_dEdX_SCE);
tree.Branch("reco_beam_calibrated_dEdX_SCE",&reco_beam_calibrated_dEdX_SCE);
tree.Branch("true_beam_incidentEnergies",&true_beam_incidentEnergies);
tree.Branch("reco_beam_incidentEnergies",&reco_beam_incidentEnergies);
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
tree.Branch("reco_beam_diffXY",&reco_beam_diffXY,"reco_beam_diffXY/D");

tree.Branch("reco_beam_dCos",&reco_beam_dCos,"reco_beam_dCos/D");
tree.Branch("reco_beam_dirCos",&reco_beam_dirCos,"reco_beam_dirCos/D");
tree.Branch("reco_beam_zEndPointCos",&reco_beam_zEndPointCos,"reco_beam_zEndPointCos/D");
tree.Branch("reco_beam_cosXZ",&reco_beam_cosXZ,"reco_beam_cosXZ/D");
tree.Branch("reco_beam_cosYZ",&reco_beam_cosYZ,"reco_beam_cosYZ/D");
tree.Branch("true_beam_traj_Z",&true_beam_traj_Z);
tree.Branch("true_beam_traj_KE",&true_beam_traj_KE);
tree.Branch("reco_beam_type",&reco_beam_type,"reco_beam_type/I");
tree.Branch("reco_daughter_allTrack_len",&reco_daughter_allTrack_len);
tree.Branch("reco_daughter_allTrack_ID",&reco_daughter_allTrack_ID);
tree.Branch("reco_daughter_PFP_true_byHits_PDG",&reco_daughter_PFP_true_byHits_PDG);
tree.Branch("reco_daughter_PFP_true_byHits_ID",&reco_daughter_PFP_true_byHits_ID);
int reco_beam_nDaughters=-1; int reco_daughter_PFP_true_nDaughters=-1;
tree.Branch("reco_beam_nDaughters",&reco_beam_nDaughters);
tree.Branch("reco_daughter_PFP_true_nDaughters",&reco_daughter_PFP_true_nDaughters);
//tree.Branch("reco_daughter_PFP_true_byHits_ID",&reco_daughter_PFP_true_byHits_ID);
tree.Branch("reco_daughter_PFP_true_byHits_process",&reco_daughter_PFP_true_byHits_process);
tree.Branch("reco_beam_calibrated_interactingEnergy",&reco_beam_calibrated_interactingEnergy);
tree.Branch("reco_beam_calibrated_incidentEnergies",&reco_beam_calibrated_incidentEnergies);

tree.Branch("true_beam_nElasticScatters",&true_beam_nElasticScatters);
tree.Branch("true_beam_elastic_costheta",&true_beam_elastic_costheta);
tree.Branch("true_beam_elastic_X",&true_beam_elastic_X);
tree.Branch("true_beam_elastic_Y",&true_beam_elastic_Y);
tree.Branch("true_beam_elastic_Z",&true_beam_elastic_Z);

std::vector<double>* reco_daughter_allTrack_dirCos=0x0;
tree.Branch("reco_daughter_allTrack_dirCos",&reco_daughter_allTrack_dirCos);
tree.Branch("reco_daughter_allTrack_startX",&reco_daughter_allTrack_startX);
tree.Branch("reco_daughter_allTrack_startY",&reco_daughter_allTrack_startY);
tree.Branch("reco_daughter_allTrack_startZ",&reco_daughter_allTrack_startZ);

tree.Branch("reco_daughter_allTrack_endX",&reco_daughter_allTrack_endX);
tree.Branch("reco_daughter_allTrack_endY",&reco_daughter_allTrack_endY);
tree.Branch("reco_daughter_allTrack_endZ",&reco_daughter_allTrack_endZ);


tree.Branch("reco_daughter_PFP_true_byHits_startX",&reco_daughter_PFP_true_byHits_startX);
tree.Branch("reco_daughter_PFP_true_byHits_startY",&reco_daughter_PFP_true_byHits_startY);
tree.Branch("reco_daughter_PFP_true_byHits_startZ",&reco_daughter_PFP_true_byHits_startZ);

tree.Branch("reco_daughter_PFP_true_byHits_endX",&reco_daughter_PFP_true_byHits_endX);
tree.Branch("reco_daughter_PFP_true_byHits_endY",&reco_daughter_PFP_true_byHits_endY);
tree.Branch("reco_daughter_PFP_true_byHits_endZ",&reco_daughter_PFP_true_byHits_endZ);
tree.Branch("true_beam_traj_incidentEnergies",&true_beam_traj_incidentEnergies);
tree.Branch("true_beam_traj_interactingEnergy",&true_beam_traj_interactingEnergy);
tree.Branch("true_beam_traj_slice_z",&true_beam_traj_slice_z);
tree.Branch("true_beam_traj_slice_index",&true_beam_traj_slice_index);
//tree.Branch("g4rw_p1",&g4rw_p1);
tree.Branch("g4rw_full_grid_kplus_weights",&g4rw_full_grid_kplus_weights);
tree.Branch("g4rw_full_grid_kplus_coeffs",&g4rw_full_grid_kplus_coeffs);   
Long64_t nbytes = 0, nb = 0;
   //nentries=10000;
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
   beam_inst_P=beam_inst_P*6.f;
   TVector3 vec(dX, dY, dZ);
   true_beam_traj_incidentEnergies->clear();
   reco_beam_calibrated_incidentEnergies->clear();
   true_beam_traj_slice_z->clear();
   true_beam_traj_slice_index->clear();
   reco_beam_calibrated_interactingEnergy=-999;
   TVector3 vhat=vec.Unit();
    reco_beam_zEndPointCos=vhat(2);
    trklen=TMath::Sqrt(dX*dX+dY*dY+dZ*dZ);
     
    reco_beam_nDaughters=reco_daughter_allTrack_ID->size();
   reco_daughter_PFP_true_nDaughters=reco_daughter_PFP_true_byHits_ID->size();
     reco_beam_daughter_sameID=0;
     true_beam_len=-99999;
     for(int daughterID:*reco_daughter_PFP_true_byHits_ID){

     if (daughterID==reco_beam_true_byE_ID) reco_beam_daughter_sameID=1;

      }
    // trklen=reco_beam_len;
    beam_inst_KE=1000.f*(TMath::Sqrt(beam_inst_P*beam_inst_P+0.493*0.493)-0.493);
//double totalEvent=0;
std::string inelProcess="kaon+Inelastic";
   if(PDG==211){
     if(true_beam_PDG!=13 && true_beam_PDG!=211) continue;
beam_inst_KE=1000.f*(TMath::Sqrt(beam_inst_P*beam_inst_P+0.1396*0.1396)-0.1396);
   inelProcess="pi+Inelastic";
   }
   else if(PDG==2212){
   if(true_beam_PDG!=2212) continue;
beam_inst_KE=1000.f*(TMath::Sqrt(beam_inst_P*beam_inst_P+0.9383*0.9383)-0.9383);
inelProcess="protonInelastic";
 
 } 
   else{
     if(true_beam_PDG!=321) continue;
   }
    reco_beam_diffXY=TMath::Sqrt(TMath::Power(beam_inst_X-reco_beam_startX,2)+TMath::Power(beam_inst_Y-reco_beam_startY,2));
    reco_beam_dCos=beam_inst_dirX*reco_beam_trackDirX+beam_inst_dirY*reco_beam_trackDirY+beam_inst_dirZ*reco_beam_trackDirZ;
    reco_beam_dirCos=reco_beam_trackEndDirX*reco_beam_trackDirX+reco_beam_trackEndDirY*reco_beam_trackDirY+reco_beam_trackEndDirZ*reco_beam_trackDirZ;
       reco_beam_cosXZ=(reco_beam_endZ-reco_beam_startZ)/(reco_beam_endX-reco_beam_startX);
       reco_beam_cosYZ=(reco_beam_endZ-reco_beam_startZ)/(reco_beam_endY-reco_beam_startY);
   
   if (beam_inst_P<4 || beam_inst_P>8) continue;
   if(reco_beam_len<-800 || reco_beam_calo_wire->size()<1 ) selection_ID=4;
  // if(reco_beam_len<-800 || reco_beam_calo_wire->size()<1) continue;
   else if (abs(reco_beam_dCos)<0.997 || reco_beam_dCos>1.0 || abs(beam_inst_Y-reco_beam_startY+0.59)>0.9 || abs(beam_inst_Z-reco_beam_startZ+0.3)>0.9 || abs(beam_inst_X-reco_beam_startX-1.65)>0.75 || abs(reco_beam_diffXY-1.78)>0.75) selection_ID=3;
   else if (reco_beam_endZ>220.0) selection_ID=2;
   else if (reco_beam_zEndPointCos<.99 && reco_beam_zEndPointCos>0.87) selection_ID=1;
   else{ selection_ID=1;

     }

   if (true_beam_endZ<0) sample_ID=3;
   else if(true_beam_endZ>226.0) sample_ID=2;
   else if(true_beam_endProcess->find("kaon+Inelastic")!=std::string::npos ){
   sample_ID=1;


   }

   else sample_ID=4;

for(long unsigned int i=0; i<reco_daughter_allTrack_startX->size();++i){
    if (abs(reco_daughter_allTrack_startX->at(i))>700) continue;
     double trklen=reco_beam_alt_len; 
    double dZ_daughter=(reco_daughter_allTrack_endZ->at(i)-reco_daughter_allTrack_startZ->at(i));
    double dX_daughter=(reco_daughter_allTrack_endX->at(i)-reco_daughter_allTrack_startX->at(i));
    double dY_daughter=(reco_daughter_allTrack_endY->at(i)-reco_daughter_allTrack_startY->at(i));
   TVector3 vecDaugh(dX_daughter, dY_daughter, dZ_daughter);
   TVector3 vhatDaugh=vec.Unit();
    double dirZ_daughter=vhatDaugh(2);
    double dirX_daughter=vhatDaugh(0);
    double dirY_daughter=vhatDaugh(1);
    reco_daughter_allTrack_dirCos->push_back(reco_beam_trackEndDirX*dirX_daughter+reco_beam_trackEndDirY*dirY_daughter+reco_beam_trackEndDirZ*dirZ_daughter);

}

double initialKE=-999;
double interactingKE=-999;


if (reco_beam_calibrated_dEdX_SCE->size()>0 && reco_beam_calo_wire->size()>0 ){
initialKE=beam_inst_KE;
interactingKE=beam_inst_KE;

double currentdEdx=BetheBloch(beam_inst_KE, 493.67);
std::vector<double> copyCalibrated=*reco_beam_calibrated_dEdX_SCE;
for(long unsigned int index=0; index<reco_beam_calibrated_dEdX_SCE->size(); ++index){

if (reco_beam_calibrated_dEdX_SCE->at(index)>5 || reco_beam_calibrated_dEdX_SCE->at(index)<0){
interactingKE=interactingKE-currentdEdx*reco_beam_TrkPitch_SCE->at(index);
if(index<reco_beam_calibrated_dEdX_SCE->size()-1) reco_beam_calibrated_incidentEnergies->push_back(interactingKE);
}
else{
 
currentdEdx=reco_beam_calibrated_dEdX_SCE->at(index);

interactingKE=interactingKE-currentdEdx*reco_beam_TrkPitch_SCE->at(index);
if(index<reco_beam_calibrated_dEdX_SCE->size()-1) reco_beam_calibrated_incidentEnergies->push_back(interactingKE);

}
}
reco_beam_calibrated_interactingEnergy=interactingKE;
//std::cout<<interactingKE<<std::endl;
}

/*if(true_beam_slices->size()){  

std::cout<<"Start point: "<<true_beam_startZ<<","<<true_beam_slices->at(0)<<std::endl;
std::cout<<"End point: "<<true_beam_endZ<<','<<true_beam_slices->at(true_beam_slices->size()-1)<<std::endl;

}*/


    double fTrajZStart=-.49375;
    double fPitch=0.47974;
    double next_slice_z = fTrajZStart;
    int next_slice_num = 0;
    true_beam_len=0;
    for (size_t j = 1; j < true_beam_traj_Z->size()-1 ; ++j) {
   double z_traj=true_beam_traj_Z->at(j-1);
   double y_traj=true_beam_traj_Y->at(j-1);
     double x_traj=true_beam_traj_X->at(j-1);
  double diffZ=true_beam_traj_Z->at(j)-z_traj;
  double diffY=true_beam_traj_Y->at(j)-y_traj;
  double diffX=true_beam_traj_X->at(j)-x_traj;
   if (z_traj>fTrajZStart && z_traj<700){
        true_beam_len=true_beam_len+TMath::Sqrt(diffZ*diffZ+diffY*diffY+diffX*diffX);



}
    }
    for (size_t j = 1; j < true_beam_traj_Z->size()-1 ; ++j) {
      double z = true_beam_traj_Z->at(j);
      double x = true_beam_traj_X->at(j);
      double y = true_beam_traj_Y->at(j);
      double ke = true_beam_traj_KE->at(j);

      if (z < fTrajZStart) {
        continue;
      }

      if (z >= next_slice_z) {
        double temp_z = true_beam_traj_Z->at(j-1);
        double temp_y = true_beam_traj_Y->at(j-1);
        double temp_x = true_beam_traj_X->at(j-1);
        double temp_e = true_beam_traj_KE->at(j-1);
        
        while (next_slice_z < z && next_slice_num < 100000) {
//          std::cout<<"Current Slice: "<<next_slice_z<<std::endl;
          double sub_z = next_slice_z - temp_z;
          double delta_e = true_beam_traj_KE->at(j-1) - ke;
          double delta_z = z - true_beam_traj_Z->at(j-1);
//     std::cout<<next_slice_z<<','<<temp_z<<std::endl;
         
          temp_e -= (sub_z/delta_z)*delta_e;
          true_beam_traj_incidentEnergies->push_back(temp_e);
          true_beam_traj_slice_index->push_back(next_slice_num);
          true_beam_traj_slice_z->push_back(next_slice_z);
          temp_z = next_slice_z;
          //temp_x= next_slice_x;
          //temp_y= next_slice_y;
          next_slice_z += fPitch;
//          std::cout<<temp_e<<std::endl;
          ++next_slice_num;
         //std::cout<<temp_e<<','<<next_slice_z<<','<<true_beam_interactingEnergy<<','<<true_beam_endZ<<std::endl;
        }
      }
    }
//if(selection_ID==1 && true_beam_ID==reco_beam_true_byE_ID)  std::cout<<true_beam_len<<","<<true_beam_endZ<<','<<reco_beam_alt_len<<','<<reco_beam_endZ<<std::endl;    
   
   true_beam_traj_interactingEnergy= sqrt(true_beam_endP*true_beam_endP*1.e6 + 493.677*493.677) - 493.677;
  // std::cout<<true_beam_traj_interactingEnergy<<std::endl;
   tree.Fill();

}
   tree.Write();
}
