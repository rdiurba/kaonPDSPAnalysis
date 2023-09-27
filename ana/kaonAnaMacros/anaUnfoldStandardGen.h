//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Aug 23 08:33:33 2021 by ROOT version 6.22/06
// from TTree ana/anaTree
// found on file: kaonAnaMacros/kaonTreeSel.root
//////////////////////////////////////////////////////////

#ifndef anaUnfoldStandardGen_h
#define anaUnfoldStandardGen_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "string"
#include "vector"
#include "vector"
#include "vector"

class anaUnfoldStandardGen {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           event;
   Int_t           subrun;
   Int_t           run;
   Int_t           reco_beam_true_byE_ID;
   Int_t           reco_beam_true_byE_PDG;
   Int_t           true_beam_ID;
   Int_t           true_beam_PDG;
   Double_t        true_beam_endP;
   Double_t        true_beam_startP;
   Int_t           selection_ID;
   Int_t           sample_ID;
   Int_t           reco_beam_daughter_sameID;
   Double_t        true_beam_interactingEnergy;
   Double_t        reco_beam_interactingEnergy;
   Double_t        reco_beam_startZ;
   Double_t        true_beam_startZ;
   Double_t        reco_beam_startY;
   Double_t        true_beam_startY;
   Double_t        reco_beam_startX;
   Double_t        true_beam_startX;
   Double_t        reco_beam_endY;
   Double_t        true_beam_endY;
   Double_t        reco_beam_endX;
   Double_t        true_beam_endX;
   Double_t        reco_beam_trackDirZ;
   Double_t        reco_beam_trackDirY;
   Double_t        reco_beam_trackDirX;
   Double_t        reco_beam_trackEndDirZ;
   Double_t        reco_beam_trackEndDirY;
   Double_t        reco_beam_trackEndDirX;
   Double_t        reco_beam_endZ;
   Double_t        true_beam_endZ;
   Double_t        reco_beam_len;
   Double_t        reco_beam_alt_len;
   string          *reco_beam_true_byE_endProcess;
   string          *true_beam_endProcess;
   vector<double>  *reco_beam_calibrated_dEdX_SCE;
   vector<double>  *true_beam_incidentEnergies;
   vector<double>  *reco_beam_incidentEnergies;
   vector<double>  *reco_beam_TrkPitch_SCE;
   vector<double>  *reco_beam_calo_wire;
   vector<double>  *reco_beam_calo_wire_z;
   vector<int>     *true_beam_slices;
   vector<double>  *true_beam_slices_dE;
   Double_t        beam_inst_KE;
   Double_t        beam_inst_P;
   Double_t        beam_inst_X;
   Double_t        beam_inst_Y;
   Double_t        beam_inst_Z;
   Double_t        beam_inst_dirX;
   Double_t        beam_inst_dirY;
   Double_t        beam_inst_dirZ;
   Double_t        reco_beam_dCos;
   Double_t        reco_beam_dirCos;
   Double_t        reco_beam_zEndPointCos;
   Double_t        reco_beam_cosXZ;
   Double_t        reco_beam_cosYZ;
   vector<double>  *true_beam_traj_Z;
   vector<double>  *true_beam_traj_KE;
   Int_t           reco_beam_type;
   vector<double>  *reco_daughter_allTrack_len;
   vector<int>     *reco_daughter_allTrack_ID;
   vector<int>     *reco_daughter_PFP_true_byHits_PDG;
   vector<int>     *reco_daughter_PFP_true_byHits_ID;
   Int_t           reco_beam_nDaughters;
   Int_t           reco_daughter_PFP_true_nDaughters;
   vector<string>  *reco_daughter_PFP_true_byHits_process;
   Double_t        reco_beam_calibrated_interactingEnergy;
   vector<double>  *reco_beam_calibrated_incidentEnergies;
   Int_t           true_beam_nElasticScatters;
   vector<double>  *true_beam_elastic_costheta;
   vector<double>  *true_beam_elastic_X;
   vector<double>  *true_beam_elastic_Y;
   vector<double>  *true_beam_elastic_Z;
   vector<double>  *reco_daughter_allTrack_dirCos;
   vector<double>  *reco_daughter_allTrack_startX;
   vector<double>  *reco_daughter_allTrack_startY;
   vector<double>  *reco_daughter_allTrack_startZ;
   vector<double>  *reco_daughter_allTrack_endX;
   vector<double>  *reco_daughter_allTrack_endY;
   vector<double>  *reco_daughter_allTrack_endZ;
   vector<double>  *reco_daughter_PFP_true_byHits_startX;
   vector<double>  *reco_daughter_PFP_true_byHits_startY;
   vector<double>  *reco_daughter_PFP_true_byHits_startZ;
   vector<double>  *reco_daughter_PFP_true_byHits_endX;
   vector<double>  *reco_daughter_PFP_true_byHits_endY;
   vector<double>  *reco_daughter_PFP_true_byHits_endZ;
   vector<double>  *true_beam_traj_incidentEnergies;
   Double_t        true_beam_traj_interactingEnergy;
   vector<double>  *true_beam_traj_slice_z;
   vector<int>     *true_beam_traj_slice_index;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_subrun;   //!
   TBranch        *b_run;   //!
   TBranch        *b_reco_beam_true_byE_ID;   //!
   TBranch        *b_reco_beam_true_byE_PDG;   //!
   TBranch        *b_true_beam_ID;   //!
   TBranch        *b_true_beam_PDG;   //!
   TBranch        *b_true_beam_endP;   //!
   TBranch        *b_true_beam_startP;   //!
   TBranch        *b_selection_ID;   //!
   TBranch        *b_sample_ID;   //!
   TBranch        *b_reco_beam_daughter_sameID;   //!
   TBranch        *b_true_beam_interactingEnergy;   //!
   TBranch        *b_reco_beam_interactingEnergy;   //!
   TBranch        *b_reco_beam_startZ;   //!
   TBranch        *b_true_beam_startZ;   //!
   TBranch        *b_reco_beam_startY;   //!
   TBranch        *b_true_beam_startY;   //!
   TBranch        *b_reco_beam_startX;   //!
   TBranch        *b_true_beam_startX;   //!
   TBranch        *b_reco_beam_endY;   //!
   TBranch        *b_true_beam_endY;   //!
   TBranch        *b_reco_beam_endX;   //!
   TBranch        *b_true_beam_endX;   //!
   TBranch        *b_reco_beam_trackDirZ;   //!
   TBranch        *b_reco_beam_trackDirY;   //!
   TBranch        *b_reco_beam_trackDirX;   //!
   TBranch        *b_reco_beam_trackEndDirZ;   //!
   TBranch        *b_reco_beam_trackEndDirY;   //!
   TBranch        *b_reco_beam_trackEndDirX;   //!
   TBranch        *b_reco_beam_endZ;   //!
   TBranch        *b_true_beam_endZ;   //!
   TBranch        *b_reco_beam_len;   //!
   TBranch        *b_reco_beam_alt_len;   //!
   TBranch        *b_reco_beam_true_byE_endProcess;   //!
   TBranch        *b_true_beam_endProcess;   //!
   TBranch        *b_reco_beam_calibrated_dEdX_SCE;   //!
   TBranch        *b_true_beam_incidentEnergies;   //!
   TBranch        *b_reco_beam_incidentEnergies;   //!
   TBranch        *b_reco_beam_TrkPitch_SCE;   //!
   TBranch        *b_reco_beam_calo_wire;   //!
   TBranch        *b_reco_beam_calo_wire_z;   //!
   TBranch        *b_true_beam_slices;   //!
   TBranch        *b_true_beam_slices_dE;   //!
   TBranch        *b_beam_inst_KE;   //!
   TBranch        *b_beam_inst_P;   //!
   TBranch        *b_beam_inst_X;   //!
   TBranch        *b_beam_inst_Y;   //!
   TBranch        *b_beam_inst_Z;   //!
   TBranch        *b_beam_inst_dirX;   //!
   TBranch        *b_beam_inst_dirY;   //!
   TBranch        *b_beam_inst_dirZ;   //!
   TBranch        *b_reco_beam_dCos;   //!
   TBranch        *b_reco_beam_dirCos;   //!
   TBranch        *b_reco_beam_zEndPointCos;   //!
   TBranch        *b_reco_beam_cosXZ;   //!
   TBranch        *b_reco_beam_cosYZ;   //!
   TBranch        *b_true_beam_traj_Z;   //!
   TBranch        *b_true_beam_traj_KE;   //!
   TBranch        *b_reco_beam_type;   //!
   TBranch        *b_reco_daughter_allTrack_len;   //!
   TBranch        *b_reco_daughter_allTrack_ID;   //!
   TBranch        *b_reco_daughter_PFP_true_byHits_PDG;   //!
   TBranch        *b_reco_daughter_PFP_true_byHits_ID;   //!
   TBranch        *b_reco_beam_nDaughters;   //!
   TBranch        *b_reco_daughter_PFP_true_nDaughters;   //!
   TBranch        *b_reco_daughter_PFP_true_byHits_process;   //!
   TBranch        *b_reco_beam_calibrated_interactingEnergy;   //!
   TBranch        *b_reco_beam_calibrated_incidentEnergies;   //!
   TBranch        *b_true_beam_nElasticScatters;   //!
   TBranch        *b_true_beam_elastic_costheta;   //!
   TBranch        *b_true_beam_elastic_X;   //!
   TBranch        *b_true_beam_elastic_Y;   //!
   TBranch        *b_true_beam_elastic_Z;   //!
   TBranch        *b_reco_daughter_allTrack_dirCos;   //!
   TBranch        *b_reco_daughter_allTrack_startX;   //!
   TBranch        *b_reco_daughter_allTrack_startY;   //!
   TBranch        *b_reco_daughter_allTrack_startZ;   //!
   TBranch        *b_reco_daughter_allTrack_endX;   //!
   TBranch        *b_reco_daughter_allTrack_endY;   //!
   TBranch        *b_reco_daughter_allTrack_endZ;   //!
   TBranch        *b_reco_daughter_PFP_true_byHits_startX;   //!
   TBranch        *b_reco_daughter_PFP_true_byHits_startY;   //!
   TBranch        *b_reco_daughter_PFP_true_byHits_startZ;   //!
   TBranch        *b_reco_daughter_PFP_true_byHits_endX;   //!
   TBranch        *b_reco_daughter_PFP_true_byHits_endY;   //!
   TBranch        *b_reco_daughter_PFP_true_byHits_endZ;   //!
   TBranch        *b_true_beam_traj_incidentEnergies;   //!
   TBranch        *b_true_beam_traj_interactingEnergy;   //!
   TBranch        *b_true_beam_traj_slice_z;   //!
   TBranch        *b_true_beam_traj_slice_index;   //!

   anaUnfoldStandardGen(TTree *tree=0);
   virtual ~anaUnfoldStandardGen();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef anaUnfoldStandardGen_cxx
anaUnfoldStandardGen::anaUnfoldStandardGen(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/dune/data/users/rdiurba/standard6GeVkaonTreeSel.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/dune/data/users/rdiurba/standard6GeVkaonTreeSel.root");
      }
      f->GetObject("ana",tree);

   }
   Init(tree);
}

anaUnfoldStandardGen::~anaUnfoldStandardGen()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t anaUnfoldStandardGen::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t anaUnfoldStandardGen::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void anaUnfoldStandardGen::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   reco_beam_true_byE_endProcess = 0;
   true_beam_endProcess = 0;
   reco_beam_calibrated_dEdX_SCE = 0;
   true_beam_incidentEnergies = 0;
   reco_beam_incidentEnergies = 0;
   reco_beam_TrkPitch_SCE = 0;
   reco_beam_calo_wire = 0;
   reco_beam_calo_wire_z = 0;
   true_beam_slices = 0;
   true_beam_slices_dE = 0;
   true_beam_traj_Z = 0;
   true_beam_traj_KE = 0;
   reco_daughter_allTrack_len = 0;
   reco_daughter_allTrack_ID = 0;
   reco_daughter_PFP_true_byHits_PDG = 0;
   reco_daughter_PFP_true_byHits_ID = 0;
   reco_daughter_PFP_true_byHits_process = 0;
   reco_beam_calibrated_incidentEnergies = 0;
   true_beam_elastic_costheta = 0;
   true_beam_elastic_X = 0;
   true_beam_elastic_Y = 0;
   true_beam_elastic_Z = 0;
   reco_daughter_allTrack_dirCos = 0;
   reco_daughter_allTrack_startX = 0;
   reco_daughter_allTrack_startY = 0;
   reco_daughter_allTrack_startZ = 0;
   reco_daughter_allTrack_endX = 0;
   reco_daughter_allTrack_endY = 0;
   reco_daughter_allTrack_endZ = 0;
   reco_daughter_PFP_true_byHits_startX = 0;
   reco_daughter_PFP_true_byHits_startY = 0;
   reco_daughter_PFP_true_byHits_startZ = 0;
   reco_daughter_PFP_true_byHits_endX = 0;
   reco_daughter_PFP_true_byHits_endY = 0;
   reco_daughter_PFP_true_byHits_endZ = 0;
   true_beam_traj_incidentEnergies = 0;
   true_beam_traj_slice_z = 0;
   true_beam_traj_slice_index = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("subrun", &subrun, &b_subrun);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("reco_beam_true_byE_ID", &reco_beam_true_byE_ID, &b_reco_beam_true_byE_ID);
   fChain->SetBranchAddress("reco_beam_true_byE_PDG", &reco_beam_true_byE_PDG, &b_reco_beam_true_byE_PDG);
   fChain->SetBranchAddress("true_beam_ID", &true_beam_ID, &b_true_beam_ID);
   fChain->SetBranchAddress("true_beam_PDG", &true_beam_PDG, &b_true_beam_PDG);
   fChain->SetBranchAddress("true_beam_endP", &true_beam_endP, &b_true_beam_endP);
   fChain->SetBranchAddress("true_beam_startP", &true_beam_startP, &b_true_beam_startP);
   fChain->SetBranchAddress("selection_ID", &selection_ID, &b_selection_ID);
   fChain->SetBranchAddress("sample_ID", &sample_ID, &b_sample_ID);
   fChain->SetBranchAddress("reco_beam_daughter_sameID", &reco_beam_daughter_sameID, &b_reco_beam_daughter_sameID);
   fChain->SetBranchAddress("true_beam_interactingEnergy", &true_beam_interactingEnergy, &b_true_beam_interactingEnergy);
   fChain->SetBranchAddress("reco_beam_interactingEnergy", &reco_beam_interactingEnergy, &b_reco_beam_interactingEnergy);
   fChain->SetBranchAddress("reco_beam_startZ", &reco_beam_startZ, &b_reco_beam_startZ);
   fChain->SetBranchAddress("true_beam_startZ", &true_beam_startZ, &b_true_beam_startZ);
   fChain->SetBranchAddress("reco_beam_startY", &reco_beam_startY, &b_reco_beam_startY);
   fChain->SetBranchAddress("true_beam_startY", &true_beam_startY, &b_true_beam_startY);
   fChain->SetBranchAddress("reco_beam_startX", &reco_beam_startX, &b_reco_beam_startX);
   fChain->SetBranchAddress("true_beam_startX", &true_beam_startX, &b_true_beam_startX);
   fChain->SetBranchAddress("reco_beam_endY", &reco_beam_endY, &b_reco_beam_endY);
   fChain->SetBranchAddress("true_beam_endY", &true_beam_endY, &b_true_beam_endY);
   fChain->SetBranchAddress("reco_beam_endX", &reco_beam_endX, &b_reco_beam_endX);
   fChain->SetBranchAddress("true_beam_endX", &true_beam_endX, &b_true_beam_endX);
   fChain->SetBranchAddress("reco_beam_trackDirZ", &reco_beam_trackDirZ, &b_reco_beam_trackDirZ);
   fChain->SetBranchAddress("reco_beam_trackDirY", &reco_beam_trackDirY, &b_reco_beam_trackDirY);
   fChain->SetBranchAddress("reco_beam_trackDirX", &reco_beam_trackDirX, &b_reco_beam_trackDirX);
   fChain->SetBranchAddress("reco_beam_trackEndDirZ", &reco_beam_trackEndDirZ, &b_reco_beam_trackEndDirZ);
   fChain->SetBranchAddress("reco_beam_trackEndDirY", &reco_beam_trackEndDirY, &b_reco_beam_trackEndDirY);
   fChain->SetBranchAddress("reco_beam_trackEndDirX", &reco_beam_trackEndDirX, &b_reco_beam_trackEndDirX);
   fChain->SetBranchAddress("reco_beam_endZ", &reco_beam_endZ, &b_reco_beam_endZ);
   fChain->SetBranchAddress("true_beam_endZ", &true_beam_endZ, &b_true_beam_endZ);
   fChain->SetBranchAddress("reco_beam_len", &reco_beam_len, &b_reco_beam_len);
   fChain->SetBranchAddress("reco_beam_alt_len", &reco_beam_alt_len, &b_reco_beam_alt_len);
   fChain->SetBranchAddress("reco_beam_true_byE_endProcess", &reco_beam_true_byE_endProcess, &b_reco_beam_true_byE_endProcess);
   fChain->SetBranchAddress("true_beam_endProcess", &true_beam_endProcess, &b_true_beam_endProcess);
   fChain->SetBranchAddress("reco_beam_calibrated_dEdX_SCE", &reco_beam_calibrated_dEdX_SCE, &b_reco_beam_calibrated_dEdX_SCE);
   fChain->SetBranchAddress("true_beam_incidentEnergies", &true_beam_incidentEnergies, &b_true_beam_incidentEnergies);
   fChain->SetBranchAddress("reco_beam_incidentEnergies", &reco_beam_incidentEnergies, &b_reco_beam_incidentEnergies);
   fChain->SetBranchAddress("reco_beam_TrkPitch_SCE", &reco_beam_TrkPitch_SCE, &b_reco_beam_TrkPitch_SCE);
   fChain->SetBranchAddress("reco_beam_calo_wire", &reco_beam_calo_wire, &b_reco_beam_calo_wire);
   fChain->SetBranchAddress("reco_beam_calo_wire_z", &reco_beam_calo_wire_z, &b_reco_beam_calo_wire_z);
   fChain->SetBranchAddress("true_beam_slices", &true_beam_slices, &b_true_beam_slices);
   fChain->SetBranchAddress("true_beam_slices_dE", &true_beam_slices_dE, &b_true_beam_slices_dE);
   fChain->SetBranchAddress("beam_inst_KE", &beam_inst_KE, &b_beam_inst_KE);
   fChain->SetBranchAddress("beam_inst_P", &beam_inst_P, &b_beam_inst_P);
   fChain->SetBranchAddress("beam_inst_X", &beam_inst_X, &b_beam_inst_X);
   fChain->SetBranchAddress("beam_inst_Y", &beam_inst_Y, &b_beam_inst_Y);
   fChain->SetBranchAddress("beam_inst_Z", &beam_inst_Z, &b_beam_inst_Z);
   fChain->SetBranchAddress("beam_inst_dirX", &beam_inst_dirX, &b_beam_inst_dirX);
   fChain->SetBranchAddress("beam_inst_dirY", &beam_inst_dirY, &b_beam_inst_dirY);
   fChain->SetBranchAddress("beam_inst_dirZ", &beam_inst_dirZ, &b_beam_inst_dirZ);
   fChain->SetBranchAddress("reco_beam_dCos", &reco_beam_dCos, &b_reco_beam_dCos);
   fChain->SetBranchAddress("reco_beam_dirCos", &reco_beam_dirCos, &b_reco_beam_dirCos);
   fChain->SetBranchAddress("reco_beam_zEndPointCos", &reco_beam_zEndPointCos, &b_reco_beam_zEndPointCos);
   fChain->SetBranchAddress("reco_beam_cosXZ", &reco_beam_cosXZ, &b_reco_beam_cosXZ);
   fChain->SetBranchAddress("reco_beam_cosYZ", &reco_beam_cosYZ, &b_reco_beam_cosYZ);
   fChain->SetBranchAddress("true_beam_traj_Z", &true_beam_traj_Z, &b_true_beam_traj_Z);
   fChain->SetBranchAddress("true_beam_traj_KE", &true_beam_traj_KE, &b_true_beam_traj_KE);
   fChain->SetBranchAddress("reco_beam_type", &reco_beam_type, &b_reco_beam_type);
   fChain->SetBranchAddress("reco_daughter_allTrack_len", &reco_daughter_allTrack_len, &b_reco_daughter_allTrack_len);
   fChain->SetBranchAddress("reco_daughter_allTrack_ID", &reco_daughter_allTrack_ID, &b_reco_daughter_allTrack_ID);
   fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_PDG", &reco_daughter_PFP_true_byHits_PDG, &b_reco_daughter_PFP_true_byHits_PDG);
   fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_ID", &reco_daughter_PFP_true_byHits_ID, &b_reco_daughter_PFP_true_byHits_ID);
   fChain->SetBranchAddress("reco_beam_nDaughters", &reco_beam_nDaughters, &b_reco_beam_nDaughters);
   fChain->SetBranchAddress("reco_daughter_PFP_true_nDaughters", &reco_daughter_PFP_true_nDaughters, &b_reco_daughter_PFP_true_nDaughters);
   fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_process", &reco_daughter_PFP_true_byHits_process, &b_reco_daughter_PFP_true_byHits_process);
   fChain->SetBranchAddress("reco_beam_calibrated_interactingEnergy", &reco_beam_calibrated_interactingEnergy, &b_reco_beam_calibrated_interactingEnergy);
   fChain->SetBranchAddress("reco_beam_calibrated_incidentEnergies", &reco_beam_calibrated_incidentEnergies, &b_reco_beam_calibrated_incidentEnergies);
   fChain->SetBranchAddress("true_beam_nElasticScatters", &true_beam_nElasticScatters, &b_true_beam_nElasticScatters);
   fChain->SetBranchAddress("true_beam_elastic_costheta", &true_beam_elastic_costheta, &b_true_beam_elastic_costheta);
   fChain->SetBranchAddress("true_beam_elastic_X", &true_beam_elastic_X, &b_true_beam_elastic_X);
   fChain->SetBranchAddress("true_beam_elastic_Y", &true_beam_elastic_Y, &b_true_beam_elastic_Y);
   fChain->SetBranchAddress("true_beam_elastic_Z", &true_beam_elastic_Z, &b_true_beam_elastic_Z);
   fChain->SetBranchAddress("reco_daughter_allTrack_dirCos", &reco_daughter_allTrack_dirCos, &b_reco_daughter_allTrack_dirCos);
   fChain->SetBranchAddress("reco_daughter_allTrack_startX", &reco_daughter_allTrack_startX, &b_reco_daughter_allTrack_startX);
   fChain->SetBranchAddress("reco_daughter_allTrack_startY", &reco_daughter_allTrack_startY, &b_reco_daughter_allTrack_startY);
   fChain->SetBranchAddress("reco_daughter_allTrack_startZ", &reco_daughter_allTrack_startZ, &b_reco_daughter_allTrack_startZ);
   fChain->SetBranchAddress("reco_daughter_allTrack_endX", &reco_daughter_allTrack_endX, &b_reco_daughter_allTrack_endX);
   fChain->SetBranchAddress("reco_daughter_allTrack_endY", &reco_daughter_allTrack_endY, &b_reco_daughter_allTrack_endY);
   fChain->SetBranchAddress("reco_daughter_allTrack_endZ", &reco_daughter_allTrack_endZ, &b_reco_daughter_allTrack_endZ);
   fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_startX", &reco_daughter_PFP_true_byHits_startX, &b_reco_daughter_PFP_true_byHits_startX);
   fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_startY", &reco_daughter_PFP_true_byHits_startY, &b_reco_daughter_PFP_true_byHits_startY);
   fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_startZ", &reco_daughter_PFP_true_byHits_startZ, &b_reco_daughter_PFP_true_byHits_startZ);
   fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_endX", &reco_daughter_PFP_true_byHits_endX, &b_reco_daughter_PFP_true_byHits_endX);
   fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_endY", &reco_daughter_PFP_true_byHits_endY, &b_reco_daughter_PFP_true_byHits_endY);
   fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_endZ", &reco_daughter_PFP_true_byHits_endZ, &b_reco_daughter_PFP_true_byHits_endZ);
   fChain->SetBranchAddress("true_beam_traj_incidentEnergies", &true_beam_traj_incidentEnergies, &b_true_beam_traj_incidentEnergies);
   fChain->SetBranchAddress("true_beam_traj_interactingEnergy", &true_beam_traj_interactingEnergy, &b_true_beam_traj_interactingEnergy);
   fChain->SetBranchAddress("true_beam_traj_slice_z", &true_beam_traj_slice_z, &b_true_beam_traj_slice_z);
   fChain->SetBranchAddress("true_beam_traj_slice_index", &true_beam_traj_slice_index, &b_true_beam_traj_slice_index);
   Notify();
}

Bool_t anaUnfoldStandardGen::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void anaUnfoldStandardGen::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t anaUnfoldStandardGen::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef anaUnfoldStandardGen_cxx
