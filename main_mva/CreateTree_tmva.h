//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul 18 14:32:24 2016 by ROOT version 5.34/18
// from TTree tree/PhysicsTools.Heppy.analyzers.core.AutoFillTreeProducer.AutoFillTreeProducer_1
// found on file: VBFSpin0ToBBbar_W_1p0_M_900_TuneCUEP8M1_13TeV_pythia8_v21.root
//////////////////////////////////////////////////////////

#ifndef CreateTree_tmva_h
#define CreateTree_tmva_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class CreateTree_tmva {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
	TString file_name;

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          lumi;
   ULong64_t       evt;
   Int_t           isData;
   Float_t         xsec;
   Float_t         puWeight;
   Float_t         nTrueInt;
   Float_t         genWeight;
   Float_t         puWeightUp;
   Float_t         puWeightDown;
   Float_t         json;
   Float_t         json_silver;
   Float_t         nPU0;
   Float_t         nPVs;
   Float_t         Vtype;
   Float_t         rho;
   Float_t         rhoN;
   Float_t         rhoCHPU;
   Float_t         rhoCentral;
   Float_t         lheNj;
   Float_t         lheNb;
   Float_t         lheNc;
   Float_t         lheNg;
   Float_t         lheNl;
   Float_t         lheV_pt;
   Float_t         lheHT;
   Float_t         met_sig;
   Float_t         met_rawpt;
   Float_t         metPuppi_pt;
   Float_t         metPuppi_phi;
   Float_t         metPuppi_rawpt;
   Float_t         metType1p2_pt;
   Float_t         bTagWeight_LFUp;
   Float_t         bTagWeight_LFStats2Down;
   Float_t         bTagWeight_LFDown;
   Float_t         bTagWeight_HFUp;
   Float_t         bTagWeight_HFStats1Up;
   Float_t         bTagWeight_cErr1Down;
   Float_t         bTagWeight_cErr2Up;
   Float_t         bTagWeight_cErr1Up;
   Float_t         bTagWeight_LFStats1Down;
   Float_t         bTagWeight_JESDown;
   Float_t         bTagWeight_LFStats1Up;
   Float_t         bTagWeight;
   Float_t         bTagWeight_HFDown;
   Float_t         bTagWeight_LFStats2Up;
   Float_t         bTagWeight_JESUp;
   Float_t         bTagWeight_HFStats2Up;
   Float_t         bTagWeight_cErr2Down;
   Float_t         bTagWeight_HFStats1Down;
   Float_t         bTagWeight_HFStats2Down;
   Int_t           HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq200_v;
   Int_t           HLT_BIT_HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq460_v;
   Int_t           HLT_BIT_HLT_DiPFJetAve40_v;
   Int_t           HLT_BIT_HLT_DiPFJetAve60_v;
   Int_t           HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq240_v;
   Int_t           HLT_BIT_HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq500_v;
   Int_t           HLT_BIT_HLT_QuadPFJet_VBF_v;
   Int_t           HLT_BIT_HLT_L1_TripleJet_VBF_v;
   Int_t           softActivityVH_njets2;
   Int_t           softActivityVH_njets5;
   Int_t           softActivityVH_njets10;
   Float_t         softActivityVH_HT;
   Float_t         met_shifted_JetResUp_pt;
   Float_t         met_shifted_JetResUp_phi;
   Float_t         met_shifted_JetResUp_sumEt;
   Float_t         met_pt;
   Float_t         met_eta;
   Float_t         met_phi;
   Float_t         met_mass;
   Float_t         met_sumEt;
   Float_t         met_rawPt;
   Float_t         met_rawPhi;
   Float_t         met_rawSumEt;
   Float_t         met_genPt;
   Float_t         met_genPhi;
   Float_t         met_genEta;
   Float_t         met_shifted_MuonEnDown_pt;
   Float_t         met_shifted_MuonEnDown_phi;
   Float_t         met_shifted_MuonEnDown_sumEt;
   Float_t         met_shifted_ElectronEnUp_pt;
   Float_t         met_shifted_ElectronEnUp_phi;
   Float_t         met_shifted_ElectronEnUp_sumEt;
   Float_t         met_shifted_ElectronEnDown_pt;
   Float_t         met_shifted_ElectronEnDown_phi;
   Float_t         met_shifted_ElectronEnDown_sumEt;
   Float_t         met_shifted_TauEnDown_pt;
   Float_t         met_shifted_TauEnDown_phi;
   Float_t         met_shifted_TauEnDown_sumEt;
   Float_t         met_shifted_TauEnUp_pt;
   Float_t         met_shifted_TauEnUp_phi;
   Float_t         met_shifted_TauEnUp_sumEt;
   Float_t         met_shifted_UnclusteredEnUp_pt;
   Float_t         met_shifted_UnclusteredEnUp_phi;
   Float_t         met_shifted_UnclusteredEnUp_sumEt;
   Float_t         met_shifted_UnclusteredEnDown_pt;
   Float_t         met_shifted_UnclusteredEnDown_phi;
   Float_t         met_shifted_UnclusteredEnDown_sumEt;
   Float_t         met_shifted_JetEnUp_pt;
   Float_t         met_shifted_JetEnUp_phi;
   Float_t         met_shifted_JetEnUp_sumEt;
   Float_t         met_shifted_JetEnDown_pt;
   Float_t         met_shifted_JetEnDown_phi;
   Float_t         met_shifted_JetEnDown_sumEt;
   Float_t         met_shifted_JetResDown_pt;
   Float_t         met_shifted_JetResDown_phi;
   Float_t         met_shifted_JetResDown_sumEt;
   Int_t           softActivity_njets2;
   Int_t           softActivity_njets5;
   Int_t           softActivity_njets10;
   Float_t         softActivity_HT;
   Float_t         met_shifted_MuonEnUp_pt;
   Float_t         met_shifted_MuonEnUp_phi;
   Float_t         met_shifted_MuonEnUp_sumEt;
   Int_t           nGenBQuarkFromHafterISR;
   Int_t           GenBQuarkFromHafterISR_pdgId[1];   //[nGenBQuarkFromHafterISR]
   Float_t         GenBQuarkFromHafterISR_pt[1];   //[nGenBQuarkFromHafterISR]
   Float_t         GenBQuarkFromHafterISR_eta[1];   //[nGenBQuarkFromHafterISR]
   Float_t         GenBQuarkFromHafterISR_phi[1];   //[nGenBQuarkFromHafterISR]
   Float_t         GenBQuarkFromHafterISR_mass[1];   //[nGenBQuarkFromHafterISR]
   Float_t         GenBQuarkFromHafterISR_charge[1];   //[nGenBQuarkFromHafterISR]
   Int_t           GenBQuarkFromHafterISR_status[1];   //[nGenBQuarkFromHafterISR]
   Int_t           nvLeptons;
   Float_t         vLeptons_dz[2];   //[nvLeptons]
   Int_t           nLHE_weights_scale;
   Int_t           LHE_weights_scale_id[1];   //[nLHE_weights_scale]
   Float_t         LHE_weights_scale_wgt[1];   //[nLHE_weights_scale]
   Int_t           nGenBQuarkFromH;
   Int_t           GenBQuarkFromH_pdgId[2];   //[nGenBQuarkFromH]
   Float_t         GenBQuarkFromH_pt[2];   //[nGenBQuarkFromH]
   Float_t         GenBQuarkFromH_eta[2];   //[nGenBQuarkFromH]
   Float_t         GenBQuarkFromH_phi[2];   //[nGenBQuarkFromH]
   Float_t         GenBQuarkFromH_mass[2];   //[nGenBQuarkFromH]
   Float_t         GenBQuarkFromH_charge[2];   //[nGenBQuarkFromH]
   Int_t           GenBQuarkFromH_status[2];   //[nGenBQuarkFromH]
   Int_t           nGenHiggsSisters;
   Int_t           GenHiggsSisters_pdgId[2];   //[nGenHiggsSisters]
   Float_t         GenHiggsSisters_pt[2];   //[nGenHiggsSisters]
   Float_t         GenHiggsSisters_eta[2];   //[nGenHiggsSisters]
   Float_t         GenHiggsSisters_phi[2];   //[nGenHiggsSisters]
   Float_t         GenHiggsSisters_mass[2];   //[nGenHiggsSisters]
   Float_t         GenHiggsSisters_charge[2];   //[nGenHiggsSisters]
   Int_t           GenHiggsSisters_status[2];   //[nGenHiggsSisters]
   Int_t           nsoftActivityVHJets;
   Float_t         softActivityVHJets_pt[5];   //[nsoftActivityVHJets]
   Float_t         softActivityVHJets_eta[5];   //[nsoftActivityVHJets]
   Float_t         softActivityVHJets_phi[5];   //[nsoftActivityVHJets]
   Float_t         softActivityVHJets_mass[5];   //[nsoftActivityVHJets]
   Int_t           nselLeptons;
   Int_t           selLeptons_charge[3];   //[nselLeptons]
   Int_t           selLeptons_tightId[3];   //[nselLeptons]
   Float_t         selLeptons_relIso03[3];   //[nselLeptons]
   Int_t           selLeptons_pdgId[3];   //[nselLeptons]
   Float_t         selLeptons_pfRelIso03[3];   //[nselLeptons]
   Int_t           nJet;
   Int_t           Jet_id[18];   //[nJet]
   Int_t           Jet_puId[18];   //[nJet]
   Float_t         Jet_btagCSV[18];   //[nJet]
   Int_t           Jet_hadronFlavour[18];   //[nJet]
   Float_t         Jet_corr_JECUp[18];   //[nJet]
   Float_t         Jet_corr_JECDown[18];   //[nJet]
   Float_t         Jet_corr[18];   //[nJet]
   Float_t         Jet_corr_JERUp[18];   //[nJet]
   Float_t         Jet_corr_JERDown[18];   //[nJet]
   Float_t         Jet_corr_JER[18];   //[nJet]
   Float_t         Jet_pt[18];   //[nJet]
   Float_t         Jet_eta[18];   //[nJet]
   Float_t         Jet_phi[18];   //[nJet]
   Float_t         Jet_mass[18];   //[nJet]
   Float_t         Jet_leadTrackPt[18];   //[nJet]
   Float_t         Jet_qgl[18];   //[nJet]
   Float_t         Jet_ptd[18];   //[nJet]
   Float_t         Jet_axis2[18];   //[nJet]
   Int_t           Jet_mult[18];   //[nJet]
   Float_t         Jet_blike_VBF[18];   //[nJet]
   Float_t         Jet_pt_reg[18];   //[nJet]
   Float_t         Jet_pt_regVBF[18];   //[nJet]
   Float_t         Jet_pt_reg_corrJECUp[18];   //[nJet]
   Float_t         Jet_pt_regVBF_corrJECUp[18];   //[nJet]
   Float_t         Jet_pt_reg_corrJECDown[18];   //[nJet]
   Float_t         Jet_pt_regVBF_corrJECDown[18];   //[nJet]
   Float_t         Jet_pt_reg_corrJERUp[18];   //[nJet]
   Float_t         Jet_pt_regVBF_corrJERUp[18];   //[nJet]
   Float_t         Jet_pt_reg_corrJERDown[18];   //[nJet]
   Float_t         Jet_pt_regVBF_corrJERDown[18];   //[nJet]
   Int_t           nLHE_weights_pdf;
   Int_t           LHE_weights_pdf_id[1];   //[nLHE_weights_pdf]
   Float_t         LHE_weights_pdf_wgt[1];   //[nLHE_weights_pdf]
   Int_t           nsoftActivityJets;
   Float_t         softActivityJets_pt[5];   //[nsoftActivityJets]
   Float_t         softActivityJets_eta[5];   //[nsoftActivityJets]
   Float_t         softActivityJets_phi[5];   //[nsoftActivityJets]
   Float_t         softActivityJets_mass[5];   //[nsoftActivityJets]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_xsec;   //!
   TBranch        *b_puWeight;   //!
   TBranch        *b_nTrueInt;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_puWeightUp;   //!
   TBranch        *b_puWeightDown;   //!
   TBranch        *b_json;   //!
   TBranch        *b_json_silver;   //!
   TBranch        *b_nPU0;   //!
   TBranch        *b_nPVs;   //!
   TBranch        *b_Vtype;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_rhoN;   //!
   TBranch        *b_rhoCHPU;   //!
   TBranch        *b_rhoCentral;   //!
   TBranch        *b_lheNj;   //!
   TBranch        *b_lheNb;   //!
   TBranch        *b_lheNc;   //!
   TBranch        *b_lheNg;   //!
   TBranch        *b_lheNl;   //!
   TBranch        *b_lheV_pt;   //!
   TBranch        *b_lheHT;   //!
   TBranch        *b_met_sig;   //!
   TBranch        *b_met_rawpt;   //!
   TBranch        *b_metPuppi_pt;   //!
   TBranch        *b_metPuppi_phi;   //!
   TBranch        *b_metPuppi_rawpt;   //!
   TBranch        *b_metType1p2_pt;   //!
   TBranch        *b_bTagWeight_LFUp;   //!
   TBranch        *b_bTagWeight_LFStats2Down;   //!
   TBranch        *b_bTagWeight_LFDown;   //!
   TBranch        *b_bTagWeight_HFUp;   //!
   TBranch        *b_bTagWeight_HFStats1Up;   //!
   TBranch        *b_bTagWeight_cErr1Down;   //!
   TBranch        *b_bTagWeight_cErr2Up;   //!
   TBranch        *b_bTagWeight_cErr1Up;   //!
   TBranch        *b_bTagWeight_LFStats1Down;   //!
   TBranch        *b_bTagWeight_JESDown;   //!
   TBranch        *b_bTagWeight_LFStats1Up;   //!
   TBranch        *b_bTagWeight;   //!
   TBranch        *b_bTagWeight_HFDown;   //!
   TBranch        *b_bTagWeight_LFStats2Up;   //!
   TBranch        *b_bTagWeight_JESUp;   //!
   TBranch        *b_bTagWeight_HFStats2Up;   //!
   TBranch        *b_bTagWeight_cErr2Down;   //!
   TBranch        *b_bTagWeight_HFStats1Down;   //!
   TBranch        *b_bTagWeight_HFStats2Down;   //!
   TBranch        *b_HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq200_v;   //!
   TBranch        *b_HLT_BIT_HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq460_v;   //!
   TBranch        *b_HLT_BIT_HLT_DiPFJetAve40_v;   //!
   TBranch        *b_HLT_BIT_HLT_DiPFJetAve60_v;   //!
   TBranch        *b_HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq240_v;   //!
   TBranch        *b_HLT_BIT_HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq500_v;   //!
   TBranch        *b_HLT_BIT_HLT_QuadPFJet_VBF_v;   //!
   TBranch        *b_HLT_BIT_HLT_L1_TripleJet_VBF_v;   //!
   TBranch        *b_softActivityVH_njets2;   //!
   TBranch        *b_softActivityVH_njets5;   //!
   TBranch        *b_softActivityVH_njets10;   //!
   TBranch        *b_softActivityVH_HT;   //!
   TBranch        *b_met_shifted_JetResUp_pt;   //!
   TBranch        *b_met_shifted_JetResUp_phi;   //!
   TBranch        *b_met_shifted_JetResUp_sumEt;   //!
   TBranch        *b_met_pt;   //!
   TBranch        *b_met_eta;   //!
   TBranch        *b_met_phi;   //!
   TBranch        *b_met_mass;   //!
   TBranch        *b_met_sumEt;   //!
   TBranch        *b_met_rawPt;   //!
   TBranch        *b_met_rawPhi;   //!
   TBranch        *b_met_rawSumEt;   //!
   TBranch        *b_met_genPt;   //!
   TBranch        *b_met_genPhi;   //!
   TBranch        *b_met_genEta;   //!
   TBranch        *b_met_shifted_MuonEnDown_pt;   //!
   TBranch        *b_met_shifted_MuonEnDown_phi;   //!
   TBranch        *b_met_shifted_MuonEnDown_sumEt;   //!
   TBranch        *b_met_shifted_ElectronEnUp_pt;   //!
   TBranch        *b_met_shifted_ElectronEnUp_phi;   //!
   TBranch        *b_met_shifted_ElectronEnUp_sumEt;   //!
   TBranch        *b_met_shifted_ElectronEnDown_pt;   //!
   TBranch        *b_met_shifted_ElectronEnDown_phi;   //!
   TBranch        *b_met_shifted_ElectronEnDown_sumEt;   //!
   TBranch        *b_met_shifted_TauEnDown_pt;   //!
   TBranch        *b_met_shifted_TauEnDown_phi;   //!
   TBranch        *b_met_shifted_TauEnDown_sumEt;   //!
   TBranch        *b_met_shifted_TauEnUp_pt;   //!
   TBranch        *b_met_shifted_TauEnUp_phi;   //!
   TBranch        *b_met_shifted_TauEnUp_sumEt;   //!
   TBranch        *b_met_shifted_UnclusteredEnUp_pt;   //!
   TBranch        *b_met_shifted_UnclusteredEnUp_phi;   //!
   TBranch        *b_met_shifted_UnclusteredEnUp_sumEt;   //!
   TBranch        *b_met_shifted_UnclusteredEnDown_pt;   //!
   TBranch        *b_met_shifted_UnclusteredEnDown_phi;   //!
   TBranch        *b_met_shifted_UnclusteredEnDown_sumEt;   //!
   TBranch        *b_met_shifted_JetEnUp_pt;   //!
   TBranch        *b_met_shifted_JetEnUp_phi;   //!
   TBranch        *b_met_shifted_JetEnUp_sumEt;   //!
   TBranch        *b_met_shifted_JetEnDown_pt;   //!
   TBranch        *b_met_shifted_JetEnDown_phi;   //!
   TBranch        *b_met_shifted_JetEnDown_sumEt;   //!
   TBranch        *b_met_shifted_JetResDown_pt;   //!
   TBranch        *b_met_shifted_JetResDown_phi;   //!
   TBranch        *b_met_shifted_JetResDown_sumEt;   //!
   TBranch        *b_softActivity_njets2;   //!
   TBranch        *b_softActivity_njets5;   //!
   TBranch        *b_softActivity_njets10;   //!
   TBranch        *b_softActivity_HT;   //!
   TBranch        *b_met_shifted_MuonEnUp_pt;   //!
   TBranch        *b_met_shifted_MuonEnUp_phi;   //!
   TBranch        *b_met_shifted_MuonEnUp_sumEt;   //!
   TBranch        *b_nGenBQuarkFromHafterISR;   //!
   TBranch        *b_GenBQuarkFromHafterISR_pdgId;   //!
   TBranch        *b_GenBQuarkFromHafterISR_pt;   //!
   TBranch        *b_GenBQuarkFromHafterISR_eta;   //!
   TBranch        *b_GenBQuarkFromHafterISR_phi;   //!
   TBranch        *b_GenBQuarkFromHafterISR_mass;   //!
   TBranch        *b_GenBQuarkFromHafterISR_charge;   //!
   TBranch        *b_GenBQuarkFromHafterISR_status;   //!
   TBranch        *b_nvLeptons;   //!
   TBranch        *b_vLeptons_dz;   //!
   TBranch        *b_nLHE_weights_scale;   //!
   TBranch        *b_LHE_weights_scale_id;   //!
   TBranch        *b_LHE_weights_scale_wgt;   //!
   TBranch        *b_nGenBQuarkFromH;   //!
   TBranch        *b_GenBQuarkFromH_pdgId;   //!
   TBranch        *b_GenBQuarkFromH_pt;   //!
   TBranch        *b_GenBQuarkFromH_eta;   //!
   TBranch        *b_GenBQuarkFromH_phi;   //!
   TBranch        *b_GenBQuarkFromH_mass;   //!
   TBranch        *b_GenBQuarkFromH_charge;   //!
   TBranch        *b_GenBQuarkFromH_status;   //!
   TBranch        *b_nGenHiggsSisters;   //!
   TBranch        *b_GenHiggsSisters_pdgId;   //!
   TBranch        *b_GenHiggsSisters_pt;   //!
   TBranch        *b_GenHiggsSisters_eta;   //!
   TBranch        *b_GenHiggsSisters_phi;   //!
   TBranch        *b_GenHiggsSisters_mass;   //!
   TBranch        *b_GenHiggsSisters_charge;   //!
   TBranch        *b_GenHiggsSisters_status;   //!
   TBranch        *b_nsoftActivityVHJets;   //!
   TBranch        *b_softActivityVHJets_pt;   //!
   TBranch        *b_softActivityVHJets_eta;   //!
   TBranch        *b_softActivityVHJets_phi;   //!
   TBranch        *b_softActivityVHJets_mass;   //!
   TBranch        *b_nselLeptons;   //!
   TBranch        *b_selLeptons_charge;   //!
   TBranch        *b_selLeptons_tightId;   //!
   TBranch        *b_selLeptons_relIso03;   //!
   TBranch        *b_selLeptons_pdgId;   //!
   TBranch        *b_selLeptons_pfRelIso03;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_Jet_id;   //!
   TBranch        *b_Jet_puId;   //!
   TBranch        *b_Jet_btagCSV;   //!
   TBranch        *b_Jet_hadronFlavour;   //!
   TBranch        *b_Jet_corr_JECUp;   //!
   TBranch        *b_Jet_corr_JECDown;   //!
   TBranch        *b_Jet_corr;   //!
   TBranch        *b_Jet_corr_JERUp;   //!
   TBranch        *b_Jet_corr_JERDown;   //!
   TBranch        *b_Jet_corr_JER;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_mass;   //!
   TBranch        *b_Jet_leadTrackPt;   //!
   TBranch        *b_Jet_qgl;   //!
   TBranch        *b_Jet_ptd;   //!
   TBranch        *b_Jet_axis2;   //!
   TBranch        *b_Jet_mult;   //!
   TBranch        *b_Jet_blike_VBF;   //!
   TBranch        *b_Jet_pt_reg;   //!
   TBranch        *b_Jet_pt_regVBF;   //!
   TBranch        *b_Jet_pt_reg_corrJECUp;   //!
   TBranch        *b_Jet_pt_regVBF_corrJECUp;   //!
   TBranch        *b_Jet_pt_reg_corrJECDown;   //!
   TBranch        *b_Jet_pt_regVBF_corrJECDown;   //!
   TBranch        *b_Jet_pt_reg_corrJERUp;   //!
   TBranch        *b_Jet_pt_regVBF_corrJERUp;   //!
   TBranch        *b_Jet_pt_reg_corrJERDown;   //!
   TBranch        *b_Jet_pt_regVBF_corrJERDown;   //!
   TBranch        *b_nLHE_weights_pdf;   //!
   TBranch        *b_LHE_weights_pdf_id;   //!
   TBranch        *b_LHE_weights_pdf_wgt;   //!
   TBranch        *b_nsoftActivityJets;   //!
   TBranch        *b_softActivityJets_pt;   //!
   TBranch        *b_softActivityJets_eta;   //!
   TBranch        *b_softActivityJets_phi;   //!
   TBranch        *b_softActivityJets_mass;   //!

   CreateTree_tmva(TTree *tree=0, TString filename="");
   virtual ~CreateTree_tmva();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TString inputfile="",TString output_dir="", int data=0);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef CreateTree_tmva_cxx
CreateTree_tmva::CreateTree_tmva(TTree *tree, TString filename) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
			file_name = filename;
      	TFile *f = TFile::Open(file_name);
      f->GetObject("tree",tree);
   }
   Init(tree);
}

CreateTree_tmva::~CreateTree_tmva()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t CreateTree_tmva::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t CreateTree_tmva::LoadTree(Long64_t entry)
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

void CreateTree_tmva::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("xsec", &xsec, &b_xsec);
   fChain->SetBranchAddress("puWeight", &puWeight, &b_puWeight);
   fChain->SetBranchAddress("nTrueInt", &nTrueInt, &b_nTrueInt);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("puWeightUp", &puWeightUp, &b_puWeightUp);
   fChain->SetBranchAddress("puWeightDown", &puWeightDown, &b_puWeightDown);
   fChain->SetBranchAddress("json", &json, &b_json);
   fChain->SetBranchAddress("json_silver", &json_silver, &b_json_silver);
   fChain->SetBranchAddress("nPU0", &nPU0, &b_nPU0);
   fChain->SetBranchAddress("nPVs", &nPVs, &b_nPVs);
   fChain->SetBranchAddress("Vtype", &Vtype, &b_Vtype);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("rhoN", &rhoN, &b_rhoN);
   fChain->SetBranchAddress("rhoCHPU", &rhoCHPU, &b_rhoCHPU);
   fChain->SetBranchAddress("rhoCentral", &rhoCentral, &b_rhoCentral);
   fChain->SetBranchAddress("lheNj", &lheNj, &b_lheNj);
   fChain->SetBranchAddress("lheNb", &lheNb, &b_lheNb);
   fChain->SetBranchAddress("lheNc", &lheNc, &b_lheNc);
   fChain->SetBranchAddress("lheNg", &lheNg, &b_lheNg);
   fChain->SetBranchAddress("lheNl", &lheNl, &b_lheNl);
   fChain->SetBranchAddress("lheV_pt", &lheV_pt, &b_lheV_pt);
   fChain->SetBranchAddress("lheHT", &lheHT, &b_lheHT);
   fChain->SetBranchAddress("met_sig", &met_sig, &b_met_sig);
   fChain->SetBranchAddress("met_rawpt", &met_rawpt, &b_met_rawpt);
   fChain->SetBranchAddress("metPuppi_pt", &metPuppi_pt, &b_metPuppi_pt);
   fChain->SetBranchAddress("metPuppi_phi", &metPuppi_phi, &b_metPuppi_phi);
   fChain->SetBranchAddress("metPuppi_rawpt", &metPuppi_rawpt, &b_metPuppi_rawpt);
   fChain->SetBranchAddress("metType1p2_pt", &metType1p2_pt, &b_metType1p2_pt);
   fChain->SetBranchAddress("bTagWeight_LFUp", &bTagWeight_LFUp, &b_bTagWeight_LFUp);
   fChain->SetBranchAddress("bTagWeight_LFStats2Down", &bTagWeight_LFStats2Down, &b_bTagWeight_LFStats2Down);
   fChain->SetBranchAddress("bTagWeight_LFDown", &bTagWeight_LFDown, &b_bTagWeight_LFDown);
   fChain->SetBranchAddress("bTagWeight_HFUp", &bTagWeight_HFUp, &b_bTagWeight_HFUp);
   fChain->SetBranchAddress("bTagWeight_HFStats1Up", &bTagWeight_HFStats1Up, &b_bTagWeight_HFStats1Up);
   fChain->SetBranchAddress("bTagWeight_cErr1Down", &bTagWeight_cErr1Down, &b_bTagWeight_cErr1Down);
   fChain->SetBranchAddress("bTagWeight_cErr2Up", &bTagWeight_cErr2Up, &b_bTagWeight_cErr2Up);
   fChain->SetBranchAddress("bTagWeight_cErr1Up", &bTagWeight_cErr1Up, &b_bTagWeight_cErr1Up);
   fChain->SetBranchAddress("bTagWeight_LFStats1Down", &bTagWeight_LFStats1Down, &b_bTagWeight_LFStats1Down);
   fChain->SetBranchAddress("bTagWeight_JESDown", &bTagWeight_JESDown, &b_bTagWeight_JESDown);
   fChain->SetBranchAddress("bTagWeight_LFStats1Up", &bTagWeight_LFStats1Up, &b_bTagWeight_LFStats1Up);
   fChain->SetBranchAddress("bTagWeight", &bTagWeight, &b_bTagWeight);
   fChain->SetBranchAddress("bTagWeight_HFDown", &bTagWeight_HFDown, &b_bTagWeight_HFDown);
   fChain->SetBranchAddress("bTagWeight_LFStats2Up", &bTagWeight_LFStats2Up, &b_bTagWeight_LFStats2Up);
   fChain->SetBranchAddress("bTagWeight_JESUp", &bTagWeight_JESUp, &b_bTagWeight_JESUp);
   fChain->SetBranchAddress("bTagWeight_HFStats2Up", &bTagWeight_HFStats2Up, &b_bTagWeight_HFStats2Up);
   fChain->SetBranchAddress("bTagWeight_cErr2Down", &bTagWeight_cErr2Down, &b_bTagWeight_cErr2Down);
   fChain->SetBranchAddress("bTagWeight_HFStats1Down", &bTagWeight_HFStats1Down, &b_bTagWeight_HFStats1Down);
   fChain->SetBranchAddress("bTagWeight_HFStats2Down", &bTagWeight_HFStats2Down, &b_bTagWeight_HFStats2Down);
   fChain->SetBranchAddress("HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq200_v", &HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq200_v, &b_HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq200_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq460_v", &HLT_BIT_HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq460_v, &b_HLT_BIT_HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq460_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_DiPFJetAve40_v", &HLT_BIT_HLT_DiPFJetAve40_v, &b_HLT_BIT_HLT_DiPFJetAve40_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_DiPFJetAve60_v", &HLT_BIT_HLT_DiPFJetAve60_v, &b_HLT_BIT_HLT_DiPFJetAve60_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq240_v", &HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq240_v, &b_HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq240_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq500_v", &HLT_BIT_HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq500_v, &b_HLT_BIT_HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq500_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_QuadPFJet_VBF_v", &HLT_BIT_HLT_QuadPFJet_VBF_v, &b_HLT_BIT_HLT_QuadPFJet_VBF_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_L1_TripleJet_VBF_v", &HLT_BIT_HLT_L1_TripleJet_VBF_v, &b_HLT_BIT_HLT_L1_TripleJet_VBF_v);
   fChain->SetBranchAddress("softActivityVH_njets2", &softActivityVH_njets2, &b_softActivityVH_njets2);
   fChain->SetBranchAddress("softActivityVH_njets5", &softActivityVH_njets5, &b_softActivityVH_njets5);
   fChain->SetBranchAddress("softActivityVH_njets10", &softActivityVH_njets10, &b_softActivityVH_njets10);
   fChain->SetBranchAddress("softActivityVH_HT", &softActivityVH_HT, &b_softActivityVH_HT);
   fChain->SetBranchAddress("met_shifted_JetResUp_pt", &met_shifted_JetResUp_pt, &b_met_shifted_JetResUp_pt);
   fChain->SetBranchAddress("met_shifted_JetResUp_phi", &met_shifted_JetResUp_phi, &b_met_shifted_JetResUp_phi);
   fChain->SetBranchAddress("met_shifted_JetResUp_sumEt", &met_shifted_JetResUp_sumEt, &b_met_shifted_JetResUp_sumEt);
   fChain->SetBranchAddress("met_pt", &met_pt, &b_met_pt);
   fChain->SetBranchAddress("met_eta", &met_eta, &b_met_eta);
   fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
   fChain->SetBranchAddress("met_mass", &met_mass, &b_met_mass);
   fChain->SetBranchAddress("met_sumEt", &met_sumEt, &b_met_sumEt);
   fChain->SetBranchAddress("met_rawPt", &met_rawPt, &b_met_rawPt);
   fChain->SetBranchAddress("met_rawPhi", &met_rawPhi, &b_met_rawPhi);
   fChain->SetBranchAddress("met_rawSumEt", &met_rawSumEt, &b_met_rawSumEt);
   fChain->SetBranchAddress("met_genPt", &met_genPt, &b_met_genPt);
   fChain->SetBranchAddress("met_genPhi", &met_genPhi, &b_met_genPhi);
   fChain->SetBranchAddress("met_genEta", &met_genEta, &b_met_genEta);
   fChain->SetBranchAddress("met_shifted_MuonEnDown_pt", &met_shifted_MuonEnDown_pt, &b_met_shifted_MuonEnDown_pt);
   fChain->SetBranchAddress("met_shifted_MuonEnDown_phi", &met_shifted_MuonEnDown_phi, &b_met_shifted_MuonEnDown_phi);
   fChain->SetBranchAddress("met_shifted_MuonEnDown_sumEt", &met_shifted_MuonEnDown_sumEt, &b_met_shifted_MuonEnDown_sumEt);
   fChain->SetBranchAddress("met_shifted_ElectronEnUp_pt", &met_shifted_ElectronEnUp_pt, &b_met_shifted_ElectronEnUp_pt);
   fChain->SetBranchAddress("met_shifted_ElectronEnUp_phi", &met_shifted_ElectronEnUp_phi, &b_met_shifted_ElectronEnUp_phi);
   fChain->SetBranchAddress("met_shifted_ElectronEnUp_sumEt", &met_shifted_ElectronEnUp_sumEt, &b_met_shifted_ElectronEnUp_sumEt);
   fChain->SetBranchAddress("met_shifted_ElectronEnDown_pt", &met_shifted_ElectronEnDown_pt, &b_met_shifted_ElectronEnDown_pt);
   fChain->SetBranchAddress("met_shifted_ElectronEnDown_phi", &met_shifted_ElectronEnDown_phi, &b_met_shifted_ElectronEnDown_phi);
   fChain->SetBranchAddress("met_shifted_ElectronEnDown_sumEt", &met_shifted_ElectronEnDown_sumEt, &b_met_shifted_ElectronEnDown_sumEt);
   fChain->SetBranchAddress("met_shifted_TauEnDown_pt", &met_shifted_TauEnDown_pt, &b_met_shifted_TauEnDown_pt);
   fChain->SetBranchAddress("met_shifted_TauEnDown_phi", &met_shifted_TauEnDown_phi, &b_met_shifted_TauEnDown_phi);
   fChain->SetBranchAddress("met_shifted_TauEnDown_sumEt", &met_shifted_TauEnDown_sumEt, &b_met_shifted_TauEnDown_sumEt);
   fChain->SetBranchAddress("met_shifted_TauEnUp_pt", &met_shifted_TauEnUp_pt, &b_met_shifted_TauEnUp_pt);
   fChain->SetBranchAddress("met_shifted_TauEnUp_phi", &met_shifted_TauEnUp_phi, &b_met_shifted_TauEnUp_phi);
   fChain->SetBranchAddress("met_shifted_TauEnUp_sumEt", &met_shifted_TauEnUp_sumEt, &b_met_shifted_TauEnUp_sumEt);
   fChain->SetBranchAddress("met_shifted_UnclusteredEnUp_pt", &met_shifted_UnclusteredEnUp_pt, &b_met_shifted_UnclusteredEnUp_pt);
   fChain->SetBranchAddress("met_shifted_UnclusteredEnUp_phi", &met_shifted_UnclusteredEnUp_phi, &b_met_shifted_UnclusteredEnUp_phi);
   fChain->SetBranchAddress("met_shifted_UnclusteredEnUp_sumEt", &met_shifted_UnclusteredEnUp_sumEt, &b_met_shifted_UnclusteredEnUp_sumEt);
   fChain->SetBranchAddress("met_shifted_UnclusteredEnDown_pt", &met_shifted_UnclusteredEnDown_pt, &b_met_shifted_UnclusteredEnDown_pt);
   fChain->SetBranchAddress("met_shifted_UnclusteredEnDown_phi", &met_shifted_UnclusteredEnDown_phi, &b_met_shifted_UnclusteredEnDown_phi);
   fChain->SetBranchAddress("met_shifted_UnclusteredEnDown_sumEt", &met_shifted_UnclusteredEnDown_sumEt, &b_met_shifted_UnclusteredEnDown_sumEt);
   fChain->SetBranchAddress("met_shifted_JetEnUp_pt", &met_shifted_JetEnUp_pt, &b_met_shifted_JetEnUp_pt);
   fChain->SetBranchAddress("met_shifted_JetEnUp_phi", &met_shifted_JetEnUp_phi, &b_met_shifted_JetEnUp_phi);
   fChain->SetBranchAddress("met_shifted_JetEnUp_sumEt", &met_shifted_JetEnUp_sumEt, &b_met_shifted_JetEnUp_sumEt);
   fChain->SetBranchAddress("met_shifted_JetEnDown_pt", &met_shifted_JetEnDown_pt, &b_met_shifted_JetEnDown_pt);
   fChain->SetBranchAddress("met_shifted_JetEnDown_phi", &met_shifted_JetEnDown_phi, &b_met_shifted_JetEnDown_phi);
   fChain->SetBranchAddress("met_shifted_JetEnDown_sumEt", &met_shifted_JetEnDown_sumEt, &b_met_shifted_JetEnDown_sumEt);
   fChain->SetBranchAddress("met_shifted_JetResDown_pt", &met_shifted_JetResDown_pt, &b_met_shifted_JetResDown_pt);
   fChain->SetBranchAddress("met_shifted_JetResDown_phi", &met_shifted_JetResDown_phi, &b_met_shifted_JetResDown_phi);
   fChain->SetBranchAddress("met_shifted_JetResDown_sumEt", &met_shifted_JetResDown_sumEt, &b_met_shifted_JetResDown_sumEt);
   fChain->SetBranchAddress("softActivity_njets2", &softActivity_njets2, &b_softActivity_njets2);
   fChain->SetBranchAddress("softActivity_njets5", &softActivity_njets5, &b_softActivity_njets5);
   fChain->SetBranchAddress("softActivity_njets10", &softActivity_njets10, &b_softActivity_njets10);
   fChain->SetBranchAddress("softActivity_HT", &softActivity_HT, &b_softActivity_HT);
   fChain->SetBranchAddress("met_shifted_MuonEnUp_pt", &met_shifted_MuonEnUp_pt, &b_met_shifted_MuonEnUp_pt);
   fChain->SetBranchAddress("met_shifted_MuonEnUp_phi", &met_shifted_MuonEnUp_phi, &b_met_shifted_MuonEnUp_phi);
   fChain->SetBranchAddress("met_shifted_MuonEnUp_sumEt", &met_shifted_MuonEnUp_sumEt, &b_met_shifted_MuonEnUp_sumEt);
   fChain->SetBranchAddress("nGenBQuarkFromHafterISR", &nGenBQuarkFromHafterISR, &b_nGenBQuarkFromHafterISR);
   fChain->SetBranchAddress("GenBQuarkFromHafterISR_pdgId", &GenBQuarkFromHafterISR_pdgId, &b_GenBQuarkFromHafterISR_pdgId);
   fChain->SetBranchAddress("GenBQuarkFromHafterISR_pt", &GenBQuarkFromHafterISR_pt, &b_GenBQuarkFromHafterISR_pt);
   fChain->SetBranchAddress("GenBQuarkFromHafterISR_eta", &GenBQuarkFromHafterISR_eta, &b_GenBQuarkFromHafterISR_eta);
   fChain->SetBranchAddress("GenBQuarkFromHafterISR_phi", &GenBQuarkFromHafterISR_phi, &b_GenBQuarkFromHafterISR_phi);
   fChain->SetBranchAddress("GenBQuarkFromHafterISR_mass", &GenBQuarkFromHafterISR_mass, &b_GenBQuarkFromHafterISR_mass);
   fChain->SetBranchAddress("GenBQuarkFromHafterISR_charge", &GenBQuarkFromHafterISR_charge, &b_GenBQuarkFromHafterISR_charge);
   fChain->SetBranchAddress("GenBQuarkFromHafterISR_status", &GenBQuarkFromHafterISR_status, &b_GenBQuarkFromHafterISR_status);
   fChain->SetBranchAddress("nvLeptons", &nvLeptons, &b_nvLeptons);
   fChain->SetBranchAddress("vLeptons_dz", vLeptons_dz, &b_vLeptons_dz);
   fChain->SetBranchAddress("nLHE_weights_scale", &nLHE_weights_scale, &b_nLHE_weights_scale);
   fChain->SetBranchAddress("LHE_weights_scale_id", &LHE_weights_scale_id, &b_LHE_weights_scale_id);
   fChain->SetBranchAddress("LHE_weights_scale_wgt", &LHE_weights_scale_wgt, &b_LHE_weights_scale_wgt);
   fChain->SetBranchAddress("nGenBQuarkFromH", &nGenBQuarkFromH, &b_nGenBQuarkFromH);
   fChain->SetBranchAddress("GenBQuarkFromH_pdgId", GenBQuarkFromH_pdgId, &b_GenBQuarkFromH_pdgId);
   fChain->SetBranchAddress("GenBQuarkFromH_pt", GenBQuarkFromH_pt, &b_GenBQuarkFromH_pt);
   fChain->SetBranchAddress("GenBQuarkFromH_eta", GenBQuarkFromH_eta, &b_GenBQuarkFromH_eta);
   fChain->SetBranchAddress("GenBQuarkFromH_phi", GenBQuarkFromH_phi, &b_GenBQuarkFromH_phi);
   fChain->SetBranchAddress("GenBQuarkFromH_mass", GenBQuarkFromH_mass, &b_GenBQuarkFromH_mass);
   fChain->SetBranchAddress("GenBQuarkFromH_charge", GenBQuarkFromH_charge, &b_GenBQuarkFromH_charge);
   fChain->SetBranchAddress("GenBQuarkFromH_status", GenBQuarkFromH_status, &b_GenBQuarkFromH_status);
   fChain->SetBranchAddress("nGenHiggsSisters", &nGenHiggsSisters, &b_nGenHiggsSisters);
   fChain->SetBranchAddress("GenHiggsSisters_pdgId", GenHiggsSisters_pdgId, &b_GenHiggsSisters_pdgId);
   fChain->SetBranchAddress("GenHiggsSisters_pt", GenHiggsSisters_pt, &b_GenHiggsSisters_pt);
   fChain->SetBranchAddress("GenHiggsSisters_eta", GenHiggsSisters_eta, &b_GenHiggsSisters_eta);
   fChain->SetBranchAddress("GenHiggsSisters_phi", GenHiggsSisters_phi, &b_GenHiggsSisters_phi);
   fChain->SetBranchAddress("GenHiggsSisters_mass", GenHiggsSisters_mass, &b_GenHiggsSisters_mass);
   fChain->SetBranchAddress("GenHiggsSisters_charge", GenHiggsSisters_charge, &b_GenHiggsSisters_charge);
   fChain->SetBranchAddress("GenHiggsSisters_status", GenHiggsSisters_status, &b_GenHiggsSisters_status);
   fChain->SetBranchAddress("nsoftActivityVHJets", &nsoftActivityVHJets, &b_nsoftActivityVHJets);
   fChain->SetBranchAddress("softActivityVHJets_pt", softActivityVHJets_pt, &b_softActivityVHJets_pt);
   fChain->SetBranchAddress("softActivityVHJets_eta", softActivityVHJets_eta, &b_softActivityVHJets_eta);
   fChain->SetBranchAddress("softActivityVHJets_phi", softActivityVHJets_phi, &b_softActivityVHJets_phi);
   fChain->SetBranchAddress("softActivityVHJets_mass", softActivityVHJets_mass, &b_softActivityVHJets_mass);
   fChain->SetBranchAddress("nselLeptons", &nselLeptons, &b_nselLeptons);
   fChain->SetBranchAddress("selLeptons_charge", selLeptons_charge, &b_selLeptons_charge);
   fChain->SetBranchAddress("selLeptons_tightId", selLeptons_tightId, &b_selLeptons_tightId);
   fChain->SetBranchAddress("selLeptons_relIso03", selLeptons_relIso03, &b_selLeptons_relIso03);
   fChain->SetBranchAddress("selLeptons_pdgId", selLeptons_pdgId, &b_selLeptons_pdgId);
   fChain->SetBranchAddress("selLeptons_pfRelIso03", selLeptons_pfRelIso03, &b_selLeptons_pfRelIso03);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("Jet_id", Jet_id, &b_Jet_id);
   fChain->SetBranchAddress("Jet_puId", Jet_puId, &b_Jet_puId);
   fChain->SetBranchAddress("Jet_btagCSV", Jet_btagCSV, &b_Jet_btagCSV);
   fChain->SetBranchAddress("Jet_hadronFlavour", Jet_hadronFlavour, &b_Jet_hadronFlavour);
   fChain->SetBranchAddress("Jet_corr_JECUp", Jet_corr_JECUp, &b_Jet_corr_JECUp);
   fChain->SetBranchAddress("Jet_corr_JECDown", Jet_corr_JECDown, &b_Jet_corr_JECDown);
   fChain->SetBranchAddress("Jet_corr", Jet_corr, &b_Jet_corr);
   fChain->SetBranchAddress("Jet_corr_JERUp", Jet_corr_JERUp, &b_Jet_corr_JERUp);
   fChain->SetBranchAddress("Jet_corr_JERDown", Jet_corr_JERDown, &b_Jet_corr_JERDown);
   fChain->SetBranchAddress("Jet_corr_JER", Jet_corr_JER, &b_Jet_corr_JER);
   fChain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_mass", Jet_mass, &b_Jet_mass);
   fChain->SetBranchAddress("Jet_leadTrackPt", Jet_leadTrackPt, &b_Jet_leadTrackPt);
   fChain->SetBranchAddress("Jet_qgl", Jet_qgl, &b_Jet_qgl);
   fChain->SetBranchAddress("Jet_ptd", Jet_ptd, &b_Jet_ptd);
   fChain->SetBranchAddress("Jet_axis2", Jet_axis2, &b_Jet_axis2);
   fChain->SetBranchAddress("Jet_mult", Jet_mult, &b_Jet_mult);
   fChain->SetBranchAddress("Jet_blike_VBF", Jet_blike_VBF, &b_Jet_blike_VBF);
   fChain->SetBranchAddress("Jet_pt_reg", Jet_pt_reg, &b_Jet_pt_reg);
   fChain->SetBranchAddress("Jet_pt_regVBF", Jet_pt_regVBF, &b_Jet_pt_regVBF);
   fChain->SetBranchAddress("Jet_pt_reg_corrJECUp", Jet_pt_reg_corrJECUp, &b_Jet_pt_reg_corrJECUp);
   fChain->SetBranchAddress("Jet_pt_regVBF_corrJECUp", Jet_pt_regVBF_corrJECUp, &b_Jet_pt_regVBF_corrJECUp);
   fChain->SetBranchAddress("Jet_pt_reg_corrJECDown", Jet_pt_reg_corrJECDown, &b_Jet_pt_reg_corrJECDown);
   fChain->SetBranchAddress("Jet_pt_regVBF_corrJECDown", Jet_pt_regVBF_corrJECDown, &b_Jet_pt_regVBF_corrJECDown);
   fChain->SetBranchAddress("Jet_pt_reg_corrJERUp", Jet_pt_reg_corrJERUp, &b_Jet_pt_reg_corrJERUp);
   fChain->SetBranchAddress("Jet_pt_regVBF_corrJERUp", Jet_pt_regVBF_corrJERUp, &b_Jet_pt_regVBF_corrJERUp);
   fChain->SetBranchAddress("Jet_pt_reg_corrJERDown", Jet_pt_reg_corrJERDown, &b_Jet_pt_reg_corrJERDown);
   fChain->SetBranchAddress("Jet_pt_regVBF_corrJERDown", Jet_pt_regVBF_corrJERDown, &b_Jet_pt_regVBF_corrJERDown);
   fChain->SetBranchAddress("nLHE_weights_pdf", &nLHE_weights_pdf, &b_nLHE_weights_pdf);
   fChain->SetBranchAddress("LHE_weights_pdf_id", &LHE_weights_pdf_id, &b_LHE_weights_pdf_id);
   fChain->SetBranchAddress("LHE_weights_pdf_wgt", &LHE_weights_pdf_wgt, &b_LHE_weights_pdf_wgt);
   fChain->SetBranchAddress("nsoftActivityJets", &nsoftActivityJets, &b_nsoftActivityJets);
   fChain->SetBranchAddress("softActivityJets_pt", softActivityJets_pt, &b_softActivityJets_pt);
   fChain->SetBranchAddress("softActivityJets_eta", softActivityJets_eta, &b_softActivityJets_eta);
   fChain->SetBranchAddress("softActivityJets_phi", softActivityJets_phi, &b_softActivityJets_phi);
   fChain->SetBranchAddress("softActivityJets_mass", softActivityJets_mass, &b_softActivityJets_mass);
   Notify();
}

Bool_t CreateTree_tmva::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void CreateTree_tmva::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t CreateTree_tmva::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef CreateTree_tmva_cxx
