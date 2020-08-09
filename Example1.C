/*
This macro shows how to access the particle-level reference for reconstructed objects.
It is also shown how to loop over the jet constituents.

root -l examples/Example1.C'("delphes_output.root")'
*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#else
class ExRootTreeReader;
class ExRootResult;
#endif


#include "TPaveText.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TH1.h"
#include "TClonesArray.h"
#include <stdio.h>
#include <iostream>

using namespace std;


//------------------------------------------------------------------------------
bool cmp(pair<int, double>& a,
         pair<int, double>& b)
{
    return a.second < b.second;
}

struct TestPlots
{

  TH1 *fParticlePT;
  TH1 *fParticleEta;
  TH1 *fParticleWM;
  TH1 *fParticleTopM;

  TH1 *fElectronDeltaPT;
  TH1 *fElectronDeltaEta;

  TH1 *fPhotonDeltaPT;
  TH1 *fPhotonDeltaEta;
  TH1 *fPhotonDeltaE;

  TH1 *fMuonDeltaPT;
  TH1 *fMuonDeltaEta;

  TH1 *fJetDeltaPT;
  TH1 *fJetPT;
  TH1 *fJetEta;
  TH1 *fNJets;
  TH1 *fNBJets;
  
  TH1 *fminDRjj;
  TH1 *fDRbb;


    
  TH1 *fminDRqq;
 
  TH1 *flquarkPT;
  TH1 *flmatchPT;
  TH1 *fbquarkPT;
  TH1 *fbmatchPT;
  TH1 *flquarkEta;
  TH1 *flmatchEta;
  TH1 *fbquarkEta;
  TH1 *fbmatchEta;


    TH1 *fmatchTop2_minDRjj;
    TH1 *fmatchTop2_minDRqq;
    TH1 *fmatchTop2_truthDRbb;
    TH1 *fmatchTop2_Jet4PT;
    TH1 *fmatchTop2_Jet5PT;
    TH1 *fmatchTop2_Jet6PT;

    TH1 *fmatchTop1_minDRjj;
    TH1 *fmatchTop1_minDRqq;
    TH1 *fmatchTop1_truthDRbb;
    TH1 *fmatchTop1_Jet4PT;
    TH1 *fmatchTop1_Jet5PT;
    TH1 *fmatchTop1_Jet6PT;

    TH1 *fmatchTop0_minDRqq;
    TH1 *fmatchTop0_minDRjj;
    TH1 *fmatchTop0_truthDRbb;
    TH1 *fmatchTop0_Jet4PT;
    TH1 *fmatchTop0_Jet5PT;
    TH1 *fmatchTop0_Jet6PT;


  TH1 *fmatchWJetM;
  TH1 *fmatchTopJetM;
  TH1 *fmatchTopJetNum;

};

//------------------------------------------------------------------------------

void BookHistograms(ExRootResult *result, TestPlots *plots)
{
  TLegend *legend;
  TPaveText *comment;

  plots->fElectronDeltaPT = result->AddHist1D(
    "electron_delta_pt", "(p_{T}^{particle} - p_{T}^{electron})/p_{T}^{particle}",
    "(p_{T}^{particle} - p_{T}^{electron})/p_{T}^{particle}", "number of electrons",
    100, -0.1, 0.1);

  plots->fElectronDeltaEta = result->AddHist1D(
    "electron_delta_eta", "(#eta^{particle} - #eta^{electron})/#eta^{particle}",
    "(#eta^{particle} - #eta^{electron})/#eta^{particle}", "number of electrons",
    100, -0.1, 0.1);

  plots->fPhotonDeltaPT = result->AddHist1D(
    "photon_delta_pt", "(p_{T}^{particle} - p_{T}^{photon})/p_{T}^{particle}",
    "(p_{T}^{particle} - p_{T}^{photon})/p_{T}^{particle}", "number of photons",
    100, -0.1, 0.1);

  plots->fPhotonDeltaEta = result->AddHist1D(
    "photon_delta_eta", "(#eta^{particle} - #eta^{photon})/#eta^{particle}",
    "(#eta^{particle} - #eta^{photon})/#eta^{particle}", "number of photons",
    100, -0.1, 0.1);

  plots->fPhotonDeltaE = result->AddHist1D(
    "photon_delta_energy", "(E^{particle} - E^{photon})/E^{particle}",
    "(E^{particle} - E^{photon})/E^{particle}", "number of photons",
    100, -0.1, 0.1);

  plots->fMuonDeltaPT = result->AddHist1D(
    "muon_delta_pt", "(p_{T}^{particle} - p_{T}^{muon})/p_{T}^{particle}",
    "(p_{T}^{particle} - p_{T}^{muon})/p_{T}^{particle}", "number of muons",
    100, -0.1, 0.1);

  plots->fMuonDeltaEta = result->AddHist1D(
    "muon_delta_eta", "(#eta^{particle} - #eta^{muon})/#eta^{particle}",
    "(#eta^{particle} - #eta^{muon})/#eta^{particle}", "number of muons",
    100, -0.1, 0.1);

  plots->fJetDeltaPT = result->AddHist1D(
    "jet_delta_pt", "(p_{T}^{jet} - p_{T}^{constituents})/p_{T}^{jet}",
    "(p_{T}^{jet} - p_{T}^{constituents})/p_{T}^{jet}", "number of jets",
    100, -1.0e-1, 1.0e-1);

  plots->fJetPT = result->AddHist1D(
      "jet_pt", "p_{T}^{jet} ",
      "p_{T}^{jet}", "number of jets",
      100, 0.0, 300.0);


    
  plots->fJetEta = result->AddHist1D(
      "jet_eta", "jet_eta",
      "eta", "number",
      40,-4,4);
    
    plots->fNJets = result->AddHist1D(
    "NJets", "Njet",
    "NJets", "number",
    15,-0.5,14.5);

    plots->fNBJets = result->AddHist1D(
    "NBJets", "NBjet",
    "NBJets", "number",
    15,-0.5,14.5);

    
    plots->fParticlePT = result->AddHist1D(
    "particle_pt", "",
    "pt [GeV]", "number",
    100,0.0,400.0);

    plots->fParticleEta = result->AddHist1D(
    "particle_eta", "",
    "eta", "number",
    40,-4,4);

    plots->fParticleWM = result->AddHist1D(
      "particle_W_M", "",
      "Mass [GeV]", "number",
      100,40.0,100.0);

    plots->fParticleTopM = result->AddHist1D(
      "particle_Top_M", "",
      "Mass [GeV]", "number",
      80,60.0,300.0);

    plots->fminDRjj = result->AddHist1D(
      "minDRjj", "",
      "#Delta R", "number",
      100,0.0,3.2);

    plots->fDRbb = result->AddHist1D(
      "DRbb", "",
      "#Delta R", "number",
      100,0.0,4.0);

    plots->fminDRqq = result->AddHist1D(
      "minDRjj", "",
      "#Delta R", "number",
      100,0.0,3.2);

    plots->fmatchTop2_minDRjj = result->AddHist1D(
      "matchTop2_minDRjj", "",
      "#Delta R", "number",
      100,0.0,3.2);

    plots->fmatchTop1_minDRjj = result->AddHist1D(
      "matchTop1_minDRjj", "",
      "#Delta R", "number",
      100,0.0,3.2);

    plots->fmatchTop0_minDRjj = result->AddHist1D(
      "matchTop0_minDRjj", "",
      "#Delta R", "number",
      100,0.0,3.2);

    plots->fmatchTop2_minDRqq = result->AddHist1D(
      "matchTop2_minDRqq", "",
      "#Delta R", "number",
      100,0.0,3.2);

    plots->fmatchTop1_minDRqq = result->AddHist1D(
      "matchTop1_minDRqq", "",
      "#Delta R", "number",
      100,0.0,3.2);
    
    plots->fmatchTop0_minDRqq = result->AddHist1D(
      "matchTop0_minDRqq", "",
      "#Delta R", "number",
      100,0.0,3.2);

    plots->fmatchTop2_truthDRbb = result->AddHist1D(
      "matchTop2_truthDRqq", "",
      "#Delta R", "number",
      100,0.0,3.2);

    plots->fmatchTop1_truthDRbb = result->AddHist1D(
      "matchTop1_truthDRqq", "",
      "#Delta R", "number",
      100,0.0,3.2);

    plots->fmatchTop0_truthDRbb = result->AddHist1D(
      "matchTop0_truthDRqq", "",
      "#Delta R", "number",
      100,0.0,3.2);


    plots->flquarkPT = result->AddHist1D(
      "lquark_PT", "",
      "P_{T} [GeV]", "number",
      100,0.0,300.0);

    plots->flmatchPT = result->AddHist1D(
      "lmatch_PT", "",
      "P_{T} [GeV]", "number",
      100,0.0,300.0);

    plots->fbquarkPT = result->AddHist1D(
      "bquark_PT", "",
      "P_{T} [GeV]", "number",
      100,0.0,300.0);

    plots->fbmatchPT = result->AddHist1D(
      "bmatch_PT", "",
      "P_{T} [GeV]", "number",
      100,0.0,300.0);

    plots->flquarkEta= result->AddHist1D(
      "lquark_Eta", "",
      "#eta", "number",
      60,-3,3);

    plots->flmatchEta= result->AddHist1D(
      "lmatch_Eta", "",
      "#eta", "number",
      60,-3,3);

    plots->fbquarkEta= result->AddHist1D(
      "bquark_Eta", "",
      "#eta", "number",
      60,-3,3);

    plots->fbmatchEta= result->AddHist1D(
      "bmatch_Eta", "",
      "#eta", "number",
      60,-3,3);

    plots->fmatchWJetM = result->AddHist1D(
      "matchJet_W_M", "",
      "Mass [GeV]", "number",
      75,0.0,150.0);

    plots->fmatchTopJetM = result->AddHist1D(
      "matchJet_Top_M", "",
      "Mass [GeV]", "number",
      80,40.0,300.0);

    plots->fmatchTopJetNum = result->AddHist1D(
    "matchJets", "",
    "matchTopJets", "number",
    4,-0.5,3.5);

    plots->fmatchTop2_Jet4PT = result->AddHist1D(
        "matchTop2_jet4_pt", "",
        "p_{T}^{4^{th} jet}", "number of jets",
        100, 0.0, 400.0);
    plots->fmatchTop2_Jet5PT = result->AddHist1D(
        "matchTop2_jet5_pt", "",
        "p_{T}^{5^{th} jet}", "number of jets",
        100, 0.0, 400.0);
    plots->fmatchTop2_Jet6PT = result->AddHist1D(
        "matchTop2_jet6_pt", "",
        "p_{T}^{6^{th} jet}", "number of jets",
        100, 0.0, 400.0);

    plots->fmatchTop1_Jet4PT = result->AddHist1D(
        "matchTop1_jet4_pt", "",
        "p_{T}^{4^{th} jet}", "number of jets",
        100, 0.0, 400.0);
    plots->fmatchTop1_Jet5PT = result->AddHist1D(
        "matchTop1_jet5_pt", "",
        "p_{T}^{5^{th} jet}", "number of jets",
        100, 0.0, 400.0);
    plots->fmatchTop1_Jet6PT = result->AddHist1D(
        "matchTop1_jet6_pt", "",
        "p_{T}^{6^{th} jet}", "number of jets",
        100, 0.0, 400.0);

    plots->fmatchTop0_Jet4PT = result->AddHist1D(
        "matchTop0_jet4_pt", "",
        "p_{T}^{4^{th} jet}", "number of jets",
        100, 0.0, 400.0);
    plots->fmatchTop0_Jet5PT = result->AddHist1D(
        "matchTop0_jet5_pt", "",
        "p_{T}^{5^{th} jet}", "number of jets",
        100, 0.0, 400.0);
    plots->fmatchTop0_Jet6PT = result->AddHist1D(
        "matchTop0_jet6_pt", "",
        "p_{T}^{6^{th} jet}", "number of jets",
        100, 0.0, 400.0);

}

//https://gitlab.cern.ch/atlas/athena/-/blob/master/PhysicsAnalysis/TopPhys/xAOD/TopPartons/Root/CalcTopPartonHistory.cxx

int findAfterFSR(int index, const TClonesArray *truthParticles) {
   bool isAfterFSR(false);

   const GenParticle* particle = (GenParticle*) truthParticles->At(index);
   const int particle_ID = particle->PID;
   int forLoop  = 0;
   while(!isAfterFSR){

           forLoop  = 0;
           for (size_t j=particle->D1; j<= particle->D2; j++ ) {
               const GenParticle* tmp_children = (GenParticle*) truthParticles->At(j);
               if (tmp_children && tmp_children->PID==particle_ID){
                   particle = tmp_children;
                   index = j;
                   forLoop++;
                   break;
               }//if
           }//for
           
           if (forLoop == 0)       isAfterFSR = true;
   }//while

   return index;
}


bool hasParticleIdenticalParent(const GenParticle* particle, const TClonesArray *truthParticles) {
   bool skipit(false);
    for (size_t i=particle->M1; i<=particle->M2; i++ ) {
        const GenParticle* parent = (GenParticle*) truthParticles->At(i);
        if (parent && parent->PID==particle->PID){
               skipit=true;
               break;
        }//if
   }//for
   return skipit;
}

//
bool topWb( const TClonesArray *truthParticles,
            int start, int& t_beforeFSR, int& t_afterFSR, int& W,
                  int& b, int& Wdecay1, int& Wdecay2){
  
  bool hasT            = false;
  bool hasW            = false;
  bool hasB            = false;
  bool hasWdecayProd1     = false;
  bool hasWdecayProd2     = false;
  

  for (size_t idx =0 ;idx <truthParticles->GetEntries();idx++){
          
      const GenParticle* particle = (GenParticle*) truthParticles->At(idx);

      if (particle->PID != start)    continue;
      if (hasParticleIdenticalParent(particle, truthParticles)) continue; // kepping only top before FSR

 
      t_beforeFSR = idx; // top before FSR
      hasT = true;
          
      // demanding the last tops after FSR
      idx = findAfterFSR(idx, truthParticles);
      t_afterFSR = idx; // top after FSR
      particle = (GenParticle*) truthParticles->At(idx);

      for (size_t k=particle->D1; k <= particle->D2; k++) {
          const GenParticle* topChildren = (GenParticle*) truthParticles->At(k);

          if (fabs(topChildren->PID) == 24){
          
              W = k;  // W boson after FSR
              hasW = true;
      
              // demanding the last W after FSR
              W = findAfterFSR(W, truthParticles);
              topChildren = (GenParticle*) truthParticles->At(W);
 
              for (size_t q = topChildren->D1; q <= topChildren->D2; ++q) {
                  const GenParticle* WChildren = (GenParticle*) truthParticles->At(q);
                  if (fabs(WChildren->PID)<17){
                      if (WChildren->PID>0){
                          Wdecay1 = q;
                          hasWdecayProd1 = true;
                      }else{
                          Wdecay2 = q;
                          hasWdecayProd2 = true;
                      }//else
                  }//if
              }//for

          } else if (fabs(topChildren->PID) == 5) {
              b = k;
              hasB = true;
          } //else if
  
      } //for (size_t k=0; k < particle->nChildren(); k++)

      if (hasT && hasW && hasB && hasWdecayProd1 && hasWdecayProd2)    return true;
  } //for (const xAOD::TruthParticle* particle : *truthParticles)
  
  return false;
                
}



//------------------------------------------------------------------------------

void AnalyseEvents(ExRootTreeReader *treeReader, TestPlots *plots, const char *outputFile, bool isGenJets, bool isDEBUG)
{
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
  TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");

  TClonesArray *branchJet = 0;
  if(isGenJets)
      branchJet = treeReader->UseBranch("GenJet");
  else
      branchJet = treeReader->UseBranch("Jet");
  cout << "JetBranch= "<< branchJet->GetName()<<endl;


  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Chain contains " << allEntries << " events" << endl;

  GenParticle *particle;
  Electron *electron;
  Photon *photon;
  Muon *muon;

  Track *track;
  Tower *tower;

  Jet *jet;
  TObject *object;

  TLorentzVector momentum;

  Float_t Eem, Ehad;
  Bool_t skip;

  Long64_t entry;

  Int_t i, j, pdgCode;

    
    //Output Tree
    TFile *fout = new TFile(outputFile,"recreate");
    TTree *event =new TTree("event","Jet-parton assignment");

    int NmatchTop;
    event->Branch("NmatchTop", &NmatchTop);
    vector<double> jet_pt, jet_eta, jet_phi, jet_mass;
    vector<int> jet_btag, jet_barcode, jet_parton_index;

    event->Branch("jet_parton_index",&jet_parton_index);
    event->Branch("jet_barcode",&jet_barcode);
    event->Branch("jet_pt",&jet_pt);
    event->Branch("jet_eta",&jet_eta);
    event->Branch("jet_phi",&jet_phi);
    event->Branch("jet_mass",&jet_mass);
    event->Branch("jet_btag",&jet_btag);

    vector<double> parton_pt, parton_eta, parton_phi, parton_mass;
    vector<int> parton_pdgid, parton_barcode, parton_jet_index;

    event->Branch("parton_jet_index",&parton_jet_index);
    event->Branch("parton_barcode",&parton_barcode);
    event->Branch("parton_pt",&parton_pt);
    event->Branch("parton_eta",&parton_eta);
    event->Branch("parton_phi",&parton_phi);
    event->Branch("parton_mass",&parton_mass);
    event->Branch("parton_pdgid",&parton_pdgid);

    gDirectory->cd();
  // Loop over all events
  //  allEntries=1000;

  int CountNJets=0, CountNBJets=0;

  for(entry = 0; entry < allEntries; ++entry)
  {

      if(entry%10000 ==0) cout <<"Process "<< entry <<endl;

    // Load selected branches with data from specified event
      treeReader->ReadEntry(entry);

      //Reset tree entries
      NmatchTop=0;
      jet_pt.clear(); jet_eta.clear(); jet_phi.clear(); jet_mass.clear();
      jet_btag.clear(); jet_barcode.clear(); jet_parton_index.clear();
      parton_pt.clear(); parton_eta.clear(); parton_phi.clear(); parton_mass.clear();
      parton_pdgid.clear(); parton_barcode.clear(); parton_jet_index.clear();

    // Get Final State Particle
 
      int t_before_idx, t_after_idx, Wp_idx, b_idx, WpDecay1_idx, WpDecay2_idx;
      int tbar_before_idx, tbar_after_idx, Wm_idx, bbar_idx, WmDecay1_idx, WmDecay2_idx;
      bool event_top    = topWb(branchParticle, 6, t_before_idx, t_after_idx, Wp_idx, b_idx, WpDecay1_idx, WpDecay2_idx);
      bool event_topbar = topWb(branchParticle,-6, tbar_before_idx, tbar_after_idx, Wm_idx, bbar_idx, WmDecay1_idx, WmDecay2_idx);

      if(!event_top || !event_topbar) {
          cout <<"Error Event "<< entry << " failed finding top truth =" << event_top <<" topbar truth=" << event_topbar<<endl;
          continue;
      }
          
      //      cout << "Top "<<t_before_idx <<" "<<t_after_idx <<" "<<Wp_idx <<" "<<b_idx <<" "<<WpDecay1_idx<<" "<< WpDecay2_idx<<endl;
      //      cout << "Top "<<tbar_before_idx <<" "<<tbar_after_idx <<" "<<Wm_idx <<" "<<bbar_idx <<" "<<WmDecay1_idx<<" "<< WmDecay2_idx<<endl;

      GenParticle *w1=0, *b1=0, *q1=0, *q2=0;
      GenParticle *w2=0, *b2=0, *q3=0, *q4=0;
      w1 = (GenParticle*) branchParticle->At(Wp_idx);
      b1 = (GenParticle*) branchParticle->At(b_idx);
      q1 = (GenParticle*) branchParticle->At(WpDecay1_idx);
      q2 = (GenParticle*) branchParticle->At(WpDecay2_idx);

      w2 = (GenParticle*) branchParticle->At(Wm_idx);
      b2 = (GenParticle*) branchParticle->At(bbar_idx);
      q3 = (GenParticle*) branchParticle->At(WmDecay1_idx);
      q4 = (GenParticle*) branchParticle->At(WmDecay2_idx);


      vector<GenParticle*> partonList;
      partonList.push_back(b1);
      partonList.push_back(q1);
      partonList.push_back(q2);
      partonList.push_back(b2);
      partonList.push_back(q3);
      partonList.push_back(q4);

      for(int i=0;i<partonList.size();i++){
          GenParticle* parton = partonList[i];
          plots->fParticlePT->Fill(parton->PT);
          plots->fParticleEta->Fill(parton->Eta);
      }
      TLorentzVector w1P4 = q1->P4()+q2->P4();
      TLorentzVector w2P4 = q3->P4()+q4->P4();
      TLorentzVector t1P4 = b1->P4()+w1P4;
      TLorentzVector t2P4 = b2->P4()+w2P4;

      plots->fParticleWM->Fill(w1P4.M());
      plots->fParticleWM->Fill(w2P4.M());
      plots->fParticleTopM->Fill(t1P4.M());
      plots->fParticleTopM->Fill(t2P4.M());


    // Loop over all electrons in event
    for(i = 0; i < branchElectron->GetEntriesFast(); ++i)
    {
      electron = (Electron*) branchElectron->At(i);
      particle = (GenParticle*) electron->Particle.GetObject();

      plots->fElectronDeltaPT->Fill((particle->PT - electron->PT)/particle->PT);
      plots->fElectronDeltaEta->Fill((particle->Eta - electron->Eta)/particle->Eta);
    }

    // Loop over all photons in event
    for(i = 0; i < branchPhoton->GetEntriesFast(); ++i)
    {
      photon = (Photon*) branchPhoton->At(i);

      // skip photons with references to multiple particles
      if(photon->Particles.GetEntriesFast() != 1) continue;

      particle = (GenParticle*) photon->Particles.At(0);
      plots->fPhotonDeltaPT->Fill((particle->PT - photon->PT)/particle->PT);
      plots->fPhotonDeltaEta->Fill((particle->Eta - photon->Eta)/particle->Eta);
      plots->fPhotonDeltaE->Fill((particle->E - photon->E)/particle->E);
          
    }

    // Loop over all muons in event
    for(i = 0; i < branchMuon->GetEntriesFast(); ++i)
    {
      muon = (Muon*) branchMuon->At(i);
      particle = (GenParticle*) muon->Particle.GetObject();

      plots->fMuonDeltaPT->Fill((particle->PT - muon->PT)/particle->PT);
      plots->fMuonDeltaEta->Fill((particle->Eta - muon->Eta)/particle->Eta);
    }

    // cout << "--  New event -- " << endl;

      //==========================
      //Jet Selection
      //==========================
    int Njet40=0;
    int Njet50=0;
    vector<Jet*> jetList;
    vector<Jet*> bjetList;  
    // Loop over all jets in event
    for(i = 0; i < branchJet->GetEntriesFast(); ++i)
    {
      jet = (Jet*) branchJet->At(i);

      plots->fJetPT->Fill(jet->PT);
      plots->fJetEta->Fill(jet->Eta);

      if(jet->PT<25 || fabs(jet->Eta)>2.4) continue;
      jetList.push_back(jet);

      if(jet->PT>40) Njet40++;
      if(jet->PT>50) Njet50++;
      
      //if(jet->BTag)
      if(jet->Flavor==5)
      	bjetList.push_back(jet);
        
      momentum.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);

      // cout<<"Looping over jet constituents. Jet pt: "<<jet->PT<<", eta: "<<jet->Eta<<", phi: "<<jet->Phi<<endl;

        
      // Loop over all jet's constituents
      for(j = 0; j < jet->Constituents.GetEntriesFast(); ++j)
      {
        object = jet->Constituents.At(j);

        // Check if the constituent is accessible
        if(object == 0) continue;

        if(object->IsA() == GenParticle::Class())
        {
          particle = (GenParticle*) object;
          // cout << "    GenPart pt: " << particle->PT << ", eta: " << particle->Eta << ", phi: " << particle->Phi << endl;
          momentum += particle->P4();
        }
        else if(object->IsA() == Track::Class())
        {
          track = (Track*) object;
          // cout << "    Track pt: " << track->PT << ", eta: " << track->Eta << ", phi: " << track->Phi << endl;
          momentum += track->P4();
        }
        else if(object->IsA() == Tower::Class())
        {
          tower = (Tower*) object;
          // cout << "    Tower pt: " << tower->ET << ", eta: " << tower->Eta << ", phi: " << tower->Phi << endl;
          momentum += tower->P4();
        }
      }
      plots->fJetDeltaPT->Fill((jet->PT - momentum.Pt())/jet->PT);
    }//jet loop
 
      //=============================================
      // Event Selection
      //=============================================

      //if(Njet50<6) continue;

      if(jetList.size()<6) continue;
      CountNJets++;
      if(bjetList.size()<2) continue;
      CountNBJets++;


      double minDRjj=100;
      for(int i=0;i<jetList.size();i++){
         TLorentzVector iP4 = jetList[i]->P4();
         for(int j=i+1;j<jetList.size();j++){
           TLorentzVector jP4 = jetList[j]->P4();
           double dR= jP4.DeltaR(iP4);
           if(dR<minDRjj)
              minDRjj= dR;
         }
      }
      //if(minDRjj<0.5) continue;


      plots->fNJets->Fill(jetList.size());
      plots->fNBJets->Fill(bjetList.size());
      plots->fminDRjj->Fill( minDRjj);

      plots->fDRbb->Fill ( bjetList[0]->P4().DeltaR(bjetList[1]->P4()) );

      //parton performance
      double minDRqq=4.0;
      for(size_t i=0;i<partonList.size();i++){
          GenParticle* iparticle = (GenParticle*) partonList[i];
          for(size_t j=i+1;j<partonList.size();j++){
               GenParticle* jparticle = (GenParticle*) partonList[j];
               double tmpDR= iparticle->P4().DeltaR(jparticle->P4());
               if(tmpDR<minDRqq) minDRqq=tmpDR;  
          } 
      }
      plots->fminDRqq->Fill( minDRqq );

      //Parton-jet matching performance
      for(size_t i=0;i<partonList.size();i++){
          GenParticle* particle = (GenParticle*) partonList[i];
          bool isMatch=false;
          for(size_t j=0;j< jetList.size();j++){
              Jet* jet = (Jet*) jetList[j];

              double tmpDR= jet->P4().DeltaR(particle->P4());
              if(tmpDR>0.4)continue;
              isMatch=true;
              break;
          }
          if(i==2 || i==5){

            plots->fbquarkPT ->Fill(particle->PT);
            plots->fbquarkEta->Fill(particle->Eta);
            if(isMatch){
              plots->fbmatchPT ->Fill(particle->PT);
              plots->fbmatchEta->Fill(particle->Eta);
            }
          }else{
            plots->flquarkPT ->Fill(particle->PT);
            plots->flquarkEta->Fill(particle->Eta);
            if(isMatch){
              plots->flmatchPT ->Fill(particle->PT);
              plots->flmatchEta->Fill(particle->Eta);
            }
          }
      }

      if(isDEBUG)
          cout << "Event " << entry <<endl;

      //=============================================
      // matching jet and partons based on dR
      //
      vector< pair<double, int> > dRList;
      int jNum= (int)jetList.size();
      int pNum= (int)partonList.size();

      for(int i=0;i<jNum;i++){
          TLorentzVector jetP4 = jetList[i]->P4();
          for(int j=0;j< pNum;j++){
              TLorentzVector parP4 = partonList[j]->P4();
              dRList.push_back( make_pair(parP4.DeltaR(jetP4), i*pNum+j ) );
        }
      }
      sort(dRList.begin(), dRList.end());
      
      //
      // matching based on minimum DeltaR
      //
      vector<vector<bool> > match( jNum , vector<bool> (pNum, false));
      for(int i=0;i< dRList.size();i++){
          if(dRList[i].first >0.3) break;
          int jdx = dRList[i].second/pNum;
          int pdx = dRList[i].second%pNum;
          
          bool hasMatched=false;
          for(int j=0;j<pNum;j++)
              if(match[jdx][j]) hasMatched=true;
          for(int j=0;j<jNum;j++)
              if(match[j][pdx]) hasMatched=true;

          if(hasMatched==false) match[jdx][pdx]=true;
           if(isDEBUG)
               cout << "DEBUG dR(jet "<< jdx <<", parton "<< pdx <<" )= "<< dRList[i].first <<" match "<< match[jdx][pdx] <<endl;
      }
      

      //Create 1D list for i) Parton to Jet matching index ii) Jet to Parton matching index
      vector<int> matchJetIndex (pNum, -1);
      vector<int> matchPartonIndex (jNum, -1);
      for(int j=0;j< pNum;j++){
          for(int i=0;i<jNum;i++){
              if(match[i][j]){
                  matchPartonIndex[i]=j;
                  matchJetIndex[j]=i;
                  if(isDEBUG)
                      cout <<"match jet "<<i <<" to parton "<< j <<endl;
              }
          }
      }

      //Jet Barcode
      //b1 q1 q2 b2 q3 q4
      int barcode[6]={0b100010, 0b101000, 0b101000, 0b010001, 0b010100, 0b010100};
      vector<int> matchBarcode (jNum, -1);
      for(int i=0;i<pNum;i++){
          int jdx= matchJetIndex[i];
          if(jdx>-1){
              matchBarcode[jdx]= barcode[i];
          }
      }
      //======================================
      // Plot matched top candidates
      //
            
      if(matchJetIndex[0]>-1 && matchJetIndex[1]>-1 && matchJetIndex[2]>-1){
          NmatchTop++;
          Jet *b1_jet = jetList[matchJetIndex[0]];
          Jet *q1_jet = jetList[matchJetIndex[1]];
          Jet *q2_jet = jetList[matchJetIndex[2]];

          TLorentzVector wP4 = q1_jet->P4()+q2_jet->P4();
          TLorentzVector tP4 = b1_jet->P4()+wP4;

          plots->fmatchWJetM->Fill(( wP4.M()));
          plots->fmatchTopJetM->Fill(( tP4.M()));
      }
      if(matchJetIndex[3]>-1 && matchJetIndex[4]>-1 && matchJetIndex[5]>-1){
          NmatchTop++;
          Jet *b2_jet = jetList[matchJetIndex[3]];
          Jet *q3_jet = jetList[matchJetIndex[4]];
          Jet *q4_jet = jetList[matchJetIndex[5]];

          TLorentzVector wP4 = q3_jet->P4()+q4_jet->P4();
          TLorentzVector tP4 = b2_jet->P4()+wP4;

          plots->fmatchWJetM->Fill(( wP4.M()));
          plots->fmatchTopJetM->Fill(( tP4.M()));
      }
      plots->fmatchTopJetNum->Fill(NmatchTop);
    
      double truthDRbb = partonList[0]->P4().DeltaR(partonList[3]->P4());
      if(NmatchTop==2){
          plots->fmatchTop2_minDRjj->Fill(minDRjj);
          plots->fmatchTop2_minDRqq->Fill(minDRqq);
          plots->fmatchTop2_truthDRbb->Fill(truthDRbb);
          plots->fmatchTop2_Jet4PT->Fill(jetList[3]->PT);
          plots->fmatchTop2_Jet5PT->Fill(jetList[4]->PT);
          plots->fmatchTop2_Jet6PT->Fill(jetList[5]->PT);
      }
      else if(NmatchTop==1){
          plots->fmatchTop1_minDRjj->Fill(minDRjj);
          plots->fmatchTop1_minDRqq->Fill(minDRqq);
          plots->fmatchTop1_truthDRbb->Fill(truthDRbb);
          plots->fmatchTop1_Jet4PT->Fill(jetList[3]->PT);
          plots->fmatchTop1_Jet5PT->Fill(jetList[4]->PT);
          plots->fmatchTop1_Jet6PT->Fill(jetList[5]->PT);
      }
      else {
          plots->fmatchTop0_minDRjj->Fill(minDRjj);
          plots->fmatchTop0_minDRqq->Fill(minDRqq);
          plots->fmatchTop0_truthDRbb->Fill(truthDRbb);
          plots->fmatchTop0_Jet4PT->Fill(jetList[3]->PT);
          plots->fmatchTop0_Jet5PT->Fill(jetList[4]->PT);
          plots->fmatchTop0_Jet6PT->Fill(jetList[5]->PT);
      }
      if(isDEBUG){
          cout <<"DEBUG NmatchTop " << NmatchTop <<endl;
          cout <<"Parton"<<endl;
          for(int i=0;i<partonList.size();i++){
              cout << i<<" "<< partonList[i]->P4().Pt() <<" "<< partonList[i]->P4().Eta()<< " "<< partonList[i]->P4().Phi() <<" "<< matchJetIndex[i] <<endl;
          }
          cout <<"Jet"<<endl;
          for(int i=0;i<jetList.size();i++){
              cout << i<<" "<< jetList[i]->P4().Pt() <<" "<< jetList[i]->P4().Eta()<< " "<< jetList[i]->P4().Phi() <<" "<< matchPartonIndex[i] << endl;
          }
      }
      
      //Fill output tree
      for(int i=0;i<partonList.size();i++){
          GenParticle* particle = (GenParticle*) partonList[i];
          parton_pt.push_back(particle->PT);
          parton_eta.push_back(particle->Eta);
          parton_phi.push_back(particle->Phi);
          parton_mass.push_back(particle->Mass);
          parton_pdgid.push_back(particle->PID);
          parton_barcode.push_back(barcode[i]);
          parton_jet_index.push_back(matchJetIndex[i]);
      }
      
      for(int i=0;i<jetList.size();i++){
          Jet* particle = (Jet*) jetList[i];
          jet_pt.push_back(particle->PT);
          jet_eta.push_back(particle->Eta);
          jet_phi.push_back(particle->Phi);
          jet_mass.push_back(particle->Mass);
          if(isGenJets)
              jet_btag.push_back(particle->Flavor);
          else
              jet_btag.push_back(particle->BTag);
          jet_barcode.push_back(matchBarcode[i]);
          jet_parton_index.push_back(matchPartonIndex[i]);
      }
      event->Fill();
  }//event loop

    plots->flmatchPT ->Divide(plots->flquarkPT );
    plots->fbmatchPT ->Divide(plots->fbquarkPT );
    plots->flmatchEta->Divide(plots->flquarkEta);
    plots->fbmatchEta->Divide(plots->fbquarkEta);

    cout << CountNJets <<" "<< CountNBJets <<" / "<< allEntries <<endl;

    //output Tree
    fout->cd();
    event->Write();
    fout->Close();
}

//------------------------------------------------------------------------------

void PrintHistograms(ExRootResult *result, TestPlots *plots)
{
  result->Print("png");
}

//------------------------------------------------------------------------------

void Example1(const char *inputFile, const char *outputFile, bool isGenJets=false, bool isDEBUG=false)
{
  gSystem->Load("libDelphes");

  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile);

  cout <<"** input file = " << inputFile <<endl;

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  ExRootResult *result = new ExRootResult();

  TestPlots *plots = new TestPlots;

  BookHistograms(result, plots);

  AnalyseEvents(treeReader, plots, outputFile, isGenJets, isDEBUG);

  PrintHistograms(result, plots);

  TString output= Form("hist_%s",outputFile); 
  result->Write(output.Data());

  cout << "** Exiting..." << endl;

  delete plots;
  delete result;
  delete treeReader;
  delete chain;
}

//------------------------------------------------------------------------------
