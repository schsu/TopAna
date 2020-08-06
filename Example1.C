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
    100,0.0,300.0);

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
      100,150.0,200.0);

    plots->fmatchWJetM = result->AddHist1D(
      "matchJet_W_M", "",
      "Mass [GeV]", "number",
      100,0.0,150.0);

    plots->fmatchTopJetM = result->AddHist1D(
      "matchJet_Top_M", "",
      "Mass [GeV]", "number",
      100,0.0,220.0);

    plots->fmatchTopJetNum = result->AddHist1D(
    "matchJets", "",
    "matchTopJets", "number",
    4,-0.5,3.5);

}

//------------------------------------------------------------------------------

void AnalyseEvents(ExRootTreeReader *treeReader, TestPlots *plots)
{
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
  TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");

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

  // Loop over all events
  //  allEntries=1000;

  int CountNJets=0, CountNBJets=0;

  for(entry = 0; entry < allEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    // Loop Particle
   GenParticle *w1=0, *b1=0, *q1=0, *q2=0;
   GenParticle *w2=0, *b2=0, *q3=0, *q4=0;
   for(i = 0; i < branchParticle->GetEntriesFast(); ++i)
    {
      particle = (GenParticle*) branchParticle->At(i);

      int pid = fabs(particle->PID);
      if(pid == 6) {
         GenParticle* dau1 = (GenParticle*) branchParticle->At(particle->D1);
         GenParticle* dau2 = (GenParticle*) branchParticle->At(particle->D2);

         while(particle->PID == dau1->PID || particle->PID==dau2->PID){
             if(particle->PID == dau1->PID) particle = dau1;
             else particle = dau2;
             dau1 = (GenParticle*) branchParticle->At(particle->D1);
             dau2 = (GenParticle*) branchParticle->At(particle->D2);
             if( fabs(dau2->PID) ==24 ){
                dau1= (GenParticle*) branchParticle->At(particle->D2);
                dau2= (GenParticle*) branchParticle->At(particle->D1);
             }
         }
         GenParticle *Wtmp= dau1;
         GenParticle *btmp= dau2;
         dau1 = (GenParticle*) branchParticle->At(Wtmp->D1);
         dau2 = (GenParticle*) branchParticle->At(Wtmp->D2);
  
         while(Wtmp->PID == dau1->PID || Wtmp->PID == dau2->PID){
             if(Wtmp->PID == dau1->PID) Wtmp = dau1;
             else Wtmp = dau2; 

             dau1 = (GenParticle*) branchParticle->At(Wtmp->D1);
             dau2 = (GenParticle*) branchParticle->At(Wtmp->D2);
         }
         
         if(!w1){
           w1 = Wtmp;
           b1 = btmp;
           q1 = dau1;
           q2 = dau2;
         }      
         else{
           w2 = Wtmp;
           b2 = btmp;
           q3 = dau1;
           q4 = dau2;
         }
      }
      //Debug
      //cout <<i <<" "<< particle->PID<<" "<< particle->Status <<" "<< particle->M1 <<" "<< particle->D1<<" "<< particle->D2 <<endl;
      //plots->fParticlePT->Fill(particle->PT);
      //plots->fParticleEta->Fill(particle->Eta);
      if(w1 && w2) break;
    }
    //Debug
    //cout << q1->PID <<" "<< q2->PID<<" "<< b1->PID<<" "<<q3->PID<<" " <<q4->PID<<" "<< b2->PID<<endl;
 
    vector<GenParticle*> partonList;
      partonList.push_back(q1);
      partonList.push_back(q2);
      partonList.push_back(b1);
      partonList.push_back(q3);
      partonList.push_back(q4);
      partonList.push_back(b2);

    for(int i=0;i<partonList.size();i++){
      GenParticle* parton = partonList[i];
      cout << parton->PT <<endl;
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
    vector<Jet*> jetList;
    vector<Jet*> bjetList;  
    // Loop over all jets in event
    for(i = 0; i < branchJet->GetEntriesFast(); ++i)
    {
      jet = (Jet*) branchJet->At(i);

      plots->fJetPT->Fill(jet->PT);
      plots->fJetEta->Fill(jet->Eta);

      if(jet->PT<25 || fabs(jet->Eta)>2.5) continue;
      jetList.push_back(jet);
      
      if(jet->BTag)
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
      plots->fNJets->Fill(jetList.size());
      plots->fNBJets->Fill(bjetList.size());
 
      //=============================================
      // Event Selection
      //=============================================
      if(jetList.size()<6) continue;
      CountNJets++;
      if(bjetList.size()<2) continue;
      CountNBJets++;

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
          if(dRList[i].first >0.5) break;
          int jdx = dRList[i].second/pNum;
          int pdx = dRList[i].second%pNum;
          
          cout << jdx <<" "<< pdx <<" "<< dRList[i].first<<endl;
          bool hasMatched=false;
          for(int j=0;j<pNum;j++)
              if(match[jdx][j]) hasMatched=true;
          for(int j=0;j<jNum;j++)
              if(match[j][pdx]) hasMatched=true;

          if(hasMatched==false) match[jdx][pdx]=true;
      }
      

      //Parton to Jet matching index
      vector<int> matchJetIndex (pNum, -1);
      for(int j=0;j< pNum;j++){
          for(int i=0;i<jNum;i++){
              if(match[i][j]){
                  matchJetIndex[j]=i;
                  cout <<"match jet "<<i <<" to parton "<< j <<endl;
              }
          }
      }
      
      
      //======================================
      // Plot matched top candidates
      //
      
      int matchTop=0;
      
      if(matchJetIndex[0]>-1 && matchJetIndex[1]>-1 && matchJetIndex[2]>-1){
          matchTop++;
          Jet *q1_jet = jetList[matchJetIndex[0]];
          Jet *q2_jet = jetList[matchJetIndex[1]];
          Jet *b1_jet = jetList[matchJetIndex[2]];

          TLorentzVector wP4 = q1_jet->P4()+q2_jet->P4();
          TLorentzVector tP4 = b1_jet->P4()+wP4;

          plots->fmatchWJetM->Fill(( wP4.M()));
          plots->fmatchTopJetM->Fill(( tP4.M()));
      }
      if(matchJetIndex[3]>-1 && matchJetIndex[4]>-1 && matchJetIndex[5]>-1){
          matchTop++;
          Jet *q3_jet = jetList[matchJetIndex[3]];
          Jet *q4_jet = jetList[matchJetIndex[4]];
          Jet *b2_jet = jetList[matchJetIndex[5]];

          TLorentzVector wP4 = q3_jet->P4()+q4_jet->P4();
          TLorentzVector tP4 = b2_jet->P4()+wP4;

          plots->fmatchWJetM->Fill(( wP4.M()));
          plots->fmatchTopJetM->Fill(( tP4.M()));
      }
      plots->fmatchTopJetNum->Fill(matchTop);
      
  }//event loop

  cout << CountNJets <<" "<< CountNBJets <<" / "<< allEntries <<endl;

}

//------------------------------------------------------------------------------

void PrintHistograms(ExRootResult *result, TestPlots *plots)
{
  result->Print("png");
}

//------------------------------------------------------------------------------

void Example1(const char *inputFile)
{
  gSystem->Load("libDelphes");

  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile);

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  ExRootResult *result = new ExRootResult();

  TestPlots *plots = new TestPlots;

  BookHistograms(result, plots);

  AnalyseEvents(treeReader, plots);

  PrintHistograms(result, plots);

  result->Write("results.root");

  cout << "** Exiting..." << endl;

  delete plots;
  delete result;
  delete treeReader;
  delete chain;
}

//------------------------------------------------------------------------------
