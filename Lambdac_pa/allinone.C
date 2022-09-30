
#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <vector>
using namespace std;

#include "Rtypes.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/TrackCovariance/VertexFit.h"
#include "Geometry.C"
#include "Reconstructor.C"
#include "FeatureClass.C"

R__ADD_LIBRARY_PATH($DELPHES)
R__LOAD_LIBRARY(libDelphes)



void allinone() {
    // formmating the input output files
    string inputFile_st;
    string outputFile_st;
    Int_t foundFiles;
    // inputFile_st = "./z_tautau_3pinu.root";
    inputFile_st = "./z_B_Lambdacp3pi_pa.root";
    outputFile_st = "./z_B_Lambdacp3pi_pa_reco.root";
    const char* inputFile = inputFile_st.c_str();
    const char* outputFile = outputFile_st.c_str();

    cout << "\nReading: " << inputFile << "\n";
    if (gSystem->AccessPathName(inputFile)) {
        cout << "inputFile: " << inputFile << " does not exist" << endl;
        return;
    }
    // Load lib, and read data
    gSystem->Load("libDelphes");
    TChain chain("Delphes");
    chain.Add(inputFile);
    ExRootTreeReader* treeReader = new ExRootTreeReader(&chain);
    Int_t numberOfEntries = treeReader->GetEntries();
    // clone the branches
    TClonesArray* branchParticle = treeReader->UseBranch("Particle");
    TClonesArray* branchTrack = treeReader->UseBranch("Track");
    TClonesArray* branchElectron = treeReader->UseBranch("Electron");
    TClonesArray* branchMuon = treeReader->UseBranch("Muon");
    TClonesArray* branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
    TClonesArray* branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");
    // book feature storing tree and file
    TFile fea(outputFile, "recreate");
    TTree tr("t", "features");
    Features* features = new Features;
    tr.Branch("features", &features);

    Int_t num_test = treeReader->GetEntries();
    cout << num_test << endl;

    
    for (Int_t i_en = 0; i_en < num_test; i_en++) {
        treeReader->ReadEntry(i_en);
        // cout << endl;
        // cout << i_en<<endl;

        Int_t iALP; 
        Int_t nParticles = branchParticle->GetEntries();
        for (Int_t ipr = 0; ipr < nParticles; ipr++) {
            GenParticle* particle1 = (GenParticle*)branchParticle->At(ipr);
            if (abs(particle1->PID) == 4900022) {
                iALP = ipr; 
                break;
            }
        }
        GenParticle* particleALP = (GenParticle*)branchParticle->At(iALP);
        // cout << particleALP->E << endl;
        TLorentzVector ALPTrue;
        ALPTrue.SetPtEtaPhiE(particleALP->PT, particleALP->Eta, particleALP->Phi, particleALP->E);


        
        iFinalStates iFS = FindFinalStatesIndex(branchTrack);
        if (iFS.foundAll == 0) continue;

        Track* Proton1Tr = (Track*)branchTrack->At(iFS.iProton1);
        Track* Proton2Tr = (Track*)branchTrack->At(iFS.iProton2);
        Track* Pion1Tr = (Track*)branchTrack->At(iFS.iPion1);
        Track* Pion2Tr = (Track*)branchTrack->At(iFS.iPion2);
        Track* Pion3Tr = (Track*)branchTrack->At(iFS.iPion3);

        TVector3 vertex;
        GenParticle* Proton1TrObj= (GenParticle*)Proton1Tr->Particle.GetObject();
        vertex.SetXYZ(Proton1TrObj->X, Proton1TrObj->Y, Proton1TrObj->Z);
        
        TLorentzVector Lambdab, Proton1, Proton2, Pion1, Pion2, Pion3, ALP;
        Proton1.SetPtEtaPhiM(Proton1Tr->PT, Proton1Tr->Eta, Proton1Tr->Phi, mProtonPDG);
        Proton2.SetPtEtaPhiM(Proton2Tr->PT, Proton2Tr->Eta, Proton2Tr->Phi, mProtonPDG);
        Pion1.SetPtEtaPhiM(Pion1Tr->PT, Pion1Tr->Eta, Pion1Tr->Phi, mPionPDG);
        Pion2.SetPtEtaPhiM(Pion2Tr->PT, Pion2Tr->Eta, Pion2Tr->Phi, mPionPDG);
        Pion3.SetPtEtaPhiM(Pion3Tr->PT, Pion3Tr->Eta, Pion3Tr->Phi, mPionPDG);

        Lambdab = reconstructB4Momentum(vertex, Proton1 + Proton2 + Pion1 + Pion2 + Pion3, branchTrack, branchEFlowPhoton, branchEFlowNeutralHadron);

        ALP = Lambdab - Proton1 - Proton2 - Pion1 - Pion2 - Pion3;
        Float_t q2 = pow(ALP.E(), 2) -pow(ALP.P(), 2);

        // if (q2 < -10) cout << " ................................ " << endl;
        features->q2 = q2;

        features->pProton1 = Proton1.P();
        features->etaProton1 = Proton1.Eta();
        features->phiProton1 = Proton1.Phi();
        features->pTProton1 = Proton1.Pt();
        features->pProton2 = Proton2.P();
        features->etaProton2 = Proton2.Eta();
        features->phiProton2 = Proton2.Phi();
        features->pTProton2 = Proton2.Pt();

        features->pPion1 = Pion1.P();
        features->etaPion1 = Pion1.Eta();
        features->phiPion1 = Pion1.Phi();
        features->pTPion1 = Pion1.Pt();
        features->pPion2 = Pion2.P();
        features->etaPion2 = Pion2.Eta();
        features->phiPion2 = Pion2.Phi();
        features->pTPion2 = Pion2.Pt();
        features->pPion3 = Pion3.P();
        features->etaPion3 = Pion3.Eta();
        features->phiPion3 = Pion3.Phi();
        features->pTPion3 = Pion3.Pt();

        features->q2True = pow(ALPTrue.M(), 2);
        tr.Fill();



        
        


    }
    tr.Write();
    fea.Close();
    cout << "Writing to: " << outputFile << "\n\n";

}
