


//
// ===================================================================
// || This code does:                                               || 
// || 1. Tag the true particles (Lambdab, p, K)                     || 
// ||    For calculating q^2_true for ML regression                 || 
// || 2. || 
// ||  || 
// ||  || 
// ||  || 
// ||  || 
// ===================================================================
//
//
//
//
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
#include "SignalBackgroundClassifier.C"

R__ADD_LIBRARY_PATH($DELPHES)
R__LOAD_LIBRARY(libDelphes)



void allinone(const string type) {
    // formmating the input output files
    string inputFile_st;
    string outputFile_st;
    Int_t foundFiles;
    if (type.at(0) == 's') inputFile_st = "./z_Lambdab_pKa.root"; outputFile_st = "./z_Lambdab_pKa_reco.root";
    if (type.at(0) == 'b') inputFile_st = "./z_Lambdab_pKnunu.root"; outputFile_st = "./z_Lambdab_pKnunu_reco.root";
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

    // num_test = 10;
    
    for (Int_t i_en = 0; i_en < num_test; i_en++) {
        treeReader->ReadEntry(i_en);

        // find the correspdoning particles in the truth
        iFinalStates iFSTrue;                               // final states particles, in truth level
        TLorentzVector LambdabTrue, ProtonTrue, KaonTrue;   // define lorentz vector for the truth level particleM
        Int_t passing = 0;
        passing = Classify(branchParticle, &LambdabTrue, &ProtonTrue, &KaonTrue, &iFSTrue);
        if (passing == 0) continue;
        TLorentzVector qTrue = LambdabTrue - ProtonTrue - KaonTrue;

        // reconstruction, search for tracks in detector level
        iFinalStates iFS = FindFinalStatesIndex(branchTrack);
        if (iFS.foundAll == 0) continue;

        TLorentzVector Lambdab, Proton, Kaon, ALP;
        Track* ProtonTr = (Track*)branchTrack->At(iFS.iProton);
        Track* KaonTr = (Track*)branchTrack->At(iFS.iKaon);
        Proton.SetPtEtaPhiM(ProtonTr->PT, ProtonTr->Eta, ProtonTr->Phi, mProtonPDG);
        Kaon.SetPtEtaPhiM(KaonTr->PT, KaonTr->Eta, KaonTr->Phi, mKPDG);
    
        // Lambda_b decay vertex 
        TVector3 vertex;
        GenParticle* ProtonTrObj= (GenParticle*)ProtonTr->Particle.GetObject();
        vertex.SetXYZ(ProtonTrObj->X, ProtonTrObj->Y, ProtonTrObj->Z);

        // reconstruct Lambda_b
        Lambdab = reconstructB4Momentum(vertex, Proton + Kaon, branchTrack, branchEFlowPhoton, branchEFlowNeutralHadron);

        // reconstruct ALP
        ALP = Lambdab - Proton - Kaon;
        Float_t q2 = pow(ALP.E(), 2) - pow(ALP.P(), 2);


        // storing features
        features->iEvt = i_en;

        features->q2 = q2;
        features->pProton = Proton.P();
        features->etaProton = Proton.Eta();
        features->phiProton = Proton.Phi();
        features->pTProton = Proton.Pt();

        features->pKaon = Kaon.P();
        features->etaKaon = Kaon.Eta();
        features->phiKaon = Kaon.Phi();
        features->pTKaon = Kaon.Pt();

        features->pLambdab = Lambdab.P();
        features->etaLambdab = Lambdab.Eta();
        features->phiLambdab = Lambdab.Phi();
        features->pTLambdab = Lambdab.Pt();

        features->Xvertex = vertex.X();
        features->Yvertex = vertex.Y();
        features->Yvertex = vertex.Z();


        features->q2True =  pow(qTrue.E(), 2) -pow(qTrue.P(), 2);

        tr.Fill();



        
        


    }
    tr.Write();
    fea.Close();
    cout << "Writing to: " << outputFile << "\n\n";

}
