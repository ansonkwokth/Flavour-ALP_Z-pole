
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
    inputFile_st = "./z_tautau_3pinu_mua.root";
    outputFile_st = "./z_tautau_3pinu_mua_reco.root";
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

    num_test = 10;
    for (Int_t i_en = 0; i_en < num_test; i_en++) {
        treeReader->ReadEntry(i_en);
        // cout << endl;
        // cout << i_en<<endl;

        Int_t iALP; 
        Int_t nParticles = branchParticle->GetEntries();
        Int_t iTau;
        Int_t iLep;
        for (Int_t ipr = 0; ipr < nParticles; ipr++) {
            
            GenParticle* particle1 = (GenParticle*)branchParticle->At(ipr);
            // cout << particle1->PID << endl;
            if (abs(particle1->PID) == 4900022) {
                iALP = ipr; 
            }
            if (particle1->PID == -15) {
                iTau = ipr; 
            }
            if (particle1->PID == -11 || particle1->PID == -13) {
            //if (particle1->PID == 11 || particle1->PID == 13) {
                iLep = ipr; 
            }
        }
        GenParticle* particleALP = (GenParticle*)branchParticle->At(iALP);
        TLorentzVector ALPTrue;
        ALPTrue.SetPtEtaPhiE(particleALP->PT, particleALP->Eta, particleALP->Phi, particleALP->E);
        GenParticle* particleTau = (GenParticle*)branchParticle->At(iTau);
        GenParticle* particleLep = (GenParticle*)branchParticle->At(iLep);


       
        iFinalStates iFS = FindFinalStatesIndex(branchTrack);
        if (iFS.foundAll == 0) continue;

        Track* Pion1Tr = (Track*)branchTrack->At(iFS.iPion1);
        Track* Pion2Tr = (Track*)branchTrack->At(iFS.iPion2);
        Track* Pion3Tr = (Track*)branchTrack->At(iFS.iPion3);
        Track* LepTr = (Track*)branchTrack->At(iFS.iLepton);

        TVector3 vertexTag;
        GenParticle* Pion1TrObj= (GenParticle*)Pion1Tr->Particle.GetObject();
        vertexTag.SetXYZ(Pion1Tr->X, Pion1Tr->Y, Pion1Tr->Z);
        
        TLorentzVector Lep, Pion1, Pion2, Pion3, ALP;
        Pion1.SetPtEtaPhiM(Pion1Tr->PT, Pion1Tr->Eta, Pion1Tr->Phi, mPionPDG);
        Pion2.SetPtEtaPhiM(Pion2Tr->PT, Pion2Tr->Eta, Pion2Tr->Phi, mPionPDG);
        Pion3.SetPtEtaPhiM(Pion3Tr->PT, Pion3Tr->Eta, Pion3Tr->Phi, mPionPDG);
        Lep.SetPtEtaPhiM(LepTr->PT, LepTr->Eta, LepTr->Phi, mPionPDG);

        Float_t LTagLep, sTag, sLep;
        TVector3 LepPoint;
        GenParticle* LepTrObj= (GenParticle*)LepTr->Particle.GetObject();
        LepPoint.SetXYZ(LepTrObj->X, LepTrObj->Y, LepTrObj->Z);


        Float_t LvertexTagLep, sVertexTag, sLeptonPoint;
        // distance_2lines(0, 0, 0, vertexTag.X(), vertexTag.Y(), vertexTag.Z(),
        distance_2lines(0, 0, 0, vertexTag.X(), vertexTag.Y(), vertexTag.Z(),
                        LepPoint.X(), LepPoint.Y(), LepPoint.Z(), particleLep->Px,particleLep->Py, particleLep->Pz,
                        &LvertexTagLep, &sVertexTag, &sLeptonPoint);
        // distance_2lines(1, 1,2, 1, 1, 0,
        //                1, 1,1, 1, 1, 0,
        //                &LvertexTagLep, &sVertexTag, &sLeptonPoint);
        //cout << LvertexTagLep << endl;

        //Float_t X = vertexTag.X() + sVertexTag *vertexTag.X();
        //Float_t Y = vertexTag.Y() + sVertexTag *vertexTag.Y();
        //Float_t Z = vertexTag.Z() + sVertexTag *vertexTag.Z();
        Float_t X = sVertexTag *vertexTag.X();
        Float_t Y = sVertexTag *vertexTag.Y();
        Float_t Z = sVertexTag *vertexTag.Z();
        TVector3 vertex(X, Y, Z);

        // cout << vertex.X() << "; " << LepPoint.X() << endl;
        cout << (particleTau->Px / particleTau->E) <<"; " << LepPoint.X()<< endl;
        cout << particleALP->X<<endl;


        //Float_t q2 = pow(ALP.E(), 2) -pow(ALP.P(), 2);

        // if (q2 < -10) cout << " ................................ " << endl;
        //features->q2 = q2;

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
