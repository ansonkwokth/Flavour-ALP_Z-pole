#pragma once
#include <iostream>
using namespace std;
#include <vector>

#include "Geometry.C"
#include "Rtypes.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "classes/DelphesClasses.h"
#include "TDatabasePDG.h"
#include "FinalStatesClass.C"

TParticlePDG* BPDG= TDatabasePDG::Instance()->GetParticle(521);
extern Float_t mBPDG;
Float_t mBPDG = BPDG->Mass();

TParticlePDG* protonPDG= TDatabasePDG::Instance()->GetParticle(2212);
extern Float_t mProtonPDG;
Float_t mProtonPDG = protonPDG->Mass();

TParticlePDG* pionPDG= TDatabasePDG::Instance()->GetParticle(211);
extern Float_t mPionPDG;
Float_t mPionPDG = pionPDG->Mass();

TParticlePDG* K0PDG = TDatabasePDG::Instance()->GetParticle(311);
extern Float_t mK0PDG;
Float_t mK0PDG = K0PDG->Mass();

iFinalStates FindFinalStatesIndex(TClonesArray* branchTrack) {
    iFinalStates iFinalStatesIndexes;
    Int_t nTracks = branchTrack->GetEntries();
    Int_t foundAll = 0; 
    Int_t foundVertex = 0;
    Int_t iProton1 = 99999;
    Int_t iProton2 = 99999;
    Int_t iPion1 = 99999;
    Int_t iPion2 = 99999;
    Int_t iPion3 = 99999;

    Int_t nProtons = 0;
    Int_t nPions = 0;
    for (Int_t itr = 0; itr < nTracks; itr++) {
        Track* track = (Track*)branchTrack->At(itr);
        if (abs(track->PID) == 2212) nProtons += 1;
        if (abs(track->PID) == 211) nPions += 1;
    }
    // cout << nProtons << "; " << nPions << endl;
    // cout << nProtons << endl;

    
    // find proton1 displaced
    for (Int_t itr1 = 0; itr1 < nTracks; itr1++) {
        Track* track1 = (Track*)branchTrack->At(itr1);
        if (abs(track1->PID) != 2212) continue;
        GenParticle* trackobj1 = (GenParticle*)track1->Particle.GetObject();
        if (Length(trackobj1->X, trackobj1->Y, trackobj1->Z) < 1e-8) continue;

        // find pion1 displaced and from the same vertex as proton
        for (Int_t itr2 = 0; itr2 < nTracks; itr2++) {
            Track* track2 = (Track*)branchTrack->At(itr2);
            if (abs(track2->PID) != 211) continue;
            GenParticle* trackobj2 = (GenParticle*)track2->Particle.GetObject();
            if (not(abs(trackobj1->X - trackobj2->X) < 1e-8 && abs(trackobj1->Y - trackobj2->Y) < 1e-8 && abs(trackobj1->Z - trackobj2->Z) < 1e-8)) continue;

            // find pion2 displaced and from the same vertex as proton
            for (Int_t itr3 = 0; itr3 < nTracks; itr3++) {
                if (itr2 == itr3) continue;
                Track* track3 = (Track*)branchTrack->At(itr3);
                if (abs(track3->PID) != 211) continue;
                GenParticle* trackobj3 = (GenParticle*)track3->Particle.GetObject();
                if (not(abs(trackobj1->X - trackobj3->X) < 1e-8 && abs(trackobj1->Y - trackobj3->Y) < 1e-8 && abs(trackobj1->Z - trackobj3->Z) < 1e-8)) continue;

                // find pion3 displaced and from the same vertex as proton
                for (Int_t itr4 = 0; itr4 < nTracks; itr4++) {
                    if (itr2 == itr4 || itr3 == itr4) continue;
                    Track* track4 = (Track*)branchTrack->At(itr4);
                    if (abs(track4->PID) != 211) continue;
                    GenParticle* trackobj4 = (GenParticle*)track4->Particle.GetObject();
                    if (not(abs(trackobj1->X - trackobj4->X) < 1e-8 && abs(trackobj1->Y - trackobj4->Y) < 1e-8 && abs(trackobj1->Z - trackobj4->Z) < 1e-8)) continue;


                    if (abs(track1->Charge + track2->Charge + track3->Charge + track4->Charge) == 2 && 
                        abs(track2->Charge + track3->Charge + track4->Charge) == 1) { 
                        iProton1 = itr1;
                        iPion1 = itr2;
                        iPion2 = itr3;
                        iPion3 = itr4;
                        foundVertex = 1;
                        break;
                    }
                    if (foundVertex == 1) break;
                }
                if (foundVertex == 1) break;
            }
            if (foundVertex == 1) break;
        }
        if (foundVertex == 1) break;
    }

    if (foundVertex == 1) {
        // find proton2 displaced unpaired
        Float_t maxPt = -99999;
        for (Int_t itr1 = 0; itr1 < nTracks; itr1++) {
            if (itr1 == iProton1) continue;
            Track* trackPn = (Track*)branchTrack->At(iProton1);
            Track* track1 = (Track*)branchTrack->At(itr1);
            if (abs(track1->PID) != 2212) continue;
            GenParticle* trackobj1 = (GenParticle*)track1->Particle.GetObject();
            if (Length(trackobj1->X, trackobj1->Y, trackobj1->Z) < 1e-8) continue;
            if (abs(track1->Charge + trackPn->Charge) != 0) continue; 

            TLorentzVector Proton1, Proton2; 
            Proton1.SetPtEtaPhiM(trackPn->PT, trackPn->Eta, trackPn->Phi, mProtonPDG);
            Proton2.SetPtEtaPhiM(track1->PT, track1->Eta, track1->Phi, mProtonPDG);
            if (Proton1.Px()*Proton2.Px() + Proton1.Py()*Proton2.Py() + Proton1.Pz()*Proton2.Pz() <= 0) continue; 
            if (maxPt < Proton2.Pt()) {
                maxPt = Proton2.Pt();
                iProton2 = itr1;
                foundAll = 1;
            }
        }
    }
    
    if (foundAll == 1) {
        iFinalStatesIndexes.foundAll = 1;
        iFinalStatesIndexes.iProton1 = iProton1;
        iFinalStatesIndexes.iProton2 = iProton2;
        iFinalStatesIndexes.iPion1 = iPion1;
        iFinalStatesIndexes.iPion2 = iPion2;
        iFinalStatesIndexes.iPion3 = iPion3;
    }
    return iFinalStatesIndexes;
}


// deduce b-hadron 4 momentum be Z boson 2 body decay
TLorentzVector reconstructB4Momentum(TVector3 v3B, TLorentzVector fromSig, TClonesArray* branchTrack, TClonesArray* branchEFlowPhoton, TClonesArray* branchEFlowNeutralHadron) {
    TLorentzVector signalHemi, recoilHemi, pRem;

    Int_t nTracks = branchTrack->GetEntries();
    Int_t nEFlowPhotons = branchEFlowPhoton->GetEntries();
    Int_t nEFlowNeutralHadrons = branchEFlowNeutralHadron->GetEntries();

    Track* track;
    Tower* eflowph;
    Tower* eflownh;

    // charged particles
    if (nTracks > 0) {
        for (int it = 0; it < nTracks; it++) {
            track = (Track*)branchTrack->At(it);
            TLorentzVector chargedPi;
            chargedPi.SetPtEtaPhiE(track->PT, track->Eta, track->Phi, track->P);
            if (chargedPi.Px() * v3B.X() + chargedPi.Py() * v3B.Y() + chargedPi.Pz() * v3B.Z() > 0) {
                signalHemi += chargedPi;
            } else {
                recoilHemi += chargedPi;
            }
        }
    }

    // photon
    if (nEFlowPhotons > 0) {
        for (int iefp = 0; iefp < nEFlowPhotons; iefp++) {
            eflowph = (Tower*)branchEFlowPhoton->At(iefp);
            TLorentzVector eflowPhoi;
            eflowPhoi.SetPtEtaPhiE(eflowph->ET, eflowph->Eta, eflowph->Phi, eflowph->E);
            if (eflowPhoi.Px() * v3B.X() + eflowPhoi.Py() * v3B.Y() + eflowPhoi.Pz() * v3B.Z() > 0) {
                signalHemi += eflowPhoi;
                pRem += eflowPhoi;
            } else {
                recoilHemi += eflowPhoi;
            }
        }
    }

    // neutral hadrons
    if (nEFlowNeutralHadrons > 0) {
        for (int iefnh = 0; iefnh < nEFlowNeutralHadrons; iefnh++) {
            eflownh = (Tower*)branchEFlowNeutralHadron->At(iefnh);
            Float_t pt = eflownh->ET * pow(eflownh->E * eflownh->E - mK0PDG * mK0PDG, 0.5) / eflownh->E;
            TLorentzVector eflowNeuHi;
            eflowNeuHi.SetPtEtaPhiE(pt, eflownh->Eta, eflownh->Phi, eflownh->E);
            if (eflowNeuHi.Px() * v3B.X() + eflowNeuHi.Py() * v3B.Y() + eflowNeuHi.Pz() * v3B.Z() > 0) {
                signalHemi += eflowNeuHi;
                pRem += eflowNeuHi;
            } else {
                recoilHemi += eflowNeuHi;
            }
        }
    }
    
    pRem = signalHemi - fromSig;

    Float_t EB = (91.0 * 91.0 + signalHemi.M() * signalHemi.M() - recoilHemi.M() * recoilHemi.M()) / (2 * 91.0) - pRem.E();

    Float_t pB = pow(EB * EB - mBPDG* mBPDG, 0.5);
    TLorentzVector B;
    Float_t ratio = pB / (pow(v3B.X() * v3B.X() + v3B.Y() * v3B.Y() + v3B.Z() * v3B.Z(), 0.5));
    B.SetPxPyPzE(v3B.X() * ratio, v3B.Y() * ratio, v3B.Z() * ratio, EB);

    return B;
}

