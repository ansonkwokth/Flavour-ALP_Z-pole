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

TParticlePDG* LambdabPDG = TDatabasePDG::Instance()->GetParticle(5122);
extern Float_t mLambdabPDG;
Float_t mLambdabPDG = LambdabPDG->Mass();

TParticlePDG* protonPDG = TDatabasePDG::Instance()->GetParticle(2212);
extern Float_t mProtonPDG;
Float_t mProtonPDG = protonPDG->Mass();

TParticlePDG* KPDG = TDatabasePDG::Instance()->GetParticle(321);
extern Float_t mKPDG;
Float_t mKPDG = KPDG->Mass();

TParticlePDG* K0PDG = TDatabasePDG::Instance()->GetParticle(311);
extern Float_t mK0PDG;
Float_t mK0PDG = K0PDG->Mass();

iFinalStates FindFinalStatesIndex(TClonesArray* branchTrack) {
    iFinalStates iFinalStatesIndexes;
    Int_t nTracks = branchTrack->GetEntries();
    Int_t foundProton = 0;
    Int_t foundKaon = 0;
    Int_t iProton = 99999;
    Int_t iKaon = 99999;
    // find proton displaced
    for (Int_t itr1 = 0; itr1 < nTracks; itr1++) {
        Track* track1 = (Track*)branchTrack->At(itr1);
        if (abs(track1->PID) != 2212) continue;
        GenParticle* trackobj1 = (GenParticle*)track1->Particle.GetObject();
        if (Length(trackobj1->X, trackobj1->Y, trackobj1->Z) < 1e-8) continue;

        // find kaon displaced and from the same vertex as proton
        for (Int_t itr2 = 0; itr2 < nTracks; itr2++) {
            Track* track2 = (Track*)branchTrack->At(itr2);
            if (abs(track2->PID) != 321) continue;
            GenParticle* trackobj2 = (GenParticle*)track2->Particle.GetObject();
            if (abs(trackobj1->X - trackobj2->X) < 1e-8 && 
                abs(trackobj1->Y - trackobj2->Y) < 1e-8 && 
                abs(trackobj1->Z - trackobj2->Z) < 1e-8 &&
                track1->Charge + track2->Charge == 0) {
                iProton = itr1;
                iKaon = itr2;
                foundKaon = 1;
                break;
                // vertex.SetXYZ(trackobj1->X, trackobj1->Y, trackobj1->Z);
            }
            if (foundKaon == 1) break;
        }
        if (foundKaon == 1) break;
    }
    if (foundKaon == 1) {
        iFinalStatesIndexes.iProton = iProton;
        iFinalStatesIndexes.iKaon = iKaon;
        iFinalStatesIndexes.foundAll = 1;
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

    Float_t pB = pow(EB * EB - mLambdabPDG* mLambdabPDG, 0.5);
    TLorentzVector B;
    Float_t ratio = pB / (pow(v3B.X() * v3B.X() + v3B.Y() * v3B.Y() + v3B.Z() * v3B.Z(), 0.5));
    B.SetPxPyPzE(v3B.X() * ratio, v3B.Y() * ratio, v3B.Z() * ratio, EB);

    return B;
}

