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



TParticlePDG* pionPDG= TDatabasePDG::Instance()->GetParticle(211);
extern Float_t mPionPDG;
Float_t mPionPDG = pionPDG->Mass();


iFinalStates FindFinalStatesIndex(TClonesArray* branchTrack) {
    iFinalStates iFinalStatesIndexes;
    Int_t nTracks = branchTrack->GetEntries();
    Int_t foundTag = 0; 
    Int_t foundAll = 0; 
    Int_t iPion1 = 99999;
    Int_t iPion2 = 99999;
    Int_t iPion3 = 99999;
    Int_t iLepton = 99999;

    

    // find pion1 displaced and from the same vertex as proton
    for (Int_t itr2 = 0; itr2 < nTracks; itr2++) {
        Track* track2 = (Track*)branchTrack->At(itr2);
        if (abs(track2->PID) != 211) continue;
        GenParticle* trackobj2 = (GenParticle*)track2->Particle.GetObject();
        if (Length(trackobj2->X, trackobj2->Y, trackobj2->Z) < 1e-8) continue;

        // find pion2 displaced and from the same vertex as proton
        for (Int_t itr3 = 0; itr3 < nTracks; itr3++) {
            if (itr2 == itr3) continue;
            Track* track3 = (Track*)branchTrack->At(itr3);
            if (abs(track3->PID) != 211) continue;
            GenParticle* trackobj3 = (GenParticle*)track3->Particle.GetObject();
            if (not(abs(trackobj2->X - trackobj3->X) < 1e-8 && abs(trackobj2->Y - trackobj3->Y) < 1e-8 && abs(trackobj2->Z - trackobj3->Z) < 1e-8)) continue;

            // find pion3 displaced and from the same vertex as proton
            for (Int_t itr4 = 0; itr4 < nTracks; itr4++) {
                if (itr2 == itr4 || itr3 == itr4) continue;
                Track* track4 = (Track*)branchTrack->At(itr4);
                if (abs(track4->PID) != 211) continue;
                GenParticle* trackobj4 = (GenParticle*)track4->Particle.GetObject();
                if (not(abs(trackobj2->X - trackobj4->X) < 1e-8 && abs(trackobj2->Y - trackobj4->Y) < 1e-8 && abs(trackobj2->Z - trackobj4->Z) < 1e-8)) continue;


                if (abs(track2->Charge + track3->Charge + track4->Charge) == 1) { 
                    iPion1 = itr2;
                    iPion2 = itr3;
                    iPion3 = itr4;
                    foundTag = 1;
                    break;
                }
                if (foundTag == 1) break;
            }
            if (foundTag == 1) break;
        }
        if (foundTag == 1) break;
    }


    if (foundTag == 1) {
        Track* track2 = (Track*)branchTrack->At(iPion1);
        Track* track3 = (Track*)branchTrack->At(iPion2);
        Track* track4 = (Track*)branchTrack->At(iPion3);
        GenParticle* trackobj2 = (GenParticle*)track2->Particle.GetObject();
        Int_t nLep = 0;
        for (Int_t itr = 0; itr < nTracks; itr++) {
            if (nLep > 1) {
                iLepton = 99999;
                foundAll = 0;
                break;
            }
            Track* track = (Track*)branchTrack->At(itr);
            if (not (abs(track->PID) == 11 || abs(track->PID) == 13)) continue;
            nLep += 1;
            iLepton = itr;
            foundAll = 1;
        }
        //cout << nLep << endl;
    }


    if (foundAll == 1) {
        iFinalStatesIndexes.foundAll = 1;
        iFinalStatesIndexes.iPion1 = iPion1;
        iFinalStatesIndexes.iPion2 = iPion2;
        iFinalStatesIndexes.iPion3 = iPion3;
        iFinalStatesIndexes.iLepton = iLepton;
    }
    return iFinalStatesIndexes;
}

