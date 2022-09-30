#include <iostream>
using namespace std;
#include "FinalStatesClass.C"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "classes/DelphesClasses.h"

// classifying singals, and storing the truth level indices
Int_t Classify(TClonesArray *branchParticle, TLorentzVector *LambdabTrue, TLorentzVector *ProtonTrue, TLorentzVector *KaonTrue, iFinalStates *iFSTrue) {
    GenParticle *particle;
    GenParticle *particleB;
    GenParticle *particleC;
    GenParticle *particle0;
    GenParticle *particle0M;
    Int_t nParticles = branchParticle->GetEntries();

    Int_t iBTrue, foundBTrue = 0;
    Int_t iProtonTrue, foundProtonTrue = 0;
    Int_t iKaonTrue, foundKaonTrue = 0;
    for (int ipB = 0; ipB < nParticles; ipB++) {
        foundBTrue = 0;
        particleB = (GenParticle *)branchParticle->At(ipB);
        if (abs(particleB->PID) == 5122) {
            iBTrue = ipB;
            foundBTrue = 1;

            Int_t iProtonTrue_i, foundProtonTrue_i = 0;
            Int_t iKaonTrue_i, foundKaonTrue_i = 0;
            for (int ip2 = 0; ip2 < nParticles; ip2++) {
                particle0 = (GenParticle *)branchParticle->At(ip2);
                if (particle0->M1 == -1) continue;
                particle0M = (GenParticle *)branchParticle->At(particle0->M1);

                if (abs(particle0->PID) == 2212 && particle0->M1 == iBTrue) {
                    iProtonTrue_i = ip2;
                    foundProtonTrue_i = 1;
                }
                if (abs(particle0->PID) == 321 && particle0->M1 == iBTrue) {
                    iKaonTrue_i = ip2;
                    foundKaonTrue_i = 1;
                }
            }
            if (foundProtonTrue_i == 1 && foundKaonTrue_i == 1) {
                iProtonTrue = iProtonTrue_i;
                foundProtonTrue = 1;
                iKaonTrue = iKaonTrue_i;
                foundKaonTrue = 1;
                break;
            }
        }
    }
    if (foundBTrue == 1 && foundProtonTrue == 1 && foundKaonTrue == 1) {
        iFinalStates iFSTrue_;
        iFSTrue_.iLambdab = iBTrue;
        iFSTrue_.iProton = iProtonTrue;
        iFSTrue_.iKaon = iKaonTrue;

        TLorentzVector BTrue_, ProtonTrue_, KaonTrue_;
        GenParticle *particleProton;
        GenParticle *particleKaon;
        particleB = (GenParticle *)branchParticle->At(iBTrue);
        BTrue_.SetPtEtaPhiE(particleB->PT, particleB->Eta, particleB->Phi, particleB->E);
        particleProton = (GenParticle *)branchParticle->At(iProtonTrue);
        ProtonTrue_.SetPtEtaPhiE(particleProton->PT, particleProton->Eta, particleProton->Phi, particleProton->E);
        particleKaon = (GenParticle *)branchParticle->At(iKaonTrue);
        KaonTrue_.SetPtEtaPhiE(particleKaon->PT, particleKaon->Eta, particleKaon->Phi, particleKaon->E);


        *iFSTrue = iFSTrue_;
        *LambdabTrue = BTrue_;
        *ProtonTrue = ProtonTrue_;
        *KaonTrue = KaonTrue_;

        return 1;
    } else {
        return 0;
    }
}
