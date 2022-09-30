// Define the class of Features, to be stored and used in BDT code

#include "Rtypes.h"
class Features {
public:
    Int_t iEvt;      // event index
    Float_t q2;

    Float_t pProton;
    Float_t etaProton;
    Float_t phiProton;
    Float_t pTProton;

    Float_t pKaon;
    Float_t etaKaon;
    Float_t phiKaon;
    Float_t pTKaon;

    Float_t pLambdab;
    Float_t etaLambdab;
    Float_t phiLambdab;
    Float_t pTLambdab;
    
    Float_t Xvertex;
    Float_t Yvertex;
    Float_t Zvertex;



    Float_t q2True;
    
};
