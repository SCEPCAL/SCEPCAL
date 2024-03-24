#ifndef DETSENSITIVE_SCEPCALSD_H
#define DETSENSITIVE_SCEPCALSD_H

#include "DD4hep/Segmentations.h"

// Geant
#include "G4THitsCollection.hh"
#include "G4VSensitiveDetector.hh"

// #include "G4SystemOfUnits.hh"
// #include "G4PhysicalConstants.hh"
// #include "G4Step.hh"
// #include "G4TouchableHistory.hh"

#include "ScepcalHit.h"
#include "SCEPCALSegmentation.h"

// namespace scepcal {
// class ScepcalHit;

// }

namespace det {

class ScepcalSD : public G4VSensitiveDetector {

public:
    // ScepcalSD();

    ScepcalSD(const std::string& aDetectorName, 
              const std::string& aReadoutName, 
              const dd4hep::Segmentation& aSeg);

    virtual ~ScepcalSD();

    virtual void Initialize(G4HCofThisEvent* aHitsCollections) final;
    virtual bool ProcessHits(G4Step* aStep, G4TouchableHistory*) final;

    // inline int testNumber() {return 9989; }

private:

    G4THitsCollection<scepcal::ScepcalHit>* m_calorimeterCollection;
    dd4hep::Segmentation m_seg;

    dd4hep::DDSegmentation::SCEPCALSegmentation* fSeg;

    std::unordered_map<int, scepcal::ScepcalHit*>* m_hitExists;

    // SiPMHitsCollection* fHitCollection;
};
}

#endif
