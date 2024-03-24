#include "ScepcalSD.h"
#include "ScepcalHit.h"
#include "SCEPCALSegmentation.h"

// DD4hep
#include "DDG4/Defs.h"
#include "DDG4/Geant4Mapping.h"
#include "DDG4/Geant4VolumeManager.h"
// #include "DD4hep/DD4hepUnits.h"

// Geant4
#include "G4SDManager.hh"

// CLHEP
#include "CLHEP/Vector/ThreeVector.h"

// #include "G4HCofThisEvent.hh"
// #include "G4ParticleDefinition.hh"
// #include "G4ParticleTypes.hh"
// #include "G4SystemOfUnits.hh"

namespace det {

// ScepcalSD::ScepcalSD() {
    // fSeg = dd4hep::DDSegmentation::SCEPCALSegmentation* ;
// }

ScepcalSD::ScepcalSD(const std::string& aDetectorName,
                     const std::string& aReadoutName,
                     const dd4hep::Segmentation& aSeg)
        : G4VSensitiveDetector(aDetectorName), m_calorimeterCollection(nullptr), m_seg(aSeg) {

    collectionName.insert(aReadoutName);

    fSeg = dynamic_cast<dd4hep::DDSegmentation::SCEPCALSegmentation*>( aSeg.segmentation() );
}

ScepcalSD::~ScepcalSD() {
    delete m_hitExists;
}

void ScepcalSD::Initialize(G4HCofThisEvent* aHitsCollections) {
    // create a collection of hits and add it to G4HCofThisEvent
    // deleted in ~G4Event
    m_calorimeterCollection = new G4THitsCollection<scepcal::ScepcalHit>(SensitiveDetectorName, collectionName[0]);
    m_hitExists = new std::unordered_map<int, scepcal::ScepcalHit*>();

    std::cout << "ScepcalSD Initialize" << std::endl;

    aHitsCollections->AddHitsCollection(G4SDManager::GetSDMpointer()->GetCollectionID(m_calorimeterCollection),
                                      m_calorimeterCollection);
}


bool ScepcalSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) {

    // check if energy was deposited
    G4double edep = aStep->GetTotalEnergyDeposit(); //*CLHEP::MeV/CLHEP::GeV;
    if (edep == 0.) return false;

    //Unpack
    // G4Track*              aTrack            = aStep->GetTrack();
    G4StepPoint*          aPreStepPoint     = aStep->GetPreStepPoint();
    G4TouchableHandle     aPreStepTouchable = aPreStepPoint->GetTouchableHandle();
    // G4StepPoint*          aPostStepPoint    = aStep->GetPostStepPoint();
    // G4ParticleDefinition* aParticle         = aTrack->GetDefinition();

    auto copyNum64  = fSeg->convertFirst32to64(aPreStepTouchable->GetCopyNumber());
    int cellID = (int)copyNum64;

    scepcal::ScepcalHit* newHit = nullptr;
    scepcal::ScepcalHit* hitMatch = nullptr;

    std::cout << "ScepcalSD Process hits" << std::endl;

    if (m_hitExists->count(cellID) == 0) { //Add if not found

        std::cout << "new hit" << std::endl;

        newHit = new scepcal::ScepcalHit();

        auto pos = fSeg->myPosition(copyNum64);
        CLHEP::Hep3Vector hitpos(pos.x()*CLHEP::millimeter,
                                 pos.y()*CLHEP::millimeter,
                                 pos.z()*CLHEP::millimeter
                                );
  
        newHit->position = hitpos;
        newHit->cellID = cellID;
        newHit->energyDeposit = 0;

        m_hitExists->emplace(cellID, newHit);
        m_calorimeterCollection->insert(newHit);
    }

    hitMatch = m_hitExists->at(cellID);
    hitMatch->energyDeposit += edep;
  
    return true;
}
}