#include "SimG4SCEPCALSteppingAction.h"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4VProcess.hh"

#include "CLHEP/Units/SystemOfUnits.h"
#include "DD4hep/DD4hepUnits.h"
#include <bitset>

namespace drc {

SimG4SCEPCALSteppingAction::SimG4SCEPCALSteppingAction()
: G4UserSteppingAction(), fPrevTower(0), fPrevTowerCherenkov(0), fDebugLevel(5), fCounter(0) {}


SimG4SCEPCALSteppingAction::~SimG4SCEPCALSteppingAction() {}

void SimG4SCEPCALSteppingAction::UserSteppingAction(const G4Step* step) {

  fCounter++;
  if (fCounter==100000000) {
    std::cout<<"100 million steps calculated"<<std::endl;
    fCounter = 0;
  }

  //Unpack once
  G4Track*              track             = step->GetTrack();
  G4StepPoint*          preStepPoint      = step->GetPreStepPoint();
  G4TouchableHandle     preStepTouchable  = preStepPoint->GetTouchableHandle();
  G4StepPoint*          postStepPoint     = step->GetPostStepPoint();
  G4ParticleDefinition* particle          = track->GetDefinition();
  
  if (preStepTouchable->GetHistoryDepth()<1) {
    //historyDepth    world
    //historyDepth-1  experimentalhall (detector)
    //historyDepth-2  scepcalAssemblyVol, assemblyEnvelopVol

    return; // skip particles in the world or assembly volume
  }

  //if (system==1)
  SteppingAction(step,
                 track,
                 preStepPoint,
                 preStepTouchable,
                 particle);

  // leakage particles
  if (postStepPoint->GetStepStatus() == fWorldBoundary) {
    saveLeakage(track,preStepPoint);
    return;
  }
}

void SimG4SCEPCALSteppingAction::SteppingAction(const G4Step*         step,
                                                G4Track*              track,
                                                G4StepPoint*          preStepPoint,
                                                G4TouchableHandle     preStepTouchable,
                                                G4ParticleDefinition* particle) {

  float  edep       = step->GetTotalEnergyDeposit()*CLHEP::MeV/CLHEP::GeV;
  auto   copyNum64  = pSeg->convertFirst32to64(preStepTouchable->GetCopyNumber());

  if ((track->GetCurrentStepNumber()==1) && particle==G4OpticalPhoton::OpticalPhotonDefinition()) {
    G4String processName=track->GetCreatorProcess()->GetProcessName();

    if (processName == "Cerenkov") {
      accumulateCherenkov(fPrevTowerCherenkov, copyNum64);
      track->SetTrackStatus(fKillTrackAndSecondaries);

/*     
      //kill very long or short wavelengths
      float photWL=1239.84187/(track->GetTotalEnergy()*CLHEP::eV/CLHEP::GeV);

      if (photWL>1000||photWL<300) {
        if (fDebugLevel<2) std::cout << "Cerenkov photon killed with WL: "<<photWL<<std::endl;
        track->SetTrackStatus(fKillTrackAndSecondaries);
      } else {
        accumulateCherenkov(fPrevTowerCherenkov, copyNum64);
        //do not propagate the photon
        track->SetTrackStatus(fKillTrackAndSecondaries);
      }
*/
    
    };

    return;
  }

/*  
    // Save every hit
    if (edep>m_thres) {
    auto simEdep3d=m_Edeps3d->create();
    auto &pos=preStepPoint->GetPosition();

    simEdep3d.setCellID(static_cast<unsigned long long>(copyNum64));
    simEdep3d.setEnergy(edep);
    simEdep3d.setPosition({static_cast<float>(pos.x()*CLHEP::millimeter),
                           static_cast<float>(pos.y()*CLHEP::millimeter),
                           static_cast<float>(pos.z()*CLHEP::millimeter)});
    }
*/
  accumulate(fPrevTower, copyNum64, edep);
}

void SimG4SCEPCALSteppingAction::accumulate(unsigned int &prev,
                                            dd4hep::DDSegmentation::CellID& copyNum64,
                                            float edep) {

  // search for the element
  bool found = false;
  edm4hep::MutableSimCalorimeterHit* thePtr = nullptr;

  if ( m_Edeps->size() > prev ) { // check previous element
    auto element = m_Edeps->at(prev);
    if ( checkId(element, copyNum64) ) {
      thePtr = &element;
      found = true;
    }
  }
  if (!found) { // fall back to loop
    for (unsigned int iElement = 0; iElement<m_Edeps->size(); iElement++) {
      auto element = m_Edeps->at(iElement);
      if ( checkId(element, copyNum64) ) {
        found = true;
        prev = iElement;
        thePtr = &element;
        break;
      }
    }
  }
  if (!found) { // create

    auto simEdep = m_Edeps->create();
    simEdep.setCellID( static_cast<unsigned long long>(copyNum64) );
    simEdep.setEnergy(0.); // added later

    auto pos = pSeg->myPosition(copyNum64);

    simEdep.setPosition( { static_cast<float>(pos.x()*CLHEP::millimeter),
                           static_cast<float>(pos.y()*CLHEP::millimeter),
                           static_cast<float>(pos.z()*CLHEP::millimeter) } );

    prev = m_Edeps->size();
    thePtr = &simEdep;
  }

  auto edepPrev = thePtr->getEnergy();
  thePtr->setEnergy( edepPrev + edep );

}

void SimG4SCEPCALSteppingAction::accumulateCherenkov(unsigned int &prev,
                                                     dd4hep::DDSegmentation::CellID&  copyNum64) {

  // search for the element
  bool found = false;
  edm4hep::MutableSimCalorimeterHit* thePtr = nullptr;

  if ( m_EdepsCherenkov->size() > prev ) { // check previous element
    auto element = m_EdepsCherenkov->at(prev);
    if ( checkId(element, copyNum64) ) {
      thePtr = &element;
      found = true;
    }
  }
  if (!found) { // fall back to loop
    for (unsigned int iElement = 0; iElement<m_EdepsCherenkov->size(); iElement++) {
      auto element = m_EdepsCherenkov->at(iElement);
      if ( checkId(element, copyNum64) ) {
        found = true;
        prev = iElement;
        thePtr = &element;
        break;
      }
    }
  }
  if (!found) { // create

    auto simEdepCherenkov = m_EdepsCherenkov->create();
    simEdepCherenkov.setCellID( static_cast<unsigned long long>(copyNum64) );
    simEdepCherenkov.setEnergy(0.); // added later

    auto pos = pSeg->myPosition(copyNum64);

    simEdepCherenkov.setPosition( { static_cast<float>(pos.x()*CLHEP::millimeter),
                                    static_cast<float>(pos.y()*CLHEP::millimeter),
                                    static_cast<float>(pos.z()*CLHEP::millimeter) } );
    prev = m_EdepsCherenkov->size();
    thePtr = &simEdepCherenkov;
  }

  auto edepPrev = thePtr->getEnergy();
  thePtr->setEnergy( edepPrev + 1 );
}

bool SimG4SCEPCALSteppingAction::checkId(edm4hep::SimCalorimeterHit       edep,
                                         dd4hep::DDSegmentation::CellID&  copyNum64) {
  return (edep.getCellID()==static_cast<unsigned long long>(copyNum64) );
}

void SimG4SCEPCALSteppingAction::saveLeakage(G4Track* track,
                                             G4StepPoint* preStepPoint) {

  auto leakage = m_Leakages->create();
  leakage.setPDG( track->GetDefinition()->GetPDGEncoding() );
  leakage.setGeneratorStatus(1); // leakages naturally belong to final states
  leakage.setCharge( track->GetDefinition()->GetPDGCharge() );
  leakage.setMass( track->GetDefinition()->GetPDGMass()*CLHEP::MeV/CLHEP::GeV );
  leakage.setMomentum( { static_cast<float>(track->GetMomentum().x()*CLHEP::MeV/CLHEP::GeV),
                         static_cast<float>(track->GetMomentum().y()*CLHEP::MeV/CLHEP::GeV),
                         static_cast<float>(track->GetMomentum().z()*CLHEP::MeV/CLHEP::GeV) } );
  leakage.setVertex( { static_cast<float>(preStepPoint->GetPosition().x()*CLHEP::millimeter),
                       static_cast<float>(preStepPoint->GetPosition().y()*CLHEP::millimeter),
                       static_cast<float>(preStepPoint->GetPosition().z()*CLHEP::millimeter) } );
}

} // namespace drc
