#include "SimG4SCEPCALEventAction.h"

#include "G4RunManager.hh"

namespace drc {
SimG4SCEPCALEventAction::SimG4SCEPCALEventAction(): G4UserEventAction() {}

SimG4SCEPCALEventAction::~SimG4SCEPCALEventAction() {}

void SimG4SCEPCALEventAction::BeginOfEventAction(const G4Event*) {

  m_Edeps = new edm4hep::SimCalorimeterHitCollection();
  m_Edeps3d = new edm4hep::SimCalorimeterHitCollection();
  m_EdepsCherenkov = new edm4hep::SimCalorimeterHitCollection();

  m_Leakages = new edm4hep::MCParticleCollection();

  pSteppingAction->setEdepsCollection(m_Edeps);
  pSteppingAction->setEdeps3dCollection(m_Edeps3d);
  pSteppingAction->setEdepsCollectionCherenkov(m_EdepsCherenkov);

  pSteppingAction->setLeakagesCollection(m_Leakages);

  return;
}

void SimG4SCEPCALEventAction::EndOfEventAction(const G4Event*) {}
} // namespace drc
