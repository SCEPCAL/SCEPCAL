#include "SimG4SCEPCALActionInitialization.h"
#include "SimG4SCEPCALSteppingAction.h"
#include "SimG4SCEPCALEventAction.h"
#include "CLHEP/Units/SystemOfUnits.h"

namespace drc {
SimG4SCEPCALActionInitialization::SimG4SCEPCALActionInitialization(): G4VUserActionInitialization() {}

SimG4SCEPCALActionInitialization::~SimG4SCEPCALActionInitialization() {}

void SimG4SCEPCALActionInitialization::setBirksConstant(const std::string scintName, const double birks) {
  m_scintName = scintName;
  m_birks = birks;
}

void SimG4SCEPCALActionInitialization::Build() const {

  SimG4SCEPCALSteppingAction* steppingAction = new SimG4SCEPCALSteppingAction(); // deleted by G4

  steppingAction->setSegmentation(pSeg);
  steppingAction->setThreshold(m_thres);

  SetUserAction(steppingAction);

  SimG4SCEPCALEventAction* eventAction = new SimG4SCEPCALEventAction(); // deleted by G4

  eventAction->setSteppingAction(steppingAction);
  SetUserAction(eventAction);

  G4Material::GetMaterial(m_scintName)->GetIonisation()->SetBirksConstant(m_birks*CLHEP::millimeter/CLHEP::MeV); // makeshift for DD4hep
}
}
