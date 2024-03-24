#ifndef SimG4SCEPCALEventAction_h
#define SimG4SCEPCALEventAction_h 1

#include "G4UserEventAction.hh"
#include "G4Event.hh"

#include "SimG4SCEPCALSteppingAction.h"

#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"

namespace drc {
class SimG4SCEPCALEventAction : public G4UserEventAction {
    
public:
  SimG4SCEPCALEventAction();
  virtual ~SimG4SCEPCALEventAction();

  virtual void BeginOfEventAction(const G4Event*) final;
  virtual void EndOfEventAction(const G4Event*) final;

  void setSteppingAction(SimG4SCEPCALSteppingAction* steppingAction) { pSteppingAction = steppingAction; }


  edm4hep::SimCalorimeterHitCollection* getEdepsCollection() { return m_Edeps; }
  edm4hep::SimCalorimeterHitCollection* getEdeps3dCollection() { return m_Edeps3d; }
  edm4hep::SimCalorimeterHitCollection* getEdepsCollectionCherenkov() { return m_EdepsCherenkov; }

  edm4hep::MCParticleCollection* getLeakagesCollection() { return m_Leakages; }

private:
  SimG4SCEPCALSteppingAction* pSteppingAction;

  // ownership of collections transferred to DataWrapper<T>
  edm4hep::SimCalorimeterHitCollection* m_Edeps;
  edm4hep::SimCalorimeterHitCollection* m_Edeps3d;
  edm4hep::SimCalorimeterHitCollection* m_EdepsCherenkov;

  edm4hep::MCParticleCollection* m_Leakages;

};

}

#endif
