#ifndef SimG4SCEPCALSteppingAction_h
#define SimG4SCEPCALSteppingAction_h 1

#include "SCEPCALSegmentation.h"

#include "G4UserSteppingAction.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"

// Data model
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"

namespace drc {
class SimG4SCEPCALSteppingAction : public G4UserSteppingAction {
public:
  SimG4SCEPCALSteppingAction();
  virtual ~SimG4SCEPCALSteppingAction();

  virtual void UserSteppingAction(const G4Step*);

  void setLeakagesCollection(edm4hep::MCParticleCollection* data) { m_Leakages = data; }

  void setSegmentation(dd4hep::DDSegmentation::SCEPCALSegmentation* seg) { pSeg = seg; }

  void setEdepsCollection(edm4hep::SimCalorimeterHitCollection* data) { m_Edeps = data; }
  void setEdeps3dCollection(edm4hep::SimCalorimeterHitCollection* data) { m_Edeps3d = data; }
  void setEdepsCollectionCherenkov(edm4hep::SimCalorimeterHitCollection* data) { m_EdepsCherenkov = data; }

  void setThreshold(const double thres) { m_thres = thres; }

private:
  void SteppingAction(const G4Step*         step,
                      G4Track*              track,
                      G4StepPoint*          preStepPoint,
                      G4TouchableHandle     preStepTouchable,
                      G4ParticleDefinition* particle);

  void accumulate(unsigned int                    &prev,
                  dd4hep::DDSegmentation::CellID& copyNum64,
                  float                           edep);
  void accumulateCherenkov(unsigned int                     &prev,
                           dd4hep::DDSegmentation::CellID&  copyNum64);

  bool checkId(edm4hep::SimCalorimeterHit edep,
               dd4hep::DDSegmentation::CellID& copyNum64);

  void saveLeakage(G4Track* track,
                   G4StepPoint* preStepPoint);

  unsigned int fPrevTower;
  unsigned int fPrevTowerCherenkov;

  dd4hep::DDSegmentation::SCEPCALSegmentation* pSeg;

  // collections owned by SimG4SCEPCALEventAction
  edm4hep::MCParticleCollection* m_Leakages;

  edm4hep::SimCalorimeterHitCollection* m_Edeps;
  edm4hep::SimCalorimeterHitCollection* m_Edeps3d;
  edm4hep::SimCalorimeterHitCollection* m_EdepsCherenkov;

  double m_thres;

  int fDebugLevel;

  int fCounter;

};
}

#endif
