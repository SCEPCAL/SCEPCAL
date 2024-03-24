#ifndef SimG4SaveSCEPCALMCTruth_h
#define SimG4SaveSCEPCALMCTruth_h 1

// Data model
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"

#include "GaudiAlg/GaudiTool.h"
#include "k4FWCore/DataHandle.h"
#include "k4Interface/ISimG4Svc.h"
#include "k4Interface/ISimG4SaveOutputTool.h"

#include "SimG4SCEPCALEventAction.h"

class IGeoSvc;

class SimG4SaveSCEPCALMCTruth : public GaudiTool, virtual public ISimG4SaveOutputTool {
public:
  explicit SimG4SaveSCEPCALMCTruth(const std::string& aType, const std::string& aName, const IInterface* aParent);
  virtual ~SimG4SaveSCEPCALMCTruth();

  virtual StatusCode initialize();
  virtual StatusCode finalize();

  virtual StatusCode saveOutput(const G4Event& aEvent) final;

private:
  DataHandle<edm4hep::SimCalorimeterHitCollection> m_Edeps{"SimCalorimeterHits", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::SimCalorimeterHitCollection> m_Edeps3d{"Sim3dCalorimeterHits", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::SimCalorimeterHitCollection> m_EdepsCherenkov{"SimCalorimeterHitsCherenkov", Gaudi::DataHandle::Writer, this};

  DataHandle<edm4hep::MCParticleCollection> m_Leakages{"Leakages", Gaudi::DataHandle::Writer, this};

  drc::SimG4SCEPCALEventAction* m_eventAction;

  ServiceHandle<ISimG4Svc> m_geantSvc;
};

#endif
