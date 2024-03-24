#ifndef SimG4SaveSCEPCALHits_h
#define SimG4SaveSCEPCALHits_h 1

#include "SiPMHit.h"

// Data model
#include "edm4hep/RawCalorimeterHitCollection.h"
#include "edm4hep/SparseVectorCollection.h"

#include "GaudiAlg/GaudiTool.h"
#include "k4FWCore/DataHandle.h"
#include "k4Interface/ISimG4SaveOutputTool.h"
#include "k4Interface/IGeoSvc.h"

class IGeoSvc;

class SimG4SaveSCEPCALHits : public GaudiTool, virtual public ISimG4SaveOutputTool {
public:
  explicit SimG4SaveSCEPCALHits(const std::string& aType, const std::string& aName, const IInterface* aParent);
  virtual ~SimG4SaveSCEPCALHits();

  virtual StatusCode initialize();
  virtual StatusCode finalize();

  virtual StatusCode saveOutput(const G4Event& aEvent) final;

private:
  /// Pointer to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc;

  Gaudi::Property<std::vector<std::string>> m_readoutNames{this, "readoutNames", { "SCEPCALreadout" },
                                                           "Name of the readouts (hits collections) to save"};

  DataHandle<edm4hep::RawCalorimeterHitCollection> mRawCaloHits{"RawCalorimeterHits", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::SparseVectorCollection> mTimeStruct{"RawTimeStructs", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::SparseVectorCollection> mWavlenStruct{"RawWavlenStructs", Gaudi::DataHandle::Writer, this};

};

#endif
