#ifndef SimG4SCEPCALActions_h
#define SimG4SCEPCALActions_h 1

#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/MCParticleCollection.h"

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/SystemOfUnits.h"

#include "k4FWCore/DataHandle.h"
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/ISimG4ActionTool.h"

#include "SimG4SCEPCALActionInitialization.h"

#include "SCEPCALSegmentation.h"

class SimG4SCEPCALActions : public AlgTool, virtual public ISimG4ActionTool {
public:
  explicit SimG4SCEPCALActions(const std::string& type, const std::string& name, const IInterface* parent);
  virtual ~SimG4SCEPCALActions();

  virtual StatusCode initialize() final;
  virtual StatusCode finalize() final;

  virtual G4VUserActionInitialization* userActionInitialization() final;

private:
  ServiceHandle<IGeoSvc> m_geoSvc;

  dd4hep::DDSegmentation::SCEPCALSegmentation* pSeg;

  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "SCEPCALreadout", "readout name of SCEPCAL"};
  Gaudi::Property<std::string> m_scintName{this, "scintName", "PbWO", "Name of the scintillators"};

  Gaudi::Property<double> m_birks{this, "birks", 0.126, "Birk's constant for the scintillators in mm/MeV"};
  Gaudi::Property<double> m_thres{this, "thres", 0.0, "Energy threshold to store 3d SimCalorimeterHits in GeV"};
};

#endif
