#include "SimG4SaveSCEPCALMCTruth.h"

#include "G4RunManager.hh"

DECLARE_COMPONENT(SimG4SaveSCEPCALMCTruth)

SimG4SaveSCEPCALMCTruth::SimG4SaveSCEPCALMCTruth(const std::string& aType, const std::string& aName, const IInterface* aParent)
: GaudiTool(aType, aName, aParent), m_geantSvc("SimG4Svc", aName) {
  declareInterface<ISimG4SaveOutputTool>(this);
}

SimG4SaveSCEPCALMCTruth::~SimG4SaveSCEPCALMCTruth() {}

StatusCode SimG4SaveSCEPCALMCTruth::initialize() {
  if (GaudiTool::initialize().isFailure())
    return StatusCode::FAILURE;

  if (!m_geantSvc) {
    error() << "Unable to locate Geant Simulation Service" << endmsg;
    return StatusCode::FAILURE;
  }

  auto* runManager = G4RunManager::GetRunManager();
  const auto* eventAction = dynamic_cast<const drc::SimG4SCEPCALEventAction*>(runManager->GetUserEventAction());

  if (!eventAction) {
    error() << "Unable to cast to SimG4SCEPCALEventAction from G4UserEventAction" << endmsg;
    return StatusCode::FAILURE;
  }

  m_eventAction = const_cast<drc::SimG4SCEPCALEventAction*>(eventAction); // HACK!!!

  return StatusCode::SUCCESS;
}

StatusCode SimG4SaveSCEPCALMCTruth::finalize() { return GaudiTool::finalize(); }

StatusCode SimG4SaveSCEPCALMCTruth::saveOutput(const G4Event&) {

  auto* edeps = m_eventAction->getEdepsCollection();
  auto* edeps3d = m_eventAction->getEdeps3dCollection();
  auto* edepsCherenkov = m_eventAction->getEdepsCollectionCherenkov();

  auto* leakages = m_eventAction->getLeakagesCollection();

  m_Edeps.put(edeps);
  m_Edeps3d.put(edeps3d);
  m_EdepsCherenkov.put(edepsCherenkov);
  m_Leakages.put(leakages);

  return StatusCode::SUCCESS;
}
