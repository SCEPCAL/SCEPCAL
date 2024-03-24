#include "SimG4SCEPCALActions.h"

DECLARE_COMPONENT(SimG4SCEPCALActions)

SimG4SCEPCALActions::SimG4SCEPCALActions(const std::string& type, const std::string& name, const IInterface* parent)
: AlgTool(type, name, parent), m_geoSvc("GeoSvc", name) {
  declareInterface<ISimG4ActionTool>(this);
}

SimG4SCEPCALActions::~SimG4SCEPCALActions() {}

StatusCode SimG4SCEPCALActions::initialize() {
  if (AlgTool::initialize().isFailure())
    return StatusCode::FAILURE;

  if (!m_geoSvc) {
    error() << "Unable to locate Geometry Service. "
            << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
    return StatusCode::FAILURE;
  }

  pSeg = dynamic_cast<dd4hep::DDSegmentation::SCEPCALSegmentation*>(m_geoSvc->getDetector()->readout(m_readoutName).segmentation().segmentation());

  return StatusCode::SUCCESS;
}

StatusCode SimG4SCEPCALActions::finalize() { return AlgTool::finalize(); }

G4VUserActionInitialization* SimG4SCEPCALActions::userActionInitialization() {
  auto* actions = new drc::SimG4SCEPCALActionInitialization();

  actions->setSegmentation(pSeg);
  actions->setBirksConstant(m_scintName,m_birks);
  actions->setThreshold(m_thres);

  return actions;
}
