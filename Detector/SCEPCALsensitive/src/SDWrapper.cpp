#include "DD4hep/Detector.h"
#include "DDG4/Factories.h"

#include "ScepcalSD.h"

namespace dd4hep {
namespace sim {


static G4VSensitiveDetector* create_scepcal_sd(const std::string& aDetectorName, dd4hep::Detector& aLcdd) {
  std::string readoutName = aLcdd.sensitiveDetector(aDetectorName).readout().name();
  return new det::ScepcalSD(
    aDetectorName,readoutName,aLcdd.sensitiveDetector(aDetectorName).readout().segmentation());

  std::cout << "Creating ScepcalSD Factory Method" << std::endl;

}

}
}

DECLARE_EXTERNAL_GEANT4SENSITIVEDETECTOR(ScepcalSD, dd4hep::sim::create_scepcal_sd)