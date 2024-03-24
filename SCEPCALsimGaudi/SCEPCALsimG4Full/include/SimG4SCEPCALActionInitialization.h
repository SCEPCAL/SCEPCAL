#ifndef SimG4SCEPCALActionInitialization_h
#define SimG4SCEPCALActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

#include "SCEPCALSegmentation.h"

namespace drc {
class SimG4SCEPCALActionInitialization : public G4VUserActionInitialization {
public:
  SimG4SCEPCALActionInitialization();
  virtual ~SimG4SCEPCALActionInitialization();

  virtual void Build() const final;

  void setSegmentation(dd4hep::DDSegmentation::SCEPCALSegmentation* seg) { pSeg = seg; }
  void setThreshold(const double thres) { m_thres = thres; }
  void setBirksConstant(const std::string scintName, const double birks);

private:
  dd4hep::DDSegmentation::SCEPCALSegmentation* pSeg;

  std::string m_scintName;
  double m_birks;
  double m_thres;
};
}

#endif
