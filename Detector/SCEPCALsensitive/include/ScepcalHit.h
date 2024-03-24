#ifndef ScepcalHit_h
#define ScepcalHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

#include "DD4hep/Objects.h"
#include "DD4hep/Segmentations.h"

// CLHEP
#include "CLHEP/Vector/ThreeVector.h"

namespace scepcal {

class ScepcalHit : public G4VHit {
public:

  ScepcalHit();
  virtual ~ScepcalHit();

  // const ScepcalHit& operator=(const ScepcalHit &right);
  // G4bool operator==(const ScepcalHit &right) const;

  G4int operator==(const ScepcalHit&) const;
  inline void *operator new(size_t);
  inline void operator delete(void* aHit);

  virtual void Draw() {};
  virtual void Print() {};

  CLHEP::Hep3Vector position;
  int cellID;
  double energyDeposit;
};

// typedef G4THitsCollection<ScepcalHit> ScepcalHitsCollection;

extern G4ThreadLocal G4Allocator<ScepcalHit>* ScepcalHitAllocator;

inline void* ScepcalHit::operator new(size_t) {
  if (!ScepcalHitAllocator) ScepcalHitAllocator = new G4Allocator<ScepcalHit>;
  return (void*)ScepcalHitAllocator->MallocSingle();
}

inline void ScepcalHit::operator delete(void*aHit) {
  ScepcalHitAllocator->FreeSingle((ScepcalHit*) aHit);
}

}

#endif
