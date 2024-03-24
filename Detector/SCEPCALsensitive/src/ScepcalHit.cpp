#include "ScepcalHit.h"

namespace scepcal {

G4ThreadLocal G4Allocator<ScepcalHit>* ScepcalHitAllocator = 0;
ScepcalHit::~ScepcalHit() {}
ScepcalHit::ScepcalHit() {}

G4int ScepcalHit::operator==(const ScepcalHit& right) const { return (this == &right) ? 1 : 0; }

}
