//==========================================================================
//  AIDA Detector description implementation 
//--------------------------------------------------------------------------
// Copyright (C) Organisation europeenne pour la Recherche nucleaire (CERN)
// All rights reserved.
//
// For the licensing terms see $DD4hepINSTALL/LICENSE.
// For the list of contributors see $DD4hepINSTALL/doc/CREDITS.
//
// Author     : M.Frank
//
//==========================================================================

// Framework include files
#include "ScepcalHit.h"
#include "ScepcalSD.h"
#include "SCEPCALSegmentation.h"

// CLHEP
#include "CLHEP/Vector/ThreeVector.h"

#include <DDG4/Geant4SensDetAction.inl>
#include <DDG4/Geant4ParticleInformation.h>
#include <DDG4/Factories.h>


namespace det {

  class ScepcalSDContainer {
  public:
    // typedef scepcal::ScepcalHit Hit;
    // If we need special data to personalize the action, be put it here
    // int mumDeposits = 0;
    // double integratedDeposit = 0;
  };
}

/// Namespace for the AIDA detector description toolkit
namespace dd4hep {

  /// Namespace for the Geant4 based simulation part of the AIDA detector description toolkit
  namespace sim   {

    using namespace det;
    
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //               Geant4SensitiveAction<ScepcalSD>
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    /** \addtogroup Geant4SDActionPlugin
     *
     * @{
     * \package ScepcalSDAction
     * \brief Sensitive detector Scepcal, will produce one hit per step
     *
     * @}
     */

    /// Define collections created by this sensitivie action object
    template <> void Geant4SensitiveAction<ScepcalSDContainer>::defineCollections()    {
      m_collectionID = declareReadoutFilteredCollection<Geant4Calorimeter::Hit>();
    }

    /// Method for generating hit(s) using the information of G4Step object.
    template <> bool 
    Geant4SensitiveAction<ScepcalSDContainer>::process(const G4Step* step,G4TouchableHistory* /*hist*/ ) {

        // G4double edep = step->GetTotalEnergyDeposit(); //*CLHEP::MeV/CLHEP::GeV;
        // if (edep == 0.) return false;
        
        dd4hep::Segmentation* _geoSeg = &m_segmentation;
        auto segmentation=dynamic_cast<dd4hep::DDSegmentation::SCEPCALSegmentation *>(_geoSeg->segmentation());

        typedef Geant4Calorimeter::Hit Hit;
        
        Geant4StepHandler    h(step);
        Geant4TouchableHandler handler(step);

        Position  prePos    = h.prePos();
        Position  postPos   = h.postPos();

        HitContribution      contrib = Hit::extractContribution(step);
        Geant4HitCollection* coll    = collection(m_collectionID);

        G4StepPoint*          aPreStepPoint     = step->GetPreStepPoint();
        G4TouchableHandle     aPreStepTouchable = aPreStepPoint->GetTouchableHandle();

        auto copyNum64  = segmentation->convertFirst32to64(aPreStepTouchable->GetCopyNumber(0));
        int cellID = (int)copyNum64;

        DDSegmentation::Vector3D pos = segmentation->myPosition(copyNum64);

        //   DDSegmentation::Vector3D pos = m_segmentation.position(cell);
      
        // printM2("Segmentation name: %s", m_segmentation.name()); //SCEPCALreadout
        // printM2("SensitiveDetector name: %s", sensitiveDetector().name()); //Scepcal

        // auto scepcalsd=dynamic_cast<det::ScepcalSD *>(sensitiveDetector());

        // printM2("%s> Segmentation params: %f %f %f %f %f %f",
                // c_name(),segmentation->getFdz(),segmentation->getRdz(),segmentation->getnomfw(),
                // segmentation->getnomth(),segmentation->getEBz(),segmentation->getRin());


        Position global(pos.x(),pos.y(),pos.z());

        // Hit* hit = nullptr;
        Hit* hit = coll->findByKey<Hit>(cellID);
        if ( !hit ) {

            hit = new Hit(global);
            hit->cellID = cellID;
            coll->add(cellID, hit);

            printM2("> CREATE hit with deposit:%e MeV  prePos:%8.2f %8.2f %8.2f  %s %d %d  ",
                    contrib.deposit,
                    global.x(),global.y(),global.z(),
                    handler.path().c_str(),
                    aPreStepTouchable->GetHistoryDepth(),
                    aPreStepTouchable->GetCopyNumber(0)
                    );

        }

        hit->truth.emplace_back(contrib);
        hit->energyDeposit += contrib.deposit;
        
            printM2("> UPDATE hit with deposit:%e MeV  prePos:%8.2f %8.2f %8.2f  %s %d %d  ",
                    contrib.deposit,
                    global.x(),global.y(),global.z(),
                    handler.path().c_str(),
                    aPreStepTouchable->GetHistoryDepth(),
                    aPreStepTouchable->GetCopyNumber(0)
                    );
        
        mark(h.track);
        return true;

    }

  }
}

//--- Factory declaration
namespace dd4hep { namespace sim {
    typedef Geant4SensitiveAction<ScepcalSDContainer> ScepcalSDAction;
  }}
DECLARE_GEANT4SENSITIVE(ScepcalSDAction)
