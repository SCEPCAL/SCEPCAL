//=========================================================================
// Author: Wonyong Chung
//=========================================================================

#define VERBOSE_LEVEL 0

#include "SCEPCALSegmentationHandle.h"

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DetectorTools.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Detector.h"
#include "DDRec/DetectorData.h"

#include "TGeoTrd2.h"
#include <bitset>

using dd4hep::Transform3D;
using dd4hep::RotationZYX;
using dd4hep::RotationY;
using ROOT::Math::RotationZ;
using dd4hep::Rotation3D;
using dd4hep::Position;
using dd4hep::rec::LayeredCalorimeterData;

static dd4hep::Ref_t
create_detector_TL(dd4hep::Detector &theDetector, xml_h xmlElement, dd4hep::SensitiveDetector sens) {

  // Import xml objects from compact xml
  xml_det_t detectorXML                   = xmlElement;
  xml_comp_t dimXML                       = detectorXML.child(_Unicode(dim));
  xml_comp_t barrelXML                    = detectorXML.child(_Unicode(barrel));
  xml_comp_t endcapXML                    = detectorXML.child(_Unicode(endcap));
  xml_comp_t crystalFXML                  = detectorXML.child(_Unicode(crystalF));
  xml_comp_t crystalRXML                  = detectorXML.child(_Unicode(crystalR));
  xml_comp_t timingTrXML                  = detectorXML.child(_Unicode(timingLayerTr));
  xml_comp_t timingLgXML                  = detectorXML.child(_Unicode(timingLayerLg));
  xml_comp_t instXML                      = detectorXML.child(_Unicode(inst));
  xml_comp_t scepcalAssemblyGlobalVisXML  = detectorXML.child(_Unicode(scepcalAssemblyGlobalVis));
  xml_comp_t scepcalAssemblyXML           = detectorXML.child(_Unicode(scepcalAssembly));
  xml_comp_t printDebugXML                = detectorXML.child(_Unicode(printDebug));

  // Material definitions
  dd4hep::Material crystalFMat = theDetector.material(crystalFXML.materialStr());
  dd4hep::Material crystalRMat = theDetector.material(crystalRXML.materialStr());
  dd4hep::Material timingTrMat = theDetector.material(timingTrXML.materialStr());
  dd4hep::Material timingLgMat = theDetector.material(timingLgXML.materialStr());
  dd4hep::Material instMat     = theDetector.material(instXML.materialStr());

  // Parse input parameters from imported xml objects
  const int  debugLevel = printDebugXML.attr<int>(_Unicode(level));

  const double Fdz      = crystalFXML.attr<double>(_Unicode(length));
  const double Rdz      = crystalRXML.attr<double>(_Unicode(length));
  const double nomfw    = dimXML.attr<double>(_Unicode(crystalFaceWidthNominal));
  const double nomth    = dimXML.attr<double>(_Unicode(crystalTimingThicknessNominal));
  const double EBz      = dimXML.attr<double>(_Unicode(barrelHalfZ));
  const double Rin      = dimXML.attr<double>(_Unicode(barrelInnerR));

  const int     PHI_SEGMENTS       = dimXML.attr<int>(_Unicode(phiSegments));
  const bool    REFLECT_ENDCAP     = endcapXML.attr<bool>(_Unicode(reflect));

  const bool    CONSTRUCT_BARREL   = barrelXML.attr<bool>(_Unicode(construct));
  const int     BARREL_PHI_START   = barrelXML.attr<int>(_Unicode(phistart));
  const int     BARREL_PHI_END     = barrelXML.attr<int>(_Unicode(phiend));

  const bool    CONSTRUCT_ENDCAP   = endcapXML.attr<bool>(_Unicode(construct));
  const int     ENDCAP_PHI_START   = endcapXML.attr<int>(_Unicode(phistart));
  const int     ENDCAP_PHI_END     = endcapXML.attr<int>(_Unicode(phiend));
  const int     ENDCAP_THETA_START = endcapXML.attr<int>(_Unicode(thetastart));

  //-----------------------------------------------------------------------------------
  // Global geometry numbers
  //-----------------------------------------------------------------------------------

  const double  D_PHI_GLOBAL    = 2*M_PI/PHI_SEGMENTS;

  // Need odd number of N_THETA_BARREL to make center slice
  double  THETA_SIZE_ENDCAP     = atan(Rin/EBz);

  int     N_THETA_BARREL        = int(floor(2*EBz/nomfw))%2==1? floor(2*EBz/nomfw) : floor(2*EBz/nomfw)-1 ;
  int     N_THETA_ENDCAP        = floor(Rin/nomfw);
  
  double  D_THETA_BARREL        = (M_PI-2*THETA_SIZE_ENDCAP)/(N_THETA_BARREL);
  double  D_THETA_ENDCAP        = THETA_SIZE_ENDCAP/N_THETA_ENDCAP;

  int     N_PHI_BARREL_CRYSTAL  = floor(2*M_PI*Rin/(PHI_SEGMENTS*nomfw));
  double  D_PHI_BARREL_CRYSTAL  = D_PHI_GLOBAL/N_PHI_BARREL_CRYSTAL;


  std::cout << std::endl;
  std::cout << "==TIMING LAYER==" << std::endl; 
  std::cout << "GEOMETRY INPUTS:" << std::endl;
  std::cout << "CONSTRUCT_BARREL:   " << CONSTRUCT_BARREL   << std::endl;
  std::cout << "BARREL_PHI_START:   " << BARREL_PHI_START   << std::endl;
  std::cout << "BARREL_PHI_END  :   " << BARREL_PHI_END     << std::endl;
  std::cout << "CONSTRUCT_ENDCAP:   " << CONSTRUCT_ENDCAP   << std::endl;
  std::cout << "ENDCAP_PHI_START:   " << ENDCAP_PHI_START   << std::endl;
  std::cout << "ENDCAP_PHI_END  :   " << ENDCAP_PHI_END     << std::endl;
  std::cout << "ENDCAP_THETA_START: " << ENDCAP_THETA_START << std::endl;
  std::cout << "REFLECT_ENDCAP:     " << REFLECT_ENDCAP     << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "GLOBAL GEOMETRY PARAMETERS:" << std::endl;
  std::cout << "PHI_SEGMENTS:         " << PHI_SEGMENTS         << std::endl;
  std::cout << "D_PHI_GLOBAL:         " << D_PHI_GLOBAL         << std::endl;
  std::cout << "N_THETA_BARREL:       " << N_THETA_BARREL       << std::endl;
  std::cout << "N_THETA_ENDCAP:       " << N_THETA_ENDCAP       << std::endl;
  std::cout << "D_THETA_BARREL:       " << D_THETA_BARREL       << std::endl;
  std::cout << "D_THETA_ENDCAP:       " << D_THETA_ENDCAP       << std::endl;
  std::cout << "N_PHI_BARREL_CRYSTAL: " << N_PHI_BARREL_CRYSTAL << std::endl;
  std::cout << "D_PHI_BARREL_CRYSTAL: " << D_PHI_BARREL_CRYSTAL << std::endl;
  std::cout << "THETA_SIZE_ENDCAP:    " << THETA_SIZE_ENDCAP    << std::endl;
  std::cout << std::endl;

  //-----------------------------------------------------------------------------------
  // Initialize detector element
  //-----------------------------------------------------------------------------------

  std::string name = detectorXML.nameStr();
  dd4hep::DetElement Scepcal(name, detectorXML.id());
  dd4hep::Volume experimentalHall = theDetector.pickMotherVolume(Scepcal);

  dd4hep::xml::Dimension sdType = detectorXML.child(_Unicode(sensitive));
  sens.setType(sdType.typeStr());

  std::cout << "Sensitive Detector Type: " << sdType.typeStr() << std::endl;

  dd4hep::Box     scepcalAssemblyShape(Rin+2*Fdz+2*Rdz, Rin+2*Fdz+2*Rdz, 2*EBz+2*Fdz+2*Rdz);
  dd4hep::Volume  scepcalAssemblyVol("scepcalAssemblyVol", scepcalAssemblyShape, theDetector.material("Vacuum"));
  scepcalAssemblyVol.setVisAttributes(theDetector, scepcalAssemblyGlobalVisXML.visStr());
  // scepcalAssemblyVol.setSensitiveDetector(sens);

  // Initialize the segmentation
  dd4hep::Readout readout = sens.readout();
  dd4hep::Segmentation geomseg = readout.segmentation();
  dd4hep::Segmentation* _geoSeg = &geomseg;
  auto segmentation=dynamic_cast<dd4hep::DDSegmentation::SCEPCALSegmentation *>(_geoSeg->segmentation());
  segmentation->setGeomParams(Fdz, Rdz, nomfw, nomth, EBz, Rin, PHI_SEGMENTS);

  //-----------------------------------------------------------------------------------
  // Barrel + Timing
  //-----------------------------------------------------------------------------------

  // Barrel envelope
  double thC_end          = THETA_SIZE_ENDCAP+D_THETA_BARREL/2;

  double r0slice_end      = Rin/sin(thC_end);
  double y0slice_end      = r0slice_end*tan(D_THETA_BARREL/2.);
  double slice_front_jut  = y0slice_end*sin(M_PI/2-thC_end);

  double z1slice          = Rin -slice_front_jut;
  double z2slice          = Rin +Fdz +Rdz +slice_front_jut;
  double zheight_slice    = (z2slice-z1slice)/2;

  double y1slice          = z1slice*tan(M_PI/2-THETA_SIZE_ENDCAP);
  double y2slice          = z2slice*tan(M_PI/2-THETA_SIZE_ENDCAP);

  // Timing layer envelope
  double  rT      = z1slice -2*nomth;
  double  wT      = rT *tan(D_PHI_GLOBAL/2);
  int     nTiles  = ceil(y1slice/wT);
  double  lT      = 2*y1slice/nTiles;
  int     nCy     = floor(lT/nomth);
  double  actY    = lT/nCy;
  double  actX    = 2*wT/nCy; 

  //-----------------------------------------------------------------------------------
  // Reco struct
  //-----------------------------------------------------------------------------------
  typedef LayeredCalorimeterData::Layer CaloCellData;

  LayeredCalorimeterData* barrelData = new LayeredCalorimeterData;
  barrelData->layoutType = LayeredCalorimeterData::BarrelLayout;
  barrelData->inner_symmetry = PHI_SEGMENTS;
  barrelData->outer_symmetry = PHI_SEGMENTS; 
  barrelData->inner_phi0 = 0.; 
  barrelData->outer_phi0 = 0.; 
  barrelData->gap0 = 0.; //FIXME
  barrelData->gap1 = 0.; //FIXME
  barrelData->gap2 = 0.; //FIXME  
  barrelData->extent[0] =  rT;
  barrelData->extent[1] =  z2slice;
  barrelData->extent[2] = -y2slice;
  barrelData->extent[3] =  y2slice;

  LayeredCalorimeterData* endcapData = new LayeredCalorimeterData;
  endcapData->layoutType = LayeredCalorimeterData::ConicalLayout;
  endcapData->inner_symmetry = PHI_SEGMENTS;
  endcapData->outer_symmetry = PHI_SEGMENTS; 
  endcapData->inner_phi0 = 0.; 
  endcapData->outer_phi0 = 0.; 
  endcapData->gap0 = 0.; //FIXME
  endcapData->gap1 = 0.; //FIXME
  endcapData->gap2 = 0.; //FIXME  
  endcapData->extent[0] =  rT;
  endcapData->extent[1] =  z2slice;
  endcapData->extent[2] = -y2slice;
  endcapData->extent[3] =  y2slice;

  for (int iPhi= CONSTRUCT_BARREL? BARREL_PHI_START:BARREL_PHI_END ; iPhi<BARREL_PHI_END ; iPhi++) {

    double phiEnvBarrel = iPhi*D_PHI_GLOBAL;
    
    //-----------------------------------------------------------------------------------
    // Timing Layer
    //-----------------------------------------------------------------------------------

    dd4hep::Box timingAssemblyShape(nomth,wT,y1slice);
    dd4hep::Volume timingAssemblyVolume("TimingAssembly", timingAssemblyShape, theDetector.material("Vacuum"));
    timingAssemblyVolume.setVisAttributes(theDetector, scepcalAssemblyXML.visStr());

    RotationZYX rotTimingAssembly(0, 0, 0);
    RotationZ rotZPhi(phiEnvBarrel);
    rotTimingAssembly = rotZPhi*rotTimingAssembly;

    double rTimingAssembly = rT+nomth;
    Position dispTimingAssembly(rTimingAssembly*cos(phiEnvBarrel),
                                rTimingAssembly*sin(phiEnvBarrel),
                                0);

    scepcalAssemblyVol.placeVolume( timingAssemblyVolume, Transform3D(rotTimingAssembly, dispTimingAssembly) );

    dd4hep::Box timingCrystalLg(nomth/2,actX/2,lT/2);
    dd4hep::Box timingCrystalTr(nomth/2,wT,actY/2);
    dd4hep::Volume timingCrystalLgVol("TimingCrystalLg", timingCrystalLg, timingLgMat);
    dd4hep::Volume timingCrystalTrVol("TimingCrystalTr", timingCrystalTr, timingTrMat);
    timingCrystalLgVol.setVisAttributes(theDetector, timingLgXML.visStr());
    timingCrystalTrVol.setVisAttributes(theDetector, timingTrXML.visStr());

    for (int nTile=0; nTile<nTiles; nTile++) {

      for (int nC=0; nC<nCy; nC++) {

        int phiEnvBarrelSign = iPhi%2==0? 1:-1;
        int sign             = nTile%2==0? 1:-1;

        RotationZYX rotTiming(0, 0, 0);

        Position dispLg(sign*phiEnvBarrelSign*(nomth/2),
                        -wT +actX/2 + nC*actX,
                        -y1slice +nTile*lT + lT/2
                        );

        Position dispTr(sign*phiEnvBarrelSign*(-nomth/2),
                        0,
                        -y1slice +nTile*lT +actY/2 +nC*actY
                        );

        auto timingLgId64=segmentation->setVolumeID(1,  nTile*nCy +nC , iPhi, 0);
        // auto timingTrId64=segmentation->setVolumeID(1, -nTile*nCy -nC , iPhi, 0);
        auto timingTrId64=segmentation->setVolumeID(1,  nTile*nCy +nC , iPhi, 3);

        int timingLgId32=segmentation->getFirst32bits(timingLgId64);
        int timingTrId32=segmentation->getFirst32bits(timingTrId64);

        // Place volumes and ID them
        dd4hep::PlacedVolume timingLgp = timingAssemblyVolume.placeVolume( timingCrystalLgVol, timingLgId32, Transform3D(rotTiming,dispLg) );
        dd4hep::PlacedVolume timingTrp = timingAssemblyVolume.placeVolume( timingCrystalTrVol, timingTrId32, Transform3D(rotTiming,dispTr) );

        if (!timingLgp.volume().isSensitive()) timingLgp.volume().setSensitiveDetector(sens);
        if (!timingTrp.volume().isSensitive()) timingTrp.volume().setSensitiveDetector(sens);

        Position dispLgglobal = dispTimingAssembly+rotZPhi*dispLg;
        Position dispTrglobal = dispTimingAssembly+rotZPhi*dispTr;

        CaloCellData timingCellLg;
        timingCellLg.cellSize0 = lT;
        timingCellLg.cellSize1 = actX;
        timingCellLg.inner_thickness     = nomth;
        timingCellLg.outer_thickness     = nomth;
        timingCellLg.sensitive_thickness = nomth;
        timingCellLg.inner_nRadiationLengths   = nomth/timingLgMat.radLength();
        timingCellLg.inner_nInteractionLengths = nomth/timingLgMat.intLength();
        timingCellLg.outer_nRadiationLengths   = nomth/timingLgMat.radLength();
        timingCellLg.outer_nInteractionLengths = nomth/timingLgMat.intLength();
        timingCellLg.distance = sqrt(dispLgglobal.mag2()); 
        barrelData->layers.push_back( timingCellLg ) ;

        CaloCellData timingCellTr;
        timingCellTr.cellSize0 = actY;
        timingCellTr.cellSize1 = 2*wT;
        timingCellTr.inner_thickness     = nomth;
        timingCellTr.outer_thickness     = nomth;
        timingCellTr.sensitive_thickness = nomth;
        timingCellTr.inner_nRadiationLengths   = nomth/timingTrMat.radLength();
        timingCellTr.inner_nInteractionLengths = nomth/timingTrMat.intLength();
        timingCellTr.outer_nRadiationLengths   = nomth/timingTrMat.radLength();
        timingCellTr.outer_nInteractionLengths = nomth/timingTrMat.intLength();
        timingCellTr.distance = sqrt(dispTrglobal.mag2()); 
        barrelData->layers.push_back( timingCellTr ) ;

        timingLgp.addPhysVolID("system", 1);
        timingLgp.addPhysVolID("eta", nTile*nCy +nC);
        timingLgp.addPhysVolID("phi", iPhi);
        timingLgp.addPhysVolID("depth", 0);

        timingTrp.addPhysVolID("system", 1);
        // timingTrp.addPhysVolID("eta", -1 -nTile*nCy -nC);
        timingTrp.addPhysVolID("eta", nTile*nCy +nC);

        timingTrp.addPhysVolID("phi", iPhi);
        timingTrp.addPhysVolID("depth", 3);
      }
    }

  }

  // Place the detector
  auto scepcalAssemblyVolId =segmentation->setVolumeID(4,0,0,0);
  int scepcalAssemblyVolId32=segmentation->getFirst32bits(scepcalAssemblyVolId);

  dd4hep::PlacedVolume ScepcalPlacedVol = experimentalHall.placeVolume(scepcalAssemblyVol,scepcalAssemblyVolId32);
  ScepcalPlacedVol.addPhysVolID("system", 4);
  Scepcal.setPlacement(ScepcalPlacedVol);

  // Scepcal.addExtension< LayeredCalorimeterData >( barrelData ) ;

  return Scepcal;
}

DECLARE_DETELEMENT(SCEPCAL_TL, create_detector_TL)
