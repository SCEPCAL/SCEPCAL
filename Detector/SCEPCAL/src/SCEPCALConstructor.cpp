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
create_detector_BE(dd4hep::Detector &theDetector, xml_h xmlElement, dd4hep::SensitiveDetector sens) {

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
  std::cout << "==BARREL AND ENDCAP==" << std::endl;
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

    // dd4hep::Box timingAssemblyShape(nomth,wT,y1slice);
    // dd4hep::Volume timingAssemblyVolume("TimingAssembly", timingAssemblyShape, theDetector.material("Vacuum"));
    // timingAssemblyVolume.setVisAttributes(theDetector, scepcalAssemblyXML.visStr());

    // RotationZYX rotTimingAssembly(0, 0, 0);
    // RotationZ rotZPhi(phiEnvBarrel);
    // rotTimingAssembly = rotZPhi*rotTimingAssembly;

    // double rTimingAssembly = rT+nomth;
    // Position dispTimingAssembly(rTimingAssembly*cos(phiEnvBarrel),
                                // rTimingAssembly*sin(phiEnvBarrel),
                                // 0);

    // scepcalAssemblyVol.placeVolume( timingAssemblyVolume, Transform3D(rotTimingAssembly, dispTimingAssembly) );

    // dd4hep::Box timingCrystalLg(nomth/2,actX/2,lT/2);
    // dd4hep::Box timingCrystalTr(nomth/2,wT,actY/2);
    // dd4hep::Volume timingCrystalLgVol("TimingCrystalLg", timingCrystalLg, timingLgMat);
    // dd4hep::Volume timingCrystalTrVol("TimingCrystalTr", timingCrystalTr, timingTrMat);
    // timingCrystalLgVol.setVisAttributes(theDetector, timingLgXML.visStr());
    // timingCrystalTrVol.setVisAttributes(theDetector, timingTrXML.visStr());

    // for (int nTile=0; nTile<nTiles; nTile++) {

    //   for (int nC=0; nC<nCy; nC++) {

    //     int phiEnvBarrelSign = iPhi%2==0? 1:-1;
    //     int sign             = nTile%2==0? 1:-1;

    //     RotationZYX rotTiming(0, 0, 0);

    //     Position dispLg(sign*phiEnvBarrelSign*(nomth/2),
    //                     -wT +actX/2 + nC*actX,
    //                     -y1slice +nTile*lT + lT/2
    //                     );

    //     Position dispTr(sign*phiEnvBarrelSign*(-nomth/2),
    //                     0,
    //                     -y1slice +nTile*lT +actY/2 +nC*actY
    //                     );

    //     auto timingLgId64=segmentation->setVolumeID(1,  nTile*nCy +nC , iPhi, 0);
    //     // auto timingTrId64=segmentation->setVolumeID(1, -nTile*nCy -nC , iPhi, 0);
    //     auto timingTrId64=segmentation->setVolumeID(1,  nTile*nCy +nC , iPhi, 3);

    //     int timingLgId32=segmentation->getFirst32bits(timingLgId64);
    //     int timingTrId32=segmentation->getFirst32bits(timingTrId64);

    //     // Place volumes and ID them
    //     dd4hep::PlacedVolume timingLgp = timingAssemblyVolume.placeVolume( timingCrystalLgVol, timingLgId32, Transform3D(rotTiming,dispLg) );
    //     dd4hep::PlacedVolume timingTrp = timingAssemblyVolume.placeVolume( timingCrystalTrVol, timingTrId32, Transform3D(rotTiming,dispTr) );

    //     if (!timingLgp.volume().isSensitive()) timingLgp.volume().setSensitiveDetector(sens);
    //     if (!timingTrp.volume().isSensitive()) timingTrp.volume().setSensitiveDetector(sens);

    //     Position dispLgglobal = dispTimingAssembly+rotZPhi*dispLg;
    //     Position dispTrglobal = dispTimingAssembly+rotZPhi*dispTr;

    //     CaloCellData timingCellLg;
    //     timingCellLg.cellSize0 = lT;
    //     timingCellLg.cellSize1 = actX;
    //     timingCellLg.inner_thickness     = nomth;
    //     timingCellLg.outer_thickness     = nomth;
    //     timingCellLg.sensitive_thickness = nomth;
    //     timingCellLg.inner_nRadiationLengths   = nomth/timingLgMat.radLength();
    //     timingCellLg.inner_nInteractionLengths = nomth/timingLgMat.intLength();
    //     timingCellLg.outer_nRadiationLengths   = nomth/timingLgMat.radLength();
    //     timingCellLg.outer_nInteractionLengths = nomth/timingLgMat.intLength();
    //     timingCellLg.distance = sqrt(dispLgglobal.mag2()); 
    //     barrelData->layers.push_back( timingCellLg ) ;

    //     CaloCellData timingCellTr;
    //     timingCellTr.cellSize0 = actY;
    //     timingCellTr.cellSize1 = 2*wT;
    //     timingCellTr.inner_thickness     = nomth;
    //     timingCellTr.outer_thickness     = nomth;
    //     timingCellTr.sensitive_thickness = nomth;
    //     timingCellTr.inner_nRadiationLengths   = nomth/timingTrMat.radLength();
    //     timingCellTr.inner_nInteractionLengths = nomth/timingTrMat.intLength();
    //     timingCellTr.outer_nRadiationLengths   = nomth/timingTrMat.radLength();
    //     timingCellTr.outer_nInteractionLengths = nomth/timingTrMat.intLength();
    //     timingCellTr.distance = sqrt(dispTrglobal.mag2()); 
    //     barrelData->layers.push_back( timingCellTr ) ;

    //     timingLgp.addPhysVolID("system", 1);
    //     timingLgp.addPhysVolID("eta", nTile*nCy +nC);
    //     timingLgp.addPhysVolID("phi", iPhi);
    //     timingLgp.addPhysVolID("depth", 0);

    //     timingTrp.addPhysVolID("system", 1);
    //     // timingTrp.addPhysVolID("eta", -1 -nTile*nCy -nC);
    //     timingTrp.addPhysVolID("eta", nTile*nCy +nC);

    //     timingTrp.addPhysVolID("phi", iPhi);
    //     timingTrp.addPhysVolID("depth", 3);
    //   }
    // }

    //-----------------------------------------------------------------------------------
    // Barrel Crystals
    //-----------------------------------------------------------------------------------

    for (int nGamma=0; nGamma<N_PHI_BARREL_CRYSTAL; nGamma++) {

      double gamma = -D_PHI_GLOBAL/2+D_PHI_BARREL_CRYSTAL/2+D_PHI_BARREL_CRYSTAL*nGamma;

      // Make assembly slice
      double x0y0le = z1slice*tan(gamma-D_PHI_BARREL_CRYSTAL/2);
      double x0y0re = z1slice*tan(gamma+D_PHI_BARREL_CRYSTAL/2);

      double x1y1le = z2slice*tan(gamma-D_PHI_BARREL_CRYSTAL/2);
      double x1y1re = z2slice*tan(gamma+D_PHI_BARREL_CRYSTAL/2);

      double verticesS[]={
                          -y1slice, x0y0le,
                          -y1slice, x0y0re,
                           y1slice, x0y0re,
                           y1slice, x0y0le,

                          -y2slice, x1y1le,
                          -y2slice, x1y1re,
                           y2slice, x1y1re,
                           y2slice, x1y1le
                         };
      
      dd4hep::EightPointSolid barrelSliceAssemblyShape(zheight_slice,verticesS);
      dd4hep::Volume barrelSliceAssemblyVolume("BarrelSliceAssembly", barrelSliceAssemblyShape, theDetector.material("Vacuum"));
      barrelSliceAssemblyVolume.setVisAttributes(theDetector, scepcalAssemblyXML.visStr());

      RotationZYX rotSlice(0, M_PI/2, 0);
      RotationZ rotZSlice(phiEnvBarrel);
      rotSlice = rotZSlice*rotSlice;

      double rSlice =(z1slice+z2slice)/2;

      Position dispSlice(rSlice*cos(phiEnvBarrel),
                          rSlice*sin(phiEnvBarrel),
                          0);

      scepcalAssemblyVol.placeVolume( barrelSliceAssemblyVolume, Transform3D(rotSlice, dispSlice) );

      for (int iTheta=0; iTheta<N_THETA_BARREL; iTheta++) {
    
        if (debugLevel>1) std::cout << "  Barrel: theta: " << iTheta << std::endl;

        double thC =THETA_SIZE_ENDCAP+D_THETA_BARREL/2+(iTheta*D_THETA_BARREL);

        double r0e=Rin/sin(thC);
        double r1e=r0e+Fdz;
        double r2e=r1e+Rdz;

        double y0e=r0e*tan(D_THETA_BARREL/2.);
        double y1e=r1e*tan(D_THETA_BARREL/2.);
        double y2e=r2e*tan(D_THETA_BARREL/2.);

        // Make crystal shapes
        // Projective trapezoids using EightPointSolids. see: https://root.cern.ch/doc/master/classTGeoArb8.html

        double x0y0 = (r0e*cos(thC) +y0e*sin(thC)) *tan(thC -D_THETA_BARREL/2.);
        double x1y0 = (r0e*cos(thC) -y0e*sin(thC)) *tan(thC +D_THETA_BARREL/2.);

        double x0y0l = x0y0*tan(gamma-D_PHI_BARREL_CRYSTAL/2);
        double x0y0r = x0y0*tan(gamma+D_PHI_BARREL_CRYSTAL/2);

        double x1y0l = x1y0*tan(gamma-D_PHI_BARREL_CRYSTAL/2);
        double x1y0r = x1y0*tan(gamma+D_PHI_BARREL_CRYSTAL/2);

        double x0y1 = (r1e*cos(thC) +y1e*sin(thC)) *tan(thC -D_THETA_BARREL/2.);
        double x1y1 = (r1e*cos(thC) -y1e*sin(thC)) *tan(thC +D_THETA_BARREL/2.);

        double x0y1l = x0y1*tan(gamma-D_PHI_BARREL_CRYSTAL/2);
        double x0y1r = x0y1*tan(gamma+D_PHI_BARREL_CRYSTAL/2);

        double x1y1l = x1y1*tan(gamma-D_PHI_BARREL_CRYSTAL/2);
        double x1y1r = x1y1*tan(gamma+D_PHI_BARREL_CRYSTAL/2);

        double x0y2 = (r2e*cos(thC) +y2e*sin(thC)) *tan(thC -D_THETA_BARREL/2.);
        double x1y2 = (r2e*cos(thC) -y2e*sin(thC)) *tan(thC +D_THETA_BARREL/2.);

        double x0y2l = x0y2*tan(gamma-D_PHI_BARREL_CRYSTAL/2);
        double x0y2r = x0y2*tan(gamma+D_PHI_BARREL_CRYSTAL/2);

        double x1y2l = x1y2*tan(gamma-D_PHI_BARREL_CRYSTAL/2);
        double x1y2r = x1y2*tan(gamma+D_PHI_BARREL_CRYSTAL/2);

        // Front (F) and Rear (R) crystal shapes
        double verticesF[]={
                            x0y0r,  y0e,
                            x1y0r, -y0e,
                            x1y0l, -y0e,
                            x0y0l,  y0e,

                            x0y1r,  y1e,
                            x1y1r, -y1e,
                            x1y1l, -y1e,
                            x0y1l,  y1e
                           };

        double verticesR[]={
                            x0y1r,  y1e,
                            x1y1r, -y1e,
                            x1y1l, -y1e,
                            x0y1l,  y1e,

                            x0y2r,  y2e,
                            x1y2r, -y2e,
                            x1y2l, -y2e,
                            x0y2l,  y2e
                           };

        dd4hep::EightPointSolid crystalFShape(Fdz/2, verticesF);
        dd4hep::EightPointSolid crystalRShape(Rdz/2, verticesR);

        // Promote shapes to volumes and set attributes
        dd4hep::Volume crystalFVol("BarrelCrystalF", crystalFShape, crystalFMat);
        dd4hep::Volume crystalRVol("BarrelCrystalR", crystalRShape, crystalRMat);

        crystalFVol.setVisAttributes(theDetector, crystalFXML.visStr());
        crystalRVol.setVisAttributes(theDetector, crystalRXML.visStr());

        RotationZYX rot(M_PI/2, -M_PI/2 +thC, 0);

        double rF=r0e+Fdz/2.;
        Position dispF(-rF*cos(thC),
                      0,
                      -(rSlice-rF*sin(thC))
                      );

        double rR=r1e+Rdz/2.;
        Position dispR(-rR*cos(thC),
                      0,
                      -(rSlice-rR*sin(thC))
                      );

        auto crystalFId64=segmentation->setVolumeID(1, N_THETA_ENDCAP+iTheta , iPhi*N_PHI_BARREL_CRYSTAL+nGamma, 1);
        auto crystalRId64=segmentation->setVolumeID(1, N_THETA_ENDCAP+iTheta , iPhi*N_PHI_BARREL_CRYSTAL+nGamma, 2);
        int crystalFId32=segmentation->getFirst32bits(crystalFId64);
        int crystalRId32=segmentation->getFirst32bits(crystalRId64);

        // Place volumes and ID them
        dd4hep::PlacedVolume crystalFp = barrelSliceAssemblyVolume.placeVolume( crystalFVol, crystalFId32, Transform3D(rot,dispF) );
        dd4hep::PlacedVolume crystalRp = barrelSliceAssemblyVolume.placeVolume( crystalRVol, crystalRId32, Transform3D(rot,dispR) );

        if (!crystalFp.volume().isSensitive()) crystalFp.volume().setSensitiveDetector(sens);
        if (!crystalRp.volume().isSensitive()) crystalRp.volume().setSensitiveDetector(sens);

        Position dispFgamma(0, rF*sin(thC)*tan(gamma), 0);
        Position dispFslice(rF*sin(thC)*cos(phiEnvBarrel),
                    rF*sin(thC)*sin(phiEnvBarrel),
                    rF*cos(thC));
        Position dispFglobal(dispFslice+rotZSlice*dispFgamma);

        Position dispRgamma(0, rR*sin(thC)*tan(gamma), 0);
        Position dispRslice(rR*sin(thC)*cos(phiEnvBarrel),
                    rR*sin(thC)*sin(phiEnvBarrel),
                    rR*cos(thC));
        Position dispRglobal(dispRslice+rotZSlice*dispRgamma);

        CaloCellData barrelCellF;
        barrelCellF.cellSize0 = 2*y0e;
        barrelCellF.cellSize1 = x0y0l+x0y0r;
        barrelCellF.inner_thickness     = Fdz;
        barrelCellF.outer_thickness     = Fdz;
        barrelCellF.sensitive_thickness = Fdz;
        barrelCellF.inner_nRadiationLengths   = Fdz/crystalFMat.radLength();
        barrelCellF.inner_nInteractionLengths = Fdz/crystalFMat.intLength();
        barrelCellF.outer_nRadiationLengths   = Fdz/crystalFMat.radLength();
        barrelCellF.outer_nInteractionLengths = Fdz/crystalFMat.intLength();
        barrelCellF.distance = sqrt(dispFglobal.mag2()); 
        barrelData->layers.push_back( barrelCellF ) ;

        CaloCellData barrelCellR;
        barrelCellR.cellSize0 = 2*y1e;
        barrelCellR.cellSize1 = x0y1l+x0y1r;
        barrelCellR.inner_thickness     = Rdz;
        barrelCellR.outer_thickness     = Rdz;
        barrelCellR.sensitive_thickness = Rdz;
        barrelCellR.inner_nRadiationLengths   = Rdz/crystalRMat.radLength();
        barrelCellR.inner_nInteractionLengths = Rdz/crystalRMat.intLength();
        barrelCellR.outer_nRadiationLengths   = Rdz/crystalRMat.radLength();
        barrelCellR.outer_nInteractionLengths = Rdz/crystalRMat.intLength();
        barrelCellR.distance = sqrt(dispRglobal.mag2()); 
        barrelData->layers.push_back( barrelCellR ) ;

        crystalFp.addPhysVolID("system", 1);
        crystalFp.addPhysVolID("eta", N_THETA_ENDCAP+iTheta);
        crystalFp.addPhysVolID("phi", iPhi*N_PHI_BARREL_CRYSTAL+nGamma);
        crystalFp.addPhysVolID("depth", 1);
        
        crystalRp.addPhysVolID("system", 1);
        crystalRp.addPhysVolID("eta", N_THETA_ENDCAP+iTheta);
        crystalRp.addPhysVolID("phi", iPhi*N_PHI_BARREL_CRYSTAL+nGamma);
        crystalRp.addPhysVolID("depth", 2);
      }
    }
  }

  //-----------------------------------------------------------------------------------
  // Endcap
  //-----------------------------------------------------------------------------------

  for (int iTheta=CONSTRUCT_ENDCAP? ENDCAP_THETA_START:N_THETA_ENDCAP ; iTheta<N_THETA_ENDCAP ; iTheta++) {

    double thC        = D_THETA_ENDCAP/2+ iTheta*D_THETA_ENDCAP;
    double RinEndcap  = EBz*tan(thC);

    int    nPhiEndcapCrystal = floor(2*M_PI*RinEndcap/(PHI_SEGMENTS*nomfw));
    double dPhiEndcapCrystal = D_PHI_GLOBAL/nPhiEndcapCrystal;

    double r0e=RinEndcap/sin(thC);
    double r1e=r0e+Fdz;
    double r2e=r1e+Rdz;

    double y0e=r0e*tan(D_THETA_ENDCAP/2.);
    double y1e=r1e*tan(D_THETA_ENDCAP/2.);
    double y2e=r2e*tan(D_THETA_ENDCAP/2.);

    // Make assembly polyhedra
    double a = r0e/cos(D_THETA_ENDCAP/2);
    double z1 = a*cos(thC+D_THETA_ENDCAP/2);
    double r1min = z1*tan(thC-D_THETA_ENDCAP/2);
    double r1max = z1*tan(thC+D_THETA_ENDCAP/2);

    double b = sqrt(r2e*r2e +y2e*y2e);
    double z2 = b*cos(thC-D_THETA_ENDCAP/2);
    double r2min = z2*tan(thC-D_THETA_ENDCAP/2);
    double r2max = z2*tan(thC+D_THETA_ENDCAP/2);

    std::vector<double> zPolyhedra = {z1,z2};
    std::vector<double> rminPolyhedra = {r1min, r2min};
    std::vector<double> rmaxPolyhedra = {r1max, r2max};

    dd4hep::Polyhedra phiRingAssemblyShape(PHI_SEGMENTS, D_PHI_GLOBAL/2, 2*M_PI, zPolyhedra, rminPolyhedra, rmaxPolyhedra);

    Position  dispCone(0,0,0);
    Position  dispCone1(0,0,0);
    RotationY rotMirror(M_PI);

    // Endcap assembly volume
    dd4hep::Volume phiRingAssemblyVolume("EndcapRingAssembly", phiRingAssemblyShape, theDetector.material("Vacuum"));
    phiRingAssemblyVolume.setVisAttributes(theDetector, scepcalAssemblyXML.visStr());
    scepcalAssemblyVol.placeVolume( phiRingAssemblyVolume, dispCone );

    dd4hep::Volume* phiRingAssemblyVolume1 = nullptr;

    if (REFLECT_ENDCAP) {
      phiRingAssemblyVolume1 = new dd4hep::Volume("EndcapRingAssembly1", phiRingAssemblyShape, theDetector.material("Vacuum"));
      phiRingAssemblyVolume1->setVisAttributes(theDetector, scepcalAssemblyXML.visStr());
      scepcalAssemblyVol.placeVolume( *phiRingAssemblyVolume1, Transform3D(rotMirror, dispCone1) );
    }

    for (int iPhi=ENDCAP_PHI_START ; iPhi<ENDCAP_PHI_END ; iPhi++) {
      
      for (int nGamma=0; nGamma<nPhiEndcapCrystal; nGamma++) {

        double gamma = -D_PHI_GLOBAL/2+dPhiEndcapCrystal/2+dPhiEndcapCrystal*nGamma;

        double x0y0 = (r0e*cos(thC) +y0e*sin(thC)) *tan(thC -D_THETA_ENDCAP/2.);
        double x1y0 = (r0e*cos(thC) -y0e*sin(thC)) *tan(thC +D_THETA_ENDCAP/2.);

        double x0y0l = x0y0*tan(gamma-dPhiEndcapCrystal/2);
        double x0y0r = x0y0*tan(gamma+dPhiEndcapCrystal/2);

        double x1y0l = x1y0*tan(gamma-dPhiEndcapCrystal/2);
        double x1y0r = x1y0*tan(gamma+dPhiEndcapCrystal/2);

        double x0y1 = (r1e*cos(thC) +y1e*sin(thC)) *tan(thC -D_THETA_ENDCAP/2.);
        double x1y1 = (r1e*cos(thC) -y1e*sin(thC)) *tan(thC +D_THETA_ENDCAP/2.);

        double x0y1l = x0y1*tan(gamma-dPhiEndcapCrystal/2);
        double x0y1r = x0y1*tan(gamma+dPhiEndcapCrystal/2);

        double x1y1l = x1y1*tan(gamma-dPhiEndcapCrystal/2);
        double x1y1r = x1y1*tan(gamma+dPhiEndcapCrystal/2);

        double x0y2 = (r2e*cos(thC) +y2e*sin(thC)) *tan(thC -D_THETA_ENDCAP/2.);
        double x1y2 = (r2e*cos(thC) -y2e*sin(thC)) *tan(thC +D_THETA_ENDCAP/2.);

        double x0y2l = x0y2*tan(gamma-dPhiEndcapCrystal/2);
        double x0y2r = x0y2*tan(gamma+dPhiEndcapCrystal/2);

        double x1y2l = x1y2*tan(gamma-dPhiEndcapCrystal/2);
        double x1y2r = x1y2*tan(gamma+dPhiEndcapCrystal/2);

        double verticesF[]={
                            x0y0r,  y0e,
                            x1y0r, -y0e,
                            x1y0l, -y0e,
                            x0y0l,  y0e,

                            x0y1r,  y1e,
                            x1y1r, -y1e,
                            x1y1l, -y1e,
                            x0y1l,  y1e
                           };

        double verticesR[]={
                            x0y1r,  y1e,
                            x1y1r, -y1e,
                            x1y1l, -y1e,
                            x0y1l,  y1e,

                            x0y2r,  y2e,
                            x1y2r, -y2e,
                            x1y2l, -y2e,
                            x0y2l,  y2e
                           };

        dd4hep::EightPointSolid crystalFShape(Fdz/2, verticesF);
        dd4hep::EightPointSolid crystalRShape(Rdz/2, verticesR);

        dd4hep::Volume crystalFVol("EndcapCrystalF", crystalFShape, crystalFMat);
        dd4hep::Volume crystalRVol("EndcapCrystalR", crystalRShape, crystalRMat);

        crystalFVol.setVisAttributes(theDetector, crystalFXML.visStr());
        crystalRVol.setVisAttributes(theDetector, crystalRXML.visStr());
        
        double phiEnvEndcap=iPhi*D_PHI_GLOBAL;

        RotationZYX rot(M_PI/2, thC, 0);
        RotationZ rotZ(phiEnvEndcap);
        rot = rotZ*rot;

        double rF=r0e+Fdz/2.;
        Position dispF(rF*sin(thC)*cos(phiEnvEndcap),
                      rF*sin(thC)*sin(phiEnvEndcap),
                      rF*cos(thC));

        double rR=r1e+Rdz/2.;
        Position dispR(rR*sin(thC)*cos(phiEnvEndcap),
                      rR*sin(thC)*sin(phiEnvEndcap),
                      rR*cos(thC));

        auto crystalFId64=segmentation->setVolumeID(1, iTheta, iPhi*nPhiEndcapCrystal+nGamma, 1);
        auto crystalRId64=segmentation->setVolumeID(1, iTheta, iPhi*nPhiEndcapCrystal+nGamma, 2);
        int crystalFId32=segmentation->getFirst32bits(crystalFId64);
        int crystalRId32=segmentation->getFirst32bits(crystalRId64);

        dd4hep::PlacedVolume crystalFp = phiRingAssemblyVolume.placeVolume( crystalFVol, crystalFId32, Transform3D(rot,dispF-dispCone) );
        dd4hep::PlacedVolume crystalRp = phiRingAssemblyVolume.placeVolume( crystalRVol, crystalRId32, Transform3D(rot,dispR-dispCone) );

        if (!crystalFp.volume().isSensitive()) crystalFp.volume().setSensitiveDetector(sens);
        if (!crystalRp.volume().isSensitive()) crystalRp.volume().setSensitiveDetector(sens);

        Position dispFgamma(0, rF*sin(thC)*tan(gamma), 0);
        Position dispFslice(rF*sin(thC)*cos(phiEnvEndcap),
                    rF*sin(thC)*sin(phiEnvEndcap),
                    rF*cos(thC));
        Position dispFglobal(dispFslice+rotZ*dispFgamma);

        Position dispRgamma(0, rR*sin(thC)*tan(gamma), 0);
        Position dispRslice(rR*sin(thC)*cos(phiEnvEndcap),
                    rR*sin(thC)*sin(phiEnvEndcap),
                    rR*cos(thC));
        Position dispRglobal(dispRslice+rotZ*dispRgamma);

        CaloCellData endcapCellF;
        endcapCellF.cellSize0 = 2*y0e;
        endcapCellF.cellSize1 = x0y0l+x0y0r;
        endcapCellF.inner_thickness     = Fdz;
        endcapCellF.outer_thickness     = Fdz;
        endcapCellF.sensitive_thickness = Fdz;
        endcapCellF.inner_nRadiationLengths   = Fdz/crystalFMat.radLength();
        endcapCellF.inner_nInteractionLengths = Fdz/crystalFMat.intLength();
        endcapCellF.outer_nRadiationLengths   = Fdz/crystalFMat.radLength();
        endcapCellF.outer_nInteractionLengths = Fdz/crystalFMat.intLength();
        endcapCellF.distance = sqrt(dispFglobal.mag2()); 
        barrelData->layers.push_back( endcapCellF ) ;

        CaloCellData endcapCellR;
        endcapCellR.cellSize0 = 2*y1e;
        endcapCellR.cellSize1 = x0y1l+x0y1r;
        endcapCellR.inner_thickness     = Rdz;
        endcapCellR.outer_thickness     = Rdz;
        endcapCellR.sensitive_thickness = Rdz;
        endcapCellR.inner_nRadiationLengths   = Rdz/crystalRMat.radLength();
        endcapCellR.inner_nInteractionLengths = Rdz/crystalRMat.intLength();
        endcapCellR.outer_nRadiationLengths   = Rdz/crystalRMat.radLength();
        endcapCellR.outer_nInteractionLengths = Rdz/crystalRMat.intLength();
        endcapCellR.distance = sqrt(dispRglobal.mag2()); 
        barrelData->layers.push_back( endcapCellR ) ;

        crystalFp.addPhysVolID("system", 1);
        crystalFp.addPhysVolID("eta", iTheta);
        crystalFp.addPhysVolID("phi", iPhi*nPhiEndcapCrystal+nGamma);
        crystalFp.addPhysVolID("depth", 1);

        crystalRp.addPhysVolID("system", 1);
        crystalRp.addPhysVolID("eta", iTheta);
        crystalRp.addPhysVolID("phi", iPhi*nPhiEndcapCrystal+nGamma);
        crystalRp.addPhysVolID("depth", 2);

        if (REFLECT_ENDCAP) {
          auto crystalFId641=segmentation->setVolumeID(1, N_THETA_ENDCAP+N_THETA_BARREL+N_THETA_ENDCAP-iTheta, iPhi*nPhiEndcapCrystal+nGamma, 1);
          auto crystalRId641=segmentation->setVolumeID(1, N_THETA_ENDCAP+N_THETA_BARREL+N_THETA_ENDCAP-iTheta, iPhi*nPhiEndcapCrystal+nGamma, 2);
          
          int crystalFId321=segmentation->getFirst32bits(crystalFId641);
          int crystalRId321=segmentation->getFirst32bits(crystalRId641);

          dd4hep::PlacedVolume crystalFp1 = phiRingAssemblyVolume1->placeVolume( crystalFVol, crystalFId321, Transform3D(rot,dispF-dispCone) );
          dd4hep::PlacedVolume crystalRp1 = phiRingAssemblyVolume1->placeVolume( crystalRVol, crystalRId321, Transform3D(rot,dispR-dispCone) );

          if (!crystalFp1.volume().isSensitive()) crystalFp1.volume().setSensitiveDetector(sens);
          if (!crystalRp1.volume().isSensitive()) crystalRp1.volume().setSensitiveDetector(sens);

          crystalFp1.addPhysVolID("system", 1);
          crystalFp1.addPhysVolID("eta", N_THETA_ENDCAP+N_THETA_BARREL+N_THETA_ENDCAP-iTheta);
          crystalFp1.addPhysVolID("phi", iPhi*nPhiEndcapCrystal+nGamma);
          crystalFp1.addPhysVolID("depth", 1);

          crystalRp1.addPhysVolID("system", 1);
          crystalRp1.addPhysVolID("eta", N_THETA_ENDCAP+N_THETA_BARREL+N_THETA_ENDCAP-iTheta);
          crystalRp1.addPhysVolID("phi", iPhi*nPhiEndcapCrystal+nGamma);
          crystalRp1.addPhysVolID("depth", 2);

          Position dispFglobal1(rotMirror*(dispFslice+rotZ*dispFgamma));
          Position dispRglobal1(rotMirror*(dispRslice+rotZ*dispRgamma));

          CaloCellData endcapCellF1;
          endcapCellF1.cellSize0 = 2*y0e;
          endcapCellF1.cellSize1 = x0y0l+x0y0r;
          endcapCellF1.inner_thickness     = Fdz;
          endcapCellF1.outer_thickness     = Fdz;
          endcapCellF1.sensitive_thickness = Fdz;
          endcapCellF1.inner_nRadiationLengths   = Fdz/crystalFMat.radLength();
          endcapCellF1.inner_nInteractionLengths = Fdz/crystalFMat.intLength();
          endcapCellF1.outer_nRadiationLengths   = Fdz/crystalFMat.radLength();
          endcapCellF1.outer_nInteractionLengths = Fdz/crystalFMat.intLength();
          endcapCellF1.distance = sqrt(dispFglobal1.mag2()); 
          barrelData->layers.push_back( endcapCellF1 ) ;

          CaloCellData endcapCellR1;
          endcapCellR1.cellSize0 = 2*y1e;
          endcapCellR1.cellSize1 = x0y1l+x0y1r;
          endcapCellR1.inner_thickness     = Rdz;
          endcapCellR1.outer_thickness     = Rdz;
          endcapCellR1.sensitive_thickness = Rdz;
          endcapCellR1.inner_nRadiationLengths   = Rdz/crystalRMat.radLength();
          endcapCellR1.inner_nInteractionLengths = Rdz/crystalRMat.intLength();
          endcapCellR1.outer_nRadiationLengths   = Rdz/crystalRMat.radLength();
          endcapCellR1.outer_nInteractionLengths = Rdz/crystalRMat.intLength();
          endcapCellR1.distance = sqrt(dispRglobal1.mag2()); 
          barrelData->layers.push_back( endcapCellR1 ) ;

        }
      }
    }
  }

  // Place the detector
  auto scepcalAssemblyVolId =segmentation->setVolumeID(3,0,0,0);
  int scepcalAssemblyVolId32=segmentation->getFirst32bits(scepcalAssemblyVolId);

  dd4hep::PlacedVolume ScepcalPlacedVol = experimentalHall.placeVolume(scepcalAssemblyVol,scepcalAssemblyVolId32);
  ScepcalPlacedVol.addPhysVolID("system", 3);
  Scepcal.setPlacement(ScepcalPlacedVol);

  // Scepcal.addExtension< LayeredCalorimeterData >( barrelData ) ;

  return Scepcal;
}

DECLARE_DETELEMENT(SCEPCAL_BE, create_detector_BE)
