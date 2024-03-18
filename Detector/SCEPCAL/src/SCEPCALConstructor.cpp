//=========================================================================
// Author: Wonyong Chung
//=========================================================================

#define VERBOSE_LEVEL 0

#include "SCEPCALSegmentationHandle.h"

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DetectorTools.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Detector.h"

#include "TGeoTrd2.h"
#include <bitset>

using dd4hep::Transform3D;
using dd4hep::RotationZYX;
using dd4hep::Rotation3D;
using dd4hep::Position;
//using dd4hep::PlacedVolume::VolIDs;

namespace ddSCEPCAL {
    static dd4hep::Ref_t
    create_detector(dd4hep::Detector &theDetector, xml_h xmlElement, dd4hep::SensitiveDetector sensDet) {

      // Initialize detector element
      xml_det_t detectorXML=xmlElement;
      std::string name=detectorXML.nameStr();
      dd4hep::DetElement drDet(name, detectorXML.id());

      // Import xml objects from compact xml
      dd4hep::xml::Dimension sensDetType=detectorXML.child(_Unicode(sensitive));
      xml_comp_t dimXML=detectorXML.child(_Unicode(dim));
      xml_comp_t crystalFXML=detectorXML.child(_Unicode(crystalF));
      xml_comp_t crystalRXML=detectorXML.child(_Unicode(crystalR));
      xml_comp_t timingTrXML=detectorXML.child(_Unicode(timingLayerTr));
      xml_comp_t timingLgXML=detectorXML.child(_Unicode(timingLayerLg));
      xml_comp_t instXML=detectorXML.child(_Unicode(inst));
      xml_comp_t scepcalAssemblyGlobalVisXML=detectorXML.child(_Unicode(scepcalAssemblyGlobalVis));
      xml_comp_t scepcalAssemblyXML=detectorXML.child(_Unicode(scepcalAssembly));

      xml_comp_t printDebugXML=detectorXML.child(_Unicode(printDebug));

      //-----------------------------------------------------------------------------------

      sensDet.setType(sensDetType.typeStr());

      // Set the segmentation class
      auto segmentation=dynamic_cast<dd4hep::DDSegmentation::SCEPCALSegmentation *>( sensDet.readout().segmentation().segmentation());

      dd4hep::Assembly experimentalHall("hall");

      const int  debugLevel = printDebugXML.attr<bool>(_Unicode(level));

      // Parse input parameters from imported xml objects
      const double Fdz    =crystalFXML.attr<double>(_Unicode(length));
      const double Rdz    =crystalRXML.attr<double>(_Unicode(length));
      const double nomfw  =dimXML.attr<double>(_Unicode(crystalFaceWidthNominal));
      const double nomth  =dimXML.attr<double>(_Unicode(crystalTimingThicknessNominal));
      const double EBz    =dimXML.attr<double>(_Unicode(barrelHalfZ));
      const double Rin    =dimXML.attr<double>(_Unicode(barrelInnerR));

      // Material definitions
      dd4hep::Material crystalFMat =theDetector.material(crystalFXML.materialStr());
      dd4hep::Material crystalRMat =theDetector.material(crystalRXML.materialStr());
      dd4hep::Material timingTrMat =theDetector.material(timingTrXML.materialStr());
      dd4hep::Material timingLgMat =theDetector.material(timingLgXML.materialStr());
      dd4hep::Material instMat     =theDetector.material(instXML.materialStr());

      // Begin geometry calculations

      // Need odd number of nThetaBarrel to make center crystal
      int nThetaBarrel  =int(floor(2*EBz/nomfw))%2==1? floor(2*EBz/nomfw) : floor(2*EBz/nomfw)-1 ;
      int nThetaEndcap  =floor(Rin/nomfw);

      double thetaSizeEndcap=atan(Rin/EBz);

      double dThetaBarrel =(M_PI-2*thetaSizeEndcap)/(nThetaBarrel);
      double dThetaEndcap =thetaSizeEndcap/nThetaEndcap;

      int    nPhiBarrel = 16;
      double dPhiBarrel = 2*M_PI/nPhiBarrel;
      
      int    nPhiBarrelCrystal   =floor(2*M_PI*Rin/(nPhiBarrel*nomfw));
      double dPhiBarrelCrystal   =dPhiBarrel/nPhiBarrelCrystal;

      dd4hep::Box scepcalAssemblyShape(Rin+2*Fdz+2*Rdz, Rin+2*Fdz+2*Rdz, 2*EBz+2*Fdz+2*Rdz);
      dd4hep::Volume scepcalAssemblyVol("scepcalAssemblyVol", scepcalAssemblyShape, theDetector.material("Vacuum"));
      scepcalAssemblyVol.setVisAttributes(theDetector, scepcalAssemblyGlobalVisXML.visStr());

      /** 
       * 
       * Barrel
       * 
       * **/

      for (int iPhi=0; iPhi<nPhiBarrel; iPhi++) {

        if (debugLevel>1) std::cout << "Barrel: phi: " << iPhi << std::endl;

        // Shared calculations for timing layer and barrel envelopes
        double thC_end = thetaSizeEndcap+dThetaBarrel/2;

        double r0slice_end =Rin/sin(thC_end);
        double y0slice_end =r0slice_end*tan(dThetaBarrel/2.);
        double slice_front_jut = y0slice_end*sin(M_PI/2-thC_end);

        double z1slice =Rin -slice_front_jut;
        double z2slice =Rin +Fdz +Rdz +slice_front_jut;
        double zheight_slice = (z2slice-z1slice)/2;

        double y1slice =z1slice*tan(M_PI/2-thetaSizeEndcap);
        double y2slice =z2slice*tan(M_PI/2-thetaSizeEndcap);

        double phi=iPhi*dPhiBarrel;


        // Timing layer

        double rT = z1slice -2*nomth;
        double w  = rT *tan(dPhiBarrel/2);
        int nTiles= ceil(y1slice/w);
        double lT = 2*y1slice/nTiles;
        int nCy = floor(lT/nomth);
        double actY = lT/nCy;
        double actX = 2*w/nCy; 

        dd4hep::Box timingAssemblyShape(nomth,w,y1slice);
        dd4hep::Volume timingAssemblyVolume("TimingAssembly", timingAssemblyShape, theDetector.material("Vacuum"));
        timingAssemblyVolume.setVisAttributes(theDetector, scepcalAssemblyXML.visStr());

        RotationZYX rotTimingAssembly(0, 0, 0);
        ROOT::Math::RotationZ rotZPhi = ROOT::Math::RotationZ(phi);
        rotTimingAssembly = rotZPhi*rotTimingAssembly;

        double rTimingAssembly = rT+nomth;
        Position dispTimingAssembly(rTimingAssembly*cos(phi),
                                    rTimingAssembly*sin(phi),
                                    0);

        scepcalAssemblyVol.placeVolume( timingAssemblyVolume, Transform3D(rotTimingAssembly, dispTimingAssembly) );

        dd4hep::Box timingCrystalTr(nomth/2,w,actY/2);
        dd4hep::Box timingCrystalLg(nomth/2,actX/2,lT/2);
        dd4hep::Volume timingCrystalTrVol("TimingCrystalTr", timingCrystalTr, timingTrMat);
        dd4hep::Volume timingCrystalLgVol("TimingCrystalLg", timingCrystalLg, timingLgMat);
        timingCrystalTrVol.setVisAttributes(theDetector, timingTrXML.visStr());
        timingCrystalLgVol.setVisAttributes(theDetector, timingLgXML.visStr());

        for (int nTile=0; nTile<nTiles; nTile++) {

          for (int nC=0; nC<nCy; nC++) {

            int phiSign =  iPhi%2==0? 1:-1;
            int sign    = nTile%2==0? 1:-1;

            RotationZYX rotTiming(0, 0, 0);

            Position dispLg(sign*phiSign*(nomth/2),
                            -w +actX/2 + nC*actX,
                            -y1slice +nTile*lT + lT/2
                            );

            Position dispTr(sign*phiSign*(-nomth/2),
                            0,
                            -y1slice +nTile*lT +actY/2 +nC*actY
                            );

            auto timingLgId64=segmentation->setVolumeID(1, nTile*nCy +nC , iPhi, 1);
            auto timingTrId64=segmentation->setVolumeID(1, -nTile*nCy -nC , iPhi, 2);
            int timingLgId32=segmentation->getFirst32bits(timingLgId64);
            int timingTrId32=segmentation->getFirst32bits(timingTrId64);

            // Place volumes and ID them
            dd4hep::PlacedVolume timingLgp = timingAssemblyVolume.placeVolume( timingCrystalLgVol, timingLgId32, Transform3D(rotTiming,dispLg) );
            dd4hep::PlacedVolume timingTrp = timingAssemblyVolume.placeVolume( timingCrystalTrVol, timingTrId32, Transform3D(rotTiming,dispTr) );

            timingLgp.addPhysVolID("eta", nTile*nCy +nC);
            timingLgp.addPhysVolID("phi", iPhi);
            timingLgp.addPhysVolID("depth", 0);
            timingLgp.addPhysVolID("system", 1);
            
            timingTrp.addPhysVolID("eta", -nTile*nCy -nC);
            timingTrp.addPhysVolID("phi", iPhi);
            timingTrp.addPhysVolID("depth", 0);
            timingTrp.addPhysVolID("system", 1);
          }
        }

        for (int nGamma=0; nGamma<nPhiBarrelCrystal; nGamma++) {

          double gamma = -dPhiBarrel/2+dPhiBarrelCrystal/2+dPhiBarrelCrystal*nGamma;

          // Make assembly slice
          double x0y0le = z1slice*tan(gamma-dPhiBarrelCrystal/2);
          double x0y0re = z1slice*tan(gamma+dPhiBarrelCrystal/2);

          double x1y1le = z2slice*tan(gamma-dPhiBarrelCrystal/2);
          double x1y1re = z2slice*tan(gamma+dPhiBarrelCrystal/2);

          double verticesS[]={
                              -y1slice,x0y0le,
                              -y1slice,x0y0re,
                              y1slice,x0y0re,
                              y1slice,x0y0le,

                              -y2slice,x1y1le,
                              -y2slice,x1y1re,
                              y2slice,x1y1re,
                              y2slice,x1y1le
                              };
          
          dd4hep::EightPointSolid barrelSliceAssemblyShape(zheight_slice,verticesS);
          dd4hep::Volume barrelSliceAssemblyVolume("BarrelSliceAssembly", barrelSliceAssemblyShape, theDetector.material("Vacuum"));
          barrelSliceAssemblyVolume.setVisAttributes(theDetector, scepcalAssemblyXML.visStr());

          RotationZYX rotSlice(0, M_PI/2, 0);
          ROOT::Math::RotationZ rotZSlice = ROOT::Math::RotationZ(phi);
          rotSlice = rotZSlice*rotSlice;

          double rSlice =(z1slice+z2slice)/2;

          Position dispSlice(rSlice*cos(phi),
                             rSlice*sin(phi),
                             0);

          scepcalAssemblyVol.placeVolume( barrelSliceAssemblyVolume, Transform3D(rotSlice, dispSlice) );

          for (int iTheta=0; iTheta<nThetaBarrel; iTheta++) {
        
            if (debugLevel>1) std::cout << "  Barrel: theta: " << iTheta << std::endl;

            double thC =thetaSizeEndcap+dThetaBarrel/2+(iTheta*dThetaBarrel);

            double r0e=Rin/sin(thC);
            double r1e=r0e+Fdz;
            double r2e=r1e+Rdz;

            double y0e=r0e*tan(dThetaBarrel/2.);
            double y1e=r1e*tan(dThetaBarrel/2.);
            double y2e=r2e*tan(dThetaBarrel/2.);

            // Make crystal shapes
            // Projective trapezoids using EightPointSolids. see: https://root.cern.ch/doc/master/classTGeoArb8.html

            double x0y0 = (r0e*cos(thC) +y0e*sin(thC)) *tan(thC -dThetaBarrel/2.);
            double x1y0 = (r0e*cos(thC) -y0e*sin(thC)) *tan(thC +dThetaBarrel/2.);

            double x0y0l = x0y0*tan(gamma-dPhiBarrelCrystal/2);
            double x0y0r = x0y0*tan(gamma+dPhiBarrelCrystal/2);

            double x1y0l = x1y0*tan(gamma-dPhiBarrelCrystal/2);
            double x1y0r = x1y0*tan(gamma+dPhiBarrelCrystal/2);

            double x0y1 = (r1e*cos(thC) +y1e*sin(thC)) *tan(thC -dThetaBarrel/2.);
            double x1y1 = (r1e*cos(thC) -y1e*sin(thC)) *tan(thC +dThetaBarrel/2.);

            double x0y1l = x0y1*tan(gamma-dPhiBarrelCrystal/2);
            double x0y1r = x0y1*tan(gamma+dPhiBarrelCrystal/2);

            double x1y1l = x1y1*tan(gamma-dPhiBarrelCrystal/2);
            double x1y1r = x1y1*tan(gamma+dPhiBarrelCrystal/2);

            double x0y2 = (r2e*cos(thC) +y2e*sin(thC)) *tan(thC -dThetaBarrel/2.);
            double x1y2 = (r2e*cos(thC) -y2e*sin(thC)) *tan(thC +dThetaBarrel/2.);

            double x0y2l = x0y2*tan(gamma-dPhiBarrelCrystal/2);
            double x0y2r = x0y2*tan(gamma+dPhiBarrelCrystal/2);

            double x1y2l = x1y2*tan(gamma-dPhiBarrelCrystal/2);
            double x1y2r = x1y2*tan(gamma+dPhiBarrelCrystal/2);

            // Front (F) and Rear (R) crystal shapes

            double verticesF[]={
                               (x0y0r), (0+y0e),
                               (x1y0r), (0-y0e),
                               ( x1y0l), (0-y0e),
                               ( x0y0l), (0+y0e),

                                (x0y1r), (0+y1e),
                                (x1y1r), (0-y1e),
                                ( x1y1l), (0-y1e),
                                ( x0y1l), (0+y1e)

                                };

            double verticesR[]={
                                (x0y1r), (0+y1e),
                                (x1y1r), (0-y1e),
                                ( x1y1l), (0-y1e),
                                ( x0y1l), (0+y1e),

                                (x0y2r), (0+y2e),
                                (x1y2r), (0-y2e),
                                ( x1y2l), (0-y2e),
                                ( x0y2l), (0+y2e)

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

            auto crystalFId64=segmentation->setVolumeID(1, nThetaEndcap+iTheta , iPhi*nPhiBarrelCrystal+nGamma, 1);
            auto crystalRId64=segmentation->setVolumeID(1, nThetaEndcap+iTheta , iPhi*nPhiBarrelCrystal+nGamma, 2);
            int crystalFId32=segmentation->getFirst32bits(crystalFId64);
            int crystalRId32=segmentation->getFirst32bits(crystalRId64);

            // Place volumes and ID them
            dd4hep::PlacedVolume crystalFp = barrelSliceAssemblyVolume.placeVolume( crystalFVol, crystalFId32, Transform3D(rot,dispF) );
            dd4hep::PlacedVolume crystalRp = barrelSliceAssemblyVolume.placeVolume( crystalRVol, crystalRId32, Transform3D(rot,dispR) );

            crystalFp.addPhysVolID("eta", nThetaEndcap+iTheta);
            crystalFp.addPhysVolID("phi", iPhi*nPhiBarrelCrystal+nGamma);
            crystalFp.addPhysVolID("depth", 1);
            crystalFp.addPhysVolID("system", 1);
            
            crystalRp.addPhysVolID("eta", nThetaEndcap+iTheta);
            crystalRp.addPhysVolID("phi", iPhi*nPhiBarrelCrystal+nGamma);
            crystalRp.addPhysVolID("depth", 2);
            crystalRp.addPhysVolID("system", 1);

            // std::bitset<10> _eta(nThetaEndcap+iTheta);
            // std::bitset<10> _phi(iPhi*nPhiBarrelCrystal+nGamma);
            // std::bitset<3>  depthF(1);
            // std::bitset<3>  depthR(2);
            // std::bitset<32> id32F(crystalFId32);
            // std::bitset<32> id32R(crystalRId32);

            //VolIDs crystalFpVID = static_cast<VolIDs> crystalFp.VolIDs();
            //VolIDs crystalRpVID = static_cast<VolIDs> crystalRp.VolIDs();

            //std::cout << "B crystalF eta: " << nThetaEndcap+iTheta << " phi: " << iPhi << " depth: " << " 1 " << std::endl;
            //std::cout << "B crystalF eta: " << _eta << " phi: " << _phi << " depth: " << depthF << std::endl;
            //std::cout << "B crystalF eta: " << segmentation->Eta(crystalFId64) << " phi: " << segmentation->Phi(crystalFId64) << " depth: " << depthF << std::endl;
            //std::cout << "B crystalFId32: " << id32F << std::endl;
            //std::cout << "B crystalF copyNum: " <<  crystalFId32 << std::endl;
            //std::cout << "B crystalF copyNum: " <<  crystalFId64 << std::endl;
            //std::cout << "B crystalR copyNum: " << crystalRp.copyNumber() << " VolIDs: " << dd4hep::detail::tools::toString( crystalRpVID ) << std::endl;

            //std::cout << "B crystalR eta: " << nThetaEndcap+iTheta << " phi: " << iPhi << " depth: " << " 2 " << std::endl;
            //std::cout << "B crystalR eta: " << _eta << " phi: " << _phi << " depth: " << depthR << std::endl;
            //std::cout << "B crystalRId32: " << id32R << std::endl;
            //std::cout << "B crystalR copyNum: " <<  crystalRId32 << std::endl;
            //std::cout << "B crystalR copyNum: " << crystalRp.copyNumber() << " VolIDs: " << dd4hep::detail::tools::toString( crystalRpVID ) << std::endl;
          }
        }
      }

      /** 
       * 
       * Endcap
       * 
       * **/

      for (int iTheta=5; iTheta<nThetaEndcap; iTheta++) {

        double thC        = dThetaEndcap/2+ iTheta*dThetaEndcap;
        double RinEndcap  = EBz*tan(thC);

        int    nPhiEndcap = nPhiBarrel;
        double dPhiEndcap = dPhiBarrel;

        int    nPhiEndcapCrystal = floor(2*M_PI*RinEndcap/(nPhiEndcap*nomfw));
        double dPhiEndcapCrystal = dPhiEndcap/nPhiEndcapCrystal;

        double r0e=RinEndcap/sin(thC);
        double r1e=r0e+Fdz;
        double r2e=r1e+Rdz;

        double y0e=r0e*tan(dThetaEndcap/2.);
        double y1e=r1e*tan(dThetaEndcap/2.);
        double y2e=r2e*tan(dThetaEndcap/2.);

        // Skip crystals at low theta (near the beampipe) if the crystal face aspect ratio is outside 14% of unity
        // double centerCrystalWidth = RinEndcap*sin(dPhiEndcapCrystal/2);
        // if (abs(1-y0e/centerCrystalWidth)>0.14) continue;

        // Make assembly polyhedra
        double a = r0e/cos(dThetaEndcap/2);
        double z1 = a*cos(thC+dThetaEndcap/2);
        double r1min = z1*tan(thC-dThetaEndcap/2);
        double r1max = z1*tan(thC+dThetaEndcap/2);

        double b = sqrt(r2e*r2e +y2e*y2e);
        double z2 = b*cos(thC-dThetaEndcap/2);
        double r2min = z2*tan(thC-dThetaEndcap/2);
        double r2max = z2*tan(thC+dThetaEndcap/2);

        std::vector<double> zPolyhedra = {z1,z2};
        std::vector<double> rminPolyhedra = {r1min, r2min};
        std::vector<double> rmaxPolyhedra = {r1max, r2max};

        dd4hep::Polyhedra phiRingAssemblyShape(nPhiEndcap, dPhiEndcap/2, 2*M_PI, zPolyhedra, rminPolyhedra, rmaxPolyhedra);

        Position dispCone(0,0,0);
        Position dispCone1(0,0,0);
        RotationZYX rotMirror(0, 0, M_PI);

        // Endcap assembly volume
        dd4hep::Volume phiRingAssemblyVolume("EndcapRingAssembly", phiRingAssemblyShape, theDetector.material("Vacuum"));
        phiRingAssemblyVolume.setVisAttributes(theDetector, scepcalAssemblyXML.visStr());
        scepcalAssemblyVol.placeVolume( phiRingAssemblyVolume, dispCone );

        // Reflected endcap
        dd4hep::Volume phiRingAssemblyVolume1("EndcapRingAssembly1", phiRingAssemblyShape, theDetector.material("Vacuum"));
        phiRingAssemblyVolume1.setVisAttributes(theDetector, scepcalAssemblyXML.visStr());
        scepcalAssemblyVol.placeVolume( phiRingAssemblyVolume1, Transform3D(rotMirror, dispCone1) );

        for (int iPhi=0; iPhi<nPhiEndcap; iPhi++) {
          
          for (int nGamma=0; nGamma<nPhiEndcapCrystal; nGamma++) {

            double gamma = -dPhiEndcap/2+dPhiEndcapCrystal/2+dPhiEndcapCrystal*nGamma;
            // Make crystal shapes

            double x0y0 = (r0e*cos(thC) +y0e*sin(thC)) *tan(thC -dThetaEndcap/2.);
            double x1y0 = (r0e*cos(thC) -y0e*sin(thC)) *tan(thC +dThetaEndcap/2.);

            double x0y0l = x0y0*tan(gamma-dPhiEndcapCrystal/2);
            double x0y0r = x0y0*tan(gamma+dPhiEndcapCrystal/2);

            double x1y0l = x1y0*tan(gamma-dPhiEndcapCrystal/2);
            double x1y0r = x1y0*tan(gamma+dPhiEndcapCrystal/2);

            double x0y1 = (r1e*cos(thC) +y1e*sin(thC)) *tan(thC -dThetaEndcap/2.);
            double x1y1 = (r1e*cos(thC) -y1e*sin(thC)) *tan(thC +dThetaEndcap/2.);

            double x0y1l = x0y1*tan(gamma-dPhiEndcapCrystal/2);
            double x0y1r = x0y1*tan(gamma+dPhiEndcapCrystal/2);

            double x1y1l = x1y1*tan(gamma-dPhiEndcapCrystal/2);
            double x1y1r = x1y1*tan(gamma+dPhiEndcapCrystal/2);

            double x0y2 = (r2e*cos(thC) +y2e*sin(thC)) *tan(thC -dThetaEndcap/2.);
            double x1y2 = (r2e*cos(thC) -y2e*sin(thC)) *tan(thC +dThetaEndcap/2.);

            double x0y2l = x0y2*tan(gamma-dPhiEndcapCrystal/2);
            double x0y2r = x0y2*tan(gamma+dPhiEndcapCrystal/2);

            double x1y2l = x1y2*tan(gamma-dPhiEndcapCrystal/2);
            double x1y2r = x1y2*tan(gamma+dPhiEndcapCrystal/2);

            double verticesF[]={
                               (x0y0r), (0+y0e),
                               (x1y0r), (0-y0e),
                               ( x1y0l), (0-y0e),
                               ( x0y0l), (0+y0e),

                                (x0y1r), (0+y1e),
                                (x1y1r), (0-y1e),
                                ( x1y1l), (0-y1e),
                                ( x0y1l), (0+y1e)

                                };

            double verticesR[]={
                                (x0y1r), (0+y1e),
                                (x1y1r), (0-y1e),
                                ( x1y1l), (0-y1e),
                                ( x0y1l), (0+y1e),

                                (x0y2r), (0+y2e),
                                (x1y2r), (0-y2e),
                                ( x1y2l), (0-y2e),
                                ( x0y2l), (0+y2e)

                                };

            dd4hep::EightPointSolid crystalFShape(Fdz/2, verticesF);
            dd4hep::EightPointSolid crystalRShape(Rdz/2, verticesR);

            dd4hep::Volume crystalFVol("EndcapCrystalF", crystalFShape, crystalFMat);
            dd4hep::Volume crystalRVol("EndcapCrystalR", crystalRShape, crystalRMat);

            crystalFVol.setVisAttributes(theDetector, crystalFXML.visStr());
            crystalRVol.setVisAttributes(theDetector, crystalRXML.visStr());
            
            double phi=iPhi*dPhiEndcap;

            RotationZYX rot(M_PI/2, thC, 0);
            ROOT::Math::RotationZ rotZ = ROOT::Math::RotationZ(phi);
            rot = rotZ*rot;

            double rF=r0e+Fdz/2.;
            Position dispF(rF*sin(thC)*cos(phi),
                          rF*sin(thC)*sin(phi),
                          rF*cos(thC));

            double rR=r1e+Rdz/2.;
            Position dispR(rR*sin(thC)*cos(phi),
                          rR*sin(thC)*sin(phi),
                          rR*cos(thC));

            auto crystalFId64=segmentation->setVolumeID(1, iTheta, iPhi*nPhiEndcapCrystal+nGamma, 1);
            auto crystalRId64=segmentation->setVolumeID(1, iTheta, iPhi*nPhiEndcapCrystal+nGamma, 2);
            int crystalFId32=segmentation->getFirst32bits(crystalFId64);
            int crystalRId32=segmentation->getFirst32bits(crystalRId64);

            dd4hep::PlacedVolume crystalFp = phiRingAssemblyVolume.placeVolume( crystalFVol, crystalFId32, Transform3D(rot,dispF-dispCone) );
            dd4hep::PlacedVolume crystalRp = phiRingAssemblyVolume.placeVolume( crystalRVol, crystalRId32, Transform3D(rot,dispR-dispCone) );

            crystalFp.addPhysVolID("eta", iTheta);
            crystalFp.addPhysVolID("phi", iPhi*nPhiEndcapCrystal+nGamma);
            crystalFp.addPhysVolID("depth", 1);
            crystalFp.addPhysVolID("system", 1);

            crystalRp.addPhysVolID("eta", iTheta);
            crystalRp.addPhysVolID("phi", iPhi*nPhiEndcapCrystal+nGamma);
            crystalRp.addPhysVolID("depth", 2);
            crystalRp.addPhysVolID("system", 1);

            // Add suffix 1 for the other endcap (mirrored)

            auto crystalFId641=segmentation->setVolumeID(1, nThetaEndcap+nThetaBarrel+nThetaEndcap-iTheta, iPhi*nPhiEndcapCrystal+nGamma, 1);
            auto crystalRId641=segmentation->setVolumeID(1, nThetaEndcap+nThetaBarrel+nThetaEndcap-iTheta, iPhi*nPhiEndcapCrystal+nGamma, 2);
            int crystalFId321=segmentation->getFirst32bits(crystalFId641);
            int crystalRId321=segmentation->getFirst32bits(crystalRId641);

            dd4hep::PlacedVolume crystalFp1 = phiRingAssemblyVolume1.placeVolume( crystalFVol, crystalFId321, Transform3D(rot,dispF-dispCone) );
            dd4hep::PlacedVolume crystalRp1 = phiRingAssemblyVolume1.placeVolume( crystalRVol, crystalRId321, Transform3D(rot,dispR-dispCone) );

            crystalFp1.addPhysVolID("eta", nThetaEndcap+nThetaBarrel+nThetaEndcap-iTheta);
            crystalFp1.addPhysVolID("phi", iPhi*nPhiEndcapCrystal+nGamma);
            crystalFp1.addPhysVolID("depth", 1);
            crystalFp1.addPhysVolID("system", 1);

            crystalRp1.addPhysVolID("eta", nThetaEndcap+nThetaBarrel+nThetaEndcap-iTheta);
            crystalRp1.addPhysVolID("phi", iPhi*nPhiEndcapCrystal+nGamma);
            crystalRp1.addPhysVolID("depth", 2);
            crystalRp1.addPhysVolID("system", 1);


            // std::bitset<10> _eta(iTheta);
            // std::bitset<10> _eta1(iTheta);
            // std::bitset<10> _phi(iPhi*nPhiEndcapCrystal+nGamma);
            // std::bitset<3>  depthF(1);
            // std::bitset<3>  depthR(2);
            // std::bitset<32> id32F(crystalFId32);
            // std::bitset<32> id32R(crystalRId32);

            // std::bitset<32> id32F1(crystalFId321);
            // std::bitset<32> id32R1(crystalRId321);

            //std::cout << "E crystalF eta: " << iTheta << " phi: " << iPhi << " depth: " << 1 << std::endl;
            //std::cout << "E crystalF eta: " << _eta << " phi: " << _phi << " depth: " << depthF << std::endl;
            //std::cout << "E crystalFId32: " << id32F << std::endl;
            //std::cout << "E crystalF copyNum: " << crystalFp.copyNumber() << " VolIDs: " << dd4hep::detail::tools::toString( crystalFp.VolIDs()) << std::endl;

            //std::cout << "E crystalR eta: " << iTheta << " phi: " << iPhi << " depth: " << 2 << std::endl;
            //std::cout << "E crystalR eta: " << _eta << " phi: " << _phi << " depth: " << depthR << std::endl;
            //std::cout << "E crystalRId32: " << id32R << std::endl;
            //std::cout << "E crystalR copyNum: " << crystalRp.copyNumber() << " VolIDs: " << dd4hep::detail::tools::toString( crystalRp.VolIDs()) << std::endl;

            //std::cout << "E crystalF1 eta: " << iTheta << " phi: " << iPhi << " depth: " << 1 << std::endl;
            //std::cout << "E crystalF1 eta: " << _eta1 << " phi: " << _phi << " depth: " << depthF << std::endl;
            //std::cout << "E crystalFId321: " << id32F1 << std::endl;
            //std::cout << "E crystalF1 copyNum: " << crystalFp1.copyNumber() << " VolIDs: " << dd4hep::detail::tools::toString(crystalFp1.VolIDs()) << std::endl;

            //std::cout << "E crystalR1 eta: " << iTheta << " phi: " << iPhi << " depth: " << 2 << std::endl;
            //std::cout << "E crystalR1 eta: " << _eta1 << " phi: " << _phi << " depth: " << depthR << std::endl;
            //std::cout << "E crystalRId321: " << id32R1 << std::endl;
            //std::cout << "E crystalR1 copyNum: " << crystalRp1.copyNumber() << " VolIDs: " << dd4hep::detail::tools::toString(crystalRp1.VolIDs()) << std::endl;
          }
        }
      }

      // Place the detector
      experimentalHall.placeVolume( scepcalAssemblyVol );

      dd4hep::PlacedVolume hallPlace=theDetector.pickMotherVolume(drDet).placeVolume(experimentalHall);

      //hallPlace.addPhysVolID("system", detectorXML.id());
      hallPlace.addPhysVolID("system", 0);

      drDet.setPlacement(hallPlace);

      return drDet;
    }
}


DECLARE_DETELEMENT(SCEPCAL, ddSCEPCAL::create_detector)
