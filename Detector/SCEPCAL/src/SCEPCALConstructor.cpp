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
      xml_comp_t timingXML=detectorXML.child(_Unicode(timingLayer));
      xml_comp_t instXML=detectorXML.child(_Unicode(inst));
      xml_comp_t towerAssemblyXML=detectorXML.child(_Unicode(towerAssembly));

      //-----------------------------------------------------------------------------------

      sensDet.setType(sensDetType.typeStr());

      // Set the segmentation class
      auto segmentation=dynamic_cast<dd4hep::DDSegmentation::SCEPCALSegmentation *>( sensDet.readout().segmentation().segmentation());

      dd4hep::Assembly experimentalHall("hall");

      // Parse input parameters from imported xml objects
      const double Fdz    =crystalFXML.attr<double>(_Unicode(length));
      const double Rdz    =crystalRXML.attr<double>(_Unicode(length));
      const double nomfw  =dimXML.attr<double>(_Unicode(crystalFaceWidthNominal));
      const double EBz    =dimXML.attr<double>(_Unicode(barrelHalfZ));
      const double Rin    =dimXML.attr<double>(_Unicode(barrelInnerR));
      const double cube   =towerAssemblyXML.attr<double>(_Unicode(cube));

      // Material definitions
      dd4hep::Material crystalFMat =theDetector.material(crystalFXML.materialStr());
      dd4hep::Material crystalRMat =theDetector.material(crystalRXML.materialStr());
      dd4hep::Material timingMat   =theDetector.material(timingXML.materialStr());
      dd4hep::Material instMat     =theDetector.material(instXML.materialStr());

      // Begin geometry calculations
      int nThetaBarrel  =floor(EBz/nomfw);
      int nThetaEndcap  =floor(Rin/nomfw);

      double thetaSizeEndcap=atan(Rin/EBz);

      double dThetaBarrel =(M_PI/2-thetaSizeEndcap)/(nThetaBarrel);
      double dThetaEndcap =thetaSizeEndcap/nThetaEndcap;

      int    nPhiBarrel   =floor(2*M_PI*Rin/nomfw);
      double dPhiBarrel   =2*M_PI/nPhiBarrel;

      // Make towers along theta (eta) and make rotations in phi each step (i.e. phi nested in theta)
      for (int iTheta=0; iTheta<2*nThetaBarrel+1; iTheta++) {
        if (iTheta == nThetaBarrel) continue;
        double thC =thetaSizeEndcap+(iTheta*dThetaBarrel);

        // Projective towers using EightPointSolids. see: https://root.cern.ch/doc/master/classTGeoArb8.html

        double r0 =Rin/sin(thC);
        double r1 =r0+Fdz;
        double r2 =r1+Rdz;

        double y0 =r0*tan(dThetaBarrel/2.);
        double y1 =r1*tan(dThetaBarrel/2.);
        double y2 =r2*tan(dThetaBarrel/2.);

        double x0y0 = (r0*cos(thC) +y0*sin(thC)) *tan(thC -dThetaBarrel/2.) *tan(dPhiBarrel/2.);
        double x1y0 = (r0*cos(thC) -y0*sin(thC)) *tan(thC +dThetaBarrel/2.) *tan(dPhiBarrel/2.);

        double x0y1 = (r1*cos(thC) +y1*sin(thC)) *tan(thC -dThetaBarrel/2.) *tan(dPhiBarrel/2.);
        double x1y1 = (r1*cos(thC) -y1*sin(thC)) *tan(thC +dThetaBarrel/2.) *tan(dPhiBarrel/2.);

        double x0y2 = (r2*cos(thC) +y2*sin(thC)) *tan(thC -dThetaBarrel/2.) *tan(dPhiBarrel/2.);
        double x1y2 = (r2*cos(thC) -y2*sin(thC)) *tan(thC +dThetaBarrel/2.) *tan(dPhiBarrel/2.);

        // Tower Assembly shape (A) holds Forward (F) and Rear (R) crystal shapes

        double verticesF[]={x0y0,y0,x1y0,-y0,-x1y0,-y0,-x0y0,y0,
                            x0y1,y1,x1y1,-y1,-x1y1,-y1,-x0y1,y1};

        double verticesR[]={x0y1,y1,x1y1,-y1,-x1y1,-y1,-x0y1,y1,
                            x0y2,y2,x1y2,-y2,-x1y2,-y2,-x0y2,y2};

        double verticesA[]={x0y0,y0,x1y0,-y0,-x1y0,-y0,-x0y0,y0,
                            x0y2,y2,x1y2,-y2,-x1y2,-y2,-x0y2,y2};

        dd4hep::EightPointSolid crystalFShape(Fdz/2, verticesF);
        dd4hep::EightPointSolid crystalRShape(Rdz/2, verticesR);
        dd4hep::EightPointSolid towerAssemblyShapeBarrel((Fdz+Rdz)/2, verticesA);

        // Promote shapes to volumes and set attributes
        dd4hep::Volume crystalFVol("BarrelCrystalR", crystalFShape, crystalFMat);
        dd4hep::Volume crystalRVol("BarrelCrystalF", crystalRShape, crystalRMat);

        crystalFVol.setVisAttributes(theDetector, crystalFXML.visStr());
        crystalRVol.setVisAttributes(theDetector, crystalRXML.visStr());

        // Rotations in phi to place volumes
        for (int iPhi=0; iPhi<nPhiBarrel; iPhi++) {

          auto crystalFId64=segmentation->setVolumeID(1, (nThetaBarrel-iTheta) , iPhi, 1);
          auto crystalRId64=segmentation->setVolumeID(1, (nThetaBarrel-iTheta) , iPhi, 2);

          int crystalFId32=segmentation->getFirst32bits(crystalFId64);
          int crystalRId32=segmentation->getFirst32bits(crystalRId64);

          dd4hep::Volume towerAssemblyVol("towerAssemblyVol", towerAssemblyShapeBarrel, theDetector.material("Vacuum"));
          towerAssemblyVol.setVisAttributes(theDetector, towerAssemblyXML.visStr());

          double rt=r0+(Fdz+Rdz)/2.;
          double phi=iPhi*dPhiBarrel;

          // Have to use the ROOT rotation class here because Euler rotations are not implemented in dd4hep
          RotationZYX rot(M_PI/2, thC, 0);
          ROOT::Math::RotationZ rotZ = ROOT::Math::RotationZ(phi);
          rot = rotZ*rot;

          Position disp(rt*sin(thC)*cos(phi),
                        rt*sin(thC)*sin(phi),
                        rt*cos(thC));

          experimentalHall.placeVolume( towerAssemblyVol, Transform3D(rot,disp) );

          // Place volumes and ID them
          dd4hep::PlacedVolume crystalFp = towerAssemblyVol.placeVolume( crystalFVol, crystalFId32, Position(0,0,-Rdz/2) );
          dd4hep::PlacedVolume crystalRp = towerAssemblyVol.placeVolume( crystalRVol, crystalRId32, Position(0,0,Fdz/2) );

          crystalFp.addPhysVolID("eta", nThetaBarrel-iTheta);
          crystalFp.addPhysVolID("phi", iPhi);
          crystalFp.addPhysVolID("depth", 1);
          crystalFp.addPhysVolID("system", 1);
          
          crystalRp.addPhysVolID("eta", nThetaBarrel-iTheta );
          crystalRp.addPhysVolID("phi", iPhi);
          crystalRp.addPhysVolID("depth", 2);
          crystalRp.addPhysVolID("system", 1);

          std::bitset<10> _eta((nThetaBarrel-iTheta) );
          std::bitset<10> _phi(iPhi);
          std::bitset<3> depthF(1);
          std::bitset<3> depthR(2);
          std::bitset<32> id32F(crystalFId32);
          std::bitset<32> id32R(crystalRId32);

          //VolIDs crystalFpVID = static_cast<VolIDs> crystalFp.VolIDs();
          //VolIDs crystalRpVID = static_cast<VolIDs> crystalRp.VolIDs();

          //std::cout << "B crystalF eta: " << ((nThetaBarrel-iTheta) ) << " phi: " << iPhi << " depth: " << " 1 " << std::endl;
          //std::cout << "B crystalF eta: " << _eta << " phi: " << _phi << " depth: " << depthF << std::endl;
          //std::cout << "B crystalF eta: " << segmentation->Eta(crystalFId64) << " phi: " << segmentation->Phi(crystalFId64) << " depth: " << depthF << std::endl;
          //std::cout << "B crystalFId32: " << id32F << std::endl;
          //std::cout << "B crystalF copyNum: " <<  crystalFId32 << std::endl;
          //std::cout << "B crystalF copyNum: " <<  crystalFId64 << std::endl;
          //std::cout << "B crystalR copyNum: " << crystalRp.copyNumber() << " VolIDs: " << dd4hep::detail::tools::toString( crystalRpVID ) << std::endl;

          //std::cout << "B crystalR eta: " << ((nThetaBarrel-iTheta) ) << " phi: " << iPhi << " depth: " << " 2 " << std::endl;
          //std::cout << "B crystalR eta: " << _eta << " phi: " << _phi << " depth: " << depthR << std::endl;
          //std::cout << "B crystalRId32: " << id32R << std::endl;
          //std::cout << "B crystalR copyNum: " <<  crystalRId32 << std::endl;
          //std::cout << "B crystalR copyNum: " << crystalRp.copyNumber() << " VolIDs: " << dd4hep::detail::tools::toString( crystalRpVID ) << std::endl;
        }
      }

      // Endcap
      for (int iTheta=1; iTheta<nThetaEndcap; iTheta++) {

        double thC        =iTheta*dThetaEndcap;
        double RinEndcap  = EBz*tan(thC);

        // Same calculations, except nPhiEndcap changes for each theta instead of being constant
        int    nPhiEndcap =floor(2*M_PI*RinEndcap/nomfw);
        double dPhiEndcap =2*M_PI/nPhiEndcap;

        double r0=RinEndcap/sin(thC);
        double r1=r0+Fdz;
        double r2=r1+Rdz;

        double y0=r0*tan(dThetaEndcap/2.);
        double y1=r1*tan(dThetaEndcap/2.);
        double y2=r2*tan(dThetaEndcap/2.);

        // Skip crystals at low theta (near the beampipe) if the crystal face aspect ratio is over 15% of unity
        double centralHalfWidthActual = RinEndcap*sin(dPhiEndcap/2);
        //std::cout << "For theta " << nThetaBarrel+nThetaEndcap-iTheta <<  ", crystal face aspect ratio " << abs(1-y0/centralHalfWidthActual) << std::endl;
        if (abs(1-y0/centralHalfWidthActual)>0.14) { // changed to 14, since aspect ratio of first tower (42/-42) is 14.9
          //std::cout << " More than 14% --> SKIPPING ... " << std::endl;
          continue;}

        double x0y0 = (r0*cos(thC) +y0*sin(thC)) *tan(thC -dThetaEndcap/2.) *tan(dPhiEndcap/2.);
        double x1y0 = (r0*cos(thC) -y0*sin(thC)) *tan(thC +dThetaEndcap/2.) *tan(dPhiEndcap/2.);

        double x0y1 = (r1*cos(thC) +y1*sin(thC)) *tan(thC -dThetaEndcap/2.) *tan(dPhiEndcap/2.);
        double x1y1 = (r1*cos(thC) -y1*sin(thC)) *tan(thC +dThetaEndcap/2.) *tan(dPhiEndcap/2.);

        double x0y2 = (r2*cos(thC) +y2*sin(thC)) *tan(thC -dThetaEndcap/2.) *tan(dPhiEndcap/2.);
        double x1y2 = (r2*cos(thC) -y2*sin(thC)) *tan(thC +dThetaEndcap/2.) *tan(dPhiEndcap/2.);

        double verticesF[]={x0y0,y0,x1y0,-y0,-x1y0,-y0,-x0y0,y0,
                            x0y1,y1,x1y1,-y1,-x1y1,-y1,-x0y1,y1};

        double verticesR[]={x0y1,y1,x1y1,-y1,-x1y1,-y1,-x0y1,y1,
                            x0y2,y2,x1y2,-y2,-x1y2,-y2,-x0y2,y2};

        double verticesA[]={x0y0,y0,x1y0,-y0,-x1y0,-y0,-x0y0,y0,
                            x0y2,y2,x1y2,-y2,-x1y2,-y2,-x0y2,y2};

        dd4hep::EightPointSolid crystalFShape(Fdz/2, verticesF);
        dd4hep::EightPointSolid crystaslRShape(Rdz/2, verticesR);
        dd4hep::EightPointSolid towerAssemblyShapeEndcap((Fdz+Rdz)/2, verticesA);

        dd4hep::Volume crystalFVol("EndcapCrystalR", crystalFShape, crystalFMat);
        dd4hep::Volume crystalRVol("EndcapCrystalF", crystaslRShape, crystalRMat);

        crystalFVol.setVisAttributes(theDetector, crystalFXML.visStr());
        crystalRVol.setVisAttributes(theDetector, crystalRXML.visStr());

        dd4hep::Box towerAssemblyBoxEndcap(cube,cube,cube);

        for (int iPhi=0; iPhi<nPhiEndcap; iPhi++) {

          auto crystalFId64=segmentation->setVolumeID(1,(nThetaBarrel+nThetaEndcap-iTheta), iPhi, 1);
          auto crystalRId64=segmentation->setVolumeID(1,(nThetaBarrel+nThetaEndcap-iTheta), iPhi, 2);
          auto crystalFId641=segmentation->setVolumeID(1, -(nThetaBarrel+nThetaEndcap-iTheta), iPhi, 1);
          auto crystalRId641=segmentation->setVolumeID(1, -(nThetaBarrel+nThetaEndcap-iTheta), iPhi, 2);

          int crystalFId32=segmentation->getFirst32bits(crystalFId64);
          int crystalRId32=segmentation->getFirst32bits(crystalRId64);
          int crystalFId321=segmentation->getFirst32bits(crystalFId641);
          int crystalRId321=segmentation->getFirst32bits(crystalRId641);

          dd4hep::Volume towerAssemblyVolEndcap("towerAssemblyVolEndcap", towerAssemblyShapeEndcap, theDetector.material("Vacuum"));
          dd4hep::Volume towerAssemblyVolEndcap1("towerAssemblyVolEndcap1", towerAssemblyShapeEndcap, theDetector.material("Vacuum"));
          towerAssemblyVolEndcap.setVisAttributes(theDetector, towerAssemblyXML.visStr());
          towerAssemblyVolEndcap1.setVisAttributes(theDetector, towerAssemblyXML.visStr());

          double rt=r0+(Fdz+Rdz)/2.;
          double phi=iPhi*dPhiEndcap;

          ROOT::Math::RotationZ rotZ = ROOT::Math::RotationZ(phi);

          RotationZYX rot(M_PI/2, thC, 0);
          rot = rotZ*rot;
          Position disp(rt*sin(thC)*cos(phi),
                        rt*sin(thC)*sin(phi),
                        rt*cos(thC));

          RotationZYX rot1(M_PI/2, thC, M_PI);
          rot1 = rotZ*rot1;
          Position disp1(rt*sin(thC)*cos(phi),
                        rt*sin(thC)*sin(phi),
                        -rt*cos(thC));

          experimentalHall.placeVolume( towerAssemblyVolEndcap, Transform3D(rot,disp) );
          experimentalHall.placeVolume( towerAssemblyVolEndcap1, Transform3D(rot1,disp1) );

          dd4hep::PlacedVolume crystalFp = towerAssemblyVolEndcap.placeVolume( crystalFVol, crystalFId32, Position(0,0,-Rdz/2) );
          dd4hep::PlacedVolume crystalRp = towerAssemblyVolEndcap.placeVolume( crystalRVol, crystalRId32, Position(0,0,Fdz/2) );

          dd4hep::PlacedVolume crystalFp1 = towerAssemblyVolEndcap1.placeVolume( crystalFVol, crystalFId321, Position(0,0,-Rdz/2) );
          dd4hep::PlacedVolume crystalRp1 = towerAssemblyVolEndcap1.placeVolume( crystalRVol, crystalRId321, Position(0,0,Fdz/2) );

          crystalFp.addPhysVolID("eta", (nThetaBarrel+nThetaEndcap-iTheta));
          crystalFp.addPhysVolID("phi", iPhi);
          crystalFp.addPhysVolID("depth", 1);
          crystalFp.addPhysVolID("system", 1);

          crystalRp.addPhysVolID("eta", (nThetaBarrel+nThetaEndcap-iTheta));
          crystalRp.addPhysVolID("phi", iPhi);
          crystalRp.addPhysVolID("depth", 2);
          crystalRp.addPhysVolID("system", 1);


          crystalFp1.addPhysVolID("eta", -(nThetaBarrel+nThetaEndcap-iTheta));
          crystalFp1.addPhysVolID("phi", iPhi);
          crystalFp1.addPhysVolID("depth", 1);
          crystalFp1.addPhysVolID("system", 1);

          crystalRp1.addPhysVolID("eta", -(nThetaBarrel+nThetaEndcap-iTheta));
          crystalRp1.addPhysVolID("phi", iPhi);
          crystalRp1.addPhysVolID("depth", 2);
          crystalRp1.addPhysVolID("system", 1);

          std::bitset<10> _eta((nThetaBarrel+nThetaEndcap-iTheta));
          std::bitset<10> _eta1(-(nThetaBarrel+nThetaEndcap-iTheta));
          std::bitset<10> _phi(iPhi);
          std::bitset<3> depthF(1);
          std::bitset<3> depthR(2);
          std::bitset<32> id32F(crystalFId32);
          std::bitset<32> id32F1(crystalFId321);
          std::bitset<32> id32R(crystalRId32);
          std::bitset<32> id32R1(crystalRId321);



          //std::cout << "E crystalF eta: " << (nThetaBarrel+nThetaEndcap-iTheta) << " phi: " << iPhi << " depth: " << 1 << std::endl;
          //std::cout << "E crystalF eta: " << _eta << " phi: " << _phi << " depth: " << depthF << std::endl;
          //std::cout << "E crystalFId32: " << id32F << std::endl;
          //std::cout << "E crystalF copyNum: " << crystalFp.copyNumber() << " VolIDs: " << dd4hep::detail::tools::toString( crystalFp.VolIDs()) << std::endl;

          //std::cout << "E crystalR eta: " << (nThetaBarrel+nThetaEndcap-iTheta) << " phi: " << iPhi << " depth: " << 2 << std::endl;
          //std::cout << "E crystalR eta: " << _eta << " phi: " << _phi << " depth: " << depthR << std::endl;
          //std::cout << "E crystalRId32: " << id32R << std::endl;
          //std::cout << "E crystalR copyNum: " << crystalRp.copyNumber() << " VolIDs: " << dd4hep::detail::tools::toString( crystalRp.VolIDs()) << std::endl;

          //std::cout << "E crystalF1 eta: " << -(nThetaBarrel+nThetaEndcap-iTheta) << " phi: " << iPhi << " depth: " << 1 << std::endl;
          //std::cout << "E crystalF1 eta: " << _eta1 << " phi: " << _phi << " depth: " << depthF << std::endl;
          //std::cout << "E crystalFId321: " << id32F1 << std::endl;
          //std::cout << "E crystalF1 copyNum: " << crystalFp1.copyNumber() << " VolIDs: " << dd4hep::detail::tools::toString(crystalFp1.VolIDs()) << std::endl;

          //std::cout << "E crystalR1 eta: " << -(nThetaBarrel+nThetaEndcap-iTheta) << " phi: " << iPhi << " depth: " << 2 << std::endl;
          //std::cout << "E crystalR1 eta: " << _eta1 << " phi: " << _phi << " depth: " << depthR << std::endl;
          //std::cout << "E crystalRId321: " << id32R1 << std::endl;
          //std::cout << "E crystalR1 copyNum: " << crystalRp1.copyNumber() << " VolIDs: " << dd4hep::detail::tools::toString(crystalRp1.VolIDs()) << std::endl;
        }


      }

      // Place the detector

      dd4hep::PlacedVolume hallPlace=theDetector.pickMotherVolume(drDet).placeVolume(experimentalHall);

      //hallPlace.addPhysVolID("system", detectorXML.id());
      hallPlace.addPhysVolID("system", 0);

      drDet.setPlacement(hallPlace);

      return drDet;
    }
}


DECLARE_DETELEMENT(SCEPCAL, ddSCEPCAL::create_detector)
