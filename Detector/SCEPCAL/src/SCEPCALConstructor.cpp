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
      const double EBz    =dimXML.attr<double>(_Unicode(barrelHalfZ));
      const double Rin    =dimXML.attr<double>(_Unicode(barrelInnerR));

      // Material definitions
      dd4hep::Material crystalFMat =theDetector.material(crystalFXML.materialStr());
      dd4hep::Material crystalRMat =theDetector.material(crystalRXML.materialStr());
      dd4hep::Material timingMat   =theDetector.material(timingXML.materialStr());
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

      // Make one theta slice then rotate in phi (i.e. theta nested in phi)

      // for (int iPhi=0; iPhi<nPhiBarrel; iPhi++) {
      for (int iPhi=1; iPhi<1; iPhi++) {

        if (debugLevel>1) std::cout << "Barrel: phi: " << iPhi << std::endl;

        for (int nGamma=0; nGamma<nPhiBarrelCrystal; nGamma++) {
          // if (nGamma%2==1) continue;
          double aBarrel = Rin*cos(dPhiBarrel/2);
          double RinBn = aBarrel/ abs(cos(dPhiBarrel/2-dPhiBarrelCrystal/2-dPhiBarrelCrystal*nGamma));

          // Make assembly slice
          double thC_end = thetaSizeEndcap+dThetaBarrel/2;
          double r0slice_end =RinBn/sin(thC_end);
          double y0slice_end =r0slice_end*tan(dThetaBarrel/2.);

          double slice_front_jut = y0slice_end*sin(M_PI/2-thC_end);

          double z1slice =RinBn -slice_front_jut;
          double z2slice =RinBn +Fdz +Rdz +slice_front_jut;

          double zheight_slice = (z2slice-z1slice)/2;

          double y1slice =z1slice*tan(M_PI/2-thetaSizeEndcap);
          double y2slice =z2slice*tan(M_PI/2-thetaSizeEndcap);

          double x1slice =z1slice*tan(dPhiBarrelCrystal/2);
          double x2slice =z2slice*tan(dPhiBarrelCrystal/2);

          double verticesS[]={y1slice,x1slice,y1slice,-x1slice,-y1slice,-x1slice,-y1slice,x1slice,
                              y2slice,x2slice,y2slice,-x2slice,-y2slice,-x2slice,-y2slice,x2slice};
          
          dd4hep::EightPointSolid barrelSliceAssemblyShape(zheight_slice,verticesS);
          dd4hep::Volume barrelSliceAssemblyVolume("BarrelSliceAssembly", barrelSliceAssemblyShape, theDetector.material("Vacuum"));

          barrelSliceAssemblyVolume.setVisAttributes(theDetector, scepcalAssemblyXML.visStr());

          double phi=iPhi*dPhiBarrel+nGamma*dPhiBarrelCrystal;
          RotationZYX rotSlice(0, M_PI/2, 0);
          ROOT::Math::RotationZ rotZSlice = ROOT::Math::RotationZ(phi);
          rotSlice = rotZSlice*rotSlice;

          double rSlice =(z1slice+z2slice)/2;
          Position dispSlice(rSlice*cos(phi),rSlice*sin(phi),0);

          scepcalAssemblyVol.placeVolume( barrelSliceAssemblyVolume, Transform3D(rotSlice, dispSlice) );


          for (int iTheta=0; iTheta<nThetaBarrel; iTheta++) {
        
            if (debugLevel>1) std::cout << "  Barrel: theta: " << iTheta << std::endl;

            double thC =thetaSizeEndcap+dThetaBarrel/2+(iTheta*dThetaBarrel);

            // Projective trapezoids using EightPointSolids. see: https://root.cern.ch/doc/master/classTGeoArb8.html

            double r0 =RinBn/sin(thC);
            double r1 =r0+Fdz;
            double r2 =r1+Rdz;

            double y0 =r0*tan(dThetaBarrel/2.);
            double y1 =r1*tan(dThetaBarrel/2.);
            double y2 =r2*tan(dThetaBarrel/2.);

            double x0y0 = (r0*cos(thC) +y0*sin(thC)) *tan(thC -dThetaBarrel/2.) *tan(dPhiBarrelCrystal/2.);
            double x1y0 = (r0*cos(thC) -y0*sin(thC)) *tan(thC +dThetaBarrel/2.) *tan(dPhiBarrelCrystal/2.);

            double x0y1 = (r1*cos(thC) +y1*sin(thC)) *tan(thC -dThetaBarrel/2.) *tan(dPhiBarrelCrystal/2.);
            double x1y1 = (r1*cos(thC) -y1*sin(thC)) *tan(thC +dThetaBarrel/2.) *tan(dPhiBarrelCrystal/2.);

            double x0y2 = (r2*cos(thC) +y2*sin(thC)) *tan(thC -dThetaBarrel/2.) *tan(dPhiBarrelCrystal/2.);
            double x1y2 = (r2*cos(thC) -y2*sin(thC)) *tan(thC +dThetaBarrel/2.) *tan(dPhiBarrelCrystal/2.);

            // Front (F) and Rear (R) crystal shapes

            double verticesF[]={x0y0,y0,x1y0,-y0,-x1y0,-y0,-x0y0,y0,
                                x0y1,y1,x1y1,-y1,-x1y1,-y1,-x0y1,y1};

            double verticesR[]={x0y1,y1,x1y1,-y1,-x1y1,-y1,-x0y1,y1,
                                x0y2,y2,x1y2,-y2,-x1y2,-y2,-x0y2,y2};

            dd4hep::EightPointSolid crystalFShape(Fdz/2, verticesF);
            dd4hep::EightPointSolid crystalRShape(Rdz/2, verticesR);

            // Promote shapes to volumes and set attributes
            dd4hep::Volume crystalFVol("BarrelCrystalR", crystalFShape, crystalFMat);
            dd4hep::Volume crystalRVol("BarrelCrystalF", crystalRShape, crystalRMat);

            crystalFVol.setVisAttributes(theDetector, crystalFXML.visStr());
            crystalRVol.setVisAttributes(theDetector, crystalRXML.visStr());

            auto crystalFId64=segmentation->setVolumeID(1, nThetaEndcap+iTheta , iPhi*nPhiBarrel+nGamma, 1);
            auto crystalRId64=segmentation->setVolumeID(1, nThetaEndcap+iTheta , iPhi*nPhiBarrel+nGamma, 2);

            int crystalFId32=segmentation->getFirst32bits(crystalFId64);
            int crystalRId32=segmentation->getFirst32bits(crystalRId64);

            // Use the ROOT rotation class here because Euler rotations are not implemented in dd4hep
            // double phiCrystal=dPhiBarrel/2-nGamma*dPhiBarrelCrystal;

            RotationZYX rot(M_PI/2, -M_PI/2 +thC, 0);
            // ROOT::Math::RotationZ rotZ = ROOT::Math::RotationZ(phiCrystal);
            // rot = rotZ*rot;

            double rF=r0+Fdz/2.;
            Position dispF(-rF*cos(thC),
                          0, //rF*sin(thC)*sin(phiCrystal),
                          -(rSlice-rF*sin(thC))

                          );

            double rR=r1+Rdz/2.;
            Position dispR(-rR*cos(thC),
                          0, //rR*sin(thC)*sin(phiCrystal),
                          -(rSlice-rR*sin(thC))
                          );

            // Place volumes and ID them
            dd4hep::PlacedVolume crystalFp = barrelSliceAssemblyVolume.placeVolume( crystalFVol, crystalFId32, Transform3D(rot,dispF) );
            dd4hep::PlacedVolume crystalRp = barrelSliceAssemblyVolume.placeVolume( crystalRVol, crystalRId32, Transform3D(rot,dispR) );

            crystalFp.addPhysVolID("eta", nThetaEndcap+iTheta);
            crystalFp.addPhysVolID("phi", iPhi*nPhiBarrel+nGamma);
            crystalFp.addPhysVolID("depth", 1);
            crystalFp.addPhysVolID("system", 1);
            
            crystalRp.addPhysVolID("eta", nThetaEndcap+iTheta );
            crystalRp.addPhysVolID("phi", iPhi*nPhiBarrel+nGamma);
            crystalRp.addPhysVolID("depth", 2);
            crystalRp.addPhysVolID("system", 1);

            std::bitset<10> _eta(nThetaEndcap+iTheta);
            std::bitset<10> _phi(iPhi*nPhiBarrel+nGamma);
            std::bitset<3>  depthF(1);
            std::bitset<3>  depthR(2);
            std::bitset<32> id32F(crystalFId32);
            std::bitset<32> id32R(crystalRId32);

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

      // for (int iTheta=1; iTheta<nThetaEndcap; iTheta++) {
      for (int iTheta=nThetaEndcap-6; iTheta<nThetaEndcap-5; iTheta++) {

        if (debugLevel>1) std::cout << "Endcap: theta: " << iTheta << std::endl;

        double thC        = dThetaEndcap/2+ iTheta*dThetaEndcap;
        double RinEndcap  = EBz*tan(thC);

        // Same calculations, except nPhiEndcap changes for each theta instead of being constant
        // int    nPhiEndcap =floor(2*M_PI*RinEndcap/nomfw);
        // double dPhiEndcap =2*M_PI/nPhiEndcap;
        int nPhiEndcap = nPhiBarrel;
        double dPhiEndcap = dPhiBarrel;

        int nPhiEndcapCrystal = floor(2*M_PI*RinEndcap/(nPhiEndcap*nomfw));
        double dPhiEndcapCrystal = dPhiEndcap/nPhiEndcapCrystal;

        double r0e=RinEndcap/sin(thC);
        double r1e=r0e+Fdz;
        double r2e=r1e+Rdz;

        double y0e=r0e*tan(dThetaEndcap/2.);
        double y1e=r1e*tan(dThetaEndcap/2.);
        double y2e=r2e*tan(dThetaEndcap/2.);

        // Skip crystals at low theta (near the beampipe) if the crystal face aspect ratio is outside 15% of unity
        // double centralHalfWidthActual = RinEndcap*sin(dPhiEndcap/2);
        //std::cout << "For theta " << nThetaBarrel+nThetaEndcap-iTheta <<  ", crystal face aspect ratio " << abs(1-y0/centralHalfWidthActual) << std::endl;
        // changed to 14, since aspect ratio of first tower (42/-42) is 14.9
        // if (abs(1-y0/centralHalfWidthActual)>0.14) {
          //std::cout << " More than 14% --> SKIPPING ... " << std::endl;
          // std::cout << "skipping nPhiEndcap: " << nPhiEndcap << std::endl;
          // continue;
        // }

        // Make assembly polyhedra
        double a = r0e/cos(dThetaEndcap/2);
        double z1 = a*cos(thC+dThetaEndcap/2);
        double r1min = z1*tan(thC-dThetaEndcap/2);
        double r1max = z1*tan(thC+dThetaEndcap/2);

        double b = sqrt(r2e*r2e +y2e*y2e);
        double z2 = b*cos(thC-dThetaEndcap/2);
        double r2min = z2*tan(thC-dThetaEndcap/2);
        double r2max = z2*tan(thC+dThetaEndcap/2);

        double zcone = (z2-z1)/2;

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
        // dd4hep::Volume phiRingAssemblyVolume1("EndcapRingAssembly1", phiRingAssemblyShape, theDetector.material("Vacuum"));
        // phiRingAssemblyVolume1.setVisAttributes(theDetector, scepcalAssemblyXML.visStr());
        // scepcalAssemblyVol.placeVolume( phiRingAssemblyVolume1, Transform3D(rotMirror, dispCone1) );

        // for (int iPhi=0; iPhi<nPhiEndcap; iPhi++) {
        for (int iPhi=3; iPhi<4; iPhi++) {

          for (int nGamma=0; nGamma<nPhiEndcapCrystal; nGamma++) {
            
            // double aEndcap = RinEndcap*cos(dPhiEndcap/2);
            double phiGamma = dPhiEndcap/2-dPhiEndcapCrystal/2-dPhiEndcapCrystal*nGamma;
            double aEndcap = RinEndcap;

            // double thCe = atan(RinBnE/EBz);

            // Make crystal shapes
            double r0adj= aEndcap*abs(tan(phiGamma));
            double r0= sqrt(r0adj*r0adj + r0e*r0e);

            double r1=r0+Fdz;
            double r2=r1+Rdz;

            double ratio0te = r0/r0e;
            double ratio1te = r1/r0e;
            double ratio2te = r2/r0e;

            double r1adj = (r1e/r0e)*r0adj;
            double r2adj = (r2e/r0e)*r0adj;

            double y0=r0*tan(dThetaEndcap/2.);
            double y1=r1*tan(dThetaEndcap/2.);
            double y2=r2*tan(dThetaEndcap/2.);


            double RinBnEl = aEndcap/ cos(phiGamma-dPhiEndcapCrystal/2);
            double RinBnEr = aEndcap/ cos(phiGamma+dPhiEndcapCrystal/2);

            double y0adjl = -(RinBnEl-aEndcap)*cos(phiGamma-dPhiEndcapCrystal/2);
            double y0adjr = -(RinBnEr-aEndcap)*cos(phiGamma+dPhiEndcapCrystal/2);

            double y1adjl = y0adjl*(r1/r0);
            double y1adjr = y0adjr*(r1/r0);

            double y2adjl = y0adjl*(r2/r0);
            double y2adjr = y0adjr*(r2/r0);

            double sf = cos(phiGamma);
            // double sf = 1;

            double alpha = dPhiEndcapCrystal/2;
            double beta = M_PI/2 -alpha;
            double gamma = abs(phiGamma);
            double epsilon = M_PI/2 -gamma;

            double rho = M_PI -(beta+epsilon);
            double rho_l = M_PI/2 -beta;

            double omega = M_PI -(epsilon+rho);
            double nu = M_PI/2 - omega;

            // double a = x0y0;
            // double b = a*cos(gamma);
            // double c = a*sin(gamma);

            // double e = c/cos(rho);
            // double e_l = c/cos(nu);

            // double f_l = e_l*cos(beta);
            // double f = e/cos(omega);

            // double g_l = e_l*cos(rho_l);
            // double g = e/cos(nu);

            double x0y0 = (r0e*cos(thC) +y0e*sin(thC)) *tan(thC -dThetaEndcap/2.) *tan(dPhiEndcapCrystal/2.);
            double x0y0corrl = ((x0y0*sin(gamma))/cos(nu))*cos(beta); //f_l
            double x0y0corrr = ((x0y0*sin(gamma))/cos(rho))/cos(omega); //f_r
            double y0lcorrx0 = ((x0y0*sin(gamma))/cos(nu))*cos(rho_l); //g_l
            double y0rcorrx0 = ((x0y0*sin(gamma))/cos(rho))/cos(nu); //g_r

            double x1y0 = (r0e*cos(thC) -y0e*sin(thC)) *tan(thC +dThetaEndcap/2.) *tan(dPhiEndcapCrystal/2.);
            double x1y0corrl = ((x1y0*sin(gamma))/cos(nu))*cos(beta); //f_l
            double x1y0corrr = ((x1y0*sin(gamma))/cos(rho))/cos(omega); //f_r
            double y0lcorrx1 = ((x1y0*sin(gamma))/cos(nu))*cos(rho_l); //g_l
            double y0rcorrx1 = ((x1y0*sin(gamma))/cos(rho))/cos(nu); //g_r

            double x0y1 = (r1e*cos(thC) +y1e*sin(thC)) *tan(thC -dThetaEndcap/2.) *tan(dPhiEndcapCrystal/2.);
            double x0y1corrl = ((x0y1*sin(gamma))/cos(nu))*cos(beta); //f_l
            double x0y1corrr = ((x0y1*sin(gamma))/cos(rho))/cos(omega); //f_r
            double y1lcorrx0 = ((x0y1*sin(gamma))/cos(nu))*cos(rho_l); //g_l
            double y1rcorrx0 = ((x0y1*sin(gamma))/cos(rho))/cos(nu); //g_r

            double x1y1 = (r1e*cos(thC) -y1e*sin(thC)) *tan(thC +dThetaEndcap/2.) *tan(dPhiEndcapCrystal/2.);
            double x1y1corrl = ((x1y1*sin(gamma))/cos(nu))*cos(beta); //f_l
            double x1y1corrr = ((x1y1*sin(gamma))/cos(rho))/cos(omega); //f_r
            double y1lcorrx1 = ((x1y1*sin(gamma))/cos(nu))*cos(rho_l); //g_l
            double y1rcorrx1 = ((x1y1*sin(gamma))/cos(rho))/cos(nu); //g_r

            double x0y2 = (r2e*cos(thC) +y2e*sin(thC)) *tan(thC -dThetaEndcap/2.) *tan(dPhiEndcapCrystal/2.);
            double x0y2corrl = ((x0y2*sin(gamma))/cos(nu))*cos(beta); //f_l
            double x0y2corrr = ((x0y2*sin(gamma))/cos(rho))/cos(omega); //f_r
            double y2lcorrx0 = ((x0y2*sin(gamma))/cos(nu))*cos(rho_l); //g_l
            double y2rcorrx0 = ((x0y2*sin(gamma))/cos(rho))/cos(nu); //g_r

            double x1y2 = (r2e*cos(thC) -y2e*sin(thC)) *tan(thC +dThetaEndcap/2.) *tan(dPhiEndcapCrystal/2.);
            double x1y2corrl = ((x1y2*sin(gamma))/cos(nu))*cos(beta); //f_l
            double x1y2corrr = ((x1y2*sin(gamma))/cos(rho))/cos(omega); //f_r
            double y2lcorrx1 = ((x1y2*sin(gamma))/cos(nu))*cos(rho_l); //g_l
            double y2rcorrx1 = ((x1y2*sin(gamma))/cos(rho))/cos(nu); //g_r

            double verticesF[]={
                                x0y0 -x0y0corrl,  y0e +y0lcorrx0,
                                x1y0 -x1y0corrl, -y0e +y0lcorrx1,
                               -x1y0 -x1y0corrr, -y0e -y0rcorrx1,
                               -x0y0 -x0y0corrr,  y0e -y0rcorrx0,

                                 x0y1 -x0y1corrl,  y1e +y1lcorrx0,
                                 x1y1 -x1y1corrl, -y1e +y1lcorrx1,
                                -x1y1 -x1y1corrr, -y1e -y1rcorrx1,
                                -x0y1 -x0y1corrr,  y1e -y1rcorrx0
                                };

            double verticesR[]={
                                 x0y1 -x0y1corrl,  y1e +y1lcorrx0,
                                 x1y1 -x1y1corrl, -y1e +y1lcorrx1,
                                -x1y1 -x1y1corrr, -y1e -y1rcorrx1,
                                -x0y1 -x0y1corrr,  y1e -y1rcorrx0,

                                 x0y2 -x0y2corrl,  y2e +y2lcorrx0,
                                 x1y2 -x1y2corrl, -y2e +y2lcorrx1,
                                -x1y2 -x1y2corrr, -y2e -y2rcorrx1,
                                -x0y2 -x0y2corrr,  y2e -y2rcorrx0
                                };

            dd4hep::EightPointSolid crystalFShape(Fdz/2, verticesF);
            dd4hep::EightPointSolid crystalRShape(Rdz/2, verticesR);

            dd4hep::Volume crystalFVol("EndcapCrystalR", crystalFShape, crystalFMat);
            dd4hep::Volume crystalRVol("EndcapCrystalF", crystalRShape, crystalRMat);

            crystalFVol.setVisAttributes(theDetector, crystalFXML.visStr());
            crystalRVol.setVisAttributes(theDetector, crystalRXML.visStr());
            
            if (debugLevel>1) std::cout << "  Endcap: phi: " << iPhi << std::endl;

            double phi=iPhi*dPhiEndcap+dPhiEndcapCrystal/2.+nGamma*dPhiEndcapCrystal +dPhiEndcap/2;
            RotationZYX rot(M_PI/2, thC, 0);
            ROOT::Math::RotationZ rotZ = ROOT::Math::RotationZ(phi);
            rot = rotZ*rot;

            double rF=r0+Fdz/2.;
            Position dispF(rF*sin(thC)*cos(phi),
                          rF*sin(thC)*sin(phi),
                          rF*cos(thC));

            double rR=r1+Rdz/2.;
            Position dispR(rR*sin(thC)*cos(phi),
                          rR*sin(thC)*sin(phi),
                          rR*cos(thC));

            auto crystalFId64=segmentation->setVolumeID(1, iTheta, iPhi*nPhiEndcap+nGamma, 1);
            auto crystalRId64=segmentation->setVolumeID(1, iTheta, iPhi*nPhiEndcap+nGamma, 2);
            int crystalFId32=segmentation->getFirst32bits(crystalFId64);
            int crystalRId32=segmentation->getFirst32bits(crystalRId64);

            dd4hep::PlacedVolume crystalFp = phiRingAssemblyVolume.placeVolume( crystalFVol, crystalFId32, Transform3D(rot,dispF-dispCone) );
            dd4hep::PlacedVolume crystalRp = phiRingAssemblyVolume.placeVolume( crystalRVol, crystalRId32, Transform3D(rot,dispR-dispCone) );

            crystalFp.addPhysVolID("eta", iTheta);
            crystalFp.addPhysVolID("phi", iPhi*nPhiEndcap+nGamma);
            crystalFp.addPhysVolID("depth", 1);
            crystalFp.addPhysVolID("system", 1);

            crystalRp.addPhysVolID("eta", iTheta);
            crystalRp.addPhysVolID("phi", iPhi*nPhiEndcap+nGamma);
            crystalRp.addPhysVolID("depth", 2);
            crystalRp.addPhysVolID("system", 1);

            // Add suffix 1 for the other endcap (mirrored)

            // auto crystalFId641=segmentation->setVolumeID(1, nThetaEndcap+nThetaBarrel+iTheta, iPhi*nPhiEndcap+nGamma, 1);
            // auto crystalRId641=segmentation->setVolumeID(1, nThetaEndcap+nThetaBarrel+iTheta, iPhi*nPhiEndcap+nGamma, 2);
            // int crystalFId321=segmentation->getFirst32bits(crystalFId641);
            // int crystalRId321=segmentation->getFirst32bits(crystalRId641);

            // dd4hep::PlacedVolume crystalFp1 = phiRingAssemblyVolume1.placeVolume( crystalFVol, crystalFId321, Transform3D(rot,dispF-dispCone) );
            // dd4hep::PlacedVolume crystalRp1 = phiRingAssemblyVolume1.placeVolume( crystalRVol, crystalRId321, Transform3D(rot,dispR-dispCone) );

            // crystalFp1.addPhysVolID("eta", nThetaEndcap+nThetaBarrel+iTheta);
            // crystalFp1.addPhysVolID("phi", iPhi*nPhiEndcap+nGamma);
            // crystalFp1.addPhysVolID("depth", 1);
            // crystalFp1.addPhysVolID("system", 1);

            // crystalRp1.addPhysVolID("eta", nThetaEndcap+nThetaBarrel+iTheta);
            // crystalRp1.addPhysVolID("phi", iPhi*nPhiEndcap+nGamma);
            // crystalRp1.addPhysVolID("depth", 2);
            // crystalRp1.addPhysVolID("system", 1);

            std::bitset<10> _eta(iTheta);
            std::bitset<10> _eta1(nThetaEndcap+nThetaBarrel+iTheta);
            std::bitset<10> _phi(iPhi*nPhiEndcap+nGamma);
            std::bitset<3>  depthF(1);
            std::bitset<3>  depthR(2);
            std::bitset<32> id32F(crystalFId32);
            std::bitset<32> id32R(crystalRId32);

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

            //std::cout << "E crystalF1 eta: " << nThetaEndcap+nThetaBarrel+iTheta << " phi: " << iPhi << " depth: " << 1 << std::endl;
            //std::cout << "E crystalF1 eta: " << _eta1 << " phi: " << _phi << " depth: " << depthF << std::endl;
            //std::cout << "E crystalFId321: " << id32F1 << std::endl;
            //std::cout << "E crystalF1 copyNum: " << crystalFp1.copyNumber() << " VolIDs: " << dd4hep::detail::tools::toString(crystalFp1.VolIDs()) << std::endl;

            //std::cout << "E crystalR1 eta: " << nThetaEndcap+nThetaBarrel+iTheta << " phi: " << iPhi << " depth: " << 2 << std::endl;
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
