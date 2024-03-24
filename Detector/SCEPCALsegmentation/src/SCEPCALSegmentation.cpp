#include "SCEPCALSegmentation.h"

#include <climits>
#include <cmath>
#include <stdexcept>

namespace dd4hep {
namespace DDSegmentation {

/// default constructor using an encoding string
SCEPCALSegmentation::SCEPCALSegmentation(const std::string& cellEncoding) : Segmentation(cellEncoding) {
    _type = "SCEPCALSegmentation";
    _description = "SCEPCAL segmentation based on side/eta/phi/depth/S/C";

    registerIdentifier("identifier_system", "Cell ID identifier for numSystem", fSystemId, "system");
    registerIdentifier("identifier_eta", "Cell ID identifier for numEta", fEtaId, "eta");
    registerIdentifier("identifier_phi", "Cell ID identifier for numPhi", fPhiId, "phi");
    registerIdentifier("identifier_depth", "Cell ID identifier for numDepth", fDepthId, "depth");
}

SCEPCALSegmentation::SCEPCALSegmentation(const BitFieldCoder* decoder) : Segmentation(decoder) {
    _type = "SCEPCALSegmentation";
    _description = "SCEPCAL segmentation based on side/eta/phi/depth/S/C";

    registerIdentifier("identifier_system", "Cell ID identifier for numSystem", fSystemId, "system");
    registerIdentifier("identifier_eta", "Cell ID identifier for Eta", fEtaId, "eta");
    registerIdentifier("identifier_phi", "Cell ID identifier for Phi", fPhiId, "phi");
    registerIdentifier("identifier_depth", "Cell ID identifier for Depth", fDepthId, "depth");
}

SCEPCALSegmentation::~SCEPCALSegmentation() {}

Vector3D SCEPCALSegmentation::position(const CellID& cID) const {
    return Vector3D(0,0,0);
};

Vector3D SCEPCALSegmentation::myPosition(const CellID& cID) {

    int copyNum = (int)cID;

    // // now in mm
    // double Fdz= f_Fdz; // getFdz();
    // double Rdz= f_Rdz; // getRdz();
    // double nomfw= f_nomfw; // getnomfw();
    // double nomth= f_nomth; // getnomth();
    // double EBz= f_EBz; // getEBz();
    // double Rin= f_Rin; // getRin();
    // int phiSegments= f_phiSegments; // getphiSegments();

    // now in mm
    double Fdz= getFdz();
    double Rdz= getRdz();
    double nomfw= getnomfw();
    double nomth= getnomth();
    double EBz= getEBz();
    double Rin= getRin();
    int phiSegments= getphiSegments();

// std::cout << "Fdz: " << Fdz << std::endl;
// std::cout << "Rdz: " << Rdz << std::endl;
// std::cout << "nomfw: " << nomfw << std::endl;
// std::cout << "nomth: " << nomth << std::endl;
// std::cout << "EBz: " << EBz << std::endl;
// std::cout << "Rin: " << Rin << std::endl;
// std::cout << "phiSegments: " << phiSegments << std::endl;


    int system=System(copyNum);
    int nEta_in=Eta(copyNum);
    int nPhi_in=Phi(copyNum);
    int nDepth_in=Depth(copyNum);

    // std::cout << "Segmentation: " << system << " " << nEta_in << " " << nPhi_in << " " << nDepth_in << std::endl;

    if (system==3) return Vector3D(0,0,0);

    // if (fPositionOf.count(copyNum) == 0) { //Add if not found

        //int system=(copyNum)&(32-1);
        //int nEta_in=(copyNum>>5)&(1024-1);
        //int nPhi_in=(copyNum>>15)&(1024-1);
        //int nDepth_in=(copyNum>>25)&(8-1);


        // Begin geometry calculations

        // Need odd number of nThetaBarrel to make center crystal
        int nThetaBarrel  =int(floor(2*EBz/nomfw))%2==1? floor(2*EBz/nomfw) : floor(2*EBz/nomfw)-1 ;
        int nThetaEndcap  =floor(Rin/nomfw);

        double thetaSizeEndcap=atan(Rin/EBz);

        double dThetaBarrel =(M_PI-2*thetaSizeEndcap)/(nThetaBarrel);
        double dThetaEndcap =thetaSizeEndcap/nThetaEndcap;

        int    nPhiBarrel = phiSegments;
        double dPhiBarrel = 2*M_PI/nPhiBarrel;

        int    nPhiEndcap = nPhiBarrel;
        double dPhiEndcap = dPhiBarrel;
        
        int    nPhiBarrelCrystal   =floor(2*M_PI*Rin/(nPhiBarrel*nomfw));
        double dPhiBarrelCrystal   =dPhiBarrel/nPhiBarrelCrystal;

        // Shared calculations for timing layer and barrel envelopes
        double thC_end = thetaSizeEndcap+dThetaBarrel/2;

        double r0slice_end =Rin/sin(thC_end);
        double y0slice_end =r0slice_end*tan(dThetaBarrel/2.);
        double slice_front_jut = y0slice_end*sin(M_PI/2-thC_end);

        double z1slice =Rin -slice_front_jut;
        double z2slice =Rin +Fdz +Rdz +slice_front_jut;

        double y1slice =z1slice*tan(M_PI/2-thetaSizeEndcap);

        // Timing layer
        if (nDepth_in==0) {

            double phi=nPhi_in*dPhiBarrel;

            double rT = z1slice -2*nomth;
            double w  = rT *tan(dPhiBarrel/2);
            int nTiles= ceil(y1slice/w);
            double lT = 2*y1slice/nTiles;
            int nCy = floor(lT/nomth);
            double actY = lT/nCy;
            double actX = 2*w/nCy; 

            double rTimingAssembly = rT+nomth;
            ROOT::Math::XYZVector dispTimingAssembly(rTimingAssembly*cos(phi),
                                        rTimingAssembly*sin(phi),
                                        0);

            int nTile = int(nEta_in/nCy);
            int nC = nEta_in - nTile*nCy;

            int phiSign =    nPhi_in%2==0? 1:-1;
            int sign    = abs(nTile)%2==0? 1:-1;

            ROOT::Math::RotationZ rotZ(phi);

            if (nTile>0) {
                ROOT::Math::XYZVector dispLg(sign*phiSign*(nomth/2),
                                -w +actX/2 + abs(nC)*actX,
                                -y1slice +abs(nTile)*lT + lT/2
                                );
                // fPositionOf.emplace(copyNum,dispTimingAssembly+rotZ*dispLg);
                return (dispTimingAssembly+rotZ*dispLg);
                
            }
            else if (nTile<0) {
                ROOT::Math::XYZVector dispTr(sign*phiSign*(-nomth/2),
                                0,
                                -y1slice +abs(nTile)*lT +actY/2 +abs(nC)*actY
                                );
                // fPositionOf.emplace(copyNum,dispTimingAssembly+rotZ*dispTr);
                return (dispTimingAssembly+rotZ*dispTr);

            }
        }

        // Endcap
        else if (nEta_in < nThetaEndcap || nEta_in > nThetaEndcap+nThetaBarrel) {
            double thC;

            if (nEta_in < nThetaEndcap) {
                thC = dThetaEndcap/2+ nEta_in*dThetaEndcap;
            }
            else if (nEta_in > nThetaEndcap+nThetaBarrel) {
                thC = dThetaEndcap/2+ (nThetaEndcap -(nEta_in%(nThetaEndcap+nThetaBarrel)) )*dThetaEndcap;
            }
            // double thC        = dThetaEndcap/2+ (nEta_in%(nThetaEndcap+nThetaBarrel))*dThetaEndcap;
            double RinEndcap  = EBz*tan(thC);

            int    nPhiEndcapCrystal = floor(2*M_PI*RinEndcap/(nPhiEndcap*nomfw));
            double dPhiEndcapCrystal = dPhiEndcap/nPhiEndcapCrystal;

            double r0e=RinEndcap/sin(thC);
            double r1e=r0e+Fdz;

            int nPhi   = int(nPhi_in/nPhiEndcapCrystal);
            int nGamma = nPhi_in%nPhiEndcapCrystal;

            double phi   = nPhi*dPhiEndcap;
            double gamma = -dPhiEndcap/2+dPhiEndcapCrystal/2+dPhiEndcapCrystal*nGamma;
            
            int mirror = nEta_in<nThetaEndcap? 0:1;

            ROOT::Math::RotationZ rotZ(phi);
            ROOT::Math::RotationY rotY(M_PI*mirror);

            if (nDepth_in==1) {
                double rF=r0e+Fdz/2.;
                ROOT::Math::XYZVector dispGamma(0, rF*sin(thC)*tan(gamma), 0);
                ROOT::Math::XYZVector dispF(rF*sin(thC)*cos(phi),
                            rF*sin(thC)*sin(phi),
                            rF*cos(thC));
                // fPositionOf.emplace(copyNum, rotY*(dispF+rotZ*dispGamma) );
                return (rotY*(dispF+rotZ*dispGamma)) ;

            }

            else if (nDepth_in==2) {
                double rR=r1e+Rdz/2.;
                ROOT::Math::XYZVector dispGamma(0, rR*sin(thC)*tan(gamma), 0);
                ROOT::Math::XYZVector dispR(rR*sin(thC)*cos(phi),
                            rR*sin(thC)*sin(phi),
                            rR*cos(thC));
                // fPositionOf.emplace(copyNum, rotY*(dispR+rotZ*dispGamma) );
                return (rotY*(dispR+rotZ*dispGamma)) ;

            }
        }

        // Barrel
        else if (nEta_in>=nThetaEndcap && nEta_in<nThetaEndcap+nThetaBarrel) {

            int nTheta = nEta_in -nThetaEndcap;
            int nPhi   = int(nPhi_in/nPhiBarrelCrystal);
            int nGamma = nPhi_in%nPhiBarrelCrystal;

            double phi   = nPhi*dPhiBarrel;
            double gamma = -dPhiBarrel/2+dPhiBarrelCrystal/2+dPhiBarrelCrystal*nGamma;

            double thC =thetaSizeEndcap+dThetaBarrel/2+(nTheta*dThetaBarrel);

            double r0e=Rin/sin(thC);
            double r1e=r0e+Fdz;

            double rSlice =(z1slice+z2slice)/2;
            Position dispSlice(rSlice*cos(phi),
                               rSlice*sin(phi),
                               0);

            ROOT::Math::RotationZ rotZ(phi);

            if (nDepth_in==1) {
                double rF=r0e+Fdz/2.;
                ROOT::Math::XYZVector dispF(
                              rF*sin(thC)-rSlice,
                              rF*sin(thC)*tan(gamma),
                              rF*cos(thC)
                          );
                // fPositionOf.emplace(copyNum, dispSlice+ rotZ*dispF );
                return (dispSlice+ rotZ*dispF) ;

            }

            else if (nDepth_in==2) {
                double rR=r1e+Rdz/2.;
                ROOT::Math::XYZVector dispR(
                              rR*sin(thC)-rSlice,
                              rR*sin(thC)*tan(gamma),
                              rR*cos(thC)
                          );
                // fPositionOf.emplace(copyNum, dispSlice+ rotZ*dispR );
                return (dispSlice+ rotZ*dispR) ;

            }
        }

/**
        // To test positions in detector constructor

        // int nEta_in = nTile*nCy +nC;

        // int nTilet = int(nEta_in/nCy);
        // int nCt = nEta_in - nTilet*nCy;

        // int phiSignt =    iPhi%2==0? 1:-1;
        // int signt    = abs(nTilet)%2==0? 1:-1;

        // ROOT::Math::RotationZ rotPhi(phi);

        // Position dispLgt(signt*phiSignt*(nomth/2),
        //                 -w +actX/2 + abs(nCt)*actX,
        //                 -y1slice +abs(nTilet)*lT + lT/2
        //                 );

        // dd4hep::Box testBox(nomth, nomth, nomth);
        // dd4hep::Volume testBoxVol("BarrelCrystalF", testBox, timingLgMat);
        // testBoxVol.setVisAttributes(theDetector, timingLgXML.visStr());
        // dd4hep::PlacedVolume timingLgp = scepcalAssemblyVol.placeVolume( testBoxVol, timingLgId32,  dispTimingAssembly+rotPhi*dispLgt );




        // Position dispGamma(0, rF*sin(thC)*tan(gamma), 0);
        // Position dispFt(rF*sin(thC)*cos(phi),
        //             rF*sin(thC)*sin(phi),
        //             rF*cos(thC));

        // ROOT::Math::RotationZ rotPhi(phi);
        // ROOT::Math::RotationY rotY(M_PI*0);

        // dd4hep::Box testBox(Fdz/2, Fdz/2, Fdz/2);
        // dd4hep::Volume testBoxVol("BarrelCrystalF", testBox, crystalFMat);
        // testBoxVol.setVisAttributes(theDetector, crystalFXML.visStr());
        // dd4hep::PlacedVolume crystalFp = scepcalAssemblyVol.placeVolume( testBoxVol, crystalFId32,  rotY*(dispFt+rotPhi*dispGamma) );

**/
    // }

    return Vector3D(0,0,0);
    // return fPositionOf.at(copyNum);
}


CellID SCEPCALSegmentation::cellID(const Vector3D& /*localPosition*/, 
                                   const Vector3D& /*globalPosition*/, 
                                   const VolumeID& vID) const {
    return setCellID(System(vID), Eta(vID), Phi(vID), Depth(vID) );
}

VolumeID SCEPCALSegmentation::setVolumeID(int System, int Eta, int Phi, int Depth) const {
    VolumeID SystemId = static_cast<VolumeID>(System);
    VolumeID EtaId = static_cast<VolumeID>(Eta);
    VolumeID PhiId = static_cast<VolumeID>(Phi);
    VolumeID DepthId = static_cast<VolumeID>(Depth);
    VolumeID vID = 0;

    //std::cout << " In setVolumeID:: " << std::endl;
    //std::cout << " EtaID:: " << EtaId <<std::endl;
    //std::cout << " PhiID:: " << PhiId <<std::endl;
    //std::cout << " DepthID:: " << DepthId <<std::endl;

    _decoder->set(vID, fSystemId, SystemId);
    _decoder->set(vID, fEtaId, EtaId);
    _decoder->set(vID, fPhiId, PhiId);
    _decoder->set(vID, fDepthId, DepthId);

    return vID;
}

CellID SCEPCALSegmentation::setCellID(int System, int Eta, int Phi, int Depth) const {
    VolumeID SystemId = static_cast<VolumeID>(System);
    VolumeID EtaId = static_cast<VolumeID>(Eta);
    VolumeID PhiId = static_cast<VolumeID>(Phi);
    VolumeID DepthId = static_cast<VolumeID>(Depth);
    VolumeID vID = 0;

    //std::cout << " In setCellID:: " << std::endl;
    //std::cout << " EtaID:: " << EtaId <<std::endl;
    //std::cout << " PhiID:: " << PhiId <<std::endl;
    //std::cout << " DepthID:: " << DepthId <<std::endl;

    _decoder->set(vID, fSystemId, SystemId);
    _decoder->set(vID, fEtaId, EtaId);
    _decoder->set(vID, fPhiId, PhiId);
    _decoder->set(vID, fDepthId, DepthId);

    return vID;
}

int SCEPCALSegmentation::System(const CellID& aCellID) const {
    VolumeID System = static_cast<VolumeID>(_decoder->get(aCellID, fSystemId));
    return static_cast<int>(System);
}

int SCEPCALSegmentation::Eta(const CellID& aCellID) const {
    VolumeID Eta = static_cast<VolumeID>(_decoder->get(aCellID, fEtaId));
    return static_cast<int>(Eta);
}

int SCEPCALSegmentation::Phi(const CellID& aCellID) const {
    VolumeID Phi = static_cast<VolumeID>(_decoder->get(aCellID, fPhiId));
    return static_cast<int>(Phi);
}

int SCEPCALSegmentation::Depth(const CellID& aCellID) const {
    VolumeID Depth = static_cast<VolumeID>(_decoder->get(aCellID, fDepthId));
    return static_cast<int>(Depth);
}

int SCEPCALSegmentation::getLast32bits(const CellID& aCellID) const {
    CellID aId64 = aCellID >> sizeof(int)*CHAR_BIT;
    int aId32 = (int)aId64;
    return aId32;
}

CellID SCEPCALSegmentation::convertLast32to64(const int aId32) const {
    CellID aId64 = (CellID)aId32;
    aId64 <<= sizeof(int)*CHAR_BIT;
    return aId64;
}

}
}
