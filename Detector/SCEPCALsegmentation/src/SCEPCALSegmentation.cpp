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

    // now in mm
    double Fdz      = getFdz();
    double Rdz      = getRdz();
    double nomfw    = getnomfw();
    double nomth    = getnomth();
    double EBz      = getEBz();
    double Rin      = getRin();
    int PHI_SEGMENTS = getphiSegments();

    int system      = System(copyNum);
    int nEta_in     = Eta(copyNum);
    int nPhi_in     = Phi(copyNum);
    int nDepth_in   = Depth(copyNum);

    if (system==3) return Vector3D(0,0,0);

    // if (fPositionOf.count(copyNum) == 0) { //Add if not found

    //-----------------------------------------------------------------------------------
    // Global geometry numbers
    //-----------------------------------------------------------------------------------

        double  D_PHI_GLOBAL = 2*M_PI/PHI_SEGMENTS;

        // Need odd number of N_THETA_BARREL to make center slice
        double  THETA_SIZE_ENDCAP     = atan(Rin/EBz);

        int     N_THETA_BARREL        = int(floor(2*EBz/nomfw))%2==1? floor(2*EBz/nomfw) : floor(2*EBz/nomfw)-1 ;
        int     N_THETA_ENDCAP        = floor(Rin/nomfw);

        double  D_THETA_BARREL        = (M_PI-2*THETA_SIZE_ENDCAP)/(N_THETA_BARREL);
        double  D_THETA_ENDCAP        = THETA_SIZE_ENDCAP/N_THETA_ENDCAP;

        int     N_PHI_BARREL_CRYSTAL  = floor(2*M_PI*Rin/(PHI_SEGMENTS*nomfw));
        double  D_PHI_BARREL_CRYSTAL  = D_PHI_GLOBAL/N_PHI_BARREL_CRYSTAL;

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

        double y1slice          = z1slice*tan(M_PI/2-THETA_SIZE_ENDCAP);

        // Timing layer
        if (nDepth_in==0 || nDepth_in==3) {

            double phiTiming=nPhi_in*D_PHI_GLOBAL;

            // Timing layer envelope
            double  rT      = z1slice -2*nomth;
            double  wT      = rT *tan(D_PHI_GLOBAL/2);
            int     nTiles  = ceil(y1slice/wT);
            double  lT      = 2*y1slice/nTiles;
            int     nCy     = floor(lT/nomth);
            double  actY    = lT/nCy;
            double  actX    = 2*wT/nCy; 

            double rTimingAssembly = rT+nomth;
            ROOT::Math::XYZVector dispTimingAssembly(rTimingAssembly*cos(phiTiming),
                                        rTimingAssembly*sin(phiTiming),
                                        0);

            int nTile = int(nEta_in/nCy);
            int nC = nEta_in - nTile*nCy;

            int phiTimingSign =    nPhi_in%2==0? 1:-1;
            int sign    = abs(nTile)%2==0? 1:-1;

            ROOT::Math::RotationZ rotZ(phiTiming);

            // if (nTile>0) {
            if (nDepth_in==0) {

                ROOT::Math::XYZVector dispLg(sign*phiTimingSign*(nomth/2),
                                -wT +actX/2 + abs(nC)*actX,
                                -y1slice +abs(nTile)*lT + lT/2
                                );
                // fPositionOf.emplace(copyNum,dispTimingAssembly+rotZ*dispLg);
                return (dispTimingAssembly+rotZ*dispLg);
                
            }
            // else if (nTile<0) {
            else if (nDepth_in==3) {

                ROOT::Math::XYZVector dispTr(sign*phiTimingSign*(-nomth/2),
                                0,
                                -y1slice +abs(nTile)*lT +actY/2 +abs(nC)*actY
                                );
                // fPositionOf.emplace(copyNum,dispTimingAssembly+rotZ*dispTr);
                return (dispTimingAssembly+rotZ*dispTr);

            }
        }

        // Endcap
        else if (nEta_in < N_THETA_ENDCAP || nEta_in > N_THETA_ENDCAP+N_THETA_BARREL) {
            double thC;

            if (nEta_in < N_THETA_ENDCAP) {
                thC = D_THETA_ENDCAP/2+ nEta_in*D_THETA_ENDCAP;
            }
            else if (nEta_in > N_THETA_ENDCAP+N_THETA_BARREL) {
                thC = M_PI -THETA_SIZE_ENDCAP +D_THETA_ENDCAP/2 + (nEta_in -(N_THETA_ENDCAP+N_THETA_BARREL)) *D_THETA_ENDCAP;
                // thC = D_THETA_ENDCAP/2+ (N_THETA_ENDCAP -(nEta_in%(N_THETA_ENDCAP+N_THETA_BARREL)) )*D_THETA_ENDCAP;
            }
            // double thC        = D_THETA_ENDCAP/2+ (nEta_in%(N_THETA_ENDCAP+N_THETA_BARREL))*D_THETA_ENDCAP;
            double RinEndcap  = EBz*tan(thC);

            int    nPhiEndcapCrystal = floor(2*M_PI*RinEndcap/(PHI_SEGMENTS*nomfw));
            double dPhiEndcapCrystal = D_PHI_GLOBAL/nPhiEndcapCrystal;

            double r0e=RinEndcap/sin(thC);
            double r1e=r0e+Fdz;

            int nPhi   = int(nPhi_in/nPhiEndcapCrystal);
            int nGamma = nPhi_in%nPhiEndcapCrystal;

            double phi   = nPhi*D_PHI_GLOBAL;
            double gamma = -D_PHI_GLOBAL/2+dPhiEndcapCrystal/2+dPhiEndcapCrystal*nGamma;
            
            int mirror = nEta_in<N_THETA_ENDCAP? 0:1;

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
        else if (nEta_in>=N_THETA_ENDCAP && nEta_in<N_THETA_ENDCAP+N_THETA_BARREL) {

            int nTheta = nEta_in -N_THETA_ENDCAP;
            int nPhi   = int(nPhi_in/N_PHI_BARREL_CRYSTAL);
            int nGamma = nPhi_in%N_PHI_BARREL_CRYSTAL;

            double phi   = nPhi*D_PHI_GLOBAL;
            double gamma = -D_PHI_GLOBAL/2+D_PHI_BARREL_CRYSTAL/2+D_PHI_BARREL_CRYSTAL*nGamma;

            double thC =THETA_SIZE_ENDCAP+D_THETA_BARREL/2+(nTheta*D_THETA_BARREL);

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
        //                 -wT +actX/2 + abs(nCt)*actX,
        //                 -y1slice +abs(nTilet)*lT + lT/2
        //                 );

        // dd4hep::Box testBox(nomth, nomth, nomth);
        // dd4hep::Volume testBoxVol("BarrelCrystalF", testBox, timingLgMat);
        // testBoxVol.setVisAttributes(theDetector, timingLgXML.visStr());
        // dd4hep::PlacedVolume timingLgp = scepcalAssemblyVol.placeVolume( testBoxVol, timingLgId32,  dispTimingAssembly+rotPhi*dispLgt );


        // ROOT::Math::RotationZ rotPhi(phi);
        // ROOT::Math::RotationY rotY(M_PI*0);
        // Position dispGamma(0, rF*sin(thC)*tan(gamma), 0);
        // Position dispFt(rF*sin(thC)*cos(phi),
        //             rF*sin(thC)*sin(phi),
        //             rF*cos(thC));

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
