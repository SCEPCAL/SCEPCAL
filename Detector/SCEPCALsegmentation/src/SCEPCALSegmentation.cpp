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
    registerIdentifier("identifier_sipm", "Cell ID identifier for numSipm", fSipmId, "sipm");
    registerIdentifier("identifier_IsCerenkov", "Cell ID identifier for IsCerenkov", fIsCerenkovId, "c");
    registerIdentifier("identifier_module", "Cell ID identifier for module", fModule, "module");
}

SCEPCALSegmentation::SCEPCALSegmentation(const BitFieldCoder* decoder) : Segmentation(decoder) {
    _type = "SCEPCALSegmentation";
    _description = "SCEPCAL segmentation based on side/eta/phi/depth/S/C";

    registerIdentifier("identifier_system", "Cell ID identifier for numSystem", fSystemId, "system");
    registerIdentifier("identifier_eta", "Cell ID identifier for Eta", fEtaId, "eta");
    registerIdentifier("identifier_phi", "Cell ID identifier for Phi", fPhiId, "phi");
    registerIdentifier("identifier_depth", "Cell ID identifier for Depth", fDepthId, "depth");
    registerIdentifier("identifier_sipm", "Cell ID identifier for Sipm", fSipmId, "sipm");
    registerIdentifier("identifier_IsCerenkov", "Cell ID identifier for IsCerenkov", fIsCerenkovId, "c");
    registerIdentifier("identifier_module", "Cell ID identifier for module", fModule, "module");
}

SCEPCALSegmentation::~SCEPCALSegmentation() {}

Vector3D SCEPCALSegmentation::position(const CellID& cID) const {
    return Vector3D(0,0,0);
};

Vector3D SCEPCALSegmentation::myPosition(const CellID& cID) {

    int copyNum = (int)cID;

    if (fPositionOf.count(copyNum) == 0) { //Add if not found
        //int system=(copyNum)&(32-1);
        //int nEta_in=(copyNum>>5)&(1024-1);
        //int nPhi_in=(copyNum>>15)&(1024-1);
        //int nDepth_in=(copyNum>>25)&(8-1);

        int system=System(copyNum);
        int nEta_in=Eta(copyNum);
        int nPhi_in=Phi(copyNum);
        int nDepth_in=Depth(copyNum);

        // now in mm
        double Fdz=50;
        double Rdz=150;
        double nomfw=50;
        double nomth=10;
        double EBz=2250;
        double Rin=2000;

        // Begin geometry calculations

        // Need odd number of nThetaBarrel to make center crystal
        int nThetaBarrel  =int(floor(2*EBz/nomfw))%2==1? floor(2*EBz/nomfw) : floor(2*EBz/nomfw)-1 ;
        int nThetaEndcap  =floor(Rin/nomfw);

        double thetaSizeEndcap=atan(Rin/EBz);

        double dThetaBarrel =(M_PI-2*thetaSizeEndcap)/(nThetaBarrel);
        double dThetaEndcap =thetaSizeEndcap/nThetaEndcap;

        int    nPhiBarrel = 16;
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
                fPositionOf.emplace(copyNum,dispTimingAssembly+rotZ*dispLg);
                
            }
            else if (nTile<0) {
                ROOT::Math::XYZVector dispTr(sign*phiSign*(-nomth/2),
                                0,
                                -y1slice +abs(nTile)*lT +actY/2 +abs(nC)*actY
                                );
                fPositionOf.emplace(copyNum,dispTimingAssembly+rotZ*dispTr);
            }
        }

        // Endcap
        else if (nEta_in < nThetaEndcap || nEta_in > nThetaEndcap+nThetaBarrel) {

            double thC        = dThetaEndcap/2+ (nEta_in%(nThetaEndcap+nThetaBarrel))*dThetaEndcap;
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
                fPositionOf.emplace(copyNum, rotY*(dispF+rotZ*dispGamma) );
            }

            else if (nDepth_in==2) {
                double rR=r1e+Rdz/2.;
                ROOT::Math::XYZVector dispGamma(0, rR*sin(thC)*tan(gamma), 0);
                ROOT::Math::XYZVector dispR(rR*sin(thC)*cos(phi),
                            rR*sin(thC)*sin(phi),
                            rR*cos(thC));
                fPositionOf.emplace(copyNum, rotY*(dispR+rotZ*dispGamma) );
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
                fPositionOf.emplace(copyNum, dispSlice+ rotZ*dispF );
            }

            else if (nDepth_in==2) {
                double rR=r1e+Rdz/2.;
                ROOT::Math::XYZVector dispR(
                              rR*sin(thC)-rSlice,
                              rR*sin(thC)*tan(gamma),
                              rR*cos(thC)
                          );
                fPositionOf.emplace(copyNum, dispSlice+ rotZ*dispR );
            }
        }

        // nTheta always positive, so that thC is going from 0 to pi.
        // [0,pi/2) --> z, eta positive ; (pi/2,pi] --> z,eta negative; pi/2 (cos(theta) = 0) is excluded
        // int nTheta = (nThetaBarrel+nThetaEndcap)-nEta_in;
        // double thC=nTheta*dTheta;
            
        //phi should run form 0 to 2pi, so that x,y are both neg and positive.
        // double phi=nPhi_in*dPhi;

        //double r0=EBz/cos(thC); // May be this is the case for the endcaps
        // for EE fixed proj along EBz of r0, while for EB proj along transv plane is fixed:
        // double r0 = abs(nEta_in) > nThetaBarrel ? EBz/abs(cos(thC)) : Rin/abs(sin(thC));

        //double r1=r0+Fdz;
        //double r2=r1+Rdz;
        // double rF=r0+Fdz/2; // r is connecting the IP to center of crystal
        // double rR=rF+Rdz/2;

        //double R=nDepth_in==1 ? (r0+r1)/2 : (r1+r2)/2;
        // double R=nDepth_in==1 ? rF : rR;
        // double x=R*sin(thC)*cos(phi);
        // double y=R*sin(thC)*sin(phi);
        // double z=R*cos(thC);

        //std::cout << "These are nEta_in :: " << nEta_in << " Phi_in ::" << nPhi_in << " Depth :: " << nDepth_in << std::endl;
        //std::cout << "So these are nTheta :: " << nTheta << "and nPhi:: "<< nPhi_in << std::endl;
        //std::cout << "rF :: " << rF << " rR:: " << rR << std::endl;
        //std::cout << "... theta :: " << thC <<  "and  R :: " << R <<  "... finally cos(theta),sin(theta) :: " << cos(thC) << "," << sin(thC) << std::endl;
        //std::cout << "... phi   :: " << phi <<  "and  R :: " << R <<  "... finally cos(phi),sin(phi) :: " << cos(phi) << "," <<  sin(phi) << std::endl;
        //std::cout << "to get  z=R*cos(theta)= " << z << " and y=R*sin(theta)*sin(phi)= " << y <<  " and x=R*sin(theta)*cos(phi)= " << x << std::endl;

        // ROOT::Math::XYZVector position(x, y, z);
        // fPositionOf.emplace(copyNum,position);


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

    }

    return fPositionOf.at(copyNum);
}

/*
ROOT::Math::XYZVector SCEPCALSegmentation::localPosition(const CellID& cID) const {
    int numx = numX(cID);
    int numy = numY(cID);
    int numz = numZ(cID);
    int x_ = x(cID);
    int y_ = y(cID);
    int z_ = z(cID);

    return localPosition(numx,numy,numz,x_,y_,z_);
}

Vector3D SCEPCALSegmentation::localPosition(int numx, int numy, int numz, int x_, int y_, int z_) const {
    float ptX = -fGridSize*static_cast<float>(numx/2) + static_cast<float>(x_)*fGridSize + ( numx%2==0 ? fGridSize/2. : 0. );
    float ptY = -fGridSize*static_cast<float>(numy/2) + static_cast<float>(y_)*fGridSize + ( numy%2==0 ? fGridSize/2. : 0. );
    float ptZ = -fGridSize*static_cast<float>(numz/2) + static_cast<float>(z_)*fGridSize + ( numz%2==0 ? fGridSize/2. : 0. );

    return Vector3D(ptX,ptY,ptZ);
}
*/

CellID SCEPCALSegmentation::cellID(const Vector3D& /*localPosition*/, const Vector3D& /*globalPosition*/, const VolumeID& vID) const {
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

    VolumeID module = 0; // Tower
    _decoder->set(vID, fModule, module);

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

    VolumeID module = 1; // Fiber, SiPM, etc.
    _decoder->set(vID, fModule, module);

    VolumeID isCeren = IsCerenkov(Eta,Phi) ? 1 : 0;
    _decoder->set(vID, fIsCerenkovId, isCeren);

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

int SCEPCALSegmentation::Sipm(const CellID& aCellID) const {
    VolumeID Sipm = static_cast<VolumeID>(_decoder->get(aCellID, fSipmId));
    return static_cast<int>(Sipm);
}

bool SCEPCALSegmentation::IsCerenkov(const CellID& aCellID) const {
    VolumeID isCeren = static_cast<VolumeID>(_decoder->get(aCellID, fIsCerenkovId));
    return static_cast<bool>(isCeren);
}

bool SCEPCALSegmentation::IsCerenkov(int eta, int phi) const {
    bool isCeren = false;
    if ( eta%2 == 1 ) { isCeren = !isCeren; }
    if ( phi%2 == 1 ) { isCeren = !isCeren; }
    return isCeren;
}

bool SCEPCALSegmentation::IsTower(const CellID& aCellID) const {
    VolumeID module = static_cast<VolumeID>(_decoder->get(aCellID, fModule));
    return module==0;
}

bool SCEPCALSegmentation::IsSiPM(const CellID& aCellID) const {
    VolumeID module = static_cast<VolumeID>(_decoder->get(aCellID, fModule));
    return module==1;
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

/*        
DRparamBase* SCEPCALSegmentation::setParamBase(int noEta) const {
    DRparamBase* paramBase = nullptr;

    if ( fParamEndcap->unsignedTowerNo(noEta) >= fParamBarrel->GetTotTowerNum() ) paramBase = static_cast<DRparamBase*>(fParamEndcap);
    else paramBase = static_cast<DRparamBase*>(fParamBarrel);

    if ( paramBase->GetCurrentTowerNum()==noEta ) return paramBase;

    // This should not be called while building detector geometry
    if (!paramBase->IsFinalized()) throw std::runtime_error("SCEPCALSegmentation::position should not be called while building detector geometry!");

    paramBase->SetDeltaThetaByTowerNo(noEta, fParamBarrel->GetTotTowerNum());
    paramBase->SetThetaOfCenterByTowerNo(noEta, fParamBarrel->GetTotTowerNum());
    paramBase->SetIsRHSByTowerNo(noEta);
    paramBase->SetCurrentTowerNum(noEta);
    paramBase->init();

    return paramBase;
}
*/

}
}
