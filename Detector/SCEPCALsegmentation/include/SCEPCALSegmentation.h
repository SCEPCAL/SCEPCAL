#ifndef SCEPCALSegmentation_h
#define SCEPCALSegmentation_h 1

#include "DDSegmentation/Segmentation.h"

#include "TVector3.h"
#include "DD4hep/DetFactoryHelper.h"

#include <vector>
#include <cmath>

namespace dd4hep {
namespace DDSegmentation {

class SCEPCALSegmentation : public Segmentation {
    public:
        SCEPCALSegmentation(const std::string& aCellEncoding);
        SCEPCALSegmentation(const BitFieldCoder* decoder);
        virtual ~SCEPCALSegmentation() override;

        virtual Vector3D position(const CellID& aCellID) const;

        virtual Vector3D myPosition(const CellID& aCellID) ;

        virtual CellID cellID(const Vector3D& aLocalPosition,
                              const Vector3D& aGlobalPosition,
                              const VolumeID& aVolumeID) const;

        VolumeID setVolumeID(int System, int Eta, int Phi, int Depth) const;
        CellID setCellID(int System, int Eta, int Phi, int Depth) const;

        int System(const CellID& aCellID) const;
        int Eta(const CellID& aCellID) const;
        int Phi(const CellID& aCellID) const;
        int Depth(const CellID& aCellID) const;

        int getFirst32bits(const CellID& aCellID) const { return (int)aCellID; }
        int getLast32bits(const CellID& aCellID) const;

        CellID convertFirst32to64(const int aId32) const { return (CellID)aId32; }
        CellID convertLast32to64(const int aId32) const;

        int System(const int& aId32) const { return System( convertFirst32to64(aId32) ); }
        int Eta(const int& aId32) const { return Eta( convertFirst32to64(aId32) ); }
        int Phi(const int& aId32) const { return Phi( convertFirst32to64(aId32) ); }
        int Depth(const int& aId32) const { return Depth( convertFirst32to64(aId32) ); }

        inline void setGeomParams(double Fdz, double Rdz,
                                  double nomfw, double nomth, 
                                  double EBz, double Rin,
                                  int phiSegments) {
            f_Fdz = Fdz;
            f_Rdz = Rdz;
            f_nomfw = nomfw;
            f_nomth = nomth;
            f_EBz = EBz;
            f_Rin = Rin;
            f_phiSegments = phiSegments;
        }

        inline double getFdz() const{ return f_Fdz; }
        inline double getRdz() const{ return f_Rdz; }
        inline double getnomfw() const{ return f_nomfw; }
        inline double getnomth() const{ return f_nomth; }
        inline double getEBz() const{ return f_EBz; }
        inline double getRin() const{ return f_Rin; }
        inline int    getphiSegments() const{ return f_phiSegments; }

        // inline const std::string& fieldNameEta() const { return fEtaId; }
        // inline const std::string& fieldNamePhi() const { return fPhiId; }
        // inline const std::string& fieldNameDepth() const { return fDepthId; }

        // inline void setFieldNameEta(const std::string& fieldName) { fEtaId = fieldName; }
        // inline void setFieldNamePhi(const std::string& fieldName) { fPhiId = fieldName; }
        // inline void setFieldNameDepth(const std::string& fieldName) { fDepthId = fieldName; }


    protected:
        std::string fSystemId;
        std::string fEtaId;
        std::string fPhiId;
        std::string fDepthId;

        double f_Fdz;
        double f_Rdz;
        double f_nomfw;
        double f_nomth;
        double f_EBz;
        double f_Rin;
        int    f_phiSegments;

    private:
        std::unordered_map<int, Vector3D> fPositionOf;


};

}
}

#endif
