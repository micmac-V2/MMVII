#ifndef CCODEDTARGETDESCRIBE_H
#define CCODEDTARGETDESCRIBE_H

#endif // CCODEDTARGETDESCRIBE_H

//#include "MMVII_Sensor.h"
#include "CodedTarget.h"
#include "MMVII_PCSens.h"
#include "cMMVII_Appli.h"
#include "MMVII_Geom3D.h"

namespace MMVII
{

extern const std::vector<cPt2dr> eSqCorners;//-> square corners coordinates from square centre
struct cCdTDetec;
class cAugCdT;

typedef cSegment<tREAL8,3> tSeg3dr;
typedef cIm2D<tU_INT1>      tIm;
typedef cDataIm2D<tU_INT1>  tDIm;


/*!
 * @brief The cExtract class stores CdT extractions
 */
struct cExtract
{
    cExtract(const cSensorCamPC* aCam, cSaveExtrEllipe aEll);
    const cSensorCamPC* mCam;
    const cSaveExtrEllipe mEll;
};

/*!
 * @brief The cAugCdT class performs CdT augmentation
 */
class cAugCdT
{
public:
    cAugCdT(std::string aName, std::shared_ptr<cFullSpecifTarget> aFSpec);
    cAugCdT();
    void                Spatialize(tREAL8 aGndInterTol=1e-2);
    void                AddExtract(cExtract aExt);
    tU_INT1             NbExtracts() const;
    void                AddData(const cAuxAr2007& anAux);
    std::string         Show() const;
    static std::string  NameFile(const cPhotogrammetricProject & aPhProj, bool Input);
    bool operator       <(const cAugCdT& aAug) const;
    std::vector<cPt3dr> GndCorners() const;
    std::string             mName;
    bool                    mOKAug;
    bool                    mOKInter;
    tREAL8                  m3DPrec;
    cSimilitud3D<tREAL8>    mRef2Gnd;

private:
    std::vector<cPt2dr>                 Corners() const;
    std::shared_ptr<cFullSpecifTarget>  mFSpec;
    std::vector<cExtract>               mVExtracts;
};

void AddData(const cAuxAr2007& anAux, cAugCdT& anEx);
}
