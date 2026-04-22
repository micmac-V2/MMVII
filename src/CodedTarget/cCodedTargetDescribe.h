#ifndef CCODEDTARGETDESCRIBE_H
#define CCODEDTARGETDESCRIBE_H

#endif // CCODEDTARGETDESCRIBE_H

//#include "MMVII_Sensor.h"
#include "CodedTarget.h"
#include "MMVII_ImageMorphoMath.h"
#include "MMVII_Interpolators.h"
#include "MMVII_PCSens.h"
#include "cMMVII_Appli.h"
#include "MMVII_Geom3D.h"
#include "MMVII_Matrix.h"

namespace MMVII
{

typedef cSegment<tREAL8,3> tSeg3dr;

class cDesCdT;
struct cDetCdT;

class cAppli_CodedTargetDescribe : public cMMVII_Appli
{
    public:
        //------
        cAppli_CodedTargetDescribe(const std::vector<std::string>& aVArgs,
                                   const cSpecMMVII_Appli& aSpec);
    private:
        //------ MMVII mandatory/usual stuff
        int Exe() override;
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override;
        cPhotogrammetricProject mPhProj;
        std::string mSpecImIn;
        bool mShow;
        std::string mFSpecName;
        std::unique_ptr<cFullSpecifTarget> mFSpec;
        //------ members
        std::vector<cDesCdT> mVDesCdT;
        //------ methods
        void AddDesCdT(std::string aName, std::unique_ptr<cFullSpecifTarget>& aSpec);
};

/** Class to store coded target refined description **/

class cDesCdT
{
    public:
    cDesCdT(std::string aName, std::unique_ptr<cFullSpecifTarget>& aSpec);
        //----- members
        std::string mName;
        std::vector<cDetCdT> mVDetects;
        //----- methods
        void AddDetect(const cSensorCamPC* aCam, cMesIm1Pt aMes, cAff2D_r aAff2D);
        void InterBitCenters(bool& aShow);
        void InterGndCorners(bool& aShow);
        void Estimate3DSimilOnCorners(bool& aShow);
        cPt2dr Gnd2CdT(cPt3dr& aPt, const cDetCdT& aDet);
        cPt3dr CdT2GndByInter(const cPt2dr& aPt, std::vector<tREAL8>* aVRes = nullptr);
        cPt3dr CdT2GndBySimil();
    private:
        //----- members
        const cOneEncoding* mEnc;
        tREAL8 mRes;
        std::vector<cPt2dr> mVBitCenters2D;
        std::vector<cPt2dr> mVCdTCorners;
        std::vector<cPt3dr> mVCdtCorners3D;
        std::vector<cPt3dr> mVGndCorners;
        cSimilitud3D<tREAL8> m3DSimil;
        //----- methods
        void Estimate3DSimil(std::vector<cPt3dr>& aVInPts, std::vector<cPt3dr>& aVOutPts, bool& aShow);
};

/** Class to store coded target detection **/

struct cDetCdT
{
    cDetCdT(const cSensorCamPC* aCam, cMesIm1Pt aMes, cAff2D_r aAff2D);
    const cSensorCamPC* mCam;
    cMesIm1Pt mMes;
    cAff2D_r mIm2Ref;
};
}
