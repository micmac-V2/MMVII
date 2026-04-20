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
        void AddDesCdT(const cOneEncoding* aEnc);
        bool DesCdtAdd(std::string& aName);
};

/** Class to store coded target refined description **/

class cDesCdT
{
    public:
        cDesCdT(const cOneEncoding* aEnc);
        std::string mName;
        std::vector<cDetCdT> mVDetects;
        void AddDetect(const cSensorCamPC* aCam, cMesIm1Pt aMes, cAff2D_r aAff2D);
        void InterBitCenters();
        void InterCorners();
        void Estimate3DSimil();
    private:
        const cOneEncoding* mEnc;
        tREAL8 mRes;
        std::vector<cPt3dr> mVBitCenters3D;
        std::vector<cPt3dr> mVCorners3D;
};

/** Class to store coded target detection **/

struct cDetCdT
{
    cDetCdT(const cSensorCamPC* aCam, cMesIm1Pt aMes, cAff2D_r aAff2D);
    const cSensorCamPC* mCam;
    cMesIm1Pt mMes;
    cAff2D_r mAff2D;
};
}
