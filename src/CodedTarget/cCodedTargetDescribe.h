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

    typedef cSegment<tREAL8,3> tSeg3dr;

    class cCdTDescr;
    struct cCdTDetec;

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
            //------ members
            std::string mFSpecName;
            std::unique_ptr<cFullSpecifTarget> mFSpec;
            std::vector<cCdTDescr> mVDescr;
            //------ methods
            void AddDescr(std::string aName, std::unique_ptr<cFullSpecifTarget>& aSpec);
    };

    /*!
     * @brief The cCdTDescr class is used to store and compute CoDed Target DEScription
     */

    class cCdTDescr
    {
        public:
        cCdTDescr(std::string aName, std::unique_ptr<cFullSpecifTarget>& aSpec);
            //----- members
            std::string mName;
            //----- methods
            void AddDetect(const cSensorCamPC* aCam, cMesIm1Pt aMes, cAff2D_r aAff2D);
            void InterBitCenters(bool& aShow);
            void InterGndCorners(bool& aShow);
            void EstimateSimil3DOnCorners(bool& aShow);
            void AddData(const cAuxAr2007 &anAux);
            static std::string NameFile(const cPhotogrammetricProject &aPhProj, std::string aName, bool Input);
        private:
            //----- members
            std::vector<cCdTDetec> mVDetec;
            const cOneEncoding* mEnc;
            tREAL8 mRes;
            std::vector<cPt2dr> mVBitCenters2D;
            std::vector<cPt2dr> mVCdTCorners;
            std::vector<cPt3dr> mVCdtCorners3D;
            std::vector<cPt3dr> mVGndCorners;
            cSimilitud3D<tREAL8> mSimil3D;
            //----- methods
            cPt2dr Gnd2CdT(cPt3dr& aPt, const cCdTDetec& aDet);
            cPt3dr CdT2GndByInter(const cPt2dr& aPt, std::vector<tREAL8>* aVRes = nullptr);
            cPt3dr CdT2GndBySimil();
            void Estimate3DSimil(std::vector<cPt3dr>& aVInPts, std::vector<cPt3dr>& aVOutPts, bool& aShow);
    };
    void AddData(const cAuxAr2007 &anAux, cSetTargetMap &anEx);

    /*!
     * @brief The cCdTDetec class is used to store CoDed Target DETection
     */

    struct cCdTDetec
    {
        cCdTDetec(const cSensorCamPC* aCam, cMesIm1Pt aMes, cAff2D_r aAff2D);
        const cSensorCamPC* mCam;
        cMesIm1Pt mMes;
        cAff2D_r mIm2Ref;
    };
}
