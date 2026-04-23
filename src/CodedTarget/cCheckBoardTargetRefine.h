#ifndef CCHECKBOARDTARGETREFINE_H
#define CCHECKBOARDTARGETREFINE_H

#endif // CCHECKBOARDTARGETREFINE_H

#include "CodedTarget.h"
#include "MMVII_PCSens.h"
#include "cMMVII_Appli.h"
#include "MMVII_Geom3D.h"
#include "cCodedTargetDescribe.h"

namespace MMVII
{
    class cCdTSampl;

    class cAppli_CheckBoardTargetRefine : public cMMVII_Appli
    {
    public:
        //------
        cAppli_CheckBoardTargetRefine(const std::vector<std::string>& aVArgs,
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
        //------ methods
    };

    /*!
     * @brief The cCdTDiscr class stores all about coded target image discretisation
     */

    class cCdTDiscr
    {
        public:
        cCdTDiscr(cCdTDescr aDes, cSensorCamPC* aCam);
            //----- members
            std::string mName;

            //----- methods
        private:
            //----- members
            std::string mNum;
            cIm2D<tU_INT1> mTrue;
            cIm2D<tU_INT1> mRef;
            cSensorCamPC* mCam;

            //----- methods
            void getSampledIm();
            void getTrueIm();
            void True2Sample();
    };
}
