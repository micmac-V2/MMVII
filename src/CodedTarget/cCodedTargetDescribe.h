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


class TargetDescript;

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
        std::vector<TargetDescript> mVTargetsDescriptions;
        //------ methods
        void LoadImMes();//→ init & fill mVTargetsDecriptions
        void ComputeBaseBundles();
};

class TargetDescript
{
    public:
        TargetDescript(const std::string& aName);
    private:
        std::string mName;
        cTargetMap mFieldMap;
        std::vector<cMesIm1Pt> mVMesIm;
        std::vector<std::string> mVBundles;
        void AddBundle(tSeg3dr aBundle);
};

}
