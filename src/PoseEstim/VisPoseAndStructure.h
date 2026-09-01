#ifndef  _MMVII_VISPOSEANDSTRUCTURE_H_
#define  _MMVII_VISPOSEANDSTRUCTURE_H_

#include "cMMVII_Appli.h"
#include "MMVII_Sensor.h"
#include "MMVII_UtiSort.h"
#include "treethread.h"
#include "MMVII_Tpl_ElemStrToVal.h"

//#include "MMVII_util_tpl.h"
//#include "MMVII_Geom3D.h"

namespace MMVII
{

/* ********************************************************** */
/*                                                            */
/*                     cAppli_VisuPoseStr3D                   */
/*                                                            */
/* ********************************************************** */

class cStaticLidar;


struct cEstimLengthPyrCam
{
    double mProp;
    double mMult;

    cEstimLengthPyrCam () :
        mProp  (0.3),
        mMult  (0.1)
    {
    }
    ARG2007_STRUCT_FIELDS (
        mProp,FieldSem({{eTA2007::AddCom,"Prop of estimation of depth"},eTA2007::HDV}),
        mMult,FieldSem({{eTA2007::AddCom,"Mult depth->lengtght pyram of cam"},eTA2007::HDV})
     //   OutDir,FieldSem({eTA2007::HDV,{eTA2007::AddCom,"Optional output dir"}}),
     //   ExportSigma,FieldSem({eTA2007::HDV,{eTA2007::AddCom,"Optional compensated sigma"}})
    )
};

class cAppli_VisuPoseStr3D : public cMMVII_Appli
{
public:
    cAppli_VisuPoseStr3D(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec);

    int Exe() override;
    cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override ;
    cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override ;

    std::vector<cOneHelpSampleCmp>  Samples() const override;


private:

    cPhotogrammetricProject   mPhProj;
    std::string               mPatImIn;
    double                    mErrProjMax;
    double                    mCamScale;
    cEstimLengthPyrCam        mEstimmLengthPyrCam;
    int                       mTSLCloudDezoom; //  0 for no point clouds
    std::string               mOutfile;
    bool                      mBinary;
    bool                      mWithRGB;
    bool                      mWithAvgRGB; //RGB values averaged over all images, slower? but less noisy

    void AddCameras(cPlyVertices& aPlyverts, cComputeMergeMulTieP * &, const std::vector<cSensorImage *>& );
    void AddPointCould(cPlyVertices& aPlyverts, cStaticLidar* aScan, int aTSLCloudDezoom, cPt3dr aScanColor, bool aOnlyEdges=false);
    double CalculateFDepth(const cPt2di&, const double&);
};

}

#endif // _MMVII_VISPOSEANDSTRUCTURE_H_
