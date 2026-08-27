#ifndef  _MMVII_VISPOSEANDSTRUCTURE_H_
#define  _MMVII_VISPOSEANDSTRUCTURE_H_

#include "cMMVII_Appli.h"
#include "MMVII_Sensor.h"
#include "MMVII_UtiSort.h"
#include "treethread.h"
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
