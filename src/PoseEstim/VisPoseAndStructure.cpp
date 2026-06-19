#include "VisPoseAndStructure.h"
#include "MMVII_StaticLidar.h"

/**
   \file VisPoseAndStructure.cpp

   \brief Visualise a set of poses and the sparse 3D structure

*/

namespace MMVII
{




cAppli_VisuPoseStr3D::cAppli_VisuPoseStr3D(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec) :
    cMMVII_Appli (aVArgs,aSpec),
    mPhProj      (*this),
    mErrProjMax  (10.0),
    mCamScale    (0.1),
    mTSLCloudDezoom(0),
    mOutfile     ("VisSFM_${ori}_${features}.ply"),
    mBinary      (true),
    mWithRGB     (true),
    mWithAvgRGB  (false)
{
}

cCollecSpecArg2007 & cAppli_VisuPoseStr3D::ArgObl(cCollecSpecArg2007 & anArgObl)
{
    return anArgObl
           << Arg2007(mPatImIn,"Pattern/file of images",{{eTA2007::MPatFile,"0"},{eTA2007::FileDirProj}})
           << mPhProj.DPOrient().ArgDirInMand("Input orientation (Sfm)")
        ;
}

cCollecSpecArg2007 & cAppli_VisuPoseStr3D::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
    return    anArgOpt
           << AOpt2007(mErrProjMax,"ErrMax","Outlier threshold",{eTA2007::HDV})
           << AOpt2007(mCamScale,"CamScale","Scale camera frustum",{eTA2007::HDV})
           << AOpt2007(mTSLCloudDezoom,"TSLcloudDezoom","Add TSL point clouds with DeZoom (0 for no points)",{eTA2007::HDV})
           << AOpt2007(mOutfile,"Outfile","Output filename",{eTA2007::HDV})
           << AOpt2007(mBinary,"Bin","Output in binary format",{eTA2007::HDV})
           << AOpt2007(mWithRGB,"RGB","Output colored pointcloud",{eTA2007::HDV})
           << AOpt2007(mWithAvgRGB,"AvgRGB","RGB values averaged over all images",{eTA2007::Tuning,eTA2007::HDV})
           << mPhProj.DPMulTieP().ArgDirInOpt("","Input features (multiple tie-points format)")
           << mPhProj.DPGndPt2D().ArgDirInOpt("","Input features (image measurements format)")
        ;
}

int cAppli_VisuPoseStr3D::Exe()
{
    mPhProj.FinishInit();

    if (!IsInit(&mOutfile))
        mOutfile =  mPhProj.DirVisuAppli() +
                     ("VisSFM_"+mPhProj.DPOrient().DirIn()+
                    (mPhProj.DPMulTieP().DirInIsInit() ? "_"+mPhProj.DPMulTieP().DirIn() :
                     mPhProj.DPGndPt2D().DirInIsInit() ? "_"+mPhProj.DPGndPt2D().DirIn() : "")  +".ply");


    // vector of all image names
    std::vector<std::string> aVNames;
    // cameras
    std::vector<cSensorImage *> aVSens ;

    // read names and cameras; skip images without an orientation
    cAutoTimerSegm  aATS(TimeSegm(),"READ CAMERAS");
    for (const auto & aIm : VectMainSet(0))
    {
        cSensorImage* aSens = mPhProj.ReadSensor(aIm,true,true);
        if (aSens)
        {
            aVNames.push_back(aIm);
            aVSens.push_back(aSens);
        }
        else StdOut() << "Image " << aIm << " has no orientation" << std::endl;
    }
    // sort images alphbetically (and aVSens accordingly) for AllocStdFromMTPFromFolder
    Sort2VectFirstOne(aVNames,aVSens);

    // read the tie points
    TimeSegm().SetIndex("READ TIE-POINTS");
    cComputeMergeMulTieP * aTPts = nullptr;
    if (mPhProj.DPMulTieP().DirInIsInit())
    {
        aTPts = AllocStdFromMTPFromFolder(
                mPhProj.DPMulTieP().DirIn(),aVNames,mPhProj,true,false,true);
    }
    else if (mPhProj.DPGndPt2D().DirInIsInit())
    {
        // aVNames is already sorted (Sort2VectFirstOne above)
        cMemoryInterfImportHom aMIIH;
        for (size_t aK1=0; aK1<aVNames.size(); aK1++)
        {
            cSetMesPtOf1Im * aSetM1 = mPhProj.RemanentLoadMeasureIm(aVNames[aK1]);
            if (!aSetM1) continue;
            for (size_t aK2=aK1+1; aK2<aVNames.size(); aK2++)
            {
                cSetMesPtOf1Im * aSetM2 = mPhProj.RemanentLoadMeasureIm(aVNames[aK2]);
                if (!aSetM2) continue;
                cSetHomogCpleIm aCpleH;
                aCpleH.AddPairSet(*aSetM1,*aSetM2);
                if (aCpleH.NbH() > 0)
                    aMIIH.Add(aCpleH,aVNames[aK1],aVNames[aK2]);
            }
        }
        aTPts = new cComputeMergeMulTieP(aVNames,&aMIIH);
    }

    // intersect in 3d
    TimeSegm().SetIndex("INTERSECT 3D TIE-POINTS");
    if (aTPts)
        for (auto & aPair : aTPts->Pts())
            MakePGround(aPair,aVSens);


    std::string aPrintRGB = (mWithRGB) ? "+RGB": "";
    std::string aPrintPT3D = "ADD 3D PTS";
    std::string aPrintCAM = "+CAMERAS";
    std::string aPrintPtsRGBCam = aPrintPT3D + aPrintRGB + aPrintCAM;
    TimeSegm().SetIndex(aPrintPtsRGBCam);
    cPlyVertices aPlyverts;

    AddCameras(aPlyverts,aTPts,aVSens);

    if (mTSLCloudDezoom>0)
        for (auto &aSens : aVSens)
        {
            cStaticLidar* aScan = dynamic_cast<cStaticLidar*>(aSens);
            if (aScan)
            {
                aScan->ReadRasters(mPhProj.DirStaticLidarRasters());
                cPt3dr aColor(RandUnif_0_1(),RandUnif_0_1(), RandUnif_0_1());
                AddPointCould(aPlyverts, aScan, mTSLCloudDezoom, aColor);
            }
        }

    aPlyverts.ToPly(mOutfile, mBinary);

    delete aTPts;

    return EXIT_SUCCESS;
}


tMMVII_UnikPApli Alloc_VisuPoseStr3D(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec)
{
    return tMMVII_UnikPApli(new cAppli_VisuPoseStr3D(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpec_VisuPoseStr3D
    (
        "VisuPose3D",
        Alloc_VisuPoseStr3D,
        "Create PLY with poses and 3D structure",
        {eApF::Ori},
        {eApDT::Ori},
        {eApDT::Orient},
        __FILE__
        );

}; // namespace MMVII
