#include "MMVII_Ptxd.h"
#include "cMMVII_Appli.h"
#include "MMVII_Geom3D.h"
#include "MMVII_PCSens.h"
#include "MMVII_Tpl_Images.h"
#include "MMVII_StaticLidar.h"

/**
   \file GCPCompare.cpp


 */

namespace MMVII
{

/*cSet2D3D

cPerspCamIntrCalib*/
//cSet2D3D

/* ==================================================== */
/*                                                      */
/*                   cAppli_CmpCalib                    */
/*                                                      */
/* ==================================================== */

class cAppli_CmpCalib : public cMMVII_Appli
{
    public :

        cAppli_CmpCalib(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli &,bool isModeLocal);
        int Exe() override;
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override;

        std::vector<cOneHelpSampleCmp> Samples() const override;

    private :

        void CmpCalib(const cPerspCamIntrCalib & aCal1,const cPerspCamIntrCalib & aCal2);

        void Add1Pt(const cPt2dr &);


        void OneRes(bool Adj);

        cPhotogrammetricProject  mPhProj;
        bool                     mIsModeLocal;

        std::string              mNameIm1;
        std::string              mOri1;

        std::string              mNameIm2;
        std::string              mOri2;

        int                      mNbPts;
        cPt2di mSzCal;
        const cPerspCamIntrCalib * mCal1;
        const cPerspCamIntrCalib * mCal2;
        cRotation3D<tREAL8>        mRot2to1;

        std::vector<cPt2dr>   mPtsIm;
        std::vector<cPt3dr>   mDirsB1;
        std::vector<cPt3dr>   mDirsB2;
        std::vector<cPt3dr>   mDirsB2Adj;

};

cAppli_CmpCalib::cAppli_CmpCalib
(
    const std::vector<std::string> &  aVArgs,
    const cSpecMMVII_Appli & aSpec,
    bool  isModeLocal
) :
    cMMVII_Appli  (aVArgs,aSpec),
    mPhProj       (*this),
    mIsModeLocal  (isModeLocal),
    mNameIm1      (MMVII_NONE),
    mOri1         (MMVII_NONE),
    mNameIm2      (MMVII_NONE),
    mOri2         (MMVII_NONE),
    mNbPts        (5000)
{
}


cCollecSpecArg2007 & cAppli_CmpCalib::ArgObl(cCollecSpecArg2007 & anArgObl)
{
    if (mIsModeLocal)
    {
        anArgObl
                << Arg2007(mNameIm1,"Name of first image",{{eTA2007::FileImage},{eTA2007::FileDirProj}})
                << mPhProj.DPOrient().ArgDirInMand()
        ;
    }
    else
    {

    }

    return anArgObl;
}

cCollecSpecArg2007 & cAppli_CmpCalib::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
    if (mIsModeLocal)
    {
        anArgOpt
              << AOpt2007(mOri2,"Ori2","Second orientation folder when != first",{{eTA2007::Input},{eTA2007::Orient}})
              << AOpt2007(mNameIm2,"Im2","Second image when != first ",{{eTA2007::Input},{eTA2007::FileImage}})
            ;
     }
    else
    {
    }

    anArgOpt
             << AOpt2007(mNbPts,"NbPts","Number of points on the grid",{{eTA2007::HDV}})
     ;


     return anArgOpt;
}

std::vector<cOneHelpSampleCmp>  cAppli_CmpCalib::Samples() const
{
    return {
            // {"MMVII ReportGCPCmp Init Final Filter='.*_'"}
    };
}


void cAppli_CmpCalib::Add1Pt(const cPt2dr & aPIm)
{
    if (
             (mCal1->DegreeVisibilityOnImFrame(aPIm)<=1.0)
          || (mCal2->DegreeVisibilityOnImFrame(aPIm)<=1.0)
     )
    {
        return;
    }

   cPt3dr aDir1=  VUnit(mCal1->DirBundle(aPIm)) ;
   cPt3dr aDir2=  VUnit(mCal2->DirBundle(aPIm)) ;

   mPtsIm.push_back(aPIm);
   mDirsB1.push_back(aDir1);
   mDirsB2.push_back(aDir2);
}


void cAppli_CmpCalib::CmpCalib(const cPerspCamIntrCalib & aCal1,const cPerspCamIntrCalib & aCal2)
{
    mCal1 = & aCal1;
    mCal2 = & aCal2;
    mSzCal = mCal1->SzPix();

    MMVII_INTERNAL_ASSERT_User_UndefE
    (
         mSzCal == mCal2->SzPix(),
         "Different number of pixels"
    );

    std::vector<cPt2dr>  aVecPtsGrid = RegularGrid(mSzCal,mNbPts);

    for (const auto & aPtIm : aVecPtsGrid )
        Add1Pt(aPtIm);



    StdOut()  << " O1=" << mOri1 << " F1=" << aCal1.F()
              << " O2=" << mOri2 << " F2=" << aCal2.F()
               <<  " Nb=" << aVecPtsGrid.size() << " => " << mPtsIm.size()
              << "\n";

    mRot2to1 =  cRotation3D<tREAL8>::StdGlobEstimate(mDirsB2,mDirsB1,nullptr, nullptr,cParamCtrlOpt::Default());

    for (const auto aB2 : mDirsB2)
        mDirsB2Adj.push_back(mRot2to1.Value(aB2));

    OneRes(false);
    OneRes(true);

}

void cAppli_CmpCalib::OneRes(bool Adj)
{
    std::vector<cPt3dr>& aDirB2 = Adj ? mDirsB2Adj : mDirsB2;

    cWeightAv<tREAL8,tREAL8> aWRes;

    for (size_t aK=0 ; aK<mPtsIm.size() ; aK++)
    {
        cPt3dr aB1 = mDirsB1.at(aK);
        cPt3dr aB2 = aDirB2.at(aK);
        aWRes.Add(1.0,Norm2(aB1-aB2));
    }

    tREAL8 aRes = aWRes.Average();
    tREAL8 aResPix =aRes * mCal1->F();
    tREAL8 aDistPts =  std::sqrt(MulCoord(mCal1->SzPix()) / mNbPts) ;

    tREAL8 anAmpl = (aDistPts  / aResPix) / 3.0;


    std::string  aNameRes= (Adj ?  "Adjusted" : "Initial");
    StdOut() << Color::title  << aNameRes << "\n" << Color::end ;
    StdOut() << Color::sub_title  << "  Residuals" << "\n" << Color::end ;
    StdOut() <<   "   " << aRes    <<  Color::descr    <<  " Radian \n" << Color::end ;
    StdOut() <<   "   " << aResPix  <<  Color::descr   <<  " Pixels \n" << Color::end ;
    StdOut() <<  Color::sub_title  << "     Amplitude visu" <<Color::end  << anAmpl  << "  \n" << Color::end ;


    cRGBImage anImage(mSzCal,cRGBImage::Black);


    for (size_t aK=0 ; aK<mPtsIm.size() ; aK++)
    {
        cPt3dr aB1 = mDirsB1.at(aK);
        cPt3dr aB2 = aDirB2.at(aK);

        cPt2dr aPIm1 = mCal1->Value(aB1);
        cPt2dr aPIm2 = mCal1->Value(aB2);

        anImage.DrawCircle(cRGBImage::Green,aPIm1,aDistPts/20.0);
        anImage.DrawLine(aPIm1,aPIm1+(aPIm2-aPIm1)*anAmpl,cRGBImage::Red,3.0);

    }

    anImage.ToFileDeZoom(mPhProj.DirVisuAppli()+"Res"+aNameRes+".jpg",3);
}


int cAppli_CmpCalib::Exe()
{
   mPhProj.FinishInit();


  if (mIsModeLocal)
  {
      mOri1 = mPhProj.DPOrient().DirIn();
      if (!IsInit(&mOri2))
          mOri2 = mOri1;

      if (!IsInit(&mNameIm2))
          mNameIm2 = mNameIm1;

      cPerspCamIntrCalib * aCal1 = mPhProj.InternalCalibFromStdName(mNameIm1);
      cPerspCamIntrCalib * aCal2 = mPhProj.InternalCalibFromFolderStdName(mOri2,mNameIm2);


      CmpCalib(*aCal1,*aCal2);


  }
  else
  {
  }

   return EXIT_SUCCESS;
}

/* ==================================================== */
/*                                                      */
/*               MMVII                                  */
/*                                                      */
/* ==================================================== */



tMMVII_UnikPApli Alloc_CmpCalib_Local(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec)
{
   return tMMVII_UnikPApli(new cAppli_CmpCalib(aVArgs,aSpec,true));
}

cSpecMMVII_Appli  TheSpec_CmpCalib_Local
(
     "OriCmpCalibLocal",
      Alloc_CmpCalib_Local,
      "Compare internal calibrations, coming from local folder",
      {eApF::GCP},
      {eApDT::ObjCoordWorld},
      {eApDT::Xml},
      __FILE__
);



}; // MMVII

