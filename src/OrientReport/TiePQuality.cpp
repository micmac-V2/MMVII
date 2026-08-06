#include "MMVII_Ptxd.h"
#include "cMMVII_Appli.h"
#include "MMVII_Geom3D.h"
#include "MMVII_PCSens.h"
#include "MMVII_Tpl_Images.h"


/**
   \file GCPQuality.cpp


 */

namespace MMVII
{

/* *********************************************************** */
/*                                                             */
/*                       cImageVectorField                     */
/*                                                             */
/* *********************************************************** */



void cImageVectorField::Save(const std::string & aName)
{
  mIm.ToFileDeZoom(aName+".tif",mDeZoom);
}



tPtrArg2007  cImageVectorField::ArgOpt(cCollecSpecArg2007 & anArgOpt,std::vector<tREAL8>  & aParamsFV)
{

    return    AOpt2007
             (
                 aParamsFV,
                 "ParamFV",
                 "Parameters for field of vectors [Ampl,Width?,Ray?,Zoom?,Jpeg?]",
                 {{eTA2007::HDV},{eTA2007::ISizeV,"[1,5]"}}
             );
}


cImageVectorField::cImageVectorField(const std::string & aNameIm,const  std::vector<tREAL8> &aVParams) :
    mNameIm       (aNameIm),
    mIm           (cRGBImage::FromFile(mNameIm,1,eForceGray::Yes)),
    mAmpl         (aVParams.at(0)),
    mWidth        (GetDef(aVParams,1,DefParam().at(1))),
    mRay          (GetDef(aVParams,2,DefParam().at(2))),
    mDeZoom       (GetDef(aVParams,3,DefParam().at(3))),
    mJPeg         (GetDef(aVParams,4,DefParam().at(4))),
    mColorArrow   (cRGBImage::Red),
    mColorCircle  (cRGBImage::Blue)
{
}

const std::vector<tREAL8> & cImageVectorField::DefParam()
{
    static std::vector<tREAL8> aRes {100.0,1.0,4.0,1.0,1.0};
    return aRes;
}

cRGBImage cImageVectorField::Im() const {return mIm;}

void cImageVectorField::DrawArrow_P1Vect(const cPt2dr & aP1,const cPt2dr& aVect)
{
    mIm.DrawCircle(mColorCircle,aP1,mRay);
    mIm.DrawLine(aP1,aP1+aVect*mAmpl,mColorArrow,mWidth);
}

void cImageVectorField::DrawArrow_P1P2(const cPt2dr & aP1,const cPt2dr& aP2)
{
    DrawArrow_P1Vect(aP1,aP2-aP1);
}

void cImageVectorField::DrawPointsRemarq(const cPt2dr& aPt,const cPt3di & aCoul,const std::vector<tREAL8>& aVRay)
{
    for (const auto & aRay : aVRay)
        mIm.DrawCircle(aCoul,aPt,aRay);
}
void cImageVectorField::DrawPointsRemarq(const cPt2dr&aPt)
{
    DrawPointsRemarq(aPt,mColorCircle,{2.0,30.0,32.0,34.0});
}

//==============================================

struct cStatRes2D
{
  public :
    cStdStatRes mStatRes;
    cWeightAv<tREAL8,cPt2dr>  mAvg2d;

    void AddResidu(const cPt2dr &aRes)
    {
        mStatRes.Add(Norm2(aRes));
        mAvg2d.Add(1.0,aRes);
    }

};

class cAppli_TiePReport : public cMMVII_Appli
{
     public :

        cAppli_TiePReport(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli &);
        int Exe() override;
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override;
        std::vector<std::string>  Samples() const override;



     private :

            cPt3dr ColourOfResidual(tREAL8 aRes);

             //  void ShowImage(const cSetIm4SparseDist<tREAL8>&,const std::string& aVal) const;
             void ProcessBySingleImage();
             void ProcessByPair();
             void ProcessPly();

             void ProcessOnePair(const std::string&,const std::string&);

             cPhotogrammetricProject  mPhProj;
             std::vector<int>         mPropStat;
             std::string              mSpecImIn;   ///  Pattern of xml file
             std::vector<std::string> mSetNames;



             void RegisterStatRes(const std::string & aFile,const cStatRes2D&, const std::string anAggreg);
             void InitReporStat(const std::string & aFile);

             cComputeMergeMulTieP *   mCMTP;
             std::string              mPrefixCSVIma;
             std::string              mPrefixCSVPairs;
             std::string              mPrefixCSVCam;
             std::string              mPrefixCSVGlob;





             std::string              mPatImRes;       // Pattern for residual images
             std::string              mPatImFV;        // Patternof images for field vector
             tREAL8                   mFactRedIm;      // Reduction factor for residual images
             bool                     mDoImResCalib;   // do we generate residual images by calib

             std::string              mPatternVFByPair;
             std::string              mPatternResByPair;

             std::vector<tREAL8>      mParamPly;
             cPlyVertices             mPlyFile;

             std::vector<tREAL8>      mParamsFV;  /// Parameters for fields of vectors
             tREAL8                   mStepHisto;
             tREAL8                   mMaxHisto;
             tHistoCumulrr            mHistoRes;

             std::string              mPrefixSave;
};

void cAppli_TiePReport::RegisterStatRes(const std::string & aFile,const cStatRes2D& aS2D, const std::string anAggreg)
{
   cPt2dr anAvg2 = aS2D.mAvg2d.Average();
   AddStdStatCSV
   (
      aFile,anAggreg,aS2D.mStatRes,mPropStat,
      {ToStr(anAvg2.x()),ToStr(anAvg2.y())}
   );
}

void cAppli_TiePReport::InitReporStat(const std::string & aFile)
{
    InitReportCSV(aFile,"csv",false);
    AddStdHeaderStatCSV(aFile,"Image",mPropStat,{"AvgX","AvgY"});
}


cAppli_TiePReport::cAppli_TiePReport
(
     const std::vector<std::string> &  aVArgs,
     const cSpecMMVII_Appli & aSpec
) :
     cMMVII_Appli  (aVArgs,aSpec),
     mPhProj       (*this),
     mPropStat     ({50,75}),
     mCMTP         (nullptr),
     mFactRedIm    (100.0),
     mDoImResCalib (true),
     mParamsFV     (cImageVectorField::DefParam()),
     mStepHisto    (1e-3),
     mMaxHisto     (1e2),
     mHistoRes     (mMaxHisto/mStepHisto)
{
}

std::vector<std::string>  cAppli_TiePReport::Samples() const
{
    return {
        "MMVII ReportTieP .*JPG V1 Adjust  ParamPly=[1,4.0]",
        "MMVII ReportTieP IMG_03.*JPG V1Dense Adjust InTieP=V1Dense  PairPatRes=IMG_035[0-3].JPG  PairPatFV=IMG_035[1-2].JPG"
    }   ;
}

cCollecSpecArg2007 & cAppli_TiePReport::ArgObl(cCollecSpecArg2007 & anArgObl)
{
      return     anArgObl
              << Arg2007(mSpecImIn,"Pattern/file for images",{{eTA2007::MPatFile,"0"},{eTA2007::FileDirProj}})
              << mPhProj.DPMulTieP().ArgDirInMand()
              << mPhProj.DPOrient().ArgDirInMand()
      ;
}


cCollecSpecArg2007 & cAppli_TiePReport::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{

    return   anArgOpt

          << cHeaderSectionArg ("Parameters for per/image stas/visu")
          << AOpt2007(mPatImRes,"PatImRes","Pattern of  images  where we generat residual (if any)")
          << AOpt2007(mPatImFV,"PatImFV","Pattern of  images  where we generat field vector (if any)")
          << AOpt2007(mDoImResCalib,"DoImResCalib","Do we generate images of residual by calibration ?",{eTA2007::HDV})

          << cHeaderSectionArg("Parameters by pairs ")
          << mPhProj.DPTieP().ArgDirInOpt("","Tie point in pair format for by pair stat/visu")
          << AOpt2007(mPatternVFByPair,"PairPatFV","Pattern for pairs generatin vector fiels by pair")
          << AOpt2007(mPatternResByPair,"PairPatRes","Pattern for pairs generating residual images")




          << cHeaderSectionArg("General pameters (for per-images and pairs)")
          << AOpt2007(mPropStat,"Perc","Percentil for stat exp",{eTA2007::HDV})
          << cImageVectorField::ArgOpt(anArgOpt,mParamsFV)
             /*
          << AOpt2007
             (
                 mParamsFV,
                 "ParamFV",
                 "Parameters for field of vectors [Ampl,Width?,Ray?,Zoom?,Jpeg?]",
                 {{eTA2007::HDV},{eTA2007::ISizeV,"[1,5]"}}
             )*/

          << cHeaderSectionArg("3D pameters")
          << AOpt2007(mParamPly,"ParamPly","Generate a 3D visualization of ply files [Colored,ExpQual]",{{eTA2007::ISizeV,"[2,2]"}})
    ;
}




void  cAppli_TiePReport::ProcessOnePair(const std::string & aName1,const std::string& aName2)
{

       cSensorImage* aSens1 = mPhProj.ReadSensor(aName1,true);
       cSensorImage* aSens2 = mPhProj.ReadSensor(aName2,true);

    //   std::string aDirHom =  mParamByPair.at(0);

       cSetHomogCpleIm aSetH;

       bool HasHom = mPhProj.GenReadHomol(aSetH,aName1,aName2);
       if (!HasHom)
           return;

       cImageVectorField * aIVF = nullptr;
       if (IsInit(&mPatternVFByPair) && (MatchRegex(aName1,mPatternVFByPair)) && (MatchRegex(aName2,mPatternVFByPair)))
          aIVF = new cImageVectorField(aName1,mParamsFV);

       cSetIm4SparseDist<tREAL8>* aImAvg = nullptr;
       if (IsInit(&mPatternResByPair)  && (MatchRegex(aName1,mPatternResByPair)) && (MatchRegex(aName2,mPatternResByPair)))
       {
           aImAvg = new cSetIm4SparseDist<tREAL8>({"N2","x","y"},aSens1->Sz(),mFactRedIm);
       }


       cStatRes2D aStat;

       for (const auto & aCple : aSetH.SetH())
       {
           cPt3dr aPG = aSens1->PInterBundle(aCple,*aSens2);
           cPt2dr aVect = (aSens1->Ground2Image(aPG)-aCple.mP1) * 2.0;

           if (aIVF)
               aIVF->DrawArrow_P1Vect(aCple.mP1,aVect);
           if (aImAvg)
               aImAvg->Add(aCple.mP1,{Norm2(aVect),aVect.x(),aVect.y()});
          aStat.AddResidu(aVect);

       }

       RegisterStatRes(mPrefixCSVPairs,aStat,aName1+"/" + aName2);
       if ( aIVF)
       {
          aIVF->Save(mPhProj.DirVisuAppli() + "PairIVF-"+aName1+"-"+aName2 + ".tif");
          delete aIVF;
       }
       if (aImAvg)
       {
           aImAvg->SaveDenseSave(mPhProj.DirVisuAppli(),"PairRes-"+aName1+"-"+aName2+".tif");
           delete aImAvg;
       }

}


void  cAppli_TiePReport::ProcessByPair()
{
    InitReporStat(mPrefixCSVPairs);

    for (size_t aKIm1 = 0 ; aKIm1<mSetNames.size() ; aKIm1++)
    {       
        for (size_t aKIm2 = aKIm1+1 ; aKIm2<mSetNames.size() ; aKIm2++)
        {
           ProcessOnePair(mSetNames.at(aKIm1),mSetNames.at(aKIm2));
        }
    }
}

cPt3dr cAppli_TiePReport::ColourOfResidual(tREAL8 aRes)
{
   tREAL8 aProp=    mHistoRes.PropCumul(aRes/mStepHisto);

   tREAL8 aQual = 1.0 - std::pow(aProp,mParamPly.at(1));

   if (mParamPly.at(0)==0.0)
   {
      return cPt3dr::PCste(1.0 -aQual);
   }
   else
   {
      static  cPt3dr aColRed (1.0,0,0);
      static  cPt3dr aColOrange(1.0,0.5,0);
      static  cPt3dr aColGreen(0,1.0,0.0);

       if (aQual<0.5)
       {
           tREAL8 aWRed = 2.0 * (0.5-aQual);
           return aColRed * aWRed + aColOrange * (1-aWRed);
       }
       tREAL8 aWGreen = (aQual-0.5) * 2.0;
       return aColGreen * aWGreen + aColOrange * (1-aWGreen);
   }
}


void cAppli_TiePReport::ProcessPly()
{

    // Compute residual & histogra!!
    std::vector<tREAL8> aVRes;
    for (const auto&  aPair : mCMTP->Pts())
    {
        const auto & [aConfig,aValue]  = aPair;
        size_t aNbPts = NbPtsMul(aPair);
        int aNbIm = Multiplicity(aPair);
        for (size_t aKPt=0 ; aKPt<aNbPts ; aKPt++)
        {
            cPt3dr aPGround = aValue.mVPGround.at(aKPt);
            tREAL8 aRes = 0.0;
            for (int aKIm = 0 ; aKIm<aNbIm ; aKIm++)
            {
                 cSensorImage * aSensI = mCMTP->VSensors().at(aConfig.at(aKIm));
                cPt2dr aMesIm = KthPt(aPair,aKIm,aKPt);
                cPt2dr aProj =aSensI->Ground2Image(aPGround);
                aRes += Norm2(aMesIm-aProj);
            }
            aRes /= aNbIm-1;
            aVRes.push_back(aRes);

            mHistoRes.AddV(aRes/mStepHisto,1.0);
        }
    }
    // make cumulative histogramm
    mHistoRes.MakeCumul();

    // generate date for ply
    size_t aKRes=0;
    for (const auto&  aPair : mCMTP->Pts())
    {
        const auto & [aConfig,aValue]  = aPair;
        size_t aNbPts = NbPtsMul(aPair);
        for (size_t aKPt=0 ; aKPt<aNbPts ; aKPt++)
        {
            cPt3dr aPGround = aValue.mVPGround.at(aKPt);

            cPt3dr aCol = ColourOfResidual(aVRes.at(aKRes++));
            mPlyFile.AddVert(aPGround,aCol);
        }
    }

    // export
    mPlyFile.ToPly(mPhProj.DirVisuAppli()+"RTP.ply",false);


}

void cAppli_TiePReport::ProcessBySingleImage()
{
   std::map<std::string,cStatRes2D>                 aMapStatByCam;
   cStatRes2D                                       aStatGlob;
   std::map<std::string,cSetIm4SparseDist<tREAL8>*>  aMapResByCam;

   InitReporStat(mPrefixCSVIma);
   InitReporStat(mPrefixCSVCam);
   InitReporStat(mPrefixCSVGlob);

   /// InitReportCSV(mPrefixCSVIma,"csv",false);
   /// AddStdHeaderStatCSV(mPrefixCSVIma,"Image",mPropStat,{"AvgX","AvgY"});

   // Parse all images
   for (size_t aKImGlob = 0 ; aKImGlob<mSetNames.size() ; aKImGlob++)
   {
       cSensorImage * aSensI = mCMTP->VSensors().at(aKImGlob);
       std::string aNameIm = aSensI->NameImage();


     //  std::string aNameCal = mCMTP->VSensors()

       // a litle check on indexe for these complexe structures
       MMVII_INTERNAL_ASSERT_tiny(mSetNames[aKImGlob]==aNameIm,"Chek names in Appli_TiePReport::MakeStatByImage");

       cStatRes2D  aStat2;


       // ---  Data for images of residual per images--------------------------
       cSetIm4SparseDist<tREAL8>* aImAvg = nullptr;
       if (IsInit(&mPatImRes) && MatchRegex(aNameIm,mPatImRes))
       {
           aImAvg = new cSetIm4SparseDist<tREAL8>({"N2","x","y"},aSensI->Sz(),mFactRedIm);
       }

       //  ---------- Data for image of vector's field ---------------
       cImageVectorField * aIVF = nullptr;
       if (IsInit(&mPatImFV) &&  MatchRegex(aNameIm,mPatImFV))
       {
          aIVF = new cImageVectorField(aNameIm,mParamsFV);
       }

      // StdOut() << "PATFV " << aNameIm << " " << aIVF << "\n";

       // ---------  Data for images of residual per camera ------------------
       std::string aNameCamera;
       cSetIm4SparseDist<tREAL8>* aImAvgCamera = nullptr;
       if (aSensI->IsSensorCamPC() && mDoImResCalib)
       {
           aNameCamera = aSensI->GetSensorCamPC()->InternalCalib()->Name();
           if (! MapBoolFind(aMapResByCam,aNameCamera))
               aMapResByCam[aNameCamera] = new cSetIm4SparseDist<tREAL8>({"N2","x","y"},aSensI->Sz(),mFactRedIm);
           aImAvgCamera =   aMapResByCam[aNameCamera];
       }

       // Parse all the configuration that contain the image
       const auto & aLInd = mCMTP->IndexeOf1Image(aKImGlob);
       for (const auto& [aKImLoc,aPairConfTiep] : aLInd)
       {
           // little check that aKImLoc is the postisition where aKImGlob is found in the configuration of tie points
           {
              // Unused in mode release
              [[maybe_unused]] const auto & aConfig = aPairConfTiep->first;
              // a litle check on indexe for these complexe structures
              MMVII_INTERNAL_ASSERT_tiny(aKImGlob==(size_t)aConfig.at(aKImLoc),"Check nums in cAppli_TiePReport::MakeStatByImage");
           }

           const cVal1ConfTPM & aVal = aPairConfTiep->second; // The data itsefl : tiep, Pt3D ? ,
           size_t aMult = Multiplicity(*aPairConfTiep);
           tREAL8 aAmplMul = aMult / tREAL8(aMult); // take into account degree of freedom

           // Parse all the points of this configuration that contain the image
           for (size_t aKp=0 ; aKp<NbPtsMul(*aPairConfTiep) ; aKp++)
           {
               cPt2dr aPtIm =  KthPt(*aPairConfTiep,aKImLoc,aKp);
               cPt3dr aPGr = aVal.mVPGround.at(aKp);
               cPt2dr aPProj = aSensI->Ground2Image(aPGr);
               cPt2dr aResidu = (aPProj - aPtIm) * aAmplMul;

               aStat2.AddResidu(aResidu);
               aStatGlob.AddResidu(aResidu);
               if (aNameCamera!="")
                  aMapStatByCam[aNameCamera].AddResidu(aResidu);


               if (aImAvg)
               {
                   aImAvg->Add(aPtIm,{Norm2(aResidu),aResidu.x(),aResidu.y()});
               }
               if (aImAvgCamera)
               {
                   aImAvgCamera->Add(aPtIm,{Norm2(aResidu),aResidu.x(),aResidu.y()});
               }
               if (aIVF)
               {
                   aIVF->DrawArrow_P1Vect(aPtIm,aResidu);
               }

           }
       }
       // Eventually write and clear  Images of resisual
       if (aImAvg)
       {
           aImAvg->SaveDenseSave(mPhProj.DirVisuAppli(),aNameIm);
           delete aImAvg;
       }
       // Eventutaly, write and clear images of vector's field
       if (aIVF)
       {
           StdOut()  << " SAVEFV "  << (mPhProj.DirVisuAppli()+"ResVF_"+LastPrefix(aNameIm)) << "\n";
            aIVF->Save(mPhProj.DirVisuAppli()+"ResVF_"+LastPrefix(aNameIm));
           delete aIVF;
       }

       // Save statistique for this image
       RegisterStatRes(mPrefixCSVIma,aStat2,aNameIm);
   }

   for (auto [aNameCam,aStat] : aMapStatByCam)
      RegisterStatRes(mPrefixCSVCam,aStat,aNameCam);

    RegisterStatRes(mPrefixCSVGlob,aStatGlob,"Glob");

    for (auto [aName,aAvgByCam] : aMapResByCam)
    {
        aAvgByCam->SaveDenseSave(mPhProj.DirVisuAppli(),aName);
        delete aAvgByCam;
    }
}


int cAppli_TiePReport::Exe()
{
   mPhProj.FinishInit();

   mPrefixSave = mPhProj.DPOrient().DirIn() + "_" + mPhProj.DPMulTieP().DirIn();
   SetReportSubDir(mPrefixSave);
   mPhProj.SetVisuSubDir(mPrefixSave);


   mSetNames = VectMainSet(0);
   // ...    bool  WithIndexPt,   bool  WithSensor,bool  WithIndexImages
   mCMTP = AllocStdFromMTP(mSetNames,mPhProj,true,true,true);
   mCMTP->SetPGround();

   mPrefixCSVIma =  "ByImages" ;
   mPrefixCSVPairs = "ByPairs";
   mPrefixCSVCam   = "ByCams";
   mPrefixCSVGlob   = "Glob";



   ProcessBySingleImage();

   if ( mPhProj.DPTieP().DirInIsInit())
       ProcessByPair();
   else
   {
       MMVII_INTERNAL_ASSERT_User_UndefE(!IsInit(&mPatternVFByPair),"Pattern field vector without TieP By Pair");
       MMVII_INTERNAL_ASSERT_User_UndefE(!IsInit(&mPatternResByPair),"Pattern residual without TieP By Pair");

     //  << mPhProj.DPTieP().ArgDirInOpt("","Tie point in pair format for by pair stat/visu")
     //  << AOpt2007(mPatternVFByPair,"PairPatFV","Pattern for pairs generatin vector fiels by pair")
     //  << AOpt2007(mPatternResByPair,"PairPatRes","Pattern for pairs generating residual images")
   }

   if (IsInit(&mParamPly))
       ProcessPly();





   delete mCMTP;
   return EXIT_SUCCESS;
}


/* ==================================================== */

tMMVII_UnikPApli Alloc_TiePReport(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec)
{
   return tMMVII_UnikPApli(new cAppli_TiePReport(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpec_TiePReport
(
     "ReportTieP",
      Alloc_TiePReport,
      "Reports on TieP projection",
      {eApF::TieP,eApF::Ori},
      {eApDT::TieP,eApDT::Orient},
      {eApDT::Image,eApDT::Xml},
      __FILE__
);



}; // MMVII

