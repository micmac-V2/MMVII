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

class cAppli_TiePReport : public cMMVII_Appli
{
     public :

        cAppli_TiePReport(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli &);
        int Exe() override;
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override;

     private :

             //  void ShowImage(const cSetIm4SparseDist<tREAL8>&,const std::string& aVal) const;
             void MakeStatByImage();
             void ProcessByPair();
             void ProcessPly();

             void ProcessOnePair(const std::string&,const std::string&);

             cPhotogrammetricProject  mPhProj;
             std::vector<int>         mPropStat;
             std::string              mSpecImIn;   ///  Pattern of xml file
             std::vector<std::string> mSetNames;
             std::map<std::string,cSetIm4SparseDist<tREAL8>*>  mMapResByCam;


             cComputeMergeMulTieP *   mCMTP;
             std::string              mPrefixCSV;
             std::string              mPrefixCSVIma;

             std::string              mPatImRes;       // Pattern for residual images
             tREAL8                   mFactRedIm;      // Reduction factor for residual images
             bool                     mDoImResCalib;   // do we generate residual images by calib

             std::vector<std::string> mParamByPair;

             bool                     mWithPly;
             cPlyVertices             mPlyFile;


};

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
     mWithPly      (false)
{
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
          << AOpt2007(mPropStat,"Perc","Percentil for stat exp",{eTA2007::HDV})
          << AOpt2007(mPatImRes,"PatImRes","Pattern of  images  where we generat residual (if any)")
          << AOpt2007(mDoImResCalib,"DoImResCalib","Do we generate images of residual by calibration ?",{eTA2007::HDV})
          << AOpt2007(mParamByPair,"ByPair","Parameter for detail by pair [TieP,Pat,Ampl?=100,SsRes?=2] ",{{eTA2007::ISizeV,"[2,4]"}})
          << AOpt2007(mWithPly,"WithPly","Generate a 3D visualization of ply files",{eTA2007::HDV})
    ;
}


class cImageVectorField
{
   public :
     cImageVectorField(const std::string & aNameIm,const std::vector<tREAL8> &aVParams);

     void DrawArrow(const cPt2dr & aP1,const cPt2dr& aP2);
   private :
     std::string   mNameIm;
     cRGBImage     mIm;
     tREAL8        mAmpl;
     tREAL8        mWidth;
     tREAL8        mRay;
     cPt3di        mColorArrow;
     cPt3di        mColorCircle;
};

cImageVectorField::cImageVectorField(const std::string & aNameIm,const  std::vector<tREAL8> &aVParams) :
    mNameIm       (aNameIm),
    mIm           (cRGBImage::FromFile(mNameIm,1,eForceGray::Yes)),
    mAmpl         (aVParams.at(0)),
    mWidth        (GetDef(aVParams,1,1.0)),
    mRay          (GetDef(aVParams,2,4.0)),
    mColorArrow   (cRGBImage::Red),
    mColorCircle  (cRGBImage::Blue)
{
}

void cImageVectorField::DrawArrow(const cPt2dr & aP1,const cPt2dr& aP2)
{
    mIm.DrawCircle(mColorCircle,aP1,mRay);
    mIm.DrawLine(aP1,aP1+(aP2-aP1)*mAmpl,mColorArrow,mWidth);
}




void  cAppli_TiePReport::ProcessOnePair(const std::string & aName1,const std::string& aName2)
{
       cSensorImage* aSens1 = mPhProj.ReadSensor(aName1,true);
       cSensorImage* aSens2 = mPhProj.ReadSensor(aName2,true);

       std::string aDirHom =  mParamByPair.at(0);
       tREAL8 anAmpl = cStrIO<tREAL8>::FromStr(GetDef(mParamByPair,2,std::string("100.0")));
       int aDeZoom = cStrIO<int>::FromStr(GetDef(mParamByPair,3,std::string("2")));


       cSetHomogCpleIm aSetH;

       bool HasHom = mPhProj.GenReadHomol(aSetH,aName1,aName2,aDirHom);
       if (!HasHom)
           return;

       tNameSelector aSelIm = AllocRegex(mParamByPair.at(1));

       if (! (aSelIm.Match(aName1) && aSelIm.Match(aName2)))
          return;

       StdOut()  << HasHom << " "<< aName1 << aSens1->Sz()  << " ;; " << aName2 << aSens2->Sz()   << " A=" << anAmpl << " Z="<< aDeZoom << "\n";

        cRGBImage anIm = cRGBImage::FromFile(aName1,1,eForceGray::Yes);
        anIm.ResetGray();

        for (const auto & aCple : aSetH.SetH())
        {
            cPt3dr aPG = aSens1->PInterBundle(aCple,*aSens2);
            cPt2dr aVect = aSens1->Ground2Image(aPG)-aCple.mP1;

           // StdOut() << "POL" << ToPolar(aVect,1000.0) << "\n";

            anIm.DrawCircle(cRGBImage::Blue,aCple.mP1,4.0);
            anIm.DrawLine(aCple.mP1,aCple.mP1+aVect*anAmpl,cRGBImage::Red,1);
        }
        anIm.ToFileDeZoom(mPhProj.DirVisuAppli()+"Res_"+LastPrefix(aName1)+"_"+LastPrefix(aName2)+".tif",aDeZoom);
}


void  cAppli_TiePReport::ProcessByPair()
{
    for (size_t aKIm1 = 0 ; aKIm1<mSetNames.size() ; aKIm1++)
    {       
        for (size_t aKIm2 = aKIm1+1 ; aKIm2<mSetNames.size() ; aKIm2++)
        {
           ProcessOnePair(mSetNames.at(aKIm1),mSetNames.at(aKIm2));
        }
    }
}


void cAppli_TiePReport::ProcessPly()
{

}


void cAppli_TiePReport::MakeStatByImage()
{
   InitReportCSV(mPrefixCSVIma,"csv",false);
   AddStdHeaderStatCSV(mPrefixCSVIma,"Image",mPropStat,{"AvgX","AvgY"});

   // const std::vector<std::list<std::pair<size_t,tPairTiePMult*>>>
   const auto & aVII  = mCMTP->IndexeOfImages();
   for (size_t aKImGlob = 0 ; aKImGlob<mSetNames.size() ; aKImGlob++)
   {
       // std::list<std::pair<size_t,tPairTiePMult*>>
       const auto & aLInd = aVII.at(aKImGlob);
       cSensorImage * aSensI = mCMTP->VSensors().at(aKImGlob);
       std::string aNameIm = aSensI->NameImage();


     //  std::string aNameCal = mCMTP->VSensors()

       // a litle check on indexe for these complexe structures
       MMVII_INTERNAL_ASSERT_tiny(mSetNames[aKImGlob]==aNameIm,"Chek names in Appli_TiePReport::MakeStatByImage");

       cWeightAv<tREAL8,cPt2dr>  aAvg2d;
       cStdStatRes               aStat;

       cSetIm4SparseDist<tREAL8>* aImAvg = nullptr;

       if (IsInit(&mPatImRes) && MatchRegex(aNameIm,mPatImRes))
       {
           aImAvg = new cSetIm4SparseDist<tREAL8>({"N2","x","y"},aSensI->Sz(),mFactRedIm);
       }

       std::string aNameCalib;
       cSetIm4SparseDist<tREAL8>* aImAvgCalib = nullptr;

     //  StdOut() <<  "PCCCCC= "  << aSensI->IsSensorCamPC() << " " <<  mDoImResCalib << "\n";
       if (aSensI->IsSensorCamPC() && mDoImResCalib)
       {
           aNameCalib = aSensI->GetSensorCamPC()->InternalCalib()->Name();
           if (! MapBoolFind(mMapResByCam,aNameCalib))
               mMapResByCam[aNameCalib] = new cSetIm4SparseDist<tREAL8>({"N2","x","y"},aSensI->Sz(),mFactRedIm);
           aImAvgCalib =   mMapResByCam[aNameCalib];
       }

       for (const auto& aPairKC : aLInd)
       {
           // std::pair<size_t,tPairTiePMult*W
           size_t aKImLoc =  aPairKC.first;
           // Unused in mode release
           [[maybe_unused]] const auto & aConfig = aPairKC.second->first;
           // a litle check on indexe for these complexe structures
           MMVII_INTERNAL_ASSERT_tiny(aKImGlob==(size_t)aConfig.at(aKImLoc),"Check nums in cAppli_TiePReport::MakeStatByImage");

           const auto & aVal = aPairKC.second->second;
           size_t aNbP = NbPtsMul(*aPairKC.second);
           // size_t aMult = Multiplicity(*aPairKC.second);

           for (size_t aKp=0 ; aKp<aNbP ; aKp++)
           {
               // cPt2dr aPt  = aVal.mVPIm.at(aKImLoc + aMult*aKp);
               cPt2dr aPt =  KthPt(*aPairKC.second,aKImLoc,aKp);
               cPt3dr aPGr = aVal.mVPGround.at(aKp);
               cPt2dr aResidu = aSensI->Ground2Image(aPGr) - aPt;

               aStat.Add(Norm2(aResidu));
               aAvg2d.Add(1.0,aResidu);

               if (aImAvg)
               {
                   aImAvg->Add(aPt,{Norm2(aResidu),aResidu.x(),aResidu.y()});
               }
               if (aImAvgCalib)
                   aImAvgCalib->Add(aPt,{Norm2(aResidu),aResidu.x(),aResidu.y()});

           }

       }
     //  StdOut() << "NIIIIIIiimmm= " << aNameIm << " " << aImAvg << "\n";
       if (aImAvg)
       {
       //    StdOut() << "aImAvgaImAvg " << mPhProj.DirVisuAppli() << " " << aNameIm << "\n";
           aImAvg->SaveDenseSave(mPhProj.DirVisuAppli(),aNameIm);
           delete aImAvg;
       }

       AddStdStatCSV
       (
          mPrefixCSVIma,aNameIm,aStat,mPropStat,
          {ToStr(aAvg2d.Average().x()),ToStr(aAvg2d.Average().y())}
       );
      // StdOut() << aNameIm  << " Avg=" << aStat.Avg() << std::endl;
   }


}


int cAppli_TiePReport::Exe()
{
   mPhProj.FinishInit();

   mSetNames = VectMainSet(0);
   // ...    bool  WithIndexPt,   bool  WithSensor,bool  WithIndexImages
   mCMTP = AllocStdFromMTP(mSetNames,mPhProj,true,true,true);
   mCMTP->SetPGround();

   mPrefixCSV =  "_Ori-"+  mPhProj.DPOrient().DirIn() +  "_Mes-"+  mPhProj.DPMulTieP().DirIn() ;
   mPrefixCSVIma =  "ByImages" + mPrefixCSV;


   MakeStatByImage();

   if (IsInit(&mParamByPair))
       ProcessByPair();

   if (mWithPly)
       ProcessPly();

   for (auto [aName,aAvgByCam] : mMapResByCam)
   {
       aAvgByCam->SaveDenseSave(mPhProj.DirVisuAppli(),aName);
       delete aAvgByCam;
   }

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

