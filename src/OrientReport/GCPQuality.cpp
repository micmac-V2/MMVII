#include "MMVII_Ptxd.h"
#include "cMMVII_Appli.h"
#include "MMVII_Geom3D.h"
#include "MMVII_PCSens.h"
#include "MMVII_Tpl_Images.h"
#include "MMVII_StaticLidar.h"

/**
   \file GCPQuality.cpp


 */

namespace MMVII
{


/* ==================================================== */
/*                                                      */
/*                  cSetIm4SparseDist                   */
/*                                                      */
/* ==================================================== */

template <class Type>
   cSetIm4SparseDist<Type>::cSetIm4SparseDist
   (
           const std::vector<std::string> & aVNames,
           const cPt2di& aSzInit,
           tREAL8 aFactRed
    ) :
      mFiltered  (false),
      mNames     (aVNames),
      mNbIm      (mNames.size()),
      mNbMeasure (0.0),
      mSumW      (0.0),
      mSzInit    (aSzInit),
      mFactRed   (aFactRed),
      mSzRed     (Pt_round_up(ToR(mSzInit)/aFactRed)),
      mImW       (mSzRed,nullptr,eModeInitImage::eMIA_Null)
{
     for (size_t aKIm=0 ; aKIm<mNbIm ; aKIm++)
         mIm2Avg.push_back(tIm(mSzRed,nullptr,eModeInitImage::eMIA_Null));
}

template <class Type> void cSetIm4SparseDist<Type>::Add(const cPt2dr &aPFull,const std::vector<tREAL8> & aVValues,tREAL8 aWeight)
{
   MMVII_INTERNAL_ASSERT_tiny(!mFiltered,"Already filtered in MakeDense");

    mNbMeasure++;
    mSumW += aWeight;
    cPt2dr aPRed = aPFull/mFactRed;

    if (!mImW.DIm().InsideBL(aPRed))
        return;
    mImW.DIm().AddVBL(aPRed,aWeight);
    MMVII_INTERNAL_ASSERT_always(aVValues.size()==mNbIm,"Bad size in cSetIm4SparseDist::Add");

    for (size_t aKIm=0 ; aKIm<mNbIm ; aKIm++)
    {
       mIm2Avg[aKIm].DIm().AddVBL(aPRed,aVValues.at(aKIm)*aWeight);
    }
}

template <class Type> void cSetIm4SparseDist<Type>::MakeDense(tREAL8 aMulSigma)
{
    MMVII_INTERNAL_ASSERT_always(!mFiltered,"Already filtered in MakeDense");
    mFiltered = true;
    //  S^2 * aNbPtsTot = Sz.x() * Sz.y()
    tREAL8 aSigma = std::sqrt(mImW.DIm().NbElem() / tREAL8(mNbMeasure)) * aMulSigma;
    aSigma = std::sqrt(Square(aSigma)+Square(1.0));
    mImW = mImW.GaussFilter(aSigma);

    SetMaxCsteInPlace(mImW.DIm(),mImW.DIm(),Type(1e-10));
    for (size_t aKIm=0 ; aKIm<mNbIm ; aKIm++)
    {
        mIm2Avg[aKIm] = mIm2Avg[aKIm].GaussFilter(aSigma);
        DivImageInPlace(mIm2Avg[aKIm].DIm(),mIm2Avg[aKIm].DIm(),mImW.DIm());
    }
}

template <class Type> void cSetIm4SparseDist<Type>::Gen1File
                           (
                               tIm anIm,
                               const std::string& aDir,
                               std::string aLayer,
                               const std::string & aPref
                           )
{
  aLayer = aLayer + std::string(mFiltered?"_Filtered_" : "_Raw_");
  std::string aName = aDir + aLayer + aPref + ".tif";

  anIm.DIm().ToFile(aName);
}



template <class Type> void cSetIm4SparseDist<Type>::GenFiles(const std::string& aDir,const std::string & aPref,bool WithLayer)
{
       Gen1File(mImW,aDir,"Weight",aPref);
       if (WithLayer)
       {
           for (size_t aKIm=0 ; aKIm<mNbIm; aKIm++)
           {
               Gen1File(mIm2Avg.at(aKIm),aDir,mNames.at(aKIm),aPref);
           }
       }
}

template <class Type> void cSetIm4SparseDist<Type>::SaveDenseSave
                           (
                                  const std::string& aDir,
                                  const std::string & aPref,
                                  bool WithLayebefore,
                                  bool WithLayerAfter
                            )
{
    GenFiles(aDir,aPref,WithLayebefore);
    MakeDense();
    GenFiles(aDir,aPref,WithLayerAfter);
}



template <class Type>  cIm2D<Type>  cSetIm4SparseDist<Type>::ImW() const {return mImW;}
template <class Type>  cIm2D<Type>  cSetIm4SparseDist<Type>::ImAvg(size_t aKth) const {return mIm2Avg.at(aKth);}
//tIm ImAvg(size_t aKTh) const;

template class cSetIm4SparseDist<tREAL4>;
template class cSetIm4SparseDist<tREAL8>;


/* ==================================================== */
/*                                                      */
/*          cAppli_CalibratedSpaceResection             */
/*                                                      */
/* ==================================================== */

class cAppli_CGPReport : public cMMVII_Appli
{
     public :

        cAppli_CGPReport(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli &);
        int Exe() override;
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override;

        std::vector<std::string>  Samples() const override;


     private :
        /** make the report  by image, for each image a cvs file with all GCP,
         * optionnaly make a visualisation of the residual fielsd for each image */
        void  MakeOneIm(const std::string & aNameIm);

        /** Make a report with an average for each GCP */
        void  ReportsByGCP();
        /** Make a visualization of residual in sensor plane*/
        void  ReportsByCam();

        std::string              mSpecImIn;   ///  Pattern of xml file
        cPhotogrammetricProject  mPhProj;

        std::vector<double>      mGeomFiedlVec;
        tREAL8                   mFactRed; ///< Reduction factor
        std::vector<int>         mPropStat;

        std::string              mPrefixReport;

        std::string              mNameReportDetail;
        std::string              mNameReportIm;
        std::string              mNameReportGCP;
        std::string              mNameReportGCP_Ground;
        std::string              mNameReportGCP_Ground_Glob;
        std::string              mNameReportCam;
        std::string              mNameReportMissed;

        double                   mMarginMiss;  ///  Margin for counting missing targets
        std::string              mSuffixReportSubDir; // additional name for report subdir
        std::string              mFilterName;  // pattern to filter names of GCP
        std::string              mFilterAdd;  // pattern to filter GCP by additional info

        std::vector<tREAL8>     mParamVF;
        std::string             mPatternParamIVF;
};

cAppli_CGPReport::cAppli_CGPReport
(
     const std::vector<std::string> &  aVArgs,
     const cSpecMMVII_Appli & aSpec
) :
     cMMVII_Appli  (aVArgs,aSpec),
     mPhProj       (*this),
     mFactRed      (100.0),
     mPropStat     ({50,75}),
     mMarginMiss   (50.0),
     mSuffixReportSubDir (""),
     mFilterName         (""),
     mFilterAdd          (""),
     mParamVF            (cImageVectorField::DefParam())
{
}



cCollecSpecArg2007 & cAppli_CGPReport::ArgObl(cCollecSpecArg2007 & anArgObl)
{
    return anArgObl
                << Arg2007(mSpecImIn,"Pattern/file for images",{{eTA2007::MPatFile,"0"},{eTA2007::FileDirProj}})
                << mPhProj.DPGndPt3D().ArgDirInMand()
                << mPhProj.DPGndPt2D().ArgDirInMand()
                << mPhProj.DPOrient().ArgDirInMand();
}

cCollecSpecArg2007 & cAppli_CGPReport::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
    return anArgOpt
                << AOpt2007(mPropStat,"Perc","Percentil for stat exp",{eTA2007::HDV})
                << AOpt2007(mGeomFiedlVec,"GFV","Geom Fiel Vect for visu [Mul,Witdh,Ray,Zoom?=2,JPeg?=1]",{{eTA2007::ISizeV,"[3,5]"}})
                << AOpt2007(mMarginMiss,"MargMiss","Margin to border for counting missed target",{eTA2007::HDV})
                << AOpt2007(mSuffixReportSubDir, "Suffix", "Suffix to report subdirectory name")
                << AOpt2007(mFilterName, "Filter", "Pattern to filter GCP by name")
                << AOpt2007(mFilterAdd, "FilterAdd", "Pattern to filter GCP by additional info")
                << AOpt2007(mFactRed, "RedFact", "Factor of reduction for sensor images (dist & bias)",{eTA2007::HDV})
                << cImageVectorField::ArgOpt(anArgOpt,mParamVF)
                <<  AOpt2007( mPatternParamIVF,"PatImFV","Pattern of  images  where we generat field vector (if any)")
            ;

}

std::vector<std::string>  cAppli_CGPReport::Samples() const
{
    return {
        "MMVII ReportGCP \"C2_0002.*JPG\" targets-out-wiSigma Targets-2D BA-Block  GFV=[100,1,1,1,0] "
    };
}


//================================================


void cAppli_CGPReport::MakeOneIm(const std::string & aNameIm)
{
    if (! ExistFile(mPhProj.NameMeasureGCPIm(aNameIm,true)) )
       return ;

    cSetMesGndPt             aSetMes;
    mPhProj.LoadGCP3D(aSetMes,nullptr,"",mFilterName,mFilterAdd);

    // cSet2D3D aSet32;
    // mSetMes.ExtractMes1Im(aSet32,aNameIm);
    cSensorImage*  aCam = mPhProj.ReadSensor(aNameIm,true,false);

    mPhProj.LoadIm(aSetMes,aNameIm,nullptr,aCam);

    const cSetMesPtOf1Im  &  aSetMesIm = aSetMes.MesImInitOfName(aNameIm);


    cImageVectorField * aIVF = nullptr;
    if (IsInit(&mPatternParamIVF) && MatchRegex(aNameIm,mPatternParamIVF))
    {
        aIVF = new cImageVectorField(aNameIm,mParamVF);
       for (const auto & aMes : aSetMesIm.Measures())
        {
            cPt2dr aP2 = aMes.mPt;
            cPt3dr aPGr = aSetMes.MesGCPOfName(aMes.mNamePt).mPt;
            cPt2dr aProj = aCam->Ground2Image(aPGr);
            if (Norm2(aProj-aP2)>10)
            {
                StdOut() << aP2 << aProj << "\n";
                 aIVF->Im().DrawLine(aP2,aProj,cRGBImage::Red);
            }
        }
    }

    // ------------ Generate the csv for residual : individual && aggregated by image -------------------

    cWeightAv<tREAL8,cPt2dr>  aAvg2d;  // Average of 2d-residual for bias
    cWeightAv<tREAL8,tREAL8>  aAvgD;   // for 3D distance measurments
    cStdStatRes               aStat;   // Statistique for residual

    for (const auto & aMes : aSetMesIm.Measures())
    {
        cPt2dr aP2 = aMes.mPt;
        cPt3dr aPGr = aSetMes.MesGCPOfName(aMes.mNamePt).mPt;
        cPt2dr aProj = aCam->Ground2Image(aPGr);
        cPt2dr  aVec = (aP2-aProj);
        std::string aResDist3DStr = "XXX";
        cStaticLidar * aStaticLidar = dynamic_cast<cStaticLidar*>(aCam);
        if (aStaticLidar)
        {
            aStaticLidar->ReadRasters(mPhProj.DirStaticLidarRasters());
            tREAL8 aMesDistance = aStaticLidar->Image2Distance(aP2);
            tREAL8 aResDist3D = aMesDistance - Norm2(aPGr-aStaticLidar->Center());
            aResDist3DStr = ToStr(aResDist3D);
            if (aMesDistance == 0.)
                aResDist3DStr = "XXX"; // no measurement
            aAvgD.Add(1.0, aResDist3D);
        }
        aAvg2d.Add(1.0,aVec);
        tREAL8 aDist = Norm2(aVec);
        aStat.Add(aDist);
        AddOneReportCSV(mNameReportDetail,{aNameIm,aMes.mNamePt,ToStr(aDist),ToStr(aVec.x()),ToStr(aVec.y()),aResDist3DStr});
    }


    for (const auto & aGCP : aSetMes.MesGCP())
    {
        if (starts_with(aGCP.mNamePt,"SL"))
        {
            AddOneReportCSV(mNameReportMissed,{aNameIm,aGCP.mNamePt,"UnObservable","XXX","XXX"});
        }
        else
        {
            if (aCam->IsVisible(aGCP.mPt))
            {
               cPt2dr aPIm = aCam->Ground2Image(aGCP.mPt);
               tREAL8 aDeg = aCam->DegreeVisibilityOnImFrame(aPIm);
               if (!aSetMesIm.NameHasMeasure(aGCP.mNamePt))
               {
                  if (aDeg> mMarginMiss)
                  {
                      AddOneReportCSV(mNameReportMissed,{aNameIm,aGCP.mNamePt,"UnDetected",ToStr(aPIm.x()),ToStr(aPIm.y())});

                      if (aIVF)
                      {
                          aIVF->DrawPointsRemarq(aPIm);
                      }
                  }
                  else
                      AddOneReportCSV(mNameReportMissed,{aNameIm,aGCP.mNamePt,"Margin",ToStr(aPIm.x()),ToStr(aPIm.y())});
               }
               else
               {
                   cPt2dr aP2 = aSetMesIm.MeasuresOfName(aGCP.mNamePt).mPt;
                   cPt2dr aProj = aCam->Ground2Image(aGCP.mPt);
                   if (aIVF)
                      aIVF->DrawArrow_P1P2(aP2,aProj);
               }
            }
            else
            {
                AddOneReportCSV(mNameReportMissed,{aNameIm,aGCP.mNamePt,"ProjOutImage","XXX","XXX"});
            }
        }

    }


    auto aMesX = (aAvg2d.SW()>0.) ? ToStr(aAvg2d.Average().x()) : "XXX";
    auto aMesY = (aAvg2d.SW()>0.) ? ToStr(aAvg2d.Average().y()) : "XXX";
    auto aMesD = (aAvgD.SW()>0.) ? ToStr(aAvgD.Average()) : "XXX";
    AddStdStatCSV
    (
       mNameReportIm,aNameIm,aStat,mPropStat,
       {aMesX, aMesY, aMesD}
    );

    if (aIVF)
    {
        aIVF->Save(mPhProj.DirVisuAppli()+LastPrefix(aNameIm));
        delete aIVF;
    }
}



void cAppli_CGPReport::ReportsByGCP()
{
   cSetMesGndPt             aSetMes;
   mPhProj.LoadGCP3D(aSetMes,nullptr,"",mFilterName,mFilterAdd);

   for (const auto & aNameIm : VectMainSet(0))
   {
       mPhProj.LoadIm(aSetMes,aNameIm,nullptr,mPhProj.ReadSensor(aNameIm,true,false),true);
   }

   const std::vector<cSensorImage*> &  aVSens =  aSetMes.VSens() ;

   InitReportCSV(mNameReportGCP,"csv",false);
   AddStdHeaderStatCSV(mNameReportGCP,"GCP",mPropStat);

   InitReportCSV(mNameReportGCP_Ground,"csv",false);
   AddOneReportCSV(mNameReportGCP_Ground,{"Name","Dx","Dy","Dz"});

   std::vector<cStdStatRes> aVStatXYZ{cStdStatRes(),cStdStatRes(),cStdStatRes()};

   for (const auto &  aMesIm :  aSetMes.MesImOfPt())
   {
        const auto & aGCP  = aSetMes.MesGCPOfMulIm(aMesIm);
        const std::vector<int> &  aVIndI = aMesIm.VImages() ;
        cStdStatRes               aStat;

        for (size_t aKIm = 0 ; aKIm<  aVIndI.size() ; aKIm++)
        {
            aStat.Add(Norm2( aMesIm.VMeasures()[aKIm]  - aVSens[aVIndI[aKIm]]->Ground2Image(aGCP.mPt)));
        }
        AddStdStatCSV(mNameReportGCP,aGCP.mNamePt,aStat,mPropStat);
    if (aVIndI.size()>1)
    {
        cPt3dr aDelta = aGCP.mPt -  aSetMes.BundleInter(aMesIm);
        AddOneReportCSV(mNameReportGCP_Ground,{aGCP.mNamePt,ToStr(aDelta.x()),ToStr(aDelta.y()),ToStr(aDelta.z())});

        for (int aKC=0 ; aKC<3 ; aKC++)
            aVStatXYZ[aKC].Add(aDelta[aKC]);
    }
    else
       AddOneReportCSV(mNameReportGCP_Ground,{aGCP.mNamePt,"xxx","yyy","zzz"});
   }

   InitReportCSV(mNameReportGCP_Ground_Glob,"csv",false);
   AddStdHeaderStatCSV(mNameReportGCP_Ground_Glob,"Coord",{});
   std::vector<std::string> aVCoord{"x","y","z"};
   for (int aKC=0 ; aKC<3 ; aKC++)
        AddStdStatCSV(mNameReportGCP_Ground_Glob,aVCoord[aKC],aVStatXYZ[aKC],{});
}

void cAppli_CGPReport::ReportsByCam()
{
   std::map<cPerspCamIntrCalib*,std::vector<cSensorCamPC*>>  aMapCam;
   cSetMesGndPt             aSetMes;
   mPhProj.LoadGCP3D(aSetMes,nullptr,"",mFilterName,mFilterAdd);

   for (const auto & aNameIm : VectMainSet(0))
   {
       cSensorCamPC *  aCam = mPhProj.ReadCamPC(aNameIm,true,true);
       if (aCam)
       {
            mPhProj.LoadIm(aSetMes,aNameIm,nullptr,aCam,true);
            aMapCam[aCam->InternalCalib()].push_back(aCam);
       }
   }

   InitReportCSV(mNameReportCam,"csv",false);
   AddStdHeaderStatCSV(mNameReportCam,"Cam",mPropStat);

   for (const auto& aPair : aMapCam)
   {
       cPerspCamIntrCalib * aCalib =  aPair.first;
       cSetIm4SparseDist<tREAL8> aImAvg({"x","y"},aCalib->SzPix(),mFactRed);
       cStdStatRes               aStat;

       for (const auto & aCam : aPair.second)
       {
           cSet2D3D aSet32;
           aSetMes.ExtractMes1Im(aSet32,aCam->NameImage(),true);

           for (const auto & aPair : aSet32.Pairs())
           {
               cPt2dr aP2 = aPair.mP2;
               cPt2dr aRes = (aCam->Ground2Image(aPair.mP3) - aP2) ;
               aStat.Add(Norm2(aRes));
               aImAvg.Add(aP2,aRes.ToStdVector());
           }
       }
       aImAvg.SaveDenseSave(mPhProj.DirVisuAppli(), aCalib->Name());


       AddStdStatCSV(mNameReportCam,aCalib->Name(),aStat,mPropStat);
   }
}





int cAppli_CGPReport::Exe()
{
   mPhProj.FinishInit();

   auto nameSubDir = mPhProj.DPOrient().DirIn() +  "_Mes-"+  mPhProj.DPGndPt3D().DirIn()
                                                +  "-"+  mPhProj.DPGndPt2D().DirIn();
   if (IsInit(&mSuffixReportSubDir))
       nameSubDir += "_" + mSuffixReportSubDir;
   SetReportSubDir(nameSubDir);

   mNameReportIm   =  "ByImage" ;
   mNameReportDetail   =  "Detail" ;
   mNameReportGCP  =  "ByGCP"   ;
   mNameReportCam   =  "ByCam"   ;

   mNameReportGCP_Ground   =  "ByGCP_3D"   ;
   mNameReportGCP_Ground_Glob   =  "ByGCP_3D_Stat"   ;

   mNameReportMissed   =  "MissedPoint"   ;

   InitReportCSV(mNameReportIm,"csv",true);
   InitReportCSV(mNameReportDetail,"csv",true);
   InitReportCSV(mNameReportMissed,"csv",true);

   if (LevelCall()==0)
   {
       AddStdHeaderStatCSV(mNameReportIm,"Image",mPropStat,{"AvgX","AvgY","AvgD"});
       AddOneReportCSV(mNameReportDetail,{"Image","GCP","Err","Dx","Dy","Ddist"});
       AddOneReportCSV(mNameReportMissed,{"Image","GCP","Natrue","XTh","YTh"});
   }

   // --------- If a pattern was used, run in // by a recall to itself  0->Param 0->Set  -------------
   if (RunMultiSet(0,0))
   {
       // After the run in // extract the result
      int aRes = ResultMultiSet();

      // Do the stat that agregate by Cam/GCP
      // N.B. this is not called when you have a single image; but it's not realy a pb (?)
      // because it would be equal to Image/Detail ...

       ReportsByGCP();
       ReportsByCam();


      return aRes;
   }


    MakeOneIm(FileOfPath(mSpecImIn,false));


   return EXIT_SUCCESS;
}

/* ==================================================== */
/*                                                      */
/*               MMVII                                  */
/*                                                      */
/* ==================================================== */



tMMVII_UnikPApli Alloc_CGPReport(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec)
{
   return tMMVII_UnikPApli(new cAppli_CGPReport(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpec_CGPReport
(
     "ReportGCP",
      Alloc_CGPReport,
      "Reports on GCP projection",
      {eApF::GCP, eApF::Ori},
      {eApDT::ObjCoordWorld, eApDT::ObjMesInstr, eApDT::Orient},
      {eApDT::Image,eApDT::Xml},
      __FILE__
);



}; // MMVII

