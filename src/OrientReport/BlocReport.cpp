#include "MMVII_BlocRig.h"

#include "MMVII_Ptxd.h"
#include "cMMVII_Appli.h"
#include "MMVII_Geom3D.h"
#include "MMVII_PCSens.h"
#include "MMVII_Tpl_Images.h"
#include "MMVII_2Include_Serial_Tpl.h"
#include "MMVII_Clino.h"
#include "MMVII_InstrumentalBlock.h"


/**
 
   \file cWire3DInit.cpp

   This file contains a command for computing the 3D position of a wire from multiple
   images.

 */

namespace MMVII
{

/* ==================================================== */
/*                                                      */
/*          cAppli_CalibratedSpaceResection             */
/*                                                      */
/* ==================================================== */

class  cStatDistPtWire
{
    public :
    /*
        cStatDistPtWire():
            mIsGlob (false),
            mDist   (-1)
        {
        }

        bool         mIsGlob;
        tREAL8       mDist;
*/

        cStdStatRes  mStat3d;  // stat on 3D dist
        cStdStatRes  mStatH;  // stat on horiz dist
        cStdStatRes  mStatV;  // stat on vert dist
};



class cAppli_ReportBlock : public cMMVII_Appli
{
     public :
        typedef std::pair<cSensorCamPC *,cMesIm1Pt> tPairCamPt;

        cAppli_ReportBlock(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli &);
        int Exe() override;
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override;

        //  std::vector<std::string>  Samples() const override;

     private :
        void ProcessOneBloc(const std::string& anIdSync,const cIrbComp_TimeS &);
        void AddStatDistWirePt(const cPt3dr& aPt,const cPt3dr& aDirLoc,const std::string & anIdSync,const std::string & aNamePt);

        void TestWire3D(const std::string& anIdSync,const std::vector<cSensorCamPC *> & aVCam);
        /// For a given Id of Sync and a bloc of cameras, compute the stat on accuracy of intersection
        void TestPoint3D(const std::string& anIdSync,const std::vector<cSensorCamPC *> & aVCam,bool DoStat);

        std::string StrDev(const cStdStatRes & aStat) const
        {
            return (aStat.NbMeasures() >= 2) ? (ToStr(1e6*aStat.UBDevStd(-1))) : "???";
        }

        std::string StrDiff(const cStdStatRes & aS1,const cStdStatRes & aRef) const
        {
            return (aS1.NbMeasures() >= 1) ? (ToStr(1e6*(aS1.Avg()-aRef.Avg()))) : "???";
        }

        cPhotogrammetricProject     mPhProj;
        bool                        mShow;  // do we print residual on terminal
        std::string                 mSpecImIn;

        std::string                  mNameBloc;
        cIrbCal_Block*               mCalBlocInstr;
        cIrbComp_Block*              mCompBlocInstr;

        std::string                  mIdRepWire;
        std::string                  mIdRepPtIndiv;
        std::string                  mIdRepDWirePt;
        std::string                  mIdRepPtGlob;
        std::string                  mPatNameGCP;

        std::string                  mStrM2T;  /// String of measure to test
        std::string                  mAddExReport;
        std::string                  mDirExReport;
        cWeightAv<tREAL8,tREAL8>     mAvgGlobRes;
        cStdStatRes                  mStatGlobPt;
        cStdStatRes                  mStatGlobWire;
        std::vector<int>             mPercStat;
        std::map<std::string,cStdStatRes>    mMapStatPair;
        std::map<std::string,cStdStatRes>    mMap1Image;
      
        bool                         mWithTargetPoint; ///< is there targeted points
        bool                         mWithCoordGround; ///< do we have the 3D coordinat of these points

        //  Pair of 3d coordinate to compute transfo  Local in bloc of cam <-> "absolute" like sphere coord
        std::vector<cPt3dr>          mCoordGround;  ///<  3D coordinate "absolute"
        std::vector<cPt3dr>          mCoordLoc;     ///< 3D coordinate in bloc of cameras repai

        //  Add a statistic results in csv-file
        void CSV_AddStat(const std::string& anId,const std::string& aMes,const cStdStatRes &) ;

        /// 3d position of the wire that may be computed (if required and enough image)
        cSegmentCompiled<tREAL8,3> *  mCurWire ;
        bool                          mWithWire;  ///< Do we require Wire computation
        std::string                   mWireFolder; ///< Wire folder, that may differ from point folder

        /// mStatWirePt[NamePt][Loc]
        std::map<std::string, std::map<std::string,cStatDistPtWire>>   mStatWirePt;

        //  Stuff  for  CERN-Sphere-Center like manip
        bool                                    mSCFreeScale;  ///< do we have a free scale
        cSetMesGndPt                            mMesGround;     ///< "Absolute" coordinate of the point
        std::string                             mPatDistWire;   ///< Pattern of point name for computing distance point <--> WIRE
        std::string                             mPatAgregateDist;  ///< Pattern for aggregatinf distance

        bool                                    mWithClino;  /// Where clino measure loade
        bool                                    mWithAgregateDist; ///< Do we agregate the measure
        cSetMeasureClino                        mMesClino; ///< store clinometers measures

        /// For a given name of point memorize the set of pair "Cam+Im Measure"
        std::map<std::string,std::vector<tPairCamPt>> mMapMatch;

};




cAppli_ReportBlock::cAppli_ReportBlock
(
     const std::vector<std::string> &  aVArgs,
     const cSpecMMVII_Appli & aSpec
) :
     cMMVII_Appli   (aVArgs,aSpec),
     mPhProj        (*this),
     mShow          (false),
     mNameBloc      (cIrbCal_Block::theDefaultName),
     mCalBlocInstr  (nullptr),
     mCompBlocInstr (nullptr),
     mIdRepWire     ("Wire"),
     mIdRepPtIndiv  ("Pt"),
     mIdRepDWirePt  ("DistWP"),
     mIdRepPtGlob   ("GlobPt"),
     mPatNameGCP    (".*"),
     mStrM2T        ("T"),
     mPercStat      {15,25,50,75,85},
     mSCFreeScale   (false),
     mPatDistWire   ("CENTRE")

{
    FakeUseIt(mSCFreeScale);
}


cCollecSpecArg2007 & cAppli_ReportBlock::ArgObl(cCollecSpecArg2007 & anArgObl)
{
      return anArgObl
             <<  Arg2007(mSpecImIn,"Pattern/file for images", {{eTA2007::MPatFile,"0"},{eTA2007::FileDirProj}}  )
             <<  mPhProj.DPOrient().ArgDirInMand()
             <<  mPhProj.DPGndPt2D().ArgDirInMand()
             <<   mPhProj.DPBlockInstr().ArgDirInMand()
           ;
}


cCollecSpecArg2007 & cAppli_ReportBlock::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{

    return      anArgOpt

             << AOpt2007(mWireFolder,"WireFolder","Folder for wire if != folder points")
             << AOpt2007(mPatNameGCP,"PatFiltGCP","Pattern to filter name of GCP",{{eTA2007::HDV}})
             << AOpt2007(mStrM2T,"M2T","Measure to test : T-arget W-ire",{{eTA2007::HDV}})
             << AOpt2007(mAddExReport,"AddExRep","Addditional Extension in Report Name")
             << AOpt2007(mDirExReport,"DirExRep","Fix globally Directory of Report Name")
             << AOpt2007(mPercStat,"PercStat","Percentils for stat in global report",{{eTA2007::HDV}})

             << mPhProj.DPGndPt3D().ArgDirInOpt("","GCP 3D coordinate for computing centre")

             << AOpt2007(mPatDistWire,"Cern-PatDistWire","Pattern of points name to computing distance point<-->wire",{{eTA2007::HDV}})
             << AOpt2007(mPatAgregateDist,"Cern-AggregateDist","Pattern of points name to computing distance point<-->wire",{{eTA2007::HDV}})
             << mPhProj.DPMeasuresClino().ArgDirInOpt("Cern-Clino")

             << AOpt2007(mSCFreeScale,"Cern-SphFreeScale","Do we have free scale for sphere",{{eTA2007::HDV}})

             << AOpt2007(mShow,"Show","Show details on results",{{eTA2007::HDV}})
    ;
}

static const std::string StrGLOB = "GLOB";

void cAppli_ReportBlock::AddStatDistWirePt
     (
            const cPt3dr& aPt,
            const cPt3dr& aVertLoc,
            const std::string & anIm0,
            const std::string & aNamePt
     )
{
    if (!mCurWire)
       return;

    std::map<std::string,cStatDistPtWire> & aMapPt = mStatWirePt[aNamePt];

     cPt3dr  aProj = mCurWire->Proj(aPt);
     cPt3dr  anEc = aProj-aPt;
     tREAL8 aD3 = Norm2(anEc);
     tREAL8 aDH = -1;
     tREAL8 aDV = -1;

     std::string anIdGlob = StrGLOB;
     std::string anIdAgreg;

     if (mWithAgregateDist)
     {
         std::string aName = ReplacePattern(mPatAgregateDist,"$1",anIm0);
         int aL = cStrIO<int>::FromStr(aName);
         aName = ToStr(aL,3);
         anIdAgreg =  "LOC-"+aName;
     }

     if (mWithClino)
     {
        cPt3dr  aCompV = aVertLoc * Scal(aVertLoc,anEc);
        cPt3dr  aCompH = anEc-aCompV;
        aDH = Norm2(aCompH);
        aDV = Norm2(aCompV);
        aMapPt[anIdGlob].mStatH.Add(aDH);
        aMapPt[anIdGlob].mStatV.Add(aDV);
        if (mWithAgregateDist)
        {
            aMapPt[anIdAgreg].mStatH.Add(aDH);
            aMapPt[anIdAgreg].mStatV.Add(aDV);
        }
     }


     aMapPt[anIdGlob].mStat3d.Add(Norm2(anEc));
     if (mWithAgregateDist)
     {
         aMapPt[anIdAgreg].mStat3d.Add(Norm2(anEc));
     }

     // Add individual statistic for each point each image
     AddOneReportCSV(mIdRepDWirePt,{anIm0,aNamePt,ToStr(aD3),ToStr(aDH),ToStr(aDV)});
}





void cAppli_ReportBlock::TestWire3D(const std::string & anIdSync,const std::vector<cSensorCamPC *> & aVCam)
{
     std::vector<cPlane3D>  aVPlane;       // vector of plane
     std::vector<cSensorCamPC *>  aVCamOk; // camera for which we have plane
     std::vector<tSeg2dr>         aVSegOk; // segmment of camera corresponding to planes

     // ideally, we shoud use the class,  cCam2Wire_2Dto3D, but there is too much intrication
     // with CorrMesSeg, so for now I accept this "code duplication" ;-((

     // [1]  compute the planes
     for (const auto & aCam : aVCam)
     {
          const std::string & aNameIm = aCam->NameImage();
          if (mPhProj.HasFileLinesFolder(mWireFolder,aNameIm))
          {
            //  ReadLinesFolder

              // cLinesAntiParal1Im   aSetL  = mPhProj.ReadLines(aNameIm);
              cLinesAntiParal1Im   aSetL  = mPhProj.ReadLinesFolder(mWireFolder,aNameIm);

              const std::vector<cOneLineAntiParal> & aVL  = 	aSetL.mLines;

              // At this step we dont handle multiple lines
              if (aVL.size()==1)
              {
                 // comput seg, and correct it if we are at this step
                 tSeg2dr aSeg = aVL.at(0).mSeg;
                // aSeg =  CorrMesSeg(aCam,aSeg);

                 //  memorize plane, seg and cam
                 aVPlane.push_back(aCam->SegImage2Ground(aSeg));
                 aVCamOk.push_back(aCam);
                 aVSegOk.push_back(aSeg);
              }
          }
     }

     int aNbPl = aVPlane.size();

   //   StdOut() << "NBBBBB " << aNbPl << " " << mWireFolder << "\n";
     // if we can compute plane
     if (aNbPl>=2)
     {
        mCurWire  = new cSegmentCompiled<tREAL8,3>(cPlane3D::InterPlane(aVPlane));

        // if we are the step where we compute the correction
      /*  if (mStepCompCCorr)
           for (size_t aKC=0 ; aKC<aVCamOk.size() ; aKC++)
               mMapCorrSyst[aVCamOk.at(aKC)].AddPairSeg(*mCurWire,aVSegOk.at(aKC));*/

        // if we have enough plane to compute residuals
        if (aNbPl>=3)
        {
            cWeightAv<tREAL8> aWGr;  // Average of ground distance
            cWeightAv<tREAL8> aWPix; // Average of pixel distance

            // Parse the camera where a seg was detected
            for (size_t aKC=0 ; aKC<aVCamOk.size() ; aKC++)
            {
                cSensorCamPC * aCam = aVCamOk[aKC];
                cPerspCamIntrCalib * aCalib = aCam->InternalCalib();

                int aNSampleSeg = 5; //  Number of sample on the seg
                for (int aKS=0 ; aKS<=aNSampleSeg ; aKS++)
                {
                      // compute a point on curve corresponding to the undistorde line
                      cPt2dr aPIm = aCalib->InterpolOnUDLine(aVSegOk[aKC],aKS/tREAL8(aNSampleSeg));
                      aWGr.Add(1.0,aCam->GroundDistBundleSeg(aPIm,*mCurWire));
                      aWPix.Add(1.0,aCam->PixDistBundleSeg(aPIm,*mCurWire));
                }
            }
            //  if we are the step where we compute

            {
                tREAL8 aRatio = aNbPl /(aNbPl-2.0); // ratio contraints, degre of freedom
                tREAL8 aDist3D =  aWGr.Average() * aRatio;
                tREAL8 aDistPix =  aWPix.Average() * aRatio;

                mStatGlobWire.Add(aDistPix);
                AddOneReportCSV(mIdRepWire,{anIdSync,ToStr(aNbPl),ToStr(aDist3D),ToStr(aDistPix)});
            }
        }
    }
}


void cAppli_ReportBlock::CSV_AddStat(const std::string& anId,const std::string& aMes,const cStdStatRes & aStat)
{
   AddStdStatCSV(anId,aMes,aStat,mPercStat);
}

std::string ToCernStr(const cPt3dr & aPt)
{
   return ToStr(aPt.x()) + " " + ToStr(aPt.y()) + " " + ToStr(aPt.z());
}

void cAppli_ReportBlock::TestPoint3D
     (
         const std::string & anIdSync,
         const std::vector<cSensorCamPC *> & aVCam,
        bool doStat
     )
{
     // for a given name of point, store  Mes+Cam , that will allow to compute bundles
     mMapMatch.clear();
     cStdStatRes  aStatRes;

     // [1]  Parse all the camera to group measur by name of point (A) load points (B) parse them to store image measure + Cam
     for (const auto & aCam : aVCam)
     {
          int aNbMesOK = 0;
         // if images measures were  computed
          if (mPhProj.HasMeasureIm(aCam->NameImage()))
          {
              cSetMesPtOf1Im  aSet = mPhProj.LoadMeasureIm(aCam->NameImage()); // (A) Load the points

              for (const auto & aMes : aSet.Measures()) // (B) parse the points
              {
                 // Dont select points if NotCodes or not selected by user-regex-filtering
                 if ((!starts_with( aMes.mNamePt,MMVII_NONE)) && MatchRegex(aMes.mNamePt,mPatNameGCP))
                 {
                    aNbMesOK++;
                    mMapMatch[aMes.mNamePt].push_back(tPairCamPt(aCam,aMes));
                  }
             }
             if (aNbMesOK==0)
             {
                 StdOut() << "NO Measure valide  for " << aCam->NameImage() << "\n";
             }
         }
          else
             StdOut() << "NO Measure file  for " << aCam->NameImage() << "\n";
     }



     // [2]  Parse the measure grouped by points
     for (const auto & [aNamePt,aVect] : mMapMatch )
     {
         int aNbPt = aVect.size();
         if (aNbPt >  2)
         {

             // [2.1]  compute the vector of bundles and their intersection
             std::vector<tSeg3dr> aVSeg;
             for (const auto & [aCam,aMes] : aVect)
             {
                 aVSeg.push_back(aCam->Image2Bundle(aMes.mPt));
             }
             cPt3dr aPLoc =   BundleInters(aVSeg);
             const cPt3dr * aPGround = nullptr;

             if (mWithCoordGround &&  mMesGround.NameIsGCP(aNamePt))
             {
                 const cMes1Gnd3D &   aMes =      mMesGround.MesGCPOfName(aNamePt) ;
                 aPGround = & aMes.mPt;
                 mCoordGround.push_back(*aPGround);
                 mCoordLoc.push_back(aPLoc);

             }
             // [2.2] compute residual and eventually memo data for correction
             cWeightAv<tREAL8> aWPix;

             for (const auto & [aCam,aMes] : aVect)
             {
                 cPt2dr aPProj = aCam->Ground2Image(aPLoc);
                 cPt2dr aResidual = aMes.mPt-aPProj;
                 aWPix.Add(1.0,Norm2(aResidual));
             }
             
             // [2.3] export residual in csv file
             if (doStat)
             {

                 tREAL8 aDistPix = aWPix.Average() * (aNbPt*2.0) / (aNbPt*2.0 -3.0);
                 AddOneReportCSV(mIdRepPtIndiv,{anIdSync,aNamePt,ToStr(aNbPt),ToStr(aDistPix)});
                 aStatRes.Add(aDistPix);
                 mStatGlobPt.Add(aDistPix);

                 //  Now make the computation by pair of camera
                 for (size_t aK1=0 ; aK1<aVect.size() ; aK1++)
                 {
                     const auto & [aCam1,aMes1] = aVect.at(aK1);
                     for (size_t aK2=aK1+1 ; aK2<aVect.size() ; aK2++)
                     {
                         const auto & [aCam2,aMes2] = aVect.at(aK2);
                         cHomogCpleIm aCple(aMes1.mPt,aMes2.mPt);
                         // The weighting is not exaclt the same for pair of cam  ??
                         tREAL8 aRes12 = aCam1->PixResInterBundle(aCple,*aCam2) * 4.0;  // 4.0 = DOF = 4 / (4-3)
                         std::string anId1 = "Cam:"+aCam1->InternalCalib()->Name();
                         std::string anId2 = "Cam:"+aCam2->InternalCalib()->Name();
                         std::string aNamePair = anId1 + "/" + anId2;
                         mMapStatPair[aNamePair].Add(aRes12);
                         mMap1Image[anId1].Add(aRes12);
                         mMap1Image[anId2].Add(aRes12);
                     }
                 }
             }
         }
     }
     //  [3]  export global resiudal
     if (doStat)
        CSV_AddStat(mIdRepPtGlob,"AVG "+anIdSync,aStatRes);
}


void cAppli_ReportBlock::ProcessOneBloc(const std::string& anIdSync,const cIrbComp_TimeS & aDataTS)
{
     mCurWire = nullptr ;
     mCoordGround.clear();
     mCoordLoc.clear();

     std::vector<cSensorCamPC *>  aVCam;
     for (const cIrbComp_Cam1 & aCompCam : aDataTS.SetCams().VCompPoses())
         aVCam.push_back(aCompCam.CamPC());




     if (mWithTargetPoint)
     {
        TestPoint3D(anIdSync,aVCam,true);
     }

     if (mWithWire)
     {
         TestWire3D(anIdSync,aVCam);
     }

     // compute the stat specific to CERN
     if (mWithClino && (mCoordGround.size() >=3) && mCurWire)
     {
        auto  aWMin =  MMVII::ExtractVerticalLoc(aDataTS,false);
        tPoseR aPose = tPoseR::RansacL1Estimate(mCoordGround,mCoordLoc,1000);
        tSim3dR aSim = tSim3dR::RansacL1Estimate(mCoordGround,mCoordLoc,1000);

        aPose = aPose.LeastSquareRefine(mCoordGround,mCoordLoc);
        aSim = aSim.LeastSquareRefine(mCoordGround,mCoordLoc);


        for (const auto &  aMes3d : mMesGround.MesGCP())
        {

            if (MatchRegex(aMes3d.mNamePt,mPatDistWire))
            {
                cPt3dr aPtLoc = mSCFreeScale                 ?
                                   aSim.Value(aMes3d.mPt)    :
                                   aPose.Value(aMes3d.mPt)   ;
                AddStatDistWirePt(aPtLoc,aWMin.IndexExtre(),aVCam[0]->NameImage(),aMes3d.mNamePt);
            }
        }
     }


     delete mCurWire;
}

int cAppli_ReportBlock::Exe()
{
    mPhProj.FinishInit();  // the final construction of  photogrammetric project manager can only be done now

    mWithTargetPoint = contains(mStrM2T,'T');
    mWithWire = contains(mStrM2T,'W') || IsInit(&mWireFolder);

    mWithClino = mPhProj.DPMeasuresClino().DirInIsInit() ;
    mWithCoordGround = mPhProj.DPGndPt3D().DirInIsInit();
    mWithAgregateDist = mWithClino && mWithClino && IsInit(&mPatAgregateDist);

    if (mWithWire && !IsInit(&mWireFolder))
        mWireFolder = mPhProj.DPGndPt2D().DirIn();

    if (mWithCoordGround)
    {
        mPhProj.LoadGCP3D(mMesGround);
    }


    std::string aDirRep =          mPhProj.DPOrient().DirIn()
                           + "-" + mPhProj.DPGndPt2D().DirIn()
                           + "-" + mPhProj.DPBlockInstr().DirIn() ;
    if (IsInit(&mAddExReport))
       aDirRep =  mAddExReport + "-" + aDirRep;
    if (IsInit(&mDirExReport))
       aDirRep = mDirExReport;
    SetReportSubDir(aDirRep);


    InitReportCSV(mIdRepWire,"csv",false,{"TimeBloc","NbPlane","Dist Ground","Dist Pix"});
    // AddOneReportCSV(mIdRepWire,{"TimeBloc","NbPlane","Dist Ground","Dist Pix"});
    
    InitReportCSV(mIdRepPtIndiv,"csv",false,{"TimeBloc","Point","Mult","Dist Pix"});

    InitReportCSV(mIdRepDWirePt,"csv",false,{"TimeStamp","NamePt","D3","DH","DV"});

    InitReportCSV(mIdRepPtGlob,"csv",false);
    AddStdHeaderStatCSV(mIdRepPtGlob,"NameAggreg",mPercStat);

    mCompBlocInstr = new cIrbComp_Block(mPhProj,mNameBloc);
    mCalBlocInstr =  & mCompBlocInstr->CalBlock();
    mCompBlocInstr->AddImagesPoses(VectMainSet(0),false,true);

    if (mWithClino )
    {
        mCompBlocInstr->SetClinoValues(true);
    }


    for (const auto & [aTimeS,aDataTS] : mCompBlocInstr->DataTS() )
    {
        ProcessOneBloc(aTimeS,aDataTS);
    }

   // Add the stat for all pairs
   for (const auto & [aNameImage,aStatImage] : mMap1Image )
       CSV_AddStat(mIdRepPtGlob,aNameImage,aStatImage);

   // Add the stat for all pairs
   for (const auto & [aNamePair,aStatPair] : mMapStatPair )
       CSV_AddStat(mIdRepPtGlob,aNamePair,aStatPair);

   // Add the stat for all the points
    CSV_AddStat(mIdRepPtGlob,"GlobAVG ",mStatGlobPt);
    StdOut() << mStatGlobPt.Show("Pt-Pix",{50,85}) << "\n";

    if (mShow)
    {

            StdOut() << mStatGlobPt.Show("Pt-Pix",{50,85}) << "\n";
            StdOut() << mStatGlobWire.Show("Wire-Pix",{50,85}) << "\n";
    }

    //  ----------------- Make the export specific to CERN : distance to wire ------------

    for ( auto & [aNamePt,aMap] : mStatWirePt ) // parse all points
    {

        std::string anIdStatWire ="StatDistWire-"+aNamePt;

        std::vector<std::string> aVHeader{"Name","Nb","3D-Avg","3D-Delta micr","3D-Dev micr"};
        if (mWithClino)
            AppendIn(aVHeader,{"H-Avg","H-Delta micr","H-Dev","V-Avg","V-Delta micr","V-Dev"});

        InitReportCSV(anIdStatWire,"csv",false,aVHeader);
        const auto& aRef = aMap[StrGLOB];

        for (const auto & [aName,aStat] : aMap ) // parse all stats
        {
           // size_t aNb = aStat.mStat3d.NbMeasures();
            std::vector<std::string> aVStat
                                     {
                                          aName,
                                          ToStr(aStat.mStat3d.NbMeasures()),
                                          aStat.mStat3d.StrAvg(),
                                          StrDiff(aStat.mStat3d,aRef.mStat3d),
                                          StrDev(aStat.mStat3d)
                                     };
            if (mWithClino)
            {
                  AppendIn(aVStat,{aStat.mStatH.StrAvg(),StrDiff(aStat.mStatH,aRef.mStatH),StrDev(aStat.mStatH)});
                  AppendIn(aVStat,{aStat.mStatV.StrAvg(),StrDiff(aStat.mStatV,aRef.mStatV),StrDev(aStat.mStatV)});
            }
            AddOneReportCSV(anIdStatWire,aVStat);
        }
    }


    delete mCompBlocInstr;
    // StdOut() << "cIrbCal_BlockcIrbCal_Block\n";
    return EXIT_SUCCESS;
}

/* ==================================================== */
/*                                                      */
/*               MMVII                                  */
/*                                                      */
/* ==================================================== */

tMMVII_UnikPApli Alloc_ReportBlock(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec)
{
   return tMMVII_UnikPApli(new cAppli_ReportBlock(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpec_BlocReport
(
     "ReportBlock",
      Alloc_ReportBlock,
     "Report different measures relative to a block of cam",
      {eApF::Ori},
      {eApDT::Orient},
      {eApDT::Xml},
      __FILE__
);

}; // MMVII

