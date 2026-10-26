#include "cMMVII_Appli.h"
#include "MMVII_DeclareCste.h"
#include "MMVII_Geom3D.h"
#include "MMVII_2Include_Serial_Tpl.h"
#include "MMVII_PointCloud.h"
#include "cColorateCloud.h"
#include "MMVII_Image2D.h"

namespace MMVII
{
cPt3dr VertSphericalDir(const cPt2dr & aSun)
{
  return   spher2cart(cPt3dr(aSun.x()*(M_PI/2.0),M_PI/2.0+aSun.y(),1.0));
}
cPt3dr VertSphericalDir(const cPt3dr & aSun)
{
    return VertSphericalDir(Proj(aSun));
}

/* =============================================== */
/*                                                 */
/*                 cAppli_MMVII_CloudColorate      */
/*                                                 */
/* =============================================== */

enum eModeCloudCol
{
    eColShade,
    eColZ,
    eColXY,
    eSimulRand
};

class cAppli_MMVII_CloudColorate : public cMMVII_Appli
{
     public :

        cAppli_MMVII_CloudColorate(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec);

     private :
        int Exe() override;
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override ;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override ;

        // --- Mandatory ----
        std::string   mNameCloudIn;
        // --- Optionnal ----
        std::string mNameCloudOut;

        cPt2dr   mPropRayLeaf;
        tREAL8   mSurResol;
        int      mNbSampS;
        cPt3dr   mSun;
        bool     mShowMsg;
        bool     mExportIm;
        bool     mProfIsZ0;
        int           mIMCol;
        eModeCloudCol mModeCol;
};

cAppli_MMVII_CloudColorate::cAppli_MMVII_CloudColorate
(
     const std::vector<std::string> & aVArgs,
     const cSpecMMVII_Appli & aSpec
) :
     cMMVII_Appli    (aVArgs,aSpec),
     mPropRayLeaf    (1.0,1.0),
     mSurResol       (2.0),
     mNbSampS        (5),
     mShowMsg        (false),
     mExportIm       (false),
     mProfIsZ0       (false),
     mIMCol          (0),
     mModeCol        (eModeCloudCol::eColShade)
{
}

cCollecSpecArg2007 & cAppli_MMVII_CloudColorate::ArgObl(cCollecSpecArg2007 & anArgObl)
{
 return anArgObl
          <<   Arg2007(mNameCloudIn,"Name of input cloud/mesh", {eTA2007::FileDirProj,eTA2007::FileDmp})
   ;
}


cCollecSpecArg2007 & cAppli_MMVII_CloudColorate::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
   return anArgOpt
          << AOpt2007(mNameCloudOut,CurOP_Out,"Name of output file, def=Colorate_+InPut")
          << AOpt2007(mPropRayLeaf,"RayLeaves","Ray of leaves (/ avg dist)",{eTA2007::HDV})
          << AOpt2007(mSurResol,"SurResol","Sur resol in computation (/ avg dist)",{eTA2007::HDV})
          << AOpt2007(mNbSampS,"NbSampS","Number of sample/face for sphere discretization",{eTA2007::HDV})
          << AOpt2007(mSun,"Sun","Sun : Dir3D=(x,y,1)  ,  Z=WEIGHT !! ")
          << AOpt2007(mShowMsg,"ShowMsg","Print detailled message at each computation",{{eTA2007::HDV},{eTA2007::Tuning}})
          << AOpt2007(mExportIm,"ExportIm","Export all individual images",{{eTA2007::HDV},{eTA2007::Tuning}})
          << AOpt2007(mProfIsZ0,"ProfIsZ0","Prof is ZInit/\"Z in Dir proj\"",{{eTA2007::HDV},{eTA2007::Tuning}})
          << AOpt2007(mIMCol,"ICol","Col mode 0-Shde 1-Z 2-XY")
  ;
}


int  cAppli_MMVII_CloudColorate::Exe()
{
   if (! IsInit(&mNameCloudOut))
      mNameCloudOut = "Colorate_"+ mNameCloudIn;

  
   mModeCol = eModeCloudCol(mIMCol);
   cAutoTimerSegm aTSRead(TimeSegm(),"Read");
   cPointCloud   aPC_In ;
   ReadFromFile(aPC_In,mNameCloudIn);

   StdOut() << "ColNB COLIN=" << aPC_In.GetNbColours() << "\n";

   cIm2D<tREAL8> aImTextSim(cPt2di(1,1));
   tREAL8 aStep = 1 / std::sqrt(aPC_In.CurStdDensity());
   cBox2dr aBox = aPC_In.Box2d();

   if (mModeCol== eModeCloudCol::eSimulRand)
   {

       cPt3dr aC = aPC_In.Centroid();
       tREAL8 aZ0 = aC.z();

       cPt2di aSz = ToI(aBox.Sz());
       cRect2 aRect (cPt2di(0,0),aSz);

       std::vector<cParamNoise> aVParamZ{{100.0,1/aStep}};
       cIm2D<tREAL8> aImZ = CreateImageNoise(aSz,aZ0,aVParamZ);

       std::vector<cParamNoise> aVParamT{{30,30},{15,20},{5,15},{2,10}};
       aImTextSim = CreateImageNoise(aSz+cPt2di(3,3),128.0,aVParamT);

       cPointCloud aNewPC(false) ;
       aNewPC.SetOffset(aC);

       cPt2dr aLastP;
       for (const auto & aPix : aRect)
       {
           cPt2dr aP2G = aBox.P0() + ToR(aPix)*aStep;
           aNewPC.AddPt(TP3z(aP2G,aImZ.DIm().GetV(aPix)));
          // StdOut() << " ZZZ " << TP3z(aP2G,aImZ.DIm().GetV(aPix)) << "\n";
           aLastP = aP2G;
       }
     //  aNewPC.ComputeCurFineDensity();
       aPC_In = aNewPC;
   }

   // generate the sz of leaves
   if (! aPC_In.LeavesIsInit())
   {
       aPC_In.SetLeavesUnit(0.05,SVP::Yes);  // fix unit step,
       for (size_t aKPt=0 ; aKPt<aPC_In.NbPts() ; aKPt++)
       {
           tREAL8  aRayLeaf  = RandInInterval(mPropRayLeaf)  * aPC_In.GroundSampling();
           aPC_In.SetSzLeaves(aKPt,aRayLeaf);
       }
   }

   cAutoTimerSegm aTSInit(TimeSegm(),"Init");
   tREAL8 aWeightInit = (mNbSampS==0);
   cProjPointCloud  aPPC(aPC_In,aWeightInit);  // Weight Init 0 if NbS ,
// cProjPointCloud::cProjPointCloud(cPointCloud& aPC,tREAL8 aWeightInit) :

    
   cAutoTimerSegm aTSProj(TimeSegm(),"1Proj");

   int aNbStd=0;
   aPC_In.SetMulDegVis(1e4);
   if (mModeCol ==eModeCloudCol::eColShade)
   {
       if (mNbSampS>0)
       {
           cSampleSphere3D aSampS(mNbSampS);
           for (int aK=0 ; aK< aSampS.NbSamples() ; aK++)
           {
               cPt3dr aDir = VUnit(aSampS.KthPt(aK));
               if (aDir.z() >= 0.2)
               {
                   std::unique_ptr<cCamOrthoC> aCam (aPPC.PPC_CamOrtho(aK,mProfIsZ0,aDir));
                   cPt3di aDirI = ToI(aDir*100.0);
                   std::string aMsg = ToStr(aDirI.x()) + "_" +  ToStr(aDirI.y()) + "_" +  ToStr(aDirI.z());
                   aPPC.ProcessOneProj(mSurResol,*aCam,1.0,false,eModeImaDepth::eZGround,aMsg,mShowMsg,mExportIm);
                   aNbStd++;
                   StdOut() << "Still " << aSampS.NbSamples() - aK << "\n";
               }
           }
        }
        if (IsInit(&mSun))
        {
            tREAL8 aW0  = mNbSampS ? aNbStd : 1.0;
            std::unique_ptr<cCamOrthoC> aCam (aPPC.PPC_CamOrtho(0,mProfIsZ0,VertSphericalDir(mSun)));
            aPPC.ProcessOneProj(mSurResol,*aCam,aW0 * mSun.z(),false,eModeImaDepth::eZGround,"",false,false);
        }
        aPPC.ColorizePC();
   }
   else
   {

       if (1)
       {
           cTplBoxOfPts<tREAL8,2> aBoxOfPt;
           for (size_t aKPt=0 ; aKPt<aPC_In.NbPts(); aKPt++)
           {
               aBoxOfPt.Add(Proj( aPC_In.KthPt(aKPt)));
           }
           cBox2dr aBox = aBoxOfPt.CurBox();
           tREAL8 aSurf1Pt = aBox.NbElem() / aPC_In.NbPts();
           tREAL8 aStep = std::sqrt(aSurf1Pt);

           StdOut() << " NB PTS=" << aPC_In.NbPts() << " Step=" << aStep << " " << aPC_In.CurStdDensity() << "\n";

       }
       cBoundVals<tREAL8> aBounds;
       for (size_t aKPt=0 ; aKPt<aPC_In.NbPts(); aKPt++)
       {
           cPt3dr aPt = aPC_In.KthPt(aKPt);
           tREAL8 aDegVis=0.0;
           if (mModeCol==eModeCloudCol::eColXY)
           {
               int Ix=  round_ni(aPt.x());
               int Iy=  round_ni(aPt.y());

              // aDegVis =   ((1+std::sin(aPt.x() /3.0 )) * (1+std::sin(aPt.y() /9.0 )))  / 4.0 ;
               aDegVis =   (Ix/ 100 + Iy/100)%2;

              // aDegVis =  ((Ix%60)<5) || ((Iy%60)<5);
              // aDegVis = aDegVis>0.5;
           }
           if (mModeCol==eModeCloudCol::eSimulRand)
           {
               cPt2di aPix = ToI((Proj(aPt)- aBox.P0()) /aStep);
               tREAL8 aText = aImTextSim.DIm().GetV(aPix);
               aDegVis = std::max(0.0,std::min(1.0,aText / 255.0));
           }


           /*
           cPt2dr aP2G = aBox.P0() + ToR(aPix)*aStep;
           aNewPC.AddPt(TP3z(aP2G,aImZ.DIm().GetV(aPix)));
          // StdOut() << " ZZZ " << TP3z(aP2G,aImZ.DIm().GetV(aPix)) << "\n";
           aLastP = aP2G;
           */
           aPC_In.SetDegVis(aKPt,aDegVis);
           aBounds.Add(aPC_In.GetDegVis(aKPt) );
       }

       int aNbT=-1;
       for (int aKT=0 ; aKT<aNbT ; aKT++)
       {
           int aKPt = (aPC_In.NbPts()*aKT) / aNbT;
           StdOut() << aPC_In.GetDegVis(aKPt) << "\n";
       }
       StdOut() << "ColNB COLOUT=" << aPC_In.GetNbColours() << " DVInt="<< aBounds.VMin() << " " << aBounds.VMax() << "\n";

   }


   SaveInFile(aPC_In,mNameCloudOut);

   return EXIT_SUCCESS;
}

     /* =============================================== */
     /*                       ::                        */
     /* =============================================== */

tMMVII_UnikPApli Alloc_MMVII_CloudColorate(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli & aSpec)
{
   return tMMVII_UnikPApli(new cAppli_MMVII_CloudColorate(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpec_MMVII_CloudColorate
(
     "CloudMMVIIColorate",
      Alloc_MMVII_CloudColorate,
      "Generate a colorate version of  MMVII-Cloud",
      {eApF::Cloud},
      {eApDT::Ply},
      {eApDT::Ply},
      __FILE__
);


};
