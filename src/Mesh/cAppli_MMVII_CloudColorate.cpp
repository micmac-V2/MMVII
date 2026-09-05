#include "cMMVII_Appli.h"
#include "MMVII_DeclareCste.h"
#include "MMVII_Geom3D.h"
#include "MMVII_Sensor.h"
#include "MMVII_2Include_Serial_Tpl.h"
#include "MMVII_Tpl_Images.h"
#include "MMVII_PointCloud.h"
#include "MMVII_Linear2DFiltering.h"
#include "MMVII_Interpolators.h"
#include "MMVII_PCSens.h"
#include "MMVII_2Include_Tiling.h"
#include "../ImagesBase/cGdalApi.h"

#include "cColorateCloud.h"

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
    eColOrtho,
    eColXY
};

/** Light object to put VISIBLE points (with their true ortho colour) in a
 *  spatial index, so occluded points can borrow the nearest visible colour.
 *  The index is 2D (on X,Y) because cTiling::GetObjAtDist only supports 2D
 *  (it uses cRect2==cPixBox<2>), but we keep the full 3D point so candidates
 *  can be ranked by 3D distance : this avoids borrowing from the canopy point
 *  located right above an occluded ground point (small XY dist, large Z). */
class cVisPtCol
{
    public :
        static constexpr int TheDim = 2;      // 2D indexing on (X,Y)
        typedef cPt2dr tPrimGeom;             // geometric primitive is a 2D point
        typedef int    tArgPG;                // unused call-back argument

        const tPrimGeom & GetPrimGeom(int =-1) const {return mP2;}
        cVisPtCol(const cPt3dr & aPt,tREAL8 aCol) : mP2(Proj(aPt)),mP3(aPt),mCol(aCol) {}
        const cPt3dr & Pt3d() const {return mP3;}
        tREAL8 Col() const {return mCol;}
    private :
        cPt2dr mP2;
        cPt3dr mP3;
        tREAL8 mCol;
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
        std::string   mOrthoFile;
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
          << AOpt2007(mOrthoFile,"OrthoFile","Ortho image to colorate the cloud",{eTA2007::HDV})
          << AOpt2007(mSun,"Sun","Sun : Dir3D=(x,y,1)  ,  Z=WEIGHT !! ")
          << AOpt2007(mShowMsg,"ShowMsg","Print detailled message at each computation",{{eTA2007::HDV},{eTA2007::Tuning}})
          << AOpt2007(mExportIm,"ExportIm","Export all individual images",{{eTA2007::HDV},{eTA2007::Tuning}})
          << AOpt2007(mProfIsZ0,"ProfIsZ0","Prof is ZInit/\"Z in Dir proj\"",{{eTA2007::HDV},{eTA2007::Tuning}})
          << AOpt2007(mIMCol,"ICol","Col mode 0-Shde 1-Z 2-Ortho 3-XY")
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
                   aPPC.ProcessOneProj(mSurResol,*aCam,1.0,false,aMsg,mShowMsg,mExportIm);
                   aNbStd++;
                   StdOut() << "Still " << aSampS.NbSamples() - aK << "\n";
               }
           }
        }
        if (IsInit(&mSun))
        {
            tREAL8 aW0  = mNbSampS ? aNbStd : 1.0;
            std::unique_ptr<cCamOrthoC> aCam (aPPC.PPC_CamOrtho(0,mProfIsZ0,VertSphericalDir(mSun)));
            aPPC.ProcessOneProj(mSurResol,*aCam,aW0 * mSun.z(),false,"",false,false);
        }
        aPPC.ColorizePC();
   }
   else if(mModeCol == eModeCloudCol::eColOrtho)
   {
        // we colorize with an ortho image, suppose one orthocamera 
        std::unique_ptr<cCamOrthoC> aCam (aPPC.PPC_CamOrtho(0,mProfIsZ0,cPt3dr(0,0,1.0)));

        aPPC.ProcessOneProj(mSurResol,*aCam,1.0,false,"",false,false);

        // ProcessOneProj only accumulates visibility inside aPPC (mSumRad);
        // ColorizePC pushes it into aPC_In as DegVis = visible-pixel fraction.
        aPPC.ColorizePC();

        // Snapshot visibility BEFORE DegVis is reused to store the ortho color.
        // DegVis==0  => point is occluded (never front-most in the nadir depth buffer).
        std::vector<bool> aVisible(aPC_In.NbPts(),false);
        for (size_t aKPt=0 ; aKPt<aPC_In.NbPts() ; aKPt++)
            aVisible[aKPt] = (aPC_In.GetDegVis(aKPt) > 0.0);

        // read ortho image and use it to colorate the cloud
        StdOut() << "Colorate in Ortho mode, using file " << mOrthoFile << "\n";

        cDataFileIm2D aFOrtho = cDataFileIm2D::Create(mOrthoFile,eForceGray::No);
        cDataIm2D<tU_INT1> aIDmOrtho(cPt2di(0,0),        
                                    aFOrtho.Sz());
        tREAL8 aTransform[6];
        std::vector<const cDataIm2D<tU_INT1>*> aVIms({&aIDmOrtho});

        cGdalApi::ReadWrite(cGdalApi::IoMode::Read,
                            aVIms,
                            aFOrtho,
                            cPt2di(0,0),
                            1.0,
                            cPixBox<2>(cPt2di(0,0),aFOrtho.Sz()),
                            aTransform);

        // fill cAffin2D
        cAffin2D<tREAL8> aTF(cPt2dr(aTransform[0],aTransform[3]),
                            cPt2dr(aTransform[1],0.0),
                            cPt2dr(0.0,aTransform[5]));

        // -------- Pass 1 : colour the VISIBLE points from their own ortho pixel.
        // A point's own ortho pixel is only its true colour when the point is on
        // the top surface. For an occluded point (e.g. ground under a tree) that
        // pixel is the canopy colour, so occluded points are handled in pass 2.
        std::vector<tREAL8> aVCol(aPC_In.NbPts(),0.0);

        // 2D spatial index of the visible points (storing their full 3D point)
        cBox2dr aBox2d = aPC_In.Box2d();
        cTiling<cVisPtCol> aTilVis(aBox2d,true,std::max<size_t>(1,aPC_In.NbPts()/20),-1);

        for(size_t aKPt=0 ; aKPt<aPC_In.NbPts() ; aKPt++)
        {
            cPt3dr aPGr = aPC_In.KthPt(aKPt);
            cPt2dr aPIm = aTF.Inverse(Proj(aPGr)); // ortho pixel
            if (aVisible[aKPt] && aIDmOrtho.InsideBL(aPIm))
            {
                aVCol[aKPt] = aIDmOrtho.GetVBL(aPIm) / 255.0;
                aTilVis.Add(cVisPtCol(aPGr,aVCol[aKPt]));
            }
        }

        // -------- Pass 2 : occluded points borrow the nearest visible colour,
        // giving a near-realistic ground colour instead of the canopy/black.
        // Search is done on the 2D index but ranked by 3D distance : since
        // 3Ddist >= XYdist, once the best 3D distance found is <= the current
        // search radius, no unseen point can be closer, so it is the 3D nearest.
        tREAL8 aGSD  = aPC_In.GroundSampling();       // ~ mean point spacing
        tREAL8 aRMax = aBox2d.Sz().x() + aBox2d.Sz().y();
        for(size_t aKPt=0 ; aKPt<aPC_In.NbPts() ; aKPt++)
        {
            if (aVisible[aKPt])
                continue;                             // already coloured in pass 1

            cPt3dr aP3  = aPC_In.KthPt(aKPt);
            cPt2dr aP2  = Proj(aP3);
            tREAL8 aBestD2 = 1e30;
            for (tREAL8 aRay=4*aGSD ; ; aRay*=2.0)
            {
                for (const cVisPtCol * aV : aTilVis.GetObjAtDist(aP2,aRay))
                {
                    tREAL8 aD2 = SqN2(aV->Pt3d()-aP3);   // full 3D distance
                    if (aD2 < aBestD2)
                    {
                        aBestD2 = aD2;
                        aVCol[aKPt] = aV->Col();
                    }
                }
                if (aBestD2 <= aRay*aRay)  // true 3D nearest guaranteed found
                    break;
                if (aRay > aRMax)          // safety : no visible point at all
                    break;
            }
        }
        for(size_t aKPt=0 ; aKPt<aPC_In.NbPts() ; aKPt++)
            aPC_In.SetDegVis(aKPt,aVCol[aKPt]);
   }
   else
   {
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
