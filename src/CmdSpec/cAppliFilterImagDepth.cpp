
// #include "include/MMVII_2Include_Serial_Tpl.h"
#include "MMVII_Matrix.h"


/** \file cMMVII_CalcSet.cpp
    \brief Command for set calculation

    This file contain the command that compute a  set of file from File/Regex
  It's also the first "real" command of MMVII, so an occasion for tuning a
  a lot of thing.

*/

namespace MMVII
{


typedef tREAL4             tValIm;
typedef cIm2D<tValIm>      tIm;
typedef cDataIm2D<tValIm>  tDIm;

typedef cIm2D<tU_INT1>     tImFlag;
typedef cDataIm2D<tU_INT1>     tDImFlag;



/**  Created initially to compute discontinuities in Corona image resulting
     from scanning artefact.

      Handle only vertical/horizontal line.  May evolve or not, not sure there is
     interest for much more complex tool (i.e polyline ...) in a command line
*/

class cAppli_FilterImDepth : public cMMVII_Appli
{
     public :


        cAppli_FilterImDepth(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli &);  ///< constructor
        int Exe() override;                                             ///< execute action
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override; ///< return spec of  mandatory args
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override; ///< return spec of optional args

     private :

        void EstimateStat();
        void MakeImFlag(tIm anIm);

        tIm  OneIterImageFilter(tIm);


        tREAL8  OneFilterVal(const cPt2di &,const tDIm&) const;


      // Read mandatory parameters
         std::string       mNameFile;    ///< Name of Image File
         tIm               mIm;
         tImFlag           mImFlag;
         cPt2di            mSzIm;


         tREAL8 mNbTestMed;

         tREAL8 mPercDev;
         tREAL8 mMulDev;
         tREAL8 mStdDev;
         tREAL8 mDivDev;
         int    mNbIter;

         
};

cCollecSpecArg2007 &  cAppli_FilterImDepth::ArgObl(cCollecSpecArg2007 & anArgObl)
{
   return
      anArgObl
          << Arg2007(mNameFile,"Name of file", {eTA2007::FileDirProj})
/*
          << Arg2007(mCenters,"Centers" )
          << Arg2007(mDP0,"Upper left offset")
          << Arg2007(mDP1,"Bottom right offset")
*/
      ;
}

cCollecSpecArg2007 & cAppli_FilterImDepth::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
   return
      anArgOpt
           << AOpt2007(mNbTestMed,"NbMed","Number of samples in median estimate",{eTA2007::HDV})
           << AOpt2007(mNbIter,"NbIter","Number of iteration",{eTA2007::HDV})
           << AOpt2007(mMulDev,"MulDev","multiplier M for W=1/(1+(R/M*Dev)^2)",{eTA2007::HDV})


      ;
}

cAppli_FilterImDepth::cAppli_FilterImDepth(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec) :
  cMMVII_Appli  (aVArgs,aSpec),
  mIm           (cPt2di(1,1)),
  mImFlag       (cPt2di(1,1)),
  mNbTestMed    (1e6),
  mPercDev      (0.5),
  mMulDev       (3.0),
  mNbIter       (2)
{

}



void cAppli_FilterImDepth::EstimateStat()
{
  //  tIm anIm = tIm::FromFile(mNameFile);
    tDIm & aDIm = mIm.DIm();

    cRect2  aRectInt = aDIm.Interior(1);
    const std::vector<cPt2di> & aV8 = Alloc8Neighbourhood();


    std::vector<tREAL8> aVDif;
    for (int aK=0 ; aK<mNbTestMed ; aK++)
    {
        cPt2di aPix1 = aRectInt.GeneratePointInside();
        cPt2di aPix2 = aPix1+aV8.at(RandUnif_N(8));
        tREAL8 aDif = std::abs(aDIm.GetV(aPix1)-aDIm.GetV(aPix2));

        aVDif.push_back(aDif);
    }

     mStdDev = NC_KthVal(aVDif,mPercDev) ;
     mDivDev = mStdDev * mMulDev;

    StdOut() << "MEDIANE*1000=" << mStdDev*1e3 << " Sz=" << aVDif.size() << "\n";
}


tIm  cAppli_FilterImDepth::OneIterImageFilter(tIm anIm)
{

    tIm aRes = anIm.Dup();

    for (const auto & aPix : aRes.DIm().Interior(1))
    {
        tREAL8 aNewV = OneFilterVal(aPix,anIm.DIm());
        aRes.DIm().SetV(aPix,aNewV);
    }

    return aRes;
}

tREAL8  cAppli_FilterImDepth::OneFilterVal(const cPt2di & aPix,const tDIm& aDIm) const
{
    cWeightAv<tREAL8> anAv;
    tValIm aV0 = aDIm.GetV(aPix);

    for (int aDx=-1 ; aDx<=1 ; aDx++)
    {
        for (int aDy=-1 ; aDy<=1 ; aDy++)
        {
             tREAL8 aV = aDIm.GetV(aPix+cPt2di(aDx,aDy));
             tREAL8 aWeigth = 1 /(1+ Square((aV-aV0)/mDivDev));
             anAv.Add(aWeigth,aV);
        }
    }

    return anAv.Average();
}

void   cAppli_FilterImDepth::MakeImFlag(tIm anIm)
{
    const std::vector<cPt2di> & aV8 = Alloc8Neighbourhood();
    const tDIm & aDIm = anIm.DIm();
     tDImFlag & aDImFlag = mImFlag.DIm();

    for (const auto & aPix1 : mImFlag.DIm())
    {
        for (int aK=0 ; aK<4 ; aK++)
        {
             cPt2di   aPix2 =  aPix1 + aV8[aK];
             if (aDIm.Inside(aPix2))
             {
                 tREAL8 aDif = aDIm.GetV(aPix1)-aDIm.GetV(aPix2);
                 if (std::abs(aDif) < mStdDev)
                 {
                     aDImFlag.GetReference_V(aPix1) |=  1<< aK;
                     aDImFlag.GetReference_V(aPix2) |=  1<< (aK+4);
                 }
             }
        }
    }

}


int cAppli_FilterImDepth::Exe()
{

    mIm = tIm::FromFile(mNameFile);
    mSzIm = mIm.DIm().Sz();
    mImFlag = tImFlag(mSzIm);

    EstimateStat();
    tIm aImF = mIm.Dup();

    for (int aK=0 ; aK<mNbIter ; aK++)
    {
        aImF = OneIterImageFilter(aImF);
    }
    MakeImFlag(aImF);

    aImF.DIm().ToFile("FilterInit_"+mNameFile);
    mImFlag.DIm().ToFile("Flage_"+mNameFile);


    return EXIT_SUCCESS;
}




tMMVII_UnikPApli Alloc_FilterImDepth(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec)
{
   return tMMVII_UnikPApli(new cAppli_FilterImDepth(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpecFilterImDepth
(
     "ImageFilterImDepth",
      Alloc_FilterImDepth,
      "Depth filtering in images",
      {eApF::ImProc},
      {eApDT::Image,eApDT::Console},
      {eApDT::Console},
      __FILE__
);

};

