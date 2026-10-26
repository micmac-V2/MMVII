#include "cMMVII_Appli.h"
#include "MMVII_PCSens.h"

namespace MMVII
{


/* ********************************************************************** */
/*                                                                        */
/*                   cExportSimCampPC                                     */
/*                                                                        */
/* ********************************************************************** */


cExportSimCampPC::cExportSimCampPC()
{
}

cExportSimCampPC::cExportSimCampPC(const cSensorCamPC * aCamPC,const cPt3dr & aPGround):
    mPGround    (aPGround),
    mPCam       (aCamPC->Pt_W2L(mPGround)),
    mPProj      (aCamPC->InternalCalib()->Dir_Proj()->Value(mPCam)),
    mPUnitDist  (aCamPC->InternalCalib()->Dir_Dist()->Value(mPProj)),
    mPIm        (aCamPC->InternalCalib()->MapPProj2Im().Value(mPUnitDist))
{
}
void AddData(const  cAuxAr2007 & anAux, cExportSimCampPC  &  anExp)
{
    AddData(cAuxAr2007("Ground",anAux),anExp.mPGround);
    AddData(cAuxAr2007("Cam",anAux),anExp.mPCam);
    AddData(cAuxAr2007("Proj",anAux),anExp.mPProj);
    AddData(cAuxAr2007("Dist",anAux),anExp.mPUnitDist);
    AddData(cAuxAr2007("Im",anAux),anExp.mPIm);
}


/* ********************************************************************** */
/*                                                                        */
/*                   cAlloc_ExportSimulMeasures                           */
/*                                                                        */
/* ********************************************************************** */


class cAlloc_ExportSimulMeasures : public cMMVII_Appli
{
    public:
        cAlloc_ExportSimulMeasures(const std::vector<std::string> & aVArgs, const cSpecMMVII_Appli & aSpec);

    private:
        // ================== Overiding cMMVII_Appli ========================
        int Exe() override;
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override;
         std::vector<cOneHelpSampleCmp>  Samples() const override;


         /// Make the export for central perspective camera
        void DoExportSimCampPC(const cSensorCamPC * aCamPC, const cSet2D3D &);
        /// Name of the file to save for a given image
        std::string NameFileSauv(const std::string & aNaemIm) const;

        //===================================================================
        //=============                    DATA PART            =============
        //===================================================================

        cPhotogrammetricProject  mPhProj;

        //=================  Mandatory args
        std::string              mSpecImIn;  ///< Patten image

        //==============  Optionnal args =======================

        eTypeSerial              mTypeSerial;   ///<  type of export (xml...)
        int                      mNbXY;         ///< number of sample in 1 direction
        int                      mNbDepth;      ///< number of layer in depth
        bool                     mOnlyPairs;    ///< if true , only generate Ground/Image even for central perspective camera
        tREAL8                   mEpsilonCheck; ///<

};

cAlloc_ExportSimulMeasures::cAlloc_ExportSimulMeasures(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec) :
    cMMVII_Appli    (aVArgs,aSpec),
    mPhProj         (*this),
    mNbXY           (6),
    mNbDepth        (2),
    mOnlyPairs      (false),
    mEpsilonCheck   (1e-3)
{
}

cCollecSpecArg2007 & cAlloc_ExportSimulMeasures::ArgObl(cCollecSpecArg2007 & anArgObl)
{
      return anArgObl
             << Arg2007(mSpecImIn,"Pattern/file for images",{{eTA2007::MPatFile,"0"},{eTA2007::FileDirProj}})
             << mPhProj.DPOrient().ArgDirInMand()
      ;
}




cCollecSpecArg2007 & cAlloc_ExportSimulMeasures::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
    return anArgOpt

               << cHeaderSectionArg("\"Standard\" args")
               << AOpt2007(mNbXY,"NbPlani","Number of samples in each dim X-Y",{eTA2007::HDV})
               << AOpt2007(mNbDepth,"NbDepth","Number of depth sampled",{eTA2007::HDV})
               << AOpt2007(mTypeSerial,"TypeS","Type of serial file (xml,json,csv..) default from user profile")

               << cHeaderSectionArg("\"Expert level\" args")
               << AOpt2007(mOnlyPairs,"OnlyGrIm","Export only Ground/Image (for Central Persp)",{eTA2007::HDV})
               << AOpt2007(mEpsilonCheck,"EpsilonCheck","Epsilon for checking with image proj with PC-Cam)",{eTA2007::HDV,eTA2007::Tuning})
            ;
}

std::vector<cOneHelpSampleCmp>  cAlloc_ExportSimulMeasures::Samples() const
{
    return
    {
        {"MMVII ExportSimulMesIm \"IMG_.*.JPG\" Adjust   TypeS=json"}
    };
}



void cAlloc_ExportSimulMeasures::DoExportSimCampPC(const cSensorCamPC * aCamPC, const cSet2D3D & aSetP32)
{
    std::vector<cExportSimCampPC> aVec;

    for (const auto & aP32 : aSetP32.Pairs() )
    {
        cExportSimCampPC anExp(aCamPC,aP32.mP3);
        /// little chack
        MMVII_INTERNAL_ASSERT_always(Norm2(anExp.mPIm - aP32.mP2)<mEpsilonCheck,"cAlloc_ExportSimulMeasures::DoExportSimCampPC");
        aVec.push_back(anExp);
    }
    SaveInFile(aVec,NameFileSauv(aCamPC->NameImage() ));
}

std::string cAlloc_ExportSimulMeasures::NameFileSauv(const std::string & aNaemIm) const
{
    return mPhProj.DirSimulAppli() + aNaemIm + "." + E2Str(mTypeSerial);
}


int cAlloc_ExportSimulMeasures::Exe()
{
    mPhProj.FinishInit();

    if (! IsInit(&mTypeSerial))
        mTypeSerial = TaggedDefSerial();

    //read the image pattern
    std::vector<std::string> aVecIm = VectMainSet(0);//interface to MainSet
    
    for (const std::string& aCImage : aVecIm)
    {
        cSensorImage * aSens = mPhProj.ReadSensor(aCImage,true,SVP::No);
        cSet2D3D  aSet32 =  aSens->SyntheticsCorresp3D2D (mNbXY,mNbDepth,1,2,true,0.05) ;

        if (aSens->IsSensorCamPC() && (!mOnlyPairs))
        {
            DoExportSimCampPC(aSens->GetSensorCamPC(),aSet32);
        }
        else
        {
            std::vector<cPair2D3D> aV32;
            for (const auto aPair : aSet32.Pairs())
                aV32.push_back(aPair);
            SaveInFile(aV32,NameFileSauv(aCImage ));
        }
    }

    return EXIT_SUCCESS;
}

tMMVII_UnikPApli Alloc_ExportSimulMeasures(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli & aSpec)
{
      return tMMVII_UnikPApli(new cAlloc_ExportSimulMeasures(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpec_ExportSimulMeasures
(
     "ExportSimulMesIm",
      Alloc_ExportSimulMeasures,
      "Export artificiel coorespondance 3D-2D (and more) for check/ reverse engeneering",
      {eApF::Ori},
      {eApDT::ObjMesInstr},
      {eApDT::Console},
      __FILE__
);

}
