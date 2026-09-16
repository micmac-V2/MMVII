#if MMVII_USE_LIBTORCH
#include "MMVII_2Include_Serial_Tpl.h"
#include "LearnDM.h"

//#include "include/MMVII_Tpl_Images.h"
//#include "include/MMVII_TplLayers3D.h"
#include <thread>
// included model cnn
#include "cCnnModelPredictor.h"
#include "MMVII_Tpl_Images.h"
#include "MMVII_TplLayers3D.h"

namespace MMVII
{

namespace  cNS_FillCubeCost
{

class cAppliFillCubeCost;
struct cOneModele
{
    public :
        typedef cIm2D<tREAL4>  tImRad;
        cOneModele
        (
            const std::string & aNameModele,
            cAppliFillCubeCost  & aAppliLearn
        );


        double ComputeCost(bool &Ok,const cPt2di & aPC1,const cPt2di & aPC2,int aDZ) const;
        void CalcCorrelExterneTerm(const cBox2di & aBoxInitIm1,int aPxMin,int aPxMax);
        void CalcCorrelExterneRecurs(const cBox2di & aBoxIm1);
        void CalcCorrelExterne();
        
        
       // ADDED METHODS FOR MVCNN
        void CalcCorrelMvCNN();
        void UseZInfZSup();

        cAppliFillCubeCost  * mAppli;
        std::string           mNameModele;
        bool                  mWithIntCorr;
        bool                  mWithExtCorr;
        bool                  mWithMMVIISimLearned;
        bool                  mWithStatModele;
        cHistoCarNDim         mModele;
        cPyr1ImLearnMatch *   mPyrL1;
        cPyr1ImLearnMatch *   mPyrL2;
        int                   mSzW;
        cPt2di                mPSzW;
        
        // instantiate a pointer to null CNN PREDICTOR
        aCnnModelPredictor * mCNNPredictor=nullptr;
        std::string mArchitecture ="";
        std::string mModelBinDir="";
        cPt2di              mCNNWin;
        
        // Networks architectures
        torch::jit::script::Module mMSNet;
        torch::jit::script::Module mDecisionNet;
        torch::jit::script::Module mMatcherNet;
        
};
static const std::string TheNameCorrel  = "MMVIICorrel";
static const std::string TheNameExtCorr = "ExternCorrel";
static const std::string TheNameSimLearned = "MMVIISimLearned";

static const std::string TheUnetMlpCubeMatcher="UnetMlp";
//.....................................................

class cAppliFillCubeCost : public cAppliLearningMatch
{
     public :
        typedef tINT2                      tElemZ;
        typedef cIm2D<tElemZ>              tImZ;
        typedef cDataIm2D<tElemZ>          tDataImZ;
        typedef cIm2D<tREAL4>              tImRad;
        typedef cDataIm2D<tREAL4>          tDataImRad;
        typedef cLayer3D<float,tElemZ>     tLayerCor;
        typedef cIm2D<tREAL4>              tImPx;
        typedef cDataIm2D<tU_INT1>         tDataImMasq;
        typedef cDataIm2D<tREAL4>          tDataImPx;
        typedef cIm2D<tREAL4>              tImFiltred;
        typedef cDataIm2D<tREAL4>          tDataImF;
        typedef cGaussianPyramid<tREAL4>   tPyr;
        typedef std::shared_ptr<tPyr>      tSP_Pyr;
        typedef cPyr1ImLearnMatch *        tPtrPyr1ILM;

        cAppliFillCubeCost(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec);
        double ComputCorrel(const cPt2di & aPI1,const cPt2dr & aPI2,int mSzW) const;
        void OneSliceSimil(
                    cPt2di aPx,
                    std::vector<std::pair<int,std::pair<torch::Tensor,int>>> * AllSimils,
                    std::vector<cOneModele*> VMods,
                    torch::Tensor & EmbedL,
                    torch::Tensor & EmbedR,
                    tDataImZ & aDZMIN,
                    tDataImZ & aDZMAX,
                    int FeatSize);
        void ExtractSurface(torch::Tensor & aSimilTensor, 
                            std::string & aNamePaxOut, 
                            tREAL8 aStepZ);
        const tDataImRad & DI1() {return *mDI1;}
        const tDataImRad & DI2() {return *mDI2;}
        
        const tImRad & IMNorm1() {return mImNorm1;}
        const tImRad & IMNorm2() {return mImNorm2;}
        const tDataImRad & NDI1() {return *mDINorm1;}
        const tDataImRad & NDI2() {return *mDINorm2;}

        cPyr1ImLearnMatch * PyrL1 () {return PyrL(mPyrL1,mBoxGlob1,mNameI1);}
        cPyr1ImLearnMatch * PyrL2 () {return PyrL(mPyrL2,mBoxGlob2,mNameI2);}
        tREAL8     StepZ() const {return mStepZ;}
        bool Ok1(int aX) const {return Ok(aX,mVOk1);}
        bool Ok2(int aX) const {return Ok(aX,mVOk2);}
        const cAimePCar & PC1(int aX) const {return mVPC1.at(aX);}
        const cAimePCar & PC2(int aX) const {return mVPC2.at(aX);}
        const cBox2di  & BoxGlob1() const {return mBoxGlob1;}  ///< Accessor
        const cBox2di  & BoxGlob2() const {return mBoxGlob2;}  ///< Accessor
        const std::string   & NameI1() const {return mNameI1;}  ///< Accessor
        const std::string   & NameI2() const {return mNameI2;}  ///< Accessor
        
        const std::string  & NameArch() const {return mModelArchitecture;} // ACCESSOR
        const std::string  & NameDirModel() const {return mModelBinaries;} // ACCESSOR

        cBox2di BoxFile1() const {return cDataFileIm2D::Create(mNameI1,eForceGray::Yes);}
        cBox2di BoxFile2() const {return cDataFileIm2D::Create(mNameI2,eForceGray::Yes);}

        int  SzW() const {return mSzW;}
        bool InterpolLearn() const {return mInterpolLearn;}
        double ExpLearn()    const {return mExpLearn; }
        double FactLearn()   const {return mFactLearn; }
        const cFilterPCar  & FPC() const {return mFPC;}  ///< Used to compute Pts

        const tImZ  & ImZMin() {return  mImZMin;}
        const tImZ  & ImZMax() {return  mImZMax;}
        void MakeNormalizedIm();
        bool mUseCuda=false;
        bool mUsePredicNet=false;

        // -------------- Internal variables -------------------
     private :

        cPyr1ImLearnMatch * PyrL (tPtrPyr1ILM & aPtrPyr,const cBox2di & aBoxI,const std::string & aNameIm)
        {
           if (aPtrPyr==nullptr)
              aPtrPyr = new cPyr1ImLearnMatch(aBoxI,aBoxI,aNameIm,*this,mFPC,false);
           return aPtrPyr;
        }
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override ;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override ;

        int Exe() override;
        //void PushCost(double aCost);
        
        bool Ok(int aX,const std::vector<bool> &  aV) const
        {
            return (aX>=0) && (aX<int(aV.size())) && (aV.at(aX)) ;
        }

        void MakeLinePC(int aYLoc,bool Im1);
        // -------------- Mandatory args -------------------
        std::string   mNameI1;
        std::string   mNameI2;
        std::string   mNamePax;
        std::string   mNameModele;
        cPt2di        mP0Z;  // Pt corresponding in Im1 to (0,0)
        cBox2di       mBoxGlob1;  // Box to Load, taking into account siwe effect
        cBox2di       mBoxGlob2;
        std::string   mNamePost;

        // -------------- Optionnal args -------------------
        tREAL8        mStepZ;
        bool          mCmpCorLearn; //  Create a comparison between Correl & Learn
        bool          mInterpolLearn; // Use interpolation mode for learned cost
        double        mExpLearn; // Exposant to adapt learned cost
        double        mFactLearn; // Factor to adapt learned cost
        std::string   mNameCmpModele;
        int           mSzW;
        tREAL4 mThresholdSim;
    
        std::string mModelBinaries;
        std::string mModelArchitecture=TheUnetMlpCubeMatcher;

        // -------------- Internal variables -------------------
        
        std::string StdName(const std::string & aPre,const std::string & aPost);

        int         mNbCmpCL;
            cIm2D<tREAL8>  mImCmp;
        std::string mNameZMin;
        std::string mNameZMax;
        std::string mNameCube;
        //cMMVII_Ofs* mFileCube;

        tImZ        mImZMin;
        tImZ        mImZMax;
        tImRad      mIm1;
        tDataImRad  *mDI1;
        tImRad      mIm2;
        tDataImRad  *mDI2;

        // Normalized images, in radiometry, avearge 0, std dev 1.0
        tImRad      mImNorm1;
        tDataImRad  *mDINorm1;
        tImRad      mImNorm2;
        tDataImRad  *mDINorm2;
        tLayerCor   mLayerCor;

        double      ToCmpCost(double aCost) const;

        cPyr1ImLearnMatch * mPyrL1;
        cPyr1ImLearnMatch * mPyrL2;
        cFilterPCar  mFPC;  ///< Used to compute Pts
        std::vector<bool>         mVOk1;
        std::vector<cAimePCar>    mVPC1;
        std::vector<bool>         mVOk2;
        std::vector<cAimePCar>    mVPC2;
};


/* *************************************************** */
/*                                                     */
/*                   cOneModele                        */
/*                                                     */
/* *************************************************** */


cOneModele::cOneModele
(
    const std::string & aNameModele,
    cAppliFillCubeCost  & aAppliLearn
    //aCnnModelPredictor & aPredicor
) :
   mAppli          (&aAppliLearn),
   //**************************,
   mNameModele     (aNameModele),
   mWithIntCorr    (mNameModele==TheNameCorrel),
   mWithExtCorr    (mNameModele==TheNameExtCorr),
   /***********************************************************************/
   mWithMMVIISimLearned (mNameModele==TheNameSimLearned),
   /***********************************************************************/
   mWithStatModele (! (mWithIntCorr || mWithExtCorr || mWithMMVIISimLearned)),
   mPyrL1          (nullptr),
   mPyrL2          (nullptr),
   mSzW            (mAppli->SzW()),
   mPSzW           (mSzW,mSzW),
   mCNNWin          (0,0)
{
    if (mWithStatModele)
    {
       ReadFromFile(mModele,mNameModele);
       mPyrL1 = mAppli->PyrL1 ();
       mPyrL2 = mAppli->PyrL2 ();
    }
    else if (mWithExtCorr)
    {
         CalcCorrelExterne();
    }
    else if (mWithMMVIISimLearned)
    {
        mArchitecture=mAppli->NameArch();
        mModelBinDir =mAppli->NameDirModel();
        
        MMVII_INTERNAL_ASSERT_strong(mArchitecture!="",
            "The network architecture should be specified :  "+ TheUnetMlpCubeMatcher+" !");
        MMVII_INTERNAL_ASSERT_strong(mModelBinDir!=""," Model params dir must be specified ! ");

        if (mArchitecture==TheUnetMlpCubeMatcher)
            {
                mCNNPredictor = new aCnnModelPredictor(TheUnetMlpCubeMatcher,
                                            mModelBinDir,
                                            mAppli->mUseCuda);
                mCNNPredictor->PopulateModelFeatures(mMSNet);
                mCNNWin=cPt2di(1,1);
                if (mAppli->mUsePredicNet)
                {
                  mCNNPredictor->PopulateModelDecision(mDecisionNet);
                }
            }
    }
}
void cOneModele::CalcCorrelExterneTerm(const cBox2di & aBoxInitIm1,int aPxMin,int aPxMax)
{
     // cBox2di aBoxDil = aBoxInitIm1.Inter
}

void cOneModele::CalcCorrelExterneRecurs(const cBox2di & aBoxIm1)
{
     cPt2di aP0Glob = mAppli->BoxGlob1().P0();
     int aMinPxMin = 1e6;
     int aMaxPxMax = -1e6;
     int aTotPx = 0;
     int aNbPx = 0;
     const cAppliFillCubeCost::tDataImZ & aDIZMin = mAppli->ImZMin().DIm();
     const cAppliFillCubeCost::tDataImZ & aDIZMax = mAppli->ImZMax().DIm();

     cRect2 aR2(aBoxIm1.P0(),aBoxIm1.P1());
     for (const auto & aP : aR2)
     {
          cPt2di aPLoc = aP-aP0Glob;
          int aPxMin = aDIZMin.GetV(aPLoc);
          int aPxMax = aDIZMax.GetV(aPLoc);
          UpdateMin(aMinPxMin,aPxMin);
          UpdateMax(aMaxPxMax,aPxMax);
          aTotPx += aPxMax - aPxMin;
          aNbPx++;
     }

     double aAvgPx = aTotPx / double(aNbPx);
     int   aIntervPx = (aMaxPxMax-aMinPxMin);

     bool   isTerminal = aIntervPx < 2 * aAvgPx;

     if ((aIntervPx*aIntervPx) > 1e8)
     {
        isTerminal = false;
     }

     {
        int aSzMin = MinAbsCoord(aBoxIm1.Sz());
        if (aSzMin<200)
           isTerminal =true;
     }


     if (isTerminal)
     {
         CalcCorrelExterneTerm(aBoxIm1,aMinPxMin,aMaxPxMax);
     }
     else
     {
        cPt2di aP0 = aBoxIm1.P0();
        cPt2di aP1 = aBoxIm1.P1();
        std::vector<int> aVx{aP0.x(),(aP0.x()+aP1.x())/2,aP1.x()};
        std::vector<int> aVy{aP0.y(),(aP0.y()+aP1.y())/2,aP1.y()};
        for (int aKx=0 ; aKx<2 ; aKx++)
        {
             for (int aKy=0 ; aKy<2 ; aKy++)
             {
                  cPt2di aQ0(aVx.at(aKx),aVy.at(aKy));
                  cPt2di aQ1(aVx.at(aKx+1),aVy.at(aKy+1));
                  CalcCorrelExterneRecurs(cBox2di(aQ0,aQ1));
             }
        }
     }
}

void cOneModele::CalcCorrelExterne()
{
   mAppli->MakeNormalizedIm();
   CalcCorrelExterneRecurs(mAppli->BoxGlob1());
}

void cOneModele::CalcCorrelMvCNN()
{
   mAppli->MakeNormalizedIm();
}


void UseZInfZSup()
{
    std::cout<<"Handling Bounding Layers "<<std::endl;
}

double cOneModele::ComputeCost(bool & Ok,const cPt2di & aPC1,const cPt2di & aPC20,int aDZ) const
{
    Ok = false;
    double aCost= 1.0;
    if (mWithStatModele)
    {
       int aX1 =  aPC1.x();
       int aX2 =  aPC20.x() + aDZ;
       aCost = 0.5;
       if (mAppli->Ok1(aX1) && mAppli->Ok2(aX2))
       {
           cVecCaracMatch aVCM(*mPyrL1,*mPyrL2,mAppli->PC1(aX1),mAppli->PC2(aX2));
           aCost = 1-mModele.HomologyLikelihood(aVCM,mAppli->InterpolLearn());
           aCost = mAppli->FactLearn() * pow(std::max(0.0,aCost),mAppli->ExpLearn());
           Ok = true;
       };
    }
    else if (mWithIntCorr)
    {
        cPt2dr aPC2Z(aPC20.x()+aDZ*mAppli->StepZ(),aPC20.y());
        double aCorrel = 0.0;

        if (WindInside4BL(mAppli->DI1(),aPC1,mPSzW) && WindInside4BL(mAppli->DI2(),aPC2Z,mPSzW))
        {
            aCorrel = mAppli->ComputCorrel(aPC1,aPC2Z,mSzW);
            Ok = true;
        }
        aCost=(1-aCorrel)/2.0;
    }
    else if (mWithExtCorr)
    {
        //COMPUTE CORREL USING STATISTICAL LEARNING
        
    }
    return aCost;
}

/* *************************************************** */
/*                                                     */
/*              cAppliFillCubeCost                     */
/*                                                     */
/* *************************************************** */

cAppliFillCubeCost::cAppliFillCubeCost(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec) :
   cAppliLearningMatch  (aVArgs,aSpec),
   mBoxGlob1            (cBox2di::Empty()),
   mBoxGlob2            (cBox2di::Empty()),
   mStepZ               (1.0),
   mCmpCorLearn         (true),
   mInterpolLearn       (true),
   mExpLearn            (0.5),
   mFactLearn           (0.33333),
   mSzW                 (3),
   mThresholdSim        (0.5),
   mNbCmpCL             (200),
   mImCmp               (cPt2di(mNbCmpCL+2,mNbCmpCL+2),nullptr,eModeInitImage::eMIA_Null),
   //mFileCube            (nullptr),
   mImZMin              (cPt2di(1,1)),
   mImZMax              (cPt2di(1,1)),
   mIm1                 (cPt2di(1,1)),
   mDI1                 (nullptr),
   mIm2                 (cPt2di(1,1)),
   mDI2                 (nullptr),
   mImNorm1             (cPt2di(1,1)),
   mDINorm1             (nullptr),
   mImNorm2             (cPt2di(1,1)),
   mDINorm2             (nullptr),
   mLayerCor            (tLayerCor::Empty()),
   mPyrL1               (nullptr),
   mPyrL2               (nullptr),
   mFPC                 (false)
{
    mFPC.FinishAC();
    mFPC.Check();
}

void cAppliFillCubeCost::MakeNormalizedIm()
{
    if (mDINorm1!= nullptr) return;

    mImNorm1 = NormalizedAvgDev(mIm1,1e-4);
    mDINorm1 = &(mImNorm1.DIm());

    mImNorm2 = NormalizedAvgDev(mIm2,1e-4);
    mDINorm2 = &(mImNorm2.DIm());

    //mLayerCor  = tLayerCor(mImZMin,mImZMax);
}


double cAppliFillCubeCost::ToCmpCost(double aCost) const
{
   return mNbCmpCL * std::max(0.0,std::min(1.0,aCost));
}

cCollecSpecArg2007 & cAppliFillCubeCost::ArgObl(cCollecSpecArg2007 & anArgObl)
{
 return
      anArgObl
          <<   Arg2007(mNameI1,"Name of first image")
          <<   Arg2007(mNameI2,"Name of second image")
          <<   Arg2007(mNamePax, "Name of Parallax image")
          <<   Arg2007(mNameModele,"Name for modele : .*dmp|MMVIICorrel|MMVIISimLearned")
          <<   Arg2007(mP0Z,"Origin in first image")
          <<   Arg2007(mBoxGlob1,"Box to read 4 Im1")
          <<   Arg2007(mBoxGlob2,"Box to read 4 Im2")
          <<   Arg2007(mNamePost,"Post fix for other names (ZMin,ZMax,Cube)")
   ;
}

cCollecSpecArg2007 & cAppliFillCubeCost::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
   return anArgOpt
          << AOpt2007(mStepZ, "StepZ","Step for paralax",{eTA2007::HDV})
          << AOpt2007(mNameCmpModele, "ModCmp","Modele for Comparison")
          << AOpt2007(mSzW, "SzW","Size for windows to match",{eTA2007::HDV})
          << AOpt2007(mModelBinaries,"CNNParams" ,"Model Directory : Contient des fichiers binaires *.bin")
          << AOpt2007(mModelArchitecture,"CNNArch" ,"Modek architecture :  " + TheUnetMlpCubeMatcher)
          << AOpt2007(mUseCuda,"UseCuda","USE CUDA TO LAUNCH MODELS")
          << AOpt2007(mUsePredicNet,"UsePredicNet","Use the mlp to compute learnt similarities")
   ;
}

std::string cAppliFillCubeCost::StdName(const std::string & aPre,const std::string & aPost)
{
        return aPre + "_" + mNamePost + "." + aPost;
}

double cAppliFillCubeCost::ComputCorrel(const cPt2di & aPCI1,const cPt2dr & aPCI2,int aSzW) const
{
   cMatIner2Var<tREAL4> aMat;

   for (int aDx=-aSzW ; aDx<=aSzW  ; aDx++)
   {
       for (int aDy=-aSzW ; aDy<=aSzW  ; aDy++)
       {
            aMat.Add
            (
                mDI1->GetV  (aPCI1+cPt2di(aDx,aDy)),
                mDI2->GetVBL(aPCI2+cPt2dr(aDx,aDy))
            );
       }
   }

   return aMat.Correl();
}

/*void cAppliFillCubeCost::PushCost(double aCost)
{
   tU_INT2 aICost = round_ni(1e4*(std::max(0.0,std::min(1.0,aCost))));
   mFileCube->Write(aICost);
}*/


void cAppliFillCubeCost::MakeLinePC(int aYLoc,bool Im1)
{
   if (mPyrL1==nullptr)
      return;

   MMVII_INTERNAL_ASSERT_strong(mStepZ==1.0,"For now do not handle StepZ!=1 with model");

   std::vector<bool>     & aVOK  = Im1 ? mVOk1      : mVOk2;
   std::vector<cAimePCar>& aVPC  = Im1 ? mVPC1      : mVPC2;
   const cBox2di         & aBox  = Im1 ? mBoxGlob1  : mBoxGlob2;
   cPyr1ImLearnMatch     & aPyrL = Im1 ? *mPyrL1    : *mPyrL2;

   aVOK.clear();
   aVPC.clear();

   for (int aX=aBox.P0().x() ; aX<=aBox.P1().x()  ; aX++)
   {
       cPt2di aPAbs (aX,aYLoc+mP0Z.y());
       cPt2di aPLoc = aPAbs - aBox.P0();
       aVOK.push_back(aPyrL.CalculAimeDesc(ToR(aPLoc)));
       aVPC.push_back(aPyrL.DupLPIm());
   }
}

void cAppliFillCubeCost::OneSliceSimil(
                                    cPt2di aPx,
                                    std::vector<std::pair<int,std::pair<torch::Tensor,int>>> * AllSimils,
                                    std::vector<cOneModele*> VMods,
                                    torch::Tensor & EmbedL,
                                    torch::Tensor & EmbedR,
                                    tDataImZ & aDZMIN,
                                    tDataImZ & aDZMAX,
                                    int FS
                                )
{
    cPt2di aPAbs = aPx + mP0Z;
    cPt2di aPC1  = aPAbs-mBoxGlob1.P0();
    cPt2di aPC20 = aPAbs-mBoxGlob2.P0();
    using namespace torch::indexing;
    auto aVecL=EmbedL.index({0,Slice(0,FS,1),aPC1.y(),aPC1.x()}).unsqueeze(0).unsqueeze(2).unsqueeze(3);
    int liminf=round_ni(aPC20.x()+aDZMIN.GetV(aPx)*this->StepZ());
    int limsup=round_ni(aPC20.x()+aDZMAX.GetV(aPx)*this->StepZ());
    if (liminf<0 && limsup>0)
    {
            liminf=0;
    }
    auto aVecR=EmbedR.index(
                {0,
                 Slice(0,FS,1),
                 aPC20.y(),
                 Slice(liminf,
                       limsup,
                       1)
                }
                );
    aVecR=aVecR.unsqueeze(0).unsqueeze(2);
    auto aRepeatedVecL=at::repeat_interleave(aVecL,aVecR.size(-1),-1);
    auto SimSlice=VMods.at(0)->mCNNPredictor->PredictDecisionNet(VMods.at(0)->mDecisionNet,aRepeatedVecL,aVecR);
    AllSimils->push_back(std::make_pair(aPx.x(),std::make_pair(SimSlice,liminf)));
}

void cAppliFillCubeCost::ExtractSurface(torch::Tensor & aSimilTensor, 
                                        std::string & aNamePaxOut, 
                                        tREAL8 aStepZ)
{
    // aSimilTensor is the correlation cube [H, W1, W2] (height, width of Im1, width of Im2).
    // For each left pixel (H,W1) the disparity is the argmax over the W2 axis (dim 2),
    // which yields a [H, W1] map matching the layout of Im1.
    auto aDispArgmax = at::argmax(aSimilTensor, 2);

    // argmax returns an int64 tensor : convert to float32, move to CPU and make it
    // contiguous so that its linear layout matches the MMVII image buffer.
    auto aDispArgmaxCPU = aDispArgmax.to(torch::kCPU)
                                     .to(torch::kFloat32)
                                     .contiguous();

    // Similarity at argmax
    auto aSimAtMax = aSimilTensor.gather(2, aDispArgmax.unsqueeze(2)).squeeze(2)
                                     .to(torch::kCPU)
                                     .to(torch::kFloat32)
                                     .contiguous();

    // Allocate an image that owns its buffer (do NOT alias the tensor data), then copy.
    cIm2D<tREAL4> aImDisp(mDI1->Sz(), nullptr, eModeInitImage::eMIA_Null);
    tDataImRad & aDID = aImDisp.DIm();

    // Allocate a similarity image that owns its buffer (do NOT alias the tensor data), then copy.
    cIm2D<tREAL4> aImSim(mDI1->Sz(), nullptr, eModeInitImage::eMIA_Null);
    tDataImRad & aDIS = aImSim.DIm();

    MMVII_INTERNAL_ASSERT_always(
        aDispArgmaxCPU.numel() == aDID.NbElem(),
        "ExtractSurface: disparity tensor size does not match Im1 size");

    const tREAL4 * aSrc = aDispArgmaxCPU.data_ptr<tREAL4>();
    tREAL4 *       aDst = aDID.RawDataLin();

    const tREAL4 * aSrcSim = aSimAtMax.data_ptr<tREAL4>();
    tREAL4 *       aDstSim = aDIS.RawDataLin();

    for (int aK=0 ; aK<aDID.NbElem() ; aK++)
    {
        // disparity is d = x2 - x1, where x2 is the index of the max similarity in Im2 for pixel x1 in Im1.
        aDst[aK] = aSrc[aK] - (aK % aDID.Sz().x());
        aDstSim[aK] = aSrcSim[aK];
    }

    // if StepZ is not 1, we need to scale the disparity values accordingly
    if (aStepZ!=1.0)
    {
        for (const auto & aP: aDID)
        {
            aDID.SetV(aP, aDID.GetV(aP) * aStepZ);
        }
    }

    cIm2D<tU_INT1> aImMask(mDI1->Sz(), nullptr, eModeInitImage::eMIA_Null);
    tDataImMasq & aDM = aImMask.DIm();

    for (const auto & aP: aDIS)
    {
        aDM.SetV(aP, (aDIS.GetV(aP) >= mThresholdSim) ? 1 : 0);
    }
    // Write the disparity image to a file
    aDID.ToFile(aNamePaxOut);
    aDM.ToFile(aNamePaxOut + "_AutoMask.tif");
    aDIS.ToFile(aNamePaxOut + "_AutoSim.tif");
}

int  cAppliFillCubeCost::Exe()
{
#ifdef _WIN32
  if (mUseCuda) LoadLibraryA("torch_cuda.dll");
#endif
   // Compute names

   
   //mNameZMin = StdName("ZMin","tif");
   //mNameZMax = StdName("ZMax","tif");
   //mNameCube = StdName("MatchingCube","data");
   //  Read images
   //mImZMin = tImZ::FromFile(mNameZMin);
   //tDataImZ & aDZMin = mImZMin.DIm();
   //mImZMax = tImZ::FromFile(mNameZMax);
   //tDataImZ & aDZMax = mImZMax.DIm();

   mIm1 = tImRad::FromFile(mNameI1);
   mDI1 = &(mIm1.DIm());

   mIm2 = tImRad::FromFile(mNameI2);
   mDI2 = &(mIm2.DIm());
   //mFileCube = new cMMVII_Ofs(mNameCube, eFileModeOut::CreateBinary);

   mCmpCorLearn = IsInit(&mNameCmpModele);
   std::vector<cOneModele*> aVMods;
   aVMods.push_back(new cOneModele(mNameModele,*this));
   if (mCmpCorLearn)
       aVMods.push_back(new cOneModele(mNameCmpModele,*this));

    if (aVMods.at(0)->mWithMMVIISimLearned)
    {

        torch::Device device(mUseCuda ? torch::kCUDA : torch::kCPU);
        cPt2di aSzL = DI1().Sz();
        cPt2di aSzR = DI2().Sz();
        //int FeatSize;
        //if (aVMods.at(0)->mArchitecture==TheUnetMlpCubeMatcher) FeatSize=64 ;
        torch::Tensor LREmbeddingsL,LREmbeddingsR;
        if (aVMods.at(0)->mArchitecture==TheUnetMlpCubeMatcher)
             {
               LREmbeddingsL=aVMods.at(0)->mCNNPredictor->PredictMSNetTileFeatures(aVMods.at(0)->mMSNet,
                                                                            this->mIm1,
                                                                            aSzL);
               LREmbeddingsR=aVMods.at(0)->mCNNPredictor->PredictMSNetTileFeatures(aVMods.at(0)->mMSNet,
                                                                            this->mIm2,
                                                                            aSzR);
                cAutoTimerSegm aTEmbeddingsCreate(TimeSegm(),
                                        "DM4FillCubeCost::Embeddings Calcultion"); 
             }

        StdOut()  <<" EMBEDDING TENSOR SIZE LEFT   "<<LREmbeddingsL.sizes()<<"\n";
        StdOut()  <<" EMBEDDING TENSOR SIZE RIGHT  "<<LREmbeddingsR.sizes()<<"\n";
        
        //cPt2di aPix;
        /* Discard Mlp prediction for now, it is not used in the current version of the code
        if (mUsePredicNet)
        {
            for (aPix.y()=0 ; aPix.y()<aSzL.y() ; aPix.y()++)
            {
                int  aMinZmin=1e8;
                int  aMaxZmax=-1e8;
                for (aPix.x()=0 ; aPix.x()<aSzL.x() ; aPix.x()++)
                {
                        if (aDZMin.GetV(aPix)<aMinZmin) aMinZmin=aDZMin.GetV(aPix);
                        if (aDZMax.GetV(aPix)>aMaxZmax) aMaxZmax=aDZMax.GetV(aPix);
                }
                torch::Device TheAvailDevice(mUseCuda ? torch::kCUDA : torch::kCPU);
                MMVII_INTERNAL_ASSERT_always(aMaxZmax-aMinZmin, "PAX INTERVAL NULL !");
                using namespace torch::indexing;

                auto aSlcL=LREmbeddingsL.index({Slice(0,FeatSize,1),
                                                Slice(aPix.y()+mP0Z.y()-mBoxGlob1.P0().y()
                                                ,aPix.y()+mP0Z.y()-mBoxGlob1.P0().y()+1,1),
                                                Slice(mP0Z.x()-mBoxGlob1.P0().x(),
                                                aSzL.x()+mP0Z.x()-mBoxGlob1.P0().x(),1)}); // FeatSize,1,W

                int Intervalle_DPAX=round_ni(aMaxZmax-aMinZmin);// /this->StepZ());
                torch::Tensor CUBE= torch::ones({2*FeatSize,Intervalle_DPAX,1,aSzL.x()},
                        torch::TensorOptions().dtype(torch::kFloat32).device(TheAvailDevice)).mul(1.0);

                MMVII_INTERNAL_ASSERT_always(Intervalle_DPAX==aMaxZmax-aMinZmin, "ISSUE WITH PAWX INTERVAL ");

                for (int dd=0;dd<Intervalle_DPAX;dd++)
                    {
                        // Get relevant right features

                        CUBE.index({Slice(0,FeatSize,1),dd,Slice(0,None,1),Slice(0,None,1)}).copy_(aSlcL);

                        int lim_inf=round_ni(mP0Z.x()-mBoxGlob2.P0().x()+dd+aMinZmin);
                        int lim_sup=round_ni(mP0Z.x()-mBoxGlob2.P0().x()+dd+aMinZmin+aSzL.x());
                        if ((lim_inf<0 && lim_sup<0) || (lim_inf>aSzR.x() && lim_sup>aSzR.x()) )
                            {

                            }
                        else if(lim_inf<0 && lim_sup>aSzR.x())
                            {
                            auto aSlcR=LREmbeddingsR.index({Slice(0,FeatSize,1),
                                                                Slice(aPix.y()+mP0Z.y()-mBoxGlob2.P0().y()
                                                                ,aPix.y()+mP0Z.y()-mBoxGlob2.P0().y()+1,1),
                                                                Slice(0,aSzR.x(),1)}); // FEatSize,1,W
                            CUBE.index({Slice(FeatSize,2*FeatSize,1),dd,
                                        Slice(0,None,1),Slice(-lim_inf,aSzR.x()-lim_inf,1)}).copy_(aSlcR);

                            }
                        else if (lim_inf<0)
                            {
                                auto aSlcR=LREmbeddingsR.index({Slice(0,FeatSize,1),
                                                                Slice(aPix.y()+mP0Z.y()-mBoxGlob2.P0().y()
                                                                ,aPix.y()+mP0Z.y()-mBoxGlob2.P0().y()+1,1),
                                                                Slice(0,lim_sup,1)}); // FEatSize,1,W

                                CUBE.index({Slice(FeatSize,2*FeatSize,1),dd,Slice(0,None,1),Slice(-lim_inf,CUBE.size(3),1)}).copy_(aSlcR);
                            }
                        else if (lim_sup>aSzR.x())
                            {
                                auto aSlcR=LREmbeddingsR.index({Slice(0,FeatSize,1),
                                                                    Slice(aPix.y()+mP0Z.y()-mBoxGlob2.P0().y()
                                                                    ,aPix.y()+mP0Z.y()-mBoxGlob2.P0().y()+1,1),
                                                                    Slice(lim_inf,aSzR.x(),1)}); // FEatSize,1,W
                                CUBE.index({Slice(FeatSize,2*FeatSize,1),dd,Slice(0,None,1),Slice(0,aSlcR.size(2),1)}).copy_(aSlcR);
                            }
                        else
                            {
                                // TAKE USUAL RANGES
                                auto aSlcR=LREmbeddingsR.index({Slice(0,FeatSize,1),
                                                                    Slice(aPix.y()+mP0Z.y()-mBoxGlob2.P0().y()
                                                                    ,aPix.y()+mP0Z.y()-mBoxGlob2.P0().y()+1,1),
                                                                    Slice(lim_inf,lim_sup,1)}); // FEatSize,1,W
                                CUBE.index({Slice(FeatSize,2*FeatSize,1),
                                            dd,
                                            Slice(0,None,1),Slice(0,None,1)}).copy_(aSlcR);
                            }
                    }

                auto aSimilCUBE=aVMods.at(0)->mCNNPredictor->PredictONCUBE(aVMods.at(0)->mDecisionNet, CUBE);
                CUBE.resize_(at::IntArrayRef{0});

                if (aSimilCUBE.dim()==4)
                {
                    aSimilCUBE.squeeze();
                }
                // Fill cube cost
                for (aPix.x()=0 ; aPix.x()<aSzL.x() ; aPix.x()++)
                {
                        cPt2di aPAbs = aPix + mP0Z;
                        cPt2di aPC1  = aPAbs-mBoxGlob1.P0();
                        cPt2di aPC20 = aPAbs-mBoxGlob2.P0();
                        for (int aDz=aDZMin.GetV(aPix) ; aDz<aDZMax.GetV(aPix) ; aDz++)
                        {
                            double aTabCost[2]={1.0,1.0};
                            bool   aTabOk[2]={false,false};
                            cPt2di aPC2Z(round_ni(aPC20.x()+aDz*this->StepZ()),aPC20.y());  // INTEG FOR NOW
                            for (int aK=0 ; aK<int(aVMods.size()) ; aK++)
                                    {
                                        bool IsInside=WindInside4BL(this->DI1(),aPC1,aVMods[aK]->mCNNWin) && WindInside4BL(this->DI2(),aPC2Z,aVMods[aK]->mCNNWin);
                                        if(IsInside)
                                        {
                                            //auto aVecR=LREmbeddingsR.slice(2,aPC2Z.y(),aPC2Z.y()+1).slice(3,aPC2Z.x(),aPC2Z.x()+1);
                                            using namespace torch::indexing;
                                            auto aSim=aSimilCUBE.index({Slice(aDz-aMinZmin,aDz-aMinZmin+1,1),
                                                                            Slice(aPix.x(),aPix.x()+1,1)});
                                            MMVII_INTERNAL_ASSERT_always(
                                                aSim.item<float>()<=1.0 && aSim.item<float>()>=0,
                                                "Similarity values issue not in bound 0 ,1 "
                                                );
                                            aTabCost[aK] =(1-(double)aSim.item<float>());
                                            aTabOk[aK]=true;
                                        }
                                    }
                            //PushCost(aTabCost[0]);
                            if (mCmpCorLearn && aTabOk[0] && aTabOk[1])
                            {
                                    double aC0 = ToCmpCost(aTabCost[0]);
                                    double aC1 = ToCmpCost(aTabCost[1]);
                                    mImCmp.DIm().AddVBL(cPt2dr(aC1,aC0),1.0);
                            }
                        }
                }

            }
        }

        */
        //else
        {
            // calculer le min sur les nappes
            // JUST IGNORE FOR NOW BOUNDING LAYERS int  aMinZmin=1e8;
            // JUST IGNORE FOR NOW BOUNDING LAYERS int  aMaxZmax=-1e8;
            // JUST IGNORE FOR NOW BOUNDING LAYERS cPt2di aPix;
            // JUST IGNORE FOR NOW BOUNDING LAYERS for (aPix.y()=0 ; aPix.y()<aSzL.y() ; aPix.y()++)
            // JUST IGNORE FOR NOW BOUNDING LAYERS {
            // JUST IGNORE FOR NOW BOUNDING LAYERS     for (aPix.x()=0 ; aPix.x()<aSzL.x() ; aPix.x()++)
            // JUST IGNORE FOR NOW BOUNDING LAYERS     {
            // JUST IGNORE FOR NOW BOUNDING LAYERS         if (aDZMin.GetV(aPix)<aMinZmin) aMinZmin=aDZMin.GetV(aPix);
            // JUST IGNORE FOR NOW BOUNDING LAYERS         if (aDZMax.GetV(aPix)>aMaxZmax) aMaxZmax=aDZMax.GetV(aPix);
            // JUST IGNORE FOR NOW BOUNDING LAYERS     }
            // JUST IGNORE FOR NOW BOUNDING LAYERS }
            // Construct a cube of features to be forwarded into the network
            torch::Device TheAvailDevice(mUseCuda ? torch::kCUDA : torch::kCPU);
            //MMVII_INTERNAL_ASSERT_always(aMaxZmax-aMinZmin, "PAX INTERVAL NULL !");

            // compute all possible correlation by dot product between feature vectors
            LREmbeddingsL=torch::nn::functional::normalize(LREmbeddingsL,F::NormalizeFuncOptions().p(2).dim(0));

            if( StepZ()<1.0)
            {
                std::vector<tREAL8> aScale={1.0,1.0/StepZ()};
                LREmbeddingsR= torch::nn::functional::interpolate(LREmbeddingsR.unsqueeze(0),
                                        F::InterpolateFuncOptions().mode(torch::kBilinear).
                                        align_corners(true).
                                        scale_factor(aScale)).
                                contiguous().squeeze();
            }

            cAutoTimerSegm aTEmbeddingsNormalizeInterpol(TimeSegm(),
                                            "DM4FillCubeCost::Embeddings Normalization/Interpol"); 


            StdOut()<<LREmbeddingsR.sizes()<<std::endl;

            LREmbeddingsR=torch::nn::functional::normalize(LREmbeddingsR,F::NormalizeFuncOptions().p(2).dim(0));

            auto aCorrelCube=torch::einsum("ijk,ijh->jkh",
                                                {LREmbeddingsL,LREmbeddingsR}).contiguous();

            
            cAutoTimerSegm aTSSimCompute(TimeSegm(),"DM4FillCubeCost::Scalar Product Via Einsum"); 


            // Extract the disparity map from the correlation cube
            ExtractSurface(aCorrelCube,mNamePax,StepZ());

            // argmax
            /*
            
            For now, no bounding layer filling for MMV1 as it is outdated
            just, argmax to extract a disparity 
            MAC: Implement Semi Global Matching algo for regularization
            */    
           /*
            StdOut()
                <<"Correl Cube size "<<
                aCorrelCube.sizes()<<std::endl;

            aCorrelCube=aCorrelCube.mul(-1.0).add(1.0).div(2.0);

            auto end = std::chrono::system_clock::now();
            std::chrono::duration<double>  elapsed_seconds = end-start;
            std::cout << "similarity computation elapsed time: " << elapsed_seconds.count() << "s"
                        << std::endl;

            // aCorrelCube: H,W1,W2 : H: height, W1: Width image1, W2: Width image2

            std::cout<<"GENERATED CUBE OF FEATURES "<<std::endl;
            // Refill tables defining layers of NAPPES WITH GENERATED CORREL
            start = std::chrono::system_clock::now();
            using namespace torch::indexing;
            int aSizeR=aCorrelCube.size(-1);
            for (aPix.y()=0 ; aPix.y()<aSzL.y() ; aPix.y()++)
            {
                for (aPix.x()=0 ; aPix.x()<aSzL.x() ; aPix.x()++)
                {
                    cPt2di aPAbs = aPix + mP0Z;
                    cPt2di aPC1  = aPAbs-mBoxGlob1.P0();
                    cPt2di aPC20 = aPAbs-mBoxGlob2.P0();
                    int aBInf=std::max(aVMods[0]->mCNNWin.x(),round_down(aPC20.x()/StepZ())+aDZMin.GetV(aPix));
                    int aBSup=std::min(round_up(aPC20.x()/StepZ())+aDZMax.GetV(aPix),aSizeR-aVMods[0]->mCNNWin.x());
                    auto aSim=aCorrelCube.index({Slice(aPix.y(),aPix.y()+1,1),
                                                    Slice(aPix.x(),aPix.x()+1,1),
                                                    Slice(aBInf,
                                                            aBSup+1)}).squeeze().contiguous();
                    std::vector<float> aVSim(aSim.data_ptr<float>(),
                                                aSim.data_ptr<float>() + aSim.numel());
                    int id_sim=0;
                    double aCost;
                    for (int aDz=aDZMin.GetV(aPix) ; aDz<aDZMax.GetV(aPix) ; aDz++)
                    {
                        aCost=1.0;
                        //cPt2di aPC2Z(round_ni(aPC20.x()+aDz*this->StepZ()),aPC20.y());  // INTEG FOR NOW
                        cPt2dr aPC2Z(aPC20.x()+aDz*this->StepZ(),aPC20.y());  // INTEG FOR NOW
                        bool IsInside=WindInside4BL(this->DI1(),aPC1,aVMods[0]->mCNNWin)
                                        && WindInside4BL(this->DI2(),aPC2Z,aVMods[0]->mCNNWin);
                        if (IsInside)
                        {
                            aCost=(double)aVSim[id_sim];
                            id_sim++;
                        }
                        PushCost(aCost);
                    }
                }
            }
            */
        }
        
    }
    else
    {
            /*cPt2di aSz = aDZMin.Sz();
            cPt2di aPix;

            int aCpt=0;

            for (aPix.y()=0 ; aPix.y()<aSz.y() ; aPix.y()++)
            {
                StdOut() << "Line " << aPix.y() << " on " << aSz.y()  << "\n";
                for (aPix.x()=0 ; aPix.x()<aSz.x() ; aPix.x()++)
                {
                        cPt2di aPAbs = aPix + mP0Z;
                        cPt2di aPC1  = aPAbs-mBoxGlob1.P0();
                        cPt2di aPC20 = aPAbs-mBoxGlob2.P0();
                        for (int aDz=aDZMin.GetV(aPix) ; aDz<aDZMax.GetV(aPix) ; aDz++)
                        {
                        double aTabCost[2];
                        bool   aTabOk[2];
                    for (int aK=0 ; aK<int(aVMods.size()) ; aK++)
                                aTabCost[aK] = aVMods[aK]->ComputeCost(aTabOk[aK],aPC1,aPC20,aDz);
                        aCpt++;
                        //PushCost(aTabCost[0]);

                    if (mCmpCorLearn && aTabOk[0] && aTabOk[1])
                    {
                            double aC0 = ToCmpCost(aTabCost[0]);
                            double aC1 = ToCmpCost(aTabCost[1]);
                            mImCmp.DIm().AddVBL(cPt2dr(aC1,aC0),1.0);
                    }
                        }
                }
            }*/
   }


   if (mCmpCorLearn)
   {
       mImCmp.DIm().ToFile("CmpCorrLearn_"+ mNamePost + ".tif");
   }

   //delete mFileCube;
   delete mPyrL1;
   delete mPyrL2;
   DeleteAllAndClear(aVMods);

   return EXIT_SUCCESS;
}



};

/* =============================================== */
/*                                                 */
/*                       ::                        */
/*                                                 */
/* =============================================== */
using namespace  cNS_FillCubeCost;

tMMVII_UnikPApli Alloc_FillCubeCost(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli & aSpec)
{
   return tMMVII_UnikPApli(new cAppliFillCubeCost(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpecFillCubeCost
(
     "DM4FillCubeCost",
      Alloc_FillCubeCost,
      "Fill a cube with matching costs",
      {eApF::Match},
      {eApDT::Image},
      {eApDT::ToDef},
      __FILE__
);

};
#endif
