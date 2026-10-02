//#include "MMVII_PCSens.h"

// Test commit

#include "MMVII_Image2D.h"
#include "MMVII_ImageMorphoMath.h"
#include "MMVII_Sensor.h"
#include "MMVII_Ptxd.h"
#include "MMVII_2Include_Serial_Tpl.h"
#include "MMVII_Tpl_ElemStrToVal.h"
#include "MMVII_Tpl_Images.h"
#include "MMVII_Tpl_ElemFilterLocImages.h"

#include "MMVII_Tpl_GraphStruct.h"
#include "MMVII_Tpl_Graph_SubGraph.h"
#include "MMVII_Tpl_GraphAlgo_SPCC.h"




namespace MMVII
{

namespace NS_FrangesDetect
{

/**  An application for  extaction curves on interference images
 *   rather very specific...
 */
typedef tREAL4            tElIm;
typedef cIm2D<tElIm >     tIm;
typedef cDataIm2D<tElIm > tDIm;
typedef cIm1D<tREAL8>     tIm1D;
typedef cDataIm1D<tREAL8>     tDIm1D;


class cConnComp
{
   public :

     cConnComp(std::vector<cPt2di> & aVPts,const tSeg2dr & aSeg,bool isHor,tREAL8 aScore,bool isOk) :
         mPts   (aVPts),
         mSeg   (aSeg),
         mIsHor (isHor),
         mScore (aScore),
         mIsOk   (isOk)
     {
     }

    std::vector<cPt2di> mPts;
    tSeg2dr mSeg;
    bool    mIsHor;
    tREAL8  mScore;
    bool    mIsOk;

};

/* =============================================== */
/*                                                 */
/*                 cAppliNewFrange                   */
/*                                                 */
/* =============================================== */


class cAppliNewFrange : public cMMVII_Appli
{
     public :
        cAppliNewFrange(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec);
     private :
        //================= typedef part ===============================


        typedef cIm2D<tU_INT1 >    tImLabel;
        typedef cDataIm2D<tU_INT1> tDImLabel;

        //================================================================
        //       METHODS DECLARATION
        //================================================================

            //------------------------- overidding cMMVII_Appli ----------------
        int Exe() override;
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override ;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override ;
         std::vector<cOneHelpSampleCmp>  Samples() const override;
        virtual ~cAppliNewFrange();

     private :

        void  DoOneImage(const std::string & aNameIm);
        std::string NameVisu(const std::string & aPref) const;

        /** Generete simulation images with "perfect" parabols, not realistic model, but
         *  to test correctnes of some computation */
        void MakeImSimul();

        /** Compute images of local maxims, this will reduce "drastically" the number of potential point
         in shortest path approach (and allow more complexe lins with "big" jumps ) */
        void MakeImageMaxLoc();


        bool NewCC(bool isHoriz,const cPt2di&);
        void ConnecCompMaxLoc(bool isHoriz);

        /** Make different tensor like computation :
             compute 2D tensor
             compute tan of Y=F(X) as tensor accumalation
             extract the center
         */
        void DoTensorProcessing();

        /// Compute the theoretical value of back ground & franges
        void ComputeRadiomCste();

        /// Estimate for a given Y, with window aSzY, how it is a good center of anti-symetry
        tREAL8 QualityCenterAntiSymetry(int aY,int aSzY);
        /// Compute the centre as "best" anti symetric local point
        int ComputeCenterByAntiSym();

        /// Integrate one way the tangent image and save it in integrale image
        void OneWayIntegrateTangent(int aDy,int aYLim);

        cPhotogrammetricProject  mPhProj;

        // ----------- Mandatory Args -----------
        std::string              mPatImage; /// Pattern of all images
        std::string              mNameIm;

        tREAL8    mZoomRed;     ///< Zoom for initial reduction
        tREAL8    mDerFactZ1;   ///< Factor of deriche gradient for zoom 1
        tREAL8    mSigmaTensZ1; ///< Sigma filter on tensor cumulated for zoom 1
        bool      mDoSimul;
        int       mDoVisu;

        tIm      mImZ1;
        tDIm*    mDImZ1;

        bool IsMax(cPt2di aP0,cPt2di aDp,int aNb) const;


        //-------- These values are related to reduced image ----------------
        tIm      mImRed;   /// reduced images
        tDIm*    mDImRed;  /// data reduced im
        cPt2di   mSzRed;     /// sz of reduced ima
        tIm      mImRedBlur;  /// im blured
        tDIm*    mDImRedBlur;  /// data image blurred


        tIm      mImTx;   /// image tensor x
        tDIm*    mDImTx;  /// data
        tIm      mImTy;   /// image tensor y
        tDIm*    mDImTy;   /// data

         cIm2D<tU_INT1> mImMaxHor;
         cDataIm2D<tU_INT1>* mDImMaxHor;
         cIm2D<tU_INT1> mImMaxVert;
         cDataIm2D<tU_INT1>* mDImMaxVert;

         cIm2D<tU_INT1> mImMaxLoc;
         cDataIm2D<tU_INT1>* mDImMaxLoc;


         tIm1D    mImTeta;    /// Image of tangent as x=F(Y)
         tDIm1D*  mDImTeta;
        tIm1D    mImTgt;    /// Image of tangent as x=F(Y)
        tDIm1D*  mDImTgt;
        int      mYC;
        tIm1D    mImIntegr;
        tDIm1D*  mDImIntegr;
        cRGBImage mImVisu;


        std::list<cConnComp>  mListCC;

        tREAL8   mRadiomBackGround;  ///< Estimation of radiometry of background
        tIm1D    mImRadFrange;       ///< Estimation of radiometry of franges, depend of Y
        tDIm1D*  mDImRadFrange;      ///< Data of mImRadFrange
};

}; // namespace NS_FrangesDetect
}; //  namespace MMVII

