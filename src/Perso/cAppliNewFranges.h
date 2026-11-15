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
typedef cIm1D<tREAL8>      tIm1D;
typedef cDataIm1D<tREAL8>  tDIm1D;
typedef cImGrad<tElIm>    tImGrad;


class cConnComp
{
   public :

     cConnComp(std::vector<cPt2di> & aVPts,const tSeg2dr & aSeg,cPt2dr aTg1,cPt2dr aTg2,bool isHor,tREAL8 aScore,bool isOk);

        std::vector<cPt2di> mPts;
        tSeg2dr mSeg;
        cPt2dr  mTg1;
        cPt2dr  mTg2;
        bool    mIsHor;
        tREAL8  mScore;
        bool    mIsOk;

};


/* *********************************************************** */
/*                                                             */
/*          Declaration of Vector-Graph                        */
/*                                                             */
/* *********************************************************** */

enum class eTypePix
{
    eTmp=0,
    eMaxVert=1,
    eMinVert=2,
    eMaxTens=3,
    eMinTens=4,
    eProlong=5,
    eTypeSmalCC=6,
    eTypeBadAvgCC=7
};



class cMaxLocSum
{
   public :
      tREAL8 mXC;
      tREAL8 mX0;
      tREAL8 mX1;
};

class cAVertReg : public cMemCheck
{
   public :
       cAVertReg(const  std::vector<cPt2di> & aVPts,int aNum);



       bool IsBorderExt() const;

       std::vector<cPt2di> mVPts;

       int mNumV;       // Num of vertice ~ num in graphe
       int mNumCC0;     // Num of connect comp it belongs to
       int mNumMaxLoc;  // Possible num of local maxima
       tREAL8    mXC;   // Absice where it eventually intersecte central line


       /// Extrem, such that P1->P2 goes Up for horizontal line
       cPt2di   mP1;
       cPt2di   mP2;
       cPt2di   mCenter;
       cBox2di  mBox;  // Box of Vpts

       bool  mUpSide;
       bool  mBottomSide;
       bool  mNoCross;
       bool  mTouchLeft;
       bool  mTouchRight;
};

class cAOrientReg
{
    public :
       cAOrientReg();
       cPt2di mPt;
};

class cASymReg
{
     public :

        cASymReg();
        tREAL8 mMaxVal;
        bool mInMST;
};

typedef cVG_Vertex<cAVertReg,cAOrientReg,cASymReg>    tVertexReg;
typedef cVG_Edge<cAVertReg,cAOrientReg,cASymReg>      tEdgeReg;
typedef cVG_Graph<cAVertReg,cAOrientReg,cASymReg>     tGraphReg;


/* =============================================== */
/*                                                 */
/*                 cAppliNewFrange                 */
/*                                                 */
/* =============================================== */


struct cSigmaTensors
{
    double mSigmaTeta = 15.0;
    double mSigmaTensX=20.0;
    double mSigmaTensY= 5.0;
    double mSigmaGradX=2.0;
    double mSigmaGradY=2.0;


    ARG2007_STRUCT_FIELDS (
        mSigmaTeta,FieldSem({{eTA2007::AddCom,"Sigma Teta"},eTA2007::HDV}),
        mSigmaTensX,FieldSem({{eTA2007::AddCom,"Sigma Tensor-X"},eTA2007::HDV}),
        mSigmaTensY,FieldSem({{eTA2007::AddCom,"Sigma Tensor-Y"},eTA2007::HDV}),
        mSigmaGradX,FieldSem({{eTA2007::AddCom,"Sigma Grad-X"},eTA2007::HDV}),
        mSigmaGradY,FieldSem({{eTA2007::AddCom,"Sigma Grad-Y"},eTA2007::HDV})
    )
};

struct cParamVisu
{
   public :
    tREAL8 mMulGauss = 0.7;
    tREAL8 mMulLAB   = 0.3;
    tREAL8 mMulSat   = 0.5;

    ARG2007_STRUCT_FIELDS (

            mMulGauss,FieldSem({{eTA2007::AddCom,"Multiplier in gaussian equalisation"},eTA2007::HDV}),
            mMulLAB  ,FieldSem({{eTA2007::AddCom,"Multiplier of AB in LAB-space-colors"},eTA2007::HDV}),
            mMulSat  ,FieldSem({{eTA2007::AddCom,"Multiplier of Saturation in HSI space"},eTA2007::HDV})

    )
};


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

        //-------------------------- FONCTION RELATED TO VISUALIZATION -----------------
        /// Generate several visualization
        void DoVisu();
        /// Generate back ground gray image for diff visu, used Im-Reduc-Blured  & Gaussian equam
        cRGBImage ImBckGround(tREAL8) const ;
        /// Put in anIm a colorizatuion depenaning of labels, aCoefAB multiply AB in LAB-space
        void SetImageLabel(cRGBImage & anIm,tREAL8 aCoefAB) const;

        /// Show extrmitys of line
        void ShowExtres(cRGBImage & anIm) const;


        void ShowCurvesIntegral(cRGBImage&);


        /// Name of visualisation image
        std::string NameVisu(const std::string & aPref) const;

        template <class TypeIm> void GenVisu(const TypeIm & anIm,const std::string & aPref)
        {
            if (MatchRegex(aPref,mPatVisu))
            {
                anIm.ToFile(NameVisu(aPref));
                mNbVisuGen++;
            }
            mVNameVisu.push_back(aPref);
        }

        /// Classical  saddle point deprecated 4 now
        bool IsSaddlePoint(const cPt2di & aPix) const;

        ///  compute next extremity of a line, by following (+-) crest-line close to previous direction
        cPt2dr NextExtreLine(const cPt2dr & aP0,const cPt2dr& aP1) const;

        /// compute direction of tensor by reading 2D-images of tensor
        cPt2dr TensorImage(cPt2di aPix) const;

        /// compute direction of tensor by teta  of y
        cPt2dr TensorTeta(cPt2di aPix) const;


        void  ProlongateExtre(std::vector<cPt2dr>& aRes,cPt2di  aPix,eTypePix aLabel) const;

        /*==========================  LABEL MANIPULATION =====================*/

        static int Lab2I(eTypePix aLabel) {return (int) aLabel;}
        static int LabTmp2I() {return   Lab2I(eTypePix::eTmp);}


              // ------------------ Fix On/Off -----------------------------
        void SetLabelOn(const cPt2di & aPix,eTypePix aLabel)
        {
            SetBitOn(*mDImFlagLabel,aPix,Lab2I(aLabel));
        }
        void SetLabel(const cPt2di & aPix,eTypePix aLabel)
        {
            SetFlagOfBit(*mDImFlagLabel,aPix,Lab2I(aLabel));
        }

        void SetLabelOn(const cPt2di & aPix,const std::vector<eTypePix> & aLabel);

        void SetLabelOff(const cPt2di & aPix,eTypePix aLabel)
        {
            SetBitOff(*mDImFlagLabel,aPix,Lab2I(aLabel));
        }
             // ----------------- extract value ------------------------------
        bool LabelIsOn(const cPt2di & aPix,eTypePix aLabel) const
        {
            return BitIsOn(*mDImFlagLabel,aPix,(int)aLabel);
        }

        const cPt3di * CoulOfLabel(const cPt2di & aPt) const;

                // ----------------- Specific to Tmp ------------------------------
        bool LabelTmpIsOff(const cPt2di & aPix) const
        {
            return ! LabelIsOn(aPix,eTypePix::eTmp);
        }
        void SetLabelTmpOff(const cPt2di & aPix)
        {
            SetLabelOff(aPix,eTypePix::eTmp);
        }



        /** Generete simulation images with "perfect" parabols, not realistic model, but
         *  to test correctnes of some computation */
        void MakeImSimul();

        /// Compute the theoretical value of back ground & franges
        void ComputeRadiomCste();



        ///
        void MakeImageRegion();

        void MakeImageMarqIn();


        void ComputeArribute(cAVertReg&,cDataIm2D<tU_INT1> &);

        void PostProcEdge(tEdgeReg&);

        void PostProcVertex(tVertexReg&);




        // =================== MaxLoc & Pix -CC  ==================

        /** Compute images of local maxims, this will reduce "drastically" the number of potential point
         in shortest path approach (and allow more complexe lins with "big" jumps ) */
        void AnalyseCC();

        /// is P0 it a maxmima in direction Dp on aNbPixels
      //  bool IsMinMax(cPt2di aP0,cPt2di aDp,int aNb,bool IsMax) const;

        /** Compute the connect component of Horiz/Vertical chain
         !! Horiz  is the direction for compuing max, so the chain with be +- orthog */
        void ConnecCompMaxLoc(eTypePix );
        /// Compute the connected component and eventually memorize it
        eTypePix NewCC(eTypePix,std::vector<cPt2di>  &,const cPt2di&);


        /// Add heuristically some horizontal CC
        void ComplemCCHoriz();

        /// Must a point be added to assured 4-Connectivity
        bool AddFor8Cx(const cDataIm2D<tU_INT1> & aImMax,const cPt2di&) const;
        /// Add all the point that assure 8-Connecitvity
        void Make8Cx(cDataIm2D<tU_INT1> & aImMax);





        /// Follow the chain, on a fix number of points stop if dead end or branch
        cPt2di BackTrackExtre( cPt2di aP0,const cDataIm2D<tU_INT1> & aImMax) const;
        cPt2di OneStepBackTrackExtre(const cPt2di & aP0,const cPt2di & aPLast,const cDataIm2D<tU_INT1> & aImMax) const;

        cPt2di NewBackTrackExtre( cPt2di aP0,const std::vector<cPt2di> & aVPts,tREAL8 aDist) const;

        // =================== "High-level" Graph (on max-loc)  ==================

         //void GraphCompute();
         //void GraphCreate();
        // void GraphCCHoriz();




        // ===================  Tensor & angle from tensor processing ==================

        /** Make different tensor like computation :
             compute 2D tensor
             compute tan of Y=F(X) as tensor accumalation
             extract the center
         */
        void DoTensorProcessing(bool isFirst,bool isLast);


        void ExtractTensorLines(bool isMax);
        void ExtractTensorLines();


        /// Is it a maxima in the tensor direction ?
        bool IsExtreLocTens(const cPt2di & aPt,const cPt2dr & aDir,tREAL8 aSignVal) const;


        /// Estimate for a given Y, with window aSzY, how it is a good center of anti-symetry
        tREAL8 QualityCenterAntiSymetry(bool First,int aY,int aSzY,const std::vector<cPt2dr>&,const std::vector<tREAL8>&);
        /// Compute the centre as "best" anti symetric local point
        int ComputeCenterByAntiSym(bool First,const std::vector<cPt2dr>&,const std::vector<tREAL8>&);

        /// Integrate one way the tangent image and save it in integrale image
        void OneWayIntegrateTangent(int aDy,int aYLim);

        /// Integrate one way the tangent image and save it in integrale image
        int ComputeZoneYCenter(int aDy,int aYLim);

        ///< compute max +- basic with radiom of center
        void ComputeMaxRadCenter();

        //---------- Version using integral cuvre, OK but get false detection with broken curve;
        ///< compute projection on x using integral
        void ComputeSumOnIntegralCurve();
        ///< Extract maxima on sum of integral curve
        void ExtractMaxSumIC();
        ///<
        bool IsMaxLocSumIC(int aY) const;



        void MakeGrapheTensor(int aX0,int aDx,int aDy, cRGBImage * = nullptr);
        void MakeGrapheTensor(const cIm2D<tU_INT1>);

        //=========================================  DATA ========================================

        cPhotogrammetricProject  mPhProj;

        // ----------- Mandatory Args -----------
        std::string              mPatImage; /// Pattern of all images
        std::string              mNameIm;


        ///  Sigmas for first iter tensor-computation
        cSigmaTensors            mSigTensors;
        /// Param for visualisation
        cParamVisu               mParamVisu;


        tREAL8    mZoomRed;     ///< Zoom for initial reduction
        tREAL8    mDerFactZ1;   ///< Factor of deriche gradient for zoom 1
        tREAL8    mSigmaTensZ1; ///< Sigma filter on tensor cumulated for zoom 1
        int       mDoSimul;
        std::string mPatVisu;        ///< Pattern for selecting visualizations
        std::vector<std::string> mVNameVisu;  ///< Memorize all visu
        int         mNbVisuGen; ///< Number of visu generated

         std::vector<cMaxLocSum> mVMaxLoc;
         tREAL8   mIntervMax;     ///< Inteval betwen consec max +- period of the bouzin
         int      mMinZCYC;       /// Begin of central zone arround YC
         int      mMaxZYC;        ///< End of central zone arround YC



        bool mHasMask;  ///< Is there a mask
        cIm2D<tU_INT1> mImMask; ///< Possible mask
        tIm      mImZ1;
        tDIm*    mDImZ1;



        //-------- These values are related to reduced image ----------------
        tIm      mImRed;       ///< reduced images
        tDIm*    mDImRed;      ///< data reduced im
        cPt2di   mSzRed;       ///< sz of reduced ima
        tIm      mImRedBlur;   ///< im blured
        tDIm*    mDImRedBlur;  ///< data image blurred


        tImGrad mImGrad;
        tImGrad mImTens;


        /// Image of flag of labels
        cIm2D<tU_INT4>        mImFlagLabel;
        cDataIm2D<tU_INT4>*   mDImFlagLabel;
        cIm2D<tINT2>          mImRegion;
        cDataIm2D<tINT2>*     mDImRegion;

        cIm2D<tU_INT1>          mImMarqIn;
        cDataIm2D<tU_INT1>*     mDImMarqIn;


         tIm1D    mImCumNormT;    /// Image of teta as x=F(Y)
         tDIm1D*  mDImCumNormT;
         tIm1D    mImTeta;    /// Image of teta as x=F(Y)
         tDIm1D*  mDImTeta;
        tIm1D    mImTgt;      /// Image of tangent as x=F(Y)
        tDIm1D*  mDImTgt;
        int      mYC;          ///<  estimation of Y center
        tIm1D    mImIntegr;        ///<  integration of tangent
        tDIm1D*  mDImIntegr;       ///<  data
        tIm      mImSmoothTgt;  ///< Image smoothed takin account tgt
        tDIm *   mDImSmoothTgt;  ///< Data
        tIm1D    mImProjInt;     ///< cumulated projection
        tDIm1D*  mDImProjInt;    ///< Data
        tREAL8   mStdDevProjI;   ///< standard deviation
        tREAL8   mMinDPI;
        tREAL8   mMaxDPI;


        tREAL8   mRadiomBackGround;  ///< Estimation of radiometry of background
        tIm1D    mImRadFrange;       ///< Estimation of radiometry of franges, depend of Y
        tDIm1D*  mDImRadFrange;      ///< Data of mImRadFrange


        tGraphReg     mGrReg;
        tVertexReg*   mVirtExtVert; ///< Virtual Exten vertex : the border

        /// Image of backround used in visu if any
        cIm2D<tU_INT1>  mImGrayBG;
        cRGBImage       mImRgbBG;

};

}; // namespace NS_FrangesDetect
}; //  namespace MMVII

