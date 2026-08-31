//#include "MMVII_PCSens.h"

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


struct cASom
{
  public :
      cASom(const cPt2di & aPt ) : mPt (aPt) {}
      cPt2di mPt;
};

struct cAEdgeSym
{
   public :
};

struct cAEdgeOriented
{
    public :
       cAEdgeOriented(const tREAL8 aCost) : mCost (aCost) {}
       tREAL8 mCost;
};

typedef cVG_Vertex<cASom,cAEdgeOriented,cAEdgeSym> tVertexFrange;
typedef cVG_Edge<cASom,cAEdgeOriented,cAEdgeSym>  tEdgeFrange;
typedef cVG_Graph<cASom,cAEdgeOriented,cAEdgeSym> tGraphFrange;

class cWeithGr : public cAlgo_ParamVG<tGraphFrange>
{
   public :
      tREAL8 WeightEdge(const    tEdgeFrange & anEdge) const  override
      {
             return anEdge.AttrOriented().mCost;
      }
};

//template <class TGraph>  class cAlgo_ParamVG


struct cCalcSym
{
    double mFactRed;
    double mWitdhCalc;
    double mDistFrange;


    cCalcSym () :
      mFactRed    (5.0),
      mWitdhCalc  (40.0),
      mDistFrange (500.0)
    {
    }

    ARG2007_STRUCT_FIELDS (
        mFactRed,FieldSem({eTA2007::AddCom,"Reduction factor"}),
        mWitdhCalc,FieldSem({eTA2007::AddCom,"Witdh for calc of symetry(at full resol)"}),
        mDistFrange,FieldSem({eTA2007::AddCom,"Typicall dist between Frange (at full resol)"})

        )
};

class cFitParabol
{
   public :
    cFitParabol (tREAL8 anYC,tREAL8 anInY,tDIm *,tREAL8 aStep=1.0);
    void SetX(tREAL8 anX,tREAL8 anXAtLim);

    tREAL8 Score(tREAL8 anX0,tREAL8 aXAtLim,tREAL8 aY0Rel,tREAL8 aY1Rel) ;

    cPt2dr PtOfDY(tREAL8 anY) const;

    //cPt2dr
   private :
       tREAL8   mYCenter;
       tREAL8   mIntervY;
       int      mNbY;
       tREAL8   mStepY;
       tDIm*    mDIm;
       tREAL8   mX0;
       tREAL8   mXAtLim;
};

cFitParabol::cFitParabol (tREAL8 anYC,tREAL8 anIntervY,tDIm * aDIm,tREAL8 aStepY) :
    mYCenter  (anYC),
    mIntervY  (anIntervY),
    mNbY      (round_ni(mIntervY/aStepY)),
    mStepY    (mIntervY/mNbY),
    mDIm      (aDIm)
{
}

void cFitParabol::SetX(tREAL8 anX0,tREAL8 anXAtLim)
{
    mX0 = anX0;
    mXAtLim = anXAtLim;
}

cPt2dr  cFitParabol::PtOfDY(tREAL8 aDY) const
{
    tREAL8 anX = mX0 + Square(aDY/tREAL8(mIntervY)) * mXAtLim;
    return cPt2dr(anX,mYCenter+aDY);
}



tREAL8 cFitParabol::Score(tREAL8 anX0,tREAL8 aXAtLim,tREAL8 aY0Rel,tREAL8 aY1Rel)
{
    SetX(anX0,aXAtLim);
    tREAL8 aSom = 0;
    for (int aDY= round_ni(mNbY*aY0Rel) ; aDY<= round_ni(mNbY*aY1Rel) ; aDY++)
    {
        cPt2dr aPt = PtOfDY(aDY*mStepY);
        if (! mDIm->InsideBL(aPt))
            return -1e10;

        aSom += mDIm->GetVBL(aPt);
    }

    return aSom;
}

/* =============================================== */
/*                                                 */
/*                 cCurvFrange                     */
/*                                                 */
/* =============================================== */

/**
 * @brief Class for creating/storing a curve in the Frange appli
 */
class cCurvFrange
{
    public :

       /// Constructor, box use to size the curve
       cCurvFrange(const cBox2di &);

       /// Create points from the average,  z is the weighting
       std::vector<cPt3dr>  ExtractCurve(const cPt2dr & anOffset) const;

       /// Add a pix for averaging
       void AddPix(const cPt2di & aPix,tREAL8 aW);

       ///  Comparator used for sorting from left to right
       bool operator < (const cCurvFrange& aC2) const;
    private :
       cBox2di              mBox;
       std::vector<tWArr>   mAvX;
};



cCurvFrange::cCurvFrange(const cBox2di & aBox) :
    mBox (aBox),
    mAvX (aBox.P1().y()+1)  // +1 => P1 is include
{
}

//  Compare on P0.x for sorting left to right
bool cCurvFrange::operator < (const cCurvFrange& aC2) const
{
   return  mBox.P0().x() < aC2.mBox.P0().x();
}

// Accumulate a pixel in averaging
void cCurvFrange::AddPix(const cPt2di & aPix,tREAL8 aW)
{
   mAvX.at(aPix.y()).Add(aW,aPix.x());
}


// Extract average point with no null weighting
std::vector<cPt3dr>  cCurvFrange::ExtractCurve(const cPt2dr & anOffset) const
{
   std::vector<cPt3dr> aRes;
   for (size_t anY =0 ; anY<mAvX.size() ; anY++)
   {
       const auto & anAv = mAvX.at(anY);
       if (anAv.SW() != 0)
       {
           aRes.push_back(cPt3dr(anAv.Average()+anOffset.x(),anY+anOffset.y(),anAv.SW()));
       }
   }
   return aRes;
}

/* =============================================== */
/*                                                 */
/*                 cExportFrange                   */
/*                                                 */
/* =============================================== */

class cExportFrange
{
   public :
       int                 mLabel;
       std::vector<cPt3dr> mVPts;
};

void AddData(const cAuxAr2007 & anAux, cExportFrange & anEF)
{
    AddData(cAuxAr2007("Label", anAux),anEF.mLabel);
    StdContAddData(cAuxAr2007("Pts", anAux),anEF.mVPts);
}

/* =============================================== */
/*                                                 */
/*                 cAppliFranges                   */
/*                                                 */
/* =============================================== */


class cAppliFranges : public cMMVII_Appli
{
     public :
        cAppliFranges(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec,int aMode);
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
        virtual ~cAppliFranges();

            //------------------------- Specific method----------------
        ///  Do all the job for one image
        void  DoOneImage(const std::string & aNameIm) ;

        /// Do the job for one connected component
        void  AnalyseOneConnectedComp(const cPt2di& aPix);
        /// Validate a connected component, decision made on bounding box for now
        bool  ConnectedCompIsValide(const std::vector<cPt2di> & aVPts,const cBox2di & aBox);
        ///  Create a new curve once the CC is valid
        void  CreateNewCurve(const std::vector<cPt2di> & aVPts,const cBox2di&);
        /// convert a gray level in a weigthing inside [0,1]
        tREAL8 WeightOfGray(tREAL8 aGray) const;

        /// Generate visualisation if required
        void MakeVisu();

        void ReadImage(const cBox2di& aBox);

        void CalcAxeSym();

        void ComputeCurvePrgDyn(int aX0);
        void HalfComputeCurvePrgDyn(int aX0,int  aSignY);
        tREAL8 Cost(const tVertexFrange * aV0,const tVertexFrange * aV1,int aSignY) const;
        tREAL8 Cost(tREAL8 aValIm) const;

        void   ComputeParabol(tREAL8 aX0) ;
        void   ComputeParabol(tREAL8 aX0,tREAL8 aY0Rel,tREAL8 aY1Rel) ;

        tREAL8 ScoreParabol(tREAL8 aX0,int aIntY,tREAL8 aXAtLim) const;


        //================================================================
        //       DATA DEFINITION
        //================================================================

        /// == Used for file naming, even if obviously not a photogrammetric context ========
        cPhotogrammetricProject  mPhProj;
        int                      mMode;

        // ----------- Mandatory Args -----------
        std::string              mPatImage; /// Pattern of all images

        // -----------  Optionnal args  -----------

        cPt2di                   mIntY;    ///<  Y Interval of image
        cPt2dr                   mSigma;   ///<  Sigma of smoothing
        tREAL8                   mSigCurv;  ///< Sigma for comutig Cuvres
        int                      mNbIter;  ///<  Number of iteration of smoothing

        tREAL8                   mWidhMinAll;    ///< Witdh Min of Connected compenent
        tREAL8                   mWidhMinBorder; ///< Idem but for border image
        tREAL8                   mMinHeightBorderRight;
        bool                     mDoVisu;        ///< Do we generate visualisation
        std::vector<std::string> mParamStack;    ///< Param for image stacking if any

        cCalcSym                 mCalcSym;  ///< Parameter for symetry calc
        int                      mWidthCS;  ///< With of calc sym


        //  -------------  Variable for Prg-Dyn approach -----------------------

        int                     mYCenter; ///< Center of sym
        tREAL8                  mFactRed;
        tREAL8                  mLowVal;
        tREAL8                  mVeryLowVal;
        tREAL8                  mHighVal;
        tIm                     mImRed;
        tDIm*                   mDImRed;
        tIm                     mImBlur;
        tDIm*                   mDImBlur;
        cPt2di                  mSzImRed;
        cRGBImage               mImVisuPrgDyn;
        tGraphFrange *          mGraph;

        //  ----------- Variables for images -----------
        std::string              mNameCurIm;   ///< Name of current image
        cPt2di                   mCurSz;       ///< Size of current image
        tIm                      mCurIm;       ///< Loaded current image
        tDIm*                    mDCurIm;      ///< Data current im
        tIm                      mImSmooth;    ///< Smoothed image
        tDIm*                    mDImSmooth;   ///< Data smoothed
        tImLabel                 mImLabel;     ///< Label image
        tDImLabel*               mDImLabel;    ///< Data label

        //  ----------- miscellaneous -----------

        cTimerSegm               mTimeSeg;   ///< Used for timing info
        int                      mMarqCurCC; ///< Counter for connected component
        tREAL8                   mThreshW;   ///< Trehshold for white
        std::vector<cCurvFrange> mCurves;    ///< Result of curves extraction
        std::string              mIdExport;  ///< Identiant for Export (CSV/Tiff ...)
};

    /* =================================================== */
    /*      Overiding of  cMMVII_Appli                     */
    /* =================================================== */


cAppliFranges::cAppliFranges(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec,int aMode) :
    cMMVII_Appli      (aVArgs,aSpec),
    mPhProj           (*this),
    mMode             (aMode),

    // -------- Defautl values of opt param ----
    mIntY             (0,10000000),
    mSigma            (30.0,5.0),
    mNbIter           (5),
    mWidhMinAll       (200.0),
    mWidhMinBorder    (300.0),
    mMinHeightBorderRight (200.0),
    mDoVisu           (mMode==1),
    mParamStack       {},
    // -----------------
    mYCenter          (-1),
    mFactRed          (0.0),
    mImRed            (cPt2di(1,1)),
    mDImRed           (nullptr),
    mImBlur           (cPt2di(1,1)),
    mDImBlur          (nullptr),
    mSzImRed          (-1,-1),
    mImVisuPrgDyn     (cPt2di(1,1)),
    mGraph            (nullptr),
    //----------  Mandatory init of images --------
    mCurIm            (cPt2di(1,1)),
    mImSmooth         (cPt2di(1,1)),
    mImLabel          (cPt2di(1,1)),
    mTimeSeg          (this)

{
    mCurves.reserve(20);  // to save memory
}

cAppliFranges::~cAppliFranges()
{
}

cCollecSpecArg2007 & cAppliFranges::ArgObl(cCollecSpecArg2007 & anArgObl)
{
      return    anArgObl
            <<  Arg2007(mPatImage,"Name of input Image", {eTA2007::FileDirProj,{eTA2007::MPatFile,"0"}})
      ;
}

cCollecSpecArg2007 & cAppliFranges::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{

     if (mMode==0)
     {
         anArgOpt
            << AOpt2007(mSigma,"Sigma","Sima x/y of initial smoothig for CC",{eTA2007::HDV})
            << AOpt2007(mSigCurv,"SigCurv","Sima for smoothig curve",{eTA2007::HDV})
            << AOpt2007(mIntY,"IntY","Interval for Y",{eTA2007::HDV})
            << AOpt2007(mNbIter,"NbIt","Number of iter initial smoothing",{eTA2007::HDV})
            << AOpt2007(mDoVisu,"DoVisu","Generate Visualisation ?",{eTA2007::HDV})
            << AOpt2007(mWidhMinAll,"MinWidth","Minimal witdh, general case",{eTA2007::HDV})
            << AOpt2007(mWidhMinBorder,"BorderMinWidth","Minimal witdh for border",{eTA2007::HDV})
            << AOpt2007(mMinHeightBorderRight,"MinHeightBorderRight","Minimal Heitgh for right border",{eTA2007::HDV})
            << AOpt2007(mParamStack,"Stack","Stacking parm [Pat,Nb,Mode] ",{{eTA2007::ISizeV,"[3,3]"}})
         ;
     }
     else if (mMode==1)
     {
         anArgOpt
            << AOpt2007(mCalcSym,"CalcSym","Paramater for symetry calc",{eTA2007::HDV})
         ;
     }

     return anArgOpt;
}


std::vector<cOneHelpSampleCmp>  cAppliFranges::Samples() const
{
   return
   {
       {"MMVII ExtractFranges Retiga_000000105.tif Sigma=[10,2] DoVisu=1"},
       {"MMVII ExtractFranges Retiga_.*.tif"}
   };
}

int cAppliFranges::Exe()
{
    mPhProj.FinishInit();

    if (RunMultiSet(0,0)) // Case several images in //
    {
       return ResultMultiSet();
    }
    DoOneImage(UniqueStr(0)); // Case 1 image (may be recalled from multiple)
    return EXIT_SUCCESS;
}

/* =================================================== */
/*              Specific functions                     */
/* =================================================== */

// Linear ressample of [Thr,255]  in [0,1]
tREAL8 cAppliFranges::WeightOfGray(tREAL8 aGray) const
{
   return  std::min(1.0,std::max(0.0,(aGray-mThreshW) /(255.0-mThreshW)));
}




bool cAppliFranges::ConnectedCompIsValide(const std::vector<cPt2di> & aVPts,const cBox2di & aBox)
{
    // CC touches  extrem left, not full -> false
    if (aBox.P0().x() <= 10) return false;

    if (aBox.Sz().x() < mWidhMinAll)  // not width enough
        return false;

    // more strict width rule for CC touching the right
    if (aBox.P1().x() >=  mCurSz.x()-10)
    {
        if  (aBox.Sz().x() < mWidhMinBorder)
            return false;
    }

      // right border, component may have been cut, relax constraint on height
     if (aBox.P1().x() >=  mCurSz.x()-10)
     {
         if (aBox.Sz().y() < mMinHeightBorderRight) return false;
     }
     else  // Do not touch the right must be full height
     {
       // CC must fill the entire hight
       if (aBox.P0().y()!= 1) return false;
       if (aBox.P1().y()!= mCurSz.y()-2) return false;
     }

    // So far, so good ...
    return true;
}

void cAppliFranges::CreateNewCurve(const std::vector<cPt2di> & aVPts,const cBox2di& aBox)
{
    mCurves.push_back(cCurvFrange(aBox));  // add an empty curve
    cCurvFrange & aCurv = mCurves.back();  // reference it

    // Parse all the pixel to add them in curve
    for (const auto & aPix : aVPts)
    {
        aCurv.AddPix(aPix, WeightOfGray(mDImSmooth->GetV(aPix)));
    }
}


void  cAppliFranges::AnalyseOneConnectedComp(const cPt2di& aPix)
{
   std::vector<cPt2di> aVPts;
   // Call the method of lib to push CC in aVPts
   ConnectedComponent
   (
       aVPts ,
       *mDImLabel,
       Alloc8Neighbourhood(),
       aPix, // seed
       255,   // CC of point=255
       mMarqCurCC // Value set in point selected
   );

   // Create new box, possible for example all point have same x => Allow Empty
   cBox2di aBox = cBox2di::FromVect(aVPts,eAllowEmpty::Yes);

   if (ConnectedCompIsValide(aVPts,aBox)) // are Vpt/Box OK
   {
       CreateNewCurve(aVPts,aBox);  // is ok, create curve
       mMarqCurCC++;       // to have different label on next CC
   }
   else
   {
       mDImLabel->VPtsSetV(aVPts,1); // put 1 as label for rejected CC
   }
}

void cAppliFranges::ReadImage(const cBox2di& aBox )
{
   // if no stack, basic read
   if (!IsInit(&mParamStack))
   {
       mCurIm = tIm::FromFile(mNameCurIm,aBox);
       return;
   }

   // --- read parameters of mParamStack = [Pat,NbIm,Mode] -------
   std::string aPat = mParamStack.at(0);
   int aNbImRequired = cStrIO<int>::FromStr(mParamStack.at(1));
   int aMode = cStrIO<int>::FromStr(mParamStack.at(2));

   // --- read Names, sort and extract index of CurName --------------
   std::vector<std::string> aVName =  GetFilesFromDir(DirProject(),AllocRegex(aPat)); // read All Names
   std::sort(aVName.begin(),aVName.end()); // sort name
   const auto & anIter =  std::find(aVName.begin(),aVName.end(),mNameCurIm); // extract index of name
   if (anIter==aVName.end())
   {
       MMVII_WARNING(" Image not found in stacking, use standard read \n");
       mCurIm = tIm::FromFile(mNameCurIm,aBox);
       return;
   }
   int aKC = anIter - aVName.begin();

   // --- Compute correct interval of index of image, warantee to be inside aVName
   int aK0 = std::max(0,aKC-aNbImRequired);
   int aK1 = std::min(int(aVName.size())-1,aKC+aNbImRequired);
   int aNbImUsed = std::min(aKC-aK0,aK1-aKC);
   aK0 = aKC-aNbImUsed;
   aK1 = aKC+aNbImUsed;

  // ------- Read Images ------------------
   std::vector<tIm> aVIm;
   for (int aKIm=aK0 ; aKIm<=aK1; aKIm++)
   {
       aVIm.push_back(tIm::FromFile(aVName.at(aKIm),aBox));
   }
   mCurIm = tIm(mCurSz);


   // ------- Compute weigthing of index ------------------
   std::vector<tREAL8> aVW;
   for (size_t aKV=0 ; aKV<aVIm.size() ; aKV++)
   {
       tREAL8 aW=1.0;
       tREAL8 aRnk = (aKV+0.5) / aVIm.size();
       if (aMode==1)
           aW =  (aKV == (aVIm.size()/2)); // Median
       else if (aMode==2)
       {
           tREAL8 aTeta = aRnk * 2 * M_PI;
           aW = 1+std::cos(aTeta +M_PI);
       }
       else if (aMode==3)
       {
           aW = 0.5 - std::abs(aRnk-0.5);
       }
       aVW.push_back(aW);
   }

   //--- Compute weighted average --------------------------
   for (const auto & aPix : mCurIm.DIm())
   {
       std::vector<tREAL8> aVVals;
       for (const auto & anIm : aVIm )
           aVVals.push_back(anIm.DIm().GetV(aPix));
       std::sort(aVVals.begin(),aVVals.end());
       cWeightAv<tREAL8,tREAL8> aWAv;

       for (size_t aKV=0 ; aKV<aVVals.size() ; aKV++)
       {
           aWAv.Add(aVW.at(aKV),aVVals.at(aKV));
       }
       mCurIm.DIm().SetV(aPix,aWAv.Average());
   }
/*
   StdOut() << " STACKKK " << aVName.size()
            << " N0=" << aVName.at(0)
            << " KF=" << aKC
            << " NBI=" << aNbImUsed
            << " MODE=" << aMode
            << " W " << aVW
            << "\n";
   getchar();
*/
}


tREAL8 cAppliFranges::Cost(tREAL8 aValIm) const
{
    tREAL8 aCoeff = (aValIm-mVeryLowVal) / (mHighVal-mVeryLowVal);

    // + to be continous and should be sufficcient, btw :  max to be 100% sure >0,
    aCoeff = std::max(1e-2,aCoeff+1e-1);

    return std::pow(1.0 / aCoeff,3.0);
}

tREAL8 cAppliFranges::Cost(const tVertexFrange * aVertex0,const tVertexFrange * aVertex1,int aSignY) const
{
   cPt2di aPix0 = aVertex0->Attr().mPt;
   cPt2di aPix1 = aVertex1->Attr().mPt;
   cPt2di aV01 = aPix1-aPix0;

   tREAL8 aV0 = mDImRed->GetV(aPix0);
   tREAL8 aV1 = mDImRed->GetV(aPix1);

   tREAL8 aCost = Cost(aV0) + Cost(aV1);

   //aCost = std::pow(aCost,3.0);
   return aCost * (std::abs(aV01.x()) * 1.0 + std::abs(aV01.y()) * 1.0);

}




void cAppliFranges::HalfComputeCurvePrgDyn(int aX0,int  aSignY)
{
    mGraph = new tGraphFrange();
    int mWitdh = 500 / mFactRed;
    int mLenght = 1200 / mFactRed;

    int aX1 = std::min(aX0+mLenght,mSzImRed.x());
    int aY0 =  mYCenter;
    int aY1 = aY0 + mWitdh * aSignY;
    aY1 = std::max(0,std::min(aY1,mSzImRed.y()));

    std::vector<std::vector<tVertexFrange *>> aMatVertex;
    tVertexFrange * aVertexInit = nullptr;
    std::vector<tVertexFrange *> aVTarget;

    for (int aX=aX0 ; aX != aX1 ; aX++)
    {
        aMatVertex.push_back(std::vector<tVertexFrange *>());
        for (int aY=aY0 ; aY!=aY1 ; aY+= aSignY )
        {
            cPt2di aPix(aX,aY);
            cASom anAttr(aPix);

            auto  aNewVertex = mGraph->NewVertex(anAttr);
            aMatVertex.back().push_back(aNewVertex);
            if ((aX==aX0) && (aY==aY0))
                aVertexInit = aNewVertex;

            if ((aX==aX1-1) || (aY==aY1-aSignY))
                aVTarget.push_back(aNewVertex);

            cPt3di aCoul = mImVisuPrgDyn.GetRGBPix(aPix);
            if (aSignY>0)
                aCoul.x() = 255 - aCoul.x();
            else
                aCoul.y() = 255-aCoul.y();
            mImVisuPrgDyn.SetRGBPixWithAlpha(aPix,aCoul,cPt3dr(0.8,0.8,0.8));
        }
   }

   int aNbX = aX1-aX0;
   int aNbY = (aY1-aY0) * aSignY;

   std::vector<cPt2di> aVNeigh{{1,0},{1,aSignY},{0,aSignY},{-1,aSignY}};
   bool GotV0  = false;
   for (int aX=aX0 ; aX != aX1 ; aX++)
   {
       for (int aY=aY0 ; aY!=aY1 ; aY+= aSignY )
       {
           int aIndX0 = aX-aX0;
           int aIndY0 = (aY-aY0)*aSignY;
           auto aVertex0 =    aMatVertex.at(aIndX0).at(aIndY0);

           if (aVertex0== aVertexInit)
               GotV0 = true;
           for (const auto aNeigh : aVNeigh)
           {
               int aIndX1 = aIndX0 + aNeigh.x();
               // we must re-invert to be coherent with previous, and separation dir-init
               int aIndY1 = aIndY0 + aNeigh.y() * aSignY;
               if ((aIndX1>=0) && (aIndX1<aNbX) && (aIndY1>=0) && (aIndY1<aNbY))
               {
                   auto aVertex1 =    aMatVertex.at(aIndX1).at(aIndY1);
                   cAEdgeSym anAttrSym;
                   cAEdgeOriented anAttrOrient_01(Cost(aVertex0,aVertex1,aSignY));
                   cAEdgeOriented anAttrOrient_10(1e5); // Cost(aVertex1,aVertex0,aSignY));


                   mGraph->AddEdge
                   (
                        *aVertex0,*aVertex1,
                        anAttrOrient_01,anAttrOrient_10,anAttrSym
                   );
               }
           }
       }
    }

    {
       cPt2di aP0(aX0,mYCenter);
         for (const auto anE : aMatVertex[1][1]->EdgesSucc() )
         {
             if (anE->IsDirInit())
             {
                StdOut() <<  anE->VertexInit().Attr().mPt - aP0
                         << " => "
                         <<    anE->Succ().Attr().mPt - aP0
                         << " \n";
             }
         }
    }


    cAlgoSP<tGraphFrange> anAlgo;
    cSubGraphOfVertices<tGraphFrange> aGrTarget(*mGraph,aVTarget);

    tVertexFrange * anEnd = anAlgo.MakeShortestPathGen
                              (
                                  *mGraph,
                                  true,
                                  {aVertexInit},
                                  cWeithGr(),
                                  aGrTarget
                              );

    StdOut() <<  "  VVVV " << anEnd << " GV0=" << GotV0 << "\n";
    std::vector<tVertexFrange*> aPath;

    if (anEnd)
        anEnd->BackTrackFathersPath(aPath);


    for (const auto & aVert : aPath)
    {
        mImVisuPrgDyn.SetRGBPix(aVert->Attr().mPt,(aSignY>0) ? cRGBImage::Blue : cRGBImage::Cyan);
    }

    StdOut()  << "SOMM0 : " << aVertexInit->Attr().mPt  << " NB=" << aVTarget.size() << " PATH=" << aPath.size() << "\n";
    delete mGraph;
    mGraph = nullptr;
}

void  cAppliFranges::ComputeCurvePrgDyn(int aX0)
{
    HalfComputeCurvePrgDyn(aX0,1);
    HalfComputeCurvePrgDyn(aX0,-1);
    StdOut() << " \n=============================================\n\n";
}




tREAL8 cAppliFranges::ScoreParabol(tREAL8 aX0,int aIntY,tREAL8 aXAtLim) const
{
    tREAL8 aSom = 0;
    for (int aDY= -aIntY ; aDY<= aIntY ; aDY++)
    {
        tREAL8 aX = aX0  + Square(aDY/tREAL8(aIntY)) * aXAtLim;
        cPt2dr aPt(aX,mYCenter+aDY);
        if (! mDImBlur->InsideBL(aPt))
            return -1e10;

        aSom += mDImBlur->GetVBL(aPt);
    }

    return aSom;
}

/*
class cFitParabol
{
   public :
    cFitParabol (tREAL8 anYC,tREAL8 anInY,tDIm *,tREAL8 aStep=1.0);
    void SetX(tREAL8 anX,tREAL8 anXAtLim);

    tREAL8 Score(tREAL8 anX0,tREAL8 aXAtLim) ;

    cPt2dr PtOfDY(tREAL8 anY) const;
    */

void   cAppliFranges::ComputeParabol(tREAL8 aX0,tREAL8 aYRel0,tREAL8 aYRel1)
{
    int aNbY = 20;

    cWhichMax<cPt2dr,tREAL8> aWMax;
    cFitParabol aFitP(mYCenter,aNbY,mDImBlur,1.0);

    for (int aDX=-4 ; aDX<=4 ; aDX+=2)
    {
        for (int aXAtLim =0 ; aXAtLim < 800 ; aXAtLim++)
        {
            cPt2dr aParam(aX0+aDX,aXAtLim);
            tREAL8 aScore = aFitP.Score(aX0+aDX,aXAtLim,aYRel0,aYRel1);
            aWMax.Add(aParam,aScore);
        }
    }

    StdOut() << "SCORE " << aX0 << " =" << aWMax.ValExtre()  << " " << aWMax.IndexExtre() << "\n";

     cPt2dr aParam =  aWMax.IndexExtre();
     aFitP.SetX(aParam.x(),aParam.y());

     for (int aY= round_ni(2*aNbY*aYRel0) ; aY<= round_ni(2*aNbY*aYRel1) ; aY++)
     {
         cPt2dr aPt = aFitP.PtOfDY(aY);
         if (mImVisuPrgDyn.InsideBL(aPt))
         {
            mImVisuPrgDyn.SetRGBPix(ToI(aPt),(aYRel0==0) ? cRGBImage::Orange :  cRGBImage::Magenta);
         }
     }
}

void   cAppliFranges::ComputeParabol(tREAL8 aX0)
{
    ComputeParabol(aX0,-1,0);
    ComputeParabol(aX0,0,1);

}


void cAppliFranges::CalcAxeSym()
{
    StdOut()  << " cAppliFranges::CalcAxeSym...\n";

    mFactRed = mCalcSym.mFactRed;

    mImRed = mCurIm.BiCubicDeZoom(mFactRed,1.0);
    mDImRed = &mImRed.DIm();
    mSzImRed = mImRed.DIm().Sz();

    mImBlur = mImRed.Dup();
    mDImBlur = & mImBlur.DIm();
    ExpFilterOfStdDev(*mDImBlur,5,2.0,4.0);
    mDImBlur->ToFile("toto-blured.tif");

    mWidthCS = std::max(3,round_up(mCalcSym.mWitdhCalc / mFactRed));
    //tREAL8 aDisFS = mCalcSym.mDistFrange /


    if (mDoVisu)
    {
        tElIm aMin,aMax;
        GetBounds(aMin,aMax,mImRed.DIm());
         mImVisuPrgDyn = cRGBImage(mSzImRed);

         for (const auto & aPix : mImRed.DIm())
         {
             mImVisuPrgDyn.SetGrayPix(aPix,mImRed.DIm().GetV(aPix) * (255.0/aMax));
         }
     }

    //  ======================== Compute Y of symetry ==================
    int aY0 = mWidthCS;
    int aY1 = mSzImRed.y()-mWidthCS;


    cWhichMin<int,tREAL8> aMinSym;
    for (int anYC=aY0 ; anYC<aY1 ; anYC++)
    {
        cSymMeasure<tREAL8>  aMesSym;
        for (int anX=0 ; anX<mSzImRed.x() ; anX++)
        {
            for (int aDY=0 ; aDY<mWidthCS ; aDY++)
            {
                cPt2di aP1(anX,anYC+aDY);
                cPt2di aP2(anX,anYC-aDY);

                tREAL8 aWeight =  std::abs(mWidthCS-aDY)  / mWidthCS;
                aMesSym.Add
                (
                    aWeight,
                    mDImRed->GetV(aP1),
                    mDImRed->GetV(aP2)
                );
            }
        }
        aMinSym.Add(anYC,aMesSym.Sym(1e-2));
    }

     mYCenter =  aMinSym.IndexExtre() ;

    StdOut() << "YSIM=" << mYCenter << "\n";

    //  ============  compute threshold  ========
    mLowVal = -1;
    mVeryLowVal = -1;
    mHighVal = -1;
    std::vector<tREAL8> aVecMed;
    std::vector<int>  aVecMax;
    {
       cWhichMax<tINT4,tREAL8>  aMaxMed;
       cRect2 aRectN = cRect2::BoxWindow(cPt2di(0,0),cPt2di(1,10));
       for (int aX=0 ; aX<mSzImRed.x() ; aX++)
       {
           cPt2di aPt(aX,mYCenter);
           tREAL8 aValMed = KthValLoc(*mDImRed,aPt,aRectN,0.5,tREAL4(-1.0));
           aVecMed.push_back(aValMed);
           aMaxMed.Add(aX,aValMed);
       }

       mLowVal = Cst_KthVal(aVecMed,0.5);
       mVeryLowVal = Cst_KthVal(aVecMed,0.2);

       mHighVal = aMaxMed.ValExtre();
       aVecMax.push_back(aMaxMed.IndexExtre());
    }

    StdOut() << " VeryLow=" << mVeryLowVal << " Low=" << mLowVal << " High=" << mHighVal << "\n";

    // ===================   Compute All centers ===================
    {
        bool doCont = true;
        tREAL8 aThrVal = (mLowVal+mHighVal) / 2.0;
        int aThrDist = (mCalcSym.mDistFrange / mFactRed) * 0.5;

        while (doCont)
        {
            cWhichMax<tINT4,tREAL8>  aBestNext;

            for (int aX=0 ; aX<mSzImRed.x() ; aX++)
            {
                if (aVecMed.at(aX) > aThrVal )
                {
                     bool isOk = true;
                     for (const auto & aMax : aVecMax)
                     {
                         if (std::abs(aMax-aX) < aThrDist)
                         {
                            isOk=false;
                         }
                     }

                    if (isOk)
                       aBestNext.Add(aX,aVecMed.at(aX));
                }
            }
            if (aBestNext.IsInit())
            {
                aVecMax.push_back(aBestNext.IndexExtre());
            }
            else
                doCont = false;
        }
    }




     if (mDoVisu)
     {
         for (int anX=0 ; anX<mSzImRed.x() ; anX++)
         {
             mImVisuPrgDyn.SetRGBPix(cPt2di(anX,mYCenter),cRGBImage::Magenta);
         }

         for (const auto & aX : aVecMax)
         {
             ComputeParabol(aX);
             ComputeCurvePrgDyn(aX);
         }

         for (const auto & aXC : aVecMax)
         {
             mImVisuPrgDyn.DrawCircle(cRGBImage::Red,cPt2dr(aXC,mYCenter),3.0);
         }




         mImVisuPrgDyn.ToFile("toto.tif");

         cIm2D<tREAL4> aImCost(mSzImRed);
         for (const auto & aPix : *mDImRed)
         {
             aImCost.DIm().SetV(aPix,Cost(mDImRed->GetV(aPix)));
         }
         aImCost.DIm().ToFile("tutu_Cost.tif");

//         aImRed.DIm().ToFile("toto.tif");
    }

}


void  cAppliFranges::DoOneImage(const std::string & aNameIm)
{
    mNameCurIm = aNameIm;
    mIdExport = "Frange-MMVII-" + Prefix(mNameCurIm);
    InitReportCSV (mIdExport,"csv",false,{"Label","X","Y","Weight"});



    //  ============== Create & read the images
    cDataFileIm2D aDataFI2D = cDataFileIm2D::Create(mNameCurIm,eForceGray::No);

    UpdateMax(mIntY.x(),0);
    UpdateMin(mIntY.y(),aDataFI2D.Sz().y());

    cBox2di aBox(cPt2di(0,mIntY.x()), cPt2di(aDataFI2D.Sz().x(),mIntY.y()));
   // mCurIm = tIm::FromFile(mNameCurIm,aBox);
    mCurSz =  aBox.Sz();
    ReadImage(aBox);
    mDCurIm = & mCurIm.DIm();


    if (mMode==1)
    {
        CalcAxeSym();
        return;
    }

    mImSmooth = mCurIm.Dup();
    mDImSmooth= &mImSmooth.DIm();
    mImLabel = tImLabel(mCurSz);
    mDImLabel = &mImLabel.DIm();

    // =========== create smoothing of images ==========================
    tREAL8 aFactExpX = FactExpFromSigma2(Square(mSigma.x())/mNbIter);
    tREAL8 aFactExpY = FactExpFromSigma2(Square(mSigma.y())/mNbIter);
    ExponentialFilter(true,*mDImSmooth,mNbIter,*mDImSmooth,aFactExpX,aFactExpY);

    // =================== init labeling by thresholding =========
    mThreshW = mDCurIm->MoyVal();
    for ( auto& aPix : *mDImLabel )
    {
        bool isOver = mDImSmooth->GetV(aPix)>mThreshW;
        mDImLabel->SetV(aPix,isOver ? 255 : 0);
    }
    mDImLabel->InitBorder(0); // Current precaution to avoid image outside

    // ============== Now Put in mDImSmooth , Convol on Y, or Original
    mDCurIm->DupIn(*mDImSmooth);
    if (IsInit(&mSigCurv))
    {
        tREAL8 aFactExpY = FactExpFromSigma2(Square(mSigCurv)/mNbIter);
        ExponentialFilter(true,*mDImSmooth,mNbIter,*mDImSmooth,0,aFactExpY);
    }



    // ====   Create the curves of connected component analysis =====
    mMarqCurCC = 2;
    for (const auto& aPix : *mDImLabel )
    {
         if (mDImLabel->GetV(aPix) == 255)
             AnalyseOneConnectedComp(aPix);
    }
    // sort the curve to have them from left to right, using operator <
    std::sort (mCurves.begin(),mCurves.end());

    // ===============  Generate the export =====================

    // Generate the csv
    std::vector<cExportFrange> aVCE;
    for (size_t aKC=0 ; aKC<mCurves.size() ; aKC++)
    {
        cExportFrange anEF;
        anEF.mLabel = aKC;
        anEF.mVPts = mCurves.at(aKC).ExtractCurve(cPt2dr(0,mIntY.x()));
        for (const auto aPt : anEF.mVPts)
        {
            AddOneReportCSV(mIdExport,{ToStr(aKC),ToStr(aPt.x()),ToStr(aPt.y()),ToStr(aPt.z())});
        }
        aVCE.push_back(anEF);
    }
    SaveInFile(aVCE,DirReport()+mIdExport+".json");

    // Generate the visualisation if required
    if (mDoVisu)
       MakeVisu();
}


void cAppliFranges::MakeVisu()
{
    std::string aPref = mPhProj.DirVisuAppli() + mIdExport ;

    mDImSmooth->ToFile(aPref+"-Smooth.tif");
    mDImLabel->ToFile(aPref+"-Label.tif");

    cRGBImage aIRGB(mCurSz);

    tIm aImGray(mCurSz);
    for (const auto& aPix : *mDImLabel)
    {
        tREAL8 aVal = 0;
        if (mDImLabel->GetV(aPix) >1)
        {
           aVal = WeightOfGray(mDImSmooth->GetV(aPix)) * 255;
        }
        aImGray.DIm().SetV(aPix,aVal);
        aIRGB.SetGrayPix(aPix,aVal);
    }
    for (const auto & aFr : mCurves)
    {
        for (const auto & aPt : aFr.ExtractCurve(cPt2dr(0,0)))
        {
             aIRGB.SetRGBPix(ToI(Proj(aPt)),cRGBImage::Red);
        }
    }

    aIRGB.ToFile(aPref+"-Curves.tif");
    aImGray.DIm().ToFile(aPref+"-GrayLab.tif");
}
};

using  namespace NS_FrangesDetect;

// ============================= Old version ===============

tMMVII_UnikPApli Alloc_AppliFranges_0(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli & aSpec)
{
      return tMMVII_UnikPApli(new cAppliFranges(aVArgs,aSpec,0));
}


cSpecMMVII_Appli  TheSpecAppliFranges_0
(
     "ExtractFranges_0",
      Alloc_AppliFranges_0,
      "Extraction of Franges, version 0 (image filter)",
      {eApF::ImProc},
      {eApDT::Image},
      {eApDT::Console},
      __FILE__
);


// ============================= New version ===============

tMMVII_UnikPApli Alloc_AppliFranges_1(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli & aSpec)
{
      return tMMVII_UnikPApli(new cAppliFranges(aVArgs,aSpec,1));
}


cSpecMMVII_Appli  TheSpecAppliFranges_1
(
     "ExtractFranges_1",
      Alloc_AppliFranges_1,
      "Extraction of Franges, version 1 (prog dyn, still unachieved)",
      {eApF::ImProc},
      {eApDT::Image},
      {eApDT::Console},
      __FILE__
);


#if (0)
#endif

};
