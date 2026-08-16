#include "MMVII_Ptxd.h"
#include "MMVII_Image2D.h"
#include "MMVII_Geom2D.h"
#include "MMVII_Geom3D.h"
#include "MMVII_Triangles.h"
#include "MMVII_Mappings.h"

#include "MMVII_ZBuffer.h"


namespace MMVII
{

bool  ZBufLabIsOk(eZBufRes aLab)
{
     return (aLab==eZBufRes::Visible) || (aLab==eZBufRes::LikelyVisible) ;
}


/* =============================================== */
/*                                                 */
/*                 cTri3DIterator                  */
/*                                                 */
/* =============================================== */

void cTri3DIterator::ResetAll()
{
    ResetTri();
    ResetPts();
}

cCountTri3DIterator * cTri3DIterator::CastCount() {return nullptr;}

/* =============================================== */
/*                                                 */
/*            cCountTri3DIterator                  */
/*                                                 */
/* =============================================== */

cCountTri3DIterator::cCountTri3DIterator(size_t aNbP,size_t aNbF) :
    mNbP  (aNbP),
    mNbF  (aNbF)
{
   ResetPts();
   ResetTri();
}

void cCountTri3DIterator::ResetTri() { mIndexF=0;}
void cCountTri3DIterator::ResetPts() { mIndexP=0;}

bool cCountTri3DIterator::GetNextPoint(cPt3dr & aP )
{
    if (mIndexP>=mNbP) return false;
    aP = KthP(mIndexP);
    mIndexP++;
    return true;
}

bool cCountTri3DIterator::GetNextTri(tTri3dr & aTri)
{
    if (mIndexF>=mNbF) return false;
    aTri = KthF(mIndexF);
    mIndexF++;
    return true;
}

cCountTri3DIterator * cCountTri3DIterator::CastCount() {return this;}

/* =============================================== */
/*                                                 */
/*              cMeshTri3DIterator                 */
/*                                                 */
/* =============================================== */

cMeshTri3DIterator::cMeshTri3DIterator(cTriangulation3D<tREAL8> * aTri) :
    cCountTri3DIterator(aTri->NbPts(),aTri->NbFace()),
    mTri (aTri)
{
}

cPt3dr  cMeshTri3DIterator::KthP(int aKP) const {return mTri->KthPts(aKP);}
tTri3dr cMeshTri3DIterator::KthF(int aKF) const {return mTri->KthTri(aKF);}

/* =============================================== */
/*                                                 */
/*              cResModeSurfD                      */
/*                                                 */
/* =============================================== */

void  AddData(const cAuxAr2007  &anAux,cResModeSurfD& aRMS )
{
     int aResult = int(aRMS.mResult);
     AddData(cAuxAr2007("Result",anAux),aResult);
     if (anAux.Input())
        aRMS.mResult = eZBufRes(aResult);
     AddData(cAuxAr2007("Resol",anAux),aRMS.mResol);
}

/* =============================================== */
/*                                                 */
/*              cZBuffer                           */
/*                                                 */
/* =============================================== */


cZBuffer::cZBuffer(cTri3DIterator & aMesh,const tSet &  aSetIn,const tMap & aMapI2O,
                   const tSet &  aSetOut,double aResolOut,bool aSameOri, bool IsBascProc,
                    cTri3DIterator * aMesh2D, const cPt2di & aForcedOutSize) :
    mIsOk       (true),
    mZF_SameOri (aSameOri),
    mMultZ      (mZF_SameOri ? -1 : 1),
    mMesh       (aMesh),
    mMesh2DDepth(*aMesh2D),
    mCountMesh  (mMesh.CastCount()),
    mMapI2O     (aMapI2O),
    mSetIn      (aSetIn),
    mSetOut     (aSetOut),
    mResolOut   (aResolOut),
    mIsBasc     (IsBascProc),
    mBoxIn      (cBox3dr::Empty()),
    mBoxOut     (cBox3dr::Empty()),
    mROut2Pix   (),
    mZBufIm     (cPt2di(1,1)),
    mImSign     (cPt2di(1,1)),
    mMaxStretching(4.0)
{
    cTplBoxOfPts<tREAL8,3> aBoxOfPtsIn;
    cTplBoxOfPts<tREAL8,3> aBoxOfPtsOut;

    //  compute the box in put and output space
     //  compute the box in put and output space
    cPt3dr aPIn;

    mMesh.ResetAll();

    if (mIsBasc)
        {
            //cPt3dr aPOutReady;
            mMesh2DDepth.ResetAll();
            while (mMesh.GetNextPoint(aPIn))
            {
                cPt3dr aPOut = mMapI2O.Value(aPIn);
                /*MMVII_INTERNAL_ASSERT_strong(Norm1(aPOut-aPOutReady)<0.000001,
                                "ready points are not correctly mapped in mMesh2DDepth");*/

                aBoxOfPtsIn.Add(aPIn);
                aBoxOfPtsOut.Add(aPOut);
            }

            mMesh2DDepth.ResetPts();
        }
    else  ///< General case no precomputed 2D + Depth
        {
            while (mMesh.GetNextPoint(aPIn))
            {
                if (mSetIn.InsideWithBox(aPIn))
                {
                    cPt3dr aPOut = mMapI2O.Value(aPIn);

                    if (mSetOut.InsideWithBox(aPOut))
                    {
                        aBoxOfPtsIn.Add(aPIn);
                        aBoxOfPtsOut.Add(aPOut);
                    }
                }
            }
        }

    mMesh.ResetPts();

    if ((aBoxOfPtsIn.NbPts()<3) || (aBoxOfPtsOut.NbPts()<3))
    {
        mIsOk = false;
        return;
    }

    mBoxIn = aBoxOfPtsIn.CurBox();
    mBoxOut = aBoxOfPtsOut.CurBox();


    if (aForcedOutSize.x()>0 && aForcedOutSize.y()>0)
    {
        mBoxOut.P0ByRef() = cPt3dr(0,0, mBoxOut.P0().z());
        mBoxOut.P1ByRef() = cPt3dr(aForcedOutSize.x(),aForcedOutSize.y(), mBoxOut.P1().z());
    }


    cPt2di aBrd(2,2);
    //   aP0/aResout + aTr -> 1,1
    cPt2dr aTr = ToR(aBrd) - Proj(mBoxOut.P0()) * (1.0/mResolOut);
    mROut2Pix = cHomot2D<tREAL8>(aTr,1.0/mResolOut);

    mSzPix =  Pt_round_up(ToPix(mBoxOut.P1())) + aBrd;


    mZBufIm = tIm(mSzPix);
    mZBufIm.DIm().InitCste(mInfty);
    mImSign = tImSign(mSzPix,nullptr,eModeInitImage::eMIA_Null);
}

cPt2dr  cZBuffer::ToPix(const cPt3dr & aPt) const {return mROut2Pix.Value(Proj(aPt));}
cZBuffer::tIm  cZBuffer::ZBufIm() const {return mZBufIm;}
cResModeSurfD&  cZBuffer::ResSurfD(size_t aK)  {return mResSurfD.at(aK);}
double  cZBuffer::MaxRSD() const {return mMaxRSD;}

std::vector<cResModeSurfD> & cZBuffer::VecResSurfD() {return mResSurfD;}

void cZBuffer::AssertIsOk() const
{
   MMVII_INTERNAL_ASSERT_tiny(mIsOk,"Non ok Buffer");
}

const cPt2di & cZBuffer::SzPix() const
{
    return mSzPix;
}

bool cZBuffer::IsOk() const {return mIsOk;}

void cZBuffer::MakeZBuf(eZBufModeIter aMode)
{

    if (aMode==eZBufModeIter::SurfDevlpt)
    {
        mResSurfD.clear();
        mMaxRSD = 0.0;
    }

    tTri3dr  aTriIn = tTri3dr::Tri000();
    int aNbTriVis = 0;
    while (mMesh.GetNextTri(aTriIn))
    {
        mLastResSurfDev = -1;
        eZBufRes aRes = eZBufRes::Undefined;
        if (!mIsOk)
        {
        }
        //  not sure this is us to test that, or the user to assure it give clean data ...
        else if (aTriIn.Regularity() <=0)
           aRes = eZBufRes::UnRegIn;
        else if (! mSetIn.InsideWithBox(aTriIn))
           aRes = eZBufRes::OutIn;
        else
        {
            tTri3dr aTriOut = mMapI2O.TriValue(aTriIn);

            if (aTriOut.Regularity() <=0)
               aRes = eZBufRes::UnRegOut;
            else if (! mSetOut.InsideWithBox(aTriOut))
               aRes = eZBufRes::OutOut;
            else
            {
               aNbTriVis++;
               aRes = MakeOneTri(aTriIn,aTriOut,aMode);
            }
        }
        if (aMode==eZBufModeIter::SurfDevlpt)
        {
           cResModeSurfD aRMS;
           aRMS.mResult = aRes;
           aRMS.mResol  = mLastResSurfDev;
           mResSurfD.push_back(aRMS);
        }
    }
    if (aNbTriVis==0)
       mIsOk=false;

    mMesh.ResetTri();
}


void cZBuffer::MakeZBufForBasc(eZBufModeIter aMode)
{

    if (aMode==eZBufModeIter::SurfDevlpt)
    {
        mResSurfD.clear();
        mMaxRSD = 0.0;
    }

    tTri3dr  aTriIn  = tTri3dr::Tri000();
    tTri3dr  aTriOut = tTri3dr::Tri000();

    double aStretchingThreshold = 1.0 ;
    int aNbTriVis = 0;

    while (mMesh.GetNextTri(aTriIn) &&
           mMesh2DDepth.GetNextTri(aTriOut))
    {
        mLastResSurfDev = -1;
        eZBufRes aRes = eZBufRes::Undefined;
        if (!mIsOk)
        {
        }
        //  not sure this is us to test that, or the user to assure it give clean data ...
        else if (aTriIn.Regularity() <=0)
            aRes = eZBufRes::UnRegIn;
        else
        {

            if (aTriOut.Regularity() <=0)
                aRes = eZBufRes::UnRegOut;
            else
            {
                aNbTriVis++;
                aRes = MakeOneTri(aTriIn,aTriOut,aMode);
            }
        }
        if (aMode==eZBufModeIter::SurfDevlpt)
        {
            // change triangle class to distorted 
            
            if (IsStretched(aTriIn,aTriOut,aStretchingThreshold))
                aRes= eZBufRes::Distorted;

            cResModeSurfD aRMS;
            aRMS.mResult = aRes;
            aRMS.mResol  = mLastResSurfDev;
            aRMS.mStretchThresh=aStretchingThreshold;
            mResSurfD.push_back(aRMS);
        }

    }

    if (aNbTriVis==0)
        mIsOk=false;

    mMesh.ResetTri();
    mMesh2DDepth.ResetTri();
}


double cZBuffer::ComputeResol(const tTri3dr & aTri3In ,const tTri3dr & aTri3Out) const
{
        // input triangle, developped isometrically on the plane
        tTri2dr aTri2In  = cIsometry3D<tREAL8>::ToPlaneZ0(0,aTri3In,true);
        // output triangle, projected on the plane
        tTri2dr aTri2Out = Proj(aTri3Out);
        // Affinity  Input-Dev -> Output proj
        cAffin2D<tREAL8> aAffI2O =  cAffin2D<tREAL8>::Tri2Tri(aTri2In,aTri2Out);

        return aAffI2O.MinResolution();
}


eZBufRes cZBuffer::MakeOneTri(const tTri3dr & aTriIn,const tTri3dr &aTri3,eZBufModeIter  aMode)
{
    eZBufRes aRes = eZBufRes::Undefined;

    //  cTriangle2DCompiled<tREAL8>  aTri2(ToPix(aTri3.Pt(0)) , ToPix(aTri3.Pt(1)) ,ToPix(aTri3.Pt(2)));
    cTriangle2DCompiled<tREAL8>  aTri2 = ImageOfTri(Proj(aTri3),mROut2Pix);

    // do no use wrapped triangles
    int aIndexLongest = aTri2.IndexLongestSeg();
    float aMaxLen = Norm2(aTri2.KVect(aIndexLongest));
    if (aMaxLen > SzPix().x()/1.5)
        return aRes;

    cPt3dr aPtZ(aTri3.Pt(0).z(),aTri3.Pt(1).z(),aTri3.Pt(2).z());

    std::vector<cPt2di> aVPix;
    std::vector<cPt3dr> aVW;


    cPt3dr aNorm = Normal(aTri3);

    int aSign = (aNorm.z() > 0) ? 1 : - 1;
     ///  the axe K of camera is in direction of view, the normal is in direction of visibility => they are opposite
     ///  the axe K of camera is in direction of view, the normal is in direction of visibility => they are opposite
    bool WellOriented =  mZF_SameOri ?  (aSign>0)  :(aSign<0);

    aTri2.PixelsInside(aVPix,1e-9,&aVW);

    tDIm & aDZImB = mZBufIm.DIm();
    int aNbVis = 0;
    for (size_t aK=0 ; aK<aVPix.size() ; aK++)
    {
       const cPt2di  & aPix = aVPix[aK];
       tElem aNewZ = mMultZ * Scal(aPtZ,aVW[aK]);
       tElem aZCur = aDZImB.GetV(aPix);
       if (aMode==eZBufModeIter::ProjInit)
       {
           if (aNewZ> aZCur)
           {
               aDZImB.SetV(aPix,aNewZ);
           }
       }
       else
       {
           if (aNewZ==aZCur)
              aNbVis++;
       }
    }

    if (aMode==eZBufModeIter::SurfDevlpt)
    {
       if (! WellOriented)
          aRes =  eZBufRes::BadOriented;
       else
       {
           bool IsVis = ((aNbVis*2)  > int(aVPix.size()));
           aRes = IsVis ? eZBufRes::Visible : eZBufRes::Hidden;
           if(!mIsBasc) // mIsBasc= true -> no need to compute resol as it is just a mapping from depth to DEM
           {
                mLastResSurfDev = ComputeResol(aTriIn,aTri3);
               if (IsVis)
               {
                   UpdateMax(mMaxRSD,mLastResSurfDev);
               }
           }

           if ((aVPix.size()<=0) && (aNbVis==0))
              aRes = eZBufRes::NoPix;
       }
    }

    return aRes;
}



bool cZBuffer::IsStretched(const tTri3dr & aTriIn, const tTri3dr & aTri3, double & aStretchingThreshold)
{

     ///< Evaluate stretching of triangle where Mapping(aTri3) = aTriIn

     cPt3dr aP01 = aTriIn.Pt(1) - aTriIn.Pt(0);
     cPt3dr aP12 = aTriIn.Pt(2) - aTriIn.Pt(0);

     cPt3dr aP01_3d = aTri3.Pt(1) - aTri3.Pt(0);
     cPt3dr aP12_3d = aTri3.Pt(2) - aTri3.Pt(0);


    /*tREAL8 anInArea =  Norm1(aP01 ^ aP12);
    tREAL8 anOutArea =  Norm1(aP01_3d ^ aP12_3d);


    aStretchingThreshold = (anInArea/anOutArea) ;*/


    // create the jacobian of the deformation matrix
    cDenseMatrix<tREAL8> aMat_Out(2,3);
    MMVII::SetCol(aMat_Out,0,aP01);
    MMVII::SetCol(aMat_Out,1,aP12);

    cDenseMatrix<tREAL8> aMat_In(2,3);
    MMVII::SetCol(aMat_In,0,aP01_3d);
    MMVII::SetCol(aMat_In,1,aP12_3d);


    auto aC = aMat_In.Transpose()*aMat_In;
    auto aG = aMat_Out.Transpose()*aMat_Out;


    double A_ = aC.Det();
    double C_ = aG.Det();

    // double trCG  = c22*g11 - 2.0*c12*g12 + c11*g22;

    double tr_CG = aC.GetElem(cPt2di(1,1)) * aG.GetElem(cPt2di(0,0)) 
                    - 2.0 * aC.GetElem(cPt2di(0,1)) * aG.GetElem(cPt2di(0,1))
                    + aC.GetElem(cPt2di(0,0)) * aG.GetElem(cPt2di(1,1));

    double B_ = -tr_CG;

    //      |A_ B_|
    // MAT =|     |
    //      |B_ C_|

    double disc  = std::sqrt(std::max(0.0, B_*B_ - 4.0*A_*C_));

    // biggest eigenvalue
    double lambda_1    = (-B_ + disc) / (2.0*A_);
    double lambda_2    = (-B_ - disc) / (2.0*A_);

    double sigma1 = std::sqrt(std::max(0.0, lambda_1));
    double sigma2 = std::sqrt(std::max(0.0, lambda_2));

    aStretchingThreshold = (sigma1/sigma2) ;

    if( (sigma2<1e-6) ||  (aStretchingThreshold > mMaxStretching))
        return true;
    else
    {
        // keep stretching to 1.0
        aStretchingThreshold = 1.0;
    }

    return false;

}

};
