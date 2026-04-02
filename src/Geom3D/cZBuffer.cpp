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
                    cTri3DIterator * aMesh2D) :
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
    mMaxStretching(5.0)
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

    cPt3dr aPtZ(aTri3.Pt(0).z(),aTri3.Pt(1).z(),aTri3.Pt(2).z());

    std::vector<cPt2di> aVPix;
    std::vector<cPt3dr> aVW;


    cPt3dr aNorm = Normal(aTri3);

    int aSign = (aNorm.z() > 0) ? 1 : - 1;
     ///  the axe K of camera is in direction of view, the normal is in direction of visibility => they are opposite
     ///  the axe K of camera is in direction of view, the normal is in direction of visibility => they are opposite
    bool WellOriented =  mZF_SameOri ?  (aSign>0)  :(aSign<0);

    aTri2.PixelsInside(aVPix,1e-8,&aVW);
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
     cPt3dr aP12 = aTriIn.Pt(2) - aTriIn.Pt(1);

     cPt3dr aP01_3d = aTri3.Pt(1) - aTri3.Pt(0);
     cPt3dr aP12_3d = aTri3.Pt(2) - aTri3.Pt(1);

     // Gram-Schmidt orthonormal basis definition 
     double aP01norm = Norm1(aP01);
     cPt3dr aP01N = aP01/aP01norm;
     double aX2 = Scal(aP12,aP01N);
     double aY2 = Norm1( aP12 - (aP01N*aX2 )) ;

    /* cDenseMatrix<tREAL8> aDm_1= M2x2FromCol(cPt2dr(aP01norm,0.0), 
                                                cPt2dr(aX2,aY2)).Inverse();*/

    
    // compute the inverse analytically 
    tREAL8 aDet = std::abs(aP01norm*aY2);

    if( aDet<=0)
        return true;

    cDenseMatrix<tREAL8> aDm_1 = M2x2FromCol(cPt2dr(aY2,-aX2)/aDet,
                                            cPt2dr(0.0,aP01norm)/aDet);

    // compute the Jacobian Matrix 
    /*cDenseMatrix<tREAL8>  aDS(2,3); // col, row
    SetCol(aDS,0, aP01_3d);
    SetCol(aDS,1, aP12_3d);*/

    // Jacobian J = aDS * aDm^-1 manually
    cPt3dr aJ01 = aP01_3d * aDm_1.GetElem(cPt2di(0,0)) + aP12_3d * aDm_1.GetElem(cPt2di(0,1));
    cPt3dr aJ12 = aP01_3d * aDm_1.GetElem(cPt2di(1,0)) + aP12_3d * aDm_1.GetElem(cPt2di(1,1));

    // JT J manually

    /*
         |aa bb|
    JTJ =|     |
         |bb cc|
    */

    tREAL8 aa = Scal(aJ01,aJ01);
    tREAL8 bb = Scal(aJ01,aJ12);
    tREAL8 cc = Scal(aJ12,aJ12);

    tREAL8 trace = aa + cc;
    tREAL8 det = aa*cc - bb*bb;
    tREAL8 temp = std::sqrt(std::max(0.0, trace*trace * 0.25 - det));
    

    // eigen values
    tREAL8 lambda1 = trace * 0.5 + temp;
    tREAL8 lambda2 = trace * 0.5 - temp;
    
    tREAL8 sigma1 = std::sqrt(lambda1);
    tREAL8 sigma2 = std::sqrt(lambda2);

    aStretchingThreshold = (sigma1/sigma2) ;

    // dirichlet energy 
    //aStretchingThreshold= (sigma1*sigma1)+ (sigma2*sigma2)+ 1/(sigma1*sigma1) + 1/(sigma2*sigma2);

    if( (sigma2<1e-6) || (aStretchingThreshold > mMaxStretching))
        return true;
    else
    {
        // keep stretching to 1.0
        aStretchingThreshold = 1.0;
    }

    return false;


    /*

    auto aJac = aDS * aDm_1;




    auto aJacTaJac = aJac.Transpose()* aJac;

    // eigen val
    const cDenseVect<tREAL8> anEV = aJacTaJac.Eigen_Decomposition().mEigenVal_R;

    //evaluate stretching level 
    double aRatio = anEV(0)/anEV(1);
    //double anArea = anEV(0)*anEV(1);
    //StdOut()<<"aRatio "<<aRatio<<std::endl;
    if( aRatio > mMaxStretching) return true;

    // symmetric dirichlet energy 

    if(anEV(1)<1e-6)
        return true ;

    if( aRatio > mMaxStretching) 
        return true;

    return false;

    */

}

};
