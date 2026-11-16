//#include "MMVII_PCSens.h"

// Test commit

#include "cAppliNewFranges.h"




namespace MMVII
{

namespace NS_FrangesDetect
{

#if (0)

/** Add heuristically some horizontal CC :
 *    - they are on center
 *    - they are max loc on the center line
 *    - they are not close to existing CC
*/
void cAppliNewFrange::ComplemCCHoriz()
{
    tREAL8 aPropLow = 0.6;
    tREAL8 aThreshold = ((aPropLow*mRadiomBackGround) + (1-aPropLow) * mDImRadFrange->GetV(mYC));
    int aNbMaxLoc = 40;
    tREAL8 aDistCC = 100.0;
    //int aNbX = mSzRed.x();


    for (int anX=aNbMaxLoc ; anX< mSzRed.x()-1-aNbMaxLoc ; anX++)
    {
        // Over a thrashold
        cPt2di aPt(anX,mYC);
        tREAL8 aVC = mDImRedBlur->GetV(aPt) ;
        bool isMaxLoc = aVC > aThreshold;

        // is a local maxima
        if (isMaxLoc)
        {
            for (int aDx=1 ; (aDx<aNbMaxLoc) && isMaxLoc ; aDx++)
            {
                tREAL8 aVPDx =  mDImRedBlur->GetV(cPt2di(anX+aDx,mYC));
                tREAL8 aVPMx =  mDImRedBlur->GetV(cPt2di(anX-aDx,mYC));

                 isMaxLoc = isMaxLoc && (aVC>=aVPMx) && (aVC>aVPDx);
            }
        }

        // is not close to existing max loc
        if (isMaxLoc)
        {
           for (size_t aKC=0 ; (aKC<mVectCC.size()) && isMaxLoc ; aKC++)
           {
               const cConnComp & aCC = *mVectCC.at(aKC);
               if (aCC.mIsHor && aCC.mIsOk)
               {
                   cSegment2DCompiled<tREAL8> aSeg(aCC.mSeg);
                   tREAL8 aD = aSeg.DistClosedSeg(ToR(aPt));
                   isMaxLoc = isMaxLoc && (aD>aDistCC);
               }
           }
        }


        // Add it by prolongatin for center as long as it is over a threshold
        if (isMaxLoc)
        {
            int aBorder = std::min(mSzRed.y()/4,aNbMaxLoc);
            int anY0 = mYC;
            while ((anY0>aBorder) &&  (mDImRedBlur->GetV(cPt2di(anX,anY0))>aThreshold) )
                anY0--;

            int anY1 = mYC;
            while (anY1<mSzRed.y() - aBorder &&  (mDImRedBlur->GetV(cPt2di(anX,anY1))>aThreshold) )
                anY1++;

            std::vector<cPt2di>  aVPt;
            for (int anY=anY0 ; anY<=anY1 ; anY++)
                aVPt.push_back(cPt2di(anX,anY));

            tSeg2dr aSeg(cPt2dr(anX,anY0),cPt2dr(anX,anY1));
            cConnComp* aCC = new cConnComp(aVPt,aSeg,cPt2dr(0,-1),cPt2dr(0,1),true,1.0,true);
            mVectCC.push_back(aCC);
        }

    }
}
#endif

};  // NS_FrangesDetect



}; //  namespace MMVII
