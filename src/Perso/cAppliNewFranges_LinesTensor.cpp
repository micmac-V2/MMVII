//#include "MMVII_PCSens.h"

// Test commit

#include "cAppliNewFranges.h"




namespace MMVII
{


namespace NS_FrangesDetect
{


bool cAppliNewFrange::IsSaddlePoint(const cPt2di & aPix) const
{
    if (mDImRedBlur->Insideness(aPix)<=1)
        return false;
   return  FlagSup8Neigh(*mDImRedBlur,aPix).NbConComp() >=4;
}


cPt2dr cAppliNewFrange::TensorImage(cPt2di aPix) const
{
   return ToR(cTensor<tElIm>::Tens2Grad(mImTens.Grad(aPix),1.0));
}

cPt2dr cAppliNewFrange::TensorTeta(cPt2di aPix) const
{
   return FromPolar(1.0,mDImTeta->GetV(aPix.y()));
}




cPt2dr cAppliNewFrange::NextExtreLine(const cPt2dr & aP1,const cPt2dr& aP2) const
{                  
    // compute once for all Direction arround P1->P2 that must be explored
    static std::vector<cPt2dr> aVDir;
    if (aVDir.empty())
    {
        tREAL8 aStep = 0.15;
        tREAL8 aTetaMax = 0.3;
        int aNbTeta = round_up(aTetaMax/aStep);
        aStep = aTetaMax / aNbTeta;
        for (int aKTeta=-aNbTeta; aKTeta<=aNbTeta ; aKTeta++)
        {
            aVDir.push_back(FromPolar(1.0,aStep*aKTeta));
        }
    }

    // extract maxima of image, near to P2, in direction close to P1->P2
    cWhichMax<cPt2dr,tREAL8> aMaxP;

    cPt2dr aDirNext = VUnit(aP2-aP1);
    for (const auto&  aDir : aVDir)
    {
        cPt2dr aNewP = aP2+aDir * aDirNext;
        tREAL8 aVal = mDImRedBlur->DefGetVBL(aNewP,-1e8);
        aMaxP.Add(aNewP,aVal);
    }

    return aMaxP.IndexExtre();
}


void cAppliNewFrange:: ProlongateExtre(std::vector<cPt2dr>  & aRes, cPt2di aPix,eTypePix aLabel) const
{
    cPt2dr aRPix = ToR(aPix);
    aRes.clear();

    // first point must be in the set
    if (! LabelIsOn(aPix,aLabel))
        return;

    // check it's an extremity by computing the number of CC in neighboor
    {
        cFlag8Neigh aFlag8N;
        for (size_t aKN=0 ; aKN<8 ; aKN++)
        {
              if (LabelIsOn(aPix+FreemanV8[aKN],aLabel))
                  aFlag8N.AddNeigh(aKN);
        }
        if (aFlag8N.NbConComp()!=2)
            return;
    }

    // we want to avoid to grow to high, else small parasite may be linked to real line
    // so we compute a threshold we must not go over
    tREAL8 aVal = mDImRedBlur->GetV(aPix);
    tREAL8 aThrs = aVal + (aVal-mRadiomBackGround);


    aRes.push_back(aRPix);
    cPt2dr aDirT = TensorImage(aPix);

    cWhichMax<cPt2dr,tREAL8> aClosTens;
    for (const auto& aDN8 : Alloc8Neighbourhood())
    {
        cPt2di aN8 = aPix + aDN8;
        if ( LabelIsOn(aN8,aLabel))
        {
            aClosTens.Add(ToR(aDN8),std::abs(Scal(aDirT,VUnit(ToR(aDN8)))));
        }
    }
    // Note the - sign as we want to go in direction exterior to the line
    cPt2dr aDelta =  VUnit(-aClosTens.IndexExtre());
    //cPt2dr aNext = aRPix + aDelta;
    // Refin poistion of firs point (as it was limited to 8 neigh)
    cPt2dr aNext =    NextExtreLine(ToR(aPix)-aDelta,ToR(aPix));

    aRes.push_back(aNext);

    for (int aK=0 ; aK<6 ; aK++)
    {
        cPt2dr aNew = NextExtreLine(aRes.at(aRes.size()-2),aRes.at(aRes.size()-1));
        aRes.push_back(aNew);

        cPt2di aINew = ToI(aNew);

        // if we are dangerously close to to border or to high this is bad
        if ((mDImLabel->Insideness(aINew)<3) ||  (mDImRedBlur->GetV(aINew) > aThrs))
        {
            aRes.clear();
            return ;
        }

        // if we are "touching" the label, we have reached our goal=connecting, so can stop
        for (const auto & aNeigh : Alloc8Neighbourhood())
        {
             if (LabelIsOn(aINew+aNeigh,aLabel))
                 return ;
        }
    }

    // to many iteration w/o reaching something : useless
    aRes.clear();
}



bool cAppliNewFrange::IsExtreLocTens(const cPt2di & aPix,const cPt2dr & aDir,tREAL8 aSignVal) const
{

    tREAL8 aV0 = mDImRedBlur->GetV(aPix)*aSignVal;

    tREAL8 aContrastMin = (mDImRadFrange->GetV(aPix.y()) - mRadiomBackGround) * 0.05 * aSignVal;
    tREAL8  aVLow = aV0 - aContrastMin; // on each side, at least one value will have to be > aVLow


    // Contrats will be ok when OK on 2 side
    bool hasEnoughContrast[2]= {false,false};

    int aNbContrast = 15;  // number of step for testing enough contrat
    int aNbExtreLoc =   5;  // number of step for max loc, not to high to avoid joining next line


    for (int aK=1 ; aK<aNbContrast ; aK++)
    {
        for (int aKSign =0; aKSign<2 ; aKSign++) // Explorate "Left & Right"
        {
            int aSignDir = (aKSign==0) ? 1 : -1;
            cPt2dr aP1 = ToR(aPix) + aDir * (tREAL8)(aK * aSignDir);

            if (mDImRedBlur->InsideBL(aP1))
            {
                tREAL8 aV =  mDImRedBlur->GetVBL(aP1) * aSignVal;
                hasEnoughContrast[aKSign] |= (aV<aVLow);

                if (aK<aNbExtreLoc) // test for being local extrema
                {
                    if ((aV>aV0) || (aK && (aV==aV0)) )
                       return false;
                }
                else
                {
                    // if we are here, local extrem is OK, and OK as soon as we get enough contrast
                    if (hasEnoughContrast[0] && hasEnoughContrast[1])
                        return true;
                }
            }
        }
    }


    return false;
}



void cAppliNewFrange::ExtractTensorLines(bool isMax)
{
     eTypePix aLabel = isMax ? eTypePix::eMaxTens : eTypePix::eMinTens;
     tREAL8 aSign = isMax ? 1.0 : -1.0;

     // tentative extraction of centers of lines as maxima of grey
     // in direction of tensor
     int aBorder = 5;
     //std::vector<cPt2dr> aVN;
     std::vector<cPt2di> aVPixMax;

     for (const auto & aPix : mDImLabel->Interior(aBorder))
     {
        // cPt2dr aDirTens = ToR(cTensor<tElIm>::Tens2Grad(mImTens.Grad(aPix),1.0));

        cPt2dr aDirTens =  isMax ? TensorTeta(aPix) : cPt2dr(0,1) ; // FromPolar(1.0,mDImTeta->GetV(aPix.y()));

         if (IsExtreLocTens(aPix,aDirTens,aSign))
         {
            aVPixMax.push_back(aPix);
            SetLabelOn(aPix,aLabel);
          }
     }

     for (const auto & aPixMax : aVPixMax)
     {
         std::vector<cPt2dr> aVProl;
         ProlongateExtre(aVProl,aPixMax,aLabel);
         for (const auto & aProl : aVProl)
         {
             SetLabelOn(ToI(aProl),aLabel);
             SetLabelOn(ToI(aProl),eTypePix::eProlong);
         }
     }

}

void cAppliNewFrange::ExtractTensorLines()
{
    ExtractTensorLines(true);
    ExtractTensorLines(false);
}





};  // NS_FrangesDetect
}; //  namespace MMVII

