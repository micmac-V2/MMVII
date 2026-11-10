#include "MMVII_ImageMorphoMath.h"
namespace MMVII
{

template <class FoncBelong,class FoncUDdate,class FoncRestore>
void Tpl_ConnectedComponent
     (
        std::vector<cPt2di> & aVPts,
        const std::vector<cPt2di> & aNeighbourhood,
        const  std::vector<cPt2di> & aVecSeed,
        const FoncBelong& aFBelong,
        const FoncUDdate& aFUpate,
        bool doRestore,
        const FoncRestore & aFRestore
     )
{

//    std::cout << "Tpl_ConnectedComponentTpl_ConnectedComponentTpl_ConnectedComponent\n"; // getchar();
    for (const auto & aPSeed : aVecSeed)
    {
        //if (aDIm.GetV(aPSeed) == aMarqInit)
       if (aFBelong(aPSeed))
        {
           aFUpate(aPSeed);
//         aDIm.SetV(aPSeed,aNewMarq);
           aVPts.push_back(aPSeed);
        }
    }
    size_t aIndBottom = 0;

    while (aIndBottom!=aVPts.size())
    {
          cPt2di aP0 = aVPts[aIndBottom];
          for (const auto & aDelta : aNeighbourhood)
          {
              cPt2di aNeigh = aP0 + aDelta;
              //if (aDIm.GetV(aNeigh)==aMarqInit)
             if (aFBelong(aNeigh))
             {
                  aFUpate(aNeigh);
                 // aDIm.SetV(aNeigh,aNewMarq);
                  aVPts.push_back(aNeigh);
             }
          }
          aIndBottom++;
    }

    if (doRestore)
    {
       for (const auto & aPix : aVPts)
           aFRestore(aPix);
    }
}


void ConnectedComponent
     (
         std::vector<cPt2di> & aVPts,
         cDataIm2D<tU_INT1>  & aDIm,
         const std::vector<cPt2di> & aNeighbourhood,
         const  std::vector<cPt2di> & aVecSeed,
         int aMarqInit,
         int aNewMarq,
         bool doRestore
     )
{

    cDataIm2D<tU_INT1>  * aPtrDIm = & aDIm;


    Tpl_ConnectedComponent
    (
            aVPts,
            aNeighbourhood,
            aVecSeed,
            [aMarqInit,aPtrDIm](const cPt2di & aPix) {return aPtrDIm->GetV(aPix)== aMarqInit; },
            [aNewMarq,aPtrDIm](const cPt2di & aPix) {aPtrDIm->SetV(aPix,aNewMarq); },
            doRestore,
            [aMarqInit,aPtrDIm](const cPt2di & aPix) { aPtrDIm->SetV(aPix,aMarqInit); }
    );
}

void ConnectedComponent
     (
         std::vector<cPt2di> & aVPts,
         cDataIm2D<tU_INT1>  & aDIm,
         const std::vector<cPt2di> & aNeighbourhood,
         const cPt2di& aPtSeed,
         int aMarqInit,
         int aNewMarq,
         bool Restore
     )
{
    std::vector<cPt2di> aVSeed {aPtSeed};
    ConnectedComponent(aVPts,aDIm,aNeighbourhood,aVSeed,aMarqInit,aNewMarq,Restore);
}


void FlagBitConnectedComponent
     (
         std::vector<cPt2di> & aVPts,
         cDataIm2D<tU_INT4>  & aDIm,
         const std::vector<cPt2di> & aNeighbourhood,
         const  std::vector<cPt2di> & aVecSeed,
         int aBitIn,
         int aBitTmp,
         bool Restore
     )
{

    cDataIm2D<tU_INT4>  * aPtrDIm = & aDIm;


    Tpl_ConnectedComponent
    (
            aVPts,
            aNeighbourhood,
            aVecSeed,
            [aBitIn,aBitTmp,aPtrDIm](const cPt2di & aPix) {auto aV=aPtrDIm->GetV(aPix); return  BitIsOn(aV,aBitIn)&&(!BitIsOn(aV,aBitTmp)); },
            [aBitTmp,aPtrDIm](const cPt2di & aPix) {SetBitOn(*aPtrDIm,aPix,aBitTmp); },
            Restore,
            [aBitTmp,aPtrDIm](const cPt2di & aPix) {SetBitOff(*aPtrDIm,aPix,aBitTmp); }

    );
}

void BenchConnectedComponent(cParamExeBench & aParam)
{
    if (! aParam.NewBench("ConnectedComponent")) return;

    cPt2di aSz (100,120);

    cIm2D<tU_INT1>  anIm(aSz,nullptr,eModeInitImage::eMIA_Null);

    cIm2D<tU_INT4> aImBit(aSz,nullptr,eModeInitImage::eMIA_Null);


    // we fill image with regular squae on perido 10x12 with a little hole
    for (const auto & aPix : anIm.DIm())
    {
        int anX = aPix.x() % 10;
        int anY = aPix.y() % 12;

        bool InsideBig = (anX>=2) && (anX<8) && (anY>=2) && (anY<10);
        bool InsideSmall = (anX>=4) && (anX<6) && (anY>=4) && (anY<8);

        int aNum = 1+(aPix.x()/10) + (aPix.y()/12);

        if (InsideBig && (!InsideSmall))
        {
            anIm.DIm().SetV(aPix,aNum);
            aImBit.DIm().SetV(aPix,1<<(aNum+1));
        }
    }


    cIm2D<tU_INT1> anIm0 = anIm.Dup();
    cIm2D<tU_INT4> aImBit0 = aImBit.Dup();

    const auto & aV8= Alloc8Neighbourhood();
    int aNbCC = 0;
    for (const auto & aPix : anIm.DIm())
    {
        int aVal = anIm.DIm().GetV(aPix);
        if (aVal && (aVal<100))
        {
            std::vector<cPt2di> aCC;
            ConnectedComponent(aCC, anIm.DIm(), aV8,{aPix}, aVal, aVal+100 );

            std::vector<cPt2di> aCCBits;
            FlagBitConnectedComponent(aCCBits, aImBit.DIm(), aV8,{aPix}, 1+aVal, 0 ,true);

           //   StdOut() << "CBBBSz=" << aCCBits.size() << "\n";

            size_t aSzCC = (6*8) -(2*4); // size expected
            MMVII_INTERNAL_ASSERT_bench( aCC.size()==aSzCC,"Bad CC sz in bench");
            MMVII_INTERNAL_ASSERT_bench( aCC.size()==aCCBits.size(),"Bad CC sz in bench");

            aNbCC++;

            for (const auto & aPixCC : aCC)
                anIm0.DIm().AddVal(aPixCC,100);
        }
    }

 //   StdOut() << " DISTTT=" << anIm.DIm().LInfDist(anIm0.DIm())  << " " << aImBit.DIm().LInfDist(aImBit0.DIm())  << "\n";
    MMVII_INTERNAL_ASSERT_bench( aNbCC==100,"Bad CC number in bench");

    MMVII_INTERNAL_ASSERT_bench(anIm.DIm().LInfDist(anIm0.DIm())==0,"Bad Dist Check in bench-CC");
    MMVII_INTERNAL_ASSERT_bench(aImBit.DIm().LInfDist(aImBit0.DIm())==0,"Bad Dist Check in bench-CC");



    aParam.EndBench();
}


void MakeImageDist(cIm2D<tU_INT1> aImIn,const std::string & aNameChamfer)
{
     // Apparently Distance image are not used for now, will write it if needed
     MMVII_INTERNAL_ERROR("No MakeImageDist");
     // const Chamfer & aChamf = Chamfer::ChamferFromName(aNameChamfer);
     // aChamf.im_dist(cMMV1_Conv<tU_INT1>::ImToMMV1(aImIn.DIm())) ;
}


};
