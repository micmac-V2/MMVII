#include "cMMVII_Appli.h"
#include "MMVII_PCSens.h"
#include "MMVII_Geom2D.h"
#include "MMVII_Geom3D.h"
#include "MMVII_AllClassDeclare.h"
#include "MMVII_DeclareCste.h"
#include "MMVII_2Include_Serial_Tpl.h"
#include "MMVII_Triangles.h"
#include "MMVII_Image2D.h"
#include "MMVII_ZBuffer.h"
#include "MeshDev.h"
#include "MMVII_Sys.h"
#include "MMVII_Radiom.h"
#include "MMVII_CloudRaster.h"
#include "ogrsf_frmts.h"
#include "ranges"

namespace MMVII
{
static std::string NameIndBoxRecal="INTERNAL_IndexBoxRecall";

class cAppliParsedBoxVirtualIm : public cAppliParseBoxIm<tREAL4>
{
    private :
        cBox2di mGlobBoxToParse;
        std::vector<cBox2di> mSetBoxes;

    public:
        cAppliParsedBoxVirtualIm(cMMVII_Appli & anAppli,
                                const cPt2di & aSzTiles,
                                const cPt2di & aSzOverlap,
                                bool ParalTiles);

        void APBI_ExecAll (bool Silence=false);

        cBox2di GlobBox2Parse(){return mGlobBoxToParse;}
        std::vector<cBox2di> SetOfBoxes(){return mSetBoxes;}

        void SetBox2Parse(cBox2di & aBox){mGlobBoxToParse=aBox;}
        void SetTheSetOfBoxes(std::vector<cBox2di> & aSet) {mSetBoxes=aSet;}
};


cAppliParsedBoxVirtualIm::cAppliParsedBoxVirtualIm(cMMVII_Appli & anAppli,
                                                    const cPt2di & aSzTiles,
                                                    const cPt2di & aSzOverlap,
                                                    bool ParalTiles):
     cAppliParseBoxIm<tREAL4>(anAppli,eForceGray::No,aSzTiles,aSzOverlap,ParalTiles),
     mGlobBoxToParse(cBox2di::Empty()),
     mSetBoxes({})
{    
}

void cAppliParsedBoxVirtualIm::APBI_ExecAll(bool Silence)
{
    if (APBI_TestMode())
     {
        LoadI(CurBoxIn());
        mAppli.ExeOnParsedBox();
        return;
     }
    AssertNotInParsing();

    // int cParseBoxInOut
    cParseBoxInOut<2> aPBIO=cParseBoxInOut<2>::CreateFromSize(cBox2di(cPt2di(2000,2000)),
                                        mSzTiles);

    MMVII_INTERNAL_ASSERT_strong(!mGlobBoxToParse.IsEmpty(),
                                "Box to parse should not be empty for tiling");

    // update mSzTiles
    std::vector<cPt2di> aSetOfPixIndexes;
    bool aCustomSetOfBoxes=false;
    if(!mSetBoxes.empty())
    {
        mSzTiles = mSetBoxes[0].Sz();
        aCustomSetOfBoxes=true;
    }

    aPBIO  = cParseBoxInOut<2>::CreateFromSize(mGlobBoxToParse,mSzTiles);

    if(!mSetBoxes.empty())
    {
        for (const auto & aBox: mSetBoxes)
        {
            aSetOfPixIndexes.push_back(
                        CByC1P(aPBIO.BoxIndex().FromNormaliseCoord(aPBIO.BoxGlob().ToNormaliseCoord(aBox.Middle())),
                        round_ni)
                    );
        }
    }

    mParseBox = & aPBIO;
    std::list<cParamCallSys>  aLComParal;

    if(aCustomSetOfBoxes)
    {
        for (const auto & aPixI : aSetOfPixIndexes)
        {
            // if a the top level of paralelization, construct the string
            // For first box, run it classically so that files are created only once
            if (TopCallParallTile() && (aPixI!=aSetOfPixIndexes[0]))
            {
            //std::string aCom = mAppli.CommandOfMain() + " " +NameIndBoxRecal + "=" + ToStr(aPixI);
                cParamCallSys aCom = mAppli.CommandOfMain();
                aCom.AddArgs(NameIndBoxRecal + "=" + ToStr(aPixI));
                aLComParal.push_back(aCom);
            }
            else
            {
                // If not in paral do all box, else do only the box indicate by recall
                if ((!mParalTiles) || (aPixI==mIndBoxRecal))
                {
                    mCurPixIndex = aPixI;
                    //LoadI(CurBoxIn());
                    mAppli.ExeOnParsedBox();
                }
            }
        }
    }
    else
    {
        for (const auto & aPixI : aPBIO.BoxIndex())
        {
            // if a the top level of paralelization, construct the string
            // For first box, run it classically so that files are created only once
            if (TopCallParallTile() && (aPixI!=cPt2di(0,0)))
            {
            //std::string aCom = mAppli.CommandOfMain() + " " +NameIndBoxRecal + "=" + ToStr(aPixI);
                cParamCallSys aCom = mAppli.CommandOfMain();
                aCom.AddArgs(NameIndBoxRecal + "=" + ToStr(aPixI));
                aLComParal.push_back(aCom);
            }
            else
            {
                // If not in paral do all box, else do only the box indicate by recall
                if ((!mParalTiles) || (aPixI==mIndBoxRecal))
                {
                    mCurPixIndex = aPixI;
                    //LoadI(CurBoxIn());
                    mAppli.ExeOnParsedBox();
                }
            }
        }
    }

    mParseBox = nullptr ;   // No longer inside parsing
    mAppli.ExeComParal(aLComParal,Silence);

}


class cAppliCloudFuser : public cMMVII_Appli,
                         public cAppliParsedBoxVirtualIm

{
    private:
    int Exe() override;
    int ExeOnParsedBox() override;
    cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override ;
    cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override ;

    std::string mSetDemsPattern;
    std::vector<std::string> mSetDemsNames;
    std::vector<cCloudRaster> mSetDems;
    cBox2di mGlobalPixBoxFusion;
    cTplBoxOfPts<tREAL8,2> mGlobalRealBoxFusion;
    cAffin2D<tREAL8> mGlobTf;
    cAffin2D<tREAL8> mCurBoxInAffTf;
    std::vector<cBox2dr> mSetRealBoxCalc;
    std::vector<cBox2di> mSetPixBoxCalc;
    tREAL8 mGlobGSD;
    std::string mTilingIndexFile;
    bool IsOutputTiled;
    std::string mNamePrefixOut;
    std::string mNameFusionDemResult;
    std::string mNameFusionCorrelResult;
    std::string mNameFusionMasqResult;


    // here add private objets 
    public:

        cAppliCloudFuser(const std::vector<std::string> & aVArgs, const cSpecMMVII_Appli & aSpec );
        static constexpr tREAL8 mInfty =  -1e10;
        void ReadIndexTilingAndMetadataFile();
        void GenTFW(const cAffin2D<tREAL8> & anAff, const std::string & aNameTFW);
        void FuseDemsByWAvg( std::vector<cIm2D<tREAL4>> & aVDems,
                                 std::vector<cIm2D<tREAL4>> & aVWeighters,
                                 cIm2D<tREAL4> & aMergedDem,
                                 cIm2D<tU_INT1> & aMergedMask,
                                 cIm2D<tU_INT1> & aMergedCorrel);
        void FuseDemsByMedian( std::vector<cIm2D<tREAL4>> & aVDems,
                                 std::vector<cIm2D<tREAL4>> & aVWeighters,
                                 cIm2D<tREAL4> & aMergedDem,
                                 cIm2D<tU_INT1> & aMergedMask,
                                 cIm2D<tU_INT1> & aMergedCorrel);
        void ExtractInterpDems(cCloudRaster & aDem,
                                const cBox2dr & aBoxInDem, 
                                const cBox2dr & aBoxInCurTileCalc,
                                const cAffin2D<tREAL8> & anAffCurBox,
                                cIm2D<tREAL4> & aDemToExtract,
                                cIm2D<tREAL4> & aWeightingIm);
};

cAppliCloudFuser::cAppliCloudFuser(const std::vector<std::string> & aVArgs, const cSpecMMVII_Appli & aSpec):
    cMMVII_Appli(aVArgs,aSpec), 
    cAppliParsedBoxVirtualIm(*this,cPt2di(2000,2000),cPt2di(50,50),true),
    mSetDemsPattern(""),
    mSetDemsNames({}),
    mSetDems({}),
    mGlobalPixBoxFusion(cBox2di::Empty()),
    mGlobalRealBoxFusion(cTplBoxOfPts<tREAL8,2>()),
    mGlobTf(cAffin2D<tREAL8>()),
    mCurBoxInAffTf(cAffin2D<tREAL8>()),
    mSetRealBoxCalc({}),
    mSetPixBoxCalc({}),
    mGlobGSD(0.0),
    mTilingIndexFile(""),
    IsOutputTiled(false),
    mNamePrefixOut("Fusion"),
    mNameFusionDemResult(""),
    mNameFusionCorrelResult(""),
    mNameFusionMasqResult("")
{
}



cCollecSpecArg2007 & cAppliCloudFuser::ArgObl(cCollecSpecArg2007 & anArgObl)
{
    return
            //APBI_ArgObl(anArgObl)
            anArgObl
           <<   Arg2007(mSetDemsPattern,
                        "Set of xml files each defines a dem, an ambiuity/correl map and a Masq of definition",
                        {{eTA2007::MPatFile,"0"}})
        ;
}

cCollecSpecArg2007 & cAppliCloudFuser::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
    return 
            APBI_ArgOpt
        (
            anArgOpt
            <<AOpt2007(mTilingIndexFile,
                        "TilingIndex",
                        "Tiling Index file, i.e a GeoJson file that defines the tiling scheme.\n"
                        "If not given, a big tif is created as the fusion result.", 
                        {eTA2007::FileTagged})
            << AOpt2007(mNamePrefixOut,"Out"," prefix of output files, default=Fusion.tif",{eTA2007::HDV})
        )
        ;
}


void cAppliCloudFuser::ReadIndexTilingAndMetadataFile()
{
    GDALAllRegister();  // register all drivers (call once at startup)

    GDALDataset* poDS = (GDALDataset*) GDALOpenEx(
        mTilingIndexFile.c_str(),
        GDAL_OF_VECTOR | GDAL_OF_READONLY,
        nullptr,   // allowed drivers (nullptr = all)
        nullptr,   // open options
        nullptr    // sibling files
    );

    if (poDS == nullptr) {
        MMVII_INTERNAL_ASSERT_strong(false,
                    "Unable to read The tiling Index "+mTilingIndexFile);
    }

    // read features and capture geometry, i.e each tile index 

    // assume one layer

    OGRLayer* poLayer = poDS->GetLayer(0);
    poLayer->ResetReading();  

    OGRFeature* poFeature;

    while((poFeature=poLayer->GetNextFeature())!=nullptr)
    {
        OGRGeometry* poGeom = poFeature->GetGeometryRef();

        if(poGeom!=nullptr)
        {
            OGRwkbGeometryType eType = wkbFlatten(poGeom->getGeometryType());
            MMVII_INTERNAL_ASSERT_User((eType==wkbPolygon) || (eType==wkbMultiPolygon),
                                        eTyUEr::eBadEnum,
                                    "Bad Geometry type in GeoJson/shp file, It Should be Polygon or MultiPolygon");

            // now read 
           /*OGRPolygon* poPoly = poGeom->toPolygon();
            OGRLinearRing* poRing = poPoly->getExteriorRing();*/

            OGREnvelope aTileEnvelope;
            poGeom->getEnvelope(&aTileEnvelope);
            // get tile extent 
            cBox2dr aBoxTile(cPt2dr(aTileEnvelope.MinX,aTileEnvelope.MinY),
                            cPt2dr(aTileEnvelope.MaxX,aTileEnvelope.MaxY));

            // get Global real Box for fusion
            mGlobalRealBoxFusion.Add(aBoxTile.P0());
            mGlobalRealBoxFusion.Add(aBoxTile.P1());

            mSetRealBoxCalc.push_back(aBoxTile);
        }
        OGRFeature::DestroyFeature(poFeature);
    }

    GDALClose(poDS);

    //sort the set of pix boxes from low to high box coordinates
    //std::sort(mSetPixBoxCalc.begin(),mSetPixBoxCalc.end());
}

void cAppliCloudFuser::GenTFW(const cAffin2D<tREAL8> & anAff, const std::string & aNameTFW)
{
    std::ofstream aFtfw(aNameTFW.c_str());
    aFtfw.precision(10);

    aFtfw << anAff.VX().x() << "\n" << anAff.VX().y() << "\n";
    aFtfw << anAff.VY().x() << "\n" << anAff.VY().y() << "\n";
    aFtfw << anAff.Tr().x() << "\n" << anAff.Tr().y() << "\n";

    aFtfw.close();
}


void cAppliCloudFuser::FuseDemsByMedian( std::vector<cIm2D<tREAL4>> & aVDems,
                                 std::vector<cIm2D<tREAL4>> & aVWeighters,
                                 cIm2D<tREAL4> & aMergedDem,
                                 cIm2D<tU_INT1> & aMergedMask,
                                 cIm2D<tU_INT1> & aMergedCorrel)

{
    std::vector<std::pair<tREAL8,tREAL8>> aVZwithWeights;
    tREAL8 aTotWeights;
    int aNb;
    tU_INT1 aCorrel;
    for (const auto & aPix: aMergedDem.DIm())
    {
        aVZwithWeights.clear() ;
        aTotWeights= 0.0 ;
        aNb=0;
        aCorrel=0.0;

        for (size_t i=0; i<aVDems.size();i++)
        {
            tREAL8 aWeight = aVWeighters[i].DIm().GetV(aPix);
            if( aWeight > 0) // masq has a 1 value
            {
                aTotWeights += aWeight ;
                aVZwithWeights.push_back({aVDems[i].DIm().GetV(aPix),aWeight});
                aNb++;
            }
        }
        aMergedDem.DIm().SetV(aPix,
            (aTotWeights==0) ? 0.0 : NonConstWMediane(aVZwithWeights,aTotWeights));
        aCorrel = ( aNb==0 ) ? 0: std::min(round_ni(255.0*(aTotWeights/aNb)),255); 

        aMergedCorrel.DIm().SetV(aPix,aCorrel);
        aMergedMask.DIm().SetV(aPix,aNb>0 ? 1 : 0);
    }



}


void cAppliCloudFuser::FuseDemsByWAvg( std::vector<cIm2D<tREAL4>> & aVDems,
                                 std::vector<cIm2D<tREAL4>> & aVWeighters,
                                 cIm2D<tREAL4> & aMergedDem,
                                 cIm2D<tU_INT1> & aMergedMask,
                                 cIm2D<tU_INT1> & aMergedCorrel)
{
    tREAL8 aVWeightedZ=0.0;
    tREAL8 aWeights =0.0;
    tU_INT1 aCorrel =0 ;
    int aNb=0;
    // assume all image sizes are equal so that at one pixel location we have multiple values of z correl and Masq
    for (const auto & aPix: aMergedDem.DIm())
    {
        aVWeightedZ=0.0;
        aWeights =0.0;
        aCorrel =0 ;
        aNb=0;  

        for (size_t i=0; i<aVDems.size();i++)
        {
            tREAL8 aWeight = aVWeighters[i].DIm().GetV(aPix);
            if( aWeight >= 0) // masq has a 1 value
            {
                aWeights += aWeight ;
                aVWeightedZ+=aWeight*aVDems[i].DIm().GetV(aPix);
                aNb++;
            }
        }
        // Weight dem values
        aMergedDem.DIm().SetV(aPix,(aWeights==0) ? 0.0 : aVWeightedZ/aWeights);
        aCorrel = ( aNb==0 ) ? 0: std::min(round_ni(255.0*(aWeights/aNb)),255); 

        aMergedCorrel.DIm().SetV(aPix,aCorrel);
        aMergedMask.DIm().SetV(aPix,aNb>0 ? 1 : 0);
    }
}


void cAppliCloudFuser::ExtractInterpDems(cCloudRaster & aDem,
                                        const cBox2dr & aBoxInDem, 
                                        const cBox2dr & aBoxInCurTileCalc,
                                        const cAffin2D<tREAL8> & anAffCurBox,
                                        cIm2D<tREAL4> & aDemToExtract,
                                        cIm2D<tREAL4> & aWeightingIm)
{  
    //initialize values 
    aDemToExtract.DIm().InitCste(mInfty);
    aWeightingIm.DIm().InitCste(-1.0);

    cAffin2D<tREAL8> aDemAff= aDem.Transform();
    
    cDataFileIm2D aDemDFile =  cDataFileIm2D::Create(aDem.NameDem(),eForceGray::No) ;
    cDataFileIm2D aMasqDFile =  cDataFileIm2D::Create(aDem.NameMasq(),eForceGray::No) ;
    cDataFileIm2D aCorrelDFile =  cDataFileIm2D::Create(aDem.NameCorrel(),eForceGray::No) ;  
    
    cBox2di aBoxIntInDem = aBoxInDem.ToI();
    cBox2di aBoxIntInCurTileCalc = aBoxInCurTileCalc.ToI();

    //check boundaries 
    aBoxIntInDem = aBoxIntInDem.Inter(cBox2di(aDem.SzPix()));
    aBoxIntInCurTileCalc = aBoxIntInCurTileCalc.Inter(CurBoxInLoc());


    // make safe boundaries of aBoxIntInCurTileCalc 
    cBox2di aBoxDemTransform = anAffCurBox.InverseBox(
                                aDemAff.ValueBox(aBoxIntInDem.ToR())
                            ).ToI();

    aBoxIntInCurTileCalc = aBoxIntInCurTileCalc.Inter(aBoxDemTransform);

    // the way back 
    /*cBox2di aBoxCurTransform = aBoxedDemAffT.InverseBox(
                                anAffCurBox.ValueBox(aBoxIntInCurTileCalc.ToR())
                            ).ToI();

    aBoxIntInDem = aBoxIntInDem.Inter(aBoxCurTransform);*/
    // compute box in dem affine transform

    cAffin2D<tREAL8> aBoxedDemAffT = cAffin2D<tREAL8>(aDemAff.Value(ToR(aBoxIntInDem.P0())),
                                    aDemAff.VX(),
                                    aDemAff.VY());

    StdOut()<<"aBoxInInDem "<<aBoxIntInDem<<" BOX IM "<<
    " aBoxIntInCurTileCalc "<<aBoxIntInCurTileCalc<<"  "<<
    " aBoxIntInDemAfterTransform "<<aBoxDemTransform<<" "<<
    aDemToExtract.DIm().Sz()<<std::endl;


    // read tiled data and not the entire dem 
    cIm2D<tREAL4> aBoxedDem(aBoxIntInDem.Sz());
    aBoxedDem.Read(aDemDFile,aBoxIntInDem.P0(),1.0);

    cIm2D<tU_INT1> aBoxedCorrel(aBoxIntInDem.Sz());
    aBoxedCorrel.Read(aCorrelDFile,aBoxIntInDem.P0(),1.0);

    cIm2D<tU_INT1> aBoxedMasq(aBoxIntInDem.Sz());
    aBoxedMasq.Read(aMasqDFile,aBoxIntInDem.P0(),1.0);

    tREAL4 aZ,aWeighter;
    for (const auto & aPix: cPixBox<2>(aBoxIntInCurTileCalc.P0(),
                                        aBoxIntInCurTileCalc.P1()))
    {
        aWeighter = -1.0;
        // perform Transforms to get loc in boxed dems
        cPt2dr aPWorld  = anAffCurBox.Value(ToR(aPix));
        cPt2dr aPInBoxDem  = aBoxedDemAffT.Inverse(aPWorld);

        if(aBoxedMasq.DIm().Inside(ToI(aPInBoxDem)))
        {
            if( aBoxedMasq.DIm().GetV(ToI(aPInBoxDem)) )
            {
                if (aBoxedDem.DIm().InsideBL(aPInBoxDem))
                {
                    aWeighter = aBoxedCorrel.DIm().GetVBL(aPInBoxDem)/255.0;
                    aZ = aBoxedDem.DIm().GetVBL(aPInBoxDem);
                    // set 
                    aDemToExtract.DIm().SetV(aPix,aZ);
                    aWeightingIm.DIm().SetV(aPix,aWeighter);
                }

            }
        }
    }
}


int cAppliCloudFuser::ExeOnParsedBox()
{
    //Init Tiles affine transforms 

    // first get Tile transform

    // output region to write in cas we want the ouput to be tiled 

    cPt2dr anOffx= cPt2dr(mGlobGSD,0.0);
    cPt2dr anOffy=cPt2dr(0.0,-mGlobGSD);
    cPt2dr anOrigin = mGlobTf.Value(ToR(CurBoxOut().P0()));
    cAffin2D<tREAL8> aLocTileAffTransform = cAffin2D<tREAL8>(anOrigin,anOffx,anOffy);


    // the curbox affine transform
    mCurBoxInAffTf = cAffin2D<tREAL8>(mGlobTf.Value(ToR(CurBoxIn().P0())),
                                    anOffx,
                                    anOffy);
    

    // fast search of dems that overlap the concerned CurBoxIn()
    std::vector<cCloudRaster> mSetDemsOverlapCurBox;
    std::vector<cBox2dr> mSetDemsOverlapBoxes;
    std::vector<cBox2dr> mSetDemWhere2WriteBoxes;
    
    std::vector<cIm2D<tREAL4>> mSetOfBoxedDems;
    std::vector<cIm2D<tREAL4>> mSetOfBoxedWeighters;

    cBox2dr aBoxCurGlob = mGlobTf.ValueBox(CurBoxIn().ToR());
                            
    for (auto & aDem : mSetDems)
    {
        cAffin2D<tREAL8> aLocAff = aDem.Transform();

        cBox2dr aBoxDem = aLocAff.ValueBox(cBox2dr(ToR(aDem.SzPix())));

        //StdOut()<<aBoxDem<<std::endl;
        cBox2dr aBoxInter = aBoxCurGlob.Inter(aBoxDem); 
        if (! aBoxInter.IsEmpty())
        {
            mSetDemsOverlapCurBox.push_back(aDem);

            cBox2dr aBoxLocInDemR = aLocAff.InverseBox(aBoxInter);

            cBox2dr aBoxInCurBoxInR = mGlobTf.InverseBox(aBoxInter);
                                
            // aBox In CurBoxIn() coordinate system
            cBox2dr aBoxLocInCurBoxInR(aBoxInCurBoxInR.P0()-ToR(CurBoxIn().P0()),
                                        aBoxInCurBoxInR.P1()-ToR(CurBoxIn().P0()));

            mSetDemsOverlapBoxes.push_back(aBoxLocInDemR);
            mSetDemWhere2WriteBoxes.push_back(aBoxLocInCurBoxInR);
        }
    }

    // check overlapping dems
    if (!mSetDemsOverlapCurBox.size())
        return EXIT_SUCCESS; 
    
    auto anItDem = mSetDemsOverlapCurBox.begin();
    auto anItBoxIndem = mSetDemsOverlapBoxes.begin();
    auto anItBoxInCurBox = mSetDemWhere2WriteBoxes.begin();

    for (; anItDem != mSetDemsOverlapCurBox.end() && 
           anItBoxIndem != mSetDemsOverlapBoxes.end() && 
           anItBoxInCurBox != mSetDemWhere2WriteBoxes.end();
            ++anItDem, ++anItBoxIndem,++anItBoxInCurBox)
    {
        auto [aDem, aBoxInDem,aBoxInCurB] = std::tie(*anItDem, *anItBoxIndem,*anItBoxInCurBox);
        cIm2D<tREAL4> aDemMap(CurBoxIn().Sz());
        cIm2D<tREAL4> aWeighterMap(CurBoxIn().Sz());

        // read dems and construct weighters 
        ExtractInterpDems(aDem,aBoxInDem,aBoxInCurB,mCurBoxInAffTf,
                            aDemMap,aWeighterMap);

        mSetOfBoxedDems.push_back(aDemMap);
        mSetOfBoxedWeighters.push_back(aWeighterMap);
    }

    cIm2D<tREAL4> aFinalDem = cIm2D<tREAL4>(CurBoxIn().Sz());
    cIm2D<tU_INT1> aFinalMask = cIm2D<tU_INT1>(CurBoxIn().Sz());
    cIm2D<tU_INT1> aFinalCorrel = cIm2D<tU_INT1>(CurBoxIn().Sz());

    FuseDemsByMedian(mSetOfBoxedDems,mSetOfBoxedWeighters,
              aFinalDem,aFinalMask, aFinalCorrel);
    // Save DEM, MASK and CORREL

    if(IsOutputTiled)
    {
        std::string aNameMergedDem =DirOfPath(mNamePrefixOut,false)+"DEM_"+ 
                                    FileOfPath(mNamePrefixOut,false)+
                                    ToStr(mIndBoxRecal.x())+"_"+ToStr(mIndBoxRecal.y());

        
        GenTFW(aLocTileAffTransform,aNameMergedDem +".tfw");                            
        
        cDataFileIm2D  aDFDem = cDataFileIm2D::Create(aNameMergedDem+".tif",
                                                eTyNums::eTN_REAL4,
                                                CurBoxOut().Sz(),
                                                1);
                                
        aFinalDem.Write(aDFDem,-CurBoxOutLoc().P0(),1.0,CurBoxOutLoc());

        //Masq
        std::string aNameMergedMask =DirOfPath(mNamePrefixOut,false)+"MASQ_"+ 
                                    FileOfPath(mNamePrefixOut,false)+
                                    ToStr(mIndBoxRecal.x())+"_"+ToStr(mIndBoxRecal.y());

        GenTFW(aLocTileAffTransform,aNameMergedMask + ".tfw");     


        cDataFileIm2D  aDFMasq = cDataFileIm2D::Create(aNameMergedMask+".tif",
                                                eTyNums::eTN_U_INT1,
                                                CurBoxOut().Sz(),
                                                1);

        aFinalMask.Write(aDFMasq,-CurBoxOutLoc().P0(),1.0,CurBoxOutLoc());


        //Correl
        std::string aNameMergedCorrel =DirOfPath(mNamePrefixOut,false)+"CORREL_"+ 
                                    FileOfPath(mNamePrefixOut,false)+
                                    ToStr(mIndBoxRecal.x())+"_"+ToStr(mIndBoxRecal.y());

        GenTFW(aLocTileAffTransform,aNameMergedCorrel+".tfw");    

        cDataFileIm2D  aDFCorrel = cDataFileIm2D::Create(aNameMergedCorrel+".tif",
                                                eTyNums::eTN_U_INT1,
                                                CurBoxOut().Sz(),
                                                1);

        aFinalCorrel.Write(aDFCorrel,-CurBoxOutLoc().P0(),1.0,CurBoxOutLoc());
    }
    else
    {

        // create if needed 
        cDataFileIm2D  aDFH = cDataFileIm2D::Create(mNameFusionDemResult+".tif",
                                              eTyNums::eTN_REAL4,
                                              GlobBox2Parse().Sz(),
                                              1);
        cDataFileIm2D  aDFC = cDataFileIm2D::Create(mNameFusionCorrelResult+".tif",
                                              eTyNums::eTN_U_INT1,
                                              GlobBox2Parse().Sz(),
                                              1);
        cDataFileIm2D  aDFM = cDataFileIm2D::Create(mNameFusionMasqResult+".tif",
                                              eTyNums::eTN_U_INT1,
                                              GlobBox2Parse().Sz(),
                                              1); 

        aFinalDem.Write(aDFH,CurP0(),1.0,CurBoxOutLoc());
        aFinalCorrel.Write(aDFC,CurP0(),1.0,CurBoxOutLoc());
        aFinalMask.Write(aDFM,CurP0(),1.0,CurBoxOutLoc());
    }

    return EXIT_SUCCESS;
}


int cAppliCloudFuser::Exe()
{
    // may be mPhProj.FinishInit();
   mSetDemsNames = VectMainSet(0);

   
   int aNB=0;
   // read serialized cloudRaster objects from bascule to get global context
    for(const auto & aCldName: mSetDemsNames)
    {
        cCloudRaster aCldRst;
        ReadFromFile_Std(aCldRst,aCldName);
        mSetDems.push_back(aCldRst);

        // increment global context
        cAffin2D<tREAL8> aTF = aCldRst.Transform();


        // compute global context
        if(!IsInit(&mTilingIndexFile))
        {
            mGlobalRealBoxFusion.Add(aTF.Value(cPt2dr(aCldRst.SzPix().x(),0)));
            mGlobalRealBoxFusion.Add(aTF.Value(cPt2dr(0,aCldRst.SzPix().y())));
        }
        // gsd
        mGlobGSD+=aTF.VX().x();
        aNB++;
    }
    // gsd of output fused dem
    mGlobGSD/=aNB;

    //  Check if Tiling Index was initialized or not 

    if(IsInit(&mTilingIndexFile))
    {
        IsOutputTiled =true;
        // Computation parallelization will be directly be computed from the tiling index itself
        // i.e Box per computation

        bool CheckExtensionReadable = (ends_with(mTilingIndexFile,"geojson") ||
                                        ends_with(mTilingIndexFile,"shp") );

        MMVII_INTERNAL_ASSERT_User(CheckExtensionReadable,
                                eTyUEr::eBadPostfix,
                                "Index File is not a common geojson or shp file");

        // read mTilingIndexFile
        ReadIndexTilingAndMetadataFile();

        // global affine transform for the resulting output image
        mGlobTf = cAffin2D<tREAL8>(cPt2dr(mGlobalRealBoxFusion.CurBox().P0().x(),
                                    mGlobalRealBoxFusion.CurBox().P1().y()),
                                    cPt2dr(mGlobGSD,0),
                                    cPt2dr(0,-mGlobGSD));

        
        for (const auto & aRBox: mSetRealBoxCalc)
        {
            mSetPixBoxCalc.push_back(mGlobTf.InverseBox(aRBox).ToI());
        }
    

        // some sorting and indexing with a mask to specify what should be compuated 
        std::sort(mSetPixBoxCalc.begin(), mSetPixBoxCalc.end(), [](const cBox2di& a, const cBox2di& b) 
        {
            if (a.P0().y() != b.P0().y()) return a.P0().y() < b.P0().y();
            return a.P0().x() < b.P0().x();                          
        });

        // run on a set of boxes in parallel
        SetTheSetOfBoxes(mSetPixBoxCalc);

        mGlobalPixBoxFusion = ToI(mGlobalRealBoxFusion.CurBox().Sz()/mGlobGSD);
        SetBox2Parse(mGlobalPixBoxFusion);
        
    }
    else
    {
        // global affine transform for the resulting output image
        mGlobTf = cAffin2D<tREAL8>(cPt2dr(mGlobalRealBoxFusion.CurBox().P0().x(),
                                        mGlobalRealBoxFusion.CurBox().P1().y()),
                                        cPt2dr(mGlobGSD,0),
                                        cPt2dr(0,-mGlobGSD));
        // we will follow the standard MMVII tiling scheme in cAppliParseBoxImVirtual
        mGlobalPixBoxFusion = cBox2di(Pt_round_up(mGlobalRealBoxFusion.CurBox().Sz()/mGlobGSD));

        SetBox2Parse(mGlobalPixBoxFusion);
        
        //init output file names
        mNameFusionDemResult = DirOfPath(mNamePrefixOut,false)+"DEM_"+ 
                                    FileOfPath(mNamePrefixOut,false);
        mNameFusionCorrelResult =DirOfPath(mNamePrefixOut,false)+"CORREL_"+ 
                                    FileOfPath(mNamePrefixOut,false);
        mNameFusionMasqResult =DirOfPath(mNamePrefixOut,false)+"MASQ_"+ 
                                    FileOfPath(mNamePrefixOut,false);

        // TFW FILES
        GenTFW(mGlobTf,mNameFusionDemResult+".tfw");
        GenTFW(mGlobTf,mNameFusionCorrelResult+".tfw");   
        GenTFW(mGlobTf,mNameFusionMasqResult+".tfw");   
        
    }

    APBI_ExecAll();

    return EXIT_SUCCESS;
}


/*  ============================================= */
/*       ALLOCATION                               */
/*  ============================================= */

tMMVII_UnikPApli Alloc_CloudFuser(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli & aSpec)
{
    return tMMVII_UnikPApli(new cAppliCloudFuser(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpecCloudFuser
    (
        "CloudMMVII_Fuse",
        Alloc_CloudFuser,
        "Fuse multiple digital elevation maps given ambiguity/correlation maps",
        {eApF::Cloud},
        {eApDT::Image},
        {eApDT::Image},
        __FILE__
        );

};