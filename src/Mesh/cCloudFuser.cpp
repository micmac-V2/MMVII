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
                                const cPt2di & aSzOverlap,
                                bool ParalTiles);

        void APBI_ExecAll (bool Silence=false);

        cBox2di GlobBox2Parse(){return mGlobBoxToParse;}
        std::vector<cBox2di> SetOfBoxes(){return mSetBoxes;}

        void SetBox2Parse(cBox2di & aBox){mGlobBoxToParse=aBox;}
        void SetTheSetOfBoxes(std::vector<cBox2di> & aSet) {mSetBoxes=aSet;}
};


cAppliParsedBoxVirtualIm::cAppliParsedBoxVirtualIm(cMMVII_Appli & anAppli,
                                                    const cPt2di & aSzOverlap,
                                                    bool ParalTiles):
     cAppliParseBoxIm<tREAL4>(anAppli,eForceGray::Yes,cPt2di(2000,2000),aSzOverlap,ParalTiles),
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

    StdOut()<<mGlobBoxToParse<<" "<<mSzTiles<<std::endl;
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
    std::vector<cBox2dr> mSetRealBoxCalc;
    std::vector<cBox2di> mSetPixBoxCalc;
    tREAL8 mGlobGSD;
    std::string mTilingIndexFile;
    std::string mNameFusionResult;


    // here add private objets 
    public:

        cAppliCloudFuser(const std::vector<std::string> & aVArgs, const cSpecMMVII_Appli & aSpec );
        static constexpr tREAL8 mInfty =  -1e10;
        void ReadIndexTilingAndMetadataFile();
        void MergeDems(std::vector<cIm2D<tREAL4>> & aVDems,
                        std::vector<cIm2D<tU_INT1>> & aVMasq,
                        std::vector<cIm2D<tU_INT1>> & aVCorrel,
                        cIm2D<tREAL4> & aMergedDem,
                        cIm2D<tU_INT1> & aMergedMask,
                        cIm2D<tU_INT1> & aMergedCorrel);
};

cAppliCloudFuser::cAppliCloudFuser(const std::vector<std::string> & aVArgs, const cSpecMMVII_Appli & aSpec):
    cMMVII_Appli(aVArgs,aSpec), 
    cAppliParsedBoxVirtualIm(*this,cPt2di(50,50),true),
    mSetDemsPattern(""),
    mSetDemsNames({}),
    mSetDems({}),
    mGlobalPixBoxFusion(cBox2di::Empty()),
    mGlobalRealBoxFusion(cTplBoxOfPts<tREAL8,2>()),
    mGlobTf(cAffin2D<tREAL8>()),
    mSetRealBoxCalc({}),
    mSetPixBoxCalc({}),
    mGlobGSD(0.0),
    mTilingIndexFile(""),
    mNameFusionResult("Fusion")
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
            << AOpt2007(mNameFusionResult,"Out"," prefix of output files, default=Fusion.tif",{eTA2007::HDV})
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

void cAppliCloudFuser::MergeDems( std::vector<cIm2D<tREAL4>> & aVDems,
                                std::vector<cIm2D<tU_INT1>> & aVMasq,
                                std::vector<cIm2D<tU_INT1>> & aVCorrel,
                                cIm2D<tREAL4> & aMergedDem,
                                cIm2D<tU_INT1> & aMergedMask,
                                cIm2D<tU_INT1> & aMergedCorrel)
{

    // assume all image sizes are equal so that at one pixel location we have multiple values of z correl and Masq
    for (const auto & aPix: aMergedDem.DIm())
    {
        tREAL8 aVWeightedZ=0.0;
        tREAL8 aWeights =0.0;
        int aNb=0;
        for (size_t i=0; i<aVDems.size();i++)
        {
            if (aVMasq[i].DIm().GetV(aPix))
            {
                tREAL8 aWeight = aVCorrel[i].DIm().GetV(aPix)/255.0;
                aWeights += aWeight ;
                aVWeightedZ+=aWeight*aVDems[i].DIm().GetV(aPix);
                aNb++;
            }
        }
        // Weight dem values
        aMergedDem.DIm().SetV(aPix,(aWeights==0) ? 0.0 : aVWeightedZ/aWeights);
        aMergedCorrel.DIm().SetV(aPix,
                            std::min(round_ni(255.0*(aWeights/aVDems.size())),255));
        aMergedMask.DIm().SetV(aPix,aNb>0 ? 1 : 0);
    }
}


int cAppliCloudFuser::ExeOnParsedBox()
{
    StdOut()<<CurBoxIn()<<std::endl;
    // fast search of dems that overlap the concerned CurBoxIn()
    std::vector<cCloudRaster> mSetDemsOverlapCurBox;
    std::vector<cBox2di> mSetDemsOverlapBoxes;
    std::vector<cBox2di> mSetDemWhere2WriteBoxes;
    
    std::vector<cIm2D<tREAL4>> mSetOfBoxedDems;
    std::vector<cIm2D<tU_INT1>> mSetOfBoxedMasqs;
    std::vector<cIm2D<tU_INT1>> mSetOfBoxedCorrels;


    cPt2dr aP0Ul = mGlobTf.Value(ToR(CurBoxIn().P0()));
    cPt2dr aP1Lr = mGlobTf.Value(ToR(CurBoxIn().P1()));

    cBox2dr aBoxCurGlob = cBox2dr(cPt2dr(aP0Ul.x(),aP1Lr.y()),
                                    cPt2dr(aP1Lr.x(),aP0Ul.y()));


    StdOut()<<aBoxCurGlob<<std::endl;

    for (auto & aDem : mSetDems)
    {
        cAffin2D<tREAL8> aLocAff = aDem.Transform();

        cPt2dr aP0UlDem = aLocAff.Value(cPt2dr(0,0));
        cPt2dr aP1LrDem = aLocAff.Value(ToR(aDem.SzPix()));

        cBox2dr aBoxDem = cBox2dr(cPt2dr(aP0UlDem.x(),aP1LrDem.y()),
                                    cPt2dr(aP1LrDem.x(),aP0UlDem.y()));

        //StdOut()<<aBoxDem<<std::endl;

        cBox2dr aBoxInter = aBoxCurGlob.Inter(aBoxDem); 
        if (! aBoxInter.IsEmpty())
        {
            mSetDemsOverlapCurBox.push_back(aDem);

            cPt2di aP0UlLocInter = ToI(aLocAff.Inverse(aBoxInter.P0()));
            cPt2di aP1LrLocInter = ToI(aLocAff.Inverse(aBoxInter.P1()));

            cBox2di aBoxLocInDem = cBox2di(cPt2di(aP0UlLocInter.x(),aP1LrLocInter.y()),
                                            cPt2di(aP1LrLocInter.x(),aP0UlLocInter.y())); 

            cPt2di aP0UlGlobInter = ToI(mGlobTf.Inverse(aBoxInter.P0()));
            cPt2di aP1LrGlobInter = ToI(mGlobTf.Inverse(aBoxInter.P1()));

            cBox2di aBoxLocInCurBoxIn = cBox2di(cPt2di(aP0UlGlobInter.x(),aP1LrGlobInter.y()),
                                            cPt2di(aP1LrGlobInter.x(),aP0UlGlobInter.y()));  

            mSetDemsOverlapBoxes.push_back(aBoxLocInDem);
            mSetDemWhere2WriteBoxes.push_back(aBoxLocInCurBoxIn);
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
        //load relevant boxes of dems, correl/ambiguity, masqs for fusion
        cPt2di aP0InCurBox = aBoxInCurB.P0()-CurBoxIn().P0();
        cPt2di aP1InCurBox = aBoxInCurB.P1()-CurBoxIn().P0();

        StdOut()<<aBoxInCurB<<" "<<aP0InCurBox<<" "<<aP1InCurBox<<std::endl;
        
        cIm2D<tREAL4> aDemMap = cIm2D<tREAL4> (CurBoxIn().Sz());
        aDemMap.DIm().InitCste(mInfty);
        StdOut()<<aDem.NameDem()<<std::endl;
        aDemMap.Read(cDataFileIm2D::Create(aDem.NameDem(),eForceGray::No),
                    aBoxInDem.P0(),1.0,cPixBox<2>(aP0InCurBox,aP1InCurBox));
        

        StdOut()<<"aDemMap "<<std::endl;
        // mask
        cIm2D<tU_INT1> aDemMask = cIm2D<tU_INT1> (CurBoxIn().Sz());  
        aDemMask.DIm().InitCste(0);
        aDemMask.Read(cDataFileIm2D::Create(aDem.NameMasq(),eForceGray::No),
                    aBoxInDem.P0(),1.0,cPixBox<2>(aP0InCurBox,aP1InCurBox));

        StdOut()<<"aDemMapMasq "<<std::endl;
        cIm2D<tU_INT1> aDemCorrel = cIm2D<tU_INT1> (CurBoxIn().Sz());  
        aDemCorrel.DIm().InitCste(0);
        aDemCorrel.Read(cDataFileIm2D::Create(aDem.NameCorrel(),eForceGray::No),
                        aBoxInDem.P0(),1.0,cPixBox<2>(aP0InCurBox,aP1InCurBox));

        StdOut()<<"aDemMapCorrel "<<std::endl;
        mSetOfBoxedDems.push_back(aDemMap);
        mSetOfBoxedMasqs.push_back(aDemMask);
        mSetOfBoxedCorrels.push_back(aDemCorrel);
    }

    cIm2D<tREAL4> aFinalDem = cIm2D<tREAL4>(CurBoxIn().Sz());
    cIm2D<tU_INT1> aFinalMask = cIm2D<tU_INT1>(CurBoxIn().Sz());
    cIm2D<tU_INT1> aFinalCorrel = cIm2D<tU_INT1>(CurBoxIn().Sz());

    MergeDems(mSetOfBoxedDems,mSetOfBoxedMasqs,mSetOfBoxedCorrels,
              aFinalDem,aFinalMask, aFinalCorrel);

    // Save DEM, MASK and CORREL

    std::string aNameMergedDem =DirOfPath(mNameFusionResult,false)+"DEM_"+ 
                                FileOfPath(mNameFusionResult,false)+
                                ToStr(mIndBoxRecal.x())+"_"+ToStr(mIndBoxRecal.y())+".tif";

    cDataFileIm2D  aDFDem = cDataFileIm2D::Create(aNameMergedDem,
                                            eTyNums::eTN_REAL4,
                                            CurBoxOut().Sz(),
                                            1);

    aFinalDem.Write(aDFDem,CurP0(),1.0,CurBoxOut());

    //Masq
    std::string aNameMergedMask =DirOfPath(mNameFusionResult,false)+"MASQ_"+ 
                                FileOfPath(mNameFusionResult,false)+
                                ToStr(mIndBoxRecal.x())+"_"+ToStr(mIndBoxRecal.y())+".tif";


    cDataFileIm2D  aDFMasq = cDataFileIm2D::Create(aNameMergedMask,
                                            eTyNums::eTN_U_INT1,
                                            CurBoxOut().Sz(),
                                            1);

    aFinalMask.Write(aDFMasq,CurP0(),1.0,CurBoxOut());


    //Correl
    std::string aNameMergedCorrel =DirOfPath(mNameFusionResult,false)+"CORREL_"+ 
                                FileOfPath(mNameFusionResult,false)+
                                ToStr(mIndBoxRecal.x())+"_"+ToStr(mIndBoxRecal.y())+".tif";


    cDataFileIm2D  aDFCorrel = cDataFileIm2D::Create(aNameMergedCorrel,
                                            eTyNums::eTN_U_INT1,
                                            CurBoxOut().Sz(),
                                            1);

    aFinalCorrel.Write(aDFCorrel,CurP0(),1.0,CurBoxOut());


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
            mGlobalRealBoxFusion.Add(aTF.Value(cPt2dr(0,0)));
            mGlobalRealBoxFusion.Add(aTF.Value(ToR(aCldRst.SzPix())));
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
            mSetPixBoxCalc.push_back(
                                cBox2di(
                                    ToI(mGlobTf.Inverse(aRBox.P0())),
                                    ToI(mGlobTf.Inverse(aRBox.P1()))
                                ));
        }

        // some sorting and indexing with a mask to specify what should be compuated 
        std::sort(mSetPixBoxCalc.begin(), mSetPixBoxCalc.end(), [](const cBox2di& a, const cBox2di& b) 
        {
            if (a.P0().y() != b.P0().y()) return a.P0().y() < b.P0().y();
            return a.P0().x() < b.P0().x();                          
        });

        // run on a set of boxes in parallel
        SetTheSetOfBoxes(mSetPixBoxCalc);

        mGlobalPixBoxFusion = cPt2di(Pt_round_up(mGlobalRealBoxFusion.CurBox().Sz()/mGlobGSD));
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
        mGlobalPixBoxFusion = cPt2di(Pt_round_up(mGlobalRealBoxFusion.CurBox().Sz()/mGlobGSD));

        StdOut()<<mGlobalPixBoxFusion<<std::endl;
        SetBox2Parse(mGlobalPixBoxFusion);

        StdOut()<<GlobBox2Parse()<<std::endl;
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
        "Fuse multiple digital elevation maps",
        {eApF::Cloud},
        {eApDT::Image},
        {eApDT::Image},
        __FILE__
        );

};