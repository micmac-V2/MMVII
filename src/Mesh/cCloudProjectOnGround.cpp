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
#include "MMVII_Sys.h"
#include "MMVII_Radiom.h"
#include "MMVII_CloudRaster.h"
#include <fstream>
#include <iostream>



namespace MMVII
{

 class cAppliCloudProjectOnGround : public cMMVII_Appli,
                            public cAppliParseBoxIm<tREAL4>
{
    private:
        int Exe() override;
        int ExeOnParsedBox() override;


        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override ;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override ;

        std::string mNameDMEpipPairResult;
        std::string mNameCloud2D_DepthIn;
        std::string mNameMasq1;
        std::string mNameCorrel;
        std::string mNameSec;

        cSensorImage * mCamGen;
        cSensorImage * mSecCamGen;
        cIm2D<tREAL4> mImPx1;
        cIm2D<tU_INT1> mImMasq1;
        cIm2D<tU_INT1> mImCorrel;
        cIm2D<tU_INT1> mImMasqOut;
        cIm2D<tU_INT1> mImCorrelOut;
        cIm2D<tREAL4> mImRed;

        std::string mNameBascOut;
        std::string mNameMasqOut;
        std::string mNameCorrelOut;

        // photogrammetric project
        cPhotogrammetricProject mPhProj;
        cTriangulation3D<tREAL8> * mTri3D;
        cTriangulation3D<tREAL8> * mTri2DDepth;
        std::string mNameResult;
        eModeGeom mModeGeom;
        tREAL8 mGSD;
        cBox2di mBoxGlobOutPix;
        cAffin2D<tREAL8> mGlobAff;
        std::string mOutDir;
        tREAL8 mNoiseZ;
        tREAL8 mThreshGrad;
        bool mBascCorrel;
        bool mOriInFatherDir;  ///< if true, search orientation in father dir, else in DirProj
        double      mMII;   ///<  Marge Inside Image
        double mStretschingThresh;
        //cBox2dr mBoxLocTarget;
        cTplBoxOfPts<tREAL8,2> mBoxLocTarget;
        cTplBoxOfPts<tREAL8,2> mBoxGlobTarget;
        std::vector<tU_INT1> mVCorrel;
        std::vector <tU_INT1> mVMasq;

    public:
        cAppliCloudProjectOnGround(const std::vector<std::string> & aVArgs, const cSpecMMVII_Appli & aSpec );
        static constexpr tREAL8 mInfty =  -1e10;
        std::pair<cPt3dr,cPt3dr> BascOnePoint(const cPt2di A, const cPt2di anOffSet, bool & oValid);
        /// Estimate the Ground Sampling Distance
        tREAL8 ComputeGSD();
        void MakeFastBasc();
        void MakeBasculeTris(cZBuffer & aZB);
        void ProcessNoPix(cZBuffer &  aZB);
        void GenTFW(const cAffin2D<tREAL8> & anAff, const std::string & aNameTFW);
        static cAffin2D<tREAL8> ReadTFW(const std::string & aNameTFW);
        cBox2di  BoxUtile(cIm2D<tU_INT1> & anImMasq);
        std::vector<cOneHelpSampleCmp> Samples() const;
        void MergeResults();
};


cAppliCloudProjectOnGround::cAppliCloudProjectOnGround(const std::vector<std::string> & aVArgs, const cSpecMMVII_Appli & aSpec ):
    cMMVII_Appli(aVArgs,aSpec),
    cAppliParseBoxIm<tREAL4>(*this,eForceGray::No,cPt2di(2000,2000),cPt2di(50,50),true),
    mCamGen(nullptr),
    mSecCamGen(nullptr),
    mImPx1(cPt2di(1,1)),
    mImMasq1(cPt2di(1,1)),
    mImCorrel(cPt2di(1,1)),
    mImMasqOut(cPt2di(1,1)),
    mImCorrelOut(cPt2di(1,1)),
    mImRed(cPt2di(1,1)),
    mNameBascOut(""),
    mNameMasqOut(""),
    mNameCorrelOut(""),
    mPhProj(*this),
    mTri3D(nullptr),
    mTri2DDepth(nullptr),
    mNameResult("Dem_"),
    mModeGeom(eModeGeom::eGEOM_EPIP),
    mGSD(0.2),
    mBoxGlobOutPix(cBox2di::Empty()),
    mGlobAff(),
    mNoiseZ(0.2),
    mThreshGrad(0.3),
    mBascCorrel(false),
    mOriInFatherDir(false),
    mMII(0.0),
    mStretschingThresh(4.0)
{
}


cCollecSpecArg2007 & cAppliCloudProjectOnGround::ArgObl(cCollecSpecArg2007 & anArgObl)
{
    // use xml file serialized by cDMEpipPairResult class in MMVII_CloudRaster.h
    return
        anArgObl
           <<   Arg2007(mNameDMEpipPairResult,"xml file containing the result of epipolar resampling",{eTA2007::FileTagged})
           <<   mPhProj.DPOrient().ArgDirInMand()
           <<   mPhProj.DPMeshDev().ArgDirOutMand()
        ;
}


cCollecSpecArg2007 & cAppliCloudProjectOnGround::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
    return
        APBI_ArgOpt
        (
            anArgOpt
                << AOpt2007(mNameResult,"NameResult","Prefix for output files")
                << AOpt2007(mGSD,"GSD", "Ground sampling distance for rasterization")
                << AOpt2007(mOutDir,"DirOut","Directory of output files, (def=VISU/"+Specs().Name()+")")
                << AOpt2007(mMII,"MII","Margin Inside Image (for triangle validation)", {eTA2007::HDV})
                << AOpt2007(mOriInFatherDir,"OriInFatherDir","If true, search orientation in father dir, else in DirProj",{eTA2007::HDV})
                << AOpt2007(mStretschingThresh,"DistortionThreshold","Level of triangle distortion to discard from bascule")
        )
        ;
}

cBox2di cAppliCloudProjectOnGround::BoxUtile(cIm2D<tU_INT1> & anImMasq)
{
    cTplBoxOfPts<int,2> aBox;
    for (const auto & aPix: anImMasq.DIm())
    {
        if(anImMasq.DIm().GetV(aPix)>0.99)
            aBox.Add(aPix);
    }

    if (!aBox.NbPts()) return cBox2di::Empty();

    return aBox.CurBox();
}


std::pair<cPt3dr,cPt3dr> cAppliCloudProjectOnGround::BascOnePoint(const cPt2di A, const cPt2di anOffSet, bool & oValid)
{
    cPt3dr aPZA, aPWithDepth;

    cPt2di  aPix1 = A+anOffSet;
    oValid = true;

    if(mModeGeom==eModeGeom::eGEOM_EPIP)  ///< entry is a disparity maps so need to compute pseudo intersection
    {
        double  aX2= aPix1.x() + mImPx1.DIm().GetV(A);
        cPt2dr aPix2=cPt2dr(aX2, aPix1.y());
        tSeg3dr aP2TerCam1 = mCamGen->Image2Bundle(ToR(aPix1));
        tSeg3dr aP2TerCam2= mSecCamGen->Image2Bundle(aPix2);

        aPZA=RobustBundleInters({aP2TerCam1,aP2TerCam2});

        tREAL4 aHgtValue = - aPZA.z();
        
        
        if (dynamic_cast<cSensorCamPC*>(mCamGen))
            aHgtValue = mCamGen->Ground2ImageAndDepth(aPZA).z();


        aPWithDepth= cPt3dr(aPix1.x(),
                             aPix1.y(),
                             aHgtValue);
    }
    else if (mModeGeom==eModeGeom::eGEOM_DEPTH) ///< entry is a depth map so just compute 3D world coordinates
    {
        tREAL8 aDepth = mImPx1.DIm().GetV(A);
        if (aDepth==0)  ///< 0 = no-data in the depth map : skip this point
        {
            oValid = false;
            return {cPt3dr(aPix1.x(),aPix1.y(),0.0),cPt3dr(0,0,0)};
        }

        aPWithDepth= cPt3dr(aPix1.x(),
                            aPix1.y(),
                            aDepth);

        MMVII_INTERNAL_ASSERT_tiny(dynamic_cast<cSensorCamPC*>(mCamGen),
            "Camera is not a cSensorCamPC in BascOnePoint Mode GEOM_DEPTH");

        aPZA = mCamGen->ImageAndDepth2Ground(aPWithDepth);
    }
    else
    {
        ///< nothing for now
    }
    mBoxLocTarget.Add(Proj(aPZA));

    return {aPWithDepth,aPZA};
}


tREAL8 cAppliCloudProjectOnGround::ComputeGSD()
{
    // Estimate the Ground Sampling Distance
    if( mModeGeom != eModeGeom::eGEOM_EPIP && mModeGeom != eModeGeom::eGEOM_DEPTH )
    {
        MMVII_INTERNAL_ASSERT_tiny(false,"ComputeGSD : Unsupported eModeGeom mode");
    }
    // read mImPx1 not necessarily the whole image, just a box as it is enough to estimate GSD
    const cPt2di aOffsetRead = cBox2di(mCamGen->Sz()).Middle()-cPt2di(50,50);
    mImPx1 = cIm2D<tREAL4>(cPt2di(100,100));
    mImPx1.Read(cDataFileIm2D::Create(DirProject()+mNameCloud2D_DepthIn,mIsGray),aOffsetRead);

    cIm2D<tU_INT1> aMaskMap = cIm2D<tU_INT1>(cPt2di(100,100));
    aMaskMap.Read(cDataFileIm2D::Create(DirProject()+mNameMasq1,mIsGray),
                  aOffsetRead);

    bool IsCamLoaded = ( mCamGen != nullptr );
    if (mModeGeom == eModeGeom::eGEOM_EPIP)
    {
        IsCamLoaded = IsCamLoaded && (mSecCamGen != nullptr);
    }

    MMVII_INTERNAL_ASSERT_tiny(IsCamLoaded,
        "ComputeGSD : camera(s) not loaded before computing GSD");

    std::vector<tREAL8> aVGSD;
    tREAL8 aDepth;
    cPt3dr aPGround;
    for (const auto aPix: aMaskMap.DIm())
    {
        if (aMaskMap.DIm().GetV(aPix) < 0.99)
            continue;

        bool aValid = true;
        auto [aP00_Depth,aP00_3D] = BascOnePoint(aPix,aOffsetRead,aValid);

        if (! aValid)
            continue;

        cWeightAv<tREAL8> aWeighD;
        aDepth  = aP00_Depth.z();
        for (const auto & aNeigh : AllocNeighbourhood<2>(1))
        {
            cPt2dr aPIm = Proj(aP00_Depth) + ToR(aNeigh);
            if(dynamic_cast<cSensorCamPC*>(mCamGen))
                aPGround = mCamGen->ImageAndDepth2Ground(TP3z(aPIm,aDepth));
            else
                aPGround = mCamGen->ImageAndZ2Ground(TP3z(aPIm,aP00_3D.z()));

            tREAL8 aD = Norm2(aPGround-aP00_3D);

            aWeighD.Add(1.0,aD);
        }
        aVGSD.push_back(aWeighD.Average());
    }

    MMVII_INTERNAL_ASSERT_tiny(aVGSD.size()>3,"No sufficient GSD values to compute median");

    // back to init of mImPx1
    mImPx1 = cIm2D<tREAL4>(cPt2di(1,1));

    // keep only three decimal places to avoid over-optimistic precision
    const tREAL8 aMedGSD = ConstMediane(aVGSD);
    const tREAL8 aTruncGSD = std::trunc(aMedGSD * 1000.0) / 1000.0;

    // a very fine GSD (< 0.001) would truncate to 0, which would later cause a
    // division by zero when sizing the output raster (Sz/mGSD) : keep the raw median
    return (aTruncGSD > 0.0) ? aTruncGSD : aMedGSD;
}


void cAppliCloudProjectOnGround::MakeFastBasc()
{
    //mImRed.DIm().Resize(CurSzIn(),eModeInitImage::eMIA_Null);
    cPt2di aP0 = CurP0();
    cPt2di aP1= CurBoxIn().P1();
    std::vector<cPt3dr> aVPts;
    std::vector<cPt3dr> aVPtsDepth;
    std::vector<cPt3di> aVFaces;

    StdOut() << "MakeFastBasc " << aP0 << " " << aP1 << "\n";
    cBox2di  aBoxUtileWithMasq = BoxUtile(mImMasq1);
    if (aBoxUtileWithMasq.IsEmpty())
        return;

    StdOut() << "BoxUtileWithMasq " << aBoxUtileWithMasq << "\n";

    ///  Init index for faces
    cIm2D<int> anImIndex(aBoxUtileWithMasq.Sz());
    cDataIm2D<int> & aDIdx= anImIndex.DIm();
    aDIdx.InitCste(-1);   ///<  -1 => no valid vertex for this pixel (e.g. depth==0)

    cAutoTimerSegm aTSMeshCreate(TimeSegm(),"ZBuffer::PrepareProjections"); 
    /// Create mesh points, keeping only those with a valid depth
    int IndBegin=0;
    for(const auto & aPix: cPixBox<2>(aBoxUtileWithMasq.P0(),
                                       aBoxUtileWithMasq.P1())
         )
    {
        bool aValid = true;
        auto [aP00_Depth,aP00_3D] = BascOnePoint(aPix,aP0,aValid);

        if (! aValid)  ///< no-data depth : do not create a vertex here
            continue;

        aDIdx.SetV(aPix-aBoxUtileWithMasq.P0(),IndBegin);
        IndBegin++;

        aVPtsDepth.push_back(aP00_Depth);
        aVPts.push_back(aP00_3D);

        if(mBascCorrel)
            mVCorrel.push_back(mImCorrel.DIm().GetV(aPix));
    }

    /// Create mesh faces (only when the 4 corners have a valid vertex)
    for(const auto & aPix: cPixBox<2>(aBoxUtileWithMasq.P0(),
                                       aBoxUtileWithMasq.P1()-cPt2di(1,1))
         )
    {
        cPt2di P00(aPix.x(),aPix.y());
        cPt2di P10(aPix.x()+1,aPix.y());
        cPt2di P01(aPix.x(),aPix.y()+1);
        cPt2di P11(aPix.x()+1,aPix.y()+1);

        int aI00 = aDIdx.GetV(P00-aBoxUtileWithMasq.P0());
        int aI10 = aDIdx.GetV(P10-aBoxUtileWithMasq.P0());
        int aI01 = aDIdx.GetV(P01-aBoxUtileWithMasq.P0());
        int aI11 = aDIdx.GetV(P11-aBoxUtileWithMasq.P0());

        if ((aI00<0)||(aI10<0)||(aI01<0)||(aI11<0))  ///< a corner is no-data => skip face
            continue;

        aVFaces.push_back(cPt3di(aI00,aI10,aI11));
        aVFaces.push_back(cPt3di(aI00,aI11,aI01));
    }

    cAutoTimerSegm aTSBufferInit(TimeSegm(),"ZBuffer::Init"); 

    /// ZBUFFER
    mTri3D = new cTriangulation3D<tREAL8>(aVPts,aVFaces);
    mTri2DDepth = new cTriangulation3D<tREAL8> (aVPtsDepth,aVFaces);

    cMeshTri3DIterator  aTriIt(mTri3D);

    cMeshTri3DIterator * aTriIT2DDepth= new cMeshTri3DIterator(mTri2DDepth);

    // Keep a real polymorphic mapping object (held by pointer) : cZBuffer stores a reference to
    // it (const tMap & mMapI2O), so it must outlive the ZBuffer and must NOT be sliced into a
    // base tIMap_R3 value (slicing would drop the Value()/Inverse() overrides and trigger the
    // default Value<->Values mutual recursion -> "Forbiden Recursive Call").
    std::unique_ptr<tIMap_R3> aMapCamDepth;

    if (dynamic_cast<cSensorCamPC*>(mCamGen))
        aMapCamDepth = std::make_unique<cSIMap_Ground2ImageAndProf>(mCamGen);
    else
        aMapCamDepth = std::make_unique<cSIMap_Ground2ImageAndZ>(mCamGen);


    cSetVisibility aSetVis(mCamGen,mMII);

    double Infty =1e20;

    cBox3dr  aBox(cPt3dr(aP0.x(),aP0.y(),-Infty),cPt3dr(aP1.x(),aP1.y(),Infty));

    cDataBoundedSet<tREAL8,3>  aSetCam(aBox);

    StdOut()<<"ZBUFFER::ZBUFFER"<<std::endl;
    cZBuffer aZBuf(aTriIt,
                   aSetVis,
                   *aMapCamDepth,
                   aSetCam,
                   1,
                   true,
                   true,
                   aTriIT2DDepth);


    cAutoTimerSegm aTSBufferProjInit(TimeSegm(),"ZBuffer::ProjInit"); 

    aZBuf.MakeZBufForBasc(eZBufModeIter::ProjInit);


    

    cAutoTimerSegm aTSBufferSurfDev(TimeSegm(),"ZBuffer::SurfDevlpt"); 

    aZBuf.MakeZBufForBasc(eZBufModeIter::SurfDevlpt);

    cAutoTimerSegm aTSBufferProcNoPix(TimeSegm(),"ZBuffer::ProcessNoPix"); 

    ProcessNoPix(aZBuf);
    
    cAutoTimerSegm aTSBufferPixellizer(TimeSegm(),"ZBuffer::Pixellization"); 
    MakeBasculeTris(aZBuf);

    delete mTri2DDepth;
    delete mTri3D;
}

void cAppliCloudProjectOnGround::GenTFW(const cAffin2D<tREAL8> & anAff, const std::string & aNameTFW)
{
    std::ofstream aFtfw(aNameTFW.c_str());
    aFtfw.precision(10);

    aFtfw << anAff.VX().x() << "\n" << anAff.VX().y() << "\n";
    aFtfw << anAff.VY().x() << "\n" << anAff.VY().y() << "\n";
    aFtfw << anAff.Tr().x() << "\n" << anAff.Tr().y() << "\n";

    aFtfw.close();
}

cAffin2D<tREAL8> cAppliCloudProjectOnGround::ReadTFW(const std::string & aNameTFW)
{
    std::ifstream aFtfw(aNameTFW.c_str());
    MMVII_INTERNAL_ASSERT_User(aFtfw.good(),eTyUEr::eOpenFile,
        "ReadTFW : cannot open TFW file " + aNameTFW);

    std::string aline;
    std::vector< std::string > acontent;
    while(std::getline(aFtfw,aline))
    {
        acontent.push_back(aline);
    }
    aFtfw.close();

    MMVII_INTERNAL_ASSERT_User(acontent.size() >= 6,eTyUEr::eUnClassedError,
        "ReadTFW : expected at least 6 lines in TFW file " + aNameTFW);

    // TFW layout matches GenTFW : VX.x VX.y VY.x VY.y Tr.x Tr.y
    tREAL8 aVXx = stod(acontent.at(0));
    tREAL8 aVXy = stod(acontent.at(1));
    tREAL8 aVYx = stod(acontent.at(2));
    tREAL8 aVYy = stod(acontent.at(3));
    tREAL8 aTrx = stod(acontent.at(4));
    tREAL8 aTry = stod(acontent.at(5));

    return cAffin2D<tREAL8>(cPt2dr(aTrx,aTry),cPt2dr(aVXx,aVXy),cPt2dr(aVYx,aVYy));
}

void cAppliCloudProjectOnGround::ProcessNoPix(cZBuffer &  aZB)
{
    // comppute dual graph to have neigbouring relation between faces
    mTri3D->MakeTopo();
    const cGraphDual &  aGrD = mTri3D->DualGr() ;

    bool  GoOn = true;
    while (GoOn)  // continue as long as we get some modification
    {
        GoOn = false;
        // parse all face
        for (size_t aKF1 = 0 ; aKF1<mTri3D->NbFace() ; aKF1++)
        {
            // check for each face labled  "NoPix" if it  has a neigboor visible (or likely)
            if (aZB.ResSurfD(aKF1).mResult == eZBufRes::NoPix)
            {
                std::vector<int> aVF2;
                aGrD.GetFacesNeighOfFace(aVF2,aKF1);
                for (const auto &  aKF2 : aVF2)
                {
                    if ( ZBufLabIsOk(aZB.ResSurfD(aKF2).mResult) )
                    {
                        // Got 1 => this face is likely , and we must prolongate the global process
                        aZB.ResSurfD(aKF1).mResult =  eZBufRes::LikelyVisible;
                        GoOn = true;
                    }
                }
            }
        }
    }
}

void cAppliCloudProjectOnGround::MakeBasculeTris(cZBuffer & aZB)
{
    cPt3di aVInd;
    cPt3di aTriCorrel;
    cPt2dr anOffX(mGSD, 0.0);
    cPt2dr anOffY(0.0,-mGSD);
    cAffin2D<tREAL8> anAffinetoTarget(cPt2dr(mBoxLocTarget.P0().x(),mBoxLocTarget.P1().y()),
                                      anOffX,
                                      anOffY);
    cAffin2D<tREAL8> anInvAfftoPixel= anAffinetoTarget.MapInverse();

    cBox2dr aBoxTarget= mBoxLocTarget.CurBox();
    //Write aBoxTarget without rounding 
    cPt2di aSzTarget = Pt_round_up(aBoxTarget.Sz()/mGSD);

    mImRed.DIm().Resize(aSzTarget,eModeInitImage::eMIA_Null);
    mImRed.DIm().InitCste(mInfty);

    mImMasqOut.DIm().Resize(aSzTarget,eModeInitImage::eMIA_Null);

    if (mBascCorrel)
        mImCorrelOut.DIm().Resize(aSzTarget,eModeInitImage::eMIA_Null);


    // iterate over all over triangles and ( add : confidence and occlusion info
    int aNInOut=0;
    int aNGood=0;
    int aNHidden=0;
    int aNNoPix=0;


    for (size_t  aKF=0; aKF<mTri3D->NbFace(); aKF++)
    {
        if(aZB.ResSurfD(aKF).mResult == eZBufRes::Distorted )
        {
            aNInOut++;
        }

        if( aZB.ResSurfD(aKF).mResult==eZBufRes::Hidden)
        {
            aNHidden++;
        }

        if( aZB.ResSurfD(aKF).mResult==eZBufRes::NoPix)
        {
            aNNoPix++;
        }

        if(ZBufLabIsOk(aZB.ResSurfD(aKF).mResult))
        {
            aNGood++;
            // apply affine transform
            tTri3dr aTri3DW = mTri3D->KthTri(aKF);

            if (mBascCorrel)
            {
                aVInd = mTri3D->KthFace(aKF);
                aTriCorrel=cPt3di(mVCorrel[aVInd[0]],mVCorrel[aVInd[1]],mVCorrel[aVInd[2]]);
            }

            cTriangle2DCompiled<tREAL8>  aTriComp(anInvAfftoPixel.Value(Proj(aTri3DW.Pt(0))),
                                                 anInvAfftoPixel.Value(Proj(aTri3DW.Pt(1))),
                                                 anInvAfftoPixel.Value(Proj(aTri3DW.Pt(2))));

            cPt3dr aElev(aTri3DW.Pt(0).z(),aTri3DW.Pt(1).z(),aTri3DW.Pt(2).z());


            std::vector<cPt2di> aVPix;
            std::vector<cPt3dr> aVW;


            aTriComp.PixelsInside(aVPix,1e-10,&aVW);

            for (size_t aK=0; aK<aVPix.size();aK++)
            {
                const cPt2di aPix = aVPix[aK];

                tREAL8 aNewZ =  Scal(aElev,aVW[aK]);
                mImRed.DIm().SetV(aPix,aNewZ);

                if( mBascCorrel)
                {
                    tREAL8 aWeightStretch = std::min(1.0,1/aZB.ResSurfD(aKF).mStretchThresh);
                    tREAL8 aWeightedCorrel = aWeightStretch * Scal(ToR(aTriCorrel),aVW[aK]);

                    // mask out pixels whose weighted correlation is near 0
                    mImMasqOut.DIm().SetV(aPix, (std::abs(aWeightedCorrel) < 1e-6) ? 0 : 1);

                    mImCorrelOut.DIm().SetV(aPix,
                                            std::min(255,round_ni(aWeightedCorrel))
                                            );
                }
                else
                {
                    mImMasqOut.DIm().SetV(aPix,1);
                }
            }
        }
    }

    StdOut()<<"NB GOOD "<<aNGood<<" NB BAD "<<aNInOut<<" All Faces "<<mTri3D->NbFace()<<std::endl;
    // write individual images

    StdOut()<<"NB HIDDEN "<<aNHidden<<" NB NO PIX "<<aNNoPix<<std::endl;


    cDataFileIm2D aDF= cDataFileIm2D::Create(mPhProj.DPMeshDev().FullDirOut()+
                                                  "BLOC-"+
                                                  ToStr(mIndBoxRecal.x())+"-"+
                                                  ToStr(mIndBoxRecal.y())+"-"+
                                                  FileOfPath(mNameResult,false),
                                                  eTyNums::eTN_REAL8,
                                                  mImRed.DIm().Sz(),
                                                  1);

    mImRed.DIm().Write(aDF,aDF.P0());

    cDataFileIm2D aDFM= cDataFileIm2D::Create(mPhProj.DPMeshDev().FullDirOut()+
                                                   "MASQ-"+
                                                   ToStr(mIndBoxRecal.x())+"-"+
                                                   ToStr(mIndBoxRecal.y())+"-"+
                                                   FileOfPath(mNameResult,false),
                                                   eTyNums::eTN_U_INT1,
                                                   mImMasqOut.DIm().Sz(),
                                                   1);
    mImMasqOut.DIm().Write(aDFM,aDFM.P0());

    // write tFW DATA
    GenTFW(anAffinetoTarget, mPhProj.DPMeshDev().FullDirOut()+
                                 "BLOC-"+
                                 ToStr(mIndBoxRecal.x())+"-"+
                                 ToStr(mIndBoxRecal.y())+"-"+
                                 ChgPostix(FileOfPath(mNameResult,false),"tfw"));

    GenTFW(anAffinetoTarget, mPhProj.DPMeshDev().FullDirOut()+
                                 "MASQ-"+
                                 ToStr(mIndBoxRecal.x())+"-"+
                                 ToStr(mIndBoxRecal.y())+"-"+
                                 ChgPostix(FileOfPath(mNameResult,false),"tfw"));


    /// CORREL
    if( mBascCorrel)
    {
        cDataFileIm2D aDFC= cDataFileIm2D::Create(mPhProj.DPMeshDev().FullDirOut()+
                                                       "CORR-"+
                                                       ToStr(mIndBoxRecal.x())+"-"+
                                                       ToStr(mIndBoxRecal.y())+"-"+
                                                       FileOfPath(mNameResult,false),
                                                       eTyNums::eTN_U_INT1,
                                                       mImCorrelOut.DIm().Sz(),
                                                       1);
        mImCorrelOut.DIm().Write(aDFC,aDFC.P0());

        GenTFW(anAffinetoTarget, mPhProj.DPMeshDev().FullDirOut()+
                                     "CORR-"+
                                     ToStr(mIndBoxRecal.x())+"-"+
                                     ToStr(mIndBoxRecal.y())+"-"+
                                     ChgPostix(FileOfPath(mNameResult,false),"tfw"));
    }

}



void cAppliCloudProjectOnGround::MergeResults()
{
    mDFI2d = cDataFileIm2D::Create(mNameIm,mIsGray);
    cParseBoxInOut<2> aPBIO =  cParseBoxInOut<2>::CreateFromSize(mDFI2d,mSzTiles);

    //cPt2dr aP0 (1e10,-1e10);
    std::vector<cBox2di> aSzLocTiles;
    std::vector<cAffin2D<tREAL8>> aLocAffOut;

    /// COMPUTE GLOBAL CONTEXT

    for( const auto & PixI: aPBIO.BoxIndex())
    {
        //StdOut()<<PixI<<std::endl;        // READ CORRESPONDING LOC AFFINE TRANSFORMATIONS

        std::string aNameTileTfW = mPhProj.DPMeshDev().FullDirOut()+
                                              "MASQ-"+
                                              ToStr(PixI.x())+"-"+
                                              ToStr(PixI.y())+"-"+
                                              ChgPostix(FileOfPath(mNameResult,false),"tfw");

        if (! MMVII::ExistFile(aNameTileTfW)) 
            continue;

        cAffin2D<tREAL8> mTrfLocBox = ReadTFW(aNameTileTfW);

        aLocAffOut.push_back(mTrfLocBox);
        // read masq files to get the extent of the number of pixels
        // parse all images and fill global content image
        std::string aNameMasq = mPhProj.DPMeshDev().FullDirOut()+
                                "MASQ-"+
                                ToStr(PixI.x())+"-"+
                                ToStr(PixI.y())+"-"+
                                FileOfPath(mNameResult,false) ;

        cDataFileIm2D mDF = cDataFileIm2D::Create(aNameMasq,eForceGray::No);

        // read mask to get useful masq info
        cIm2D<tU_INT1> aMasqDalle(mDF.Sz()); 
        cBox2di aBoxUsefullMasq = BoxUtile(aMasqDalle);

        aSzLocTiles.push_back(aBoxUsefullMasq);

        cPt2dr aPUL = mTrfLocBox.Value(ToR(aBoxUsefullMasq.P0()));
        cPt2dr aPLR = mTrfLocBox.Value(ToR(aBoxUsefullMasq.P1()));

        mBoxGlobTarget.Add(cPt2dr(aPUL.x(),aPLR.y()));
        mBoxGlobTarget.Add(cPt2dr(aPLR.x(),aPUL.y()));
    }


    mGlobAff=cAffin2D<tREAL8>(cPt2dr(mBoxGlobTarget.CurBox().P0().x(),
                                     mBoxGlobTarget.CurBox().P1().y()),
                              cPt2dr(mGSD,0),
                              cPt2dr(0,-mGSD));



    StdOut()<<"BOX GLOB TARGET "<<mBoxGlobTarget.CurBox()<<std::endl;
    mBoxGlobOutPix=cPt2di(Pt_round_up(mBoxGlobTarget.CurBox().Sz()/mGSD));
    StdOut()<<"BOX GLOB OUT PIX "<<mBoxGlobOutPix<<std::endl;

    mNameBascOut = mOutDir+"Prof_"+ 
                                FileOfPath(mNameResult,false);

    mNameMasqOut = mOutDir+"Masq_"+ 
                                FileOfPath(mNameResult,false);

    std::string aNameTFWProfGlb =  mOutDir+"Prof_"+ 
                                    ChgPostix(FileOfPath(mNameResult,false),"tfw");
    std::string aNameTFWMasqGlb =  mOutDir+"Masq_"+ 
                                    ChgPostix(FileOfPath(mNameResult,false),"tfw");


    GenTFW(mGlobAff,aNameTFWProfGlb);
    GenTFW(mGlobAff,aNameTFWMasqGlb);

    if(mBascCorrel)
    {
        std::string aNameTFWCorrelGlb =  mOutDir+"Correl_"+ 
                                            ChgPostix(FileOfPath(mNameResult,false),"tfw");
        GenTFW(mGlobAff,aNameTFWCorrelGlb);
    }

    cDataFileIm2D  aDF = cDataFileIm2D::Create(mNameBascOut,
                                              eTyNums::eTN_REAL4,
                                              mBoxGlobOutPix.Sz(),
                                              1);

    cIm2D<tREAL4> aGlobIm(aDF.Sz(),aDF);
    aGlobIm.DIm().InitCste(-1e9);


    cDataFileIm2D  aDFM = cDataFileIm2D::Create(mNameMasqOut,
                                              eTyNums::eTN_U_INT1,
                                              mBoxGlobOutPix.Sz(),
                                              1);
    cIm2D<tU_INT1> aGlobMasqIm(aDFM.Sz(),aDFM);

    cIm2D<tU_INT1> aGlobCorrelIm(cPt2di(1,1));
    cDataFileIm2D aDFC=cDataFileIm2D::Empty();

    if (mBascCorrel)
    {
        mNameCorrelOut= mOutDir+"Correl_"+ FileOfPath(mNameResult,false);
        aDFC = cDataFileIm2D::Create(mNameCorrelOut,
                                                   eTyNums::eTN_U_INT1,
                                                   mBoxGlobOutPix.Sz(),
                                                   1);
        aGlobCorrelIm=cIm2D<tU_INT1>(aDFC.Sz(),aDFC);
    }

    int aK=0;
    cPt2di aP0G,aP1G;
    for( const auto & PixI: aPBIO.BoxIndex())
    {
        StdOut()<<PixI<<std::endl;
        std::string aNameProfDalle = mPhProj.DPMeshDev().FullDirOut()+
                                "BLOC-"+
                                ToStr(PixI.x())+"-"+
                                ToStr(PixI.y())+"-"+
                                FileOfPath(mNameResult,false) ;
        if (! MMVII::ExistFile(aNameProfDalle))
            continue;

        std::string aNameMasqDalle = mPhProj.DPMeshDev().FullDirOut()+
                                "MASQ-"+
                                ToStr(PixI.x())+"-"+
                                ToStr(PixI.y())+"-"+
                                FileOfPath(mNameResult,false) ;

        aP0G = ToI(mGlobAff.Inverse(aLocAffOut[aK].Value(ToR(aSzLocTiles[aK].P0()))));
        aP1G = ToI(mGlobAff.Inverse(aLocAffOut[aK].Value(ToR(aSzLocTiles[aK].P1()))));

        cBox2di aBoxLocInGlob(aP0G,aP1G);
        StdOut()<<"LOC "<<aBoxLocInGlob<<" "<<"GLOB "<<mBoxGlobOutPix<<std::endl;

        cIm2D<tREAL4>  aImProfDalle(aSzLocTiles[aK].Sz());
        cIm2D<tU_INT1> aImMasqDalle(aSzLocTiles[aK].Sz());

        aImProfDalle.Read(cDataFileIm2D::Create(aNameProfDalle,eForceGray::No),
                            aSzLocTiles[aK].P0()) ;
        aImMasqDalle.Read(cDataFileIm2D::Create(aNameMasqDalle,eForceGray::No),
                            aSzLocTiles[aK].P0()) ;

        //aImMasqDalle.DIm().Dilate(-mSzOverlap);

        cIm2D<tU_INT1> aImCorrelDalle(cPt2di(1,1));

        if (mBascCorrel)
        {
            std::string aNameCorrelDalle = mPhProj.DPMeshDev().FullDirOut()+
                                         "CORR-"+
                                         ToStr(PixI.x())+"-"+
                                         ToStr(PixI.y())+"-"+
                                         FileOfPath(mNameResult,false) ;
            aImCorrelDalle=cIm2D<tU_INT1>(aSzLocTiles[aK].Sz());
            aImCorrelDalle.Read(cDataFileIm2D::Create(aNameCorrelDalle,eForceGray::No),
                                aSzLocTiles[aK].P0()) ;
        }

        // only select masked depths  ??????
        cDataIm2D<tU_INT1> & aDIMm = aImMasqDalle.DIm();
        for (const auto & aPix: aDIMm)
        {
            if (aDIMm.GetV(aPix))
            {
                tREAL8  aSavedZ = aGlobIm.DIm().GetV(aBoxLocInGlob.P0()+aPix);
                tREAL8 aCurUpdateZ= aImProfDalle.DIm().GetV(aPix);

                if (aCurUpdateZ>aSavedZ)
                {
                    aGlobIm.DIm().SetV(aBoxLocInGlob.P0()+aPix,aCurUpdateZ);
                    aGlobMasqIm.DIm().SetV(aBoxLocInGlob.P0()+aPix,1);

                    if( mBascCorrel)
                    {
                        aGlobCorrelIm.DIm().SetV(aBoxLocInGlob.P0()+aPix,
                                                 aImCorrelDalle.DIm().GetV(aPix));
                    }
                }

            }
        }

        aK++;
    }

    // Write images

    aGlobIm.Write(aDF,cPt2di(0,0));
    aGlobMasqIm.Write(aDFM,cPt2di(0,0));

    if( mBascCorrel)
    {
        aGlobCorrelIm.Write(aDFC,cPt2di(0,0));
    }
}


int cAppliCloudProjectOnGround::ExeOnParsedBox()
{
    StdOut()<<"CURBX "<<CurBoxIn()<<std::endl;
    mImPx1= APBI_ReadIm<tREAL4>(DirProject() + mNameCloud2D_DepthIn);
    mImMasq1 = ReadMasqWithDef(CurBoxIn(),DirProject() + mNameMasq1);

    // check if mNameCorrel is provided, if so read the correlation image and set mBascCorrel to true
    if (!mNameCorrel.empty())
    {
        mBascCorrel=true;
        mImCorrel = APBI_ReadIm<tU_INT1>(DirProject() + mNameCorrel);
    }

    MakeFastBasc();


    return EXIT_SUCCESS;
}



std::vector<cOneHelpSampleCmp> cAppliCloudProjectOnGround::Samples() const
{

    return 
    {
        cOneHelpSampleCmp::Header(
            "Project a disparity map to ground given an xml file of this structure\n"
            "    <Data>\n"
            "        <ModeGeom>DEPTH</ModeGeom>\n"
            "        <NameEpipIm1>Epip1.png</NameEpipIm1>\n"
            "        <NameEpipIm2>Epip2.png</NameEpipIm2>\n"
            "        <NameDisparity>DenseMatchResult/disparity.tif</NameDisparity>\n"
            "        <NameCorrel>DenseMatchResult/Correlation.tif</NameCorrel>\n"
            "        <NameMasq>DenseMatchResult/occlusion_mask.tif</NameMasq>\n"
            "    </Data>"),
            {"MMVII CloudProjectOnGround /path/to/DMEpipPairResult.xml ori BascTmpFiles GSD=0.1"}
    };
}

int cAppliCloudProjectOnGround::Exe()
{
    // In father-dir mode the orientation folder lives at the project root, not under the
    // epipolar-resampling sub-folder that DirProject() points to. Prevent FinishInit() from
    // validating the orientation input dir against the (wrong) standard project-relative path.
    if (mOriInFatherDir)
        mPhProj.DPOrient().SetAllowDirInEmpty();

    mPhProj.FinishInit();

    // Now that all dirs are finalized, redirect the orientation input to the father dir.
    if (mOriInFatherDir)
    {
        std::string aFatherDir = DirProject();
        size_t aPos = aFatherDir.find(MMVII_DirPhp);
        if (aPos != std::string::npos)
            aFatherDir = aFatherDir.substr(0,aPos);

        std::string aOriFullDir =   aFatherDir
                              + mPhProj.DPOrient().DirLocOfMode()   // "MMVII-PhgrProj/Ori/"
                              + mPhProj.DPOrient().DirIn()
                              + StringDirSeparator();

        mPhProj.DPOrient().SetFullDirIn(aOriFullDir);
    }


    // Read mNameDMEpipPairResult
    cDMEpipPairResult aDMEpipPairResult;
    ReadFromFile_Std(aDMEpipPairResult,
            mNameDMEpipPairResult);

    // Instantiate all strings from the DMEpipPairResult
    mModeGeom = Str2E<eModeGeom>(aDMEpipPairResult.ModeGeom());
    mNameIm =  aDMEpipPairResult.NameEpipIm1();
    mNameSec =  aDMEpipPairResult.NameEpipIm2();
    mNameCloud2D_DepthIn = aDMEpipPairResult.NameDisparity();
    mNameMasq1  =  aDMEpipPairResult.NameMasq();
    mNameCorrel =  aDMEpipPairResult.NameCorrel();

    if (! IsInit(&mOutDir))
    {
        mOutDir = mPhProj.DirVisuAppli();;
    }
    if (! mOutDir.empty())
    {
        mOutDir += StringDirSeparator();
    }

    CreateDirectories(mOutDir);

    // check if mNameResult is Init, if not use default name
    if (! IsInit(&mNameResult))
    {
        mNameResult = "Dem_" + APBI_NameIm();
        if (Postfix(mNameResult) != "tif" && Postfix(mNameResult) != "tiff")
        {
            mNameResult += ".tif";
        }
    }

    // read camera
    mCamGen = mPhProj.ReadSensor(APBI_NameIm(),true);

    if (mModeGeom== eModeGeom::eGEOM_EPIP)
        {
            MMVII_INTERNAL_ASSERT_strong( (!mNameSec.empty() ),"should provide secondary image ");
            mSecCamGen = mPhProj.ReadSensor(mNameSec,true);
            mNameSec = DirProject() + mNameSec;
        }

    // restore name of input image, as it is used in the APBI_ExecAll() function
    mNameIm = DirProject() + mNameIm;

    //Compute GSD 
    
    cAutoTimerSegm aTGSDCompute(TimeSegm(),"ComputeGSD"); 
    if ( !IsInit(&mGSD) )
    {
        mGSD = ComputeGSD();
        mBoxLocTarget = cTplBoxOfPts<tREAL8,2>();
        // Share the auto-estimated GSD with the parallel tile subprocesses on their command line,
        // so each recall reuses it (IsInit true) instead of recomputing it.
        mArgsRecall.push_back("GSD=" + ToStr(mGSD));
    }
    cAutoTimerSegm aTBascAll(TimeSegm(),"Projection");


    APBI_ExecAll();

    // Merge all results of bascul
    if (!InsideParalRecall())
    {

        MergeResults();

        // REMOVE NON NECESSARY INDIVIDUAL TILES
        RemovePatternFile(mPhProj.DPMeshDev().FullDirOut()+"BLOC.*",false);
        RemovePatternFile(mPhProj.DPMeshDev().FullDirOut()+"MASQ.*",false);
        RemovePatternFile(mPhProj.DPMeshDev().FullDirOut()+"CORR.*",false);
    }

    eTypeSerial aTsOut= eTypeSerial::exml;
    std::string aNameXmlOut= mNameBascOut + "." + E2Str(aTsOut);
    // save Serialized Info about bascule
    cCloudRaster aCldRaster(mNameBascOut,
                          mNameCorrelOut,
                          mNameMasqOut,
                          mBoxGlobOutPix.Sz(),
                          mGlobAff);
    SaveInFile(aCldRaster, aNameXmlOut);
    return EXIT_SUCCESS;
}



/*  ============================================= */
/*       ALLOCATION                               */
/*  ============================================= */

tMMVII_UnikPApli Alloc_CloudProjectOnGround(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli & aSpec)
{
    return tMMVII_UnikPApli(new cAppliCloudProjectOnGround(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpecCloudProjectOnGround
    (
        "CloudProjectOnGround",
        Alloc_CloudProjectOnGround,
        "Project of Pax/Depth maps to ground",
        {eApF::Cloud},
        {eApDT::Image},
        {eApDT::Image},
        __FILE__
        );

void BenchCloudEpipProjectOnGround(cParamExeBench & aParam)
{
    if (! aParam.NewBench("CloudEpipProjectOnGround")) return;
    
    cMMVII_Appli & anAp = cMMVII_Appli::CurrentAppli();
    const std::string & aInDir = anAp.InputDirTestMMVII() + "/EPIP/Synthetic/";
    const std::string & aTmpDir = anAp.TmpDirTestMMVII();

    const std::string aNameDMEpipPairResult = aInDir + "DenseMatchResult/DMEpipPairResult.xml";
    const std::string aOutDir = aTmpDir + "CloudEpipProjectOnGround/";
    const std::string aFullOriDir =aInDir + "MMVII-PhgrProj/Ori/Synthetic";  
    // read aNameDMEpipPairResult
    cDMEpipPairResult aDMPRes;
    ReadFromFile_Std(aDMPRes,aNameDMEpipPairResult);

    const eModeGeom aModeGeom = Str2E<eModeGeom>(aDMPRes.ModeGeom());
    const std::string aNameIm =  aDMPRes.NameEpipIm1();
    const std::string aNameSec =  aDMPRes.NameEpipIm2();
    const std::string aNameCloud2D_DepthIn = aDMPRes.NameDisparity();
    const std::string aNameMasq1  =  aDMPRes.NameMasq();
    const std::string aNameCorrel =  aDMPRes.NameCorrel();

    MMVII_INTERNAL_ASSERT_bench(aModeGeom==eModeGeom::eGEOM_EPIP,
        "When using epipolair pairs with disparity, the eModeGeom should be eGEOM_EPIP");

    // Read cameras  and project one triangle, get relevant ground truth height and compare

    const std::string aFullNameOriIm= aFullOriDir + StringDirSeparator() + cSensorCamPC::NameOri_From_Image(aNameIm);
    const std::string aFullNameOriSec = aFullOriDir + StringDirSeparator() + cSensorCamPC::NameOri_From_Image(aNameSec);
    const std::string aFullNamePax = aInDir + aNameCloud2D_DepthIn;
    const std::string aFullNameGnd = aInDir + "GroundTruth_Dem_Epip1.tif";
    const std::string aFullNameGndTfw = aInDir + "GroundTruth_Dem_Epip1.tfw";

    // read camers 
    cSensorCamPC * aCamEpip1 = cSensorCamPC::FromFile(aFullNameOriIm, true);
    cSensorCamPC * aCamEpip2 = cSensorCamPC::FromFile(aFullNameOriSec,true);

    cMMVII_Appli::AddObj2DelAtEnd(aCamEpip1);
    cMMVII_Appli::AddObj2DelAtEnd(aCamEpip2);


    // read the disparity map
    cIm2D<tREAL4> aImPax = cIm2D<tREAL4>(aCamEpip1->SzPix());
    aImPax.Read(cDataFileIm2D::Create(aFullNamePax,
                                    eTyNums::eTN_REAL4,
                                    aCamEpip1->SzPix()),cPt2di(0,0));

        
    //read ground truth dem 
    cIm2D<tREAL4> aImGnd = cIm2D<tREAL4>::FromFile(aFullNameGnd);
    cAffin2D<tREAL8> anAffGT = cAppliCloudProjectOnGround::ReadTFW(aFullNameGndTfw);

    std::list<cPt2di> aVSetOfPoints = {cPt2di(856,372),cPt2di(494,369),cPt2di(488,648)};

    // project some triangles, get relevant ground truth height and compare
    cPt2di Up[2]={cPt2di(1,0),cPt2di(1,1)};
    cPt2di Down[2]={cPt2di(1,1),cPt2di(0,1)};
    for (const auto aTriOff : {Up,Down})
    {
        for (const auto aP1: aVSetOfPoints)
        {
            cPt2di aP2= aP1+aTriOff[0];
            cPt2di aP3= aP1+aTriOff[1];

            std::vector<cPt2di> aVTri = {aP1,aP2,aP3} ;
            std::vector<cPt3dr> aVTriGround;

            for (auto aP: aVTri)
            {
                tREAL4 aPaxV = aImPax.DIm().GetV(aP);
                cPt2dr aPInEp2 = cPt2dr(aPaxV,0) + ToR(aP) ;

                tSeg3dr aS3OfP = aCamEpip1->Image2Bundle(ToR(aP));
                tSeg3dr aS3OfPInEp2 =  aCamEpip2->Image2Bundle(aPInEp2);

                cPt3dr aP3OfP = RobustBundleInters({aS3OfP,aS3OfPInEp2});
                
                aVTriGround.push_back(aP3OfP) ;
            }   
            tTri3dr aTriC(aVTriGround[0],aVTriGround[1],aVTriGround[2]);



            cTriangle2DCompiled<tREAL8> aTri2D(anAffGT.Inverse(Proj(aTriC.Pt(0))),
                                            anAffGT.Inverse(Proj(aTriC.Pt(1))),
                                            anAffGT.Inverse(Proj(aTriC.Pt(2))));

            cPt3dr aElev(aTriC.Pt(0).z(),aTriC.Pt(1).z(),aTriC.Pt(2).z());

            std::vector<cPt2di> aVPix;
            std::vector<cPt3dr> aVW;

            tREAL8 aScale = 1.0;
            for (int aKp=0 ; aKp<3 ; aKp++)
            {
                aScale = std::max(aScale,std::abs(aTri2D.Pt(aKp).x()));
                aScale = std::max(aScale,std::abs(aTri2D.Pt(aKp).y()));
            }
            tREAL8 aTol = -8.0 * std::numeric_limits<tREAL8>::epsilon() * aScale;

            aTri2D.PixelsInside(aVPix,aTol,&aVW); 

            //MMVII_INTERNAL_ASSERT_bench(aVPix.size()>0, "cTriangle2DCompiled::PixelsInside yields no pixels inside Triangle");
            for (int aK=0; aK<(int)aVPix.size(); aK++)
            {
                MMVII_INTERNAL_ASSERT_bench(aImGnd.DIm().Inside(aVPix[aK]), 
                                "Pixel inside triangle is outside ground truth image");
                // get the ground truth height
                tREAL4 aGTZ = aImGnd.DIm().GetV(aVPix[aK]);
                tREAL4 aBascZ = Scal(aElev,aVW[aK]);
                std::cout<<"diff between zgt and zbasc "<<aGTZ<<" "<<aBascZ<<std::endl;
                MMVII_INTERNAL_ASSERT_bench(std::abs(aGTZ-aBascZ)<0.01, 
                            "Pixel inside triangle height does not match ground truth");
            }

        }
    }
    aParam.EndBench();
};


void BenchRunCloudProjectOnGround(cParamExeBench & aParam)
{
    if (! aParam.NewBench("RunCloudProjectOnGround")) return;

    cMMVII_Appli & anAp = cMMVII_Appli::CurrentAppli();
    const std::string & aInDir = anAp.InputDirTestMMVII() + "/EPIP/Synthetic/";
    const std::string & aTmpDir = anAp.TmpDirTestMMVII();

    const std::string aNameDMEpipPairResult = "DenseMatchResult/DMEpipPairResult.xml";
    const std::string aOutDir = aTmpDir + "RunCloudProjectOnGround/";
    const std::string aFullOriDir =aInDir + "MMVII-PhgrProj/Ori/Synthetic";
    tREAL8 aNoData = -1e9;


    // run the appli 
    anAp.ExeCallMMVII
    (
        TheSpecCloudProjectOnGround,
        anAp.StrObl() << aInDir + StringDirSeparator() + aNameDMEpipPairResult << "Synthetic" << "TmpTiles",
        anAp.StrOpt() << std::make_pair("GSD","0.01") << std::make_pair("DirProj",aInDir) <<std::make_pair("DirOut",aOutDir)
    );

    // compare projected dem with ground truth 
    const std::string aFullNameGnd = aInDir + "GroundTruth_Dem_Epip1.tif";
    const std::string aFullNameGndTfw = aInDir + "GroundTruth_Dem_Epip1.tfw";

    const std::string aFullNamePredicted = aOutDir+"Prof_Dem_Epip1.png.tif";
    const std::string aFullNamePredictedTfw = aOutDir + "Prof_Dem_Epip1.png.tfw";

    MMVII_INTERNAL_ASSERT_bench((ExistFile(aFullNamePredicted) && ExistFile(aFullNameGndTfw)),
                    "CloudProjectOnGround did not write the predicted projected dem to "+aOutDir);

    cIm2D<tREAL4> aImGnd  = cIm2D<tREAL4>::FromFile(aFullNameGnd) ;
    cIm2D<tREAL4> aImPred = cIm2D<tREAL4>::FromFile(aFullNamePredicted);

    // TFW 
    cAffin2D<tREAL8> anAffGt   = cAppliCloudProjectOnGround::ReadTFW(aFullNameGndTfw);
    cAffin2D<tREAL8> anAffPred = cAppliCloudProjectOnGround::ReadTFW(aFullNamePredictedTfw);

    tREAL8 aSumAbsDiff = 0.0 ;
    tREAL8 aWeighter = 0.0;

    for (const auto aP : aImPred.DIm())
    {
        tREAL4 aValPred = aImPred.DIm().GetV(aP);

        if (aValPred != aNoData)
        {
            cPt2dr aPGlobPred = anAffPred.Value(ToR(aP));
            cPt2dr aPLocInGt  = anAffGt.Inverse(aPGlobPred);
            tREAL4 aValGt     = aImGnd.DIm().GetVBL(aPLocInGt);
            if (aValGt>0)
            {
                //std::cout<<aValPred<<" "<<aValGt<<std::endl;
                aSumAbsDiff+=std::abs(aValPred-aValGt);
                aWeighter+=1.0;
            }

        }
    }

    aSumAbsDiff/=aWeighter;
    MMVII_INTERNAL_ASSERT_bench(aSumAbsDiff<0.015, 
            "Average Height difference between computed Dem and ground truth dem is above 1.5 cm");

    aParam.EndBench();        
};


};// MMVII namespace

