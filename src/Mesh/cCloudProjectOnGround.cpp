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

        std::string mNameCloud2D_DepthIn;
        std::string mDirOri;
        std::string mNameMasq1;
        std::string mNameCorrel;
        std::string mNameSec;

        cSensorCamPC * mCamPC;
        cSensorCamPC * mSecCamPC;
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
        bool mZF_SameOri;
        bool mBascCorrel;
        bool mOriInFatherDir;  ///< if true, search orientation in father dir, else in DirProj
        int  mMultZ;
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
        std::pair<cPt3dr,cPt3dr> BascOnePoint(cPt2di A,  cPt2di anOffSet, bool & oValid);
        /// Read a camera whose orientation folder is stored at the project root (father dir),
        /// and not under DirProj which here points to the epipolar-resampling sub-folder.
        cSensorCamPC * ReadCamFromFatherDir(const std::string & aNameIm);
        //void MakeBasc();
        void MakeFastBasc();
        void MakeBasculeTris(cZBuffer & aZB);
        void ProcessNoPix(cZBuffer &  aZB);
        void GenTFW(const cAffin2D<tREAL8> & anAff, const std::string & aNameTFW);
        cAffin2D<tREAL8> ReadTFW(const std::string & aNameTFW);
        cBox2di  BoxUtile(cIm2D<tU_INT1> & anImMasq);
        void MergeResults();
};


cAppliCloudProjectOnGround::cAppliCloudProjectOnGround(const std::vector<std::string> & aVArgs, const cSpecMMVII_Appli & aSpec ):
    cMMVII_Appli(aVArgs,aSpec),
    cAppliParseBoxIm<tREAL4>(*this,eForceGray::No,cPt2di(2000,2000),cPt2di(50,50),true),
    mCamPC(nullptr),
    mSecCamPC(nullptr),
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
    mZF_SameOri(true),
    mBascCorrel(false),
    mOriInFatherDir(false),
    mMultZ(mZF_SameOri ? 1 : -1),
    mMII(0.0),
    mStretschingThresh(4.0)
{
}


cCollecSpecArg2007 & cAppliCloudProjectOnGround::ArgObl(cCollecSpecArg2007 & anArgObl)
{
    return
            APBI_ArgObl(anArgObl)
           <<   Arg2007(mNameCloud2D_DepthIn,"Name of input depth map", {eTA2007::FileImage} )
           <<   Arg2007(mDirOri,"Mandatory directory of orientation files (stored under the project root)")
           <<   mPhProj.DPMeshDev().ArgDirOutMand()
        ;
}


cCollecSpecArg2007 & cAppliCloudProjectOnGround::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
    return
        APBI_ArgOpt
        (
            anArgOpt
                << AOpt2007(mModeGeom,"ModeGeom","Either epipolar geometry, def=geometry of input depth",{AC_ListVal<eModeGeom>()})
                << AOpt2007(mNameResult,"Out"," prefix of output files",{eTA2007::HDV})
                << AOpt2007(mNameMasq1, "Masq1","Masq of first image if any",{eTA2007::HDV})
                << AOpt2007(mNameSec, "Im2","Secondary image if geom is Epip",{eTA2007::HDV})
                << AOpt2007(mNameCorrel,"ImCorrel","Name of correlation or confidence image")
                << AOpt2007(mGSD,"GroundResolution", "Ground sampling distance of bascule")
                << AOpt2007(mOutDir,"DirOut","Directory of output files, (def=VISU/"+Specs().Name()+")")
                << AOpt2007(mMII,"MII","Margin Inside Image (for triangle validation)", {eTA2007::HDV})
                << AOpt2007(mOriInFatherDir,"OriInFatherDir","If true, search orientation in father dir, else in DirProj",{eTA2007::HDV})
                << AOpt2007(mStretschingThresh,"ThresholdDistortion","Level of triangle distortion to discard from bascule")
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


std::pair<cPt3dr,cPt3dr> cAppliCloudProjectOnGround::BascOnePoint(cPt2di A,  cPt2di anOffSet, bool & oValid)
{
    cPt3dr aPZA, aPWithDepth;

    cPt2di  aPix1 = A+anOffSet;
    oValid = true;

    if(mModeGeom==eModeGeom::eGEOM_EPIP)  ///< entry is a disparity maps so need to compute pseudo intersection
    {
        double  aX2= aPix1.x() + mImPx1.DIm().GetV(A);
        cPt2dr aPix2=cPt2dr(aX2, aPix1.y());
        tSeg3dr aP2TerCam1 = mCamPC->Image2Bundle(ToR(aPix1));
        tSeg3dr aP2TerCam2= mSecCamPC->Image2Bundle(aPix2);

        aPZA=RobustBundleInters({aP2TerCam1,aP2TerCam2});

        // compute 2D + Depth and no need to perform Value of mapI2O Value or TriValue as it is already computed

        // aPWithDepth = Value (aPZA)

        aPWithDepth= cPt3dr(aPix1.x(),
                             aPix1.y(),
                             mCamPC->Pose().Inverse(aPZA).z());
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

        aPZA = mCamPC->ImageAndDepth2Ground(aPWithDepth);
    }
    else
    {
        ///< nothing for now
    }
    mBoxLocTarget.Add(Proj(aPZA));

    return {aPWithDepth,aPZA};
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

    cSIMap_Ground2ImageAndProf aMapCamDepth(mCamPC);

    cSetVisibility aSetVis(mCamPC,mMII);

    double Infty =1e20;

    cBox3dr  aBox(cPt3dr(aP0.x(),aP0.y(),-Infty),cPt3dr(aP1.x(),aP1.y(),Infty));

    cDataBoundedSet<tREAL8,3>  aSetCam(aBox);

    StdOut()<<"ZBUFFER::ZBUFFER"<<std::endl;
    cZBuffer aZBuf(aTriIt,
                   aSetVis,
                   aMapCamDepth,
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
    
    cAutoTimerSegm aTSBufferPixellizer(TimeSegm(),"rasterization-pixellization"); 
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
    std::string aline;
    std::vector< std::string > acontent;
    while(std::getline(aFtfw,aline))
    {
        acontent.push_back(aline);
    }
    aFtfw.close();

    tREAL8 gsd_x=stof(acontent.at(0));
    tREAL8 gsd_y=stof(acontent.at(3));
    tREAL8 x_ul=stof(acontent.at(4));
    tREAL8 y_ul=stof(acontent.at(5));

    acontent.clear();

    return cAffin2D<tREAL8>(cPt2dr(x_ul,y_ul),cPt2dr(gsd_x,0),cPt2dr(0,gsd_y));
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
            //StdOut()<<"hidden"<<std::endl;
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

                tREAL8 aNewZ = mMultZ * Scal(aElev,aVW[aK]);
                mImRed.DIm().SetV(aPix,aNewZ);
                mImMasqOut.DIm().SetV(aPix,1);

                if( mBascCorrel)
                {
                    tREAL8 aWeightStretch = std::min(1.0,1/aZB.ResSurfD(aKF).mStretchThresh);
                    mImCorrelOut.DIm().SetV(aPix,
                                            std::min(255,round_ni(aWeightStretch * Scal(ToR(aTriCorrel),aVW[aK])))
                                            );
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



    mBoxGlobOutPix=cPt2di(Pt_round_up(mBoxGlobTarget.CurBox().Sz()/mGSD));

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

cSensorCamPC * cAppliCloudProjectOnGround::ReadCamFromFatherDir(const std::string & aNameIm)
{
    // When mOriInFatherDir is set, DirProj points to the epipolar-resampling working sub-folder
    // (where the images live), but the orientation folder is stored at the project root (father
    // dir). So we anchor the orientation folder to the father dir : everything in DirProject()
    // before "MMVII-PhgrProj/". Otherwise we keep looking under DirProj (standard behaviour).
    std::string aFatherDir = DirProject();
    if (mOriInFatherDir)
    {
        size_t aPos = aFatherDir.find(MMVII_DirPhp);
        if (aPos != std::string::npos)
            aFatherDir = aFatherDir.substr(0,aPos);
    }

    std::string aOriFullDir =   aFatherDir
                              + mPhProj.DPOrient().DirLocOfMode()   // "MMVII-PhgrProj/Ori/"
                              + mDirOri
                              + StringDirSeparator();

    std::string aNameCam = aOriFullDir + cSensorCamPC::NameOri_From_Image(aNameIm);

    cSensorCamPC * aCamPC = cSensorCamPC::FromFile(aNameCam,true);
    cMMVII_Appli::AddObj2DelAtEnd(aCamPC);

    return aCamPC;
}

int cAppliCloudProjectOnGround::ExeOnParsedBox()
{
    StdOut()<<"CURBX "<<CurBoxIn()<<std::endl;
    mImPx1= APBI_ReadIm<tREAL4>(DirProject() + mNameCloud2D_DepthIn);
    mImMasq1 = ReadMasqWithDef(CurBoxIn(),DirProject() + mNameMasq1);

    if (IsInit(&mNameCorrel))
    {
        mBascCorrel=true;
        mImCorrel = APBI_ReadIm<tU_INT1>(DirProject() + mNameCorrel);
    }

    MakeFastBasc();


    return EXIT_SUCCESS;
}


int cAppliCloudProjectOnGround::Exe()
{
    mPhProj.FinishInit();

    if (! IsInit(&mOutDir))
    {
        mOutDir = mPhProj.DirVisuAppli();;
    }
    if (! mOutDir.empty())
    {
        mOutDir += "/";
    }

    CreateDirectories(mOutDir);

    // check if mNameResult is Init, if not use default name
    if (! IsInit(&mNameResult))
        mNameResult = "Dem_" + APBI_NameIm();
    // read camera
    mCamPC = ReadCamFromFatherDir(APBI_NameIm());

    if (mModeGeom== eModeGeom::eGEOM_EPIP)
        {
            MMVII_INTERNAL_ASSERT_strong(IsInit(&mNameSec),"should provide secondary image ");
            mSecCamPC = ReadCamFromFatherDir(mNameSec);
        }

    // restore name of input image, as it is used in the APBI_ExecAll() function
    mNameIm = DirProject() + mNameIm;

    APBI_ExecAll();

    // Merge all results of bascule
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


};
