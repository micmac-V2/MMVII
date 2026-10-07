#include "cEpipolarRectification.h"
#include "MMVII_Error.h"
#include "MMVII_Sensor.h"
#include "MMVII_PCSens.h"
#include "MMVII_Bench.h"
#include "cMMVII_Appli.h"
#include "MMVII_DeclareAllCmd.h"

namespace MMVII {

// ============================================================
//  cEpipMappingPC
// ============================================================

cEpipMappingPC::cEpipMappingPC(const cSensorCamPC & aCam,const cRotation3D<tREAL8> & aEpipRot,tREAL8 aFocal,const cPt2dr & aZInterval)
    : cEpipolarMapping(aZInterval)
    , mCam   (&aCam)
    , mRot   (aEpipRot)
    , mFocal (aFocal)
{
}

cEpipMappingPC::cEpipMappingPC()
    : cEpipolarMapping(cPt2dr(0,0))
    , mCam   (nullptr)
    , mRot   (cRotation3D<tREAL8>::Identity())
    , mFocal (1.0)
{
}

const cSensorCamPC & cEpipMappingPC::Cam() const
{
    MMVII_INTERNAL_ASSERT_always(mCam!=nullptr,"cEpipMappingPC: source camera not bound (SetSourceSensor)");
    return *mCam;
}

void cEpipMappingPC::SetSourceSensor(const cSensorImage & aSensor)
{
    mCam = aSensor.GetSensorCamPC();
    MMVII_INTERNAL_ASSERT_User(mCam!=nullptr,eTyUEr::eUnClassedError,"The epipolar model needs a central perspective camera for image " + aSensor.NameImage());
}

void cEpipMappingPC::AddData(const cAuxAr2007 &anAux)
{
    AddDataBase(anAux);
    MMVII::AddData(cAuxAr2007("Focal",anAux),mFocal);
    cPt3dr aAxes[3];
    for (int aK=0 ; aK<3 ; aK++)
        GetCol(aAxes[aK],mRot.Mat(),aK);
    MMVII::AddData(cAuxAr2007("AxisX",anAux),aAxes[0]);
    MMVII::AddData(cAuxAr2007("AxisY",anAux),aAxes[1]);
    MMVII::AddData(cAuxAr2007("AxisZ",anAux),aAxes[2]);
    if (anAux.Input())
        mRot = cRotation3D<tREAL8>(aAxes[0],aAxes[1],aAxes[2],true);   // text precision may leave a slight non-orthogonality
}

// Mapping of a point where the distortion inversion of the source camera is valid
cPt2dr cEpipMappingPC::ValueInDomain(const cPt2dr& aPt) const
{
    // direction of the bundle in the frame of the virtual camera, then pinhole projection
    const cPt3dr aDir = mRot.Inverse(Cam().Image2Bundle(aPt).V12());
    MMVII_INTERNAL_ASSERT_User(aDir.z()>0,eTyUEr::eUnClassedError,"Epipolar rectification: a bundle is behind the virtual camera");
    return cPt2dr(aDir.x(),aDir.y()) * (mFocal/aDir.z()) - ToR(mEpipImFrame.P0());
}

cPt2dr cEpipMappingPC::Value(const cPt2dr& aPt) const
{
    // Far from the image the distortion cannot be inverted : first order extrapolation from the nearest valid point
    const cBox2dr aBox = Cam().PixelDomain().Box();
    const cBox2dr aBoxOk = aBox.Dilate(0.15*Cam().InternalCalib()->F());
    if (aBoxOk.Inside(aPt))
        return ValueInDomain(aPt);

    const cPt2dr aP0(std::min(std::max(aPt.x(),aBox.P0().x()+2.0),aBox.P1().x()-2.0),std::min(std::max(aPt.y(),aBox.P0().y()+2.0),aBox.P1().y()-2.0));
    const cPt2dr aV0 = ValueInDomain(aP0);
    const cPt2dr aDx = ValueInDomain(aP0+cPt2dr(1,0)) - aV0;
    const cPt2dr aDy = ValueInDomain(aP0+cPt2dr(0,1)) - aV0;
    const cPt2dr aDelta = aPt - aP0;
    return aV0 + aDx*aDelta.x() + aDy*aDelta.y();
}

cPt2dr cEpipMappingPC::Inverse(const cPt2dr& aPt) const
{
    const cPt2dr aPV = aPt + ToR(mEpipImFrame.P0());
    return Cam().Ground2Image(Cam().Center() + mRot.Value(cPt3dr(aPV.x(),aPV.y(),mFocal)));
}

cSensorCamPC * cEpipMappingPC::VirtualCamera(const std::string & aName) const
{
    const cPt2dr aPP = -ToR(mEpipImFrame.P0());
    auto * aCalib = cPerspCamIntrCalib::SimpleCalib(aName+"-Calib",eProjPC::eStenope,EpipImSz(),cPt3dr(aPP.x(),aPP.y(),mFocal),cPt3di(0,0,0));
    cMMVII_Appli::AddObj2DelAtEnd(aCalib);  // a sensor never owns its calib
    return new cSensorCamPC(aName,cIsometry3D<tREAL8>(Cam().Center(),mRot),aCalib);
}

// ============================================================
//  cEpipolarRectificationPC
// ============================================================

cEpipolarRectificationPC::cEpipolarRectificationPC(const cSensorCamPC & aCam1,const cSensorCamPC & aCam2,const cParams & aParams)
    : mCam1   (aCam1)
    , mCam2   (aCam2)
    , mParams (aParams)
{
}

cEpipModelPC cEpipolarRectificationPC::Compute()
{
    // A 360 degrees image has bundles in every direction, no pinhole virtual camera can receive them
    for (const cSensorCamPC * aCam : {&mCam1,&mCam2})
        MMVII_INTERNAL_ASSERT_User(aCam->InternalCalib()->TypeProj()!=eProjPC::eEquiRect,eTyUEr::eUnClassedError,
            "Closed form epipolar rectification is not implemented for EquiRect (360 degrees) cameras, image " + aCam->NameImage());

    // X axis : baseline
    MMVII_INTERNAL_ASSERT_User(Norm2(mCam2.Center()-mCam1.Center())>1e-12,eTyUEr::eUnClassedError,"Epipolar rectification: the two cameras have the same centre");
    const cPt3dr aAxisX = VUnit(mCam2.Center()-mCam1.Center());

    // Z axis : mean viewing direction, made orthogonal to the baseline
    const cPt3dr aMeanK = VUnit(mCam1.AxeK()+mCam2.AxeK());
    const cPt3dr aZ = aMeanK - aAxisX * Scal(aMeanK,aAxisX);
    MMVII_INTERNAL_ASSERT_User(SqN2(aZ)>1e-12,eTyUEr::eUnClassedError,"Epipolar rectification: viewing direction parallel to the baseline");
    const cPt3dr aAxisZ = VUnit(aZ);
    const cPt3dr aAxisY = aAxisZ ^ aAxisX;   // right-handed (X,Y,Z), already unit

    const cRotation3D<tREAL8> aRot(aAxisX,aAxisY,aAxisZ,false);
    const tREAL8 aFocal = 0.5 * (mCam1.InternalCalib()->F() + mCam2.InternalCalib()->F());

    // Z interval of each image (user, tie points or sensor), information for the crop of the second image
    std::optional<cPt2dr> aCache;
    const cPt2dr aZ1 = EpipEffectiveZInterval(mParams,mCam1,mCam1,mCam2,aCache);
    const cPt2dr aZ2 = EpipEffectiveZInterval(mParams,mCam2,mCam1,mCam2,aCache);

    cEpipModelPC aModel(std::make_unique<cEpipMappingPC>(mCam1,aRot,aFocal,aZ1),std::make_unique<cEpipMappingPC>(mCam2,aRot,aFocal,aZ2));
    aModel.ComputeCommonFraming(mCam1.PixelDomain().Box(),mCam2.PixelDomain().Box(),mParams.mEpipFrm,mParams.mMargin);
    return aModel;
}

// ============================================================
//  BenchEpipolarPC : synthetic pair with distortion
//  Same ground point => same epipolar row in both images ; mapping round trip ; virtual camera == mapping.
// ============================================================

void BenchEpipolarPCBody(cParamExeBench & aParam)
{

    for (int aKTest=0 ; aKTest<6 ; aKTest++)
    {
        // Same calibration for both cameras (distortion degree depends on the test), second camera shifted and slightly rotated
        auto * aCalib = cPerspCamIntrCalib::RandomCalib(eProjPC::eStenope,aKTest%4);
        const cPt2di aSz = aCalib->SzPix();
        const tREAL8 aDepth = 10.0;
        const cIsometry3D<tREAL8> aPose1(cPt3dr(0,0,0),cRotation3D<tREAL8>::RandomSmallElem(0.05));
        const cIsometry3D<tREAL8> aPose2(cPt3dr(aDepth*0.1,aDepth*0.02*RandUnif_C(),aDepth*0.02*RandUnif_C()),cRotation3D<tREAL8>::RandomSmallElem(0.1));
        cSensorCamPC aCam1("BenchEpip1",aPose1,aCalib);
        cSensorCamPC aCam2("BenchEpip2",aPose2,aCalib);

        cEpipolarRectificationPC::cParams aParams;
        aParams.mZIntv = cPt2dr(0.5*aDepth,1.5*aDepth);
        cEpipolarRectificationPC aRectif(aCam1,aCam2,aParams);
        auto aModel = aRectif.Compute();
        const cEpipolarMapping & aMap1 = aModel.EpipMap1();
        const cEpipolarMapping & aMap2 = aModel.EpipMap2();
        std::unique_ptr<cSensorCamPC> aEpipCam1(static_cast<const cEpipMappingPC&>(aMap1).VirtualCamera("EpipBench1"));
        std::unique_ptr<cSensorCamPC> aEpipCam2(static_cast<const cEpipMappingPC&>(aMap2).VirtualCamera("EpipBench2"));

        MMVII_INTERNAL_ASSERT_bench((aMap1.ZInterval()==aParams.mZIntv.value()) && (aMap2.ZInterval()==aParams.mZIntv.value()),"EpipolarPC : Z interval");
        MMVII_INTERNAL_ASSERT_bench((aMap1.EpipImSz().x()>0) && (aMap1.EpipImSz().y()>0),"EpipolarPC : empty frame");
        MMVII_INTERNAL_ASSERT_bench(aMap1.EpipFrame().P0().y()==aMap2.EpipFrame().P0().y(),"EpipolarPC : frames not aligned on rows");

        int aNbTested = 0;
        for (int aKPt=0 ; aKPt<100 ; aKPt++)
        {
            // ground point seen by the central part of image 1 (distortion inverse is approximate near the border)
            const cPt2dr aP1(aSz.x()*(0.2+0.6*RandUnif_0_1()),aSz.y()*(0.2+0.6*RandUnif_0_1()));
            const cPt3dr aPG = aCam1.ImageAndDepth2Ground(cPt3dr(aP1.x(),aP1.y(),aDepth*(0.7+0.6*RandUnif_0_1())));
            const cPt2dr aP2 = aCam2.Ground2Image(aPG);
            if ((aP2.x()<0.1*aSz.x()) || (aP2.x()>0.9*aSz.x()) || (aP2.y()<0.1*aSz.y()) || (aP2.y()>0.9*aSz.y()))
                continue;
            aNbTested++;

            const cPt2dr aE1 = aMap1.Value(aP1);
            const cPt2dr aE2 = aMap2.Value(aP2);
            MMVII_INTERNAL_ASSERT_bench(std::abs(aE1.y()-aE2.y())<1e-3,"EpipolarPC : rows differ");
            MMVII_INTERNAL_ASSERT_bench(Norm2(aMap1.Inverse(aE1)-aP1)<1e-3,"EpipolarPC : mapping 1 round trip");
            MMVII_INTERNAL_ASSERT_bench(Norm2(aMap2.Inverse(aE2)-aP2)<1e-3,"EpipolarPC : mapping 2 round trip");
            MMVII_INTERNAL_ASSERT_bench(Norm2(aEpipCam1->Ground2Image(aPG)-aE1)<1e-3,"EpipolarPC : virtual camera 1 differs from the mapping");
            MMVII_INTERNAL_ASSERT_bench(Norm2(aEpipCam2->Ground2Image(aPG)-aE2)<1e-3,"EpipolarPC : virtual camera 2 differs from the mapping");
            MMVII_INTERNAL_ASSERT_bench(std::abs(aEpipCam1->Ground2Image(aPG).y()-aEpipCam2->Ground2Image(aPG).y())<1e-3,"EpipolarPC : virtual cameras rows differ");
        }
        MMVII_INTERNAL_ASSERT_bench(aNbTested>30,"EpipolarPC : too few points tested");

        // Pair model on disk : mappings rebuilt from the file, bound to their cameras, give the same values
        const std::string aModelFile = cMMVII_Appli::CurrentAppli().TmpDirTestMMVII() + "EpipPCModel." + GlobTaggedNameDefSerial();
        cEpipPairModel(aMap1.Clone(),aMap2.Clone(),"OriTest","BenchEpip1","BenchEpip2").ToFile(aModelFile);
        auto aReloaded = cEpipPairModel::FromFile(aModelFile);
        RemoveFile(aModelFile,false);
        MMVII_INTERNAL_ASSERT_bench((aReloaded.Map(1).TypeName()=="PC") && (aReloaded.Map(2).TypeName()=="PC"),"EpipolarPC : model type");
        aReloaded.BindSensors(aCam1,aCam2);
        for (int aKPt=0 ; aKPt<20 ; aKPt++)
        {
            const cPt2dr aP(aSz.x()*(0.1+0.8*RandUnif_0_1()),aSz.y()*(0.1+0.8*RandUnif_0_1()));
            MMVII_INTERNAL_ASSERT_bench(Norm2(aReloaded.Map(1).Value(aP)-aMap1.Value(aP))<1e-6,"EpipolarPC : model round trip (Value)");
            MMVII_INTERNAL_ASSERT_bench(Norm2(aReloaded.Map(2).Inverse(aMap2.Value(aP))-aP)<1e-3,"EpipolarPC : model round trip (Inverse)");
        }
        MMVII_INTERNAL_ASSERT_bench(aReloaded.Map(1).EpipFrame().P0()==aMap1.EpipFrame().P0(),"EpipolarPC : model round trip (frame)");

        // EpipSensor gives the virtual camera
        std::unique_ptr<cSensorImage> aSI(EpipSensor(aCam1,aMap1,"EpipBenchSensor"));
        MMVII_INTERNAL_ASSERT_bench(aSI->GetSensorCamPC()!=nullptr,"EpipolarPC : EpipSensor is not a camera");
        delete aCalib;
    }

}

// ============================================================
//  BenchEpipolarPCCmd : EpipRectification then EpipResampling (crop) on real central perspective cameras.
//  The virtual cameras written by the commands must satisfy the epipolar property ; a crop is a shift of the full frame.
// ============================================================

// Camera read from a file, deleted at end of application (as the project does)
static cSensorCamPC * ReadCamAutoDel(const std::string & aName)
{
    cSensorCamPC * aCam = cSensorCamPC::FromFile(aName);
    cMMVII_Appli::AddObj2DelAtEnd(aCam);
    return aCam;
}

void BenchEpipolarPCCmd(cParamExeBench & aParam)
{
    if (! aParam.NewBench("EpipolarPCCmd")) return;

    cMMVII_Appli & anAp = cMMVII_Appli::CurrentAppli();
    const std::string aInDir = anAp.InputDirTestMMVII() + "/Saisies-MMV1/";
    const std::string aProj = anAp.TmpDirTestMMVII() + "EpipPCProj/";
    const std::string aOriDir = aProj + "MMVII-PhgrProj/Ori/toto/";
    const std::string aOutDir = aProj + "Out/";
    RemoveRecurs(aProj,false,SVP::Yes);
    CreateDirectories(aOriDir,SVP::No);
    const std::string aNameIm1 = "IMGP4167.JPG", aNameIm2 = "IMGP4168.JPG";
    for (const std::string & aFile : {std::string("Calib-PerspCentral-Foc-28000_Cam-PENTAX_K5.xml"),"Ori-PerspCentral-" + aNameIm1 + ".xml","Ori-PerspCentral-" + aNameIm2 + ".xml"})
        CopyFile(aInDir + "MMVII-PhgrProj/Ori/toto/" + aFile,aOriDir + aFile);

    // Rectification (closed form : two cameras), model and sensors only
    anAp.ExeCallMMVII
        (
            TheSpec_EpipRectification,
            anAp.StrObl() << aInDir+aNameIm1 << aInDir+aNameIm2 << "toto",
            anAp.StrOpt() << std::make_pair("DirProj",aProj) << std::make_pair("ZIntv","[1.4,1.6]") << std::make_pair("SaveModel","true")
                          << std::make_pair("OutDir",aOutDir) << std::make_pair("NoImage","true") << std::make_pair("StdOut","0"+aProj+"Rectification.txt")
        );

    const std::string aBase1 = "Epip_IMGP4167_IMGP4168.tif", aBase2 = "Epip_IMGP4168_IMGP4167.tif";
    const std::string aModelFile = aOutDir + "Epip_IMGP4167_IMGP4168.EpipModel." + GlobTaggedNameDefSerial();
    MMVII_INTERNAL_ASSERT_bench(ExistFile(aModelFile),"EpipolarPCCmd : no model file");
    auto aModel = cEpipPairModel::FromFile(aModelFile);
    MMVII_INTERNAL_ASSERT_bench((aModel.Map(1).TypeName()=="PC") && (aModel.Map(2).TypeName()=="PC"),"EpipolarPCCmd : model is not closed form");
    MMVII_INTERNAL_ASSERT_bench(aModel.Map(1).ZInterval()==cPt2dr(1.4,1.6),"EpipolarPCCmd : Z interval of the model");

    cSensorCamPC * aCam1 = ReadCamAutoDel(aOriDir + "Ori-PerspCentral-" + aNameIm1 + ".xml");
    cSensorCamPC * aCam2 = ReadCamAutoDel(aOriDir + "Ori-PerspCentral-" + aNameIm2 + ".xml");
    aModel.BindSensors(*aCam1,*aCam2);
    const std::string aOri1 = aOutDir + cSensorCamPC::NameOri_From_Image(aBase1), aOri2 = aOutDir + cSensorCamPC::NameOri_From_Image(aBase2);
    MMVII_INTERNAL_ASSERT_bench(ExistFile(aOri1) && ExistFile(aOri2),"EpipolarPCCmd : no epipolar orientation written");
    cSensorCamPC * aEpip1 = ReadCamAutoDel(aOri1);
    cSensorCamPC * aEpip2 = ReadCamAutoDel(aOri2);

    // Ground points seen by both images : same row in the epipolar orientations and the model, virtual camera == mapping
    int aNbTested = 0;
    const cPt2di aSz = aCam1->SzPix();
    for (int aKPt=0 ; aKPt<400 ; aKPt++)
    {
        const cPt2dr aP1(aSz.x()*(0.15+0.7*RandUnif_0_1()),aSz.y()*(0.15+0.7*RandUnif_0_1()));
        const cPt3dr aPG = aCam1->ImageAndDepth2Ground(cPt3dr(aP1.x(),aP1.y(),4.0+8.0*RandUnif_0_1()));
        const cPt2dr aP2 = aCam2->Ground2Image(aPG);
        if ((aCam2->DegreeVisibility(aPG)<=0) || (aCam2->PixelDomain().Insideness(aP2)<=0.1*aSz.y()))
            continue;
        aNbTested++;
        const cPt2dr aE1 = aEpip1->Ground2Image(aPG), aE2 = aEpip2->Ground2Image(aPG);
        MMVII_INTERNAL_ASSERT_bench(std::abs(aE1.y()-aE2.y())<1e-2,"EpipolarPCCmd : rows of the written orientations differ");
        MMVII_INTERNAL_ASSERT_bench(Norm2(aE1-aModel.Map(1).Value(aP1))<1e-2,"EpipolarPCCmd : orientation 1 differs from the model");
        MMVII_INTERNAL_ASSERT_bench(Norm2(aE2-aModel.Map(2).Value(aP2))<1e-2,"EpipolarPCCmd : orientation 2 differs from the model");
    }
    MMVII_INTERNAL_ASSERT_bench(aNbTested>20,"EpipolarPCCmd : too few points tested");

    // Resampling of a crop : sizes, and the cropped orientation is the full one shifted by the crop origin
    const cPt2di aCrop0(300,300), aCrop1(700,600);
    anAp.ExeCallMMVII
        (
            TheSpec_EpipResampling,
            anAp.StrObl() << aModelFile,
            anAp.StrOpt() << std::make_pair("DirProj",aProj) << std::make_pair("OutDir",aProj+"Crop/") << std::make_pair("StdOut","0"+aProj+"Resampling.txt")
                          << std::make_pair("CropP0",cStrIO<cPt2di>::ToStr(aCrop0)) << std::make_pair("CropP1",cStrIO<cPt2di>::ToStr(aCrop1))
        );
    MMVII_INTERNAL_ASSERT_bench(cDataFileIm2D::Create(aProj+"Crop/"+aBase1,eForceGray::No).Sz()==aCrop1-aCrop0,"EpipolarPCCmd : size of the cropped image");
    cSensorCamPC * aCropCam1 = ReadCamAutoDel(aProj + "Crop/" + cSensorCamPC::NameOri_From_Image(aBase1));
    for (int aKPt=0 ; aKPt<20 ; aKPt++)
    {
        const cPt3dr aPG = aCam1->ImageAndDepth2Ground(cPt3dr(aSz.x()*RandUnif_0_1(),aSz.y()*RandUnif_0_1(),4.0+8.0*RandUnif_0_1()));
        MMVII_INTERNAL_ASSERT_bench(Norm2(aCropCam1->Ground2Image(aPG)-(aEpip1->Ground2Image(aPG)-ToR(aCrop0)))<1e-6,"EpipolarPCCmd : cropped orientation");
    }

    RemoveRecurs(aProj,false,SVP::Yes);
    aParam.EndBench();
}

}; // MMVII
