#include "MMVII_PlyExport.h"
#include "MMVII_Sensor.h"
#include "MMVII_Image2D.h"
#include "MMVII_PCSens.h"

namespace MMVII
{

/* ******************************************** */
/*                                              */
/*            cPlyExportCamGeom                 */
/*                                              */
/* ******************************************** */

cPlyExportCamGeom::cPlyExportCamGeom(tREAL8 aCamScale,int aNbSteps) :
    mCamScale(aCamScale),
    mNbSteps(aNbSteps)
{
    MMVII_INTERNAL_ASSERT_tiny(mNbSteps>0,"cPlyExportCamGeom : NbSteps must be >0");
}

void cPlyExportCamGeom::AddCenters(cPlyVertices & aPlyverts, const std::vector<cSensorImage *>& aVSens) const
{
    for (auto aCam : aVSens)
    {
        // pushbroom sensor
        /* if (aCam->CenterOfPC() == nullptr)
        {
            // Two bundles at the same scan-line Y but opposite X ends converge at the
            // perspective center of that line (unlike cross-track pairs which are near-parallel).
            double aY = aCam->Sz().y() / 2.0;
            tSeg3dr aBund0 = aCam->Image2Bundle(cPt2dr(0,              aY));
            tSeg3dr aBund1 = aCam->Image2Bundle(cPt2dr(aCam->Sz().x(), aY));
            aCenter = BundleInters(aBund0, aBund1);
            aCenter = aCam->PseudoCenterOfProj();
        } // perspective sensor
        else aCenter = aCam->PseudoCenterOfProj(); */

        cPt3dr aCenter = aCam->PseudoCenterOfProj();
        aPlyverts.AddVert(aCenter, {1.,0.,0.});

    }
}

void cPlyExportCamGeom::AddImagePlanes(cPlyVertices & aPlyverts, const std::vector<cSensorImage *>& aVSens) const
{
    // add points in image plane
    int aSteps = 50;
    for (auto aSens : aVSens)
    {
        cPt2di aSz = aSens->Sz();
        cPt2dr aImStepSz(aSz[0]/aSteps,aSz[1]/aSteps);
        // pushbroom sensor
        if (aSens->CenterOfPC() == nullptr)
        {
            // Image plane at scene level: use midpoint of each bundle (within valid RPC altitude range)
            auto BundMid = [&](cPt2dr aPt) {
                tSeg3dr aB = aSens->Image2Bundle(aPt);
                return (aB.P1() + aB.P2()) * 0.5;
            };

            for (int aS=0; aS<=aSteps; aS++)
            {
                std::vector<cPt3dr> aImVPts;
                double aDX = aImStepSz[0]*aS;
                double aDY = aImStepSz[1]*aS;

                aImVPts.push_back(BundMid(cPt2dr(0,       aDY   )));
                aImVPts.push_back(BundMid(cPt2dr(aDX,     0     )));
                aImVPts.push_back(BundMid(cPt2dr(aSz[0],  aDY   )));
                aImVPts.push_back(BundMid(cPt2dr(aDX,     aSz[1])));

                for (auto aP : aImVPts)
                    aPlyverts.AddVert(aP, {1.,0.,0.});
            }

            // Satellite trajectory: perspective centers sampled along the along-track direction
            for (int aS=0; aS<=aSteps; aS++)
            {
                double aDY = aImStepSz[1]*aS;
                tSeg3dr aBund0 = aSens->Image2Bundle(cPt2dr(0,       aDY));
                tSeg3dr aBund1 = aSens->Image2Bundle(cPt2dr(aSz[0],  aDY));
                cPt3dr aSatPos = BundleInters(aBund0, aBund1);

                aPlyverts.AddVert(aSatPos, {1.,0.,0.});
            }
        }
        else
        {
            cSensorCamPC *  aCamPC = aSens->GetSensorCamPC();

            double aFPix = aCamPC->InternalCalib()->F();
            double aF = CalculateFDepth(aSz,aFPix);


            // add image border
            std::vector<cPt3dr> aImVPts;
            for (int aS=0; aS<=aSteps; aS++)
            {

                double aDX = aImStepSz[0]*aS-1;
                double aDY = aImStepSz[1]*aS-1;

                // image plane border
                aImVPts.push_back(aSens->ImageAndDepth2Ground( cPt3dr(0,aDY,aF) ));
                aImVPts.push_back(aSens->ImageAndDepth2Ground( cPt3dr(aDX,0,aF) ));
                aImVPts.push_back(aSens->ImageAndDepth2Ground( cPt3dr(aSz[0],aDY,aF) ));
                aImVPts.push_back(aSens->ImageAndDepth2Ground( cPt3dr(aDX,aSz[1],aF) ));

            }
            // add image plane grid for non-pinhole
            if (aSens->IsSensorCamPC() )
            {
                if (aSens->GetSensorCamPC()->InternalCalib()->TypeProj()!=eProjPC::eStenope)
                {
                    int aSmallSteps = aSteps/2;
                    cPt2dr aImSmallStepSz(aSz[0]/aSmallSteps,aSz[1]/aSmallSteps);

                    for (int aSX=0; aSX<=aSmallSteps; aSX++)
                    {
                        double aDX = aImSmallStepSz[0]*aSX-1;

                        for (int aSY=0; aSY<=aSmallSteps; aSY++)
                        {
                            double aDY = aImSmallStepSz[1]*aSY-1;

                            aImVPts.push_back(aSens->ImageAndDepth2Ground( cPt3dr(aDX,aDY,aF) ));
                            //StdOut() << cPt2dr(aDX,aDY) << std::endl;
                        }
                    }
                }
            }

            for (auto aP : aImVPts)
                aPlyverts.AddVert(aP, {1.,0.,0.});
        }
    }
}

void cPlyExportCamGeom::AddFrustums(cPlyVertices & aPlyverts, const std::vector<cSensorImage *>& aVSens) const
{
    // add points on the frustum
    int aSteps = 50;
    for (auto aSens : aVSens)
    {
        // ignore for pushbroom camera
        if (aSens->CenterOfPC() == nullptr) continue;
        else
        {
            cSensorCamPC *  aCamPC = aSens->GetSensorCamPC();

            cPt2di aSz = aSens->Sz();
            double aFPix = aCamPC->InternalCalib()->F();
            double aF = CalculateFDepth(aSz,aFPix);

            double aFStepSz = aF/aSteps;


            for (int aS=0; aS<=aSteps; aS++)
            {
                std::vector<cPt3dr> aImVPts;
                double aDepth = aFStepSz*aS;
                aImVPts.push_back(aSens->ImageAndDepth2Ground( cPt3dr(0,0,aDepth) ));
                aImVPts.push_back(aSens->ImageAndDepth2Ground( cPt3dr(0,aSz[1],aDepth) ));
                aImVPts.push_back(aSens->ImageAndDepth2Ground( cPt3dr(aSz[0],0,aDepth) ));
                aImVPts.push_back(aSens->ImageAndDepth2Ground( cPt3dr(aSz[0],aSz[1],aDepth) ));


                for (auto aP : aImVPts)
                {
                    aPlyverts.AddVert(aP, {0.,1.,0.});
                }
            }
        }
    }
}

tREAL8 cPlyExportCamGeom::CalculateFDepth(const cPt2di & aSz,tREAL8 aF) const
{
    //    double aDiag = std::sqrt(std::pow(aSz[0],2)+std::pow(aSz[1],2));
    double aDiag = Norm2(aSz) ;  // MPDER : simpler
    double aRatioDiagF = aDiag/aF ;

    return aRatioDiagF*mCamScale;
}

void cPlyExportCamGeom::AddCameras(cPlyVertices& aPlyVerts, const std::vector<cSensorImage *>& aVSens) const
{
    AddCenters(aPlyVerts,aVSens);
    AddImagePlanes(aPlyVerts,aVSens);
    AddFrustums(aPlyVerts,aVSens);

}

/* ******************************************** */
/*                                              */
/*            cPlyExportTiePoints               */
/*                                              */
/* ******************************************** */

cPlyExportTiePoints::cPlyExportTiePoints(tREAL8 aErrProjMax,bool aWithRGB,bool aWithAvgRGB) :
    mErrProjMax(aErrProjMax),
    mWithRGB(aWithRGB),
    mWithAvgRGB(aWithAvgRGB)
{}

void cPlyExportTiePoints::AddTiePoints(cPlyVertices &                      aPlyVerts,
                                       cComputeMergeMulTieP &              aTPts,
                                       const std::vector<cSensorImage*> &  aVSens) const
{

    // Three-step processing for better memory management
    // 1- collect all 3D points and valid camera observations
    // 2- load each image once, collect colors, the free image
    // 3- save colored points to ply

    size_t aNbCam=aVSens.size();

    // Collect 3D and 2D points
    struct tObs { size_t mPtIdx; cPt2dr mPIm; }; //global point index, image obs
    std::vector<cPt3dr> aAllPts3D;
    std::vector<std::vector<tObs>> aObsByCam(aNbCam);

    for (auto& aAllConfigs : aTPts.Pts())
    {
        const auto & aConfig = aAllConfigs.first;
        auto & aVals = aAllConfigs.second;
        size_t aNbIm = aConfig.size();
        size_t aNbPts = aVals.mVIdPts.empty() ? NbPtsMul(aAllConfigs) : aVals.mVIdPts.size();

        for (size_t aKPts=0; aKPts<aNbPts; aKPts++)
        {
            std::vector<std::pair<size_t,cPt2dr>> aValidPtObs;

            const cPt3dr & aP3D = aVals.mVPGround.at(aKPts);

            for (size_t aKIm=0; aKIm<aNbIm; aKIm++)
            {
                size_t aKImSorted = aConfig.at(aKIm);
                const cPt2dr aPIm = aVals.mVPIm.at(aKPts*aNbIm+aKIm);
                cSensorImage* aCam = aVSens.at(aKImSorted);

                if (aCam->IsVisibleOnImFrame(aPIm) && aCam->IsVisible(aP3D))
                {
                    double aResidual = Norm2(aPIm - aCam->Ground2Image(aP3D));
                    if (aResidual<mErrProjMax)
                        aValidPtObs.push_back({aKImSorted,aPIm});

                }
            }

            if (aValidPtObs.size()>1)
            {
                size_t aGPtIdx = aAllPts3D.size(); // global index follows the order of pts in aAllPts3D
                aAllPts3D.push_back(aP3D);
                for (auto & [aCamIdx,aPIm] : aValidPtObs)
                {
                    aObsByCam.at(aCamIdx).push_back({aGPtIdx,aPIm});
                    if (!mWithAvgRGB) break; //no averaging so no need to collect more RGB vals
                }
            }
        }
    }

    // read an image at a time, accumulate colors, free the image
    std::vector<cPt3dr> aSumRGB(aAllPts3D.size(),{1.,1.,1.});
    std::vector<int>    aNbRGB(aAllPts3D.size(),0);

    if (mWithRGB)
    {
        // RGB values averaged over all images
        // neater results but slower?
        if (mWithAvgRGB)
        {
            for (size_t aKCam=0; aKCam<aNbCam; aKCam++)
            {
                StdOut() << "(" << aKCam+1 << "/" << aNbCam << ") "
                         << aVSens[aKCam]->NameImage() << (aObsByCam[aKCam].empty() ? " : no features" : "") << std::endl;

                if (aObsByCam[aKCam].empty()) continue;

                cRGBImage aImRGB = cRGBImage::FromFile(aVSens[aKCam]->NameImage());

                for (auto& aObs : aObsByCam[aKCam])
                {
                    if (aImRGB.InsideBL(aObs.mPIm))
                    {
                        aSumRGB[aObs.mPtIdx] += ToR(aImRGB.GetRGBPixBL(aObs.mPIm));
                        aNbRGB[aObs.mPtIdx] += 1;
                    }
                }
            }
        }
        else
        {
            //cMemManager::SetActiveMemoryCount(false);
            //#pragma omp parallel for schedule(dynamic)
            for (size_t aKCam = 0; aKCam < aVSens.size(); aKCam++)
            {
                //StdOutLock::lock();
                StdOut() << "(" << aKCam+1 << "/" << aNbCam << ") "
                         << aVSens[aKCam]->NameImage() << (aObsByCam[aKCam].empty() ? " : no features" : "") << std::endl;
                //StdOutLock::unlock();

                if (aObsByCam[aKCam].empty()) continue;

                cRGBImage aIm = cRGBImage::FromFile(aVSens[aKCam]->NameImage());

                for (auto & aObs : aObsByCam[aKCam])
                {
                    if (aIm.InsideBL(aObs.mPIm))
                    {
                        aSumRGB[aObs.mPtIdx] = ToR(aIm.GetRGBPixBL(aObs.mPIm));
                        aNbRGB[aObs.mPtIdx] = 1;
                    }
                }
            }
            //cMemManager::SetActiveMemoryCount(true);
        }
    }

    // send points to the ply pointcloud
    cPt3dr aDefRGB (1.,1.,1.);
    for (size_t aK=0; aK<aAllPts3D.size(); aK++)
    {
        cPt3dr aRGBVal = (mWithRGB && aNbRGB[aK]>0) ? (aSumRGB[aK] / (255.0 * aNbRGB[aK])) : aDefRGB ;
        aPlyVerts.AddVert(aAllPts3D[aK],aRGBVal);
    }

}
}; // MMVII