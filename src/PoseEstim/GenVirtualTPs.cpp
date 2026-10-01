/**
   \file GenVirtualTPs.cpp

   \brief Generation of virtual tie points from a 3D point cloud, shared by the
   pair and triplet relative-orientation commands

*/

#include "MMVII_PoseRel.h"
#include "MMVII_Geom3D.h"
#include "MMVII_enums.h"

namespace MMVII
{

bool GenerateVirtualPts
    (
        const std::vector<cPt3dr> & aVPts,    ///< 3D points
        const std::vector<tREAL8> & aVW,      ///< weights, same size as aVPts
        const std::vector<tPoseR> & aVPoses,  ///< camera poses
        const std::vector<const cPerspCamIntrCalib*> & aVCalib,  ///< internal calibrations
        eDistrVirTPs                       aDistrib, ///< distribution (5, 9, 27 pts, ...)
        std::vector<std::vector<cPt2dr>> & aVVProj   ///< result : [virtual pt][camera]
        )
{
    MMVII_INTERNAL_ASSERT_tiny(aVPts.size()==aVW.size(),"GenerateVirtualPts: points/weights size mismatch");
    MMVII_INTERNAL_ASSERT_tiny(aVPoses.size()==aVCalib.size(),"GenerateVirtualPts: poses/calibs size mismatch");
    MMVII_INTERNAL_ASSERT_tiny(!aVPoses.empty(),"GenerateVirtualPts: no camera");

    const size_t aNbCam       = aVPoses.size();
    const int    NbSplits     = 10;   ///< number of bisection steps on the scale
    const cPt2dr aSz          = ToR(aVCalib.at(0)->SzPix());
    const tREAL8 TheMarginPix = 0.05 * std::min(aSz.x(),aSz.y());

    aVVProj.clear();

    // not enough points
    if (aVPts.size()<10) return false;

    // ignore points too far from the first camera
    std::vector<tREAL8> aVDist;
    for (const auto & aP : aVPts)
        aVDist.push_back(Norm2(aP-aVPoses.at(0).Tr()));
    std::vector<tREAL8> aVDistSort = aVDist;
    const tREAL8 aDMax = 5 * NC_KthVal(aVDistSort,0.75);

    // weighted covariance of the 3D points
    cStrStat2<tREAL8> aCovMat(3);
    for (size_t aK=0; aK<aVPts.size(); aK++)
    {
        if (aVDist.at(aK)<=aDMax)
            aCovMat.WeightedAdd(aVPts.at(aK).ToVect(),aVW.at(aK));
    }
    aCovMat.Normalise(true);
    cResulSymEigenValue<tREAL8> aVp = aCovMat.DoEigen();

    // degenerate cloud
    if (aVp.EigenValues()(0) <= 1e-12 * aVp.EigenValues()(2)) return false;

    cGenGauss3D aG3D(aVp.EigenVectors(),aVp.EigenValues(),aCovMat.Moy());

    // project the virtual points in all cameras; false if any falls outside the margin
    auto Project = [&](const std::vector<cPt3dr> & aVVirt,
                       std::vector<std::vector<cPt2dr>> & aVProj) -> bool
    {
        aVProj.clear();
        for (const auto & aP : aVVirt)
        {
            std::vector<cPt2dr> aVIm;
            for (size_t aKC=0; aKC<aNbCam; aKC++)
            {
                const cPt3dr aPL = aVPoses.at(aKC).Inverse(aP);
                if (aVCalib.at(aKC)->DegreeVisibility(aPL) <= TheMarginPix) return false;
                aVIm.push_back(aVCalib.at(aKC)->Value(aPL));
            }
            aVProj.push_back(aVIm);
        }
        return true;
    };

    std::vector<cPt3dr> aVVirtPts;
    std::vector<std::vector<cPt2dr>> aVProj;
    tREAL8 aLow  = 0.0;
    tREAL8 aHigh = 1.0;

    // generate virtual points
    aG3D.GetDistribNPts(aVVirtPts,aDistrib,1.0);

    // find the scale reprojecting all points in cameras:
    //     largest scale in ]0,1] with all points visible in all cameras
    if (Project(aVVirtPts,aVVProj))
        aLow = 1.0;
    else
    {
        for (int aK=0; aK<NbSplits; aK++)
        {
            tREAL8 aMid = (aLow+aHigh)/2.0;
            aG3D.GetDistribNPts(aVVirtPts,aDistrib,aMid);
            if (Project(aVVirtPts,aVProj))
            {
                aLow = aMid;
                aVVProj = aVProj;
            }
            else
                aHigh = aMid;
        }
    }

    // all points or none
    if (aLow <= 0.0)
    {
        aVVProj.clear();
        return false;
    }

    return true;
}

}; // namespace MMVII
