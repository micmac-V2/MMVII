#include "MMVII_StaticLidar.h"

#include <functional>
#include <fstream>

#include "../Mesh/happly.h"
#include "E57SimpleReader.h"
#include "MMVII_TplGradImFilter.h"
#include "MMVII_Tpl_Images.h"
#include "MMVII_ImageInfoExtract.h"
#include "../SymbDerGen/Formulas_CentralProj.h"
#include "MMVII_Interpolators.h"
#include "MMVII_util.h"
#include "MMVII_QuadTree.h"


namespace MMVII
{

cPt3dr cart2spher(const cPt3dr & aPtCart)
{
    tREAL8 dist = Norm2(aPtCart);
    tREAL8 theta =  atan2(aPtCart.y(),aPtCart.x());
    tREAL8 distxy = Norm2(Proj(aPtCart));
    tREAL8 phi =  atan2(aPtCart.z(),distxy);
    return {theta, phi, dist};
}

cPt3dr spher2cart(const cPt3dr & aPtspher)
{
    tREAL8 dhz = aPtspher.z()*cos(aPtspher.y());
    tREAL8 x = dhz* cos(aPtspher.x());
    tREAL8 y = dhz* sin(aPtspher.x());
    tREAL8 z = aPtspher.z()*sin(aPtspher.y());
    return {x, y, z};
}

tREAL8 toMinusPiPlusPi(tREAL8 aAng, tREAL8 aOffset)
{
    if (!std::isfinite(aAng))
        return aAng;
    while (aAng-aOffset<-M_PI)  aAng += 2*M_PI;
    while (aAng-aOffset>M_PI)  aAng -= 2*M_PI;
    return aAng;
}

cStaticLidarImporter::cStaticLidarImporter() :
    mHasRowCol(false), mNoMiss(false), mAllPointsReturn(false), mIsStrucured(false),
    mReadPose(), mDistMinToExist(1e-5),
    mNbCol            (0),
    mNbLine           (0),
    mThetaStart        (NAN),
    mThetaStep         (NAN),
    mPhiStart          (NAN),
    mPhiStep           (NAN),
    mRotInput2TSL      (cRotation3D<tREAL8>::Identity()),
    mVertRot           (cRotation3D<tREAL8>::Identity())
{

}

using  namespace happly;

void cStaticLidarImporter::readPlyPoints(std::string aPlyFileName, bool aForceGreenAsIntensity)
{
    StdOut() << "Read ply file " << aPlyFileName << "..." << std::endl;
    mVectPtsXYZ.clear();
    mVectPtsTPD.clear();
    mVectPtsIntens.clear();
    mVectPtsCol.clear();
    mVectPtsLine.clear();
    try
    {
        PLYData  aPlyF(aPlyFileName,false);
        auto aElementsNames = aPlyF.getElementNames();
        // Read points
        {
            std::vector<std::array<double, 3>> aVecPts = aPlyF.getVertexPositions();
            mVectPtsXYZ.resize(aVecPts.size());
            for (size_t i=0; i<aVecPts.size(); ++i)
            {
                mVectPtsXYZ[i] = cPt3dr(aVecPts[i][0],aVecPts[i][1],aVecPts[i][2]);
            }
        }

        mHasCartesian = true;
        mHasIntensity = false;
        mHasSpherical = false;
        mHasRowCol = false;

        // try to fill points attribute with any "i", "I" or "*intens*" property
        auto aVertProps = aPlyF.getElement("vertex").getPropertyNames();
        auto aPropIntensityName = std::find_if(aVertProps.begin(), aVertProps.end(), [](const std::string &s){
            return (ToLower(s)=="i") || (ToLower(s).find("intens") != std::string::npos);
        });
        auto aPropGreenName = std::find_if(aVertProps.begin(), aVertProps.end(), [](const std::string &s){
            return (ToLower(s)=="g") || (ToLower(s).find("green") != std::string::npos);
        });
        if (aForceGreenAsIntensity)
        {
            if (aPropGreenName!= aVertProps.end())
            {
                mHasIntensity = true;
                mVectPtsIntens = aPlyF.getElement("vertex").getProperty<tREAL8>(*aPropGreenName);
                for (size_t i=0; i<mVectPtsIntens.size(); ++i)
                {
                    mVectPtsIntens[i] = mVectPtsIntens[i]/255.;
                }
            }
        } else {
            if (aPropIntensityName!= aVertProps.end())
            {
                mHasIntensity = true;
                mVectPtsIntens = aPlyF.getElement("vertex").getProperty<tREAL8>(*aPropIntensityName);
            }
        }

    }
    catch (const std::runtime_error &e)
    {
        MMVII_UserError(eTyUEr::eReadFile, std::string("Error reading PLY file \"") + aPlyFileName + "\": " + e.what());
    }
}

void cStaticLidarImporter::readE57Points(std::string aE57FileName, bool aForceGreenAsIntensity)
{
    StdOut() << "Read e57 file " << aE57FileName << "..." << std::endl;
    mVectPtsXYZ.clear();
    mVectPtsTPD.clear();
    mVectPtsIntens.clear();
    mVectPtsCol.clear();
    mVectPtsLine.clear();
    try
    {
        e57::Reader reader( aE57FileName, {});
        MMVII_INTERNAL_ASSERT_tiny(reader.IsOpen(), "Error: unable to open file " + aE57FileName)
        StdOut() << "Image2DCount: " << reader.GetImage2DCount() << "\n";
        StdOut() << "Data3DCount: " << reader.GetData3DCount() << "\n";
        MMVII_INTERNAL_ASSERT_tiny(reader.GetData3DCount()==1, "Error: File should have exactly 1 scan for now")
        e57::E57Root fileHeader;
        reader.GetE57Root( fileHeader );
        /*StdOut() << fileHeader.formatName << " =? " << "ASTM E57 3D Imaging Data File" << std::endl;
        StdOut() << fileHeader.versionMajor << " =? " << 1 << std::endl;
        StdOut() << fileHeader.versionMinor << " =? " << 0 << std::endl;
        StdOut() << fileHeader.guid << " =? " << "Zero Points GUID" << std::endl;*/
        e57::Data3D data3DHeader;
        reader.ReadData3D( 0, data3DHeader );
        // data3DHeader.indexBounds is not correct
        const uint64_t cNumPoints = data3DHeader.pointCount;
        e57::Data3DPointsFloat pointsData( data3DHeader );
        auto vectorReader = reader.SetUpData3DPointsData( 0, cNumPoints, pointsData );
        const uint64_t cNumRead = vectorReader.read();
        MMVII_INTERNAL_ASSERT_tiny(cNumPoints==cNumRead, "Error: cNumPoints!=cNumRead")

        mHasCartesian = pointsData.cartesianX && pointsData.cartesianY && pointsData.cartesianZ;
        mHasIntensity = aForceGreenAsIntensity ? (pointsData.colorGreen!=nullptr) : (pointsData.intensity!=nullptr);
        mHasSpherical = pointsData.sphericalAzimuth && pointsData.sphericalElevation && pointsData.sphericalRange;
        mHasRowCol = pointsData.columnIndex && pointsData.rowIndex;

        if (mHasCartesian){
            mVectPtsXYZ.resize(cNumRead);
            for (uint64_t i=0;i<cNumRead;++i)
                mVectPtsXYZ[i] = {pointsData.cartesianX[i], pointsData.cartesianY[i], pointsData.cartesianZ[i]};
        }
        if (mHasSpherical){
            mVectPtsTPD.resize(cNumRead);
            for (uint64_t i=0;i<cNumRead;++i)
                mVectPtsTPD[i] = {pointsData.sphericalAzimuth[i], pointsData.sphericalElevation[i], pointsData.sphericalRange[i]};
        }
        if (mHasIntensity){
            mVectPtsIntens.resize(cNumRead);
            if (aForceGreenAsIntensity)
            {
                for (uint64_t i=0;i<cNumRead;++i)
                    mVectPtsIntens[i] = pointsData.colorGreen[i]/255.;
            } else {
                for (uint64_t i=0;i<cNumRead;++i)
                    mVectPtsIntens[i] = pointsData.intensity[i];
            }
        }
        if (mHasRowCol){
            mVectPtsLine.resize(cNumRead);
            mVectPtsCol.resize(cNumRead);
            mNbLine = 0;
            for (uint64_t i=0;i<cNumRead;++i)
            {
                mVectPtsLine[i] = pointsData.rowIndex[i];
                if (pointsData.rowIndex[i] >= mNbLine)
                    mNbLine = pointsData.rowIndex[i] + 1;
            }
            mNbCol = 0;
            for (uint64_t i=0;i<cNumRead;++i)
            {
                mVectPtsCol[i] = pointsData.columnIndex[i];
                if (pointsData.columnIndex[i] >= mNbCol)
                    mNbCol = pointsData.columnIndex[i] + 1;
            }
        }

        vectorReader.close();
        cRotation3D<tREAL8> aRotTSL2MM = cRotation3D<tREAL8>::RotFromCanonicalAxes("k-i-j");
        cRotation3D<tREAL8> aRot = cRotation3D(Quat2MatrRot<tREAL8>({data3DHeader.pose.rotation.w,data3DHeader.pose.rotation.x,
                                                 data3DHeader.pose.rotation.y,data3DHeader.pose.rotation.z}).Transpose(), false);

        mReadPose.emplace(cPt3dr(data3DHeader.pose.translation.x,data3DHeader.pose.translation.y,data3DHeader.pose.translation.z),
                          (aRotTSL2MM*aRot).MapInverse());
    }
    catch (const std::runtime_error &e)
    {
        MMVII_UserError(eTyUEr::eReadFile, std::string("Error reading E57 file \"") + aE57FileName + "\": " + e.what());
    }
}


void cStaticLidarImporter::readPtxPoints(std::string aPtxFileName, bool aForceGreenAsIntensity)
{
    StdOut() << "Read PTX file " << aPtxFileName << "..." << std::endl;
    mVectPtsXYZ.clear();
    mVectPtsTPD.clear();
    mVectPtsIntens.clear();
    mVectPtsCol.clear();
    mVectPtsLine.clear();
    try
    {
        std::ifstream  aPtxFile(aPtxFileName);
        if(!aPtxFile.is_open())
            throw std::runtime_error("Unable to open file "+aPtxFileName);
        mHasCartesian = true;
        mHasIntensity = true; // to check
        mHasSpherical = false;
        mHasRowCol = true; // directly computed
        mIsStrucured = true;
        mNoMiss = false;
        mAllPointsReturn = false;

        aPtxFile >> mNbCol;
        aPtxFile >> mNbLine;
        tREAL8 aTx, aTy, aTz;
        aPtxFile >> aTx >> aTy >> aTz;
        tREAL8 aR11, aR12, aR13;
        aPtxFile >> aR11 >> aR12 >> aR13;
        tREAL8 aR21, aR22, aR23;
        aPtxFile >> aR21 >> aR22 >> aR23;
        tREAL8 aR31, aR32, aR33;
        aPtxFile >> aR31 >> aR32 >> aR33;
        mReadPose.emplace(cPt3dr(aTx, aTy, aTz),
                          cRotation3D<tREAL8>({aR11, aR12, aR13}, {aR21, aR22, aR23}, {aR31, aR32, aR33}, false));
        char tmp[200];
        aPtxFile.getline(tmp, 200); // for now just skip transformation matrix
        aPtxFile.getline(tmp, 200);
        aPtxFile.getline(tmp, 200);
        aPtxFile.getline(tmp, 200);

        mVectPtsXYZ.resize(mNbCol*mNbLine);
        mVectPtsCol.resize(mNbCol*mNbLine);
        mVectPtsLine.resize(mNbCol*mNbLine);

        long i =0 ;
        tREAL8 aX, aY, aZ, aI;
        // for 1st line we test if there is intensity
        aPtxFile.getline(tmp, 200);
        std::istringstream iss(tmp);
        iss >> aX >> aY >> aZ >> aI;
        mVectPtsXYZ[i] = {aX, aY, aZ};
        if (iss.bad())
        {
            mHasIntensity = false;
        } else {
            mVectPtsIntens.resize(mNbCol*mNbLine);
            mVectPtsIntens[i] = aI;
        }

        for (long aCol = 0; aCol < mNbCol; aCol++)
        {
            for (long aRow = 0; aRow < mNbLine; aRow++)
            {
                ++i;
                aPtxFile.getline(tmp, 200);
                std::istringstream iss(tmp);
                iss >> aX >> aY >> aZ;
                mVectPtsXYZ[i] = {aX, aY, aZ};
                mVectPtsCol[i] = aCol;
                mVectPtsLine[i] = aRow;
                if (mHasIntensity)
                {
                    iss >> aI;
                    mVectPtsIntens[i] = aI;
                }
            }
        }


    }
    catch (const std::runtime_error &e)
    {
        MMVII_UserError(eTyUEr::eReadFile, std::string("Error reading PTX file \"") + aPtxFileName + "\": " + e.what());
    }
}

bool cStaticLidarImporter::read(const std::string & aName, bool OkNone,
                                bool aForceStructured, std::string aStrInput2TSL, bool aForceGreenAsIntensity)
{
    std::string aPost = LastPostfix(aName);
    if (UCaseEqual(aPost,"ply"))
       readPlyPoints(aName, aForceGreenAsIntensity);
    else if (UCaseEqual(aPost,"e57"))
       readE57Points(aName, aForceGreenAsIntensity);
    else if (UCaseEqual(aPost,"ptx"))
        readPtxPoints(aName, aForceGreenAsIntensity);
    else
    {
        if (! OkNone)
        {
           MMVII_UnclasseUsEr("Cannot read cloud for " + aName);
        }
        return false;
    }
    if (aForceStructured)
        mIsStrucured = true;

    StdOut() << "Read data: ";
    if (HasCartesian())
        StdOut() << mVectPtsXYZ.size() << " cartesian points";
    if (HasSpherical())
        StdOut() << mVectPtsTPD.size() << " spherical points";
    if (HasIntensity())
        StdOut() << " with intensity";
    if (HasRowCol())
        StdOut() << " with row-col";
    StdOut() << "\n";

    if (HasCartesian() && !HasSpherical())
    {
        mRotInput2TSL = cRotation3D<tREAL8>::RotFromCanonicalAxes(aStrInput2TSL);
        // apply rotframe to original points
        for (auto & aPtXYZ : mVectPtsXYZ)
        {
            aPtXYZ = mRotInput2TSL.Value(aPtXYZ);
        }
        convertToThetaPhiDist();
        // go back to original xyz
        for (auto & aPtXYZ : mVectPtsXYZ)
        {
            aPtXYZ = mRotInput2TSL.Inverse(aPtXYZ);
        }
    } else if (!HasCartesian() && HasSpherical()) // mTransfoIJK not used if spherical
    {
        convertToXYZ();
    }


    return true;
}



void cStaticLidarImporter::convertToThetaPhiDist()
{
    StdOut() << "convertToThetaPhiDist\n";
    mNoMiss = false;
    mAllPointsReturn = false;
    mVectPtsTPD.resize(mVectPtsXYZ.size()); // all points in theta-phi-dist
    size_t aNbPtsNul = 0;
    for (size_t i=0; i<mVectPtsTPD.size(); ++i)
    {
        mVectPtsTPD[i] = cart2spher(mVectPtsXYZ[i]);
        if (mVectPtsTPD[i].z()<mDistMinToExist)
            aNbPtsNul++;
    }
    StdOut() << aNbPtsNul << " null points\n";
    mNoMiss = mIsStrucured || (aNbPtsNul>0);
    mAllPointsReturn = (aNbPtsNul==0);
}


void cStaticLidarImporter::convertToXYZ()
{
    StdOut() << "convertToXYZ\n";
    mNoMiss = false;
    mAllPointsReturn = false;
    mVectPtsXYZ.resize(mVectPtsTPD.size());
    size_t aNbPtsNul = 0;
    for (size_t i=0; i<mVectPtsTPD.size(); ++i)
    {
        mVectPtsXYZ[i] = spher2cart(mVectPtsTPD[i]);
        if (mVectPtsTPD[i].z()<mDistMinToExist)
            aNbPtsNul++;
    }
    StdOut() << aNbPtsNul << " null points\n";
    mNoMiss = mIsStrucured || (aNbPtsNul>0);
    mAllPointsReturn = (aNbPtsNul==0);
}

void cStaticLidarImporter::ComputeRotInput2Raster(std::string aTransfoIJK)
{
    MMVII_INTERNAL_ASSERT_tiny(HasRowCol(),"Error: ComputeRotInput2Raster needs row/col");

    cRotation3D<tREAL8> aRotInput2TSL = cRotation3D<tREAL8>::RotFromCanonicalAxes(aTransfoIJK);
    //  scanner     to      raster
    //  z|   y
    //   | /
    //   + --- x            + --- z
    //                      | `x
    //                     y|
    cRotation3D<tREAL8> aRotZvert2view = cRotation3D<tREAL8>::RotFromCanonicalAxes("k-i-j");

    std::cout<<"aRotZvert2view:\n"<<aRotZvert2view.AxeI()<<"\n"
              <<aRotZvert2view.AxeJ()<<"\n"<<aRotZvert2view.AxeK()<<std::endl;
    std::cout<<"aRotFrame:\n"<<aRotInput2TSL.AxeI()<<"\n"
              <<aRotInput2TSL.AxeJ()<<"\n"<<aRotInput2TSL.AxeK()<<std::endl;

    mRotInput2Raster = aRotZvert2view * aRotInput2TSL; // TODO: use aRotFrame
    std::cout<<"mRotInput2Raster:\n"<<mRotInput2Raster.AxeI()<<"\n"
              <<mRotInput2Raster.AxeJ()<<"\n"<<mRotInput2Raster.AxeK()<<std::endl;
}


bool cStaticLidarImporter::checkLineCol()
{
    auto [aColMinIt, aColMaxIt] = std::minmax_element(mVectPtsCol.begin(), mVectPtsCol.end());
    auto [aLineMinIt, aLineMaxIt] = std::minmax_element(mVectPtsLine.begin(), mVectPtsLine.end());
    bool isOk = true;
    if ((*aColMinIt<0) || (*aColMaxIt>mNbCol-1))
        isOk = false;
    if ((*aLineMinIt<0) || (*aLineMaxIt>mNbLine-1))
        isOk = false;
    StdOut() << "checkLineCol: " << *aColMinIt << " " << *aColMaxIt
             << " " << *aLineMinIt << " " << *aLineMaxIt
             << " " << mNbCol << " " << mNbLine << "\n";
    MMVII_INTERNAL_ASSERT_tiny(isOk,"Error: checkLineCol");
    mNbCol = *aColMaxIt + 1;
    mNbLine = *aLineMaxIt + 1;
    return isOk;
}


float cStaticLidarImporter::ColToLocalThetaApprox(float aCol) const
{
    return mThetaStart + aCol * mThetaStep;
}

float cStaticLidarImporter::LineToLocalPhiApprox(float aLine) const
{
    return mPhiStart + aLine * mPhiStep;
}

float cStaticLidarImporter::LocalThetaToColApprox(float aTheta) const
{
    tREAL8 aNbCol2pi = 2 * M_PI / fabs(mThetaStep);
    tREAL8 aCol = (aTheta - mThetaStart) / mThetaStep;
    // try to return to [-aNbCol2pi/2 : aNbCol2pi/2]
    if (aCol<-aNbCol2pi/2)
        return aCol+aNbCol2pi;
    if (aCol>aNbCol2pi*3/2)
        return aCol-aNbCol2pi;
    return aCol;
}

float cStaticLidarImporter::LocalPhiToLineApprox(float aPhi) const
{
    return (aPhi - mPhiStart) / mPhiStep;
}

void cStaticLidarImporter::decimXY(const cPt2di & aDecimXY)
{
    if (aDecimXY==cPt2di(1,1))
        return;
    StdOut() << "decimXY\n";
    MMVII_INTERNAL_ASSERT_tiny(HasRowCol(),"Error: ComputeRotInput2Raster needs row/col");

    size_t aNewNbPts = mVectPtsTPD.size() / aDecimXY.x() / aDecimXY.y() + 1; // ???
    std::vector<int> aNewVectPtsLine;
    aNewVectPtsLine.reserve(aNewNbPts);
    std::vector<int> aNewVectPtsCol;
    aNewVectPtsCol.reserve(aNewNbPts);
    std::vector<cPt3dr> aNewVectPtsXYZ;
    aNewVectPtsXYZ.reserve(aNewNbPts);
    std::vector<tREAL8> aNewVectPtsIntens;
    aNewVectPtsIntens.reserve(aNewNbPts);
    std::vector<cPt3dr> aNewVectPtsTPD;
    aNewVectPtsTPD.reserve(aNewNbPts);
    size_t j = 0;
    cPt2di aDecimMid(aDecimXY.x() / 2,aDecimXY.y() / 2);

    for (size_t i=0; i<mVectPtsTPD.size(); ++i)
    {
        if (((mVectPtsLine[i] % aDecimXY.y()) == aDecimMid.y())
            && ((mVectPtsCol[i] % aDecimXY.x()) == aDecimMid.x()))
        {
            aNewVectPtsLine.push_back(mVectPtsLine[i] / aDecimXY.y());
            aNewVectPtsCol.push_back(mVectPtsCol[i] / aDecimXY.x());
            aNewVectPtsXYZ.push_back(mVectPtsXYZ[i]);
            aNewVectPtsIntens.push_back(mVectPtsIntens[i]);
            aNewVectPtsTPD.push_back(mVectPtsTPD[i]);

            /*if (mVectPtsCol[i]==0)
            {
                std::cout<<i<<" "<<j<<" "<<aNewVectPtsLine[j]<<" "<<mVectPtsLine[i]<<" "<<aDecimXY.y()
                          <<" "<<aNewVectPtsCol[j]<<" "<<mVectPtsCol[i]<<" "<<aDecimXY.x()
                          <<aNewVectPtsXYZ[j]<<" "<<aNewVectPtsTPD[j]<<"\n";
            }*/
            ++j;
        }
    }
    StdOut() << "decim " <<  mVectPtsTPD.size() << " " << aNewNbPts << " " << j << "\n";
    StdOut() << "before " <<  mThetaStep << " " << mPhiStep << " " << "\n";
    StdOut() << "size " <<  mNbCol << " " << mNbLine << " " << "\n";
    //MMVII_INTERNAL_ASSERT_tiny(aNewNbPts==j,"Error in decimation");

    mThetaStep *= aDecimXY.x();
    mPhiStep *= aDecimXY.y();

    mNbCol = (mNbCol + aDecimXY.x() - 1) / aDecimXY.x(); // if original size is 6991, decim 2 => 3495, but last was 6990 => 3495
    mNbLine = (mNbLine + aDecimXY.y() - 1) / aDecimXY.y();

    std::swap(aNewVectPtsLine, mVectPtsLine);
    std::swap(aNewVectPtsCol, mVectPtsCol);
    std::swap(aNewVectPtsXYZ, mVectPtsXYZ);
    std::swap(aNewVectPtsIntens, mVectPtsIntens);
    std::swap(aNewVectPtsTPD, mVectPtsTPD);

    StdOut() << "after " <<  mThetaStep << " " << mPhiStep << " " << "\n";
    StdOut() << "size " <<  mNbCol << " " << mNbLine << " " << "\n";
    StdOut() << "=> " <<  mVectPtsTPD.size() << "\n";
}


void cStaticLidarImporter::ComputeAgregatedAngles()
{
    mVectPhisCol.resize(mNbLine, 0.);
    mVectThetasLine.resize(mNbCol, 0.);
    std::vector<int> aNbMesPhisCol(mNbLine,0);
    std::vector<int> aNbMesThetasLine(mNbCol,0);
    for (size_t i=0; i<mVectPtsTPD.size(); ++i)
    {
        if (mVectPtsTPD[i].z()<mDistMinToExist)
            continue;
        mVectPhisCol.at(mVectPtsLine[i]) += mVectPtsTPD[i].y();
        aNbMesPhisCol.at(mVectPtsLine[i])++;
        mVectThetasLine.at(mVectPtsCol[i]) += mVectPtsTPD[i].x();
        aNbMesThetasLine.at(mVectPtsCol[i])++;
    }
    // compute average (NAN for no data)
    for (int i=0; i<mNbLine; ++i)
    {
        mVectPhisCol[i] /= aNbMesPhisCol[i];
    }
    for (int i=0; i<mNbCol; ++i)
    {
        mVectThetasLine[i] /= aNbMesThetasLine[i];
    }
    // TODO: check that phis are constant among first and last lines

    // export mean thetas per col
    std::fstream file_thetas;
    file_thetas.open("thetas.txt", std::ios_base::out);
    for (auto aTheta:mVectThetasLine)
    {
        file_thetas << aTheta << "\n";
    }
    //file_thetas << "\n";
    file_thetas.close();
}

float cStaticLidarImporter::LocalPhiToLinePrecise(float aPhi) const
{
    MMVII_INTERNAL_ASSERT_tiny(!mVectPhisCol.empty(),"Error: run ComputeAgregatedAngles() before LocalPhiToLinePrecise()");

    float aLineApprox = LocalPhiToLineApprox(aPhi);
    if ((aLineApprox<0) || (aLineApprox>=mNbLine))
        return aLineApprox;
    for (int i=0; i<5;++i)
    {
        int aLineBefore = (int)aLineApprox;
        int aLineAfter = (int)aLineApprox + 1;
        if ((aLineBefore<0) || (aLineBefore>=mNbLine))
            return aLineApprox;
        if ((aLineAfter<0) || (aLineAfter>=mNbLine))
            return aLineApprox;
        float aPhiBefore = mVectPhisCol[aLineBefore];
        float aPhiAfter = mVectPhisCol[aLineAfter];
        //std::cout<<"iter "<<aLineApprox<<"\n";
        aLineApprox = aLineBefore + (aPhi-aPhiBefore)/(aPhiAfter-aPhiBefore)*(aLineAfter-aLineBefore);
    }
    return aLineApprox;
}

float cStaticLidarImporter::LocalThetaToColPrecise(float aTheta) const
{
    MMVII_INTERNAL_ASSERT_tiny(!mVectThetasLine.empty(),"Error: run ComputeAgregatedAngles() before LocalThetaToColPrecise()");

    float aColApprox = LocalThetaToColApprox(aTheta);
    if ((aColApprox<0) || (aColApprox>=mNbCol))
        return aColApprox;
    for (int i=0; i<5;++i)
    {
        int aColBefore = (int)aColApprox;
        int aColAfter = (int)aColApprox + 1;
        if ((aColBefore<0) || (aColBefore>=mNbCol))
            return aColApprox;
        if ((aColAfter<0) || (aColAfter>=mNbCol))
            return aColApprox;
        float aThetaBefore = mVectThetasLine[aColBefore];
        float aThetaAfter = mVectThetasLine[aColAfter];
        aColApprox = aColBefore + (aTheta-aThetaBefore)/(aThetaAfter-aThetaBefore)*(aColAfter-aColBefore);
    }
    return aColApprox;
}


cPt2dr cStaticLidarImporter::Input3DtoRasterAngle(const cPt3dr &aPt3DInput) const
{
    //std::cout<<"Input3DtoRasterAngle: "<<aPt3DInput<<"\n";
    cPt3dr aPt3DInputNorm = aPt3DInput/Norm2(aPt3DInput);
    cPt3dr aPt3DRaster = RotInput2Raster().Value(aPt3DInputNorm);
    cProj_EquiRect aProjEquiRect(M_PI);
    auto aThetaPhi = cPt2dr::FromStdVector(aProjEquiRect.Proj(aPt3DRaster.ToStdVector()));
    return aThetaPhi;
}
/*    cPt2dr aP2d_approx;
    aP2d_approx.x() = LocalThetaToColApprox(aThetaPhi.x());
    aP2d_approx.y() = LocalPhiToLineApprox(aThetaPhi.y());
    //TODO: iterative search around aP2d_approx in X/Z and Y/Z rasters
    cPt2dr aP2d = aP2d_approx;
    std::cout<<"  => "<<aThetaPhi<<"  => "<<aP2d<<"\n";
    return aP2d;
}*/

std::string cStaticLidarImporter::DefaultPoseName(const std::string & aDirStaticLidarRasters, const std::string & aLidarId)
{
    return aDirStaticLidarRasters + "PoseFromCloudFile-" + aLidarId + ".xml";
}

void cStaticLidarImporter::MakeIdImage(const std::string & aNameFile) const
{

    const int aTargetWidth = 800;
    int aFactor = NbCol()/aTargetWidth;
    if (aFactor<1)
        aFactor = 1;
    auto aIdIm = cIm2D<tU_INT1>(cPt2di((NbCol()+aFactor-1)/aFactor, (NbLine()+aFactor-1)/aFactor+1), 0, eModeInitImage::eMIA_Null);
    auto & aIdImData = aIdIm.DIm();

    tREAL8 aDistMin = INFINITY;
    tREAL8 aDistMax = 0.;
    tREAL8 aIntensMin = INFINITY;
    tREAL8 aIntensMax = 0.;
    for (size_t i=0; i<mVectPtsTPD.size(); ++i)
    {
        tREAL8 aDist = mVectPtsTPD[i].z();
        if (aDist==0)
            continue;

        if (aDist<aDistMin)
            aDistMin = aDist;
        if (aDist>aDistMax)
            aDistMax = aDist;

        tREAL8 aIntens = mVectPtsIntens[i];
        if (aIntens<aIntensMin)
            aIntensMin = aIntens;
        if (aIntens>aIntensMax)
            aIntensMax = aIntens;
    }
    //StdOut() <<"D: "<<aDistMin<<"-"<<aDistMax<<", I:"<<aIntensMin<<"-"<<aIntensMax<<"\n";

    for (size_t i=0; i<mVectPtsTPD.size(); ++i)
    {
        if (((mVectPtsLine[i] % aFactor) == 0)
            && ((mVectPtsCol[i] % aFactor) == 0))
        {
            cPt2di aPcl(mVectPtsCol[i],mVectPtsLine[i]);

            float aValDist = 127 * (1+cos(2*M_PI*10*(mVectPtsTPD[i].z()-aDistMin)/(aDistMax-aDistMin)));
            float aVal = HasIntensity() ? (mVectPtsIntens[i]-aIntensMin)/(aIntensMax-aIntensMin) :  aValDist;
            if (aVal<0)
                aVal =0;
            if (aVal>255)
                aVal = 255;
            aIdImData.SetV(aPcl/aFactor, aVal);
        }
    }

    aIdImData.ToFile(aNameFile);
}



//-----------------------------------------------------------------------

cStaticLidar::cStaticLidar(const std::string & aNameFile, const std::string & aStationName,
                           const std::string & aScanName, const tPose & aPose, cPerspCamIntrCalib * aCalib,
                           cRotation3D<tREAL8> aRotInput2Raster, tREAL8 aSigma) :
    cSensorCamPC(aNameFile, aPose, aCalib),
    mStationName(aStationName),
    mScanName(aScanName),
    mAreRastersReady(false),
    mSigma(aSigma),
    mRotInput2Raster(aRotInput2Raster),
    mTriangulation(nullptr),
    mEqDistColinearityDist(nullptr)
{

}

cStaticLidar::~cStaticLidar()
{
    if (mTriangulation)
        delete mTriangulation;
    for (auto & [aDistMin, aInterp] : mVInterpN)
        delete aInterp;
}

std::string cStaticLidar::GetIdSuffix()
{
    return ".tsl.gif";
}

std::string cStaticLidar::GetIdSuffixRegex()
{
    return "\\.tsl\\.gif";
}

bool cStaticLidar::DoAddCalibToUk() const
{
    return false; // F and PP fixed on lidar for now
}

std::tuple<double, double, cPt3dr> cStaticLidar::getDistSigmaNormalPlane(cPt2dr aCenter, const cPixBox<2> & aPixBox) const
{
    //static std::pair<cPlane3D,tREAL8> LSQEstimate(const std::vector<cPt3dr> & aP0,const std::vector<tREAL8>* =nullptr);
    std::vector<cPt3dr> aVPtCam;
    for (const auto & aPt : aPixBox)
    {
        if (IsValidPoint(aPt))
            aVPtCam.push_back(Image2Camera3D(aPt));
    }
    if (aVPtCam.empty() || (!IsValidPoint(aCenter)))
        return std::make_tuple(NAN,NAN,cPt3dr::Dummy());

    if (IsValidPoint(aCenter) && (aVPtCam.size()<3))
    {
        return std::make_tuple(Image2Distance(aCenter), Sigma(),
                               Pose().Rot().Value(Image2NormalInstr(aCenter))); // which interpolator to use for Target extract?
    }

    std::cout<<"Box "<<aPixBox<<" => "<<aPixBox.Sz().x()*aPixBox.Sz().y()
              <<" usable "<<aVPtCam.size()<<" "<<(100.*aVPtCam.size())/(aPixBox.Sz().x()*aPixBox.Sz().y())<<"%\n";
    auto [aPlane3D, aRes] = cPlane3D::LSQEstimate(aVPtCam,nullptr);
    auto anInter = aPlane3D.Inter({0,0,0}, Image2Camera3D(aCenter));
    auto aDist = Norm2(anInter);
    auto aNormalGnd = Pose().Rot().Value(aPlane3D.AxeK());

    return std::make_tuple(aDist, (aRes+Sigma())/2., aNormalGnd);
}

cDiffInterpolator1D * cStaticLidar::getNormalInterpolator(const cPt2dr & aPt) const
{
    if (!getRasterDistance().InsideBL(aPt))
    {
        return mVInterpN.front().second; // if impossible to get dist, return any interpolator
    }
    double aDist = Image2Distance(aPt);
    size_t i = 0;
    while ((i<mVInterpN.size()) && (aDist>mVInterpN[i].first))
        ++i;
    return mVInterpN.at(i-1).second;
}

void cStaticLidar::Show() const
{
    cSensorCamPC::Show();
    StdOut()<<"   Station: " << mStationName << std::endl;
    StdOut()<<"   Scan   : " << mScanName << std::endl;
}

void cStaticLidar::initInterpolators()
{
    // init interpolators
    double aNormalPrec = 5. * M_PI / 180.;
    std::vector< std::pair<double,double>> aV_Start_Mid_dist_class = {{0,1.25},{2,4},{6,12}}; // distances classes for normal interpolator (start and middle)
    //std::vector<std::string> aParamInterpolN = {"Scale","3","1000","Linear"};
    std::vector<std::string> aParamInterpolN = {"Scale","3","1000","Cubic","-0.5"};
    for (auto & [aStart,aMid]: aV_Start_Mid_dist_class)
    {
        int aSz = pow(mSigma/(tan(aNormalPrec)*aMid*1./InternalCalib()->F()),2./3.)+1;
        aParamInterpolN[1] = std::to_string(aSz);
        mVInterpN.emplace_back(aStart,cDiffInterpolator1D::AllocFromNames(aParamInterpolN));
    }
}

cStaticLidar * cStaticLidar::FromFile(const std::string & aNameScanOriFile, bool aSVP)
{
    if ( aSVP && (!ExistFile(aNameScanOriFile) || IsDirectory(aNameScanOriFile)) )
    {
        return nullptr;
    }

    cStaticLidar * aRes = new cStaticLidar("NONE","?","?",tPoseR::Identity(),nullptr,tRotR::Identity(),NAN);
    ReadFromFile(*aRes, aNameScanOriFile);

    cPerspCamIntrCalib* aCalib = cPerspCamIntrCalib::FromFile(DirOfPath(aNameScanOriFile)
                                                              + aRes->mTmpNameCalib + "." + GlobTaggedNameDefSerial());
    aRes->mInternalCalib = aCalib;

    aRes->initInterpolators();

    return aRes;
}

bool cStaticLidar::ReadRasters(const std::string & aDataDir)
{
    if (mAreRastersReady)
        return true;
    StdOut() << "Read rasters " << NameImage() <<"..."<<std::endl;
    mRasterDistance = std::make_unique<cIm2D<tREAL4>>(cIm2D<tREAL4>::FromFile(aDataDir+"/"+mRasterDistancePath));
    if (mRasterIntensityPath.size()>0)
        mRasterIntensity = std::make_unique<cIm2D<tU_INT1>>(cIm2D<tU_INT1>::FromFile(aDataDir+"/"+mRasterIntensityPath));
    mRasterMask = std::make_unique<cIm2D<tU_INT1>>(cIm2D<tU_INT1>::FromFile(aDataDir+"/"+mRasterMaskPath));
    mRasterX = std::make_unique<cIm2D<tREAL4>>(cIm2D<tREAL4>::FromFile(aDataDir+"/"+mRasterXPath));
    mRasterY = std::make_unique<cIm2D<tREAL4>>(cIm2D<tREAL4>::FromFile(aDataDir+"/"+mRasterYPath));
    mRasterZ = std::make_unique<cIm2D<tREAL4>>(cIm2D<tREAL4>::FromFile(aDataDir+"/"+mRasterZPath));
    mAreRastersReady = true;
    return true;
}

cPt3dr cStaticLidar::Image2InputXYZ(cPt2di aRasterPxI) const
{
    MMVII_INTERNAL_ASSERT_tiny(mAreRastersReady, "Error: rasters not ready");
    cPt2dr aRasterPx = ToR(aRasterPxI);
    InternalCalib()->FixLoopPixelsInImage(aRasterPx);
    aRasterPxI = ToI(aRasterPx);
    auto & aRasterXData = mRasterX->DIm();
    auto & aRasterYData = mRasterY->DIm();
    auto & aRasterZData = mRasterZ->DIm();
    if (!aRasterXData.Inside(aRasterPxI))
        return {0.,0.,0.};
    return cPt3dr{
        aRasterXData.GetV(aRasterPxI),
        aRasterYData.GetV(aRasterPxI),
        aRasterZData.GetV(aRasterPxI),
    };
}

cPt3dr cStaticLidar::Image2InputXYZ(cPt2dr aRasterPx) const
{
    MMVII_INTERNAL_ASSERT_tiny(mAreRastersReady, "Error: rasters not ready");
    InternalCalib()->FixLoopPixelsInImage(aRasterPx);
    auto & aRasterXData = mRasterX->DIm();
    auto & aRasterYData = mRasterY->DIm();
    auto & aRasterZData = mRasterZ->DIm();
    if (!aRasterXData.InsideBL(aRasterPx))
        return {0.,0.,0.};
    return cPt3dr{
        aRasterXData.GetVBL(aRasterPx),
        aRasterYData.GetVBL(aRasterPx),
        aRasterZData.GetVBL(aRasterPx),
    };
}

cPt3dr cStaticLidar::Image2InputXYZ_InterpoleDist(cPt2dr aRasterPx,
                                                  const cInterpolator1D *anInterpol) const
{
    cPt3dr aPtXYZ = Image2InputXYZ(aRasterPx);
    tREAL8 aDistFactor = Image2DistanceInterpol(aRasterPx, anInterpol) / Norm2(aPtXYZ);
    aPtXYZ = aPtXYZ * aDistFactor;
    return aPtXYZ;
}

cPt3dr cStaticLidar::Image2Ground(const cPt2di & aRasterPxI) const
{
    cPt2dr aRasterPx = ToR(aRasterPxI);
    return Image2Ground(aRasterPx);
}

cPt3dr cStaticLidar::Image2Ground(cPt2dr aRasterPx) const
{
    InternalCalib()->FixLoopPixelsInImage(aRasterPx);
    cPt3dr aCam3DPt = Image2Camera3D(aRasterPx);
    return Pose().Value(aCam3DPt);
}

cPt3dr cStaticLidar::Image2Ground_InterpoleDist(cPt2dr aRasterPx,
                                                const cInterpolator1D *anInterpol) const
{
    InternalCalib()->FixLoopPixelsInImage(aRasterPx);
    cPt3dr aCam3DPt = Image2Camera3D_InterpoleDist(aRasterPx, anInterpol);
    return Pose().Value(aCam3DPt);
}

tREAL4 cStaticLidar::Image2Distance(cPt2dr aRasterPx) const
{
    InternalCalib()->FixLoopPixelsInImage(aRasterPx);
    return getRasterDistance().GetVBL(aRasterPx);
}

tREAL4 cStaticLidar::Image2DistanceInterpol(cPt2dr aRasterPx, const cInterpolator1D *anInterpol) const
{
    InternalCalib()->FixLoopPixelsInImage(aRasterPx);
    return getRasterDistance().GetValueInterpol(*anInterpol, aRasterPx);
}

cPt3dr cStaticLidar::ImageAndDepth2Ground(const cPt3dr & aPIm3) const
{
    cPt3dr aCam3DPt = Image2Camera3D(Proj(aPIm3));
    auto aPtCamNorm = Norm2(aCam3DPt);
    if (aPtCamNorm>0)
    {
        // can use rasters
        cPt3dr aCam3DPtDist = aCam3DPt/aPtCamNorm * aPIm3.z();
        return Pose().Value(aCam3DPtDist);
    } else {
        return cSensorCamPC::ImageAndDepth2Ground(aPIm3);
    }
}

cPt3dr cStaticLidar::Image2Camera3D_InterpoleDist(const cPt2dr & aRasterPx,
                                                  const cInterpolator1D *anInterpol) const
{
    cPt3dr aPtInput3D = Image2InputXYZ_InterpoleDist(aRasterPx, anInterpol);
    cPt3dr aPtCam3D = mRotInput2Raster.Value(aPtInput3D);
    return aPtCam3D;
}


std::tuple<tREAL8,tREAL8,tREAL8> cStaticLidar::AvgDistNbValidAndNbNotMasked() const
{
    // take mean squared or cubed dist?
    tREAL16 aAvg = 0.;
    int aNb = 0;
    int aNbNotMasked = 0;
    for (int l = 0 ; l < PixelDomain().Sz().y(); l++)
        for (int c = 0 ; c < PixelDomain().Sz().x(); c++)
        {
            auto aDist = getRasterDistance().GetV(cPt2di(c, l));
            if (aDist>0)
            {
                aAvg+=aDist*aDist;//*aDist;
                ++aNb;
            }
            if (!IsMaskedPoint(cPt2dr(c, l)))
                ++aNbNotMasked;
        }
    //return {pow(aAvg/aNb,1./2.), aNb, aNbNotMasked};
    return {sqrt(aAvg/aNb), aNb, aNbNotMasked};
}


/// normal in sensor frame
cPt3dr cStaticLidar::Image2NormalInstr(const cPt2dr & aRasterPx) const
{
    //std::cout<<"normal for im pt "<<aRasterPx<<"\n";
    //std::cout<<"gnd pt "<<Image2Ground(aRasterPx)<<"\n";
    //std::cout<<"instr vect "<<Image2Camera3D(aRasterPx)<<"\n";
    //std::cout<<"TPD vect "<<Image2ThetaPhiDist(aRasterPx)<<"\n";

    // all computation is done in Raster frame
    cDataGenUnTypedIm<2> & aGenDImDist = getRasterDistance();
    const auto aInterpN = getNormalInterpolator( aRasterPx);
    auto [aDist, aDistGr] = aGenDImDist.GetValueAndGradInterpol(*aInterpN,aRasterPx);
    auto aTPD = Image2ThetaPhiDist(aRasterPx);
    // differencial of cartesian point regarding theta and phi
    // x = d cos(phi) sin(theta)
    // y = d sin(phi)
    // z = d cos(phi) cos(theta)

    //std::cout<<"recalcul instr "<<aTPD.z() * cos(aTPD.y()) * sin(aTPD.x())
    //    <<" "<<aTPD.z() * sin(aTPD.y())<<" "<<aTPD.z() * cos(aTPD.y()) * cos(aTPD.x())<<"\n";


    tREAL8 aDiffDistTheta = aDistGr.x() * InternalCalib()->F();
    tREAL8 aDiffDistPhi   = aDistGr.y() * InternalCalib()->F();
    cPt3dr aDiffPtTheta(
        aTPD.z() * cos(aTPD.y()) * cos(aTPD.x()) + aDiffDistTheta * cos(aTPD.y()) * sin(aTPD.x()),
        aDiffDistTheta * sin(aTPD.y()),
        -aTPD.z() * cos(aTPD.y()) * sin(aTPD.x()) + aDiffDistTheta * cos(aTPD.y()) * cos(aTPD.x())
        );
    cPt3dr aDiffPtPhi(
        -aTPD.z() * sin(aTPD.y()) * sin(aTPD.x()) + aDiffDistPhi * cos(aTPD.y()) * sin(aTPD.x()),
        aTPD.z() * cos(aTPD.y()) + aDiffDistPhi * sin(aTPD.y()),
        -aTPD.z() * sin(aTPD.y()) * cos(aTPD.x()) + aDiffDistPhi * cos(aTPD.y()) * cos(aTPD.x())
        );
    cPt3dr aN = aDiffPtTheta ^ aDiffPtPhi;
    tREAL8 aNNorm = Norm2(aN);
    if (aNNorm>0.)
        aN = aN / aNNorm;
    //std::cout<<"Dist v dx dy "<<aDist<<" "<<aDiffDistTheta<<" "<<aDiffDistPhi<<"\n";
    //std::cout<<"Normal in Raster frame: "<<aN<<"\n";
    //std::cout<<"Normal in Ground frame: "<<Pose().Rot().Value(aN)<<"\n";


    // test with std::pair<cPlane3D,tREAL8> LSQEstimate(co
    return aN;
}



void cStaticLidar::TriangulateRegular(const std::string & aVisuPath, int aFactor)
{
    if (mTriangulation)
        delete mTriangulation;

    //ScopedTimer aTimer("ToTriangulation3DRegular");
    tREAL8 aLimitCosTriangles = 0.999;
    tREAL8 aLimitLenOnMinDist = 50/InternalCalib()->F();
    //TODO remove triangles too much in view direction?

    //TODO make triangulation in instrument frame, add pose later?
    std::vector<cPt3dr> aVPt3D;
    std::vector<cPt2di> aVPt2D;
    std::vector<bool>   aVPtOk;
    std::vector<cPt3di> aVFace;
    MMVII_INTERNAL_ASSERT_tiny(aFactor>0, "Error factor triangulation")

    int aNewW = int((PixelDomain().Sz().x()-1)/aFactor) + 1;
    int aNewH = int((PixelDomain().Sz().y()-1)/aFactor) + 1;
    aVPt3D.reserve(aNewH*aNewW);
    aVPt2D.reserve(aNewH*aNewW);
    aVPtOk.reserve(aNewH*aNewW);
    aVFace.reserve(aVPt3D.capacity()*2);

    for (int l = 0 ; l < PixelDomain().Sz().y(); l+=aFactor)
        for (int c = 0 ; c < PixelDomain().Sz().x(); c+=aFactor)
        {
            aVPt3D.push_back(Image2Ground(cPt2di(c, l))); // or Image2Camera3D
            aVPt2D.push_back(cPt2di(c, l));
            aVPtOk.push_back(IsValidPoint(cPt2dr(c, l)));
        }

    // TODO do not make triangle if any point inside triangle is masked
    for (int l = 0 ; l < aNewH-1; ++l)
        for (int c = 0 ; c < aNewW-1; ++c)
        {
            int aK0 = l*aNewW+c;
            const std::array<std::array<int,3>,2> aVKs = {{{aK0, aK0+aNewW, aK0+1}, {aK0+1, aK0+aNewW, aK0+1+aNewW}}};
            for (auto & [aKa, aKb, aKc]: aVKs)
            {
                if (aVPtOk[aKa] && aVPtOk[aKb] && aVPtOk[aKc])
                {
                    // do not use elongated triangles
                    auto aAB = aVPt3D[aKb]-aVPt3D[aKa];
                    auto aBC = aVPt3D[aKc]-aVPt3D[aKb];
                    auto aCA = aVPt3D[aKa]-aVPt3D[aKc];
                    tREAL4 aLenAB = Norm2(aAB);
                    tREAL4 aLenBC = Norm2(aBC);
                    tREAL4 aLenCA = Norm2(aCA);
                    tREAL4 aABdotBC = Scal(aAB,aBC);
                    tREAL4 aBCdotCA = Scal(aBC,aCA);
                    tREAL4 aCAdotAB = Scal(aCA,aAB);
                    if (fabs(aABdotBC) / (aLenAB*aLenBC) > aLimitCosTriangles)
                        continue;
                    if (fabs(aBCdotCA) / (aLenBC*aLenCA) > aLimitCosTriangles)
                        continue;
                    if (fabs(aCAdotAB) / (aLenCA*aLenAB) > aLimitCosTriangles)
                        continue;
                    // check size vs distance
                    tREAL4 aDistA = getRasterDistance().GetV(cPt2di(aVPt2D[aKa].x(), aVPt2D[aKa].y()));
                    tREAL4 aDistB = getRasterDistance().GetV(cPt2di(aVPt2D[aKb].x(), aVPt2D[aKb].y()));
                    tREAL4 aDistC = getRasterDistance().GetV(cPt2di(aVPt2D[aKc].x(), aVPt2D[aKc].y()));
                    tREAL4 aMaxLen = std::max({aLenAB, aLenBC, aLenCA});
                    tREAL4 aMinDist = std::min({aDistA, aDistB, aDistC});
                    //std::cout<<aMaxLen <<" / "<<aMinDist << " = " << aMaxLen/aMinDist << " > " << aLimitLenOnMinDist << ": " << (aMaxLen/aMinDist > aLimitLenOnMinDist) <<"\n";
                    if (aMaxLen/aMinDist > aLimitLenOnMinDist)
                        continue;
                    aVFace.push_back(cPt3di(aKa, aKb, aKc));
                }
            }
        }

    mTriangulation = new cTriangulation3D<tREAL8>(aVPt3D, aVFace);
    StdOut() <<"Scan triangulation regular "<< mStationName+"_"+mScanName <<": "<<aVPt3D.size()<<" pts, "<<aVFace.size()<<" faces\n";
    mTriangulation->WriteFile(aVisuPath + mStationName+"_"+mScanName+"_regular.ply",true);
}


void cStaticLidar::Triangulate(const std::string & aVisuPath, int aFactor)
{
    if (mTriangulation)
        delete mTriangulation;

    //ScopedTimer aTimer("ToTriangulation3D");
    tREAL8 aLimitCosTriangles = 0.9999;
    tREAL8 aLimitLenOnMinDist = 0.2;
    //TODO make triangulation in instrument frame, add pose later?
    std::vector<cPt2dr> aVPt2D; // get regular points in 2D
    std::vector<cPt3dr> aVPt3D; // swap points 2D for 3D
    std::vector<cPt3di> aVFaceFiltered; // keep good 3D triangles
    MMVII_INTERNAL_ASSERT_tiny(aFactor>0, "Error factor triangulation")

    int aNewW = int((PixelDomain().Sz().x()-1)/aFactor) + 1;
    int aNewH = int((PixelDomain().Sz().y()-1)/aFactor) + 1;
    // maximum sizes
    aVPt3D.reserve(aNewH*aNewW);
    aVPt2D.reserve(aNewH*aNewW);
    aVFaceFiltered.reserve(aVPt3D.capacity()*2);

    // get 2d points
    for (int l = 0 ; l < PixelDomain().Sz().y(); l+=aFactor)
        for (int c = 0 ; c < PixelDomain().Sz().x(); c+=aFactor)
        {
            if (getRasterDistance().GetV(cPt2di(c, l))>1e-4)
                aVPt2D.push_back(cPt2dr(c, l));
        }

    // triangulation
    cTriangulation2D<tREAL8> aTriRaster(aVPt2D);
    aTriRaster.MakeDelaunay();

    // get 3D corresponding points
    for (auto & aPt2d: aVPt2D)
        aVPt3D.push_back(Image2Ground(cPt2di(aPt2d.x(),aPt2d.y())));

    for (const auto & aFace: aTriRaster.VFaces())
    {
        const int & aKa = aFace.x();
        const int & aKb = aFace.y();
        const int & aKc = aFace.z();
        auto aAB = aVPt3D[aKb]-aVPt3D[aKa];
        auto aBC = aVPt3D[aKc]-aVPt3D[aKb];
        auto aCA = aVPt3D[aKa]-aVPt3D[aKc];

        // check angles
        tREAL4 aLenAB = Norm2(aAB);
        tREAL4 aLenBC = Norm2(aBC);
        tREAL4 aLenCA = Norm2(aCA);
        tREAL4 aABdotBC = Scal(aAB,aBC);
        tREAL4 aBCdotCA = Scal(aBC,aCA);
        tREAL4 aCAdotAB = Scal(aCA,aAB);
        if (fabs(aABdotBC) / (aLenAB*aLenBC) > aLimitCosTriangles)
            continue;
        if (fabs(aBCdotCA) / (aLenBC*aLenCA) > aLimitCosTriangles)
            continue;
        if (fabs(aCAdotAB) / (aLenCA*aLenAB) > aLimitCosTriangles)
            continue;

        // check size vs distance
        tREAL4 aDistA = getRasterDistance().GetV(cPt2di(aVPt2D[aKa].x(), aVPt2D[aKa].y()));
        tREAL4 aDistB = getRasterDistance().GetV(cPt2di(aVPt2D[aKb].x(), aVPt2D[aKb].y()));
        tREAL4 aDistC = getRasterDistance().GetV(cPt2di(aVPt2D[aKc].x(), aVPt2D[aKc].y()));
        tREAL4 aMaxLen = std::max({aLenAB, aLenBC, aLenCA});
        tREAL4 aMinDist = std::min({aDistA, aDistB, aDistC});
        if (aMaxLen/aMinDist > aLimitLenOnMinDist)
            continue;
        aVFaceFiltered.push_back(aFace);
    }

    mTriangulation = new cTriangulation3D<tREAL8>(aVPt3D, aVFaceFiltered);
    //StdOut() <<"Scan triangulation "<< mStationName+"_"+mScanName <<": "<<aVPt3D.size()<<" pts, "<<aVFaceFiltered.size()<<" faces\n";
    mTriangulation->WriteFile(aVisuPath + mStationName+"_"+mScanName+".ply",true);

}

cTriangulation3D<tREAL8> * cStaticLidar::getTriangulation() const
{
    return mTriangulation;
}

std::string  cStaticLidar::V_PrefixName() const { return PrefixName() ; }
std::string  cStaticLidar::PrefixName()  { return "Scan";}

std::string  cStaticLidar::Pat2Sup(const std::string & aPatSelect)
{
    return "Ori-" + cStaticLidar::PrefixName() + "-" + aPatSelect + "\\." + GlobTaggedNameDefSerial();
}

cDataIm2D<tREAL4> & cStaticLidar::getRasterDistance() const
{
    MMVII_INTERNAL_ASSERT_tiny(mAreRastersReady, "Error: rasters not ready");
    return mRasterDistance.get()->DIm();
}

bool cStaticLidar::IsValidPoint(const cPt2dr &aRasterPx) const
{
    MMVII_INTERNAL_ASSERT_tiny(mRasterDistance, "Error: mRasterMask must be computed first");
    auto & aDistanceImData = mRasterDistance->DIm();
    return aDistanceImData.InsideBL(aRasterPx)
           && (aDistanceImData.GetV(cPt2di(aRasterPx.x()+0.5,aRasterPx.y()+0.5))>0);
}

bool cStaticLidar::IsValidPoint(const cPt2di &aRasterPx) const
{
    MMVII_INTERNAL_ASSERT_tiny(mRasterDistance, "Error: mRasterMask must be computed first");
    auto & aDistanceImData = mRasterDistance->DIm();
    return aDistanceImData.Inside(aRasterPx)
           && (aDistanceImData.GetV(aRasterPx)>0);
}


bool cStaticLidar::IsMaskedPoint(const cPt2dr &aRasterPx) const
{
    MMVII_INTERNAL_ASSERT_tiny(mRasterMask, "Error: mRasterMask must be computed first");
    auto & aMaskImData = mRasterMask->DIm();
    return aMaskImData.InsideBL(aRasterPx)
           && (aMaskImData.GetV(cPt2di(aRasterPx.x()+0.5,aRasterPx.y()+0.5))==0);
}


tREAL8 cStaticLidar::Sigma() const
{
    return mSigma;
}

const std::vector<cPt2di> & cStaticLidar::PatchCenters() const
{
    return mPatchCenters;
}

void cStaticLidar::FixPtPxLoopAroundPP(cPt2dr &aPtPx) const
{
    tREAL8 aW2piInPixels = InternalCalib()->MapPProj2Im().F()*2*M_PI;
    auto aDX = aPtPx.x() - InternalCalib()->MapPProj2Im().PP().x();
    if (aDX>aW2piInPixels/2)
        aPtPx.x() -=  aW2piInPixels;
    if (aDX<-aW2piInPixels/2)
        aPtPx.x() +=  aW2piInPixels;
}

cPt2dr cStaticLidar::Ground2ImagePrecise(const cPt3dr & aGroundPt) const
{
    MMVII_INTERNAL_ASSERT_tiny(mAreRastersReady, "Error: rasters not ready");
    //std::cout<<"  Ground2ImagePrecise for point "<<aGroundPt<<"\n";
    cPt2dr aDirCam3DTheoretical = InternalCalib()->Dir_Proj()->Value(Pose().Inverse(aGroundPt));
    //std::cout<<"  UV th: "<<aDirCam3DTheoretical<<"\n";
    cPt3dr aPtRasterApprox =cSensorCamPC::Ground2ImageAndDepth(aGroundPt);
    cPt2dr aPtRaster = {aPtRasterApprox.x(), aPtRasterApprox.y()};

    // test if int value
    cPt2di aPtRasterRounded(round(aPtRaster.x()),round(aPtRaster.y()));
    if (Norm2(aPtRaster - cPt2dr(aPtRasterRounded.x(),aPtRasterRounded.y()))< 1e-5)
    {
        // in this case Image2Camera3D will not use GetVBL => works on first and last columns
        cPt2dr aDirTest = InternalCalib()->Dir_Proj()->Value(Image2Camera3D(aPtRasterRounded));
        if (Norm2(aDirTest - aDirCam3DTheoretical)< 1e-5)
        {
            //std::cout<<"  skip iter\n";
            return aPtRaster;
        }
    }

    // if approx is sufficient, no need for iter
    cPt3dr aPtCam3D = Image2Camera3D(aPtRaster);
    if (IsNull(aPtCam3D))
    {
        // we have no data for this point => use default Ground2Image
        return cSensorCamPC::Ground2Image(aGroundPt);
    }

    cPt2dr aDirTest = InternalCalib()->Dir_Proj()->Value(aPtCam3D);
    if (Norm2(aDirTest - aDirCam3DTheoretical)< 1e-5)
    {
        //std::cout<<"  skip iter\n";
        return aPtRaster;
    }

    for (int i = 0; i<3; ++i)
    {
        //std::cout<<"   raster: "<<aPtRaster<<"\n";
        cPt2di aPtRasterUL((int)aPtRaster.x(), (int)aPtRaster.y());
        cPt2di aPtRasterLR((int)aPtRaster.x()+1, (int)aPtRaster.y()+1);
        cPt3dr aPtCam3DUL = Image2Camera3D(aPtRasterUL);
        if (IsNull(aPtCam3DUL))
            return aPtRaster; // impossible to continue
        cPt2dr aDirUL = InternalCalib()->Dir_Proj()->Value(aPtCam3DUL);
        cPt3dr aPtCam3DLR = Image2Camera3D(aPtRasterLR);
        if (IsNull(aPtCam3DLR))
            return aPtRaster; // impossible to continue
        cPt2dr aDirLR = InternalCalib()->Dir_Proj()->Value(aPtCam3DLR);
        //std::cout<<"   Dirs: "<<aDirUL<<" "<<aDirLR<<"\n";
        float aDiffDirX = aDirLR.x()-aDirUL.x();
        float aDiffDirY = aDirLR.y()-aDirUL.y();
        float aBetterX =  (aDiffDirX!=0) ?
                            aPtRasterUL.x() + (aDirCam3DTheoretical.x()-aDirUL.x())/aDiffDirX
                                               *(aPtRasterLR.x()-aPtRasterUL.x())
                                        : aPtRasterUL.x();
        float aBetterY = (aDiffDirY!=0) ?
                            aPtRasterUL.y() + (aDirCam3DTheoretical.y()-aDirUL.y())/aDiffDirY
                                               *(aPtRasterLR.y()-aPtRasterUL.y())
                                        : aPtRasterUL.y();
        aPtRaster = {aBetterX, aBetterY};
    }

    return aPtRaster;
}

cPt2dr cStaticLidar::Ground2Image(const cPt3dr & aGroundPt) const
{
    return Ground2ImagePrecise(aGroundPt);
}

cPt3dr cStaticLidar::Ground2ImageAndDepth(const cPt3dr & aGroundPt) const
{
    cPt2dr aImPt = Ground2Image(aGroundPt);
    tREAL8 aDist = Image2Distance(aImPt);
    return {aImPt.x(), aImPt.y(), aDist};
}


void cStaticLidar::ToFile(const std::string & aNameFile) const
{
    SaveInFile((*this),aNameFile);
    std::string aNameCalib = DirOfPath(aNameFile) + mInternalCalib->Name() + "." + GlobTaggedNameDefSerial();
    mInternalCalib->ToFile(aNameCalib);
}


void cStaticLidar::ToPly(const std::string & aName,bool useMask) const
{
    std::vector<cPt4dr> aSelectionXYZI;
    aSelectionXYZI.reserve(SzPix().x()*SzPix().y());
    MMVII_INTERNAL_ASSERT_tiny(mRasterMask, "Error: mRasterMask must be computed first");
    auto & aMaskImData = mRasterMask->DIm();
    cDataIm2D<tU_INT1> * aRasterIntensityData = mRasterIntensity ? &mRasterIntensity->DIm() : nullptr;
    auto & aRasterXData = mRasterX->DIm();
    auto & aRasterYData = mRasterY->DIm();
    auto & aRasterZData = mRasterZ->DIm();
    for (int l = 0 ; l < SzPix().y(); ++l)
    {
        for (int c = 0 ; c < SzPix().x(); ++c)
        {
            cPt2di aPt(c, l);
            if (aMaskImData.GetV(aPt))
                aSelectionXYZI.push_back(cPt4dr(aRasterXData.GetV(aPt),aRasterYData.GetV(aPt),
                                                aRasterZData.GetV(aPt),
                                                aRasterIntensityData?aRasterIntensityData->GetV(aPt):255));
        }
    }

    cMMVII_Ofs anOfs(aName,eFileModeOut::CreateText);

    bool  aMode8 =  false;

    std::string aSpecCoord = aMode8 ? "float64" : "float32";
    anOfs.Ofs() <<  "ply\n";
    anOfs.Ofs() <<  "format ascii 1.0\n";
    anOfs.Ofs() <<  "comment Generated by MMVVI\n";
    anOfs.Ofs() <<  "element vertex " << aSelectionXYZI.size() << "\n";
    anOfs.Ofs() <<  "property " <<  aSpecCoord  <<" x\n";
    anOfs.Ofs() <<  "property " <<  aSpecCoord  <<" y\n";
    anOfs.Ofs() <<  "property " <<  aSpecCoord  <<" z\n";
    anOfs.Ofs() <<  "property uchar intensity\n";
    anOfs.Ofs() <<  "element face 0\n";
    anOfs.Ofs() <<  "end_header\n";


    for (auto& aPt : aSelectionXYZI)
    {
        anOfs.Ofs() << aPt.x() << " " << aPt.y() << " " << aPt.z();
        anOfs.Ofs() << " " << aPt.t();
        anOfs.Ofs() << "\n";
    }
}

template <typename TYPE> void cStaticLidar::fillRaster(const cStaticLidarImporter & aSL_importer,
                              std::function<TYPE (int)> func, std::unique_ptr<cIm2D<TYPE> > & aIm)
{
    MMVII_INTERNAL_ASSERT_tiny(aSL_importer.mVectPtsCol.size()==aSL_importer.mVectPtsXYZ.size(), "Error: Compute line/col numbers before fill raster");

    aIm.reset(new cIm2D<TYPE>(cPt2di(aSL_importer.NbCol(), aSL_importer.NbLine()), 0, eModeInitImage::eMIA_Null));
    auto & aRasterData = aIm->DIm();
    for (size_t i=0; i<aSL_importer.mVectPtsTPD.size(); ++i)
    {
        cPt2di aPcl = {aSL_importer.mVectPtsCol[i], aSL_importer.mVectPtsLine[i]};
        aRasterData.SetV(aPcl, func(i));
    }
}

std::string cStaticLidar::RasterIntensityPath(const std::string & aImName)
{
    return aImName + "_intensity.tif";
}

std::string cStaticLidar::RasterIntensityPath(const cPhotogrammetricProject & aPhProj, const std::string & aImIDName)
{
    return RasterIntensityPath(aPhProj.DirStaticLidarRasters()+cStaticLidar::NameFromId(aImIDName, false));
}

void cStaticLidar::FillRasters(const cStaticLidarImporter & aSL_importer)
{
    fillRaster<tU_INT1>(aSL_importer,
                        [&aSL_importer](int i)
                        {
                            auto aPtAng = aSL_importer.mVectPtsTPD[i];
                            return (aPtAng.z()<aSL_importer.DistMinToExist())?0:255;
                        }, mRasterMask);

    if (aSL_importer.HasIntensity())
        fillRaster<tU_INT1>(aSL_importer,
                            [&aSL_importer](int i){return aSL_importer.mVectPtsIntens[i]*255;}, mRasterIntensity);

    fillRaster<tREAL4>(aSL_importer,
                      [&aSL_importer](int i){auto aPtAng = aSL_importer.mVectPtsTPD[i];return aPtAng.z();},
                       mRasterDistance);

    fillRaster<tREAL4>(aSL_importer, [&aSL_importer](int i){auto aPtXYZ = aSL_importer.mVectPtsXYZ[i];return aPtXYZ.x();}, mRasterX);
    fillRaster<tREAL4>(aSL_importer, [&aSL_importer](int i){auto aPtXYZ = aSL_importer.mVectPtsXYZ[i];return aPtXYZ.y();}, mRasterY);
    fillRaster<tREAL4>(aSL_importer, [&aSL_importer](int i){auto aPtXYZ = aSL_importer.mVectPtsXYZ[i];return aPtXYZ.z();}, mRasterZ);

    /*fillRaster<tREAL4>(aSL_importer, aPhProjDirOut, mRasterThetaPath, [&aSL_importer](int i){auto aPtAng = aSL_importer.mVectPtsTPD[i];return aPtAng.x();}, saveRasters );
    fillRaster<tREAL4>(aSL_importer, aPhProjDirOut, mRasterPhiPath, [&aSL_importer](int i){auto aPtAng = aSL_importer.mVectPtsTPD[i];return aPtAng.y();}, saveRasters );
    fillRaster<tREAL4>(aSL_importer, aPhProjDirOut, mRasterThetaErrPath, [&aSL_importer](int i)
                      {
                          auto aPtAng = aSL_importer.mVectPtsTPD[i];
                          tREAL8 aThetaCol = aSL_importer.ThetaStart() + aSL_importer.ThetaStep() * aSL_importer.mVectPtsCol[i];
                          aThetaCol = toMinusPiPlusPi(aThetaCol);
                          return aPtAng.x()-aThetaCol;
                      }, saveRasters );
    fillRaster<tREAL4>(aSL_importer, aPhProjDirOut, mRasterPhiErrPath, [&aSL_importer](int i)
                      {
                          auto aPtAng = aSL_importer.mVectPtsTPD[i];
                          tREAL8 aPhiLine = aSL_importer.PhiStart() + aSL_importer.PhiStep() * aSL_importer.mVectPtsLine[i];
                          aPhiLine = toMinusPiPlusPi(aPhiLine);
                          return aPtAng.y()-aPhiLine;
                      }, saveRasters );
    */
    mRasterScore.reset(new cIm2D<tREAL4>(cPt2di(aSL_importer.NbCol()+1, aSL_importer.NbLine()+1), 0, eModeInitImage::eMIA_Null));

    mAreRastersReady = true;
}


void cStaticLidar::SaveRasters(const cStaticLidarImporter & aSL_importer, const std::string &aPhProjDirOut)
{
    mRasterDistancePath = NameImage() + "_distance.tif";
    if (aSL_importer.HasIntensity())
        mRasterIntensityPath = RasterIntensityPath(NameImage());
    else
        mRasterIntensityPath = "";
    mRasterMaskPath = NameImage() + "_mask.tif";
    mRasterXPath = NameImage() + "_X.tif";
    mRasterYPath = NameImage() + "_Y.tif";
    mRasterZPath = NameImage()+ "_Z.tif";

    //mRasterThetaPath = NameImage() + "_Theta.tif";
    //mRasterPhiPath = NameImage() + "_Phi.tif";
    //mRasterThetaErrPath = NameImage() + "_ThetaErr.tif";
    //mRasterPhiErrPath = NameImage() + "_PhiErr.tif";

    mRasterDistance->DIm().ToFile(aPhProjDirOut + mRasterDistancePath);

    // do not save intensity raster, it should have been done before, before decimation
    //mRasterIntensity->DIm().ToFile(aPhProjDirOut + mRasterIntensityPath);

    mRasterMask->DIm().ToFile(aPhProjDirOut + mRasterMaskPath);
    mRasterX->DIm().ToFile(aPhProjDirOut + mRasterXPath);
    mRasterY->DIm().ToFile(aPhProjDirOut + mRasterYPath);
    mRasterZ->DIm().ToFile(aPhProjDirOut + mRasterZPath);
}


std::string cStaticLidar::NameFromId(const std::string &aIdName, bool getOriName)
{
    if (!IsNameTSL(aIdName))
        return MMVII_NONE;
    //MMVII_INTERNAL_ASSERT_User(ends_with(aIdName,GetIdSuffix()),eTyUEr::eBadFileRelName,"Error, Scan Id image must end in "+GetIdSuffix());
    std::string aNameImage = aIdName;
    if (getOriName)
        return NameOri_From_PrefixAndImage(PrefixName(),aNameImage);
    else
        return aNameImage;
}

bool cStaticLidar::IsNameTSL(const std::string &aImageName)
{
    return ends_with(aImageName,GetIdSuffix());
}

cCalculator<double> * cStaticLidar::CreateEqColinearity(bool WithDerives, int aSzBuf, bool ReUse)
{
    return EqTSL_GCP(WithDerives,aSzBuf,ReUse);
}

cCalculator<double> * cStaticLidar::CreateEqColinearityDist(bool WithDerives, int aSzBuf, bool ReUse)
{
    return EqTSL_GCPD(WithDerives,aSzBuf,ReUse);
}

cCalculator<double> * cStaticLidar::GetEqColinearityDist()
{
    if (!mEqDistColinearityDist)
        mEqDistColinearityDist  = CreateEqColinearityDist(true,10,true);
    return mEqDistColinearityDist;
}

void cStaticLidar::cStaticLidar::PushOwnObsColinearity(std::vector<double> & aVObs, const cPt3dr &)
{
    aVObs.push_back(InternalCalib()->F());
    aVObs.push_back(InternalCalib()->PP().x());
    aVObs.push_back(InternalCalib()->PP().y());
    mPose_WU.PushObs(aVObs,true);
}


void cStaticLidar::cStaticLidar::PushOwnObsColinearityDistance(std::vector<double> & aVObs, tREAL4 aMesDistance)
{
    aVObs.push_back(aMesDistance);
    PushOwnObsColinearity(aVObs,{});
}

void cStaticLidar::FilterIntensity(const cStaticLidarImporter &aSL_importer, tREAL8 aLowest, tREAL8 aHighest)
{
    if (!aSL_importer.HasIntensity())
        return;
    MMVII_INTERNAL_ASSERT_tiny(mRasterMask, "Error: mRasterMask must be computed first");
    auto & aMaskImData = mRasterMask->DIm();
    auto & aRasterScoreData = mRasterScore->DIm();
    tREAL8 aMiddle = (aLowest + aHighest) / 2.;
    for (size_t i=0; i<aSL_importer.mVectPtsTPD.size(); ++i)
    {
        cPt2di aPcl = {aSL_importer.mVectPtsCol[i], aSL_importer.mVectPtsLine[i]};
        if ((aSL_importer.mVectPtsIntens[i]<aLowest) || (aSL_importer.mVectPtsIntens[i]>aHighest))
            aMaskImData.SetV(aPcl, 0);
        aRasterScoreData.SetV(aPcl, aRasterScoreData.GetV(aPcl) + fabs(aSL_importer.mVectPtsIntens[i]-aMiddle));
    }
    //aMaskImData.ToFile("MaskIntens.png");
}

void cStaticLidar::FilterIncidence(const cStaticLidarImporter &aSL_importer, tREAL8 aAngMax)
{
    MMVII_INTERNAL_ASSERT_tiny(mRasterMask, "Error: mRasterMask must be computed first");
    auto & aMaskImData = mRasterMask->DIm();
    auto & aRasterScoreData = mRasterScore->DIm();

    // TODO: use im.InitCste()
    cIm2D<tREAL4> aImDistGrX(cPt2di(aSL_importer.NbCol()+1, aSL_importer.NbLine()+1), 0, eModeInitImage::eMIA_Null);
    auto & aImDistGrXData = aImDistGrX.DIm();
    cIm2D<tREAL4> aImDistGrY(cPt2di(aSL_importer.NbCol()+1, aSL_importer.NbLine()+1), 0, eModeInitImage::eMIA_Null);
    auto & aImDistGrYData = aImDistGrY.DIm();

    tREAL4 aTanAngMax = tan(aAngMax);

    // gaussian blur of masked image: blur image and mask, for valid pixels, result = blured_im/blured_mask
    auto aRasterDistGauss = mRasterDistance->Dup();
    auto & aRasterDistGaussData = aRasterDistGauss.DIm();
    ExpFilterOfStdDev(aRasterDistGaussData, 2, 2.);

    //mRasterMask->DIm().ToFile("Mask.tif");
    auto aRasterMaskGauss = Convert((float*)nullptr, mRasterMask->DIm()) * (1./255.);
    auto & aRasterMaskGaussData = aRasterMaskGauss.DIm();
    ExpFilterOfStdDev(aRasterMaskGaussData, 2, 2.);

    //aRasterDistGaussData.ToFile("DistGaussData.tif");
    //aRasterMaskGaussData.ToFile("MaskGaussData.tif");

    cImGrad<tREAL4> aDistGradIm(aRasterDistGauss);
    ComputeSobel<tREAL4,tREAL4>(*aDistGradIm.mDGx, *aDistGradIm.mDGy, aRasterDistGaussData);
    for (int l = 0 ; l < aSL_importer.NbLine(); ++l)
    {
        //tREAL4 phi = lToPhiApprox(l, aSL_importer.PhiStart(), aSL_importer.PhiStep());
        tREAL4 phi = InternalCalib()->DirBundle({0.,(double)l}).y();
        tREAL4 aStepThetaFix = aSL_importer.ThetaStep()*cos(phi);
        for (int c = 0 ; c < aSL_importer.NbCol(); ++c)
        {
            cPt2di aPt(c, l);
            tREAL4 aDist = mRasterDistance->DIm().GetV(aPt);
            tREAL4 aValDistGradX = aDistGradIm.Gx(aPt);
            tREAL4 aValDistGradY = aDistGradIm.Gy(aPt);
            tREAL4 aValGaussMask = aRasterMaskGaussData.GetV(aPt);
            tREAL4 aValMask = mRasterMask->DIm().GetV(aPt);
            if (! aValMask)
                continue;
            tREAL4 aTanIncidX = aValDistGradX / (aStepThetaFix * aDist) / aValGaussMask;
            aImDistGrXData.SetV(aPt, aTanIncidX);
            tREAL4 aTanIncidY = aValDistGradY / (aSL_importer.PhiStep() * aDist) / aValGaussMask;
            aImDistGrYData.SetV(aPt, aTanIncidY);
            if (fabs(aTanIncidX*aTanIncidX+aTanIncidY*aTanIncidY)>aTanAngMax*aTanAngMax)
                aMaskImData.SetV(aPt, 0);
            aRasterScoreData.SetV(aPt, aRasterScoreData.GetV(aPt) + 10.*fabs(aTanIncidX*aTanIncidX+aTanIncidY*aTanIncidY));
        }
    }
    //aImDistGrXData.ToFile("DistGrXData.tif");
    //aImDistGrYData.ToFile("DistGrYData.tif");
    //aMaskImData.ToFile("MaskIncidence.png");
}

void cStaticLidar::FilterDistance(tREAL8 aDistMin, tREAL8 aDistMax)
{
    MMVII_INTERNAL_ASSERT_tiny(mRasterMask, "Error: mRasterMask must be computed first");
    auto & aMaskImData = mRasterMask->DIm();
    //auto & aRasterScoreData = mRasterScore->DIm(); // add something to mRasterScore?
    auto & aRasterDistData = mRasterDistance->DIm();
    for (int l = 0 ; l < SzPix().y(); ++l)
    {
        for (int c = 0 ; c < SzPix().x(); ++c)
        {
            cPt2di aPt(c, l);
            tREAL4 aDist = aRasterDistData.GetV(aPt);
            if ((aDist<aDistMin)||(aDist>aDistMax))
            {
                aMaskImData.SetV(aPt, 0);
            }
        }
    }
}

void cStaticLidar::MaskBuffer(const cStaticLidarImporter &aSL_importer, tREAL8 aAngBuffer, const std::string &aPhProjDirOut)
{
    if (aAngBuffer==0)
        return;
    StdOut() << "Computing Mask buffer..."<<std::endl;
    MMVII_INTERNAL_ASSERT_tiny(mRasterMask, "Error: mRasterMask must be computed first");
    auto & aMaskImData = mRasterMask->DIm();
    mRasterMaskBuffer.reset( new cIm2D<tU_INT1>(cPt2di(aSL_importer.NbCol(), aSL_importer.NbLine()), 0, eModeInitImage::eMIA_NoInit));
    auto & aMaskBufImData = mRasterMaskBuffer->DIm();

    auto & aRasterScoreData = mRasterScore->DIm();

    bool aHzLoop = false;
    if (fabs(fabs(aSL_importer.ThetaStep()) * (aSL_importer.NbCol()+1) - 2 * M_PI) < 2 * fabs(aSL_importer.ThetaStep()))
        aHzLoop = true;

    tREAL8 aRadPx = aAngBuffer/aSL_importer.PhiStep();
    aMaskBufImData.InitCste(255);

    std::vector<bool> aLinesFull(aSL_importer.NbLine()+1, false); // record lignes completely masked to pass them next time
    for (int l = 0 ; l < aSL_importer.NbLine(); ++l)
    {
        for (int c = 0 ; c < aSL_importer.NbCol(); ++c)
        {
            auto aMaskVal = aMaskImData.GetV(cPt2di(c, l));
            if (aMaskVal==0)
            {
                for (int il = l - aRadPx; il <= l + aRadPx; ++il)
                {
                    if ((il<0) || (il>=aSL_importer.NbLine())) continue;
                    if (aLinesFull[il]) continue;
                    //tREAL8 phi = InternalCalib()->DirBundle({0.,(double)l}).y();
                    tREAL8 w = fabs(sqrt(aRadPx*aRadPx - (il-l)*(il-l)) ); // is  /cos(phi)); useful?
                    if (w>aSL_importer.NbCol())
                    {
                        w=aSL_importer.NbCol();
                        aLinesFull[il] = true;
                        // TODO: fill line and continue
                    }
                    for (int ic = c - w; ic <= c + w; ++ic)
                    {
                        int icc = ic; // working copy
                        if (aHzLoop)
                        {
                            if (icc<0)
                                icc += (aSL_importer.NbCol()+1);
                            if (icc>aSL_importer.NbCol())
                                icc -= (aSL_importer.NbCol()+1);
                        }
                        if ((icc<0)||(icc>=aSL_importer.NbCol()))
                            continue;
                        aMaskBufImData.SetV(cPt2di(icc, il), 0);
                    }
                }
            }
        }
    }
    for (int l = 0 ; l < aSL_importer.NbLine(); ++l)
        for (int c = 0 ; c < aSL_importer.NbCol(); ++c)
        {
            if (aMaskBufImData.GetV(cPt2di(c, l))==0)
            {
                aRasterScoreData.SetV(cPt2di(c, l), 1000.);
                aMaskImData.SetV(cPt2di(c, l), 0);
            }
        }
}

void cStaticLidar::SelectPatchCenters1(int aNbPatches)
{
    MMVII_INTERNAL_ASSERT_tiny(mAreRastersReady, "Error: rasters not ready");
    mPatchCenters.clear();
    float aNbPatchesFactor = 2.; // a priori search for aNbPatches * aNbPatchesFactor
    auto & aRasterScoreData = mRasterScore->DIm();
    cResultExtremum aRes;
    double aRadius = sqrt(aRasterScoreData.SzX()*aRasterScoreData.SzY()/(M_PI*aNbPatches))/aNbPatchesFactor;
    ExtractExtremum1(aRasterScoreData, aRes, aRadius);
    mPatchCenters = aRes.mPtsMin;
    StdOut() << "Nb patches: " << mPatchCenters.size() <<"\n";
    //std::fstream file1;
    //file1.open("centers.txt", std::ios_base::out);
    //for (auto & aCenter : mPatchCenters)
    //    file1 << aCenter.x() << " " << -aCenter.y() <<"\n";
    //aRasterScoreData.ToFile("Score.tif");
}

void cStaticLidar::SelectPatchCenters2(int aNbPatches, cDataIm2D<tU_INT1> * aSupMaskDIm)
{
    MMVII_INTERNAL_ASSERT_tiny(mAreRastersReady, "Error: rasters not ready");
    if (aSupMaskDIm)
        MMVII_INTERNAL_ASSERT_tiny(
            aSupMaskDIm->Sz() == InternalCalib()->SzPix(),
            "Error: Sup mask must have the same size as TSL mask");
    mPatchCenters.clear();
    auto & aRasterMaskData = mRasterMask->DIm();
    /*cResultExtremum aRes;
    double aRadius = sqrt(aRasterMaskBufferData.SzX()*aRasterMaskBufferData.SzY()/(M_PI*nbPatches));
    ExtractExtremum1(aRasterMaskBufferData, aRes, aRadius);
    mPatchCenters = aRes.mPtsMax;*/

    // regular grid
    auto [aAvgDist, aNbValid, aNbNotMasked] = AvgDistNbValidAndNbNotMasked();
    MMVII_INTERNAL_ASSERT_tiny(
        aNbNotMasked>0,
        "Error: all the scan is masked!");

    auto & aRasterDistData = mRasterDistance->DIm();
    aNbPatches = aNbPatches * sqrt(2); // ??
    float aXYratio=((float)aRasterMaskData.SzX())/aRasterMaskData.SzY();
    int aNbPatchesX = sqrt((double)aNbPatches)*sqrt(aXYratio)+1;
    int aNbPatchesY = sqrt((double)aNbPatches)/sqrt(aXYratio)+1;
    float aNbPatchesFactor = // a priori search for aNbPatches * aNbPatchesFactor, to adjust for no return
        sqrt(PixelDomain().Sz().x()*PixelDomain().Sz().y()/ ((aNbValid+aNbNotMasked)/2.));
    float aX;
    float aY = float(aRasterMaskData.SzY()) / aNbPatchesY / 3.;
    float aYmax = PixelDomain().Sz().y() - aY - 1;
    float aXStep;
    float aYStep = float(aRasterMaskData.SzY()) / aNbPatchesY / aNbPatchesFactor;
    if (aYStep<1.)
        aYStep = 1.;
    int aLineCounter = 0;
    float aXdecal = float(aRasterMaskData.SzX()) / aNbPatchesX;

    std::cout<<"aAvgDist="<<aAvgDist<<" aNbValid="<<aNbValid<<
               " aNbPatchesFactor="<<aNbPatchesFactor<<" aYStep="<<aYStep<<"\n";
    while (aY<aYmax)
    {
        aX = aXdecal * ((aLineCounter%2)?1./3.:2./3.);
        auto aPhi = (aY - InternalCalib()->PP().y()) / InternalCalib()->F();
        float aYStepCurr = aYStep;
        double aLineAvgDist = 0;
        int aLineNbPts = 0;
        while (aX<aRasterMaskData.SzX()-aXdecal*1./3.)
        {
            // take lat/long proj into account
            aXStep = fabs(((float)aRasterMaskData.SzX()) / aNbPatchesX / aNbPatchesFactor / cos(aPhi))+1;
            auto aPt = cPt2di(aX, aY);
            if (aRasterMaskData.GetV(aPt) && ( (!aSupMaskDIm) || aSupMaskDIm->GetV(aPt)))
            {
                mPatchCenters.push_back(aPt);
                auto aDist = aRasterDistData.GetV(aPt);
                aXStep *= aAvgDist/aDist; // take depth into account
                aLineAvgDist += aDist;//*aDist;
                aLineNbPts++;
            } else
                aXStep /= 9.; // if this pixel has no response, search next closer than normal step
            if (aXStep<2.)
                aXStep = 2.;
            aX += aXStep;
        }
        if (aLineNbPts>0)
            aYStepCurr *= aAvgDist/(aLineAvgDist/aLineNbPts);

        aY += aYStepCurr;
        aLineCounter++;
    }

    StdOut() << "Nb patches: " << mPatchCenters.size() <<"\n";
    //std::fstream file1;
    //file1.open("centers.txt", std::ios_base::out);
    //for (auto & aCenter : mPatchCenters)
    //    file1 << aCenter.x() << " " << -aCenter.y() <<"\n";
}

// set a regular grid, subdivide parts if deep, keep best score in cell
void cStaticLidar::SelectPatchCenters3(int aNbPatches, cDataIm2D<tU_INT1> * aSupMaskDIm)
{
    MMVII_INTERNAL_ASSERT_tiny(mAreRastersReady, "Error: rasters not ready");

    if (aSupMaskDIm)
        MMVII_INTERNAL_ASSERT_tiny(
            aSupMaskDIm->Sz() == InternalCalib()->SzPix(),
            "Error: Sup mask must have the same size as TSL mask");

    //auto [aAvgDist, aNbValid, aNbNotMasked] = AvgDistNbValidAndNbNotMasked();
    //std::cout<<"aAvgDist="<<aAvgDist<<" aNbValid="<<aNbValid<<
    //    " aNbNotMasked="<<aNbNotMasked<<"\n";

    auto & aRasterMaskData = mRasterMask->DIm();
    auto & aRasterDistData = mRasterDistance->DIm();

    // increase the number of cells to correspond to empty areas
    float aNbPatchesFactor = 1.1;//sqrt((aRasterDistData.SzX()*aRasterDistData.SzY())/(float(aNbNotMasked+1)));

    //fix dist image : d <= d*cos(aPhi)
    cIm2D<tREAL4> aImDistFix(aRasterDistData.Sz(),0,eModeInitImage::eMIA_Null);
    auto & aImDistFixData = aImDistFix.DIm();
    for (int l = 0 ; l < aRasterDistData.Sz().y(); ++l)
    {
        tREAL4 aPhi = InternalCalib()->DirBundle({0.,(double)l}).y();
        for (int c = 0 ; c <aRasterDistData.Sz().x(); ++c)
        {
            aImDistFixData.GetLine(l)[c] = aRasterDistData.GetLine(l)[c]*aRasterMaskData.GetLine(l)[c]*cos(aPhi);
        }
    }
    mPatchCenters.clear();

    cQuadTree aQuadTree(&aImDistFixData);
    aQuadTree.Split(aNbPatches*aNbPatchesFactor);

    std::cout<<"QuadTree: "<<aQuadTree.GetCurNbCell()<<"\n";

    // make patches, good points on each leafs
    for (auto& aLeaf: aQuadTree.GetVLeafs())
    {
        auto aLeafCenter = (aLeaf->GetArea().P0() + aLeaf->GetArea().P1())/2;
        cPt2di aBestPt = aLeafCenter;
        float aBestScore = 0;
        for (auto &aPx: aLeaf->GetArea())
        {
            float aCurScore = aRasterMaskData.GetV(aPx) * aRasterDistData.GetV(aPx) / (1. + Norm2(aPx-aLeafCenter));
            if (aSupMaskDIm)
                aCurScore *= aSupMaskDIm->GetV(aPx);
            if (aCurScore > aBestScore)
            {
                aBestScore = aCurScore;
                aBestPt = aPx;
            }
        }
        if (aBestScore>0)
            mPatchCenters.push_back(aBestPt);
    }
    StdOut() << "Nb patches: " << mPatchCenters.size() <<"\n";
}

void cStaticLidar::MakeVisu(const cPhotogrammetricProject & aPhProj) const
{
    MMVII_INTERNAL_ASSERT_tiny(mAreRastersReady, "Error: rasters not ready");
    auto & aRasterDistData = mRasterDistance->DIm();
    double aDistMax = 0.;
    int aPtSize = 1 + mRasterDistance->DIm().SzX()/1000;
    for (auto & aPt :  aRasterDistData)
    {
        if (aRasterDistData.GetV(aPt)>aDistMax)
            aDistMax=aRasterDistData.GetV(aPt);
    }
    cRGBImage  aImDist8b = mRasterIntensity?
        RGBImFromGray(mRasterIntensity->DIm(), 1., 1) :
        RGBImFromGray(aRasterDistData, 255./aDistMax,1);
    for (auto & aCenter : mPatchCenters)
    {
        //aImDist8b.SetRGBPoint(ToR(aCenter)-cPt2dr(0.,0.), cRGBImage::Red);
        //aImDist8b.RawSetPoint(aCenter, cRGBImage::Red);
        aImDist8b.FillRectangle(cRGBImage::Red, aCenter-cPt2di(aPtSize,aPtSize), aCenter+cPt2di(aPtSize,aPtSize), {0.,0.,0.});
    }
    std::string aPath = aPhProj.DirVisuAppli() + mStationName + "_" + mScanName + "_patches.png";
    aImDist8b.ToFile(aPath);

}

bool cStaticLidar::isInsideNormalInterpolator(cPt2dr & aPt) const
{
    const auto aInterp = getNormalInterpolator( aPt );
    return getRasterDistance().InsideInterpolator(*aInterp,aPt,1.0);
}


void cStaticLidar::MakePatches
    (std::list<cLidarRasterPatch> &aLPatches,
     const std::vector<cSensorCamPC *> & aVCam,
     int    aNbPointByPatch,
     int    aSzMin
     ) const
{
    //StdOut() << "MakePatches\n";
    MMVII_INTERNAL_ASSERT_tiny(mAreRastersReady, "Error: rasters not ready");
    auto & aRasterDistData = mRasterDistance->DIm();
    auto & aRasterMaskData = mRasterMask->DIm();

    // shortcut if only 1 point needed: just get the centers
    if (aNbPointByPatch==1)
    {
        for (size_t i=0; i<mPatchCenters.size(); ++i)
        {
            auto & aCenter = mPatchCenters[i];
            MMVII_INTERNAL_ASSERT_tiny(IsMaskedPoint(ToR(aCenter))==false, "Error: patch " + ToStr(aCenter) + " is on a masked area");
            auto aCenterR = cPt2dr(aCenter.x(),aCenter.y());
            if (isInsideNormalInterpolator(aCenterR))  // is it sufficiently inside
            {
                auto aN = Image2NormalInstr(aCenterR);
                aLPatches.push_back({i, {aCenter}, {}, aN});
                //std::cout<<"make patch "<<aCenter<<" N="<<aN<<"\n";
            }
        }
        return;
    }
//#define NUMMAKEPATCHDEBUG 0
    std::vector<tREAL8> aVectGndPixelSize;
    aVectGndPixelSize.resize(aVCam.size());
    // parse center points
    for (size_t i=0; i<mPatchCenters.size(); ++i)
    {
        auto & aCenter = mPatchCenters[i];
        auto aCenterR = cPt2dr(aCenter.x(),aCenter.y());
        if (!isInsideNormalInterpolator(aCenterR))  // is it sufficiently inside
            continue;
        //search for average GndPixelSize
        aVectGndPixelSize.clear();
        cPt3dr aGndCenter = Image2Ground(aCenter);
        int aNumCamVisib = 0;
        for (const auto & aCam: aVCam)
        {
            if (aCam->IsVisible(aGndCenter))
            {
                ++aNumCamVisib;
                tREAL8 aDist = Norm2(aCam->Center()-aGndCenter);
                cPt2dr aImCenter = aCam->Ground2Image(aGndCenter);

                aVectGndPixelSize.push_back(Norm2(aCam->ImageAndDepth2Ground(cPt3dr{aImCenter.x(), aImCenter.y(), aDist})
                                                  - aCam->ImageAndDepth2Ground(cPt3dr{aImCenter.x()+1., aImCenter.y(), aDist})));

            }
        }
        if (aNumCamVisib<2) continue;
        tREAL8 aGndPixelSize = NonConstMediane(aVectGndPixelSize);
    #ifdef NUMMAKEPATCHDEBUG
        if (i==NUMMAKEPATCHDEBUG)
            StdOut() << "GndPixelSize: " << aGndPixelSize << "\n";
    #endif
        // compute raster step to get aNbPointByPatch separated by aGndPixelSize
        tREAL4 aMeanDepth = aRasterDistData.GetV(aCenter);

        cPt3dr aCenterThetaPhiDist = Image2ThetaPhiDist(aCenter);
        //StdOut() << "CenterThetaPhiDist: " << aCenterThetaPhiDist << "\n";
        tREAL4 aThetaStep = NAN;
        if (aCenter.x() < mRasterX->DIm().SzX()/2) // TODO: strange way to compute steps?
            aThetaStep = aCenterThetaPhiDist.x() - Image2ThetaPhiDist(aCenter+cPt2di(1, 0)).x();
        else
            aThetaStep = Image2ThetaPhiDist(aCenter-cPt2di(1, 0)).x() - aCenterThetaPhiDist.x();

        tREAL4 aProjColFactor = 1/cos(aCenterThetaPhiDist.y());
        tREAL4 aNbStepRadius = sqrt(aNbPointByPatch+2) / 2.;
        tREAL4 aRasterPxGndW = fabs(aThetaStep) * aMeanDepth * aProjColFactor;
        tREAL4 aRasterPxGndH = fabs(aThetaStep) * aMeanDepth;
        tREAL4 aRasterStepPixelsY = aGndPixelSize / aRasterPxGndH;
        tREAL4 aRasterStepPixelsX = aGndPixelSize / aRasterPxGndW;
    #ifdef NUMMAKEPATCHDEBUG
        if (i==NUMMAKEPATCHDEBUG)
            StdOut() << "RasterStepPixels: " << aRasterStepPixelsX << " " << aRasterStepPixelsY << "\n";
    #endif
        // have a least one scan step of difference between patch points
        if (aRasterStepPixelsX < 1.)
            aRasterStepPixelsX = 1.;
        if (aRasterStepPixelsY < 1.)
            aRasterStepPixelsY = 1.;

        std::vector<cPt2di> aPatchPts = {aCenter}; // convention: center is at first
#ifdef NUMMAKEPATCHDEBUG
        if (i==NUMMAKEPATCHDEBUG)
        {
            for (auto const & aPt:aPatchPts)
                std::cout<<aPt<<" ";
            std::cout<<"\n";
        }
#endif
        for (int aJ = -aNbStepRadius; aJ<=aNbStepRadius; ++aJ)
            for (int aI = -aNbStepRadius; aI<=aNbStepRadius; ++aI)
            {
                if ((aI==0)&&(aJ==0))
                    continue; // do not add center twice
                cPt2di aPt = aCenter + cPt2di(aI*aRasterStepPixelsX,aJ*aRasterStepPixelsY);
                if (aRasterMaskData.Inside(aPt) && aRasterMaskData.GetV(aPt))
                    aPatchPts.push_back(aPt);
            }

        // some requirement on minimal size
        if ((int)aPatchPts.size() > aSzMin)
        {
            auto aN = Image2NormalInstr(aCenterR);
            aLPatches.push_back({i, aPatchPts, {}, aN});
            //std::cout<<"make patch "<<aCenter<<" N="<<aN<<"\n";
        #ifdef NUMMAKEPATCHDEBUG
            if (i==NUMMAKEPATCHDEBUG)
            {
                for (auto const & aPt:aPatchPts)
                    std::cout<<aPt<<" ";
                std::cout<<"\n";
            }
        #endif
        }
    }
}

cIm2D<tU_INT1> cStaticLidar::projectIntensityFrom(const cStaticLidar& aFrom) const
{
    StdOut() << "Reproject " << aFrom.NameImage() << " on " << NameImage() << "\n";
    cIm2D<tU_INT1> aProj(Sz(),nullptr,eModeInitImage::eMIA_Null);
    if (!mRasterIntensity)
    {
        MMVII_UserError(eTyUEr::eUnClassedError, std::string("Error, no intensity for scan \"") + NameImage() + "\": ");
        return aProj;
    }

    auto & aProjDIm = aProj.DIm();
    auto & aFromDIm = aFrom.mRasterIntensity->DIm();
    for (const auto & aP : aProjDIm)
    {
        auto aPgnd = Image2Ground(aP);
        auto aPfrom = ToI(aFrom.Ground2Image(aPgnd));
        if (aFromDIm.Inside(aPfrom))
            aProjDIm.SetV(aP, aFromDIm.GetV(aPfrom));
    }
    return aProj;
}

void cStaticLidar::AddData(const  cAuxAr2007 & anAux)
{
    cSensorCamPC::AddData(anAux);
    MMVII::AddData(cAuxAr2007("StationName",anAux),mStationName);
    MMVII::AddData(cAuxAr2007("ScanName",anAux),mScanName);

    MMVII::AddData(cAuxAr2007("RasterDistance",anAux),mRasterDistancePath);
    MMVII::AddData(cAuxAr2007("RasterIntensity",anAux),mRasterIntensityPath);
    MMVII::AddData(cAuxAr2007("RasterMask",anAux),mRasterMaskPath);
    MMVII::AddData(cAuxAr2007("RasterX",anAux),mRasterXPath);
    MMVII::AddData(cAuxAr2007("RasterY",anAux),mRasterYPath);
    MMVII::AddData(cAuxAr2007("RasterZ",anAux),mRasterZPath);

    //MMVII::AddData(cAuxAr2007("RasterTheta",anAux),mRasterThetaPath);
    //MMVII::AddData(cAuxAr2007("RasterPhi",anAux),mRasterPhiPath);
    //MMVII::AddData(cAuxAr2007("RasterThetaErr",anAux),mRasterThetaErrPath);
    //MMVII::AddData(cAuxAr2007("RasterPhiErr",anAux),mRasterPhiErrPath);
    MMVII::AddData(cAuxAr2007("Sigma",anAux),mSigma);

    MMVII::AddData(cAuxAr2007("RotInput2Raster",anAux),mRotInput2Raster);

    MMVII::AddData(cAuxAr2007("PatchCenters",anAux),mPatchCenters);
}

void AddData(const  cAuxAr2007 & anAux,cStaticLidar & aSL)
{
   aSL.AddData(anAux);
}


template<typename TYPE>
void TestRaster2Gnd2Raster(const std::vector<TYPE> &aVectPtsTest, cStaticLidar * aScan)
{
    float aPrecision = 1e-2;
    if constexpr (std::is_same_v<TYPE, cPt2di>)
        aPrecision = 1e-3;
    long i=0;
    for (auto & aPIm: aVectPtsTest)
    {
        //std::cout<<"Test " << i << ": "<<aPIm<<"\n";
        auto aPgnd = aScan->Image2Ground(aPIm);
        auto aPImtest = aScan->Ground2Image(aPgnd);
        //std::cout<<"Result: "<<aPIm<<" -> "<<aPgnd<<" -> "<<aPImtest<<"\n";
        ++i;
        MMVII_INTERNAL_ASSERT_bench(Norm2(cPt2dr(aPIm.x(), aPIm.y())-aPImtest)<aPrecision ,"TestRaster2Gnd2Raster: " + std::to_string(i));
    }
}

/// tests the scans of a cube, where summit is {0,0,-8.66} in ground coords
void TestPose(const std::string & aInPath, const std::string & aScanName, const cPt2dr& aSummitPx)
{
    cStaticLidar * aScan =  cStaticLidar::FromFile(aInPath + aScanName, false);
    aScan->ReadRasters(aInPath);
    auto aRasterPx = aScan->Ground2Image({0,0,-8.66});
    //std::cout<<"Result: "<<aRasterPx<<" - theoritical "<<aSummitPx<<" -> error "<<Norm2(aRasterPx-aSummitPx)<<"\n";
    MMVII_INTERNAL_ASSERT_bench(Norm2(aRasterPx-aSummitPx)<1e-3 ,"TestPose " + aScanName);
    delete aScan;
}

void BenchTSL(cParamExeBench & aParam)
{
    if (! aParam.NewBench("TSL")) return;

    const std::string & aInPath = cMMVII_Appli::CurrentAppli().InputDirTestMMVII() + "/TSL/Scan1/";

    // test with scan pose = Id
    cStaticLidar * aScan =  cStaticLidar::FromFile(aInPath + "Scan-St1-Sc1.xml", false);
    aScan->ReadRasters(aInPath);

    //aScan->ToPly(cMMVII_Appli::CurrentAppli().TmpDirTestMMVII() + "/TSL.ply");

    auto & pp = aScan->InternalCalib()->PP();
    cPt2di ppInt = cPt2di(round(pp.x()), round(pp.y()));
    auto & sz = aScan->InternalCalib()->SzPix();
    std::vector<cPt2di> aVectPtsTest1 = {ppInt, {0, ppInt.y()}, {sz.x()-2, ppInt.y()}, {0, 0}, {sz.x()-2, sz.y()-2}};
    TestRaster2Gnd2Raster(aVectPtsTest1, aScan);

    std::vector<cPt2dr> aVectPtsTest2;
    for (int i = 0; i<10; ++i)
        aVectPtsTest2.push_back( pp + cPt2dr(pp.x()/10. * i * cos(2*M_PI*i/10),
                                            pp.y()/10. * i * sin(2*M_PI*i/10)));
    TestRaster2Gnd2Raster(aVectPtsTest2, aScan);

    // set a bad F, to make Ground2ImagePrecise different from simple projection
    tREAL8 aOriginalF = aScan->InternalCalib()->MapPProj2Im().F();
    aScan->InternalCalib()->MapPProj2Im().F()*=0.9;
    // TestRaster2Gnd2Raster(aVectPtsTest1, aScan); // TODO: what to do with points out of raster?
    TestRaster2Gnd2Raster(aVectPtsTest2, aScan);
    // restore original For next iteration
    aScan->InternalCalib()->MapPProj2Im().F() = aOriginalF;

    delete aScan;

    // tests with scan translation
    TestPose(aInPath, "Scan-St2-Sc1.xml", {35.5508,50});

    // tests with scan translation + rotation
    TestPose(aInPath, "Scan-St3-Sc1.xml", {40.7816,64.542});

    // just rot x
    TestPose(aInPath, "Scan-St4-Sc1.xml", {66.9119,60.7288});

    // just rot xyz
    TestPose(aInPath, "Scan-St5-Sc1.xml", {67.6836,71.4344});

    //std::cout<<"Bench TSL finished."<<std::endl;
    aParam.EndBench();
    return;
}


};

