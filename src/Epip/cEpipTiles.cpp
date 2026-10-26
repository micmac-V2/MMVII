#include "cEpipolarRectification.h"
#include "MMVII_2Include_Serial_Tpl.h"
#include <algorithm>
#include <limits>

namespace MMVII {

std::string EpipNameWithExtension(const std::string & aPattern)
{
    const size_t aPos = aPattern.rfind('.');
    if ((aPos!=std::string::npos) && (aPos+1<aPattern.size()) && (aPattern.find_first_of("%$/",aPos)==std::string::npos))
        return aPattern;
    return aPattern + ((aPos!=std::string::npos) && (aPos+1==aPattern.size()) ? "tif" : ".tif");
}

std::string EpipTileName(const std::string & aPattern,int aRow,int aCol,int aNbRow,int aNbCol)
{
    const size_t aWidth = std::max<size_t>(2,ToStr(std::max(aNbRow,aNbCol)-1).size());
    auto Pad = [&](int aK)
    {
        const std::string aStr = ToStr(aK);
        return std::string(aWidth>aStr.size() ? aWidth-aStr.size() : 0,'0') + aStr;
    };
    const std::string aTag = "_t" + Pad(aRow) + "_" + Pad(aCol);
    const size_t aPos = aPattern.rfind('.');
    return (aPos==std::string::npos) ? aPattern + aTag : aPattern.substr(0,aPos) + aTag + aPattern.substr(aPos);
}

void cEpipTileEntry::AddData(const cAuxAr2007 &anAux)
{
    MMVII::AddData(cAuxAr2007("Row", anAux), mRow);
    MMVII::AddData(cAuxAr2007("Col", anAux), mCol);
    MMVII::AddData(cAuxAr2007("InfoFile", anAux), mInfoFile);
}

void AddData(const cAuxAr2007 &anAux, cEpipTileEntry &anEntry)
{
    anEntry.AddData(anAux);
}

void cEpipTilesInfo::AddData(const cAuxAr2007 &anAux)
{
    MMVII::AddData(cAuxAr2007("ModelName", anAux), mNameModel);
    MMVII::AddData(cAuxAr2007("Master", anAux), mMaster);
    MMVII::AddData(cAuxAr2007("SzTiles", anAux), mSzTiles);
    MMVII::AddData(cAuxAr2007("SzOverL", anAux), mSzOverL);
    MMVII::AddData(cAuxAr2007("RegionP0", anAux), mRegion0);
    MMVII::AddData(cAuxAr2007("RegionP1", anAux), mRegion1);
    MMVII::AddData(cAuxAr2007("DispRangeSlaveMinusMaster", anAux), mDispRange);
    MMVII::AddData(cAuxAr2007("Tiles", anAux), mTiles);
}

void AddData(const cAuxAr2007 &anAux, cEpipTilesInfo &anInfo)
{
    anInfo.AddData(anAux);
}

void cEpipTilesInfo::ToFile(const std::string &aNameFile) const
{
    // Same formatting precautions as the model : no truncation of small values
    PushPrecTxtSerial(std::numeric_limits<tREAL8>::max_digits10);
    SetFixedFloatTxtSerial(false);
    SaveInFile(*this, aNameFile);
    SetFixedFloatTxtSerial(true);
    PopPrecTxtSerial();
}

cEpipTilesInfo cEpipTilesInfo::FromFile(const std::string &aNameFile)
{
    cEpipTilesInfo anInfo;
    ReadFromFile(anInfo, aNameFile);
    return anInfo;
}

} // namespace MMVII
