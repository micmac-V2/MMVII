#include "cMMVII_Appli.h"
#include "MMVII_Geom3D.h"
#include "MMVII_Mappings.h"

#include "lasreader.hpp"


#include <memory>

namespace MMVII {

  /* ********************************************************** */
  /*                                                            */
  /*                 cTriangulation3D                           */
  /*                                                            */
  /* ********************************************************** */



  struct ClassificationTags
  {
    const int8_t Unclassified=1;
    const int8_t Ground=2;
    const int8_t Low_Vegetation=3;
    const int8_t Medium_Vegetation=4;
    const int8_t High_Vegetation=5;
    const int8_t Building=6;
    const int8_t Water=9;
    const std::string DSMMarker="dsm_marker";
    const std::string DTMMarker="dtm_marker";
  };

  enum class eLabelIm_MASQ : tU_INT1
  {
     eFree,     // Mode MicMac V1
     eReached,  // Mode filled
     eNbVals
  };


template <class Type>  cTriangulation3DLas<Type>::cTriangulation3DLas(const std::string & aName):
        cTriangulation<Type,3>(std::vector<tPt>())
{
    if (UCaseEqual(LastPostfix(aName),"las") || UCaseEqual(LastPostfix(aName),"laz") )
    {
       LasInit(aName);
    }
    else
    {
       MMVII_UserError(eTyUEr::eBadPostfix,"Unknown postfix in cTriangulation3D");
    }
}

template <class Type> unsigned int cTriangulation3DLas<Type>::NbPts()
  {
    return this->mVPts.size();
  }

template <class Type> const char* cTriangulation3DLas<Type>::ProjStr()
  {
    return mProjStr;
  }

template <class Type> bool cTriangulation3DLas<Type>::HasTime()
  {
    return mHasTime;
  }

template <class Type> bool cTriangulation3DLas<Type>::HasColor()
  {
    return mHasColor;
  }



template <class Type> void cTriangulation3DLas<Type>::SamplePts(const bool & targetted,const Type & aStep)
  {
    /// < Sample points either by targetted or by random sampling
    mVSelectedIds.clear();
   if (targetted)
     {
       // sample points in a grid

       /*
        * |----- |  -----  | -------|
        * |----- |  -----  | -------|
        * |----- |  -----  | -------|
        */
       // Random Ordering of points
       //this->mVPts=RandomOrder(this->mVPts);

       // Empty grid of points to sample with size bbox/aStep
       cDataTypedIm<tU_INT1,2> aD_Grid(cPt2di(0,0),
                                     cPt2di(Pt_round_down(mDelimitBox.Sz()/aStep))
                                     );
       aD_Grid.InitCste(0);

       std::cout<<"Grid Size "<<aD_Grid.Sz()<<std::endl;
       // fill grid
       size_t it=0;
       int AreallCellsFilled=aD_Grid.NbElem();
       std::cout<<"Nb Elem "<<AreallCellsFilled<<std::endl;
       std::cout<<mDelimitBox.P0()<<"  "<<mDelimitBox.P1()<<std::endl;
       while((it<this->mVPts.size()) && AreallCellsFilled)
         {
            tPt aPt=this->mVPts.at(it);
            //std::cout<<aPt<<std::endl;
            cPt2di aPix=cPt2di((aPt.x()-mDelimitBox.P0().x())/aStep,
                               (aPt.y()-mDelimitBox.P0().y())/aStep);
            //std::cout<<aPix<<std::endl;
            if(aD_Grid.Inside(aPix))
              {
                if (aD_Grid.VI_GetV(aPix)==tU_INT1(eLabelIm_MASQ::eFree))
                  {
                    //std::cout<<"  "<<" cell "<<allCellsFilled<<" c"<<eLabelIm_MASQ::eFree<<std::endl;
                    aD_Grid.VI_SetV(aPix,1);
                    mVSelectedIds.push_back(it);
                    AreallCellsFilled--;
                  }
              }
            it++;
         }
     }
   else
     {
     }
}

  // ============================================================
  //  Constants for LAS Point Data Format IDs
  // ============================================================
  static constexpr int kLasFmt_HasTime[]  = {1, 3, 5, 6, 7, 8, 9, 10};
  static constexpr int kLasFmt_HasColor[] = {2, 3, 5, 7, 8, 10};

  // ============================================================
  //  Constants for LAS VLR record IDs and user IDs
  // ============================================================
  static constexpr U16        kVlrRecordId_WKT        = 2112;
  static constexpr U16        kVlrRecordId_ExtraBytes  = 4;
  static constexpr const char kVlrUserId_Projection[] = "LASF_Projection";
  static constexpr const char kVlrUserId_Spec[]       = "LASF_Spec";

  // ============================================================
  //  Constants for Extra Bytes types (LAS 1.4 spec Table 24)
  // ============================================================
  static constexpr U8  kExtraByteType_UInt8  = 1;
  static constexpr U8  kExtraByteType_Int8   = 2;
  static constexpr U8  kExtraByteType_UInt16 = 3;
  static constexpr U8  kExtraByteType_Int16  = 4;
  static constexpr U8  kExtraByteType_UInt32 = 5;
  static constexpr U8  kExtraByteType_Int32  = 6;
  static constexpr U8  kExtraByteType_UInt64 = 7;
  static constexpr U8  kExtraByteType_Int64  = 8;
  static constexpr U8  kExtraByteType_Float  = 9;
  static constexpr U8  kExtraByteType_Double = 10;

  // Size in bytes for each Extra Byte type (index = type id)
  static constexpr int kExtraByteTypeSizes[] = {1, 1, 1, 2, 2, 4, 4, 8, 8, 4, 8};

  // Size of one Extra Bytes VLR record entry (LAS 1.4 spec)
  static constexpr int kExtraByteRecordSize  = 192;
  // Offset of the 'type' field within one Extra Bytes record entry
  static constexpr int kExtraByteTypeOffset  = 2;
  // Offset of the 'name' field within one Extra Bytes record entry
  static constexpr int kExtraByteNameOffset  = 4;
  // Max length of the name field
  static constexpr int kExtraByteNameLength  = 32;

  // ============================================================
  //  Structure holding DSM marker extra-byte metadata
  // ============================================================
  struct cDsmMarkerInfo
  {
      bool mHasMarker  = false;
      int  mByteOffset = -1;  ///< Byte offset within point extra_bytes buffer
      int  mByteSize   =  0;  ///< Size in bytes of the marker field
  };

  // ============================================================
  //  Compute byte size of an Extra Bytes field from its type id
  // ============================================================
  static int ExtraByteSize(U8 aType)
  {
      if (aType < sizeof(kExtraByteTypeSizes) / sizeof(kExtraByteTypeSizes[0]))
          return kExtraByteTypeSizes[aType];
      return 1; // fallback: 1 byte for unknown types
  }

  // ============================================================
  //  Read the WKT projection string from LAS VLRs
  // ============================================================
  static std::string ReadProjStr(const LASheader & aHeader)
  {
      for (int i = 0; i < aHeader.number_of_variable_length_records; i++)
      {
          const auto & vlr = aHeader.vlrs[i];
          if (std::string(vlr.user_id) == kVlrUserId_Projection
              && vlr.record_id == kVlrRecordId_WKT)
          {
              return std::string(reinterpret_cast<const char*>(vlr.data),
                                 vlr.record_length_after_header);
          }
      }
      return "";
  }

  // ============================================================
  //  Parse VLRs to find DSM marker extra-byte offset and size
  // ============================================================
  static cDsmMarkerInfo ReadDsmMarkerInfo(const LASheader & aHeader)
  {
      cDsmMarkerInfo aInfo;

      for (int i = 0; i < aHeader.number_of_variable_length_records; i++)
      {
          const auto & vlr = aHeader.vlrs[i];
          if (std::string(vlr.user_id) != kVlrUserId_Spec
              || vlr.record_id != kVlrRecordId_ExtraBytes)
              continue;

          int nExtra      = vlr.record_length_after_header / kExtraByteRecordSize;
          int aByteOffset = 0;

          for (int e = 0; e < nExtra; e++)
          {
              const U8* aRec  = vlr.data + e * kExtraByteRecordSize;
              U8        aType = aRec[kExtraByteTypeOffset];

              std::string aName(reinterpret_cast<const char*>(aRec + kExtraByteNameOffset),
                                kExtraByteNameLength);
              aName = aName.c_str(); // trim trailing nulls

              int aSize = ExtraByteSize(aType);

              if (aName == ClassificationTags().DSMMarker)
              {
                  aInfo.mHasMarker  = true;
                  aInfo.mByteOffset = aByteOffset;
                  aInfo.mByteSize   = aSize;
                  return aInfo;
              }
              aByteOffset += aSize;
          }
          break; // only one Extra Bytes VLR expected
      }
      return aInfo;
  }

  // ============================================================
  //  Read the DSM marker value for a single point
  // ============================================================
  static bool IsForDSM(const LASpoint & aPoint, const cDsmMarkerInfo & aInfo)
  {
      if (!aInfo.mHasMarker || aInfo.mByteOffset < 0 || aPoint.extra_bytes == nullptr)
          return false;

      const U8* eb = aPoint.extra_bytes + aInfo.mByteOffset;

      switch (aInfo.mByteSize)
      {
      case 1:  return (int)(*eb)        != 0;
      case 2:  return (int)(*(U16*)eb)  != 0;
      case 4:  return (int)(*(I32*)eb)  != 0;
      default: return false;
      }
  }

  // ============================================================
  //  Main LAS initialisation function
  // ============================================================
  template <class Type>
  void cTriangulation3DLas<Type>::LasInit(const std::string & aNameFile)
  {
      LASreadOpener lasreadopener;
      lasreadopener.set_file_name(aNameFile.c_str());

      LASreader* lasreader = lasreadopener.open();
      if (!lasreader)
      {
          std::cerr << "ERROR: cannot open file " << aNameFile << std::endl;
          return;
      }

      const LASheader & aHeader = lasreader->header;

      std::cout << "POINT VIEW SIZE " << aHeader.number_of_point_records << std::endl;

      // ---- Header metadata ----
      {
          mNbPoints = aHeader.number_of_point_records;
          // LAS 1.4 may store counts > 2^32 in the extended field
          if (mNbPoints == 0 && aHeader.extended_number_of_point_records > 0)
              mNbPoints = (I64)aHeader.extended_number_of_point_records;

          mProjStr = ReadProjStr(aHeader);

          // Detect GPS time and RGB colour from point data format
          int fmt   = aHeader.point_data_format;
          mHasTime  = (fmt==1||fmt==3||fmt==5||fmt==6||fmt==7||fmt==8||fmt==9||fmt==10);
          mHasColor = (fmt==2||fmt==3||fmt==5||fmt==7||fmt==8||fmt==10);

          mDelimitBox = cTplBox<Type,2>(
              cPtxd<Type,2>(aHeader.min_x, aHeader.min_y),
              cPtxd<Type,2>(aHeader.max_x, aHeader.max_y));
      }

      // ---- DSM marker extra-byte descriptor ----
      const cDsmMarkerInfo aDsmInfo = ReadDsmMarkerInfo(aHeader);

      // ---- Read points ----
      if (aDsmInfo.mHasMarker)
      {
          std::cout << "HAS DSM MARKER" << std::endl;
          while (lasreader->read_point())
          {
              const LASpoint & p       = lasreader->point;
              int              Classif = (int)p.get_classification();

              bool IsVeg = ((Classif == ClassificationTags().Low_Vegetation)   ||
                            (Classif == ClassificationTags().Medium_Vegetation) ||
                            (Classif == ClassificationTags().High_Vegetation));

              if (IsForDSM(p, aDsmInfo) && !IsVeg)
              {
                  tPt aP(p.get_x(), p.get_y(), p.get_z());
                  this->mVPts.push_back(aP);
              }
          }
      }
      else
      {
          // No DSM marker: fall back to classification-based filtering
          // Keep Ground, Buildings, and Unclassified points
          while (lasreader->read_point())
          {
              const LASpoint & p       = lasreader->point;
              int              Classif = (int)p.get_classification();

              bool IsBuilding     = (Classif == ClassificationTags().Building);
              bool IsGround       = (Classif == ClassificationTags().Ground);
              bool IsUnclassified = (Classif == ClassificationTags().Unclassified);

              if (IsBuilding || IsGround || IsUnclassified)
              {
                  tPt aP(p.get_x(), p.get_y(), p.get_z());
                  this->mVPts.push_back(aP);
              }
          }
      }

      lasreader->close();
      delete lasreader;
  }

#if 0
template <class Type> void cTriangulation3DLas<Type>::LasInit(const std::string & aNameFile)
{
  pdal::Option las_opt("filename", aNameFile);
  pdal::Options las_opts;
  las_opts.add(las_opt);
  pdal::PointTable table;
  pdal::LasReader las_reader;
  las_reader.setOptions(las_opts);
  las_reader.prepare(table);
  pdal::PointViewSet point_view_set = las_reader.execute(table);
  pdal::PointViewPtr point_view = *point_view_set.begin();
  //pdal::Dimension::IdList dims = point_view->dims();
  pdal::LasHeader las_header = las_reader.header();

  std::cout<<"POINT VIEW SIZE "<<point_view->size()<<std::endl;

  {
    /* Point cloud properties */

    mNbPoints= las_header.pointCount();
    /** @brief spatial reference
     */
    mProjStr = table.spatialReference().getWKT().c_str();

    mHasTime = las_header.hasTime();
    mHasColor = las_header.hasColor();
    //pdal::Dimension::IdList dims = point_view->dims()
    //bounding box
    mDelimitBox=cTplBox<Type,2>(cPtxd<Type,2>(las_header.minX(),las_header.minY()),
                        cPtxd<Type,2>(las_header.maxX(),las_header.maxY()));
  }
  // Read Points and classification


  //std::cout<<" GET ALL DIMENSIONS "<<table.layout()->dims()<<std::endl;
  auto aDsmMarkerDim = table.layout()->findProprietaryDim(ClassificationTags().DSMMarker);
  bool HasDsmMarker=table.layout()->hasDim(aDsmMarkerDim);

  if (HasDsmMarker) // read points tagged as useful for DSM generation and not on trees
    {
      std::cout<<"HAS DSM MARKER "<<std::endl;
      for (pdal::PointId idx = 0; idx < point_view->size(); ++idx)
      {
         using namespace pdal::Dimension;
         auto IsForDsm=point_view->getFieldAs<int>(aDsmMarkerDim,idx);
         auto Classif=point_view->getFieldAs<int>(Id::Classification,idx);
         bool IsVeg=((Classif==ClassificationTags().Low_Vegetation) ||
                     (Classif==ClassificationTags().Medium_Vegetation) ||
                     (Classif==ClassificationTags().High_Vegetation)
                     );

         if (IsForDsm && !IsVeg)
           {
             tPt aP(point_view->getFieldAs<tREAL8>(Id::X, idx),
                    point_view->getFieldAs<tREAL8>(Id::Y, idx),
                    point_view->getFieldAs<tREAL8>(Id::Z, idx));
             this->mVPts.push_back(aP);
           }
      }
    }
  else  // assume point cloud is classified -> if there is not a tag dsm marker get points in GROUND, BUILDINGS
    {

        using namespace pdal::Dimension;

      for (pdal::PointId idx = 0; idx < point_view->size(); ++idx)
      {
        auto Classif=point_view->getFieldAs<int>(Id::Classification, idx);
        bool IsBuilding=(Classif==ClassificationTags().Building);
        bool IsGround=(Classif==ClassificationTags().Ground);
        bool IsUnclassified=(Classif==ClassificationTags().Unclassified);

        if ( IsBuilding || IsGround || IsUnclassified)
          {
            tPt aP(point_view->getFieldAs<tREAL8>(Id::X, idx),
                   point_view->getFieldAs<tREAL8>(Id::Y, idx),
                   point_view->getFieldAs<tREAL8>(Id::Z, idx));
            this->mVPts.push_back(aP);
          }
       }
    }

  // Read faces
  /*{
      std::vector<std::vector<size_t>> aVFace =   aLazF.getFaceIndices<size_t>();
      for (const auto & aFace : aVFace)
      {
          MMVII_INTERNAL_ASSERT_tiny(aFace.size()==3,"Bad face");
          this->AddFace(cPt3di(aFace[0],aFace[1],aFace[2]));
      }
  }*/

}

#endif
/* ========================== */
/*     INSTANTIATION          */
/* ========================== */

#define INSTANTIATE_TRI3DLAS(TYPE)\
template class cTriangulation3DLas<TYPE>;

INSTANTIATE_TRI3DLAS(tREAL4)
INSTANTIATE_TRI3DLAS(tREAL8)
INSTANTIATE_TRI3DLAS(tREAL16)

};
