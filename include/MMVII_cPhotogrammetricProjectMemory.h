#ifndef  _MMVII_cPhotogrammetricProjectMemory_H_
#define  _MMVII_cPhotogrammetricProjectMemory_H_

#include "MMVII_Sensor.h"


using namespace NS_SymbolicDerivative;

namespace MMVII
{

/// Mode for handling intrinsic calibration in cPhotogrammetricProjectMemory
enum class eModeBenchPhMI
           {
              eDupl,      ///<  case duplication of each intrinsic camera
          //     eTL,        ///< case we use a thread local static map
              eTLS        ///< case thread local static

           };


/**  In-memory implementation of cIPhProj.
 *   Calibrations, orientations and homologous points are stored
 *   in maps keyed by image name rather than read from disk.
 *
 *   Populate with Add*() methods before use; the Read*() methods then return
 *   the stored objects.  The Save*() methods are not meaningful in memory mode
 *   and will trigger an error if called.
 */
class cPhotogrammetricProjectMemory : public cIPhProj
{
    public :

        cPhotogrammetricProjectMemory();
        ~cPhotogrammetricProjectMemory();

        // === Population methods ===

        /// Register calibration for a given image name (does NOT take ownership)
        void  AddCalib(const std::string & aNameIm, cPerspCamIntrCalib *);
        /// Register sensor for a given image name (does NOT take ownership)
   //  Unused ?     void  AddSensor(const std::string & aNameIm, cSensorCamPC *);
        /// Store homologous points for an image pair (copied)
        void  AddHomol(const std::string & aNameIm1, const std::string & aNameIm2,
                       const cSetHomogCpleIm &);
        /// Store multiple tie-points for one image (copied)
        void  AddMulTieP(const std::string & aNameIm, const cVecTiePMul &);

        // === cIPhProj interface ===

        cPerspCamIntrCalib *  InternalCalibFromStdName(const std::string aNameIm, bool isRemanent=true) const override;
        cPerspCamIntrCalib *  InternalCalibFromImage(const std::string & aNameIm) const override;

        cSensorImage *  ReadSensor(const std::string & aNameIm, bool ToDeleteAutom, bool SVP=false) const override;
        cSensorCamPC *  ReadCamPC(const std::string &, bool ToDeleteAutom, bool SVP=false) const override;
        tPoseR          ReadPoseCamPC(const std::string & aNameIm, bool * SVP=nullptr) const override;
        void            SaveSensor(const cSensorImage &) const override;
        void            SaveCamPC(const cSensorCamPC &) const override;
        void            SaveCalibPC(const cPerspCamIntrCalib &) const override;

        void  ReadHomol(cSetHomogCpleIm &,
                        const std::string & aNameIm1,
                        const std::string & aNameIm2,
                        const std::string & aDir="") const override;
        void  SaveHomol(const cSetHomogCpleIm &,
                        const std::string & aNameIm1,
                        const std::string & aNameIm2,
                        const std::string & aDir="") const override;

        void  ReadMultipleTiePFromFolder(const std::string & aFolder,
                                         cVecTiePMul &,
                                         const std::string & aNameIm,
                                         bool SVP=false) const override;
        std::string  MulTiePDirIn() const override;  ///< Returns ""

        const std::map<std::string, cSensorCamPC *> & SensorMap() const { return mSensorMap; }  ///< Accessor

    private :
        cPhotogrammetricProjectMemory(const cPhotogrammetricProjectMemory &) = delete;

        std::map<std::string, cPerspCamIntrCalib *>                           mGLOB_CalibMap;
        // thread_local   std::map<std::string, cPerspCamIntrCalib *>            mTL_CalibMap;
        thread_local static   std::map<std::string, cPerspCamIntrCalib *>     mTLS_CalibMap;

        eModeBenchPhMI                                                 mMode;
        mutable std::map<std::string, cSensorCamPC *>                  mSensorMap;
        mutable std::map<std::string, cSensorCamPC *>                  mOwnedSensorMap;  ///< owns sensors written via SaveCamPC (destructor deletes)
        std::map<std::pair<std::string,std::string>, cSetHomogCpleIm>  mHomolMap;
        std::map<std::string, cVecTiePMul>                             mMulTiePMap;
};


};

#endif  //  _MMVII_SENSOR_H_
