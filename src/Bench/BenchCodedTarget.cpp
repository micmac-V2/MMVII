#include "MMVII_DeclareAllCmd.h"
#include "MMVII_DeclareCste.h"
#include "MMVII_Image2D.h"

#include <filesystem>


namespace MMVII
{
namespace fs = std::filesystem;
class ScopedChdir 
{
    fs::path old_;
public:
    explicit ScopedChdir(const fs::path& dir) : old_(fs::current_path()) {
        fs::current_path(dir);
    }
    ~ScopedChdir() {
        std::error_code ec;
        fs::current_path(old_, ec);          // pas d'exception dans un destructeur
    }
    ScopedChdir(const ScopedChdir&) = delete;
    ScopedChdir& operator=(const ScopedChdir&) = delete;
};


static const std::string  PrefixIm = "ImSynt_";

void BenchTargetGenerateImage(int aNbIm,const cPt2di&  aSzIm)
{
    for (int aKIm=0 ; aKIm<aNbIm ; aKIm++)
    {
        cParamCallSys aCom(TheSpec_cAppliGenRandomImage,true);
        aCom.AddArgs(ToStr(aSzIm));

        aCom.AddMMVIIArgsOpt("Type","U_INT1");
        aCom.AddMMVIIArgsOpt(CurOP_Out,PrefixIm+ToStr(aKIm) + ".tif");
        std::vector<std::vector<double>> aVecNoise{{100.0,70.0},{5.0,10.0}};

        for (const auto & aGaussN : aVecNoise)
            aCom.AddMMVIIArgsOpt("GaussNoise",ToStr(aGaussN));

        aCom.Execute(aKIm==0);
     }
}


std::string BenchTargetGenerate_Specif(eTyCodeTarget aType,int aNBB)
{

    //  ----------------- Generate the encoding -----------------------------------------------

      std::string aNameEncoding = "Encoding_" + E2Str(aType) + "_" + ToStr(aNBB) + ".xml";
      {
          cParamCallSys aComEncoding(TheSpecGenerateEncoding,true);
          aComEncoding.AddArgs(E2Str(aType));
          aComEncoding.AddArgs(ToStr(aNBB));
          aComEncoding.AddMMVIIArgsOpt(CurOP_Out,aNameEncoding);
          aComEncoding.Execute(false);
      }

      //  ----------------- Generate the full spec, with geometry  -----------------------------------------------

      std::string aNameSpec = "FullSpec_" + E2Str(aType) + "_" + ToStr(aNBB) + ".xml";
      {
          cParamCallSys aComGenerate(TheSpecGenCodedTarget,true);
          aComGenerate.AddArgs(aNameEncoding);
          aComGenerate.AddMMVIIArgsOpt(CurOP_Out,aNameSpec);
          // For (still) unexplained reason, there is a memory leak if we do  Execute(false)
          // (do it in the same process) by the way, the command if it goes untill the end is probably
          // memory correct as its has been widely used . So, for now, run in a separate process
          aComGenerate.Execute(true);
      }

      return aNameSpec;
}


void BenchTarget_GenerateSimul(const std::string & aNameFullSpec,const std::string & aPatNames)
{
    // MMVII CodedTargetSimul ImSynt_.*tif FullSpec_IGNIndoor_14.xml Radius=[30,60] Ratio=[0.5,1.0] NoiseAmpl=[0,0] PropLinBias=[0,0]

    cParamCallSys aComSimul(TheSpecSimulCodedTarget,true);

    aComSimul.AddArgs(PrefixIm+".*tif");
    aComSimul.AddArgs(aNameFullSpec);

    aComSimul.AddMMVIIArgsOpt("PatNames",aPatNames);

    aComSimul.AddMMVIIArgsOpt("Radius",ToStr(cPt2dr(40,60)));
    aComSimul.AddMMVIIArgsOpt("Ratio",ToStr(cPt2dr(0.7,1.0)));
    aComSimul.AddMMVIIArgsOpt("NoiseAmpl",ToStr(cPt2dr(0,0)));
    aComSimul.AddMMVIIArgsOpt("PropLinBias",ToStr(cPt2dr(0,0)));
    // aComSimul.AddMMVIIArgsOpt("Show","false");

   aComSimul.Execute(true);

    // StdOut() << " COM= " << aComSimul.Com() << "\n";
}



void BenchCodedTarget(cParamExeBench & aParam)
{
    if (! aParam.NewBench("CodedTarget")) return;

    std::string aDir = cMMVII_Appli::TmpDirTestMMVII() ;

    StdOut() << "PROFIL NAME=" << cMMVII_Appli::ProfileName() << "\n";
    if (UserIsMPD())
    {
        aDir = "/home/MPierrot-Deseilligny/MMVII/MMVII-TestDir/Tmp/";
        StdOut() << " DIR=" << aDir << "\n";
    }

    // push to the tmp directory
    ScopedChdir aScopeChD(aDir);


    if (true)
    {
       BenchTargetGenerateImage(3,cPt2di(3000,2000));
    }
    std::string aNameFullSpec = BenchTargetGenerate_Specif(eTyCodeTarget::eIGNIndoor,14);

    BenchTarget_GenerateSimul(aNameFullSpec,".*[05]");

    aParam.EndBench();
}





};
