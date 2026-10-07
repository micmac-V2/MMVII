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


/*



MMVII ImageGenRandom  [3000,2000]   GaussNoise=[100,200] GaussNoise=[5,10]  Out="/home/MPierrot-Deseilligny/MMVII/MMVII-TestDir/Tmp/ImSynt_1.tif" SeedRand=1
MMVII ImageGenRandom  [3000,2000]   GaussNoise=[100,200] GaussNoise=[5,10]  Out="/home/MPierrot-Deseilligny/MMVII/MMVII-TestDir/Tmp/ImSynt_2.tif" SeedRand=2



MMVII CodedTargetGenerateEncoding IGNIndoor 14 Out=/home/MPierrot-Deseilligny/MMVII/MMVII-TestDir/Tmp/Encoding.xml

MMVII CodedTargetGenerate  /home/MPierrot-Deseilligny/MMVII/MMVII-TestDir/Tmp/Encoding.xml Out=/home/MPierrot-Deseilligny/MMVII/MMVII-TestDir/Tmp/FullSpec.xml


MMVII CodedTargetSimul /home/MPierrot-Deseilligny/MMVII/MMVII-TestDir/Tmp/ImSynt_1.tif  /home/MPierrot-Deseilligny/MMVII/MMVII-TestDir/Tmp/FullSpec.xml


  /home/MPierrot-Deseilligny/MMVII/MMVII-TestDir/Tmp/Encoding.xml Out=/home/MPierrot-Deseilligny/MMVII/MMVII-TestDir/Tmp/FullSpec.xml




*/

void BenchCodedTarget(cParamExeBench & aParam)
{
    if (! aParam.NewBench("CodedTarget")) return;

    std::string aDir = cMMVII_Appli::TmpDirTestMMVII() ;
//    std::string aDir = "/home/MPierrot-Deseilligny/MMVII/MMVII-TestDir/Tmp/";
    StdOut() << " TMP=" << aDir << "\n";


    ScopedChdir aScopeChD(aDir);

    cIm2D<tU_INT1> aIm(cPt2di(200,200));
    aIm.DIm().ToFile("toto.tif");
    getchar();


    int aNbIm = 3;
    cPt2di aSzIm(3000,2000);
    for (int aKIm=0 ; aKIm<aNbIm ; aKIm++)
    {
        cParamCallSys aCom(TheSpec_cAppliGenRandomImage,true);
        aCom.AddArgs(ToStr(aSzIm));

        aCom.AddMMVIIArgsOpt(CurOP_Out,"ImSynt_"+ToStr(aKIm) + ".tif");
        std::vector<std::vector<double>> aVecNoise{{100.0,200.0},{5.0,10.0}};

        for (const auto & aGaussN : aVecNoise)
            aCom.AddMMVIIArgsOpt("GaussNoise",ToStr(aGaussN));

        ///aCom.AddArgs("MMVII");

        ///aCom.AddArgs(TheSpec_cAppliGenRandomImage.Name());


       // aCom.AddArgs("GaussNoise=[100,200]");
       //  aCom.AddArgs("GaussNoise=[5,10]");


       //  StdOut() << "COM=" << aCom.Com() << "\n";  getchar();

        aCom.Execute(true);
     }



    /*
    std::string aDirTmp = "/home/MPierrot-Deseilligny/MMVII/MMVII-TestDir/Tmp/";
    ScopedChdir aSChd(aDirTmp);


    std::string aCmd1 = "MMVII ImageGenRandom  [3000,2000]   GaussNoise=[100,200] GaussNoise=[5,10]  Out=ImSynt_1.tif SeedRand=1";*/

    aParam.EndBench();
}





};
