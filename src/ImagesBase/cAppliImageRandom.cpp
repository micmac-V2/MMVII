#include "MMVII_Tpl_Images.h"
#include "MMVII_Linear2DFiltering.h"
#include "MMVII_DeclareCste.h"
#include "MMVII_Geom2D.h"
#include "MMVII_Interpolators.h"

namespace MMVII
{

class cAppliGenRandomImage : public cMMVII_Appli
{
     public :
        cAppliGenRandomImage(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec);
     private :
        //================= typedef part ===============================


        typedef cIm2D<tU_INT1 >    tImLabel;
        typedef cDataIm2D<tU_INT1> tDImLabel;

        //================================================================
        //       METHODS DECLARATION
        //================================================================

            //------------------------- overidding cMMVII_Appli ----------------
        int Exe() override;
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override ;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override ;
        std::vector<cOneHelpSampleCmp>  Samples() const override;
        virtual ~cAppliGenRandomImage();
};  

}; // namespace MMVII

