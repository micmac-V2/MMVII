#if MMVII_USE_LIBTORCH
// Tool for calculating disparity between two Tiles using a CNN Trained Model

#include <torch/torch.h>
#include <torch/script.h>
#include <ATen/ATen.h>
#include "cCnnModelPredictor.h"

namespace F = torch::nn::functional;

namespace MMVII
{


/**************************************************************************/
aCnnModelPredictor::aCnnModelPredictor(std::string anArchitecture, std::string aModelBinDir, bool Cuda):
    mArchitecture(anArchitecture),IsCuda(Cuda)
{
    mSetModelBinaries = ToVect(SetNameFromPat(aModelBinDir,true));
}

/***********************************************************************/
void aCnnModelPredictor::PopulateModelMSNetHead(/*MSNetHead Network*/ torch::jit::script::Module & Network)
{
    //StdOut()<<"TO LOAD MODEL "<<"\n";
    std::string aModel= mSetModelBinaries.at(0);
    Network=torch::jit::load(aModel);
    auto cuda_available = torch::cuda::is_available();
    torch::Device device(cuda_available ? torch::kCUDA : torch::kCPU);
    Network.to(device);
    StdOut()<<"TORCH LOAD  "<<"\n";
}
/***********************************************************************/
void aCnnModelPredictor::PopulateModelFeatures(torch::jit::script::Module & Network)
{
#ifdef _WIN32
  if (IsCuda) LoadLibraryA("torch_cuda.dll");
#endif

    // add a convention on Model Name TAKE FOR EXAMPLES FEATURES AS A KEY FOR THE FEATURE MODULE
    std::string aModel;
    for (unsigned int i=0;i<mSetModelBinaries.size();i++)
    {
        if (mSetModelBinaries.at(i).find("FEATURES") != std::string::npos)
        {
            aModel=mSetModelBinaries.at(i);
            std::cout<<"Models checked "<<mSetModelBinaries.at(i)<<std::endl;
            break;
        }
    }
    StdOut()<<"Model Name "<<aModel<<"\n";
    torch::Device device(IsCuda ? torch::kCUDA : torch::kCPU);
    Network=torch::jit::load(aModel);
    Network.to(device);
    StdOut()<<"Loaded Model Feature Learning !!!!!! "<<"\n";
}

/***********************************************************************/
void aCnnModelPredictor::PopulateModelFeatures(torch::jit::script::Module & Network,bool DeviceCuda)
{
    // add a convention on Model Name TAKE FOR EXAMPLES FEATURES AS A KEY FOR THE FEATURE MODULE
    std::string aModel;
    for (unsigned int i=0;i<mSetModelBinaries.size();i++)
    {
        if (mSetModelBinaries.at(i).find("FEATURES") != std::string::npos)
        {
            aModel=mSetModelBinaries.at(i);
            std::cout<<"Models checked "<<mSetModelBinaries.at(i)<<std::endl;
            break;
        }
    }
    //<$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$>
    //auto cuda_available = WhichDevice;  //torch::cuda::is_available();
    //<$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$>
    StdOut()<<"Model Name "<<aModel<<"\n";
    //auto cuda_available=false;
    torch::Device device(DeviceCuda ? torch::kCUDA : torch::kCPU);
    try
    {
            Network= torch::jit::load(aModel);
    }
    catch (std::exception& e)
    {
            std::cout << e.what() << std::endl;
    }
    //Network=torch::jit::load(aModel);
    Network.to(device);
    StdOut()<<"MODEL FEATURES LOADED !!!!!! "<<"\n";
}
/***************************************************************************************/

void aCnnModelPredictor::PopulateModelDecision(torch::jit::script::Module & Network)
{
    // add a convention on Model Name TAKE FOR EXAMPLES FEATURES AS A KEY FOR THE FEATURE MODULE
    std::string aModel;
    for (unsigned int i=0;i<mSetModelBinaries.size();i++)
    {
        if (mSetModelBinaries.at(i).find("DECISION_NET") != std::string::npos)
        {
            std::cout<<"Models checked for the decision NEtwork"<<mSetModelBinaries.at(i)<<std::endl;
            aModel=mSetModelBinaries.at(i);
            break;
        }
    }
    //auto cuda_available = torch::cuda::is_available();
    torch::Device device(IsCuda ? torch::kCUDA : torch::kCPU);
    Network=torch::jit::load(aModel);
    Network.to(device);
    StdOut()<<"MODEL DECISION LOADED !!  "<<"\n";
}

/***************************************************************************************/

void aCnnModelPredictor::PopulateModelDecision(torch::jit::script::Module & Network,bool DeviceCuda)
{
    // add a convention on Model Name TAKE FOR EXAMPLES FEATURES AS A KEY FOR THE FEATURE MODULE
    std::string aModel;
    for (unsigned int i=0;i<mSetModelBinaries.size();i++)
    {
        if (mSetModelBinaries.at(i).find("DECISION_NET") != std::string::npos)
        {
            std::cout<<"Models checked for the decision NEtwork"<<mSetModelBinaries.at(i)<<std::endl;
            aModel=mSetModelBinaries.at(i);
            break;
        }
    }
    torch::Device device(DeviceCuda ? torch::kCUDA : torch::kCPU);
    Network=torch::jit::load(aModel);
    Network.to(device);
    StdOut()<<"MODEL DECISION LOADED !!  "<<"\n";
}
/***************************************************************************************/

void aCnnModelPredictor::PopulateModelMatcher(torch::jit::script::Module & Network)
{
    // add a convention on Model Name TAKE FOR EXAMPLES FEATURES AS A KEY FOR THE FEATURE MODULE
    std::string aModel;
    for (unsigned int i=0;i<mSetModelBinaries.size();i++)
    {
        if (mSetModelBinaries.at(i).find("MATCHER_NET") != std::string::npos)
        {
            std::cout<<"Models checked for the decision NEtwork"<<mSetModelBinaries.at(i)<<std::endl;
            aModel=mSetModelBinaries.at(i);
            break;
        }
    }
    auto cuda_available = torch::cuda::is_available();
    torch::Device device(cuda_available ? torch::kCUDA : torch::kCPU);
    Network=torch::jit::load(aModel);
    Network.to(device);
    StdOut()<<"MODEL MATCHER LOADED !!  "<<"\n";
}

/**********************************************************************************************************************/
torch::Tensor aCnnModelPredictor::PredictUNetWDecision(torch::jit::script::Module mNet, 
                                                        std::vector<tTImV2> aMasterP,
                                                        std::vector<tTImV2> aPatchLV, 
                                                        cPt2di aPSz)
{
    auto cuda_available=false;
    torch::Device device(cuda_available ? torch::kCUDA : torch::kCPU);
    torch::NoGradGuard no_grad;
    mNet.eval();
    torch::Tensor aPAllMasters=torch::empty({(int) aMasterP.size(),aPSz.y(),aPSz.x()},
                                            torch::TensorOptions().dtype(torch::kFloat32));
    torch::Tensor aPAllSlaves=torch::empty({(int) aPatchLV.size(),aPSz.y(),aPSz.x()},
                                           torch::TensorOptions().dtype(torch::kFloat32));
    for (int cc=0;cc<(int) aMasterP.size();cc++)
    {
        tREAL4 ** mPatchLData=aMasterP.at(cc).DIm().ExtractRawData2D();
        torch::Tensor aPL=torch::from_blob((*mPatchLData), {1,aPSz.y(),aPSz.x()},
                                           torch::TensorOptions().dtype(torch::kFloat32));
        //normalize apl
        //auto std=aPL.std();
        //aPL=aPL.sub(aPL.mean());
        //aPL=aPL.div(std.add(1e-12));
        //std::cout<<"  PATCH CONTENT "<<aPL<<std::endl;
        aPL=aPL.div(255.0);
        aPL=(aPL.sub(0.4353755468)).div(0.19367880);
        aPAllMasters.index_put_({cc},aPL);
    }
    //StdOut()<<"master "<<aPAllMasters.sizes()<<"\n";
    for (int cc=0;cc<(int) aPatchLV.size();cc++)
    {
        tREAL4 ** mPatchLData=aPatchLV.at(cc).DIm().ExtractRawData2D();
        torch::Tensor aPL=torch::from_blob((*mPatchLData), {1,aPSz.y(),aPSz.x()},
                                           torch::TensorOptions().dtype(torch::kFloat32));
        //auto std=aPL.std();
        //aPL=aPL.sub(aPL.mean());
        //aPL=aPL.div(std.add(1e-12));
        aPL=aPL.div(255.0);
        aPL=(aPL.sub(0.4353755468)).div(0.19367880);
        aPAllSlaves.index_put_({cc},aPL);
    }
    auto aPAll=torch::cat({aPAllMasters.unsqueeze(0),aPAllSlaves.unsqueeze(0)},0).to(device); 
    //StdOut()<<"Patches "<<aPAll.sizes()<<"\n";
    torch::jit::IValue inp(aPAll);
    std::vector<torch::jit::IValue> allinp={inp};
    auto out=mNet.forward(allinp);
    auto output=out.toTensor().squeeze();
    return output.to(torch::kCPU);
}

/**********************************************************************************************************************/
torch::Tensor aCnnModelPredictor::PredictUnetFeaturesOnly(torch::jit::script::Module mNet,
                                                            std::vector<tTImV2> aPatchLV, 
                                                            cPt2di aPSz)
{
    torch::Device device(IsCuda ? torch::kCUDA : torch::kCPU);
    torch::NoGradGuard no_grad;
    mNet.eval();
    torch::Tensor aPAllSlaves=torch::empty({(int) aPatchLV.size(),aPSz.y(),aPSz.x()},
                                           torch::TensorOptions().dtype(torch::kFloat32)).to(device);
    for (int cc=0;cc<(int) aPatchLV.size();cc++)
    {
        tREAL4 ** mPatchLData=aPatchLV.at(cc).DIm().ExtractRawData2D();
        torch::Tensor aPL=torch::from_blob((*mPatchLData), {1,aPSz.y(),aPSz.x()},
                                           torch::TensorOptions().dtype(torch::kFloat32));
        aPL=((aPL.div(255.0)).mul(2.0)).sub(1.0);
        aPAllSlaves.index_put_({cc},aPL.to(device));
    }
    // rotate by 90°
    //aPAllSlaves=aPAllSlaves.rot90(1,{1,2});
    torch::jit::IValue inp(aPAllSlaves.unsqueeze(0));
    std::vector<torch::jit::IValue> allinp={inp};
    auto out=mNet.forward(allinp);
    auto output=out.toTensor().squeeze();
    return output;
}

/**********************************************************************************************************************/
torch::Tensor aCnnModelPredictor::PredictUnetFeaturesOnly(torch::jit::script::Module mNet,
                                                            torch::Tensor aPAllSlaves)
{
    //std::cout<<"Cuda is available ? "<<cuda_available<<std::endl;
    torch::Device device(IsCuda ? torch::kCUDA : torch::kCPU);
    torch::NoGradGuard no_grad;
    mNet.eval();
    aPAllSlaves=aPAllSlaves.to(device);
    torch::jit::IValue inp(aPAllSlaves.unsqueeze(0).unsqueeze(0));
    std::vector<torch::jit::IValue> allinp={inp};
    auto out=mNet.forward(allinp);
    auto output=out.toTensor().squeeze();
    return output;//.to(torch::kCPU);
}
/**********************************************************************************************************************/
torch::Tensor aCnnModelPredictor::PredictMSNetTile(torch::jit::script::Module mNet, tTImV2 aPatchLV, cPt2di aPSz)
{
    auto cuda_available = torch::cuda::is_available();
    //std::cout<<"Cuda is available ? "<<cuda_available<<std::endl;
    torch::Device device(cuda_available ? torch::kCPU : torch::kCPU);
    torch::NoGradGuard no_grad;
    mNet.eval();
    tREAL4 ** mPatchLData=aPatchLV.DIm().ExtractRawData2D();
    torch::Tensor aPL=torch::from_blob((*mPatchLData), {1,1,aPSz.y(),aPSz.x()},
                                         torch::TensorOptions().dtype(torch::kFloat32)).to(device);
    torch::jit::IValue inp(aPL);
    std::vector<torch::jit::IValue> allinp={inp};
    auto out=mNet.forward(allinp);
    auto output=out.toTensor();
    return output;
}
/**********************************************************************************************************************/
/**********************************************************************************************************************/
torch::Tensor aCnnModelPredictor::PredictMSNetTileFeatures(torch::jit::script::Module mNet, 
                                                            tTImV2 aPatchLV, 
                                                            cPt2di aPSz)
{
    #ifdef _WIN32
    if(IsCuda)
      {
        LoadLibrary("torch_cuda.dll");
      }
    #endif
    //auto cuda_available = torch::cuda::is_available();
    //std::cout<<"Cuda is available ? "<<cuda_available<<std::endl;
    torch::Device device(IsCuda ? torch::kCUDA : torch::kCPU);
    torch::NoGradGuard no_grad;
    mNet.eval();
    tREAL4 ** mPatchLData=aPatchLV.DIm().ExtractRawData2D();
    torch::Tensor aPL=torch::from_blob((*mPatchLData),
                                         {1,1,aPSz.y(),aPSz.x()},
                                         torch::TensorOptions().dtype(torch::kFloat32)).to(device);
    aPL=((aPL.div(255.0)).mul(2.0)).sub(1.0);//.sub(0.5);
    torch::jit::IValue inp(aPL);
    std::vector<torch::jit::IValue> allinp={inp};
    //std::cout<<"IVALUE CREATED "<<std::endl;
    auto out=mNet.forward(allinp);
    auto output=out.toTensor().squeeze();
    return output;  
}
/**********************************************************************************************************************/
/*********************************************************************************************************************/
torch::Tensor aCnnModelPredictor::PredictDecisionNet(torch::jit::script::Module mNet, 
                                                    torch::Tensor Left, 
                                                    torch::Tensor Right)
{
    auto cuda_available = torch::cuda::is_available();
    torch::Device device(cuda_available ? torch::kCUDA : torch::kCPU);
    torch::NoGradGuard no_grad;
    mNet.eval();
    auto CatTensor=torch::cat({Left,Right},1); // to get a size of {1,FeatsSIZE}
    //CatTensor=CatTensor.squeeze(0).unsqueeze(3);
    torch::jit::IValue inp(CatTensor);
    std::vector<torch::jit::IValue> allinp={inp};
    torch::Tensor OutSim=mNet.forward(allinp).toTensor().squeeze();
    return torch::sigmoid(OutSim).to(torch::kCPU);
}

/*********************************************************************************************************************/
torch::Tensor aCnnModelPredictor::PredictONCUBE(torch::jit::script::Module mMlp, torch::Tensor & aCube)
{
    //torch::Device device(torch::kCUDA);
    torch::NoGradGuard no_grad;
    mMlp.eval();
    torch::jit::IValue inp(aCube);
    std::vector<torch::jit::IValue> allinp={inp};
    torch::Tensor OutSimBrut=mMlp.forward(allinp).toTensor().sigmoid();
    return OutSimBrut.to(torch::kCPU);
}
/**********************************************************************************************************************/
torch::Tensor aCnnModelPredictor::PredictMSNetAtt(MSNet_Attention mNet, std::vector<tTImV2> aPatchLV, cPt2di aPSz)
{
        torch::Device device(torch::kCPU);
        torch::NoGradGuard no_grad;
        mNet->eval();
    torch::Tensor aPAllScales=torch::empty({4,aPSz.y(),aPSz.x()}, torch::TensorOptions().dtype(torch::kFloat32));
    //std::cout<<"SIZE OF MSCALE TILES "<<aPatchLV.size()<<std::endl;
    for (int cc=0;cc<(int) aPatchLV.size();cc++)
    {
        StdOut()<<"Size of tile Mul Scale is "<<aPatchLV.at(cc).DIm().Sz()<<"\n";
        tREAL4 ** mPatchLData=aPatchLV.at(cc).DIm().ExtractRawData2D();
        torch::Tensor aPL=torch::from_blob((*mPatchLData), {1,aPSz.y(),aPSz.x()}, 
                                            torch::TensorOptions().dtype(torch::kFloat32));
        aPAllScales.index_put_({cc},aPL);
    }
    aPAllScales=aPAllScales.unsqueeze(0);

    auto output=mNet->forward(aPAllScales).squeeze();
    return output;
}
/**********************************************************************************************************************/
torch::Tensor aCnnModelPredictor::PredictMSNetHead(torch::jit::script::Module mNet, 
                                                std::vector<tTImV2> aPatchLV, 
                                                cPt2di aPSz)
{
    auto cuda_available = torch::cuda::is_available();
    std::cout<<"Cuda is available ? "<<cuda_available<<std::endl;
    torch::Device device(cuda_available ? torch::kCUDA : torch::kCPU);
        torch::NoGradGuard no_grad;
        mNet.eval();
    torch::Tensor aPAllScales=torch::empty({(int) aPatchLV.size(),aPSz.y(),aPSz.x()},
                                         torch::TensorOptions().dtype(torch::kFloat32));;
    for (int cc=0;cc<(int) aPatchLV.size();cc++)
    {
        tREAL4 ** mPatchLData=aPatchLV.at(cc).DIm().ExtractRawData2D();
        torch::Tensor aPL=torch::from_blob((*mPatchLData), {1,aPSz.y(),aPSz.x()}, 
                                            torch::TensorOptions().dtype(torch::kFloat32));
        aPAllScales.index_put_({cc},aPL);
    }
    aPAllScales=aPAllScales.unsqueeze(0).to(device);
    StdOut()<<"Patches "<<aPAllScales.sizes()<<"\n";
    torch::jit::IValue inp(aPAllScales);
    std::vector<torch::jit::IValue> allinp={inp};
    auto out=mNet.forward(allinp);
    auto output=out.toTensor().squeeze();
    return output.to(torch::kCPU);
}

};
#endif
