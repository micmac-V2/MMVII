#include "cMMVII_Appli.h"
#include "MMVII_PCSens.h"
#include "MMVII_Geom2D.h"
#include "MMVII_Geom3D.h"
#include "MMVII_AllClassDeclare.h"
#include "MMVII_DeclareCste.h"
#include "MMVII_2Include_Serial_Tpl.h"
#include "MMVII_Triangles.h"
#include "MMVII_Image2D.h"
#include "MMVII_ZBuffer.h"
#include "MeshDev.h"
#include "MMVII_Sys.h"
#include "MMVII_Radiom.h"
#include "MMVII_CloudRaster.h"
#include "ogrsf_frmts.h"
#include "MMVII_Tpl_GraphStruct.h"
#include "MMVII_Tpl_Graph_SubGraph.h"
#include "ranges"

namespace MMVII
{

static std::string NameIndBoxRecal="INTERNAL_IndexBoxRecall";


// classes of Vertex Attribute, OrientedEdge Attribute and non-Oriented Edge Attribute
template <class Type> class cGG_AttrVert;
template <class Type> class cGG_AttrEdgeOr;
template <class Type> class cGG_AttrEdgeSym;
template <class Type> class cAlgo_ChambollePockParams;
template<class tGraph, class Type> class cAlgo_ChambollePock;

/* Class Vertex Attributes */ 
template <class Type> class cGG_AttrVert
{
    public:
        cGG_AttrVert(cPt2di aPt, Type aZ):
            mPt(aPt),
            mAbsCurv(aZ),
            mValOptInit(0.0),
            mValToOpt(0.0),
            mValDbar(0.0),
            mWeightsVert({}),
            mWeightPerNode(0.0)
            {
            }
        cPt2di mPt;
        Type mAbsCurv;
        Type mValOptInit; // initial value of the vertex in the optimization variable, used for debug
        Type mValToOpt; // value of the vertex in the optimization variable, used for debug
        Type mValDbar; // value of the vertex in the over-relaxation step, used for debug
        std::vector<Type> mWeightsVert;
        Type mWeightPerNode;

    private:
};

/* Class of Oriented Edge Attributes */
template <class Type>  class cGG_AttrEdgeOr
{
    public:
        cGG_AttrEdgeOr(Type aDeltaAc):
            mDeltaAC(aDeltaAc)
            {

            }
        Type mDeltaAC;
    private:
};

/* Class of NON Oriented (SYMMETRIC) Edge Attributes */
template <class Type> class cGG_AttrEdgeSym
{
    public:
        cGG_AttrEdgeSym(Type aDist) :
            mDist     (aDist) ,
            mRanCost  (RandUnif_C()),
            mIsOk     (true),
            mCptPath0 (0),
            mCptPath1 (0),
            mHCode0   (0),
            mHCode1   (0)
           {}
        Type mDist;      // euclidean distance, not really use 4 now
        Type mRanCost;   // use to modify the weighting
        bool   mIsOk;      // use to creat some specific subgraph
        int    mCptPath0;  // count the number of path
        int    mCptPath1;  // use to "back up" path counter for check
        size_t mHCode0;    // h-code of pathes going throuh it
        size_t mHCode1;    // h-code of pathes going throuh it
    private:
};



// class for Grid Graph where vertexes coincide with pixels 
/* this class will be used to merge multiple dems with confidence maps using 
the Chambolle-Pock algorithm */
template <class Type> class cAlgo_ChambollePockParams
{
    public:
        cAlgo_ChambollePockParams(Type aLambda, 
                                    Type aSigmaR, 
                                    int aNbIter, 
                                    Type aTolerance, 
                                    Type aRho) :
            mLambda(aLambda),
            mSigmaR(aSigmaR),
            mNbIter(aNbIter),
            mTolerance(aTolerance),
            mRho(aRho)          
        {}
        Type mLambda; // regularization parameter
        Type mSigmaR; // sigmas for edge weight, for example in exp(-|D_i^0 - D_j^0|^2 / (2*sigmaR^2))
        int mNbIter; // number of iterations for the optimization
        Type mTolerance; // tolerance for the optimization        
        Type mRho; // Step Size for the Chambolle-Pock algorithm, can be tuned for convergence speed

    private:
};  




template <class tGraph, class Type> class cAlgo_ChambollePock
{
    /* The functional to optimize has the following form
    
    E(D) = sum_i [ w_i * (D_i - D_i^0)^2 ] + lambda * sum_{i,j} [ w_{i,j} * ||D_i - D_j||^2 ] for the moment use quadratic regularization, 
    but it can be easily changed to other type of regularization (e.g. TV) by changing the second term 

    where: 
        - D is the variable to optimize
        - D^0 is the initial value (data term)
        - w_i is the confidence of the data term
        - w_{i,j} is the weight of the regularization term
        - lambda is the regularization parameter

    For digital elevation model fusion, 
        - D would be the fused DEM
        - D^0 would be the input DEMs
        - w_i would be the confidence maps
        - w_{i,j} would be the weights for the regularization term (e.g exp(-|D_i^0 - D_j^0|) to preserve edges)
    */
    public:
        //typedef typename cAlgo_ParamVG<tGraph> tParamA;
        /*typedef typename tGraph::tVertex tVertex;
        typedef typename tGraph::tEdge tEdge;
        typedef typename std::vector<tVertex*> tVecVert;*/

        void PrepareInitTerms(tGraph & aGraph,
                            const std::vector<cIm2D<Type>> & aVMes, 
                            const std::vector<cIm2D<Type>> & aVConf,
                            const cAlgo_ChambollePockParams<Type> & aParams)

        {
            // 3-  precompute effective weights of the functional:
            //  mWeightPerNode = sum_i w_i and mValToOpt = sum_i w_i * D_i^0 / sum_i w_i for each vertex
            for(const auto & aVertex : aGraph.AllVertices())
            {
                cPt2di aPt = aVertex->Attr().mPt;
                Type aSumConf = 0.0;
                Type aValInitToOpt = 0.0;
                for (size_t k=0; k<aVConf.size(); k++)                
                {   
                    Type conf = aVConf[k].DIm().GetV(aPt);
                    if( conf >= 0.0)
                    {
                        Type mes = aVMes[k].DIm().GetV(aPt);
                        MMVII_INTERNAL_ASSERT_tiny(mes>=0, "all heights must be positive, use a mask to ignore negative heights"+ToStr(mes));
                        aValInitToOpt += conf * mes; // use confidence as weight for data term
                        aSumConf += conf;
                    }
                }

                aVertex->Attr().mWeightPerNode = aSumConf; //w_total per node
                aVertex->Attr().mValToOpt = aValInitToOpt / (aSumConf+1e-6);  // D : value of the vertex in the optimization variable
                aVertex->Attr().mValDbar = aVertex->Attr().mValToOpt; // D̄ : value of the vertex in the over-relaxation step, initialized to D
                aVertex->Attr().mValOptInit = aValInitToOpt; // d_bar : initial value of the vertex in the optimization variable, used for debug
            }

            // 4- initialize edge weights for the regularization term, for example using the gradient of initialized Value for each Vertex of Node 
            for (const auto & anE : aGraph.AllEdges())
            {  
                // for example, use the difference of initial values as weight to preserve edges
                Type aDeltaVal = std::abs(anE->VertexInit().Attr().mValToOpt - anE->Succ().Attr().mValToOpt);
                Type weight = std::exp(-aDeltaVal*aDeltaVal / (2 * aParams.mSigmaR * aParams.mSigmaR)); 
                anE->AttrSym().mRanCost = weight; // use symmetric edge attribute to store weight
            }
        }

        void SolveChambollePock(tGraph & aGraph, 
                        Type & aHighestEigenVal,
                        cAlgo_ChambollePockParams<Type> & aParamsAlgo,
                        std::vector<cIm2D<Type>> & aVMes, 
                        std::vector<cIm2D<Type>> & aVConf)
        {
            // 1- first initialize the terms of the functional
            StdOut() << " PrepareInitTerms \n"; 
            PrepareInitTerms(aGraph, aVMes, aVConf,aParamsAlgo);
            StdOut() << " PrepareInitTerms done \n"; 
            //2- Estimate Step Sizes for Primal and Dual Updates
            //HighestEigenVal(aGraph, aHighestEigenVal);
            StdOut() << " Highest Eigen Value = " << aHighestEigenVal << "\n";
            // 3- Use ||K|| = sqrt(lambda_max(B^T W^2 B)) for Chambolle-Pock step sizes.
            Type aNormK = std::sqrt(std::max<Type>(aHighestEigenVal,0.0));
            Type tau = aParamsAlgo.mRho / (aNormK + 1e-6); // step size for primal update
            Type sigma = aParamsAlgo.mRho / (aNormK + 1e-6); // step size for dual update

            // 4 - Init Dual Var
            for (const auto & anEdge : aGraph.AllEdges())
            {
                anEdge->AttrSym().mDist = 0.0 ; // use mDist for DualVar p_e = 0.0   
            }
            // 5- Iteration loop
            StdOut() << " Start Iteration Loop \n";
            for(int anIter =0 ; anIter< aParamsAlgo.mNbIter; anIter++)
                {
                    if (anIter % 10 == 0)
                        StdOut() << " Iteration " << anIter << "\n";

                    //1- Dual Update on all Edges using Incidence Matrix
                    // (logging removed - was called every iteration)
                    //#pragma omp parallel for
                    for(const auto & anEdge : aGraph.AllEdges())
                    {
                        Type aVOptIn = anEdge->VertexInit().Attr().mValDbar; // D̄bar_u
                        Type aVOptSucc = anEdge->Succ().Attr().mValDbar; // D̄bar_v
                        Type aWghtEdg = anEdge->AttrSym().mRanCost;

                        // p_e  +  σ · w_e · (D̄_u − D̄_v)
                        Type anUpdate_P_E = anEdge->AttrSym().mDist + sigma * aWghtEdg * (aVOptIn - aVOptSucc);
                        // clamp between −λ·w_e,  +λ·w_e
                        anEdge->AttrSym().mDist = std::clamp(anUpdate_P_E, - aParamsAlgo.mLambda * aWghtEdg,aParamsAlgo.mLambda * aWghtEdg);
                    }

                    // 2-  Primal Update on all Vertices
                    // for each node/vertex  
                    // div_u  =  Σ_{e leaving u} p_e · w_e  −  Σ_{e entering u} p_e · w_e from mIncidenceMatrix
                    // z_u    =  D_u  −  τ · div_u
                    // D_new_u    =  (z_u  +  τ · (sum_k w_k·D_k^0)) / (1 + τ · w_u) where w_u is total confidence weight and D_k^0 are initial values

                    Type aSqChange = 0.0;
                    Type aSqD = 0.0;

                    // Primal Update on all Vertices (logging removed)
                    //#pragma omp parallel for reduction(+:aSqChange,aSqD)
                    for(const auto & aVertex : aGraph.AllVertices())
                    {
                        Type w_u = aVertex->Attr().mWeightPerNode; // total confidence weight for the data term at vertex u
                        // compute div_u using incidence matrix
                        Type div_u = 0.0;
                        // compute div_u using edge weights
                        for (const auto & anEdge : aVertex->EdgesSucc()) // for each edge incident to u
                        {
                            Type w_e = anEdge->AttrSym().mRanCost;
                            if( anEdge->VertexInit().Attr().mAbsCurv == aVertex->Attr().mAbsCurv) // edge leaving u
                            {
                                div_u -= anEdge->AttrSym().mDist * w_e; // p_e · w_e
                            }
                            else if (anEdge->Succ().Attr().mAbsCurv == aVertex->Attr().mAbsCurv) // edge entering u
                            {
                                div_u += anEdge->AttrSym().mDist * w_e; // - p_e · w_e
                            } 
                        }
                        // compute z_u  
                        Type z_u = aVertex->Attr().mValToOpt - tau * div_u; // z_u = D_u - τ · div_u
                        // D_new_u = (z_u + τ · (sum_k w_k·D_k^0)) / (1 + τ · w_u)
                        Type D_new_u = (z_u + tau * aVertex->Attr().mValOptInit) / (1 + tau * w_u); 

                        // OVER-RELAXATION

                        aVertex->Attr().mValDbar = 2 * D_new_u - aVertex->Attr().mValToOpt; 

                        // CONVERGENCE CHECK
                        aSqChange += (D_new_u - aVertex->Attr().mValToOpt) * (D_new_u - aVertex->Attr().mValToOpt);
                        aSqD      += D_new_u * D_new_u;
    
                        aVertex->Attr().mValToOpt = D_new_u; // update primal variable
                    }
                    Type aRelChange = std::sqrt(aSqChange / (aSqD + 1e-6));

                    StdOut() << " Iteration " << anIter << ", Squared Change = " << aSqChange << ", Relative Change = " << aRelChange << "\n";
                    if (aRelChange < aParamsAlgo.mTolerance)
                    {
                        StdOut() << "Convergence reached at iteration " << anIter << " with relative change " << aRelChange << std::endl;
                        break;
                    }

                }
        }

        /*void HighestEigenVal(tGraph& aGraph,
                            Type & aHighestEigenVal)
        {
            // Power iteration on L = B^T W^2 B (symmetric positive semi-definite).
            const size_t aNbV = aGraph.AllVertices().size();
            if (aNbV == 0)
            {
                aHighestEigenVal = 0.0;
                return;
            }

            auto ApplyL = [&aGraph,aNbV](const std::vector<Type> & aX,std::vector<Type> & aY)
            {
                aY.assign(aNbV,0.0);
                for (const auto & anEdge : aGraph.AllEdges())
                {
                    int j1 = (int)anEdge->VertexInit().Attr().mAbsCurv;
                    int j2 = (int)anEdge->Succ().Attr().mAbsCurv;
                    Type w_e = anEdge->AttrSym().mRanCost;
                    Type coeff = (w_e*w_e) * (aX.at(j1) - aX.at(j2));
                    aY.at(j1) += coeff;
                    aY.at(j2) -= coeff;
                }
            };

            std::vector<Type> aX(aNbV,1.0);
            std::vector<Type> aY;

            // Normalize the initial vector.
            Type aNormX = std::sqrt((Type)aNbV);
            for (auto & aVal : aX)
                aVal /= (aNormX + 1e-6);

            Type aPrevLambda = -1.0;
            Type aLambda = 0.0;

            for (int aK=0; aK<100; aK++)
            {
                ApplyL(aX,aY);

                Type aNormY2 = 0.0;
                for (const auto & aVal : aY)
                    aNormY2 += aVal*aVal;
                Type aNormY = std::sqrt(aNormY2);
                if (aNormY < 1e-14)
                {
                    aHighestEigenVal = 0.0;
                    return;
                }

                for (size_t aV=0; aV<aNbV; aV++)
                    aX[aV] = aY[aV] / aNormY;

                // Rayleigh quotient lambda = (x^T L x)/(x^T x).
                ApplyL(aX,aY);
                Type aNum = 0.0;
                Type aDen = 0.0;
                for (size_t aV=0; aV<aNbV; aV++)
                {
                    aNum += aX[aV]*aY[aV];
                    aDen += aX[aV]*aX[aV];
                }
                aLambda = aNum / (aDen + 1e-12);

                if ((aPrevLambda > 0.0) && (std::abs(aLambda-aPrevLambda) / (aPrevLambda+1e-12) < 1e-5))
                    break;

                aPrevLambda = aLambda;
            }

            aHighestEigenVal = std::max<Type>(aLambda,0.0);
        }*/
    private:
};


template <class Type>
class cGG_Graph : public cVG_Graph<cGG_AttrVert<Type>,cGG_AttrEdgeOr<Type>,cGG_AttrEdgeSym<Type>>
{
    public:
        /* Typedef part*/
        typedef cVG_Graph<cGG_AttrVert<Type>,cGG_AttrEdgeOr<Type>,cGG_AttrEdgeSym<Type>> tGraph;
        typedef typename tGraph::tVertex tVertex;
        typedef typename tGraph::tEdge tEdge;
        typedef typename std::vector<tVertex*> tVecVert;
        typedef tVertex* tPtrVert;

        // add other typdefs later

        cGG_Graph(cPt2di aSzGrid, bool is8Connex=true);


        tPtrVert & VertOfPt(const cPt2di & aPt)
        {
            return mGridVertices.at(aPt.y()).at(aPt.x());
        }

        tEdge&  EdgeOfPts(const cPt2di & aP1,const cPt2di & aP2)
        {
            tVertex & aV1 = *VertOfPt(aP1);
            tVertex & aV2 = *VertOfPt(aP2);
            return *aV1.EdgeOfSucc(aV2)->EdgeInitOr();
        }

        tPtrVert & VertOfAbsCurv(int aAbsCurv)
        {
            return FindIf(this->AllVertices(),[aAbsCurv](const auto & aV) {return aV->Attr().mAbsCurv==aAbsCurv;});
        }

        int NumQuad(const cPt2di & aPt,int aNb)
        {
            cPt2dr aSzR = ToR(mSzGrid) / tREAL8(aNb-0.5);
            return   round_ni(aPt.x()/aSzR.x())  + round_ni((aPt.y()/aSzR.y())) * aNb;
        }

        /// use to compute a unique id up to a circular permutation
        //size_t PathHashCode(const std::vector<tEdge *>&) const;

        void AddEdge(const cPt2di & aP1,const cPt2di & aP2);

         class cNeigh_4_Connex : public  cAlgo_ParamVG<cGG_Graph>
         {
            public :
              // this formula validate the edge iff  |x1-x2|+|y1-y2| <= 1
                   bool   InsideEdge(const    tEdge & anE) const override
                   {
                     return Norm1(anE.VertexInit().Attr().mPt-anE.Succ().Attr().mPt)<=1;
                   }
         };

        class cNeigh_8_Connex : public  cAlgo_ParamVG<cGG_Graph>
         {
            public :
              // this formula validate the edge iff  |x1-x2| <= 1 and |y1-y2| <= 1
                   bool   InsideEdge(const tEdge & anE) const override
                   {
                     cPt2di aD = (anE.VertexInit().Attr().mPt-anE.Succ().Attr().mPt);
                     return (abs(aD.x())<=1) && (abs(aD.y())<=1);
                   }
         };

        bool Is8Connex() const { return mIs8Connex; }

    private:
        cPt2di mSzGrid;
        cRect2  mBox;
        bool mIs8Connex; // flag to indicate if the graph is 8 connexity or 4 connexity, used for edge creation and validation
        std::vector<std::vector<tPtrVert>> mGridVertices; // 2D vector of pointers to vertex, for direct access to vertex from pixel coordinates

    };


template <class Type>
cGG_Graph<Type>::cGG_Graph(cPt2di aSzGrid, bool is8Connex) :
    cVG_Graph<cGG_AttrVert<Type>,cGG_AttrEdgeOr<Type>,cGG_AttrEdgeSym<Type>>(),
    mSzGrid(aSzGrid),
    mBox(mSzGrid),
    mIs8Connex(is8Connex),

    // Allocate the vector of vector
    mGridVertices(mSzGrid.y(),std::vector<tPtrVert>(mSzGrid.x(),nullptr))
{
    // 1 - create vertexes

    for (const auto & aP : mBox) // parse the box, for each pixel create a vertex
    {
         // Formula to create curvilinear absisca 
        cPt2di aPix4Num
                (
                (aP.y()%2) ? (mSzGrid.x()-aP.x()-1) : aP.x(),
                aP.y()
                );

        Type anAbsicsa = mBox.IndexeLinear(aPix4Num);
        VertOfPt(aP) = this->NewVertex(cGG_AttrVert<Type>(aP,anAbsicsa));
    }

    // 2 - create edges
    for (const auto & aP : mBox) // parse the box, for each pixel create edge with right and down neighbour
    {
        if (mIs8Connex)
        {
            AddEdge(aP,aP+cPt2di(1,0)); // right neighbour
            AddEdge(aP,aP+cPt2di(0,1)); // down neighbour
            AddEdge(aP,aP+cPt2di(1,1)); // right down neighbour
            AddEdge(aP,aP+cPt2di(-1,1)); // left down neighbour
        }
        else
        {
            AddEdge(aP,aP+cPt2di(1,0)); // right neighbour
            AddEdge(aP,aP+cPt2di(0,1)); // down neighbour
        }
    }
}


template <class Type>
void cGG_Graph<Type>::AddEdge(const cPt2di & aP1,const cPt2di & aP2)
{
    if (! mBox.Inside(aP1) || ! mBox.Inside(aP2)) // if one of the point is outside the box, do not create the edge
        return;

    tVertex * aV1 = VertOfPt(aP1);
    tVertex * aV2 = VertOfPt(aP2);

    Type aDist = Norm2(aV1->Attr().mPt-aV2->Attr().mPt);
    Type aDiffAbs12= aV2->Attr().mAbsCurv - aV1->Attr().mAbsCurv;

    cGG_AttrEdgeOr<Type> aAttrOr12(aDiffAbs12);
    cGG_AttrEdgeOr<Type> aAttrOr21(-aDiffAbs12);

    cGG_AttrEdgeSym<Type> aAttrSym(aDist);

    tGraph::AddEdge(*aV1,*aV2,aAttrOr12,aAttrOr21,aAttrSym);
}



// class cImageImpainting 
template <class Type> class cImageImpainting
{
    public: 

    cImageImpainting (cIm2D<tU_INT1> & aImMasqInit,
                       cIm2D<tU_INT1> & aImMasqFinal,
                       cIm2D<Type> & aImToFill,
                       tREAL4 aDyn):
        mImMasqInit(aImMasqInit),
        mImMasqFinal(aImMasqFinal),
        mImToFill(aImToFill),
        mBuffIm(cIm2D<Type>(aImToFill.DIm().Sz())),
        mSetPixToFill({}),
        mSz(aImToFill.DIm().Sz()),
        mDyn(aDyn),
        mNbPixToFill(0)
    {

        for(const auto & aP: cPixBox<2>(cPt2di(1,1), mSz-cPt2di(1,1)))
         {
            if (mImMasqInit.DIm().GetV(aP) &&  (!mImMasqFinal.DIm().GetV(aP)))
            {
                mSetPixToFill.push_back(aP);
            }
         }
         StdOut() << "Number of pixels to fill: " << mSetPixToFill.size() << "\n";
         mNbPixToFill = mSetPixToFill.size();
    }

    void FillMissingPixels()
    {
        std::vector<cPt2di> aNeighPts={cPt2di(-1,0),cPt2di(1,0),cPt2di(0,-1),cPt2di(0,1),
                                cPt2di(-1,-1),cPt2di(-1,1),cPt2di(1,-1),cPt2di(1,1)};

        for (const auto & aP: mImToFill.DIm())
        {
            Type aVal = mImMasqFinal.DIm().GetV(aP) ? mImToFill.DIm().GetV(aP)*mDyn : 1e10;
            mBuffIm.DIm().SetV(aP,aVal);
        }

        for (int aKIter=0 ; aKIter< 6 ; aKIter++)
        {

            bool Pair= ((aKIter%2)==0);

            int IndDeb = Pair ? 0       : (mNbPixToFill-1);
            int IndOut = Pair ? mNbPixToFill  : (-1)      ;
            int Incr   = Pair ? 1       : (-1)      ;

            for (int Ind=IndDeb ; Ind!=IndOut ; Ind+=Incr)
            {
                    cPt2di aP2Cur = mSetPixToFill[Ind];
                    Type aValOfZMin = (Type)mImToFill.DIm().GetV(aP2Cur);
                    Type aZMin = mBuffIm.DIm().GetV(aP2Cur);
                    for (int aKV = 0 ; aKV<8 ; aKV++)
                    {
                        cPt2di aPVois = aP2Cur + aNeighPts[aKV];
                        if (mImMasqFinal.DIm().GetV(aPVois) || mImMasqInit.DIm().GetV(aPVois))
                        {
                            Type aZAugm = mBuffIm.DIm().GetV(aPVois) + ((aKV%2) ? 3 : 2);
                            if (aZAugm < aZMin)
                            {
                                aZMin = aZAugm;
                                aValOfZMin = (Type)mImToFill.DIm().GetV(aPVois);
                            }
                        }
                    }
                    //StdOut() << "Filling pixel " << aP2Cur << " with value " << aValOfZMin << "\n";
                    mImToFill.DIm().SetV(aP2Cur,aValOfZMin);
                    mBuffIm.DIm().SetV(aP2Cur,aZMin);
            }
        }
    }

    cIm2D<tU_INT1> & mImMasqInit;
    cIm2D<tU_INT1> & mImMasqFinal;
    cIm2D<Type> & mImToFill;
    cIm2D<Type> mBuffIm;
    std::vector<cPt2di> mSetPixToFill;
    cPt2di mSz;
    tREAL4 mDyn;
    size_t mNbPixToFill;

    private:

};

class cAppliParsedBoxVirtualIm : public cAppliParseBoxIm<tREAL4>
{
    private :
        cBox2di mGlobBoxToParse;
        std::vector<cBox2di> mSetBoxes;

    public:
        cAppliParsedBoxVirtualIm(cMMVII_Appli & anAppli,
                                const cPt2di & aSzTiles,
                                const cPt2di & aSzOverlap,
                                bool ParalTiles);

        void APBI_ExecAll (bool Silence=false);

        cBox2di GlobBox2Parse(){return mGlobBoxToParse;}
        std::vector<cBox2di> SetOfBoxes(){return mSetBoxes;}

        void SetBox2Parse(cBox2di & aBox){mGlobBoxToParse=aBox;}
        void SetTheSetOfBoxes(std::vector<cBox2di> & aSet) {mSetBoxes=aSet;}
};


cAppliParsedBoxVirtualIm::cAppliParsedBoxVirtualIm(cMMVII_Appli & anAppli,
                                                    const cPt2di & aSzTiles,
                                                    const cPt2di & aSzOverlap,
                                                    bool ParalTiles):
     cAppliParseBoxIm<tREAL4>(anAppli,eForceGray::No,aSzTiles,aSzOverlap,ParalTiles),
     mGlobBoxToParse(cBox2di::Empty()),
     mSetBoxes({})
{    
}

void cAppliParsedBoxVirtualIm::APBI_ExecAll(bool Silence)
{
    if (APBI_TestMode())
     {
        LoadI(CurBoxIn());
        mAppli.ExeOnParsedBox();
        return;
     }
    AssertNotInParsing();

    // int cParseBoxInOut
    cParseBoxInOut<2> aPBIO=cParseBoxInOut<2>::CreateFromSize(cBox2di(cPt2di(2000,2000)),
                                        mSzTiles);

    MMVII_INTERNAL_ASSERT_strong(!mGlobBoxToParse.IsEmpty(),
                                "Box to parse should not be empty for tiling");

    // update mSzTiles
    std::vector<cPt2di> aSetOfPixIndexes;
    bool aCustomSetOfBoxes=false;
    if(!mSetBoxes.empty())
    {
        mSzTiles = mSetBoxes[0].Sz();
        aCustomSetOfBoxes=true;
    }

    aPBIO  = cParseBoxInOut<2>::CreateFromSize(mGlobBoxToParse,mSzTiles);

    if(!mSetBoxes.empty())
    {
        for (const auto & aBox: mSetBoxes)
        {

            //StdOut()<<"box glob to Index coord "<<aPBIO.BoxIndex().FromNormaliseCoord(aPBIO.BoxGlob().ToNormaliseCoord(aBox.P0()))<<"\n";
            cPt2di BoxPixIndex = CByC1P(aPBIO.BoxIndex().FromNormaliseCoord(aPBIO.BoxGlob().ToNormaliseCoord(aBox.P0())),
                        round_ni);
            aSetOfPixIndexes.push_back(BoxPixIndex);
            //StdOut() << "Box to parse: " << aBox << " with pixel index: " << BoxPixIndex << "\n";
        }
    }

    mParseBox = & aPBIO;
    std::list<cParamCallSys>  aLComParal;

    if(aCustomSetOfBoxes)
    {
        for (const auto & aPixI : aSetOfPixIndexes)
        {
            // if a the top level of paralelization, construct the string
            // For first box, run it classically so that files are created only once
            if (TopCallParallTile() /*&& (aPixI!=aSetOfPixIndexes[0])*/)
            {
                cParamCallSys aCom = mAppli.CommandOfMain();
                aCom.AddArgs(NameIndBoxRecal + "=" + ToStr(aPixI));
                aLComParal.push_back(aCom);
            }
            else
            {
                // If not in paral do all box, else do only the box indicate by recall

                if ((!mParalTiles) || (aPixI==mIndBoxRecal))
                {
                    mCurPixIndex = aPixI;
                    //StdOut() << "Processing box with pixel index: " << aPixI << "\n";
                    //LoadI(CurBoxIn());
                    mAppli.ExeOnParsedBox();
                }
            }
        }
    }
    else
    {
        for (const auto & aPixI : aPBIO.BoxIndex())
        {
            // if a the top level of paralelization, construct the string
            // For first box, run it classically so that files are created only once
            if (TopCallParallTile() && (aPixI!=cPt2di(0,0)))
            {
            //std::string aCom = mAppli.CommandOfMain() + " " +NameIndBoxRecal + "=" + ToStr(aPixI);
                cParamCallSys aCom = mAppli.CommandOfMain();
                aCom.AddArgs(NameIndBoxRecal + "=" + ToStr(aPixI));
                aLComParal.push_back(aCom);
            }
            else
            {
                // If not in paral do all box, else do only the box indicate by recall
                if ((!mParalTiles) || (aPixI==mIndBoxRecal))
                {
                    mCurPixIndex = aPixI;
                    //LoadI(CurBoxIn());
                    mAppli.ExeOnParsedBox();
                }
            }
        }
    }

    mParseBox = nullptr ;   // No longer inside parsing

    //StdOut() << "Number of boxes to process in parallel: " << aLComParal.size() << "\n";
    mAppli.ExeComParal(aLComParal,Silence);

}


class cAppliCloudFuser : public cMMVII_Appli,
                         public cAppliParsedBoxVirtualIm

{
    private:
    int Exe() override;
    int ExeOnParsedBox() override;
    cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override ;
    cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override ;

    std::string mSetDemsPattern;
    std::vector<std::string> mSetDemsNames;
    std::vector<cCloudRaster> mSetDems;
    cBox2di mGlobalPixBoxFusion;
    cTplBoxOfPts<tREAL8,2> mGlobalRealBoxFusion;
    cAffin2D<tREAL8> mGlobTf;
    cAffin2D<tREAL8> mCurBoxInAffTf;
    std::vector<cBox2dr> mSetRealBoxCalc;
    std::vector<cBox2di> mSetPixBoxCalc;
    tREAL8 mGlobGSD;
    std::string mTilingIndexFile;
    bool IsOutputTiled;
    std::string mNamePrefixOut;
    std::string mNameFusionDemResult;
    std::string mNameFusionCorrelResult;
    std::string mNameFusionMasqResult;


    // here add private objets 
    public:

        cAppliCloudFuser(const std::vector<std::string> & aVArgs, const cSpecMMVII_Appli & aSpec );
        static constexpr tREAL8 mInfty =  -1e10;
        void ReadIndexTilingAndMetadataFile();
        void GenTFW(const cAffin2D<tREAL8> & anAff, const std::string & aNameTFW);
        void FuseDemsByWAvg( std::vector<cIm2D<tREAL4>> & aVDems,
                                 std::vector<cIm2D<tREAL4>> & aVWeighters,
                                 cIm2D<tREAL4> & aMergedDem,
                                 cIm2D<tU_INT1> & aMergedMask,
                                 cIm2D<tU_INT1> & aMergedCorrel);
        void FuseDemsByMedian( std::vector<cIm2D<tREAL4>> & aVDems,
                                 std::vector<cIm2D<tREAL4>> & aVWeighters,
                                 cIm2D<tREAL4> & aMergedDem,
                                 cIm2D<tU_INT1> & aMergedMask,
                                 cIm2D<tU_INT1> & aMergedCorrel);
        void FuseDemsByChP( std::vector<cIm2D<tREAL4>> & aVDems,
                                 std::vector<cIm2D<tREAL4>> & aVWeighters,
                                 cIm2D<tREAL4> & aMergedDem,
                                 cIm2D<tU_INT1> & aMergedMask,
                                 cIm2D<tU_INT1> & aMergedCorrel);
        void FillVoids(cIm2D<tREAL4> & aDem, cIm2D<tU_INT1> & aMask);
        void ExtractInterpDems(cCloudRaster & aDem,
                                const cBox2dr & aBoxInDem, 
                                const cBox2dr & aBoxInCurTileCalc,
                                const cAffin2D<tREAL8> & anAffCurBox,
                                cIm2D<tREAL4> & aDemToExtract,
                                cIm2D<tREAL4> & aWeightingIm);
};

cAppliCloudFuser::cAppliCloudFuser(const std::vector<std::string> & aVArgs, const cSpecMMVII_Appli & aSpec):
    cMMVII_Appli(aVArgs,aSpec), 
    cAppliParsedBoxVirtualIm(*this,cPt2di(2000,2000),cPt2di(50,50),true),
    mSetDemsPattern(""),
    mSetDemsNames({}),
    mSetDems({}),
    mGlobalPixBoxFusion(cBox2di::Empty()),
    mGlobalRealBoxFusion(cTplBoxOfPts<tREAL8,2>()),
    mGlobTf(cAffin2D<tREAL8>()),
    mCurBoxInAffTf(cAffin2D<tREAL8>()),
    mSetRealBoxCalc({}),
    mSetPixBoxCalc({}),
    mGlobGSD(0.0),
    mTilingIndexFile(""),
    IsOutputTiled(false),
    mNamePrefixOut("Fusion"),
    mNameFusionDemResult(""),
    mNameFusionCorrelResult(""),
    mNameFusionMasqResult("")
{
}



cCollecSpecArg2007 & cAppliCloudFuser::ArgObl(cCollecSpecArg2007 & anArgObl)
{
    return
            //APBI_ArgObl(anArgObl)
            anArgObl
           <<   Arg2007(mSetDemsPattern,
                        "Set of xml files each defines a dem, an ambiuity/correl map and a Masq of definition",
                        {{eTA2007::MPatFile,"0"}})
        ;
}

cCollecSpecArg2007 & cAppliCloudFuser::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
    return 
            APBI_ArgOpt
        (
            anArgOpt
            <<AOpt2007(mTilingIndexFile,
                        "TilingIndex",
                        "Tiling Index file, i.e a GeoJson file that defines the tiling scheme.\n"
                        "If not given, a big tif is created as the fusion result.", 
                        {eTA2007::FileTagged})
            << AOpt2007(mNamePrefixOut,"Out"," prefix of output files, default=Fusion.tif",{eTA2007::HDV})
        )
        ;
}


void cAppliCloudFuser::ReadIndexTilingAndMetadataFile()
{
    GDALAllRegister();  // register all drivers (call once at startup)

    GDALDataset* poDS = (GDALDataset*) GDALOpenEx(
        mTilingIndexFile.c_str(),
        GDAL_OF_VECTOR | GDAL_OF_READONLY,
        nullptr,   // allowed drivers (nullptr = all)
        nullptr,   // open options
        nullptr    // sibling files
    );

    if (poDS == nullptr) {
        MMVII_INTERNAL_ASSERT_strong(false,
                    "Unable to read The tiling Index "+mTilingIndexFile);
    }

    // read features and capture geometry, i.e each tile index 

    // assume one layer

    OGRLayer* poLayer = poDS->GetLayer(0);
    poLayer->ResetReading();  

    OGRFeature* poFeature;

    while((poFeature=poLayer->GetNextFeature())!=nullptr)
    {
        OGRGeometry* poGeom = poFeature->GetGeometryRef();

        if(poGeom!=nullptr)
        {
            OGRwkbGeometryType eType = wkbFlatten(poGeom->getGeometryType());
            MMVII_INTERNAL_ASSERT_User((eType==wkbPolygon) || (eType==wkbMultiPolygon),
                                        eTyUEr::eBadEnum,
                                    "Bad Geometry type in GeoJson/shp file, It Should be Polygon or MultiPolygon");

            // now read 
           /*OGRPolygon* poPoly = poGeom->toPolygon();
            OGRLinearRing* poRing = poPoly->getExteriorRing();*/

            OGREnvelope aTileEnvelope;
            poGeom->getEnvelope(&aTileEnvelope);
            // get tile extent 
            cBox2dr aBoxTile(cPt2dr(aTileEnvelope.MinX,aTileEnvelope.MinY),
                            cPt2dr(aTileEnvelope.MaxX,aTileEnvelope.MaxY));

            // get Global real Box for fusion
            mGlobalRealBoxFusion.Add(aBoxTile.P0());
            mGlobalRealBoxFusion.Add(aBoxTile.P1());

            mSetRealBoxCalc.push_back(aBoxTile);
        }
        OGRFeature::DestroyFeature(poFeature);
    }

    GDALClose(poDS);

    //sort the set of pix boxes from low to high box coordinates
    //std::sort(mSetPixBoxCalc.begin(),mSetPixBoxCalc.end());
}

void cAppliCloudFuser::GenTFW(const cAffin2D<tREAL8> & anAff, const std::string & aNameTFW)
{
    std::ofstream aFtfw(aNameTFW.c_str());
    aFtfw.precision(10);

    aFtfw << anAff.VX().x() << "\n" << anAff.VX().y() << "\n";
    aFtfw << anAff.VY().x() << "\n" << anAff.VY().y() << "\n";
    aFtfw << anAff.Tr().x() << "\n" << anAff.Tr().y() << "\n";

    aFtfw.close();
}


void cAppliCloudFuser::FuseDemsByMedian( std::vector<cIm2D<tREAL4>> & aVDems,
                                 std::vector<cIm2D<tREAL4>> & aVWeighters,
                                 cIm2D<tREAL4> & aMergedDem,
                                 cIm2D<tU_INT1> & aMergedMask,
                                 cIm2D<tU_INT1> & aMergedCorrel)

{
    std::vector<std::pair<tREAL8,tREAL8>> aVZwithWeights;
    tREAL8 aTotWeights;
    int aNb;
    tU_INT1 aCorrel;
    for (const auto & aPix: aMergedDem.DIm())
    {
        aVZwithWeights.clear() ;
        aTotWeights= 0.0 ;
        aNb=0;
        aCorrel=0.0;

        for (size_t i=0; i<aVDems.size();i++)
        {
            tREAL8 aWeight = aVWeighters[i].DIm().GetV(aPix);
            if( aWeight > 0) // masq has a 1 value
            {
                aTotWeights += aWeight ;
                aVZwithWeights.push_back({aVDems[i].DIm().GetV(aPix),aWeight});
                aNb++;
            }
        }
        aMergedDem.DIm().SetV(aPix,
            (aTotWeights==0) ? 0.0 : NonConstWMediane(aVZwithWeights,aTotWeights));
        aCorrel = ( aNb==0 ) ? 0: std::min(round_ni(255.0*(aTotWeights/aNb)),255); 

        aMergedCorrel.DIm().SetV(aPix,aCorrel);
        aMergedMask.DIm().SetV(aPix,aNb>0 ? 1 : 0);
    }
}


void cAppliCloudFuser::FuseDemsByChP( std::vector<cIm2D<tREAL4>> & aVDems,
                                 std::vector<cIm2D<tREAL4>> & aVWeighters,
                                 cIm2D<tREAL4> & aMergedDem,
                                 cIm2D<tU_INT1> & aMergedMask,
                                 cIm2D<tU_INT1> & aMergedCorrel)

{
    // initialize the graph with 8 connexity
    cGG_Graph<tREAL4> aGraph(aMergedDem.DIm().Sz(),false);
    cAlgo_ChambollePockParams<tREAL4> aParamsAlgo(0.05,2.0,60,1e-5,0.99);
    cAlgo_ChambollePock<cGG_Graph<tREAL4>,tREAL4> mOptimizer; 

    tREAL4 aHighestEigenVal=2*std::sqrt(2.0); // for 8 connexity, the highest eigenvalue of the incidence matrix is 2*sqrt(2)

    if (aGraph.Is8Connex())
    {
        aHighestEigenVal=4.0; // for 4 connexity, the highest eigenvalue of the incidence matrix is 4
    }
    // if connexity is 4, the highest eigenvalue of the incidence matrix is 2 sqrt(2)
    // if connexity is 8, the highest eigenvalue of the incidence matrix is 4
    mOptimizer.SolveChambollePock(aGraph,
                            aHighestEigenVal,
                            aParamsAlgo,
                            aVDems,     
                            aVWeighters);

    // Save the optimized variable back to aMergedDem
    for(const auto & aVert : aGraph.AllVertices())
    {
        aMergedDem.DIm().SetV(aVert->Attr().mPt,aVert->Attr().mValToOpt);
        aMergedMask.DIm().SetV(aVert->Attr().mPt,(aVert->Attr().mWeightPerNode > 0.0) ? 1 : 0);
            // set correl value based on the weight, for example, you can set it to
        MMVII_INTERNAL_ASSERT_tiny(aVert->Attr().mWeightPerNode >= 0.0, "Weight per node should be non-negative");
        aMergedCorrel.DIm().SetV(aVert->Attr().mPt, 
                                    std::min(
                                    round_ni(255.0*aVert->Attr().mWeightPerNode/aVDems.size())
                                    ,255)
                                    );    
    }
}


void cAppliCloudFuser::FillVoids(cIm2D<tREAL4> & aDem, cIm2D<tU_INT1> & aMask)
{
    cIm2D<tU_INT1> aMasqInit(aMask.DIm().Sz());
    aMasqInit.DIm().InitCste(1);
    cImageImpainting<tREAL4> aImpaint(aMasqInit,aMask,aDem,1.0);
    aImpaint.FillMissingPixels();
}


void cAppliCloudFuser::FuseDemsByWAvg( std::vector<cIm2D<tREAL4>> & aVDems,
                                 std::vector<cIm2D<tREAL4>> & aVWeighters,
                                 cIm2D<tREAL4> & aMergedDem,
                                 cIm2D<tU_INT1> & aMergedMask,
                                 cIm2D<tU_INT1> & aMergedCorrel)
{
    tREAL8 aVWeightedZ=0.0;
    tREAL8 aWeights =0.0;
    tU_INT1 aCorrel =0 ;
    int aNb=0;
    // assume all image sizes are equal so that at one pixel location we have multiple values of z correl and Masq
    for (const auto & aPix: aMergedDem.DIm())
    {
        aVWeightedZ=0.0;
        aWeights =0.0;
        aCorrel =0 ;
        aNb=0;  

        for (size_t i=0; i<aVDems.size();i++)
        {
            tREAL8 aWeight = aVWeighters[i].DIm().GetV(aPix);
            if( aWeight >= 0) // masq has a 1 value
            {
                aWeights += aWeight ;
                aVWeightedZ+=aWeight*aVDems[i].DIm().GetV(aPix);
                aNb++;
            }
        }
        // Weight dem values
        aMergedDem.DIm().SetV(aPix,(aWeights==0) ? 0.0 : aVWeightedZ/aWeights);
        aCorrel = ( aNb==0 ) ? 0: std::min(round_ni(255.0*(aWeights/aNb)),255); 

        aMergedCorrel.DIm().SetV(aPix,aCorrel);
        aMergedMask.DIm().SetV(aPix,aNb>0 ? 1 : 0);
    }
}


void cAppliCloudFuser::ExtractInterpDems(cCloudRaster & aDem,
                                        const cBox2dr & aBoxInDem, 
                                        const cBox2dr & aBoxInCurTileCalc,
                                        const cAffin2D<tREAL8> & anAffCurBox,
                                        cIm2D<tREAL4> & aDemToExtract,
                                        cIm2D<tREAL4> & aWeightingIm)
{  
    //initialize values 
    aDemToExtract.DIm().InitCste(mInfty);
    aWeightingIm.DIm().InitCste(-1.0);

    cAffin2D<tREAL8> aDemAff= aDem.Transform();
    
    cDataFileIm2D aDemDFile =  cDataFileIm2D::Create(aDem.NameDem(),eForceGray::No) ;
    cDataFileIm2D aMasqDFile =  cDataFileIm2D::Create(aDem.NameMasq(),eForceGray::No) ;
    cDataFileIm2D aCorrelDFile =  cDataFileIm2D::Create(aDem.NameCorrel(),eForceGray::No) ;  
    
    cBox2di aBoxIntInDem = aBoxInDem.ToI();
    cBox2di aBoxIntInCurTileCalc = aBoxInCurTileCalc.ToI();

    //check boundaries 
    aBoxIntInDem = aBoxIntInDem.Inter(cBox2di(aDem.SzPix()));
    aBoxIntInCurTileCalc = aBoxIntInCurTileCalc.Inter(CurBoxInLoc());


    // make safe boundaries of aBoxIntInCurTileCalc 
    cBox2di aBoxDemTransform = anAffCurBox.InverseBox(
                                aDemAff.ValueBox(aBoxIntInDem.ToR())
                            ).ToI();

    aBoxIntInCurTileCalc = aBoxIntInCurTileCalc.Inter(aBoxDemTransform);

    // the way back 
    /*cBox2di aBoxCurTransform = aBoxedDemAffT.InverseBox(
                                anAffCurBox.ValueBox(aBoxIntInCurTileCalc.ToR())
                            ).ToI();

    aBoxIntInDem = aBoxIntInDem.Inter(aBoxCurTransform);*/
    // compute box in dem affine transform

    cAffin2D<tREAL8> aBoxedDemAffT = cAffin2D<tREAL8>(aDemAff.Value(ToR(aBoxIntInDem.P0())),
                                    aDemAff.VX(),
                                    aDemAff.VY());

    StdOut()<<"aBoxInInDem "<<aBoxIntInDem<<" BOX IM "<<
    " aBoxIntInCurTileCalc "<<aBoxIntInCurTileCalc<<"  "<<
    " aBoxIntInDemAfterTransform "<<aBoxDemTransform<<" "<<
    aDemToExtract.DIm().Sz()<<std::endl;


    // read tiled data and not the entire dem 
    cIm2D<tREAL4> aBoxedDem(aBoxIntInDem.Sz());
    aBoxedDem.Read(aDemDFile,aBoxIntInDem.P0(),1.0);

    cIm2D<tU_INT1> aBoxedCorrel(aBoxIntInDem.Sz());
    aBoxedCorrel.Read(aCorrelDFile,aBoxIntInDem.P0(),1.0);

    cIm2D<tU_INT1> aBoxedMasq(aBoxIntInDem.Sz());
    aBoxedMasq.Read(aMasqDFile,aBoxIntInDem.P0(),1.0);

    tREAL4 aZ,aWeighter;
    for (const auto & aPix: cPixBox<2>(aBoxIntInCurTileCalc.P0(),
                                        aBoxIntInCurTileCalc.P1()))
    {
        aWeighter = -1.0;
        aZ = mInfty;
        // perform Transforms to get loc in boxed dems
        cPt2dr aPWorld  = anAffCurBox.Value(ToR(aPix));
        cPt2dr aPInBoxDem  = aBoxedDemAffT.Inverse(aPWorld);

        cPt2di aPInBoxDemI = ToI(aPInBoxDem);
        cPt2di aPInBoxDemUl = Pt_round_down(aPInBoxDem);

        if(aBoxedMasq.DIm().InsideBL(aPInBoxDem))
        {
            if (aBoxedDem.DIm().InsideBL(aPInBoxDem))
            {
                // avoid weighting and computing z with nodata in the neighboouthood
                auto & aBoxMasq = aBoxedMasq.DIm();
                if(
                    aBoxMasq.Inside(aPInBoxDemUl+cPt2di(1,0)) && 
                    aBoxMasq.Inside(aPInBoxDemUl+cPt2di(0,1))
                    )
                {
                    bool isAllNeighborDefined = aBoxMasq.GetV(aPInBoxDemUl) &&
                                                aBoxMasq.GetV(aPInBoxDemUl+cPt2di(1,0)) &&
                                                aBoxMasq.GetV(aPInBoxDemUl+cPt2di(0,1)) &&
                                                aBoxMasq.GetV(aPInBoxDemUl+cPt2di(1,1));
                    if (isAllNeighborDefined)
                    {

                        //bilinear interpolation
                        aWeighter = aBoxedCorrel.DIm().GetVBL(aPInBoxDem)/255.0;
                        aZ = aBoxedDem.DIm().GetVBL(aPInBoxDem);
                        /*MMVII_INTERNAL_ASSERT_tiny(aZ>=0, "postive dem value should be expected"+ToStr(aZ)+" Weighter "+ToStr(aWeighter)
                                + ToStr(aBoxedDem.DIm().GetV(aPInBoxDemUl))+" "+ToStr(aBoxedDem.DIm().GetV(aPInBoxDemUl+cPt2di(1,0)))+" "+
                                ToStr(aBoxedDem.DIm().GetV(aPInBoxDemUl+cPt2di(0,1)))+" "+ToStr(aBoxedDem.DIm().GetV(aPInBoxDemUl+cPt2di(1,1)))
                                );*/
                    }
                }
                else
                {
                    // nearest neighbor
                    aWeighter=aBoxedCorrel.DIm().GetV(aPInBoxDemI)/255.0;
                    aZ=aBoxedDem.DIm().GetV(aPInBoxDemI);
                }
                // set 
                aDemToExtract.DIm().SetV(aPix,aZ);
                aWeightingIm.DIm().SetV(aPix,aWeighter);
            }

        }
    }
}


int cAppliCloudFuser::ExeOnParsedBox()
{
    //Init Tiles affine transforms 

    // first get Tile transform

    // output region to write in cas we want the ouput to be tiled 

    cPt2dr anOffx= cPt2dr(mGlobGSD,0.0);
    cPt2dr anOffy=cPt2dr(0.0,-mGlobGSD);
    cPt2dr anOrigin = mGlobTf.Value(ToR(CurBoxOut().P0()));
    cAffin2D<tREAL8> aLocTileAffTransform = cAffin2D<tREAL8>(anOrigin,anOffx,anOffy);


    // the curbox affine transform
    mCurBoxInAffTf = cAffin2D<tREAL8>(mGlobTf.Value(ToR(CurBoxIn().P0())),
                                    anOffx,
                                    anOffy);
    

    // fast search of dems that overlap the concerned CurBoxIn()
    std::vector<cCloudRaster> mSetDemsOverlapCurBox;
    std::vector<cBox2dr> mSetDemsOverlapBoxes;
    std::vector<cBox2dr> mSetDemWhere2WriteBoxes;
    
    std::vector<cIm2D<tREAL4>> mSetOfBoxedDems;
    std::vector<cIm2D<tREAL4>> mSetOfBoxedWeighters;

    cBox2dr aBoxCurGlob = mGlobTf.ValueBox(CurBoxIn().ToR());
                            
    for (auto & aDem : mSetDems)
    {
        cAffin2D<tREAL8> aLocAff = aDem.Transform();

        cBox2dr aBoxDem = aLocAff.ValueBox(cBox2dr(ToR(aDem.SzPix())));

        //StdOut()<<aBoxDem<<std::endl;
        cBox2dr aBoxInter = aBoxCurGlob.Inter(aBoxDem); 
        if (! aBoxInter.IsEmpty())
        {
            mSetDemsOverlapCurBox.push_back(aDem);

            cBox2dr aBoxLocInDemR = aLocAff.InverseBox(aBoxInter);

            cBox2dr aBoxInCurBoxInR = mGlobTf.InverseBox(aBoxInter);
                                
            // aBox In CurBoxIn() coordinate system
            cBox2dr aBoxLocInCurBoxInR(aBoxInCurBoxInR.P0()-ToR(CurBoxIn().P0()),
                                        aBoxInCurBoxInR.P1()-ToR(CurBoxIn().P0()));

            mSetDemsOverlapBoxes.push_back(aBoxLocInDemR);
            mSetDemWhere2WriteBoxes.push_back(aBoxLocInCurBoxInR);
        }
    }

    // check overlapping dems
    if (!mSetDemsOverlapCurBox.size())
        return EXIT_SUCCESS; 
    
    auto anItDem = mSetDemsOverlapCurBox.begin();
    auto anItBoxIndem = mSetDemsOverlapBoxes.begin();
    auto anItBoxInCurBox = mSetDemWhere2WriteBoxes.begin();

    for (; anItDem != mSetDemsOverlapCurBox.end() && 
           anItBoxIndem != mSetDemsOverlapBoxes.end() && 
           anItBoxInCurBox != mSetDemWhere2WriteBoxes.end();
            ++anItDem, ++anItBoxIndem,++anItBoxInCurBox)
    {
        auto [aDem, aBoxInDem,aBoxInCurB] = std::tie(*anItDem, *anItBoxIndem,*anItBoxInCurBox);
        cIm2D<tREAL4> aDemMap(CurBoxIn().Sz());
        cIm2D<tREAL4> aWeighterMap(CurBoxIn().Sz());

        // read dems and construct weighters 
        ExtractInterpDems(aDem,aBoxInDem,aBoxInCurB,mCurBoxInAffTf,
                            aDemMap,aWeighterMap);

        mSetOfBoxedDems.push_back(aDemMap);
        mSetOfBoxedWeighters.push_back(aWeighterMap);
    }

    cIm2D<tREAL4> aFinalDem = cIm2D<tREAL4>(CurBoxIn().Sz());
    cIm2D<tU_INT1> aFinalMask = cIm2D<tU_INT1>(CurBoxIn().Sz());
    cIm2D<tU_INT1> aFinalCorrel = cIm2D<tU_INT1>(CurBoxIn().Sz());

    cAutoTimerSegm aTSChambolle(TimeSegm(),"TotalVariationFusion"); 

    FuseDemsByChP(mSetOfBoxedDems,mSetOfBoxedWeighters,
              aFinalDem,aFinalMask, aFinalCorrel);

    FillVoids(aFinalDem,aFinalMask);

    // Save DEM, MASK and CORREL
    cAutoTimerSegm aTSFinalWrite(TimeSegm(),"FinalWriting"); 
    if(IsOutputTiled)
    {
        std::string aNameMergedDem =DirOfPath(mNamePrefixOut,false)+"DEM_"+ 
                                    FileOfPath(mNamePrefixOut,false)+
                                    ToStr(mIndBoxRecal.x())+"_"+ToStr(mIndBoxRecal.y());

        
        GenTFW(aLocTileAffTransform,aNameMergedDem +".tfw");                            
        
        cDataFileIm2D  aDFDem = cDataFileIm2D::Create(aNameMergedDem+".tif",
                                                eTyNums::eTN_REAL4,
                                                CurBoxOut().Sz(),
                                                1);
                                
        aFinalDem.Write(aDFDem,-CurBoxOutLoc().P0(),1.0,CurBoxOutLoc());

        //Masq
        std::string aNameMergedMask =DirOfPath(mNamePrefixOut,false)+"MASQ_"+ 
                                    FileOfPath(mNamePrefixOut,false)+
                                    ToStr(mIndBoxRecal.x())+"_"+ToStr(mIndBoxRecal.y());

        GenTFW(aLocTileAffTransform,aNameMergedMask + ".tfw");     


        cDataFileIm2D  aDFMasq = cDataFileIm2D::Create(aNameMergedMask+".tif",
                                                eTyNums::eTN_U_INT1,
                                                CurBoxOut().Sz(),
                                                1);

        aFinalMask.Write(aDFMasq,-CurBoxOutLoc().P0(),1.0,CurBoxOutLoc());


        //Correl
        std::string aNameMergedCorrel =DirOfPath(mNamePrefixOut,false)+"CORREL_"+ 
                                    FileOfPath(mNamePrefixOut,false)+
                                    ToStr(mIndBoxRecal.x())+"_"+ToStr(mIndBoxRecal.y());

        GenTFW(aLocTileAffTransform,aNameMergedCorrel+".tfw");    

        cDataFileIm2D  aDFCorrel = cDataFileIm2D::Create(aNameMergedCorrel+".tif",
                                                eTyNums::eTN_U_INT1,
                                                CurBoxOut().Sz(),
                                                1);

        aFinalCorrel.Write(aDFCorrel,-CurBoxOutLoc().P0(),1.0,CurBoxOutLoc());
    }
    else
    {

        // create if needed 
        cDataFileIm2D  aDFH = cDataFileIm2D::Create(mNameFusionDemResult+".tif",
                                              eTyNums::eTN_REAL4,
                                              GlobBox2Parse().Sz(),
                                              1);
        cDataFileIm2D  aDFC = cDataFileIm2D::Create(mNameFusionCorrelResult+".tif",
                                              eTyNums::eTN_U_INT1,
                                              GlobBox2Parse().Sz(),
                                              1);
        cDataFileIm2D  aDFM = cDataFileIm2D::Create(mNameFusionMasqResult+".tif",
                                              eTyNums::eTN_U_INT1,
                                              GlobBox2Parse().Sz(),
                                              1); 

        aFinalDem.Write(aDFH,CurP0(),1.0,CurBoxOutLoc());
        aFinalCorrel.Write(aDFC,CurP0(),1.0,CurBoxOutLoc());
        aFinalMask.Write(aDFM,CurP0(),1.0,CurBoxOutLoc());
    }

    return EXIT_SUCCESS;
}


int cAppliCloudFuser::Exe()
{
    // may be mPhProj.FinishInit();
   mSetDemsNames = VectMainSet(0);

   
   int aNB=0;
   // read serialized cloudRaster objects from bascule to get global context
    for(const auto & aCldName: mSetDemsNames)
    {
        cCloudRaster aCldRst;
        ReadFromFile_Std(aCldRst,aCldName);
        mSetDems.push_back(aCldRst);

        // increment global context
        cAffin2D<tREAL8> aTF = aCldRst.Transform();


        // compute global context
        if(!IsInit(&mTilingIndexFile))
        {
            mGlobalRealBoxFusion.Add(aTF.Value(cPt2dr(aCldRst.SzPix().x(),0)));
            mGlobalRealBoxFusion.Add(aTF.Value(cPt2dr(0,aCldRst.SzPix().y())));
        }
        // gsd
        mGlobGSD+=aTF.VX().x();
        aNB++;
    }
    // gsd of output fused dem
    mGlobGSD/=aNB;

    //  Check if Tiling Index was initialized or not 

    if(IsInit(&mTilingIndexFile))
    {
        IsOutputTiled =true;
        // Computation parallelization will be directly be computed from the tiling index itself
        // i.e Box per computation

        bool CheckExtensionReadable = (ends_with(mTilingIndexFile,"geojson") ||
                                        ends_with(mTilingIndexFile,"shp") );

        MMVII_INTERNAL_ASSERT_User(CheckExtensionReadable,
                                eTyUEr::eBadPostfix,
                                "Index File is not a common geojson or shp file");

        // read mTilingIndexFile
        ReadIndexTilingAndMetadataFile();

        //StdOut()<<"mSetRealBoxCalc.size() "<<mSetRealBoxCalc.size()<<std::endl;

        // global affine transform for the resulting output image
        mGlobTf = cAffin2D<tREAL8>(cPt2dr(mGlobalRealBoxFusion.CurBox().P0().x(),
                                    mGlobalRealBoxFusion.CurBox().P1().y()),
                                    cPt2dr(mGlobGSD,0),
                                    cPt2dr(0,-mGlobGSD));

        
        for (const auto & aRBox: mSetRealBoxCalc)
        {
            mSetPixBoxCalc.push_back(mGlobTf.InverseBox(aRBox).ToI());
        }
    

        // some sorting and indexing with a mask to specify what should be computed 
        std::sort(mSetPixBoxCalc.begin(), mSetPixBoxCalc.end(), [](const cBox2di& a, const cBox2di& b) 
        {
            if (a.P0().y() != b.P0().y()) return a.P0().y() < b.P0().y();
            return a.P0().x() < b.P0().x();                          
        });

        // run on a set of boxes in parallel
        SetTheSetOfBoxes(mSetPixBoxCalc);

        mGlobalPixBoxFusion = ToI(mGlobalRealBoxFusion.CurBox().Sz()/mGlobGSD);
        SetBox2Parse(mGlobalPixBoxFusion);
        
    }
    else
    {
        // global affine transform for the resulting output image
        mGlobTf = cAffin2D<tREAL8>(cPt2dr(mGlobalRealBoxFusion.CurBox().P0().x(),
                                        mGlobalRealBoxFusion.CurBox().P1().y()),
                                        cPt2dr(mGlobGSD,0),
                                        cPt2dr(0,-mGlobGSD));
        // we will follow the standard MMVII tiling scheme in cAppliParseBoxImVirtual
        mGlobalPixBoxFusion = cBox2di(Pt_round_up(mGlobalRealBoxFusion.CurBox().Sz()/mGlobGSD));

        SetBox2Parse(mGlobalPixBoxFusion);
        
        //init output file names
        mNameFusionDemResult = DirOfPath(mNamePrefixOut,false)+"DEM_"+ 
                                    FileOfPath(mNamePrefixOut,false);
        mNameFusionCorrelResult =DirOfPath(mNamePrefixOut,false)+"CORREL_"+ 
                                    FileOfPath(mNamePrefixOut,false);
        mNameFusionMasqResult =DirOfPath(mNamePrefixOut,false)+"MASQ_"+ 
                                    FileOfPath(mNamePrefixOut,false);

        // TFW FILES
        GenTFW(mGlobTf,mNameFusionDemResult+".tfw");
        GenTFW(mGlobTf,mNameFusionCorrelResult+".tfw");   
        GenTFW(mGlobTf,mNameFusionMasqResult+".tfw");   
        
    }

    APBI_ExecAll();

    return EXIT_SUCCESS;
}


/*  ============================================= */
/*       ALLOCATION                               */
/*  ============================================= */

tMMVII_UnikPApli Alloc_CloudFuser(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli & aSpec)
{
    return tMMVII_UnikPApli(new cAppliCloudFuser(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpecCloudFuser
    (
        "CloudMMVII_Fuse",
        Alloc_CloudFuser,
        "Fuse multiple digital elevation maps given ambiguity/correlation maps",
        {eApF::Cloud},
        {eApDT::Image},
        {eApDT::Image},
        __FILE__
        );

};