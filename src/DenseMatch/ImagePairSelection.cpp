#include "cMMVII_Appli.h"
#include "MMVII_PCSens.h"
#include "MMVII_Geom2D.h"
#include "MMVII_Geom3D.h"
#include "MMVII_AllClassDeclare.h"
#include "MMVII_DeclareCste.h"
#include "MMVII_Interpolators.h"
#include "MMVII_Tpl_GraphAlgo_SPCC.h"
#include "../src/PoseEstim/VisPoseAndStructure.h"

#include <bit>

namespace MMVII
{
    #if __cplusplus >= 202002L
        #include <bit>
        #define POPCOUNT64(x) std::popcount(x)
    #else
        #define POPCOUNT64(x) __builtin_popcountll(x)
    #endif



    // Occupancy Level with general image aspect ratio
    class OccupancyLevel 
    {
        public:
            int w = 0, h = 0;
            std::vector<tU_INT8> bits;
        
            void Init(int width, int height) {
                w = std::max(1, width);
                h = std::max(1, height);
                // ceil(w*h / 64) words: enough 64-bit words to hold one bit per cell,
                // rounding up since a partial word still needs a full word allocated.
                size_t num_words = (static_cast<size_t>(w) * h + 63) / 64;
                bits.assign(num_words, 0ull);
            }
        
            // u, v must be normalized image coordinates in [0,1)
            void InsetTieP(cPt2dr aPt) {
                // clamp on BOTH sides : tie points falling outside the image frame
                // give negative (or >=1) normalized coords, which would otherwise
                // produce an out-of-range index and an out-of-bounds write.
                int cx = std::min(std::max(static_cast<int>(aPt.x() * w), 0), w - 1);
                int cy = std::min(std::max(static_cast<int>(aPt.y() * h), 0), h - 1);
                std::size_t idx = static_cast<size_t>(cy) * w + cx;
                bits[idx / 64] |= (tU_INT8(1) << (idx % 64));
            }
        
            int OccupiedCount() const {
                int s = 0;
                for (tU_INT8 word : bits) s += POPCOUNT64(word);
                return s;
            }
        
            int TotalCells() const { return w * h; }
    };


    // Occupancy Pyramid 
    class OccupancyPyramid 
    {
        public:
            OccupancyPyramid(cPt2di aSz):
            mSzIm(aSz)
            {
            }
            std::vector<OccupancyLevel> levels;
            cPt2di mSzIm;
        
            // image_width / image_height: pixel dimensions of the source image
            // L: number of pyramid levels (levels 0 .. L-1)
            // base_cells: total cell budget at level 0; doubles each subsequent level
            void build(int L, int base_cells = 512) {
                levels.clear();
                levels.resize(L);
        
                tREAL8 aAspectR = static_cast<tREAL8>(mSzIm.x()) / static_cast<tREAL8>(mSzIm.y());
        
                for (int l = 0; l < L; ++l) {
                    int cells = base_cells << l;   // total cell target for this level, doubling per level
                    int gw = std::max(1, static_cast<int>(std::round(std::sqrt(cells * aAspectR))));
                    int gh = std::max(1, static_cast<int>(std::round(std::sqrt(cells / aAspectR))));
                    levels[l].Init(gw, gh);
                }
            }
        
            // insert a tie point given in pixel coordinates (px, py) within the
            // original image_width x image_height used in build()
            void InsertInPyrStruct(cPt2dr aPt) {
                aPt= aPt/ToR(mSzIm);
                for (auto & lvl : levels) {
                    lvl.InsetTieP(aPt);
                }
            }
        
            // weighted coverage score: finer levels contribute more, same spirit
            // as the SfM next-best-view scoring function
            tREAL8 HierarchicalScore() const {
                tREAL8 s = 0.0;
                for (std::size_t l = 0; l < levels.size(); ++l) {
                    s += static_cast<tREAL8>(tU_INT8(1) << l) * levels[l].OccupiedCount();
                }
                return s;
            }
    };

    // Graph definition 

    /// vertex : image name, perspective center
    class cDIPS_AV
    {
        public:
           cDIPS_AV(cSensorCamPC * aSens) :
               mSens(aSens)   
           {}
           cSensorCamPC * mSens = nullptr;
    };

    /// Class for oriented attribute : nothing
    class cDIPS_Or
    {
    public :

    };

   /// Class for  symetric attribute : nothing
    class cDIPS_Sym
    {
    public :
        cSetHomogCpleIm mCpleH;
        cPt2dr          mWeightHomol;
        tREAL8          mWeightGeom;
    };
 
    typedef cVG_Graph<cDIPS_AV,cDIPS_Or,cDIPS_Sym>  tGrDIPS;
    typedef tGrDIPS::tVertex  tVGrDIPS;
    typedef tGrDIPS::tEdge    tEGrDIPS;



    class cAppliImagePairSelector: public cMMVII_Appli
    {
        public :
            cAppliImagePairSelector(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec);

        private :
            // photogrammetric project
            cPhotogrammetricProject mPhProj;
            tGrDIPS mGraph;
            std::string mPattenIm;
            std::string mNameOut;
            bool mIsWithTieP; 
            int mNbMinHomol; // minimum number of tie points to consider a pair
            int mCellSize; // cell size for tie point distribution analysis
            int mNbMaxHomol; // maximum number of tie points to keep per edge
            int mKBestImages;
            cPt2dr mTriangulationAngleOpt;

            cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override ;
            cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override ; 

            tREAL8 WeightTriangulationAngle
                (
                        tREAL8 aAngle,                             
                        tREAL8 aAlphaOpt,
                        tREAL8 aHalfWidth 
                );       
                
            bool IsValidPair(tREAL8 & aWeight,const cSensorCamPC * aCam1, const cSensorCamPC * aCam2);
            void AddCamera(cPlyVertices & aPlyV, cSensorCamPC * aCam);
            void AddConnexion(cPlyVertices & aPlyV, const tEGrDIPS * anEdge);
            int Exe() override;
    };



    cAppliImagePairSelector::cAppliImagePairSelector(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec):
        cMMVII_Appli(aVArgs,aSpec),
        mPhProj(*this),
        mPattenIm(""),
        mNameOut(""),
        mIsWithTieP(false),
        mNbMinHomol(10),
        mCellSize(10),
        mNbMaxHomol(1000),
        mKBestImages(10),
        mTriangulationAngleOpt(cPt2dr(12.0, 12.0))  // optimal triangulation angle range in degrees
    {
    }

    cCollecSpecArg2007 & cAppliImagePairSelector::ArgObl(cCollecSpecArg2007 & anArgObl)
    {
        return anArgObl
            << Arg2007(mPattenIm,"Pattern for images",{{eTA2007::MPatFile,"0"}})
            << mPhProj.DPOrient().ArgDirInMand()
            << Arg2007(mNameOut,"Name of output xml file ",{{eTA2007::FileTxt,"1"}})

            ;
    }   

    cCollecSpecArg2007 & cAppliImagePairSelector::ArgOpt(cCollecSpecArg2007 & anArgOpt)
    {
        return anArgOpt
            <<  mPhProj.DPTieP().ArgDirInOpt()
            <<  mPhProj.DPMulTieP().ArgDirInOpt()
            << AOpt2007(mIsWithTieP,"WithTieP","If true, only keep pairs with tie points",{{eTA2007::HDV}})
            << AOpt2007(mNbMinHomol,"NbMinHomol","Minimum number of tie points to consider a pair",{{eTA2007::HDV}})
            << AOpt2007(mCellSize,"CellSize","Cell size for tie point distribution analysis",{{eTA2007::HDV}})
            << AOpt2007(mNbMaxHomol,"NbMaxHomol","Maximum number of tie points to keep per edge",{{eTA2007::HDV}})
            << AOpt2007(mKBestImages,"KBestImages","Number of best images to keep for each image",{{eTA2007::HDV}})
            << AOpt2007(mTriangulationAngleOpt,"TriangulationAngleOpt","Optimal triangulation angle range in degrees",{{eTA2007::HDV}})
            ; 
    }       

     bool cAppliImagePairSelector::IsValidPair(tREAL8 & aWeight, 
                                    const cSensorCamPC * aCam1, 
                                    const cSensorCamPC * aCam2) 
     {
        return true; // --- IGNORE ---

        cPt3dr aBase = aCam2->Center() - aCam1->Center();
        cPt3dr aBaseN = aBase/ Norm2(aBase);

        cPt3dr aV1_Dir = aCam1->Image2Bundle(ToR(aCam1->Sz()/2)).V12();
        aV1_Dir /= Norm2(aV1_Dir);

        cPt3dr aV2_Dir = aCam2->Image2Bundle(ToR(aCam2->Sz()/2)).V12();
        aV2_Dir /=Norm2(aV2_Dir);

        //  1- Base infinitesimal
        if(Norm2(aBase)<1e-2)
            return false;


        // specific case of nadir viewing (aerial and drone imagery) where the cameras are looking straight down, the triangulation angle is small, and the rays are nearly parallel.
        cPt3dr aNadirView = cPt3dr(0,0,-1); // Assuming nadir direction is along negative Z-axis
        
        tREAL8 anOffNadir1 = std::acos(Scal(aV1_Dir,aNadirView));
        tREAL8 anOffNadir2 = std::acos(Scal(aV2_Dir,aNadirView));
        tREAL8 aNadirThreshold = 10 * (M_PI / 180.0); // 10 degrees in radians

        StdOut()<<aCam1->NameImage()<<"  "<<aCam2->NameImage()<<std::endl;

        StdOut()<<" near parallel views "<<Norm2(aV1_Dir ^ aV2_Dir) << " << "<<  std::sin(1 * (M_PI / 180.0)) <<  "\n";
        
        bool isNearParallelViews = Norm2(aV1_Dir ^ aV2_Dir) < std::sin(1 * (M_PI / 180.0)); // Check if the cross product is near zero
    

        if(isNearParallelViews)
        {
            if(anOffNadir1 < aNadirThreshold && 
            anOffNadir2 < aNadirThreshold)
            {
                // Both cameras are looking nearly straight down and their viewing directions are nearly parallel
                // this is the case of nadir aerial imagery, we should have a proxy for depth to estimate overlap or b/h ratio
                return true;
                // later if camera have a proxy for depth, we can compute the overlap and b/h ratio to decide if the pair is valid
            }
            return false; // Discard pairs with near-parallel viewing directions
        }


        // 3- Base not included in either one of the camera frustums 
        if( (Scal(aBaseN,aV1_Dir)<=std::cos(90 * (M_PI / 180.0))) ||
            (Scal(-aBaseN,aV2_Dir)<=std::cos(90 * (M_PI / 180.0))) 
            ) // cameras not looking at same region or heavily skewed, discard
            return false;
        
        // 4- divergent rays
        tREAL8 d = Scal(aV1_Dir,aV2_Dir);
        tREAL8 D = 1 - d*d;

        if (std::abs(D)<1e-6)
            return false;

        tREAL8 e = Scal(aV1_Dir,-aBase);
        tREAL8 f = Scal(aV2_Dir,-aBase);

        tREAL8 t1 = (d*f-e)/(D);
        tREAL8 t2 = (f- d*e)/(D);

        if ((t1<=0) || (t2<=0))
            return false;

        cPt3dr aP1 =  aCam1->Center() + t1 * aV1_Dir;
        cPt3dr aP2 =  aCam2->Center() + t2 * aV2_Dir;

        cPt3dr aPMedian = ( aP1 + aP2 ) / 2.0 ;
        cPt3dr aV1Ray = aCam1->Center() - aPMedian;
        cPt3dr aV2Ray = aCam2->Center() - aPMedian;

        aV1Ray/=Norm2(aV1Ray);
        aV2Ray/=Norm2(aV2Ray);


        tREAL8 aAngleTri = std::acos(Scal(aV1Ray,aV2Ray));
        aWeight = WeightTriangulationAngle(aAngleTri,
                                    mTriangulationAngleOpt.x(),
                                    mTriangulationAngleOpt.y());

        // discard pairs with triangulation angle outside the range of 1.0° to 70°
        if  ( ( aAngleTri < (1.0 * (M_PI / 180.0))) || 
             (aAngleTri > (70 * (M_PI / 180.0))) 
            ) // angle not between 1.0° and 70°
            return false;

        return true;
    }

    tREAL8 cAppliImagePairSelector::WeightTriangulationAngle
           (
                tREAL8 aAngle,
                tREAL8 aAlphaOpt,
                tREAL8 aHalfWidth
           )
                /**  Weight of a triangulation angle for image-pair / view selection.
     *
     *   Maps a triangulation angle to a quality weight in [0,1] using MMVII's
     *   compact-support cubic bump (CubAppGaussVal) :
     *
     *        w(alpha) = CubAppGaussVal( (alpha - aAlphaOpt) / aHalfWidth )
     *
     *   - w == 1 at the optimal angle aAlphaOpt (best depth/matchability tradeoff)
     *   - w decreases smoothly and reaches 0 at aAlphaOpt +/- aHalfWidth
     *   - w == 0 outside that band, so near-parallel pairs (unstable depth) and
     *     very wide pairs (hard to match) are naturally discarded.
     *
     *   All angles are expressed in radians.
     */
    {
        // to radians
        aAlphaOpt *= (M_PI / 180.0);
        aHalfWidth *= (M_PI / 180.0);
        return CubAppGaussVal((aAngle - aAlphaOpt) / aHalfWidth);
    }

    void cAppliImagePairSelector::AddCamera(cPlyVertices & aPlyV, cSensorCamPC * aCam)
    {
        aPlyV.AddVert(aCam->Center(),cPt3dr(1,0,0));

        //create planes to visualize camera frustum
        
        int aSzStep =20;
        cPt2di aSz = aCam->Sz();
        cPt2dr aImStepSz(aSz[0]/aSzStep,aSz[1]/aSzStep);

        double aFPix = aCam->InternalCalib()->F();
        double aDiag = std::sqrt(std::pow(aSz[0],2)+std::pow(aSz[1],2));
        double aF = 0.2 * aDiag / aFPix;


        // add image border
        std::vector<cPt3dr> aImVPts;
        for (int aS=0; aS<=aSzStep; aS++)
        {

            double aDX = aImStepSz[0]*aS-1;
            double aDY = aImStepSz[1]*aS-1;

            // image plane border
            aImVPts.push_back(aCam->ImageAndDepth2Ground( cPt3dr(0,aDY,aF) ));
            aImVPts.push_back(aCam->ImageAndDepth2Ground( cPt3dr(aDX,0,aF) ));
            aImVPts.push_back(aCam->ImageAndDepth2Ground( cPt3dr(aSz[0],aDY,aF) ));
            aImVPts.push_back(aCam->ImageAndDepth2Ground( cPt3dr(aDX,aSz[1],aF) ));

        }

        for (auto & aP : aImVPts)
        {
            aPlyV.AddVert(aP,cPt3dr(0,1,0));
        }
    }

    void cAppliImagePairSelector::AddConnexion(cPlyVertices & aPlyV, const tEGrDIPS * anEdge)
    {
        cSensorCamPC * aCam1 = anEdge->VertexInit().Attr().mSens;
        cSensorCamPC * aCam2 = anEdge->Succ().Attr().mSens;

        cPt3dr aC1 = aCam1->Center();
        cPt3dr aC2 = aCam2->Center();

        cPt3dr aColor= cPt3dr(RandUnif_0_1(),RandUnif_0_1(),RandUnif_0_1());

        //aPlyV.AddVert(aC1,aColor);
        //aPlyV.AddVert(aC2,aColor);

        for (int aK=1; aK<99; aK++)
        {
            cPt3dr aP = aC1 + (aC2 - aC1) * (aK / 100.0);
            aPlyV.AddVert(aP,aColor);
        }
    }

    int cAppliImagePairSelector::Exe()
    {
        mPhProj.FinishInit();
        std::vector<std::string> aVIm = VectMainSet(0);
        // vertexes are cameras
        for ( const auto & aNameIm : aVIm)
        {
            cSensorCamPC * aSensor = mPhProj.ReadCamPC(aNameIm,true);
            cDIPS_AV aAV(aSensor);
            mGraph.NewVertex(aAV);
        }


        // logic to create edges between images (a criterion should be defined)
        for(size_t aK1=0; aK1<mGraph.NbVertex();aK1++)
        {   
            tVGrDIPS & aV1 = mGraph.VertexOfNum(aK1);
            StdOut() << "Finding Neighbors for : " << aV1.Attr().mSens->NameImage() << "\n";

            for( size_t aK2=aK1+1; aK2<mGraph.NbVertex();aK2++)
            {
                tVGrDIPS & aV2 = mGraph.VertexOfNum(aK2);
                cDIPS_Or aOr;
                cDIPS_Sym aSym;
                aSym.mWeightHomol = cPt2dr(0.0,0.0); // Assign a weight of 1.
                aSym.mWeightGeom  = -1.0;
                bool isToAdd = true;

                // initial geometry prior check 
                cSensorCamPC * aCam1 = aV1.Attr().mSens;
                cSensorCamPC * aCam2 = aV2.Attr().mSens;

                if (!IsValidPair(aSym.mWeightGeom, aCam1, aCam2))
                {
                    isToAdd = false;
                }

                if(mIsWithTieP)
                {
                    int aNbInit = 0;
                    mPhProj.ReadHomolMultiSrce(aNbInit,
                                            aSym.mCpleH,
                                            aV1.Attr().mSens->NameImage(),
                                            aV2.Attr().mSens->NameImage());
                    
                    if (aNbInit==0)
                    {
                        MMVII_UnclasseUsEr("No homologous source indicated");
                    }

                    isToAdd = isToAdd && ((int)aSym.mCpleH.NbH() > mNbMinHomol);
                }

                if (isToAdd)
                    {
                        mGraph.AddEdge(aV1, aV2, aOr, aOr, aSym);
                        StdOut() << " Possible neighbor: number of homol " << aSym.mCpleH.NbH() << " for " << aV2.Attr().mSens->NameImage() << "\n";
                        
                    }

            }
        }



        StdOut() << "Graph built with " << mGraph.NbVertex() << " vertices and " << mGraph.AllEdges().size() << " edges\n";

        // Edge pruning based on scoring the distribution of tie points                   
        for( auto & anEdge: mGraph.AllEdges())
        {   
            // add info on distribution of tie points if mIsWithTieP is true
            if( mIsWithTieP)
            {
                OccupancyPyramid aV1Occ(anEdge->VertexInit().Attr().mSens->Sz()); 
                OccupancyPyramid aV2Occ(anEdge->Succ().Attr().mSens->Sz());
                aV1Occ.build(5); // 5 pyramid levels
                aV2Occ.build(5);

                //StdOut()<<"Occupancy pyramid built for edge between " << anEdge->VertexInit().Attr().mSens->NameImage() << " and " << anEdge->Succ().Attr().mSens->NameImage() << "\n";
                // compute a proxy for the triangulation angle using tie points 
                std::vector<tREAL8> aSetOfWeightsGeom;
                aSetOfWeightsGeom.reserve(anEdge->AttrSym().mCpleH.NbH());

                cSetHomogCpleIm aFilteredCpleH = anEdge->AttrSym().mCpleH.SelectOnSpatialCriteria(mNbMaxHomol);

                for( const auto & aCple : aFilteredCpleH.SetH())
                {
                    cPt2dr aP1 = aCple.Pt(true);
                    cPt2dr aP2 = aCple.Pt(false);
                    aV1Occ.InsertInPyrStruct(aP1);
                    aV2Occ.InsertInPyrStruct(aP2);

                    // triangulate
                    tSeg3dr aC1Bundle = anEdge->VertexInit().Attr().mSens->Image2Bundle(aP1);
                    tSeg3dr aC2Bundle= anEdge->Succ().Attr().mSens->Image2Bundle(aP2);
                    cPt3dr aP3D=RobustBundleInters({aC1Bundle,aC2Bundle});

                    cPt3dr aV1Ray = anEdge->VertexInit().Attr().mSens->Center() - aP3D;
                    cPt3dr aV2Ray= anEdge->Succ().Attr().mSens->Center() - aP3D;

                    aV1Ray = aV1Ray/Norm2(aV1Ray);
                    aV2Ray = aV2Ray/Norm2(aV2Ray);

                    tREAL8 aAngleTri = std::acos(Scal(aV1Ray,aV2Ray));
                    aSetOfWeightsGeom.push_back(WeightTriangulationAngle(aAngleTri,
                                                    mTriangulationAngleOpt.x(),
                                                    mTriangulationAngleOpt.y()));

                }
                anEdge->AttrSym().mWeightHomol = cPt2dr(aV1Occ.HierarchicalScore(), aV2Occ.HierarchicalScore());
                anEdge->AttrSym().mWeightGeom = aSetOfWeightsGeom.empty() ? 0 : ConstMediane(aSetOfWeightsGeom);
            }
        }

        // Output the selected pairs to an XML file
        tNameRel aSetOfPairs;
        cPlyVertices aPlyVertices;
        // Now select based on mWeighHomol and mWeightGeom, keep only the best K edges for each vertex
        for(size_t aK1=0; aK1<mGraph.NbVertex();aK1++)
        {

            // Add Camera
            AddCamera(aPlyVertices, mGraph.VertexOfNum(aK1).Attr().mSens);


            tVGrDIPS & aV1 = mGraph.VertexOfNum(aK1);
            std::vector<tEGrDIPS*> aEdges = aV1.EdgesSucc();

            // sort edges based on a combined score of mWeightHomol and mWeightGeom
            std::sort(aEdges.begin(), aEdges.end(), [](tEGrDIPS* e1, tEGrDIPS* e2) {
                tREAL8 score1 = e1->AttrSym().mWeightHomol.x() + e1->AttrSym().mWeightHomol.y() ; //+ e1->AttrSym().mWeightGeom;
                tREAL8 score2 = e2->AttrSym().mWeightHomol.x() + e2->AttrSym().mWeightHomol.y() ; //+ e2->AttrSym().mWeightGeom;
                return score1 > score2; // descending order
            });

            // display sorted edges for debugging
            /*StdOut() << "Sorted edges for " << aV1.Attr().mSens->NameImage() << ":\n";
            for (const auto& edge : aEdges)
            {
                StdOut() << "  Neighbor: " << edge->Succ().Attr().mSens->NameImage() 
                         << ", Score: " << edge->AttrSym().mWeightHomol.x() + edge->AttrSym().mWeightHomol.y() 
                         << ", Geom Weight: " << edge->AttrSym().mWeightGeom << "\n";
            }
            */
            // keep only the best K edges

            for (size_t aK2=0; aK2<std::min(aEdges.size(), static_cast<size_t>(mKBestImages)); aK2++)
            {
                tEGrDIPS* anEdge = aEdges[aK2];
                std::string aName1= anEdge->VertexInit().Attr().mSens->NameImage();
                std::string aName2= anEdge->Succ().Attr().mSens->NameImage();
                aSetOfPairs.Add(tNamePair(aName1, aName2));

                // Add Edge
                AddConnexion(aPlyVertices, anEdge);
            }

        }

        // save aSetOfPairs to XML
        SaveInFile(aSetOfPairs, mNameOut);

        // Save to ply file for visualization
        std::string aNamePly = mNameOut + ".ply";
        aPlyVertices.ToPly(aNamePly,false);

        return EXIT_SUCCESS;
    }


tMMVII_UnikPApli Alloc_DMSelectBestPairs(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec)
{
    return tMMVII_UnikPApli(new cAppliImagePairSelector(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpec_DMSelectBestPairs
    (
        "DMSelectBestPairs",
        Alloc_DMSelectBestPairs,
        "Select the best image pairs based on geometric and tie point distribution criteria",
        {eApF::Ori},
        {eApDT::TieP},
        {eApDT::Orient},
        __FILE__
    );


}; // namespace MMVII
