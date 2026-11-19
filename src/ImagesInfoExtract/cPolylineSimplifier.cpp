#include "MMVII_Ptxd.h"
#include "MMVII_Bench.h"
#include "cMMVII_Appli.h"
#include "MMVII_Matrix.h"

#include <vector>
#include <cmath>
#include <algorithm>

// Suppose que cPt2dr (avec x() et y()) est déjà déclarée (ex: include MMVII).

namespace MMVII {


class cApproxPoly
{
    public :
       cApproxPoly(const std::vector<cPt2dr> &);

     private :
         cPt2dr                    mCenter;
         size_t                     mNbPts;
         std::vector<cPt2dr>        mVPts;
         std::vector<tMatIn2r>      mVMat;


};


cApproxPoly::cApproxPoly(const std::vector<cPt2dr> & aVInit) :
    mCenter (Centroid(aVInit)),
    mNbPts  (aVInit.size())
{
    tMatIn2r aMatCum;
    for (const auto & aPInit : aVInit)
    {
        mVPts.push_back(aPInit-mCenter);
    }



}





namespace ClaudeAppoxPoly {


/**
 * Approximation polygonale d'une polyligne par un sous-ensemble de ses sommets
 * (approche Imai-Iri + programmation dynamique).
 *
 *  - Toute arête i->j avec j-i <= aNMax est autorisée; son coût est
 *    Err(i,j) + aPenalty, où Err(i,j) = sqrt(somme des carrés des distances des points
 *    i..j à la corde [Pi,Pj]) : homogène à une distance, et qui croît avec la longueur
 *    du segment (une somme, pas une moyenne).
 *  - On minimise la somme des coûts des arêtes du chemin 0 -> n-1. aPenalty règle le
 *    compromis : grand => peu de points, petit => polyligne fidèle (beaucoup de points).
 *    Une courbe très chahutée peut donc conserver plus de points que n/aNMax.
 *  - Err d'une corde est obtenu en O(1) grâce aux moments d'ordre 2 cumulés.
 *    Complexité totale : O(n * aNMax).
 */
class cPolylineSimplifier
{
    public :
        /// Précalcule les moments cumulés; les points doivent rester inchangés pendant l'usage
        explicit cPolylineSimplifier(const std::vector<cPt2dr> & aVPts);

        /// Indices (croissants) des sommets conservés; contient toujours 0 et n-1 (si n>0).
        /// aPenalty : coût (en unité de distance) ajouté à chaque segment, règle le compromis
        /// fidélité / nombre de points. aNMax : taille maximale d'un saut entre indices conservés.
        std::vector<int> Simplify(double aPenalty,int aNMax,tREAL8 aWNMin=0,tREAL8 aNMin=3.0) const;

        /// sqrt(somme des carrés des distances des points i..j (inclus) à la corde [Pi,Pj])  (0 <= aI < aJ < n)
        double ErrToChord(int aI,int aJ) const;

    private :
        /// Moments d'ordre 0, 1 et 2 d'un ensemble de points
        struct cMoments
        {
            long double mN  = 0;
            long double mX  = 0;
            long double mY  = 0;
            long double mXX = 0;
            long double mXY = 0;
            long double mYY = 0;
        };

        /// Moments des points d'indices aI..aJ inclus, exprimés dans le repère centré sur le point aI
        cMoments MomentsCenteredOnFirst(int aI,int aJ) const;

        int                   mNbPts;
        double                mX0;      ///< Origine de recentrage (1er point), limite les erreurs numériques
        double                mY0;
        std::vector<cMoments> mCumul;   ///< mCumul[k] = somme des moments des points 0..k-1
};

/* ================================ Implémentation ================================ */

inline cPolylineSimplifier::cPolylineSimplifier(const std::vector<cPt2dr> & aVPts) :
    mNbPts (int(aVPts.size())),
    mX0    (0.0),
    mY0    (0.0),
    mCumul (aVPts.size()+1)
{
    if (mNbPts==0)
       return;

    //  Modif MPD : centroid is more stable to
    //mX0 = aVPts[0].x();
    //mY0 = aVPts[0].y();


    cPt2dr aC0 = Centroid(aVPts);
    mX0 = aC0.x();
    mY0 = aC0.y();

    for (int aK=0 ; aK<mNbPts ; aK++)
    {
        const long double aX = aVPts[aK].x() - mX0;
        const long double aY = aVPts[aK].y() - mY0;

        cMoments aM = mCumul[aK];   // cumul précédent, auquel on ajoute le point K
        aM.mN  += 1;
        aM.mX  += aX;
        aM.mY  += aY;
        aM.mXX += aX*aX;
        aM.mXY += aX*aY;
        aM.mYY += aY*aY;
        mCumul[aK+1] = aM;
    }
}

inline cPolylineSimplifier::cMoments cPolylineSimplifier::MomentsCenteredOnFirst(int aI,int aJ) const
{
    // 1) Moments des points I..J par différence de cumuls, dans le repère d'origine P0
    const cMoments & aB = mCumul[aJ+1];
    const cMoments & aA = mCumul[aI];
    cMoments aM;
    aM.mN  = aB.mN  - aA.mN;
    aM.mX  = aB.mX  - aA.mX;
    aM.mY  = aB.mY  - aA.mY;
    aM.mXX = aB.mXX - aA.mXX;
    aM.mXY = aB.mXY - aA.mXY;
    aM.mYY = aB.mYY - aA.mYY;

    // 2) Changement d'origine : P0 -> Pi  (formules de translation des moments)
    //    Pi est le point d'indice aI, de coordonnées (cX,cY) dans le repère P0 :
    //    il se retrouve par différence de cumuls de l'intervalle [I,I]
    const long double aCx = mCumul[aI+1].mX - mCumul[aI].mX;
    const long double aCy = mCumul[aI+1].mY - mCumul[aI].mY;

    cMoments aRes;
    aRes.mN  = aM.mN;
    aRes.mX  = aM.mX  - aM.mN*aCx;
    aRes.mY  = aM.mY  - aM.mN*aCy;
    aRes.mXX = aM.mXX - 2*aCx*aM.mX + aM.mN*aCx*aCx;
    aRes.mXY = aM.mXY - aCx*aM.mY - aCy*aM.mX + aM.mN*aCx*aCy;
    aRes.mYY = aM.mYY - 2*aCy*aM.mY + aM.mN*aCy*aCy;
    return aRes;
}

inline double cPolylineSimplifier::ErrToChord(int aI,int aJ) const
{
    const cMoments aM = MomentsCenteredOnFirst(aI,aJ);

    // Vecteur directeur de la corde Pi->Pj (dans n'importe quel repère de même orientation)
    const long double aDx = (mCumul[aJ+1].mX - mCumul[aJ].mX) - (mCumul[aI+1].mX - mCumul[aI].mX);
    const long double aDy = (mCumul[aJ+1].mY - mCumul[aJ].mY) - (mCumul[aI+1].mY - mCumul[aI].mY);
    const long double aL2 = aDx*aDx + aDy*aDy;

    long double aSomD2;   // somme des carrés des distances à la corde

    if (aL2 == 0)
    {
        // Cas dégénéré (Pi = Pj) : distance au point Pi, qui est l'origine du repère courant
        aSomD2 = aM.mXX + aM.mYY;
    }
    else
    {
        // Distance à la droite de normale unitaire n = (-Dy,Dx)/L :  d = n.(p-Pi)
        // Somme des d^2 = n^T M n  avec M la matrice d'inertie (non normalisée) centrée sur Pi
        aSomD2 = (aDy*aDy*aM.mXX - 2*aDx*aDy*aM.mXY + aDx*aDx*aM.mYY) / aL2;
    }

    return double(std::sqrt(std::max<long double>(0,aSomD2)));   // max : protège des valeurs négatives dues aux arrondis
}

inline std::vector<int> cPolylineSimplifier::Simplify(double aPenalty,int aNMax,tREAL8 aWNMin,tREAL8 aNMin) const
{
    if (mNbPts==0) return {};
    if (mNbPts==1) return {0};

    aNMax    = std::max(1,aNMax);        // un saut de 1 est toujours possible => un chemin existe toujours
    aPenalty = std::max(0.0,aPenalty);

    // aBestCost[j] = coût minimal d'un chemin 0 -> j, chaque arête i->j coûtant ErrToChord(i,j) + aPenalty
    std::vector<double> aBestCost (mNbPts,0.0);
    std::vector<int>    aPrev     (mNbPts,-1);

    for (int aJ=1 ; aJ<mNbPts ; aJ++)
    {
        double aMinCost = 0.0;
        int    aMinPrev = -1;

        for (int aI=std::max(0,aJ-aNMax) ; aI<aJ ; aI++)
        {
            int aDelta = aJ-aI;
            // Entre voisins il n'y a aucun point intermédiaire : erreur nulle
            const double aErr  = (aDelta==1) ? 0.0 : ErrToChord(aI,aJ);

            const double aCost = aBestCost[aI] + aErr + aPenalty * (1+aWNMin * std::clamp((aNMin-aDelta)/aNMin,0.0,1.0)) ;

            if ((aMinPrev<0) || (aCost<aMinCost))
            {
                aMinCost = aCost;
                aMinPrev = aI;
            }
        }

        aBestCost[aJ] = aMinCost;
        aPrev[aJ]     = aMinPrev;
    }

    // Remontée du chemin de n-1 vers 0, puis retournement
    std::vector<int> aRes;
    for (int aK=mNbPts-1 ; aK>=0 ; aK=aPrev[aK])
        aRes.push_back(aK);
    std::reverse(aRes.begin(),aRes.end());
    return aRes;
}


/* Exemple d'usage :

   std::vector<cPt2dr> aVPts = ...;
   cPolylineSimplifier aSimp(aVPts);
   std::vector<int> aVIndex = aSimp.Simplify(0.5,200);   // pénalité 0.5 par segment, sauts <= 200 points
*/





}; // namespace ClaudeAppoxPoly

using namespace ClaudeAppoxPoly;


std::vector<int> PolynSimplify(const std::vector<cPt2dr> & aVPts, double aPenalty,int aNMax,tREAL8 aWNmin,tREAL8 aNMin)
{
    cPolylineSimplifier aPS(aVPts);

    return aPS.Simplify(aPenalty,aNMax,aWNmin,aNMin);
}

std::vector<int> PolynSimplify(const std::vector<cPt2di> & aVPtsI, double aPenalty,int aNMax,tREAL8 aWNmin,tREAL8 aNMin)
{
   std::vector<cPt2dr> aVPtsR;
   for (const auto & aPix : aVPtsI)
       aVPtsR.push_back(ToR(aPix));

   return PolynSimplify(aVPtsR,aPenalty,aNMax,aWNmin,aNMin);
}

template <class Type,const int Dim> int IndexNearesPoint
    (const std::vector<cPtxd<Type,Dim>> & aVPts,const cPtxd<Type,Dim>& aPt,int aDefInd0=-1,int aDefIndMax=-1)
{
    size_t anInd0 = aDefInd0>=0 ? aDefInd0 : 0;
    size_t anInd1 = aDefInd0>=1 ? aDefInd0 : aVPts.size();

    cWhichMin<int,tREAL8> aWhichMin(-1,1e30);
    for (size_t anInd = anInd0 ; anInd<anInd1 ; anInd++)
        aWhichMin.Add(anInd,SqN2(aPt-aVPts.at(anInd)));

    return aWhichMin.IndexExtre();

}

#define INSTANTIATE_NEAREST_POINT(TYPE,DIM)\
template  int IndexNearesPoint(const std::vector<cPtxd<TYPE,DIM>> &,const cPtxd<TYPE,DIM>&,int,int);

INSTANTIATE_NEAREST_POINT(tREAL8,2)
INSTANTIATE_NEAREST_POINT(tINT4,2)


/**
 * @brief The cOneBenchApproxPoly class
 *
 *  Generate point for test aprox poly :
 *
 *     * the node are on a sinusoid (to be sure that the solution by line is OK)
 *     * inside node, point are pefectlly allined
 */

class  cOneBenchApproxPoly
{
   public :

       cOneBenchApproxPoly();


       tREAL8 mPer;
       tREAL8 mOmega;
       int    mNbNode;
       int    mNbMaxPerInterv;

       cPt2dr PtOnCurve(const tREAL8 anAbs) const;

       std::vector<cPt2dr>  mVPts;
       std::vector<int>     mVIndexNodes;


       void AddPts(const cPt2dr & aPt) {mVPts.push_back(aPt);}
       void AddNode(const cPt2dr & aPt)
       {
           mVIndexNodes.push_back(mVPts.size());
           AddPts(aPt);
       }

       void Show(const  std::vector<int> & );



};


cPt2dr cOneBenchApproxPoly::PtOnCurve(const tREAL8 anX) const
{
    return cPt2dr(anX,std::sin(mOmega*anX)/mOmega);
}

void cOneBenchApproxPoly::Show(const  std::vector<int> &  aVI)
{
    cPolylineSimplifier  aPLS(mVPts);

    for (size_t aK=1 ; aK<aVI.size() ; aK++)
    {
        int aI = aVI.at(aK-1);
        int aJ = aVI.at(aK);

        StdOut() << " " << aI << " " << aJ << " " << aPLS.ErrToChord(aI,aJ) << "\n";
    }
}


cOneBenchApproxPoly::cOneBenchApproxPoly() :
    mPer             (RandInInterval(20,40)),
    mOmega           (2*M_PI/mPer),
    mNbNode          (2), // (RandInInterval(5,20)),
    mNbMaxPerInterv  (RandInInterval(5,15))
{
    tREAL8 aCurAbs = RandInInterval(-mPer,mPer);
    cPt2dr aCurNode = PtOnCurve(aCurAbs);
    AddNode(aCurNode);

    for (int aK=1 ; aK<= mNbNode ; aK++)
    {
        tREAL8 aNextAbs = aCurAbs + RandInInterval(mPer*0.1,mPer*0.5);
        cPt2dr aNexNode = PtOnCurve(aNextAbs);
        int aNbPts = RandInInterval(3,mNbMaxPerInterv);

        for (int aK=1 ; aK<aNbPts ; aK++)
        {
            tREAL8 aW1 = aK / (tREAL8) aNbPts;
            AddPts(aCurNode*(1-aW1) + aNexNode*aW1);
        }

        aCurAbs = aNextAbs;
        aCurNode = aNexNode;
        AddNode(aCurNode);
    }


    std::vector<int>   aVInd = PolynSimplify(mVPts, 1e-4,mNbMaxPerInterv+1,0,1.0);

    MMVII_INTERNAL_ASSERT_bench(aVInd==mVIndexNodes,"Error on PolynSimplify");


}



void BenchApproxPoly(cParamExeBench & aParam)
{
      if (! aParam.NewBench("ApproxPoly")) return;



      for (int aK=0 ; aK<100; aK++)
          cOneBenchApproxPoly();

      aParam.EndBench();
}




}; // namespace MMVII


/* Exemple d'usage :

   std::vector<cPt2dr> aVPts = ...;
   cPolylineSimplifier aSimp(aVPts);
   std::vector<int> aVIndex = aSimp.Simplify(0.5,200);   // RMS <= 0.5, sauts <= 200 points
*/
