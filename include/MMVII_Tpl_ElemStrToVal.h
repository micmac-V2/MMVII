#include <map>
#include "MMVII_enums.h"
#include "MMVII_Stringifier.h"
#include "cMMVII_Appli.h"

/** \file uti_e2string.cpp
    \brief Implementation enum <=> string conversion

    Probably sooner or later this file will be generated
automatically. Waiting for that, I do it quick and dirty
by hand.
*/


namespace MMVII
{
    /* =========================================== */
    /*                                             */
    /*              cE2Str<TypeEnum>               */
    /*                                             */
    /* =========================================== */

/// Those will never be automatized


template <class TypeEnum> class cE2Str
{
   public :
     static const std::string & E2Str(const TypeEnum & anE)
     {
         /// In this direction no need to create
         typename tMapE2Str::iterator anIt = mE2S.find(anE);
         // Enum to string is not user error (user do not create enum)
         if (anIt == mE2S.end())
            MMVII_INTERNAL_ASSERT_always(false,"E2Str for enum : " + ToStr(int(anE)) + ", for type: " + cStrIO<TypeEnum>::msNameType);
         return anIt->second;
     }

     //static const TypeEnum &  Str2E(const std::string & aStr,bool WithDef)
     static TypeEnum   Str2E(const std::string & aStr,bool WithDef)
     {
         /// If first time we create mS2E by inverting the  mE2S
         if (mS2E==0)
         {
            // mS2E = new tMapStr2E;
            mS2E.reset(new tMapStr2E);
            for (const auto & it : mE2S)
                (*mS2E)[it.second] = it.first;
         }
         typename tMapStr2E::iterator anIt = mS2E->find(aStr);
         // String to enum is probably a user error (programm create enum)
         if (anIt == mS2E->end())
         {
            if (WithDef)
                return TypeEnum::eNbVals;
            MMVII_UserError(eTyUEr::eBadEnum,"Str2E for : "+aStr+" ; valids are : "+ StrAllVal() );
         }
         return anIt->second;
     }

     static const std::string   StrAllVal()
     {
         std::string aRes;
         for (const auto & it : mE2S)
         {
             if (aRes!="") aRes += ",";
             aRes += it.second;
         }
         return aRes;
     }

     static std::vector<TypeEnum> VecOfPat(const std::string & aPat,bool AcceptEmpy)
     {
          std::vector<TypeEnum> aRes;
          tNameSelector  aSel =  AllocRegex(aPat);
          for (const auto & it : mE2S)
          {
              if (aSel.Match(it.second))
              {
                 aRes.push_back(it.first);
              }
          }
          if ((!AcceptEmpy) && aRes.empty())
          {
             MMVII_UserError
             (
                eTyUEr::eEmptyPattern,
                "No value for enum, allowed are :"+StrAllVall<TypeEnum>()
             );

          }

          return aRes;
     }

     static std::vector<TypeEnum> AllVals() { return VecOfPat(".*",true); }

     static std::vector<bool> VecBoolOfPat(const std::string & aPat,bool AcceptEmpy)
     {
         std::vector<TypeEnum>  aVEnum = VecOfPat(aPat,AcceptEmpy);
         std::vector<bool> aResult(size_t(TypeEnum::eNbVals)+1,false);

         for (const auto & aLab : aVEnum)
                 aResult.at(size_t(aLab)) = true;
         return aResult;
     }

   private :
     typedef std::map<TypeEnum,std::string> tMapE2Str;
     typedef std::map<std::string,TypeEnum> tMapStr2E;

     static tMapE2Str                   mE2S;
     static std::unique_ptr<tMapStr2E > mS2E;
};

#define TPL_ENUM_2_STRING(TypeEnum)\
template<> std::unique_ptr<cE2Str<TypeEnum>::tMapStr2E > cE2Str<TypeEnum>::mS2E = nullptr;\
const std::string & E2Str(const TypeEnum & anOp)\
{\
   return cE2Str<TypeEnum>::E2Str(anOp);\
}\
template <> TypeEnum  Str2E<TypeEnum>(const std::string & aName,bool WithDef)\
{\
   return cE2Str<TypeEnum>::Str2E(aName,WithDef);\
}\
template <> std::string   StrAllVall<TypeEnum>()\
{\
   return cE2Str<TypeEnum>::StrAllVal();\
}\
template <> std::vector<TypeEnum> SubOfPat<TypeEnum>(const std::string & aPat,bool AcceptEmpty)\
{\
   return cE2Str<TypeEnum>::VecOfPat(aPat,AcceptEmpty);\
}\
template <> std::vector<TypeEnum> AllEnumValues<TypeEnum>()\
{\
   return cE2Str<TypeEnum>::AllVals();\
}\
template <> tSemA2007  AC_ListVal<TypeEnum>()\
{\
   return {eTA2007::AllowedValues,StrAllVall<TypeEnum>()};\
}\
template <> std::vector<bool> VBoolOfPat<TypeEnum>(const std::string & aPat,bool AcceptEmpty)\
{\
   return cE2Str<TypeEnum>::VecBoolOfPat(aPat,AcceptEmpty);\
}\


/* ========================== */
/*          cEnumAttr         */
/* ========================== */

template <class TypeEnum> cEnumAttr<TypeEnum>::cEnumAttr(TypeEnum aType,const std::string & anAux) :
   mType (aType),
   mAux  (anAux)
{
}
template <class TypeEnum> cEnumAttr<TypeEnum>::cEnumAttr(TypeEnum aType) :
   cEnumAttr<TypeEnum>(aType,"")
{
}
template <class TypeEnum> TypeEnum            cEnumAttr<TypeEnum>::Type() const {return mType;}
template <class TypeEnum> const std::string & cEnumAttr<TypeEnum>::Aux()  const {return mAux;}


/* ========================== */
/*    cES_PropertyList        */
/* ========================== */

template <class TypeEnum> cES_PropertyList<TypeEnum>::cES_PropertyList(const tAllPairs & aAllPairs) :
   mAllPairs (aAllPairs)
{
}

template <class TypeEnum> const typename  cES_PropertyList<TypeEnum>::tAllPairs & cES_PropertyList<TypeEnum>::AllPairs() const
{
   return mAllPairs;
}


template class cEnumAttr<eTA2007>;
template class cES_PropertyList<eTA2007>;




/* ============================================== */
/*                                                */
/*       cInstReadOneArgCL2007                    */
/*                                                */
/* ============================================== */

template<typename T> struct is_vector : public std::false_type {};

template<typename T, typename A>
struct is_vector<std::vector<T, A>> : public std::true_type {};


template <class Type> void  GlobCheckSize(const Type & ,const std::string & anArg)
{
    MMVII_INTERNAL_ASSERT_always(false,"Check size vect for non vect arg");
}

template <class Type> void  GlobCheckSize(const std::vector<Type> & aVal,const std::string & anArg)
{
    cPt2di aSz = cStrIO<cPt2di>::FromStr(anArg);
    if ((int(aVal.size()) < aSz.x()) || ((int(aVal.size()) > aSz.y())))
    {
       MMVII_UserError(eTyUEr::eBadSize4Vect,"IntervalOk=" + anArg + " Got=" + ToStr(int(aVal.size())));
    }
}


template <class Type> class cInstReadOneArgCL2007 : public cSpecOneArg2007
{
    public :

       void  CheckSize(const std::string & anArg) const override
       {
               GlobCheckSize(mVal,anArg);
       }

       bool IsVector() const override
       {
           return  is_vector<Type>::value;
       }



        void V_InitParam(const std::string & aStr) override
        {
            mVal = cStrIO<Type>::FromStr(aStr);
        }
        cInstReadOneArgCL2007 (Type & aVal,const std::string & aName,const std::string & aCom,const tAllSemPL & aVSem) :
              cSpecOneArg2007(aName,aCom,aVSem),
              mVal         (aVal)
        {
            if constexpr (std::is_enum_v<Type>)
            {
                auto tmp = aVSem;
                tmp.push_back(AC_ListVal<Type>());
                mSemPL = tmp;
            }
        }

        const std::string & NameType() const override
        {
            return  cStrIO<Type>::msNameType;
        }
        void * AdrParam() override {return &mVal;}
        std::string NameValue() const override {return ToStr(mVal);}

    private :
        Type &          mVal;
};


template <class Type> tPtrArg2007 Arg2007(Type & aVal, const std::string & aCom,const cSpecOneArg2007::tAllSemPL & aVSem )
{
   return tPtrArg2007(new cInstReadOneArgCL2007<Type>(aVal,"",aCom,aVSem));
}




template <class Type> tPtrArg2007 AOpt2007(Type & aVal,const std::string & aName, const std::string &aCom,const cSpecOneArg2007::tAllSemPL & aVSem)
{
   return  tPtrArg2007(new cInstReadOneArgCL2007<Type>(aVal,aName,aCom,aVSem));
}

#define MACRO_INSTANTIATE_ARG2007(Type)\
template tPtrArg2007 Arg2007<Type>(Type &, const std::string & aCom,const cSpecOneArg2007::tAllSemPL & aVSem);\
template tPtrArg2007 AOpt2007<Type>(Type &,const std::string & aName, const std::string & aCom,const cSpecOneArg2007::tAllSemPL & aVSem);



/**
    This file contains the implemenation of conversion between strings and
   atomic and some non atomic object
*/
/* ==================================== */
/*                                      */
/*          Enumerated type             */
/*    eOpAff,                           */
/*                                      */
/* ==================================== */

#define MACRO_DECLARE_STRIO_ENUM(ETYPE) template<> const std::string cStrIO<ETYPE>::msNameType;


#define MACRO_INSTANTIATE_STRIO_ENUM(ETYPE,ENAME)\
MACRO_INSTANTIATE_ARG2007(ETYPE)\
TPL_ENUM_2_STRING(ETYPE)\
template <>  std::string cStrIO<ETYPE>::ToStr(const ETYPE & anEnum) { return  E2Str(anEnum); }\
template <>  ETYPE cStrIO<ETYPE>::FromStr(const std::string & aStr) { return Str2E<ETYPE>(aStr); }\
template <>  const std::string cStrIO<ETYPE>::msNameType = ENAME;


} // namespace MMVII
