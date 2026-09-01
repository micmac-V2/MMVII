#ifndef  _MMVII_Tpl_ElemStrToVal_H_
#define  _MMVII_Tpl_ElemStrToVal_H_

#include <map>
#include "MMVII_enums.h"
#include "MMVII_Stringifier.h"
#include "cMMVII_Appli.h"
//#include "MMVII_Tpl_ElemStrToVal.h"
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
             MMVII_INTERNAL_ASSERT_always(false,"E2Str for enum : " + ToStr(int(anE)) + ", for type: " + cStrIO<TypeEnum>::msNameType());
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
            MMVII_UserError(eTyUEr::eBadEnum,"Str2E for : "+aStr+" ; valids are : " + StrAllVal() );
         }
         return anIt->second;
     }

     static const std::string   StrAllVal()
     {
         std::string aRes;
         char aSep = '[';
         for (const auto & it : mE2S)
         {
             aRes += aSep + it.second;
             aSep=',';
         }
         return aRes + ']';
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
struct is_vector<std::vector<T, A>> : public std::true_type { using value_type = T; };

template<typename T, typename = void> struct is_struct_args : std::false_type {};
template<typename T>
struct is_struct_args<T,std::void_t<decltype(std::declval<T>().ToStr()),
                                     decltype(std::declval<T>().FromStr(std::declval<const std::string&>(), std::declval<bool>())),
                                    decltype(std::declval<T>().NameType()),
                                    decltype(std::declval<T>().Arg2007Fields())
                                     >> : std::true_type
{
};

template<class T> inline constexpr bool is_vector_v = is_vector<T>::value;
template<class T> inline constexpr bool is_struct_args_v = is_struct_args<T>::value;

// Number of leading '[' expected when parsing a fully-bracketed value of type T:
// one level per std::vector<> nesting, plus one more if the innermost type is
// itself a struct-args (its generated FromStr also consumes its own bracket).
template <class T>
struct bracket_depth
{
    static constexpr int value = is_struct_args_v<T> ? 1 : 0;
    using value_type = T;
};

template <class U, class Alloc>
struct bracket_depth<std::vector<U, Alloc>>
{
    static constexpr int value = 1 + bracket_depth<U>::value;
    using value_type = typename bracket_depth<U>::value_type;
};

template<typename T, typename = void> struct is_vector_of_struct_args : std::false_type {};
template<typename T>
struct is_vector_of_struct_args<T,std::enable_if_t<is_vector_v<T>>>
    : is_struct_args<typename T::value_type>
{
};
template<class T> inline constexpr bool is_vector_of_struct_args_v = is_vector_of_struct_args<T>::value;


template <class Type> class cInstReadOneArgCL2007 : public cSpecOneArg2007
{
    public :

       void  CheckSize(const std::string & anArg) const override
       {
           if constexpr(is_vector_v<Type>) {
               cPt2di aSz = cStrIO<cPt2di>::FromStr(anArg);
               if ((int(mVal.size()) < aSz.x()) || ((int(mVal.size()) > aSz.y())))
               {
                   MMVII_UserError(eTyUEr::eBadSize4Vect,"IntervalOk=" + anArg + " Got=" + ToStr(int(mVal.size())));
               }
           } else {
               MMVII_INTERNAL_ASSERT_always(false,"Check size vect for non vect arg");
           }
       }

       bool IsVector() const override
       {
           return  is_vector_v<Type>;
       }

        void V_InitParam(const std::string & aStr, bool aFirstInit) override
        {
            if constexpr (is_vector_v<Type>) {
                constexpr int aBracketDepth = bracket_depth<Type>::value;   // expected leading '[' count for a fully-bracketed value

                size_t n=0;
                while ( n< aStr.size() && aStr[n] == '[') n++;
                //  A repeatable arg is given one element at a time, so it is accepted with one
                //  level of brackets less than the full vector. Reserved to CanRepeat : for an
                //  ordinary vector the value is the whole vector, and taking it as one element
                //  would add it to the default value instead of replacing it
                if (HasType(eTA2007::CanRepeat) && (aBracketDepth == n+1)) {    // singleton
                    if (aFirstInit)                                     // else the default value is kept
                       mVal.clear();
                    mVal.push_back(cStrIO<typename is_vector<Type>::value_type>::FromStr(aStr));
                } else {                                                // else: push_back all values from the string into the vector
                    auto aVals = cStrIO<Type>::FromStr(aStr);
                    if (aFirstInit || !HasType(eTA2007::CanRepeat)) {
                        mVal = aVals;
                    } else {
                        mVal.insert(mVal.end(), aVals.begin(), aVals.end());
                    }
                }
            } else {                                                    // not a vector
                mVal = cStrIO<Type>::FromStr(aStr);
            }
        }
        cInstReadOneArgCL2007 (Type & aVal,const std::string & aName,const std::string & aCom,const tAllSemPL & aVSem) :
            cSpecOneArg2007(aName,aCom,aVSem),
            mVal        (aVal),
            mDefaultVal (aVal)
        {

            if constexpr (is_struct_args_v<Type>)
            {
                mStructFields = std::move(aVal.Arg2007Fields());
            }
            else if constexpr (is_vector_of_struct_args_v<Type>)
            {
                typename Type::value_type aValue{};
                mStructFields = std::move(aValue.Arg2007Fields());
            }
            else if constexpr (std::is_enum_v<Type>)
            {
                auto tmp = aVSem;
                tmp.push_back(AC_ListVal<Type>());
                mSemPL = tmp;
            }
        }

        std::string NameType() const override
        {
            return  cStrIO<Type>::msNameType();
        }
        void * AdrParam() override {return &mVal;}
        std::string DefaultNameValue() const override {return ToS(mDefaultVal);}
        const tVecArg2007 & StructFields() const override { return mStructFields; }

    private :
        Type &       mVal;
        Type         mDefaultVal;
        tVecArg2007  mStructFields;
};


template <class Type> tPtrArg2007 Arg2007(Type & aVal, const std::string & aCom,const cSpecOneArg2007::tAllSemPL & aVSem )
{

   return tPtrArg2007(new cInstReadOneArgCL2007<Type>(aVal,"",aCom,aVSem));
}


template <class Type> tPtrArg2007 AOpt2007(Type & aVal,const std::string & aName, const std::string &aCom,const cSpecOneArg2007::tAllSemPL & aVSem)
{
   return  tPtrArg2007(new cInstReadOneArgCL2007<Type>(aVal,aName,aCom,aVSem));
}

/*
int * NullIntPtr() ;
tPtrArg2007   ArgPureComment(const std::string & aCom)
{
    return tPtrArg2007(new cInstReadOneArgCL2007<int>(*NullIntPtr(),"",aCom,{}));
}

int * NullIntPtr() {return nullptr;}
*/



#define MACRO_INSTANTIATE_ARG2007(Type)\
template tPtrArg2007 Arg2007<Type>(Type &, const std::string & aCom,const cSpecOneArg2007::tAllSemPL & aVSem);\
template tPtrArg2007 AOpt2007<Type>(Type &,const std::string & aName, const std::string & aCom,const cSpecOneArg2007::tAllSemPL & aVSem);


namespace StructuredArg
{

template<class T> inline constexpr bool is_field_semantics_v =
    std::is_same_v<std::decay_t<T>,cArg2007FieldSemantics>;


/// Split the field names of ARG2007_STRUCT_FIELDS, ignoring the commas inside brackets
std::vector<std::string> FieldNames(const std::string & aNames);

///  Apply to the spec of a field the semantics declared by the FieldSem that follows it.
void ApplyFieldSemantics(const tPtrArg2007 & aSpec,const cArg2007FieldSemantics & aFieldSem);


template<typename T>
void ReadOne
(
    const std::vector<std::string> & aValues,
    size_t & aValueIndex,
    const tVecArg2007 & aFields,
    T & aField,
    bool ExceptionOnError
)
{
    if constexpr (!is_field_semantics_v<T>)
    {
        if (aValueIndex<aValues.size())
        {
            if (aValues[aValueIndex].empty())
            {
                if (! aFields[aValueIndex]->HasType(eTA2007::HDV))
                {
                    ExceptionOrError(ExceptionOnError,eTyUEr::eBadSize4Vect,
                                     "Empty argument [" + aFields[aValueIndex]->Name() + "]");
                }
            } else {
                using tField = std::decay_t<T>;
                aField = cStrIO<tField>::FromStr(aValues[aValueIndex],ExceptionOnError);
            }
        }
        else if (! aFields[aValueIndex]->HasType(eTA2007::HDV))
        {
            ExceptionOrError(ExceptionOnError,eTyUEr::eBadSize4Vect,
                             "Missing argument [" + aFields[aValueIndex]->Name() + "]");
        }
        ++aValueIndex;
    }
}

template<typename ...Ts>
void FromStr
(
    const std::string & aStr,
    bool ExceptionOnError,
    const tVecArg2007 & aFields,
    Ts && ... args
)
{
    const auto aValues = cStrIO<std::vector<std::string>>::FromStr(aStr,ExceptionOnError);

    if (aValues.size() > aFields.size())
    {
        ExceptionOrError(ExceptionOnError,eTyUEr::eBadSize4Vect,
                         "Too many arguments [" + aStr+ "]: got "
                             + std::to_string(aValues.size()) + ", expected a maximum of"
                             + std::to_string(aFields.size()));
    }

    size_t aValueIndex = 0;
    (ReadOne(aValues,aValueIndex,aFields,args,ExceptionOnError),...);
}

template<typename ...Ts>
std::string ToStr(const Ts& ...args)
{
    std::string aRes = "[";
    bool first = true;
    auto AddOne = [&aRes,&first](const auto & anArg)
    {
        using tArg = std::decay_t<decltype(anArg)>;
        if constexpr (!is_field_semantics_v<tArg>)
        {
            aRes += first ? "" : ",";
            first = false;
            aRes += ToS(anArg);
        }
    };
    (AddOne(args),...);
    aRes += "]";
    return aRes;
}

template<typename ...Ts>
std::vector<std::string> ListOfTypes(const Ts& ...args)
{
    std::vector<std::string> aV;
    auto AddOne = [&aV](const auto & anArg)
    {
        using tArg = std::decay_t<decltype(anArg)>;
        if constexpr (is_field_semantics_v<tArg>)
        {
            MMVII_INTERNAL_ASSERT_always(!aV.empty(),"FieldSem must follow a field in ARG2007_STRUCT_FIELDS");
            const bool isOptional = std::any_of
            (
                anArg.mValues.begin(),anArg.mValues.end(),
                [](const tSemA2007 & aSem) { return aSem.Type()==eTA2007::HDV; }
            );
            if (isOptional && (aV.back().empty() || (aV.back().back()!='?')))
                aV.back() += "?";
        }
        else
        {
            aV.push_back(cStrIO<tArg>::msNameType());
        }
    };
    (AddOne(args),...);
    return aV;
}



template<typename ...Ts> tVecArg2007 Fields
(
    const std::string & aFieldNames,
    Ts && ... aValues
)
{
    const auto aNames = FieldNames(aFieldNames);
    MMVII_INTERNAL_ASSERT_always(aNames.size()==sizeof...(Ts),"Bad ARG2007_STRUCT_FIELDS tokenization");

    size_t anIndex = 0;
    tPtrArg2007 aSpec;
    tVecArg2007 aResult;

    auto AddOne = [&aResult,&aNames,&anIndex,&aSpec](auto & anArg)
    {
        using tArg = std::decay_t<decltype(anArg)>;
        if constexpr (! is_field_semantics_v<tArg>)
        {
            if (aSpec)
                aResult.push_back(aSpec);
            aSpec = std::make_shared<cInstReadOneArgCL2007<tArg>>(anArg, aNames.at(anIndex), "", cSpecOneArg2007::TheEmptySem);
        }
        else
        {
            ApplyFieldSemantics(aSpec,anArg);
            aResult.push_back(std::move(aSpec));
        }
        ++anIndex;
    };
    (AddOne(aValues),...);
    if (aSpec)
        aResult.push_back(aSpec);
    return aResult;
}

} // namespace StructuredArg

template<class T>
struct cStrIO<T, std::enable_if_t<is_struct_args_v<T>>>
{
    static std::string ToStr(const T& aObj)
    {
        return aObj.ToStr();
    }
    static T FromStr(const std::string & aStr, bool ExceptionOnError = false)
    {
        T t;
        t.FromStr(aStr, ExceptionOnError);
        return t;
    }
    static std::string msNameType() {T t; return t.NameType();}
};



#define ARG2007_STRUCT_FIELDS(...) \
std::string ToStr() const { return StructuredArg::ToStr(__VA_ARGS__) ; } \
void FromStr(const std::string& str, bool ExceptionOnError=false) { \
    StructuredArg::FromStr(str,ExceptionOnError,Arg2007Fields(),__VA_ARGS__); \
} \
std::string NameType() { return ToS(StructuredArg::ListOfTypes(__VA_ARGS__)); } \
tVecArg2007 Arg2007Fields() { return StructuredArg::Fields(#__VA_ARGS__,__VA_ARGS__); }



#define MACRO_INSTANTIATE_STRIO_ENUM(ETYPE,ENAME)\
MACRO_INSTANTIATE_ARG2007(ETYPE)\
TPL_ENUM_2_STRING(ETYPE)\
template <>  std::string cStrIO<ETYPE>::ToStr(const ETYPE & anEnum) { return  E2Str(anEnum); }\
template <>  ETYPE cStrIO<ETYPE>::FromStr(const std::string & aStr, bool ExceptionOnError) { \
    auto val = Str2E<ETYPE>(aStr,true);\
    if (val == ETYPE::eNbVals) {\
        ExceptionOrError(ExceptionOnError,eTyUEr::eBadEnum,\
            "Str2E for : "+aStr+" ; valids are : "+ StrAllVall<ETYPE>() );\
        return ETYPE::eNbVals;\
    }\
    return val;\
}\
template <>  std::string cStrIO<ETYPE>::msNameType() { return ENAME ;}




} // namespace MMVII

#endif // _MMVII_Tpl_ElemStrToVal_H_
