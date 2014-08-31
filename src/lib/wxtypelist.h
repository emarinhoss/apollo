#ifndef __wxtypelist__h__
#define __wxtypelist__h__

// WarpX includes

// std includes

class WxNullType; // does nothing, just serves as a sentinel

/**
 * WxTypeList class provides a means of defining a list of types
 * (hence the name "typelist"). These are very useful for generic
 * programming, for example can be used to create abstract factories
 * using the WxGenScatterHier template class. The macros below make it
 * easy to generate typelist of a given length quite easily.
 */
template<class T, class U>
struct WxTypeList 
{
    typedef T Head;
    typedef U Tail;
};

// The following macros were generated automatically. They make
// constructing new typelists quite easy. For example, a four element
// typelist can be created using
//
// WX_TYPELIST_4(int, usigned int, char, unsigned char)
//
#define WX_TYPELIST_1(t1)                       \
    WxTypeList<t1, WxNullType>
#define WX_TYPELIST_2(t1, t2)                   \
    WxTypeList<t1, WX_TYPELIST_1(t2) >
#define WX_TYPELIST_3(t1, t2, t3)               \
    WxTypeList<t1, WX_TYPELIST_2(t2, t3) >
#define WX_TYPELIST_4(t1, t2, t3, t4)           \
    WxTypeList<t1, WX_TYPELIST_3(t2, t3, t4) >
#define WX_TYPELIST_5(t1, t2, t3, t4, t5)               \
    WxTypeList<t1, WX_TYPELIST_4(t2, t3, t4, t5) >
#define WX_TYPELIST_6(t1, t2, t3, t4, t5, t6)           \
    WxTypeList<t1, WX_TYPELIST_5(t2, t3, t4, t5, t6) >
#define WX_TYPELIST_7(t1, t2, t3, t4, t5, t6, t7)               \
    WxTypeList<t1, WX_TYPELIST_6(t2, t3, t4, t5, t6, t7) >
#define WX_TYPELIST_8(t1, t2, t3, t4, t5, t6, t7, t8)           \
    WxTypeList<t1, WX_TYPELIST_7(t2, t3, t4, t5, t6, t7, t8) >
#define WX_TYPELIST_9(t1, t2, t3, t4, t5, t6, t7, t8, t9)               \
    WxTypeList<t1, WX_TYPELIST_8(t2, t3, t4, t5, t6, t7, t8, t9) >
#define WX_TYPELIST_10(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10)         \
    WxTypeList<t1, WX_TYPELIST_9(t2, t3, t4, t5, t6, t7, t8, t9, t10) >
#define WX_TYPELIST_11(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11)    \
    WxTypeList<t1, WX_TYPELIST_10(t2, t3, t4, t5, t6, t7, t8, t9, t10, t11) >
#define WX_TYPELIST_12(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12) \
    WxTypeList<t1, WX_TYPELIST_11(t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12) >
#define WX_TYPELIST_13(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13) \
    WxTypeList<t1, WX_TYPELIST_12(t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13) >
#define WX_TYPELIST_14(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14) \
    WxTypeList<t1, WX_TYPELIST_13(t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14) >
#define WX_TYPELIST_15(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15) \
    WxTypeList<t1, WX_TYPELIST_14(t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15) >
#define WX_TYPELIST_16(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16) \
    WxTypeList<t1, WX_TYPELIST_15(t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16) >
#define WX_TYPELIST_17(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17) \
    WxTypeList<t1, WX_TYPELIST_16(t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17) >
#define WX_TYPELIST_18(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18) \
    WxTypeList<t1, WX_TYPELIST_17(t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18) >
#define WX_TYPELIST_19(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19) \
    WxTypeList<t1, WX_TYPELIST_18(t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19) >
#define WX_TYPELIST_20(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20) \
    WxTypeList<t1, WX_TYPELIST_19(t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20) >

/**
 * WxTypeMap can be used to generate a whole class hierachy at compile
 * time. It is used in conjunction with a typelist and a class Unit
 * which takes a single template argument. An object of the WxTypeMap
 * inherits from all classes Unit<T> where T is a member of the
 * supplied typelist.
 */
template<class TList, template <class> class Unit> class WxTypeMap;

/**
 * Specialization 1: Inherit from WxTypeMap generated from the
 * elements of the typelist.
 */
template <class T1, class T2, template <class> class Unit>
class WxTypeMap<WxTypeList<T1, T2>, Unit> :
    public WxTypeMap<T1, Unit> , public WxTypeMap<T2, Unit> 
{
  public:
    template <typename T> struct Rebind 
    {
        typedef Unit<T> Result;
    };
};

/**
 * Specialization 2: Inherit from Unit<AtomicType>
 */
template <class AtomicType, template <class> class Unit>
class WxTypeMap : public Unit<AtomicType> 
{
  public:
    template <typename T> struct Rebind 
    {
        typedef Unit<T> Result;
    };
};

/**
 * Specialization 3: For WxNullType do nothing
 */
template <template <class> class Unit>
class WxTypeMap<WxNullType, Unit> 
{
  public:
    template <typename T> struct Rebind 
    {
        typedef Unit<T> Result;
    };
};

/**
 * Say one has created WxTypeMap from a typelist and class Unit. Now,
 * given a type T, the WxScatterHierField function returns a reference
 * to the class Unit<T> portion of the obj. Obj is of type WxTypeMap.
 */
template <class T, class H>
typename H::template Rebind<T>::Result& wxTypeMapExtract(H& obj)
{
    return obj;
}

template <class T, class H>
const
typename H::template Rebind<T>::Result& wxTypeMapExtract(const H& obj)
{
    return obj;
}

#endif //  __wxtypelist__h__
