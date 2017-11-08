#ifndef FFTPACK_CACHE_INC
#define FFTPACK_CACHE_INC

#include <map>
#include <cassert>

// Forwards
class FactorsKey;
class FactorsArrays;
class FactorsCache;
class FactorsLocker;

class FactorsKey
{
private:
        size_t mCookie;
      
public:
        FactorsKey() : mCookie(0) {}
        FactorsKey(size_t N, bool dblPrecision, bool realFFT) : mCookie(N << 2 | (dblPrecision ? 1 : 0) | (realFFT ? 0 : 2)){}

        bool   isValid() const { return mCookie != 0; }
        bool   isSingle() const { return (mCookie & 1) == 0; }
        bool   isDouble() const { return (mCookie & 1) == 1; }
        bool   isComplexFFT() const { return (mCookie & 2) == 2; }
        bool   isRealFFT() const { return (mCookie & 2) == 0; }
        size_t getSize() const { return mCookie >> 2; }

        bool operator==(const FactorsKey& key) const { return mCookie == key.mCookie; }
        bool operator!=(const FactorsKey& key) const { return mCookie != key.mCookie; }
        bool operator<(const FactorsKey& key) const { return mCookie < key.mCookie; }

        size_t getCost() const { return getSize() << ((isComplexFFT() ? 1 : 0) + (isDouble() ? 1 : 0)); }
        size_t getArraySize() const { return getSize() << (isComplexFFT() ? 1 : 0); }
};

class FactorsArrays
{
        friend class FactorsCache;

private:
        // No copy
        FactorsArrays(const FactorsArrays&);
        FactorsArrays& operator=(const FactorsArrays&);

        FactorsKey      mKey;
        FactorsArrays*  mMoreRU;
        FactorsArrays*  mLessRU;
        FactorsCache*   mCache;
        int             mFacSave[15];
        void*           mWASave;
        int             mLockCount;

        FactorsArrays(const FactorsKey& key);
        ~FactorsArrays();

        void init();

public:
        const int*    getFacSave() const { return mFacSave; }
        const float*  getSingleWASave() const { assert(mKey.isSingle()); return reinterpret_cast<const float*>(mWASave); }
        const double* getDoubleWASave() const { assert(mKey.isDouble()); return reinterpret_cast<const double*>(mWASave); }
        
        FactorsKey getKey() const { return mKey; }
        bool       isSingle() const { return mKey.isSingle(); }
        bool       isDouble() const { return mKey.isDouble(); }
        bool       isRealFFT() const { return mKey.isRealFFT(); }
        bool       isComplexFFT() const { return mKey.isComplexFFT(); }
        size_t     getSize() const { return mKey.getSize(); }
        size_t     getArraySize() const { return mKey.getArraySize(); }
        size_t     getCost() const { return mKey.getCost(); }

        FactorsCache* getCache() const { return mCache; }
};

/**
 * Cache for saved arrays.
 * Note: To make everything thread-safe, protect calls to this class by a mutex. The class itself
 * is already designed to be thread-safe then.
 **/
class FactorsCache
{
private:
        // No copy
        FactorsCache(const FactorsCache&);
        FactorsCache& operator=(const FactorsCache&);

        size_t          mMaxCost;
        size_t          mCurrentCost;
        FactorsArrays*  mMostRU;
        FactorsArrays*  mLastRU;

        std::map<FactorsKey, FactorsArrays*> mMap;

        void addMRU(FactorsArrays* node);
        void remove(FactorsArrays* node);
        void destroy(FactorsArrays* node);

public:
        FactorsCache(size_t maxCost);
        ~FactorsCache();

        /** remove all unlocked entries */
        void clear();

        /** get entry for size N */
        FactorsArrays* lock(size_t N, bool dblPrecision, bool realFFT) { return lock(FactorsKey(N, dblPrecision, realFFT)); }
        FactorsArrays* lock(const FactorsKey& key);

        /** unlock entry */
        void unlock(FactorsArrays* array);
};

class FactorsLocker
{
private:
        // No copy
        FactorsLocker(const FactorsLocker&);
        FactorsLocker& operator=(const FactorsLocker&);

        FactorsArrays* mArrays;

public:
        FactorsLocker(FactorsArrays* arrays) : mArrays(arrays) {}
        ~FactorsLocker() { mArrays->getCache()->unlock(mArrays); }

        const int*     getFacSave() const { return mArrays->getFacSave(); }
        const float*   getSingleWASave() const { return mArrays->getSingleWASave(); }
        const double*  getDoubleWASave() const { return mArrays->getDoubleWASave(); }
        FactorsArrays* getFactorArrays() const { return mArrays; }
        FactorsKey     getKey() const { return mArrays->getKey(); }

        template <class T> const T* getWASave() const { assert(0); return 0; }
        template <> const float*  getWASave<float>() const  { return getSingleWASave(); }
        template <> const double* getWASave<double>() const { return getDoubleWASave(); }
};

#endif // FFTPACK_CACHE_INC