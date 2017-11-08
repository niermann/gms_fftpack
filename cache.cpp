#include "cache.h"
#include "fftpack.h"

using namespace std;

FactorsArrays::FactorsArrays(const FactorsKey& key)
: mKey(key), mCache(0), mWASave(0), mLockCount(0)
{
        init();
}

FactorsArrays::~FactorsArrays()
{
        assert(mCache == 0);
        assert(mLockCount == 0);

        if (mKey.isSingle())
                delete[] reinterpret_cast<float*>(mWASave);
        else
                delete[] reinterpret_cast<double*>(mWASave);
};

void FactorsArrays::init()
{
        assert(mKey.isValid());

        if (mKey.isComplexFFT()) {
                if (mKey.isSingle()) {
                        float *wasave = new float[mKey.getSize() * 2];
                        mWASave = wasave;
                        ::cffti(mKey.getSize(), wasave, mFacSave);
                } else {
                        double *wasave = new double[mKey.getSize() * 2];
                        mWASave = wasave;
                        ::dcffti(mKey.getSize(), wasave, mFacSave);
                }
        } else {
                if (mKey.isSingle()) {
                        float *wasave = new float[mKey.getSize()];
                        mWASave = wasave;
                        ::rffti(mKey.getSize(), wasave, mFacSave);
                } else {
                        double *wasave = new double[mKey.getSize()];
                        mWASave = wasave;
                        ::drffti(mKey.getSize(), wasave, mFacSave);
                }
        }
}

void FactorsCache::addMRU(FactorsArrays* node)
{
        node->mMoreRU = 0;
        node->mLessRU = mMostRU;
        mMostRU = node;
        if (node->mLessRU)
                node->mLessRU->mMoreRU = node;
        else
                mLastRU = node;
}

void FactorsCache::remove(FactorsArrays* node)
{
        // Remove from list
        if (node->mMoreRU)
                node->mMoreRU->mLessRU = node->mLessRU;
        else
                mMostRU = node->mLessRU;
        if (node->mLessRU)
                node->mLessRU->mMoreRU = node->mMoreRU;
        else
                mLastRU = node->mMoreRU;
}

void FactorsCache::destroy(FactorsArrays* node)
{
        assert(node->mLockCount == 0);

        mMap.erase(node->getKey());
        remove(node);

        mCurrentCost -= node->getCost();
        node->mCache = 0;
        delete node;
}

FactorsCache::FactorsCache(size_t maxCost)
: mMaxCost(maxCost), mCurrentCost(0), mMostRU(0), mLastRU(0)
{}

FactorsCache::~FactorsCache()
{
        clear();
}

void FactorsCache::clear()
{
        while (mMostRU)
                destroy(mMostRU);
}

FactorsArrays* FactorsCache::lock(const FactorsKey& key)
{
        // Lookup
        FactorsArrays* node;
        map<FactorsKey, FactorsArrays*>::const_iterator it = mMap.find(key);
        if (it != mMap.end()) {
                node = it->second;
                if (!node->mLockCount)
                        remove(node);
        } else {
                // First reduce overall costs
                size_t cost = key.getCost();
                while (mLastRU && (mMaxCost - mCurrentCost) < cost)
                        destroy(mLastRU);

                // then create new node
                node = new FactorsArrays(key);
                node->mCache = this;
                mCurrentCost += cost;

                // Add to map
                mMap.insert(pair<FactorsKey,FactorsArrays*>(key, node));
        }

        node->mLockCount++;
        return node;
}

void FactorsCache::unlock(FactorsArrays* node)
{
        assert(node->mLockCount > 0);
        node->mLockCount--;
        if (!node->mLockCount)
                addMRU(node);
}
