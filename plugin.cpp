#include "plugin.h"
#include "fftpack.h"
#include <vector>

#define _GATANPLUGIN_USE_CLASS_PLUGINMAIN
#include "DMPlugInMain.h"

using namespace Gatan;
using namespace std;

/**
 * Returns plugin version as string.
 **/
static DM_StringToken_1Ref fftpack_version(void)
{
	DM::String str;
	
	PLUG_IN_ENTRY
	
	str = DM::String(FFTPACK_PLUGIN_VERSION);

	PLUG_IN_EXIT

	return str.release();
}

// Workaround as template typedefs are illegal
template <class T> struct fft { 
        typedef void (*func_t)(int N, T* c, T* ch, const T* RESTRICT wasave, const int* RESTRICT facsave);
};

// Does CFFT over axis 0
template <class T>
static void do_fft(typename fft<T>::func_t func, char* destData, const char* srcData, int rank, const size_t* dims, const size_t* strides, T normalization)
{
        FactorsKey fftKey(dims[0], sizeof(T) == sizeof(double), false);
        FactorsLocker fftFactors(g_FactorsCache->lock(fftKey));
        vector<T> workArrays(dims[0]*4);

        int n = rank - 1;
        size_t counts[MAX_RANK];
        counts[n] = dims[n];

        while (true) {
                if (n == 0) {
                        // Copy to lower half of workArrays
                        T* ptr = &workArrays.at(0);
                        for (size_t n = dims[0]; n > 0; --n) {
                                *ptr++ = *reinterpret_cast<const T*>(srcData);
                                *ptr++ = *reinterpret_cast<const T*>(srcData + sizeof(T));
                                srcData += strides[0];
                        }

                        func(dims[0], &workArrays.at(0), &workArrays.at(dims[0]*2),
                                fftFactors.getWASave<T>(), fftFactors.getFacSave());

                        // Copy from lower half workArrays
                        ptr = &workArrays.at(0);
                        for (size_t n = dims[0]; n > 0; --n) {
                                *reinterpret_cast<T*>(destData)             = *ptr++ * normalization;
                                *reinterpret_cast<T*>(destData + sizeof(T)) = *ptr++ * normalization;
                                destData += strides[0];
                        }

                        counts[0] = 0;
                } else if (counts[n]) {
                        --counts[n];
                        --n;
                        counts[n] = dims[n];
                        continue;
                }

                ++n;
                if (n >= rank)
                        break;
                srcData  += strides[n] - dims[n-1] * strides[n-1];
                destData += strides[n] - dims[n-1] * strides[n-1];
        }
}

// Does RFFT over axis 0
// srcData is assumed to be real
// destData is assumed to be complex, on return the lower half is set to the FFT, the upper half
// is the complex conjugate
template <class T>
static void do_rfft(typename fft<T>::func_t func, char* destData, const char* srcData, int rank, const size_t* dims, const size_t* destStrides, const size_t* srcStrides, T normalization)
{
        FactorsKey fftKey(dims[0], sizeof(T) == sizeof(double), true);
        FactorsLocker fftFactors(g_FactorsCache->lock(fftKey));
        vector<T> workArrays(dims[0]*2);

        int n = rank - 1;
        size_t counts[MAX_RANK];
        counts[n] = dims[n];

        while (true) {
                if (n == 0) {
                        // Copy to lower half of workArrays
                        T* ptr = &workArrays.at(0);
                        for (size_t n = dims[0]; n > 0; --n) {
                                *ptr++ = *reinterpret_cast<const T*>(srcData);
                                srcData += srcStrides[0];
                        }

                        func(dims[0], &workArrays.at(0), &workArrays.at(dims[0]),
                                fftFactors.getWASave<T>(), fftFactors.getFacSave());

                        // Copy from lower half workArrays
                        ptr = &workArrays.at(0);

                        // Index 0
                        *reinterpret_cast<T*>(destData)             = *ptr++ * normalization;
                        *reinterpret_cast<T*>(destData + sizeof(T)) = 0;
                        destData += destStrides[0];

                        size_t half = (dims[0] + 1) / 2;
                        for (size_t n = 1; n < half; ++n) {
                                *reinterpret_cast<T*>(destData)             = *ptr++ * normalization;
                                *reinterpret_cast<T*>(destData + sizeof(T)) = *ptr++ * normalization;
                                destData += destStrides[0];
                        }

                        // Index N/2 (only for even N)
                        if ((dims[0] & 1) == 0) {
                                *reinterpret_cast<T*>(destData)             = *ptr++ * normalization;
                                *reinterpret_cast<T*>(destData + sizeof(T)) = 0;
                        }
                        destData += destStrides[0] * (dims[0] - half);

                        counts[0] = 0;
                } else if (counts[n]) {
                        --counts[n];
                        --n;
                        counts[n] = dims[n];
                        continue;
                }

                ++n;
                if (n >= rank)
                        break;
                srcData  += srcStrides[n] - dims[n-1] * srcStrides[n-1];
                destData += destStrides[n] - dims[n-1] * destStrides[n-1];
        }
}

// Does RIFFT over axis 0
// srcData is assumed to be complex, only the lower half is used
// destData is assumed to be real
template <class T>
static void do_rifft(typename fft<T>::func_t func, char* destData, const char* srcData, int rank, const size_t* dims, const size_t* destStrides, const size_t* srcStrides, T normalization)
{
        FactorsKey fftKey(dims[0], sizeof(T) == sizeof(double), true);
        FactorsLocker fftFactors(g_FactorsCache->lock(fftKey));
        vector<T> workArrays(dims[0]*2);

        int n = rank - 1;
        size_t counts[MAX_RANK];
        counts[n] = dims[n];

        while (true) {
                if (n == 0) {
                        // Copy to lower half of workArrays
                        T* ptr = &workArrays.at(0);
                        *ptr++ = *reinterpret_cast<const T*>(srcData);
                        srcData += srcStrides[0];

                        size_t half = (dims[0] + 1) / 2;
                        for (size_t n = 1; n < half; ++n) {
                                *ptr++ = *reinterpret_cast<const T*>(srcData);
                                *ptr++ = *reinterpret_cast<const T*>(srcData + sizeof(T));
                                srcData += srcStrides[0];
                        }

                        // Index N/2 (only for even N)
                        if ((dims[0] & 1) == 0)
                                *ptr++ = *reinterpret_cast<const T*>(srcData);
                        srcData += srcStrides[0] * (dims[0] - half);
                        
                        func(dims[0], &workArrays.at(0), &workArrays.at(dims[0]),
                                fftFactors.getWASave<T>(), fftFactors.getFacSave());

                        // Copy from lower half workArrays
                        ptr = &workArrays.at(0);
                        for (size_t n = dims[0]; n > 0; --n) {
                                *reinterpret_cast<T*>(destData) = *ptr++ * normalization;
                                destData += destStrides[0];
                        }

                        counts[0] = 0;
                } else if (counts[n]) {
                        --counts[n];
                        --n;
                        counts[n] = dims[n];
                        continue;
                }

                ++n;
                if (n >= rank)
                        break;
                srcData  += srcStrides[n] - dims[n-1] * srcStrides[n-1];
                destData += destStrides[n] - dims[n-1] * destStrides[n-1];
        }
}

// Unpacks hermite data over axis 0
// Data expected to be complex
template <class T>
static void unpack_hermite1(char* data, int rank, const size_t* dims, const size_t* strides)
{
        int n = rank - 1;
        size_t counts[MAX_RANK];
        counts[n] = dims[n];

        while (true) {
                if (n == 0) {
                        size_t half = (dims[0] + 1) / 2;

                        const char* srcPtr = data + strides[0];
                        char* destPtr = data + (dims[0] - 1) * strides[0];

                        for (size_t m = 1; m < half; ++m) {
                                T re = *reinterpret_cast<const T*>(srcPtr);
                                T im = *reinterpret_cast<const T*>(srcPtr + sizeof(T));
                                *reinterpret_cast<T*>(destPtr)             = re;
                                *reinterpret_cast<T*>(destPtr + sizeof(T)) = -im;
                                srcPtr  += strides[0];
                                destPtr -= strides[0];
                        }

                        data += strides[0] * dims[0];
                        counts[0] = 0;
                } else if (counts[n]) {
                        --counts[n];
                        --n;
                        counts[n] = dims[n];
                        continue;
                }

                ++n;
                if (n >= rank)
                        break;
                data  += strides[n] - dims[n-1] * strides[n-1];
        }
}

// Unpacks hermite data over axis 0 and 1
// Data expected to be complex
template <class T>
static void unpack_hermite2(char* data, int rank, const size_t* dims, const size_t* strides)
{
        assert(rank >= 2);

        int n = rank - 1;
        size_t counts[MAX_RANK];
        counts[n] = dims[n];

        while (true) {
                if (n == 1) {
                        size_t half = (dims[0] + 1) / 2;

                        for (size_t s = dims[1]; s > 0; --s) {
                                const char* srcPtr = data + strides[0];
                                char* destPtr = data + (dims[0] - 1) * strides[0];

                                if (s > 1) {
                                        srcPtr  += strides[1] * (s - 1);
                                        destPtr += (dims[1] - s + 1) * strides[1];
                                }

                                for (size_t m = 1; m < half; ++m) {
                                        T re = *reinterpret_cast<const T*>(srcPtr);
                                        T im = *reinterpret_cast<const T*>(srcPtr + sizeof(T));
                                        *reinterpret_cast<T*>(destPtr)             = re;
                                        *reinterpret_cast<T*>(destPtr + sizeof(T)) = -im;
                                        srcPtr  += strides[0];
                                        destPtr -= strides[0];
                                }
                        }

                        data += strides[1] * dims[1];
                        counts[1] = 0;
                } else if (counts[n]) {
                        --counts[n];
                        --n;
                        counts[n] = dims[n];
                        continue;
                }

                ++n;
                if (n >= rank)
                        break;
                data  += strides[n] - dims[n-1] * strides[n-1];
        }
}

// Unpacks hermite data over axis 0, 1 and 2
// Data expected to be complex
template <class T>
static void unpack_hermite3(char* data, int rank, const size_t* dims, const size_t* strides)
{
        assert(rank >= 3);

        int n = rank - 1;
        size_t counts[MAX_RANK];
        counts[n] = dims[n];

        while (true) {
                if (n == 2) {
                        size_t half = (dims[0] + 1) / 2;

                        for (size_t t = dims[2]; t > 0; --t) {
                                for (size_t s = dims[1]; s > 0; --s) {
                                        const char* srcPtr = data + strides[0];
                                        char* destPtr = data + (dims[0] - 1) * strides[0];

                                        if (s > 1) {
                                                srcPtr  += strides[1] * (s - 1);
                                                destPtr += (dims[1] - s + 1) * strides[1];
                                        }
                                        if (t > 1) {
                                                srcPtr  += strides[2] * (t - 1);
                                                destPtr += (dims[2] - t + 1) * strides[2];
                                        }

                                        for (size_t m = 1; m < half; ++m) {
                                                T re = *reinterpret_cast<const T*>(srcPtr);
                                                T im = *reinterpret_cast<const T*>(srcPtr + sizeof(T));
                                                *reinterpret_cast<T*>(destPtr)             = re;
                                                *reinterpret_cast<T*>(destPtr + sizeof(T)) = -im;
                                                srcPtr  += strides[0];
                                                destPtr -= strides[0];
                                        }
                                }
                        }

                        data += strides[2] * dims[2];
                        counts[2] = 0;
                } else if (counts[n]) {
                        --counts[n];
                        --n;
                        counts[n] = dims[n];
                        continue;
                }

                ++n;
                if (n >= rank)
                        break;
                data  += strides[n] - dims[n-1] * strides[n-1];
        }
}

/**
 * Does a forward FFT along the axis @a axis (0=X, 1=Y, ...).
 * @param src Source image (must be complex)
 * @param axis Axis along the which the FFT is done (0=X, 1=Y, ...).
 * @returns Complex image containing the FFT
 **/
DM_ImageToken_1Ref fftpack_fft(DM_ImageToken srcToken, long axis)
{
        DM::Image destImage;

        PLUG_IN_ENTRY

        // Test input parameters
        DM::Image srcImage(srcToken);
        int dataType = DM::ImageGetDataType(srcImage);
        if (dataType != ImageData::COMPLEX8_DATA && dataType != ImageData::COMPLEX16_DATA)
                ThrowString("Source image has invalid data type.");
        int rank = DM::ImageGetNumDimensions(srcImage);
        if (rank < 1 || rank > MAX_RANK)
                ThrowString("Source image has too many or too few dimensions.");
        if (axis < 0)
                axis += rank;
        if (axis < 0 || axis >= rank)
                ThrowString("Invalid axes selected.");
        size_t dim[MAX_RANK];
        for (int n = 0; n < rank; n++)
                dim[n] = DM::ImageGetDimensionSize(srcImage, n);

	// Create destination image
        switch (rank) {
        case 1: destImage = DM::NewImage("", dataType, dim[0]); break;
        case 2: destImage = DM::NewImage("", dataType, dim[0], dim[1]); break;
        case 3: destImage = DM::NewImage("", dataType, dim[0], dim[1], dim[2]); break;
        case 4: destImage = DM::NewImage("", dataType, dim[0], dim[1], dim[2], dim[3]); break;
        };

        // Prepare (FFT is carried out over work_dim[0])
        size_t stride = DM::ImageGetDataElementByteSize(srcImage);
        size_t work_strides[MAX_RANK];
        size_t work_dim[MAX_RANK];
        for (int n = 0, m = 1; n < rank; n++) {
                if (n == axis) {
                        work_strides[0] = stride;
                        work_dim[0]     = dim[n];
                } else {
                        work_strides[m] = stride;
                        work_dim[m]     = dim[n];
                        ++m;
                }
        	stride *= dim[n];
        }

        
	// prepare to manipulate the data..
        {
                PlugIn::ImageDataLocker srcLocker(srcImage);
                PlugIn::ImageDataLocker destLocker(destImage);
	        const char *srcData = (const char*)srcLocker.get();
                char *destData      = (char*)destLocker.get();

                if (dataType == ImageData::COMPLEX8_DATA)
                        do_fft<float>(::cfftf, destData, srcData, rank, work_dim, work_strides, 1.0);
                else
                        do_fft<double>(::dcfftf, destData, srcData, rank, work_dim, work_strides, 1.0);

		// tell DigitalMicrograph the image is changed...
		destImage.DataChanged();
	}
	
	PLUG_IN_EXIT

	return destImage.release();	
}

DM_ImageToken_1Ref fftpack_fft_axis0(DM_ImageToken srcToken)
{
        return fftpack_fft(srcToken, 0);
}

/**
 * Does a backward FFT
 * @param src Source image (must be complex)
 * @param axis Axis along the which the FFT is done (0=X, 1=Y, ...).
 * @returns Complex image containing the FFT
 **/
DM_ImageToken_1Ref fftpack_ifft(DM_ImageToken srcToken, long axis)
{
        DM::Image destImage;

        PLUG_IN_ENTRY

        // Test input parameters
        DM::Image srcImage(srcToken);
        int dataType = DM::ImageGetDataType(srcImage);
        if (dataType != ImageData::COMPLEX8_DATA && dataType != ImageData::COMPLEX16_DATA)
                ThrowString("Source image has invalid data type.");
        int rank = DM::ImageGetNumDimensions(srcImage);
        if (rank < 1 || rank > MAX_RANK)
                ThrowString("Source image has too many or too few dimensions.");
        if (axis < 0)
                axis += rank;
        if (axis < 0 || axis >= rank)
                ThrowString("Invalid axes selected.");
        size_t dim[MAX_RANK];
        for (int n = 0; n < rank; n++)
                dim[n] = DM::ImageGetDimensionSize(srcImage, n);

	// Create destination image
        switch (rank) {
        case 1: destImage = DM::NewImage("", dataType, dim[0]); break;
        case 2: destImage = DM::NewImage("", dataType, dim[0], dim[1]); break;
        case 3: destImage = DM::NewImage("", dataType, dim[0], dim[1], dim[2]); break;
        case 4: destImage = DM::NewImage("", dataType, dim[0], dim[1], dim[2], dim[3]); break;
        };

        // Prepare (FFT is carried out over work_dim[0])
        size_t stride = DM::ImageGetDataElementByteSize(srcImage);
        size_t work_strides[MAX_RANK];
        size_t work_dim[MAX_RANK];
        for (int n = 0, m = 1; n < rank; n++) {
                if (n == axis) {
                        work_strides[0] = stride;
                        work_dim[0]     = dim[n];
                } else {
                        work_strides[m] = stride;
                        work_dim[m]     = dim[n];
                        ++m;
                }
        	stride *= dim[n];
        }
        
	// prepare to manipulate the data..
        {
                PlugIn::ImageDataLocker srcLocker(srcImage);
                PlugIn::ImageDataLocker destLocker(destImage);
                const char *srcData = (const char*)srcLocker.get();
                char *destData      = (char*)destLocker.get();

                if (dataType == ImageData::COMPLEX8_DATA)
                        do_fft<float>(::cfftb, destData, srcData, rank, work_dim, work_strides, 1.0f / work_dim[0]);
                else
                        do_fft<double>(::dcfftb, destData, srcData, rank, work_dim, work_strides, 1.0 / work_dim[0]);

		// tell DigitalMicrograph the image is changed...
		destImage.DataChanged();
	}
	
	PLUG_IN_EXIT

	return destImage.release();	
}

DM_ImageToken_1Ref fftpack_ifft_axis0(DM_ImageToken srcToken)
{
        return fftpack_ifft(srcToken, 0);
}

/**
 * Does a 2D forward FFT along (always along X and Y).
 * @param src Source image (must be complex)
 * @returns Complex image containing the FFT
 **/
DM_ImageToken_1Ref fftpack_fft2(DM_ImageToken srcToken)
{
        DM::Image destImage;

        PLUG_IN_ENTRY

        // Test input parameters
        DM::Image srcImage(srcToken);
        int dataType = DM::ImageGetDataType(srcImage);
        if (dataType != ImageData::COMPLEX8_DATA && dataType != ImageData::COMPLEX16_DATA)
                ThrowString("Source image has invalid data type.");
        int rank = DM::ImageGetNumDimensions(srcImage);
        if (rank < 2 || rank > MAX_RANK)
                ThrowString("Source image has too many or too few dimensions.");
        size_t dim[MAX_RANK];
        for (int n = 0; n < rank; n++)
                dim[n] = DM::ImageGetDimensionSize(srcImage, n);

	// Create destination image
        switch (rank) {
        case 1: destImage = DM::NewImage("", dataType, dim[0]); break;
        case 2: destImage = DM::NewImage("", dataType, dim[0], dim[1]); break;
        case 3: destImage = DM::NewImage("", dataType, dim[0], dim[1], dim[2]); break;
        case 4: destImage = DM::NewImage("", dataType, dim[0], dim[1], dim[2], dim[3]); break;
        };

 	// Prepare work array
	size_t stride = DM::ImageGetDataElementByteSize(srcImage);
	size_t work_strides[MAX_RANK];
	size_t work_dim[MAX_RANK];
	for (int n = 0; n < rank; n++) {
		work_strides[n] = stride;
		work_dim[n]     = dim[n];
        	stride *= dim[n];
	}

	// prepare to manipulate the data..
        {
                PlugIn::ImageDataLocker srcLocker(srcImage);
                PlugIn::ImageDataLocker destLocker(destImage);
	        const char *srcData = (const char*)srcLocker.get();
                char *destData      = (char*)destLocker.get();

		// Do FFT over X
                if (dataType == ImageData::COMPLEX8_DATA)
                        do_fft<float>(::cfftf, destData, srcData, rank, work_dim, work_strides, 1.0f);
                else
                        do_fft<double>(::dcfftf, destData, srcData, rank, work_dim, work_strides, 1.0);

		// Switch dim 0 and 1
		swap(work_strides[0], work_strides[1]);
		swap(work_dim[0], work_dim[1]);

		// Do FFT over Y
                if (dataType == ImageData::COMPLEX8_DATA)
                        do_fft<float>(::cfftf, destData, destData, rank, work_dim, work_strides, 1.0f);
                else
                        do_fft<double>(::dcfftf, destData, destData, rank, work_dim, work_strides, 1.0);

		// tell DigitalMicrograph the image is changed...
		destImage.DataChanged();
	}
	
	PLUG_IN_EXIT

	return destImage.release();	
}

/**
 * Does a 2D backward FFT along (always along X and Y).
 * @param src Source image (must be complex)
 * @returns Complex image containing the FFT
 **/
DM_ImageToken_1Ref fftpack_ifft2(DM_ImageToken srcToken)
{
        DM::Image destImage;

        PLUG_IN_ENTRY

        // Test input parameters
        DM::Image srcImage(srcToken);
        int dataType = DM::ImageGetDataType(srcImage);
        if (dataType != ImageData::COMPLEX8_DATA && dataType != ImageData::COMPLEX16_DATA)
                ThrowString("Source image has invalid data type.");
        int rank = DM::ImageGetNumDimensions(srcImage);
        if (rank < 2 || rank > MAX_RANK)
                ThrowString("Source image has too many or too few dimensions.");
        size_t dim[MAX_RANK];
        for (int n = 0; n < rank; n++)
                dim[n] = DM::ImageGetDimensionSize(srcImage, n);

	// Create destination image
        switch (rank) {
        case 1: destImage = DM::NewImage("", dataType, dim[0]); break;
        case 2: destImage = DM::NewImage("", dataType, dim[0], dim[1]); break;
        case 3: destImage = DM::NewImage("", dataType, dim[0], dim[1], dim[2]); break;
        case 4: destImage = DM::NewImage("", dataType, dim[0], dim[1], dim[2], dim[3]); break;
        };

 	// Prepare work array
	size_t stride = DM::ImageGetDataElementByteSize(srcImage);
	size_t work_strides[MAX_RANK];
	size_t work_dim[MAX_RANK];
	for (int n = 0; n < rank; n++) {
		work_strides[n] = stride;
		work_dim[n]     = dim[n];
        	stride *= dim[n];
	}

	// prepare to manipulate the data..
        {
                PlugIn::ImageDataLocker srcLocker(srcImage);
                PlugIn::ImageDataLocker destLocker(destImage);
	        const char *srcData = (const char*)srcLocker.get();
                char *destData      = (char*)destLocker.get();

		// Do FFT over X
                if (dataType == ImageData::COMPLEX8_DATA)
                        do_fft<float>(::cfftb, destData, srcData, rank, work_dim, work_strides, 1.0f);
                else
                        do_fft<double>(::dcfftb, destData, srcData, rank, work_dim, work_strides, 1.0);

		// Switch dim 0 and 1
		swap(work_strides[0], work_strides[1]);
		swap(work_dim[0], work_dim[1]);

		// Do FFT over Y
                if (dataType == ImageData::COMPLEX8_DATA)
			do_fft<float>(::cfftb, destData, destData, rank, work_dim, work_strides, 1.0f / work_dim[0] / work_dim[1]);
                else
                        do_fft<double>(::dcfftb, destData, destData, rank, work_dim, work_strides, 1.0 / work_dim[0] / work_dim[1]);

		// tell DigitalMicrograph the image is changed...
		destImage.DataChanged();
	}
	
	PLUG_IN_EXIT

	return destImage.release();	
}

/**
 * Does a 3D forward FFT along (always along X, Y and Z).
 * @param src Source image (must be complex)
 * @returns Complex image containing the FFT
 **/
DM_ImageToken_1Ref fftpack_fft3(DM_ImageToken srcToken)
{
        DM::Image destImage;

        PLUG_IN_ENTRY

        // Test input parameters
        DM::Image srcImage(srcToken);
        int dataType = DM::ImageGetDataType(srcImage);
        if (dataType != ImageData::COMPLEX8_DATA && dataType != ImageData::COMPLEX16_DATA)
                ThrowString("Source image has invalid data type.");
        int rank = DM::ImageGetNumDimensions(srcImage);
        if (rank < 3 || rank > MAX_RANK)
                ThrowString("Source image has too many or too few dimensions.");
        size_t dim[MAX_RANK];
        for (int n = 0; n < rank; n++)
                dim[n] = DM::ImageGetDimensionSize(srcImage, n);

	// Create destination image
        switch (rank) {
        case 1: destImage = DM::NewImage("", dataType, dim[0]); break;
        case 2: destImage = DM::NewImage("", dataType, dim[0], dim[1]); break;
        case 3: destImage = DM::NewImage("", dataType, dim[0], dim[1], dim[2]); break;
        case 4: destImage = DM::NewImage("", dataType, dim[0], dim[1], dim[2], dim[3]); break;
        };

 	// Prepare work array
	size_t stride = DM::ImageGetDataElementByteSize(srcImage);
	size_t work_strides[MAX_RANK];
	size_t work_dim[MAX_RANK];
	for (int n = 0; n < rank; n++) {
		work_strides[n] = stride;
		work_dim[n]     = dim[n];
        	stride *= dim[n];
	}

	// prepare to manipulate the data..
        {
                PlugIn::ImageDataLocker srcLocker(srcImage);
                PlugIn::ImageDataLocker destLocker(destImage);
	        const char *srcData = (const char*)srcLocker.get();
                char *destData      = (char*)destLocker.get();

		// Do FFT over X
                if (dataType == ImageData::COMPLEX8_DATA)
                        do_fft<float>(::cfftf, destData, srcData, rank, work_dim, work_strides, 1.0f);
                else
                        do_fft<double>(::dcfftf, destData, srcData, rank, work_dim, work_strides, 1.0);

		// Switch dim 0 and 1
		swap(work_strides[0], work_strides[1]);
		swap(work_dim[0], work_dim[1]);

		// Do FFT over Y
                if (dataType == ImageData::COMPLEX8_DATA)
                        do_fft<float>(::cfftf, destData, destData, rank, work_dim, work_strides, 1.0f);
                else
                        do_fft<double>(::dcfftf, destData, destData, rank, work_dim, work_strides, 1.0);

		// Switch dim 0 and 2
		swap(work_strides[0], work_strides[2]);
		swap(work_dim[0], work_dim[2]);

		// Do FFT over Z
                if (dataType == ImageData::COMPLEX8_DATA)
                        do_fft<float>(::cfftf, destData, destData, rank, work_dim, work_strides, 1.0f);
                else
                        do_fft<double>(::dcfftf, destData, destData, rank, work_dim, work_strides, 1.0);

		// tell DigitalMicrograph the image is changed...
		destImage.DataChanged();
	}
	
	PLUG_IN_EXIT

	return destImage.release();	
}

/**
 * Does a 3D backward FFT along (always along X, Y and Z).
 * @param src Source image (must be complex)
 * @returns Complex image containing the FFT
 **/
DM_ImageToken_1Ref fftpack_ifft3(DM_ImageToken srcToken)
{
        DM::Image destImage;

        PLUG_IN_ENTRY

        // Test input parameters
        DM::Image srcImage(srcToken);
        int dataType = DM::ImageGetDataType(srcImage);
        if (dataType != ImageData::COMPLEX8_DATA && dataType != ImageData::COMPLEX16_DATA)
                ThrowString("Source image has invalid data type.");
        int rank = DM::ImageGetNumDimensions(srcImage);
        if (rank < 3 || rank > MAX_RANK)
                ThrowString("Source image has too many or too few dimensions.");
        size_t dim[MAX_RANK];
        for (int n = 0; n < rank; n++)
                dim[n] = DM::ImageGetDimensionSize(srcImage, n);

	// Create destination image
        switch (rank) {
        case 1: destImage = DM::NewImage("", dataType, dim[0]); break;
        case 2: destImage = DM::NewImage("", dataType, dim[0], dim[1]); break;
        case 3: destImage = DM::NewImage("", dataType, dim[0], dim[1], dim[2]); break;
        case 4: destImage = DM::NewImage("", dataType, dim[0], dim[1], dim[2], dim[3]); break;
        };

 	// Prepare work array
	size_t stride = DM::ImageGetDataElementByteSize(srcImage);
	size_t work_strides[MAX_RANK];
	size_t work_dim[MAX_RANK];
	for (int n = 0; n < rank; n++) {
		work_strides[n] = stride;
		work_dim[n]     = dim[n];
        	stride *= dim[n];
	}

	// prepare to manipulate the data..
        {
                PlugIn::ImageDataLocker srcLocker(srcImage);
                PlugIn::ImageDataLocker destLocker(destImage);
	        const char *srcData = (const char*)srcLocker.get();
                char *destData      = (char*)destLocker.get();

		// Do FFT over X
                if (dataType == ImageData::COMPLEX8_DATA)
                        do_fft<float>(::cfftb, destData, srcData, rank, work_dim, work_strides, 1.0f);
                else
                        do_fft<double>(::dcfftb, destData, srcData, rank, work_dim, work_strides, 1.0);

		// Switch dim 0 and 1
		swap(work_strides[0], work_strides[1]);
		swap(work_dim[0], work_dim[1]);

		// Do FFT over Y
                if (dataType == ImageData::COMPLEX8_DATA)
			do_fft<float>(::cfftb, destData, destData, rank, work_dim, work_strides, 1.0f);
                else
                        do_fft<double>(::dcfftb, destData, destData, rank, work_dim, work_strides, 1.0);

		// Switch dim 0 and 2
		swap(work_strides[0], work_strides[2]);
		swap(work_dim[0], work_dim[2]);

		// Do FFT over Z
                if (dataType == ImageData::COMPLEX8_DATA)
			do_fft<float>(::cfftb, destData, destData, rank, work_dim, work_strides, 1.0f / work_dim[0] / work_dim[1] / work_dim[2]);
                else
                        do_fft<double>(::dcfftb, destData, destData, rank, work_dim, work_strides, 1.0 / work_dim[0] / work_dim[1] / work_dim[2]);

		// tell DigitalMicrograph the image is changed...
		destImage.DataChanged();
	}
	
	PLUG_IN_EXIT

	return destImage.release();	
}

/**
 * Does a forward RFFT along the axis @a axis (0=X, 1=Y, ...).
 * 
 * @param src Source image (must be floating point type)
 * @param axis Axis along the which the FFT is done (0=X, 1=Y, ...).
 * @returns Complex image containing the FFT
 **/
DM_ImageToken_1Ref fftpack_rfft(DM_ImageToken srcToken, long axis)
{
        DM::Image destImage;

        PLUG_IN_ENTRY

        // Test input parameters
        DM::Image srcImage(srcToken);
        int dataType = DM::ImageGetDataType(srcImage);
        int destType;
        switch (dataType) {
        case ImageData::REAL4_DATA:
                destType = ImageData::COMPLEX8_DATA;
                break;
        case ImageData::REAL8_DATA:
                destType = ImageData::COMPLEX16_DATA;
                break;
        default:
                ThrowString("Source image has invalid data type.");
                break;
        }
        int rank = DM::ImageGetNumDimensions(srcImage);
        if (rank < 1 || rank > MAX_RANK)
                ThrowString("Source image has too many or too few dimensions.");
        if (axis < 0)
                axis += rank;
        if (axis < 0 || axis >= rank)
                ThrowString("Invalid axes selected.");
        size_t dim[MAX_RANK];
        for (int n = 0; n < rank; n++)
                dim[n] = DM::ImageGetDimensionSize(srcImage, n);

	// Create destination image
        switch (rank) {
        case 1: destImage = DM::NewImage("", destType, dim[0]); break;
        case 2: destImage = DM::NewImage("", destType, dim[0], dim[1]); break;
        case 3: destImage = DM::NewImage("", destType, dim[0], dim[1], dim[2]); break;
        case 4: destImage = DM::NewImage("", destType, dim[0], dim[1], dim[2], dim[3]); break;
        };

        // Prepare (FFT is carried out over work_dim[0])
        size_t srcStride  = DM::ImageGetDataElementByteSize(srcImage);
        size_t destStride = DM::ImageGetDataElementByteSize(destImage);
        size_t srcStrides[MAX_RANK];
        size_t destStrides[MAX_RANK];
        size_t work_dim[MAX_RANK];
        for (int n = 0, m = 1; n < rank; n++) {
                if (n == axis) {
                        srcStrides[0]  = srcStride;
                        destStrides[0] = destStride;
                        work_dim[0]    = dim[n];
                } else {
                        srcStrides[m]  = srcStride;
                        destStrides[m] = destStride;
                        work_dim[m]    = dim[n];
                        ++m;
                }
        	srcStride  *= dim[n];
                destStride *= dim[n];
        }

	// prepare to manipulate the data..
        {
                PlugIn::ImageDataLocker srcLocker(srcImage);
                PlugIn::ImageDataLocker destLocker(destImage);
	        const char *srcData = (const char*)srcLocker.get();
                char *destData      = (char*)destLocker.get();

                if (dataType == ImageData::REAL4_DATA) {
                        do_rfft<float>(::rfftf, destData, srcData, rank, work_dim, destStrides, srcStrides, 1.0);
                        unpack_hermite1<float>(destData, rank, work_dim, destStrides);
                } else {
                        do_rfft<double>(::drfftf, destData, srcData, rank, work_dim, destStrides, srcStrides, 1.0);
                        unpack_hermite1<double>(destData, rank, work_dim, destStrides);
                }

		// tell DigitalMicrograph the image is changed...
		destImage.DataChanged();
	}
	
	PLUG_IN_EXIT

	return destImage.release();	
}

DM_ImageToken_1Ref fftpack_rfft_axis0(DM_ImageToken srcToken)
{
        return fftpack_rfft(srcToken, 0);
}

/**
 * Does a 2D forward FFT from real data (always along X and Y).
 * @param src Source image (must be real)
 * @returns Complex image containing the FFT
 **/
DM_ImageToken_1Ref fftpack_rfft2(DM_ImageToken srcToken)
{
        DM::Image destImage;

        PLUG_IN_ENTRY

        // Test input parameters
        DM::Image srcImage(srcToken);
        int dataType = DM::ImageGetDataType(srcImage);
        int destType;
        switch (dataType) {
        case ImageData::REAL4_DATA:
                destType = ImageData::COMPLEX8_DATA;
                break;
        case ImageData::REAL8_DATA:
                destType = ImageData::COMPLEX16_DATA;
                break;
        default:
                ThrowString("Source image has invalid data type.");
                break;
        }
        int rank = DM::ImageGetNumDimensions(srcImage);
        if (rank < 2 || rank > MAX_RANK)
                ThrowString("Source image has too many or too few dimensions.");
        size_t dim[MAX_RANK];
        for (int n = 0; n < rank; n++)
                dim[n] = DM::ImageGetDimensionSize(srcImage, n);

	// Create destination image
        switch (rank) {
        case 1: destImage = DM::NewImage("", destType, dim[0]); break;
        case 2: destImage = DM::NewImage("", destType, dim[0], dim[1]); break;
        case 3: destImage = DM::NewImage("", destType, dim[0], dim[1], dim[2]); break;
        case 4: destImage = DM::NewImage("", destType, dim[0], dim[1], dim[2], dim[3]); break;
        };

 	// Prepare work array
	size_t srcStride = DM::ImageGetDataElementByteSize(srcImage);
        size_t destStride = DM::ImageGetDataElementByteSize(destImage);
	size_t destStrides[MAX_RANK];
        size_t srcStrides[MAX_RANK];
	size_t work_dim[MAX_RANK];
	for (int n = 0; n < rank; n++) {
		srcStrides[n]  = srcStride;
                destStrides[n] = destStride;
		work_dim[n]     = dim[n];
        	srcStride  *= dim[n];
                destStride *= dim[n];
	}

	// prepare to manipulate the data..
        {
                PlugIn::ImageDataLocker srcLocker(srcImage);
                PlugIn::ImageDataLocker destLocker(destImage);
	        const char *srcData = (const char*)srcLocker.get();
                char *destData      = (char*)destLocker.get();

		// Do RFFT over X
                if (dataType == ImageData::REAL4_DATA)
                        do_rfft<float>(::rfftf, destData, srcData, rank, work_dim, destStrides, srcStrides, 1.0f);
                else
                        do_rfft<double>(::drfftf, destData, srcData, rank, work_dim, destStrides, srcStrides, 1.0);

		// Switch dim 0 and 1, and do only half of the original X-axis
                work_dim[0] = work_dim[0] / 2 + 1;
                swap(srcStrides[0], srcStrides[1]);
                swap(destStrides[0], destStrides[1]);
		swap(work_dim[0], work_dim[1]);
                
		// Do FFT over Y
                if (destType == ImageData::COMPLEX8_DATA)
                        do_fft<float>(::cfftf, destData, destData, rank, work_dim, destStrides, 1.0f);
                else
                        do_fft<double>(::dcfftf, destData, destData, rank, work_dim, destStrides, 1.0);

		// Restore the original axis 0 and 1
                swap(srcStrides[0], srcStrides[1]);
                swap(destStrides[0], destStrides[1]);
		swap(work_dim[0], work_dim[1]);
                work_dim[0] = dim[0];

                if (dataType == ImageData::REAL4_DATA)
                        unpack_hermite2<float>(destData, rank, work_dim, destStrides);
                else
                        unpack_hermite2<double>(destData, rank, work_dim, destStrides);

		// tell DigitalMicrograph the image is changed...
		destImage.DataChanged();
	}
	
	PLUG_IN_EXIT

	return destImage.release();	
}

/**
 * Does a 3D forward FFT from real data (always along X/Y/Z).
 * @param src Source image (must be real)
 * @returns Complex image containing the FFT
 **/
DM_ImageToken_1Ref fftpack_rfft3(DM_ImageToken srcToken)
{
        DM::Image destImage;

        PLUG_IN_ENTRY

        // Test input parameters
        DM::Image srcImage(srcToken);
        int dataType = DM::ImageGetDataType(srcImage);
        int destType;
        switch (dataType) {
        case ImageData::REAL4_DATA:
                destType = ImageData::COMPLEX8_DATA;
                break;
        case ImageData::REAL8_DATA:
                destType = ImageData::COMPLEX16_DATA;
                break;
        default:
                ThrowString("Source image has invalid data type.");
                break;
        }
        int rank = DM::ImageGetNumDimensions(srcImage);
        if (rank < 2 || rank > MAX_RANK)
                ThrowString("Source image has too many or too few dimensions.");
        size_t dim[MAX_RANK];
        for (int n = 0; n < rank; n++)
                dim[n] = DM::ImageGetDimensionSize(srcImage, n);

	// Create destination image
        switch (rank) {
        case 1: destImage = DM::NewImage("", destType, dim[0]); break;
        case 2: destImage = DM::NewImage("", destType, dim[0], dim[1]); break;
        case 3: destImage = DM::NewImage("", destType, dim[0], dim[1], dim[2]); break;
        case 4: destImage = DM::NewImage("", destType, dim[0], dim[1], dim[2], dim[3]); break;
        };

 	// Prepare work array
	size_t srcStride = DM::ImageGetDataElementByteSize(srcImage);
        size_t destStride = DM::ImageGetDataElementByteSize(destImage);
	size_t destStrides[MAX_RANK];
        size_t srcStrides[MAX_RANK];
	size_t work_dim[MAX_RANK];
	for (int n = 0; n < rank; n++) {
		srcStrides[n]  = srcStride;
                destStrides[n] = destStride;
		work_dim[n]     = dim[n];
        	srcStride  *= dim[n];
                destStride *= dim[n];
	}

	// prepare to manipulate the data..
        {
                PlugIn::ImageDataLocker srcLocker(srcImage);
                PlugIn::ImageDataLocker destLocker(destImage);
	        const char *srcData = (const char*)srcLocker.get();
                char *destData      = (char*)destLocker.get();

		// Do RFFT over X
                if (dataType == ImageData::REAL4_DATA)
                        do_rfft<float>(::rfftf, destData, srcData, rank, work_dim, destStrides, srcStrides, 1.0f);
                else
                        do_rfft<double>(::drfftf, destData, srcData, rank, work_dim, destStrides, srcStrides, 1.0);

		// Switch dim 0 and 1, and do only half of the original X-axis
                work_dim[0] = work_dim[0] / 2 + 1;
                swap(srcStrides[0], srcStrides[1]);
                swap(destStrides[0], destStrides[1]);
		swap(work_dim[0], work_dim[1]);
                
		// Do FFT over Y
                if (destType == ImageData::COMPLEX8_DATA)
                        do_fft<float>(::cfftf, destData, destData, rank, work_dim, destStrides, 1.0f);
                else
                        do_fft<double>(::dcfftf, destData, destData, rank, work_dim, destStrides, 1.0);

		// Switch dim 0 and 2
                swap(srcStrides[0], srcStrides[2]);
                swap(destStrides[0], destStrides[2]);
		swap(work_dim[0], work_dim[2]);

		// Do FFT over Z
                if (destType == ImageData::COMPLEX8_DATA)
			do_fft<float>(::cfftf, destData, destData, rank, work_dim, destStrides, 1.0f);
                else
                        do_fft<double>(::dcfftf, destData, destData, rank, work_dim, destStrides, 1.0);

		// Restore the original axes
                swap(srcStrides[0], srcStrides[2]);
                swap(destStrides[0], destStrides[2]);
                swap(srcStrides[0], srcStrides[1]);
                swap(destStrides[0], destStrides[1]);
                swap(work_dim[0], work_dim[2]);
		swap(work_dim[0], work_dim[1]);
                work_dim[0] = dim[0];

                if (dataType == ImageData::REAL4_DATA)
                        unpack_hermite3<float>(destData, rank, work_dim, destStrides);
                else
                        unpack_hermite3<double>(destData, rank, work_dim, destStrides);

		// tell DigitalMicrograph the image is changed...
		destImage.DataChanged();
	}
	
	PLUG_IN_EXIT

	return destImage.release();	
}


/**
 * Does a backward RFFT along the axis @a axis (0=X, 1=Y, ...).
 * 
 * @param src Source image (must be floating point type)
 * @param axis Axis along the which the FFT is done (0=X, 1=Y, ...).
 * @returns Complex image containing the FFT
 **/
DM_ImageToken_1Ref fftpack_rifft(DM_ImageToken srcToken, long axis)
{
        DM::Image destImage;

        PLUG_IN_ENTRY

        // Test input parameters
        DM::Image srcImage(srcToken);
        int dataType = DM::ImageGetDataType(srcImage);
        int destType;
        switch (dataType) {
        case ImageData::COMPLEX8_DATA:
                destType = ImageData::REAL4_DATA;
                break;
        case ImageData::COMPLEX16_DATA:
                destType = ImageData::REAL8_DATA;
                break;
        default:
                ThrowString("Source image has invalid data type.");
                break;
        }
        int rank = DM::ImageGetNumDimensions(srcImage);
        if (rank < 1 || rank > MAX_RANK)
                ThrowString("Source image has too many or too few dimensions.");
        if (axis < 0)
                axis += rank;
        if (axis < 0 || axis >= rank)
                ThrowString("Invalid axes selected.");
        size_t dim[MAX_RANK];
        for (int n = 0; n < rank; n++)
                dim[n] = DM::ImageGetDimensionSize(srcImage, n);

	// Create destination image
        switch (rank) {
        case 1: destImage = DM::NewImage("", destType, dim[0]); break;
        case 2: destImage = DM::NewImage("", destType, dim[0], dim[1]); break;
        case 3: destImage = DM::NewImage("", destType, dim[0], dim[1], dim[2]); break;
        case 4: destImage = DM::NewImage("", destType, dim[0], dim[1], dim[2], dim[3]); break;
        };

        // Prepare (FFT is carried out over work_dim[0])
        size_t srcStride  = DM::ImageGetDataElementByteSize(srcImage);
        size_t destStride = DM::ImageGetDataElementByteSize(destImage);
        size_t srcStrides[MAX_RANK];
        size_t destStrides[MAX_RANK];
        size_t work_dim[MAX_RANK];
        for (int n = 0, m = 1; n < rank; n++) {
                if (n == axis) {
                        srcStrides[0]  = srcStride;
                        destStrides[0] = destStride;
                        work_dim[0]    = dim[n];
                } else {
                        srcStrides[m]  = srcStride;
                        destStrides[m] = destStride;
                        work_dim[m]    = dim[n];
                        ++m;
                }
        	srcStride  *= dim[n];
                destStride *= dim[n];
        }

	// prepare to manipulate the data..
        {
                PlugIn::ImageDataLocker srcLocker(srcImage);
                PlugIn::ImageDataLocker destLocker(destImage);
	        const char *srcData = (const char*)srcLocker.get();
                char *destData      = (char*)destLocker.get();

                if (destType == ImageData::REAL4_DATA)
                        do_rifft<float>(::rfftb, destData, srcData, rank, work_dim, destStrides, srcStrides, 1.0f / work_dim[0]);
                else
                        do_rifft<double>(::drfftb, destData, srcData, rank, work_dim, destStrides, srcStrides, 1.0 / work_dim[0]);

		// tell DigitalMicrograph the image is changed...
		destImage.DataChanged();
	}
	
	PLUG_IN_EXIT

	return destImage.release();	
}

DM_ImageToken_1Ref fftpack_rifft_axis0(DM_ImageToken srcToken)
{
        return fftpack_rifft(srcToken, 0);
}

void fftpackPlugin::Start()
{
        g_FactorsCache = new FactorsCache(1 << 20);

        AddFunction("dm_string fftpack_version()", &fftpack_version);

        AddFunction("ImageRef fftpack_fft(Image*, long)", &fftpack_fft);
        AddFunction("ImageRef fftpack_fft(Image*)", &fftpack_fft_axis0);
        AddFunction("ImageRef fftpack_ifft(Image*, long)", &fftpack_ifft);
        AddFunction("ImageRef fftpack_ifft(Image*)", &fftpack_ifft_axis0);
	AddFunction("ImageRef fftpack_fft2(Image*)", &fftpack_fft2);
        AddFunction("ImageRef fftpack_ifft2(Image*)", &fftpack_ifft2);
	AddFunction("ImageRef fftpack_fft3(Image*)", &fftpack_fft3);
        AddFunction("ImageRef fftpack_ifft3(Image*)", &fftpack_ifft3);

        AddFunction("ImageRef fftpack_rfft(Image*, long)", &fftpack_rfft);
        AddFunction("ImageRef fftpack_rfft(Image*)", &fftpack_rfft_axis0);
	AddFunction("ImageRef fftpack_rfft2(Image*)", &fftpack_rfft2);
        AddFunction("ImageRef fftpack_rfft3(Image*)", &fftpack_rfft3);
        AddFunction("ImageRef fftpack_rifft(Image*, long)", &fftpack_rifft);
        AddFunction("ImageRef fftpack_rifft(Image*)", &fftpack_rifft_axis0);
}

void fftpackPlugin::Run()
{
}

void fftpackPlugin::Cleanup()
{
        g_FactorsCache->clear();
}

void fftpackPlugin::End()
{
        delete g_FactorsCache;
        g_FactorsCache = 0;
}

fftpackPlugin g_fftpackPlugin;
FactorsCache* g_FactorsCache = 0;
