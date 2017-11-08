#ifndef FFTPACK_PLUGIN_INC
#define FFTPACK_PLUGIN_INC

#define _GATANPLUGIN_USES_LIBRARY_VERSION 2
#define _GATAN_USE_STL_STRING
#include "DMPlugInBasic.h"

#include <string.h>
#include <boost/cstdint.hpp>
#include "cache.h"

// Plugin version string
#define FFTPACK_PLUGIN_VERSION    "1.1"

// Max supported image rank
#define MAX_RANK 4

//----------------------------------------------------------------------------------------
// Plugin (plugin.cpp)

class fftpackPlugin : public GatanPlugIn::PlugInMain
{
private:
        virtual void Start();
	virtual void Run();
	virtual void Cleanup();
	virtual void End();
};

//----------------------------------------------------------------------------------------
// Globals

extern fftpackPlugin g_fftpackPlugin;
extern FactorsCache* g_FactorsCache;

#endif // FFTPACK_PLUGIN_INC
