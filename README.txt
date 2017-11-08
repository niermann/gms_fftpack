1. Introduction
======================================================================

The plugin contains FFT functions for Digital Micrograph:

    * with non-power of two array shape
    * with single and double precision
    
The functions are based on the fftpack.c implementation provided with 
the numpy package (see www.numpy.org). 

2. Building
======================================================================

    For the DMSDKs contact Gatan or look on their website www.gatan.com.
    
    It is assumed, that the DMSDKs and their dependencies (e.g. boost) 
    were already properly installed.

    GMS-2.X:
    * Open vc2008/fftpack.vcproj with Visual Studio 2008
    * Add your DMSDK include and library directories to the VC++ Directories 
      in Visual Studio (in Tools/Options/Projects and Solutions/VC++ Directories). The
      DMSDK is usally installed to C:\ProgramData\Gatan\DMSDK
    * Build the Release build
    * You should find your plugin as vc2008/Release/fftpack_GMS2X_x86.dll

    GMS-1.X:
    * Open vc2003/fftpack.vcproj with Visual Studio 2003
    * Add your DMSDK include and library directories to the VC++ Directories 
      in Visual Studio (in Tools/Options/Projects/VC++ Directories). 
    * Build the Release build
    * You should find your plugin as vc2003/Release/fftpack_GMS1X.dll

3. Install
======================================================================

Simply copy the fftpack.dll to the plugin directory of your Digital
Micrograph installation. By default these are for GMS-2.X:
    c:\ProgramData\Gatan\Plugins
and for GMS-1.X:
    c:\Program Files\Gatan\DigitalMicrograph\Plugins

4. Documentation
======================================================================

The documentation uses the Sphinx (sphinx.pocoo.org) documentation system,
which in turn requires a working Python (www.python.org) installation. 
To build the documentation simply go into the documentation subdirectory
(doc/) and invoke the sphinx system (e.g. by calling make.bat html)
