.. fftpack documentation master file, created by
   sphinx-quickstart on Thu Jul 05 15:19:02 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. default-domain:: cpp

Welcome to fftpack's documentation!
===================================

The fftpack plugin for DigitalMicrograph is merely a wrapper around the
Fast Fourier Transform (FFT) routines of the numpy package (www.numpy.org).
Those are based on algorithms of `FFTPACK <http://www.netlib.org/fftpack/>`_ 
(originally written by Paul N. Schwartztrauber from the National Center for Atmospheric
Research, Boulder, CO).

The FFTPACK contains routines do fourier transforms of arbitrary shaped 
one dimensional arrays. The plugin also contains routines to do multidimensional
FFTs. 

Different from DigitalMicrograph the FFTPACK routines uses the usual
FFT conventions (see below). DigitalMicrograph uses the conventions of FFTPACK's
forward FFT as inverse FFT and FFTPACK's backward FFT as forward FFT. Also
the routines in the plugin don't move the origin of the fourier spectrum. This 
means the origin of the fourier transform is found at the pixel indices (0,0).

Different from the original FFTPACK routines the plugin normalizes the inverse
transformation. So calling :cpp:func:`fftpack_fft()` and immediately 
:cpp:func:`fftpack_ifft()` on the result will return the original array.

What the forward FFT calculates is (input and output arrays have *N* 
elements, starting at index 0, *i* is sqrt(-1))::

    output[j] = sum(0 <= k < N) of
                input[k]*exp(-i*j*k*2*pi/N)

The backward fft calculates::

    output[j] = 1/N sum(0 <= k < N) of
                input[k]*exp(+i*j*k*2*pi/N)
                
The multidimensional FFTs just subsequently perform a FFT over the different
dimensions.
                
The routines respect the precision of the argument. Thus if the input
has double precision (i.e. Complex16) the FFT is done with double precision
and the result is also double precision. If input is only of single precision 
(i.e. Complex8) the whole operation is done with single precision only.
                
Functions of the plugin
=======================

.. cpp:function:: string fftpack_version()

    Returns version string for the plugin.

.. cpp:function:: image fftpack_fft(image input, number axis = 0)

    Calculates one-dimensional forward FFT of *input*. The array *input*
    must be of complex type. 
    
    *axis* indicates over which axis the FFT is calculated. A value of 0 
    (default) means the X-axis, 1 the Y-Axis, etc. Example: on a 2D input 
    array, with an *axis* value of 0, every row is fourier transformed, for an value of 1
    every column is fourier transformed.

.. cpp:function:: image fftpack_ifft(image input, number axis = 0)

    Calculates one-dimensional backward FFT of *input*. The array *input*
    must be of complex type.
    
    *axis* indicates over which axis the FFT is calculated. A value of 0 
    (default) means the X-axis, 1 the Y-Axis, etc. Example: on a 2D input 
    array, with an *axis* value of 0, every row is fourier transformed, for an value of 1
    every column is fourier transformed.

.. cpp:function:: image fftpack_fft2(image input)

    Calculates the two-dimensional forward FFT of *input*. The FFT is calculated
    over the X and Y axes. For a 3D array this means every "slice" is transformed
    individually.
    
.. cpp:function:: image fftpack_ifft2(image input)

    Calculates the two-dimensional backward FFT of *input*. The FFT is calculated
    over the X and Y axes. For a 3D array this means every "slice" is transformed
    individually.

.. function:: image fftpack_fft3(image input)

    Calculates the three-dimensional forward FFT of *input*. The FFT is calculated
    over the X, Y, and Z axes. 
    
.. function:: image fftpack_ifft3(image input)

    Calculates the three-dimensional backward FFT of *input*. The FFT is calculated
    over the X, Y, and Z axes. 

.. cpp:function:: image fftpack_rfft(image input, number axis = 0)

    Calculates one-dimensional forward FFT of the real data in *input*. Returns
    complex image (same layout as result of :cpp:func:`fftpack_fft`). The array *input*
    must be of real type. 
    
    *axis* indicates over which axis the FFT is calculated. A value of 0 
    (default) means the X-axis, 1 the Y-Axis, etc. Example: on a 2D input 
    array, with an *axis* value of 0, every row is fourier transformed, for an value of 1
    every column is fourier transformed.

.. cpp:function:: image fftpack_rifft(image input, number axis = 0)

    Calculates one-dimensional backward FFT of the complex data in *input*. 
    The returned image is real. This function assumes, that the data in *input*
    represents the Fourier transform of real data (thus is hermite symmetric).
    Effectively only the lower half of the *input* data is used.
    The array *input* must be of complex type.
    
    *axis* indicates over which axis the FFT is calculated. A value of 0 
    (default) means the X-axis, 1 the Y-Axis, etc. Example: on a 2D input 
    array, with an *axis* value of 0, every row is fourier transformed, for an value of 1
    every column is fourier transformed.

.. cpp:function:: image fftpack_rfft2(image input)

    Calculates the two-dimensional forward FFT of real data in *input*. 
    The FFT is calculated over the X and Y axes. For a 3D array this means 
    every "slice" is transformed individually. *input* is expected to
    be of real type, the result is of complex type. 

.. function:: image fftpack_rfft3(image input)

    Calculates the three-dimensional forward FFT of real data in *input*. 
    The FFT is calculated over the X, Y, and Z axes. *input* is expected to
    be of real type, the result is of complex type.  

License
=======

| Copyright (c) 2005-2012, Tore Niermann and Numpy Developers
| Contact: niermann (at) physik.tu-berlin.de
| All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the author nor the names of any contributors may 
      be used to endorse or promote products derived from this software 
      without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF 
THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
