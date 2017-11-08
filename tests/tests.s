Image make_random_complex1D(number byteSize, number dim0)
{
	Image result := ComplexImage("", byteSize, dim0)
	result = complex(random(), random());
	return result
}

Image make_random_complex2D(number byteSize, number dim0, number dim1)
{
	Image result := ComplexImage("", byteSize, dim0, dim1)
	result = complex(random(), random());
	return result
}

Image make_random_complex3D(number byteSize, number dim0, number dim1, number dim2)
{
	Image result := ComplexImage("", byteSize, dim0, dim1, dim2)
	result = complex(random(), random());
	return result
}

Image make_random_real1D(number byteSize, number dim0)
{
	Image result := RealImage("", byteSize, dim0)
	result = random()
	return result
}

Image make_random_real2D(number byteSize, number dim0, number dim1)
{
	Image result := RealImage("", byteSize, dim0, dim1)
	result = random()
	return result
}

Image make_random_real3D(number byteSize, number dim0, number dim1, number dim2)
{
	Image result := RealImage("", byteSize, dim0, dim1, dim2)
	result = random()
	return result
}

void test_fft_InvalidType()
{
	Result("test_fft_InvalidType: ")
	
	Image test := RealImage("", 4, 128, 128)
	try {
		fftpack_fft(test)
	} catch {
		Result("PASSED\n")
		return
		break
	} 

	Result("FAILED\n")
}

void test_fft_InvalidAxis()
{
	Result("test_fft_InvalidAxis: ")
	
	Image test := ComplexImage("", 8, 128, 128)
	try {
		fftpack_fft(test, 2)
	} catch {
		Result("PASSED\n")
		return
		break
	} 

	Result("FAILED\n")
}

void test_fft_single_axis0()
{
	Result("test_fft_single_axis0: ")
	
	ComplexImage testLine := make_random_complex1D(8, 128)
	ComplexImage test := ComplexImage("", 8, 128, 128)
	test = testLine[icol,0]
	
	// DM Has forward and backward FFTs swapped
	ShiftCenter(testLine)
	testLine := IFFT(testLine) * 128.0
	test := fftpack_fft(test)

	ComplexImage cmp := ComplexImage("", 8, 128, 128)
	cmp = testLine[icol,0]
	
	RealImage diff = abs(cmp - test)
	if (mean(diff) < 1e-6)
		Result("PASSED\n")
	else
		Result("FAILED: " + mean(diff) + "\n")
}

void test_fft_single_axis1()
{
	Result("test_fft_single_axis1: ")
	
	ComplexImage testLine := make_random_complex1D(8, 128)
	ComplexImage test := ComplexImage("", 8, 128, 128)
	test = testLine[irow, 0]
	
	// DM Has forward and backward FFTs swapped
	ShiftCenter(testLine)
	testLine := IFFT(testLine) * 128.0
	test := fftpack_fft(test, 1)

	ComplexImage cmp := ComplexImage("", 8, 128, 128)
	cmp = testLine[irow, 0]
	
	RealImage diff = abs(cmp - test)
	if (mean(diff) < 1e-6)
		Result("PASSED\n")
	else
		Result("FAILED: " + mean(diff) + "\n")
}

void test_ifft_InvalidType()
{
	Result("test_ifft_InvalidType: ")

	Image test := RealImage("", 4, 128, 128)
	try {
		fftpack_ifft(test)
	} catch {
		Result("PASSED\n")
		return
		break
	} 
	
	Result("FAILED\n")
}

void test_ifft_InvalidAxis()
{
	Result("test_ifft_InvalidAxis: ")

	Image test := ComplexImage("", 8, 128, 128)
	try {
		fftpack_ifft(test, -3)
	} catch {
		Result("PASSED\n")
		return
		break
	} 
	
	Result("FAILED\n")
}

void test_ifft_single_axis0()
{
	Result("test_ifft_single_axis0: ")
	
	ComplexImage testLine := make_random_complex1D(8, 128)
	ComplexImage test := ComplexImage("", 8, 128, 128)
	test = testLine[icol,0]
	
	// DM Has forward and backward FFTs swapped
	testLine := FFT(testLine) / 128.0
	ShiftCenter(testLine)
	test := fftpack_ifft(test)

	ComplexImage cmp := ComplexImage("", 8, 128, 128)
	cmp = testLine[icol,0]
	
	
	RealImage diff = abs(cmp - test)
	if (mean(diff) < 1e-6)
		Result("PASSED\n")
	else
		Result("FAILED: " + mean(diff) + "\n")
}

void test_ifft_single_axis1()
{
	Result("test_ifft_single_axis1: ")
	
	ComplexImage testLine := make_random_complex1D(8, 128)
	ComplexImage test := ComplexImage("", 8, 128, 128)
	test = testLine[irow, 0]
	
	// DM Has forward and backward FFTs swapped
	testLine := FFT(testLine) / 128.0
	ShiftCenter(testLine)
	test := fftpack_ifft(test, 1)

	ComplexImage cmp := ComplexImage("", 8, 128, 128)
	cmp = testLine[irow, 0]
	
	RealImage diff = abs(cmp - test)
	if (mean(diff) < 1e-6)
		Result("PASSED\n")
	else
		Result("FAILED\n")
}

void test_single_fft_forward_backward_even()
{
	Result("test_single_fft_forward_backward_even: ")
	
	ComplexImage test := make_random_complex2D(8, 128, 128)
	ComplexImage cmp := fftpack_fft(test, 0)
	cmp := fftpack_ifft(cmp, 0)

	RealImage diff = abs(cmp - test)
	if (mean(diff) < 1e-6)
		Result("PASSED\n")
	else
		Result("FAILED: " + mean(diff) + "\n")
}

void test_single_fft_forward_backward_odd()
{
	Result("test_single_fft_forward_backward_odd: ")
	
	ComplexImage test := make_random_complex2D(8, 2*3*4*5, 2*3*4*5)
	ComplexImage cmp := fftpack_fft(test, 0)
	cmp := fftpack_ifft(cmp, 0)

	RealImage diff = abs(cmp - test)
	if (mean(diff) < 1e-6)
		Result("PASSED\n")
	else
		Result("FAILED: " + mean(diff) + "\n")
}

void test_double_fft_forward_backward_even()
{
	Result("test_double_fft_forward_backward_even: ")
	
	ComplexImage test := make_random_complex2D(16, 128, 128)
	ComplexImage cmp := fftpack_fft(test, 0)
	if (cmp.ImageGetDataElementByteSize() != 16) {
		Result("FAILED: Not double.\n")
		return
	}
	cmp := fftpack_ifft(cmp, 0)
	if (cmp.ImageGetDataElementByteSize() != 16) {
		Result("FAILED: Not double.\n")
		return
	}

	RealImage diff = abs(cmp - test)
	if (mean(diff) < 1e-12)
		Result("PASSED\n")
	else
		Result("FAILED: " + mean(diff) + "\n")
}

void test_double_fft_forward_backward_odd()
{
	Result("test_double_fft_forward_backward_odd: ")
	
	ComplexImage test := make_random_complex2D(16, 2*3*4*5, 2*3*4*5)
	ComplexImage cmp := fftpack_fft(test, 0)
	if (cmp.ImageGetDataElementByteSize() != 16) {
		Result("FAILED: Not double.\n")
		return
	}
	cmp := fftpack_ifft(cmp, 0)
	if (cmp.ImageGetDataElementByteSize() != 16) {
		Result("FAILED: Not double.\n")
		return
	}

	RealImage diff = abs(cmp - test)
	if (mean(diff) < 1e-12)
		Result("PASSED\n")
	else
		Result("FAILED: " + mean(diff) + "\n")
}

void test_fft2_InvalidRank()
{
	Result("test_fft2_InvalidRank: ")
	
	Image test := ComplexImage("", 8, 128)
	try {
		fftpack_fft2(test)
	} catch {
		Result("PASSED\n")
		return
		break
	} 

	Result("FAILED\n")
}

void test_fft2_single()
{
	Result("test_fft2_single: ")
	
	ComplexImage test := make_random_complex2D(8, 128, 128)
	ComplexImage cmp = test
	
	// DM Has forward and backward FFTs swapped
	ShiftCenter(cmp)
	cmp := IFFT(cmp) * 128.0 * 128.0
	test := fftpack_fft2(test)

	RealImage diff = abs(cmp - test)
	if (mean(diff) < 5e-5)
		Result("PASSED\n")
	else
		Result("FAILED: " + mean(diff) + "\n")
}

void test_ifft2_single()
{
	Result("test_ifft2_single: ")
	
	ComplexImage test := make_random_complex2D(8, 128, 128)
	ComplexImage cmp = test
	
	// DM Has forward and backward FFTs swapped
	cmp := FFT(cmp) / 128.0 / 128.0
	ShiftCenter(cmp)
	test := fftpack_ifft2(test)

	RealImage diff = abs(cmp - test)
	if (mean(diff) < 5e-5)
		Result("PASSED\n")
	else
		Result("FAILED: " + mean(diff) + "\n")
}

void test_double_fft2_forward_backward_even()
{
	Result("test_double_fft2_forward_backward_even: ")
	
	ComplexImage test := make_random_complex2D(16, 128, 256)
	ComplexImage cmp := fftpack_fft2(test)
	if (cmp.ImageGetDataElementByteSize() != 16) {
		Result("FAILED: Not double.\n")
		return
	}
	cmp := fftpack_ifft2(cmp)
	if (cmp.ImageGetDataElementByteSize() != 16) {
		Result("FAILED: Not double.\n")
		return
	}

	RealImage diff = abs(cmp - test)
	if (mean(diff) < 1e-12)
		Result("PASSED\n")
	else
		Result("FAILED: " + mean(diff) + "\n")
}

void test_double_fft2_forward_backward_odd()
{
	Result("test_double_fft2_forward_backward_odd: ")
	
	ComplexImage test := make_random_complex2D(16, 192, 2*3*4*5)
	ComplexImage cmp := fftpack_fft2(test)
	if (cmp.ImageGetDataElementByteSize() != 16) {
		Result("FAILED: Not double.\n")
		return
	}
	cmp := fftpack_ifft2(cmp)
	if (cmp.ImageGetDataElementByteSize() != 16) {
		Result("FAILED: Not double.\n")
		return
	}

	RealImage diff = abs(cmp - test)
	if (mean(diff) < 1e-12)
		Result("PASSED\n")
	else
		Result("FAILED: " + mean(diff) + "\n")
}

void test_single_fft2_forward_backward_even()
{
	Result("test_single_fft2_forward_backward_even: ")
	
	ComplexImage test := make_random_complex2D(16, 128, 256)
	ComplexImage cmp := fftpack_fft2(test)
	cmp := fftpack_ifft2(cmp)

	RealImage diff = abs(cmp - test)
	if (mean(diff) < 1e-6)
		Result("PASSED\n")
	else
		Result("FAILED: " + mean(diff) + "\n")
}

void test_single_fft2_forward_backward_odd()
{
	Result("test_single_fft2_forward_backward_odd: ")
	
	ComplexImage test := make_random_complex2D(16, 192, 2*3*4*5)
	ComplexImage cmp := fftpack_fft2(test)
	cmp := fftpack_ifft2(cmp)

	RealImage diff = abs(cmp - test)
	if (mean(diff) < 1e-6)
		Result("PASSED\n")
	else
		Result("FAILED: " + mean(diff) + "\n")
}

void test_fft3_InvalidRank()
{
	Result("test_fft3_InvalidRank: ")
	
	Image test := ComplexImage("", 8, 128, 128)
	try {
		fftpack_fft3(test)
	} catch {
		Result("PASSED\n")
		return
		break
	} 

	Result("FAILED\n")
}

void test_single_fft3_forward_backward_even()
{
	Result("test_single_fft3_forward_backward_even: ")
	
	ComplexImage test := make_random_complex3D(16, 128, 256, 64)
	ComplexImage cmp := fftpack_fft3(test)
	cmp := fftpack_ifft3(cmp)

	RealImage diff = abs(cmp - test)
	if (mean(diff) < 1e-6)
		Result("PASSED\n")
	else
		Result("FAILED: " + mean(diff) + "\n")
}

void test_single_fft3_forward_backward_odd()
{
	Result("test_single_fft3_forward_backward_odd: ")
	
	ComplexImage test := make_random_complex3D(16, 192, 2*3*4*5, 60)
	ComplexImage cmp := fftpack_fft3(test)
	cmp := fftpack_ifft3(cmp)

	RealImage diff = abs(cmp - test)
	if (mean(diff) < 1e-6)
		Result("PASSED\n")
	else
		Result("FAILED: " + mean(diff) + "\n")
}

void test_rfft_InvalidType()
{
	Result("test_rfft_InvalidType: ")
	
	Image test := ComplexImage("", 8, 128, 128)
	try {
		fftpack_rfft(test)
	} catch {
		Result("PASSED\n")
		return
		break
	} 

	Result("FAILED\n")
}

void test_rfft_InvalidAxis()
{
	Result("test_rfft_InvalidAxis: ")
	
	Image test := RealImage("", 4, 128, 128)
	try {
		fftpack_rfft(test, 2)
	} catch {
		Result("PASSED\n")
		return
		break
	} 

	Result("FAILED\n")
}

void test_rifft_InvalidType()
{
	Result("test_rifft_InvalidType: ")

	Image test := RealImage("", 4, 128, 128)
	try {
		fftpack_rifft(test)
	} catch {
		Result("PASSED\n")
		return
		break
	} 
	
	Result("FAILED\n")
}

void test_rifft_InvalidAxis()
{
	Result("test_rifft_InvalidAxis: ")

	Image test := ComplexImage("", 8, 128, 128)
	try {
		fftpack_rifft(test, -3)
	} catch {
		Result("PASSED\n")
		return
		break
	} 
	
	Result("FAILED\n")
}

void test_rfft_single_axis0()
{
	Result("test_rfft_single_axis0: ")
	
	Image testLine := make_random_real1D(4, 128)
	Image test := RealImage("", 4, 128, 128)
	test = testLine[icol,0]
	
	test := fftpack_rfft(test)

	ComplexImage cmp := ComplexImage("", 8, 128, 128)
	cmp = testLine[icol,0]
	cmp := fftpack_fft(cmp)
	
	RealImage diff = abs(cmp - test)
	if (mean(diff) < 1e-6)
		Result("PASSED\n")
	else
		Result("FAILED: " + mean(diff) + "\n")
}

void test_rfft_single_axis1()
{
	Result("test_rfft_single_axis1: ")
	
	Image testLine := make_random_real1D(4, 128)
	Image test := RealImage("", 4, 128, 128)
	test = testLine[irow,0]
	
	test := fftpack_rfft(test, 1)

	ComplexImage cmp := ComplexImage("", 8, 128, 128)
	cmp = testLine[irow,0]
	cmp := fftpack_fft(cmp, 1)
	
	RealImage diff = abs(cmp - test)
	if (mean(diff) < 1e-6)
		Result("PASSED\n")
	else
		Result("FAILED: " + mean(diff) + "\n")
}

void test_single_rfft_forward_backward_even()
{
	Result("test_single_rfft_forward_backward_even: ")
	
	Image test := make_random_real2D(4, 128, 128)
	ComplexImage tmp := fftpack_rfft(test, 0)
	Image cmp := fftpack_rifft(tmp, 0)

	RealImage diff = abs(cmp - test)
	if (mean(diff) < 1e-6)
		Result("PASSED\n")
	else
		Result("FAILED: " + mean(diff) + "\n")
}

void test_single_rfft_forward_backward_odd()
{
	Result("test_single_rfft_forward_backward_odd: ")
	
	Image test := make_random_real2D(8, 125, 125)
	ComplexImage tmp := fftpack_rfft(test, 0)
	Image cmp := fftpack_rifft(tmp, 0)

	RealImage diff = abs(cmp - test)
	if (mean(diff) < 1e-6)
		Result("PASSED\n")
	else
		Result("FAILED: " + mean(diff) + "\n")
}

void test_double_rfft_forward_backward_even()
{
	Result("test_double_rfft_forward_backward_even: ")
	
	Image test := make_random_real2D(8, 128, 128)
	ComplexImage tmp := fftpack_rfft(test, 0)
	if (tmp.ImageGetDataElementByteSize() != 16) {
		Result("FAILED: Not double.\n")
		return
	}
	Image cmp := fftpack_rifft(tmp, 0)
	if (cmp.ImageGetDataElementByteSize() != 8) {
		Result("FAILED: Not double.\n")
		return
	}

	RealImage diff = abs(cmp - test)
	if (mean(diff) < 1e-12)
		Result("PASSED\n")
	else
		Result("FAILED: " + mean(diff) + "\n")
}

void test_double_rfft_forward_backward_odd()
{
	Result("test_double_rfft_forward_backward_odd: ")
	
	Image test := make_random_real2D(8, 125, 125)
	ComplexImage tmp := fftpack_rfft(test, 0)
	if (tmp.ImageGetDataElementByteSize() != 16) {
		Result("FAILED: Not double.\n")
		return
	}
	Image cmp := fftpack_rifft(tmp, 0)
	if (cmp.ImageGetDataElementByteSize() != 8) {
		Result("FAILED: Not double.\n")
		return
	}

	RealImage diff = abs(cmp - test)
	if (mean(diff) < 1e-12)
		Result("PASSED\n")
	else
		Result("FAILED: " + mean(diff) + "\n")
}

void test_rfft2_InvalidRank()
{
	Result("test_rfft2_InvalidRank: ")
	
	Image test := RealImage("", 8, 128)
	try {
		fftpack_rfft2(test)
	} catch {
		Result("PASSED\n")
		return
		break
	} 

	Result("FAILED\n")
}

void test_rfft2_single_even()
{
	Result("test_rfft2_single_even: ")
	
	Image test := make_random_real2D(4, 128, 128)

	ComplexImage cmp := ComplexImage("", 8, 128, 128)
	cmp[0,0,128,128] = test[0,0,128,128]
	cmp := fftpack_fft2(cmp)
	
	// DM Has forward and backward FFTs swapped
	test := fftpack_rfft2(test)

	RealImage diff = abs(cmp - test)
	if (mean(diff) < 5e-5)
		Result("PASSED\n")
	else
		Result("FAILED: " + mean(diff) + "\n")
}

void test_rfft2_single_odd()
{
	Result("test_rfft2_single_odd: ")
	
	Image test := make_random_real2D(4, 125, 125)

	ComplexImage cmp := ComplexImage("", 8, 125, 125)
	cmp = test
	cmp := fftpack_fft2(cmp)
	
	// DM Has forward and backward FFTs swapped
	test := fftpack_rfft2(test)

	RealImage diff = abs(cmp - test)
	if (mean(diff) < 5e-5)
		Result("PASSED\n")
	else
		Result("FAILED: " + mean(diff) + "\n")
}

void test_rfft3_single_even()
{
	Result("test_rfft3_single_even: ")
	
	Image test := make_random_real3D(4, 128, 128, 128)

	ComplexImage cmp := ComplexImage("", 8, 128, 128, 128)
	cmp  = test
	cmp := fftpack_fft3(cmp)
	
	// DM Has forward and backward FFTs swapped
	test := fftpack_rfft3(test)

	RealImage diff = abs(cmp - test)
	if (mean(diff) < 2e-4)
		Result("PASSED\n")
	else
		Result("FAILED: " + mean(diff) + "\n")
}

void test_rfft3_single_odd()
{
	Result("test_rfft3_single_odd: ")
	
	Image test := make_random_real3D(4, 125, 125, 125)

	ComplexImage cmp := ComplexImage("", 8, 125, 125, 125)
	cmp = test
	cmp := fftpack_fft3(cmp)
	
	// DM Has forward and backward FFTs swapped
	test := fftpack_rfft3(test)

	RealImage diff = abs(cmp - test)
	if (mean(diff) < 2e-4)
		Result("PASSED\n")
	else
		Result("FAILED: " + mean(diff) + "\n")
}

Result("Testing fftpack version " + fftpack_version() + "\n");

test_fft_InvalidType()
test_fft_InvalidAxis()
test_fft_single_axis0()
test_fft_single_axis1()
test_ifft_InvalidType()
test_ifft_InvalidAxis()
test_ifft_single_axis0()
test_ifft_single_axis1()
test_single_fft_forward_backward_even()
test_single_fft_forward_backward_odd()
test_double_fft_forward_backward_even()
test_double_fft_forward_backward_odd()
test_fft2_InvalidRank()
test_fft2_single()
test_ifft2_single()
test_double_fft2_forward_backward_even()
test_double_fft2_forward_backward_odd()
test_single_fft2_forward_backward_even()
test_single_fft2_forward_backward_odd()
test_fft3_InvalidRank()
test_single_fft3_forward_backward_even()
test_single_fft3_forward_backward_odd()

test_rfft_InvalidType()
test_rfft_InvalidAxis()
test_rifft_InvalidType()
test_rifft_InvalidAxis()
test_rfft_single_axis0()
test_rfft_single_axis1()
test_single_rfft_forward_backward_even()
test_single_rfft_forward_backward_odd()
test_double_rfft_forward_backward_even()
test_double_rfft_forward_backward_odd()
test_rfft2_InvalidRank()
test_rfft2_single_even()
test_rfft2_single_odd()
test_rfft3_single_even()
test_rfft3_single_odd()