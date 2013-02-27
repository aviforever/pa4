#include "ImageCleaner.h"
#include <math.h>
#include <omp.h>
#include <sys/time.h>
#include <stdio.h>

#define PI	3.14159265

void cpu_fftx(float *real_image, float *imag_image, int size_x, int size_y)
{
  // Create some space for storing temporary values
  float *realOutBuffer = new float[size_x];
  float *imagOutBuffer = new float[size_x];

  // Local values
  //float *fft_real = new float[size_y];
  //float *fft_imag = new float[size_y];
  //#pragma omp parallel 
  {
    //#pragma omp for
  for(unsigned int x = 0; x < size_x; x++)
  {
    for(unsigned int y = 0; y < size_y; y++)
    {
      float realOut = 0;
      float imagOut = 0;
      // Compute the value for this index
      realOutBuffer[y] = 0.0f;
      imagOutBuffer[y] = 0.0f;

      // Compute the frequencies for this index
      for(unsigned int n = 0; n < size_y; n++)
      {
	float term = -2 * PI * y * n / size_y;
	float fft_real = 0;
	float fft_imag = 0;

	fft_real = cos(term);
	fft_imag = sin(term);

	//	realOutBuffer[y] += (real_image[x*size_x + n] * fft_real) - (imag_image[x*size_x + n] * fft_imag);
	//	imagOutBuffer[y] += (imag_image[x*size_x + n] * fft_real)
	//	+ (real_image[x*size_x + n] * fft_imag);
	realOut += (real_image[x*size_x + n] * fft_real) - (imag_image[x*size_x + n] * fft_imag);
	imagOut += (imag_image[x*size_x + n] * fft_real) + (real_image[x*size_x + n] * fft_imag);

      }

      realOutBuffer[y] = realOut;
      imagOutBuffer[y] = imagOut;

    }
    // Write the buffer back to were the original values were
    for(unsigned int y = 0; y < size_y; y++)
    {
      real_image[x*size_x + y] = realOutBuffer[y];
      imag_image[x*size_x + y] = imagOutBuffer[y];
    }
  }
  }

  // Reclaim some memory
  delete [] realOutBuffer;
  delete [] imagOutBuffer;
  //  delete [] fft_real;
  //  delete [] fft_imag;
}

// This is the same as the thing above, except it has a scaling factor added to it
void cpu_ifftx(float *real_image, float *imag_image, int size_x, int size_y)
{
  // Create some space for storing temporary values
  float *realOutBuffer = new float[size_x];
  float *imagOutBuffer = new float[size_x];
  //  float *fft_real = new float[size_y];
  //  float *fft_imag = new float[size_y];
  //#pragma omp parallel 
  {
    //#pragma omp for
  for(unsigned int x = 0; x < size_x; x++)
  {
    for(unsigned int y = 0; y < size_y; y++)
    {
      float realOut = 0;
      float imagOut = 0;
      // Compute the value for this index
      realOutBuffer[y] = 0.0f;
      imagOutBuffer[y] = 0.0f;

      for(unsigned int n = 0; n < size_y; n++)
      {
        // Compute the frequencies for this index
	float term = 2 * PI * y * n / size_y;
	float fft_real = 0;
	float fft_imag = 0;

	fft_real = cos(term);
	fft_imag = sin(term);

	realOut += (real_image[x*size_x + n] * fft_real) - (imag_image[x*size_x + n] * fft_imag);
	imagOut += (imag_image[x*size_x + n] * fft_real) + (real_image[x*size_x + n] * fft_imag);
      }
      realOutBuffer[y] = realOut;
      imagOutBuffer[y] = imagOut;

      // Incoporate the scaling factor here
      realOutBuffer[y] /= size_y;
      imagOutBuffer[y] /= size_y;

    }

    // Write the buffer back to were the original values were
    for(unsigned int y = 0; y < size_y; y++)
    {
      real_image[x*size_x + y] = realOutBuffer[y];
      imag_image[x*size_x + y] = imagOutBuffer[y];
    }
  }
  }
  // Reclaim some memory
  delete [] realOutBuffer;
  delete [] imagOutBuffer;
  //  delete [] fft_real;
  //  delete [] fft_imag;
}

void cpu_ffty(float *real_image, float *imag_image, int size_x, int size_y)
{
  // Allocate some space for temporary values
  float *realOutBuffer = new float[size_y];
  float *imagOutBuffer = new float[size_y];
  //  float *fft_real = new float[size_x];
  //  float *fft_imag = new float[size_x];
  //#pragma omp parallel 
  {
    //#pragma omp for
  for(unsigned int y = 0; y < size_y; y++)
  {
    for(unsigned int x = 0; x < size_x; x++)
    {
      float realOut = 0;
      float imagOut = 0;
      // Compute the value for this index
      realOutBuffer[x] = 0.0f;
      imagOutBuffer[x] = 0.0f;
      // Compute the frequencies for this index
      for(unsigned int n = 0; n < size_y; n++)
      {
	float term = -2 * PI * x * n / size_x;
	float fft_real = 0;
	float fft_imag = 0;

	fft_real = cos(term);
	fft_imag = sin(term);

	realOut += (real_image[n*size_x + y] * fft_real) - (imag_image[n*size_x + y] * fft_imag);
	imagOut += (imag_image[n*size_x + y] * fft_real) + (real_image[n*size_x + y] * fft_imag);
      }
      realOutBuffer[x] = realOut;
      imagOutBuffer[x] = imagOut;
    }
    // Write the buffer back to were the original values were
    for(unsigned int x = 0; x < size_x; x++)
    {
      real_image[x*size_x + y] = realOutBuffer[x];
      imag_image[x*size_x + y] = imagOutBuffer[x];
      }
  }
  }
  // Reclaim some memory
  delete [] realOutBuffer;
  delete [] imagOutBuffer;
  //  delete [] fft_real;
  //  delete [] fft_imag;
}

// This is the same as the thing about it, but it includes a scaling factor
void cpu_iffty(float *real_image, float *imag_image, int size_x, int size_y)
{
  // Create some space for storing temporary values
  float *realOutBuffer = new float[size_y];
  float *imagOutBuffer = new float[size_y];
  //  float *fft_real = new float[size_x];
  //  float *fft_imag = new float[size_x];
  //#pragma omp parallel 
  {
    //#pragma omp for
  for(unsigned int y = 0; y < size_y; y++)
  {
    for(unsigned int x = 0; x < size_x; x++)
    {
      float realOut = 0;
      float imagOut = 0;
      // Compute the value for this index
      realOutBuffer[x] = 0.0f;
      imagOutBuffer[x] = 0.0f;
      // Compute the frequencies for this index
      for(unsigned int n = 0; n < size_y; n++)
      {
	// Note that the negative sign goes away for the term
	float term = 2 * PI * x * n / size_x;
	float fft_real = 0;
	float fft_imag = 0;

	fft_real = cos(term);
	fft_imag = sin(term);

	realOut += (real_image[n*size_x + y] * fft_real) - (imag_image[n*size_x + y] * fft_imag);
	imagOut += (imag_image[n*size_x + y] * fft_real) + (real_image[n*size_x + y] * fft_imag);
      }
      realOutBuffer[x] = realOut;
      imagOutBuffer[x] = imagOut;

      // Incorporate the scaling factor here
      realOutBuffer[x] /= size_x;
      imagOutBuffer[x] /= size_x;
    }

    // Write the buffer back to were the original values were
    for(unsigned int x = 0; x < size_x; x++)
    {
      real_image[x*size_x + y] = realOutBuffer[x];
      imag_image[x*size_x + y] = imagOutBuffer[x];
    }
  }
  }
  // Reclaim some memory
  delete [] realOutBuffer;
  delete [] imagOutBuffer;
  //  delete [] fft_real;
  //  delete [] fft_imag;
}

void cpu_filter(float *real_image, float *imag_image, int size_x, int size_y)
{
  int eightX = size_x/8;
  int eight7X = size_x - eightX;
  int eightY = size_y/8;
  int eight7Y = size_y - eightY;
  for(unsigned int x = 0; x < size_x; x++)
  {
    for(unsigned int y = 0; y < size_y; y++)
    {
      if(!(x < eightX && y < eightY) &&
	 !(x < eightX && y >= eight7Y) &&
	 !(x >= eight7X && y < eightY) &&
	 !(x >= eight7X && y >= eight7Y))
      {
	// Zero out these values
	real_image[y*size_x + x] = 0;
	imag_image[y*size_x + x] = 0;
      }
    }
  }
}

float imageCleaner(float *real_image, float *imag_image, int size_x, int size_y)
{
  // These are used for timing
  struct timeval tv1, tv2;
  struct timezone tz1, tz2;

  printf (" size_x %d size_y %d\n", size_x, size_y);
  // Start timing
  gettimeofday(&tv1,&tz1);

  // Perform fft with respect to the x direction
  cpu_fftx(real_image, imag_image, size_x, size_y);
  
  printf (" fftx done\n");

  // Perform fft with respect to the y direction
  cpu_ffty(real_image, imag_image, size_x, size_y);

  printf (" ffty done\n");
  // Filter the transformed image
  cpu_filter(real_image, imag_image, size_x, size_y);
  printf (" cpu_filter done\n");
  // Perform an inverse fft with respect to the x direction
  cpu_ifftx(real_image, imag_image, size_x, size_y);
  printf (" ifftx done\n");
  // Perform an inverse fft with respect to the y direction
  cpu_iffty(real_image, imag_image, size_x, size_y);
  printf (" iffty done\n");
  // End timing
  gettimeofday(&tv2,&tz2);

  // Compute the time difference in micro-seconds
  float execution = ((tv2.tv_sec-tv1.tv_sec)*1000000+(tv2.tv_usec-tv1.tv_usec));
  // Convert to milli-seconds
  execution /= 1000;
  // Print some output
  printf("OPTIMIZED IMPLEMENTATION STATISTICS:\n");
  printf("  Optimized Kernel Execution Time: %f ms\n\n", execution);
  return execution;
}
