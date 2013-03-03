#include "ImageCleaner.h"
#include <math.h>
#include <omp.h>
#include <sys/time.h>
#include <stdio.h>
#include <string.h>

#define PI	3.14159265
unsigned int chunks = 0;

void cpu_fftx(float *real_image, float *imag_image, int size_x, int size_y)
{
  float *realInBuffer0 = new float[size_x];
  float *imagInBuffer0 = new float[size_x];
  float *realInBuffer1 = new float[size_x];
  float *imagInBuffer1 = new float[size_x];
  float *realInBuffer2 = new float[size_x];
  float *imagInBuffer2 = new float[size_x];
  float *realInBuffer3 = new float[size_x];
  float *imagInBuffer3 = new float[size_x];
  float *realInBuffer4 = new float[size_x];
  float *imagInBuffer4 = new float[size_x];
  float *realInBuffer5 = new float[size_x];
  float *imagInBuffer5 = new float[size_x];
  float *realInBuffer6 = new float[size_x];
  float *imagInBuffer6 = new float[size_x];
  float *realInBuffer7 = new float[size_x];
  float *imagInBuffer7 = new float[size_x];


  // Local values
  for(unsigned int x = 0; x < size_x; x+=8)
  {
    // copy from real image to shared buffer so threads
    // can work in parallel and exploit shared memory
    for (unsigned int i = 0; i < size_y; i++)
      {
	realInBuffer0[i] = real_image[x*size_x + i];
	imagInBuffer0[i] = imag_image[x*size_x + i];
	
	realInBuffer1[i] = real_image[(x+1)*size_x + i];
	imagInBuffer1[i] = imag_image[(x+1)*size_x + i];

	realInBuffer2[i] = real_image[(x+2)*size_x + i];
	imagInBuffer2[i] = imag_image[(x+2)*size_x + i];

	realInBuffer3[i] = real_image[(x+3)*size_x + i];
	imagInBuffer3[i] = imag_image[(x+3)*size_x + i];

	realInBuffer4[i] = real_image[(x+4)*size_x + i];
	imagInBuffer4[i] = imag_image[(x+4)*size_x + i];
	
	realInBuffer5[i] = real_image[(x+5)*size_x + i];
	imagInBuffer5[i] = imag_image[(x+5)*size_x + i];

	realInBuffer6[i] = real_image[(x+6)*size_x + i];
	imagInBuffer6[i] = imag_image[(x+6)*size_x + i];

	realInBuffer7[i] = real_image[(x+7)*size_x + i];
	imagInBuffer7[i] = imag_image[(x+7)*size_x + i];

      }

#pragma omp parallel 
  {
#pragma omp for schedule(dynamic, chunks)
    for(unsigned int y = 0; y < size_y; y++)
    {
      float realOut1 = 0;
      float imagOut1 = 0;
      float realOut0 = 0;
      float imagOut0 = 0;
      float realOut2 = 0;
      float imagOut2 = 0;
      float realOut3 = 0;
      float imagOut3 = 0;
      float realOut4 = 0;
      float imagOut4 = 0;
      float realOut5 = 0;
      float imagOut5 = 0;
      float realOut6 = 0;
      float imagOut6 = 0;
      float realOut7 = 0;
      float imagOut7 = 0;


      // Compute the value for this index
      // Compute the frequencies for this index
      for(unsigned int n = 0; n < size_y; n++)
      {
	float term = -2 * PI * y * n / size_y;
	float fft_real = 0;
	float fft_imag = 0;

	fft_real = cos(term);
	fft_imag = sin(term);

	realOut0 += (realInBuffer0[n] * fft_real) - (imagInBuffer0[n] * fft_imag);
	imagOut0 += (imagInBuffer0[n] * fft_real) + (realInBuffer0[n] * fft_imag);

	realOut1 += (realInBuffer1[n] * fft_real) - (imagInBuffer1[n] * fft_imag);
	imagOut1 += (imagInBuffer1[n] * fft_real) + (realInBuffer1[n] * fft_imag);

	realOut2 += (realInBuffer2[n] * fft_real) - (imagInBuffer2[n] * fft_imag);
	imagOut2 += (imagInBuffer2[n] * fft_real) + (realInBuffer2[n] * fft_imag);

	realOut3 += (realInBuffer3[n] * fft_real) - (imagInBuffer3[n] * fft_imag);
	imagOut3 += (imagInBuffer3[n] * fft_real) + (realInBuffer3[n] * fft_imag);

	realOut4 += (realInBuffer4[n] * fft_real) - (imagInBuffer4[n] * fft_imag);
	imagOut4 += (imagInBuffer4[n] * fft_real) + (realInBuffer4[n] * fft_imag);

	realOut5 += (realInBuffer5[n] * fft_real) - (imagInBuffer5[n] * fft_imag);
	imagOut5 += (imagInBuffer5[n] * fft_real) + (realInBuffer5[n] * fft_imag);

	realOut6 += (realInBuffer6[n] * fft_real) - (imagInBuffer6[n] * fft_imag);
	imagOut6 += (imagInBuffer6[n] * fft_real) + (realInBuffer6[n] * fft_imag);

	realOut7 += (realInBuffer7[n] * fft_real) - (imagInBuffer7[n] * fft_imag);
	imagOut7 += (imagInBuffer7[n] * fft_real) + (realInBuffer7[n] * fft_imag);


      }

      // Write the buffer back after computing
      real_image[x*size_x + y] = realOut0;
      imag_image[x*size_x + y] = imagOut0;

      real_image[(x+1)*size_x + y] = realOut1;
      imag_image[(x+1)*size_x + y] = imagOut1;

      real_image[(x+2)*size_x + y] = realOut2;
      imag_image[(x+2)*size_x + y] = imagOut2;

      real_image[(x+3)*size_x + y] = realOut3;
      imag_image[(x+3)*size_x + y] = imagOut3;

      real_image[(x+4)*size_x + y] = realOut4;
      imag_image[(x+4)*size_x + y] = imagOut4;

      real_image[(x+5)*size_x + y] = realOut5;
      imag_image[(x+5)*size_x + y] = imagOut5;

      real_image[(x+6)*size_x + y] = realOut6;
      imag_image[(x+6)*size_x + y] = imagOut6;

      real_image[(x+7)*size_x + y] = realOut7;
      imag_image[(x+7)*size_x + y] = imagOut7;

    }
  }
  }

  // Reclaim some memory
  delete [] realInBuffer0;
  delete [] imagInBuffer0;
  delete [] realInBuffer1;
  delete [] imagInBuffer1;
  delete [] realInBuffer2;
  delete [] imagInBuffer2;
  delete [] realInBuffer3;
  delete [] imagInBuffer3;
  delete [] realInBuffer4;
  delete [] imagInBuffer4;
  delete [] realInBuffer5;
  delete [] imagInBuffer5;
  delete [] realInBuffer6;
  delete [] imagInBuffer6;
  delete [] realInBuffer7;
  delete [] imagInBuffer7;
}

// This is the same as the thing above, except it has a scaling factor added to it
void cpu_ifftx(float *real_image, float *imag_image, int size_x, int size_y)
{
  // Create some space for storing temporary values
  float *realInBuffer0 = new float[size_x];
  float *imagInBuffer0 = new float[size_x];

  float *realInBuffer1 = new float[size_x];
  float *imagInBuffer1 = new float[size_x];

  float *realInBuffer2 = new float[size_x];
  float *imagInBuffer2 = new float[size_x];

  float *realInBuffer3 = new float[size_x];
  float *imagInBuffer3 = new float[size_x];

  float *realInBuffer4 = new float[size_x];
  float *imagInBuffer4 = new float[size_x];

  float *realInBuffer5 = new float[size_x];
  float *imagInBuffer5 = new float[size_x];

  float *realInBuffer6 = new float[size_x];
  float *imagInBuffer6 = new float[size_x];

  float *realInBuffer7 = new float[size_x];
  float *imagInBuffer7 = new float[size_x];

  for(unsigned int x = 0; x < size_x; x+=8)
  {

    // copy from real image to shared buffer so threads
    // can work in parallel and exploit shared memory
    for (unsigned int i = 0; i < size_y; i++)
      {
	realInBuffer0[i] = real_image[x*size_x + i];
	imagInBuffer0[i] = imag_image[x*size_x + i];

	realInBuffer1[i] = real_image[(x+1)*size_x + i];
	imagInBuffer1[i] = imag_image[(x+1)*size_x + i];

	realInBuffer2[i] = real_image[(x+2)*size_x + i];
	imagInBuffer2[i] = imag_image[(x+2)*size_x + i];

	realInBuffer3[i] = real_image[(x+3)*size_x + i];
	imagInBuffer3[i] = imag_image[(x+3)*size_x + i];

	realInBuffer4[i] = real_image[(x+4)*size_x + i];
	imagInBuffer4[i] = imag_image[(x+4)*size_x + i];

	realInBuffer5[i] = real_image[(x+5)*size_x + i];
	imagInBuffer5[i] = imag_image[(x+5)*size_x + i];

	realInBuffer6[i] = real_image[(x+6)*size_x + i];
	imagInBuffer6[i] = imag_image[(x+6)*size_x + i];

	realInBuffer7[i] = real_image[(x+7)*size_x + i];
	imagInBuffer7[i] = imag_image[(x+7)*size_x + i];

      }
  #pragma omp parallel
  {
    #pragma omp for schedule(dynamic, chunks)
    for(unsigned int y = 0; y < size_y; y++)
    {
      float realOut0 = 0;
      float imagOut0 = 0;
      float realOut1 = 0;
      float imagOut1 = 0;
      float realOut2 = 0;
      float imagOut2 = 0;
      float realOut3 = 0;
      float imagOut3 = 0;
      float realOut4 = 0;
      float imagOut4 = 0;
      float realOut5 = 0;
      float imagOut5 = 0;
      float realOut6 = 0;
      float imagOut6 = 0;
      float realOut7 = 0;
      float imagOut7 = 0;
   
      // Compute the value for this index
      for(unsigned int n = 0; n < size_y; n++)
      {
        // Compute the frequencies for this index
	float term = 2 * PI * y * n / size_y;
	float fft_real = 0;
	float fft_imag = 0;

	fft_real = cos(term);
	fft_imag = sin(term);
	realOut0 += (realInBuffer0[n] * fft_real) - (imagInBuffer0[n] * fft_imag);
	imagOut0 += (imagInBuffer0[n] * fft_real) + (realInBuffer0[n] * fft_imag);

	realOut1 += (realInBuffer1[n] * fft_real) - (imagInBuffer1[n] * fft_imag);
	imagOut1 += (imagInBuffer1[n] * fft_real) + (realInBuffer1[n] * fft_imag);

	realOut2 += (realInBuffer2[n] * fft_real) - (imagInBuffer2[n] * fft_imag);
	imagOut2 += (imagInBuffer2[n] * fft_real) + (realInBuffer2[n] * fft_imag);

	realOut3 += (realInBuffer3[n] * fft_real) - (imagInBuffer3[n] * fft_imag);
	imagOut3 += (imagInBuffer3[n] * fft_real) + (realInBuffer3[n] * fft_imag);

	realOut4 += (realInBuffer4[n] * fft_real) - (imagInBuffer4[n] * fft_imag);
	imagOut4 += (imagInBuffer4[n] * fft_real) + (realInBuffer4[n] * fft_imag);

	realOut5 += (realInBuffer5[n] * fft_real) - (imagInBuffer5[n] * fft_imag);
	imagOut5 += (imagInBuffer5[n] * fft_real) + (realInBuffer5[n] * fft_imag);

	realOut6 += (realInBuffer6[n] * fft_real) - (imagInBuffer6[n] * fft_imag);
	imagOut6 += (imagInBuffer6[n] * fft_real) + (realInBuffer6[n] * fft_imag);

	realOut7 += (realInBuffer7[n] * fft_real) - (imagInBuffer7[n] * fft_imag);
	imagOut7 += (imagInBuffer7[n] * fft_real) + (realInBuffer7[n] * fft_imag);


      }

      realOut0 /= size_y;
      imagOut0 /= size_y;

      realOut1 /= size_y;
      imagOut1 /= size_y;

      realOut2 /= size_y;
      imagOut2 /= size_y;

      realOut3 /= size_y;
      imagOut3 /= size_y;

      realOut4 /= size_y;
      imagOut4 /= size_y;

      realOut5 /= size_y;
      imagOut5 /= size_y;

      realOut6 /= size_y;
      imagOut6 /= size_y;

      realOut7 /= size_y;
      imagOut7 /= size_y;
      

      real_image[x*size_x + y] = realOut0;
      imag_image[x*size_x + y] = imagOut0;

      real_image[(x+1)*size_x + y] = realOut1;
      imag_image[(x+1)*size_x + y] = imagOut1;

      real_image[(x+2)*size_x + y] = realOut2;
      imag_image[(x+2)*size_x + y] = imagOut2;

      real_image[(x+3)*size_x + y] = realOut3;
      imag_image[(x+3)*size_x + y] = imagOut3;

      real_image[(x+4)*size_x + y] = realOut4;
      imag_image[(x+4)*size_x + y] = imagOut4;

      real_image[(x+5)*size_x + y] = realOut5;
      imag_image[(x+5)*size_x + y] = imagOut5;

      real_image[(x+6)*size_x + y] = realOut6;
      imag_image[(x+6)*size_x + y] = imagOut6;

      real_image[(x+7)*size_x + y] = realOut7;
      imag_image[(x+7)*size_x + y] = imagOut7;

    }
  }
  }

  // Reclaim some memory
  delete [] realInBuffer0;
  delete [] imagInBuffer0;
  delete [] realInBuffer1;
  delete [] imagInBuffer1;
  delete [] realInBuffer2;
  delete [] imagInBuffer2;
  delete [] realInBuffer3;
  delete [] imagInBuffer3;
  delete [] realInBuffer4;
  delete [] imagInBuffer4;
  delete [] realInBuffer5;
  delete [] imagInBuffer5;
  delete [] realInBuffer6;
  delete [] imagInBuffer6;
  delete [] realInBuffer7;
  delete [] imagInBuffer7;

}

void cpu_ffty(float *real_image, float *imag_image, int size_x, int size_y)
{
  // Allocate some space for temporary values
  float *realInBuffer0 = new float[size_y];
  float *imagInBuffer0 = new float[size_y];
  float *realInBuffer1 = new float[size_y];
  float *imagInBuffer1 = new float[size_y];
  float *realInBuffer2 = new float[size_y];
  float *imagInBuffer2 = new float[size_y];
  float *realInBuffer3 = new float[size_y];
  float *imagInBuffer3 = new float[size_y];
  float *realInBuffer4 = new float[size_y];
  float *imagInBuffer4 = new float[size_y];
  float *realInBuffer5 = new float[size_y];
  float *imagInBuffer5 = new float[size_y];
  float *realInBuffer6 = new float[size_y];
  float *imagInBuffer6 = new float[size_y];
  float *realInBuffer7 = new float[size_y];
  float *imagInBuffer7 = new float[size_y];

  for(unsigned int y = 0; y < size_y; y+=8)
  {
    // copy from real image to shared buffer so threads
    // can work in parallel and exploit shared memory
    for (unsigned int i = 0 ;  i < size_x ; i++) 
      {
	realInBuffer0[i] = real_image[i*size_x + y];
	imagInBuffer0[i] = imag_image[i*size_x + y];

	realInBuffer1[i] = real_image[i*size_x + (y+1)];
	imagInBuffer1[i] = imag_image[i*size_x + (y+1)];

	realInBuffer2[i] = real_image[i*size_x + (y+2)];
	imagInBuffer2[i] = imag_image[i*size_x + (y+2)];

	realInBuffer3[i] = real_image[i*size_x + (y+3)];
	imagInBuffer3[i] = imag_image[i*size_x + (y+3)];

	realInBuffer4[i] = real_image[i*size_x + (y+4)];
	imagInBuffer4[i] = imag_image[i*size_x + (y+4)];

	realInBuffer5[i] = real_image[i*size_x + (y+5)];
	imagInBuffer5[i] = imag_image[i*size_x + (y+5)];

	realInBuffer6[i] = real_image[i*size_x + (y+6)];
	imagInBuffer6[i] = imag_image[i*size_x + (y+6)];

	realInBuffer7[i] = real_image[i*size_x + (y+7)];
	imagInBuffer7[i] = imag_image[i*size_x + (y+7)];

      }
  #pragma omp parallel
  {
    #pragma omp for schedule(dynamic, chunks)
    for(unsigned int x = 0; x < size_x; x++)
    {
      float realOut0 = 0;
      float imagOut0 = 0;
      float realOut1 = 0;
      float imagOut1 = 0;
      float realOut2 = 0;
      float imagOut2 = 0;
      float realOut3 = 0;
      float imagOut3 = 0;
      float realOut4 = 0;
      float imagOut4 = 0;
      float realOut5 = 0;
      float imagOut5 = 0;
      float realOut6 = 0;
      float imagOut6 = 0;
      float realOut7 = 0;
      float imagOut7 = 0;


      // Compute the value for this index
      // Compute the frequencies for this index
      for(unsigned int n = 0; n < size_y; n++)
      {
	float term = -2 * PI * x * n / size_x;
	float fft_real = 0;
	float fft_imag = 0;

	fft_real = cos(term);
	fft_imag = sin(term);

	realOut0 += (realInBuffer0[n] * fft_real) - (imagInBuffer0[n] * fft_imag);
	imagOut0 += (imagInBuffer0[n] * fft_real) + (realInBuffer0[n] * fft_imag);

	realOut1 += (realInBuffer1[n] * fft_real) - (imagInBuffer1[n] * fft_imag);
	imagOut1 += (imagInBuffer1[n] * fft_real) +  (realInBuffer1[n] * fft_imag);

	realOut2 += (realInBuffer2[n] * fft_real) - (imagInBuffer2[n] * fft_imag);
	imagOut2 += (imagInBuffer2[n] * fft_real) +  (realInBuffer2[n] * fft_imag);

	realOut3 += (realInBuffer3[n] * fft_real) - (imagInBuffer3[n] * fft_imag);
	imagOut3 += (imagInBuffer3[n] * fft_real) +  (realInBuffer3[n] * fft_imag);

	realOut4 += (realInBuffer4[n] * fft_real) - (imagInBuffer4[n] * fft_imag);
	imagOut4 += (imagInBuffer4[n] * fft_real) +  (realInBuffer4[n] * fft_imag);

	realOut5 += (realInBuffer5[n] * fft_real) - (imagInBuffer5[n] * fft_imag);
	imagOut5 += (imagInBuffer5[n] * fft_real) +  (realInBuffer5[n] * fft_imag);

	realOut6 += (realInBuffer6[n] * fft_real) - (imagInBuffer6[n] * fft_imag);
	imagOut6 += (imagInBuffer6[n] * fft_real) +  (realInBuffer6[n] * fft_imag);

	realOut7 += (realInBuffer7[n] * fft_real) - (imagInBuffer7[n] * fft_imag);
	imagOut7 += (imagInBuffer7[n] * fft_real) +  (realInBuffer7[n] * fft_imag);

      }

      real_image[x*size_x + y] = realOut0;
      imag_image[x*size_x + y] = imagOut0;

      real_image[x*size_x + (y+1)] = realOut1;
      imag_image[x*size_x + (y+1)] = imagOut1;

      real_image[x*size_x + (y+2)] = realOut2;
      imag_image[x*size_x + (y+2)] = imagOut2;

      real_image[x*size_x + (y+3)] = realOut3;
      imag_image[x*size_x + (y+3)] = imagOut3;

      real_image[x*size_x + (y+4)] = realOut4;
      imag_image[x*size_x + (y+4)] = imagOut4;

      real_image[x*size_x + (y+5)] = realOut5;
      imag_image[x*size_x + (y+5)] = imagOut5;

      real_image[x*size_x + (y+6)] = realOut6;
      imag_image[x*size_x + (y+6)] = imagOut6;

      real_image[x*size_x + (y+7)] = realOut7;
      imag_image[x*size_x + (y+7)] = imagOut7;

    }
  }
  }
  // Reclaim some memory
  delete [] realInBuffer0;
  delete [] imagInBuffer0;
  delete [] realInBuffer1;
  delete [] imagInBuffer1;
  delete [] realInBuffer2;
  delete [] imagInBuffer2;
  delete [] realInBuffer3;
  delete [] imagInBuffer3;
  delete [] realInBuffer4;
  delete [] imagInBuffer4;
  delete [] realInBuffer5;
  delete [] imagInBuffer5;
  delete [] realInBuffer6;
  delete [] imagInBuffer6;
  delete [] realInBuffer7;
  delete [] imagInBuffer7;

}

// This is the same as the thing about it, but it includes a scaling factor
void cpu_iffty(float *real_image, float *imag_image, int size_x, int size_y)
{
  // Create some space for storing temporary values
  float *realInBuffer0 = new float[size_y];
  float *imagInBuffer0 = new float[size_y];
  float *realInBuffer1 = new float[size_y];
  float *imagInBuffer1 = new float[size_y];
  float *realInBuffer2 = new float[size_y];
  float *imagInBuffer2 = new float[size_y];
  float *realInBuffer3 = new float[size_y];
  float *imagInBuffer3 = new float[size_y];
  float *realInBuffer4 = new float[size_y];
  float *imagInBuffer4 = new float[size_y];
  float *realInBuffer5 = new float[size_y];
  float *imagInBuffer5 = new float[size_y];
  float *realInBuffer6 = new float[size_y];
  float *imagInBuffer6 = new float[size_y];
  float *realInBuffer7 = new float[size_y];
  float *imagInBuffer7 = new float[size_y];


  for(unsigned int y = 0; y < size_y; y+=8)
  {
    for (unsigned int i = 0; i < size_x; i++)
      {
	realInBuffer0[i] = real_image[i*size_x + y];
	imagInBuffer0[i] = imag_image[i*size_x + y];

	realInBuffer1[i] = real_image[i*size_x + (y+1)];
	imagInBuffer1[i] = imag_image[i*size_x + (y+1)];

	realInBuffer2[i] = real_image[i*size_x + (y+2)];
	imagInBuffer2[i] = imag_image[i*size_x + (y+2)];

	realInBuffer3[i] = real_image[i*size_x + (y+3)];
	imagInBuffer3[i] = imag_image[i*size_x + (y+3)];

	realInBuffer4[i] = real_image[i*size_x + (y+4)];
	imagInBuffer4[i] = imag_image[i*size_x + (y+4)];

	realInBuffer5[i] = real_image[i*size_x + (y+5)];
	imagInBuffer5[i] = imag_image[i*size_x + (y+5)];

	realInBuffer6[i] = real_image[i*size_x + (y+6)];
	imagInBuffer6[i] = imag_image[i*size_x + (y+6)];

	realInBuffer7[i] = real_image[i*size_x + (y+7)];
	imagInBuffer7[i] = imag_image[i*size_x + (y+7)];

      }
  #pragma omp parallel
  {
    #pragma omp for schedule(dynamic, chunks)
    for(unsigned int x = 0; x < size_x; x++)
    {
      float realOut0 = 0;
      float imagOut0 = 0;
      float realOut1 = 0;
      float imagOut1 = 0;
      float realOut2 = 0;
      float imagOut2 = 0;
      float realOut3 = 0;
      float imagOut3 = 0;
      float realOut4 = 0;
      float imagOut4 = 0;
      float realOut5 = 0;
      float imagOut5 = 0;
      float realOut6 = 0;
      float imagOut6 = 0;
      float realOut7 = 0;
      float imagOut7 = 0;

      // Compute the value for this index
      // Compute the frequencies for this index
      for(unsigned int n = 0; n < size_y; n++)
      {
	// Note that the negative sign goes away for the term
	float term = 2 * PI * x * n / size_x;
	float fft_real = 0;
	float fft_imag = 0;

	fft_real = cos(term);
	fft_imag = sin(term);

	realOut0 += (realInBuffer0[n] * fft_real) - (imagInBuffer0[n] * fft_imag);
	imagOut0 += (imagInBuffer0[n] * fft_real) + (realInBuffer0[n] * fft_imag);

	realOut1 += (realInBuffer1[n] * fft_real) - (imagInBuffer1[n] * fft_imag);
	imagOut1 += (imagInBuffer1[n] * fft_real) + (realInBuffer1[n] * fft_imag);

	realOut2 += (realInBuffer2[n] * fft_real) - (imagInBuffer2[n] * fft_imag);
	imagOut2 += (imagInBuffer2[n] * fft_real) + (realInBuffer2[n] * fft_imag);

	realOut3 += (realInBuffer3[n] * fft_real) - (imagInBuffer3[n] * fft_imag);
	imagOut3 += (imagInBuffer3[n] * fft_real) + (realInBuffer3[n] * fft_imag);

	realOut4 += (realInBuffer4[n] * fft_real) - (imagInBuffer4[n] * fft_imag);
	imagOut4 += (imagInBuffer4[n] * fft_real) + (realInBuffer4[n] * fft_imag);

	realOut5 += (realInBuffer5[n] * fft_real) - (imagInBuffer5[n] * fft_imag);
	imagOut5 += (imagInBuffer5[n] * fft_real) + (realInBuffer5[n] * fft_imag);

	realOut6 += (realInBuffer6[n] * fft_real) - (imagInBuffer6[n] * fft_imag);
	imagOut6 += (imagInBuffer6[n] * fft_real) + (realInBuffer6[n] * fft_imag);

	realOut7 += (realInBuffer7[n] * fft_real) - (imagInBuffer7[n] * fft_imag);
	imagOut7 += (imagInBuffer7[n] * fft_real) + (realInBuffer7[n] * fft_imag);


      }

      realOut0 /= size_x;
      imagOut0 /= size_x;

      realOut1 /= size_x;
      imagOut1 /= size_x;

      realOut2 /= size_x;
      imagOut2 /= size_x;

      realOut3 /= size_x;
      imagOut3 /= size_x;

      realOut4 /= size_x;
      imagOut4 /= size_x;

      realOut5 /= size_x;
      imagOut5 /= size_x;

      realOut6 /= size_x;
      imagOut6 /= size_x;

      realOut7 /= size_x;
      imagOut7 /= size_x;



      real_image[x*size_x + y] = realOut0;
      imag_image[x*size_x + y] = imagOut0;

      real_image[x*size_x + (y+1)] = realOut1;
      imag_image[x*size_x + (y+1)] = imagOut1;

      real_image[x*size_x + (y+2)] = realOut2;
      imag_image[x*size_x + (y+2)] = imagOut2;

      real_image[x*size_x + (y+3)] = realOut3;
      imag_image[x*size_x + (y+3)] = imagOut3;

      real_image[x*size_x + (y+4)] = realOut4;
      imag_image[x*size_x + (y+4)] = imagOut4;

      real_image[x*size_x + (y+5)] = realOut5;
      imag_image[x*size_x + (y+5)] = imagOut5;

      real_image[x*size_x + (y+6)] = realOut6;
      imag_image[x*size_x + (y+6)] = imagOut6;

      real_image[x*size_x + (y+7)] = realOut7;
      imag_image[x*size_x + (y+7)] = imagOut7;

    }
  }
  }
  // Reclaim some memory
  delete [] realInBuffer0;
  delete [] imagInBuffer0;
  delete [] realInBuffer1;
  delete [] imagInBuffer1;
  delete [] realInBuffer2;
  delete [] imagInBuffer2;
  delete [] realInBuffer3;
  delete [] imagInBuffer3;
  delete [] realInBuffer4;
  delete [] imagInBuffer4;
  delete [] realInBuffer5;
  delete [] imagInBuffer5;
  delete [] realInBuffer6;
  delete [] imagInBuffer6;
  delete [] realInBuffer7;
  delete [] imagInBuffer7;

}

void cpu_filter(float *real_image, float *imag_image, int size_x, int size_y)
{
  int eightX = size_x/8;
  int eight7X = size_x - eightX;
  int eightY = size_y/8;
  int eight7Y = size_y - eightY;

  for(unsigned int x = 0; x < size_x; x+=4)
  {
#pragma omp parallel
  {
#pragma omp for schedule(dynamic, chunks)
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
      if(!((x+1) < eightX && y < eightY) &&
	 !((x+1) < eightX && y >= eight7Y) &&
	 !((x+1) >= eight7X && y < eightY) &&
	 !((x+1) >= eight7X && y >= eight7Y))
      {
	// Zero out these values
	real_image[y*size_x + x+1] = 0;
	imag_image[y*size_x + x+1] = 0;
      }
      if(!((x+2) < eightX && y < eightY) &&
	 !((x+2) < eightX && y >= eight7Y) &&
	 !((x+2) >= eight7X && y < eightY) &&
	 !((x+2) >= eight7X && y >= eight7Y))
      {
	// Zero out these values
	real_image[y*size_x + x+2] = 0;
	imag_image[y*size_x + x+2] = 0;
      }

      if(!((x+3) < eightX && y < eightY) &&
	 !((x+3) < eightX && y >= eight7Y) &&
	 !((x+3) >= eight7X && y < eightY) &&
	 !((x+3) >= eight7X && y >= eight7Y))
      {
	// Zero out these values
	real_image[y*size_x + x+3] = 0;
	imag_image[y*size_x + x+3] = 0;
      }
    }
  }
  }
}

float imageCleaner(float *real_image, float *imag_image, int size_x, int size_y)
{
  // These are used for timing
  struct timeval tv1, tv2;
  struct timezone tz1, tz2;
  int maxt = 8;

  maxt = omp_get_max_threads();
  chunks = size_x/maxt;

  // Start timing
  gettimeofday(&tv1,&tz1);

  // Perform fft with respect to the x direction
  cpu_fftx(real_image, imag_image, size_x, size_y);
  
  // Perform fft with respect to the y direction
  cpu_ffty(real_image, imag_image, size_x, size_y);

  // Filter the transformed image
  cpu_filter(real_image, imag_image, size_x, size_y);

  // Perform an inverse fft with respect to the x direction
  cpu_ifftx(real_image, imag_image, size_x, size_y);

  // Perform an inverse fft with respect to the y direction
  cpu_iffty(real_image, imag_image, size_x, size_y);

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
