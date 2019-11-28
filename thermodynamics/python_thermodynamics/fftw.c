/*
 * main.c
 *
 *  Created on: Jul 10, 2012
 *      Author: toto
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <fftw3.h>

int main( int argc, char** argv )
{
	fftwf_complex *data;
	fftwf_complex *fdata;
	fftwf_complex *idata;
	int ndata, nfft, i;
	unsigned int seed = 1234;

	if(argc<2)
		return(1);
	sscanf(argv[1], "%i", &ndata);
	sscanf(argv[2], "%i", &nfft);

	//allocate memory for data
	data = ( fftwf_complex* ) fftwf_malloc( sizeof( fftwf_complex ) * ndata );

	// CREATE INPUT DATA
	srand ( seed );
	for( i=0; i<ndata; i++ ) {
		data[i][0] = (double) (rand()/12345678);	//real data
		data[i][1] = 0.0;		//imaginer data
	}
	//print input data
	printf("Input Data \n");
	for( i = 0 ; i < ndata ; i++ )
		printf( "%i.  %2.5f  %2.5f \n", i, data[i][0], data[i][1]);

	// FFT PROCESS
	fdata = fftwf_data(data, ndata, nfft);
	printf("\nFFT Result \n");
	for( i = 0 ; i < nfft ; i++ )
		printf( "%i.  %2.5f  %2.5f \n", i, fdata[i][0], fdata[i][1]);

	// INVERS FFT PROCESS
	idata = ifftwf_data(fdata, ndata, nfft);
	printf("\nIFFT Result \n");
	for( i = 0 ; i < ndata ; i++ )
		printf( "%i.  %2.5f  %2.5f \n", i, idata[i][0], idata[i][1]);
	printf("\n");

	fftwf_free( data );
	fftwf_free( fdata );
	fftwf_free( idata );

	return 0;
}
