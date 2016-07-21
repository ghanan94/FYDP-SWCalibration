#include <stdio.h>
#include <stdint.h>
/*
	Rishi says:
	Module sourced from: https://sestevenson.wordpress.com/implementation-of-fir-filtering-in-c-part-1/
	Applies an FIR filter when "ApplyFIR" method is called.
	Generated coeffs using MATLAB, but it has it's own attenuation nevertheless.
*/


//////////////////////////////////////////////////////////////
//  Filter Code Definitions
//////////////////////////////////////////////////////////////
 
// maximum number of inputs that can be handled
// in one function call
#define MAX_INPUT_LEN   81
// maximum length of filter than can be handled
#define MAX_FLT_LEN     64
// buffer to hold all of the input samples
#define BUFFER_LEN      (MAX_FLT_LEN - 1 + MAX_INPUT_LEN)
 
// array to hold input samples
double insamp[ BUFFER_LEN ];
 
// FIR init
void firFloatInit( void )
{
    memset( insamp, 0, sizeof( insamp ) );
}
 
// the FIR filter function
void firFloat( double *coeffs, double *input, double *output,
      int length, int filterLength )
{
    double acc;     // accumulator for MACs
    double *coeffp; // pointer to coefficients
    double *inputp; // pointer to input samples
    int n;
    int k;
 
    // put the new samples at the high end of the buffer
    memcpy( &insamp[filterLength - 1], input,
            length * sizeof(double) );
 
    // apply the filter to each input sample
    for ( n = 0; n < length; n++ ) {
        // calculate output n
        coeffp = coeffs;
        inputp = &insamp[filterLength - 1 + n];
        acc = 0;
        for ( k = 0; k < filterLength; k++ ) {
            acc += (*coeffp++) * (*inputp--);
        }
        output[n] = acc;
    }
    // shift input samples back in time for next time
    memmove( &insamp[0], &insamp[length],
            (filterLength - 1) * sizeof(double) );
 
}
 
/*
% All frequency values are in Hz.
Fs = 8000;  % Sampling Frequency

N = 64;       % Order
Fc1 = 900;      % First Cutoff Frequency
Fc2 = 1100;     % Second Cutoff Frequency
flag = 'scale';  % Sampling Flag
% Create the window vector for the design algorithm.
win = hamming(N + 1);

% Calculate the coefficients using the FIR1 function.
b = fir1(N, [Fc1 Fc2] / (Fs / 2), 'bandpass', win, flag);
Hd = dfilt.dffir(b);

*/
#define FILTER_LEN  64
double coeffs[ FILTER_LEN ] =
{
	0.001271415,0.000613399,-0.000743051,-0.002237912,-0.002839066,-0.001499995,0.001907744,0.005809714,0.007246582,0.003696763,-0.004493201,-0.013012618,-0.015411974,-0.007469019,0.008637426,0.023849303,0.026991773,0.01252796,-0.013905401,-0.036925862,-0.040266172,-0.018036982,0.019350233,0.049730696,0.052544774,0.022829115,-0.023775485,-0.05936261,-0.060973286,-0.02576614,0.02611099,0.063457233,0.063457233,0.02611099,-0.02576614,-0.060973286,-0.05936261,-0.023775485,0.022829115,0.052544774,0.049730696,0.019350233,-0.018036982,-0.040266172,-0.036925862,-0.013905401,0.01252796,0.026991773,0.023849303,0.008637426,-0.007469019,-0.015411974,-0.013012618,-0.004493201,0.003696763,0.007246582,0.005809714,0.001907744,-0.001499995,-0.002839066,-0.002237912,-0.000743051,0.000613399,0.001271415
};
 
// number of samples to read per loop
#define SAMPLES   81

/*
	@signal: the input signal that needs to be filtered
	@signal_length: the length of the input and output signals
	@output_signal: memory has been allocated before calling this function
		NOTE: output_signal could be the same as signal (for in-place filtering)
	
*/ 
void ApplyFIR( double * signal, unsigned int signal_length, double * output_signal )
{
    int size;
    double floatInput[SAMPLES];
    double floatOutput[SAMPLES];
 
    // initialize the filter
    firFloatInit();
 
    // process all of the samples
	unsigned int i;
	for(i = 0; i < signal_length; i+= SAMPLES)
	{
		size = (signal_length - i < SAMPLES) ? signal_length - i : SAMPLES;
		memcpy(floatInput, & signal[i], size * sizeof(double));
				
		// perform the filtering
        firFloat( coeffs, floatInput, floatOutput, size, FILTER_LEN);
        
		//Place filtered signal in the output signal
		memcpy(& output_signal[i], floatOutput, size * sizeof(double));
	}
	
}