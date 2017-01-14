//
//Definitions: 
	//Strategy document: https://docs.google.com/document/d/10Bu1yvGFMuSSRcf0y2kXYmD8a38s4vKvYzKdGpoMt9Q/edit
	//MOCK: Simplifies a hardware/embedded problem to the best that can be done to model it using command line sw.
	//STUB: can't do anything in the command line sw to actually represent this action.
/*Pending TODOs:
	-remove mallocs and statically allocate and reuse via memset
	-change doubles to floats
	
	-think about how the results will be passed to the antivoice algorithm	
		revise these these 2 points:
		-it should be passed comma separated by the debug monitor to the computer
		-then with the calibraiton results, the antivoice algorithm is rebuilt and flashed
	
	-do not need to store test_signal as you can calculate each index
	-find a way to perform obs signal analysis without storing the entire signal in memory (lower memory consumption)
		-perhaps send it via the serial debug monitor to the pc?
*/

#include "CalibrationSetup.h"


/*
	Helper function for getting frequency index from Hertz
	params:
		f is the frequency in Hertz
	returns:
		i, index into frequencies array
*/
unsigned int GetFrequenciesArrayIndex(double f)
{
	unsigned int i = 0; //i is the index
	while (i < freq_count_pos)
	{
		//incase there is a small round error see if approximately the same
		double abs_difference = MAX(f,frequencies[i]) - MIN(f,frequencies[i]);
		
		if (f == frequencies[i] || abs_difference < 0.01)
		{
			break;
		}
		i++;
	}
	return i;
}

/*
	params: 
		f: frequency
		fs: sample rate for DAC, thus signal resolution
		ampltiude: amplitude of sinusoid
		signal_array_length: length of test signal to be generated
	returns:
		an array of points where each point represents the time domain height of the signal
*/
double * TestSignalGenerator(double f, unsigned int fs, double amplitude, unsigned int signal_array_length)
{
	double * test_signal = (double *)malloc(signal_array_length * sizeof(double)); //allocate memory

	unsigned int i; //the sample number
	for (i = 0; i < signal_array_length; i++)
	{
		//create the test signal point for each sample i 
		test_signal[i] = amplitude * cos(2 * PI * f * (i / (double) fs));
	}

	return test_signal;
}

/*
	Helper function
	Apply @num_samples delay to @signal in-place
*/
void ApplySampleDelay(double * signal, unsigned int num_samples, unsigned int signal_array_length)
{
	if (num_samples < signal_array_length)
	{
		//Apply specified delay
		unsigned int delay_bytes = (signal_array_length - num_samples) * sizeof(double);
		memcpy(&signal[num_samples], signal, delay_bytes);

		//zero out the delayed amount
		memset(signal, 0, num_samples * sizeof(double));
	}
	else
	{
		//zero out the entire signal and show warning
		memset(signal, 0, signal_array_length * sizeof(double));
		printf("\n ****** \n [WARNING] delay is greater than signal length (duration) \n ****** \n");
	}
}

/*
	Applies magnitude and group delay transfer function to test signal at frequency f
	parameters:
		f: frequency in Hz
		test_signal: array of points, with each point representing signal height at that sample number
	returns:
		changed_signal, with transfer function applied
*/
double * ApplyMicTransducerTransferFunction(double f, unsigned int fs, double * test_signal, unsigned int signal_array_length)
{
	unsigned int i;

	double * changed_signal = (double *)malloc(signal_array_length * sizeof(double));
	memcpy(changed_signal, test_signal, signal_array_length * sizeof(double)); //perform deep copy of test_signal to changed_signal

	//Get gain and delay amounts for this frequency f
	double gain_amount = system_gain[GetFrequenciesArrayIndex(f)];
	double seconds_delay = system_delay[GetFrequenciesArrayIndex(f)];
	unsigned int samples_delay = (unsigned int)round(seconds_delay * fs);

	for (i = 0; i < signal_array_length; i++)
	{
		//Apply gain
		changed_signal[i] *= gain_amount;
	}

	ApplySampleDelay(changed_signal, samples_delay, signal_array_length);
	
	return changed_signal;
}

/* 
	Generates additive white Gaussian Noise samples with zero mean and a standard deviation of @d. 
	This method is invoked once per sample, and the value it returns should be added to the signal.
	Author: https://www.embeddedrelated.com/showcode/311.php
*/
double AWGN_generator(double d)
{
	double temp1;
	double temp2;
	double result;
	int p;

	p = 1;

	while (p > 0)
	{
		//rand() function generates an integer between 0 and RAND_MAX,
			//which is defined in stdlib.h.
		temp2 = (rand() / ((double)RAND_MAX));

		if (temp2 == 0)
		{// temp2 is >= (RAND_MAX / 2)
			p = 1;
		}// end if
		else
		{// temp2 is < (RAND_MAX / 2)
			p = -1;
		}// end else

	}// end while()

	temp1 = cos((2.0 * (double)PI) * rand() / ((double)RAND_MAX));
	result = sqrt(-2.0 * log(temp2)) * temp1;

	return (result * d);
}

/*
	Applies channel noise to signal and air time
parameters:
	in_signal: array of points, withe each point representing signal height at that sample number
	ch_delay: time in seconds the channel delays by
returns:
	out_signal, with channel noise and delay applied
*/
double * ApplyChannel(double * in_signal, unsigned int samples_delay, unsigned int fs, unsigned int signal_array_length)
{
	unsigned int i;

	double * out_signal = in_signal; //no deep copy required since in_signal is intermediate and can be readily modified

	//Apply channel delay
	ApplySampleDelay(out_signal, samples_delay, signal_array_length);

	//Apply channel noise 
	for (i = 0; i < signal_array_length; i++)
	{
		out_signal[i] += AWGN_generator(WHITE_NOISE_STD_DEVIATION);
	}

	return out_signal;
}

/*
	Returns the amplitude of exactly periodic cosine single-tone signal, but amplitude is affected by noise.  
	params:
		h: signal level at each sample
		f: the frequency in Hz
		fs: sample frequency in Hz
	returns:
		{amplitude of sinusoid, signal delay}
*/
AMPDELAY_t DetectSinusoidalAmplitudeAndDelay(double * h, const double f, const unsigned int fs, const unsigned int signal_array_length)
{
	AMPDELAY_t ret = { 0, 0};
	const int N = (unsigned int)ceil((double)fs / f); //Samples in period

	if (signal_array_length < N)
	{
		printf("\n ****** \n [ERROR] signal length (duration) is less than period! \n ****** \n");
		return ret;
	}

	//Get the sums of the amplitude to accurately determine the signal offset to the first cosine max
	double * amplitudeSums = (double *)malloc(sizeof(double) * N);
	double * sumsCount = (double *)malloc(sizeof(double) * N);
	for (int i = 0; i < signal_array_length; i++)
	{
		if (i < N)
		{
			amplitudeSums[i] = h[i];
			sumsCount[i] = 1;
		}
		else
		{
			amplitudeSums[i % N] += h[i];
			sumsCount[i % N]++;
		}
	}

	//determine the index at which max of sinusoids occur across 1 cycle (% N)
	unsigned int minindex = 0;
	unsigned int maxindex = 0;
	for (int i = 1; i < N; i++)
	{
		if (amplitudeSums[i] > amplitudeSums[maxindex])
		{
			maxindex = i;
		}
	}

	unsigned int lastmaxima;
	for (int i = maxindex; i < signal_array_length; i+=N)
	{
		lastmaxima = i;
	}

	//Note: due to delay in signal, the first several cycles may be noise and no signal
	//Identify the number of samples delayed for start of signal
	for (int i = maxindex ; i < signal_array_length; i+= N)
	{
		//If at least half of the last maxima
		if (h[i] > 0.5 * h[lastmaxima])
		{
			//We found the first true signal maxima
			maxindex = i;
			break;
		}
		else
		{
			amplitudeSums[maxindex] -= h[i];
			sumsCount[maxindex]--;
		}
	}

	//now we do the same for corresponding minimas
	double minsum = 0;
	double minsumcount = 0;
	for (int i = maxindex + N/2; i < signal_array_length; i += N)
	{
		minsum += h[i];
		minsumcount++;
	}

	unsigned int delay_samples = maxindex;
	double delay_seconds = (double)delay_samples / (double)fs;
	ret.delay = delay_seconds;

	//Find average of min and max sums, to calculate amplitude
	double minavg = minsum / minsumcount;
	double maxavg = amplitudeSums[maxindex % N] / sumsCount[maxindex % N];
	ret.amplitude = (fabs(maxavg) + fabs(minavg)) / 2; //fabs is abs of floating number

	free(amplitudeSums);
	free(sumsCount);

	return ret;
}

/*
	Performs calibration routine
	parameters:
		f: the frequency at which to calibrate in Hertz
		test_amplitude: the amplitude of the signal being generated
	returns:
		CALIB_t containing calibration results
*/
CALIB_t GetCalibration(double f, double test_amplitude, unsigned int fs, unsigned int channel_delay_samples)
{
	//Generate test signal
	

	unsigned int signal_array_length = (unsigned int)ceil(duration * fs); //total number of samples

	double * test_signal = TestSignalGenerator(f, fs, test_amplitude, signal_array_length);

	//
	//STUB: Use DAC to send to Speaker driver
	//

	//MOCK: We model the time it takes for sound to travel from transducer to observation microphone
		//v = d/t => t = d/v;  
	
	//
	//STUB: ADC of observation microphone to MCU
	//

	//MOCKS (transducer and mic transfer functions unified for this simulation):
		//MOCK: transducer magnitude response and group delay
		//MOCK: observation magnitude response and group delay
		//transducer and obs mic are bundled together as 1 transfer function
	
	double * changed_signal = ApplyMicTransducerTransferFunction(f, fs, test_signal, signal_array_length);


	//MOCK: Channel (air / thermal) noise, model as gaussian noise
	double * obs_signal = ApplyChannel(changed_signal, channel_delay_samples, fs, signal_array_length);

	/*Not performing filtering since accurate enough when signal contains over 50 periods, assuming in a quiet room*/
		//Analyze the obs signal relative to test signal to generate calibration outputs needed for antivoice algorithm
			//filter noise from obs signal by using BPF centered at f
		//double * filtered_signal = (double *)malloc(signal_array_length * sizeof(double)); //allocate memory
		//ApplyFIR(obs_signal, signal_array_length, filtered_signal);
	//TODO: terrible assumption.  Use FFT to remove all other frequencies not of interest so that no need to hardcode linear filters
	
	//determine amplitude and seconds delay of signal relative to test
	AMPDELAY_t obs_ampdelay = DetectSinusoidalAmplitudeAndDelay(obs_signal, f, fs, signal_array_length);
	
	//Determine results
	CALIB_t result;
	result.frequency = f;
	result.amplitude = test_amplitude;
	result.eqgain = test_amplitude / obs_ampdelay.amplitude;
	result.rtt = obs_ampdelay.delay;
	result.internaltime = result.rtt - channel_delay_seconds;

	//free memory
	free(test_signal);
	free(obs_signal);
	//free(filtered_signal);
	return result;
}		

void PrintCalib_t(CALIB_t c)
{
	printf("%.2f | %.3f      | %.5f | %.5f | %.5f | %.9f | %.9f |\n", 
		c.frequency, c.amplitude, 
		system_gain[GetFrequenciesArrayIndex(c.frequency)], system_delay[GetFrequenciesArrayIndex(c.frequency)],
		c.eqgain, c.rtt, c.internaltime);
}

int main(int argc, char ** argv)
{
	//TODO: make frequencies be passed as an array by command line?

	
	//Calculate channel delay in terms of samples
	unsigned int channel_delay_samples = (unsigned int)round(channel_delay_seconds * fs);
	printf("\n[Warning] channel delay of %.6f seconds rounded to %.6f seconds due to nearest sample.",
		channel_delay_seconds, channel_delay_samples / (double) fs);
	channel_delay_seconds = channel_delay_samples / (double)fs;
	
	CALIB_t calib_results[freq_count_pos * amplitude_count]; //stores the calibration results

	printf("\nFreq | Amplitude  | Sysgain | SysDelay| eqGain  |  RTT        | Int. Delay  | \n");

	//iterate over all frequencies and collect results for each frequency
	unsigned int i,j;
	for (i = 0; i < freq_count_pos; i++)
	{
		for (j = 0; j < amplitude_count; j++)
		{
			calib_results[i + j] = GetCalibration(frequencies[i], amplitudes[j], fs, channel_delay_samples);
			PrintCalib_t(calib_results[i + j]);
		}
	}

    return 0;
}

