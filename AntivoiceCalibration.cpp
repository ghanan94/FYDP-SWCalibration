//
//Definitions: 
	//Strategy document: https://docs.google.com/document/d/10Bu1yvGFMuSSRcf0y2kXYmD8a38s4vKvYzKdGpoMt9Q/edit
	//MOCK: Simplifies a hardware/embedded problem to the best that can be done to model it using command line sw.
	//STUB: can't do anything in the command line sw to actually represent this action.
/*Pending TODOs:
	-remove mallocs and statically allocate and reuse via memset
	-change doubles to floats
	-add #ifdef SIM
	
	-how we going to have modes of operation be able to change them at bootup
	-think about how calibration data is fed into antivoice algorithm	
		-then with the calibraiton results, the antivoice algorithm is rebuilt and flashed
	
	-do not need to store test_signal as you can calculate each index
	-find a way to perform obs signal analysis without storing the entire signal in memory (lower memory consumption)
		-perhaps send it via the serial debug monitor to the pc?

	- make a realtime version of this program
*/

#include "calibration_setup.h"
#include "physical_mocking.h"


/*
	params: 
		f: frequency
		fs: sample rate for DAC, thus signal resolution
		ampltiude: amplitude of sinusoid
	returns:
		an array of points where each point represents the time domain height of the signal
*/
double * TestSignalGenerator(double f, unsigned int fs, double amplitude)
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

void playTestToneSample(double f, unsigned int fs, double amplitude, unsigned int sample_number)
{
	//create the test signal point for each sample i 
	double output_sample = amplitude * cos(2 * PI * f * (sample_number / (double)fs));

#ifdef SIM
	
#else 

#endif
	//TODO: send to speaker... (there is only 1 in this case)
		//in simulation actually send it a statically allocated array
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
AMPDELAY_t DetectSinusoidalAmplitudeAndDelay(double * h, const double f, const unsigned int fs)
{
	AMPDELAY_t ret = { 0, 0};
	const unsigned int N = (unsigned int)ceil((double)fs / f); //Samples in period

	if (signal_array_length < N)
	{
		printf("\n ****** \n [ERROR] signal length (duration) is less than period! \n ****** \n");
		return ret;
	}

	//Get the sums of the amplitude to accurately determine the signal offset to the first cosine max
	double * amplitudeSums = (double *)malloc(sizeof(double) * N);
	double * sumsCount = (double *)malloc(sizeof(double) * N);
	for (unsigned int i = 0; i < signal_array_length; i++)
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
	for (unsigned int i = 1; i < N; i++)
	{
		if (amplitudeSums[i] > amplitudeSums[maxindex])
		{
			maxindex = i;
		}
	}

	unsigned int lastmaxima;
	for (unsigned int i = maxindex; i < signal_array_length; i+=N)
	{
		lastmaxima = i;
	}

	//Note: due to delay in signal, the first several cycles may be noise and no signal
	//Identify the number of samples delayed for start of signal
	for (unsigned int i = maxindex ; i < signal_array_length; i+= N)
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
	for (unsigned int i = maxindex + N/2; i < signal_array_length; i += N)
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
CALIB_t GetCalibration(double f, double test_amplitude, unsigned int fs)
{
	//Generate test signal
	

	

	double * test_signal = TestSignalGenerator(f, fs, test_amplitude);

	//Need to alterate between sending and receiving...?
	unsigned int samples_sent = 0;
	while (samples_sent < signal_array_length)
	{
		//send tone to speaker
		playTestToneSample(f, fs, test_amplitude, samples_sent);

	}

	//instead of generating entire array of test signal ... we should have a function with a static variable i
		//it should create one audio sample at a time
		//there should be playTone function which in simulation simply copies it to a statically allocated global array
			//which then get mic sample gets it from that array..

	//MOCK: We model the time it takes for sound to travel from transducer to observation microphone
		//v = d/t => t = d/v;  
	
	//
	//STUB: ADC of observation microphone to MCU
	//

	//MOCKS (transducer and mic transfer functions unified for this simulation):
		//MOCK: transducer magnitude response and group delay
		//MOCK: observation magnitude response and group delay
		//transducer and obs mic are bundled together as 1 transfer function
	
	//TODO: move this to get mic sample
	double * changed_signal = ApplyMicTransducerTransferFunction(f, fs, test_signal);

	//TODO move this to get mic sample
	//MOCK: Channel (air / thermal) noise, model as gaussian noise
	double * obs_signal = ApplyChannel(changed_signal, fs);

	/*Not performing filtering since accurate enough when signal contains over 50 periods, assuming in a quiet room*/
		//Analyze the obs signal relative to test signal to generate calibration outputs needed for antivoice algorithm
			//filter noise from obs signal by using BPF centered at f
		//double * filtered_signal = (double *)malloc(signal_array_length * sizeof(double)); //allocate memory
		//ApplyFIR(obs_signal, signal_array_length, filtered_signal);
	//TODO: terrible assumption.  Use FFT to remove all other frequencies not of interest so that no need to hardcode linear filters
	
	//determine amplitude and seconds delay of signal relative to test
	AMPDELAY_t obs_ampdelay = DetectSinusoidalAmplitudeAndDelay(obs_signal, f, fs);
	
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


int main(int argc, char ** argv)
{
	//TODO: make frequencies be passed as an array by command line?

	/* MOCKING */
	//Calculate channel delay in terms of samples
	unsigned int channel_delay_samples = (unsigned int)round(channel_delay_seconds * MIC_FS);
	printf("\n[Info] channel delay of %.6f seconds rounded to %.6f seconds due to nearest sample.",
		channel_delay_seconds, channel_delay_samples / (double)MIC_FS);
	//Update the seconds value to the actual number of samples we are using in the system
	channel_delay_seconds = channel_delay_samples / (double)MIC_FS;
	/* END MOCKING */
	
	CALIB_t calib_results[freq_count_pos * amplitude_count]; //stores the calibration results

	printf("\nFreq | Amplitude  | Sysgain | SysDelay| eqGain  |  RTT        | Int. Delay  | \n");

	//iterate over all frequencies and collect results for each frequency
	unsigned int i,j;
	for (i = 0; i < freq_count_pos; i++)
	{
		for (j = 0; j < amplitude_count; j++)
		{
			calib_results[i + j] = GetCalibration(frequencies[i], amplitudes[j], MIC_FS);

			//MOCK
			PrintCalib_t(calib_results[i + j]);
		}
	}

    return 0;
}

