
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <limits.h>

#include <pthread.h>
#include <alsa/asoundlib.h>

#include "mailbox.h"
#include "gpu_fft.h"

void writeToFile(char* reading);
void spectrumAnalyzer();
int fft_compute_forward(float * input, int log2_N, float * output, double sampling_interval);
void *playTone(void *f);


float peakAndIntensity[32768][2]; // Can be optimized ?
float foundPeakFrequency[32768][2];
float BPF[12000][2]; // 12Khz max freq
float sourceDistance = 2.0;
#define airTemperature 20
float soundVelocity = 331 + (0.6 * airTemperature); // v = 331m/s + 0.6m/s/C * T

int globalPeakCount = 0;
int frequency = 100;
char* speakerID = "hw:0,0";
char* micID = "hw:1,0";






main(int argc, char *argv[])
{
	// pthread is a multithreading library
	// In this case playTone() function is runned in a seperate thread
	pthread_t toneThread;
	
	pthread_create(&toneThread,NULL,playTone,(void *)frequency);
	spectrumAnalyzer();
	
	pthread_exit(NULL);
}


void spectrumAnalyzer(){
	
	int16_t *recordBuffer;
	float* fRecordBuffer;
	int recordBufferLenght = 32768;
	unsigned int recordRate = 44100;
	
	int err;

	// Initialize audio in
	snd_pcm_t *capture_handle;
	snd_pcm_hw_params_t *hw_params;
	snd_pcm_format_t format = SND_PCM_FORMAT_S16_LE;
	
	snd_pcm_open(&capture_handle, micID, SND_PCM_STREAM_CAPTURE, 0);
	snd_pcm_hw_params_malloc(&hw_params);
	snd_pcm_hw_params_any(capture_handle, hw_params);
	snd_pcm_hw_params_set_access(capture_handle, hw_params, SND_PCM_ACCESS_RW_INTERLEAVED);
	snd_pcm_hw_params_set_format(capture_handle, hw_params, format);
	snd_pcm_hw_params_set_rate_near(capture_handle, hw_params, &recordRate, 0);
	snd_pcm_hw_params_set_channels(capture_handle, hw_params, 1);
	snd_pcm_hw_params(capture_handle, hw_params);

	snd_pcm_hw_params_free(hw_params);
	snd_pcm_prepare(capture_handle);

	//allocate buffer of 16bit ints, as specified in PCM_FORMAT
	recordBuffer = malloc(recordBufferLenght * snd_pcm_format_width(format) / 8 * 2); // 
	fRecordBuffer = malloc(recordBufferLenght*sizeof(float));
	fprintf(stdout, "recordBuffer allocated\n");
		
	float *FFTData;	
	int NLOG2 = 15;
	unsigned long N = 1 << NLOG2; // Used to determine maximum frequency count of FFT, "Bitwise left operation"
	FFTData = malloc(N * sizeof(float));	
	float freqResolution = recordRate / recordBufferLenght;
	int i,j;
	
	// Array to store all peak frequencies
	int peakFreqIndex[recordBufferLenght];
	
	int internalPeakCount = 0;
	
	// Phase calculation in angles for each frequency from 1 Hz to 12KHz
	for(i = 1; i <= 12000; i++){
		
		/* 
		i = f
		T = 1 / i
		Phase = 180 - ( ( (1-(d/vT)) / T) * 360)
		*/
		
		double numerator = (1 - (sourceDistance/(soundVelocity * 1/i)));
		
		BPF[i][1] = fmod(180 - ( (numerator/ (1.00/i)) * 360),360.0);
		printf("Phase: %.10f\n", BPF[i][1]);
	}
	
	//@TODO Reallocate BPF initialization
	BPF[200][0] = 1;
	BPF[220][0] = 1;
	BPF[300][0] = 1;
	BPF[500][0] = 1;
	

	
	while(1){
	
		internalPeakCount = 0;
		
		//Record audio to recordBuffer
		snd_pcm_readi(capture_handle, recordBuffer, recordBufferLenght);

		for (i = 0; i < recordBufferLenght; i++){
			fRecordBuffer[i] = (float)recordBuffer[i] * (0.54 - (0.46 * (cos(2*M_PI*(i/(recordBufferLenght - 1)))))); // Convert to float, apply hammin window			
		}
		
		// Apply FFT to fRecordBuffer
		fft_compute_forward(fRecordBuffer, NLOG2, FFTData, 1);	
		fprintf(stdout, "asd\n");
		// Find frequency with the highest intensity(Amplitude)
		int peakIntensityIndex = 0;
		float peakIntensity = 0;
		
		for(i = 0; i <= recordBufferLenght/2; i ++){
			if(FFTData[i] > peakIntensity){
				peakIntensity = FFTData[i];
				peakIntensityIndex = i;
			}
		}

		// Peak finding
		// This algorithm tries to find all frequency peaks based on their intensities(Amplitudes)

		for(i = 1; i< recordBufferLenght/2; i ++){ // @TODO index BUG?
			if(FFTData[i] > FFTData[i-1] && FFTData[i] > FFTData[i+1]){
				peakAndIntensity[0][i] = (((i * freqResolution) * 1.3459689603076) - 2);
				peakAndIntensity[1][i] = (FFTData[i] / peakIntensity) * 10000;
				if(peakAndIntensity[0][i] > 100 && peakAndIntensity[1][i] > 500){
					printf("Peak found. [%d]: F: %f, Intensity: %f RAW DATA: %d\n",i, peakAndIntensity[0][i] , peakAndIntensity[1][i], FFTData[i]);
					foundPeakFrequency[0][internalPeakCount] = peakAndIntensity[0][i];
					foundPeakFrequency[1][internalPeakCount] = peakAndIntensity[1][i];
					internalPeakCount++;
				}
			}
		}
		
		globalPeakCount = internalPeakCount;
}// End of while loop


	snd_pcm_close(capture_handle);
}

int fft_compute_forward(float * input, int log2_N, float * output, double sampling_interval)
{
	int mb = mbox_open();
	int jobs = 1, i;
	unsigned long N = 1 << log2_N;

	float FloatN = (float)N;
	float HalfN = FloatN / 2.0;
	struct GPU_FFT_COMPLEX * base;
	struct GPU_FFT * fft;
	struct GPU_FFT_COMPLEX *DataIn, *DataOut;
	int ret = gpu_fft_prepare(mb, log2_N, GPU_FFT_FWD, jobs, &fft);
	base = fft->in;
	
	for (i = 0; i<N; i++)
	{
		base[i].re = input[i];
		base[i].im = 0.0;
	}

	gpu_fft_execute(fft);
	base = fft->out;

	for (i = 0; i<N; i++){
		output[i] = (base[i].re * base[i].re) + (base[i].im * base[i].im);
	}
	gpu_fft_release(fft);
	mbox_close(mb);
}



// assume -1.0 <= x <= 1.0
int16_t double_to_int16_t (double x)
{
    return trunc (x * 32767); // rounds toward zero
}


void *playTone(void *f)
{
		int playbackRate = 44100;
		int playbackBufferLenght = 11025;
		
		uint16_t playbackBuffer[playbackBufferLenght];
		int i;
		int j;

		
		snd_pcm_t *phandle;
		snd_pcm_hw_params_t *hw_params;

		snd_pcm_open(&phandle, speakerID, SND_PCM_STREAM_PLAYBACK, 0);
		snd_pcm_hw_params_malloc(&hw_params);
		snd_pcm_hw_params_any(phandle, hw_params);
		snd_pcm_hw_params_set_access(phandle, hw_params, SND_PCM_ACCESS_RW_INTERLEAVED);
		snd_pcm_hw_params_set_format(phandle, hw_params, SND_PCM_FORMAT_S16_LE);
		snd_pcm_hw_params_set_rate_near(phandle, hw_params, &playbackRate, 0);
		snd_pcm_hw_params_set_channels(phandle, hw_params, 1);
		snd_pcm_hw_params(phandle, hw_params);
		snd_pcm_hw_params_free(hw_params);
		snd_pcm_prepare(phandle);

		static int internalPhase[12000];
		
		int phase;
		
		
		
		
		while(1){
		
		// Generate a sine wave
		float volume = 0.1;
		
		
		// Clear playbackBuffer
		for (i = 0; i < playbackBufferLenght; i++){
			playbackBuffer[i] = 0;
		}
		
			
		for (i = 0; i < globalPeakCount; i++){
			// Here we check wheter a certain peak frequency is listed in the Bandwidth Pass Filter
			// BPF 2D-array is structured so, that each row's index represents a corresponding frequency 
			// First column (index 0) contains numeric value 0 or 1. If value is set to 1, frequency is allowed to pass through.
			// e.g. BPF[100][0] = 1; 100Hz frequency is passed 
			//		BPF[123][0] = 0; 123Hz frequency is not passed
			
			//Here we check if BPF contains 1 or 0
			// (int)floor(foundPeakFrequency[0][i]) returns a found frequency.
			// if 1 is not found, For look runs to break; and found frequency peak is not added to the output sound signal
			if( (BPF[ (int)floor(foundPeakFrequency[0][i])][0]) != 1){
				break;
			}
			
			phase = internalPhase[(int)floor(foundPeakFrequency[0][i])];
			
			for (j = 0; j < playbackBufferLenght; j++)
			{	
				playbackBuffer[j] += htole16(double_to_int16_t(volume * ((double)  (sin(foundPeakFrequency[0][i] * (2 * M_PI) * phase / playbackRate)) )));
				phase++;
				internalPhase[(int)floor(foundPeakFrequency[0][i])] = phase;
			}
		}
		
		//Write playbackBuffer to audio driver (Play generated tone)
		for (i = 0; i < 1; i++) {
			snd_pcm_writei(phandle, playbackBuffer, playbackBufferLenght);
		}

}
	snd_pcm_close(phandle);

    return 0;
}

// function for saving values to .txt file
void writeToFile(char* reading){

	char *filePath = "/home/pi/Documents/AANC/test.txt";
	FILE *file = fopen(filePath, "a"); 
	
	if (file != NULL)
	{
		
		fputs(reading,file);	
		fputs("\n",file);
		
		fclose(file);
	}
}


