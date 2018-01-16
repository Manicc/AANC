/*
Capturing from ALSA
From Paul David's tutorial : http://equalarea.com/paul/alsa-audio.html

sudo apt-get install libasound2-dev
gcc -o alsa-record-example -lasound alsa-record-example.c && ./alsa-record-example hw:0

GPU_FFT thanks to Andrew Holme

other chunks of ccode copied wholesale from PeterO of raspberrypi.org, many thanks!
*/

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
//#include "tone.h"

void writeToFile(char* reading);
void spectrumAnalyzer();
int fft_compute_forward(float * input, int log2_N, float * output, double sampling_interval);
void *playTone(void *f);

int peakFrequencies[32768];
float peakAndIntensity[32768][2]; // Can be optimized ?
float playedPeakAndIntensity[32768][2];
float BPF[32768][3];

int peakCount = 0;
int frequency = 100;
char* speakerID = "hw:0,0";
char* micID = "hw:1,0";






main(int argc, char *argv[])
{
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
	unsigned long N = 1 << NLOG2;
	printf("N = %6.4d\n", N);
	FFTData = malloc(N * sizeof(float));	
	
	static struct timeval tm1,tm2;
	float freqResolution = recordRate / recordBufferLenght;
	int freq = 0;
	int i,j;
	int counter = 0;
	int peakFreqIndex[recordBufferLenght];
	//int intensityPeakFrequencies[recordBufferLenght];
	
	int internalPeakCount = 0;
	
	
	//@TODO Reallocate BPF initialization
	BPF[220][0] = 1;
	BPF[500][0] = 1;
	

	
	
	while(1){
	
		internalPeakCount = 0;
		gettimeofday(&tm1,NULL);
		
		//read from audio device into buffer
		snd_pcm_readi(capture_handle, recordBuffer, recordBufferLenght);

		for (i = 0; i < recordBufferLenght; i++){
			fRecordBuffer[i] = (float)recordBuffer[i] * (0.54 - (0.46 * (cos(2*M_PI*(i/(recordBufferLenght - 1)))))); // Convert to float, apply hammin window
			peakFreqIndex[i] = 0; // Initialize to 0
			peakFrequencies[i] = 0;
			
		}

	
		fft_compute_forward(fRecordBuffer, NLOG2, FFTData, 1);	
		
		// Calculate peak Hz
		int peakIndex = 0;
		float peakValue = 0;
		
		char reading[100];
		
		for(i = 0; i <= recordBufferLenght/2; i ++){
			if(FFTData[i] > peakValue){
				peakValue = FFTData[i];
				peakIndex = i;
			}
		}
		
		// Peak finding test
		int intensityIndex = 0;
		for(i = 1; i<(recordBufferLenght)/2; i ++){ // index BUG?
			if(FFTData[i] > FFTData[i-1] && FFTData[i] > FFTData[i+1]){
				peakFreqIndex[intensityIndex] = i;
				intensityIndex++;
				peakAndIntensity[0][i] = (((i * freqResolution) * 1.3459689603076) - 2);
				peakAndIntensity[1][i] = (FFTData[i] / peakValue) * 10000;
				if(peakAndIntensity[0][i] > 100 && peakAndIntensity[1][i] > 500){
					printf("Peak found. [%d]: F: %f, Intensity: %f RAW DATA: %d\n",i, peakAndIntensity[0][i] , peakAndIntensity[1][i], FFTData[i]);
					playedPeakAndIntensity[0][internalPeakCount] = peakAndIntensity[0][i];
					playedPeakAndIntensity[1][internalPeakCount] = peakAndIntensity[1][i];
					internalPeakCount++;
				}
			}
		}
		peakCount = internalPeakCount;


		gettimeofday(&tm2,NULL);
		unsigned long long t = (tm2.tv_usec - tm1.tv_usec);
		//printf("FFT TIME: %llu ms\n\n",t/1000);
}
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

		static struct timeval tm1,tm2;
		
		static int internalPhase[12000];
		
		int phase;
		
		
		
		
		while(1){
		
		//printf("OUTPUT frequency: %d\n", frequency);
		
		// Generate a sine wave
		float volume = 0.1;
		
		
		// Clear
		for (i = 0; i < playbackBufferLenght; i++){
			playbackBuffer[i] = 0;
		}
		printf("PeakCount: %d\n", peakCount);
			
		for (i = 0; i < peakCount; i++)
		{
			int freq = (int)floor(playedPeakAndIntensity[0][i]);
			printf("F: %d\n", freq);
			if( (BPF[freq][0]) != 1){
				break;
			}
			//printf("F: %d   PHASE: %d\n", (int)floor(playedPeakAndIntensity[0][i]), phase);
			phase = internalPhase[(int)floor(playedPeakAndIntensity[0][i])];
			
			for (j = 0; j < playbackBufferLenght; j++)
			{	
				playbackBuffer[j] += htole16(double_to_int16_t(volume * ((double)  (sin(playedPeakAndIntensity[0][i] * (2 * M_PI) * phase / playbackRate)) )));
				phase++;
				internalPhase[(int)floor(playedPeakAndIntensity[0][i])] = phase;
			}
		}
		

		
		//printf("TONE FREQUENCY: %d",peakFrequencies[0]);
		for (i = 0; i < 1; i++) {
			gettimeofday(&tm1,NULL);
			
			snd_pcm_writei(phandle, playbackBuffer, playbackBufferLenght);
			
			
			gettimeofday(&tm2,NULL);
			unsigned long long t = (tm2.tv_usec - tm1.tv_usec);
			//printf("TONE GENERATOR TIME: %llu ms\n",t/1000);
			printf("\n");
			
		}

}
	snd_pcm_close(phandle);

    return 0;
}

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


