#include <alsa/asoundlib.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <float.h>


void *amplitude();
void *playTone();
void *timerToggle();

int frequency = 100;
char* speakerID = "hw:0,0";
char* micID = "hw:1,0";

float newDecibel = 0;

int timeToggle = 0;


main(int argc, char *argv[]){
	// pthread is a multithreading library
	// In this case playTone() function is runned in a seperate thread
	pthread_t amplitudeThread;
	pthread_t toneThread;
	pthread_t timerThread;
	
	pthread_create(&amplitudeThread,NULL,amplitude,(void *)frequency);
	pthread_create(&toneThread,NULL,playTone,(void *)frequency);
	//pthread_create(&timerThread,NULL,timerToggle,(void *)frequency);
	
	pthread_exit(NULL);
}

void *timerToggle(void *f){

	while(1){
		usleep(60000);
		timeToggle ^= 1;
		printf("TimeToggle: %d\n", timeToggle);
		
	}
}


// assume -1.0 <= x <= 1.0
int16_t double_to_int16_t (double x)
{
    return trunc (x * 32767); // rounds toward zero
}

void *playTone(void *f){
		int playbackRate = 44100;
		int playbackBufferLenght = 882;
		
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

		float freq = 300;
		float T = 1.0 / freq;
		float timeResolution = 1.0 / playbackRate;
		
		int steps = T / timeResolution;
		
		static float attuenationArray[147][1];
		int stepIndex = 0;
		
		int bufferIndexTimeOffset = 0;
		float speakerOffset = 0;
		static float oldDecibel = 0;
		int lowFound = 0;
		
	
		static int phase = 0;
		while(stepIndex < 148){
		
			bufferIndexTimeOffset = 0;
			
			// Generate a sine wave
			float volume = 0.1;
			
			// Clear playbackBuffer
			for (i = 0; i < playbackBufferLenght; i++){
				playbackBuffer[i] = 0;
			}
			

			
			for (j = 0; j < playbackBufferLenght; j++){	
				
				playbackBuffer[j] += htole16(double_to_int16_t(volume * ((double)  (sin((freq + speakerOffset) * (2 * M_PI) * phase / playbackRate)) )));
				phase++;
			}
			
			//Write playbackBuffer to audio driver (Play generated tone)
			for (i = 0; i < 1; i++) {
				snd_pcm_writei(phandle, playbackBuffer, playbackBufferLenght);
			}
			//printf("Phase: %d\n", phase);
			attuenationArray[stepIndex][0] = newDecibel;
			attuenationArray[stepIndex][1] = phase;
			stepIndex++;
			phase += 1;
	}
	float lowestDecibel = 200;
	for(i = 0; i < 147; i++){
		if(attuenationArray[i][0] < lowestDecibel && attuenationArray[i][0] > 10 ){
			lowestDecibel = attuenationArray[i][0];
			//phase = (i * 2) + 50;
			attuenationArray[stepIndex][1];
		}
	}
	printf("LowestdB: %f, Phase: %d \n", lowestDecibel, phase);
		
	
	while(1){
		
			bufferIndexTimeOffset = 0;
			
			// Generate a sine wave
			float volume = 0.1;
			
			// Clear playbackBuffer
			for (i = 0; i < playbackBufferLenght; i++){
				playbackBuffer[i] = 0;
			}
			

			
			for (j = 0; j < playbackBufferLenght; j++){	
				
				playbackBuffer[j] += htole16(double_to_int16_t(volume * ((double)  (sin((freq + speakerOffset) * (2 * M_PI) * phase / playbackRate)) )));
				phase++;

			}
			
			//Write playbackBuffer to audio driver (Play generated tone)
			for (i = 0; i < 1; i++) {
				snd_pcm_writei(phandle, playbackBuffer, playbackBufferLenght);
			}
			printf("Phase: %d\n", phase);
	}
	
	//for(i = 0; i < 147; i++){
		//printf("StepIndex: %d Decibel: %f\n",i, attuenationArray[i][0]);
	//}
	
	snd_pcm_close(phandle);

    return 0;
}

void *amplitude(void *f){
	
	int16_t *recordBuffer;
	int recordBufferLenght = 441;
	unsigned int recordRate = 44100;
	char* micID = "hw:1,0";
	
	int err;
	int i,j;

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
	
	fprintf(stdout, "recordBuffer allocated\n");

	while(1){
		
		//Record audio to recordBuffer
		snd_pcm_readi(capture_handle, recordBuffer, recordBufferLenght);
		
		int totalAmplitude = 0;
		for(i = 0; i < recordBufferLenght; i++){
			totalAmplitude += abs(recordBuffer[i]);
		}
		totalAmplitude = totalAmplitude / recordBufferLenght;
		
		float decibel = 20 * log10f(totalAmplitude);
		newDecibel = decibel;
		
		int lines = (decibel * 20) / 60;
		
		//printf("Decibel: %f", decibel);
		for(i = 0; i <= lines; i++){
		//	printf("-");
		}
		//printf("\n");
		//printf("Total Amplitude: %f\n", totalAmplitude);
		//printf("Decibel: %f\n", decibel);
	}
	
	snd_pcm_close(capture_handle);
}
