#include <alsa/asoundlib.h>

int main(int argc, char *argv[]){
	
	int16_t *recordBuffer;
	int recordBufferLenght = 2205;
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
		
		int lines = (decibel * 20) / 60;
		
		printf("Decibel: %f", decibel);
		for(i = 0; i <= lines; i++){
			printf("-");
		}
		printf("\n");
		//printf("Total Amplitude: %f\n", totalAmplitude);
		//printf("Decibel: %f\n", decibel);
	}
	
	snd_pcm_close(capture_handle);
}
