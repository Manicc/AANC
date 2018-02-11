



main(int argc, char *argv[])
{
	// pthread is a multithreading library
	// In this case playTone() function is runned in a seperate thread
	pthread_t toneThread;
	
	pthread_create(&toneThread,NULL,playTone,(void *)frequency);
	spectrumAnalyzer();
	
	pthread_exit(NULL);
}
