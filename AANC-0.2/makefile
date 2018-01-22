C = sa.c mailbox.c gpu_fft.c gpu_fft_base.c gpu_fft_twiddles.c gpu_fft_shaders.c

H = gpu_fft.h mailbox.h 

F =   -lrt -lm -ldl -pthread -lasound

FILES = $(shell sudo amixer -c 0 cset numid=3 1)


SA: $(C) $(H)
	gcc -o SA $(F) $(C)

fix: 
	echo $(FILES)

