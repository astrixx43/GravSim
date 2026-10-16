all:
	gcc bodies.c -c bodies.o
	gcc Wavelength_Rgb.c -c Wavelength_Rgb.o
	gcc main.c -g src/SDL3_gfx/*.c -I src/include -I src/SDL3_gfx -L src/lib -o main -lSDL3 *.o -lm