all:
	gcc main.c -g src/SDL3_gfx/*.c src/*.c -I src/include -I src/SDL3_gfx -L src/lib -o main -lSDL3 -lSDL3_ttf -lm