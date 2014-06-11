
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <SDL.h>

#include "analogtv.h"

static void sdl_panic()
{
	fprintf(stderr, "SDL: %s\n", SDL_GetError());
	exit(EXIT_FAILURE);
}


int main(int argc, char** argv)
{
	if(SDL_Init(SDL_INIT_VIDEO) != 0) sdl_panic();

	int width = 400 * 2;
	int height = 300 * 2;

	SDL_Window* sdl_window = SDL_CreateWindow(
		"PTKS",
		SDL_WINDOWPOS_UNDEFINED,
		SDL_WINDOWPOS_UNDEFINED,
		width,
		height,
		0);

	SDL_Renderer* sdl_renderer = SDL_CreateRenderer(sdl_window, -1, 0);

	SDL_Texture* sdl_texture = SDL_CreateTexture(
		sdl_renderer,
		SDL_PIXELFORMAT_ARGB8888,
		SDL_TEXTUREACCESS_STREAMING,
		width,
		height);

	analogtv* tv = analogtv_allocate(width, height);
	analogtv_input* tvi = analogtv_input_allocate();

	analogtv_reception tvr;
	bzero(&tvr, sizeof(analogtv_reception));
	tvr.input = tvi;
	tvr.level = 1.40;
	tvr.multipath = 0.1;
	tvr.hfloss = 0;

	int desync = 0;
	int exiting = 0;
	int frame = 0;
	int sync = 1;
	while (!exiting) {

		SDL_Event e;
		while (SDL_PollEvent(&e)) {
			if (e.type == SDL_QUIT) {
				exiting = 1;
			}
			if (e.type == SDL_KEYDOWN) {
				if (e.key.keysym.sym == SDLK_ESCAPE || e.key.keysym.sym == SDLK_q) {
					exiting = 1;
				}
				if (e.key.keysym.sym == SDLK_SPACE) {
					tvr.level = 0.13f;
				}
				if (e.key.keysym.sym == SDLK_s) {
					desync = 1;
				}
				if (e.key.keysym.sym == SDLK_x) {
					sync = 0;
				}
			}
			if (e.type == SDL_KEYUP) {
				if (e.key.keysym.sym == SDLK_SPACE) {
					tvr.level = 1.5f;
				}
				if (e.key.keysym.sym == SDLK_s) {
					desync = 0;
				}
				if (e.key.keysym.sym == SDLK_x) {
					sync = 1;
				}
			}
		}

		tv->flutter_horiz_desync = 1;

		if (sync) {
			analogtv_setup_synclevel(tvi, 1, -18 + (desync ? 20 : 0));
		}
		analogtv_setup_frame(tv);

		int ntsc[4];
		analogtv_rgb_to_ntsc(0, 100, 255, ntsc);
		int x = 100;
		int y = 100;
		analogtv_draw_solid(tvi, x * 4, x * 4 + 64, y, y + 32, ntsc);
		analogtv_rgb_to_ntsc(255, 150, 0, ntsc);
		x += 50;
		y += 50;
		analogtv_draw_solid(tvi, x * 4, x * 4 + 64, y, y + 32, ntsc);

		analogtv_reception_update(&tvr);
		{
			const analogtv_reception* tvrp = &tvr;
			analogtv_draw(tv, 0.09, &tvrp, 1);
		}

		//#define RENDER_TO_DISK
		#ifdef RENDER_TO_DISK
		char filename[1024];
		sprintf(filename, "frame%05d.ppm", frame);
		FILE* out = fopen(filename, "w");
		fprintf(out, "P3\n%d %d\n255\n", width, height);
		int x0, y0;
		for (y0 = 0; y0 < height; y0++) {
			for (x0 = 0; x0 < width; x0++) {
				uint8_t* pixel = tv->image + y0 * tv->bytes_per_line + x0 * tv->bytes_per_pixel;
				fprintf(out, "%d %d %d\n", pixel[2], pixel[1], pixel[0]);
			}
		}
		fclose(out);
		#endif

		SDL_UpdateTexture(sdl_texture, NULL, tv->image, width * sizeof(uint32_t));
		//SDL_RenderClear(sdl_renderer);
		SDL_RenderCopy(sdl_renderer, sdl_texture, NULL, NULL);
		SDL_RenderPresent(sdl_renderer);

		frame++;

		//printf("frame: %d\n", frame);
	}

	return EXIT_SUCCESS;
}


