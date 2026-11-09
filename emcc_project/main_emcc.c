#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#endif

// Give me a pencil and eternity, and I will write eventually.

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <SDL3/SDL.h>
#include <SDL3_ttf/SDL_ttf.h>

#include "SDL3_gfxPrimitives.h"
#include "SDL3_framerate.h"

#include "../src/GravSimViewBodies.h"




/* Below Are all values that may be adjusted for your simulation. */

/* Simmulation Constants*/
const double EPSILON = 0.01; // Simmulation step size in second units. i.e EPSILON == 1 => every sim. update will tick the sim clock 1 second
const double C =  1920.0 + 0 * 299792458.0; // The Constant of the speed of light
double G =  6.67430 * 0.000000000001; // The Gravitational Constant

const int NUM_BODIES = 2; // The number of massive body to be simmulated.
const double MIN_SPEED = 0.15; // The minimum magnitude of the velocity of a body as a fraction of C. MUST BE << C_speed
const double MAX_SPEED = 0.20; // The maximum magnitude of the velocity of a body as a fraction of C. MUST BE  > MIN_SPEED

double MIN_RADIUS = 10; // Minimum radius of a massive body
double MAX_RADIUS = 50; // Maximum radius of a massive body

#define GRAV_TYPE 0 // 0 => Wave-Like Gravitons, 1 => Particled Gravitons, 2 => Newtonian Gravity. More to come...


/* Visualizor Constants*/
bool SHOWTEXT = true; // Display ANY text on screen or not. For best preformence turn all text false
bool SHOWSPEED = false; // Display C normalized velocity values for each body. For best preformence turn all text false
bool SHOWFPS = true; // Display true frames per second and sim's C value. For best preformence turn all text false

bool SHOWGRAVITON = true; // Display gravtionso or not. Safe to be true for non-graviton/Newtonian gravity type. For best preformence turn false
#define MAX_GRAVITONS_DRAWN 10 // Max number of gravitons to be drawn for each massive body. For particles, set low, < 70. For waves set 80. DOES NOT EFFECT GRAVITON EMMISION
#define GRAVITON_INTERVAL 1 // Number of gravitons to skip in drawing gravitons. For particles, set high, > 30. For waves set between 0 - 20. DOES NOT EFFECT GRAVITON EMMISION

const int NUM_TRACES = 10; // This value controls the total number of the tails of fast moving bodies. The higher the longer.
#define TRACE_INTERVAL 3  // This value controles the length between each update of the tail of fast moving bodies. The higher the longer.

const uint64_t FRAMERATE = 60; // Framterate of simmulation. Keep in mind, the simulation is running in real time.


/* You may Ignore The Rest ... */
#define TRACE_MAX_ALPHA 120  
const uint32_t WALL_COLOUR = 0X0F00F000;
const int WALL_THICKNESS = 10;

int WIDTH = 1920/2;
int HEIGHT = 1080/2;
const double ONE_C = 1/C;

#define NUM_SEGMENTS 100 // How smooth gravtion's edges should be. The lower, the better for proformance
double TRIGTABLE[NUM_SEGMENTS][2]; // A table of NUM_SEGMENTS many of [sin, cos] values form 0 to 2pi. Commpute once.
#define NUM_FRAMES_TO_AVERAGE_OVER 60.0
#define WINDOW_TYPE SDL_WINDOW_RESIZABLE

int kk = 0;

SDL_Event windowEvent;



// double stepsPerSecond = 1.0/(EPSILON);
const double framePerSecond = FRAMERATE;
double nanosecondsPerFrame = 1000000000/framePerSecond;


Uint64 current_time;
Uint64 prev_time;
Uint64 elapsed_time = 0;


double frame_time_accumilator = 0.0;
double EP_NS = EPSILON * 1000000000;

double sumActualTimeBetweenFrames = 0.0;
int firstMemberOfSum = 0;
int fps_avg_index = 0;
bool ready_to_sub_avg_fps = 0;
double actualTimeBeteenFrames[(int)(NUM_FRAMES_TO_AVERAGE_OVER)];
double averageFPS = 0;

GravSimModelManager* model_manager;
viewBodiesManager* view_body_manager;
TTF_TextEngine* text_engine;
SDL_Renderer* renderer;
SDL_Window* window;

void renderBodies(){
    Uint32 c = 0;
    
    TTF_Text* body_text;
    char body_str[40];
    uint32_t fontSize;

    viewBody* bodies = view_body_manager->bodies;
    for (int i = 0; i < NUM_BODIES; i++){
        c = *bodies[i].relative_colour;

        if (!filledCircleRGBA(renderer, bodies[i].pos[0], bodies[i].pos[1], (Sint16)*bodies[i].radius, (c >> 24) & 0xFF, (c >> 16) & 0xFF, (c >> 8) & 0xFF, c & 0xFF)){
            printf("FAILED TO DRAW BODY");
        };
                
        if (SHOWSPEED && SHOWTEXT){

            snprintf(body_str, sizeof(body_str), "(%.2fC, %.2fC)", bodies[i].velocity[0]*ONE_C, bodies[i].velocity[1]*ONE_C);
            
            fontSize = TTF_GetFontSize(bodies[i].font);
            body_text = TTF_CreateText(text_engine, bodies[i].font, body_str, 0); 
            
            TTF_SetTextColor(body_text, 0XFF - (((*bodies[i].relative_colour) >> 24) & 0X00FF), 0XFF - (((*bodies[i].relative_colour) >> 16) & 0X00FF), 0XFF - (((*bodies[i].relative_colour) >> 8)  & 0X00FF), (*bodies[i].relative_colour) & 0X000000FF );
            
            if (!TTF_DrawRendererText(body_text, bodies[i].pos[0] - (fontSize + 6)* 1.3333 - *bodies[i].radius/2, bodies[i].pos[1] - (fontSize - 2)* 1.3333 +(strlen(body_str)/2))){
                // printf("FAILED TO DRAW TEXT: %s\n", SDL_GetError());
            }
            TTF_DestroyText(body_text);
        }

        
    }
}

void renderTraces(){
    uint32_t c = 0;
    for (int i = 0; i < NUM_BODIES; i++){
        for (int age = view_body_manager->num_traces - 1; age >= 0; age--){

            int j = (view_body_manager->indexOfTrace - 1 - age + view_body_manager->num_traces) % view_body_manager->num_traces;

            double scale = 1.0 - (double)(age + 1) / (3.0 * (view_body_manager->num_traces + 1));

            c = view_body_manager->traces[i][j].colour;

            if (!filledCircleRGBA(renderer, view_body_manager->traces[i][j].pos[0], view_body_manager->traces[i][j].pos[1], (int)(scale * view_body_manager->traces[i][j].radius), (c >> 24) & 0xFF, (c >> 16) & 0xFF, (c >> 8) & 0xFF, (uint8_t)(scale * (c & 0xFF)))){
                printf("FAILED TO DRAW TRACE BODY");
            };

            // SDL_SetRenderDrawColor(renderer, (c >> 24) & 0xFF, (c >> 16) & 0xFF, (c >> 8) & 0xFF, (uint8_t)(scale * (c & 0xFF)));

            // SDL_FPoint points[NUM_SEGMENTS];

            // for (int l = 0; l<NUM_SEGMENTS; l++){
            //     points[l].x =  (float)(view_body_manager->traces[i][j].pos[0] + TRIGTABLE[l][0]*view_body_manager->traces[i][j].radius);
            //     points[l].y =  (float)(view_body_manager->traces[i][j].pos[1] + TRIGTABLE[l][1]*view_body_manager->traces[i][j].radius);
            // }

            // SDL_RenderLines(renderer, points, NUM_SEGMENTS);

        }
    }
}

void renderGravitons(){   

    for (int i = 0; i< view_body_manager->num_bodies; i ++){
        gravitonNode* graviton = ((gravitonList*) *view_body_manager->bodies[i].gravitonList)->last->prev;
        int j = 0;
        
        double num_graviton = 0;
        while (j < GRAVITON_INTERVAL && graviton != NULL && graviton->node != NULL){
            graviton = graviton->prev;
            j += 1;
        }

        num_graviton += j;

        int i_graviton = 0;
        int num_graviton_drawn = 0;

        
        SDL_FPoint points[NUM_SEGMENTS];
        
        while(num_graviton_drawn < MAX_GRAVITONS_DRAWN && graviton != NULL && graviton->node != NULL){
            j = 0;
            
            double scale = num_graviton/((gravitonList*) view_body_manager->bodies[i].gravitonList)->count;

            SDL_SetRenderDrawColor(renderer, 0XFF, 0XFF, 0XFF, 0XFF/(1-scale));
            
            for (int l = 0; l< NUM_SEGMENTS; l++){
                num_graviton  += 1;
                
                points[l+ (i_graviton * NUM_SEGMENTS)].x =  (float)(graviton->node->pos[0] + (TRIGTABLE[l][0]*(graviton->node->radius)));
                points[l+ (i_graviton * NUM_SEGMENTS)].y =  (float)(graviton->node->pos[1] + (TRIGTABLE[l][1]*(graviton->node->radius)));
            }

            num_graviton_drawn += 1;
          
            graviton = graviton->prev;
            
            while (j < GRAVITON_INTERVAL && graviton != NULL && graviton->node != NULL){
                graviton = graviton->prev;
                j += 1;
            }
            num_graviton += j;

            SDL_RenderLines(renderer, points, NUM_SEGMENTS);
        }

    }
}

int drawFrame(){
    // Draw bacground colours
    if (!(SDL_SetRenderDrawColor(renderer, 20, 30, 60, 255) && SDL_RenderClear(renderer))){
        printf("SDL_SetRenderDrawColor || SDL_RenderClear didn't work: %s\n", SDL_GetError());
        return -1; 
    }

    renderTraces();
    renderBodies();
    if (SHOWGRAVITON) renderGravitons();
    
    // Draw Walls
    if (!(thickLineColor(renderer, 0, 0, WIDTH + WALL_THICKNESS/2, 0, WALL_THICKNESS, WALL_COLOUR) + thickLineColor(renderer, WIDTH, 0, WIDTH, HEIGHT+ WALL_THICKNESS/2, WALL_THICKNESS, WALL_COLOUR) + thickLineColor(renderer, 0, HEIGHT, WIDTH + WALL_THICKNESS/2, HEIGHT, WALL_THICKNESS, WALL_COLOUR) + thickLineColor(renderer, 0, 0, 0, HEIGHT+WALL_THICKNESS/2, WALL_THICKNESS, WALL_COLOUR) == 4)){
        printf("SDL_RenderLine didn't work: %s\n", SDL_GetError());
        return -1;
    }        

    return 0;
}

int drawTextScreen(TTF_Font* font, double fps){
    char body_str[40];
    if(!snprintf(body_str, sizeof(body_str), "C = %.1fm/s\n FPS: %.2f", C, fps)){
        printf("FAILED TO WRITE TEXT TO body IN 'drawTextScreen'\n");
        return -1;
    }

    // int fontSize = TTF_GetFontSize(font);
    
    TTF_Text* body_text = TTF_CreateText(text_engine, font, body_str, 0);
    if(body_text == NULL){
        printf("FAILED TO MAKE 'body_text' IN 'drawTextScreen': %s\n", SDL_GetError());
        return -1;
    } 
    if(!TTF_SetTextColor(body_text, 0xFF, 0x0F, 0xFF, 0xFF)){
        printf("FAILED TO MAKE 'TTF_SetTextColor' IN 'drawTextScreen': %s\n", SDL_GetError());
        return -1;
    }
    
    if (!TTF_DrawRendererText(body_text, 0, 0)){
        printf("FAILED TO DRAW TEXT: %s\n", SDL_GetError());
        return -1;
    }

    TTF_DestroyText(body_text);
    return 1;
}

void setUpTrigTable(){
    double step_size = 2.0*M_PI/(NUM_SEGMENTS - 1);
    for(int i = 0; i < NUM_SEGMENTS; i++){
        TRIGTABLE[i][0] = sin(i*step_size);
        TRIGTABLE[i][1] = cos(i*step_size);
    }
}

void notifyChange(){


}

void main_end(void){
        
    GravSim_clean(model_manager);
    GravSimViewBodies_end_bodyView(view_body_manager);
    
    if (SHOWTEXT){
        TTF_DestroyRendererTextEngine(text_engine);
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    TTF_Quit();
    SDL_Quit();
    #ifdef __EMSCRIPTEN__
    emscripten_cancel_main_loop();
    #else
    exit(0);
    #endif
}


void main_loop(void) {
    // poll events, update, render — one frame's worth

    prev_time = current_time;
    current_time = SDL_GetTicksNS();

    elapsed_time = ((current_time - prev_time) > 250000000 )? 250000000 : (current_time - prev_time);
    
    frame_time_accumilator += elapsed_time;
    
    sumActualTimeBetweenFrames += elapsed_time;
    
    actualTimeBeteenFrames[fps_avg_index] = elapsed_time;
    fps_avg_index = (fps_avg_index >= NUM_FRAMES_TO_AVERAGE_OVER)? 0 : fps_avg_index + 1;;

    if (SDL_PollEvent(&windowEvent) && windowEvent.type == SDL_EVENT_QUIT){
        main_end();
        return;
    }else if(windowEvent.type == SDL_EVENT_WINDOW_RESIZED) {
        int new_width = windowEvent.window.data1;
        int new_height = windowEvent.window.data2;
        GravSim_update_box_size(model_manager, new_width, new_height);
        WIDTH = new_width;
        HEIGHT = new_height;
    
    }
    // Update your viewport, camera matrices, or UI layout here
    
    while (frame_time_accumilator > EP_NS){
        
        if (++kk > TRACE_INTERVAL){
            // 1. Update Traces
            GravSimViewBodies_updateTraces(view_body_manager);
            kk = 0;
        }
        
        // 2. Update Bodies
        GravSim_model_take_step(model_manager);
        

        frame_time_accumilator -= EP_NS;
    
    }
    averageFPS = NUM_FRAMES_TO_AVERAGE_OVER * 1000000000/sumActualTimeBetweenFrames;

    
    // 3. Draw the frame

    if (drawFrame() != 0){
        
        TTF_DestroyRendererTextEngine(text_engine);
        SDL_DestroyRenderer(renderer);
        SDL_DestroyWindow(window);
        SDL_Quit();
        GravSim_clean(model_manager);
        GravSimViewBodies_end_bodyView(view_body_manager);
        main_end();
    }
    
    if (SHOWTEXT && SHOWFPS){
        if (drawTextScreen(view_body_manager->menuFont, averageFPS) < 0){
            TTF_DestroyRendererTextEngine(text_engine);
            SDL_DestroyRenderer(renderer);
            SDL_DestroyWindow(window);
            SDL_Quit();
            GravSim_clean(model_manager);
            GravSimViewBodies_end_bodyView(view_body_manager);
            main_end();
        }

    }

    if(ready_to_sub_avg_fps || fps_avg_index >= (NUM_FRAMES_TO_AVERAGE_OVER)){
        sumActualTimeBetweenFrames -= actualTimeBeteenFrames[firstMemberOfSum];
        firstMemberOfSum = ((firstMemberOfSum) >= NUM_FRAMES_TO_AVERAGE_OVER)? 0 : firstMemberOfSum + 1;
        ready_to_sub_avg_fps = 1;
    }

    
    // if (ready_to_sub_avg_fps && averageFPS < FRAMERATE/1.5){
    //     GravSim_change_grav_func(model_manager, 1);
    // }else if (fps_avg_index >= (NUM_FRAMES_TO_AVERAGE_OVER) && GRAV_TYPE != 1 && GravSim_get_grav_type(model_manager) == 1 && averageFPS > FRAMERATE/1.5){
    //     GravSim_change_grav_func(model_manager, 0);
    // }


    
    SDL_RenderPresent(renderer);

    SDL_DelayPrecise((nanosecondsPerFrame - (SDL_GetTicksNS() - current_time) > 0 )? (nanosecondsPerFrame - (SDL_GetTicksNS() - current_time)) : 0);

}


int main(int argc, char *argv[]){

    if (!SDL_Init(SDL_INIT_VIDEO)) {
        printf("SDL_Init failed: %s\n", SDL_GetError());
        return -1;
    }
    window = SDL_CreateWindow("GravSim Visualizer", WIDTH, HEIGHT, WINDOW_TYPE);

    if (window == NULL){
        printf("SDL_CreateWindow Window didn't work: %s\n", SDL_GetError());
        return -1;
    }

    renderer = SDL_CreateRenderer(window, NULL);

    if (renderer == NULL){
        printf("SDL_CreateRenderer Renderer didn't work: %s\n", SDL_GetError());
        
        SDL_DestroyWindow(window);
        SDL_Quit();
        
        return -1;
    }

    if (SHOWTEXT && !TTF_Init()) {
        printf("TTF_Init failed: %s\n", SDL_GetError());
        SDL_DestroyRenderer(renderer);
        SDL_DestroyWindow(window);
        SDL_Quit();
        return -1;
    }

    text_engine = NULL;
    // TTF_Font** fonts = NULL;

    if (SHOWTEXT){
        text_engine = TTF_CreateRendererTextEngine(renderer);
    }
    
    if (SHOWTEXT && text_engine == NULL){
        printf("TTF_CreateRendererTextEngine Renderer didn't work: %s\n", SDL_GetError());
        SDL_DestroyRenderer(renderer);
        SDL_DestroyWindow(window);
        SDL_Quit();
        
        return -1;
    }
    

    model_manager =  initalize_Sim(NUM_BODIES, EPSILON, C, G, WIDTH, HEIGHT, WALL_THICKNESS, GRAV_TYPE, MIN_SPEED, MAX_SPEED, MIN_RADIUS, MAX_RADIUS, notifyChange);

    view_body_manager = viewBodiesManage_init(model_manager, NUM_TRACES, SHOWTEXT);
    setUpTrigTable();


#ifdef __EMSCRIPTEN__
    emscripten_set_main_loop(main_loop, 0, 1);  // 0 = use browser's refresh rate
#else
    while (true) main_loop();
#endif
    // cleanup
    return 0;
}