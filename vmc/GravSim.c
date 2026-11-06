// Give me a pencil and eternity, and I will write eventually.

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <SDL3/SDL.h>
#include <SDL3_ttf/SDL_ttf.h>

#include "SDL3_gfxPrimitives.h"
#include "SDL3_framerate.h"

#include "src/GravSimViewBodies.h"

#define TRACE_MAX_ALPHA 120  

const double C =  .01*1920.0 + 0 * 299792458.0;
const double ONE_C = 1/C;
double G =  6.67430 * 0.000000000001;

const int NUM_BODIES = 2;
const int NUM_TRACES = 10; // Set KK equal to this and set this to 4-6 give good lookinfg results

const int WALL_THICKNESS = 10;
const uint32_t WALL_COLOUR = 0X0F00F000;

const int WIDTH = 1920;
const int HEIGHT = 1080;
double EPSILON = 0.0001; // simmulation step size. 1 == a second

bool SHOWTEXT = true;

bool SHOWSPEED = false;

bool SHOWFPS = true;

bool SHOWGRAVITON = false;

uint64_t FRAMERATE = 60;

#define GRAV_TYPE 2

#define TRACE_INTERVAL 7
#define MAX_GRAVITONS_DRAWN 3000
#define GRAVITON_INTERVAL 10
#define NUM_SEGMENTS 50

#define NUM_FRAMES_TO_AVERAGE 60.0

double TRIGTABLE[NUM_SEGMENTS][2];


void renderBodies(SDL_Renderer* renderer, TTF_TextEngine *text_engine, viewBodiesManager* VBM){
    Uint32 c = 0;
    
    TTF_Text* body_text;
    char body_str[40];
    uint32_t fontSize;

    viewBody* bodies = VBM->bodies;
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

void renderTraces(SDL_Renderer* renderer, viewBodiesManager* VBM){
    uint32_t c = 0;
    for (int i = 0; i < NUM_BODIES; i++){
        for (int age = VBM->num_traces - 1; age >= 0; age--){

            int j = (VBM->indexOfTrace - 1 - age + VBM->num_traces) % VBM->num_traces;

            double scale = 1.0 - (double)(age + 1) / (3.0 * (VBM->num_traces + 1));

            c = VBM->traces[i][j].colour;

            if (!filledCircleRGBA(renderer, VBM->traces[i][j].pos[0], VBM->traces[i][j].pos[1], (int)(scale * VBM->traces[i][j].radius), (c >> 24) & 0xFF, (c >> 16) & 0xFF, (c >> 8) & 0xFF, (uint8_t)(scale * (c & 0xFF)))){
                printf("FAILED TO DRAW TRACE BODY");
            };

            // SDL_SetRenderDrawColor(renderer, (c >> 24) & 0xFF, (c >> 16) & 0xFF, (c >> 8) & 0xFF, (uint8_t)(scale * (c & 0xFF)));

            // SDL_FPoint points[NUM_SEGMENTS];

            // for (int l = 0; l<NUM_SEGMENTS; l++){
            //     points[l].x =  (float)(VBM->traces[i][j].pos[0] + TRIGTABLE[l][0]*VBM->traces[i][j].radius);
            //     points[l].y =  (float)(VBM->traces[i][j].pos[1] + TRIGTABLE[l][1]*VBM->traces[i][j].radius);
            // }

            // SDL_RenderLines(renderer, points, NUM_SEGMENTS);

        }
    }
}

void renderGravitons(SDL_Renderer* renderer, viewBodiesManager* VBM){   

    for (int i = 0; i< VBM->num_bodies; i ++){
        gravitonNode* graviton = ((gravitonList*) *VBM->bodies[i].gravitonList)->last->prev;
        int j = 0;
        
        while (j < GRAVITON_INTERVAL && graviton != NULL && graviton->node != NULL){
            graviton = graviton->prev;
            j += 1;
        }

        int i_graviton = 0;
        int num_graviton_drawn = 0;

        
        SDL_FPoint points[NUM_SEGMENTS];
        while(num_graviton_drawn < MAX_GRAVITONS_DRAWN && graviton != NULL && graviton->node != NULL){
            j = 0;
            SDL_SetRenderDrawColor(renderer, 0XFF/(EPSILON+graviton->node->radius), 0XFF/(EPSILON+graviton->node->radius), 0XFF/(EPSILON+graviton->node->radius), 0XFF/(EPSILON+graviton->node->radius));
            for (int l = 0; l< NUM_SEGMENTS; l++){
                points[l+ (i_graviton * NUM_SEGMENTS)].x =  (float)(graviton->node->pos[0] + (TRIGTABLE[l][0]*(graviton->node->radius+1)));
                points[l+ (i_graviton * NUM_SEGMENTS)].y =  (float)(graviton->node->pos[1] + (TRIGTABLE[l][1]*(graviton->node->radius+1)));
            }

            num_graviton_drawn += 1;
          
            graviton = graviton->prev;
            while (j < GRAVITON_INTERVAL && graviton != NULL && graviton->node != NULL){
                graviton = graviton->prev;
                j += 1;
            }

            SDL_RenderLines(renderer, points, NUM_SEGMENTS);
        }

    }
}

int drawFrame(SDL_Renderer * renderer, TTF_TextEngine* text_engine, viewBodiesManager* view_body_manager){
    // Draw bacground colours
    if (!(SDL_SetRenderDrawColor(renderer, 20, 30, 60, 255) && SDL_RenderClear(renderer))){
        printf("SDL_SetRenderDrawColor || SDL_RenderClear didn't work: %s\n", SDL_GetError());
        return -1; 
    }

    renderTraces(renderer, view_body_manager);
    renderBodies(renderer, text_engine, view_body_manager);
    if (SHOWGRAVITON) renderGravitons(renderer, view_body_manager);
    
    // Draw Walls
    if (!(thickLineColor(renderer, 0, 0, WIDTH + WALL_THICKNESS/2, 0, WALL_THICKNESS, WALL_COLOUR) + thickLineColor(renderer, WIDTH, 0, WIDTH, HEIGHT+ WALL_THICKNESS/2, WALL_THICKNESS, WALL_COLOUR) + thickLineColor(renderer, 0, HEIGHT, WIDTH + WALL_THICKNESS/2, HEIGHT, WALL_THICKNESS, WALL_COLOUR) + thickLineColor(renderer, 0, 0, 0, HEIGHT+WALL_THICKNESS/2, WALL_THICKNESS, WALL_COLOUR) == 4)){
        printf("SDL_RenderLine didn't work: %s\n", SDL_GetError());
        return -1;
    }        



    return 0;
}

int drawTextScreen(TTF_TextEngine* text_engine, TTF_Font* font, double fps){
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

int main(){

    int kk = 0;

    if (!SDL_Init(SDL_INIT_VIDEO)) {
        printf("SDL_Init failed: %s\n", SDL_GetError());
        return -1;
    }
    SDL_Window* window = SDL_CreateWindow("Hello SDL World", WIDTH, HEIGHT, SDL_WINDOW_BORDERLESS);

    if (window == NULL){
        printf("SDL_CreateWindow Window didn't work: %s\n", SDL_GetError());
        return -1;
    }

    SDL_Renderer* renderer = SDL_CreateRenderer(window, NULL);

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

    TTF_TextEngine *text_engine = NULL;
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
    

    GravSimModelManager* model_manager =  initalize_Sim(NUM_BODIES, EPSILON, C, G, WIDTH, HEIGHT, WALL_THICKNESS, GRAV_TYPE , notifyChange);

    viewBodiesManager* view_body_manager = viewBodiesManage_init(model_manager, NUM_TRACES, SHOWTEXT);
    setUpTrigTable();


    SDL_Event windowEvent;


    
    // double stepsPerSecond = 1.0/(EPSILON);
    double framePerSecond = FRAMERATE;


    double nanosecondsPerFrame = 1000000000/framePerSecond;

    
    Uint64 current_time = SDL_GetTicksNS();
    Uint64 prev_time = SDL_GetTicksNS();
    Uint64 elapsed_time = 0;


    double frame_time_accumilator = 0.0;
    double EP_NS = EPSILON * 1000000000;
    
    double sumActualTimeBetweenFrames = 0.0;
    int firstMemberOfSum = 0;
    int fps_avg_index = 0;
    bool ready_to_sub_avg_fps = 0;
    double actualTimeBeteenFrames[(int)(NUM_FRAMES_TO_AVERAGE)];
    double averageFPS = 0;
    while (true)
    {
        prev_time = current_time;
        current_time = SDL_GetTicksNS();

        elapsed_time = ((current_time - prev_time) > 250000000 )? 250000000 : (current_time - prev_time);
        
        frame_time_accumilator += elapsed_time;
        
        sumActualTimeBetweenFrames += elapsed_time;
        
        actualTimeBeteenFrames[fps_avg_index] = elapsed_time;
        fps_avg_index = (fps_avg_index >= NUM_FRAMES_TO_AVERAGE)? 0 : fps_avg_index + 1;;

        if (SDL_PollEvent(&windowEvent) && windowEvent.type == SDL_EVENT_QUIT){
            break;
        }
        
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
        averageFPS = NUM_FRAMES_TO_AVERAGE * 1000000000/sumActualTimeBetweenFrames;

        
        // 3. Draw the frame

        if (drawFrame(renderer, text_engine, view_body_manager) != 0){
            
            TTF_DestroyRendererTextEngine(text_engine);
            SDL_DestroyRenderer(renderer);
            SDL_DestroyWindow(window);
            SDL_Quit();
            GravSim_clean(model_manager);
            GravSimViewBodies_end_bodyView(view_body_manager);
            return -1;
        }
        
        if (SHOWTEXT && SHOWFPS){
            if (drawTextScreen(text_engine, view_body_manager->menuFont, averageFPS) < 0){
                TTF_DestroyRendererTextEngine(text_engine);
                SDL_DestroyRenderer(renderer);
                SDL_DestroyWindow(window);
                SDL_Quit();
                GravSim_clean(model_manager);
                GravSimViewBodies_end_bodyView(view_body_manager);
                return -1;
            }

        }

        if(ready_to_sub_avg_fps || fps_avg_index >= (NUM_FRAMES_TO_AVERAGE)){
            sumActualTimeBetweenFrames -= actualTimeBeteenFrames[firstMemberOfSum];
            firstMemberOfSum = ((firstMemberOfSum) >= NUM_FRAMES_TO_AVERAGE)? 0 : firstMemberOfSum + 1;
            ready_to_sub_avg_fps = 1;
        }

        
        // if (ready_to_sub_avg_fps && averageFPS < FRAMERATE/1.5){
        //     GravSim_change_grav_func(model_manager, 1);
        // }else if (fps_avg_index >= (NUM_FRAMES_TO_AVERAGE) && GRAV_TYPE != 1 && GravSim_get_grav_type(model_manager) == 1 && averageFPS > FRAMERATE/1.5){
        //     GravSim_change_grav_func(model_manager, 0);
        // }


        
        SDL_RenderPresent(renderer);

        SDL_DelayPrecise((nanosecondsPerFrame - (SDL_GetTicksNS() - current_time) > 0 )? (nanosecondsPerFrame - (SDL_GetTicksNS() - current_time)) : 0);
    }
    
    GravSim_clean(model_manager);
    GravSimViewBodies_end_bodyView(view_body_manager);
    
    if (SHOWTEXT){
        TTF_DestroyRendererTextEngine(text_engine);
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    TTF_Quit();
    SDL_Quit();
    
    return 0;
}