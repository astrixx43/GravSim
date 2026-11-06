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

const double C =  0.1*1920.0 + 0 * 299792458.0;
const double ONE_C = 1/C;
double G = 100 * 6.67430 * 0.000000000001;

const int NUM_BODIES = 3;
const int NUM_TRACES = 10; // Set KK equal to this and set this to 4-6 give good lookinfg results

const double MIN_RADIUS = 50;

const int WALL_THICKNESS = 10;
const uint32_t WALL_COLOUR = 0X0F00F000;


const int WIDTH = 1920;
const int HEIGHT = 1080;
double EPSILON = 0.001; // simmulation step size. 1 == a second

bool SHOWTEXT = false;

bool SHOWGRAVITON = true;

uint64_t FRAMERATE = 60;

#define TRACE_INTERVAL 100
#define NUM_SEGMENTS 75
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
                
        if (SHOWTEXT){
            snprintf(body_str, sizeof(body_str), "(%.2fC, %.2fC)\0", bodies[i].velocity[0]*ONE_C, bodies[i].velocity[1]*ONE_C);
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
    int NUM_RINGS_SKIP = 0;
    for (int i = 0; i< VBM->num_bodies; i ++){
        gravitonNode* graviton = ((gravitonList*) *VBM->bodies[i].gravitonList)->first;
        int j = 0;
        
        while (j < NUM_RINGS_SKIP && graviton != NULL && graviton->node != NULL){
            graviton = graviton->next;
            j += 1;
        }

        int i_graviton = 0;
        int num_graviton_drawn = 0;
        SDL_SetRenderDrawColor(renderer, 00, 00, 00, 0XFF);
        
        SDL_FPoint points[NUM_SEGMENTS];
        while(graviton != NULL && graviton->node != NULL){
            j = 0;
            // if (!circleRGBA(renderer, graviton->node->pos[0], graviton->node->pos[1], (Sint16)graviton->node->radius+MIN_RADIUS*EPSILON, 00, 00, 00, 0XFF)){
            //     printf("FAILED TO DRAW GRAVITON");
            //     break;
            // };  
            
            for (int l = 0; l< NUM_SEGMENTS; l++){
                points[l+ (i_graviton * NUM_SEGMENTS)].x =  (float)(graviton->node->pos[0] + TRIGTABLE[l][0]*(Sint16)graviton->node->radius+MIN_RADIUS*EPSILON);
                points[l+ (i_graviton * NUM_SEGMENTS)].y =  (float)(graviton->node->pos[1] + TRIGTABLE[l][1]*(Sint16)graviton->node->radius+MIN_RADIUS*EPSILON);
            }

            num_graviton_drawn += 1;
          
            graviton = graviton->next;
            while (j < NUM_RINGS_SKIP && graviton != NULL && graviton->node != NULL){
                graviton = graviton->next;
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

int drawTextScreen(TTF_TextEngine* text_engine, TTF_Font* font){
    char body_str[20];
    snprintf(body_str, sizeof(body_str), "C = %.3f\0", C);
    // int fontSize = TTF_GetFontSize(font);
    TTF_Text* body_text = TTF_CreateText(text_engine, font, body_str, 0); 
    TTF_SetTextColor(body_text, 0xFF, 0x0F, 0xFF, 0xFF);
    
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
    

    GravSimModelManager* model_manager =  initalize_Sim(NUM_BODIES, EPSILON, C, G, WIDTH, HEIGHT, WALL_THICKNESS, MIN_RADIUS, 0);

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
    
    while (true)
    {
        prev_time = current_time;
        current_time = SDL_GetTicksNS();

        elapsed_time = ((current_time - prev_time) > 250000000 )? 250000000 : (current_time - prev_time);
        
        frame_time_accumilator += elapsed_time;
       
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
            GraveSim_model_take_step(model_manager);
            

            frame_time_accumilator -= EP_NS;
            
            kk += 1;
        }

        
        // 3. Draw the frame

        if (drawFrame(renderer, text_engine, view_body_manager) != 0){
            
            TTF_DestroyRendererTextEngine(text_engine);
            SDL_DestroyRenderer(renderer);
            SDL_DestroyWindow(window);
            SDL_Quit();
            GraveSim_model_Sim(model_manager);
            GravSimViewBodies_end_bodyView(view_body_manager);
            return -1;
        }

        // if ( SHOWTEXT && !drawTextScreen(text_engine, fonts[NUM_BODIES])){
        //     closeFonts(fonts);
        //     TTF_DestroyRendererTextEngine(text_engine);
        //     SDL_DestroyRenderer(renderer);
        //     SDL_DestroyWindow(window);
        //     SDL_Quit();
        //     end_Sim(MMbodies);
        //     end_bodyView(traces);
        //     return -1;
        // }

        SDL_RenderPresent(renderer);

        SDL_DelayPrecise((nanosecondsPerFrame - (SDL_GetTicksNS() - current_time) > 0 )? (nanosecondsPerFrame - (SDL_GetTicksNS() - current_time)) : 0);
    }

    GraveSim_model_Sim(model_manager);
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