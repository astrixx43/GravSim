#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <SDL3/SDL.h>
#include <SDL3_ttf/SDL_ttf.h>

#include "SDL3_gfxPrimitives.h"
#include "SDL3_framerate.h"
#include "src/bodies.h"
#include "src/Wavelength_Rgb.h"

#define TRACE_MAX_ALPHA 120  

const double C =  1920.0 + 0 * 299792458.0;
const double C_SQUARED = C*C;
double G = 6.67430 * 0.000000000001;

const double MIN_SPEED = 0.001; // The minimum speed of body i of C. MUST BE << C
const double MAX_SPEED = 0.10; // The minimum speed of body i of C. MUST BE << C


const int NUM_BODIES = 30;
const int NUM_TRACES = 10; // Set KK equal to this and set this to 4-6 give good lookinfg results
int indexOfTrace = NUM_TRACES - 1;

const int WALL_THICKNESS = 10;
const uint32_t WALL_COLOUR = 0X0F00F000;


int WIDTH = 1920;
int HEIGHT = 1080;
double EPSILON = 0.001; // simmulation step size. 1 == a second

bool SHOWTEXT = true;

uint64_t FRAMERATE = 100;

double Q_rsqrt( double number )
{
    union {
        double f;
        uint64_t i;
    } conv;
    
    double x2;
    const double threehalfs = 1.5;

    x2 = number * 0.5;
    conv.f = number;
    
    conv.i = 0x5fe6ec85e7de30da - (conv.i >> 1); // what the fuck?
    

    conv.f = conv.f * (threehalfs - (x2 * conv.f * conv.f)); // 1st iteration
    
    conv.f = conv.f * (threehalfs - (x2 * conv.f * conv.f));

    return conv.f;
}

double calculateGamma(double speed){
    return Q_rsqrt((double)( 1 - ((speed*speed)/C_SQUARED)));
}

double calculateGamma_alt(double speed_squared){
    return Q_rsqrt((double)( 1 - ((speed_squared)/C_SQUARED)));
}

void relativisticDopplerEffect(body* body){
    int rest_wl = hexToWavelength(body->rest_colour);
    if (rest_wl < 0) {                 // grey/dark colour: no single wavelength
        body->relative_colour = body->rest_colour;
        return;
    }

    // Vector from the observer (screen centre) to the body
    double r_x = body->pos[0] - (WIDTH / 2.0);
    double r_y = body->pos[1] - (HEIGHT / 2.0);

    // v · r : positive = moving away from the centre, negative = moving toward it
    double radial = body->velocity[0] * r_x + body->velocity[1] * r_y;

    double speed = body->speed;
    double beta = speed / C;
    if (beta >= 1.0) beta = 0.999999;  // guard against division by zero

    double one_r_len = Q_rsqrt(r_x * r_x + r_y * r_y);
    double v_r = (r_x != r_y || r_y != 0.0) ? radial * one_r_len : 0.0;   // velocity component away from centre
    double factor = (1.0 + v_r / C) * Q_rsqrt(1.0 - beta * beta);

    body->relative_colour = wavelengthToHex(rest_wl * factor);
}

void change_body_velocity(body* body, double v_x, double v_y){
    body->velocity[0] = v_x;
    body->velocity[1] = v_y;
    body->speed = pow((body->velocity[0]*body->velocity[0]) + (body->velocity[1]*body->velocity[1]), 1.0/2.0);
    relativisticDopplerEffect(body);
}

void boost(double *E, double* p_x, double* p_y, double V_x, double V_y, double V_squared, double gamma_V) {

    if (V_squared == 0.0) return;

    double Vp = V_x * (*p_x) + V_y * (*p_y);
    double k  = (gamma_V - 1.0) * Vp / V_squared - gamma_V * (*E) / (C_SQUARED);

    *E = gamma_V * (*E - Vp);

    *p_x += k * V_x;
    *p_y += k * V_y;
}

void realtivisticCollison(body* b1, body* b2){
    double gamma_1 = calculateGamma(b1->speed);
    double gamma_2 = calculateGamma(b2->speed);

    double E_1 = gamma_1*b1->mass*C_SQUARED;
    double E_2 = gamma_2*b2->mass*C_SQUARED;

    double p_1_x = gamma_1*b1->mass*b1->velocity[0];
    double p_1_y = gamma_1*b1->mass*b1->velocity[1];

    double p_2_x = gamma_2*b2->mass*b2->velocity[0];
    double p_2_y = gamma_2*b2->mass*b2->velocity[1];

    double One_E_total = 1 / (E_1 + E_2);

    double V_x = (p_1_x + p_2_x)*C_SQUARED*One_E_total;
    double V_y = (p_1_y + p_2_y)*C_SQUARED*One_E_total;
    double V_sqaured = V_x*V_x + V_y*V_y;
    double gamma_V = calculateGamma_alt((V_x*V_x)+(V_y*V_y));

    boost(&E_1, &p_1_x, &p_1_y, V_x, V_y, V_sqaured, gamma_V);
    boost(&E_2, &p_2_x, &p_2_y, V_x, V_y, V_sqaured, gamma_V);

    double one_n = Q_rsqrt((b1->pos[0] - b2->pos[0]) * (b1->pos[0] - b2->pos[0]) + (b1->pos[1] - b2->pos[1]) * (b1->pos[1] - b2->pos[1]));
    double n_x = one_n * (b2->pos[0] - b1->pos[0]);
    double n_y = one_n * (b2->pos[1] - b1->pos[1]);
    double pn = p_1_x* n_x + p_1_y*n_y;

    if (pn <= 0.0) return;  

    p_1_x -= 2 * pn * n_x;
    p_1_y -= 2 * pn * n_y;
    
    p_2_x += 2 * pn * n_x;  
    p_2_y += 2 * pn * n_y;


    boost(&E_1, &p_1_x, &p_1_y, -1.0 * V_x, -1.0 * V_y, V_sqaured, gamma_V);
    boost(&E_2, &p_2_x, &p_2_y, -1.0 * V_x, -1.0 * V_y, V_sqaured, gamma_V);

    change_body_velocity(b1, p_1_x * C_SQUARED / E_1,  p_1_y * C_SQUARED / E_1);
    change_body_velocity(b2, p_2_x * C_SQUARED / E_2, p_2_y * C_SQUARED / E_2);


}

int checkBodyCollision(body* b1, body*b2){
    double distance = pow((b1->pos[0] - b2->pos[0])*(b1->pos[0] - b2->pos[0]) + (b1->pos[1] - b2->pos[1])*(b1->pos[1] - b2->pos[1]), 1.0/2.0);
    return distance <= b1->radius + b2->radius;
}

int checkWallCollision(body* b){
    int hit = 0;
    if ((b->pos[0] - b->radius <= 0 + WALL_THICKNESS/2 && b->velocity[0] < 0) ||
        (b->pos[0] + b->radius >= WIDTH - WALL_THICKNESS/2 && b->velocity[0] > 0)) {
        change_body_velocity(b, -b->velocity[0], b->velocity[1]);
        hit = 1;
    }
    if ((b->pos[1] - b->radius <= 0 + WALL_THICKNESS/2 && b->velocity[1] < 0) ||
        (b->pos[1] + b->radius >= HEIGHT - WALL_THICKNESS/2 && b->velocity[1] > 0)) {
        change_body_velocity(b, b->velocity[0], -b->velocity[1]);
        hit = 1;
    }
    return hit;
}

void handelCollision(body* bodies, int num){


    for (int i = 0; i < num; i++){
        if (!checkWallCollision(&bodies[i])){
            for (int j = i + 1; j < num; j++){
                if (checkBodyCollision(&bodies[i], &bodies[j])){
                    realtivisticCollison(&bodies[i], &bodies[j]);
                }
            }
        }
    }

}

void makeNewBody(body* body){
    int radius = (rand() % 100) + 15;
    
    double x_pos = ((double)rand()/(double)RAND_MAX) * (double)(WIDTH - radius) + (double)radius;
    double y_pos = ((double)rand()/(double)RAND_MAX) * (double)(HEIGHT - radius) + (double)radius;
    double theta = ((double)rand()/(double)RAND_MAX) * 2.0f * (double)(M_PI);

    double v_total = (((double)rand()/(double)RAND_MAX) * (MAX_SPEED- MIN_SPEED) +  MIN_SPEED) * C;
    double v_x =  v_total * cos(theta);
    double v_y = v_total * sin(theta);

    uint64_t colour = wavelengthToHex((uint64_t)380 + (uint64_t)rand() % (700- 380));
    // uint32_t colour = wavelengthToHex(580);
    

    body_init(body, 1000000000000 * M_PI*(double)(radius*radius), x_pos, y_pos, v_x, v_y, radius, colour);
}

void setupBodies(body* bodies, int num){
    for (int i = 0; i < num; i++){
        makeNewBody(&bodies[i]);

        for (int j = 0; j < i; j++){
            if (checkBodyCollision(&bodies[i], &bodies[j])){
                makeNewBody(&bodies[i]);
                j = -1;
            }
        }

        relativisticDopplerEffect(&bodies[i]);
    }

}

void setupBodyTraces(traceBody** traces, body* bodies){
    for (int i = 0; i < NUM_BODIES; i++){
        traces[i] = malloc(NUM_TRACES*sizeof(traceBody));

        for (int j = 0; j < NUM_TRACES; j++){
            trace_init(&traces[i][j], &bodies[i]);
        }

    }
}

void newtonianGravity(body* bodies){
    for (int i = 0; i < NUM_BODIES; i++){

        double a_x = 0;
        double a_y = 0;
        for (int j = 0; j < NUM_BODIES; j++){
            if (i != j){
                double masses_G = bodies[j].mass*G;
                double r_normal_factor = Q_rsqrt((bodies[j].pos[0] - bodies[i].pos[0])*(bodies[j].pos[0] - bodies[i].pos[0]) + (bodies[j].pos[1] - bodies[i].pos[1])*(bodies[j].pos[1] - bodies[i].pos[1]));
                a_x += r_normal_factor * masses_G *(bodies[j].pos[0] - bodies[i].pos[0])/((bodies[j].pos[0] - bodies[i].pos[0])*(bodies[j].pos[0] - bodies[i].pos[0]) + (bodies[j].pos[1] - bodies[i].pos[1])*(bodies[j].pos[1] - bodies[i].pos[1]));
                a_y += r_normal_factor * masses_G *(bodies[j].pos[1] - bodies[i].pos[1])/((bodies[j].pos[0] - bodies[i].pos[0])*(bodies[j].pos[0] - bodies[i].pos[0]) + (bodies[j].pos[1] - bodies[i].pos[1])*(bodies[j].pos[1] - bodies[i].pos[1]));

            }
        }
    
        
        bodies[i].velocity[0] = bodies[i].velocity[0] + EPSILON * a_x;
        bodies[i].velocity[1] = bodies[i].velocity[1] + EPSILON * a_y;
        bodies[i].speed = pow((bodies[i].velocity[0]*bodies[i].velocity[0]) + (bodies[i].velocity[1] * bodies[i].velocity[1]), 1.0/2.0);
    }


}

void updateBodies(body* bodies, int num){
    
    handelCollision(bodies, num);
    newtonianGravity(bodies);

    for (int i = 0; i < num; i++){
        bodies[i].pos[0] = bodies[i].pos[0] + EPSILON*bodies[i].velocity[0];
        bodies[i].pos[1] = bodies[i].pos[1] + EPSILON*bodies[i].velocity[1];

        relativisticDopplerEffect(&bodies[i]);
    }
}

void updateTraces(traceBody** traces, body* bodies){
    uint32_t traceAlpha;

    for (int i = 0; i < NUM_BODIES; i++){
        trace_init(&traces[i][indexOfTrace], &bodies[i]);
    }

    indexOfTrace = (indexOfTrace + 1) % NUM_TRACES;
}

void renderBodies(SDL_Renderer* renderer, TTF_TextEngine *text_engine, TTF_Font** fonts, body* bodies){
    Uint32 c = 0;
    
    TTF_Text* body_text;
    char body_str[40];
    uint32_t fontSize;

    for (int i = 0; i < NUM_BODIES; i++){
        c = bodies[i].relative_colour;


        if (!filledCircleRGBA(renderer, bodies[i].pos[0], bodies[i].pos[1], bodies[i].radius, (c >> 24) & 0xFF, (c >> 16) & 0xFF, (c >> 8) & 0xFF, c & 0xFF)){
            printf("FAILED TO DRAW BODY");
        };

        
        
        if (SHOWTEXT){
            snprintf(body_str, sizeof(body_str), "(%.2fC, %.2fC)\0", bodies[i].velocity[0]/C, bodies[i].velocity[1]/C);
            fontSize = TTF_GetFontSize(fonts[i]);
            body_text = TTF_CreateText(text_engine, fonts[i], body_str, 0); 
            TTF_SetTextColor(body_text, 0XFF - (bodies[i].relative_colour >> 24) & 0X00FF, 0XFF - (bodies[i].relative_colour >> 16) & 0X00FF, 0XFF - (bodies[i].relative_colour >> 8)  & 0X00FF, (bodies[i].relative_colour) & 0X000000FF );
            
            if (!TTF_DrawRendererText(body_text, bodies[i].pos[0] - (fontSize + 6)* 1.3333 - bodies[i].radius/2, bodies[i].pos[1] - (fontSize - 2)* 1.3333 +(strlen(body_str)/2))){
                // printf("FAILED TO DRAW TEXT: %s\n", SDL_GetError());
            }
            TTF_DestroyText(body_text);
        }

        
    }
}

void renderTraces(SDL_Renderer* renderer, traceBody** traces){
    uint32_t c = 0;
    for (int i = 0; i < NUM_BODIES; i++){
        for (int age = NUM_TRACES - 1; age >= 0; age--){
            int j = (indexOfTrace - 1 - age + NUM_TRACES) % NUM_TRACES;
            double scale = 1.0 - (double)(age + 1) / (3.0 * (NUM_TRACES + 1));
            c = traces[i][j].colour;

            if (!filledCircleRGBA(renderer, traces[i][j].pos[0], traces[i][j].pos[1], (int)(scale * traces[i][j].radius), (c >> 24) & 0xFF, (c >> 16) & 0xFF, (c >> 8) & 0xFF, (uint8_t)(scale * (c & 0xFF)))){
                printf("FAILED TO DRAW TRACE BODY");
            };
        }
    }
}

int drawFrame(SDL_Renderer * renderer, TTF_TextEngine* text_engine, TTF_Font** fonts, body* bodies, traceBody** traces){
    // Draw bacground colours
    if (!(SDL_SetRenderDrawColor(renderer, 20, 30, 60, 255) && SDL_RenderClear(renderer))){
        printf("SDL_SetRenderDrawColor || SDL_RenderClear didn't work: %s\n", SDL_GetError());
        return -1; 
    }

    renderTraces(renderer, traces);
    renderBodies(renderer, text_engine, fonts, bodies);
    
    // Draw Walls
    if (!(thickLineColor(renderer, 0, 0, WIDTH + WALL_THICKNESS/2, 0, WALL_THICKNESS, WALL_COLOUR) + thickLineColor(renderer, WIDTH, 0, WIDTH, HEIGHT+ WALL_THICKNESS/2, WALL_THICKNESS, WALL_COLOUR) + thickLineColor(renderer, 0, HEIGHT, WIDTH + WALL_THICKNESS/2, HEIGHT, WALL_THICKNESS, WALL_COLOUR) + thickLineColor(renderer, 0, 0, 0, HEIGHT+WALL_THICKNESS/2, WALL_THICKNESS, WALL_COLOUR) == 4)){
        printf("SDL_RenderLine didn't work: %s\n", SDL_GetError());
        return -1;
    }        



    return 0;
}

void cleanBodies(body* bodies){
    free(bodies);
}

void cleanTraces(traceBody** traces){
    for (int i = 0; i < NUM_BODIES; i++){
        free(traces[i]);
    }
    free(traces);
}

void closeFonts(TTF_Font** fonts){
    for(int f = 0; f< NUM_BODIES; f++){
        TTF_CloseFont(fonts[f]);
    }
    free(fonts);
}

int main(int argc, char* argv[]){
    int kk = 0;
    int jk = 0 ;
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
    TTF_Font** fonts = NULL;

    if (SHOWTEXT){
        text_engine = TTF_CreateRendererTextEngine(renderer);
        fonts = malloc(NUM_BODIES*sizeof(TTF_Font*));
    }
    
    if (SHOWTEXT && text_engine == NULL){
        printf("TTF_CreateRendererTextEngine Renderer didn't work: %s\n", SDL_GetError());
        SDL_DestroyRenderer(renderer);
        SDL_DestroyWindow(window);
        SDL_Quit();
        
        return -1;
    }
    if (SHOWTEXT && fonts == NULL){
        printf("Malloc did not work for Fonts\n");

        TTF_DestroyRendererTextEngine(text_engine);
        SDL_DestroyRenderer(renderer);
        SDL_DestroyWindow(window);
        SDL_Quit();
        
        return -1;
    }

    body* bodies = malloc(NUM_BODIES * sizeof(body));
    setupBodies(bodies, NUM_BODIES);

    traceBody** traces = malloc(NUM_BODIES * sizeof(*traces));
    setupBodyTraces(traces, bodies);


    for (int f = 0; (f < NUM_BODIES && SHOWTEXT); f++){
        fonts[f] = TTF_OpenFont("./src/fonts/dejavu-sans-condensed-webfont.ttf", 4 * bodies[f].radius/15);
        TTF_SetFontStyle(fonts[f], TTF_STYLE_BOLD);
    } 

 

    SDL_Event windowEvent;

    uint32_t framteRate = ((uint32_t)(EPSILON*1000.0) > 0)? (uint32_t)(EPSILON*1000.0) : 1 ;

    while (true)
    {
        if(SDL_PollEvent(&windowEvent)){
            if (windowEvent.type == SDL_EVENT_QUIT){
                break;
            }
        }
        for (int jjk = 0; jjk < ((((1.0/FRAMERATE) )/EPSILON) + 1); jjk++){
            if (kk == NUM_TRACES){
            // 1. Update Traces
                updateTraces(traces, bodies);
                kk = 0;
            }

            kk += 1;
            // 2. Update Bodies
            updateBodies(bodies, NUM_BODIES);

        }
        // 3. Draw the frame
        if (drawFrame(renderer, text_engine, fonts, bodies, traces) != 0){
            
            closeFonts(fonts);
            TTF_DestroyRendererTextEngine(text_engine);
            SDL_DestroyRenderer(renderer);
            SDL_DestroyWindow(window);
            SDL_Quit();
            cleanBodies(bodies);
            cleanTraces(traces);
            return -1;
        }

        SDL_RenderPresent(renderer);
        SDL_Delay((int)((1.0/FRAMERATE) * 1000) - 10);
    }

    // const Uint64 FRAME_NS      = SDL_NS_PER_SECOND / FRAMERATE;  // target time per rendered frame
    // const double MAX_FRAME_DT  = 0.25;   // seconds; cap to avoid a death spiral (see below)
    // const int    TRACE_INTERVAL = 100;   // physics steps between trace snapshots

    // Uint64 previous    = SDL_GetTicksNS();
    // double accumulator = 0.0;
    // int    steps_since_trace = 0;
    // bool   running     = true;

    // while (running) {
    //     Uint64 frame_start = SDL_GetTicksNS();

    //     // Handle every pending event, not just one
    //     while (SDL_PollEvent(&windowEvent)) {
    //         if (windowEvent.type == SDL_EVENT_QUIT) running = false;
    //     }

    //     // How much real time passed since the last frame?
    //     double elapsed = (double)(frame_start - previous) / SDL_NS_PER_SECOND;
    //     previous = frame_start;
    //     if (elapsed > MAX_FRAME_DT) elapsed = MAX_FRAME_DT;
    //     accumulator += elapsed;

    //     // Run exactly as many fixed steps as that real time covers
    //     while (accumulator >= EPSILON) {
    //         updateBodies(bodies, NUM_BODIES);

    //         if (++steps_since_trace >= TRACE_INTERVAL) {
    //             updateTraces(traces, bodies);
    //             steps_since_trace = 0;
    //         }
    //         accumulator -= EPSILON;
    //     }

    //     if (drawFrame(renderer, text_engine, fonts, bodies, traces) != 0) {
    //         running = false;           // fall through to the normal cleanup below
    //         break;
    //     }
    //     SDL_RenderPresent(renderer);

    //     // Sleep only for whatever is left of this frame's time budget
    //     Uint64 frame_time = SDL_GetTicksNS() - frame_start;
    //     if (frame_time < FRAME_NS) {
    //         SDL_DelayNS(FRAME_NS - frame_time);
    //     }
    // }
    cleanBodies(bodies);
    cleanTraces(traces);

    closeFonts(fonts);
    TTF_DestroyRendererTextEngine(text_engine);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    TTF_Quit();
    SDL_Quit();
    
    return SDL_EVENT_QUIT;
}