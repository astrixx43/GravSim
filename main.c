#include <stdlib.h>
#include <stdio.h>
#include <SDL3/SDL.h>
#include "SDL3_gfxPrimitives.h"
#include "SDL3_framerate.h"
#include "bodies.h"
#include "math.h"
#include "Wavelength_Rgb.h"

int WIDTH = 1600;
int HEIGHT = 1200;
double EPSILON = (double)0.01;
const double C = 100.0;
const double C_SQUARED = C*C;

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
    if ((b->pos[0] - b->radius <= 0     && b->velocity[0] < 0) ||
        (b->pos[0] + b->radius >= WIDTH && b->velocity[0] > 0)) {
        change_body_velocity(b, -b->velocity[0], b->velocity[1]);
        hit = 1;
    }
    if ((b->pos[1] - b->radius <= 0      && b->velocity[1] < 0) ||
        (b->pos[1] + b->radius >= HEIGHT && b->velocity[1] > 0)) {
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
    double x_pos = ((double)rand()/(double)RAND_MAX) * (double)WIDTH;
    double y_pos = ((double)rand()/(double)RAND_MAX) * (double)HEIGHT;
    double theta = ((double)rand()/(double)RAND_MAX) * 2.0f * (double)(M_PI);
    double v_total = ((double)rand()/(double)RAND_MAX) * 0.90 *C + 0.01 *C ;
    // v_total = 0.5 * C;
    double v_x =  v_total * cos(theta);
    double v_y = v_total * sin(theta);

    uint64_t colour = wavelengthToHex((uint64_t)380 + (uint64_t)rand() % (700- 380));
    colour = 0XFF000000;
    //    uint64_t colour = wavelengthToHex(380.0 + ((double)rand()/(double)RAND_MAX)*(700.0 - 380.0));


    body_init(body, 100, x_pos, y_pos, v_x, v_y, 20, colour);
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

void updateBodies(body* bodies, int num){
    
    handelCollision(bodies, num);

    for (int i = 0; i < num; i++){
        bodies[i].pos[0] = bodies[i].pos[0] + EPSILON*bodies[i].velocity[0];
        bodies[i].pos[1] = bodies[i].pos[1] + EPSILON*bodies[i].velocity[1];
        relativisticDopplerEffect(&bodies[i]);
    }
}

void renderBodies(SDL_Renderer* renderer, body* bodies, int num){
    Uint32 c = 0;
    for (int i = 0; i < num; i++){
        c = bodies[i].relative_colour;
        filledCircleRGBA(renderer, bodies[i].pos[0], bodies[i].pos[1], bodies[i].radius, (c >> 24) & 0xFF, (c >> 16) & 0xFF, (c >> 8) & 0xFF, c & 0xFF);
    }
}

int main(int argc, char* argv[]){
    
    if (!SDL_Init(SDL_INIT_VIDEO)) {
        printf("SDL_Init failed: %s\n", SDL_GetError());
        return 1;
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

    int num_bodies = 10;

    body* bodies = malloc(num_bodies * sizeof(body));
    setupBodies(bodies, num_bodies);


    SDL_Event windowEvent;

    
    while (true)
    {
        if(SDL_PollEvent(&windowEvent)){
            if (windowEvent.type == SDL_EVENT_QUIT){
                break;
            }
        }
        // 1. Background: dark blue (R, G, B, A — each 0 to 255)
        SDL_SetRenderDrawColor(renderer, 20, 30, 60, 255);
        SDL_RenderClear(renderer);

        // 2. Update Bodies
        updateBodies(bodies, num_bodies);


        // 3. Show the frame
        renderBodies(renderer, bodies, num_bodies);
        SDL_RenderPresent(renderer);
        SDL_Delay(1);
    }

    SDL_DestroyWindow(window);
    SDL_Quit();
    
    return SDL_EVENT_QUIT;
}