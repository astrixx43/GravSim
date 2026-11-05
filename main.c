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

const double C =  0.1*1920.0 + 0 * 299792458.0;
const double C_SQUARED = C*C;
const double ONE_C = 1/C;
const double ONE_C_SQUARED = 1/C_SQUARED;

double one_squareroot_C;

double G =  10.0 * 6.67430 * 0.000000000001;

const double MIN_SPEED = 0.1; // The minimum speed of body i of C. MUST BE << C
const double MAX_SPEED = 0.2; // The minimum speed of body i of C. MUST BE << C


const int NUM_BODIES = 3;
const int NUM_TRACES = 10; // Set KK equal to this and set this to 4-6 give good lookinfg results

const double MIN_RADIUS = 10;
const double MIN_RADIUS_SQUARED = MIN_RADIUS*MIN_RADIUS;

int indexOfTrace = NUM_TRACES - 1;
void (*grav_function)(body*);
const int WALL_THICKNESS = 10;
const uint32_t WALL_COLOUR = 0X0F00F000;


const int WIDTH = 1920;
const int HEIGHT = 1080;
double DIAGONAL_SQUARED = (WIDTH*WIDTH + HEIGHT*HEIGHT);
double EPSILON = 0.001; // simmulation step size. 1 == a second

bool SHOWTEXT = false;

uint64_t FRAMERATE = 60;


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
    return Q_rsqrt((double)( 1 - ((speed*speed)*ONE_C_SQUARED)));
}

double calculateGamma_alt(double speed_squared){
    return Q_rsqrt((double)( 1 - ((speed_squared)*ONE_C_SQUARED)));
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
    
    body->momentum[0] = v_x * body->mass;
    body->momentum[1] = v_y * body->mass;

    body->speed = pow((body->velocity[0]*body->velocity[0]) + (body->velocity[1]*body->velocity[1]), 1.0/2.0);

}

void boost(double *E, double* p_x, double* p_y, double V_x, double V_y, double V_squared, double gamma_V) {

    if (V_squared == 0.0) return;

    double Vp = V_x * (*p_x) + V_y * (*p_y);
    double k  = (gamma_V - 1.0) * Vp / V_squared - gamma_V * (*E) / (C_SQUARED);

    *E = gamma_V * (*E - Vp);

    *p_x += k * V_x;
    *p_y += k * V_y;
}

void relativisticMassiveCollison(body* b1, body* b2){
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

void relativisticGravitonOnMassiveCollison(body* b1, body* masslessBody){

    double gamma_1 = calculateGamma(b1->speed);

    double E_1 = gamma_1*b1->mass*C_SQUARED;
    double E_2 = pow((masslessBody->momentum[0] * masslessBody->momentum[0]) + (masslessBody->momentum[1]*masslessBody->momentum[1]), 1.0/2.0) * C;

    double p_1_x = gamma_1*b1->mass*b1->velocity[0];
    double p_1_y = gamma_1*b1->mass*b1->velocity[1];

    double p_2_x = masslessBody->momentum[0];
    double p_2_y = masslessBody->momentum[1];

    double One_E_total = 1 / (E_1 + E_2);
    
    double gamma_V = gamma_1 * (1 + (E_2/E_1));

    double V_x = (p_1_x + p_2_x)/(gamma_V * b1->mass);
    double V_y = (p_1_y + p_2_y)/(gamma_V * b1->mass);

    change_body_velocity(b1, V_x, V_y);

    // boost(&E_1, &p_1_x, &p_1_y, V_x, V_y, V_sqaured, gamma_V);
    // // boost(&E_2, &p_2_x, &p_2_y, V_x, V_y, V_sqaured, gamma_V);

    // double one_n = Q_rsqrt((b1->pos[0] - masslessBody->pos[0]) * (b1->pos[0] - masslessBody->pos[0]) + (b1->pos[1] - masslessBody->pos[1]) * (b1->pos[1] - masslessBody->pos[1]));
    // double n_x = one_n * (masslessBody->pos[0] - b1->pos[0]);
    // double n_y = one_n * (masslessBody->pos[1] - b1->pos[1]);
    // double pn = p_1_x* n_x + p_1_y*n_y;

    // // if (pn <= 0.0) return;  

    // p_1_x -= 2 * pn * n_x;
    // p_1_y -= 2 * pn * n_y;
    
    // // p_2_x += 2 * pn * n_x;  
    // // p_2_y += 2 * pn * n_y;


    // boost(&E_1, &p_1_x, &p_1_y, -1.0 * V_x, -1.0 * V_y, V_sqaured, gamma_V);
    // // boost(&E_2, &p_2_x, &p_2_y, -1.0 * V_x, -1.0 * V_y, V_sqaured, gamma_V);

    // change_body_velocity(b1, p_1_x * C_SQUARED / E_1,  p_1_y * C_SQUARED / E_1);
    // // change_body_velocity(masslessBody, p_2_x * C_SQUARED / E_2, p_2_y * C_SQUARED / E_2);


}

int checkBodyCollision(body* b1, body*b2){
    double distance = pow((b1->pos[0] - b2->pos[0])*(b1->pos[0] - b2->pos[0]) + (b1->pos[1] - b2->pos[1])*(b1->pos[1] - b2->pos[1]), 1.0/2.0);
    return distance <= (b1->radius + b2->radius);
}

int checkWaveBodyCollision(body* wave, body* b1){
    double distance = pow((b1->pos[0] - wave->pos[0])*(b1->pos[0] - wave->pos[0]) + (b1->pos[1] - wave->pos[1])*(b1->pos[1] - wave->pos[1]), 1.0/2.0);
    return (distance < (b1->radius + wave->radius + C*EPSILON)) && distance > (b1->radius + wave->radius - C*EPSILON);
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

void handelCollision(body* bodies){
    for (int i = 0; i < NUM_BODIES; i++){
        if (!checkWallCollision(&bodies[i])){
            for (int j = i + 1; j < NUM_BODIES; j++){
                if (checkBodyCollision(&bodies[i], &bodies[j])){
                    relativisticMassiveCollison(&bodies[i], &bodies[j]);
                }
            }
        }
    }

}

void makeNewBody(body* body){
    double radius = ((((double)rand())/RAND_MAX) * 10) + MIN_RADIUS;
    
    double x_pos = ((double)rand()/(double)RAND_MAX) * (double)(WIDTH - radius) + (double)radius;
    double y_pos = ((double)rand()/(double)RAND_MAX) * (double)(HEIGHT - radius) + (double)radius;
    double theta = ((double)rand()/(double)RAND_MAX) * 2.0f * (double)(M_PI);

    double v_total = 0; //(((double)rand()/(double)RAND_MAX) * (MAX_SPEED- MIN_SPEED) +  MIN_SPEED) * C;
    double v_x =  v_total * cos(theta);
    double v_y = v_total * sin(theta);

    uint64_t colour = wavelengthToHex((uint64_t)380 + (uint64_t)rand() % (700- 380));
    // uint32_t colour = wavelengthToHex(580);
    

    body_init(body, 100000000000.0 * M_PI*(radius*radius), x_pos, y_pos, v_x, v_y, radius, colour);
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

void shootGravitons(body* b_emmitor, body* b_reciver){
    gravitonList* gList = (gravitonList*)b_emmitor->gravitonList;


    gravitonNode* graviton = (gList->last == NULL)? gList->first : gList->last;
    // double delT = 0;
    double one_over_dist = Q_rsqrt((b_emmitor->pos[0] - b_reciver->pos[0])*(b_emmitor->pos[0] - b_reciver->pos[0]) + (b_emmitor->pos[1] - b_reciver->pos[1])*(b_emmitor->pos[1] - b_reciver->pos[1]));
    
    // double dist_squared = (b_emmitor->pos[0] - b_reciver->pos[0])*(b_emmitor->pos[0] - b_reciver->pos[0]) + (b_emmitor->pos[1] - b_reciver->pos[1])*(b_emmitor->pos[1] - b_reciver->pos[1]);
    
    graviton->node = malloc(sizeof(body));

    double x_comp = (b_reciver->pos[0]  - b_emmitor->pos[0]) * one_over_dist;
    double y_comp =  (b_reciver->pos[1] - b_emmitor->pos[1]) * one_over_dist;
    
    body_init(graviton->node, 0, b_emmitor->pos[0] + (b_emmitor->radius * x_comp), b_emmitor->pos[1] + (b_emmitor->radius * y_comp), C*x_comp, C*y_comp, 0, 0);
    
    graviton->node->momentum[0] = b_emmitor->pos[0];
    graviton->node->momentum[1] = b_emmitor->pos[1];

    graviton->next = malloc(sizeof(gravitonNode));
    gravitonNode_init(graviton->next);

    graviton = graviton->next;


    
    double delta_width = C;

    double angle_gamma = atan2((delta_width + b_emmitor->radius), 1/one_over_dist);

    double phi = atan2((b_reciver->pos[0]  - b_emmitor->pos[0]), (b_reciver->pos[1]  - b_emmitor->pos[1]));
    
    // double del_theta = (theta_step > .5)? theta_step : min_theta;
    // del_theta = del_theta/2;

    // printf("%f\n", theta_step);

    double del_theta = 0.001;
    for (double theta = -angle_gamma; theta < 0 ; theta += .05-del_theta){

        graviton->node = malloc(sizeof(body));

        double x_comp = sin(phi - theta); 
        double y_comp = cos(phi - theta);       

        body_init(graviton->node, 0, b_emmitor->pos[0] + (b_emmitor->radius * x_comp), b_emmitor->pos[1] + (b_emmitor->radius * y_comp), C*x_comp, C*y_comp, 0, 0);
        
        graviton->node->momentum[0] = b_emmitor->pos[0];
        graviton->node->momentum[1] = b_emmitor->pos[1];
        
        graviton->next = malloc(sizeof(gravitonNode));
        gravitonNode_init(graviton->next);
        graviton = graviton->next;
        del_theta = del_theta * 3.0/4.0;
    }
    
    double old_del_theta = .05+del_theta;
    del_theta = 0.;
    for (double theta = 0; theta < angle_gamma ; theta += old_del_theta-del_theta){

        graviton->node = malloc(sizeof(body));

        double x_comp = sin(phi - theta); 
        double y_comp = cos(phi - theta);       

        body_init(graviton->node, 0, b_emmitor->pos[0] + (b_emmitor->radius * x_comp), b_emmitor->pos[1] + (b_emmitor->radius * y_comp), C*x_comp, C*y_comp, 0, 0);
        
        graviton->node->momentum[0] = b_emmitor->pos[0];
        graviton->node->momentum[1] = b_emmitor->pos[1];
        
        graviton->next = malloc(sizeof(gravitonNode));
        gravitonNode_init(graviton->next);
        graviton = graviton->next;
        del_theta = del_theta * 4.0/3.0;
    }

    gList->last = graviton;

    // for (double theta = 0; theta < 2 * M_PI; theta+= .1){
   
    //     graviton->node = malloc(sizeof(body));

    //     double x_comp = sin(theta);
    //     double y_comp =  cos(theta);
        
    //     body_init(graviton->node, 0, b_emmitor->pos[0] + (b_emmitor->radius * x_comp), b_emmitor->pos[1] + (b_emmitor->radius * y_comp), C*x_comp, C*y_comp, 0, 0);
        
    //     graviton->node->momentum[0] = b_emmitor->pos[0];
    //     graviton->node->momentum[1] = b_emmitor->pos[1];

    //     graviton->next = malloc(sizeof(gravitonNode));
    //     gravitonNode_init(graviton->next);
    //     graviton = graviton->next;
    // }

   

    // graviton->next = malloc(sizeof(gravitonNode));
    // gravitonNode_init(graviton->next);
    // graviton = graviton->next;


    // direct_theta = (direct_theta > 0) ? direct_theta : acos(y_comp);
    // printf("%f\n", direct_theta);
    
    // double min_theta = acos(1 - (MIN_RADIUS_SQUARED*2.0/DIAGONAL_SQUARED));
}

void emmitGravitons(body* bodies){
    for(int j = 0; j < NUM_BODIES; j++){
        body* b_emmitor = &bodies[j];
    
        for(int i = 0; i < NUM_BODIES; i++){
            if (&bodies[i] != b_emmitor){
                shootGravitons(b_emmitor, &bodies[i]);
            }
        }
    }
    // double step_size = M_PI/1000;
    // for (double theta = 0; theta < 2 * M_PI; theta +=step_size){
    //     graviton->node = malloc(sizeof(body));

    //     double x_comp = sin(theta);
    //     double y_comp = cos(theta);
        
    //     body_init(graviton->node, 0, b_emmitor->pos[0] + (b_emmitor->radius * x_comp), b_emmitor->pos[1] + (b_emmitor->radius * y_comp), C*x_comp, C*y_comp, 0, 0);

    //     graviton->next = malloc(sizeof(gravitonNode));
    //     gravitonNode_init(graviton->next);
    //     graviton = graviton->next;
    // }

}

void updateGravitons(body* bodies){

    for (int i = 0; i< NUM_BODIES; i++){
        
        gravitonNode* graviton = ((gravitonList*)bodies[i].gravitonList)->first;

        while (graviton->node != NULL)
        {
            if (checkWallCollision(graviton->node)){
                free(graviton->node);
                
                if (((gravitonNode*)graviton->next)->node != NULL){
                    if(graviton->prev == NULL){
                        ((gravitonList*)bodies[i].gravitonList)->first = graviton->next;
                    }else{
                        ((gravitonNode*)graviton->prev)->next = graviton->next;
                    }
                    ((gravitonNode*)graviton->next)->prev = graviton->prev;
                    gravitonNode* old_node = graviton;
                    
                    graviton = graviton->next;
                    free(old_node);


                }else{

                    graviton->node = NULL;
                    free(graviton->next);
                    graviton->next = NULL;

                    ((gravitonList*)bodies[i].gravitonList)->last = graviton;
                }

            }else{

                graviton->node->pos[0] = graviton->node->pos[0] + EPSILON*graviton->node->velocity[0];
                graviton->node->pos[1] = graviton->node->pos[1] + EPSILON*graviton->node->velocity[1];
                graviton = graviton->next;
            }
        }

    }

}

void handleGravitonBodyCollision(body* b1, body* b_emmiter){

    gravitonNode* gravitons = ((gravitonList*)b_emmiter->gravitonList)->first;

    int collided = 0;
    while (gravitons->node != NULL){
        if (checkBodyCollision(b1, gravitons->node)){
            if (!collided){
                body* graviton_body = gravitons->node;
                
                double distance_sqrd = ((graviton_body->momentum[0]-b1->pos[0])*(graviton_body->momentum[0]-b1->pos[0])) + ((graviton_body->momentum[1] - b1->pos[1])*(graviton_body->momentum[1] - b1->pos[1]));
                
                double gravitational_momentum_coefficeint = G * b1->mass * b_emmiter->mass/ (distance_sqrd * C);

                graviton_body->momentum[0] = gravitational_momentum_coefficeint * (graviton_body->momentum[0]-b1->pos[0]);
                graviton_body->momentum[1] = gravitational_momentum_coefficeint * (graviton_body->momentum[1]-b1->pos[1]);

                relativisticGravitonOnMassiveCollison(b1, graviton_body);
                collided = 1;
                
            }
            
            free(gravitons->node);
                

            if (((gravitonNode*)gravitons->next)->node != NULL){

                gravitons->node = ((gravitonNode*)(gravitons->next))->node;
                gravitons->next = ((gravitonNode*)(gravitons->next))->next;

            }else{
                gravitons->node = NULL;

                free(gravitons->next);
                gravitons->next = NULL;
                ((gravitonList*)b_emmiter->gravitonList)->last = gravitons;
            }

        }else{
            gravitons = gravitons->next;
        }
    
    }
}

void particleGravitonsWithRelevCollision(body* bodies){
    for (int i = 0; i < NUM_BODIES; i++){
        for(int j = 0; j < NUM_BODIES; j++){
            if (i != j){
                handleGravitonBodyCollision(&bodies[i], &bodies[j]);
            }
        }    
    }
    updateGravitons(bodies);
}

void emmitWaveGravitons(body* bodies){
    double corner[2];

    for(int i = 0; i < NUM_BODIES; i++){
        body* b_emmitor = &bodies[i];

        gravitonNode* gravitons;
        if (((gravitonList*)b_emmitor->gravitonList)->last == NULL){
            gravitons = (gravitonNode*)((gravitonList*)b_emmitor->gravitonList)->first; 
        }else{
            gravitons = (gravitonNode*)((gravitonList*)b_emmitor->gravitonList)->last; 
        }
        
        gravitons->node = malloc(sizeof(body));

        body_init(gravitons->node, 0, bodies[i].pos[0], bodies[i].pos[1], 0, 0, bodies[i].radius, 0X000000FF); 
        


        if ((gravitons->node->pos[0] - WIDTH/2.0) < 0 ){
            corner[0] = WIDTH;
        }else{
            corner[0] = 0;

        }
        
        if ((gravitons->node->pos[1] - HEIGHT/2.0) < 0){                
            corner[1] = HEIGHT;
        }else{
            corner[1] = 0;
        }

        gravitons->max_distance = pow((corner[0] - gravitons->node->pos[0])*(corner[0] - gravitons->node->pos[0]) + (corner[1] - gravitons->node->pos[1])*(corner[1] - gravitons->node->pos[1]), 1/2.0);
        

        gravitons->next = malloc(sizeof(gravitonNode));
        
        gravitonNode_init(gravitons->next);
        ((gravitonNode*)gravitons->next)->prev = (void*)gravitons;

        ((gravitonList*)b_emmitor->gravitonList)->last = gravitons->next;
        
    }

}

void handleWaveGravitonBodyCollision(body* b1, body* b_emmiter){

    gravitonNode* gravitons = ((gravitonList*)b_emmiter->gravitonList)->first;


    while (gravitons->node != NULL){
        if (checkWaveBodyCollision(b1, gravitons->node)){

            body* graviton_body = gravitons->node;
            
            double distance_sqrd = ((graviton_body->pos[0]-b1->pos[0])*(graviton_body->pos[0]-b1->pos[0])) + ((graviton_body->pos[1] - b1->pos[1])*(graviton_body->pos[1] - b1->pos[1]));
            
            double gravitational_momentum_coefficeint = G * b1->mass * b_emmiter->mass/ (distance_sqrd * C);

            graviton_body->momentum[0] = gravitational_momentum_coefficeint * (graviton_body->pos[0]-b1->pos[0]);
            graviton_body->momentum[1] = gravitational_momentum_coefficeint * (graviton_body->pos[1]-b1->pos[1]);

            relativisticGravitonOnMassiveCollison(b1, graviton_body);

            break;
        }else{
            gravitons = gravitons->next;
        }
    
    }
}

void updateWaveGravitons(body* bodies){

    for (int i = 0; i< NUM_BODIES; i++){
        
        gravitonNode* graviton = ((gravitonList*)bodies[i].gravitonList)->first;

        while (graviton != NULL && graviton->node != NULL)
        {

            if (graviton->node->radius >= graviton->max_distance){
                free(graviton->node);
                
                if (((gravitonNode*)graviton->next)->node != NULL){
                    
                    if(graviton->prev == NULL){
                        ((gravitonList*)bodies[i].gravitonList)->first = graviton->next;
                    }else{
                        ((gravitonNode*)graviton->prev)->next = graviton->next;
                    }

                    ((gravitonNode*)graviton->next)->prev = graviton->prev;
                    gravitonNode* old_node = graviton;
                    graviton = graviton->next;
                    free(old_node);

                }else{
                    graviton->node = NULL;

                    free(graviton->next);

                    graviton->next = NULL;
                    graviton->max_distance = 0;
                    ((gravitonList*)bodies[i].gravitonList)->last = graviton;

                }

            }else{

                graviton->node->radius = graviton->node->radius + EPSILON*C;
                graviton = graviton->next;
            }
        }

    }

}

void waveGravitonsWithRelevCollision(body* bodies){
    for (int i = 0; i < NUM_BODIES; i++){
        for(int j = 0; j < NUM_BODIES; j++){
            if (i != j){
                handleWaveGravitonBodyCollision(&bodies[i], &bodies[j]);
            }
        }    
    }
    updateWaveGravitons(bodies);
}

void relativisticEIHGravity(body* bodies){
    double r_s[NUM_BODIES*NUM_BODIES];
    double one_r_s_sqrd[NUM_BODIES*NUM_BODIES];

    for (int i = 0; i<NUM_BODIES; i++){
        for (int j = i + 1; j< NUM_BODIES; j++){
            r_s[i] = pow((bodies[i].pos[1] - bodies[j].pos[1])*(bodies[i].pos[1] - bodies[j].pos[1]) + (bodies[i].pos[0] - bodies[j].pos[0])*(bodies[i].pos[0] - bodies[j].pos[0]), 1.0/2.0);
            r_s[i+ j * NUM_BODIES] = r_s[i];
            one_r_s_sqrd[i] = 1.0/r_s[i];
            one_r_s_sqrd[i+ j * NUM_BODIES] = 1.0/r_s[i];
        }
    }

    for (int i = 0; i<NUM_BODIES; i++){



    }
    
}

void updateBodies(body* bodies){
    
    handelCollision(bodies);

    grav_function(bodies);


    for (int i = 0; i < NUM_BODIES; i++){
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

        if (!filledCircleRGBA(renderer, bodies[i].pos[0], bodies[i].pos[1], (Sint16)bodies[i].radius, (c >> 24) & 0xFF, (c >> 16) & 0xFF, (c >> 8) & 0xFF, c & 0xFF)){
            printf("FAILED TO DRAW BODY");
        };
                
        if (SHOWTEXT){
            snprintf(body_str, sizeof(body_str), "(%.2fC, %.2fC)\0", bodies[i].velocity[0]*ONE_C, bodies[i].velocity[1]*ONE_C);
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

void renderGravitons(SDL_Renderer* renderer, body* bodies){
    int NUM_RINGS_SKIP = 100;
    for (int i = 0; i< NUM_BODIES; i ++){
        gravitonNode* graviton = ((gravitonList*) bodies[i].gravitonList)->first;
        int j = 0;
        
        while (j < NUM_RINGS_SKIP && graviton != NULL && graviton->node != NULL){
            graviton = graviton->next;
            j += 1;
        }
        
        while(graviton != NULL && graviton->node != NULL){
            j = 0;

            if (!circleRGBA(renderer, graviton->node->pos[0], graviton->node->pos[1], (Sint16)graviton->node->radius+MIN_RADIUS*EPSILON, 00, 00, 00, 0XFF)){
                printf("FAILED TO DRAW GRAVITON");
                break;
            };
            while (j < NUM_RINGS_SKIP && graviton != NULL && graviton->node != NULL){
                graviton = graviton->next;
                j += 1;
            }
            //graviton = graviton->next;
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
    // renderGravitons(renderer, bodies);
    
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
    int fontSize = TTF_GetFontSize(font);
    TTF_Text* body_text = TTF_CreateText(text_engine, font, body_str, 0); 
    TTF_SetTextColor(body_text, 0xFF, 0x0F, 0xFF, 0xFF);
    
    if (!TTF_DrawRendererText(body_text, 0, 0)){
        printf("FAILED TO DRAW TEXT: %s\n", SDL_GetError());
        return -1;
    }
    TTF_DestroyText(body_text);
    return 1;
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
    for(int f = 0; f< NUM_BODIES + 1; f++){
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
        fonts = malloc((NUM_BODIES+1)*sizeof(TTF_Font*));
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

    grav_function = particleGravitonsWithRelevCollision;
    // grav_function = waveGravitonsWithRelevCollision;
    one_squareroot_C = Q_rsqrt(C);

    body* bodies = malloc(NUM_BODIES * sizeof(body));
    setupBodies(bodies, NUM_BODIES);

    traceBody** traces = malloc(NUM_BODIES * sizeof(*traces));
    setupBodyTraces(traces, bodies);


    for (int f = 0; (f < NUM_BODIES && SHOWTEXT); f++){
        fonts[f] = TTF_OpenFont("./src/fonts/dejavu-sans-condensed-webfont.ttf", 4 * bodies[f].radius/15);
        TTF_SetFontStyle(fonts[f], TTF_STYLE_BOLD);
    } 

    if (SHOWTEXT){
        fonts[NUM_BODIES] = TTF_OpenFont("./src/fonts/dejavu-sans-condensed-webfont.ttf", 20);
        TTF_SetFontStyle(fonts[NUM_BODIES], TTF_STYLE_BOLD);
    }
 

    SDL_Event windowEvent;

    uint32_t framteRate = ((uint32_t)(EPSILON*1000.0) > 0)? (uint32_t)(EPSILON*1000.0) : 1 ;
    
    double stepsPerSecond = 1.0/(EPSILON);
    double framePerSecond = FRAMERATE;

    double stepPerFrame = stepsPerSecond/(framePerSecond);

    double nanosecondsPerFrame = 1000000000/framePerSecond;
    
    int mm = 0;
    int ll = 0;
    int numberOfStepsPassedSinceLastFrame = 0;
    Uint64 time_since_last_draw = SDL_GetTicksNS();
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

        // while(SDL_PollEvent(&windowEvent)){
        //     if (windowEvent.type == SDL_EVENT_QUIT){
        //         break;
        //     }
        // }
        
        if (SDL_PollEvent(&windowEvent) && windowEvent.type == SDL_EVENT_QUIT){
            break;
        }
        
        while (frame_time_accumilator > EP_NS){
            if (kk == NUM_TRACES){
            // 1. Update Traces
                updateTraces(traces, bodies);
                kk = 0;
            }

            kk += 1;
            // 2. Update Bodies
            updateBodies(bodies);
            
            mm += 1;
            if (grav_function == particleGravitonsWithRelevCollision && mm >= (MIN_RADIUS/(C*EPSILON))){
                emmitGravitons(bodies);
                mm = 0;
            }else if (grav_function == waveGravitonsWithRelevCollision && mm >= (MIN_RADIUS/(C*EPSILON))){
                emmitWaveGravitons(bodies);
                mm = 0;
            }
            numberOfStepsPassedSinceLastFrame += 1;
            frame_time_accumilator -= EP_NS;
        }

        // 3. Draw the frame

        ll += 1;

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

        // if ( SHOWTEXT && !drawTextScreen(text_engine, fonts[NUM_BODIES])){
        //     closeFonts(fonts);
        //     TTF_DestroyRendererTextEngine(text_engine);
        //     SDL_DestroyRenderer(renderer);
        //     SDL_DestroyWindow(window);
        //     SDL_Quit();
        //     cleanBodies(bodies);
        //     cleanTraces(traces);
        //     return -1;
        // }

        SDL_RenderPresent(renderer);

        SDL_DelayPrecise((nanosecondsPerFrame - (SDL_GetTicksNS() - current_time) > 0 )? (nanosecondsPerFrame - (SDL_GetTicksNS() - current_time)) : 0);
        numberOfStepsPassedSinceLastFrame = 0;
    }

    cleanBodies(bodies);
    cleanTraces(traces);
    if (SHOWTEXT){
        closeFonts(fonts);
        TTF_DestroyRendererTextEngine(text_engine);
    }
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    TTF_Quit();
    SDL_Quit();
    
    return 0;
}