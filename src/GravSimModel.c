#include <stdlib.h>
#include <stdio.h>
#include <xmmintrin.h>
#include <string.h>
#include <math.h>

#include "GravSimModel.h"
#include "Wavelength_Rgb.h"

double C_speed;
double C_SQUARED;
double ONE_C_SQUARED;
double GRAVITATIONAL_CONSTANT;
double EPSILON_STEP;
int BOX_WALL_THICKNESS;
int WIDTH_OF_BOX;
int HEIGHT_OF_BOX;
double GRAVITON_PULSE_PERIOD;


double _MIN_SPEED; // The minimum speed of body i of C_speed. MUST BE << C_speed
double _MAX_SPEED;// The maximum speed of body i of C_speed. MUST BE > _MIN_SPEED
double _MIN_RADIUS;
double _MAX_RADIUS;


float fast_inverse_sqrt(float val) {
    // Pack the float into an XMM register (__m128)
    __m128 input = _mm_set_ss(val);
    
    // Execute rsqrtss approximation on the low 32 bits
    __m128 result = _mm_rsqrt_ss(input);
    
    // Extract the float result from the XMM register
    return _mm_cvtss_f32(result);
}


double Q_rsqrt( double number )
{   
    return fast_inverse_sqrt(number);
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
    double r_x = body->pos[0] - (WIDTH_OF_BOX / 2.0);
    double r_y = body->pos[1] - (HEIGHT_OF_BOX / 2.0);

    // v · r : positive = moving away from the centre, negative = moving toward it
    double radial = body->velocity[0] * r_x + body->velocity[1] * r_y;

    double speed = body->speed;
    double beta = speed / C_speed;
    if (beta >= 1.0) beta = 0.999999;  // guard against division by zero

    double one_r_len = Q_rsqrt(r_x * r_x + r_y * r_y);
    double v_r = (r_x != r_y || r_y != 0.0) ? radial * one_r_len : 0.0;   // velocity component away from centre
    double factor = (1.0 + v_r / C_speed) * Q_rsqrt(1.0 - beta * beta);

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
    double E_2 = pow((masslessBody->momentum[0] * masslessBody->momentum[0]) + (masslessBody->momentum[1]*masslessBody->momentum[1]), 1.0/2.0) * C_speed;

    double p_1_x = gamma_1*b1->mass*b1->velocity[0];
    double p_1_y = gamma_1*b1->mass*b1->velocity[1];

    double p_2_x = masslessBody->momentum[0];
    double p_2_y = masslessBody->momentum[1];
    
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
    double distance_sqrd = (b1->pos[0] - b2->pos[0])*(b1->pos[0] - b2->pos[0]) + (b1->pos[1] - b2->pos[1])*(b1->pos[1] - b2->pos[1]);
    return distance_sqrd <= ((b1->radius + b2->radius) * (b1->radius + b2->radius));
}

int checkWaveBodyCollision(body* wave, body* b1){
    double distance_sqrd = (b1->pos[0] - wave->pos[0])*(b1->pos[0] - wave->pos[0]) + (b1->pos[1] - wave->pos[1])*(b1->pos[1] - wave->pos[1]);
    return (distance_sqrd < ((b1->radius + wave->radius + GRAVITON_PULSE_PERIOD*C_speed)*(b1->radius + wave->radius + GRAVITON_PULSE_PERIOD*C_speed))) && distance_sqrd > ((b1->radius + wave->radius - GRAVITON_PULSE_PERIOD*C_speed) * (b1->radius + wave->radius - GRAVITON_PULSE_PERIOD*C_speed));
}

int checkWallCollision(body* b){
    int hit = 0;
    if ((b->pos[0] - b->radius <= 0 + BOX_WALL_THICKNESS/2 && b->velocity[0] < 0) ||
        (b->pos[0] + b->radius >= WIDTH_OF_BOX - BOX_WALL_THICKNESS/2 && b->velocity[0] > 0)) {
        change_body_velocity(b, -b->velocity[0], b->velocity[1]);
        hit = 1;
    }
    if ((b->pos[1] - b->radius <= 0 + BOX_WALL_THICKNESS/2 && b->velocity[1] < 0) ||
        (b->pos[1] + b->radius >= HEIGHT_OF_BOX - BOX_WALL_THICKNESS/2 && b->velocity[1] > 0)) {
        change_body_velocity(b, b->velocity[0], -b->velocity[1]);
        hit = 1;
    }

    //     int hit = 0;
    // if ((b->pos[0] - b->radius <= 0 + BOX_WALL_THICKNESS/2 && b->velocity[0] < 0) ||
    //     (b->pos[0] + b->radius >= WIDTH_OF_BOX - BOX_WALL_THICKNESS/2 && b->velocity[0] > 0)) {
    //     b->velocity[0] = -1* b->velocity[0];
    //     hit = 1;
    // }
    // if ((b->pos[1] - b->radius <= 0 + BOX_WALL_THICKNESS/2 && b->velocity[1] < 0) ||
    //     (b->pos[1] + b->radius >= HEIGHT_OF_BOX - BOX_WALL_THICKNESS/2 && b->velocity[1] > 0)) {
    //     b->velocity[1] = -1 * b->velocity[1];
    //     hit = 1;
    // }
    
    // if(hit) change_body_velocity(b, b->velocity[0], b->velocity[1]);
    
    // return hit;
    return hit;
}

int checkAndHandleWaveWallCollision(body* b){

    int source_wall = (b->pos[0] == (0 + BOX_WALL_THICKNESS/2)) || (b->pos[0] == (WIDTH_OF_BOX - BOX_WALL_THICKNESS/2)) || (b->pos[1] == (0 + BOX_WALL_THICKNESS/2)) || (b->pos[1] == (HEIGHT_OF_BOX - BOX_WALL_THICKNESS/2));
  
    int hit = 0;
    if (!source_wall && (((b->pos[0] - b->radius) <= (0 + BOX_WALL_THICKNESS/2)) || ((b->pos[0] + b->radius) >= (WIDTH_OF_BOX - BOX_WALL_THICKNESS/2)))){
        // b->pos[0] = ((b->pos[0] - b->radius) <= (0 + BOX_WALL_THICKNESS/2))? BOX_WALL_THICKNESS/2 : WIDTH_OF_BOX - BOX_WALL_THICKNESS/2;
        // b->radius = 0;
        hit = 1;
    }
    
    if (!source_wall && ((((b->pos[1] - b->radius) <= (0 + BOX_WALL_THICKNESS/2)) || ((b->pos[1] + b->radius) >= (HEIGHT_OF_BOX - BOX_WALL_THICKNESS/2))))){
        // b->pos[1] = ((b->pos[1] - b->radius) <= (0 + BOX_WALL_THICKNESS/2))? BOX_WALL_THICKNESS/2 : HEIGHT_OF_BOX - BOX_WALL_THICKNESS/2;
        // b->radius = 0;
        hit = 1;
    }
    

    return hit;
}

void handelCollision(GravSimModelManager* modelManager){
    for (int i = 0; i < modelManager->num_bodies; i++){
        if (!checkWallCollision(&modelManager->bodies[i])){
            for (int j = i + 1; j < modelManager->num_bodies; j++){
                if (checkBodyCollision(&modelManager->bodies[i], &modelManager->bodies[j])){
                    relativisticMassiveCollison(&modelManager->bodies[i], &modelManager->bodies[j]);
                }
            }
        }
    }

}

void makeNewBody(body* body, GravSimModelManager* modelManager){
    double radius = ((((double)rand())/RAND_MAX) * (_MAX_RADIUS - _MIN_RADIUS)) + modelManager->MIN_RADIUS;
    
    double x_pos = ((double)rand()/(double)RAND_MAX) * (double)(WIDTH_OF_BOX - radius) + (double)radius;
    double y_pos = ((double)rand()/(double)RAND_MAX) * (double)(HEIGHT_OF_BOX - radius) + (double)radius;
    double theta = ((double)rand()/(double)RAND_MAX) * 2.0f * (double)(M_PI);

    double v_total = (((double)rand()/(double)RAND_MAX) * (_MAX_SPEED- _MIN_SPEED) +  _MIN_SPEED) * C_speed;
    double v_x =  v_total * cos(theta);
    double v_y = v_total * sin(theta);

    uint64_t colour = wavelengthToHex((uint64_t)380 + (uint64_t)rand() % (700- 380));
    // uint32_t colour = wavelengthToHex(580);
    

    body_init(body, 100000000000.0 * M_PI*(radius*radius), x_pos, y_pos, v_x, v_y, radius, colour);
}

void setupBodies(GravSimModelManager* modelManager){
    for (int i = 0; i < modelManager->num_bodies; i++){
        makeNewBody(&modelManager->bodies[i], modelManager);

        for (int j = 0; j < i; j++){
            if (checkBodyCollision(&modelManager->bodies[i], &modelManager->bodies[j])){
                makeNewBody(&modelManager->bodies[i], modelManager);
                j = -1;
            }
        }

        relativisticDopplerEffect(&modelManager->bodies[i]);
    }

}

void newtonianGravity(body* bodies, int num_bodies){
    for (int i = 0; i < num_bodies; i++){

        double a_x = 0;
        double a_y = 0;
        for (int j = 0; j < num_bodies; j++){
            if (i != j){
                double masses_G = bodies[j].mass*GRAVITATIONAL_CONSTANT;
                double r_normal_factor = Q_rsqrt((bodies[j].pos[0] - bodies[i].pos[0])*(bodies[j].pos[0] - bodies[i].pos[0]) + (bodies[j].pos[1] - bodies[i].pos[1])*(bodies[j].pos[1] - bodies[i].pos[1]));
                a_x += r_normal_factor * masses_G *(bodies[j].pos[0] - bodies[i].pos[0])/((bodies[j].pos[0] - bodies[i].pos[0])*(bodies[j].pos[0] - bodies[i].pos[0]) + (bodies[j].pos[1] - bodies[i].pos[1])*(bodies[j].pos[1] - bodies[i].pos[1]));
                a_y += r_normal_factor * masses_G *(bodies[j].pos[1] - bodies[i].pos[1])/((bodies[j].pos[0] - bodies[i].pos[0])*(bodies[j].pos[0] - bodies[i].pos[0]) + (bodies[j].pos[1] - bodies[i].pos[1])*(bodies[j].pos[1] - bodies[i].pos[1]));

            }
        }
    
        
        bodies[i].velocity[0] = bodies[i].velocity[0] + EPSILON_STEP * a_x;
        bodies[i].velocity[1] = bodies[i].velocity[1] + EPSILON_STEP * a_y;
        bodies[i].speed = pow((bodies[i].velocity[0]*bodies[i].velocity[0]) + (bodies[i].velocity[1] * bodies[i].velocity[1]), 1.0/2.0);
    }


}

void shootGravitons(body* b_emmitor, body* b_reciver){
    gravitonList* gList = (gravitonList*)b_emmitor->gravitonList;

    body* graviton = malloc(sizeof(body));

    double one_over_dist = Q_rsqrt((b_emmitor->pos[0] - b_reciver->pos[0])*(b_emmitor->pos[0] - b_reciver->pos[0]) + (b_emmitor->pos[1] - b_reciver->pos[1])*(b_emmitor->pos[1] - b_reciver->pos[1]));

    double x_comp = (b_reciver->pos[0]  - b_emmitor->pos[0]) * one_over_dist;
    double y_comp =  (b_reciver->pos[1] - b_emmitor->pos[1]) * one_over_dist;
            
    particleGraviton_init(graviton, b_emmitor->pos[0] + (b_emmitor->radius * x_comp), b_emmitor->pos[1] + (b_emmitor->radius * y_comp), C_speed*x_comp, C_speed*y_comp, C_speed);
    
    graviton->momentum[0] = b_emmitor->pos[0];
    graviton->momentum[1] = b_emmitor->pos[1];

    push_graviton_to_list(gList, graviton, WIDTH_OF_BOX);
    
    double delta_width = C_speed;

    double angle_gamma = atan2((delta_width + b_emmitor->radius), 1/one_over_dist);

    double phi = atan2((b_reciver->pos[0]  - b_emmitor->pos[0]), (b_reciver->pos[1]  - b_emmitor->pos[1]));

    double del_theta = 0.001;
    for (double theta = -angle_gamma; theta < 0 ; theta += .05-del_theta){

        
        double x_comp = sin(phi - theta); 
        double y_comp = cos(phi - theta);

        graviton = malloc(sizeof(body));
            
        particleGraviton_init(graviton, b_emmitor->pos[0] + (b_emmitor->radius * x_comp), b_emmitor->pos[1] + (b_emmitor->radius * y_comp), C_speed*x_comp, C_speed*y_comp, C_speed);
            
        graviton->momentum[0] = b_emmitor->pos[0];
        graviton->momentum[1] = b_emmitor->pos[1];
        
        push_graviton_to_list(gList, graviton, WIDTH_OF_BOX);

        del_theta = del_theta * 3.0/4.0;
    }
    
    double old_del_theta = .05+del_theta;
    del_theta = 0.;
    for (double theta = 0; theta < angle_gamma ; theta += old_del_theta-del_theta){


        double x_comp = sin(phi - theta); 
        double y_comp = cos(phi - theta);       

        graviton = malloc(sizeof(body));
            
        particleGraviton_init(graviton, b_emmitor->pos[0] + (b_emmitor->radius * x_comp), b_emmitor->pos[1] + (b_emmitor->radius * y_comp), C_speed*x_comp, C_speed*y_comp, C_speed);
    
        graviton->momentum[0] = b_emmitor->pos[0];
        graviton->momentum[1] = b_emmitor->pos[1];
        
        push_graviton_to_list(gList, graviton, WIDTH_OF_BOX);

        del_theta = del_theta * 4.0/3.0;
    }

    // for (double theta = 0; theta < 2 * M_PI; theta+= .1){
   
    //     graviton->node = malloc(sizeof(body));

    //     double x_comp = sin(theta);
    //     double y_comp =  cos(theta);
        
    //     body_init(graviton->node, 0, b_emmitor->pos[0] + (b_emmitor->radius * x_comp), b_emmitor->pos[1] + (b_emmitor->radius * y_comp), C_speed*x_comp, C_speed*y_comp, 0, 0);
        
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
}

void emmitParticleGravitons(body* bodies, int num_bodies){
    for(int j = 0; j < num_bodies; j++){
        body* b_emmitor = &bodies[j];
    
        for(int i = 0; i < num_bodies; i++){
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
        
    //     body_init(graviton->node, 0, b_emmitor->pos[0] + (b_emmitor->radius * x_comp), b_emmitor->pos[1] + (b_emmitor->radius * y_comp), C_speed*x_comp, C_speed*y_comp, 0, 0);

    //     graviton->next = malloc(sizeof(gravitonNode));
    //     gravitonNode_init(graviton->next);
    //     graviton = graviton->next;
    // }

}

void handleParticleGravitonBodyCollision(body* b1, body* b_emmiter){

    gravitonNode* gravitons = ((gravitonList*)b_emmiter->gravitonList)->first;

    int collided = 0;
    while (gravitons->node != NULL){
        if (checkBodyCollision(b1, gravitons->node)){
            if (!collided){
                body* graviton_body = gravitons->node;
                
                double distance_sqrd = ((graviton_body->momentum[0]-b1->pos[0])*(graviton_body->momentum[0]-b1->pos[0])) + ((graviton_body->momentum[1] - b1->pos[1])*(graviton_body->momentum[1] - b1->pos[1]));
                
                double gravitational_momentum_coefficeint = GRAVITATIONAL_CONSTANT * b1->mass * b_emmiter->mass/ (distance_sqrd * C_speed);

                graviton_body->momentum[0] = gravitational_momentum_coefficeint * (graviton_body->momentum[0]-b1->pos[0]);
                graviton_body->momentum[1] = gravitational_momentum_coefficeint * (graviton_body->momentum[1]-b1->pos[1]);

                relativisticGravitonOnMassiveCollison(b1, graviton_body);
                collided = 1;
                
            }

            gravitons = remove_graviton_from_list(b_emmiter->gravitonList, gravitons);

        }else{
            gravitons = gravitons->next;
        }
    
    }
}

void handleGravitonBodyCollision(body* b1, body* b_emmiter){

    gravitonNode* gravitons = ((gravitonList*)b_emmiter->gravitonList)->first;

    int collided = 0;
    int waved = 0;

    while (gravitons->node != NULL){
        if (!waved && gravitons->node->speed == 0 && checkWaveBodyCollision(b1, gravitons->node)){
            
            body* graviton_body = gravitons->node;
            
            double distance_sqrd = ((graviton_body->pos[0]-b1->pos[0])*(graviton_body->pos[0]-b1->pos[0])) + ((graviton_body->pos[1] - b1->pos[1])*(graviton_body->pos[1] - b1->pos[1]));
            
            double gravitational_momentum_coefficeint = GRAVITATIONAL_CONSTANT * b1->mass * b_emmiter->mass/(distance_sqrd * C_speed);

            graviton_body->momentum[0] = gravitational_momentum_coefficeint * (graviton_body->pos[0]-b1->pos[0]);
            graviton_body->momentum[1] = gravitational_momentum_coefficeint * (graviton_body->pos[1]-b1->pos[1]);
            waved = 1;
            relativisticGravitonOnMassiveCollison(b1, graviton_body);

        }else if(gravitons->node->speed == C_speed && checkBodyCollision(b1, gravitons->node)){
                if (!collided){
                    body* graviton_body = gravitons->node;
                    
                    double distance_sqrd = ((graviton_body->momentum[0]-b1->pos[0])*(graviton_body->momentum[0]-b1->pos[0])) + ((graviton_body->momentum[1] - b1->pos[1])*(graviton_body->momentum[1] - b1->pos[1]));
                    
                    double gravitational_momentum_coefficeint = GRAVITATIONAL_CONSTANT * b1->mass * b_emmiter->mass/ (distance_sqrd * C_speed);

                    graviton_body->momentum[0] = gravitational_momentum_coefficeint * (graviton_body->momentum[0]-b1->pos[0]);
                    graviton_body->momentum[1] = gravitational_momentum_coefficeint * (graviton_body->momentum[1]-b1->pos[1]);

                    relativisticGravitonOnMassiveCollison(b1, graviton_body);
                    collided = 1;
                }

                gravitons = remove_graviton_from_list(b_emmiter->gravitonList, gravitons);
    
        }else{
            gravitons = gravitons->next;
        }
    
    }
}

void particleGravitonsWithRelevCollision(body* bodies, int num_bodies){
    for (int i = 0; i < num_bodies; i++){
        for(int j = 0; j < num_bodies; j++){
            if (i != j){
                handleParticleGravitonBodyCollision(&bodies[i], &bodies[j]);
            }
        }    
    }
}

void emmitWaveGravitons(body* bodies, int num_bodies){
    double corner[2];


    for(int i = 0; i < num_bodies; i++){
        body* b_emmitor = &bodies[i];

        body* graviton = malloc(sizeof(body));
        waveGraviton_init(graviton, bodies[i].pos[0], bodies[i].pos[1], bodies[i].radius); 
        
        if ((graviton->pos[0] - WIDTH_OF_BOX/2.0) < 0 ){
            corner[0] = WIDTH_OF_BOX;
        }else{
            corner[0] = 0;
        }
        
        if ((graviton->pos[1] - HEIGHT_OF_BOX/2.0) < 0){                
            corner[1] = HEIGHT_OF_BOX;
        }else{
            corner[1] = 0;
        }

        double max_distance = (corner[0] - graviton->pos[0])*(corner[0] - graviton->pos[0]) + (corner[1] - graviton->pos[1])*(corner[1] - graviton->pos[1]);       
        push_graviton_to_list(b_emmitor->gravitonList, graviton, max_distance);
        
    }

}

void handleWaveGravitonBodyCollision(body* b1, body* b_emmiter){

    gravitonNode* gravitons = ((gravitonList*)b_emmiter->gravitonList)->first;


    while (gravitons->node != NULL){
        if (checkWaveBodyCollision(b1, gravitons->node)){

            body* graviton_body = gravitons->node;
            
            double distance_sqrd = ((graviton_body->pos[0]-b1->pos[0])*(graviton_body->pos[0]-b1->pos[0])) + ((graviton_body->pos[1] - b1->pos[1])*(graviton_body->pos[1] - b1->pos[1]));
            
            double gravitational_momentum_coefficeint = GRAVITATIONAL_CONSTANT * b1->mass * b_emmiter->mass/ (distance_sqrd * C_speed);

            graviton_body->momentum[0] = gravitational_momentum_coefficeint * (graviton_body->pos[0]-b1->pos[0]);
            graviton_body->momentum[1] = gravitational_momentum_coefficeint * (graviton_body->pos[1]-b1->pos[1]);

            relativisticGravitonOnMassiveCollison(b1, graviton_body);

            break;
        }else{
            gravitons = gravitons->next;
        }
    
    }
}

void waveGravitonsWithRelevCollision(body* bodies, int num_bodies){
    for (int i = 0; i < num_bodies; i++){
        for(int j = 0; j < num_bodies; j++){
            if (i != j){
                handleWaveGravitonBodyCollision(&bodies[i], &bodies[j]);
            }
        }    
    }
}

void gravitionRelativeCollision(body* bodies, int num_bodies){
    for (int i = 0; i < num_bodies; i++){
        for(int j = 0; j < num_bodies; j++){
            if (i != j){
                handleGravitonBodyCollision(&bodies[i], &bodies[j]);
            }
        }    
    }
}

void relativisticEIHGravity(GravSimModelManager* modelManager){
    
    body* bodies = modelManager->bodies;

    double r_s[modelManager->num_bodies*modelManager->num_bodies];
    double one_r_s_sqrd[modelManager->num_bodies*modelManager->num_bodies];

    for (int i = 0; i<modelManager->num_bodies; i++){
        for (int j = i + 1; j< modelManager->num_bodies; j++){
            r_s[i] = pow((bodies[i].pos[1] - bodies[j].pos[1])*(bodies[i].pos[1] - bodies[j].pos[1]) + (bodies[i].pos[0] - bodies[j].pos[0])*(bodies[i].pos[0] - bodies[j].pos[0]), 1.0/2.0);
            r_s[i+ j * modelManager->num_bodies] = r_s[i];
            one_r_s_sqrd[i] = 1.0/r_s[i];
            one_r_s_sqrd[i+ j * modelManager->num_bodies] = 1.0/r_s[i];
        }
    }

    for (int i = 0; i<modelManager->num_bodies; i++){



    }
    
}

void updateBodies(GravSimModelManager* modelManager){
    
    body* bodies = modelManager->bodies;

    handelCollision(modelManager);

    for (int i = 0; i < modelManager->num_bodies; i++){
        bodies[i].pos[0] = bodies[i].pos[0] + EPSILON_STEP*bodies[i].velocity[0];
        bodies[i].pos[1] = bodies[i].pos[1] + EPSILON_STEP*bodies[i].velocity[1];
        relativisticDopplerEffect(&modelManager->bodies[i]);
    }
}

void updateParticleGravitons(body* bodies, int num_bodies, double C){

    C;
    for (int i = 0; i< num_bodies; i++){
        
        gravitonNode* gNode = ((gravitonList*)bodies[i].gravitonList)->first;

        while (gNode->node != NULL)
        {
            if (checkWallCollision(gNode->node)){
                gNode = remove_graviton_from_list((gravitonList*)bodies[i].gravitonList, gNode); // Remove this node and set it to the next node

            }else{

                gNode->node->pos[0] = gNode->node->pos[0] + EPSILON_STEP*gNode->node->velocity[0];
                gNode->node->pos[1] = gNode->node->pos[1] + EPSILON_STEP*gNode->node->velocity[1];
                gNode = gNode->next;
            }
        }

    }

}

void updateWaveGravitons(body* bodies, int num_bodies, double C){

    for (int i = 0; i< num_bodies; i++){
        
        gravitonNode* gNode = ((gravitonList*)bodies[i].gravitonList)->first;
        
        while (gNode->node != NULL)
        {
            if ((gNode->node->radius * gNode->node->radius) >= gNode->max_distance){ //TODO: Clean up max_distance Logic, currently it is equally ture max distance squared
                gNode = remove_graviton_from_list((gravitonList*)bodies[i].gravitonList, gNode);          

            // }else if (checkAndHandleWaveWallCollision(gNode->node) == 1){
            
            //     gNode = remove_graviton_from_list((gravitonList*)bodies[i].gravitonList, gNode);
            
            }else{
                gNode->node->radius = gNode->node->radius + EPSILON_STEP*C;
                gNode = gNode->next;
            }
        }

    }

}

void updategraviton(body* bodies, int num_bodies, double C){
    
    for (int i = 0; i< num_bodies; i++){
        
        gravitonNode* gNode = ((gravitonList*)bodies[i].gravitonList)->first;

        while (gNode->node != NULL)
        {
            if( gNode->node->speed == C_speed){
                if (checkWallCollision(gNode->node)){
                    gNode = remove_graviton_from_list((gravitonList*)bodies[i].gravitonList, gNode); // Remove this node and set it to the next node

                }else{

                    gNode->node->pos[0] = gNode->node->pos[0] + EPSILON_STEP*gNode->node->velocity[0];
                    gNode->node->pos[1] = gNode->node->pos[1] + EPSILON_STEP*gNode->node->velocity[1];
                    gNode = gNode->next;
                }
            }else{

                if ((gNode->node->radius * gNode->node->radius) >= gNode->max_distance){ //TODO: Clean up max_distance Logic, currently it is equally ture max distance squared
                    gNode = remove_graviton_from_list((gravitonList*)bodies[i].gravitonList, gNode);          

                    // }else if (checkAndHandleWaveWallCollision(gNode->node) == 1){
                    
                    //     gNode = remove_graviton_from_list((gravitonList*)bodies[i].gravitonList, gNode);
                
                }else{
                    gNode->node->radius = gNode->node->radius + EPSILON_STEP*C;
                    gNode = gNode->next;
                }
            }
        }

    }

}

body* getBodies(GravSimModelManager* modelManager){
    return modelManager->bodies;
}

void cleanGravitonList(gravitonList* gravitonList){
    gravitonNode* gNode = gravitonList->first;
    
    for (int i = 0; i<gravitonList->count; i++){
        free(gNode->node);
        free(gNode->prev);
        gNode = gNode->next;
    }

    free(gNode);
}

void cleanBodies(GravSimModelManager* modelManager){
    for (int i = 0; i<modelManager->num_bodies; i++){
        cleanGravitonList(modelManager->bodies[i].gravitonList);
        
    }
    free(modelManager->bodies);
}

void GravSim_change_grav_func(GravSimModelManager* modelManager, int grav_type){
    
    if (grav_type == 0){
        modelManager->gravity_update_function = gravitionRelativeCollision;
        modelManager->graviton_update_function = updategraviton;
        modelManager->graviton_emmision_function = emmitWaveGravitons;

    }else if (grav_type == 1){
        modelManager->gravity_update_function = gravitionRelativeCollision;
        modelManager->graviton_update_function = updategraviton;
        modelManager->graviton_emmision_function = emmitParticleGravitons;
    }else{
        modelManager->gravity_update_function = newtonianGravity;
        modelManager->graviton_update_function = updategraviton;
        modelManager->graviton_emmision_function = NULL;
    }
}

int GravSim_get_grav_type(GravSimModelManager* modelManager){
    if (modelManager->graviton_emmision_function == emmitWaveGravitons){
        return 0;
    }else if (modelManager->graviton_emmision_function == emmitParticleGravitons){
        return 1;
    }else{
        return 2;
    }
}

void GravSim_update_box_size(GravSimModelManager* modelManager, double new_width, double new_height){
    WIDTH_OF_BOX = new_width;
    HEIGHT_OF_BOX = new_height;
    modelManager->WIDTH= new_width;
    modelManager->HEIGHT = new_width;
}

GravSimModelManager* initalize_Sim(int num_bodies, double epsilon,  double C, double G, double width, double height, double thickness, int grav_type, double min_speed, double max_speed, double min_radius, double max_radius,void (*notify_observers_function)){
    GravSimModelManager* modelManager = malloc(sizeof(GravSimModelManager));

    GravSim_update_box_size(modelManager, width, height);

    BOX_WALL_THICKNESS = thickness;

    
    EPSILON_STEP = epsilon; // simmulation step size. 1 == a second
    modelManager->EPSILON = epsilon;

    C_speed = C;
    C_SQUARED = C*C;
    ONE_C_SQUARED = 1/C_SQUARED;

    _MIN_SPEED = min_speed; // The minimum speed of body i of C_speed. MUST BE << C_speed
    _MAX_SPEED = max_speed;// The maximum speed of body i of C_speed. MUST BE > _MIN_SPEED
    _MIN_RADIUS = min_radius;
    _MAX_RADIUS = max_radius;
    
    modelManager->notify_observers_function = notify_observers_function;
   
    modelManager->C = C;
    modelManager->G = G;

    GRAVITATIONAL_CONSTANT = G;

    GravSim_change_grav_func(modelManager, grav_type);
    modelManager->MIN_RADIUS = _MIN_RADIUS;


    modelManager->num_bodies = num_bodies;
    
    modelManager->bodies = malloc(modelManager->num_bodies * sizeof(body));
    setupBodies(modelManager);

    modelManager->GRAVITON_PULSE_PERIOD = _MIN_RADIUS/(C_speed);
    
    GRAVITON_PULSE_PERIOD = 1;//modelManager->_MIN_RADIUS/(C_speed);

    modelManager->time_since_last_pulse = 0;



    return modelManager;

}

void GravSim_model_take_step(GravSimModelManager* modelManager){
    
    modelManager->gravity_update_function(modelManager->bodies, modelManager->num_bodies);

    updateBodies(modelManager);  
     
    modelManager->graviton_update_function(modelManager->bodies, modelManager->num_bodies, C_speed);

    if (modelManager->graviton_emmision_function != NULL && modelManager->time_since_last_pulse >= modelManager->GRAVITON_PULSE_PERIOD){
        
        modelManager->graviton_emmision_function(modelManager->bodies, modelManager->num_bodies);

        modelManager->time_since_last_pulse = 0;
    
    }
    modelManager->time_since_last_pulse += modelManager->EPSILON;



}


void GravSim_clean(GravSimModelManager* MM){
    cleanBodies(MM);
    free(MM);
}
