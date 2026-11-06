
#include <stdlib.h>
#include <stdio.h>
#include <xmmintrin.h>
#include <string.h>
#include <math.h>

#include "GravSimModel.h"
#include "Wavelength_Rgb.h"


void relativisticGravitonImpulseOnMassiveCollison(body* b1, body* masslessBody){

    double gamma_1 = calculateGamma(b1->speed);

    double p_m_1_x_i = gamma_1*b1->velocity[0];
    double p_m_1_y_i = gamma_1*b1->velocity[1];

    double impulse_m_x = GRAVITON_PULSE_PERIOD * masslessBody->momentum[0]/b1->mass;
    double impulse_m_y = GRAVITON_PULSE_PERIOD * masslessBody->momentum[1]/b1->mass;
    
    double p_1_x_f = impulse_m_x - p_m_1_x_i;
    double p_1_y_f = impulse_m_y - p_m_1_y_i;

    double devider = Q_rsqrt((C_speed*C_speed) + ((p_1_x_f*p_1_x_f)+(p_1_y_f*p_1_y_f)));
    
    double v_x_f;
    double v_y_f;

    if (p_1_x_f > 0) {
        v_x_f = C_speed*p_1_x_f*devider;
    }else{
        v_x_f = -1 *C_speed*p_1_x_f*devider;
    }

    if (p_1_y_f > 0) {
        v_y_f = C_speed*p_1_y_f*devider;
    }else{
        v_y_f = -1 * C_speed*p_1_y_f*devider;
    }
    
    change_body_velocity(b1, v_x_f, v_y_f);

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
    

        // double step_size = M_PI/10;
        // body* graviton;

        // for (double theta = 0; theta < 2 * M_PI; theta +=step_size){
        //     graviton = malloc(sizeof(body));

        //     double x_comp = sin(theta);
        //     double y_comp = cos(theta);
            
        //     body_init(graviton, 0, b_emmitor->pos[0] + (b_emmitor->radius * x_comp), b_emmitor->pos[1] + (b_emmitor->radius * y_comp), C_speed*x_comp, C_speed*y_comp, 0, 0);

        //     push_graviton_to_list(b_emmitor->gravitonList, graviton, -1);
        // }
    }

}

void handleParticleGravitonBodyCollision(body* b1, body* b_emmiter){

    gravitonNode* gravitons = ((gravitonList*)b_emmiter->gravitonList)->first;

    int collided = 0;
    while (gravitons->node != NULL){
        if (checkBodyCollision(b1, gravitons->node)){
            if (!collided){
                body* graviton_body = gravitons->node;
                
                double distance_sqrd = ((graviton_body->momentum[0]-b1->pos[0])*(graviton_body->momentum[0]-b1->pos[0])) + ((graviton_body->momentum[1] - b1->pos[1])*(graviton_body->momentum[1] - b1->pos[1]));

                double one_distance = fast_inverse_sqrt(distance_sqrd);

                double gravitational_momentum_coefficeint = GRAVITATIONAL_CONSTANT * b1->mass * b_emmiter->mass * one_distance/ (distance_sqrd);

                graviton_body->momentum[0] = gravitational_momentum_coefficeint * (graviton_body->momentum[0]-b1->pos[0]);
                graviton_body->momentum[1] = gravitational_momentum_coefficeint * (graviton_body->momentum[1]-b1->pos[1]);

                relativisticGravitonOnMassiveCollison(b1, graviton_body);
                collided = 1;
                
            }

            gravitons = remove_graviton_from_list(b_emmiter->gravitonList, gravitons);
            // gravitons = gravitons->next;
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
        if (gravitons->node->speed == 0 && checkWaveBodyCollision(b1, gravitons->node)){

            body* graviton_body = gravitons->node;
            
            double distance_sqrd = ((graviton_body->pos[0]-b1->pos[0])*(graviton_body->pos[0]-b1->pos[0])) + ((graviton_body->pos[1] - b1->pos[1])*(graviton_body->pos[1] - b1->pos[1]));

            double one_distance = fast_inverse_sqrt(distance_sqrd);

            double gravitational_momentum_coefficeint = GRAVITATIONAL_CONSTANT * b1->mass * b_emmiter->mass * one_distance/ (distance_sqrd * C_speed);

            graviton_body->momentum[0] = gravitational_momentum_coefficeint * (graviton_body->pos[0]-b1->pos[0]);
            graviton_body->momentum[1] = gravitational_momentum_coefficeint * (graviton_body->pos[1]-b1->pos[1]);

            relativisticGravitonOnMassiveCollison(b1, graviton_body);

            break;
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
            double one_distance = 1; //Q_rsqrt(distance_sqrd);

            double gravitational_momentum_coefficeint = GRAVITATIONAL_CONSTANT * b1->mass * b_emmiter->mass * one_distance/ (distance_sqrd);

            graviton_body->momentum[0] = gravitational_momentum_coefficeint * (graviton_body->pos[0]-b1->pos[0]);
            graviton_body->momentum[1] = gravitational_momentum_coefficeint * (graviton_body->pos[1]-b1->pos[1]);

            relativisticGravitonOnMassiveCollison(b1, graviton_body);
            waved = 1;

        }else if(gravitons->node->speed == C_speed && checkBodyCollision(b1, gravitons->node)){
            if (!collided){
                body* graviton_body = gravitons->node;
                
                double distance_sqrd = ((graviton_body->momentum[0]-b1->pos[0])*(graviton_body->momentum[0]-b1->pos[0])) + ((graviton_body->momentum[1] - b1->pos[1])*(graviton_body->momentum[1] - b1->pos[1]));
                double one_distance = 1; //Q_rsqrt(distance_sqrd);
            
                double gravitational_momentum_coefficeint = GRAVITATIONAL_CONSTANT * b1->mass * b_emmiter->mass * one_distance/ (distance_sqrd);

                graviton_body->momentum[0] = gravitational_momentum_coefficeint * (graviton_body->momentum[0]-b1->pos[0]);
                graviton_body->momentum[1] = gravitational_momentum_coefficeint * (graviton_body->momentum[1]-b1->pos[1]);

                relativisticGravitonOnMassiveCollison(b1, graviton_body);
                collided = 1;
                
            }

            gravitons = remove_graviton_from_list(b_emmiter->gravitonList, gravitons);
            // gravitons = gravitons->next;
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

void relativisticUpdateBodies(GravSimModelManager* modelManager){
    
    body* bodies = modelManager->bodies;

    for (int i = 0; i < modelManager->num_bodies; i++){
        double time_prime = calculateGamma(bodies[i].speed);
       
        double step = 0;
        double graviton_pulse_period = C_speed * time_prime/MIN_RADIUS;
        double delta_dis = C_speed * EPSILON_STEP;
        double distance_since_last_pulse = 0;
        while(step <= (time_prime + modelManager->bodies[i].residual_step)){

            if (!checkWallCollision(&modelManager->bodies[i])){

                for (int j = i + 1; j < modelManager->num_bodies; j++){

                    if (checkBodyCollision(&modelManager->bodies[i], &modelManager->bodies[j])){
                        relativisticMassiveCollison(&modelManager->bodies[i], &modelManager->bodies[j]);
                    }

                    handleGravitonBodyCollision(&modelManager->bodies[i], &modelManager->bodies[j]);

                }
            }
            
            bodies[i].pos[0] = bodies[i].pos[0] + EPSILON_STEP*bodies[i].velocity[0];
            bodies[i].pos[1] = bodies[i].pos[1] + EPSILON_STEP*bodies[i].velocity[1]; 

            step += 1;
            
            
            if (distance_since_last_pulse >= graviton_pulse_period){
                
                modelManager->graviton_emmision_function(modelManager->bodies, modelManager->num_bodies);

                distance_since_last_pulse = 0;
            
            }
            
            distance_since_last_pulse += delta_dis;

        }

        modelManager->bodies[i].residual_step += (time_prime - step);
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

body* getBodies(GravSimModelManager* modelManager){
    return modelManager->bodies;
}

void GravSim_model_take_step(GravSimModelManager* modelManager){
    
    modelManager->gravity_update_function(modelManager->bodies, modelManager->num_bodies);

    relativisticUpdateBodies(modelManager);  
     
    modelManager->graviton_update_function(modelManager->bodies, modelManager->num_bodies, C_speed);
    

    if (modelManager->graviton_emmision_function != NULL && modelManager->time_since_last_pulse >= modelManager->GRAVITON_PULSE_PERIOD){
        
        modelManager->graviton_emmision_function(modelManager->bodies, modelManager->num_bodies);

        modelManager->time_since_last_pulse = 0;
    
    }
    modelManager->time_since_last_pulse += modelManager->EPSILON;



}


