#include "bodies.h"
#include "math.h"
#include <stdlib.h>
#include <stdio.h>

void gravitonNode_init(gravitonNode* node){
    node->node = NULL;
    node->next = NULL;
    node->prev = NULL;
    node->max_distance = -1.0;
}

void setupGravitonList(gravitonList* gravitonList){
    gravitonList->first = malloc(sizeof(gravitonNode));
    gravitonList->last = malloc(sizeof(gravitonNode));
    
    gravitonNode_init(gravitonList->first);
    gravitonNode_init(gravitonList->last);
    
    gravitonList->first->next = gravitonList->last;
    gravitonList->last->prev = gravitonList->first;
    gravitonList->count = 0;
}

void particleGraviton_init(body* body, double x, double y, double v_x, double v_y,double C){
    body->mass = 0;

    body->pos[0] = x;
    body->pos[1] = y;

    body->velocity[0] = v_x;
    body->velocity[1] = v_y;
    body->speed = C;


    body->radius = 0; 
    body->rest_colour = 0;
    body->relative_colour = 0;
    body->gravitonList = NULL;
}

void waveGraviton_init(body* body, double x, double y, double radius){
    body->mass = 0;

    body->pos[0] = x;
    body->pos[1] = y;

    body->velocity[0] = 0;
    body->velocity[1] = 0;


    body->speed = 0;


    body->radius = radius; 

    body->rest_colour = 0;
    body->relative_colour = 0;
    body->gravitonList = NULL;
}

void body_init(body* body, double mass, double x, double y, double v_x, double v_y, double radius, uint64_t colour){ 
    body->mass = mass;

    body->pos[0] = x;
    body->pos[1] = y;

    body->velocity[0] = v_x;
    body->velocity[1] = v_y;
    body->speed = pow((v_x * v_x) + (v_y*v_y), 1.0/2.0);

    body->radius = radius; 

    body->rest_colour = colour;
    body->relative_colour = colour;

    body->speed = pow((v_x*v_x) + (v_y*v_y), 1.0/2.0);
    body->gravitonList = malloc(sizeof(gravitonList));
    setupGravitonList(body->gravitonList);

    body->residual_step = 0.0;
}

void destroy_body(body* body){
    free(body->gravitonList);
    free(body);
}

int push_graviton_to_list(gravitonList* gravitonList, body* graviton, double max_distance){

    gravitonList->count += 1;
    if (gravitonList->first->node == NULL){

        gravitonList->first->node = graviton;
        gravitonList->first->max_distance = max_distance;   

    }else{
      
        gravitonList->last->node = graviton;
        gravitonList->last->max_distance = max_distance;

        gravitonList->last->next = malloc(sizeof(gravitonNode));
        
        if (gravitonList->last->next == NULL){
            printf("FAILED TO ALLOCATE LAST NODE FOR GList\n");
            return -1;
        }
        
        gravitonNode_init(gravitonList->last->next);

        ((gravitonNode*)gravitonList->last->next)->prev = gravitonList->last;
        
        gravitonList->last = gravitonList->last->next;

    }

    return 1;
}

gravitonNode* remove_graviton_from_list(gravitonList* gravitonList, gravitonNode* gNode){
    gravitonNode* retrun_node = NULL;
    destroy_body(gNode->node);

    if (gNode == gravitonList->first){

        if (((gravitonNode*)gravitonList->first->next)->node == NULL){
            gravitonList->first->node = NULL;
        }else{
            gravitonNode* old_node = gravitonList->first;
            gravitonList->first = gravitonList->first->next;
            free(old_node);
            gravitonList->first->prev = NULL;
        }

        retrun_node =  gravitonList->first;

    }else{
        
        retrun_node = gNode->next;
        ((gravitonNode *)gNode->prev)->next = gNode->next;
        ((gravitonNode *)gNode->next)->prev = gNode->prev;
        free(gNode);
    }
    gravitonList->count -= 1;
    return retrun_node;
}
