#include "bodies.h"
#include "math.h"


void body_init(body* body, double mass, double x, double y, double v_x, double v_y, int radius, uint64_t colour){ 
    body->mass = mass;

    body->pos[0] = x;
    body->pos[1] = y;

    body->velocity[0] = v_x;
    body->velocity[1] = v_y;
    body->speed = pow((v_x*v_x) + (v_y*v_y), 1.0/2.0);

    body->radius = radius; 

    body->rest_colour = colour;
    body->relative_colour = colour;
}

void trace_init(traceBody* trace, body* body){ 

    trace->pos[0] = body->pos[0];
    trace->pos[1] = body->pos[1];

    trace->radius = body->radius; 

    
    uint32_t traceAlpha = body->relative_colour & 0XFF;

    trace->colour = (body->relative_colour & 0XFFFFFF00) + traceAlpha;

}