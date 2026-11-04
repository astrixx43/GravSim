#include <stdint.h>
typedef struct{
    double mass; // Mass in kg
    double pos[2]; // Pos in x
    double velocity[2]; // in m/s where [0] is the x-compnet [1] is y
    double speed;
    uint64_t rest_colour;
    uint64_t relative_colour;
    int radius;
    double momentum[2];
    void* gravitonList;

}body;


typedef struct{
    double pos[2]; // Pos in x
    uint64_t colour;
    int radius;

}traceBody;

typedef struct gravitonList{
    body* head;
    void* next;
}gravitonList;

void setupGravitonList(gravitonList* gravitonList);

void body_init(body* body, double mass, double x, double y, double v_x, double v_y, int radius, uint64_t colour);

void trace_init(traceBody* trace, body* body);