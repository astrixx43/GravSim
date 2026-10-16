#include <stdint.h>
typedef struct{
    double mass; // Mass in kg
    double pos[2]; // Pos in x
    double velocity[2]; // in m/s where [0] is the x-compnet [1] is y
    double speed;
    uint64_t rest_colour;
    uint64_t relative_colour;
    int radius;

}body;

void body_init(body* body, int mass, double x, double y, double v_x, double v_y, int radius, uint64_t colour);