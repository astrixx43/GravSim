#include <stdint.h>

typedef struct{
    double mass; // Mass in kg
    double pos[2]; // Pos in x
    double velocity[2]; // in m/s where [0] is the x-compnet [1] is y
    double speed;
    uint64_t rest_colour;
    uint64_t relative_colour;
    double radius;
    double momentum[2];
    void* gravitonList;

}body;

typedef struct gravitonNode{
    body* node;
    void* next;
    void* prev;
    double max_distance;
} gravitonNode;

typedef struct gravitonList{
    gravitonNode* first;
    gravitonNode* last;
    uint64_t count;
}gravitonList;

void body_init(body* body, double mass, double x, double y, double v_x, double v_y, double radius, uint64_t colour);

int push_graviton_to_list(gravitonList* gravitonList, body* graviton, double max_distance);

gravitonNode* remove_graviton_from_list(gravitonList* gravitonList, gravitonNode* gNode);