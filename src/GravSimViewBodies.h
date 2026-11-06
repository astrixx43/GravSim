#include <stdint.h>
#include <SDL3_ttf/SDL_ttf.h>
#include "GravSimModel.h"

typedef struct{
    double* pos; // Pos in x
    double* velocity;
    uint64_t* relative_colour;
    double* radius;
    void** gravitonList;
    TTF_Font* font;
}viewBody;

typedef struct{
    double pos[2]; // Pos in x
    uint64_t colour;
    double radius;

}traceBody;


typedef struct{
    int num_bodies;
    viewBody* bodies;

    int num_traces;
    traceBody** traces;

    int indexOfTrace;

    TTF_Font* menuFont;

}viewBodiesManager;


viewBodiesManager* viewBodiesManage_init(GravSimModelManager* MM, int num_traces, int SHOWTEXT);
void GravSim_change_grav_func(GravSimModelManager* modelManager, int grav_type);
int GravSim_get_grav_type(GravSimModelManager* modelManager);
void GravSimViewBodies_updateTraces(viewBodiesManager* VBM);
void GravSimViewBodies_end_bodyView(viewBodiesManager* VBM);
