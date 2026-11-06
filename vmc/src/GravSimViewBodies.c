#include <stdint.h>
#include "GravSimViewBodies.h"
#include <stdlib.h>

void viewBody_init(viewBody* viewBody, body* modelBody){
    viewBody->pos = modelBody->pos; // Pos in x
    viewBody->velocity = modelBody->velocity;

    viewBody->relative_colour = &modelBody->relative_colour;
    viewBody->radius = &modelBody->radius;
    viewBody->gravitonList = &modelBody->gravitonList;
    
    viewBody->font = NULL;
}

void setupViewBodies(viewBodiesManager* VBM, GravSimModelManager* MM){
    VBM->bodies = malloc(VBM->num_bodies * sizeof(viewBody));

    for (int i = 0; i < MM->num_bodies; i++){
        viewBody_init(&VBM->bodies[i], &MM->bodies[i]);
    }
}

void trace_init(traceBody* trace, double positon[2], double radius, uint64_t colourRGBA){ 

    trace->pos[0] = positon[0];
    trace->pos[1] = positon[1];

    trace->radius = radius; 
    
    uint32_t traceAlpha = colourRGBA & 0XFF;

    trace->colour = (colourRGBA & 0XFFFFFF00) + traceAlpha;

}

void setupBodyTraces(viewBodiesManager* VBM){
    if (VBM->num_traces == 0) return;
    VBM->traces = malloc(VBM->num_bodies*sizeof(traceBody*));

    for (int i = 0; i < VBM->num_bodies; i++){
        VBM->traces[i] = malloc(VBM->num_traces*sizeof(traceBody));

    }
}

void GravSimViewBodies_updateTraces(viewBodiesManager* VBM){
    if (VBM->num_traces == 0) return;
    for (int i = 0; i < VBM->num_bodies; i++){
        trace_init(&VBM->traces[i][VBM->indexOfTrace], VBM->bodies[i].pos, *VBM->bodies[i].radius, *VBM->bodies[i].relative_colour);
    }

    VBM->indexOfTrace = (VBM->indexOfTrace + 1) % VBM->num_traces;
}

viewBodiesManager* viewBodiesManage_init(GravSimModelManager* MM, int num_traces, int SHOWTEXT){
    viewBodiesManager* VBM = malloc(sizeof(viewBodiesManager));

    VBM->num_bodies = MM->num_bodies;
    setupViewBodies(VBM, MM);

    VBM->num_traces = num_traces;
    setupBodyTraces(VBM);

    VBM->indexOfTrace = 0;

    if (SHOWTEXT){
        for (int i = 0; i< VBM->num_bodies; i++){
            VBM->bodies[i].font = TTF_OpenFont("./src/fonts/dejavu-sans-condensed-webfont.ttf",  4 * (*(VBM->bodies[i]).radius)/15.0);
            TTF_SetFontStyle(VBM->bodies[i].font, TTF_STYLE_BOLD);
        }
        VBM->menuFont = TTF_OpenFont("./src/fonts/dejavu-sans-condensed-webfont.ttf", 20);
        TTF_SetFontStyle( VBM->menuFont, TTF_STYLE_BOLD);
    }

    return VBM;

}

void closeViewBodies(viewBodiesManager* VBM){
    if (VBM->num_bodies == 0) return;
    for (int i = 0; i<VBM->num_bodies; i++){
        TTF_CloseFont(VBM->bodies[i].font);
        
    }
    free(VBM->bodies);

}

void closeTraces(viewBodiesManager* VBM){
    if (VBM->num_traces == 0) return;
    for (int i = 0; i < VBM->num_bodies; i++){
        free(VBM->traces[i]);
    }
    free(VBM->traces);

}

void closeMenueFonts(viewBodiesManager* VBM){
    TTF_CloseFont(VBM->menuFont);
}

void GravSimViewBodies_end_bodyView(viewBodiesManager* VBM){
    closeViewBodies(VBM);
    closeTraces(VBM);
    closeMenueFonts(VBM);
    free(VBM);
}