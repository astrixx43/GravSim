#include "bodies.h"

typedef struct{
    double EPSILON;
    double WIDTH;
    double HEIGHT;

    int num_bodies;
    body* bodies;
    
    double C;
    double G;
    double MIN_SPEED;
    double MAX_SPEED;
    double MIN_RADIUS;

    double GRAVITON_PULSE_PERIOD;
    double time_since_last_pulse;
    
    void (*gravity_update_function)(body*, int);
    void (*graviton_update_function)(body*, int, double);
    void (*graviton_emmision_function)(body*, int);

}GravSimModelManager;

GravSimModelManager* initalize_Sim(int num_bodies, double epsilon,  double C, double G, double width, double height, double thickness, double min_radius, int grav_type);
void GraveSim_model_Sim(GravSimModelManager* MM);
void GraveSim_model_take_step(GravSimModelManager* modelManager);
double Q_rsqrt( double number );
