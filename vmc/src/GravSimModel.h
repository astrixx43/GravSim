#include "bodies.h"

typedef struct{
    double EPSILON;
    double WIDTH;
    double HEIGHT;

    void (*notify_observers_function)();

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

GravSimModelManager* initalize_Sim(int num_bodies, double epsilon,  double C, double G, double width, double height, double thickness, int grav_type, void (*notify_observers_function));
void GravSim_clean(GravSimModelManager* MM);
void GravSim_model_take_step(GravSimModelManager* modelManager);
void GravSim_change_grav_func(GravSimModelManager* modelManager, int grav_type);
double Q_rsqrt( double number );
