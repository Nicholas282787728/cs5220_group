#ifndef LEAPFROG_H
#define LEAPFROG_H

#include "state.h"

void leapfrog_start(sim_state_t* s,proc_info* pInfo, double dt);
void leapfrog_step(sim_state_t* s, proc_info* pInfo,double dt);

#endif /* LEAPFROG_H */
