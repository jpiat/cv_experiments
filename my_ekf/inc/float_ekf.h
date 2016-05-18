#include <stdlib.h>

#include "mat_types.h"

#ifndef FLOAT_EKF_H
#define FLOAT_EKF_H

struct float_ekf_state {
	unsigned int state_size;
	unsigned int state_observable;
	float * x; //state vector
	struct float_matrix * P; //prediction error covariance
	struct float_matrix * Q; //process noise covariance
	struct float_matrix * R; //measurement error covariance
	struct float_matrix * K; //kalman gains
	struct float_matrix * F; //Jacobian of process model, state transition matrix
	struct float_matrix * H; //Jacobian of measurement model

	struct float_matrix * Pp; //updated P


	struct float_matrix * PpHt; //n lines, m column
	struct float_matrix * HPp; //m lines, n columns
	struct float_matrix * HPpHt; // m lines, m colums
	struct float_matrix * HPpHt_inv; // m lines, m colums
	struct float_matrix * FP; //n lines, ncolumns

	struct float_matrix * KH; // n lines, m columns, allocated in same space as PpHt
	struct float_matrix * ImKHP; //n lines, m columns, allocated in same space as HPp


	float * hx ;
	float * fx ;
};

int float_ekf_init(struct float_ekf_state * e, unsigned int state_size,
		unsigned int observed_size);

int float_ekf_predict(struct float_ekf_state * e, void * private_data);

int float_ekf_update(struct float_ekf_state * e, float * z);

void float_ekf_close(struct float_ekf_state * e);

extern void float_ekf_model(struct float_ekf_state * ekf, void * private_data);


#endif
