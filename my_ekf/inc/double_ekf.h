#include <stdlib.h>

#include "mat_types.h"

#ifndef DOUBLE_EKF_H
#define DOUBLE_EKF_H


struct double_ekf_state {
	unsigned int state_size;
	unsigned int state_observable;
	double * x; //state vector
	struct double_matrix * P; //prediction error covariance
	struct double_matrix * Q; //process noise covariance
	struct double_matrix * R; //measurement error covariance
	struct double_matrix * K; //kalman gains
	struct double_matrix * F; //Jacobian of process model, state transition matrix
	struct double_matrix * H; //Jacobian of measurement model

	struct double_matrix * Pp; //updated P


	struct double_matrix * PpHt; //n lines, m column
	struct double_matrix * HPp; //m lines, n columns
	struct double_matrix * HPpHt; // m lines, m colums
	struct double_matrix * HPpHt_inv; // m lines, m colums
	struct double_matrix * FP; //n lines, ncolumns

	struct double_matrix * KH; // n lines, m columns, allocated in same space as PpHt
	struct double_matrix * ImKHP; //n lines, m columns, allocated in same space as HPp


	double * hx ;
	double * fx ;
};

int double_ekf_init(struct double_ekf_state * e, unsigned int state_size,
		unsigned int observed_size);

int double_ekf_predict(struct double_ekf_state * e, void * private_data);

int double_ekf_update(struct double_ekf_state * e, double * z);

void double_ekf_close(struct double_ekf_state * e);

extern void double_ekf_model(struct double_ekf_state * ekf, void * private_data) ;

#endif
