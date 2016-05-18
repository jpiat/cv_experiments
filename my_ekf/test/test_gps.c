#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "float_ekf.h"
#include "float_mat_ops.h"

#include "double_ekf.h"
#include "double_mat_ops.h"

#include "mat_types.h"


#ifdef FLOAT_EKF
#define EKF_TYPE float
#else
#define EKF_TYPE double
#endif


const EKF_TYPE T = 1;


#ifdef FLOAT_EKF
static void blkfill(struct float_ekf_state * ekf, const EKF_TYPE * a, int off) {
#else
static void blkfill(struct double_ekf_state * ekf, const EKF_TYPE * a, int off) {
#endif
	off *= 2;
	MAT_PTR_ELT_AT((ekf->Q), off, off) = a[0];
	MAT_PTR_ELT_AT((ekf->Q), off, (off + 1)) = a[1];
	MAT_PTR_ELT_AT((ekf->Q), (off + 1), off) = a[2];
	MAT_PTR_ELT_AT((ekf->Q), (off + 1), (off + 1)) = a[3];
}

#ifdef FLOAT_EKF
void double_ekf_model(struct double_ekf_state * ekf, void * private_data) {

}
#else
void float_ekf_model(struct float_ekf_state * ekf, void * private_data) {

}
#endif

//EKF_TYPE SV[4][3]
#ifdef FLOAT_EKF
void float_ekf_model(struct float_ekf_state * ekf, void * private_data) {
#else
void double_ekf_model(struct double_ekf_state * ekf, void * private_data) {
#endif
	int i, j;
	EKF_TYPE * SV = (EKF_TYPE*) private_data;

	//Update of state based on previous state and Time step
	for (j = 0; j < 8; j += 2) {
		ekf->fx[j] = ekf->x[j] + T * ekf->x[j + 1];
		ekf->fx[j + 1] = ekf->x[j + 1];
	}

	//Computation of F based on knowledge of process model
	//Jacobian process F is partial derivative of state evolution (in this case T is the derivative)
	for (j = 0; j < 8; ++j)
		MAT_PTR_ELT_AT(ekf->F,j, j) = 1;

	for (j = 0; j < 4; ++j)
		MAT_PTR_ELT_AT(ekf->F, (2 * j), (2 * j + 1)) = T;

	EKF_TYPE dx[4][3];

	for (i = 0; i < 4; ++i) {
		ekf->hx[i] = 0;
		for (j = 0; j < 3; ++j) {
			EKF_TYPE d = ekf->fx[j * 2] - SV[i*3 + j];
			dx[i][j] = d;
			ekf->hx[i] += d * d;
		}
		ekf->hx[i] = pow(ekf->hx[i], 0.5) + ekf->fx[6];
	}

	for (i = 0; i < 4; ++i) {
		for (j = 0; j < 3; ++j)
			MAT_PTR_ELT_AT(ekf->H, i, j * 2) = dx[i][j] / ekf->hx[i];
		MAT_PTR_ELT_AT(ekf->H, i, 6) = 1;
	}
}

static void readline(char * line, FILE * fp) {
	line = fgets(line, 1000, fp);
}

static void readdata(FILE * fp, EKF_TYPE SV_Pos[4][3], EKF_TYPE SV_Rho[4]) {
	char line[1000];

	readline(line, fp);

	char * p = strtok(line, ",");

	int i, j;

	for (i = 0; i < 4; ++i)
		for (j = 0; j < 3; ++j) {
			SV_Pos[i][j] = atof(p);
			p = strtok(NULL, ",");
		}

	for (j = 0; j < 4; ++j) {
		SV_Rho[j] = atof(p);
		p = strtok(NULL, ",");
	}
}

static void skipline(FILE * fp) {
	char line[1000];
	readline(line, fp);
}

void error(const char * msg) {
	fprintf(stderr, "%s\n", msg);
}
#ifdef FLOAT_EKF
int init_ekf_gps(struct float_ekf_state *state) {
#else
int init_ekf_gps(struct double_ekf_state *state) {
#endif
#ifdef FLOAT_EKF
	if(float_ekf_init(state, 8, 4) == 0) {
#else
	if (double_ekf_init(state, 8, 4) == 0) {
#endif
		printf("Cannot initialize EKF. Likely a memory allocation problem \n");
		return -1;
	}
	const EKF_TYPE Sf = 36;
	const EKF_TYPE Sg = 0.01;
	const EKF_TYPE sigma = 5;         // state transition variance
	const EKF_TYPE Qb[4] = { Sf * T + Sg * T * T * T / 3, Sg * T * T / 2, Sg * T
			* T / 2, Sg * T };
	const EKF_TYPE Qxyz[4] = { sigma * sigma * T * T * T / 3, sigma * sigma * T
			* T / 2, sigma * sigma * T * T / 2, sigma * sigma * T };

	blkfill(state, Qxyz, 0);
	blkfill(state, Qxyz, 1);
	blkfill(state, Qxyz, 2);
	blkfill(state, Qb, 3);
	EKF_TYPE P0 = 10;
	EKF_TYPE R0 = 36;

	int i;

	for (i = 0; i < 8; ++i){
		MAT_PTR_ELT_AT(state->P, i,i) = P0;
		MAT_PTR_ELT_AT(state->Pp, i,i) = P0;
		//Pp and P must contain same initial values since prediction will be based on Pp
	}
	for (i = 0; i < 4; ++i)
		MAT_PTR_ELT_AT(state->R, i,i) = R0;

	state->x[0] = -2.168816181271560e+006;
	state->x[2] = 4.386648549091666e+006;
	state->x[4] = 4.077161596428751e+006;

// velocity
	state->x[1] = 0;
	state->x[3] = 0;
	state->x[5] = 0;

// clock bias
	state->x[6] = 3.575261153706439e+006;

// clock drift
	state->x[7] = 4.549246345845814e+001;

	return 0;
}

int main(void) {
// Do generic EKF initialization
#ifdef FLOAT_EKF
	struct float_ekf_state ekf;
#else
	struct double_ekf_state ekf;
#endif
// Do local initialization
	init_ekf_gps(&ekf);

// Open input data file
	FILE * ifp = fopen("gps.csv", "r");

// Skip CSV header
	skipline(ifp);

// Make a place to store the data from the file and the output of the EKF
	EKF_TYPE SV_Pos[4][3];
	EKF_TYPE SV_Rho[4];
	EKF_TYPE Pos_KF[25][3];

// Open output CSV file and write header
	const char * OUTFILE = "ekf.csv";
	FILE * ofp = fopen(OUTFILE, "w");
	fprintf(ofp, "X,Y,Z\n");

	int j, k;

// Loop till no more data
	for (j = 0; j < 25; ++j) {
		readdata(ifp, SV_Pos, SV_Rho);
		//model(&ekf, SV_Pos);
#ifdef FLOAT_EKF
		float_ekf_predict(&ekf, SV_Pos);
		float_ekf_update(&ekf, SV_Rho);
#else
		double_ekf_predict(&ekf, SV_Pos);
		if(double_ekf_update(&ekf, SV_Rho) < 0){
			printf("EKF error \n");
			exit(-1);
		}
#endif

		// grab positions, ignoring velocities
		for (k = 0; k < 3; ++k)
			Pos_KF[j][k] = ekf.x[2 * k];
	}

// Compute means of filtered positions
	EKF_TYPE mean_Pos_KF[3] = { 0, 0, 0 };
	for (j = 0; j < 25; ++j)
		for (k = 0; k < 3; ++k)
			mean_Pos_KF[k] += Pos_KF[j][k];
	for (k = 0; k < 3; ++k)
		mean_Pos_KF[k] /= 25.0;

// Dump filtered positions minus their means
	for (j = 0; j < 25; ++j) {
		fprintf(ofp, "%lf,%lf,%lf\n", Pos_KF[j][0] - mean_Pos_KF[0],
				Pos_KF[j][1] - mean_Pos_KF[1], Pos_KF[j][2] - mean_Pos_KF[2]);
		printf("%lf %lf %lf\n", Pos_KF[j][0], Pos_KF[j][1], Pos_KF[j][2]);
	}

// Done!
	fclose(ifp);
	fclose(ofp);
	printf("Wrote file %s\n", OUTFILE);
	return 0;

}

