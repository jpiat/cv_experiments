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
	EKF_TYPE * gyro = (EKF_TYPE*) private_data;

	EKF_TYPE dx[4][3];


}

static void readline(char * line, FILE * fp) {
	line = fgets(line, 1000, fp);
}

static void readdata(FILE * fp, double * timestamp, EKF_TYPE acc[3], EKF_TYPE gyro[3]) {
	char line[1000];

	readline(line, fp);

	char * p = strtok(line, ",");

	int i, j;

	(*timestamp) = acc[i] = atof(p);
	p = strtok(line, ",");

	for (i = 0; i < 3; ++i) {
		acc[i] = atof(p);
		p = strtok(NULL, ",");
	}

	for (j = 0; j < 3; ++j) {
		gyro[j] = atof(p);
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
int init_ekf_imu(struct float_ekf_state *state) {
#else
int init_ekf_imu(struct double_ekf_state *state) {
#endif
#ifdef FLOAT_EKF
	if(float_ekf_init(state, 7, 4) == 0) { //quaternion + angular velocities, only quaternion is observed through accelerometer
#else
	if (double_ekf_init(state, 7, 4) == 0) {
#endif
		printf("Cannot initialize EKF. Likely a memory allocation problem \n");
		return -1;
	}

// quaternion
	state->x[0] = 0.0;
	state->x[1] = 0.0;
	state->x[2] = 0.0;
	state->x[3] = 0.0;

// velocities
	state->x[4] = 0;
	state->x[5] = 0;
	state->x[6] = 0;

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
	init_ekf_imu(&ekf);

// Open input data file
	FILE * ifp = fopen("imu.csv", "r");

// Skip CSV header
	skipline(ifp);

// Make a place to store the data from the file and the output of the EKF
	EKF_TYPE Accelerometer_data[3];
	EKF_TYPE Gyroscope_data[3];

// Open output CSV file and write header
	const char * OUTFILE = "ekf.csv";
	FILE * ofp = fopen(OUTFILE, "w");
	fprintf(ofp, "X,Y,Z\n");

	int j, k;

// Loop till no more data
	for (j = 0; j < 25; ++j) {
		double time ;
		EKF_TYPE quaternion [4] ;
		readdata(ifp, &time, Accelerometer_data, Gyroscope_data);
		//model(&ekf, SV_Pos);
#ifdef FLOAT_EKF
		float_ekf_predict(&ekf, Gyroscope_data);
		//TODO: estimate quaternion from accelerometer data
		float_ekf_update(&ekf, quaternion);
#else
		double_ekf_predict(&ekf, Gyroscope_data);
		//TODO: estimate quaternion from accelerometer data
		double_ekf_update(&ekf, quaternion);
#endif


	}
// Done!
	fclose(ifp);
	fclose(ofp);
	printf("Wrote file %s\n", OUTFILE);
	return 0;

}

