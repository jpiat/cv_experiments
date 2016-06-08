#include "dlt_plane.hpp"
#include <cmath>

// #define NORMALIZE

int normalize_for_H(Matrix<double, 3, 1> * obs1, Matrix<double, 3, 1> * obs1n,
		Matrix<double, 3, 3> & H, unsigned int nb_pts) {
	unsigned int i;
	double scale, tx, ty;
	Matrix<double, 3, 1> mean;
	double dist = 0;
	mean(0, 0) = 0.;
	mean(1, 0) = 0.;

	for (i = 0; i < nb_pts; i++) {
		mean += obs1[i];
	}
	mean /= nb_pts;
	mean(2, 0) = 0.;

	for (i = 0; i < nb_pts; i++) {
		obs1n[i] = obs1[i] - mean;
		dist += obs1n[i].block<2, 1>(0, 0).norm();
	}
	double meandist = dist / nb_pts;
	scale = sqrt(2.) / meandist;

	for (i = 0; i < nb_pts; i++) {
		obs1n[i] = obs1n[i] * (scale);
		obs1n[i](2, 0) = obs1[i](2, 0); //force homogeneous
		//obs1n[i](2, 0) = 1.;
	}

	tx = -(scale) * mean(0);
	ty = -(scale) * mean(1);

	H << scale, 0, tx, 0, scale, ty, 0, 0, 1;

	return 1;
}

int denormalize_H(Matrix<double, 3, 3> H, Matrix<double, 3, 3> H_norm1,
		Matrix<double, 3, 3> H_norm2, Matrix<double, 3, 3> & H_denorm) {
	H_denorm = H_norm1.inverse() * H * H_norm2;
	return 1;
}

/**
 * Estimate homography based on 4-points observations
 * Observations must be provided in homogeneous coordinates
 * This DLT assumes Hobs = H * obs
 */
int estimate_homo_dlt(Matrix<double, 3, 1> * obs, Matrix<double, 3, 1> * Hobs,
		Matrix<double, 3, 3> & H) {
	int i;

	Matrix<double, 8, 8> A;
	Matrix<double, 8, 8> V;
	Matrix<double, 8, 1> b;
	Matrix<double, 8, 1> X;

	//if not doing the following, the homography is defined up to a scale factor ?
	for (i = 0; i < 4; i++) {
		Hobs[i] = Hobs[i]/Hobs[i](2,0);
	}
#ifdef NORMALIZE
	Matrix<double, 3, 3> H_norm1, H_norm2;
	Matrix<double, 3, 1> obs1n[4];
	Matrix<double, 3, 1> obs2n[4];
	normalize_for_H(Hobs, obs1n, H_norm1, 4);
	normalize_for_H(obs, obs2n, H_norm2, 4);
	Hobs = &obs1n[0];
	obs = &obs2n[0];
#endif
	//b = Matrix<double, 8, 1>::Zero();
	A = Matrix<double, 8, 8>::Zero();
	//Constructing matrix
	for (i = 0; i < 4; i++) {
		/*cout <<"Hobs = "<< Hobs[i] << endl ;
		cout <<"obs = "<< obs[i] << endl ;*/
		//A = [ -y1, -y2, -1,   0,   0,  0, x1*y1, x1*y2]
		//    [   0,   0,  0, -y1, -y2, -1, x2*y1, x2*y2]
		A(i * 2, 0) = -(obs[i](0, 0));
		A(i * 2, 1) = -(obs[i](1, 0));
		A(i * 2, 2) = -1;
		A(i * 2, 6) = (Hobs[i](0, 0)) * (obs[i](0, 0));
		A(i * 2, 7) = (obs[i](1, 0)) * (Hobs[i](0, 0));

		A(i * 2 + 1, 3) = -(obs[i](0, 0));
		A(i * 2 + 1, 4) = -(obs[i](1, 0));
		A(i * 2 + 1, 5) = -1;
		A(i * 2 + 1, 6) = (Hobs[i](1, 0)) * (obs[i](0, 0));
		A(i * 2 + 1, 7) = (Hobs[i](1, 0)) * (obs[i](1, 0));
		//A(i * 2 + 1, 8) = Hobs[i](1, 0);

		// not sure of the following. Hobs[i](2, 0) can either be used to dehomogeneize or as a b product ...
		b(i * 2, 0) = -(Hobs[i](0, 0));
		b(i * 2 + 1, 0) = -(Hobs[i](1, 0));
	}

	//cout << "A = " << A << endl;
	//cout << A.jacobiSvd(ComputeFullV | ComputeFullU).solve(b) << endl ;
	X = A.jacobiSvd(ComputeFullV | ComputeFullU).solve(b);
#ifdef NORMALIZE
	Matrix<double, 3, 3> H_normed;
	H_normed(0, 0) = X(0, 0);
	H_normed(0, 1) = X(1, 0);
	H_normed(0, 2) = X(2, 0);
	H_normed(1, 0) = X(3, 0);
	H_normed(1, 1) = X(4, 0);
	H_normed(1, 2) = X(5, 0);
	H_normed(2, 0) = X(6, 0);
	H_normed(2, 1) = X(7, 0);
	H_normed(2, 2) = 1;
	denormalize_H(H_normed, H_norm1, H_norm2, H);
	H = H / H(2, 2);
	return 1;
#endif

	//Since we performed the DLT on the 8x8 matrix, the matrix is computed for a homogeneous case
	H(0, 0) = X(0, 0);
	H(0, 1) = X(1, 0);
	H(0, 2) = X(2, 0);
	H(1, 0) = X(3, 0);
	H(1, 1) = X(4, 0);
	H(1, 2) = X(5, 0);
	H(2, 0) = X(6, 0);
	H(2, 1) = X(7, 0);
	H(2, 2) = 1;

	return 1;
}
