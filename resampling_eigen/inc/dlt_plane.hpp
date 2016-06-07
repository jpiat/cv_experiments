#include <Eigen/Dense>
#include <Eigen/SVD>
#include <iostream>

using namespace Eigen;
using namespace std;


#ifndef __DLT_PLANE_H__
#define __DLT_PLANE_H__

/**
 * Estimate homography based on 4-points observations
 * Observations must be provided in homogeneous coordinates
 */
int estimate_homo_dlt(Matrix<double, 3, 1> * obs, Matrix<double, 3, 1> * Hobs,
		Matrix<double, 3, 3> & H);

#endif
