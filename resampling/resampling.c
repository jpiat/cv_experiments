#include <stdio.h>
#include "opencv2/highgui/highgui_c.h"
#include "opencv2/core/core_c.h"
#include "opencv2/core/types_c.h"
#include "opencv2/imgproc/imgproc_c.h"

#include "mat_types.h"
#include "double_mat_ops.h"

/*
 * 627.484211 0.000000 301.426124
 0.000000 626.583268 249.654026
 0.000000 0.000000 1.000000
 */

double poly[] = { 0.095880, -0.275158, 0.000727, -0.007979, 0.000000 };
//double poly []= {-0.1, 0., 0., 0., 0.};

double K[9] = { 627.484211, 0, 301.426124, 0, 626.583268, 249.654026, 0, 0, 1. };

double H_translate[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 0.5 };

double max_diff = 0.0;

double cam_1_mat[16] = { -1.141065701958907e+03, 3.421466553556007e+02,
		-1.835992900918554e+03, 5.436898725066867e+05, 7.996462861460456e+02,
		2.013181974660253e+03, -2.017580480078388e+02, 1.094990499809007e+05,
		7.010124343170585e-01, -1.711311329622356e-01, -6.923118533319601e-01,
		1.066915525986577e+03, 0, 0, 0, 1 };
double cam_2_mat[16] = { -1.018527339744041e+03, -6.185246750592945e+02,
		-1.835754473429942e+03, 5.052096340773842e+05, -9.546670514338936e+02,
		1.944093951756833e+03, -2.053087519543848e+02, 2.403050209771340e+05,
		5.980063934130443e-01, 4.034995629891167e-01, -6.925145890916056e-01,
		1.091155473093363e+03, 0, 0, 0, 1 };

double cam_1_mk[16] = { -6.654672833894143e-01, 1.920394188556554e-01,
		-7.213003232661797e-01, 4.813190380466447e+01, 2.563881075073194e-01,
		9.663513837816159e-01, 2.073985034110462e-02, -1.258121798534581e+02,
		7.010124343170585e-01, -1.711311329622356e-01, -6.923118533319601e-01,
		1.066915525986577e+03, 0, 0, 0, 1 };

double cam_2_mk[16] = { -5.886806295004532e-01, -3.652356943184633e-01,
		-7.211504725413832e-01, 2.557518376821454e+01, -5.439149472923810e-01,
		8.389225174372792e-01, 1.911909643849757e-02, -6.888416876240703e+01,
		5.980063934130443e-01, 4.034995629891167e-01, -6.925145890916056e-01,
		1.091155473093363e+03, 0, 0, 0, 1 };

double cam_bertrand_K[9] = { 2.149349810656288e+03, 0, 4.126264577230862e+02, 0,
		2.146276553276586e+03, 3.557233656020202e+02, 0, 0,
		1.000000000000000e+00 };

double distort_plumb_bob(double xn, double yn, double * xd, double * yd,
		double * poly) {
	double r_square = xn * xn + yn * yn;
	double k = 1.0 + r_square * poly[0] + (r_square * r_square) * poly[1]
			+ (r_square * r_square) * r_square * poly[4];
	double tx = 2 * poly[2] * (xn * yn)
			+ poly[3] * (r_square + (2 * (xn * xn)));
	double ty = 2 * poly[3] * (xn * yn)
			+ poly[2] * (r_square + (2 * (yn * yn)));
	(*xd) = k * xn + tx;
	(*yd) = k * yn + ty;
	return sqrt(pow((*xd) - xn, 2) + pow((*yd) - yn, 2));
}

void undistort_image(char * img_data, int w, int h, char * img_u, int wu,
		int hu) {
	int u, v;
	for (u = 0; u < wu; u++) {
		for (v = 0; v < hu; v++) {

			double xn = (u - K[2]) / K[0]; //vers plan image normalisé avec origine en axe principale et unité métrique
			double yn = (v - K[5]) / K[4];
			double xd, yd;
			double diff;
			diff = distort_plumb_bob(xn, yn, &xd, &yd, poly);
			diff = diff * K[0];
			int udi = (int) ((xd * K[0] + K[2]) + 0.5); //vers le plan image pixel de l'image distordu, avec rounding
			int vdi = (int) ((yd * K[4] + K[5]) + 0.5);
			max_diff = (diff > max_diff) ? diff : max_diff;
			//img_u[u + (v * wu)] = ((int)(diff * (256.0/20.0))) ;
			if (udi >= 0 && udi < w && vdi >= 0 && vdi < h) {

				img_u[u + (v * wu)] = img_data[udi + (vdi * w)];
			} else {
				img_u[u + (v * wu)] = 0;
			}
		}
	}
}

void homography(char * image_data, int w, int h, char * img_h, double * H, double * K1, double * K2) {
	int u, v;
	memset(img_h, 0, w * h);
	//p_a = Hba*p_b -> dans le plan image normalisé
	for (u = 0; u < w; u++) {
		for (v = 0; v < h; v++) {
			double xn = (u - K1[2]) / K1[0]; //vers plan image normalisé avec origine en axe principale et unité métrique
			double yn = (v - K1[5]) / K1[4];

			double xnp = xn * H[0] + yn * H[1] + H[2];
			double ynp = xn * H[3] + yn * H[4] + H[5];
			double wnp = xn * H[6] + yn * H[7] + H[8];

			double xnp_norm = xnp / wnp;
			double ynp_norm = ynp / wnp;

			double xnp_norm_img = xnp_norm * K2[0] + K2[2]; //vers plan image normalisé avec origine en axe principale et unité métrique
			double ynp_norm_img = ynp_norm * K2[4] + K2[5];
			int hu = round(xnp_norm_img + 0.5);
			int hv = round(ynp_norm_img + 0.5);

			if (hu >= 0 && hv >= 0 && hu < w && hv < h) {
				img_h[u + (v * w)] = image_data[hu + (hv * w)];

				//should add interpolation instead of nearest neighbor ...
			}

		}
	}
}

/**
 * Poses are expressed in world frame whose origin is center of the scene on the ground plane
 * n is plane normal in world frame, d is plane distance to origin in world frame
 */
void compute_homography_from_cam_cam_pos(double * C1_pose, double * C2_pose,
		double * n, double d, double * K1, double * K2, double *H) {
	// poses are express as a camera pose |R t| a 4x4 matrix
	//                                    |0 1|

	int i, j;
	struct double_matrix R_C1, R_C2, tC1, tC2, tC2_C2C1, R_C2C1, t_C2_C2C1_nC1t;
	struct double_matrix C1_pose_struct, C2_pose_struct, C1_pose_struct_inv,
			C2C1_pose;
	struct double_matrix nC1, n_struct, dC1_mat, H_struct;
	struct double_matrix K1_struct, K2_struct, K1H, K2_inv;
	double dC1;

	if (alloc_double_matrix(&R_C1, 3, 3, NULL) < 0)
		printf("Allocation problem \n");
	if (alloc_double_matrix(&R_C2, 3, 3, NULL) < 0)
		printf("Allocation problem \n");
	if (alloc_double_matrix(&R_C2C1, 3, 3, NULL) < 0)
		printf("Allocation problem \n");
	if (alloc_double_matrix(&tC2_C2C1, 1, 3, NULL) < 0)
		printf("Allocation problem \n");
	if (alloc_double_matrix(&nC1, 1, 4, NULL) < 0)
		printf("Allocation problem \n");
	if (alloc_double_matrix(&tC1, 1, 3, NULL) < 0)
		printf("Allocation problem \n");
	if (alloc_double_matrix(&tC2, 1, 3, NULL) < 0)
		printf("Allocation problem \n");
	if (alloc_double_matrix(&dC1_mat, 1, 1, NULL) < 0)
		printf("Allocation problem \n");
	if (alloc_double_matrix(&t_C2_C2C1_nC1t, 3, 3, NULL) < 0)
		printf("Allocation problem \n");
	if (alloc_double_matrix(&C1_pose_struct_inv, 4, 4, NULL) < 0)
		printf("Allocation problem \n");
	if (alloc_double_matrix(&C2C1_pose, 4, 4, NULL) < 0)
		printf("Allocation problem \n");

	if (alloc_double_matrix(&K1H, 3, 3, NULL) < 0)
		printf("Allocation problem \n");
	if (alloc_double_matrix(&K2_inv, 3, 3, NULL) < 0)
		printf("Allocation problem \n");
	alloc_double_matrix(&C1_pose_struct, 4, 4, C1_pose);
	alloc_double_matrix(&C2_pose_struct, 4, 4, C2_pose);
	alloc_double_matrix(&K1_struct, 3, 3, K1);
	alloc_double_matrix(&K2_struct, 3, 3, K2);
	alloc_double_matrix(&n_struct, 1, 4, n);
	alloc_double_matrix(&H_struct, 3, 3, H);

	get_double_sub_matrix(&C1_pose_struct, &R_C1, 0, 0, 3, 3);
	get_double_sub_matrix(&C2_pose_struct, &R_C2, 0, 0, 3, 3);
	get_double_sub_matrix(&C1_pose_struct, &tC1, 3, 0, 1, 3);
	get_double_sub_matrix(&C2_pose_struct, &tC2, 3, 0, 1, 3);

	print_double_matrix(C1_pose_struct, "C1");
	print_double_matrix(C2_pose_struct, "C2");
	double_mat_inv(&C1_pose_struct, &C1_pose_struct_inv);
	print_double_matrix(C1_pose_struct_inv, "C1_inv");
	double_mat_product(&C2_pose_struct, &C1_pose_struct_inv, &C2C1_pose, 0, 0);
	print_double_matrix(C2C1_pose, "C2C1_pose");

	double intra = sqrt(
			pow(MAT_ELT_AT(C2C1_pose, 0, 3), 2)
					+ pow(MAT_ELT_AT(C2C1_pose, 1, 3), 2)
					+ pow(MAT_ELT_AT(C2C1_pose, 2, 3), 2));
	printf("%lf \n", intra);

	get_double_sub_matrix(&C2C1_pose, &R_C2C1, 0, 0, 3, 3);
	get_double_sub_matrix(&C2C1_pose, &tC2_C2C1, 3, 0, 1, 3);

	print_double_matrix(R_C2C1, "R_C2_C1");
	print_double_matrix(tC2_C2C1, "tC2_C2C1");

	print_double_matrix(R_C1, "RC1");
	print_double_matrix(n_struct, "n");
	get_double_sub_matrix(&n_struct, NULL, 0, 0, 1, 3);
	get_double_sub_matrix(&nC1, NULL, 0, 0, 1, 3);
	double_mat_product(&R_C1, &n_struct, &nC1, 0, 0);
	print_double_matrix(nC1, "nC1");
	print_double_matrix(tC1, "tC1");
	get_double_sub_matrix(&nC1, NULL, 0, 0, 1, 3);
	double_mat_product(&nC1, &tC1, &dC1_mat, 1, 0);
	print_double_matrix(dC1_mat, "dC1");
	double_mat_product(&tC2_C2C1, &nC1, &t_C2_C2C1_nC1t, 0, 1);
	print_double_matrix(t_C2_C2C1_nC1t, "t_C2_C2C1_nC1t");
	dC1 = -MAT_ELT_AT(dC1_mat, 0, 0);
	for (i = 0; i < t_C2_C2C1_nC1t.nbr; i++) {
		for (j = 0; j < t_C2_C2C1_nC1t.nbc; j++) {
			MAT_ELT_AT(t_C2_C2C1_nC1t, i , j) /= dC1;
		}
	}

	double_mat_sub(&R_C2C1, &t_C2_C2C1_nC1t, &H_struct);
/*	double_mat_product(&K1_struct, &H_struct, &K1H, 0, 0);
	double_mat_inv(&K2_struct, &K2_inv);
	double_mat_product(&K1H, &K2_inv, &H_struct, 0, 0);*/
	print_double_matrix(H_struct, "H");

	free_double_matrix(&R_C1);
	free_double_matrix(&R_C2);
	free_double_matrix(&R_C2C1);
	free_double_matrix(&tC2_C2C1);
	free_double_matrix(&nC1);
	free_double_matrix(&tC1);
	free_double_matrix(&tC2);
	free_double_matrix(&dC1_mat);
	free_double_matrix(&t_C2_C2C1_nC1t);

}

//Function to compute the projection of a ground pixel in camera 1
/**
 * Camera
 */
int compute_ground_pixel_projection(double * uv_c1, double * uv_c2, double * n,
		double d, double * C1_matrix, double * C2_matrix) {
	// matrices are expressed as a camera |R t| a 4x4 matrix
	//                                    |0 1|
	double u = uv_c1[0];
	double v = uv_c1[1];
	struct double_matrix matA_struct, matA_inv_struct, matX_struct;
	struct double_matrix matB_struct, C2_mat_struct, C2_pixel;
	if (alloc_double_matrix(&matA_struct, 3, 3, NULL) < 0)
		printf("Allocation problem \n");
	if (alloc_double_matrix(&matA_inv_struct, 3, 3, NULL) < 0)
		printf("Allocation problem \n");
	if (alloc_double_matrix(&matB_struct, 1, 3, NULL) < 0)
		printf("Allocation problem \n");
	if (alloc_double_matrix(&matX_struct, 1, 4, NULL) < 0)
		printf("Allocation problem \n");
	if (alloc_double_matrix(&C2_mat_struct, 4, 4, C2_matrix) < 0)
		printf("Allocation problem \n");
	if (alloc_double_matrix(&C2_pixel, 1, 4, NULL) < 0)
		printf("Allocation problem \n");
	MAT_ELT_AT(matA_struct, 0, 0) = (C1_matrix[0 + (0 * 4)])
			- (C1_matrix[0 + (2 * 4)] * u);
	MAT_ELT_AT(matA_struct, 0, 1) = (C1_matrix[1 + (0 * 4)])
			- (C1_matrix[1 + (2 * 4)] * u);
	MAT_ELT_AT(matA_struct, 0, 2) = (C1_matrix[2 + (0 * 4)])
			- (C1_matrix[2 + (2 * 4)] * u);

	MAT_ELT_AT(matB_struct,0, 0) = (C1_matrix[3 + (2 * 4)] * u)
			- (C1_matrix[3 + (0 * 4)]);

	MAT_ELT_AT(matA_struct, 1, 0) = (C1_matrix[0 + (1 * 4)])
			- (C1_matrix[0 + (2 * 4)] * v);
	MAT_ELT_AT(matA_struct, 1, 1) = (C1_matrix[1 + (1 * 4)])
			- (C1_matrix[1 + (2 * 4)] * v);
	MAT_ELT_AT(matA_struct, 1, 2) = (C1_matrix[2 + (1 * 4)])
			- (C1_matrix[2 + (2 * 4)] * v);

	MAT_ELT_AT(matB_struct,0, 1) = (C1_matrix[3 + (2 * 4)] * v)
			- (C1_matrix[3 + (1 * 4)]);

	MAT_ELT_AT(matA_struct, 2, 0) = n[0];
	MAT_ELT_AT(matA_struct, 2, 1) = n[1];
	MAT_ELT_AT(matA_struct, 2, 2) = n[2];

	MAT_ELT_AT(matB_struct,0, 2) = -d;

	MAT_ELT_AT(matX_struct, 0, 3) = 1.;

	get_double_sub_matrix(&matX_struct, NULL, 0, 0, 1, 3);
	//print_double_matrix(matA_struct, "matA");

	//should compute A_inv directly ...

	double_mat_inv(&matA_struct, &matA_inv_struct);
	//print_double_matrix(matA_inv_struct, "matA_inv");
	//print_double_matrix(matB_struct, "matB");
	//if(double_mat_inv(&matA_struct, &matA_inv_struct) <0) return -1;
	if (double_mat_product(&matA_inv_struct, &matB_struct, &matX_struct, 0, 0)
			< 0)
		return -1;

	get_double_sub_matrix(&matX_struct, NULL, 0, 0, 1, 4);
	MAT_ELT_AT(matX_struct, 3, 0) = 1.0;
	//print_double_matrix(matX_struct, "matX");

	//print_double_matrix(C2_mat_struct, "matC2");
	double_mat_product(&C2_mat_struct, &matX_struct, &C2_pixel, 0, 0);
	MAT_ELT_AT(C2_pixel, 0, 0) = MAT_ELT_AT(C2_pixel, 0,
			0) / MAT_ELT_AT(C2_pixel, 2, 0);
	MAT_ELT_AT(C2_pixel, 1, 0) = MAT_ELT_AT(C2_pixel, 1,
			0) / MAT_ELT_AT(C2_pixel, 2, 0);
	MAT_ELT_AT(C2_pixel, 2, 0) = 1.;

	uv_c2[0] = MAT_ELT_AT(C2_pixel, 0, 0);
	uv_c2[1] = MAT_ELT_AT(C2_pixel, 1, 0);

	//print_double_matrix(C2_pixel, "matC2_pixel");

	free_double_matrix(&matA_struct);
	free_double_matrix(&matA_inv_struct);
	free_double_matrix(&matB_struct);
	free_double_matrix(&matX_struct);
	free_double_matrix(&C2_pixel);
	return 1 ;
}

int main(int argc, char ** argv) {
	int u, v;
	IplImage * img_a = cvLoadImage("/home/jpiat/Pictures/imaobj2_rect.jpg",
			CV_LOAD_IMAGE_GRAYSCALE);
	IplImage * img_compo = cvCreateImage(cvSize(img_a->width, img_a->height),
	IPL_DEPTH_8U, 1);
	IplImage * img_compo_2 = cvCreateImage(cvSize(img_a->width, img_a->height),
	IPL_DEPTH_8U, 1);

	double uvc1[2], uvc2[2], H[3 * 3];
	double n[3] = { 0, 0, 1 };
	double d = 0.; //-104.;
	compute_homography_from_cam_cam_pos(cam_1_mk, cam_2_mk, n, d, cam_bertrand_K, cam_bertrand_K, H);
	homography(img_a->imageData, img_compo_2->width, img_compo_2->height, img_compo_2->imageData, H, cam_bertrand_K, cam_bertrand_K);
	/*for (u = 0; u < img_a->width; u++) {
		for (v = 0; v < img_a->height; v++) {
			uvc1[0] = u;
			uvc1[1] = v;
			compute_ground_pixel_projection(uvc1, uvc2, n, d, cam_1_mat,
					cam_2_mat);
			if (uvc2[0] >= 0 && uvc2[0] < img_compo->width && uvc2[1] >= 0
					&& uvc2[1] < img_compo->height) {
				img_compo->imageData[((int) uvc1[1]) * img_compo->widthStep
						+ ((int) uvc1[0])] = img_a->imageData[((int) uvc2[1])
						* img_a->widthStep + ((int) uvc2[0])];
			}
		}
	}*/
	cvShowImage("orig", img_a);
	cvShowImage("compo", img_compo);
	cvShowImage("compo_2", img_compo_2);
	cvWaitKey(0);
	return 0;

	IplImage * img_orig = cvLoadImage(argv[1], CV_LOAD_IMAGE_GRAYSCALE);
	IplImage * img_und = cvCreateImage(
			cvSize(img_orig->width, img_orig->height), IPL_DEPTH_8U, 1);
	IplImage * img_homo = cvCreateImage(
			cvSize(img_orig->width, img_orig->height), IPL_DEPTH_8U, 1);

	IplImage * sub = cvCreateImage(cvSize(img_orig->width, img_orig->height),
	IPL_DEPTH_8U, 1);

	undistort_image(img_orig->imageData, img_orig->width, img_orig->height,
			img_und->imageData, img_und->width, img_und->height);

	printf("%lf \n", max_diff);

	//compute horizontal and vertical FOV
	float fovx = 2 * atan2((0.5 * img_orig->width), K[0]);
	float fovy = 2 * atan2(0.5 * img_orig->height, K[4]);
	printf("fov X : %f\n", fovx * (180.0 / M_PI));
	printf("fov Y : %f\n", fovy * (180.0 / M_PI));
	cvSub(img_orig, img_und, sub, NULL);
	cvShowImage("distort", img_orig);
	cvShowImage("undistort", img_und);
	cvShowImage("sub", sub);
	double theta = 0;
	while (1) {
		theta += M_PI / 100;
		double H_rotate_inc[9] = { 1, 0, 0, 0, cos(theta), -sin(theta), 0, sin(
				theta), cos(theta) };
		homography(img_und->imageData, img_orig->width, img_orig->height,
				img_homo->imageData, H_rotate_inc, K, K);
		if (theta >= M_PI / 4) {
			theta = 0;
		}
		usleep(10000);
		cvShowImage("homo", img_homo);
		cvWaitKey(1);
	}

	return 1;
}
