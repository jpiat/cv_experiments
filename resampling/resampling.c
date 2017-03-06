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

double K_slam[9] = { 674.148669, 0, 312.668285, 0, 674.148669, 223.832270, 0, 0,
		1. };

double K_virt[9] = { 320, 0, 320, 0, 240, 240, 0, 0, 1. };

// à redéfinir avec la poser correcte
double cam_1_C[12] = { 1, 0, 0, 0, 0, 0, -1, 1., 0, 1, 0, 0 };
double cam_1_Ctilde[12];

#define THETA_CAM_BIRD -(M_PI)
//double cam_bird_Rtw[16] = {};

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

void homography(char * image_data, int w, int h, char * img_h, double * H,
		double * K1, double * K2) {
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
			int hu = floor(xnp_norm_img);
			int hv = floor(ynp_norm_img);

			if (hu >= 0 && hv >= 0 && hu < w && hv < h) {

				float pixel_val1, pixel_val2, pixel_val;
				float pixa, pixb, pixc, pixd;

				/*if (xnp_norm_img < hu)
				 cout << "problem" << endl;
				 if (ynp_norm_img < hv)
				 cout << "problem" << endl;*/
				pixa = image_data[(hu) + (hv * w)];
				pixb = image_data[(hu + 1) + (hv * w)];

				pixc = image_data[(hu) + ((hv + 1) * w)];
				pixd = image_data[(hu + 1) + ((hv + 1) * w)];

				pixel_val1 = pixa * (1 - (xnp_norm_img - hu));
				pixel_val1 += pixb * ((xnp_norm_img - hu));

				pixel_val2 = pixc * (1 - (xnp_norm_img - hu));
				pixel_val2 += pixd * ((xnp_norm_img - hu));

				pixel_val = pixel_val1 * (1 - (ynp_norm_img - hv));
				pixel_val += pixel_val2 * (ynp_norm_img - hv);

				img_h[u + (v * w)] = round(pixel_val);
			}
		}
	}
}

void calc_ct_and_H(double * world_to_cam, double * K, double *Ct, double * H) {
	int i;
	struct double_matrix world_to_cam_struct, K_struct, Ct_struct;
	alloc_double_matrix(&world_to_cam_struct, 4, 3, world_to_cam);
	alloc_double_matrix(&K_struct, 3, 3, K);
	alloc_double_matrix(&Ct_struct, 4, 3, Ct);

	double_mat_product(&K_struct, &world_to_cam_struct, &Ct_struct, 0, 0);

	print_double_matrix(Ct_struct, "Ct");
	H[0] = Ct[0];
	H[1] = Ct[1];
	H[2] = Ct[3];
	H[3] = Ct[4];
	H[4] = Ct[5];
	H[5] = Ct[7];
	H[6] = Ct[8];
	H[7] = Ct[9];
	H[8] = Ct[11];
	printf("H =");
	for(i = 0 ; i < 9 ; i ++) printf("%lf, ", H[i]);
}
void pixel_to_ground_plane(double * Ct, double u, double v, double * x,
		double * y) {
	struct double_matrix Ct_struct;
	alloc_double_matrix(&Ct_struct, 4, 3, Ct);

	double a1 = MAT_ELT_AT(Ct_struct, 0, 0) - u * MAT_ELT_AT(Ct_struct, 2, 0);
	double b1 = MAT_ELT_AT(Ct_struct, 0, 1) - u * MAT_ELT_AT(Ct_struct, 2, 1);
	double c1 = MAT_ELT_AT(Ct_struct, 0, 2) - u * MAT_ELT_AT(Ct_struct, 2, 2);

	double a2 = MAT_ELT_AT(Ct_struct, 1, 0) - v * MAT_ELT_AT(Ct_struct, 2, 0);
	double b2 = MAT_ELT_AT(Ct_struct, 1, 1) - v * MAT_ELT_AT(Ct_struct, 2, 1);
	double c2 = MAT_ELT_AT(Ct_struct, 1, 2) - v * MAT_ELT_AT(Ct_struct, 2, 2);

	double d1 = -(MAT_ELT_AT(Ct_struct, 0, 3) - u * MAT_ELT_AT(Ct_struct, 2, 3));
	double d2 = -(MAT_ELT_AT(Ct_struct, 1, 3) - v * MAT_ELT_AT(Ct_struct, 2, 3));

	double b3 = a2 * b1 - a1 * b2;
	double c3 = a2 * c1 - a1 * c2;
	double d3 = a2 * d1 - a1 * d2;

	(*y) = d3 / b3;
	(*x) = (d1 - ((*y) * b1)) / a1;

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
	return 1;
}

void map_from_homography(char * img_in, unsigned int w_in, unsigned int h_in, char * img_out, unsigned int w_out, unsigned int h_out, double * H, double scale){
unsigned int u, v ;
	for (u = 0; u < w_out; u++) {
			for (v = 0; v < h_out; v++) {

				//double xn = 0.997535, yn = 2.060746 ;
				double xn =(v-(h_out/2.))/scale ;
				double yn = u/scale ;



				double xnp = xn * H[0] + yn * H[1] + H[2];
				double ynp = xn * H[3] + yn * H[4] + H[5];
				double wnp = xn * H[6] + yn * H[7] + H[8];

				double xnp_norm = xnp / wnp;
				double ynp_norm = ynp / wnp;

				int hu = round(xnp_norm);
				int hv = round(ynp_norm);

				if (hu >= 0 && hv >= 0 && hu <  w_in && hv < h_in) {
					unsigned char pixel_val = img_in[(hu) + ((hv) *  w_in)];
					/*unsigned char pixel_val = ((hu/50)%2 ^ (hv)%2)*255;
					parking_image->imageData[(hu) + ((hv) *  parking_image->width)] = 0 ;*/
					img_out[u + (v * w_out)] = pixel_val;
				}else{
					img_out[u + (v * w_out)] = 128 ;
				}
			}
		}
}


int main(int argc, char ** argv) {
	int u, v;
	IplImage * img_a = cvLoadImage("/home/jpiat/Pictures/imaobj2_rect.jpg",
			CV_LOAD_IMAGE_GRAYSCALE);

	IplImage * parking_image = cvLoadImage(
			"/home/jpiat/Pictures/parking_laas/in01/pgm/image_0008.pgm",
			CV_LOAD_IMAGE_GRAYSCALE);
	IplImage * img_map = cvCreateImage(cvSize(3000, 3000), IPL_DEPTH_8U, 1);
	IplImage * img_compo = cvCreateImage(cvSize(img_a->width, img_a->height),
	IPL_DEPTH_8U, 1);
	IplImage * img_compo_2 = cvCreateImage(cvSize(img_a->width, img_a->height),
	IPL_DEPTH_8U, 1);

	double uvc1[2], uvc2[2], H[3 * 3];
	double n[3] = { 0, 0, 1 };
	double H_bvdp[9];
	double Ct_bvdp[12];
	double d = 0.; //-104.;
	/*struct double_matrix cam_1_Rtw_struct, cam_bird_Rtw_struct;
	 struct double_matrix cam_1_Rtc_struct, cam_bird_Rtc_struct;
	 alloc_double_matrix(&cam_1_Rtw_struct, 4, 4, cam_1_Rtw);
	 alloc_double_matrix(&cam_bird_Rtw_struct, 4, 4, cam_bird_Rtw);
	 alloc_double_matrix(&cam_1_Rtc_struct, 4, 4, cam_1_Rtc);
	 alloc_double_matrix(&cam_bird_Rtc_struct, 4, 4, cam_bird_Rtc);

	 print_double_matrix(cam_1_Rtw_struct, "Cam1_Rtw");
	 double_mat_inv(&cam_1_Rtw_struct, &cam_1_Rtc_struct);
	 print_double_matrix(cam_1_Rtc_struct, "Cam1_Rtc");
	 double_mat_inv(&cam_bird_Rtw_struct, &cam_bird_Rtc_struct);*/
	double x, y;

	calc_ct_and_H(cam_1_C, K_slam, Ct_bvdp, H_bvdp);

	/*intersect_ground_plane(cam_1_C, K_virt, 0, 0, &x, &y);
	 printf("%lf, %lf \n", x, y);*/
	pixel_to_ground_plane(Ct_bvdp, 0, 479, &x, &y);
	printf("%lf, %lf \n", x, y);
	pixel_to_ground_plane(Ct_bvdp, 639, 479, &x, &y);
	printf("%lf, %lf \n", x, y);
	/*intersect_ground_plane(cam_1_C, K_virt, 639, 0, &x, &y);
	 printf("%lf, %lf \n", x, y);*/

	pixel_to_ground_plane(Ct_bvdp, 0, 478, &x, &y);
	printf("%lf, %lf \n", x, y);
	pixel_to_ground_plane(Ct_bvdp, 639, 478, &x, &y);
	printf("%lf, %lf \n", x, y);
	//memset(parking_image->imageData, 255, parking_image->width*parking_image->height);
	map_from_homography(parking_image->imageData,parking_image->width, parking_image->height, img_map->imageData, img_map->width, img_map->height, H_bvdp, 100.);

		cvShowImage("img", parking_image);
		cvShowImage("map", img_map);
		cvSaveImage("./test_map.png", img_map, NULL);
		//cvWaitKey(0);
	/*compute_homography_from_cam_cam_pos(cam_1_Rtw, cam_bird_Rtw, n, d, K_slam,
	 K_slam, H);*/

	return 0;

	compute_homography_from_cam_cam_pos(cam_1_mk, cam_2_mk, n, d,
			cam_bertrand_K, cam_bertrand_K, H);
	homography(img_a->imageData, img_compo_2->width, img_compo_2->height,
			img_compo_2->imageData, H, cam_bertrand_K, cam_bertrand_K);
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
