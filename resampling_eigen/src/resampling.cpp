#include <stdio.h>
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/core/core.hpp"
#include "opencv2/core/types.hpp"
#include "opencv2/imgproc/imgproc.hpp"

#include <iostream>
#include <Eigen/Dense>

#include "dlt_plane.hpp"

extern "C" {
#include "image/fast/fast.h"
#include "image/brief.h"
}

using namespace Eigen;
using namespace cv;
using namespace std;

/*
 * 627.484211 0.000000 301.426124
 0.000000 626.583268 249.654026
 0.000000 0.000000 1.000000
 */
Matrix<double, 3, 3> K;
Matrix<double, 4, 4> cam_1_mat;
Matrix<double, 4, 4> cam_2_mat;
Matrix<double, 4, 4> cam_2_mk;
Matrix<double, 4, 4> cam_1_mk;

void to_normalized_image_plane(double u, double v, Matrix<double, 3, 3> K,
		Matrix<double, 3, 1> & normed_coord) {
	normed_coord(0, 0) = (u - K(0, 2)) / K(0, 0); //vers plan image normalisé avec origine en axe principale et unité métrique
	normed_coord(1, 0) = (v - K(1, 2)) / K(1, 1);
	normed_coord(2, 0) = 1.;
}

void homography(unsigned char * image_data, int w, int h, unsigned char * img_h,
		Matrix<double, 3, 3> H, Matrix<double, 3, 3> K1,
		Matrix<double, 3, 3> K2) {
	int u, v;
	memset(img_h, 0, w * h);
	//p_a = Hba*p_b -> dans le plan image normalisé
	for (u = 0; u < w; u++) {
		for (v = 0; v < h; v++) {
			Matrix<double, 3, 1> pixel;
			Matrix<double, 3, 1> pixelp;
			to_normalized_image_plane(u, v, K1, pixel);
			/*double xn = (u - K1(0, 2)) / K1(0, 0); //vers plan image normalisé avec origine en axe principale et unité métrique
			 double yn = (v - K1(1, 2)) / K1(1, 1);
			 */
			pixelp = H * pixel;
			/*double xnp = xn * H(0, 0) + yn * H(0, 1) + H(0, 2);
			 double ynp = xn * H(1, 0) + yn * H(1, 1) + H(1, 2);
			 double wnp = xn * H(2, 0) + yn * H(2, 1) + H(2, 2);

			 double xnp_norm = xnp / wnp;
			 double ynp_norm = ynp / wnp;*/
			pixelp = pixelp / pixelp(2, 0);

			/*double xnp_norm_img = xnp_norm * K2(0, 0) + K2(0, 2); //vers plan image normalisé avec origine en axe principale et unité métrique
			 double ynp_norm_img = ynp_norm * K2(1, 1) + K2(1, 2);*/
			double xnp_norm_img = pixelp(0, 0) * K2(0, 0) + K2(0, 2);
			double ynp_norm_img = pixelp(1, 0) * K2(1, 1) + K2(1, 2);

			/*int hu = round(xnp_norm_img + 0.5);
			 int hv = round(ynp_norm_img + 0.5);*/

			int hu = floor(xnp_norm_img);
			int hv = floor(ynp_norm_img);

			if (hu >= 0 && hv >= 0 && hu < w && hv < h) {

				float pixel_val1, pixel_val2, pixel_val;
				float pixa, pixb, pixc, pixd;

				if (xnp_norm_img < hu)
					cout << "problem" << endl;
				if (ynp_norm_img < hv)
					cout << "problem" << endl;
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

				//img_h[u + (v * w)] = image_data[hu + (hv * w)];

				//should add interpolation instead of nearest neighbor ...
			}

		}
	}
}

/**
 * Poses are expressed in world frame whose origin is center of the scene on the ground plane
 * n is plane normal in world frame, d is plane distance to origin in world frame
 */
void compute_homography_from_cam_cam_pos(Matrix<double, 4, 4> C1_pose,
		Matrix<double, 4, 4> C2_pose, Matrix<double, 3, 1> n, double d,
		Matrix<double, 3, 3> K1, Matrix<double, 3, 3> K2,
		Matrix<double, 3, 3> & H) {
	// poses are express as a camera pose |R t| a 4x4 matrix
	//                                    |0 1|

	Matrix<double, 3, 3> R_C1, R_C2, R_C2C1;
	Matrix<double, 3, 1> tC1, tC2, tC2_C2C1;
	Matrix<double, 3, 3> t_C2_C2C1_nC1t;
	Matrix<double, 4, 4> C1_pose_struct;

	Matrix<double, 4, 4> C1_pose_inv;
	Matrix<double, 4, 4> C2C1_pose;
	Matrix<double, 3, 1> nC1;
	Matrix<double, 1, 1> dC1_mat;
	double dC1;

	//cout << C1_pose << endl;
	C1_pose_inv = C1_pose.inverse();
	//cout << C1_pose_inv << endl;

	C2C1_pose = C2_pose * C1_pose_inv;
	//cout << C2C1_pose << endl;

	R_C1 = C1_pose.topLeftCorner(3, 3);
	tC1 = C1_pose.topRightCorner(3, 1);
	R_C2C1 = C2C1_pose.topLeftCorner(3, 3);
	//cout << "R_C2C1 = " << R_C2C1 << endl;
	tC2_C2C1 = C2C1_pose.topRightCorner(3, 1);
	//cout << "tC2_C2C1 = " << tC2_C2C1 << endl;
	nC1 = R_C1 * n;
	//cout << "nC1 = " << nC1 << endl;
	dC1_mat = nC1.transpose() * tC1;
	//cout << "dC1_mat = " << dC1_mat << endl;
	t_C2_C2C1_nC1t = tC2_C2C1 * nC1.transpose();
	dC1 = -dC1_mat(0, 0);
	t_C2_C2C1_nC1t = t_C2_C2C1_nC1t / dC1;
	H = R_C2C1 - t_C2_C2C1_nC1t;

	H = H / H(2, 2); //
	cout << "H = " << H << endl;

}

void detect_and_draw_corners(Mat & img) {
	xy * corners;
	int nb_corners, i;
	corners = fast9_detect_nonmax(img.data, img.cols, img.rows, img.cols, 80,
			&nb_corners);
	for (i = 0; i < nb_corners; i++) {

		circle(img, Point(corners[i].x, corners[i].y), 3, Scalar(255), 1, 8, 0);
	}
}

void get_descriptors(Mat & img, xy * corners, brief_descriptor ** descs,
		unsigned int nb_corners) {
	int i;
	gray_image img_gray;
	img_gray.width = img.cols;
	img_gray.height = img.rows;
	img_gray.width_step = img.cols;
	img_gray.imageData = (unsigned char *) img.data;
	(*descs) = (brief_descriptor *) malloc(
			nb_corners * sizeof(brief_descriptor));
	for (i = 0; i < nb_corners; i++) {
		computeBriefDescriptor(&img_gray, corners[i].x, corners[i].y,
				(unsigned char *) ((*descs)[i]));
	}
}

void match_descriptors(brief_descriptor * desc_1, brief_descriptor * desc_2,
		unsigned int nb_corners_1, unsigned int nb_corners_2,
		int ** match_list) {
	unsigned int i, j;
	(*match_list) = (int *) malloc(nb_corners_1 * sizeof(int));
	for (i = 0; i < nb_corners_1; i++) {
		(*match_list)[i] = -1; //no match so far
		unsigned int min_dist = 0xFFFFFFFF;
		for (j = 0; j < nb_corners_2; j++) {
			unsigned int dist = hammingDist(desc_1[i], desc_2[j]);
			if (dist < min_dist) {
				min_dist = dist;
				(*match_list)[i] = j;
			}
		}
	}

}

int main(int argc, char ** argv) {

	cam_1_mat << -1.141065701958907e+03, 3.421466553556007e+02, -1.835992900918554e+03, 5.436898725066867e+05, 7.996462861460456e+02, 2.013181974660253e+03, -2.017580480078388e+02, 1.094990499809007e+05, 7.010124343170585e-01, -1.711311329622356e-01, -6.923118533319601e-01, 1.066915525986577e+03, 0, 0, 0, 1;
	cam_2_mat << -1.018527339744041e+03, -6.185246750592945e+02, -1.835754473429942e+03, 5.052096340773842e+05, -9.546670514338936e+02, 1.944093951756833e+03, -2.053087519543848e+02, 2.403050209771340e+05, 5.980063934130443e-01, 4.034995629891167e-01, -6.925145890916056e-01, 1.091155473093363e+03, 0, 0, 0, 1;

	cam_1_mk << -6.654672833894143e-01, 1.920394188556554e-01, -7.213003232661797e-01, 4.813190380466447e+01, 2.563881075073194e-01, 9.663513837816159e-01, 2.073985034110462e-02, -1.258121798534581e+02, 7.010124343170585e-01, -1.711311329622356e-01, -6.923118533319601e-01, 1.066915525986577e+03, 0, 0, 0, 1;

	cam_2_mk << -5.886806295004532e-01, -3.652356943184633e-01, -7.211504725413832e-01, 2.557518376821454e+01, -5.439149472923810e-01, 8.389225174372792e-01, 1.911909643849757e-02, -6.888416876240703e+01, 5.980063934130443e-01, 4.034995629891167e-01, -6.925145890916056e-01, 1.091155473093363e+03, 0, 0, 0, 1;

	K << 2.149349810656288e+03, 0, 4.126264577230862e+02, 0, 2.146276553276586e+03, 3.557233656020202e+02, 0, 0, 1.000000000000000e+00;

	Mat img_a;
	img_a = imread("/home/jpiat/Pictures/imaobj2_rect.jpg", IMREAD_GRAYSCALE);
	Mat img_compo(img_a.rows, img_a.cols,
	CV_8UC1);
	Mat img_compo_2(img_a.rows, img_a.cols,
	CV_8UC1);

	Matrix<double, 3, 3> H;
	Matrix<double, 3, 1> n;
	brief_descriptor * desc_1, *desc_2;
	xy * corners_1, corners_2;
	int nb_corners_1, nb_corners_2;
	int * match_list;
	n << 0, 0, 1;
	double d = -104.;
	compute_homography_from_cam_cam_pos(cam_1_mk, cam_2_mk, n, d, K, K, H);

	homography((unsigned char *) img_a.data, img_compo_2.cols, img_compo_2.rows,
			(unsigned char *) img_compo_2.data, H, K, K);
	corners_1 = fast9_detect_nonmax(img_a.data, img_a.cols, img_a.rows,
			img_a.cols, 80, &nb_corners_1);

	//detect_and_draw_corners(img_a);
	get_descriptors(img_a, corners_1, &desc_1, nb_corners_1);
	get_descriptors(img_a, corners_1, &desc_2, nb_corners_1);
	match_descriptors(desc_1, desc_2, nb_corners_1, nb_corners_1, &match_list);
	imshow("input", img_a);
	imshow("result", img_compo_2);
	waitKey(0);

	cout << "test DLT : " << endl;

	Matrix<double, 3, 1> obs1[4];
	Matrix<double, 3, 1> obs2[4];
	Matrix<double, 3, 3> H_calc;

	obs2[0] << 0, 0, 1;
	obs2[1] << 1, 0, 1;
	obs2[2] << 0, 1, 1;
	obs2[3] << 1, 1, 1;

	obs1[0] = H * obs2[0];
	obs1[1] = H * obs2[1];
	obs1[2] = H * obs2[2];
	obs1[3] = H * obs2[3];

	cout << "obs[0] = " << obs1[0] << endl;
	cout << "obs[1] = " << obs1[1] << endl;
	cout << "obs[2] = " << obs1[2] << endl;
	cout << "obs[3] = " << obs1[3] << endl;
	/*	cout << "obs 1 :" << obs1[0] << endl;
	 cout << "obs 2 :" << obs1[1] << endl;
	 cout << "obs 3 :" << obs1[2] << endl;
	 cout << "obs 4 :" << obs1[3] << endl;
	 */
	estimate_homo_dlt(obs2, obs1, H_calc);

	cout << "applying transformation to dataset" << endl;

	cout << "H*obs[0] = " << (H_calc * obs2[0]) << endl;
	cout << "H*obs[1] = " << (H_calc * obs2[1]) << endl;
	cout << "H*obs[2] = " << (H_calc * obs2[2]) << endl;
	cout << "H*obs[3] = " << (H_calc * obs2[3]) << endl;

	cout << "H_calc  = " << H_calc << endl;

	return 1;
}
