#include <string.h>
#include <stdio.h>

#include "double_mat_ops.h"
#include "double_ekf.h"

int double_ekf_init(struct double_ekf_state * e, unsigned int state_size,
		unsigned int observed_size) {
	e->state_size = state_size;
	e->state_observable = observed_size;

	if ((e->x = malloc(e->state_size * sizeof(double))) == NULL)
		return -1;
	if ((e->fx = malloc(e->state_size * sizeof(double))) == NULL)
		return -1;
	if ((e->hx = malloc(e->state_size * sizeof(double))) == NULL)
		return -1;

	if ((e->P = malloc(sizeof(struct double_matrix))) == NULL)
		return -1;
	if ((alloc_double_matrix(e->P, e->state_size, e->state_size, NULL)) < 0)
		return -1;
	if ((e->Q = malloc(sizeof(struct double_matrix))) == NULL)
		return -1;
	if ((alloc_double_matrix(e->Q, e->state_size, e->state_size, NULL)) < 0)
		return -1;
	if ((e->R = malloc(sizeof(struct double_matrix))) == NULL)
		return -1;
	if ((alloc_double_matrix(e->R, e->state_observable, e->state_observable,
	NULL)) < 0)
		return -1;
	if ((e->K = malloc(sizeof(struct double_matrix))) == NULL)
		return -1;
	if ((alloc_double_matrix(e->K, e->state_observable, e->state_size, NULL))
			< 0)
		return -1;
	if ((e->F = malloc(sizeof(struct double_matrix))) == NULL)
		return -1;
	if ((alloc_double_matrix(e->F, e->state_size, e->state_size, NULL) < 0))
		return -1;
	if ((e->H = malloc(sizeof(struct double_matrix))) == NULL)
		return -1;
	if ((alloc_double_matrix(e->H, e->state_size, e->state_observable, NULL))
			< 0)
		return -1;
	if ((e->Pp = malloc(sizeof(struct double_matrix))) == NULL)
		return -1;
	if ((alloc_double_matrix(e->Pp, e->state_size, e->state_size, NULL)) < 0)
		return -1;

	if ((e->PpHt = malloc(sizeof(struct double_matrix))) == NULL)
		return -1;
	if ((alloc_double_matrix(e->PpHt, e->state_observable, e->state_size, NULL))
			< 0)
		return -1;

	if ((e->HPp = malloc(sizeof(struct double_matrix))) == NULL)
		return -1;
	if ((alloc_double_matrix(e->HPp, e->state_size, e->state_observable, NULL))
			< 0)
		return -1;

	if ((e->HPpHt = malloc(sizeof(struct double_matrix))) == NULL)
		return -1;
	if ((alloc_double_matrix(e->HPpHt, e->state_observable, e->state_observable,
	NULL)) < 0)
		return -1;

	if ((e->HPpHt_inv = malloc(sizeof(struct double_matrix))) == NULL)
		return -1;
	if ((alloc_double_matrix(e->HPpHt_inv, e->state_observable,
			e->state_observable,
			NULL)) < 0)
		return -1;

	if ((e->FP = malloc(sizeof(struct double_matrix))) == NULL)
		return -1;
	if ((alloc_double_matrix(e->FP, e->state_size, e->state_size,
	NULL)) < 0)
		return -1;

	if ((e->KH = malloc(sizeof(struct double_matrix))) == NULL)
		return -1;
	if ((alloc_double_matrix(e->KH, e->state_size, e->state_size,
			NULL)) < 0)
		return -1;

	if ((e->ImKHP = malloc(sizeof(struct double_matrix))) == NULL)
		return -1;
	if ((alloc_double_matrix(e->ImKHP, e->state_size, e->state_size,
			NULL)) < 0)
		return -1;

	zeros_double_matrix(e->F);
	zeros_double_matrix(e->Q);
	zeros_double_matrix(e->P);
	zeros_double_matrix(e->R);
	zeros_double_matrix(e->H);
	zeros_double_matrix(e->Pp);
	zeros_double_matrix(e->K);
	memset(e->x, 0, e->state_size * sizeof(double));
	/*memset(e->fx, 0, e->state_size * sizeof(double));
	 memset(e->hx, 0, e->state_observable * sizeof(double));*/

	return 1;
}

void double_ekf_close(struct double_ekf_state * e) {
	free_double_matrix(e->F);
	free_double_matrix(e->H);
	free_double_matrix(e->P);
	free_double_matrix(e->Q);
	free_double_matrix(e->R);
	free_double_matrix(e->K);
	free_double_matrix(e->Pp);
	free_double_matrix(e->HPpHt);
	free_double_matrix(e->PpHt);
	free_double_matrix(e->HPp);
	free_double_matrix(e->HPpHt_inv);
	free_double_matrix(e->KH);
	free_double_matrix(e->ImKHP);
	free_double_matrix(e->FP);

	free(e->x);
	free(e->fx);
	free(e->hx);
	free(e);
}

//user is responsible for providing F and fx in the matrix to run prediction
int double_ekf_predict(struct double_ekf_state * e, void * private_data) {
	//copy_double_mat(&F, e->F, 0);

	//Update state using State transition Jacobian
	//x_{k} = Fx_{k-1}

	//double_mat_vect_product(e->F, e->x, e->fx);

	double_ekf_model(e, private_data);

	//Update state using command
	/* P_k = F_{k-1} P_{k-1} F^T_{k-1} + Q_{k-1} */
	/*Pp initialy contains a copy of P*/
	double_mat_product(e->F, e->Pp, e->FP, 0, 0); //tmp = FP
	double_mat_product(e->FP, e->F, e->Pp, 0, 1); //Pp = (FP)Ft
	double_mat_add(e->Pp, e->Q, e->Pp); //Pp = Pp + Q
	//need to validate with existing code
	return 0;
}

int double_ekf_update(struct double_ekf_state * e, double * z) {

	double * Kv, *v, *z_minus_hx, *K_z_minus_hx;

	Kv = malloc(e->state_size * sizeof(double));
	v = malloc(e->state_size * sizeof(double));
	z_minus_hx = malloc(e->state_observable * sizeof(double));
	K_z_minus_hx = malloc(e->state_size * sizeof(double));

	/* G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1} */
	double_mat_product(e->Pp, e->H, e->PpHt, 0, 1); //tmp = PpHt
	double_mat_product(e->H, e->Pp, e->HPp, 0, 0); //tmp2 = PpH
	double_mat_product(e->HPp, e->H, e->HPpHt, 0, 1); //(PpH)Ht
	double_mat_add(e->HPpHt, e->R, e->HPpHt);
	if (double_mat_inv(e->HPpHt, e->HPpHt_inv) < 0) {
		printf("Inversion problem ... \n");
		return -1; // problem in filter
	}
	//checked up to this point with Matlab
	double_mat_product(e->PpHt, e->HPpHt_inv, e->K, 0, 0); // gain update

	//double_mat_vect_product(e->H, e->fx, e->hx); //measurement prediction based on predicted state

	double_vec_sub(z, e->hx, v, e->state_observable); //measurement prediction error tmp5 = z - (Hx)
	double_mat_vect_product(e->K, v, Kv); //Kv
	double_vec_add(e->fx, Kv, e->x, e->state_size); //need size of vector

	//checked numerically with Scilab to this point
	/*P = (I - GH)P */
	/* P_k = Pk - GHPk */
	double_mat_product(e->K, e->H, e->KH, 0, 0 ); // result should be symmetric
	double_mat_eye_sub(e->KH, e->KH); //in place substraction
	double_mat_product(e->KH, e->Pp, e->P, 0, 1 );

	copy_double_matrix(e->P, e->Pp,0); //this will allow to run multiple prediction before any update
	//we should enforce symmetry of P to avoid numerical problems ...

	free(Kv);
	free(v);
	free(z_minus_hx);
	free(K_z_minus_hx);
	return 1;
}

