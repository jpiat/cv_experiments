#include "float_ekf.h"

#include <string.h>
#include <stdio.h>

#include "float_mat_ops.h"


int float_ekf_init(struct float_ekf_state * e, unsigned int state_size,
		unsigned int observed_size) {
	e->state_size = state_size;
	e->state_observable = observed_size;

	if ((e->x = malloc(e->state_size * sizeof(float))) == NULL)
		return -1;
	if ((e->fx = malloc(e->state_size * sizeof(float))) == NULL)
		return -1;
	if ((e->hx = malloc(e->state_size * sizeof(float))) == NULL)
		return -1;

	if ((e->P = malloc(sizeof(struct float_matrix))) == NULL)
		return -1;
	if ((alloc_float_matrix(e->P, e->state_size, e->state_size, NULL)) < 0)
		return -1;
	if ((e->Q = malloc(sizeof(struct float_matrix))) == NULL)
		return -1;
	if ((alloc_float_matrix(e->Q, e->state_size, e->state_size, NULL)) < 0)
		return -1;
	if ((e->R = malloc(sizeof(struct float_matrix))) == NULL)
		return -1;
	if ((alloc_float_matrix(e->R, e->state_observable, e->state_observable,
	NULL)) < 0)
		return -1;
	if ((e->K = malloc(sizeof(struct float_matrix))) == NULL)
		return -1;
	if ((alloc_float_matrix(e->K, e->state_observable, e->state_size, NULL))
			< 0)
		return -1;
	if ((e->F = malloc(sizeof(struct float_matrix))) == NULL)
		return -1;
	if ((alloc_float_matrix(e->F, e->state_size, e->state_size, NULL) < 0))
		return -1;
	if ((e->H = malloc(sizeof(struct float_matrix))) == NULL)
		return -1;
	if ((alloc_float_matrix(e->H, e->state_size, e->state_observable, NULL))
			< 0)
		return -1;
	if ((e->Pp = malloc(sizeof(struct float_matrix))) == NULL)
		return -1;
	if ((alloc_float_matrix(e->Pp, e->state_size, e->state_size, NULL)) < 0)
		return -1;

	if ((e->PpHt = malloc(sizeof(struct float_matrix))) == NULL)
		return -1;
	if ((alloc_float_matrix(e->PpHt, e->state_observable, e->state_size, NULL))
			< 0)
		return -1;

	if ((e->HPp = malloc(sizeof(struct float_matrix))) == NULL)
		return -1;
	if ((alloc_float_matrix(e->HPp, e->state_size, e->state_observable, NULL))
			< 0)
		return -1;

	if ((e->HPpHt = malloc(sizeof(struct float_matrix))) == NULL)
		return -1;
	if ((alloc_float_matrix(e->HPpHt, e->state_observable, e->state_observable,
	NULL)) < 0)
		return -1;

	if ((e->HPpHt_inv = malloc(sizeof(struct float_matrix))) == NULL)
		return -1;
	if ((alloc_float_matrix(e->HPpHt_inv, e->state_observable,
			e->state_observable,
			NULL)) < 0)
		return -1;

	if ((e->FP = malloc(sizeof(struct float_matrix))) == NULL)
		return -1;
	if ((alloc_float_matrix(e->FP, e->state_size, e->state_size,
	NULL)) < 0)
		return -1;

	if ((e->KH = malloc(sizeof(struct float_matrix))) == NULL)
		return -1;
	if ((alloc_float_matrix(e->KH, e->state_size, e->state_size,
			NULL)) < 0)
		return -1;

	if ((e->ImKHP = malloc(sizeof(struct float_matrix))) == NULL)
		return -1;
	if ((alloc_float_matrix(e->ImKHP, e->state_size, e->state_size,
			NULL)) < 0)
		return -1;

	zeros_float_matrix(e->F);
	zeros_float_matrix(e->Q);
	zeros_float_matrix(e->P);
	zeros_float_matrix(e->R);
	zeros_float_matrix(e->H);
	zeros_float_matrix(e->Pp);
	zeros_float_matrix(e->K);
	memset(e->x, 0, e->state_size * sizeof(float));
	/*memset(e->fx, 0, e->state_size * sizeof(float));
	 memset(e->hx, 0, e->state_observable * sizeof(float));*/

	return 1;
}

void float_ekf_close(struct float_ekf_state * e) {
	free_float_matrix(e->F);
	free_float_matrix(e->H);
	free_float_matrix(e->P);
	free_float_matrix(e->Q);
	free_float_matrix(e->R);
	free_float_matrix(e->K);
	free_float_matrix(e->Pp);
	free_float_matrix(e->HPpHt);
	free_float_matrix(e->PpHt);
	free_float_matrix(e->HPp);
	free_float_matrix(e->HPpHt_inv);
	free_float_matrix(e->KH);
	free_float_matrix(e->ImKHP);
	free_float_matrix(e->FP);

	free(e->x);
	free(e->fx);
	free(e->hx);
	free(e);
}

//user is responsible for providing F and fx in the matrix to run prediction
int float_ekf_predict(struct float_ekf_state * e, void * private_data) {
	//copy_float_mat(&F, e->F, 0);

	//Update state using State transition Jacobian
	//x_{k} = Fx_{k-1}

	//float_mat_vect_product(e->F, e->x, e->fx);

	//Update state using command
	/* P_k = F_{k-1} P_{k-1} F^T_{k-1} + Q_{k-1} */
	/*Pp initialy contains a copy of P*/
	float_ekf_model(e, private_data);

	float_mat_product(e->F, e->Pp, e->FP, 0, 0); // FP = F*P
	float_mat_product(e->FP, e->F, e->Pp, 0, 1); //Pp = (FP)Ft
	float_mat_add(e->Pp, e->Q, e->Pp); //Pp = Pp + Q
	return 0;
}

int float_ekf_update(struct float_ekf_state * e, float * z) {

	float * Kv, *v, *z_minus_hx, *K_z_minus_hx;

	Kv = malloc(e->state_size * sizeof(float));
	v = malloc(e->state_size * sizeof(float));
	z_minus_hx = malloc(e->state_observable * sizeof(float));
	K_z_minus_hx = malloc(e->state_size * sizeof(float));

	/* G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1} */
	float_mat_product(e->Pp, e->H, e->PpHt, 0, 1); //PpHt = Pp *Ht
	float_mat_product(e->H, e->Pp, e->HPp, 0, 0); //HPp = H * Pp
	float_mat_product(e->HPp, e->H, e->HPpHt, 0, 1); //HPpHt = (HPp) * Ht
	float_mat_add(e->HPpHt, e->R, e->HPpHt); // HPpHt = HPpHt + R
	if (float_mat_inv(e->HPpHt, e->HPpHt_inv) < 0) { //HPpHt_inv = HPpHt{-1}
		printf("Inversion problem ... \n");
		return -1; // problem in filter
	}
	float_mat_product(e->PpHt, e->HPpHt_inv, e->K, 0, 0); // gain update

	//float_mat_vect_product(e->H, e->fx, e->hx); //measurement prediction based on predicted state

	float_vec_sub(z, e->hx, v, e->state_observable); //measurement prediction error tmp5 = z - (Hx)
	float_mat_vect_product(e->K, v, Kv); //Kv
	float_vec_add(e->fx, Kv, e->x, e->state_size); //need size of vector

	//checked numerically with Scilab to this point
	/*P = (I - GH)P */
	/* P_k = Pk - GHPk */


	float_mat_product(e->K, e->H, e->KH, 0, 0);
	float_mat_eye_sub(e->KH, e->KH); //in place substraction
	float_mat_product(e->KH, e->Pp, e->P, 0, 1); //Pp and P are symmetric ...



	copy_float_matrix(e->P, e->Pp,0); //this will allow to run multiple prediction before any update
	//we should enforce symmetry of P to avoid numerical problems ...

	free(Kv);
	free(v);
	free(z_minus_hx);
	free(K_z_minus_hx);
	return 1;
}

