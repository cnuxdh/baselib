

#include "baselib.h"



#include "5point.h"
#include "triangulate.h"
#include "fmatrix.h"
#include "util.h"
#include "sfm.h"
#include "qsort.h"
#include "matrix/matrix.h"
#include "horn.h"
#include "homography.h"
#include "fmatrix.h"




DLL_EXPORT void dll_refine_fmatrix_nonlinear_matches(int num_pts, v3_t *r_pts, v3_t *l_pts, 
	double *F0, double *Fout)
{
	refine_fmatrix_nonlinear_matches(num_pts, r_pts, l_pts, F0, Fout);
}

DLL_EXPORT double dll_fmatrix_compute_residual(double *F, v3_t r, v3_t l)
{
	return fmatrix_compute_residual(F, r, l);
}

DLL_EXPORT int dll_estimate_fmatrix_ransac_matches(int num_pts, v3_t *a_pts, v3_t *b_pts, 
	int num_trials, double threshold, 
	double success_ratio, 
	int essential, double *F)
{
	return estimate_fmatrix_ransac_matches(num_pts, a_pts, b_pts, num_trials, 
		threshold, success_ratio, essential, F);
}



int dll_compute_pose_ransac(int n, v2_t *r_pts, v2_t *l_pts, 
	double *K1, double *K2, 
	double ransac_threshold, int ransac_rounds, 
	double *R_out, double *t_out)
{
	return compute_pose_ransac(n, r_pts, l_pts, K1, K2, ransac_threshold, ransac_rounds, R_out, t_out);
}

int dll_compute_pose_ransac_pano(int n, v3_t *l_pts, v3_t *r_pts, double radius, 
	double ransac_threshold, int ransac_rounds, double* em,
	double *R_out, double *t_out)
{
	return compute_pose_ransac_pano(n, l_pts, r_pts, radius, ransac_threshold, ransac_rounds, em, R_out, t_out);
}


int dll_find_projection_3x4_ransac(int num_pts, v3_t *points, v2_t *projs, 
	double *P, int ransac_rounds, double ransac_threshold)
{
	return find_projection_3x4_ransac(num_pts, points, projs, P, ransac_rounds, ransac_threshold);
}

double dll_fmatrix_compute_residual_pano(double *F, v3_t l, v3_t r, double radius)
{
	return fmatrix_compute_residual_pano(F, l, r, radius);
}

v3_t dll_triangulate(v2_t p, v2_t q, 
	double *R0, double *t0, 
	double *R1, double *t1, double *error)
{
	return triangulate(p, q, R0, t0, R1, t1, error);
}

v3_t dll_triangulate_n(int num_points, 
	v2_t *p, double *R, double *t, double *error_out)
{
	return triangulate_n(num_points, p, R, t, error_out);
}

int dll_iround(double x)
{
	return iround(x);
}

v2_t dll_sfm_project_final(camera_params_t *params, v3_t pt, 
	int explicit_camera_centers, int undistort)
{
	return sfm_project_final(params, pt, explicit_camera_centers, undistort);
}

double dll_kth_element_copy(int n, int k, double *arr)
{
	return kth_element_copy(n, k, arr);
}

v3_t dll_triangulate_n_refine(v3_t pt, int num_points, 
	v2_t *p, double *R, double *t, double *error_out)
{
	return triangulate_n_refine(pt, num_points, p, R, t, error_out);
}

double dll_median_copy(int n, double *arr)
{
	return median_copy(n, arr);
}



/////////////////////////////// matrix /////////////////////////////////////////
DLL_EXPORT void dll_matrix_invert(int n, double *A, double *Ainv)
{
	matrix_invert(n, A, Ainv);
}

DLL_EXPORT void dll_matrix_product(int Am, int An, int Bm, int Bn, 
	const double *A, const double *B, double *R)
{
	matrix_product(Am, An, Bm, Bn, A, B, R);
}

DLL_EXPORT void dll_matrix_scale(int m, int n, double *A, double s, double *R)
{
	matrix_scale(m, n, A, s, R);
}

DLL_EXPORT void dll_matrix_product331(double *A, double *b, double *r)
{
	matrix_product331(A, b, r);
}

DLL_EXPORT void dll_matrix_transpose(int m, int n, double *A, double *AT)
{
	matrix_transpose(m, n, A, AT);
}

DLL_EXPORT void dll_matrix_sum(int Am, int An, int Bm, int Bn, 
	double *A, double *B, double *R)
{
	matrix_sum(Am, An, Bm, Bn, A, B, R);
}

DLL_EXPORT void dll_matrix_diff(int Am, int An, int Bm, int Bn, double *A, double *B, double *R)
{
	matrix_diff(Am,An,Bm,Bn,A,B,R);
}

DLL_EXPORT void dll_matrix_transpose_product(int Am, int An, int Bm, int Bn, double *A, double *B, double *R)
{
	matrix_transpose_product(Am, An, Bm,Bn, A, B, R);
}

DLL_EXPORT void dll_dgerqf_driver(int m, int n, double *A, double *R, double *Q)
{
	dgerqf_driver(m,n,A,R,Q);
}
DLL_EXPORT void dll_dgelsy_driver(double *A, double *b, double *x, int m, int n, int nrhs)
{
	dgelsy_driver(A,b,x,m,n,nrhs);
}

DLL_EXPORT void dll_matrix_print(int m, int n, double *A)
{
	matrix_print(m,n,A);
}

DLL_EXPORT double dll_matrix_norm(int m, int n, double *A)
{
	return matrix_norm(m,n,A);
}

DLL_EXPORT void dll_matrix_ident(int n, double *A)
{
	matrix_ident(n,A);
}

DLL_EXPORT v2_t dll_v2_scale(double s, v2_t v)
{
	return v2_scale(s, v);
}

int  dll_dgesvd_driver(int m, int n, double *A, double *U, double *S, double *VT)
{
	return dgesvd_driver(m,n,A,U,S,VT);
}

////////////////////////////////////////////////////////////////////////////////
DLL_EXPORT double dll_align_horn(int n, v3_t *right_pts, v3_t *left_pts, 
	double *R, double *T, double *Tout, 
	double *scale, double *weight)
{
	return align_horn(n, right_pts, left_pts, R, T, Tout, scale, weight);
}

DLL_EXPORT void dll_align_homography(int num_pts, v3_t *r_pts, v3_t *l_pts, 
	double *Tout, int refine)
{
	align_homography(num_pts, r_pts, l_pts, Tout, refine);
}


////////////////////////////////////////////////////////////////////////////////
DLL_EXPORT void dll_run_sfm(int num_pts, int num_cameras, int ncons,
	char *vmask,
	double *projections,
	int est_focal_length,
	int const_focal_length,
	int undistort,
	int explicit_camera_centers,
	camera_params_t *init_camera_params,
	v3_t *init_pts, 
	int use_constraints, 
	int use_point_constraints,
	v3_t *points_constraints,
	double point_constraint_weight,
	int fix_points,
	int optimize_for_fisheye, 
	double eps2,
	double *Vout,
	double *Sout,
	double *Uout, double *Wout)
{
	run_sfm(num_pts, num_cameras, ncons, vmask, projections, est_focal_length,
		const_focal_length,  undistort, explicit_camera_centers, 
		init_camera_params, init_pts, use_constraints, use_point_constraints, 
		points_constraints, point_constraint_weight, fix_points, optimize_for_fisheye,
		eps2, Vout, Sout, Uout,Wout);
}


DLL_EXPORT void dll_sfm_project_rd(camera_params_t *init, double *K, double *k,
	double *R, double *dt, double *b, double *p,
	int undistort, int explicit_camera_centers)
{
	sfm_project_rd(init, K, K,R, dt, b, p, undistort, explicit_camera_centers);
}

void dll_camera_refine(int num_points, v3_t *points, v2_t *projs, 
	camera_params_t *params, int adjust_focal,
	int estimate_distortion)
{
	camera_refine(num_points, points, projs, params, adjust_focal, estimate_distortion);
}

v2_t dll_v2_new(double x, double y)
{
	return v2_new(x,y);
}

v3_t dll_v3_new(double x, double y, double z)
{
	return v3_new(x,y,z);
}
