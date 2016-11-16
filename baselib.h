
#ifndef BASELIB_H
#define BASELIB_H


#ifdef _WIN32
# define DLL_EXPORT __declspec( dllexport )
#else
# define DLL_EXPORT
#endif


#include "vector.h"
#include "sfm.h"



//DLL_EXPORT v3_t dll_v3_new(double x, double y, double z); 

DLL_EXPORT void dll_refine_fmatrix_nonlinear_matches(int num_pts, v3_t *r_pts, v3_t *l_pts, 
	double *F0, double *Fout);

DLL_EXPORT double dll_fmatrix_compute_residual(double *F, v3_t r, v3_t l);

DLL_EXPORT int dll_estimate_fmatrix_ransac_matches(int num_pts, v3_t *a_pts, v3_t *b_pts, 
	int num_trials, double threshold, 
	double success_ratio, 
	int essential, double *F);


DLL_EXPORT int dll_compute_pose_ransac(int n, v2_t *r_pts, v2_t *l_pts, 
	double *K1, double *K2, 
	double ransac_threshold, int ransac_rounds, 
	double *R_out, double *t_out);


DLL_EXPORT int dll_compute_pose_ransac_pano(int n, v3_t *l_pts, v3_t *r_pts, double radius, 
	double ransac_threshold, int ransac_rounds, double* em,
	double *R_out, double *t_out);


DLL_EXPORT int dll_find_projection_3x4_ransac(int num_pts, v3_t *points, v2_t *projs, 
	double *P, int ransac_rounds, double ransac_threshold);


DLL_EXPORT double dll_fmatrix_compute_residual_pano(double *F, v3_t l, v3_t r, double radius);


DLL_EXPORT v3_t dll_triangulate(v2_t p, v2_t q, 
	double *R0, double *t0, 
	double *R1, double *t1, double *error);


DLL_EXPORT v3_t dll_triangulate_n(int num_points, 
	v2_t *p, double *R, double *t, double *error_out);

DLL_EXPORT int dll_iround(double x);

DLL_EXPORT v2_t dll_sfm_project_final(camera_params_t *params, v3_t pt, 
	int explicit_camera_centers, int undistort);


DLL_EXPORT double dll_kth_element_copy(int n, int k, double *arr);

DLL_EXPORT v3_t dll_triangulate_n_refine(v3_t pt, int num_points, 
	v2_t *p, double *R, double *t, double *error_out);

DLL_EXPORT double dll_median_copy(int n, double *arr);


/////////////////////////////// matrix /////////////////////////////////////////
DLL_EXPORT void dll_matrix_invert(int n, double *A, double *Ainv);

DLL_EXPORT void dll_matrix_product(int Am, int An, int Bm, int Bn, 
	const double *A, const double *B, double *R);

DLL_EXPORT void dll_matrix_scale(int m, int n, double *A, double s, double *R);

DLL_EXPORT void dll_matrix_product331(double *A, double *b, double *r);

DLL_EXPORT void dll_matrix_transpose(int m, int n, double *A, double *AT);

DLL_EXPORT void dll_matrix_sum(int Am, int An, int Bm, int Bn, 
	double *A, double *B, double *R);

DLL_EXPORT void dll_matrix_diff(int Am, int An, int Bm, int Bn, double *A, double *B, double *R);

DLL_EXPORT void dll_matrix_transpose_product(int Am, int An, int Bm, int Bn, double *A, double *B, double *R);

DLL_EXPORT void dll_dgerqf_driver(int m, int n, double *A, double *R, double *Q);
DLL_EXPORT void dll_dgelsy_driver(double *A, double *b, double *x, int m, int n, int nrhs);

DLL_EXPORT void dll_matrix_print(int m, int n, double *A);

DLL_EXPORT double dll_matrix_norm(int m, int n, double *A);

DLL_EXPORT void dll_matrix_ident(int n, double *A);

DLL_EXPORT v2_t dll_v2_scale(double s, v2_t v);

DLL_EXPORT int  dll_dgesvd_driver(int m, int n, double *A, double *U, double *S, double *VT);

////////////////////////////////////////////////////////////////////////////////
DLL_EXPORT double dll_align_horn(int n, v3_t *right_pts, v3_t *left_pts, 
	double *R, double *T, double *Tout, 
	double *scale, double *weight);

DLL_EXPORT void dll_align_homography(int num_pts, v3_t *r_pts, v3_t *l_pts, 
	double *Tout, int refine);


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
	double *Uout, double *Wout);


DLL_EXPORT void dll_sfm_project_rd(camera_params_t *init, double *K, double *k,
	double *R, double *dt, double *b, double *p,
	int undistort, int explicit_camera_centers);


DLL_EXPORT void dll_camera_refine(int num_points, v3_t *points, v2_t *projs, 
	camera_params_t *params, int adjust_focal,
	int estimate_distortion);

DLL_EXPORT v2_t dll_v2_new(double x, double y);

DLL_EXPORT v3_t dll_v3_new(double x, double y, double z);



#endif