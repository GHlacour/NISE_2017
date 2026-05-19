#ifndef _MCFRET4_ /* ignore */
#define _MCFRET4_

typedef struct {
    float *rho0_D;
    float *J;
    float *integrated_response_tw;
    float *U_re_t1_array;
    float *U_im_t1_array;
    float *U_re_tw_array;
    float *U_im_tw_array;
    float *U_re_t2_array;
    float *U_im_t2_array;
    float *diagram_1;
    float *diagram_2;
    float *diagram_3;
    float *diagram_4;

    int N_A;
    int N_D;
    int N_t1;
    int N_tw;
    int N_t2;
    int times_N2;
    int s_D;
    int s_A;
    int N_segments;
    int largest_segment_size;
    
    t_non *non;
} fourth_order_params;

typedef struct {
    // the four complex intermediate products that can be reused in the triple time loop
    float *intermediate_product_1_re; float *intermediate_product_1_im;
    float *intermediate_product_2_re; float *intermediate_product_2_im;
    float *intermediate_product_3_re; float *intermediate_product_3_im;
    float *intermediate_product_4_re; float *intermediate_product_4_im;

    // t1 dependent products with the density matrix
    float *UDh_rho_J_UA_re_t1; float *UDh_rho_J_UA_im_t1;
    float *UAh_Jh_rho_UD_re_t1; float *UAh_Jh_rho_UD_im_t1;

    // tw dependent matrix products
    float *Jh_UDh_tw_re;  float *Jh_UDh_tw_im;
    float *Jh_UD_tw_re;   float *Jh_UD_tw_im;
    float *UD_tw_re;      float *UD_tw_im;
    float *UA_Jh_tw_re;   float *UA_Jh_tw_im;
    float *UAh_Jh_tw_re;  float *UAh_Jh_tw_im;
    float *J_UAh_Jh_tw_re; float *J_UAh_Jh_tw_im;
    float *UA_tw_re;      float *UA_tw_im;

    float *J_zeros;
} fourth_order_workspace;

void compute_all_traces_4th_order(float *rho_0,float *J_full,t_non *non);
void compute_UJJU(float *UJJU_re, float *JJ, float *U_re, float *U_im, int N_i,int sj);
void compute_JrhoJ(float *Jij_rho_jj_Jji, float* Jij, float *rho_0_sj, int N_i, int N_j, int sj);
void compute_rhoJJ(float *rho_ii_JijJji, float *JijJji, float* Jij, float *rho_0_si, int N_i, int N_j, int sj);

void write_propagator_to_big_array(float *big_array, float *propagator, int sample_length, int si, int N_site_si, int largest_segment_size, int ti);
void read_propagator_from_big_array(float *big_array, float *propagator, int sample_length, int si, int N_site_si, int largest_segment_size, int ti);

void compute_UDh_rho_J_UA_t1(float *UDh_rho_J_UA_re_t1,
                             float *UDh_rho_J_UA_im_t1,
                             float *UAh_Jh_rho_UD_re_t1,
                             float *UAh_Jh_rho_UD_im_t1,
                             fourth_order_params *p);
void compute_rate_from_4th_response(float *responses_4th_tw, float *rate_matrix_4th_I, float *rate_matrix_4th_II, float *rate_matrix_2nd, int N_segments, int N_tw, t_non *non, float prefactor, int samples);
void full_4th_order_main(float *rho_0,float *J_full,t_non *non);
void fourth_order_response_1_sample(fourth_order_params *p);
void add_diagrams(float *diagram_1, float *diagram_2, float *diagram_3, float *diagram_4, float *diagram_1_s01, float *diagram_2_s01, float *diagram_3_s01, float *diagram_4_s01, int N_t1, int N_t2);
void compute_4_intermediate_products(fourth_order_workspace *ws, fourth_order_params *p);



#endif /* _MCFRET4_ */