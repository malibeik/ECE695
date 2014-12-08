/* Include files */

#include <stddef.h>
#include "blas.h"
#include "gen_rec_sfun.h"
#include "c2_gen_rec.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "gen_rec_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static real_T _sfTime_;
static const char * c2_debug_family_names[40] = { "lamq1", "lamq2", "lamq3",
  "lamd1", "lamd2", "lamd3", "lamdf", "imq", "imd", "iq1", "iq2", "iq3", "id1",
  "id2", "id3", "plamq1", "plamq2", "plamq3", "plamd1", "plamd2", "plamd3",
  "plamdf", "lamMq", "lamMd", "wr", "lamqS", "lamdS", "nargin", "nargout", "vfd",
  "wrm", "iqs", "ids", "Iamr", "P", "vqs", "vds", "Te", "ifd", "pIamr" };

/* Function Declarations */
static void initialize_c2_gen_rec(SFc2_gen_recInstanceStruct *chartInstance);
static void initialize_params_c2_gen_rec(SFc2_gen_recInstanceStruct
  *chartInstance);
static void enable_c2_gen_rec(SFc2_gen_recInstanceStruct *chartInstance);
static void disable_c2_gen_rec(SFc2_gen_recInstanceStruct *chartInstance);
static void c2_update_debugger_state_c2_gen_rec(SFc2_gen_recInstanceStruct
  *chartInstance);
static const mxArray *get_sim_state_c2_gen_rec(SFc2_gen_recInstanceStruct
  *chartInstance);
static void set_sim_state_c2_gen_rec(SFc2_gen_recInstanceStruct *chartInstance,
  const mxArray *c2_st);
static void finalize_c2_gen_rec(SFc2_gen_recInstanceStruct *chartInstance);
static void sf_gateway_c2_gen_rec(SFc2_gen_recInstanceStruct *chartInstance);
static void c2_chartstep_c2_gen_rec(SFc2_gen_recInstanceStruct *chartInstance);
static void initSimStructsc2_gen_rec(SFc2_gen_recInstanceStruct *chartInstance);
static void init_script_number_translation(uint32_T c2_machineNumber, uint32_T
  c2_chartNumber, uint32_T c2_instanceNumber);
static const mxArray *c2_sf_marshallOut(void *chartInstanceVoid, void *c2_inData);
static void c2_emlrt_marshallIn(SFc2_gen_recInstanceStruct *chartInstance, const
  mxArray *c2_pIamr, const char_T *c2_identifier, real_T c2_y[7]);
static void c2_b_emlrt_marshallIn(SFc2_gen_recInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[7]);
static void c2_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static const mxArray *c2_b_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static real_T c2_c_emlrt_marshallIn(SFc2_gen_recInstanceStruct *chartInstance,
  const mxArray *c2_ifd, const char_T *c2_identifier);
static real_T c2_d_emlrt_marshallIn(SFc2_gen_recInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId);
static void c2_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static const mxArray *c2_c_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static void c2_e_emlrt_marshallIn(SFc2_gen_recInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  c2_struct_U3Zfiu1q58uXf1A0Qkd5GF *c2_y);
static void c2_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static void c2_info_helper(const mxArray **c2_info);
static const mxArray *c2_emlrt_marshallOut(const char * c2_u);
static const mxArray *c2_b_emlrt_marshallOut(const uint32_T c2_u);
static void c2_rdivide(SFc2_gen_recInstanceStruct *chartInstance, real_T c2_x[3],
  real_T c2_y[3], real_T c2_z[3]);
static real_T c2_sum(SFc2_gen_recInstanceStruct *chartInstance, real_T c2_x[3]);
static void c2_b_rdivide(SFc2_gen_recInstanceStruct *chartInstance, real_T c2_x,
  real_T c2_y[3], real_T c2_z[3]);
static const mxArray *c2_d_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static int32_T c2_f_emlrt_marshallIn(SFc2_gen_recInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId);
static void c2_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static uint8_T c2_g_emlrt_marshallIn(SFc2_gen_recInstanceStruct *chartInstance,
  const mxArray *c2_b_is_active_c2_gen_rec, const char_T *c2_identifier);
static uint8_T c2_h_emlrt_marshallIn(SFc2_gen_recInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId);
static void init_dsm_address_info(SFc2_gen_recInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c2_gen_rec(SFc2_gen_recInstanceStruct *chartInstance)
{
  chartInstance->c2_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c2_is_active_c2_gen_rec = 0U;
}

static void initialize_params_c2_gen_rec(SFc2_gen_recInstanceStruct
  *chartInstance)
{
  const mxArray *c2_m0 = NULL;
  const mxArray *c2_mxField;
  c2_struct_U3Zfiu1q58uXf1A0Qkd5GF c2_r0;
  c2_m0 = sf_mex_get_sfun_param(chartInstance->S, 0, 1);
  c2_mxField = sf_mex_getfield(c2_m0, "Ll", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c2_mxField), &c2_r0.Ll, 1, 0, 0U, 0, 0U, 0);
  c2_mxField = sf_mex_getfield(c2_m0, "rl", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c2_mxField), &c2_r0.rl, 1, 0, 0U, 0, 0U, 0);
  c2_mxField = sf_mex_getfield(c2_m0, "rs", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c2_mxField), &c2_r0.rs, 1, 0, 0U, 0, 0U, 0);
  c2_mxField = sf_mex_getfield(c2_m0, "Lls", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c2_mxField), &c2_r0.Lls, 1, 0, 0U, 0, 0U,
                      0);
  c2_mxField = sf_mex_getfield(c2_m0, "Lmq", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c2_mxField), &c2_r0.Lmq, 1, 0, 0U, 0, 0U,
                      0);
  c2_mxField = sf_mex_getfield(c2_m0, "Lmd", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c2_mxField), &c2_r0.Lmd, 1, 0, 0U, 0, 0U,
                      0);
  c2_mxField = sf_mex_getfield(c2_m0, "rd1", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c2_mxField), &c2_r0.rd1, 1, 0, 0U, 0, 0U,
                      0);
  c2_mxField = sf_mex_getfield(c2_m0, "rd2", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c2_mxField), &c2_r0.rd2, 1, 0, 0U, 0, 0U,
                      0);
  c2_mxField = sf_mex_getfield(c2_m0, "rd3", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c2_mxField), &c2_r0.rd3, 1, 0, 0U, 0, 0U,
                      0);
  c2_mxField = sf_mex_getfield(c2_m0, "rfd", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c2_mxField), &c2_r0.rfd, 1, 0, 0U, 0, 0U,
                      0);
  c2_mxField = sf_mex_getfield(c2_m0, "Lld1", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c2_mxField), &c2_r0.Lld1, 1, 0, 0U, 0, 0U,
                      0);
  c2_mxField = sf_mex_getfield(c2_m0, "Lld2", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c2_mxField), &c2_r0.Lld2, 1, 0, 0U, 0, 0U,
                      0);
  c2_mxField = sf_mex_getfield(c2_m0, "Lld3", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c2_mxField), &c2_r0.Lld3, 1, 0, 0U, 0, 0U,
                      0);
  c2_mxField = sf_mex_getfield(c2_m0, "Llfd", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c2_mxField), &c2_r0.Llfd, 1, 0, 0U, 0, 0U,
                      0);
  c2_mxField = sf_mex_getfield(c2_m0, "rq1", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c2_mxField), &c2_r0.rq1, 1, 0, 0U, 0, 0U,
                      0);
  c2_mxField = sf_mex_getfield(c2_m0, "rq2", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c2_mxField), &c2_r0.rq2, 1, 0, 0U, 0, 0U,
                      0);
  c2_mxField = sf_mex_getfield(c2_m0, "rq3", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c2_mxField), &c2_r0.rq3, 1, 0, 0U, 0, 0U,
                      0);
  c2_mxField = sf_mex_getfield(c2_m0, "pole", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c2_mxField), &c2_r0.pole, 1, 0, 0U, 0, 0U,
                      0);
  c2_mxField = sf_mex_getfield(c2_m0, "Llq1", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c2_mxField), &c2_r0.Llq1, 1, 0, 0U, 0, 0U,
                      0);
  c2_mxField = sf_mex_getfield(c2_m0, "Llq2", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c2_mxField), &c2_r0.Llq2, 1, 0, 0U, 0, 0U,
                      0);
  c2_mxField = sf_mex_getfield(c2_m0, "Llq3", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c2_mxField), &c2_r0.Llq3, 1, 0, 0U, 0, 0U,
                      0);
  sf_mex_destroy(&c2_m0);
  chartInstance->c2_P = c2_r0;
}

static void enable_c2_gen_rec(SFc2_gen_recInstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c2_gen_rec(SFc2_gen_recInstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c2_update_debugger_state_c2_gen_rec(SFc2_gen_recInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c2_gen_rec(SFc2_gen_recInstanceStruct
  *chartInstance)
{
  const mxArray *c2_st;
  const mxArray *c2_y = NULL;
  real_T c2_hoistedGlobal;
  real_T c2_u;
  const mxArray *c2_b_y = NULL;
  real_T c2_b_hoistedGlobal;
  real_T c2_b_u;
  const mxArray *c2_c_y = NULL;
  int32_T c2_i0;
  real_T c2_c_u[7];
  const mxArray *c2_d_y = NULL;
  real_T c2_c_hoistedGlobal;
  real_T c2_d_u;
  const mxArray *c2_e_y = NULL;
  real_T c2_d_hoistedGlobal;
  real_T c2_e_u;
  const mxArray *c2_f_y = NULL;
  uint8_T c2_e_hoistedGlobal;
  uint8_T c2_f_u;
  const mxArray *c2_g_y = NULL;
  real_T *c2_Te;
  real_T *c2_ifd;
  real_T *c2_vds;
  real_T *c2_vqs;
  real_T (*c2_pIamr)[7];
  c2_pIamr = (real_T (*)[7])ssGetOutputPortSignal(chartInstance->S, 5);
  c2_ifd = (real_T *)ssGetOutputPortSignal(chartInstance->S, 4);
  c2_Te = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c2_vds = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c2_vqs = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c2_st = NULL;
  c2_st = NULL;
  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_createcellmatrix(6, 1), false);
  c2_hoistedGlobal = *c2_Te;
  c2_u = c2_hoistedGlobal;
  c2_b_y = NULL;
  sf_mex_assign(&c2_b_y, sf_mex_create("y", &c2_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c2_y, 0, c2_b_y);
  c2_b_hoistedGlobal = *c2_ifd;
  c2_b_u = c2_b_hoistedGlobal;
  c2_c_y = NULL;
  sf_mex_assign(&c2_c_y, sf_mex_create("y", &c2_b_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c2_y, 1, c2_c_y);
  for (c2_i0 = 0; c2_i0 < 7; c2_i0++) {
    c2_c_u[c2_i0] = (*c2_pIamr)[c2_i0];
  }

  c2_d_y = NULL;
  sf_mex_assign(&c2_d_y, sf_mex_create("y", c2_c_u, 0, 0U, 1U, 0U, 1, 7), false);
  sf_mex_setcell(c2_y, 2, c2_d_y);
  c2_c_hoistedGlobal = *c2_vds;
  c2_d_u = c2_c_hoistedGlobal;
  c2_e_y = NULL;
  sf_mex_assign(&c2_e_y, sf_mex_create("y", &c2_d_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c2_y, 3, c2_e_y);
  c2_d_hoistedGlobal = *c2_vqs;
  c2_e_u = c2_d_hoistedGlobal;
  c2_f_y = NULL;
  sf_mex_assign(&c2_f_y, sf_mex_create("y", &c2_e_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c2_y, 4, c2_f_y);
  c2_e_hoistedGlobal = chartInstance->c2_is_active_c2_gen_rec;
  c2_f_u = c2_e_hoistedGlobal;
  c2_g_y = NULL;
  sf_mex_assign(&c2_g_y, sf_mex_create("y", &c2_f_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c2_y, 5, c2_g_y);
  sf_mex_assign(&c2_st, c2_y, false);
  return c2_st;
}

static void set_sim_state_c2_gen_rec(SFc2_gen_recInstanceStruct *chartInstance,
  const mxArray *c2_st)
{
  const mxArray *c2_u;
  real_T c2_dv0[7];
  int32_T c2_i1;
  real_T *c2_Te;
  real_T *c2_ifd;
  real_T *c2_vds;
  real_T *c2_vqs;
  real_T (*c2_pIamr)[7];
  c2_pIamr = (real_T (*)[7])ssGetOutputPortSignal(chartInstance->S, 5);
  c2_ifd = (real_T *)ssGetOutputPortSignal(chartInstance->S, 4);
  c2_Te = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c2_vds = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c2_vqs = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c2_doneDoubleBufferReInit = true;
  c2_u = sf_mex_dup(c2_st);
  *c2_Te = c2_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u,
    0)), "Te");
  *c2_ifd = c2_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u,
    1)), "ifd");
  c2_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u, 2)),
                      "pIamr", c2_dv0);
  for (c2_i1 = 0; c2_i1 < 7; c2_i1++) {
    (*c2_pIamr)[c2_i1] = c2_dv0[c2_i1];
  }

  *c2_vds = c2_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u,
    3)), "vds");
  *c2_vqs = c2_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u,
    4)), "vqs");
  chartInstance->c2_is_active_c2_gen_rec = c2_g_emlrt_marshallIn(chartInstance,
    sf_mex_dup(sf_mex_getcell(c2_u, 5)), "is_active_c2_gen_rec");
  sf_mex_destroy(&c2_u);
  c2_update_debugger_state_c2_gen_rec(chartInstance);
  sf_mex_destroy(&c2_st);
}

static void finalize_c2_gen_rec(SFc2_gen_recInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c2_gen_rec(SFc2_gen_recInstanceStruct *chartInstance)
{
  int32_T c2_i2;
  int32_T c2_i3;
  real_T *c2_vfd;
  real_T *c2_vqs;
  real_T *c2_wrm;
  real_T *c2_iqs;
  real_T *c2_ids;
  real_T *c2_vds;
  real_T *c2_Te;
  real_T *c2_ifd;
  real_T (*c2_pIamr)[7];
  real_T (*c2_Iamr)[7];
  c2_pIamr = (real_T (*)[7])ssGetOutputPortSignal(chartInstance->S, 5);
  c2_ifd = (real_T *)ssGetOutputPortSignal(chartInstance->S, 4);
  c2_Te = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c2_vds = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c2_Iamr = (real_T (*)[7])ssGetInputPortSignal(chartInstance->S, 4);
  c2_ids = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
  c2_iqs = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c2_wrm = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c2_vqs = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c2_vfd = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 1U, chartInstance->c2_sfEvent);
  _SFD_DATA_RANGE_CHECK(*c2_vfd, 0U);
  chartInstance->c2_sfEvent = CALL_EVENT;
  c2_chartstep_c2_gen_rec(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_gen_recMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  _SFD_DATA_RANGE_CHECK(*c2_vqs, 1U);
  _SFD_DATA_RANGE_CHECK(*c2_wrm, 2U);
  _SFD_DATA_RANGE_CHECK(*c2_iqs, 3U);
  _SFD_DATA_RANGE_CHECK(*c2_ids, 4U);
  for (c2_i2 = 0; c2_i2 < 7; c2_i2++) {
    _SFD_DATA_RANGE_CHECK((*c2_Iamr)[c2_i2], 5U);
  }

  _SFD_DATA_RANGE_CHECK(*c2_vds, 6U);
  _SFD_DATA_RANGE_CHECK(*c2_Te, 7U);
  _SFD_DATA_RANGE_CHECK(*c2_ifd, 8U);
  for (c2_i3 = 0; c2_i3 < 7; c2_i3++) {
    _SFD_DATA_RANGE_CHECK((*c2_pIamr)[c2_i3], 9U);
  }
}

static void c2_chartstep_c2_gen_rec(SFc2_gen_recInstanceStruct *chartInstance)
{
  real_T c2_hoistedGlobal;
  real_T c2_b_hoistedGlobal;
  real_T c2_c_hoistedGlobal;
  real_T c2_d_hoistedGlobal;
  real_T c2_vfd;
  real_T c2_wrm;
  real_T c2_iqs;
  real_T c2_ids;
  int32_T c2_i4;
  real_T c2_Iamr[7];
  c2_struct_U3Zfiu1q58uXf1A0Qkd5GF c2_b_P;
  uint32_T c2_debug_family_var_map[40];
  real_T c2_lamq1;
  real_T c2_lamq2;
  real_T c2_lamq3;
  real_T c2_lamd1;
  real_T c2_lamd2;
  real_T c2_lamd3;
  real_T c2_lamdf;
  real_T c2_imq;
  real_T c2_imd;
  real_T c2_iq1;
  real_T c2_iq2;
  real_T c2_iq3;
  real_T c2_id1;
  real_T c2_id2;
  real_T c2_id3;
  real_T c2_plamq1;
  real_T c2_plamq2;
  real_T c2_plamq3;
  real_T c2_plamd1;
  real_T c2_plamd2;
  real_T c2_plamd3;
  real_T c2_plamdf;
  real_T c2_lamMq;
  real_T c2_lamMd;
  real_T c2_wr;
  real_T c2_lamqS;
  real_T c2_lamdS;
  real_T c2_nargin = 6.0;
  real_T c2_nargout = 5.0;
  real_T c2_vqs;
  real_T c2_vds;
  real_T c2_Te;
  real_T c2_ifd;
  real_T c2_pIamr[7];
  real_T c2_b_lamq1[3];
  real_T c2_c_P[3];
  real_T c2_dv1[3];
  int32_T c2_i5;
  real_T c2_dv2[3];
  real_T c2_A;
  real_T c2_d_P[3];
  int32_T c2_i6;
  real_T c2_dv3[3];
  real_T c2_B;
  real_T c2_x;
  real_T c2_y;
  real_T c2_b_x;
  real_T c2_b_y;
  real_T c2_c_x;
  real_T c2_c_y;
  real_T c2_b_A;
  real_T c2_b_B;
  real_T c2_d_x;
  real_T c2_d_y;
  real_T c2_e_x;
  real_T c2_e_y;
  real_T c2_f_x;
  real_T c2_f_y;
  real_T c2_g_y;
  real_T c2_c_B;
  real_T c2_h_y;
  real_T c2_i_y;
  real_T c2_j_y;
  real_T c2_k_y;
  real_T c2_l_y[3];
  int32_T c2_i7;
  real_T c2_b_lamd1[3];
  real_T c2_e_P[3];
  int32_T c2_i8;
  real_T c2_dv4[3];
  real_T c2_c_A;
  int32_T c2_i9;
  real_T c2_m_y[3];
  real_T c2_d_B;
  real_T c2_g_x;
  real_T c2_n_y;
  real_T c2_h_x;
  real_T c2_o_y;
  real_T c2_i_x;
  real_T c2_p_y;
  real_T c2_d_A;
  real_T c2_e_B;
  real_T c2_j_x;
  real_T c2_q_y;
  real_T c2_k_x;
  real_T c2_r_y;
  real_T c2_l_x;
  real_T c2_s_y;
  real_T c2_e_A;
  real_T c2_f_B;
  real_T c2_m_x;
  real_T c2_t_y;
  real_T c2_n_x;
  real_T c2_u_y;
  real_T c2_o_x;
  real_T c2_v_y;
  real_T c2_f_A;
  real_T c2_g_B;
  real_T c2_p_x;
  real_T c2_w_y;
  real_T c2_q_x;
  real_T c2_x_y;
  real_T c2_r_x;
  real_T c2_y_y;
  real_T c2_g_A;
  real_T c2_h_B;
  real_T c2_s_x;
  real_T c2_ab_y;
  real_T c2_t_x;
  real_T c2_bb_y;
  real_T c2_u_x;
  real_T c2_cb_y;
  real_T c2_h_A;
  real_T c2_i_B;
  real_T c2_v_x;
  real_T c2_db_y;
  real_T c2_w_x;
  real_T c2_eb_y;
  real_T c2_x_x;
  real_T c2_fb_y;
  real_T c2_i_A;
  real_T c2_j_B;
  real_T c2_y_x;
  real_T c2_gb_y;
  real_T c2_ab_x;
  real_T c2_hb_y;
  real_T c2_bb_x;
  real_T c2_ib_y;
  real_T c2_j_A;
  real_T c2_k_B;
  real_T c2_cb_x;
  real_T c2_jb_y;
  real_T c2_db_x;
  real_T c2_kb_y;
  real_T c2_eb_x;
  real_T c2_lb_y;
  real_T c2_b_plamq1[7];
  int32_T c2_i10;
  real_T c2_k_A;
  real_T c2_fb_x;
  real_T c2_gb_x;
  real_T c2_hb_x;
  real_T c2_mb_y;
  real_T c2_l_A;
  real_T c2_ib_x;
  real_T c2_jb_x;
  real_T c2_kb_x;
  real_T c2_nb_y;
  int32_T c2_i11;
  real_T *c2_b_ifd;
  real_T *c2_b_Te;
  real_T *c2_b_vds;
  real_T *c2_b_vqs;
  real_T *c2_b_ids;
  real_T *c2_b_iqs;
  real_T *c2_b_wrm;
  real_T *c2_b_vfd;
  real_T (*c2_b_pIamr)[7];
  real_T (*c2_b_Iamr)[7];
  c2_b_pIamr = (real_T (*)[7])ssGetOutputPortSignal(chartInstance->S, 5);
  c2_b_ifd = (real_T *)ssGetOutputPortSignal(chartInstance->S, 4);
  c2_b_Te = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c2_b_vds = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c2_b_Iamr = (real_T (*)[7])ssGetInputPortSignal(chartInstance->S, 4);
  c2_b_ids = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
  c2_b_iqs = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c2_b_wrm = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c2_b_vqs = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c2_b_vfd = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 1U, chartInstance->c2_sfEvent);
  c2_hoistedGlobal = *c2_b_vfd;
  c2_b_hoistedGlobal = *c2_b_wrm;
  c2_c_hoistedGlobal = *c2_b_iqs;
  c2_d_hoistedGlobal = *c2_b_ids;
  c2_vfd = c2_hoistedGlobal;
  c2_wrm = c2_b_hoistedGlobal;
  c2_iqs = c2_c_hoistedGlobal;
  c2_ids = c2_d_hoistedGlobal;
  for (c2_i4 = 0; c2_i4 < 7; c2_i4++) {
    c2_Iamr[c2_i4] = (*c2_b_Iamr)[c2_i4];
  }

  c2_b_P = chartInstance->c2_P;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 40U, 40U, c2_debug_family_names,
    c2_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_lamq1, 0U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_lamq2, 1U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_lamq3, 2U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_lamd1, 3U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_lamd2, 4U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_lamd3, 5U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_lamdf, 6U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_imq, 7U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_imd, 8U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_iq1, 9U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_iq2, 10U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_iq3, 11U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_id1, 12U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_id2, 13U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_id3, 14U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_plamq1, 15U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_plamq2, 16U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_plamq3, 17U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_plamd1, 18U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_plamd2, 19U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_plamd3, 20U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_plamdf, 21U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_lamMq, 22U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_lamMd, 23U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_wr, 24U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_lamqS, 25U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_lamdS, 26U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargin, 27U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargout, 28U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_vfd, 29U, c2_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_wrm, 30U, c2_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_iqs, 31U, c2_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_ids, 32U, c2_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c2_Iamr, 33U, c2_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_b_P, 34U, c2_c_sf_marshallOut,
    c2_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_vqs, 35U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_vds, 36U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_Te, 37U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_ifd, 38U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_pIamr, 39U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 44);
  c2_lamq1 = c2_Iamr[0];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 45);
  c2_lamq2 = c2_Iamr[1];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 46);
  c2_lamq3 = c2_Iamr[2];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 47);
  c2_lamd1 = c2_Iamr[3];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 48);
  c2_lamd2 = c2_Iamr[4];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 49);
  c2_lamd3 = c2_Iamr[5];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 50);
  c2_lamdf = c2_Iamr[6];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 53);
  c2_b_lamq1[0] = c2_lamq1;
  c2_b_lamq1[1] = c2_lamq2;
  c2_b_lamq1[2] = c2_lamq3;
  c2_c_P[0] = c2_b_P.Llq1;
  c2_c_P[1] = c2_b_P.Llq2;
  c2_c_P[2] = c2_b_P.Llq3;
  c2_rdivide(chartInstance, c2_b_lamq1, c2_c_P, c2_dv1);
  for (c2_i5 = 0; c2_i5 < 3; c2_i5++) {
    c2_dv2[c2_i5] = c2_dv1[c2_i5];
  }

  c2_A = c2_iqs + c2_sum(chartInstance, c2_dv2);
  c2_d_P[0] = c2_b_P.Llq1;
  c2_d_P[1] = c2_b_P.Llq2;
  c2_d_P[2] = c2_b_P.Llq3;
  c2_b_rdivide(chartInstance, 1.0, c2_d_P, c2_dv1);
  for (c2_i6 = 0; c2_i6 < 3; c2_i6++) {
    c2_dv3[c2_i6] = c2_dv1[c2_i6];
  }

  c2_B = 1.0 + c2_b_P.Lmq * c2_sum(chartInstance, c2_dv3);
  c2_x = c2_A;
  c2_y = c2_B;
  c2_b_x = c2_x;
  c2_b_y = c2_y;
  c2_c_x = c2_b_x;
  c2_c_y = c2_b_y;
  c2_imq = c2_c_x / c2_c_y;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 54);
  c2_b_A = c2_lamdf;
  c2_b_B = c2_b_P.Llfd;
  c2_d_x = c2_b_A;
  c2_d_y = c2_b_B;
  c2_e_x = c2_d_x;
  c2_e_y = c2_d_y;
  c2_f_x = c2_e_x;
  c2_f_y = c2_e_y;
  c2_g_y = c2_f_x / c2_f_y;
  c2_c_B = c2_b_P.Llfd;
  c2_h_y = c2_c_B;
  c2_i_y = c2_h_y;
  c2_j_y = c2_i_y;
  c2_k_y = 1.0 / c2_j_y;
  c2_l_y[0] = c2_b_P.Lld1;
  c2_l_y[1] = c2_b_P.Lld2;
  c2_l_y[2] = c2_b_P.Lld3;
  for (c2_i7 = 0; c2_i7 < 3; c2_i7++) {
    c2_l_y[c2_i7] = 1.0 / c2_l_y[c2_i7];
  }

  c2_b_lamd1[0] = c2_lamd1;
  c2_b_lamd1[1] = c2_lamd2;
  c2_b_lamd1[2] = c2_lamd3;
  c2_e_P[0] = c2_b_P.Lld1;
  c2_e_P[1] = c2_b_P.Lld2;
  c2_e_P[2] = c2_b_P.Lld3;
  c2_rdivide(chartInstance, c2_b_lamd1, c2_e_P, c2_dv1);
  for (c2_i8 = 0; c2_i8 < 3; c2_i8++) {
    c2_dv4[c2_i8] = c2_dv1[c2_i8];
  }

  c2_c_A = (c2_ids + c2_g_y) + c2_sum(chartInstance, c2_dv4);
  for (c2_i9 = 0; c2_i9 < 3; c2_i9++) {
    c2_m_y[c2_i9] = c2_l_y[c2_i9];
  }

  c2_d_B = 1.0 + c2_b_P.Lmd * (c2_k_y + c2_sum(chartInstance, c2_m_y));
  c2_g_x = c2_c_A;
  c2_n_y = c2_d_B;
  c2_h_x = c2_g_x;
  c2_o_y = c2_n_y;
  c2_i_x = c2_h_x;
  c2_p_y = c2_o_y;
  c2_imd = c2_i_x / c2_p_y;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 55);
  c2_d_A = c2_lamq1 - c2_b_P.Lmq * c2_imq;
  c2_e_B = c2_b_P.Llq1;
  c2_j_x = c2_d_A;
  c2_q_y = c2_e_B;
  c2_k_x = c2_j_x;
  c2_r_y = c2_q_y;
  c2_l_x = c2_k_x;
  c2_s_y = c2_r_y;
  c2_iq1 = c2_l_x / c2_s_y;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 56);
  c2_e_A = c2_lamq2 - c2_b_P.Lmq * c2_imq;
  c2_f_B = c2_b_P.Llq2;
  c2_m_x = c2_e_A;
  c2_t_y = c2_f_B;
  c2_n_x = c2_m_x;
  c2_u_y = c2_t_y;
  c2_o_x = c2_n_x;
  c2_v_y = c2_u_y;
  c2_iq2 = c2_o_x / c2_v_y;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 57);
  c2_f_A = c2_lamq3 - c2_b_P.Lmq * c2_imq;
  c2_g_B = c2_b_P.Llq3;
  c2_p_x = c2_f_A;
  c2_w_y = c2_g_B;
  c2_q_x = c2_p_x;
  c2_x_y = c2_w_y;
  c2_r_x = c2_q_x;
  c2_y_y = c2_x_y;
  c2_iq3 = c2_r_x / c2_y_y;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 58);
  c2_g_A = c2_lamd1 - c2_b_P.Lmd * c2_imd;
  c2_h_B = c2_b_P.Lld1;
  c2_s_x = c2_g_A;
  c2_ab_y = c2_h_B;
  c2_t_x = c2_s_x;
  c2_bb_y = c2_ab_y;
  c2_u_x = c2_t_x;
  c2_cb_y = c2_bb_y;
  c2_id1 = c2_u_x / c2_cb_y;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 59);
  c2_h_A = c2_lamd2 - c2_b_P.Lmd * c2_imd;
  c2_i_B = c2_b_P.Lld2;
  c2_v_x = c2_h_A;
  c2_db_y = c2_i_B;
  c2_w_x = c2_v_x;
  c2_eb_y = c2_db_y;
  c2_x_x = c2_w_x;
  c2_fb_y = c2_eb_y;
  c2_id2 = c2_x_x / c2_fb_y;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 60);
  c2_i_A = c2_lamd3 - c2_b_P.Lmd * c2_imd;
  c2_j_B = c2_b_P.Lld3;
  c2_y_x = c2_i_A;
  c2_gb_y = c2_j_B;
  c2_ab_x = c2_y_x;
  c2_hb_y = c2_gb_y;
  c2_bb_x = c2_ab_x;
  c2_ib_y = c2_hb_y;
  c2_id3 = c2_bb_x / c2_ib_y;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 61);
  c2_j_A = c2_lamdf - c2_b_P.Lmd * c2_imd;
  c2_k_B = c2_b_P.Llfd;
  c2_cb_x = c2_j_A;
  c2_jb_y = c2_k_B;
  c2_db_x = c2_cb_x;
  c2_kb_y = c2_jb_y;
  c2_eb_x = c2_db_x;
  c2_lb_y = c2_kb_y;
  c2_ifd = c2_eb_x / c2_lb_y;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 64);
  c2_plamq1 = -c2_b_P.rq1 * c2_iq1;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 65);
  c2_plamq2 = -c2_b_P.rq2 * c2_iq2;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 66);
  c2_plamq3 = -c2_b_P.rq3 * c2_iq3;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 69);
  c2_plamd1 = -c2_b_P.rd1 * c2_id1;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 70);
  c2_plamd2 = -c2_b_P.rd2 * c2_id2;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 71);
  c2_plamd3 = -c2_b_P.rd3 * c2_id3;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 74);
  c2_plamdf = c2_vfd - c2_b_P.rfd * c2_ifd;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 77);
  c2_b_plamq1[0] = c2_plamq1;
  c2_b_plamq1[1] = c2_plamq2;
  c2_b_plamq1[2] = c2_plamq3;
  c2_b_plamq1[3] = c2_plamd1;
  c2_b_plamq1[4] = c2_plamd2;
  c2_b_plamq1[5] = c2_plamd3;
  c2_b_plamq1[6] = c2_plamdf;
  for (c2_i10 = 0; c2_i10 < 7; c2_i10++) {
    c2_pIamr[c2_i10] = c2_b_plamq1[c2_i10];
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 80);
  c2_lamMq = c2_b_P.Lmq * c2_imq;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 81);
  c2_lamMd = c2_b_P.Lmd * c2_imd;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 84);
  c2_k_A = c2_b_P.pole;
  c2_fb_x = c2_k_A;
  c2_gb_x = c2_fb_x;
  c2_hb_x = c2_gb_x;
  c2_mb_y = c2_hb_x / 2.0;
  c2_wr = c2_mb_y * c2_wrm;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 87);
  c2_l_A = c2_b_P.pole;
  c2_ib_x = c2_l_A;
  c2_jb_x = c2_ib_x;
  c2_kb_x = c2_jb_x;
  c2_nb_y = c2_kb_x / 2.0;
  c2_Te = 1.5 * c2_nb_y * (c2_lamMd * c2_imq - c2_lamMq * c2_imd);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 89);
  c2_lamqS = c2_b_P.Lls * c2_iqs + c2_b_P.Lmq * c2_imq;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 90);
  c2_lamdS = c2_b_P.Lls * c2_ids + c2_b_P.Lmd * c2_imd;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 91);
  c2_vqs = c2_b_P.rs * c2_iqs + c2_wr * c2_lamdS;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 92);
  c2_vds = c2_b_P.rs * c2_ids - c2_wr * c2_lamqS;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, -92);
  _SFD_SYMBOL_SCOPE_POP();
  *c2_b_vqs = c2_vqs;
  *c2_b_vds = c2_vds;
  *c2_b_Te = c2_Te;
  *c2_b_ifd = c2_ifd;
  for (c2_i11 = 0; c2_i11 < 7; c2_i11++) {
    (*c2_b_pIamr)[c2_i11] = c2_pIamr[c2_i11];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 1U, chartInstance->c2_sfEvent);
}

static void initSimStructsc2_gen_rec(SFc2_gen_recInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void init_script_number_translation(uint32_T c2_machineNumber, uint32_T
  c2_chartNumber, uint32_T c2_instanceNumber)
{
  (void)c2_machineNumber;
  (void)c2_chartNumber;
  (void)c2_instanceNumber;
}

static const mxArray *c2_sf_marshallOut(void *chartInstanceVoid, void *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  int32_T c2_i12;
  real_T c2_b_inData[7];
  int32_T c2_i13;
  real_T c2_u[7];
  const mxArray *c2_y = NULL;
  SFc2_gen_recInstanceStruct *chartInstance;
  chartInstance = (SFc2_gen_recInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  for (c2_i12 = 0; c2_i12 < 7; c2_i12++) {
    c2_b_inData[c2_i12] = (*(real_T (*)[7])c2_inData)[c2_i12];
  }

  for (c2_i13 = 0; c2_i13 < 7; c2_i13++) {
    c2_u[c2_i13] = c2_b_inData[c2_i13];
  }

  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_u, 0, 0U, 1U, 0U, 1, 7), false);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, false);
  return c2_mxArrayOutData;
}

static void c2_emlrt_marshallIn(SFc2_gen_recInstanceStruct *chartInstance, const
  mxArray *c2_pIamr, const char_T *c2_identifier, real_T c2_y[7])
{
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_pIamr), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_pIamr);
}

static void c2_b_emlrt_marshallIn(SFc2_gen_recInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[7])
{
  real_T c2_dv5[7];
  int32_T c2_i14;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv5, 1, 0, 0U, 1, 0U, 1, 7);
  for (c2_i14 = 0; c2_i14 < 7; c2_i14++) {
    c2_y[c2_i14] = c2_dv5[c2_i14];
  }

  sf_mex_destroy(&c2_u);
}

static void c2_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_pIamr;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  real_T c2_y[7];
  int32_T c2_i15;
  SFc2_gen_recInstanceStruct *chartInstance;
  chartInstance = (SFc2_gen_recInstanceStruct *)chartInstanceVoid;
  c2_pIamr = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_pIamr), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_pIamr);
  for (c2_i15 = 0; c2_i15 < 7; c2_i15++) {
    (*(real_T (*)[7])c2_outData)[c2_i15] = c2_y[c2_i15];
  }

  sf_mex_destroy(&c2_mxArrayInData);
}

static const mxArray *c2_b_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  real_T c2_u;
  const mxArray *c2_y = NULL;
  SFc2_gen_recInstanceStruct *chartInstance;
  chartInstance = (SFc2_gen_recInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_u = *(real_T *)c2_inData;
  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", &c2_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, false);
  return c2_mxArrayOutData;
}

static real_T c2_c_emlrt_marshallIn(SFc2_gen_recInstanceStruct *chartInstance,
  const mxArray *c2_ifd, const char_T *c2_identifier)
{
  real_T c2_y;
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_y = c2_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_ifd), &c2_thisId);
  sf_mex_destroy(&c2_ifd);
  return c2_y;
}

static real_T c2_d_emlrt_marshallIn(SFc2_gen_recInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId)
{
  real_T c2_y;
  real_T c2_d0;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), &c2_d0, 1, 0, 0U, 0, 0U, 0);
  c2_y = c2_d0;
  sf_mex_destroy(&c2_u);
  return c2_y;
}

static void c2_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_ifd;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  real_T c2_y;
  SFc2_gen_recInstanceStruct *chartInstance;
  chartInstance = (SFc2_gen_recInstanceStruct *)chartInstanceVoid;
  c2_ifd = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_y = c2_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_ifd), &c2_thisId);
  sf_mex_destroy(&c2_ifd);
  *(real_T *)c2_outData = c2_y;
  sf_mex_destroy(&c2_mxArrayInData);
}

static const mxArray *c2_c_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  c2_struct_U3Zfiu1q58uXf1A0Qkd5GF c2_u;
  const mxArray *c2_y = NULL;
  real_T c2_b_u;
  const mxArray *c2_b_y = NULL;
  real_T c2_c_u;
  const mxArray *c2_c_y = NULL;
  real_T c2_d_u;
  const mxArray *c2_d_y = NULL;
  real_T c2_e_u;
  const mxArray *c2_e_y = NULL;
  real_T c2_f_u;
  const mxArray *c2_f_y = NULL;
  real_T c2_g_u;
  const mxArray *c2_g_y = NULL;
  real_T c2_h_u;
  const mxArray *c2_h_y = NULL;
  real_T c2_i_u;
  const mxArray *c2_i_y = NULL;
  real_T c2_j_u;
  const mxArray *c2_j_y = NULL;
  real_T c2_k_u;
  const mxArray *c2_k_y = NULL;
  real_T c2_l_u;
  const mxArray *c2_l_y = NULL;
  real_T c2_m_u;
  const mxArray *c2_m_y = NULL;
  real_T c2_n_u;
  const mxArray *c2_n_y = NULL;
  real_T c2_o_u;
  const mxArray *c2_o_y = NULL;
  real_T c2_p_u;
  const mxArray *c2_p_y = NULL;
  real_T c2_q_u;
  const mxArray *c2_q_y = NULL;
  real_T c2_r_u;
  const mxArray *c2_r_y = NULL;
  real_T c2_s_u;
  const mxArray *c2_s_y = NULL;
  real_T c2_t_u;
  const mxArray *c2_t_y = NULL;
  real_T c2_u_u;
  const mxArray *c2_u_y = NULL;
  real_T c2_v_u;
  const mxArray *c2_v_y = NULL;
  SFc2_gen_recInstanceStruct *chartInstance;
  chartInstance = (SFc2_gen_recInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_u = *(c2_struct_U3Zfiu1q58uXf1A0Qkd5GF *)c2_inData;
  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_createstruct("structure", 2, 1, 1), false);
  c2_b_u = c2_u.Ll;
  c2_b_y = NULL;
  sf_mex_assign(&c2_b_y, sf_mex_create("y", &c2_b_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_y, c2_b_y, "Ll", "Ll", 0);
  c2_c_u = c2_u.rl;
  c2_c_y = NULL;
  sf_mex_assign(&c2_c_y, sf_mex_create("y", &c2_c_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_y, c2_c_y, "rl", "rl", 0);
  c2_d_u = c2_u.rs;
  c2_d_y = NULL;
  sf_mex_assign(&c2_d_y, sf_mex_create("y", &c2_d_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_y, c2_d_y, "rs", "rs", 0);
  c2_e_u = c2_u.Lls;
  c2_e_y = NULL;
  sf_mex_assign(&c2_e_y, sf_mex_create("y", &c2_e_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_y, c2_e_y, "Lls", "Lls", 0);
  c2_f_u = c2_u.Lmq;
  c2_f_y = NULL;
  sf_mex_assign(&c2_f_y, sf_mex_create("y", &c2_f_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_y, c2_f_y, "Lmq", "Lmq", 0);
  c2_g_u = c2_u.Lmd;
  c2_g_y = NULL;
  sf_mex_assign(&c2_g_y, sf_mex_create("y", &c2_g_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_y, c2_g_y, "Lmd", "Lmd", 0);
  c2_h_u = c2_u.rd1;
  c2_h_y = NULL;
  sf_mex_assign(&c2_h_y, sf_mex_create("y", &c2_h_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_y, c2_h_y, "rd1", "rd1", 0);
  c2_i_u = c2_u.rd2;
  c2_i_y = NULL;
  sf_mex_assign(&c2_i_y, sf_mex_create("y", &c2_i_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_y, c2_i_y, "rd2", "rd2", 0);
  c2_j_u = c2_u.rd3;
  c2_j_y = NULL;
  sf_mex_assign(&c2_j_y, sf_mex_create("y", &c2_j_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_y, c2_j_y, "rd3", "rd3", 0);
  c2_k_u = c2_u.rfd;
  c2_k_y = NULL;
  sf_mex_assign(&c2_k_y, sf_mex_create("y", &c2_k_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_y, c2_k_y, "rfd", "rfd", 0);
  c2_l_u = c2_u.Lld1;
  c2_l_y = NULL;
  sf_mex_assign(&c2_l_y, sf_mex_create("y", &c2_l_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_y, c2_l_y, "Lld1", "Lld1", 0);
  c2_m_u = c2_u.Lld2;
  c2_m_y = NULL;
  sf_mex_assign(&c2_m_y, sf_mex_create("y", &c2_m_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_y, c2_m_y, "Lld2", "Lld2", 0);
  c2_n_u = c2_u.Lld3;
  c2_n_y = NULL;
  sf_mex_assign(&c2_n_y, sf_mex_create("y", &c2_n_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_y, c2_n_y, "Lld3", "Lld3", 0);
  c2_o_u = c2_u.Llfd;
  c2_o_y = NULL;
  sf_mex_assign(&c2_o_y, sf_mex_create("y", &c2_o_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_y, c2_o_y, "Llfd", "Llfd", 0);
  c2_p_u = c2_u.rq1;
  c2_p_y = NULL;
  sf_mex_assign(&c2_p_y, sf_mex_create("y", &c2_p_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_y, c2_p_y, "rq1", "rq1", 0);
  c2_q_u = c2_u.rq2;
  c2_q_y = NULL;
  sf_mex_assign(&c2_q_y, sf_mex_create("y", &c2_q_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_y, c2_q_y, "rq2", "rq2", 0);
  c2_r_u = c2_u.rq3;
  c2_r_y = NULL;
  sf_mex_assign(&c2_r_y, sf_mex_create("y", &c2_r_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_y, c2_r_y, "rq3", "rq3", 0);
  c2_s_u = c2_u.pole;
  c2_s_y = NULL;
  sf_mex_assign(&c2_s_y, sf_mex_create("y", &c2_s_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_y, c2_s_y, "pole", "pole", 0);
  c2_t_u = c2_u.Llq1;
  c2_t_y = NULL;
  sf_mex_assign(&c2_t_y, sf_mex_create("y", &c2_t_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_y, c2_t_y, "Llq1", "Llq1", 0);
  c2_u_u = c2_u.Llq2;
  c2_u_y = NULL;
  sf_mex_assign(&c2_u_y, sf_mex_create("y", &c2_u_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_y, c2_u_y, "Llq2", "Llq2", 0);
  c2_v_u = c2_u.Llq3;
  c2_v_y = NULL;
  sf_mex_assign(&c2_v_y, sf_mex_create("y", &c2_v_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_y, c2_v_y, "Llq3", "Llq3", 0);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, false);
  return c2_mxArrayOutData;
}

static void c2_e_emlrt_marshallIn(SFc2_gen_recInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  c2_struct_U3Zfiu1q58uXf1A0Qkd5GF *c2_y)
{
  emlrtMsgIdentifier c2_thisId;
  static const char * c2_fieldNames[21] = { "Ll", "rl", "rs", "Lls", "Lmq",
    "Lmd", "rd1", "rd2", "rd3", "rfd", "Lld1", "Lld2", "Lld3", "Llfd", "rq1",
    "rq2", "rq3", "pole", "Llq1", "Llq2", "Llq3" };

  c2_thisId.fParent = c2_parentId;
  sf_mex_check_struct(c2_parentId, c2_u, 21, c2_fieldNames, 0U, NULL);
  c2_thisId.fIdentifier = "Ll";
  c2_y->Ll = c2_d_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c2_u, "Ll", "Ll", 0)), &c2_thisId);
  c2_thisId.fIdentifier = "rl";
  c2_y->rl = c2_d_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c2_u, "rl", "rl", 0)), &c2_thisId);
  c2_thisId.fIdentifier = "rs";
  c2_y->rs = c2_d_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c2_u, "rs", "rs", 0)), &c2_thisId);
  c2_thisId.fIdentifier = "Lls";
  c2_y->Lls = c2_d_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c2_u, "Lls", "Lls", 0)), &c2_thisId);
  c2_thisId.fIdentifier = "Lmq";
  c2_y->Lmq = c2_d_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c2_u, "Lmq", "Lmq", 0)), &c2_thisId);
  c2_thisId.fIdentifier = "Lmd";
  c2_y->Lmd = c2_d_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c2_u, "Lmd", "Lmd", 0)), &c2_thisId);
  c2_thisId.fIdentifier = "rd1";
  c2_y->rd1 = c2_d_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c2_u, "rd1", "rd1", 0)), &c2_thisId);
  c2_thisId.fIdentifier = "rd2";
  c2_y->rd2 = c2_d_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c2_u, "rd2", "rd2", 0)), &c2_thisId);
  c2_thisId.fIdentifier = "rd3";
  c2_y->rd3 = c2_d_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c2_u, "rd3", "rd3", 0)), &c2_thisId);
  c2_thisId.fIdentifier = "rfd";
  c2_y->rfd = c2_d_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c2_u, "rfd", "rfd", 0)), &c2_thisId);
  c2_thisId.fIdentifier = "Lld1";
  c2_y->Lld1 = c2_d_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c2_u, "Lld1", "Lld1", 0)), &c2_thisId);
  c2_thisId.fIdentifier = "Lld2";
  c2_y->Lld2 = c2_d_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c2_u, "Lld2", "Lld2", 0)), &c2_thisId);
  c2_thisId.fIdentifier = "Lld3";
  c2_y->Lld3 = c2_d_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c2_u, "Lld3", "Lld3", 0)), &c2_thisId);
  c2_thisId.fIdentifier = "Llfd";
  c2_y->Llfd = c2_d_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c2_u, "Llfd", "Llfd", 0)), &c2_thisId);
  c2_thisId.fIdentifier = "rq1";
  c2_y->rq1 = c2_d_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c2_u, "rq1", "rq1", 0)), &c2_thisId);
  c2_thisId.fIdentifier = "rq2";
  c2_y->rq2 = c2_d_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c2_u, "rq2", "rq2", 0)), &c2_thisId);
  c2_thisId.fIdentifier = "rq3";
  c2_y->rq3 = c2_d_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c2_u, "rq3", "rq3", 0)), &c2_thisId);
  c2_thisId.fIdentifier = "pole";
  c2_y->pole = c2_d_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c2_u, "pole", "pole", 0)), &c2_thisId);
  c2_thisId.fIdentifier = "Llq1";
  c2_y->Llq1 = c2_d_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c2_u, "Llq1", "Llq1", 0)), &c2_thisId);
  c2_thisId.fIdentifier = "Llq2";
  c2_y->Llq2 = c2_d_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c2_u, "Llq2", "Llq2", 0)), &c2_thisId);
  c2_thisId.fIdentifier = "Llq3";
  c2_y->Llq3 = c2_d_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c2_u, "Llq3", "Llq3", 0)), &c2_thisId);
  sf_mex_destroy(&c2_u);
}

static void c2_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_b_P;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  c2_struct_U3Zfiu1q58uXf1A0Qkd5GF c2_y;
  SFc2_gen_recInstanceStruct *chartInstance;
  chartInstance = (SFc2_gen_recInstanceStruct *)chartInstanceVoid;
  c2_b_P = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_P), &c2_thisId, &c2_y);
  sf_mex_destroy(&c2_b_P);
  *(c2_struct_U3Zfiu1q58uXf1A0Qkd5GF *)c2_outData = c2_y;
  sf_mex_destroy(&c2_mxArrayInData);
}

const mxArray *sf_c2_gen_rec_get_eml_resolved_functions_info(void)
{
  const mxArray *c2_nameCaptureInfo = NULL;
  c2_nameCaptureInfo = NULL;
  sf_mex_assign(&c2_nameCaptureInfo, sf_mex_createstruct("structure", 2, 23, 1),
                false);
  c2_info_helper(&c2_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c2_nameCaptureInfo);
  return c2_nameCaptureInfo;
}

static void c2_info_helper(const mxArray **c2_info)
{
  const mxArray *c2_rhs0 = NULL;
  const mxArray *c2_lhs0 = NULL;
  const mxArray *c2_rhs1 = NULL;
  const mxArray *c2_lhs1 = NULL;
  const mxArray *c2_rhs2 = NULL;
  const mxArray *c2_lhs2 = NULL;
  const mxArray *c2_rhs3 = NULL;
  const mxArray *c2_lhs3 = NULL;
  const mxArray *c2_rhs4 = NULL;
  const mxArray *c2_lhs4 = NULL;
  const mxArray *c2_rhs5 = NULL;
  const mxArray *c2_lhs5 = NULL;
  const mxArray *c2_rhs6 = NULL;
  const mxArray *c2_lhs6 = NULL;
  const mxArray *c2_rhs7 = NULL;
  const mxArray *c2_lhs7 = NULL;
  const mxArray *c2_rhs8 = NULL;
  const mxArray *c2_lhs8 = NULL;
  const mxArray *c2_rhs9 = NULL;
  const mxArray *c2_lhs9 = NULL;
  const mxArray *c2_rhs10 = NULL;
  const mxArray *c2_lhs10 = NULL;
  const mxArray *c2_rhs11 = NULL;
  const mxArray *c2_lhs11 = NULL;
  const mxArray *c2_rhs12 = NULL;
  const mxArray *c2_lhs12 = NULL;
  const mxArray *c2_rhs13 = NULL;
  const mxArray *c2_lhs13 = NULL;
  const mxArray *c2_rhs14 = NULL;
  const mxArray *c2_lhs14 = NULL;
  const mxArray *c2_rhs15 = NULL;
  const mxArray *c2_lhs15 = NULL;
  const mxArray *c2_rhs16 = NULL;
  const mxArray *c2_lhs16 = NULL;
  const mxArray *c2_rhs17 = NULL;
  const mxArray *c2_lhs17 = NULL;
  const mxArray *c2_rhs18 = NULL;
  const mxArray *c2_lhs18 = NULL;
  const mxArray *c2_rhs19 = NULL;
  const mxArray *c2_lhs19 = NULL;
  const mxArray *c2_rhs20 = NULL;
  const mxArray *c2_lhs20 = NULL;
  const mxArray *c2_rhs21 = NULL;
  const mxArray *c2_lhs21 = NULL;
  const mxArray *c2_rhs22 = NULL;
  const mxArray *c2_lhs22 = NULL;
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("rdivide"), "name", "name", 0);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c2_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs0), "rhs", "rhs", 0);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs0), "lhs", "lhs", 0);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 1);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 1);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c2_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs1), "rhs", "rhs", 1);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs1), "lhs", "lhs", 1);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 2);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 2);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 2);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1286825996U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c2_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs2), "rhs", "rhs", 2);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs2), "lhs", "lhs", 2);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "context", "context", 3);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("isequal"), "name", "name", 3);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m"), "resolved",
                  "resolved", 3);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1286825958U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c2_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs3), "rhs", "rhs", 3);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs3), "lhs", "lhs", 3);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m"), "context",
                  "context", 4);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_isequal_core"), "name",
                  "name", 4);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isequal_core.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1286825986U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c2_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs4), "rhs", "rhs", 4);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs4), "lhs", "lhs", 4);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isequal_core.m!isequal_scalar"),
                  "context", "context", 5);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("isnan"), "name", "name", 5);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 5);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1363717458U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c2_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs5), "rhs", "rhs", 5);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs5), "lhs", "lhs", 5);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "context",
                  "context", 6);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 6);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 6);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c2_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs6), "rhs", "rhs", 6);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs6), "lhs", "lhs", 6);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 7);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_div"), "name", "name", 7);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 7);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c2_rhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs7), "rhs", "rhs", 7);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs7), "lhs", "lhs", 7);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "context",
                  "context", 8);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("coder.internal.div"), "name",
                  "name", 8);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p"), "resolved",
                  "resolved", 8);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c2_rhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs8), "rhs", "rhs", 8);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs8), "lhs", "lhs", 8);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "context", "context", 9);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("sum"), "name", "name", 9);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/sum.m"), "resolved",
                  "resolved", 9);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1363717458U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c2_rhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs9), "rhs", "rhs", 9);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs9), "lhs", "lhs", 9);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/sum.m"), "context",
                  "context", 10);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 10);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 10);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c2_rhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/sum.m"), "context",
                  "context", 11);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("isequal"), "name", "name", 11);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m"), "resolved",
                  "resolved", 11);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1286825958U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c2_rhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/sum.m"), "context",
                  "context", 12);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_const_nonsingleton_dim"),
                  "name", "name", 12);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_const_nonsingleton_dim.m"),
                  "resolved", "resolved", 12);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c2_rhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_const_nonsingleton_dim.m"),
                  "context", "context", 13);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "coder.internal.constNonSingletonDim"), "name", "name", 13);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/constNonSingletonDim.m"),
                  "resolved", "resolved", 13);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c2_rhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/sum.m"), "context",
                  "context", 14);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 14);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 14);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c2_rhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 15);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 15);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 15);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c2_rhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/sum.m"), "context",
                  "context", 16);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 16);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 16);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c2_rhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/sum.m"), "context",
                  "context", 17);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 17);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 17);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c2_rhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 18);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("intmax"), "name", "name", 18);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 18);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c2_rhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "context",
                  "context", 19);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 19);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 19);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c2_rhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "context", "context", 20);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("mrdivide"), "name", "name", 20);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 20);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1388463696U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1370017086U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c2_rhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 21);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 21);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 21);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c2_rhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 22);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("rdivide"), "name", "name", 22);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 22);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c2_rhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs22), "lhs", "lhs",
                  22);
  sf_mex_destroy(&c2_rhs0);
  sf_mex_destroy(&c2_lhs0);
  sf_mex_destroy(&c2_rhs1);
  sf_mex_destroy(&c2_lhs1);
  sf_mex_destroy(&c2_rhs2);
  sf_mex_destroy(&c2_lhs2);
  sf_mex_destroy(&c2_rhs3);
  sf_mex_destroy(&c2_lhs3);
  sf_mex_destroy(&c2_rhs4);
  sf_mex_destroy(&c2_lhs4);
  sf_mex_destroy(&c2_rhs5);
  sf_mex_destroy(&c2_lhs5);
  sf_mex_destroy(&c2_rhs6);
  sf_mex_destroy(&c2_lhs6);
  sf_mex_destroy(&c2_rhs7);
  sf_mex_destroy(&c2_lhs7);
  sf_mex_destroy(&c2_rhs8);
  sf_mex_destroy(&c2_lhs8);
  sf_mex_destroy(&c2_rhs9);
  sf_mex_destroy(&c2_lhs9);
  sf_mex_destroy(&c2_rhs10);
  sf_mex_destroy(&c2_lhs10);
  sf_mex_destroy(&c2_rhs11);
  sf_mex_destroy(&c2_lhs11);
  sf_mex_destroy(&c2_rhs12);
  sf_mex_destroy(&c2_lhs12);
  sf_mex_destroy(&c2_rhs13);
  sf_mex_destroy(&c2_lhs13);
  sf_mex_destroy(&c2_rhs14);
  sf_mex_destroy(&c2_lhs14);
  sf_mex_destroy(&c2_rhs15);
  sf_mex_destroy(&c2_lhs15);
  sf_mex_destroy(&c2_rhs16);
  sf_mex_destroy(&c2_lhs16);
  sf_mex_destroy(&c2_rhs17);
  sf_mex_destroy(&c2_lhs17);
  sf_mex_destroy(&c2_rhs18);
  sf_mex_destroy(&c2_lhs18);
  sf_mex_destroy(&c2_rhs19);
  sf_mex_destroy(&c2_lhs19);
  sf_mex_destroy(&c2_rhs20);
  sf_mex_destroy(&c2_lhs20);
  sf_mex_destroy(&c2_rhs21);
  sf_mex_destroy(&c2_lhs21);
  sf_mex_destroy(&c2_rhs22);
  sf_mex_destroy(&c2_lhs22);
}

static const mxArray *c2_emlrt_marshallOut(const char * c2_u)
{
  const mxArray *c2_y = NULL;
  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c2_u)), false);
  return c2_y;
}

static const mxArray *c2_b_emlrt_marshallOut(const uint32_T c2_u)
{
  const mxArray *c2_y = NULL;
  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", &c2_u, 7, 0U, 0U, 0U, 0), false);
  return c2_y;
}

static void c2_rdivide(SFc2_gen_recInstanceStruct *chartInstance, real_T c2_x[3],
  real_T c2_y[3], real_T c2_z[3])
{
  int32_T c2_i16;
  (void)chartInstance;
  for (c2_i16 = 0; c2_i16 < 3; c2_i16++) {
    c2_z[c2_i16] = c2_x[c2_i16] / c2_y[c2_i16];
  }
}

static real_T c2_sum(SFc2_gen_recInstanceStruct *chartInstance, real_T c2_x[3])
{
  real_T c2_y;
  int32_T c2_k;
  int32_T c2_b_k;
  (void)chartInstance;
  c2_y = c2_x[0];
  for (c2_k = 2; c2_k < 4; c2_k++) {
    c2_b_k = c2_k;
    c2_y += c2_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c2_b_k), 1, 3, 1, 0) - 1];
  }

  return c2_y;
}

static void c2_b_rdivide(SFc2_gen_recInstanceStruct *chartInstance, real_T c2_x,
  real_T c2_y[3], real_T c2_z[3])
{
  real_T c2_b_x;
  real_T c2_c_x;
  int32_T c2_i17;
  (void)chartInstance;
  c2_b_x = c2_x;
  c2_c_x = c2_b_x;
  for (c2_i17 = 0; c2_i17 < 3; c2_i17++) {
    c2_z[c2_i17] = c2_c_x / c2_y[c2_i17];
  }
}

static const mxArray *c2_d_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  int32_T c2_u;
  const mxArray *c2_y = NULL;
  SFc2_gen_recInstanceStruct *chartInstance;
  chartInstance = (SFc2_gen_recInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_u = *(int32_T *)c2_inData;
  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", &c2_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, false);
  return c2_mxArrayOutData;
}

static int32_T c2_f_emlrt_marshallIn(SFc2_gen_recInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId)
{
  int32_T c2_y;
  int32_T c2_i18;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), &c2_i18, 1, 6, 0U, 0, 0U, 0);
  c2_y = c2_i18;
  sf_mex_destroy(&c2_u);
  return c2_y;
}

static void c2_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_b_sfEvent;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  int32_T c2_y;
  SFc2_gen_recInstanceStruct *chartInstance;
  chartInstance = (SFc2_gen_recInstanceStruct *)chartInstanceVoid;
  c2_b_sfEvent = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_y = c2_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_sfEvent),
    &c2_thisId);
  sf_mex_destroy(&c2_b_sfEvent);
  *(int32_T *)c2_outData = c2_y;
  sf_mex_destroy(&c2_mxArrayInData);
}

static uint8_T c2_g_emlrt_marshallIn(SFc2_gen_recInstanceStruct *chartInstance,
  const mxArray *c2_b_is_active_c2_gen_rec, const char_T *c2_identifier)
{
  uint8_T c2_y;
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_y = c2_h_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c2_b_is_active_c2_gen_rec), &c2_thisId);
  sf_mex_destroy(&c2_b_is_active_c2_gen_rec);
  return c2_y;
}

static uint8_T c2_h_emlrt_marshallIn(SFc2_gen_recInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId)
{
  uint8_T c2_y;
  uint8_T c2_u0;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), &c2_u0, 1, 3, 0U, 0, 0U, 0);
  c2_y = c2_u0;
  sf_mex_destroy(&c2_u);
  return c2_y;
}

static void init_dsm_address_info(SFc2_gen_recInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

/* SFunction Glue Code */
#ifdef utFree
#undef utFree
#endif

#ifdef utMalloc
#undef utMalloc
#endif

#ifdef __cplusplus

extern "C" void *utMalloc(size_t size);
extern "C" void utFree(void*);

#else

extern void *utMalloc(size_t size);
extern void utFree(void*);

#endif

void sf_c2_gen_rec_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(1298447359U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(1290660304U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(816786862U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(2438715312U);
}

mxArray *sf_c2_gen_rec_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("Y1jp5TkprDx1IpJCURejZD");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,5,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,3,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,3,"type",mxType);
    }

    mxSetField(mxData,3,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(7);
      pr[1] = (double)(1);
      mxSetField(mxData,4,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,4,"type",mxType);
    }

    mxSetField(mxData,4,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,1,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(13));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxData);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,5,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,3,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,3,"type",mxType);
    }

    mxSetField(mxData,3,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(7);
      pr[1] = (double)(1);
      mxSetField(mxData,4,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,4,"type",mxType);
    }

    mxSetField(mxData,4,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c2_gen_rec_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c2_gen_rec_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c2_gen_rec(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x6'type','srcId','name','auxInfo'{{M[1],M[11],T\"Te\",},{M[1],M[12],T\"ifd\",},{M[1],M[13],T\"pIamr\",},{M[1],M[10],T\"vds\",},{M[1],M[5],T\"vqs\",},{M[8],M[0],T\"is_active_c2_gen_rec\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 6, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c2_gen_rec_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc2_gen_recInstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc2_gen_recInstanceStruct *) chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _gen_recMachineNumber_,
           2,
           1,
           1,
           0,
           11,
           0,
           0,
           0,
           0,
           0,
           &(chartInstance->chartNumber),
           &(chartInstance->instanceNumber),
           (void *)S);

        /* Each instance must initialize ist own list of scripts */
        init_script_number_translation(_gen_recMachineNumber_,
          chartInstance->chartNumber,chartInstance->instanceNumber);
        if (chartAlreadyPresent==0) {
          /* this is the first instance */
          sf_debug_set_chart_disable_implicit_casting
            (sfGlobalDebugInstanceStruct,_gen_recMachineNumber_,
             chartInstance->chartNumber,1);
          sf_debug_set_chart_event_thresholds(sfGlobalDebugInstanceStruct,
            _gen_recMachineNumber_,
            chartInstance->chartNumber,
            0,
            0,
            0);
          _SFD_SET_DATA_PROPS(0,1,1,0,"vfd");
          _SFD_SET_DATA_PROPS(1,2,0,1,"vqs");
          _SFD_SET_DATA_PROPS(2,1,1,0,"wrm");
          _SFD_SET_DATA_PROPS(3,1,1,0,"iqs");
          _SFD_SET_DATA_PROPS(4,1,1,0,"ids");
          _SFD_SET_DATA_PROPS(5,1,1,0,"Iamr");
          _SFD_SET_DATA_PROPS(6,2,0,1,"vds");
          _SFD_SET_DATA_PROPS(7,2,0,1,"Te");
          _SFD_SET_DATA_PROPS(8,2,0,1,"ifd");
          _SFD_SET_DATA_PROPS(9,2,0,1,"pIamr");
          _SFD_SET_DATA_PROPS(10,10,0,0,"P");
          _SFD_STATE_INFO(0,0,2);
          _SFD_CH_SUBSTATE_COUNT(0);
          _SFD_CH_SUBSTATE_DECOMP(0);
        }

        _SFD_CV_INIT_CHART(0,0,0,0);

        {
          _SFD_CV_INIT_STATE(0,0,0,0,0,0,NULL,NULL);
        }

        _SFD_CV_INIT_TRANS(0,0,NULL,NULL,0,NULL);

        /* Initialization of MATLAB Function Model Coverage */
        _SFD_CV_INIT_EML(0,1,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,3328);
        _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c2_b_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c2_b_sf_marshallOut,(MexInFcnForType)c2_b_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c2_b_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c2_b_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c2_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 7;
          _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c2_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(6,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c2_b_sf_marshallOut,(MexInFcnForType)c2_b_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(7,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c2_b_sf_marshallOut,(MexInFcnForType)c2_b_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(8,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c2_b_sf_marshallOut,(MexInFcnForType)c2_b_sf_marshallIn);

        {
          unsigned int dimVector[1];
          dimVector[0]= 7;
          _SFD_SET_DATA_COMPILED_PROPS(9,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c2_sf_marshallOut,(MexInFcnForType)
            c2_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(10,SF_STRUCT,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c2_c_sf_marshallOut,(MexInFcnForType)c2_c_sf_marshallIn);

        {
          real_T *c2_vfd;
          real_T *c2_vqs;
          real_T *c2_wrm;
          real_T *c2_iqs;
          real_T *c2_ids;
          real_T *c2_vds;
          real_T *c2_Te;
          real_T *c2_ifd;
          real_T (*c2_Iamr)[7];
          real_T (*c2_pIamr)[7];
          c2_pIamr = (real_T (*)[7])ssGetOutputPortSignal(chartInstance->S, 5);
          c2_ifd = (real_T *)ssGetOutputPortSignal(chartInstance->S, 4);
          c2_Te = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
          c2_vds = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
          c2_Iamr = (real_T (*)[7])ssGetInputPortSignal(chartInstance->S, 4);
          c2_ids = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
          c2_iqs = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
          c2_wrm = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
          c2_vqs = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
          c2_vfd = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, c2_vfd);
          _SFD_SET_DATA_VALUE_PTR(1U, c2_vqs);
          _SFD_SET_DATA_VALUE_PTR(2U, c2_wrm);
          _SFD_SET_DATA_VALUE_PTR(3U, c2_iqs);
          _SFD_SET_DATA_VALUE_PTR(4U, c2_ids);
          _SFD_SET_DATA_VALUE_PTR(5U, *c2_Iamr);
          _SFD_SET_DATA_VALUE_PTR(6U, c2_vds);
          _SFD_SET_DATA_VALUE_PTR(7U, c2_Te);
          _SFD_SET_DATA_VALUE_PTR(8U, c2_ifd);
          _SFD_SET_DATA_VALUE_PTR(9U, *c2_pIamr);
          _SFD_SET_DATA_VALUE_PTR(10U, &chartInstance->c2_P);
        }
      }
    } else {
      sf_debug_reset_current_state_configuration(sfGlobalDebugInstanceStruct,
        _gen_recMachineNumber_,chartInstance->chartNumber,
        chartInstance->instanceNumber);
    }
  }
}

static const char* sf_get_instance_specialization(void)
{
  return "rYHEH12dlPHmVzGdIAvTNE";
}

static void sf_opaque_initialize_c2_gen_rec(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc2_gen_recInstanceStruct*) chartInstanceVar)->S,
    0);
  initialize_params_c2_gen_rec((SFc2_gen_recInstanceStruct*) chartInstanceVar);
  initialize_c2_gen_rec((SFc2_gen_recInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c2_gen_rec(void *chartInstanceVar)
{
  enable_c2_gen_rec((SFc2_gen_recInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c2_gen_rec(void *chartInstanceVar)
{
  disable_c2_gen_rec((SFc2_gen_recInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c2_gen_rec(void *chartInstanceVar)
{
  sf_gateway_c2_gen_rec((SFc2_gen_recInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c2_gen_rec(SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c2_gen_rec((SFc2_gen_recInstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c2_gen_rec();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_raw2high'.\n");
  }

  return plhs[0];
}

extern void sf_internal_set_sim_state_c2_gen_rec(SimStruct* S, const mxArray *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[3];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxDuplicateArray(st);      /* high level simctx */
  prhs[2] = (mxArray*) sf_get_sim_state_info_c2_gen_rec();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 3, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c2_gen_rec((SFc2_gen_recInstanceStruct*)chartInfo->chartInstance,
    mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c2_gen_rec(SimStruct* S)
{
  return sf_internal_get_sim_state_c2_gen_rec(S);
}

static void sf_opaque_set_sim_state_c2_gen_rec(SimStruct* S, const mxArray *st)
{
  sf_internal_set_sim_state_c2_gen_rec(S, st);
}

static void sf_opaque_terminate_c2_gen_rec(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc2_gen_recInstanceStruct*) chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_gen_rec_optimization_info();
    }

    finalize_c2_gen_rec((SFc2_gen_recInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    if (crtInfo != NULL) {
      utFree((void *)crtInfo);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc2_gen_rec((SFc2_gen_recInstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c2_gen_rec(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    initialize_params_c2_gen_rec((SFc2_gen_recInstanceStruct*)
      (chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c2_gen_rec(SimStruct *S)
{
  /* Actual parameters from chart:
     P
   */
  const char_T *rtParamNames[] = { "P" };

  ssSetNumRunTimeParams(S,ssGetSFcnParamsCount(S));
  ssRegDlgParamAsRunTimeParam(S, 0, 0, rtParamNames[0],
    sf_get_param_data_type_id(S,0));
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_gen_rec_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,2);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,2,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,2,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,2);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 4, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,2,5);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,2,5);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=5; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 5; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,2);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(3353885895U));
  ssSetChecksum1(S,(377425455U));
  ssSetChecksum2(S,(349194156U));
  ssSetChecksum3(S,(93618966U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c2_gen_rec(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c2_gen_rec(SimStruct *S)
{
  SFc2_gen_recInstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc2_gen_recInstanceStruct *)utMalloc(sizeof
    (SFc2_gen_recInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc2_gen_recInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c2_gen_rec;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c2_gen_rec;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c2_gen_rec;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c2_gen_rec;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c2_gen_rec;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c2_gen_rec;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c2_gen_rec;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c2_gen_rec;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c2_gen_rec;
  chartInstance->chartInfo.mdlStart = mdlStart_c2_gen_rec;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c2_gen_rec;
  chartInstance->chartInfo.extModeExec = NULL;
  chartInstance->chartInfo.restoreLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.restoreBeforeLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.storeCurrentConfiguration = NULL;
  chartInstance->chartInfo.debugInstance = sfGlobalDebugInstanceStruct;
  chartInstance->S = S;
  crtInfo->instanceInfo = (&(chartInstance->chartInfo));
  crtInfo->isJITEnabled = false;
  ssSetUserData(S,(void *)(crtInfo));  /* register the chart instance with simstruct */
  init_dsm_address_info(chartInstance);
  if (!sim_mode_is_rtw_gen(S)) {
  }

  sf_opaque_init_subchart_simstructs(chartInstance->chartInfo.chartInstance);
  chart_debug_initialization(S,1);
}

void c2_gen_rec_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c2_gen_rec(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c2_gen_rec(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c2_gen_rec(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c2_gen_rec_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
