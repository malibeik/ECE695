/* Include files */

#include <stddef.h>
#include "blas.h"
#include "gen_rec_sfun.h"
#include "c1_gen_rec.h"
#include "mwmathutil.h"
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
static const char * c1_debug_family_names[67] = { "lamq1", "lamq2", "lamq3",
  "lamd1", "lamd2", "lamd3", "lamdf", "wr", "LdDP", "LqDP", "lamdDP", "lamqDP",
  "x", "epsilon", "delta", "max_iter", "record", "fx", "n", "fp", "d", "betamin",
  "eqDP", "edDP", "phi", "alpha", "Epk", "Lc", "Lt", "K", "iqsCond", "idsCond",
  "theta", "ias1", "iqs1", "ids1", "ias2", "iqs2", "ids2", "ias3", "iqs3",
  "ids3", "ias4", "iqs4", "ids4", "ias5", "iqs5", "ids5", "iqs_C", "ids_C",
  "iqsCommu", "idsCommu", "nargin", "nargout", "wrm", "betac", "ed", "Iamr",
  "idc", "betaminprev", "uprev", "P", "iqs", "ids", "pidc", "beta", "u" };

/* Function Declarations */
static void initialize_c1_gen_rec(SFc1_gen_recInstanceStruct *chartInstance);
static void initialize_params_c1_gen_rec(SFc1_gen_recInstanceStruct
  *chartInstance);
static void enable_c1_gen_rec(SFc1_gen_recInstanceStruct *chartInstance);
static void disable_c1_gen_rec(SFc1_gen_recInstanceStruct *chartInstance);
static void c1_update_debugger_state_c1_gen_rec(SFc1_gen_recInstanceStruct
  *chartInstance);
static const mxArray *get_sim_state_c1_gen_rec(SFc1_gen_recInstanceStruct
  *chartInstance);
static void set_sim_state_c1_gen_rec(SFc1_gen_recInstanceStruct *chartInstance,
  const mxArray *c1_st);
static void finalize_c1_gen_rec(SFc1_gen_recInstanceStruct *chartInstance);
static void sf_gateway_c1_gen_rec(SFc1_gen_recInstanceStruct *chartInstance);
static void c1_chartstep_c1_gen_rec(SFc1_gen_recInstanceStruct *chartInstance);
static void initSimStructsc1_gen_rec(SFc1_gen_recInstanceStruct *chartInstance);
static void init_script_number_translation(uint32_T c1_machineNumber, uint32_T
  c1_chartNumber, uint32_T c1_instanceNumber);
static const mxArray *c1_sf_marshallOut(void *chartInstanceVoid, void *c1_inData);
static real_T c1_emlrt_marshallIn(SFc1_gen_recInstanceStruct *chartInstance,
  const mxArray *c1_u, const char_T *c1_identifier);
static real_T c1_b_emlrt_marshallIn(SFc1_gen_recInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_b_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_c_emlrt_marshallIn(SFc1_gen_recInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  c1_struct_U3Zfiu1q58uXf1A0Qkd5GF *c1_y);
static void c1_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_c_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static const mxArray *c1_d_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_d_emlrt_marshallIn(SFc1_gen_recInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[5]);
static void c1_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_e_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_e_emlrt_marshallIn(SFc1_gen_recInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[100]);
static void c1_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static void c1_info_helper(const mxArray **c1_info);
static const mxArray *c1_emlrt_marshallOut(const char * c1_u);
static const mxArray *c1_b_emlrt_marshallOut(const uint32_T c1_u);
static void c1_rdivide(SFc1_gen_recInstanceStruct *chartInstance, real_T c1_x,
  real_T c1_y[3], real_T c1_z[3]);
static real_T c1_sum(SFc1_gen_recInstanceStruct *chartInstance, real_T c1_x[3]);
static void c1_b_rdivide(SFc1_gen_recInstanceStruct *chartInstance, real_T c1_x
  [3], real_T c1_y[3], real_T c1_z[3]);
static real_T c1_sqrt(SFc1_gen_recInstanceStruct *chartInstance, real_T c1_x);
static void c1_eml_error(SFc1_gen_recInstanceStruct *chartInstance);
static real_T c1_abs(SFc1_gen_recInstanceStruct *chartInstance, real_T c1_x);
static real_T c1_angle(SFc1_gen_recInstanceStruct *chartInstance, creal_T c1_x);
static real_T c1_mpower(SFc1_gen_recInstanceStruct *chartInstance, real_T c1_a);
static void c1_eml_scalar_eg(SFc1_gen_recInstanceStruct *chartInstance);
static const mxArray *c1_f_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static int32_T c1_f_emlrt_marshallIn(SFc1_gen_recInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static uint8_T c1_g_emlrt_marshallIn(SFc1_gen_recInstanceStruct *chartInstance,
  const mxArray *c1_b_is_active_c1_gen_rec, const char_T *c1_identifier);
static uint8_T c1_h_emlrt_marshallIn(SFc1_gen_recInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_b_sqrt(SFc1_gen_recInstanceStruct *chartInstance, real_T *c1_x);
static void init_dsm_address_info(SFc1_gen_recInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c1_gen_rec(SFc1_gen_recInstanceStruct *chartInstance)
{
  chartInstance->c1_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c1_is_active_c1_gen_rec = 0U;
}

static void initialize_params_c1_gen_rec(SFc1_gen_recInstanceStruct
  *chartInstance)
{
  const mxArray *c1_m0 = NULL;
  const mxArray *c1_mxField;
  c1_struct_U3Zfiu1q58uXf1A0Qkd5GF c1_r0;
  c1_m0 = sf_mex_get_sfun_param(chartInstance->S, 0, 1);
  c1_mxField = sf_mex_getfield(c1_m0, "Ll", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c1_mxField), &c1_r0.Ll, 1, 0, 0U, 0, 0U, 0);
  c1_mxField = sf_mex_getfield(c1_m0, "rl", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c1_mxField), &c1_r0.rl, 1, 0, 0U, 0, 0U, 0);
  c1_mxField = sf_mex_getfield(c1_m0, "rs", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c1_mxField), &c1_r0.rs, 1, 0, 0U, 0, 0U, 0);
  c1_mxField = sf_mex_getfield(c1_m0, "Lls", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c1_mxField), &c1_r0.Lls, 1, 0, 0U, 0, 0U,
                      0);
  c1_mxField = sf_mex_getfield(c1_m0, "Lmq", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c1_mxField), &c1_r0.Lmq, 1, 0, 0U, 0, 0U,
                      0);
  c1_mxField = sf_mex_getfield(c1_m0, "Lmd", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c1_mxField), &c1_r0.Lmd, 1, 0, 0U, 0, 0U,
                      0);
  c1_mxField = sf_mex_getfield(c1_m0, "rd1", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c1_mxField), &c1_r0.rd1, 1, 0, 0U, 0, 0U,
                      0);
  c1_mxField = sf_mex_getfield(c1_m0, "rd2", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c1_mxField), &c1_r0.rd2, 1, 0, 0U, 0, 0U,
                      0);
  c1_mxField = sf_mex_getfield(c1_m0, "rd3", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c1_mxField), &c1_r0.rd3, 1, 0, 0U, 0, 0U,
                      0);
  c1_mxField = sf_mex_getfield(c1_m0, "rfd", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c1_mxField), &c1_r0.rfd, 1, 0, 0U, 0, 0U,
                      0);
  c1_mxField = sf_mex_getfield(c1_m0, "Lld1", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c1_mxField), &c1_r0.Lld1, 1, 0, 0U, 0, 0U,
                      0);
  c1_mxField = sf_mex_getfield(c1_m0, "Lld2", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c1_mxField), &c1_r0.Lld2, 1, 0, 0U, 0, 0U,
                      0);
  c1_mxField = sf_mex_getfield(c1_m0, "Lld3", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c1_mxField), &c1_r0.Lld3, 1, 0, 0U, 0, 0U,
                      0);
  c1_mxField = sf_mex_getfield(c1_m0, "Llfd", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c1_mxField), &c1_r0.Llfd, 1, 0, 0U, 0, 0U,
                      0);
  c1_mxField = sf_mex_getfield(c1_m0, "rq1", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c1_mxField), &c1_r0.rq1, 1, 0, 0U, 0, 0U,
                      0);
  c1_mxField = sf_mex_getfield(c1_m0, "rq2", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c1_mxField), &c1_r0.rq2, 1, 0, 0U, 0, 0U,
                      0);
  c1_mxField = sf_mex_getfield(c1_m0, "rq3", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c1_mxField), &c1_r0.rq3, 1, 0, 0U, 0, 0U,
                      0);
  c1_mxField = sf_mex_getfield(c1_m0, "pole", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c1_mxField), &c1_r0.pole, 1, 0, 0U, 0, 0U,
                      0);
  c1_mxField = sf_mex_getfield(c1_m0, "Llq1", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c1_mxField), &c1_r0.Llq1, 1, 0, 0U, 0, 0U,
                      0);
  c1_mxField = sf_mex_getfield(c1_m0, "Llq2", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c1_mxField), &c1_r0.Llq2, 1, 0, 0U, 0, 0U,
                      0);
  c1_mxField = sf_mex_getfield(c1_m0, "Llq3", "P", 0);
  sf_mex_import_named("P", sf_mex_dup(c1_mxField), &c1_r0.Llq3, 1, 0, 0U, 0, 0U,
                      0);
  sf_mex_destroy(&c1_m0);
  chartInstance->c1_P = c1_r0;
}

static void enable_c1_gen_rec(SFc1_gen_recInstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c1_gen_rec(SFc1_gen_recInstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c1_update_debugger_state_c1_gen_rec(SFc1_gen_recInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c1_gen_rec(SFc1_gen_recInstanceStruct
  *chartInstance)
{
  const mxArray *c1_st;
  const mxArray *c1_y = NULL;
  real_T c1_hoistedGlobal;
  real_T c1_u;
  const mxArray *c1_b_y = NULL;
  real_T c1_b_hoistedGlobal;
  real_T c1_b_u;
  const mxArray *c1_c_y = NULL;
  real_T c1_c_hoistedGlobal;
  real_T c1_c_u;
  const mxArray *c1_d_y = NULL;
  real_T c1_d_hoistedGlobal;
  real_T c1_d_u;
  const mxArray *c1_e_y = NULL;
  real_T c1_e_hoistedGlobal;
  real_T c1_e_u;
  const mxArray *c1_f_y = NULL;
  uint8_T c1_f_hoistedGlobal;
  uint8_T c1_f_u;
  const mxArray *c1_g_y = NULL;
  real_T *c1_beta;
  real_T *c1_ids;
  real_T *c1_iqs;
  real_T *c1_pidc;
  real_T *c1_g_u;
  c1_g_u = (real_T *)ssGetOutputPortSignal(chartInstance->S, 5);
  c1_beta = (real_T *)ssGetOutputPortSignal(chartInstance->S, 4);
  c1_pidc = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c1_ids = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c1_iqs = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c1_st = NULL;
  c1_st = NULL;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_createcellmatrix(6, 1), false);
  c1_hoistedGlobal = *c1_beta;
  c1_u = c1_hoistedGlobal;
  c1_b_y = NULL;
  sf_mex_assign(&c1_b_y, sf_mex_create("y", &c1_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c1_y, 0, c1_b_y);
  c1_b_hoistedGlobal = *c1_ids;
  c1_b_u = c1_b_hoistedGlobal;
  c1_c_y = NULL;
  sf_mex_assign(&c1_c_y, sf_mex_create("y", &c1_b_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c1_y, 1, c1_c_y);
  c1_c_hoistedGlobal = *c1_iqs;
  c1_c_u = c1_c_hoistedGlobal;
  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", &c1_c_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c1_y, 2, c1_d_y);
  c1_d_hoistedGlobal = *c1_pidc;
  c1_d_u = c1_d_hoistedGlobal;
  c1_e_y = NULL;
  sf_mex_assign(&c1_e_y, sf_mex_create("y", &c1_d_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c1_y, 3, c1_e_y);
  c1_e_hoistedGlobal = *c1_g_u;
  c1_e_u = c1_e_hoistedGlobal;
  c1_f_y = NULL;
  sf_mex_assign(&c1_f_y, sf_mex_create("y", &c1_e_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c1_y, 4, c1_f_y);
  c1_f_hoistedGlobal = chartInstance->c1_is_active_c1_gen_rec;
  c1_f_u = c1_f_hoistedGlobal;
  c1_g_y = NULL;
  sf_mex_assign(&c1_g_y, sf_mex_create("y", &c1_f_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c1_y, 5, c1_g_y);
  sf_mex_assign(&c1_st, c1_y, false);
  return c1_st;
}

static void set_sim_state_c1_gen_rec(SFc1_gen_recInstanceStruct *chartInstance,
  const mxArray *c1_st)
{
  const mxArray *c1_u;
  real_T *c1_beta;
  real_T *c1_ids;
  real_T *c1_iqs;
  real_T *c1_pidc;
  real_T *c1_b_u;
  c1_b_u = (real_T *)ssGetOutputPortSignal(chartInstance->S, 5);
  c1_beta = (real_T *)ssGetOutputPortSignal(chartInstance->S, 4);
  c1_pidc = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c1_ids = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c1_iqs = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c1_doneDoubleBufferReInit = true;
  c1_u = sf_mex_dup(c1_st);
  *c1_beta = c1_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u,
    0)), "beta");
  *c1_ids = c1_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 1)),
    "ids");
  *c1_iqs = c1_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 2)),
    "iqs");
  *c1_pidc = c1_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u,
    3)), "pidc");
  *c1_b_u = c1_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 4)),
    "u");
  chartInstance->c1_is_active_c1_gen_rec = c1_g_emlrt_marshallIn(chartInstance,
    sf_mex_dup(sf_mex_getcell(c1_u, 5)), "is_active_c1_gen_rec");
  sf_mex_destroy(&c1_u);
  c1_update_debugger_state_c1_gen_rec(chartInstance);
  sf_mex_destroy(&c1_st);
}

static void finalize_c1_gen_rec(SFc1_gen_recInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c1_gen_rec(SFc1_gen_recInstanceStruct *chartInstance)
{
  int32_T c1_i0;
  real_T *c1_iqs;
  real_T *c1_wrm;
  real_T *c1_betac;
  real_T *c1_ed;
  real_T *c1_idc;
  real_T *c1_betaminprev;
  real_T *c1_uprev;
  real_T *c1_ids;
  real_T *c1_pidc;
  real_T *c1_beta;
  real_T *c1_u;
  real_T (*c1_Iamr)[7];
  c1_u = (real_T *)ssGetOutputPortSignal(chartInstance->S, 5);
  c1_beta = (real_T *)ssGetOutputPortSignal(chartInstance->S, 4);
  c1_pidc = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c1_ids = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c1_uprev = (real_T *)ssGetInputPortSignal(chartInstance->S, 6);
  c1_betaminprev = (real_T *)ssGetInputPortSignal(chartInstance->S, 5);
  c1_idc = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
  c1_Iamr = (real_T (*)[7])ssGetInputPortSignal(chartInstance->S, 3);
  c1_ed = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c1_betac = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c1_wrm = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
  c1_iqs = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 0U, chartInstance->c1_sfEvent);
  chartInstance->c1_sfEvent = CALL_EVENT;
  c1_chartstep_c1_gen_rec(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_gen_recMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  _SFD_DATA_RANGE_CHECK(*c1_iqs, 0U);
  _SFD_DATA_RANGE_CHECK(*c1_wrm, 1U);
  _SFD_DATA_RANGE_CHECK(*c1_betac, 2U);
  _SFD_DATA_RANGE_CHECK(*c1_ed, 3U);
  for (c1_i0 = 0; c1_i0 < 7; c1_i0++) {
    _SFD_DATA_RANGE_CHECK((*c1_Iamr)[c1_i0], 4U);
  }

  _SFD_DATA_RANGE_CHECK(*c1_idc, 5U);
  _SFD_DATA_RANGE_CHECK(*c1_betaminprev, 6U);
  _SFD_DATA_RANGE_CHECK(*c1_uprev, 7U);
  _SFD_DATA_RANGE_CHECK(*c1_ids, 8U);
  _SFD_DATA_RANGE_CHECK(*c1_pidc, 9U);
  _SFD_DATA_RANGE_CHECK(*c1_beta, 10U);
  _SFD_DATA_RANGE_CHECK(*c1_u, 11U);
}

static void c1_chartstep_c1_gen_rec(SFc1_gen_recInstanceStruct *chartInstance)
{
  real_T c1_hoistedGlobal;
  real_T c1_b_hoistedGlobal;
  real_T c1_c_hoistedGlobal;
  real_T c1_d_hoistedGlobal;
  real_T c1_e_hoistedGlobal;
  real_T c1_f_hoistedGlobal;
  real_T c1_wrm;
  real_T c1_betac;
  real_T c1_ed;
  int32_T c1_i1;
  real_T c1_Iamr[7];
  real_T c1_idc;
  real_T c1_betaminprev;
  real_T c1_uprev;
  c1_struct_U3Zfiu1q58uXf1A0Qkd5GF c1_b_P;
  uint32_T c1_debug_family_var_map[67];
  real_T c1_lamq1;
  real_T c1_lamq2;
  real_T c1_lamq3;
  real_T c1_lamd1;
  real_T c1_lamd2;
  real_T c1_lamd3;
  real_T c1_lamdf;
  real_T c1_wr;
  real_T c1_LdDP;
  real_T c1_LqDP;
  real_T c1_lamdDP;
  real_T c1_lamqDP;
  real_T c1_x;
  real_T c1_epsilon;
  real_T c1_delta;
  real_T c1_max_iter;
  real_T c1_record[100];
  real_T c1_fx;
  real_T c1_n;
  real_T c1_fp;
  real_T c1_d;
  real_T c1_betamin;
  real_T c1_eqDP;
  real_T c1_edDP;
  real_T c1_phi;
  real_T c1_alpha;
  real_T c1_Epk;
  real_T c1_Lc;
  real_T c1_Lt;
  real_T c1_K;
  real_T c1_iqsCond;
  real_T c1_idsCond;
  real_T c1_theta;
  real_T c1_ias1;
  real_T c1_iqs1;
  real_T c1_ids1;
  real_T c1_ias2;
  real_T c1_iqs2;
  real_T c1_ids2;
  real_T c1_ias3;
  real_T c1_iqs3;
  real_T c1_ids3;
  real_T c1_ias4;
  real_T c1_iqs4;
  real_T c1_ids4;
  real_T c1_ias5;
  real_T c1_iqs5;
  real_T c1_ids5;
  real_T c1_iqs_C[5];
  real_T c1_ids_C[5];
  real_T c1_iqsCommu;
  real_T c1_idsCommu;
  real_T c1_nargin = 8.0;
  real_T c1_nargout = 5.0;
  real_T c1_iqs;
  real_T c1_ids;
  real_T c1_pidc;
  real_T c1_beta;
  real_T c1_u;
  real_T c1_A;
  real_T c1_b_x;
  real_T c1_c_x;
  real_T c1_d_x;
  real_T c1_y;
  real_T c1_B;
  real_T c1_b_y;
  real_T c1_c_y;
  real_T c1_d_y;
  real_T c1_e_y;
  real_T c1_b_B;
  real_T c1_f_y;
  real_T c1_g_y;
  real_T c1_h_y;
  real_T c1_i_y;
  real_T c1_c_P[3];
  real_T c1_dv0[3];
  int32_T c1_i2;
  real_T c1_dv1[3];
  real_T c1_c_B;
  real_T c1_j_y;
  real_T c1_k_y;
  real_T c1_l_y;
  real_T c1_m_y;
  real_T c1_d_B;
  real_T c1_n_y;
  real_T c1_o_y;
  real_T c1_p_y;
  real_T c1_q_y;
  real_T c1_d_P[3];
  int32_T c1_i3;
  real_T c1_dv2[3];
  real_T c1_e_B;
  real_T c1_r_y;
  real_T c1_s_y;
  real_T c1_t_y;
  real_T c1_u_y;
  real_T c1_b_A;
  real_T c1_f_B;
  real_T c1_e_x;
  real_T c1_v_y;
  real_T c1_f_x;
  real_T c1_w_y;
  real_T c1_g_x;
  real_T c1_x_y;
  real_T c1_y_y;
  real_T c1_b_lamd1[3];
  real_T c1_e_P[3];
  int32_T c1_i4;
  real_T c1_dv3[3];
  real_T c1_b_lamq1[3];
  real_T c1_f_P[3];
  int32_T c1_i5;
  real_T c1_dv4[3];
  int32_T c1_i6;
  real_T c1_h_x;
  real_T c1_i_x;
  real_T c1_j_x;
  real_T c1_k_x;
  real_T c1_l_x;
  real_T c1_m_x;
  int32_T c1_b_n;
  real_T c1_n_x;
  real_T c1_o_x;
  real_T c1_p_x;
  real_T c1_q_x;
  real_T c1_r_x;
  real_T c1_s_x;
  int32_T c1_i7;
  static char_T c1_cv0[43] = { 'd', 'e', 'r', 'i', 'v', 'a', 't', 'i', 'v', 'e',
    ' ', 't', 'o', 'o', ' ', 's', 'm', 'a', 'l', 'l', ' ', '(', 'N', 'e', 'w',
    't', 'o', 'n', ' ', 'M', 'e', 't', 'h', 'o', 'd', ' ', '-', ' ', 'B', 'e',
    't', 'a', ')' };

  char_T c1_b_u[43];
  const mxArray *c1_ab_y = NULL;
  real_T c1_c_A;
  real_T c1_g_B;
  real_T c1_t_x;
  real_T c1_bb_y;
  real_T c1_u_x;
  real_T c1_cb_y;
  real_T c1_v_x;
  real_T c1_db_y;
  real_T c1_w_x;
  real_T c1_x_x;
  real_T c1_y_x;
  real_T c1_ab_x;
  real_T c1_bb_x;
  real_T c1_cb_x;
  real_T c1_b;
  static creal_T c1_dc0 = { 0.0, 1.0 };

  creal_T c1_eb_y;
  creal_T c1_b_eqDP;
  real_T c1_db_x;
  real_T c1_eb_x;
  real_T c1_fb_x;
  real_T c1_gb_x;
  real_T c1_hb_x;
  real_T c1_ib_x;
  real_T c1_jb_x;
  real_T c1_kb_x;
  real_T c1_lb_x;
  real_T c1_mb_x;
  int32_T c1_i8;
  real_T c1_nb_x;
  real_T c1_ob_x;
  real_T c1_pb_x;
  real_T c1_qb_x;
  real_T c1_rb_x;
  real_T c1_sb_x;
  int32_T c1_c_n;
  real_T c1_tb_x;
  real_T c1_ub_x;
  real_T c1_vb_x;
  real_T c1_wb_x;
  real_T c1_xb_x;
  real_T c1_yb_x;
  int32_T c1_i9;
  static char_T c1_cv1[40] = { 'd', 'e', 'r', 'i', 'v', 'a', 't', 'i', 'v', 'e',
    ' ', 't', 'o', 'o', ' ', 's', 'm', 'a', 'l', 'l', ' ', '(', 'N', 'e', 'w',
    't', 'o', 'n', ' ', 'M', 'e', 't', 'h', 'o', 'd', ' ', '-', ' ', 'U', ')' };

  char_T c1_c_u[40];
  const mxArray *c1_fb_y = NULL;
  real_T c1_d_A;
  real_T c1_h_B;
  real_T c1_ac_x;
  real_T c1_gb_y;
  real_T c1_bc_x;
  real_T c1_hb_y;
  real_T c1_cc_x;
  real_T c1_ib_y;
  real_T c1_dc_x;
  real_T c1_ec_x;
  real_T c1_fc_x;
  real_T c1_gc_x;
  real_T c1_hc_x;
  real_T c1_ic_x;
  real_T c1_jc_x;
  real_T c1_kc_x;
  real_T c1_lc_x;
  real_T c1_mc_x;
  real_T c1_nc_x;
  real_T c1_oc_x;
  real_T c1_pc_x;
  real_T c1_qc_x;
  real_T c1_rc_x;
  real_T c1_sc_x;
  real_T c1_tc_x;
  real_T c1_uc_x;
  real_T c1_vc_x;
  real_T c1_wc_x;
  real_T c1_e_A;
  real_T c1_i_B;
  real_T c1_xc_x;
  real_T c1_jb_y;
  real_T c1_yc_x;
  real_T c1_kb_y;
  real_T c1_ad_x;
  real_T c1_lb_y;
  real_T c1_mb_y;
  real_T c1_bd_x;
  real_T c1_cd_x;
  real_T c1_dd_x;
  real_T c1_ed_x;
  real_T c1_f_A;
  real_T c1_j_B;
  real_T c1_fd_x;
  real_T c1_nb_y;
  real_T c1_gd_x;
  real_T c1_ob_y;
  real_T c1_hd_x;
  real_T c1_pb_y;
  real_T c1_qb_y;
  real_T c1_id_x;
  real_T c1_jd_x;
  real_T c1_kd_x;
  real_T c1_ld_x;
  real_T c1_md_x;
  real_T c1_nd_x;
  real_T c1_od_x;
  real_T c1_pd_x;
  real_T c1_g_A;
  real_T c1_qd_x;
  real_T c1_rd_x;
  real_T c1_sd_x;
  real_T c1_td_x;
  real_T c1_ud_x;
  real_T c1_vd_x;
  real_T c1_wd_x;
  real_T c1_xd_x;
  real_T c1_yd_x;
  real_T c1_h_A;
  real_T c1_k_B;
  real_T c1_ae_x;
  real_T c1_rb_y;
  real_T c1_be_x;
  real_T c1_sb_y;
  real_T c1_ce_x;
  real_T c1_tb_y;
  real_T c1_ub_y;
  real_T c1_de_x;
  real_T c1_ee_x;
  real_T c1_fe_x;
  real_T c1_ge_x;
  real_T c1_i_A;
  real_T c1_l_B;
  real_T c1_he_x;
  real_T c1_vb_y;
  real_T c1_ie_x;
  real_T c1_wb_y;
  real_T c1_je_x;
  real_T c1_xb_y;
  real_T c1_yb_y;
  real_T c1_ke_x;
  real_T c1_le_x;
  real_T c1_me_x;
  real_T c1_ne_x;
  real_T c1_oe_x;
  real_T c1_pe_x;
  real_T c1_qe_x;
  real_T c1_re_x;
  real_T c1_j_A;
  real_T c1_se_x;
  real_T c1_te_x;
  real_T c1_ue_x;
  real_T c1_ve_x;
  real_T c1_we_x;
  real_T c1_xe_x;
  real_T c1_ye_x;
  real_T c1_af_x;
  real_T c1_bf_x;
  real_T c1_k_A;
  real_T c1_m_B;
  real_T c1_cf_x;
  real_T c1_ac_y;
  real_T c1_df_x;
  real_T c1_bc_y;
  real_T c1_ef_x;
  real_T c1_cc_y;
  real_T c1_dc_y;
  real_T c1_ff_x;
  real_T c1_gf_x;
  real_T c1_hf_x;
  real_T c1_if_x;
  real_T c1_l_A;
  real_T c1_n_B;
  real_T c1_jf_x;
  real_T c1_ec_y;
  real_T c1_kf_x;
  real_T c1_fc_y;
  real_T c1_lf_x;
  real_T c1_gc_y;
  real_T c1_hc_y;
  real_T c1_mf_x;
  real_T c1_nf_x;
  real_T c1_of_x;
  real_T c1_pf_x;
  real_T c1_qf_x;
  real_T c1_rf_x;
  real_T c1_sf_x;
  real_T c1_tf_x;
  real_T c1_m_A;
  real_T c1_uf_x;
  real_T c1_vf_x;
  real_T c1_wf_x;
  real_T c1_xf_x;
  real_T c1_yf_x;
  real_T c1_ag_x;
  real_T c1_bg_x;
  real_T c1_cg_x;
  real_T c1_dg_x;
  real_T c1_n_A;
  real_T c1_o_B;
  real_T c1_eg_x;
  real_T c1_ic_y;
  real_T c1_fg_x;
  real_T c1_jc_y;
  real_T c1_gg_x;
  real_T c1_kc_y;
  real_T c1_lc_y;
  real_T c1_hg_x;
  real_T c1_ig_x;
  real_T c1_jg_x;
  real_T c1_kg_x;
  real_T c1_o_A;
  real_T c1_p_B;
  real_T c1_lg_x;
  real_T c1_mc_y;
  real_T c1_mg_x;
  real_T c1_nc_y;
  real_T c1_ng_x;
  real_T c1_oc_y;
  real_T c1_pc_y;
  real_T c1_og_x;
  real_T c1_pg_x;
  real_T c1_qg_x;
  real_T c1_rg_x;
  real_T c1_sg_x;
  real_T c1_tg_x;
  real_T c1_ug_x;
  real_T c1_vg_x;
  real_T c1_wg_x;
  real_T c1_xg_x;
  real_T c1_yg_x;
  real_T c1_ah_x;
  real_T c1_bh_x;
  real_T c1_ch_x;
  real_T c1_p_A;
  real_T c1_q_B;
  real_T c1_dh_x;
  real_T c1_qc_y;
  real_T c1_eh_x;
  real_T c1_rc_y;
  real_T c1_fh_x;
  real_T c1_sc_y;
  real_T c1_tc_y;
  real_T c1_gh_x;
  real_T c1_hh_x;
  real_T c1_ih_x;
  real_T c1_jh_x;
  real_T c1_q_A;
  real_T c1_r_B;
  real_T c1_kh_x;
  real_T c1_uc_y;
  real_T c1_lh_x;
  real_T c1_vc_y;
  real_T c1_mh_x;
  real_T c1_wc_y;
  real_T c1_xc_y;
  real_T c1_nh_x;
  real_T c1_oh_x;
  real_T c1_ph_x;
  real_T c1_qh_x;
  real_T c1_rh_x;
  real_T c1_sh_x;
  real_T c1_th_x;
  real_T c1_uh_x;
  real_T c1_b_iqs1[5];
  int32_T c1_i10;
  real_T c1_b_ids1[5];
  int32_T c1_i11;
  real_T c1_r_A;
  real_T c1_vh_x;
  real_T c1_wh_x;
  real_T c1_xh_x;
  real_T c1_yc_y;
  real_T c1_s_A;
  real_T c1_yh_x;
  real_T c1_ai_x;
  real_T c1_bi_x;
  real_T c1_ad_y;
  real_T c1_ci_x;
  real_T c1_di_x;
  real_T c1_t_A;
  real_T c1_s_B;
  real_T c1_ei_x;
  real_T c1_bd_y;
  real_T c1_fi_x;
  real_T c1_cd_y;
  real_T c1_gi_x;
  real_T c1_dd_y;
  real_T *c1_d_u;
  real_T *c1_b_beta;
  real_T *c1_b_pidc;
  real_T *c1_b_ids;
  real_T *c1_b_iqs;
  real_T *c1_b_uprev;
  real_T *c1_b_betaminprev;
  real_T *c1_b_idc;
  real_T *c1_b_ed;
  real_T *c1_b_betac;
  real_T *c1_b_wrm;
  real_T (*c1_b_Iamr)[7];
  int32_T exitg1;
  int32_T exitg2;
  c1_d_u = (real_T *)ssGetOutputPortSignal(chartInstance->S, 5);
  c1_b_beta = (real_T *)ssGetOutputPortSignal(chartInstance->S, 4);
  c1_b_pidc = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c1_b_ids = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c1_b_uprev = (real_T *)ssGetInputPortSignal(chartInstance->S, 6);
  c1_b_betaminprev = (real_T *)ssGetInputPortSignal(chartInstance->S, 5);
  c1_b_idc = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
  c1_b_Iamr = (real_T (*)[7])ssGetInputPortSignal(chartInstance->S, 3);
  c1_b_ed = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c1_b_betac = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c1_b_wrm = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
  c1_b_iqs = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 0U, chartInstance->c1_sfEvent);
  c1_hoistedGlobal = *c1_b_wrm;
  c1_b_hoistedGlobal = *c1_b_betac;
  c1_c_hoistedGlobal = *c1_b_ed;
  c1_d_hoistedGlobal = *c1_b_idc;
  c1_e_hoistedGlobal = *c1_b_betaminprev;
  c1_f_hoistedGlobal = *c1_b_uprev;
  c1_wrm = c1_hoistedGlobal;
  c1_betac = c1_b_hoistedGlobal;
  c1_ed = c1_c_hoistedGlobal;
  for (c1_i1 = 0; c1_i1 < 7; c1_i1++) {
    c1_Iamr[c1_i1] = (*c1_b_Iamr)[c1_i1];
  }

  c1_idc = c1_d_hoistedGlobal;
  c1_betaminprev = c1_e_hoistedGlobal;
  c1_uprev = c1_f_hoistedGlobal;
  c1_b_P = chartInstance->c1_P;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 67U, 67U, c1_debug_family_names,
    c1_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_lamq1, 0U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_lamq2, 1U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_lamq3, 2U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_lamd1, 3U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_lamd2, 4U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_lamd3, 5U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_lamdf, 6U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_wr, 7U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_LdDP, 8U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_LqDP, 9U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_lamdDP, 10U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_lamqDP, 11U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_x, 12U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_epsilon, 13U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_delta, 14U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c1_max_iter, 15U, c1_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_record, 16U, c1_e_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_fx, 17U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_n, 18U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_fp, 19U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_d, 20U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_betamin, 21U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_eqDP, 22U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_edDP, 23U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_phi, 24U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_alpha, 25U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_Epk, 26U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_Lc, 27U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_Lt, 28U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_K, 29U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_iqsCond, 30U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_idsCond, 31U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c1_theta, 32U, c1_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_ias1, 33U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_iqs1, 34U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_ids1, 35U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_ias2, 36U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_iqs2, 37U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_ids2, 38U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_ias3, 39U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_iqs3, 40U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_ids3, 41U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_ias4, 42U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_iqs4, 43U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_ids4, 44U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_ias5, 45U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_iqs5, 46U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_ids5, 47U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_iqs_C, 48U, c1_d_sf_marshallOut,
    c1_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_ids_C, 49U, c1_d_sf_marshallOut,
    c1_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_iqsCommu, 50U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_idsCommu, 51U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargin, 52U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargout, 53U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c1_wrm, 54U, c1_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c1_betac, 55U, c1_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c1_ed, 56U, c1_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_Iamr, 57U, c1_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c1_idc, 58U, c1_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c1_betaminprev, 59U, c1_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c1_uprev, 60U, c1_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_b_P, 61U, c1_b_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_iqs, 62U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_ids, 63U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_pidc, 64U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_beta, 65U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_u, 66U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 43);
  c1_lamq1 = c1_Iamr[0];
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 44);
  c1_lamq2 = c1_Iamr[1];
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 45);
  c1_lamq3 = c1_Iamr[2];
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 46);
  c1_lamd1 = c1_Iamr[3];
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 47);
  c1_lamd2 = c1_Iamr[4];
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 48);
  c1_lamd3 = c1_Iamr[5];
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 49);
  c1_lamdf = c1_Iamr[6];
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 52);
  c1_A = c1_b_P.pole;
  c1_b_x = c1_A;
  c1_c_x = c1_b_x;
  c1_d_x = c1_c_x;
  c1_y = c1_d_x / 2.0;
  c1_wr = c1_y * c1_wrm;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 56);
  c1_B = c1_b_P.Lmd;
  c1_b_y = c1_B;
  c1_c_y = c1_b_y;
  c1_d_y = c1_c_y;
  c1_e_y = 1.0 / c1_d_y;
  c1_b_B = c1_b_P.Llfd;
  c1_f_y = c1_b_B;
  c1_g_y = c1_f_y;
  c1_h_y = c1_g_y;
  c1_i_y = 1.0 / c1_h_y;
  c1_c_P[0] = c1_b_P.Lld1;
  c1_c_P[1] = c1_b_P.Lld2;
  c1_c_P[2] = c1_b_P.Lld3;
  c1_rdivide(chartInstance, 1.0, c1_c_P, c1_dv0);
  for (c1_i2 = 0; c1_i2 < 3; c1_i2++) {
    c1_dv1[c1_i2] = c1_dv0[c1_i2];
  }

  c1_c_B = (c1_e_y + c1_i_y) + c1_sum(chartInstance, c1_dv1);
  c1_j_y = c1_c_B;
  c1_k_y = c1_j_y;
  c1_l_y = c1_k_y;
  c1_m_y = 1.0 / c1_l_y;
  c1_LdDP = c1_b_P.Lls + c1_m_y;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 57);
  c1_d_B = c1_b_P.Lmq;
  c1_n_y = c1_d_B;
  c1_o_y = c1_n_y;
  c1_p_y = c1_o_y;
  c1_q_y = 1.0 / c1_p_y;
  c1_d_P[0] = c1_b_P.Lld1;
  c1_d_P[1] = c1_b_P.Lld2;
  c1_d_P[2] = c1_b_P.Lld3;
  c1_rdivide(chartInstance, 1.0, c1_d_P, c1_dv0);
  for (c1_i3 = 0; c1_i3 < 3; c1_i3++) {
    c1_dv2[c1_i3] = c1_dv0[c1_i3];
  }

  c1_e_B = c1_q_y + c1_sum(chartInstance, c1_dv2);
  c1_r_y = c1_e_B;
  c1_s_y = c1_r_y;
  c1_t_y = c1_s_y;
  c1_u_y = 1.0 / c1_t_y;
  c1_LqDP = c1_b_P.Lls + c1_u_y;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 60);
  c1_b_A = c1_lamdf;
  c1_f_B = c1_b_P.Llfd;
  c1_e_x = c1_b_A;
  c1_v_y = c1_f_B;
  c1_f_x = c1_e_x;
  c1_w_y = c1_v_y;
  c1_g_x = c1_f_x;
  c1_x_y = c1_w_y;
  c1_y_y = c1_g_x / c1_x_y;
  c1_b_lamd1[0] = c1_lamd1;
  c1_b_lamd1[1] = c1_lamd2;
  c1_b_lamd1[2] = c1_lamd3;
  c1_e_P[0] = c1_b_P.Lld1;
  c1_e_P[1] = c1_b_P.Lld2;
  c1_e_P[2] = c1_b_P.Lld3;
  c1_b_rdivide(chartInstance, c1_b_lamd1, c1_e_P, c1_dv0);
  for (c1_i4 = 0; c1_i4 < 3; c1_i4++) {
    c1_dv3[c1_i4] = c1_dv0[c1_i4];
  }

  c1_lamdDP = (c1_LdDP - c1_b_P.Lls) * (c1_y_y + c1_sum(chartInstance, c1_dv3));
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 61);
  c1_b_lamq1[0] = c1_lamq1;
  c1_b_lamq1[1] = c1_lamq2;
  c1_b_lamq1[2] = c1_lamq3;
  c1_f_P[0] = c1_b_P.Llq1;
  c1_f_P[1] = c1_b_P.Llq2;
  c1_f_P[2] = c1_b_P.Llq3;
  c1_b_rdivide(chartInstance, c1_b_lamq1, c1_f_P, c1_dv0);
  for (c1_i5 = 0; c1_i5 < 3; c1_i5++) {
    c1_dv4[c1_i5] = c1_dv0[c1_i5];
  }

  c1_lamqDP = (c1_LqDP - c1_b_P.Lls) * c1_sum(chartInstance, c1_dv4);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 65);
  c1_x = c1_betaminprev;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 66);
  c1_epsilon = 1.0E-5;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 67);
  c1_delta = 1.0E-10;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 68);
  c1_max_iter = 100.0;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 69);
  for (c1_i6 = 0; c1_i6 < 100; c1_i6++) {
    c1_record[c1_i6] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 71);
  c1_h_x = c1_x;
  c1_i_x = c1_h_x;
  c1_i_x = muDoubleScalarCos(c1_i_x);
  c1_j_x = c1_x;
  c1_k_x = c1_j_x;
  c1_k_x = muDoubleScalarSin(c1_k_x);
  c1_l_x = 2.0 * c1_x - 1.0471975511965976;
  c1_m_x = c1_l_x;
  c1_m_x = muDoubleScalarSin(c1_m_x);
  c1_fx = (1.7320508075688772 * c1_lamqDP * c1_i_x + 1.7320508075688772 *
           c1_lamdDP * c1_k_x) + 2.0 * c1_idc * (c1_LqDP - c1_LdDP) * c1_m_x;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 73);
  c1_n = 1.0;
  c1_b_n = 0;
  do {
    exitg2 = 0;
    if (c1_b_n < 100) {
      c1_n = 1.0 + (real_T)c1_b_n;
      CV_EML_FOR(0, 1, 0, 1);
      _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 75);
      c1_n_x = c1_x;
      c1_o_x = c1_n_x;
      c1_o_x = muDoubleScalarSin(c1_o_x);
      c1_p_x = c1_x;
      c1_q_x = c1_p_x;
      c1_q_x = muDoubleScalarCos(c1_q_x);
      c1_r_x = 2.0 * c1_x - 1.0471975511965976;
      c1_s_x = c1_r_x;
      c1_s_x = muDoubleScalarCos(c1_s_x);
      c1_fp = (-1.7320508075688772 * c1_lamqDP * c1_o_x + 1.7320508075688772 *
               c1_lamdDP * c1_q_x) + 4.0 * c1_idc * (c1_LqDP - c1_LdDP) * c1_s_x;
      _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 76);
      if (CV_EML_IF(0, 1, 0, c1_fp < c1_delta)) {
        _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 77);
        for (c1_i7 = 0; c1_i7 < 43; c1_i7++) {
          c1_b_u[c1_i7] = c1_cv0[c1_i7];
        }

        c1_ab_y = NULL;
        sf_mex_assign(&c1_ab_y, sf_mex_create("y", c1_b_u, 10, 0U, 1U, 0U, 2, 1,
          43), false);
        sf_mex_call_debug(sfGlobalDebugInstanceStruct, "disp", 0U, 1U, 14,
                          c1_ab_y);
        _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 78);
        exitg2 = 1;
      } else {
        _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 81);
        c1_c_A = c1_fx;
        c1_g_B = c1_fp;
        c1_t_x = c1_c_A;
        c1_bb_y = c1_g_B;
        c1_u_x = c1_t_x;
        c1_cb_y = c1_bb_y;
        c1_v_x = c1_u_x;
        c1_db_y = c1_cb_y;
        c1_d = c1_v_x / c1_db_y;
        _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 82);
        c1_x -= c1_d;
        _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 84);
        c1_w_x = c1_x;
        c1_x_x = c1_w_x;
        c1_x_x = muDoubleScalarCos(c1_x_x);
        c1_y_x = c1_x;
        c1_ab_x = c1_y_x;
        c1_ab_x = muDoubleScalarSin(c1_ab_x);
        c1_bb_x = 2.0 * c1_x - 1.0471975511965976;
        c1_cb_x = c1_bb_x;
        c1_cb_x = muDoubleScalarSin(c1_cb_x);
        c1_fx = (1.7320508075688772 * c1_lamqDP * c1_x_x + 1.7320508075688772 *
                 c1_lamdDP * c1_ab_x) + 2.0 * c1_idc * (c1_LqDP - c1_LdDP) *
          c1_cb_x;
        _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 85);
        if (CV_EML_IF(0, 1, 1, c1_abs(chartInstance, c1_d) < c1_epsilon)) {
          _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 87);
          exitg2 = 1;
        } else {
          _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 89);
          c1_record[_SFD_EML_ARRAY_BOUNDS_CHECK("record", (int32_T)
            _SFD_INTEGER_CHECK("n", c1_n), 1, 100, 1, 0) - 1] = c1_x;
          c1_b_n++;
          _SF_MEX_LISTEN_FOR_CTRL_C(chartInstance->S);
        }
      }
    } else {
      CV_EML_FOR(0, 1, 0, 0);
      exitg2 = 1;
    }
  } while (exitg2 == 0);

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 94);
  c1_betamin = c1_x;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 95);
  c1_beta = c1_betamin;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 101);
  c1_eqDP = c1_wr * c1_lamdDP;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 102);
  c1_edDP = -c1_wr * c1_lamqDP;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 103);
  c1_b = c1_edDP;
  c1_eb_y.re = c1_b * c1_dc0.re;
  c1_eb_y.im = c1_b * c1_dc0.im;
  c1_b_eqDP.re = c1_eqDP - c1_eb_y.re;
  c1_b_eqDP.im = 0.0 - c1_eb_y.im;
  c1_phi = c1_angle(chartInstance, c1_b_eqDP);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 104);
  c1_alpha = c1_beta + c1_phi;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 105);
  c1_Epk = c1_mpower(chartInstance, c1_eqDP) + c1_mpower(chartInstance, c1_edDP);
  c1_b_sqrt(chartInstance, &c1_Epk);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 106);
  c1_db_x = 2.0 * c1_beta + 0.52359877559829882;
  c1_eb_x = c1_db_x;
  c1_eb_x = muDoubleScalarSin(c1_eb_x);
  c1_Lc = 0.5 * (c1_LqDP + c1_LdDP) + (c1_LdDP - c1_LqDP) * c1_eb_x;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 107);
  c1_fb_x = 2.0 * c1_beta - 0.52359877559829882;
  c1_gb_x = c1_fb_x;
  c1_gb_x = muDoubleScalarSin(c1_gb_x);
  c1_Lt = (c1_LqDP + c1_LdDP) + (c1_LdDP - c1_LqDP) * c1_gb_x;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 108);
  c1_hb_x = c1_beta;
  c1_ib_x = c1_hb_x;
  c1_ib_x = muDoubleScalarSin(c1_ib_x);
  c1_jb_x = c1_beta;
  c1_kb_x = c1_jb_x;
  c1_kb_x = muDoubleScalarCos(c1_kb_x);
  c1_lb_x = 2.0 * c1_beta + 2.0943951023931953;
  c1_mb_x = c1_lb_x;
  c1_mb_x = muDoubleScalarCos(c1_mb_x);
  c1_K = 1.7320508075688772 * (-c1_lamqDP * c1_ib_x + c1_lamqDP * c1_kb_x) +
    ((c1_LdDP - c1_LqDP) * c1_mb_x - 0.5 * (c1_LqDP + c1_LdDP)) * c1_idc;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 113);
  c1_x = c1_uprev;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 114);
  c1_epsilon = 1.0E-5;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 115);
  c1_delta = 1.0E-10;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 116);
  c1_max_iter = 100.0;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 117);
  for (c1_i8 = 0; c1_i8 < 100; c1_i8++) {
    c1_record[c1_i8] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 119);
  c1_nb_x = c1_beta + c1_x;
  c1_ob_x = c1_nb_x;
  c1_ob_x = muDoubleScalarSin(c1_ob_x);
  c1_pb_x = c1_beta + c1_x;
  c1_qb_x = c1_pb_x;
  c1_qb_x = muDoubleScalarCos(c1_qb_x);
  c1_rb_x = (2.0 * c1_beta + 2.0 * c1_x) - 2.0943951023931953;
  c1_sb_x = c1_rb_x;
  c1_sb_x = muDoubleScalarSin(c1_sb_x);
  c1_fx = (c1_K - 1.7320508075688772 * (-c1_lamqDP * c1_ob_x + c1_lamdDP *
            c1_qb_x)) - ((c1_LqDP - c1_LdDP) * c1_sb_x + 0.5 * (c1_LqDP +
    c1_LdDP)) * c1_idc;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 121);
  c1_n = 1.0;
  c1_c_n = 0;
  do {
    exitg1 = 0;
    if (c1_c_n < 100) {
      c1_n = 1.0 + (real_T)c1_c_n;
      CV_EML_FOR(0, 1, 1, 1);
      _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 123);
      c1_tb_x = c1_beta + c1_x;
      c1_ub_x = c1_tb_x;
      c1_ub_x = muDoubleScalarCos(c1_ub_x);
      c1_vb_x = c1_beta + c1_x;
      c1_wb_x = c1_vb_x;
      c1_wb_x = muDoubleScalarSin(c1_wb_x);
      c1_xb_x = (2.0 * c1_beta + 2.0 * c1_x) - 2.0943951023931953;
      c1_yb_x = c1_xb_x;
      c1_yb_x = muDoubleScalarCos(c1_yb_x);
      c1_fp = (1.7320508075688772 * c1_lamqDP * c1_ub_x + 1.7320508075688772 *
               c1_lamdDP * c1_wb_x) - 2.0 * (c1_LqDP - c1_LdDP) * c1_yb_x *
        c1_idc;
      _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 124);
      if (CV_EML_IF(0, 1, 2, c1_fp < c1_delta)) {
        _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 125);
        for (c1_i9 = 0; c1_i9 < 40; c1_i9++) {
          c1_c_u[c1_i9] = c1_cv1[c1_i9];
        }

        c1_fb_y = NULL;
        sf_mex_assign(&c1_fb_y, sf_mex_create("y", c1_c_u, 10, 0U, 1U, 0U, 2, 1,
          40), false);
        sf_mex_call_debug(sfGlobalDebugInstanceStruct, "disp", 0U, 1U, 14,
                          c1_fb_y);
        _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 126);
        exitg1 = 1;
      } else {
        _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 129U);
        c1_d_A = c1_fx;
        c1_h_B = c1_fp;
        c1_ac_x = c1_d_A;
        c1_gb_y = c1_h_B;
        c1_bc_x = c1_ac_x;
        c1_hb_y = c1_gb_y;
        c1_cc_x = c1_bc_x;
        c1_ib_y = c1_hb_y;
        c1_d = c1_cc_x / c1_ib_y;
        _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 130U);
        c1_x -= c1_d;
        _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 132U);
        c1_dc_x = c1_beta + c1_x;
        c1_ec_x = c1_dc_x;
        c1_ec_x = muDoubleScalarSin(c1_ec_x);
        c1_fc_x = c1_beta + c1_x;
        c1_gc_x = c1_fc_x;
        c1_gc_x = muDoubleScalarCos(c1_gc_x);
        c1_hc_x = (2.0 * c1_beta + 2.0 * c1_x) - 2.0943951023931953;
        c1_ic_x = c1_hc_x;
        c1_ic_x = muDoubleScalarSin(c1_ic_x);
        c1_fx = (c1_K - 1.7320508075688772 * (-c1_lamqDP * c1_ec_x + c1_lamdDP *
                  c1_gc_x)) - ((c1_LqDP - c1_LdDP) * c1_ic_x + 0.5 * (c1_LqDP +
          c1_LdDP)) * c1_idc;
        _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 133U);
        if (CV_EML_IF(0, 1, 3, c1_abs(chartInstance, c1_d) < c1_epsilon)) {
          _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 135U);
          exitg1 = 1;
        } else {
          _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 137U);
          c1_record[_SFD_EML_ARRAY_BOUNDS_CHECK("record", (int32_T)
            _SFD_INTEGER_CHECK("n", c1_n), 1, 100, 1, 0) - 1] = c1_x;
          c1_c_n++;
          _SF_MEX_LISTEN_FOR_CTRL_C(chartInstance->S);
        }
      }
    } else {
      CV_EML_FOR(0, 1, 1, 0);
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 141U);
  c1_u = c1_x;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 146U);
  c1_jc_x = c1_beta + 2.0943951023931953;
  c1_kc_x = c1_jc_x;
  c1_kc_x = muDoubleScalarCos(c1_kc_x);
  c1_lc_x = (c1_beta + c1_u) + 1.0471975511965976;
  c1_mc_x = c1_lc_x;
  c1_mc_x = muDoubleScalarCos(c1_mc_x);
  c1_iqsCond = 1.1026577908435842 * c1_idc * (c1_kc_x - c1_mc_x);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 147U);
  c1_nc_x = c1_beta + 2.0943951023931953;
  c1_oc_x = c1_nc_x;
  c1_oc_x = muDoubleScalarSin(c1_oc_x);
  c1_pc_x = (c1_beta + c1_u) + 1.0471975511965976;
  c1_qc_x = c1_pc_x;
  c1_qc_x = muDoubleScalarSin(c1_qc_x);
  c1_idsCond = 1.1026577908435842 * c1_idc * (c1_oc_x - c1_qc_x);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 150U);
  c1_theta = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 151U);
  c1_rc_x = c1_theta + c1_beta;
  c1_sc_x = c1_rc_x;
  c1_sc_x = muDoubleScalarSin(c1_sc_x);
  c1_tc_x = c1_theta + c1_beta;
  c1_uc_x = c1_tc_x;
  c1_uc_x = muDoubleScalarCos(c1_uc_x);
  c1_vc_x = 2.0 * c1_beta;
  c1_wc_x = c1_vc_x;
  c1_wc_x = muDoubleScalarCos(c1_wc_x);
  c1_e_A = 1.7320508075688772 * (c1_lamqDP * c1_sc_x - c1_lamdDP * c1_uc_x);
  c1_i_B = (c1_LqDP + c1_LdDP) - (c1_LqDP - c1_LdDP) * c1_wc_x;
  c1_xc_x = c1_e_A;
  c1_jb_y = c1_i_B;
  c1_yc_x = c1_xc_x;
  c1_kb_y = c1_jb_y;
  c1_ad_x = c1_yc_x;
  c1_lb_y = c1_kb_y;
  c1_mb_y = c1_ad_x / c1_lb_y;
  c1_bd_x = 2.0 * c1_beta - 2.0943951023931953;
  c1_cd_x = c1_bd_x;
  c1_cd_x = muDoubleScalarCos(c1_cd_x);
  c1_dd_x = 2.0 * c1_beta;
  c1_ed_x = c1_dd_x;
  c1_ed_x = muDoubleScalarCos(c1_ed_x);
  c1_f_A = ((c1_LqDP - c1_LdDP) * c1_cd_x + 0.5 * (c1_LqDP + c1_LdDP)) * c1_idc;
  c1_j_B = (c1_LqDP + c1_LdDP) - (c1_LqDP - c1_LdDP) * c1_ed_x;
  c1_fd_x = c1_f_A;
  c1_nb_y = c1_j_B;
  c1_gd_x = c1_fd_x;
  c1_ob_y = c1_nb_y;
  c1_hd_x = c1_gd_x;
  c1_pb_y = c1_ob_y;
  c1_qb_y = c1_hd_x / c1_pb_y;
  c1_ias1 = (c1_K + c1_mb_y) - c1_qb_y;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 152U);
  c1_id_x = c1_theta + c1_beta;
  c1_jd_x = c1_id_x;
  c1_jd_x = muDoubleScalarSin(c1_jd_x);
  c1_kd_x = (c1_theta + c1_beta) + 1.0471975511965976;
  c1_ld_x = c1_kd_x;
  c1_ld_x = muDoubleScalarSin(c1_ld_x);
  c1_iqs1 = 1.1547005383792515 * (-c1_ias1 * c1_jd_x - c1_idc * c1_ld_x);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 153U);
  c1_md_x = c1_theta + c1_beta;
  c1_nd_x = c1_md_x;
  c1_nd_x = muDoubleScalarCos(c1_nd_x);
  c1_od_x = (c1_theta + c1_beta) + 1.0471975511965976;
  c1_pd_x = c1_od_x;
  c1_pd_x = muDoubleScalarCos(c1_pd_x);
  c1_ids1 = 1.1547005383792515 * (c1_ias1 * c1_nd_x + c1_idc * c1_pd_x);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 155U);
  c1_g_A = c1_u;
  c1_qd_x = c1_g_A;
  c1_rd_x = c1_qd_x;
  c1_sd_x = c1_rd_x;
  c1_theta = c1_sd_x / 4.0;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 156U);
  c1_td_x = c1_theta + c1_beta;
  c1_ud_x = c1_td_x;
  c1_ud_x = muDoubleScalarSin(c1_ud_x);
  c1_vd_x = c1_theta + c1_beta;
  c1_wd_x = c1_vd_x;
  c1_wd_x = muDoubleScalarCos(c1_wd_x);
  c1_xd_x = 2.0 * c1_theta + 2.0 * c1_beta;
  c1_yd_x = c1_xd_x;
  c1_yd_x = muDoubleScalarCos(c1_yd_x);
  c1_h_A = 1.7320508075688772 * (c1_lamqDP * c1_ud_x - c1_lamdDP * c1_wd_x);
  c1_k_B = (c1_LqDP + c1_LdDP) - (c1_LqDP - c1_LdDP) * c1_yd_x;
  c1_ae_x = c1_h_A;
  c1_rb_y = c1_k_B;
  c1_be_x = c1_ae_x;
  c1_sb_y = c1_rb_y;
  c1_ce_x = c1_be_x;
  c1_tb_y = c1_sb_y;
  c1_ub_y = c1_ce_x / c1_tb_y;
  c1_de_x = (2.0 * c1_theta + 2.0 * c1_beta) - 2.0943951023931953;
  c1_ee_x = c1_de_x;
  c1_ee_x = muDoubleScalarCos(c1_ee_x);
  c1_fe_x = 2.0 * c1_theta + 2.0 * c1_beta;
  c1_ge_x = c1_fe_x;
  c1_ge_x = muDoubleScalarCos(c1_ge_x);
  c1_i_A = ((c1_LqDP - c1_LdDP) * c1_ee_x + 0.5 * (c1_LqDP + c1_LdDP)) * c1_idc;
  c1_l_B = (c1_LqDP + c1_LdDP) - (c1_LqDP - c1_LdDP) * c1_ge_x;
  c1_he_x = c1_i_A;
  c1_vb_y = c1_l_B;
  c1_ie_x = c1_he_x;
  c1_wb_y = c1_vb_y;
  c1_je_x = c1_ie_x;
  c1_xb_y = c1_wb_y;
  c1_yb_y = c1_je_x / c1_xb_y;
  c1_ias2 = (c1_K + c1_ub_y) - c1_yb_y;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 157U);
  c1_ke_x = c1_theta + c1_beta;
  c1_le_x = c1_ke_x;
  c1_le_x = muDoubleScalarSin(c1_le_x);
  c1_me_x = (c1_theta + c1_beta) + 1.0471975511965976;
  c1_ne_x = c1_me_x;
  c1_ne_x = muDoubleScalarSin(c1_ne_x);
  c1_iqs2 = 1.1547005383792515 * (-c1_ias2 * c1_le_x - c1_idc * c1_ne_x);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 158U);
  c1_oe_x = c1_theta + c1_beta;
  c1_pe_x = c1_oe_x;
  c1_pe_x = muDoubleScalarCos(c1_pe_x);
  c1_qe_x = (c1_theta + c1_beta) + 1.0471975511965976;
  c1_re_x = c1_qe_x;
  c1_re_x = muDoubleScalarCos(c1_re_x);
  c1_ids2 = 1.1547005383792515 * (c1_ias2 * c1_pe_x + c1_idc * c1_re_x);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 160U);
  c1_j_A = c1_u;
  c1_se_x = c1_j_A;
  c1_te_x = c1_se_x;
  c1_ue_x = c1_te_x;
  c1_theta = c1_ue_x / 2.0;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 161U);
  c1_ve_x = c1_theta + c1_beta;
  c1_we_x = c1_ve_x;
  c1_we_x = muDoubleScalarSin(c1_we_x);
  c1_xe_x = c1_theta + c1_beta;
  c1_ye_x = c1_xe_x;
  c1_ye_x = muDoubleScalarCos(c1_ye_x);
  c1_af_x = 2.0 * c1_theta + 2.0 * c1_beta;
  c1_bf_x = c1_af_x;
  c1_bf_x = muDoubleScalarCos(c1_bf_x);
  c1_k_A = 1.7320508075688772 * (c1_lamqDP * c1_we_x - c1_lamdDP * c1_ye_x);
  c1_m_B = (c1_LqDP + c1_LdDP) - (c1_LqDP - c1_LdDP) * c1_bf_x;
  c1_cf_x = c1_k_A;
  c1_ac_y = c1_m_B;
  c1_df_x = c1_cf_x;
  c1_bc_y = c1_ac_y;
  c1_ef_x = c1_df_x;
  c1_cc_y = c1_bc_y;
  c1_dc_y = c1_ef_x / c1_cc_y;
  c1_ff_x = (2.0 * c1_theta + 2.0 * c1_beta) - 2.0943951023931953;
  c1_gf_x = c1_ff_x;
  c1_gf_x = muDoubleScalarCos(c1_gf_x);
  c1_hf_x = 2.0 * c1_theta + 2.0 * c1_beta;
  c1_if_x = c1_hf_x;
  c1_if_x = muDoubleScalarCos(c1_if_x);
  c1_l_A = ((c1_LqDP - c1_LdDP) * c1_gf_x + 0.5 * (c1_LqDP + c1_LdDP)) * c1_idc;
  c1_n_B = (c1_LqDP + c1_LdDP) - (c1_LqDP - c1_LdDP) * c1_if_x;
  c1_jf_x = c1_l_A;
  c1_ec_y = c1_n_B;
  c1_kf_x = c1_jf_x;
  c1_fc_y = c1_ec_y;
  c1_lf_x = c1_kf_x;
  c1_gc_y = c1_fc_y;
  c1_hc_y = c1_lf_x / c1_gc_y;
  c1_ias3 = (c1_K + c1_dc_y) - c1_hc_y;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 162U);
  c1_mf_x = c1_theta + c1_beta;
  c1_nf_x = c1_mf_x;
  c1_nf_x = muDoubleScalarSin(c1_nf_x);
  c1_of_x = (c1_theta + c1_beta) + 1.0471975511965976;
  c1_pf_x = c1_of_x;
  c1_pf_x = muDoubleScalarSin(c1_pf_x);
  c1_iqs3 = 1.1547005383792515 * (-c1_ias3 * c1_nf_x - c1_idc * c1_pf_x);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 163U);
  c1_qf_x = c1_theta + c1_beta;
  c1_rf_x = c1_qf_x;
  c1_rf_x = muDoubleScalarCos(c1_rf_x);
  c1_sf_x = (c1_theta + c1_beta) + 1.0471975511965976;
  c1_tf_x = c1_sf_x;
  c1_tf_x = muDoubleScalarCos(c1_tf_x);
  c1_ids3 = 1.1547005383792515 * (c1_ias3 * c1_rf_x + c1_idc * c1_tf_x);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 165U);
  c1_m_A = 3.0 * c1_u;
  c1_uf_x = c1_m_A;
  c1_vf_x = c1_uf_x;
  c1_wf_x = c1_vf_x;
  c1_theta = c1_wf_x / 4.0;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 166U);
  c1_xf_x = c1_theta + c1_beta;
  c1_yf_x = c1_xf_x;
  c1_yf_x = muDoubleScalarSin(c1_yf_x);
  c1_ag_x = c1_theta + c1_beta;
  c1_bg_x = c1_ag_x;
  c1_bg_x = muDoubleScalarCos(c1_bg_x);
  c1_cg_x = 2.0 * c1_theta + 2.0 * c1_beta;
  c1_dg_x = c1_cg_x;
  c1_dg_x = muDoubleScalarCos(c1_dg_x);
  c1_n_A = 1.7320508075688772 * (c1_lamqDP * c1_yf_x - c1_lamdDP * c1_bg_x);
  c1_o_B = (c1_LqDP + c1_LdDP) - (c1_LqDP - c1_LdDP) * c1_dg_x;
  c1_eg_x = c1_n_A;
  c1_ic_y = c1_o_B;
  c1_fg_x = c1_eg_x;
  c1_jc_y = c1_ic_y;
  c1_gg_x = c1_fg_x;
  c1_kc_y = c1_jc_y;
  c1_lc_y = c1_gg_x / c1_kc_y;
  c1_hg_x = (2.0 * c1_theta + 2.0 * c1_beta) - 2.0943951023931953;
  c1_ig_x = c1_hg_x;
  c1_ig_x = muDoubleScalarCos(c1_ig_x);
  c1_jg_x = 2.0 * c1_theta + 2.0 * c1_beta;
  c1_kg_x = c1_jg_x;
  c1_kg_x = muDoubleScalarCos(c1_kg_x);
  c1_o_A = ((c1_LqDP - c1_LdDP) * c1_ig_x + 0.5 * (c1_LqDP + c1_LdDP)) * c1_idc;
  c1_p_B = (c1_LqDP + c1_LdDP) - (c1_LqDP - c1_LdDP) * c1_kg_x;
  c1_lg_x = c1_o_A;
  c1_mc_y = c1_p_B;
  c1_mg_x = c1_lg_x;
  c1_nc_y = c1_mc_y;
  c1_ng_x = c1_mg_x;
  c1_oc_y = c1_nc_y;
  c1_pc_y = c1_ng_x / c1_oc_y;
  c1_ias4 = (c1_K + c1_lc_y) - c1_pc_y;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 167U);
  c1_og_x = c1_theta + c1_beta;
  c1_pg_x = c1_og_x;
  c1_pg_x = muDoubleScalarSin(c1_pg_x);
  c1_qg_x = (c1_theta + c1_beta) + 1.0471975511965976;
  c1_rg_x = c1_qg_x;
  c1_rg_x = muDoubleScalarSin(c1_rg_x);
  c1_iqs4 = 1.1547005383792515 * (-c1_ias4 * c1_pg_x - c1_idc * c1_rg_x);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 168U);
  c1_sg_x = c1_theta + c1_beta;
  c1_tg_x = c1_sg_x;
  c1_tg_x = muDoubleScalarCos(c1_tg_x);
  c1_ug_x = (c1_theta + c1_beta) + 1.0471975511965976;
  c1_vg_x = c1_ug_x;
  c1_vg_x = muDoubleScalarCos(c1_vg_x);
  c1_ids4 = 1.1547005383792515 * (c1_ias4 * c1_tg_x + c1_idc * c1_vg_x);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 170U);
  c1_theta = c1_u;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 171U);
  c1_wg_x = c1_theta + c1_beta;
  c1_xg_x = c1_wg_x;
  c1_xg_x = muDoubleScalarSin(c1_xg_x);
  c1_yg_x = c1_theta + c1_beta;
  c1_ah_x = c1_yg_x;
  c1_ah_x = muDoubleScalarCos(c1_ah_x);
  c1_bh_x = 2.0 * c1_theta + 2.0 * c1_beta;
  c1_ch_x = c1_bh_x;
  c1_ch_x = muDoubleScalarCos(c1_ch_x);
  c1_p_A = 1.7320508075688772 * (c1_lamqDP * c1_xg_x - c1_lamdDP * c1_ah_x);
  c1_q_B = (c1_LqDP + c1_LdDP) - (c1_LqDP - c1_LdDP) * c1_ch_x;
  c1_dh_x = c1_p_A;
  c1_qc_y = c1_q_B;
  c1_eh_x = c1_dh_x;
  c1_rc_y = c1_qc_y;
  c1_fh_x = c1_eh_x;
  c1_sc_y = c1_rc_y;
  c1_tc_y = c1_fh_x / c1_sc_y;
  c1_gh_x = (2.0 * c1_theta + 2.0 * c1_beta) - 2.0943951023931953;
  c1_hh_x = c1_gh_x;
  c1_hh_x = muDoubleScalarCos(c1_hh_x);
  c1_ih_x = 2.0 * c1_theta + 2.0 * c1_beta;
  c1_jh_x = c1_ih_x;
  c1_jh_x = muDoubleScalarCos(c1_jh_x);
  c1_q_A = ((c1_LqDP - c1_LdDP) * c1_hh_x + 0.5 * (c1_LqDP + c1_LdDP)) * c1_idc;
  c1_r_B = (c1_LqDP + c1_LdDP) - (c1_LqDP - c1_LdDP) * c1_jh_x;
  c1_kh_x = c1_q_A;
  c1_uc_y = c1_r_B;
  c1_lh_x = c1_kh_x;
  c1_vc_y = c1_uc_y;
  c1_mh_x = c1_lh_x;
  c1_wc_y = c1_vc_y;
  c1_xc_y = c1_mh_x / c1_wc_y;
  c1_ias5 = (c1_K + c1_tc_y) - c1_xc_y;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 172U);
  c1_nh_x = c1_theta + c1_beta;
  c1_oh_x = c1_nh_x;
  c1_oh_x = muDoubleScalarSin(c1_oh_x);
  c1_ph_x = (c1_theta + c1_beta) + 1.0471975511965976;
  c1_qh_x = c1_ph_x;
  c1_qh_x = muDoubleScalarSin(c1_qh_x);
  c1_iqs5 = 1.1547005383792515 * (-c1_ias5 * c1_oh_x - c1_idc * c1_qh_x);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 173U);
  c1_rh_x = c1_theta + c1_beta;
  c1_sh_x = c1_rh_x;
  c1_sh_x = muDoubleScalarCos(c1_sh_x);
  c1_th_x = (c1_theta + c1_beta) + 1.0471975511965976;
  c1_uh_x = c1_th_x;
  c1_uh_x = muDoubleScalarCos(c1_uh_x);
  c1_ids5 = 1.1547005383792515 * (c1_ias5 * c1_sh_x + c1_idc * c1_uh_x);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 176U);
  c1_b_iqs1[0] = c1_iqs1;
  c1_b_iqs1[1] = c1_iqs2;
  c1_b_iqs1[2] = c1_iqs3;
  c1_b_iqs1[3] = c1_iqs4;
  c1_b_iqs1[4] = c1_iqs5;
  for (c1_i10 = 0; c1_i10 < 5; c1_i10++) {
    c1_iqs_C[c1_i10] = c1_b_iqs1[c1_i10];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 177U);
  c1_b_ids1[0] = c1_ids1;
  c1_b_ids1[1] = c1_ids2;
  c1_b_ids1[2] = c1_ids3;
  c1_b_ids1[3] = c1_ids4;
  c1_b_ids1[4] = c1_ids5;
  for (c1_i11 = 0; c1_i11 < 5; c1_i11++) {
    c1_ids_C[c1_i11] = c1_b_ids1[c1_i11];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 178U);
  c1_r_A = c1_u;
  c1_vh_x = c1_r_A;
  c1_wh_x = c1_vh_x;
  c1_xh_x = c1_wh_x;
  c1_yc_y = c1_xh_x / 4.0;
  c1_iqsCommu = c1_yc_y * 3.1415926535897931 * ((((c1_iqs_C[0] + 4.0 * c1_iqs_C
    [1]) + 2.0 * c1_iqs_C[2]) + 4.0 * c1_iqs_C[3]) + c1_iqs_C[4]);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 179U);
  c1_s_A = c1_u;
  c1_yh_x = c1_s_A;
  c1_ai_x = c1_yh_x;
  c1_bi_x = c1_ai_x;
  c1_ad_y = c1_bi_x / 4.0;
  c1_idsCommu = c1_ad_y * 3.1415926535897931 * ((((c1_ids_C[0] + 4.0 * c1_ids_C
    [1]) + 2.0 * c1_ids_C[2]) + 4.0 * c1_ids_C[3]) + c1_ids_C[4]);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 181U);
  c1_iqs = c1_iqsCond + c1_iqsCommu;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 182U);
  c1_ids = c1_idsCond + c1_idsCommu;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 185U);
  c1_ci_x = c1_alpha;
  c1_di_x = c1_ci_x;
  c1_di_x = muDoubleScalarCos(c1_di_x);
  c1_t_A = (1.6539866862653763 * c1_Epk * c1_di_x - c1_ed) - (c1_b_P.rl +
    0.954929658551372 * c1_wr * c1_Lc) * c1_idc;
  c1_s_B = c1_b_P.Ll + c1_Lt;
  c1_ei_x = c1_t_A;
  c1_bd_y = c1_s_B;
  c1_fi_x = c1_ei_x;
  c1_cd_y = c1_bd_y;
  c1_gi_x = c1_fi_x;
  c1_dd_y = c1_cd_y;
  c1_pidc = c1_gi_x / c1_dd_y;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, -185);
  _SFD_SYMBOL_SCOPE_POP();
  *c1_b_iqs = c1_iqs;
  *c1_b_ids = c1_ids;
  *c1_b_pidc = c1_pidc;
  *c1_b_beta = c1_beta;
  *c1_d_u = c1_u;
  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 0U, chartInstance->c1_sfEvent);
}

static void initSimStructsc1_gen_rec(SFc1_gen_recInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void init_script_number_translation(uint32_T c1_machineNumber, uint32_T
  c1_chartNumber, uint32_T c1_instanceNumber)
{
  (void)c1_machineNumber;
  (void)c1_chartNumber;
  (void)c1_instanceNumber;
}

static const mxArray *c1_sf_marshallOut(void *chartInstanceVoid, void *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  real_T c1_u;
  const mxArray *c1_y = NULL;
  SFc1_gen_recInstanceStruct *chartInstance;
  chartInstance = (SFc1_gen_recInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_u = *(real_T *)c1_inData;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", &c1_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static real_T c1_emlrt_marshallIn(SFc1_gen_recInstanceStruct *chartInstance,
  const mxArray *c1_u, const char_T *c1_identifier)
{
  real_T c1_y;
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_u), &c1_thisId);
  sf_mex_destroy(&c1_u);
  return c1_y;
}

static real_T c1_b_emlrt_marshallIn(SFc1_gen_recInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  real_T c1_y;
  real_T c1_d0;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_d0, 1, 0, 0U, 0, 0U, 0);
  c1_y = c1_d0;
  sf_mex_destroy(&c1_u);
  return c1_y;
}

static void c1_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_u;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y;
  SFc1_gen_recInstanceStruct *chartInstance;
  chartInstance = (SFc1_gen_recInstanceStruct *)chartInstanceVoid;
  c1_u = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_u), &c1_thisId);
  sf_mex_destroy(&c1_u);
  *(real_T *)c1_outData = c1_y;
  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_b_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  c1_struct_U3Zfiu1q58uXf1A0Qkd5GF c1_u;
  const mxArray *c1_y = NULL;
  real_T c1_b_u;
  const mxArray *c1_b_y = NULL;
  real_T c1_c_u;
  const mxArray *c1_c_y = NULL;
  real_T c1_d_u;
  const mxArray *c1_d_y = NULL;
  real_T c1_e_u;
  const mxArray *c1_e_y = NULL;
  real_T c1_f_u;
  const mxArray *c1_f_y = NULL;
  real_T c1_g_u;
  const mxArray *c1_g_y = NULL;
  real_T c1_h_u;
  const mxArray *c1_h_y = NULL;
  real_T c1_i_u;
  const mxArray *c1_i_y = NULL;
  real_T c1_j_u;
  const mxArray *c1_j_y = NULL;
  real_T c1_k_u;
  const mxArray *c1_k_y = NULL;
  real_T c1_l_u;
  const mxArray *c1_l_y = NULL;
  real_T c1_m_u;
  const mxArray *c1_m_y = NULL;
  real_T c1_n_u;
  const mxArray *c1_n_y = NULL;
  real_T c1_o_u;
  const mxArray *c1_o_y = NULL;
  real_T c1_p_u;
  const mxArray *c1_p_y = NULL;
  real_T c1_q_u;
  const mxArray *c1_q_y = NULL;
  real_T c1_r_u;
  const mxArray *c1_r_y = NULL;
  real_T c1_s_u;
  const mxArray *c1_s_y = NULL;
  real_T c1_t_u;
  const mxArray *c1_t_y = NULL;
  real_T c1_u_u;
  const mxArray *c1_u_y = NULL;
  real_T c1_v_u;
  const mxArray *c1_v_y = NULL;
  SFc1_gen_recInstanceStruct *chartInstance;
  chartInstance = (SFc1_gen_recInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_u = *(c1_struct_U3Zfiu1q58uXf1A0Qkd5GF *)c1_inData;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_createstruct("structure", 2, 1, 1), false);
  c1_b_u = c1_u.Ll;
  c1_b_y = NULL;
  sf_mex_assign(&c1_b_y, sf_mex_create("y", &c1_b_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c1_y, c1_b_y, "Ll", "Ll", 0);
  c1_c_u = c1_u.rl;
  c1_c_y = NULL;
  sf_mex_assign(&c1_c_y, sf_mex_create("y", &c1_c_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c1_y, c1_c_y, "rl", "rl", 0);
  c1_d_u = c1_u.rs;
  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", &c1_d_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c1_y, c1_d_y, "rs", "rs", 0);
  c1_e_u = c1_u.Lls;
  c1_e_y = NULL;
  sf_mex_assign(&c1_e_y, sf_mex_create("y", &c1_e_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c1_y, c1_e_y, "Lls", "Lls", 0);
  c1_f_u = c1_u.Lmq;
  c1_f_y = NULL;
  sf_mex_assign(&c1_f_y, sf_mex_create("y", &c1_f_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c1_y, c1_f_y, "Lmq", "Lmq", 0);
  c1_g_u = c1_u.Lmd;
  c1_g_y = NULL;
  sf_mex_assign(&c1_g_y, sf_mex_create("y", &c1_g_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c1_y, c1_g_y, "Lmd", "Lmd", 0);
  c1_h_u = c1_u.rd1;
  c1_h_y = NULL;
  sf_mex_assign(&c1_h_y, sf_mex_create("y", &c1_h_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c1_y, c1_h_y, "rd1", "rd1", 0);
  c1_i_u = c1_u.rd2;
  c1_i_y = NULL;
  sf_mex_assign(&c1_i_y, sf_mex_create("y", &c1_i_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c1_y, c1_i_y, "rd2", "rd2", 0);
  c1_j_u = c1_u.rd3;
  c1_j_y = NULL;
  sf_mex_assign(&c1_j_y, sf_mex_create("y", &c1_j_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c1_y, c1_j_y, "rd3", "rd3", 0);
  c1_k_u = c1_u.rfd;
  c1_k_y = NULL;
  sf_mex_assign(&c1_k_y, sf_mex_create("y", &c1_k_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c1_y, c1_k_y, "rfd", "rfd", 0);
  c1_l_u = c1_u.Lld1;
  c1_l_y = NULL;
  sf_mex_assign(&c1_l_y, sf_mex_create("y", &c1_l_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c1_y, c1_l_y, "Lld1", "Lld1", 0);
  c1_m_u = c1_u.Lld2;
  c1_m_y = NULL;
  sf_mex_assign(&c1_m_y, sf_mex_create("y", &c1_m_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c1_y, c1_m_y, "Lld2", "Lld2", 0);
  c1_n_u = c1_u.Lld3;
  c1_n_y = NULL;
  sf_mex_assign(&c1_n_y, sf_mex_create("y", &c1_n_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c1_y, c1_n_y, "Lld3", "Lld3", 0);
  c1_o_u = c1_u.Llfd;
  c1_o_y = NULL;
  sf_mex_assign(&c1_o_y, sf_mex_create("y", &c1_o_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c1_y, c1_o_y, "Llfd", "Llfd", 0);
  c1_p_u = c1_u.rq1;
  c1_p_y = NULL;
  sf_mex_assign(&c1_p_y, sf_mex_create("y", &c1_p_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c1_y, c1_p_y, "rq1", "rq1", 0);
  c1_q_u = c1_u.rq2;
  c1_q_y = NULL;
  sf_mex_assign(&c1_q_y, sf_mex_create("y", &c1_q_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c1_y, c1_q_y, "rq2", "rq2", 0);
  c1_r_u = c1_u.rq3;
  c1_r_y = NULL;
  sf_mex_assign(&c1_r_y, sf_mex_create("y", &c1_r_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c1_y, c1_r_y, "rq3", "rq3", 0);
  c1_s_u = c1_u.pole;
  c1_s_y = NULL;
  sf_mex_assign(&c1_s_y, sf_mex_create("y", &c1_s_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c1_y, c1_s_y, "pole", "pole", 0);
  c1_t_u = c1_u.Llq1;
  c1_t_y = NULL;
  sf_mex_assign(&c1_t_y, sf_mex_create("y", &c1_t_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c1_y, c1_t_y, "Llq1", "Llq1", 0);
  c1_u_u = c1_u.Llq2;
  c1_u_y = NULL;
  sf_mex_assign(&c1_u_y, sf_mex_create("y", &c1_u_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c1_y, c1_u_y, "Llq2", "Llq2", 0);
  c1_v_u = c1_u.Llq3;
  c1_v_y = NULL;
  sf_mex_assign(&c1_v_y, sf_mex_create("y", &c1_v_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c1_y, c1_v_y, "Llq3", "Llq3", 0);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static void c1_c_emlrt_marshallIn(SFc1_gen_recInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  c1_struct_U3Zfiu1q58uXf1A0Qkd5GF *c1_y)
{
  emlrtMsgIdentifier c1_thisId;
  static const char * c1_fieldNames[21] = { "Ll", "rl", "rs", "Lls", "Lmq",
    "Lmd", "rd1", "rd2", "rd3", "rfd", "Lld1", "Lld2", "Lld3", "Llfd", "rq1",
    "rq2", "rq3", "pole", "Llq1", "Llq2", "Llq3" };

  c1_thisId.fParent = c1_parentId;
  sf_mex_check_struct(c1_parentId, c1_u, 21, c1_fieldNames, 0U, NULL);
  c1_thisId.fIdentifier = "Ll";
  c1_y->Ll = c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c1_u, "Ll", "Ll", 0)), &c1_thisId);
  c1_thisId.fIdentifier = "rl";
  c1_y->rl = c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c1_u, "rl", "rl", 0)), &c1_thisId);
  c1_thisId.fIdentifier = "rs";
  c1_y->rs = c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c1_u, "rs", "rs", 0)), &c1_thisId);
  c1_thisId.fIdentifier = "Lls";
  c1_y->Lls = c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c1_u, "Lls", "Lls", 0)), &c1_thisId);
  c1_thisId.fIdentifier = "Lmq";
  c1_y->Lmq = c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c1_u, "Lmq", "Lmq", 0)), &c1_thisId);
  c1_thisId.fIdentifier = "Lmd";
  c1_y->Lmd = c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c1_u, "Lmd", "Lmd", 0)), &c1_thisId);
  c1_thisId.fIdentifier = "rd1";
  c1_y->rd1 = c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c1_u, "rd1", "rd1", 0)), &c1_thisId);
  c1_thisId.fIdentifier = "rd2";
  c1_y->rd2 = c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c1_u, "rd2", "rd2", 0)), &c1_thisId);
  c1_thisId.fIdentifier = "rd3";
  c1_y->rd3 = c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c1_u, "rd3", "rd3", 0)), &c1_thisId);
  c1_thisId.fIdentifier = "rfd";
  c1_y->rfd = c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c1_u, "rfd", "rfd", 0)), &c1_thisId);
  c1_thisId.fIdentifier = "Lld1";
  c1_y->Lld1 = c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c1_u, "Lld1", "Lld1", 0)), &c1_thisId);
  c1_thisId.fIdentifier = "Lld2";
  c1_y->Lld2 = c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c1_u, "Lld2", "Lld2", 0)), &c1_thisId);
  c1_thisId.fIdentifier = "Lld3";
  c1_y->Lld3 = c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c1_u, "Lld3", "Lld3", 0)), &c1_thisId);
  c1_thisId.fIdentifier = "Llfd";
  c1_y->Llfd = c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c1_u, "Llfd", "Llfd", 0)), &c1_thisId);
  c1_thisId.fIdentifier = "rq1";
  c1_y->rq1 = c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c1_u, "rq1", "rq1", 0)), &c1_thisId);
  c1_thisId.fIdentifier = "rq2";
  c1_y->rq2 = c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c1_u, "rq2", "rq2", 0)), &c1_thisId);
  c1_thisId.fIdentifier = "rq3";
  c1_y->rq3 = c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c1_u, "rq3", "rq3", 0)), &c1_thisId);
  c1_thisId.fIdentifier = "pole";
  c1_y->pole = c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c1_u, "pole", "pole", 0)), &c1_thisId);
  c1_thisId.fIdentifier = "Llq1";
  c1_y->Llq1 = c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c1_u, "Llq1", "Llq1", 0)), &c1_thisId);
  c1_thisId.fIdentifier = "Llq2";
  c1_y->Llq2 = c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c1_u, "Llq2", "Llq2", 0)), &c1_thisId);
  c1_thisId.fIdentifier = "Llq3";
  c1_y->Llq3 = c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c1_u, "Llq3", "Llq3", 0)), &c1_thisId);
  sf_mex_destroy(&c1_u);
}

static void c1_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_b_P;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  c1_struct_U3Zfiu1q58uXf1A0Qkd5GF c1_y;
  SFc1_gen_recInstanceStruct *chartInstance;
  chartInstance = (SFc1_gen_recInstanceStruct *)chartInstanceVoid;
  c1_b_P = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_P), &c1_thisId, &c1_y);
  sf_mex_destroy(&c1_b_P);
  *(c1_struct_U3Zfiu1q58uXf1A0Qkd5GF *)c1_outData = c1_y;
  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_c_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_i12;
  real_T c1_b_inData[7];
  int32_T c1_i13;
  real_T c1_u[7];
  const mxArray *c1_y = NULL;
  SFc1_gen_recInstanceStruct *chartInstance;
  chartInstance = (SFc1_gen_recInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  for (c1_i12 = 0; c1_i12 < 7; c1_i12++) {
    c1_b_inData[c1_i12] = (*(real_T (*)[7])c1_inData)[c1_i12];
  }

  for (c1_i13 = 0; c1_i13 < 7; c1_i13++) {
    c1_u[c1_i13] = c1_b_inData[c1_i13];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 1, 7), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static const mxArray *c1_d_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_i14;
  real_T c1_b_inData[5];
  int32_T c1_i15;
  real_T c1_u[5];
  const mxArray *c1_y = NULL;
  SFc1_gen_recInstanceStruct *chartInstance;
  chartInstance = (SFc1_gen_recInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  for (c1_i14 = 0; c1_i14 < 5; c1_i14++) {
    c1_b_inData[c1_i14] = (*(real_T (*)[5])c1_inData)[c1_i14];
  }

  for (c1_i15 = 0; c1_i15 < 5; c1_i15++) {
    c1_u[c1_i15] = c1_b_inData[c1_i15];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 1, 5), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static void c1_d_emlrt_marshallIn(SFc1_gen_recInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[5])
{
  real_T c1_dv5[5];
  int32_T c1_i16;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv5, 1, 0, 0U, 1, 0U, 1, 5);
  for (c1_i16 = 0; c1_i16 < 5; c1_i16++) {
    c1_y[c1_i16] = c1_dv5[c1_i16];
  }

  sf_mex_destroy(&c1_u);
}

static void c1_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_ids_C;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y[5];
  int32_T c1_i17;
  SFc1_gen_recInstanceStruct *chartInstance;
  chartInstance = (SFc1_gen_recInstanceStruct *)chartInstanceVoid;
  c1_ids_C = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_ids_C), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_ids_C);
  for (c1_i17 = 0; c1_i17 < 5; c1_i17++) {
    (*(real_T (*)[5])c1_outData)[c1_i17] = c1_y[c1_i17];
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_e_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_i18;
  real_T c1_b_inData[100];
  int32_T c1_i19;
  real_T c1_u[100];
  const mxArray *c1_y = NULL;
  SFc1_gen_recInstanceStruct *chartInstance;
  chartInstance = (SFc1_gen_recInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  for (c1_i18 = 0; c1_i18 < 100; c1_i18++) {
    c1_b_inData[c1_i18] = (*(real_T (*)[100])c1_inData)[c1_i18];
  }

  for (c1_i19 = 0; c1_i19 < 100; c1_i19++) {
    c1_u[c1_i19] = c1_b_inData[c1_i19];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 2, 1, 100), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static void c1_e_emlrt_marshallIn(SFc1_gen_recInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[100])
{
  real_T c1_dv6[100];
  int32_T c1_i20;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv6, 1, 0, 0U, 1, 0U, 2, 1,
                100);
  for (c1_i20 = 0; c1_i20 < 100; c1_i20++) {
    c1_y[c1_i20] = c1_dv6[c1_i20];
  }

  sf_mex_destroy(&c1_u);
}

static void c1_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_record;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y[100];
  int32_T c1_i21;
  SFc1_gen_recInstanceStruct *chartInstance;
  chartInstance = (SFc1_gen_recInstanceStruct *)chartInstanceVoid;
  c1_record = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_record), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_record);
  for (c1_i21 = 0; c1_i21 < 100; c1_i21++) {
    (*(real_T (*)[100])c1_outData)[c1_i21] = c1_y[c1_i21];
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

const mxArray *sf_c1_gen_rec_get_eml_resolved_functions_info(void)
{
  const mxArray *c1_nameCaptureInfo = NULL;
  c1_nameCaptureInfo = NULL;
  sf_mex_assign(&c1_nameCaptureInfo, sf_mex_createstruct("structure", 2, 51, 1),
                false);
  c1_info_helper(&c1_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c1_nameCaptureInfo);
  return c1_nameCaptureInfo;
}

static void c1_info_helper(const mxArray **c1_info)
{
  const mxArray *c1_rhs0 = NULL;
  const mxArray *c1_lhs0 = NULL;
  const mxArray *c1_rhs1 = NULL;
  const mxArray *c1_lhs1 = NULL;
  const mxArray *c1_rhs2 = NULL;
  const mxArray *c1_lhs2 = NULL;
  const mxArray *c1_rhs3 = NULL;
  const mxArray *c1_lhs3 = NULL;
  const mxArray *c1_rhs4 = NULL;
  const mxArray *c1_lhs4 = NULL;
  const mxArray *c1_rhs5 = NULL;
  const mxArray *c1_lhs5 = NULL;
  const mxArray *c1_rhs6 = NULL;
  const mxArray *c1_lhs6 = NULL;
  const mxArray *c1_rhs7 = NULL;
  const mxArray *c1_lhs7 = NULL;
  const mxArray *c1_rhs8 = NULL;
  const mxArray *c1_lhs8 = NULL;
  const mxArray *c1_rhs9 = NULL;
  const mxArray *c1_lhs9 = NULL;
  const mxArray *c1_rhs10 = NULL;
  const mxArray *c1_lhs10 = NULL;
  const mxArray *c1_rhs11 = NULL;
  const mxArray *c1_lhs11 = NULL;
  const mxArray *c1_rhs12 = NULL;
  const mxArray *c1_lhs12 = NULL;
  const mxArray *c1_rhs13 = NULL;
  const mxArray *c1_lhs13 = NULL;
  const mxArray *c1_rhs14 = NULL;
  const mxArray *c1_lhs14 = NULL;
  const mxArray *c1_rhs15 = NULL;
  const mxArray *c1_lhs15 = NULL;
  const mxArray *c1_rhs16 = NULL;
  const mxArray *c1_lhs16 = NULL;
  const mxArray *c1_rhs17 = NULL;
  const mxArray *c1_lhs17 = NULL;
  const mxArray *c1_rhs18 = NULL;
  const mxArray *c1_lhs18 = NULL;
  const mxArray *c1_rhs19 = NULL;
  const mxArray *c1_lhs19 = NULL;
  const mxArray *c1_rhs20 = NULL;
  const mxArray *c1_lhs20 = NULL;
  const mxArray *c1_rhs21 = NULL;
  const mxArray *c1_lhs21 = NULL;
  const mxArray *c1_rhs22 = NULL;
  const mxArray *c1_lhs22 = NULL;
  const mxArray *c1_rhs23 = NULL;
  const mxArray *c1_lhs23 = NULL;
  const mxArray *c1_rhs24 = NULL;
  const mxArray *c1_lhs24 = NULL;
  const mxArray *c1_rhs25 = NULL;
  const mxArray *c1_lhs25 = NULL;
  const mxArray *c1_rhs26 = NULL;
  const mxArray *c1_lhs26 = NULL;
  const mxArray *c1_rhs27 = NULL;
  const mxArray *c1_lhs27 = NULL;
  const mxArray *c1_rhs28 = NULL;
  const mxArray *c1_lhs28 = NULL;
  const mxArray *c1_rhs29 = NULL;
  const mxArray *c1_lhs29 = NULL;
  const mxArray *c1_rhs30 = NULL;
  const mxArray *c1_lhs30 = NULL;
  const mxArray *c1_rhs31 = NULL;
  const mxArray *c1_lhs31 = NULL;
  const mxArray *c1_rhs32 = NULL;
  const mxArray *c1_lhs32 = NULL;
  const mxArray *c1_rhs33 = NULL;
  const mxArray *c1_lhs33 = NULL;
  const mxArray *c1_rhs34 = NULL;
  const mxArray *c1_lhs34 = NULL;
  const mxArray *c1_rhs35 = NULL;
  const mxArray *c1_lhs35 = NULL;
  const mxArray *c1_rhs36 = NULL;
  const mxArray *c1_lhs36 = NULL;
  const mxArray *c1_rhs37 = NULL;
  const mxArray *c1_lhs37 = NULL;
  const mxArray *c1_rhs38 = NULL;
  const mxArray *c1_lhs38 = NULL;
  const mxArray *c1_rhs39 = NULL;
  const mxArray *c1_lhs39 = NULL;
  const mxArray *c1_rhs40 = NULL;
  const mxArray *c1_lhs40 = NULL;
  const mxArray *c1_rhs41 = NULL;
  const mxArray *c1_lhs41 = NULL;
  const mxArray *c1_rhs42 = NULL;
  const mxArray *c1_lhs42 = NULL;
  const mxArray *c1_rhs43 = NULL;
  const mxArray *c1_lhs43 = NULL;
  const mxArray *c1_rhs44 = NULL;
  const mxArray *c1_lhs44 = NULL;
  const mxArray *c1_rhs45 = NULL;
  const mxArray *c1_lhs45 = NULL;
  const mxArray *c1_rhs46 = NULL;
  const mxArray *c1_lhs46 = NULL;
  const mxArray *c1_rhs47 = NULL;
  const mxArray *c1_lhs47 = NULL;
  const mxArray *c1_rhs48 = NULL;
  const mxArray *c1_lhs48 = NULL;
  const mxArray *c1_rhs49 = NULL;
  const mxArray *c1_lhs49 = NULL;
  const mxArray *c1_rhs50 = NULL;
  const mxArray *c1_lhs50 = NULL;
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("mrdivide"), "name", "name", 0);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1388463696U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1370017086U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c1_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs0), "rhs", "rhs", 0);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs0), "lhs", "lhs", 0);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 1);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 1);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c1_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs1), "rhs", "rhs", 1);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs1), "lhs", "lhs", 1);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 2);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("rdivide"), "name", "name", 2);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 2);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c1_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs2), "rhs", "rhs", 2);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs2), "lhs", "lhs", 2);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 3);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 3);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 3);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c1_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs3), "rhs", "rhs", 3);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs3), "lhs", "lhs", 3);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 4);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 4);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825996U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c1_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs4), "rhs", "rhs", 4);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs4), "lhs", "lhs", 4);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 5);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_div"), "name", "name", 5);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 5);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c1_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs5), "rhs", "rhs", 5);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs5), "lhs", "lhs", 5);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "context",
                  "context", 6);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.div"), "name",
                  "name", 6);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p"), "resolved",
                  "resolved", 6);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c1_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs6), "rhs", "rhs", 6);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs6), "lhs", "lhs", 6);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 7);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("rdivide"), "name", "name", 7);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 7);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c1_rhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs7), "rhs", "rhs", 7);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs7), "lhs", "lhs", 7);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 8);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("sum"), "name", "name", 8);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/sum.m"), "resolved",
                  "resolved", 8);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717458U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c1_rhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs8), "rhs", "rhs", 8);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs8), "lhs", "lhs", 8);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/sum.m"), "context",
                  "context", 9);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 9);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 9);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c1_rhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs9), "rhs", "rhs", 9);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs9), "lhs", "lhs", 9);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/sum.m"), "context",
                  "context", 10);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("isequal"), "name", "name", 10);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m"), "resolved",
                  "resolved", 10);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825958U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c1_rhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m"), "context",
                  "context", 11);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_isequal_core"), "name",
                  "name", 11);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isequal_core.m"),
                  "resolved", "resolved", 11);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825986U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c1_rhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/sum.m"), "context",
                  "context", 12);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_const_nonsingleton_dim"),
                  "name", "name", 12);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_const_nonsingleton_dim.m"),
                  "resolved", "resolved", 12);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c1_rhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_const_nonsingleton_dim.m"),
                  "context", "context", 13);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.constNonSingletonDim"), "name", "name", 13);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/constNonSingletonDim.m"),
                  "resolved", "resolved", 13);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c1_rhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/sum.m"), "context",
                  "context", 14);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 14);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 14);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c1_rhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 15);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 15);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 15);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c1_rhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/sum.m"), "context",
                  "context", 16);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 16);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 16);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c1_rhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/sum.m"), "context",
                  "context", 17);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 17);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 17);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c1_rhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 18);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("intmax"), "name", "name", 18);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 18);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c1_rhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "context",
                  "context", 19);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 19);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 19);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c1_rhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "context", "context", 20);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("isequal"), "name", "name", 20);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m"), "resolved",
                  "resolved", 20);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825958U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c1_rhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isequal_core.m!isequal_scalar"),
                  "context", "context", 21);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("isnan"), "name", "name", 21);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 21);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717458U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c1_rhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "context",
                  "context", 22);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 22);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 22);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c1_rhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs22), "lhs", "lhs",
                  22);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 23);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("sqrt"), "name", "name", 23);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 23);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "resolved",
                  "resolved", 23);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1343837586U), "fileTimeLo",
                  "fileTimeLo", 23);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 23);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 23);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 23);
  sf_mex_assign(&c1_rhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs23), "rhs", "rhs",
                  23);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs23), "lhs", "lhs",
                  23);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "context",
                  "context", 24);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_error"), "name", "name",
                  24);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 24);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m"), "resolved",
                  "resolved", 24);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1343837558U), "fileTimeLo",
                  "fileTimeLo", 24);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 24);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 24);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 24);
  sf_mex_assign(&c1_rhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs24), "rhs", "rhs",
                  24);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs24), "lhs", "lhs",
                  24);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "context",
                  "context", 25);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_sqrt"), "name",
                  "name", 25);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m"),
                  "resolved", "resolved", 25);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825938U), "fileTimeLo",
                  "fileTimeLo", 25);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 25);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 25);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 25);
  sf_mex_assign(&c1_rhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs25), "rhs", "rhs",
                  25);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs25), "lhs", "lhs",
                  25);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 26);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("cos"), "name", "name", 26);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 26);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "resolved",
                  "resolved", 26);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1343837572U), "fileTimeLo",
                  "fileTimeLo", 26);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 26);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 26);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 26);
  sf_mex_assign(&c1_rhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs26), "rhs", "rhs",
                  26);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs26), "lhs", "lhs",
                  26);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "context",
                  "context", 27);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_cos"), "name",
                  "name", 27);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 27);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m"),
                  "resolved", "resolved", 27);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825922U), "fileTimeLo",
                  "fileTimeLo", 27);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 27);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 27);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 27);
  sf_mex_assign(&c1_rhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs27), "rhs", "rhs",
                  27);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs27), "lhs", "lhs",
                  27);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 28);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("sin"), "name", "name", 28);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 28);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "resolved",
                  "resolved", 28);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1343837586U), "fileTimeLo",
                  "fileTimeLo", 28);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 28);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 28);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 28);
  sf_mex_assign(&c1_rhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs28), "rhs", "rhs",
                  28);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs28), "lhs", "lhs",
                  28);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "context",
                  "context", 29);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_sin"), "name",
                  "name", 29);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 29);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m"),
                  "resolved", "resolved", 29);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825936U), "fileTimeLo",
                  "fileTimeLo", 29);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 29);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 29);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 29);
  sf_mex_assign(&c1_rhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs29), "rhs", "rhs",
                  29);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs29), "lhs", "lhs",
                  29);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 30);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("disp"), "name", "name", 30);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 30);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXMB]$matlabroot$/toolbox/matlab/lang/disp"), "resolved", "resolved", 30);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(MAX_uint32_T), "fileTimeLo",
                  "fileTimeLo", 30);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(MAX_uint32_T), "fileTimeHi",
                  "fileTimeHi", 30);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(MAX_uint32_T), "mFileTimeLo",
                  "mFileTimeLo", 30);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(MAX_uint32_T), "mFileTimeHi",
                  "mFileTimeHi", 30);
  sf_mex_assign(&c1_rhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs30), "rhs", "rhs",
                  30);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs30), "lhs", "lhs",
                  30);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 31);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("abs"), "name", "name", 31);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 31);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 31);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 31);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 31);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 31);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 31);
  sf_mex_assign(&c1_rhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs31), "rhs", "rhs",
                  31);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs31), "lhs", "lhs",
                  31);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 32);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 32);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 32);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 32);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 32);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 32);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 32);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 32);
  sf_mex_assign(&c1_rhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs32), "rhs", "rhs",
                  32);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs32), "lhs", "lhs",
                  32);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 33);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 33);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 33);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 33);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825912U), "fileTimeLo",
                  "fileTimeLo", 33);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 33);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 33);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 33);
  sf_mex_assign(&c1_rhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs33), "rhs", "rhs",
                  33);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs33), "lhs", "lhs",
                  33);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 34);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 34);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 34);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 34);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 34);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 34);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 34);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 34);
  sf_mex_assign(&c1_rhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs34), "rhs", "rhs",
                  34);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs34), "lhs", "lhs",
                  34);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m!common_checks"),
                  "context", "context", 35);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 35);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 35);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 35);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 35);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 35);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 35);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 35);
  sf_mex_assign(&c1_rhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs35), "rhs", "rhs",
                  35);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs35), "lhs", "lhs",
                  35);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 36);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("angle"), "name", "name", 36);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 36);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/angle.m"), "resolved",
                  "resolved", 36);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1343837570U), "fileTimeLo",
                  "fileTimeLo", 36);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 36);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 36);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 36);
  sf_mex_assign(&c1_rhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs36), "rhs", "rhs",
                  36);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs36), "lhs", "lhs",
                  36);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/angle.m"), "context",
                  "context", 37);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_angle"), "name",
                  "name", 37);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 37);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_angle.m"),
                  "resolved", "resolved", 37);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825916U), "fileTimeLo",
                  "fileTimeLo", 37);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 37);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 37);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 37);
  sf_mex_assign(&c1_rhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs37), "rhs", "rhs",
                  37);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs37), "lhs", "lhs",
                  37);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_angle.m"),
                  "context", "context", 38);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_atan2"), "name",
                  "name", 38);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 38);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_atan2.m"),
                  "resolved", "resolved", 38);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825920U), "fileTimeLo",
                  "fileTimeLo", 38);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 38);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 38);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 38);
  sf_mex_assign(&c1_rhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs38), "rhs", "rhs",
                  38);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs38), "lhs", "lhs",
                  38);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 39);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("mpower"), "name", "name", 39);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 39);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "resolved",
                  "resolved", 39);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 39);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 39);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 39);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 39);
  sf_mex_assign(&c1_rhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs39), "rhs", "rhs",
                  39);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs39), "lhs", "lhs",
                  39);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 40);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 40);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 40);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 40);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 40);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 40);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 40);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 40);
  sf_mex_assign(&c1_rhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs40), "rhs", "rhs",
                  40);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs40), "lhs", "lhs",
                  40);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 41);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("ismatrix"), "name", "name", 41);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 41);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/ismatrix.m"), "resolved",
                  "resolved", 41);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1331308458U), "fileTimeLo",
                  "fileTimeLo", 41);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 41);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 41);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 41);
  sf_mex_assign(&c1_rhs41, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs41, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs41), "rhs", "rhs",
                  41);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs41), "lhs", "lhs",
                  41);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 42);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("power"), "name", "name", 42);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 42);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "resolved",
                  "resolved", 42);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 42);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 42);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 42);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 42);
  sf_mex_assign(&c1_rhs42, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs42, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs42), "rhs", "rhs",
                  42);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs42), "lhs", "lhs",
                  42);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "context",
                  "context", 43);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 43);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 43);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 43);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 43);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 43);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 43);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 43);
  sf_mex_assign(&c1_rhs43, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs43, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs43), "rhs", "rhs",
                  43);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs43), "lhs", "lhs",
                  43);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 44);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 44);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 44);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 44);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 44);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 44);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 44);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 44);
  sf_mex_assign(&c1_rhs44, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs44, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs44), "rhs", "rhs",
                  44);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs44), "lhs", "lhs",
                  44);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 45);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 45);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 45);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 45);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 45);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 45);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 45);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 45);
  sf_mex_assign(&c1_rhs45, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs45, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs45), "rhs", "rhs",
                  45);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs45), "lhs", "lhs",
                  45);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "context", "context", 46);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.scalexpAlloc"),
                  "name", "name", 46);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 46);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalexpAlloc.p"),
                  "resolved", "resolved", 46);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 46);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 46);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 46);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 46);
  sf_mex_assign(&c1_rhs46, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs46, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs46), "rhs", "rhs",
                  46);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs46), "lhs", "lhs",
                  46);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 47);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("floor"), "name", "name", 47);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 47);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "resolved",
                  "resolved", 47);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717454U), "fileTimeLo",
                  "fileTimeLo", 47);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 47);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 47);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 47);
  sf_mex_assign(&c1_rhs47, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs47, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs47), "rhs", "rhs",
                  47);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs47), "lhs", "lhs",
                  47);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 48);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 48);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 48);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 48);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 48);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 48);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 48);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 48);
  sf_mex_assign(&c1_rhs48, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs48, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs48), "rhs", "rhs",
                  48);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs48), "lhs", "lhs",
                  48);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 49);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_floor"), "name",
                  "name", 49);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 49);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m"),
                  "resolved", "resolved", 49);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825926U), "fileTimeLo",
                  "fileTimeLo", 49);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 49);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 49);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 49);
  sf_mex_assign(&c1_rhs49, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs49, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs49), "rhs", "rhs",
                  49);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs49), "lhs", "lhs",
                  49);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 50);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 50);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 50);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 50);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 50);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 50);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 50);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 50);
  sf_mex_assign(&c1_rhs50, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs50, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs50), "rhs", "rhs",
                  50);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs50), "lhs", "lhs",
                  50);
  sf_mex_destroy(&c1_rhs0);
  sf_mex_destroy(&c1_lhs0);
  sf_mex_destroy(&c1_rhs1);
  sf_mex_destroy(&c1_lhs1);
  sf_mex_destroy(&c1_rhs2);
  sf_mex_destroy(&c1_lhs2);
  sf_mex_destroy(&c1_rhs3);
  sf_mex_destroy(&c1_lhs3);
  sf_mex_destroy(&c1_rhs4);
  sf_mex_destroy(&c1_lhs4);
  sf_mex_destroy(&c1_rhs5);
  sf_mex_destroy(&c1_lhs5);
  sf_mex_destroy(&c1_rhs6);
  sf_mex_destroy(&c1_lhs6);
  sf_mex_destroy(&c1_rhs7);
  sf_mex_destroy(&c1_lhs7);
  sf_mex_destroy(&c1_rhs8);
  sf_mex_destroy(&c1_lhs8);
  sf_mex_destroy(&c1_rhs9);
  sf_mex_destroy(&c1_lhs9);
  sf_mex_destroy(&c1_rhs10);
  sf_mex_destroy(&c1_lhs10);
  sf_mex_destroy(&c1_rhs11);
  sf_mex_destroy(&c1_lhs11);
  sf_mex_destroy(&c1_rhs12);
  sf_mex_destroy(&c1_lhs12);
  sf_mex_destroy(&c1_rhs13);
  sf_mex_destroy(&c1_lhs13);
  sf_mex_destroy(&c1_rhs14);
  sf_mex_destroy(&c1_lhs14);
  sf_mex_destroy(&c1_rhs15);
  sf_mex_destroy(&c1_lhs15);
  sf_mex_destroy(&c1_rhs16);
  sf_mex_destroy(&c1_lhs16);
  sf_mex_destroy(&c1_rhs17);
  sf_mex_destroy(&c1_lhs17);
  sf_mex_destroy(&c1_rhs18);
  sf_mex_destroy(&c1_lhs18);
  sf_mex_destroy(&c1_rhs19);
  sf_mex_destroy(&c1_lhs19);
  sf_mex_destroy(&c1_rhs20);
  sf_mex_destroy(&c1_lhs20);
  sf_mex_destroy(&c1_rhs21);
  sf_mex_destroy(&c1_lhs21);
  sf_mex_destroy(&c1_rhs22);
  sf_mex_destroy(&c1_lhs22);
  sf_mex_destroy(&c1_rhs23);
  sf_mex_destroy(&c1_lhs23);
  sf_mex_destroy(&c1_rhs24);
  sf_mex_destroy(&c1_lhs24);
  sf_mex_destroy(&c1_rhs25);
  sf_mex_destroy(&c1_lhs25);
  sf_mex_destroy(&c1_rhs26);
  sf_mex_destroy(&c1_lhs26);
  sf_mex_destroy(&c1_rhs27);
  sf_mex_destroy(&c1_lhs27);
  sf_mex_destroy(&c1_rhs28);
  sf_mex_destroy(&c1_lhs28);
  sf_mex_destroy(&c1_rhs29);
  sf_mex_destroy(&c1_lhs29);
  sf_mex_destroy(&c1_rhs30);
  sf_mex_destroy(&c1_lhs30);
  sf_mex_destroy(&c1_rhs31);
  sf_mex_destroy(&c1_lhs31);
  sf_mex_destroy(&c1_rhs32);
  sf_mex_destroy(&c1_lhs32);
  sf_mex_destroy(&c1_rhs33);
  sf_mex_destroy(&c1_lhs33);
  sf_mex_destroy(&c1_rhs34);
  sf_mex_destroy(&c1_lhs34);
  sf_mex_destroy(&c1_rhs35);
  sf_mex_destroy(&c1_lhs35);
  sf_mex_destroy(&c1_rhs36);
  sf_mex_destroy(&c1_lhs36);
  sf_mex_destroy(&c1_rhs37);
  sf_mex_destroy(&c1_lhs37);
  sf_mex_destroy(&c1_rhs38);
  sf_mex_destroy(&c1_lhs38);
  sf_mex_destroy(&c1_rhs39);
  sf_mex_destroy(&c1_lhs39);
  sf_mex_destroy(&c1_rhs40);
  sf_mex_destroy(&c1_lhs40);
  sf_mex_destroy(&c1_rhs41);
  sf_mex_destroy(&c1_lhs41);
  sf_mex_destroy(&c1_rhs42);
  sf_mex_destroy(&c1_lhs42);
  sf_mex_destroy(&c1_rhs43);
  sf_mex_destroy(&c1_lhs43);
  sf_mex_destroy(&c1_rhs44);
  sf_mex_destroy(&c1_lhs44);
  sf_mex_destroy(&c1_rhs45);
  sf_mex_destroy(&c1_lhs45);
  sf_mex_destroy(&c1_rhs46);
  sf_mex_destroy(&c1_lhs46);
  sf_mex_destroy(&c1_rhs47);
  sf_mex_destroy(&c1_lhs47);
  sf_mex_destroy(&c1_rhs48);
  sf_mex_destroy(&c1_lhs48);
  sf_mex_destroy(&c1_rhs49);
  sf_mex_destroy(&c1_lhs49);
  sf_mex_destroy(&c1_rhs50);
  sf_mex_destroy(&c1_lhs50);
}

static const mxArray *c1_emlrt_marshallOut(const char * c1_u)
{
  const mxArray *c1_y = NULL;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c1_u)), false);
  return c1_y;
}

static const mxArray *c1_b_emlrt_marshallOut(const uint32_T c1_u)
{
  const mxArray *c1_y = NULL;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", &c1_u, 7, 0U, 0U, 0U, 0), false);
  return c1_y;
}

static void c1_rdivide(SFc1_gen_recInstanceStruct *chartInstance, real_T c1_x,
  real_T c1_y[3], real_T c1_z[3])
{
  real_T c1_b_x;
  real_T c1_c_x;
  int32_T c1_i22;
  (void)chartInstance;
  c1_b_x = c1_x;
  c1_c_x = c1_b_x;
  for (c1_i22 = 0; c1_i22 < 3; c1_i22++) {
    c1_z[c1_i22] = c1_c_x / c1_y[c1_i22];
  }
}

static real_T c1_sum(SFc1_gen_recInstanceStruct *chartInstance, real_T c1_x[3])
{
  real_T c1_y;
  int32_T c1_k;
  int32_T c1_b_k;
  (void)chartInstance;
  c1_y = c1_x[0];
  for (c1_k = 2; c1_k < 4; c1_k++) {
    c1_b_k = c1_k;
    c1_y += c1_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c1_b_k), 1, 3, 1, 0) - 1];
  }

  return c1_y;
}

static void c1_b_rdivide(SFc1_gen_recInstanceStruct *chartInstance, real_T c1_x
  [3], real_T c1_y[3], real_T c1_z[3])
{
  int32_T c1_i23;
  (void)chartInstance;
  for (c1_i23 = 0; c1_i23 < 3; c1_i23++) {
    c1_z[c1_i23] = c1_x[c1_i23] / c1_y[c1_i23];
  }
}

static real_T c1_sqrt(SFc1_gen_recInstanceStruct *chartInstance, real_T c1_x)
{
  real_T c1_b_x;
  c1_b_x = c1_x;
  c1_b_sqrt(chartInstance, &c1_b_x);
  return c1_b_x;
}

static void c1_eml_error(SFc1_gen_recInstanceStruct *chartInstance)
{
  int32_T c1_i24;
  static char_T c1_cv2[30] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'E', 'l', 'F', 'u', 'n', 'D', 'o', 'm', 'a', 'i', 'n',
    'E', 'r', 'r', 'o', 'r' };

  char_T c1_u[30];
  const mxArray *c1_y = NULL;
  int32_T c1_i25;
  static char_T c1_cv3[4] = { 's', 'q', 'r', 't' };

  char_T c1_b_u[4];
  const mxArray *c1_b_y = NULL;
  (void)chartInstance;
  for (c1_i24 = 0; c1_i24 < 30; c1_i24++) {
    c1_u[c1_i24] = c1_cv2[c1_i24];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 10, 0U, 1U, 0U, 2, 1, 30), false);
  for (c1_i25 = 0; c1_i25 < 4; c1_i25++) {
    c1_b_u[c1_i25] = c1_cv3[c1_i25];
  }

  c1_b_y = NULL;
  sf_mex_assign(&c1_b_y, sf_mex_create("y", c1_b_u, 10, 0U, 1U, 0U, 2, 1, 4),
                false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    2U, 14, c1_y, 14, c1_b_y));
}

static real_T c1_abs(SFc1_gen_recInstanceStruct *chartInstance, real_T c1_x)
{
  real_T c1_b_x;
  (void)chartInstance;
  c1_b_x = c1_x;
  return muDoubleScalarAbs(c1_b_x);
}

static real_T c1_angle(SFc1_gen_recInstanceStruct *chartInstance, creal_T c1_x)
{
  real_T c1_b_y;
  real_T c1_b_x;
  (void)chartInstance;
  c1_b_y = c1_x.im;
  c1_b_x = c1_x.re;
  return muDoubleScalarAtan2(c1_b_y, c1_b_x);
}

static real_T c1_mpower(SFc1_gen_recInstanceStruct *chartInstance, real_T c1_a)
{
  real_T c1_b_a;
  real_T c1_c_a;
  real_T c1_ak;
  real_T c1_d_a;
  c1_b_a = c1_a;
  c1_c_a = c1_b_a;
  c1_eml_scalar_eg(chartInstance);
  c1_ak = c1_c_a;
  c1_d_a = c1_ak;
  c1_eml_scalar_eg(chartInstance);
  return c1_d_a * c1_d_a;
}

static void c1_eml_scalar_eg(SFc1_gen_recInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *c1_f_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_u;
  const mxArray *c1_y = NULL;
  SFc1_gen_recInstanceStruct *chartInstance;
  chartInstance = (SFc1_gen_recInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_u = *(int32_T *)c1_inData;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", &c1_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static int32_T c1_f_emlrt_marshallIn(SFc1_gen_recInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  int32_T c1_y;
  int32_T c1_i26;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_i26, 1, 6, 0U, 0, 0U, 0);
  c1_y = c1_i26;
  sf_mex_destroy(&c1_u);
  return c1_y;
}

static void c1_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_b_sfEvent;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  int32_T c1_y;
  SFc1_gen_recInstanceStruct *chartInstance;
  chartInstance = (SFc1_gen_recInstanceStruct *)chartInstanceVoid;
  c1_b_sfEvent = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_sfEvent),
    &c1_thisId);
  sf_mex_destroy(&c1_b_sfEvent);
  *(int32_T *)c1_outData = c1_y;
  sf_mex_destroy(&c1_mxArrayInData);
}

static uint8_T c1_g_emlrt_marshallIn(SFc1_gen_recInstanceStruct *chartInstance,
  const mxArray *c1_b_is_active_c1_gen_rec, const char_T *c1_identifier)
{
  uint8_T c1_y;
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_h_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c1_b_is_active_c1_gen_rec), &c1_thisId);
  sf_mex_destroy(&c1_b_is_active_c1_gen_rec);
  return c1_y;
}

static uint8_T c1_h_emlrt_marshallIn(SFc1_gen_recInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  uint8_T c1_y;
  uint8_T c1_u0;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_u0, 1, 3, 0U, 0, 0U, 0);
  c1_y = c1_u0;
  sf_mex_destroy(&c1_u);
  return c1_y;
}

static void c1_b_sqrt(SFc1_gen_recInstanceStruct *chartInstance, real_T *c1_x)
{
  if (*c1_x < 0.0) {
    c1_eml_error(chartInstance);
  }

  *c1_x = muDoubleScalarSqrt(*c1_x);
}

static void init_dsm_address_info(SFc1_gen_recInstanceStruct *chartInstance)
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

void sf_c1_gen_rec_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(1978226281U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(2011028337U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3736234588U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(1423830567U);
}

mxArray *sf_c1_gen_rec_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("GzbU7uqtQatHTcMvXulHeH");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,7,3,dataFields);

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
      pr[0] = (double)(7);
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
      pr[0] = (double)(1);
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

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,5,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,5,"type",mxType);
    }

    mxSetField(mxData,5,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,6,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,6,"type",mxType);
    }

    mxSetField(mxData,6,"complexity",mxCreateDoubleScalar(0));
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
      pr[0] = (double)(1);
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

mxArray *sf_c1_gen_rec_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c1_gen_rec_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c1_gen_rec(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x6'type','srcId','name','auxInfo'{{M[1],M[15],T\"beta\",},{M[1],M[13],T\"ids\",},{M[1],M[5],T\"iqs\",},{M[1],M[14],T\"pidc\",},{M[1],M[4],T\"u\",},{M[8],M[0],T\"is_active_c1_gen_rec\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 6, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c1_gen_rec_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc1_gen_recInstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc1_gen_recInstanceStruct *) chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _gen_recMachineNumber_,
           1,
           1,
           1,
           0,
           13,
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
          _SFD_SET_DATA_PROPS(0,2,0,1,"iqs");
          _SFD_SET_DATA_PROPS(1,1,1,0,"wrm");
          _SFD_SET_DATA_PROPS(2,1,1,0,"betac");
          _SFD_SET_DATA_PROPS(3,1,1,0,"ed");
          _SFD_SET_DATA_PROPS(4,1,1,0,"Iamr");
          _SFD_SET_DATA_PROPS(5,1,1,0,"idc");
          _SFD_SET_DATA_PROPS(6,1,1,0,"betaminprev");
          _SFD_SET_DATA_PROPS(7,1,1,0,"uprev");
          _SFD_SET_DATA_PROPS(8,2,0,1,"ids");
          _SFD_SET_DATA_PROPS(9,2,0,1,"pidc");
          _SFD_SET_DATA_PROPS(10,2,0,1,"beta");
          _SFD_SET_DATA_PROPS(11,2,0,1,"u");
          _SFD_SET_DATA_PROPS(12,10,0,0,"P");
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
        _SFD_CV_INIT_EML(0,1,1,4,0,0,0,2,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,6730);
        _SFD_CV_INIT_EML_IF(0,1,0,2781,2792,3006,3055);
        _SFD_CV_INIT_EML_IF(0,1,1,3006,3023,-1,-2);
        _SFD_CV_INIT_EML_IF(0,1,2,3921,3932,4175,4217);
        _SFD_CV_INIT_EML_IF(0,1,3,4175,4192,-1,-2);
        _SFD_CV_INIT_EML_FOR(0,1,0,2668,2685,3083);
        _SFD_CV_INIT_EML_FOR(0,1,1,3788,3805,4245);
        _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)c1_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 7;
          _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_c_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(6,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(7,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(8,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)c1_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(9,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)c1_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(10,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)c1_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(11,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)c1_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(12,SF_STRUCT,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_b_sf_marshallOut,(MexInFcnForType)c1_b_sf_marshallIn);

        {
          real_T *c1_iqs;
          real_T *c1_wrm;
          real_T *c1_betac;
          real_T *c1_ed;
          real_T *c1_idc;
          real_T *c1_betaminprev;
          real_T *c1_uprev;
          real_T *c1_ids;
          real_T *c1_pidc;
          real_T *c1_beta;
          real_T *c1_u;
          real_T (*c1_Iamr)[7];
          c1_u = (real_T *)ssGetOutputPortSignal(chartInstance->S, 5);
          c1_beta = (real_T *)ssGetOutputPortSignal(chartInstance->S, 4);
          c1_pidc = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
          c1_ids = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
          c1_uprev = (real_T *)ssGetInputPortSignal(chartInstance->S, 6);
          c1_betaminprev = (real_T *)ssGetInputPortSignal(chartInstance->S, 5);
          c1_idc = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
          c1_Iamr = (real_T (*)[7])ssGetInputPortSignal(chartInstance->S, 3);
          c1_ed = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
          c1_betac = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
          c1_wrm = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
          c1_iqs = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
          _SFD_SET_DATA_VALUE_PTR(0U, c1_iqs);
          _SFD_SET_DATA_VALUE_PTR(1U, c1_wrm);
          _SFD_SET_DATA_VALUE_PTR(2U, c1_betac);
          _SFD_SET_DATA_VALUE_PTR(3U, c1_ed);
          _SFD_SET_DATA_VALUE_PTR(4U, *c1_Iamr);
          _SFD_SET_DATA_VALUE_PTR(5U, c1_idc);
          _SFD_SET_DATA_VALUE_PTR(6U, c1_betaminprev);
          _SFD_SET_DATA_VALUE_PTR(7U, c1_uprev);
          _SFD_SET_DATA_VALUE_PTR(8U, c1_ids);
          _SFD_SET_DATA_VALUE_PTR(9U, c1_pidc);
          _SFD_SET_DATA_VALUE_PTR(10U, c1_beta);
          _SFD_SET_DATA_VALUE_PTR(11U, c1_u);
          _SFD_SET_DATA_VALUE_PTR(12U, &chartInstance->c1_P);
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
  return "XSQCKlHavrd5fJ2siC9vQF";
}

static void sf_opaque_initialize_c1_gen_rec(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc1_gen_recInstanceStruct*) chartInstanceVar)->S,
    0);
  initialize_params_c1_gen_rec((SFc1_gen_recInstanceStruct*) chartInstanceVar);
  initialize_c1_gen_rec((SFc1_gen_recInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c1_gen_rec(void *chartInstanceVar)
{
  enable_c1_gen_rec((SFc1_gen_recInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c1_gen_rec(void *chartInstanceVar)
{
  disable_c1_gen_rec((SFc1_gen_recInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c1_gen_rec(void *chartInstanceVar)
{
  sf_gateway_c1_gen_rec((SFc1_gen_recInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c1_gen_rec(SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c1_gen_rec((SFc1_gen_recInstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c1_gen_rec();/* state var info */
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

extern void sf_internal_set_sim_state_c1_gen_rec(SimStruct* S, const mxArray *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[3];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxDuplicateArray(st);      /* high level simctx */
  prhs[2] = (mxArray*) sf_get_sim_state_info_c1_gen_rec();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 3, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c1_gen_rec((SFc1_gen_recInstanceStruct*)chartInfo->chartInstance,
    mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c1_gen_rec(SimStruct* S)
{
  return sf_internal_get_sim_state_c1_gen_rec(S);
}

static void sf_opaque_set_sim_state_c1_gen_rec(SimStruct* S, const mxArray *st)
{
  sf_internal_set_sim_state_c1_gen_rec(S, st);
}

static void sf_opaque_terminate_c1_gen_rec(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc1_gen_recInstanceStruct*) chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_gen_rec_optimization_info();
    }

    finalize_c1_gen_rec((SFc1_gen_recInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    if (crtInfo != NULL) {
      utFree((void *)crtInfo);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc1_gen_rec((SFc1_gen_recInstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c1_gen_rec(SimStruct *S)
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
    initialize_params_c1_gen_rec((SFc1_gen_recInstanceStruct*)
      (chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c1_gen_rec(SimStruct *S)
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
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,1);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,1,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,1,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,1);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 4, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 5, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 6, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,1,7);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,1,5);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=5; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 7; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,1);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(28344723U));
  ssSetChecksum1(S,(3163408129U));
  ssSetChecksum2(S,(2676470032U));
  ssSetChecksum3(S,(1863141359U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c1_gen_rec(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c1_gen_rec(SimStruct *S)
{
  SFc1_gen_recInstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc1_gen_recInstanceStruct *)utMalloc(sizeof
    (SFc1_gen_recInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc1_gen_recInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c1_gen_rec;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c1_gen_rec;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c1_gen_rec;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c1_gen_rec;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c1_gen_rec;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c1_gen_rec;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c1_gen_rec;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c1_gen_rec;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c1_gen_rec;
  chartInstance->chartInfo.mdlStart = mdlStart_c1_gen_rec;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c1_gen_rec;
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

void c1_gen_rec_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c1_gen_rec(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c1_gen_rec(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c1_gen_rec(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c1_gen_rec_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
