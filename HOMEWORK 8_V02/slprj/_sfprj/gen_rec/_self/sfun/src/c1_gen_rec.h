#ifndef __c1_gen_rec_h__
#define __c1_gen_rec_h__

/* Include files */
#include "sf_runtime/sfc_sf.h"
#include "sf_runtime/sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef struct_struct_U3Zfiu1q58uXf1A0Qkd5GF_tag
#define struct_struct_U3Zfiu1q58uXf1A0Qkd5GF_tag

struct struct_U3Zfiu1q58uXf1A0Qkd5GF_tag
{
  real_T Ll;
  real_T rl;
  real_T rs;
  real_T Lls;
  real_T Lmq;
  real_T Lmd;
  real_T rd1;
  real_T rd2;
  real_T rd3;
  real_T rfd;
  real_T Lld1;
  real_T Lld2;
  real_T Lld3;
  real_T Llfd;
  real_T rq1;
  real_T rq2;
  real_T rq3;
  real_T pole;
  real_T Llq1;
  real_T Llq2;
  real_T Llq3;
};

#endif                                 /*struct_struct_U3Zfiu1q58uXf1A0Qkd5GF_tag*/

#ifndef typedef_c1_struct_U3Zfiu1q58uXf1A0Qkd5GF
#define typedef_c1_struct_U3Zfiu1q58uXf1A0Qkd5GF

typedef struct struct_U3Zfiu1q58uXf1A0Qkd5GF_tag
  c1_struct_U3Zfiu1q58uXf1A0Qkd5GF;

#endif                                 /*typedef_c1_struct_U3Zfiu1q58uXf1A0Qkd5GF*/

#ifndef typedef_SFc1_gen_recInstanceStruct
#define typedef_SFc1_gen_recInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c1_sfEvent;
  boolean_T c1_isStable;
  boolean_T c1_doneDoubleBufferReInit;
  uint8_T c1_is_active_c1_gen_rec;
  c1_struct_U3Zfiu1q58uXf1A0Qkd5GF c1_P;
} SFc1_gen_recInstanceStruct;

#endif                                 /*typedef_SFc1_gen_recInstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c1_gen_rec_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c1_gen_rec_get_check_sum(mxArray *plhs[]);
extern void c1_gen_rec_method_dispatcher(SimStruct *S, int_T method, void *data);

#endif
