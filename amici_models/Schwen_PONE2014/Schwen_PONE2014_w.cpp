#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>
#include <algorithm>

#include "Schwen_PONE2014_x.h"
#include "Schwen_PONE2014_p.h"
#include "Schwen_PONE2014_k.h"
#include "Schwen_PONE2014_w.h"

namespace amici {
namespace model_Schwen_PONE2014 {

void w_Schwen_PONE2014(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    flux_v1_ka1 = 1.0*Ins*Rec1*ka1;  // w[0]
    flux_v2_ka2fold = 1.0*Ins*Rec2*ka1*ka2fold;  // w[1]
    flux_v3_v_2 = 1.0*Ins*kon_unspec;  // w[2]
    flux_v4_v_3 = 1.0*BoundUnspec*koff_unspec;  // w[3]
    flux_v5_kd1 = 1.0*IR1*kd1;  // w[4]
    flux_v6_kd2fold = 1.0*IR2*kd1*kd2fold;  // w[5]
    flux_v7_v_6 = 1.0*IR1*kin;  // w[6]
    flux_v8_v_7 = 1.0*IR2*kin2;  // w[7]
    flux_v9_v_8 = 1.0*IR1in*kout;  // w[8]
    flux_v10_v_9 = 1.0*IR2in*kout2;  // w[9]
    flux_v11_v_10 = 1.0*IR1in*kout_frag;  // w[10]
    flux_v12_v_11 = 1.0*IR2in*kout_frag;  // w[11]
    flux_v13_v_12 = -1.0*IR1*kd1 + 1.0*Ins*Rec1*ka1;  // w[12]
    flux_v14_v_13 = -1.0*IR2*kd1*kd2fold + 1.0*Ins*Rec2*ka1*ka2fold;  // w[13]
}

} // namespace model_Schwen_PONE2014
} // namespace amici
