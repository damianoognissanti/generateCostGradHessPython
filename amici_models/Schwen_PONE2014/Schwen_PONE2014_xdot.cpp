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
#include "Schwen_PONE2014_xdot.h"

namespace amici {
namespace model_Schwen_PONE2014 {

void xdot_Schwen_PONE2014(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    xdot0 = -1.0*flux_v1_ka1 - 1.0*flux_v2_ka2fold - 1.0*flux_v3_v_2 + 1.0*flux_v4_v_3 + 1.0*flux_v5_kd1 + 1.0*flux_v6_kd2fold;  // xdot[0]
    xdot1 = 1.0*flux_v11_v_10 - 1.0*flux_v1_ka1 + 1.0*flux_v5_kd1;  // xdot[1]
    xdot2 = 1.0*flux_v12_v_11 - 1.0*flux_v2_ka2fold + 1.0*flux_v6_kd2fold;  // xdot[2]
    xdot3 = 1.0*flux_v1_ka1 - 1.0*flux_v5_kd1 - 1.0*flux_v7_v_6 + 1.0*flux_v9_v_8;  // xdot[3]
    xdot4 = 1.0*flux_v10_v_9 + 1.0*flux_v2_ka2fold - 1.0*flux_v6_kd2fold - 1.0*flux_v8_v_7;  // xdot[4]
    xdot5 = -1.0*flux_v11_v_10 + 1.0*flux_v7_v_6 - 1.0*flux_v9_v_8;  // xdot[5]
    xdot6 = -1.0*flux_v10_v_9 - 1.0*flux_v12_v_11 + 1.0*flux_v8_v_7;  // xdot[6]
    xdot7 = 1.0*flux_v13_v_12;  // xdot[7]
    xdot8 = 1.0*flux_v14_v_13;  // xdot[8]
    xdot9 = 1.0*flux_v11_v_10 + 1.0*flux_v12_v_11;  // xdot[9]
    xdot10 = 1.0*flux_v3_v_2 - 1.0*flux_v4_v_3;  // xdot[10]
}

} // namespace model_Schwen_PONE2014
} // namespace amici
