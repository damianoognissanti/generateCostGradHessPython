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
#include "Schwen_PONE2014_dxdotdw.h"

namespace amici {
namespace model_Schwen_PONE2014 {

void dxdotdw_Schwen_PONE2014(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dxdot0_dflux_v1_ka1 = -1.0;  // dxdotdw[0]
    dxdot1_dflux_v1_ka1 = -1.0;  // dxdotdw[1]
    dxdot3_dflux_v1_ka1 = 1.0;  // dxdotdw[2]
    dxdot0_dflux_v2_ka2fold = -1.0;  // dxdotdw[3]
    dxdot2_dflux_v2_ka2fold = -1.0;  // dxdotdw[4]
    dxdot4_dflux_v2_ka2fold = 1.0;  // dxdotdw[5]
    dxdot0_dflux_v3_v_2 = -1.0;  // dxdotdw[6]
    dxdot10_dflux_v3_v_2 = 1.0;  // dxdotdw[7]
    dxdot0_dflux_v4_v_3 = 1.0;  // dxdotdw[8]
    dxdot10_dflux_v4_v_3 = -1.0;  // dxdotdw[9]
    dxdot0_dflux_v5_kd1 = 1.0;  // dxdotdw[10]
    dxdot1_dflux_v5_kd1 = 1.0;  // dxdotdw[11]
    dxdot3_dflux_v5_kd1 = -1.0;  // dxdotdw[12]
    dxdot0_dflux_v6_kd2fold = 1.0;  // dxdotdw[13]
    dxdot2_dflux_v6_kd2fold = 1.0;  // dxdotdw[14]
    dxdot4_dflux_v6_kd2fold = -1.0;  // dxdotdw[15]
    dxdot3_dflux_v7_v_6 = -1.0;  // dxdotdw[16]
    dxdot5_dflux_v7_v_6 = 1.0;  // dxdotdw[17]
    dxdot4_dflux_v8_v_7 = -1.0;  // dxdotdw[18]
    dxdot6_dflux_v8_v_7 = 1.0;  // dxdotdw[19]
    dxdot3_dflux_v9_v_8 = 1.0;  // dxdotdw[20]
    dxdot5_dflux_v9_v_8 = -1.0;  // dxdotdw[21]
    dxdot4_dflux_v10_v_9 = 1.0;  // dxdotdw[22]
    dxdot6_dflux_v10_v_9 = -1.0;  // dxdotdw[23]
    dxdot1_dflux_v11_v_10 = 1.0;  // dxdotdw[24]
    dxdot5_dflux_v11_v_10 = -1.0;  // dxdotdw[25]
    dxdot9_dflux_v11_v_10 = 1.0;  // dxdotdw[26]
    dxdot2_dflux_v12_v_11 = 1.0;  // dxdotdw[27]
    dxdot6_dflux_v12_v_11 = -1.0;  // dxdotdw[28]
    dxdot9_dflux_v12_v_11 = 1.0;  // dxdotdw[29]
    dxdot7_dflux_v13_v_12 = 1.0;  // dxdotdw[30]
    dxdot8_dflux_v14_v_13 = 1.0;  // dxdotdw[31]
}

} // namespace model_Schwen_PONE2014
} // namespace amici
