#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>
#include <algorithm>

#include "Brannmark_JBC2010_x.h"
#include "Brannmark_JBC2010_p.h"
#include "Brannmark_JBC2010_k.h"
#include "Brannmark_JBC2010_h.h"
#include "Brannmark_JBC2010_w.h"
#include "Brannmark_JBC2010_dxdotdw.h"

namespace amici {
namespace model_Brannmark_JBC2010 {

void dxdotdw_Brannmark_JBC2010(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dxdot0_dflux_v1_v_0 = -1.0;  // dxdotdw[0]
    dxdot1_dflux_v1_v_0 = 1.0;  // dxdotdw[1]
    dxdot0_dflux_v2_v_1 = 1.0;  // dxdotdw[2]
    dxdot1_dflux_v2_v_1 = -1.0;  // dxdotdw[3]
    dxdot1_dflux_v3_v_2 = -1.0;  // dxdotdw[4]
    dxdot2_dflux_v3_v_2 = 1.0;  // dxdotdw[5]
    dxdot2_dflux_v4_v_3 = -1.0;  // dxdotdw[6]
    dxdot3_dflux_v4_v_3 = 1.0;  // dxdotdw[7]
    dxdot3_dflux_v5_v_4 = -1.0;  // dxdotdw[8]
    dxdot4_dflux_v5_v_4 = 1.0;  // dxdotdw[9]
    dxdot0_dflux_v6_v_5 = 1.0;  // dxdotdw[10]
    dxdot2_dflux_v6_v_5 = -1.0;  // dxdotdw[11]
    dxdot0_dflux_v7_v_6 = 1.0;  // dxdotdw[12]
    dxdot4_dflux_v7_v_6 = -1.0;  // dxdotdw[13]
    dxdot5_dflux_v8_v_7 = -1.0;  // dxdotdw[14]
    dxdot6_dflux_v8_v_7 = 1.0;  // dxdotdw[15]
    dxdot5_dflux_v9_v_8 = 1.0;  // dxdotdw[16]
    dxdot6_dflux_v9_v_8 = -1.0;  // dxdotdw[17]
    dxdot7_dflux_v10_v_9 = -1.0;  // dxdotdw[18]
    dxdot8_dflux_v10_v_9 = 1.0;  // dxdotdw[19]
    dxdot7_dflux_v11_v_10 = 1.0;  // dxdotdw[20]
    dxdot8_dflux_v11_v_10 = -1.0;  // dxdotdw[21]
}

} // namespace model_Brannmark_JBC2010
} // namespace amici
