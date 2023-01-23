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
#include "Brannmark_JBC2010_xdot.h"

namespace amici {
namespace model_Brannmark_JBC2010 {

void xdot_Brannmark_JBC2010(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    xdot0 = -1.0*flux_v1_v_0 + 1.0*flux_v2_v_1 + 1.0*flux_v6_v_5 + 1.0*flux_v7_v_6;  // xdot[0]
    xdot1 = 1.0*flux_v1_v_0 - 1.0*flux_v2_v_1 - 1.0*flux_v3_v_2;  // xdot[1]
    xdot2 = 1.0*flux_v3_v_2 - 1.0*flux_v4_v_3 - 1.0*flux_v6_v_5;  // xdot[2]
    xdot3 = 1.0*flux_v4_v_3 - 1.0*flux_v5_v_4;  // xdot[3]
    xdot4 = 1.0*flux_v5_v_4 - 1.0*flux_v7_v_6;  // xdot[4]
    xdot5 = -1.0*flux_v8_v_7 + 1.0*flux_v9_v_8;  // xdot[5]
    xdot6 = 1.0*flux_v8_v_7 - 1.0*flux_v9_v_8;  // xdot[6]
    xdot7 = -1.0*flux_v10_v_9 + 1.0*flux_v11_v_10;  // xdot[7]
    xdot8 = 1.0*flux_v10_v_9 - 1.0*flux_v11_v_10;  // xdot[8]
}

} // namespace model_Brannmark_JBC2010
} // namespace amici
