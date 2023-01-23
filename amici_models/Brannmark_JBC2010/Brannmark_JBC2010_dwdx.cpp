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
#include "Brannmark_JBC2010_dwdx.h"

namespace amici {
namespace model_Brannmark_JBC2010 {

void dwdx_Brannmark_JBC2010(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dflux_v1_v_0_dIR = 1.0*insulin*k1a + 1.0*k1aBasic;  // dwdx[0]
    dflux_v2_v_1_dIRins = 1.0*k1b;  // dwdx[1]
    dflux_v3_v_2_dIRins = 1.0*k1c;  // dwdx[2]
    dflux_v4_v_3_dIRp = 1.0*k1d;  // dwdx[3]
    dflux_v6_v_5_dIRp = 1.0*k1g;  // dwdx[4]
    dflux_v8_v_7_dIRp = 1.0*IRS*k21;  // dwdx[5]
    dflux_v5_v_4_dIRiP = 1.0*Xp*k1f/(Xp + 1) + 1.0*k1e;  // dwdx[6]
    dflux_v8_v_7_dIRiP = 1.0*IRS*k21*k22;  // dwdx[7]
    dflux_v7_v_6_dIRi = 1.0*k1r;  // dwdx[8]
    dflux_v8_v_7_dIRS = 1.0*k21*(IRiP*k22 + IRp);  // dwdx[9]
    dflux_v9_v_8_dIRSiP = 1.0*km2;  // dwdx[10]
    dflux_v10_v_9_dIRSiP = 1.0*X*k3;  // dwdx[11]
    dflux_v10_v_9_dX = 1.0*IRSiP*k3;  // dwdx[12]
    dflux_v5_v_4_dXp = 1.0*IRiP*(-Xp*k1f/std::pow(Xp + 1, 2) + k1f/(Xp + 1));  // dwdx[13]
    dflux_v11_v_10_dXp = 1.0*km3;  // dwdx[14]
}

} // namespace model_Brannmark_JBC2010
} // namespace amici
