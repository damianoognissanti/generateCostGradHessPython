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
#include "Brannmark_JBC2010_dwdp.h"

namespace amici {
namespace model_Brannmark_JBC2010 {

void dwdp_Brannmark_JBC2010(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp){
    dflux_v1_v_0_dk1a = 1.0*IR*insulin;  // dwdp[0]
    dflux_v1_v_0_dk1aBasic = 1.0*IR;  // dwdp[1]
    dflux_v2_v_1_dk1b = 1.0*IRins;  // dwdp[2]
    dflux_v3_v_2_dk1c = 1.0*IRins;  // dwdp[3]
    dflux_v4_v_3_dk1d = 1.0*IRp;  // dwdp[4]
    dflux_v5_v_4_dk1e = 1.0*IRiP;  // dwdp[5]
    dflux_v5_v_4_dk1f = 1.0*IRiP*Xp/(Xp + 1);  // dwdp[6]
    dflux_v6_v_5_dk1g = 1.0*IRp;  // dwdp[7]
    dflux_v7_v_6_dk1r = 1.0*IRi;  // dwdp[8]
    dflux_v8_v_7_dk21 = 1.0*IRS*(IRiP*k22 + IRp);  // dwdp[9]
    dflux_v8_v_7_dk22 = 1.0*IRS*IRiP*k21;  // dwdp[10]
    dflux_v10_v_9_dk3 = 1.0*IRSiP*X;  // dwdp[11]
    dflux_v9_v_8_dkm2 = 1.0*IRSiP;  // dwdp[12]
    dflux_v11_v_10_dkm3 = 1.0*Xp;  // dwdp[13]
}

} // namespace model_Brannmark_JBC2010
} // namespace amici
