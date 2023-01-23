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

namespace amici {
namespace model_Brannmark_JBC2010 {

void w_Brannmark_JBC2010(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    insulin = 1.0*Heaviside_0*insulin_dose_1 + 1.0*Heaviside_2*insulin_dose_2;  // w[0]
    flux_v1_v_0 = 1.0*IR*insulin*k1a + 1.0*IR*k1aBasic;  // w[1]
    flux_v2_v_1 = 1.0*IRins*k1b;  // w[2]
    flux_v3_v_2 = 1.0*IRins*k1c;  // w[3]
    flux_v4_v_3 = 1.0*IRp*k1d;  // w[4]
    flux_v5_v_4 = 1.0*IRiP*(Xp*k1f/(Xp + 1) + k1e);  // w[5]
    flux_v6_v_5 = 1.0*IRp*k1g;  // w[6]
    flux_v7_v_6 = 1.0*IRi*k1r;  // w[7]
    flux_v8_v_7 = 1.0*IRS*k21*(IRiP*k22 + IRp);  // w[8]
    flux_v9_v_8 = 1.0*IRSiP*km2;  // w[9]
    flux_v10_v_9 = 1.0*IRSiP*X*k3;  // w[10]
    flux_v11_v_10 = 1.0*Xp*km3;  // w[11]
}

} // namespace model_Brannmark_JBC2010
} // namespace amici
