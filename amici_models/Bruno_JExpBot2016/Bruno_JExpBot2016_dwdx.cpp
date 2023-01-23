#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>
#include <algorithm>

#include "Bruno_JExpBot2016_x.h"
#include "Bruno_JExpBot2016_p.h"
#include "Bruno_JExpBot2016_w.h"
#include "Bruno_JExpBot2016_dwdx.h"

namespace amici {
namespace model_Bruno_JExpBot2016 {

void dwdx_Bruno_JExpBot2016(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dflux_v1_ReactionName_dbcar = 1.0*kb1*kb1_multiplier;  // dwdx[0]
    dflux_v3_ReactionName_dbcry = 1.0*kc1*kc1_multiplier;  // dwdx[1]
    dflux_v4_ReactionName_dbcry = 1.0*kc2*kc2_multiplier;  // dwdx[2]
    dflux_v2_ReactionName_db10 = 1.0*kb2*kb2_multiplier;  // dwdx[3]
    dflux_v5_ReactionName_dohb10 = 1.0*kc4*kc4_multiplier;  // dwdx[4]
    dflux_v6_ReactionName_dzea = 1.0*k5*k5_multiplier;  // dwdx[5]
}

} // namespace model_Bruno_JExpBot2016
} // namespace amici
