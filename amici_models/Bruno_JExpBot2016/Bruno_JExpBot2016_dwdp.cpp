#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>
#include <algorithm>

#include "Bruno_JExpBot2016_x.h"
#include "Bruno_JExpBot2016_p.h"
#include "Bruno_JExpBot2016_w.h"
#include "Bruno_JExpBot2016_dwdp.h"

namespace amici {
namespace model_Bruno_JExpBot2016 {

void dwdp_Bruno_JExpBot2016(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp){
    dflux_v6_ReactionName_dk5 = 1.0*k5_multiplier*zea;  // dwdp[0]
    dflux_v1_ReactionName_dkb1 = 1.0*bcar*kb1_multiplier;  // dwdp[1]
    dflux_v2_ReactionName_dkb2 = 1.0*b10*kb2_multiplier;  // dwdp[2]
    dflux_v3_ReactionName_dkc1 = 1.0*bcry*kc1_multiplier;  // dwdp[3]
    dflux_v4_ReactionName_dkc2 = 1.0*bcry*kc2_multiplier;  // dwdp[4]
    dflux_v5_ReactionName_dkc4 = 1.0*kc4_multiplier*ohb10;  // dwdp[5]
    dflux_v6_ReactionName_dk5_multiplier = 1.0*k5*zea;  // dwdp[6]
    dflux_v1_ReactionName_dkb1_multiplier = 1.0*bcar*kb1;  // dwdp[7]
    dflux_v2_ReactionName_dkb2_multiplier = 1.0*b10*kb2;  // dwdp[8]
    dflux_v3_ReactionName_dkc1_multiplier = 1.0*bcry*kc1;  // dwdp[9]
    dflux_v4_ReactionName_dkc2_multiplier = 1.0*bcry*kc2;  // dwdp[10]
    dflux_v5_ReactionName_dkc4_multiplier = 1.0*kc4*ohb10;  // dwdp[11]
}

} // namespace model_Bruno_JExpBot2016
} // namespace amici
