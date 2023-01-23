#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>
#include <algorithm>

#include "Bruno_JExpBot2016_x.h"
#include "Bruno_JExpBot2016_p.h"
#include "Bruno_JExpBot2016_w.h"

namespace amici {
namespace model_Bruno_JExpBot2016 {

void w_Bruno_JExpBot2016(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    flux_v1_ReactionName = 1.0*bcar*kb1*kb1_multiplier;  // w[0]
    flux_v2_ReactionName = 1.0*b10*kb2*kb2_multiplier;  // w[1]
    flux_v3_ReactionName = 1.0*bcry*kc1*kc1_multiplier;  // w[2]
    flux_v4_ReactionName = 1.0*bcry*kc2*kc2_multiplier;  // w[3]
    flux_v5_ReactionName = 1.0*kc4*kc4_multiplier*ohb10;  // w[4]
    flux_v6_ReactionName = 1.0*k5*k5_multiplier*zea;  // w[5]
}

} // namespace model_Bruno_JExpBot2016
} // namespace amici
