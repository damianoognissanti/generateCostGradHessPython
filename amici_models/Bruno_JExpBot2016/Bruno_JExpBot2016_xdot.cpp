#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>
#include <algorithm>

#include "Bruno_JExpBot2016_x.h"
#include "Bruno_JExpBot2016_p.h"
#include "Bruno_JExpBot2016_w.h"
#include "Bruno_JExpBot2016_xdot.h"

namespace amici {
namespace model_Bruno_JExpBot2016 {

void xdot_Bruno_JExpBot2016(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    xdot0 = -1.0*flux_v1_ReactionName;  // xdot[0]
    xdot1 = -1.0*flux_v3_ReactionName - 1.0*flux_v4_ReactionName;  // xdot[1]
    xdot2 = 1.0*flux_v1_ReactionName - 1.0*flux_v2_ReactionName + 1.0*flux_v3_ReactionName;  // xdot[2]
    xdot3 = 1.0*flux_v1_ReactionName + 1.0*flux_v2_ReactionName + 1.0*flux_v4_ReactionName;  // xdot[3]
    xdot4 = 1.0*flux_v4_ReactionName - 1.0*flux_v5_ReactionName + 1.0*flux_v6_ReactionName;  // xdot[4]
    xdot5 = 1.0*flux_v3_ReactionName + 1.0*flux_v5_ReactionName + 1.0*flux_v6_ReactionName;  // xdot[5]
    xdot6 = -1.0*flux_v6_ReactionName;  // xdot[6]
}

} // namespace model_Bruno_JExpBot2016
} // namespace amici
