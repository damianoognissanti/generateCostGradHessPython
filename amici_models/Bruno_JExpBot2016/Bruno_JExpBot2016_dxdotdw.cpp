#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>
#include <algorithm>

#include "Bruno_JExpBot2016_x.h"
#include "Bruno_JExpBot2016_p.h"
#include "Bruno_JExpBot2016_w.h"
#include "Bruno_JExpBot2016_dxdotdw.h"

namespace amici {
namespace model_Bruno_JExpBot2016 {

void dxdotdw_Bruno_JExpBot2016(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dxdot0_dflux_v1_ReactionName = -1.0;  // dxdotdw[0]
    dxdot2_dflux_v1_ReactionName = 1.0;  // dxdotdw[1]
    dxdot3_dflux_v1_ReactionName = 1.0;  // dxdotdw[2]
    dxdot2_dflux_v2_ReactionName = -1.0;  // dxdotdw[3]
    dxdot3_dflux_v2_ReactionName = 1.0;  // dxdotdw[4]
    dxdot1_dflux_v3_ReactionName = -1.0;  // dxdotdw[5]
    dxdot2_dflux_v3_ReactionName = 1.0;  // dxdotdw[6]
    dxdot5_dflux_v3_ReactionName = 1.0;  // dxdotdw[7]
    dxdot1_dflux_v4_ReactionName = -1.0;  // dxdotdw[8]
    dxdot3_dflux_v4_ReactionName = 1.0;  // dxdotdw[9]
    dxdot4_dflux_v4_ReactionName = 1.0;  // dxdotdw[10]
    dxdot4_dflux_v5_ReactionName = -1.0;  // dxdotdw[11]
    dxdot5_dflux_v5_ReactionName = 1.0;  // dxdotdw[12]
    dxdot4_dflux_v6_ReactionName = 1.0;  // dxdotdw[13]
    dxdot5_dflux_v6_ReactionName = 1.0;  // dxdotdw[14]
    dxdot6_dflux_v6_ReactionName = -1.0;  // dxdotdw[15]
}

} // namespace model_Bruno_JExpBot2016
} // namespace amici
