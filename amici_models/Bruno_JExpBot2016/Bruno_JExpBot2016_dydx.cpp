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

void dydx_Bruno_JExpBot2016(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[1] = 1;
    dydx[8] = 1;
    dydx[12] = 1;
    dydx[21] = 1;
    dydx[28] = 1;
    dydx[41] = 1;
}

} // namespace model_Bruno_JExpBot2016
} // namespace amici
