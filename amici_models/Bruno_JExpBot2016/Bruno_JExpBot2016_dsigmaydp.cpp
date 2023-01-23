#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>
#include <algorithm>

#include "Bruno_JExpBot2016_p.h"
#include "Bruno_JExpBot2016_y.h"

namespace amici {
namespace model_Bruno_JExpBot2016 {

void dsigmaydp_Bruno_JExpBot2016(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip){
    switch(ip) {
        case 17:
            dsigmaydp[0] = 1;
            break;
        case 18:
            dsigmaydp[1] = 1;
            break;
        case 19:
            dsigmaydp[2] = 1;
            break;
        case 20:
            dsigmaydp[3] = 1;
            break;
        case 21:
            dsigmaydp[4] = 1;
            break;
        case 22:
            dsigmaydp[5] = 1;
            break;
    }
}

} // namespace model_Bruno_JExpBot2016
} // namespace amici
