#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>
#include <algorithm>

#include "Bruno_JExpBot2016_x.h"
#include "Bruno_JExpBot2016_p.h"

namespace amici {
namespace model_Bruno_JExpBot2016 {

void sx0_Bruno_JExpBot2016(realtype *sx0, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip){
    switch(ip) {
        case 0:
            sx0[2] = 1;
            break;
        case 1:
            sx0[0] = 1;
            break;
        case 2:
            sx0[1] = 1;
            break;
        case 3:
            sx0[4] = 1;
            break;
        case 4:
            sx0[6] = 1;
            break;
    }
}

} // namespace model_Bruno_JExpBot2016
} // namespace amici
