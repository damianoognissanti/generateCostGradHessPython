#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>
#include <algorithm>

#include "Bruno_JExpBot2016_p.h"
#include "Bruno_JExpBot2016_y.h"
#include "Bruno_JExpBot2016_sigmay.h"
#include "Bruno_JExpBot2016_my.h"
#include "Bruno_JExpBot2016_dJydy.h"

namespace amici {
namespace model_Bruno_JExpBot2016 {

void dJydy_Bruno_JExpBot2016(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*mob10 + 1.0*ob10)/std::pow(sigma_ob10, 2);
            break;
        case 1:
            dJydy[0] = (-1.0*mobcar + 1.0*obcar)/std::pow(sigma_obcar, 2);
            break;
        case 2:
            dJydy[0] = (-1.0*mobcry + 1.0*obcry)/std::pow(sigma_obcry, 2);
            break;
        case 3:
            dJydy[0] = (-1.0*mobio + 1.0*obio)/std::pow(sigma_obio, 2);
            break;
        case 4:
            dJydy[0] = (-1.0*moohb10 + 1.0*oohb10)/std::pow(sigma_oohb10, 2);
            break;
        case 5:
            dJydy[0] = (-1.0*mozea + 1.0*ozea)/std::pow(sigma_ozea, 2);
            break;
    }
}

} // namespace model_Bruno_JExpBot2016
} // namespace amici
