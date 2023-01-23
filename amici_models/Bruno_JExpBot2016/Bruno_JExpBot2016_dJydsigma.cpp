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

namespace amici {
namespace model_Bruno_JExpBot2016 {

void dJydsigma_Bruno_JExpBot2016(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_ob10 - 1.0*std::pow(-mob10 + ob10, 2)/std::pow(sigma_ob10, 3);
            break;
        case 1:
            dJydsigma[1] = 1.0/sigma_obcar - 1.0*std::pow(-mobcar + obcar, 2)/std::pow(sigma_obcar, 3);
            break;
        case 2:
            dJydsigma[2] = 1.0/sigma_obcry - 1.0*std::pow(-mobcry + obcry, 2)/std::pow(sigma_obcry, 3);
            break;
        case 3:
            dJydsigma[3] = 1.0/sigma_obio - 1.0*std::pow(-mobio + obio, 2)/std::pow(sigma_obio, 3);
            break;
        case 4:
            dJydsigma[4] = 1.0/sigma_oohb10 - 1.0*std::pow(-moohb10 + oohb10, 2)/std::pow(sigma_oohb10, 3);
            break;
        case 5:
            dJydsigma[5] = 1.0/sigma_ozea - 1.0*std::pow(-mozea + ozea, 2)/std::pow(sigma_ozea, 3);
            break;
    }
}

} // namespace model_Bruno_JExpBot2016
} // namespace amici
