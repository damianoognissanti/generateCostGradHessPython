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

void Jy_Bruno_JExpBot2016(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ob10, 2)) + 0.5*std::pow(-mob10 + ob10, 2)/std::pow(sigma_ob10, 2);
            break;
        case 1:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obcar, 2)) + 0.5*std::pow(-mobcar + obcar, 2)/std::pow(sigma_obcar, 2);
            break;
        case 2:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obcry, 2)) + 0.5*std::pow(-mobcry + obcry, 2)/std::pow(sigma_obcry, 2);
            break;
        case 3:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obio, 2)) + 0.5*std::pow(-mobio + obio, 2)/std::pow(sigma_obio, 2);
            break;
        case 4:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_oohb10, 2)) + 0.5*std::pow(-moohb10 + oohb10, 2)/std::pow(sigma_oohb10, 2);
            break;
        case 5:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ozea, 2)) + 0.5*std::pow(-mozea + ozea, 2)/std::pow(sigma_ozea, 2);
            break;
    }
}

} // namespace model_Bruno_JExpBot2016
} // namespace amici
