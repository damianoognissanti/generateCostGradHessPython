#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>
#include <algorithm>

#include "Bruno_JExpBot2016_p.h"
#include "Bruno_JExpBot2016_y.h"
#include "Bruno_JExpBot2016_sigmay.h"

namespace amici {
namespace model_Bruno_JExpBot2016 {

void sigmay_Bruno_JExpBot2016(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_ob10 = noiseParameter1_ob10;  // sigmay[0]
    sigma_obcar = noiseParameter1_obcar;  // sigmay[1]
    sigma_obcry = noiseParameter1_obcry;  // sigmay[2]
    sigma_obio = noiseParameter1_obio;  // sigmay[3]
    sigma_oohb10 = noiseParameter1_oohb10;  // sigmay[4]
    sigma_ozea = noiseParameter1_ozea;  // sigmay[5]
}

} // namespace model_Bruno_JExpBot2016
} // namespace amici
