#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>
#include <algorithm>

#include "Bruno_JExpBot2016_p.h"

namespace amici {
namespace model_Bruno_JExpBot2016 {

void x0_Bruno_JExpBot2016(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = init_bcar;
    x0[1] = init_bcry;
    x0[2] = init_b10;
    x0[4] = init_ohb10;
    x0[6] = init_zea;
}

} // namespace model_Bruno_JExpBot2016
} // namespace amici
