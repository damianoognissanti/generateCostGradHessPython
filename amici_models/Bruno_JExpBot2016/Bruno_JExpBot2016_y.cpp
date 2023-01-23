#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>
#include <algorithm>

#include "Bruno_JExpBot2016_x.h"
#include "Bruno_JExpBot2016_p.h"
#include "Bruno_JExpBot2016_w.h"

namespace amici {
namespace model_Bruno_JExpBot2016 {

void y_Bruno_JExpBot2016(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = b10;
    y[1] = bcar;
    y[2] = bcry;
    y[3] = bio;
    y[4] = ohb10;
    y[5] = zea;
}

} // namespace model_Bruno_JExpBot2016
} // namespace amici
