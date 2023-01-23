#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>
#include <algorithm>

#include "Bruno_JExpBot2016_x_rdata.h"

namespace amici {
namespace model_Bruno_JExpBot2016 {

void x_solver_Bruno_JExpBot2016(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = bcar;
    x_solver[1] = bcry;
    x_solver[2] = b10;
    x_solver[3] = bio;
    x_solver[4] = ohb10;
    x_solver[5] = ohbio;
    x_solver[6] = zea;
}

} // namespace model_Bruno_JExpBot2016
} // namespace amici
