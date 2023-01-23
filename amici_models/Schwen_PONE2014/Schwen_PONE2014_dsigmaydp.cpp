#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>
#include <algorithm>

#include "Schwen_PONE2014_p.h"
#include "Schwen_PONE2014_k.h"
#include "Schwen_PONE2014_y.h"

namespace amici {
namespace model_Schwen_PONE2014 {

void dsigmaydp_Schwen_PONE2014(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip){
    switch(ip) {
        case 23:
            dsigmaydp[0] = 1;
            break;
        case 24:
            dsigmaydp[1] = 1;
            break;
        case 25:
            dsigmaydp[2] = 1;
            break;
        case 26:
            dsigmaydp[3] = 1;
            break;
    }
}

} // namespace model_Schwen_PONE2014
} // namespace amici
