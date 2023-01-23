#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>
#include <algorithm>

#include "Schwen_PONE2014_x.h"
#include "Schwen_PONE2014_p.h"
#include "Schwen_PONE2014_k.h"

namespace amici {
namespace model_Schwen_PONE2014 {

void sx0_Schwen_PONE2014(realtype *sx0, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip){
    switch(ip) {
        case 0:
            sx0[1] = 1;
            sx0[2] = ini_R2fold;
            break;
        case 1:
            sx0[2] = ini_R1;
            break;
    }
}

} // namespace model_Schwen_PONE2014
} // namespace amici
