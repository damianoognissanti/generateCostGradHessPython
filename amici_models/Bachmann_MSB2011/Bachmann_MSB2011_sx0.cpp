#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "Bachmann_MSB2011_x.h"
#include "Bachmann_MSB2011_p.h"
#include "Bachmann_MSB2011_k.h"

namespace amici {
namespace model_Bachmann_MSB2011 {

void sx0_Bachmann_MSB2011(realtype *sx0, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip){
    switch(ip) {
        case 0:
            sx0[17] = CISEqcOE*init_CIS_multiplier;
            break;
        case 5:
            sx0[17] = CISEqc*init_CIS_multiplier;
            break;
        case 13:
            sx0[6] = init_SHP1*init_SHP1_multiplier;
            break;
        case 14:
            sx0[24] = SOCS3Eqc*init_SOCS3_multiplier;
            break;
        case 15:
            sx0[24] = SOCS3EqcOE*init_SOCS3_multiplier;
            break;
        case 24:
            sx0[0] = 1;
            break;
        case 25:
            sx0[6] = SHP1ProOE*init_SHP1_multiplier + 1;
            break;
        case 26:
            sx0[8] = 1;
            break;
    }
}

} // namespace model_Bachmann_MSB2011
} // namespace amici
