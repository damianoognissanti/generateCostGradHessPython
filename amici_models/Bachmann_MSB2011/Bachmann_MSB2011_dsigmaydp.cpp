#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "Bachmann_MSB2011_p.h"
#include "Bachmann_MSB2011_k.h"
#include "Bachmann_MSB2011_y.h"

namespace amici {
namespace model_Bachmann_MSB2011 {

void dsigmaydp_Bachmann_MSB2011(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip){
    switch(ip) {
        case 48:
            dsigmaydp[0] = 1;
            break;
        case 49:
            dsigmaydp[1] = 1;
            break;
        case 50:
            dsigmaydp[2] = 1;
            break;
        case 51:
            dsigmaydp[3] = 1;
            break;
        case 52:
            dsigmaydp[4] = 1;
            break;
        case 53:
            dsigmaydp[5] = 1;
            break;
        case 54:
            dsigmaydp[6] = 1;
            break;
        case 55:
            dsigmaydp[7] = 1;
            break;
        case 56:
            dsigmaydp[8] = 1;
            break;
        case 57:
            dsigmaydp[9] = 1;
            break;
        case 58:
            dsigmaydp[10] = 1;
            break;
        case 59:
            dsigmaydp[11] = 1;
            break;
        case 60:
            dsigmaydp[12] = 1;
            break;
        case 61:
            dsigmaydp[13] = 1;
            break;
        case 62:
            dsigmaydp[14] = 1;
            break;
        case 63:
            dsigmaydp[15] = 1;
            break;
        case 64:
            dsigmaydp[16] = 1;
            break;
        case 65:
            dsigmaydp[17] = 1;
            break;
        case 66:
            dsigmaydp[17] = 1;
            break;
        case 67:
            dsigmaydp[18] = 1;
            break;
        case 68:
            dsigmaydp[19] = 1;
            break;
    }
}

} // namespace model_Bachmann_MSB2011
} // namespace amici
