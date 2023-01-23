#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>
#include <algorithm>

#include "Brannmark_JBC2010_p.h"
#include "Brannmark_JBC2010_k.h"
#include "Brannmark_JBC2010_y.h"

namespace amici {
namespace model_Brannmark_JBC2010 {

void dsigmaydp_Brannmark_JBC2010(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip){
    switch(ip) {
        case 18:
            dsigmaydp[0] = 1;
            break;
        case 19:
            dsigmaydp[1] = 1;
            break;
        case 20:
            dsigmaydp[2] = 1;
            break;
    }
}

} // namespace model_Brannmark_JBC2010
} // namespace amici
