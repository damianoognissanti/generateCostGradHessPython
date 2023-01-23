#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "Boehm_JProteomeRes2014_p.h"
#include "Boehm_JProteomeRes2014_k.h"
#include "Boehm_JProteomeRes2014_y.h"

namespace amici {
namespace model_Boehm_JProteomeRes2014 {

void dsigmaydp_Boehm_JProteomeRes2014(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip){
    switch(ip) {
        case 6:
            dsigmaydp[0] = 1;
            break;
        case 7:
            dsigmaydp[1] = 1;
            break;
        case 8:
            dsigmaydp[2] = 1;
            break;
    }
}

} // namespace model_Boehm_JProteomeRes2014
} // namespace amici
