#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>
#include <algorithm>

#include "Schwen_PONE2014_p.h"
#include "Schwen_PONE2014_k.h"

namespace amici {
namespace model_Schwen_PONE2014 {

void sx0_fixedParameters_Schwen_PONE2014(realtype *sx0_fixedParameters, const realtype t, const realtype *x0, const realtype *p, const realtype *k, const int ip, gsl::span<const int> reinitialization_state_idxs){
    static const std::array<int, 1> _x0_fixedParameters_idxs = {
        0
    };
    for(auto idx: reinitialization_state_idxs) {
        if(std::find(_x0_fixedParameters_idxs.cbegin(), _x0_fixedParameters_idxs.cend(), idx) != _x0_fixedParameters_idxs.cend())
            sx0_fixedParameters[idx] = 0.0;
    }
}

} // namespace model_Schwen_PONE2014
} // namespace amici
