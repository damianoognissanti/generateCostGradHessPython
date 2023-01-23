#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "Boehm_JProteomeRes2014_p.h"
#include "Boehm_JProteomeRes2014_k.h"

namespace amici {
namespace model_Boehm_JProteomeRes2014 {

void sx0_fixedParameters_Boehm_JProteomeRes2014(realtype *sx0_fixedParameters, const realtype t, const realtype *x0, const realtype *p, const realtype *k, const int ip, gsl::span<const int> reinitialization_state_idxs){
    static const std::array<int, 2> _x0_fixedParameters_idxs = {
        0, 1
    };
    for(auto idx: reinitialization_state_idxs) {
        if(std::find(_x0_fixedParameters_idxs.cbegin(), _x0_fixedParameters_idxs.cend(), idx) != _x0_fixedParameters_idxs.cend())
            sx0_fixedParameters[idx] = 0.0;
    }
}

} // namespace model_Boehm_JProteomeRes2014
} // namespace amici
