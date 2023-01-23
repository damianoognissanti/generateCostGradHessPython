#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "Bachmann_MSB2011_p.h"
#include "Bachmann_MSB2011_k.h"

namespace amici {
namespace model_Bachmann_MSB2011 {

void sx0_fixedParameters_Bachmann_MSB2011(realtype *sx0_fixedParameters, const realtype t, const realtype *x0, const realtype *p, const realtype *k, const int ip, gsl::span<const int> reinitialization_state_idxs){
    static const std::array<int, 4> _x0_fixedParameters_idxs = {
        5, 6, 17, 24
    };
    for(auto idx: reinitialization_state_idxs) {
        if(std::find(_x0_fixedParameters_idxs.cbegin(), _x0_fixedParameters_idxs.cend(), idx) != _x0_fixedParameters_idxs.cend())
            sx0_fixedParameters[idx] = 0.0;
    }
    switch(ip) {
        case 0:
            if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 17) != reinitialization_state_idxs.cend())
                sx0_fixedParameters[17] = CISEqcOE*init_CIS_multiplier;
            break;
        case 5:
            if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 17) != reinitialization_state_idxs.cend())
                sx0_fixedParameters[17] = CISEqc*init_CIS_multiplier;
            break;
        case 13:
            if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 6) != reinitialization_state_idxs.cend())
                sx0_fixedParameters[6] = init_SHP1*init_SHP1_multiplier;
            break;
        case 14:
            if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 24) != reinitialization_state_idxs.cend())
                sx0_fixedParameters[24] = SOCS3Eqc*init_SOCS3_multiplier;
            break;
        case 15:
            if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 24) != reinitialization_state_idxs.cend())
                sx0_fixedParameters[24] = SOCS3EqcOE*init_SOCS3_multiplier;
            break;
        case 25:
            if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 6) != reinitialization_state_idxs.cend())
                sx0_fixedParameters[6] = SHP1ProOE*init_SHP1_multiplier + 1;
            break;
    }
}

} // namespace model_Bachmann_MSB2011
} // namespace amici
