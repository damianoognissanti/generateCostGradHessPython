#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "Bachmann_MSB2011_p.h"
#include "Bachmann_MSB2011_k.h"

namespace amici {
namespace model_Bachmann_MSB2011 {

void x0_fixedParameters_Bachmann_MSB2011(realtype *x0_fixedParameters, const realtype t, const realtype *p, const realtype *k, gsl::span<const int> reinitialization_state_idxs){
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 5) != reinitialization_state_idxs.cend())
        x0_fixedParameters[5] = init_EpoRJAK2_CIS;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 6) != reinitialization_state_idxs.cend())
        x0_fixedParameters[6] = init_SHP1*(SHP1ProOE*init_SHP1_multiplier + 1);
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 17) != reinitialization_state_idxs.cend())
        x0_fixedParameters[17] = CISEqc*CISEqcOE*init_CIS_multiplier;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 24) != reinitialization_state_idxs.cend())
        x0_fixedParameters[24] = SOCS3Eqc*SOCS3EqcOE*init_SOCS3_multiplier;
}

} // namespace model_Bachmann_MSB2011
} // namespace amici
