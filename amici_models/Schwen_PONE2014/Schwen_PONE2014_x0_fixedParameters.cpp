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

void x0_fixedParameters_Schwen_PONE2014(realtype *x0_fixedParameters, const realtype t, const realtype *p, const realtype *k, gsl::span<const int> reinitialization_state_idxs){
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 0) != reinitialization_state_idxs.cend())
        x0_fixedParameters[0] = init_Ins;
}

} // namespace model_Schwen_PONE2014
} // namespace amici
