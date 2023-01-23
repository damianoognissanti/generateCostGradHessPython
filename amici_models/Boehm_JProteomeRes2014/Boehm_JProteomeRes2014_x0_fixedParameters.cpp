#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "Boehm_JProteomeRes2014_p.h"
#include "Boehm_JProteomeRes2014_k.h"

namespace amici {
namespace model_Boehm_JProteomeRes2014 {

void x0_fixedParameters_Boehm_JProteomeRes2014(realtype *x0_fixedParameters, const realtype t, const realtype *p, const realtype *k, gsl::span<const int> reinitialization_state_idxs){
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 0) != reinitialization_state_idxs.cend())
        x0_fixedParameters[0] = 207.59999999999999*ratio;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 1) != reinitialization_state_idxs.cend())
        x0_fixedParameters[1] = 207.59999999999999 - 207.59999999999999*ratio;
}

} // namespace model_Boehm_JProteomeRes2014
} // namespace amici
