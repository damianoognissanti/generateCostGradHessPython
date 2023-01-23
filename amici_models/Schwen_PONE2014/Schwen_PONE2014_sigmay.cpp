#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>
#include <algorithm>

#include "Schwen_PONE2014_p.h"
#include "Schwen_PONE2014_k.h"
#include "Schwen_PONE2014_y.h"
#include "Schwen_PONE2014_sigmay.h"

namespace amici {
namespace model_Schwen_PONE2014 {

void sigmay_Schwen_PONE2014(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_observable_IR1 = noiseParameter1_observable_IR1;  // sigmay[0]
    sigma_observable_IR2 = noiseParameter1_observable_IR2;  // sigmay[1]
    sigma_observable_IRsum = noiseParameter1_observable_IRsum;  // sigmay[2]
    sigma_observable_Insulin = noiseParameter1_observable_Insulin;  // sigmay[3]
}

} // namespace model_Schwen_PONE2014
} // namespace amici
