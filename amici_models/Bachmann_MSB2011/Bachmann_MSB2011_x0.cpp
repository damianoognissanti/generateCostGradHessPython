#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "Bachmann_MSB2011_p.h"
#include "Bachmann_MSB2011_k.h"

namespace amici {
namespace model_Bachmann_MSB2011 {

void x0_Bachmann_MSB2011(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = init_EpoRJAK2;
    x0[5] = init_EpoRJAK2_CIS;
    x0[6] = init_SHP1*(SHP1ProOE*init_SHP1_multiplier + 1);
    x0[8] = init_STAT5;
    x0[17] = CISEqc*CISEqcOE*init_CIS_multiplier;
    x0[24] = SOCS3Eqc*SOCS3EqcOE*init_SOCS3_multiplier;
}

} // namespace model_Bachmann_MSB2011
} // namespace amici
