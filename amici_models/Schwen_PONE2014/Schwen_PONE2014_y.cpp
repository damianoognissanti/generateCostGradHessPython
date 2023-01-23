#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>
#include <algorithm>

#include "Schwen_PONE2014_x.h"
#include "Schwen_PONE2014_p.h"
#include "Schwen_PONE2014_k.h"
#include "Schwen_PONE2014_w.h"

namespace amici {
namespace model_Schwen_PONE2014 {

void y_Schwen_PONE2014(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = observableParameter2_observable_IR1*(IR1 + IR1in + observableParameter1_observable_IR1);
    y[1] = observableParameter2_observable_IR2*(IR2 + IR2in + observableParameter1_observable_IR2);
    y[2] = observableParameter2_observable_IRsum*(0.60499999999999998*IR1 + 0.60499999999999998*IR1in + 0.39500000000000002*IR2 + 0.39500000000000002*IR2in + observableParameter1_observable_IRsum);
    y[3] = observableParameter1_observable_Insulin + observableParameter2_observable_Insulin*(Ins + InsulinFragments*observableParameter3_observable_Insulin)/(1 + (Ins + InsulinFragments*observableParameter3_observable_Insulin)/observableParameter4_observable_Insulin);
}

} // namespace model_Schwen_PONE2014
} // namespace amici
