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
#include "Schwen_PONE2014_dwdx.h"

namespace amici {
namespace model_Schwen_PONE2014 {

void dydx_Schwen_PONE2014(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[3] = observableParameter2_observable_Insulin/(1 + (Ins + InsulinFragments*observableParameter3_observable_Insulin)/observableParameter4_observable_Insulin) - observableParameter2_observable_Insulin*(Ins + InsulinFragments*observableParameter3_observable_Insulin)/(observableParameter4_observable_Insulin*std::pow(1 + (Ins + InsulinFragments*observableParameter3_observable_Insulin)/observableParameter4_observable_Insulin, 2));
    dydx[12] = observableParameter2_observable_IR1;
    dydx[14] = 0.60499999999999998*observableParameter2_observable_IRsum;
    dydx[17] = observableParameter2_observable_IR2;
    dydx[18] = 0.39500000000000002*observableParameter2_observable_IRsum;
    dydx[20] = observableParameter2_observable_IR1;
    dydx[22] = 0.60499999999999998*observableParameter2_observable_IRsum;
    dydx[25] = observableParameter2_observable_IR2;
    dydx[26] = 0.39500000000000002*observableParameter2_observable_IRsum;
    dydx[39] = observableParameter2_observable_Insulin*observableParameter3_observable_Insulin/(1 + (Ins + InsulinFragments*observableParameter3_observable_Insulin)/observableParameter4_observable_Insulin) - observableParameter2_observable_Insulin*observableParameter3_observable_Insulin*(Ins + InsulinFragments*observableParameter3_observable_Insulin)/(observableParameter4_observable_Insulin*std::pow(1 + (Ins + InsulinFragments*observableParameter3_observable_Insulin)/observableParameter4_observable_Insulin, 2));
}

} // namespace model_Schwen_PONE2014
} // namespace amici
