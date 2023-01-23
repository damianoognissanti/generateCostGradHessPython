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

void dydp_Schwen_PONE2014(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *tcl, const realtype *dtcldp){
    switch(ip) {
        case 13:
            dydp[0] = observableParameter2_observable_IR1;
            break;
        case 14:
            dydp[0] = IR1 + IR1in + observableParameter1_observable_IR1;
            break;
        case 15:
            dydp[1] = observableParameter2_observable_IR2;
            break;
        case 16:
            dydp[1] = IR2 + IR2in + observableParameter1_observable_IR2;
            break;
        case 17:
            dydp[2] = observableParameter2_observable_IRsum;
            break;
        case 18:
            dydp[2] = 0.60499999999999998*IR1 + 0.60499999999999998*IR1in + 0.39500000000000002*IR2 + 0.39500000000000002*IR2in + observableParameter1_observable_IRsum;
            break;
        case 19:
            dydp[3] = 1;
            break;
        case 20:
            dydp[3] = (Ins + InsulinFragments*observableParameter3_observable_Insulin)/(1 + (Ins + InsulinFragments*observableParameter3_observable_Insulin)/observableParameter4_observable_Insulin);
            break;
        case 21:
            dydp[3] = InsulinFragments*observableParameter2_observable_Insulin/(1 + (Ins + InsulinFragments*observableParameter3_observable_Insulin)/observableParameter4_observable_Insulin) - InsulinFragments*observableParameter2_observable_Insulin*(Ins + InsulinFragments*observableParameter3_observable_Insulin)/(observableParameter4_observable_Insulin*std::pow(1 + (Ins + InsulinFragments*observableParameter3_observable_Insulin)/observableParameter4_observable_Insulin, 2));
            break;
        case 22:
            dydp[3] = observableParameter2_observable_Insulin*std::pow(Ins + InsulinFragments*observableParameter3_observable_Insulin, 2)/(std::pow(observableParameter4_observable_Insulin, 2)*std::pow(1 + (Ins + InsulinFragments*observableParameter3_observable_Insulin)/observableParameter4_observable_Insulin, 2));
            break;
    }
}

} // namespace model_Schwen_PONE2014
} // namespace amici
