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
#include "Schwen_PONE2014_my.h"

namespace amici {
namespace model_Schwen_PONE2014 {

void dJydsigma_Schwen_PONE2014(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_observable_IR1 - 1.0*std::pow(-std::log(mobservable_IR1)/M_LN10 + std::log(observable_IR1)/M_LN10, 2)/std::pow(sigma_observable_IR1, 3);
            break;
        case 1:
            dJydsigma[1] = 1.0/sigma_observable_IR2 - 1.0*std::pow(-std::log(mobservable_IR2)/M_LN10 + std::log(observable_IR2)/M_LN10, 2)/std::pow(sigma_observable_IR2, 3);
            break;
        case 2:
            dJydsigma[2] = 1.0/sigma_observable_IRsum - 1.0*std::pow(-std::log(mobservable_IRsum)/M_LN10 + std::log(observable_IRsum)/M_LN10, 2)/std::pow(sigma_observable_IRsum, 3);
            break;
        case 3:
            dJydsigma[3] = 1.0/sigma_observable_Insulin - 1.0*std::pow(-std::log(mobservable_Insulin)/M_LN10 + std::log(observable_Insulin)/M_LN10, 2)/std::pow(sigma_observable_Insulin, 3);
            break;
    }
}

} // namespace model_Schwen_PONE2014
} // namespace amici
