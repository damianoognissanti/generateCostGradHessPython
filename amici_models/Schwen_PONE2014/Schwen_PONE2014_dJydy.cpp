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
#include "Schwen_PONE2014_dJydy.h"

namespace amici {
namespace model_Schwen_PONE2014 {

void dJydy_Schwen_PONE2014(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*std::log(mobservable_IR1)/M_LN10 + 1.0*std::log(observable_IR1)/M_LN10)/(observable_IR1*std::pow(sigma_observable_IR1, 2)*M_LN10);
            break;
        case 1:
            dJydy[0] = (-1.0*std::log(mobservable_IR2)/M_LN10 + 1.0*std::log(observable_IR2)/M_LN10)/(observable_IR2*std::pow(sigma_observable_IR2, 2)*M_LN10);
            break;
        case 2:
            dJydy[0] = (-1.0*std::log(mobservable_IRsum)/M_LN10 + 1.0*std::log(observable_IRsum)/M_LN10)/(observable_IRsum*std::pow(sigma_observable_IRsum, 2)*M_LN10);
            break;
        case 3:
            dJydy[0] = (-1.0*std::log(mobservable_Insulin)/M_LN10 + 1.0*std::log(observable_Insulin)/M_LN10)/(observable_Insulin*std::pow(sigma_observable_Insulin, 2)*M_LN10);
            break;
    }
}

} // namespace model_Schwen_PONE2014
} // namespace amici
