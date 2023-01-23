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

void Jy_Schwen_PONE2014(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_IR1, 2)*std::pow(sigma_observable_IR1, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(-std::log(mobservable_IR1)/M_LN10 + std::log(observable_IR1)/M_LN10, 2)/std::pow(sigma_observable_IR1, 2);
            break;
        case 1:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_IR2, 2)*std::pow(sigma_observable_IR2, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(-std::log(mobservable_IR2)/M_LN10 + std::log(observable_IR2)/M_LN10, 2)/std::pow(sigma_observable_IR2, 2);
            break;
        case 2:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_IRsum, 2)*std::pow(sigma_observable_IRsum, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(-std::log(mobservable_IRsum)/M_LN10 + std::log(observable_IRsum)/M_LN10, 2)/std::pow(sigma_observable_IRsum, 2);
            break;
        case 3:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_Insulin, 2)*std::pow(sigma_observable_Insulin, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(-std::log(mobservable_Insulin)/M_LN10 + std::log(observable_Insulin)/M_LN10, 2)/std::pow(sigma_observable_Insulin, 2);
            break;
    }
}

} // namespace model_Schwen_PONE2014
} // namespace amici
