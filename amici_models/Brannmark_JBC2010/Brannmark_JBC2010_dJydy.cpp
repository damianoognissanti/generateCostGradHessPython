#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>
#include <algorithm>

#include "Brannmark_JBC2010_p.h"
#include "Brannmark_JBC2010_k.h"
#include "Brannmark_JBC2010_y.h"
#include "Brannmark_JBC2010_sigmay.h"
#include "Brannmark_JBC2010_my.h"
#include "Brannmark_JBC2010_dJydy.h"

namespace amici {
namespace model_Brannmark_JBC2010 {

void dJydy_Brannmark_JBC2010(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (1.0*IR1_P - 1.0*mIR1_P)/std::pow(sigma_IR1_P, 2);
            break;
        case 1:
            dJydy[0] = (1.0*IRS1_P - 1.0*mIRS1_P)/std::pow(sigma_IRS1_P, 2);
            break;
        case 2:
            dJydy[0] = (1.0*IRS1_P_DosR - 1.0*mIRS1_P_DosR)/std::pow(sigma_IRS1_P_DosR, 2);
            break;
    }
}

} // namespace model_Brannmark_JBC2010
} // namespace amici
