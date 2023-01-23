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

namespace amici {
namespace model_Brannmark_JBC2010 {

void sigmay_Brannmark_JBC2010(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_IR1_P = noiseParameter1_IR1_P;  // sigmay[0]
    sigma_IRS1_P = noiseParameter1_IRS1_P;  // sigmay[1]
    sigma_IRS1_P_DosR = noiseParameter1_IRS1_P_DosR;  // sigmay[2]
}

} // namespace model_Brannmark_JBC2010
} // namespace amici
