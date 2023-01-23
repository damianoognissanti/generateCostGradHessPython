#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>
#include <algorithm>

#include "Brannmark_JBC2010_p.h"
#include "Brannmark_JBC2010_k.h"

namespace amici {
namespace model_Brannmark_JBC2010 {

void x0_Brannmark_JBC2010(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 9.9495764278756909;
    x0[1] = 0.017397222172539299;
    x0[2] = 1.7629010620181e-9;
    x0[3] = 1.1159002615229599e-5;
    x0[4] = 0.033015189186268099;
    x0[5] = 9.8669934870136693;
    x0[6] = 0.13300651298633601;
    x0[7] = 9.9998419948735098;
    x0[8] = 0.000158005126497888;
}

} // namespace model_Brannmark_JBC2010
} // namespace amici
