#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>
#include <algorithm>

#include "Brannmark_JBC2010_x.h"
#include "Brannmark_JBC2010_p.h"
#include "Brannmark_JBC2010_k.h"
#include "Brannmark_JBC2010_h.h"
#include "Brannmark_JBC2010_w.h"

namespace amici {
namespace model_Brannmark_JBC2010 {

void y_Brannmark_JBC2010(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = observableParameter1_IR1_P*(IRiP + IRp);
    y[1] = IRSiP*observableParameter1_IRS1_P;
    y[2] = IRSiP*observableParameter1_IRS1_P_DosR;
}

} // namespace model_Brannmark_JBC2010
} // namespace amici
