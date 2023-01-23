#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>
#include <algorithm>

#include "Brannmark_JBC2010_x.h"
#include "Brannmark_JBC2010_p.h"
#include "Brannmark_JBC2010_k.h"

namespace amici {
namespace model_Brannmark_JBC2010 {

void x_rdata_Brannmark_JBC2010(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = IR;
    x_rdata[1] = IRins;
    x_rdata[2] = IRp;
    x_rdata[3] = IRiP;
    x_rdata[4] = IRi;
    x_rdata[5] = IRS;
    x_rdata[6] = IRSiP;
    x_rdata[7] = X;
    x_rdata[8] = Xp;
}

} // namespace model_Brannmark_JBC2010
} // namespace amici
