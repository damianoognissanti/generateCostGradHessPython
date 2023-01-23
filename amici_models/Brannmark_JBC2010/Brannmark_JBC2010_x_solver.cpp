#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>
#include <algorithm>

#include "Brannmark_JBC2010_x_rdata.h"

namespace amici {
namespace model_Brannmark_JBC2010 {

void x_solver_Brannmark_JBC2010(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = IR;
    x_solver[1] = IRins;
    x_solver[2] = IRp;
    x_solver[3] = IRiP;
    x_solver[4] = IRi;
    x_solver[5] = IRS;
    x_solver[6] = IRSiP;
    x_solver[7] = X;
    x_solver[8] = Xp;
}

} // namespace model_Brannmark_JBC2010
} // namespace amici
