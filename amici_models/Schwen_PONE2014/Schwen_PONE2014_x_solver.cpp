#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>
#include <algorithm>

#include "Schwen_PONE2014_x_rdata.h"

namespace amici {
namespace model_Schwen_PONE2014 {

void x_solver_Schwen_PONE2014(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = Ins;
    x_solver[1] = Rec1;
    x_solver[2] = Rec2;
    x_solver[3] = IR1;
    x_solver[4] = IR2;
    x_solver[5] = IR1in;
    x_solver[6] = IR2in;
    x_solver[7] = Uptake1;
    x_solver[8] = Uptake2;
    x_solver[9] = InsulinFragments;
    x_solver[10] = BoundUnspec;
}

} // namespace model_Schwen_PONE2014
} // namespace amici
