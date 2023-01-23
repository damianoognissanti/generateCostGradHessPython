#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>
#include <algorithm>

#include "Schwen_PONE2014_x.h"
#include "Schwen_PONE2014_p.h"
#include "Schwen_PONE2014_k.h"

namespace amici {
namespace model_Schwen_PONE2014 {

void x_rdata_Schwen_PONE2014(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = Ins;
    x_rdata[1] = Rec1;
    x_rdata[2] = Rec2;
    x_rdata[3] = IR1;
    x_rdata[4] = IR2;
    x_rdata[5] = IR1in;
    x_rdata[6] = IR2in;
    x_rdata[7] = Uptake1;
    x_rdata[8] = Uptake2;
    x_rdata[9] = InsulinFragments;
    x_rdata[10] = BoundUnspec;
}

} // namespace model_Schwen_PONE2014
} // namespace amici
