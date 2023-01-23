#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "Bachmann_MSB2011_x_rdata.h"

namespace amici {
namespace model_Bachmann_MSB2011 {

void x_solver_Bachmann_MSB2011(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = EpoRJAK2;
    x_solver[1] = EpoRpJAK2;
    x_solver[2] = p1EpoRpJAK2;
    x_solver[3] = p2EpoRpJAK2;
    x_solver[4] = p12EpoRpJAK2;
    x_solver[5] = EpoRJAK2_CIS;
    x_solver[6] = SHP1;
    x_solver[7] = SHP1Act;
    x_solver[8] = STAT5;
    x_solver[9] = pSTAT5;
    x_solver[10] = npSTAT5;
    x_solver[11] = CISnRNA1;
    x_solver[12] = CISnRNA2;
    x_solver[13] = CISnRNA3;
    x_solver[14] = CISnRNA4;
    x_solver[15] = CISnRNA5;
    x_solver[16] = CISRNA;
    x_solver[17] = CIS;
    x_solver[18] = SOCS3nRNA1;
    x_solver[19] = SOCS3nRNA2;
    x_solver[20] = SOCS3nRNA3;
    x_solver[21] = SOCS3nRNA4;
    x_solver[22] = SOCS3nRNA5;
    x_solver[23] = SOCS3RNA;
    x_solver[24] = SOCS3;
}

} // namespace model_Bachmann_MSB2011
} // namespace amici
