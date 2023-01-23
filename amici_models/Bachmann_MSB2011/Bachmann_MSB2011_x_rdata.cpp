#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "Bachmann_MSB2011_x.h"
#include "Bachmann_MSB2011_p.h"
#include "Bachmann_MSB2011_k.h"

namespace amici {
namespace model_Bachmann_MSB2011 {

void x_rdata_Bachmann_MSB2011(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = EpoRJAK2;
    x_rdata[1] = EpoRpJAK2;
    x_rdata[2] = p1EpoRpJAK2;
    x_rdata[3] = p2EpoRpJAK2;
    x_rdata[4] = p12EpoRpJAK2;
    x_rdata[5] = EpoRJAK2_CIS;
    x_rdata[6] = SHP1;
    x_rdata[7] = SHP1Act;
    x_rdata[8] = STAT5;
    x_rdata[9] = pSTAT5;
    x_rdata[10] = npSTAT5;
    x_rdata[11] = CISnRNA1;
    x_rdata[12] = CISnRNA2;
    x_rdata[13] = CISnRNA3;
    x_rdata[14] = CISnRNA4;
    x_rdata[15] = CISnRNA5;
    x_rdata[16] = CISRNA;
    x_rdata[17] = CIS;
    x_rdata[18] = SOCS3nRNA1;
    x_rdata[19] = SOCS3nRNA2;
    x_rdata[20] = SOCS3nRNA3;
    x_rdata[21] = SOCS3nRNA4;
    x_rdata[22] = SOCS3nRNA5;
    x_rdata[23] = SOCS3RNA;
    x_rdata[24] = SOCS3;
}

} // namespace model_Bachmann_MSB2011
} // namespace amici
