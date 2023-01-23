#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "Bachmann_MSB2011_x.h"
#include "Bachmann_MSB2011_p.h"
#include "Bachmann_MSB2011_k.h"
#include "Bachmann_MSB2011_w.h"

namespace amici {
namespace model_Bachmann_MSB2011 {

void dydp_Bachmann_MSB2011(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *tcl, const realtype *dtcldp){
    switch(ip) {
        case 0:
            dydp[4] = -CIS*observableParameter2_observable_CIS_au/std::pow(CISEqc, 2);
            dydp[5] = -CIS*observableParameter1_observable_CIS_au1/std::pow(CISEqc, 2);
            dydp[6] = -CIS*observableParameter1_observable_CIS_au2/std::pow(CISEqc, 2);
            break;
        case 15:
            dydp[12] = -SOCS3*observableParameter2_observable_SOCS3_au/std::pow(SOCS3Eqc, 2);
            break;
        case 24:
            dydp[14] = -observableParameter2_observable_pEpoR_au*(16*p12EpoRpJAK2 + 16*p1EpoRpJAK2 + 16*p2EpoRpJAK2)/std::pow(init_EpoRJAK2, 2);
            dydp[15] = -observableParameter2_observable_pJAK2_au*(2*EpoRpJAK2 + 2*p12EpoRpJAK2 + 2*p1EpoRpJAK2 + 2*p2EpoRpJAK2)/std::pow(init_EpoRJAK2, 2);
            break;
        case 25:
            dydp[18] = -observableParameter1_observable_tSHP1_au*(SHP1 + SHP1Act)/std::pow(init_SHP1, 2);
            break;
        case 26:
            dydp[17] = -observableParameter2_observable_pSTAT5_au*pSTAT5/std::pow(init_STAT5, 2);
            dydp[19] = -observableParameter1_observable_tSTAT5_au*(STAT5 + pSTAT5)/std::pow(init_STAT5, 2);
            break;
        case 27:
            dydp[0] = CISRNA/CISRNAEqc;
            break;
        case 28:
            dydp[1] = CISRNA/CISRNAEqc;
            break;
        case 29:
            dydp[2] = CISRNA/CISRNAEqc;
            break;
        case 30:
            dydp[4] = 1;
            break;
        case 31:
            dydp[4] = CIS/CISEqc;
            break;
        case 32:
            dydp[5] = CIS/CISEqc;
            break;
        case 33:
            dydp[6] = CIS/CISEqc;
            break;
        case 34:
            dydp[8] = SOCS3RNA/SOCS3RNAEqc;
            break;
        case 35:
            dydp[9] = SOCS3RNA/SOCS3RNAEqc;
            break;
        case 36:
            dydp[10] = SOCS3RNA/SOCS3RNAEqc;
            break;
        case 37:
            dydp[12] = 1;
            break;
        case 38:
            dydp[12] = SOCS3/SOCS3Eqc;
            break;
        case 39:
            dydp[14] = 1;
            break;
        case 40:
            dydp[14] = (16*p12EpoRpJAK2 + 16*p1EpoRpJAK2 + 16*p2EpoRpJAK2)/init_EpoRJAK2;
            break;
        case 41:
            dydp[15] = 1;
            break;
        case 42:
            dydp[15] = (2*EpoRpJAK2 + 2*p12EpoRpJAK2 + 2*p1EpoRpJAK2 + 2*p2EpoRpJAK2)/init_EpoRJAK2;
            break;
        case 43:
            dydp[16] = 1;
            break;
        case 44:
            dydp[17] = 1;
            break;
        case 45:
            dydp[17] = pSTAT5/init_STAT5;
            break;
        case 46:
            dydp[18] = (SHP1 + SHP1Act)/init_SHP1;
            break;
        case 47:
            dydp[19] = (STAT5 + pSTAT5)/init_STAT5;
            break;
    }
}

} // namespace model_Bachmann_MSB2011
} // namespace amici
