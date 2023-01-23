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

void y_Bachmann_MSB2011(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = CISRNA*observableParameter1_observable_CISRNA_foldA/CISRNAEqc + 1;
    y[1] = CISRNA*observableParameter1_observable_CISRNA_foldB/CISRNAEqc + 1;
    y[2] = CISRNA*observableParameter1_observable_CISRNA_foldC/CISRNAEqc + 1;
    y[3] = CIS;
    y[4] = CIS*observableParameter2_observable_CIS_au/CISEqc + observableParameter1_observable_CIS_au;
    y[5] = CIS*observableParameter1_observable_CIS_au1/CISEqc;
    y[6] = CIS*observableParameter1_observable_CIS_au2/CISEqc;
    y[7] = SHP1 + SHP1Act;
    y[8] = SOCS3RNA*observableParameter1_observable_SOCS3RNA_foldA/SOCS3RNAEqc + 1;
    y[9] = SOCS3RNA*observableParameter1_observable_SOCS3RNA_foldB/SOCS3RNAEqc + 1;
    y[10] = SOCS3RNA*observableParameter1_observable_SOCS3RNA_foldC/SOCS3RNAEqc + 1;
    y[11] = SOCS3;
    y[12] = SOCS3*observableParameter2_observable_SOCS3_au/SOCS3Eqc + observableParameter1_observable_SOCS3_au;
    y[13] = STAT5;
    y[14] = observableParameter1_observable_pEpoR_au + observableParameter2_observable_pEpoR_au*(16*p12EpoRpJAK2 + 16*p1EpoRpJAK2 + 16*p2EpoRpJAK2)/init_EpoRJAK2;
    y[15] = observableParameter1_observable_pJAK2_au + observableParameter2_observable_pJAK2_au*(2*EpoRpJAK2 + 2*p12EpoRpJAK2 + 2*p1EpoRpJAK2 + 2*p2EpoRpJAK2)/init_EpoRJAK2;
    y[16] = observableParameter1_observable_pSTAT5B_rel + 100*pSTAT5/(STAT5 + pSTAT5);
    y[17] = observableParameter1_observable_pSTAT5_au + observableParameter2_observable_pSTAT5_au*pSTAT5/init_STAT5;
    y[18] = observableParameter1_observable_tSHP1_au*(SHP1 + SHP1Act)/init_SHP1;
    y[19] = observableParameter1_observable_tSTAT5_au*(STAT5 + pSTAT5)/init_STAT5;
}

} // namespace model_Bachmann_MSB2011
} // namespace amici
