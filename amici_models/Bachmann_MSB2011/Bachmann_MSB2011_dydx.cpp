#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "Bachmann_MSB2011_x.h"
#include "Bachmann_MSB2011_p.h"
#include "Bachmann_MSB2011_k.h"
#include "Bachmann_MSB2011_w.h"
#include "Bachmann_MSB2011_dwdx.h"

namespace amici {
namespace model_Bachmann_MSB2011 {

void dydx_Bachmann_MSB2011(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[35] = 2*observableParameter2_observable_pJAK2_au/init_EpoRJAK2;
    dydx[54] = 16*observableParameter2_observable_pEpoR_au/init_EpoRJAK2;
    dydx[55] = 2*observableParameter2_observable_pJAK2_au/init_EpoRJAK2;
    dydx[74] = 16*observableParameter2_observable_pEpoR_au/init_EpoRJAK2;
    dydx[75] = 2*observableParameter2_observable_pJAK2_au/init_EpoRJAK2;
    dydx[94] = 16*observableParameter2_observable_pEpoR_au/init_EpoRJAK2;
    dydx[95] = 2*observableParameter2_observable_pJAK2_au/init_EpoRJAK2;
    dydx[127] = 1;
    dydx[138] = observableParameter1_observable_tSHP1_au/init_SHP1;
    dydx[147] = 1;
    dydx[158] = observableParameter1_observable_tSHP1_au/init_SHP1;
    dydx[173] = 1;
    dydx[176] = -100*pSTAT5/std::pow(STAT5 + pSTAT5, 2);
    dydx[179] = observableParameter1_observable_tSTAT5_au/init_STAT5;
    dydx[196] = -100*pSTAT5/std::pow(STAT5 + pSTAT5, 2) + 100/(STAT5 + pSTAT5);
    dydx[197] = observableParameter2_observable_pSTAT5_au/init_STAT5;
    dydx[199] = observableParameter1_observable_tSTAT5_au/init_STAT5;
    dydx[320] = observableParameter1_observable_CISRNA_foldA/CISRNAEqc;
    dydx[321] = observableParameter1_observable_CISRNA_foldB/CISRNAEqc;
    dydx[322] = observableParameter1_observable_CISRNA_foldC/CISRNAEqc;
    dydx[343] = 1;
    dydx[344] = observableParameter2_observable_CIS_au/CISEqc;
    dydx[345] = observableParameter1_observable_CIS_au1/CISEqc;
    dydx[346] = observableParameter1_observable_CIS_au2/CISEqc;
    dydx[468] = observableParameter1_observable_SOCS3RNA_foldA/SOCS3RNAEqc;
    dydx[469] = observableParameter1_observable_SOCS3RNA_foldB/SOCS3RNAEqc;
    dydx[470] = observableParameter1_observable_SOCS3RNA_foldC/SOCS3RNAEqc;
    dydx[491] = 1;
    dydx[492] = observableParameter2_observable_SOCS3_au/SOCS3Eqc;
}

} // namespace model_Bachmann_MSB2011
} // namespace amici
