#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "Bachmann_MSB2011_p.h"
#include "Bachmann_MSB2011_k.h"
#include "Bachmann_MSB2011_y.h"
#include "Bachmann_MSB2011_sigmay.h"

namespace amici {
namespace model_Bachmann_MSB2011 {

void sigmay_Bachmann_MSB2011(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_observable_CISRNA_foldA = noiseParameter1_observable_CISRNA_foldA;  // sigmay[0]
    sigma_observable_CISRNA_foldB = noiseParameter1_observable_CISRNA_foldB;  // sigmay[1]
    sigma_observable_CISRNA_foldC = noiseParameter1_observable_CISRNA_foldC;  // sigmay[2]
    sigma_observable_CIS_abs = noiseParameter1_observable_CIS_abs;  // sigmay[3]
    sigma_observable_CIS_au = noiseParameter1_observable_CIS_au;  // sigmay[4]
    sigma_observable_CIS_au1 = noiseParameter1_observable_CIS_au1;  // sigmay[5]
    sigma_observable_CIS_au2 = noiseParameter1_observable_CIS_au2;  // sigmay[6]
    sigma_observable_SHP1_abs = noiseParameter1_observable_SHP1_abs;  // sigmay[7]
    sigma_observable_SOCS3RNA_foldA = noiseParameter1_observable_SOCS3RNA_foldA;  // sigmay[8]
    sigma_observable_SOCS3RNA_foldB = noiseParameter1_observable_SOCS3RNA_foldB;  // sigmay[9]
    sigma_observable_SOCS3RNA_foldC = noiseParameter1_observable_SOCS3RNA_foldC;  // sigmay[10]
    sigma_observable_SOCS3_abs = noiseParameter1_observable_SOCS3_abs;  // sigmay[11]
    sigma_observable_SOCS3_au = noiseParameter1_observable_SOCS3_au;  // sigmay[12]
    sigma_observable_STAT5_abs = noiseParameter1_observable_STAT5_abs;  // sigmay[13]
    sigma_observable_pEpoR_au = noiseParameter1_observable_pEpoR_au;  // sigmay[14]
    sigma_observable_pJAK2_au = noiseParameter1_observable_pJAK2_au;  // sigmay[15]
    sigma_observable_pSTAT5B_rel = noiseParameter1_observable_pSTAT5B_rel;  // sigmay[16]
    sigma_observable_pSTAT5_au = noiseParameter1_observable_pSTAT5_au + noiseParameter2_observable_pSTAT5_au;  // sigmay[17]
    sigma_observable_tSHP1_au = noiseParameter1_observable_tSHP1_au;  // sigmay[18]
    sigma_observable_tSTAT5_au = noiseParameter1_observable_tSTAT5_au;  // sigmay[19]
}

} // namespace model_Bachmann_MSB2011
} // namespace amici
