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

void w_Bachmann_MSB2011(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    flux_v1_v_0 = 0.40000000000000002*Epo*EpoRJAK2*JAK2ActEpo/(SOCS3*SOCS3Inh/SOCS3Eqc + 1);  // w[0]
    flux_v2_v_1 = 0.40000000000000002*EpoRpJAK2*JAK2EpoRDeaSHP1*SHP1Act/init_SHP1;  // w[1]
    flux_v3_v_2 = 0.40000000000000002*EpoRActJAK2*EpoRpJAK2/(SOCS3*SOCS3Inh/SOCS3Eqc + 1);  // w[2]
    flux_v4_v_3 = 1.2000000000000002*EpoRActJAK2*EpoRpJAK2/((EpoRCISInh*EpoRJAK2_CIS + 1)*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));  // w[3]
    flux_v5_v_4 = 1.2000000000000002*EpoRActJAK2*p1EpoRpJAK2/((EpoRCISInh*EpoRJAK2_CIS + 1)*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));  // w[4]
    flux_v6_v_5 = 0.40000000000000002*EpoRActJAK2*p2EpoRpJAK2/(SOCS3*SOCS3Inh/SOCS3Eqc + 1);  // w[5]
    flux_v7_v_6 = 0.40000000000000002*JAK2EpoRDeaSHP1*SHP1Act*p1EpoRpJAK2/init_SHP1;  // w[6]
    flux_v8_v_7 = 0.40000000000000002*JAK2EpoRDeaSHP1*SHP1Act*p2EpoRpJAK2/init_SHP1;  // w[7]
    flux_v9_v_8 = 0.40000000000000002*JAK2EpoRDeaSHP1*SHP1Act*p12EpoRpJAK2/init_SHP1;  // w[8]
    flux_v10_v_9 = 0.40000000000000002*EpoRCISRemove*EpoRJAK2_CIS*(p12EpoRpJAK2 + p1EpoRpJAK2)/init_EpoRJAK2;  // w[9]
    flux_v11_v_10 = 0.40000000000000002*SHP1*SHP1ActEpoR*(EpoRpJAK2 + p12EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2)/init_EpoRJAK2;  // w[10]
    flux_v12_v_11 = 0.40000000000000002*SHP1Act*SHP1Dea;  // w[11]
    flux_v13_v_12 = 0.40000000000000002*STAT5*STAT5ActJAK2*(EpoRpJAK2 + p12EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2)/(init_EpoRJAK2*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));  // w[12]
    flux_v14_v_13 = 0.40000000000000002*STAT5*STAT5ActEpoR*std::pow(p12EpoRpJAK2 + p1EpoRpJAK2, 2)/(std::pow(init_EpoRJAK2, 2)*(CIS*CISInh/CISEqc + 1)*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));  // w[13]
    flux_v15_v_14 = 0.40000000000000002*STAT5Imp*pSTAT5;  // w[14]
    flux_v16_v_15 = 0.27500000000000002*STAT5Exp*npSTAT5;  // w[15]
    flux_v17_v_16 = 0.27500000000000002*ActD*CISRNAEqc*CISRNATurn*npSTAT5/init_STAT5;  // w[16]
    flux_v18_v_17 = 0.27500000000000002*CISRNADelay*CISnRNA1;  // w[17]
    flux_v19_v_18 = 0.27500000000000002*CISRNADelay*CISnRNA2;  // w[18]
    flux_v20_v_19 = 0.27500000000000002*CISRNADelay*CISnRNA3;  // w[19]
    flux_v21_v_20 = 0.27500000000000002*CISRNADelay*CISnRNA4;  // w[20]
    flux_v22_v_21 = 0.27500000000000002*CISRNADelay*CISnRNA5;  // w[21]
    flux_v23_v_22 = 0.40000000000000002*CISRNA*CISRNATurn;  // w[22]
    flux_v24_v_23 = 0.40000000000000002*CISEqc*CISRNA*CISTurn/CISRNAEqc;  // w[23]
    flux_v25_v_24 = 0.40000000000000002*CIS*CISTurn;  // w[24]
    flux_v26_v_25 = 0.40000000000000002*CISEqc*CISEqcOE*CISTurn*CISoe;  // w[25]
    flux_v27_v_26 = 0.27500000000000002*ActD*SOCS3RNAEqc*SOCS3RNATurn*npSTAT5/init_STAT5;  // w[26]
    flux_v28_v_27 = 0.27500000000000002*SOCS3RNADelay*SOCS3nRNA1;  // w[27]
    flux_v29_v_28 = 0.27500000000000002*SOCS3RNADelay*SOCS3nRNA2;  // w[28]
    flux_v30_v_29 = 0.27500000000000002*SOCS3RNADelay*SOCS3nRNA3;  // w[29]
    flux_v31_v_30 = 0.27500000000000002*SOCS3RNADelay*SOCS3nRNA4;  // w[30]
    flux_v32_v_31 = 0.27500000000000002*SOCS3RNADelay*SOCS3nRNA5;  // w[31]
    flux_v33_v_32 = 0.40000000000000002*SOCS3RNA*SOCS3RNATurn;  // w[32]
    flux_v34_v_33 = 0.40000000000000002*SOCS3Eqc*SOCS3RNA*SOCS3Turn/SOCS3RNAEqc;  // w[33]
    flux_v35_v_34 = 0.40000000000000002*SOCS3*SOCS3Turn;  // w[34]
    flux_v36_v_35 = 0.40000000000000002*SOCS3Eqc*SOCS3EqcOE*SOCS3Turn*SOCS3oe;  // w[35]
}

} // namespace model_Bachmann_MSB2011
} // namespace amici
