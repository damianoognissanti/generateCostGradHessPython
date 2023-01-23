#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "Bachmann_MSB2011_x.h"
#include "Bachmann_MSB2011_p.h"
#include "Bachmann_MSB2011_k.h"
#include "Bachmann_MSB2011_w.h"
#include "Bachmann_MSB2011_dwdp.h"

namespace amici {
namespace model_Bachmann_MSB2011 {

void dwdp_Bachmann_MSB2011(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp){
    dflux_v14_v_13_dCISEqc = 0.40000000000000002*CIS*CISInh*STAT5*STAT5ActEpoR*std::pow(p12EpoRpJAK2 + p1EpoRpJAK2, 2)/(std::pow(CISEqc, 2)*std::pow(init_EpoRJAK2, 2)*std::pow(CIS*CISInh/CISEqc + 1, 2)*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));  // dwdp[0]
    dflux_v24_v_23_dCISEqc = 0.40000000000000002*CISRNA*CISTurn/CISRNAEqc;  // dwdp[1]
    dflux_v26_v_25_dCISEqc = 0.40000000000000002*CISEqcOE*CISTurn*CISoe;  // dwdp[2]
    dflux_v14_v_13_dCISInh = -0.40000000000000002*CIS*STAT5*STAT5ActEpoR*std::pow(p12EpoRpJAK2 + p1EpoRpJAK2, 2)/(CISEqc*std::pow(init_EpoRJAK2, 2)*std::pow(CIS*CISInh/CISEqc + 1, 2)*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));  // dwdp[3]
    dflux_v18_v_17_dCISRNADelay = 0.27500000000000002*CISnRNA1;  // dwdp[4]
    dflux_v19_v_18_dCISRNADelay = 0.27500000000000002*CISnRNA2;  // dwdp[5]
    dflux_v20_v_19_dCISRNADelay = 0.27500000000000002*CISnRNA3;  // dwdp[6]
    dflux_v21_v_20_dCISRNADelay = 0.27500000000000002*CISnRNA4;  // dwdp[7]
    dflux_v22_v_21_dCISRNADelay = 0.27500000000000002*CISnRNA5;  // dwdp[8]
    dflux_v17_v_16_dCISRNATurn = 0.27500000000000002*ActD*CISRNAEqc*npSTAT5/init_STAT5;  // dwdp[9]
    dflux_v23_v_22_dCISRNATurn = 0.40000000000000002*CISRNA;  // dwdp[10]
    dflux_v24_v_23_dCISTurn = 0.40000000000000002*CISEqc*CISRNA/CISRNAEqc;  // dwdp[11]
    dflux_v25_v_24_dCISTurn = 0.40000000000000002*CIS;  // dwdp[12]
    dflux_v26_v_25_dCISTurn = 0.40000000000000002*CISEqc*CISEqcOE*CISoe;  // dwdp[13]
    dflux_v26_v_25_dCISEqcOE = 0.40000000000000002*CISEqc*CISTurn*CISoe;  // dwdp[14]
    dflux_v3_v_2_dEpoRActJAK2 = 0.40000000000000002*EpoRpJAK2/(SOCS3*SOCS3Inh/SOCS3Eqc + 1);  // dwdp[15]
    dflux_v4_v_3_dEpoRActJAK2 = 1.2000000000000002*EpoRpJAK2/((EpoRCISInh*EpoRJAK2_CIS + 1)*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));  // dwdp[16]
    dflux_v5_v_4_dEpoRActJAK2 = 1.2000000000000002*p1EpoRpJAK2/((EpoRCISInh*EpoRJAK2_CIS + 1)*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));  // dwdp[17]
    dflux_v6_v_5_dEpoRActJAK2 = 0.40000000000000002*p2EpoRpJAK2/(SOCS3*SOCS3Inh/SOCS3Eqc + 1);  // dwdp[18]
    dflux_v4_v_3_dEpoRCISInh = -1.2000000000000002*EpoRActJAK2*EpoRJAK2_CIS*EpoRpJAK2/(std::pow(EpoRCISInh*EpoRJAK2_CIS + 1, 2)*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));  // dwdp[19]
    dflux_v5_v_4_dEpoRCISInh = -1.2000000000000002*EpoRActJAK2*EpoRJAK2_CIS*p1EpoRpJAK2/(std::pow(EpoRCISInh*EpoRJAK2_CIS + 1, 2)*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));  // dwdp[20]
    dflux_v10_v_9_dEpoRCISRemove = 0.40000000000000002*EpoRJAK2_CIS*(p12EpoRpJAK2 + p1EpoRpJAK2)/init_EpoRJAK2;  // dwdp[21]
    dflux_v1_v_0_dJAK2ActEpo = 0.40000000000000002*Epo*EpoRJAK2/(SOCS3*SOCS3Inh/SOCS3Eqc + 1);  // dwdp[22]
    dflux_v2_v_1_dJAK2EpoRDeaSHP1 = 0.40000000000000002*EpoRpJAK2*SHP1Act/init_SHP1;  // dwdp[23]
    dflux_v7_v_6_dJAK2EpoRDeaSHP1 = 0.40000000000000002*SHP1Act*p1EpoRpJAK2/init_SHP1;  // dwdp[24]
    dflux_v8_v_7_dJAK2EpoRDeaSHP1 = 0.40000000000000002*SHP1Act*p2EpoRpJAK2/init_SHP1;  // dwdp[25]
    dflux_v9_v_8_dJAK2EpoRDeaSHP1 = 0.40000000000000002*SHP1Act*p12EpoRpJAK2/init_SHP1;  // dwdp[26]
    dflux_v11_v_10_dSHP1ActEpoR = 0.40000000000000002*SHP1*(EpoRpJAK2 + p12EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2)/init_EpoRJAK2;  // dwdp[27]
    dflux_v12_v_11_dSHP1Dea = 0.40000000000000002*SHP1Act;  // dwdp[28]
    dflux_v36_v_35_dSOCS3EqcOE = 0.40000000000000002*SOCS3Eqc*SOCS3Turn*SOCS3oe;  // dwdp[29]
    dflux_v1_v_0_dSOCS3Eqc = 0.40000000000000002*Epo*EpoRJAK2*JAK2ActEpo*SOCS3*SOCS3Inh/(std::pow(SOCS3Eqc, 2)*std::pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));  // dwdp[30]
    dflux_v3_v_2_dSOCS3Eqc = 0.40000000000000002*EpoRActJAK2*EpoRpJAK2*SOCS3*SOCS3Inh/(std::pow(SOCS3Eqc, 2)*std::pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));  // dwdp[31]
    dflux_v4_v_3_dSOCS3Eqc = 1.2000000000000002*EpoRActJAK2*EpoRpJAK2*SOCS3*SOCS3Inh/(std::pow(SOCS3Eqc, 2)*(EpoRCISInh*EpoRJAK2_CIS + 1)*std::pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));  // dwdp[32]
    dflux_v5_v_4_dSOCS3Eqc = 1.2000000000000002*EpoRActJAK2*SOCS3*SOCS3Inh*p1EpoRpJAK2/(std::pow(SOCS3Eqc, 2)*(EpoRCISInh*EpoRJAK2_CIS + 1)*std::pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));  // dwdp[33]
    dflux_v6_v_5_dSOCS3Eqc = 0.40000000000000002*EpoRActJAK2*SOCS3*SOCS3Inh*p2EpoRpJAK2/(std::pow(SOCS3Eqc, 2)*std::pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));  // dwdp[34]
    dflux_v13_v_12_dSOCS3Eqc = 0.40000000000000002*SOCS3*SOCS3Inh*STAT5*STAT5ActJAK2*(EpoRpJAK2 + p12EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2)/(std::pow(SOCS3Eqc, 2)*init_EpoRJAK2*std::pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));  // dwdp[35]
    dflux_v14_v_13_dSOCS3Eqc = 0.40000000000000002*SOCS3*SOCS3Inh*STAT5*STAT5ActEpoR*std::pow(p12EpoRpJAK2 + p1EpoRpJAK2, 2)/(std::pow(SOCS3Eqc, 2)*std::pow(init_EpoRJAK2, 2)*(CIS*CISInh/CISEqc + 1)*std::pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));  // dwdp[36]
    dflux_v34_v_33_dSOCS3Eqc = 0.40000000000000002*SOCS3RNA*SOCS3Turn/SOCS3RNAEqc;  // dwdp[37]
    dflux_v36_v_35_dSOCS3Eqc = 0.40000000000000002*SOCS3EqcOE*SOCS3Turn*SOCS3oe;  // dwdp[38]
    dflux_v1_v_0_dSOCS3Inh = -0.40000000000000002*Epo*EpoRJAK2*JAK2ActEpo*SOCS3/(SOCS3Eqc*std::pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));  // dwdp[39]
    dflux_v3_v_2_dSOCS3Inh = -0.40000000000000002*EpoRActJAK2*EpoRpJAK2*SOCS3/(SOCS3Eqc*std::pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));  // dwdp[40]
    dflux_v4_v_3_dSOCS3Inh = -1.2000000000000002*EpoRActJAK2*EpoRpJAK2*SOCS3/(SOCS3Eqc*(EpoRCISInh*EpoRJAK2_CIS + 1)*std::pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));  // dwdp[41]
    dflux_v5_v_4_dSOCS3Inh = -1.2000000000000002*EpoRActJAK2*SOCS3*p1EpoRpJAK2/(SOCS3Eqc*(EpoRCISInh*EpoRJAK2_CIS + 1)*std::pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));  // dwdp[42]
    dflux_v6_v_5_dSOCS3Inh = -0.40000000000000002*EpoRActJAK2*SOCS3*p2EpoRpJAK2/(SOCS3Eqc*std::pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));  // dwdp[43]
    dflux_v13_v_12_dSOCS3Inh = -0.40000000000000002*SOCS3*STAT5*STAT5ActJAK2*(EpoRpJAK2 + p12EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2)/(SOCS3Eqc*init_EpoRJAK2*std::pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));  // dwdp[44]
    dflux_v14_v_13_dSOCS3Inh = -0.40000000000000002*SOCS3*STAT5*STAT5ActEpoR*std::pow(p12EpoRpJAK2 + p1EpoRpJAK2, 2)/(SOCS3Eqc*std::pow(init_EpoRJAK2, 2)*(CIS*CISInh/CISEqc + 1)*std::pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));  // dwdp[45]
    dflux_v28_v_27_dSOCS3RNADelay = 0.27500000000000002*SOCS3nRNA1;  // dwdp[46]
    dflux_v29_v_28_dSOCS3RNADelay = 0.27500000000000002*SOCS3nRNA2;  // dwdp[47]
    dflux_v30_v_29_dSOCS3RNADelay = 0.27500000000000002*SOCS3nRNA3;  // dwdp[48]
    dflux_v31_v_30_dSOCS3RNADelay = 0.27500000000000002*SOCS3nRNA4;  // dwdp[49]
    dflux_v32_v_31_dSOCS3RNADelay = 0.27500000000000002*SOCS3nRNA5;  // dwdp[50]
    dflux_v27_v_26_dSOCS3RNATurn = 0.27500000000000002*ActD*SOCS3RNAEqc*npSTAT5/init_STAT5;  // dwdp[51]
    dflux_v33_v_32_dSOCS3RNATurn = 0.40000000000000002*SOCS3RNA;  // dwdp[52]
    dflux_v34_v_33_dSOCS3Turn = 0.40000000000000002*SOCS3Eqc*SOCS3RNA/SOCS3RNAEqc;  // dwdp[53]
    dflux_v35_v_34_dSOCS3Turn = 0.40000000000000002*SOCS3;  // dwdp[54]
    dflux_v36_v_35_dSOCS3Turn = 0.40000000000000002*SOCS3Eqc*SOCS3EqcOE*SOCS3oe;  // dwdp[55]
    dflux_v14_v_13_dSTAT5ActEpoR = 0.40000000000000002*STAT5*std::pow(p12EpoRpJAK2 + p1EpoRpJAK2, 2)/(std::pow(init_EpoRJAK2, 2)*(CIS*CISInh/CISEqc + 1)*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));  // dwdp[56]
    dflux_v13_v_12_dSTAT5ActJAK2 = 0.40000000000000002*STAT5*(EpoRpJAK2 + p12EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2)/(init_EpoRJAK2*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));  // dwdp[57]
    dflux_v16_v_15_dSTAT5Exp = 0.27500000000000002*npSTAT5;  // dwdp[58]
    dflux_v15_v_14_dSTAT5Imp = 0.40000000000000002*pSTAT5;  // dwdp[59]
    dflux_v10_v_9_dinit_EpoRJAK2 = -0.40000000000000002*EpoRCISRemove*EpoRJAK2_CIS*(p12EpoRpJAK2 + p1EpoRpJAK2)/std::pow(init_EpoRJAK2, 2);  // dwdp[60]
    dflux_v11_v_10_dinit_EpoRJAK2 = -0.40000000000000002*SHP1*SHP1ActEpoR*(EpoRpJAK2 + p12EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2)/std::pow(init_EpoRJAK2, 2);  // dwdp[61]
    dflux_v13_v_12_dinit_EpoRJAK2 = -0.40000000000000002*STAT5*STAT5ActJAK2*(EpoRpJAK2 + p12EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2)/(std::pow(init_EpoRJAK2, 2)*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));  // dwdp[62]
    dflux_v14_v_13_dinit_EpoRJAK2 = -0.80000000000000004*STAT5*STAT5ActEpoR*std::pow(p12EpoRpJAK2 + p1EpoRpJAK2, 2)/(std::pow(init_EpoRJAK2, 3)*(CIS*CISInh/CISEqc + 1)*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));  // dwdp[63]
    dflux_v2_v_1_dinit_SHP1 = -0.40000000000000002*EpoRpJAK2*JAK2EpoRDeaSHP1*SHP1Act/std::pow(init_SHP1, 2);  // dwdp[64]
    dflux_v7_v_6_dinit_SHP1 = -0.40000000000000002*JAK2EpoRDeaSHP1*SHP1Act*p1EpoRpJAK2/std::pow(init_SHP1, 2);  // dwdp[65]
    dflux_v8_v_7_dinit_SHP1 = -0.40000000000000002*JAK2EpoRDeaSHP1*SHP1Act*p2EpoRpJAK2/std::pow(init_SHP1, 2);  // dwdp[66]
    dflux_v9_v_8_dinit_SHP1 = -0.40000000000000002*JAK2EpoRDeaSHP1*SHP1Act*p12EpoRpJAK2/std::pow(init_SHP1, 2);  // dwdp[67]
    dflux_v17_v_16_dinit_STAT5 = -0.27500000000000002*ActD*CISRNAEqc*CISRNATurn*npSTAT5/std::pow(init_STAT5, 2);  // dwdp[68]
    dflux_v27_v_26_dinit_STAT5 = -0.27500000000000002*ActD*SOCS3RNAEqc*SOCS3RNATurn*npSTAT5/std::pow(init_STAT5, 2);  // dwdp[69]
}

} // namespace model_Bachmann_MSB2011
} // namespace amici
