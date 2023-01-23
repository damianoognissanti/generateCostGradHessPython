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

void dwdx_Bachmann_MSB2011(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dflux_v1_v_0_dEpoRJAK2 = 0.40000000000000002*Epo*JAK2ActEpo/(SOCS3*SOCS3Inh/SOCS3Eqc + 1);  // dwdx[0]
    dflux_v2_v_1_dEpoRpJAK2 = 0.40000000000000002*JAK2EpoRDeaSHP1*SHP1Act/init_SHP1;  // dwdx[1]
    dflux_v3_v_2_dEpoRpJAK2 = 0.40000000000000002*EpoRActJAK2/(SOCS3*SOCS3Inh/SOCS3Eqc + 1);  // dwdx[2]
    dflux_v4_v_3_dEpoRpJAK2 = 1.2000000000000002*EpoRActJAK2/((EpoRCISInh*EpoRJAK2_CIS + 1)*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));  // dwdx[3]
    dflux_v11_v_10_dEpoRpJAK2 = 0.40000000000000002*SHP1*SHP1ActEpoR/init_EpoRJAK2;  // dwdx[4]
    dflux_v13_v_12_dEpoRpJAK2 = 0.40000000000000002*STAT5*STAT5ActJAK2/(init_EpoRJAK2*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));  // dwdx[5]
    dflux_v5_v_4_dp1EpoRpJAK2 = 1.2000000000000002*EpoRActJAK2/((EpoRCISInh*EpoRJAK2_CIS + 1)*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));  // dwdx[6]
    dflux_v7_v_6_dp1EpoRpJAK2 = 0.40000000000000002*JAK2EpoRDeaSHP1*SHP1Act/init_SHP1;  // dwdx[7]
    dflux_v10_v_9_dp1EpoRpJAK2 = 0.40000000000000002*EpoRCISRemove*EpoRJAK2_CIS/init_EpoRJAK2;  // dwdx[8]
    dflux_v11_v_10_dp1EpoRpJAK2 = 0.40000000000000002*SHP1*SHP1ActEpoR/init_EpoRJAK2;  // dwdx[9]
    dflux_v13_v_12_dp1EpoRpJAK2 = 0.40000000000000002*STAT5*STAT5ActJAK2/(init_EpoRJAK2*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));  // dwdx[10]
    dflux_v14_v_13_dp1EpoRpJAK2 = 0.40000000000000002*STAT5*STAT5ActEpoR*(2*p12EpoRpJAK2 + 2*p1EpoRpJAK2)/(std::pow(init_EpoRJAK2, 2)*(CIS*CISInh/CISEqc + 1)*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));  // dwdx[11]
    dflux_v6_v_5_dp2EpoRpJAK2 = 0.40000000000000002*EpoRActJAK2/(SOCS3*SOCS3Inh/SOCS3Eqc + 1);  // dwdx[12]
    dflux_v8_v_7_dp2EpoRpJAK2 = 0.40000000000000002*JAK2EpoRDeaSHP1*SHP1Act/init_SHP1;  // dwdx[13]
    dflux_v11_v_10_dp2EpoRpJAK2 = 0.40000000000000002*SHP1*SHP1ActEpoR/init_EpoRJAK2;  // dwdx[14]
    dflux_v13_v_12_dp2EpoRpJAK2 = 0.40000000000000002*STAT5*STAT5ActJAK2/(init_EpoRJAK2*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));  // dwdx[15]
    dflux_v9_v_8_dp12EpoRpJAK2 = 0.40000000000000002*JAK2EpoRDeaSHP1*SHP1Act/init_SHP1;  // dwdx[16]
    dflux_v10_v_9_dp12EpoRpJAK2 = 0.40000000000000002*EpoRCISRemove*EpoRJAK2_CIS/init_EpoRJAK2;  // dwdx[17]
    dflux_v11_v_10_dp12EpoRpJAK2 = 0.40000000000000002*SHP1*SHP1ActEpoR/init_EpoRJAK2;  // dwdx[18]
    dflux_v13_v_12_dp12EpoRpJAK2 = 0.40000000000000002*STAT5*STAT5ActJAK2/(init_EpoRJAK2*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));  // dwdx[19]
    dflux_v14_v_13_dp12EpoRpJAK2 = 0.40000000000000002*STAT5*STAT5ActEpoR*(2*p12EpoRpJAK2 + 2*p1EpoRpJAK2)/(std::pow(init_EpoRJAK2, 2)*(CIS*CISInh/CISEqc + 1)*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));  // dwdx[20]
    dflux_v4_v_3_dEpoRJAK2_CIS = -1.2000000000000002*EpoRActJAK2*EpoRCISInh*EpoRpJAK2/(std::pow(EpoRCISInh*EpoRJAK2_CIS + 1, 2)*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));  // dwdx[21]
    dflux_v5_v_4_dEpoRJAK2_CIS = -1.2000000000000002*EpoRActJAK2*EpoRCISInh*p1EpoRpJAK2/(std::pow(EpoRCISInh*EpoRJAK2_CIS + 1, 2)*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));  // dwdx[22]
    dflux_v10_v_9_dEpoRJAK2_CIS = 0.40000000000000002*EpoRCISRemove*(p12EpoRpJAK2 + p1EpoRpJAK2)/init_EpoRJAK2;  // dwdx[23]
    dflux_v11_v_10_dSHP1 = 0.40000000000000002*SHP1ActEpoR*(EpoRpJAK2 + p12EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2)/init_EpoRJAK2;  // dwdx[24]
    dflux_v2_v_1_dSHP1Act = 0.40000000000000002*EpoRpJAK2*JAK2EpoRDeaSHP1/init_SHP1;  // dwdx[25]
    dflux_v7_v_6_dSHP1Act = 0.40000000000000002*JAK2EpoRDeaSHP1*p1EpoRpJAK2/init_SHP1;  // dwdx[26]
    dflux_v8_v_7_dSHP1Act = 0.40000000000000002*JAK2EpoRDeaSHP1*p2EpoRpJAK2/init_SHP1;  // dwdx[27]
    dflux_v9_v_8_dSHP1Act = 0.40000000000000002*JAK2EpoRDeaSHP1*p12EpoRpJAK2/init_SHP1;  // dwdx[28]
    dflux_v12_v_11_dSHP1Act = 0.40000000000000002*SHP1Dea;  // dwdx[29]
    dflux_v13_v_12_dSTAT5 = 0.40000000000000002*STAT5ActJAK2*(EpoRpJAK2 + p12EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2)/(init_EpoRJAK2*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));  // dwdx[30]
    dflux_v14_v_13_dSTAT5 = 0.40000000000000002*STAT5ActEpoR*std::pow(p12EpoRpJAK2 + p1EpoRpJAK2, 2)/(std::pow(init_EpoRJAK2, 2)*(CIS*CISInh/CISEqc + 1)*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));  // dwdx[31]
    dflux_v15_v_14_dpSTAT5 = 0.40000000000000002*STAT5Imp;  // dwdx[32]
    dflux_v16_v_15_dnpSTAT5 = 0.27500000000000002*STAT5Exp;  // dwdx[33]
    dflux_v17_v_16_dnpSTAT5 = 0.27500000000000002*ActD*CISRNAEqc*CISRNATurn/init_STAT5;  // dwdx[34]
    dflux_v27_v_26_dnpSTAT5 = 0.27500000000000002*ActD*SOCS3RNAEqc*SOCS3RNATurn/init_STAT5;  // dwdx[35]
    dflux_v18_v_17_dCISnRNA1 = 0.27500000000000002*CISRNADelay;  // dwdx[36]
    dflux_v19_v_18_dCISnRNA2 = 0.27500000000000002*CISRNADelay;  // dwdx[37]
    dflux_v20_v_19_dCISnRNA3 = 0.27500000000000002*CISRNADelay;  // dwdx[38]
    dflux_v21_v_20_dCISnRNA4 = 0.27500000000000002*CISRNADelay;  // dwdx[39]
    dflux_v22_v_21_dCISnRNA5 = 0.27500000000000002*CISRNADelay;  // dwdx[40]
    dflux_v23_v_22_dCISRNA = 0.40000000000000002*CISRNATurn;  // dwdx[41]
    dflux_v24_v_23_dCISRNA = 0.40000000000000002*CISEqc*CISTurn/CISRNAEqc;  // dwdx[42]
    dflux_v14_v_13_dCIS = -0.40000000000000002*CISInh*STAT5*STAT5ActEpoR*std::pow(p12EpoRpJAK2 + p1EpoRpJAK2, 2)/(CISEqc*std::pow(init_EpoRJAK2, 2)*std::pow(CIS*CISInh/CISEqc + 1, 2)*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));  // dwdx[43]
    dflux_v25_v_24_dCIS = 0.40000000000000002*CISTurn;  // dwdx[44]
    dflux_v28_v_27_dSOCS3nRNA1 = 0.27500000000000002*SOCS3RNADelay;  // dwdx[45]
    dflux_v29_v_28_dSOCS3nRNA2 = 0.27500000000000002*SOCS3RNADelay;  // dwdx[46]
    dflux_v30_v_29_dSOCS3nRNA3 = 0.27500000000000002*SOCS3RNADelay;  // dwdx[47]
    dflux_v31_v_30_dSOCS3nRNA4 = 0.27500000000000002*SOCS3RNADelay;  // dwdx[48]
    dflux_v32_v_31_dSOCS3nRNA5 = 0.27500000000000002*SOCS3RNADelay;  // dwdx[49]
    dflux_v33_v_32_dSOCS3RNA = 0.40000000000000002*SOCS3RNATurn;  // dwdx[50]
    dflux_v34_v_33_dSOCS3RNA = 0.40000000000000002*SOCS3Eqc*SOCS3Turn/SOCS3RNAEqc;  // dwdx[51]
    dflux_v1_v_0_dSOCS3 = -0.40000000000000002*Epo*EpoRJAK2*JAK2ActEpo*SOCS3Inh/(SOCS3Eqc*std::pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));  // dwdx[52]
    dflux_v3_v_2_dSOCS3 = -0.40000000000000002*EpoRActJAK2*EpoRpJAK2*SOCS3Inh/(SOCS3Eqc*std::pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));  // dwdx[53]
    dflux_v4_v_3_dSOCS3 = -1.2000000000000002*EpoRActJAK2*EpoRpJAK2*SOCS3Inh/(SOCS3Eqc*(EpoRCISInh*EpoRJAK2_CIS + 1)*std::pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));  // dwdx[54]
    dflux_v5_v_4_dSOCS3 = -1.2000000000000002*EpoRActJAK2*SOCS3Inh*p1EpoRpJAK2/(SOCS3Eqc*(EpoRCISInh*EpoRJAK2_CIS + 1)*std::pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));  // dwdx[55]
    dflux_v6_v_5_dSOCS3 = -0.40000000000000002*EpoRActJAK2*SOCS3Inh*p2EpoRpJAK2/(SOCS3Eqc*std::pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));  // dwdx[56]
    dflux_v13_v_12_dSOCS3 = -0.40000000000000002*SOCS3Inh*STAT5*STAT5ActJAK2*(EpoRpJAK2 + p12EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2)/(SOCS3Eqc*init_EpoRJAK2*std::pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));  // dwdx[57]
    dflux_v14_v_13_dSOCS3 = -0.40000000000000002*SOCS3Inh*STAT5*STAT5ActEpoR*std::pow(p12EpoRpJAK2 + p1EpoRpJAK2, 2)/(SOCS3Eqc*std::pow(init_EpoRJAK2, 2)*(CIS*CISInh/CISEqc + 1)*std::pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));  // dwdx[58]
    dflux_v35_v_34_dSOCS3 = 0.40000000000000002*SOCS3Turn;  // dwdx[59]
}

} // namespace model_Bachmann_MSB2011
} // namespace amici
