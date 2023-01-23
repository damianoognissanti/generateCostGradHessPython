#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "Bachmann_MSB2011_x.h"
#include "Bachmann_MSB2011_p.h"
#include "Bachmann_MSB2011_k.h"
#include "Bachmann_MSB2011_w.h"
#include "Bachmann_MSB2011_dxdotdw.h"

namespace amici {
namespace model_Bachmann_MSB2011 {

void dxdotdw_Bachmann_MSB2011(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dxdot0_dflux_v1_v_0 = -2.5;  // dxdotdw[0]
    dxdot1_dflux_v1_v_0 = 2.5;  // dxdotdw[1]
    dxdot0_dflux_v2_v_1 = 2.5;  // dxdotdw[2]
    dxdot1_dflux_v2_v_1 = -2.5;  // dxdotdw[3]
    dxdot1_dflux_v3_v_2 = -2.5;  // dxdotdw[4]
    dxdot2_dflux_v3_v_2 = 2.5;  // dxdotdw[5]
    dxdot1_dflux_v4_v_3 = -2.5;  // dxdotdw[6]
    dxdot3_dflux_v4_v_3 = 2.5;  // dxdotdw[7]
    dxdot2_dflux_v5_v_4 = -2.5;  // dxdotdw[8]
    dxdot4_dflux_v5_v_4 = 2.5;  // dxdotdw[9]
    dxdot3_dflux_v6_v_5 = -2.5;  // dxdotdw[10]
    dxdot4_dflux_v6_v_5 = 2.5;  // dxdotdw[11]
    dxdot0_dflux_v7_v_6 = 2.5;  // dxdotdw[12]
    dxdot2_dflux_v7_v_6 = -2.5;  // dxdotdw[13]
    dxdot0_dflux_v8_v_7 = 2.5;  // dxdotdw[14]
    dxdot3_dflux_v8_v_7 = -2.5;  // dxdotdw[15]
    dxdot0_dflux_v9_v_8 = 2.5;  // dxdotdw[16]
    dxdot4_dflux_v9_v_8 = -2.5;  // dxdotdw[17]
    dxdot5_dflux_v10_v_9 = -2.5;  // dxdotdw[18]
    dxdot6_dflux_v11_v_10 = -2.5;  // dxdotdw[19]
    dxdot7_dflux_v11_v_10 = 2.5;  // dxdotdw[20]
    dxdot6_dflux_v12_v_11 = 2.5;  // dxdotdw[21]
    dxdot7_dflux_v12_v_11 = -2.5;  // dxdotdw[22]
    dxdot8_dflux_v13_v_12 = -2.5;  // dxdotdw[23]
    dxdot9_dflux_v13_v_12 = 2.5;  // dxdotdw[24]
    dxdot8_dflux_v14_v_13 = -2.5;  // dxdotdw[25]
    dxdot9_dflux_v14_v_13 = 2.5;  // dxdotdw[26]
    dxdot9_dflux_v15_v_14 = -2.5;  // dxdotdw[27]
    dxdot10_dflux_v15_v_14 = 3.6363636363636362;  // dxdotdw[28]
    dxdot8_dflux_v16_v_15 = 2.5;  // dxdotdw[29]
    dxdot10_dflux_v16_v_15 = -3.6363636363636362;  // dxdotdw[30]
    dxdot11_dflux_v17_v_16 = 3.6363636363636362;  // dxdotdw[31]
    dxdot11_dflux_v18_v_17 = -3.6363636363636362;  // dxdotdw[32]
    dxdot12_dflux_v18_v_17 = 3.6363636363636362;  // dxdotdw[33]
    dxdot12_dflux_v19_v_18 = -3.6363636363636362;  // dxdotdw[34]
    dxdot13_dflux_v19_v_18 = 3.6363636363636362;  // dxdotdw[35]
    dxdot13_dflux_v20_v_19 = -3.6363636363636362;  // dxdotdw[36]
    dxdot14_dflux_v20_v_19 = 3.6363636363636362;  // dxdotdw[37]
    dxdot14_dflux_v21_v_20 = -3.6363636363636362;  // dxdotdw[38]
    dxdot15_dflux_v21_v_20 = 3.6363636363636362;  // dxdotdw[39]
    dxdot15_dflux_v22_v_21 = -3.6363636363636362;  // dxdotdw[40]
    dxdot16_dflux_v22_v_21 = 2.5;  // dxdotdw[41]
    dxdot16_dflux_v23_v_22 = -2.5;  // dxdotdw[42]
    dxdot17_dflux_v24_v_23 = 2.5;  // dxdotdw[43]
    dxdot17_dflux_v25_v_24 = -2.5;  // dxdotdw[44]
    dxdot17_dflux_v26_v_25 = 2.5;  // dxdotdw[45]
    dxdot18_dflux_v27_v_26 = 3.6363636363636362;  // dxdotdw[46]
    dxdot18_dflux_v28_v_27 = -3.6363636363636362;  // dxdotdw[47]
    dxdot19_dflux_v28_v_27 = 3.6363636363636362;  // dxdotdw[48]
    dxdot19_dflux_v29_v_28 = -3.6363636363636362;  // dxdotdw[49]
    dxdot20_dflux_v29_v_28 = 3.6363636363636362;  // dxdotdw[50]
    dxdot20_dflux_v30_v_29 = -3.6363636363636362;  // dxdotdw[51]
    dxdot21_dflux_v30_v_29 = 3.6363636363636362;  // dxdotdw[52]
    dxdot21_dflux_v31_v_30 = -3.6363636363636362;  // dxdotdw[53]
    dxdot22_dflux_v31_v_30 = 3.6363636363636362;  // dxdotdw[54]
    dxdot22_dflux_v32_v_31 = -3.6363636363636362;  // dxdotdw[55]
    dxdot23_dflux_v32_v_31 = 2.5;  // dxdotdw[56]
    dxdot23_dflux_v33_v_32 = -2.5;  // dxdotdw[57]
    dxdot24_dflux_v34_v_33 = 2.5;  // dxdotdw[58]
    dxdot24_dflux_v35_v_34 = -2.5;  // dxdotdw[59]
    dxdot24_dflux_v36_v_35 = 2.5;  // dxdotdw[60]
}

} // namespace model_Bachmann_MSB2011
} // namespace amici
