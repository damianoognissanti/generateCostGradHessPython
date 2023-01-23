#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "Bachmann_MSB2011_x.h"
#include "Bachmann_MSB2011_p.h"
#include "Bachmann_MSB2011_k.h"
#include "Bachmann_MSB2011_w.h"
#include "Bachmann_MSB2011_xdot.h"

namespace amici {
namespace model_Bachmann_MSB2011 {

void xdot_Bachmann_MSB2011(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    xdot0 = -2.5*flux_v1_v_0 + 2.5*flux_v2_v_1 + 2.5*flux_v7_v_6 + 2.5*flux_v8_v_7 + 2.5*flux_v9_v_8;  // xdot[0]
    xdot1 = 2.5*flux_v1_v_0 - 2.5*flux_v2_v_1 - 2.5*flux_v3_v_2 - 2.5*flux_v4_v_3;  // xdot[1]
    xdot2 = 2.5*flux_v3_v_2 - 2.5*flux_v5_v_4 - 2.5*flux_v7_v_6;  // xdot[2]
    xdot3 = 2.5*flux_v4_v_3 - 2.5*flux_v6_v_5 - 2.5*flux_v8_v_7;  // xdot[3]
    xdot4 = 2.5*flux_v5_v_4 + 2.5*flux_v6_v_5 - 2.5*flux_v9_v_8;  // xdot[4]
    xdot5 = -2.5*flux_v10_v_9;  // xdot[5]
    xdot6 = -2.5*flux_v11_v_10 + 2.5*flux_v12_v_11;  // xdot[6]
    xdot7 = 2.5*flux_v11_v_10 - 2.5*flux_v12_v_11;  // xdot[7]
    xdot8 = -2.5*flux_v13_v_12 - 2.5*flux_v14_v_13 + 2.5*flux_v16_v_15;  // xdot[8]
    xdot9 = 2.5*flux_v13_v_12 + 2.5*flux_v14_v_13 - 2.5*flux_v15_v_14;  // xdot[9]
    xdot10 = 3.6363636363636362*flux_v15_v_14 - 3.6363636363636362*flux_v16_v_15;  // xdot[10]
    xdot11 = 3.6363636363636362*flux_v17_v_16 - 3.6363636363636362*flux_v18_v_17;  // xdot[11]
    xdot12 = 3.6363636363636362*flux_v18_v_17 - 3.6363636363636362*flux_v19_v_18;  // xdot[12]
    xdot13 = 3.6363636363636362*flux_v19_v_18 - 3.6363636363636362*flux_v20_v_19;  // xdot[13]
    xdot14 = 3.6363636363636362*flux_v20_v_19 - 3.6363636363636362*flux_v21_v_20;  // xdot[14]
    xdot15 = 3.6363636363636362*flux_v21_v_20 - 3.6363636363636362*flux_v22_v_21;  // xdot[15]
    xdot16 = 2.5*flux_v22_v_21 - 2.5*flux_v23_v_22;  // xdot[16]
    xdot17 = 2.5*flux_v24_v_23 - 2.5*flux_v25_v_24 + 2.5*flux_v26_v_25;  // xdot[17]
    xdot18 = 3.6363636363636362*flux_v27_v_26 - 3.6363636363636362*flux_v28_v_27;  // xdot[18]
    xdot19 = 3.6363636363636362*flux_v28_v_27 - 3.6363636363636362*flux_v29_v_28;  // xdot[19]
    xdot20 = 3.6363636363636362*flux_v29_v_28 - 3.6363636363636362*flux_v30_v_29;  // xdot[20]
    xdot21 = 3.6363636363636362*flux_v30_v_29 - 3.6363636363636362*flux_v31_v_30;  // xdot[21]
    xdot22 = 3.6363636363636362*flux_v31_v_30 - 3.6363636363636362*flux_v32_v_31;  // xdot[22]
    xdot23 = 2.5*flux_v32_v_31 - 2.5*flux_v33_v_32;  // xdot[23]
    xdot24 = 2.5*flux_v34_v_33 - 2.5*flux_v35_v_34 + 2.5*flux_v36_v_35;  // xdot[24]
}

} // namespace model_Bachmann_MSB2011
} // namespace amici
