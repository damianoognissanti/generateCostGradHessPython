#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>
#include <algorithm>

#include "Schwen_PONE2014_x.h"
#include "Schwen_PONE2014_p.h"
#include "Schwen_PONE2014_k.h"
#include "Schwen_PONE2014_w.h"
#include "Schwen_PONE2014_dwdp.h"

namespace amici {
namespace model_Schwen_PONE2014 {

void dwdp_Schwen_PONE2014(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp){
    dflux_v1_ka1_dka1 = 1.0*Ins*Rec1;  // dwdp[0]
    dflux_v2_ka2fold_dka1 = 1.0*Ins*Rec2*ka2fold;  // dwdp[1]
    dflux_v13_v_12_dka1 = 1.0*Ins*Rec1;  // dwdp[2]
    dflux_v14_v_13_dka1 = 1.0*Ins*Rec2*ka2fold;  // dwdp[3]
    dflux_v2_ka2fold_dka2fold = 1.0*Ins*Rec2*ka1;  // dwdp[4]
    dflux_v14_v_13_dka2fold = 1.0*Ins*Rec2*ka1;  // dwdp[5]
    dflux_v5_kd1_dkd1 = 1.0*IR1;  // dwdp[6]
    dflux_v6_kd2fold_dkd1 = 1.0*IR2*kd2fold;  // dwdp[7]
    dflux_v13_v_12_dkd1 = -1.0*IR1;  // dwdp[8]
    dflux_v14_v_13_dkd1 = -1.0*IR2*kd2fold;  // dwdp[9]
    dflux_v6_kd2fold_dkd2fold = 1.0*IR2*kd1;  // dwdp[10]
    dflux_v14_v_13_dkd2fold = -1.0*IR2*kd1;  // dwdp[11]
    dflux_v7_v_6_dkin = 1.0*IR1;  // dwdp[12]
    dflux_v8_v_7_dkin2 = 1.0*IR2;  // dwdp[13]
    dflux_v4_v_3_dkoff_unspec = 1.0*BoundUnspec;  // dwdp[14]
    dflux_v3_v_2_dkon_unspec = 1.0*Ins;  // dwdp[15]
    dflux_v9_v_8_dkout = 1.0*IR1in;  // dwdp[16]
    dflux_v10_v_9_dkout2 = 1.0*IR2in;  // dwdp[17]
    dflux_v11_v_10_dkout_frag = 1.0*IR1in;  // dwdp[18]
    dflux_v12_v_11_dkout_frag = 1.0*IR2in;  // dwdp[19]
}

} // namespace model_Schwen_PONE2014
} // namespace amici
