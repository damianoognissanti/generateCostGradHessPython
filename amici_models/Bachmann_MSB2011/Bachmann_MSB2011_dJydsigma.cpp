#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "Bachmann_MSB2011_p.h"
#include "Bachmann_MSB2011_k.h"
#include "Bachmann_MSB2011_y.h"
#include "Bachmann_MSB2011_sigmay.h"
#include "Bachmann_MSB2011_my.h"

namespace amici {
namespace model_Bachmann_MSB2011 {

void dJydsigma_Bachmann_MSB2011(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_observable_CISRNA_foldA - 1.0*std::pow(-std::log(mobservable_CISRNA_foldA)/M_LN10 + std::log(observable_CISRNA_foldA)/M_LN10, 2)/std::pow(sigma_observable_CISRNA_foldA, 3);
            break;
        case 1:
            dJydsigma[1] = 1.0/sigma_observable_CISRNA_foldB - 1.0*std::pow(-std::log(mobservable_CISRNA_foldB)/M_LN10 + std::log(observable_CISRNA_foldB)/M_LN10, 2)/std::pow(sigma_observable_CISRNA_foldB, 3);
            break;
        case 2:
            dJydsigma[2] = 1.0/sigma_observable_CISRNA_foldC - 1.0*std::pow(-std::log(mobservable_CISRNA_foldC)/M_LN10 + std::log(observable_CISRNA_foldC)/M_LN10, 2)/std::pow(sigma_observable_CISRNA_foldC, 3);
            break;
        case 3:
            dJydsigma[3] = 1.0/sigma_observable_CIS_abs - 1.0*std::pow(-std::log(mobservable_CIS_abs)/M_LN10 + std::log(observable_CIS_abs)/M_LN10, 2)/std::pow(sigma_observable_CIS_abs, 3);
            break;
        case 4:
            dJydsigma[4] = 1.0/sigma_observable_CIS_au - 1.0*std::pow(-std::log(mobservable_CIS_au)/M_LN10 + std::log(observable_CIS_au)/M_LN10, 2)/std::pow(sigma_observable_CIS_au, 3);
            break;
        case 5:
            dJydsigma[5] = 1.0/sigma_observable_CIS_au1 - 1.0*std::pow(-std::log(mobservable_CIS_au1)/M_LN10 + std::log(observable_CIS_au1)/M_LN10, 2)/std::pow(sigma_observable_CIS_au1, 3);
            break;
        case 6:
            dJydsigma[6] = 1.0/sigma_observable_CIS_au2 - 1.0*std::pow(-std::log(mobservable_CIS_au2)/M_LN10 + std::log(observable_CIS_au2)/M_LN10, 2)/std::pow(sigma_observable_CIS_au2, 3);
            break;
        case 7:
            dJydsigma[7] = 1.0/sigma_observable_SHP1_abs - 1.0*std::pow(-std::log(mobservable_SHP1_abs)/M_LN10 + std::log(observable_SHP1_abs)/M_LN10, 2)/std::pow(sigma_observable_SHP1_abs, 3);
            break;
        case 8:
            dJydsigma[8] = 1.0/sigma_observable_SOCS3RNA_foldA - 1.0*std::pow(-std::log(mobservable_SOCS3RNA_foldA)/M_LN10 + std::log(observable_SOCS3RNA_foldA)/M_LN10, 2)/std::pow(sigma_observable_SOCS3RNA_foldA, 3);
            break;
        case 9:
            dJydsigma[9] = 1.0/sigma_observable_SOCS3RNA_foldB - 1.0*std::pow(-std::log(mobservable_SOCS3RNA_foldB)/M_LN10 + std::log(observable_SOCS3RNA_foldB)/M_LN10, 2)/std::pow(sigma_observable_SOCS3RNA_foldB, 3);
            break;
        case 10:
            dJydsigma[10] = 1.0/sigma_observable_SOCS3RNA_foldC - 1.0*std::pow(-std::log(mobservable_SOCS3RNA_foldC)/M_LN10 + std::log(observable_SOCS3RNA_foldC)/M_LN10, 2)/std::pow(sigma_observable_SOCS3RNA_foldC, 3);
            break;
        case 11:
            dJydsigma[11] = 1.0/sigma_observable_SOCS3_abs - 1.0*std::pow(-std::log(mobservable_SOCS3_abs)/M_LN10 + std::log(observable_SOCS3_abs)/M_LN10, 2)/std::pow(sigma_observable_SOCS3_abs, 3);
            break;
        case 12:
            dJydsigma[12] = 1.0/sigma_observable_SOCS3_au - 1.0*std::pow(-std::log(mobservable_SOCS3_au)/M_LN10 + std::log(observable_SOCS3_au)/M_LN10, 2)/std::pow(sigma_observable_SOCS3_au, 3);
            break;
        case 13:
            dJydsigma[13] = 1.0/sigma_observable_STAT5_abs - 1.0*std::pow(-std::log(mobservable_STAT5_abs)/M_LN10 + std::log(observable_STAT5_abs)/M_LN10, 2)/std::pow(sigma_observable_STAT5_abs, 3);
            break;
        case 14:
            dJydsigma[14] = 1.0/sigma_observable_pEpoR_au - 1.0*std::pow(-std::log(mobservable_pEpoR_au)/M_LN10 + std::log(observable_pEpoR_au)/M_LN10, 2)/std::pow(sigma_observable_pEpoR_au, 3);
            break;
        case 15:
            dJydsigma[15] = 1.0/sigma_observable_pJAK2_au - 1.0*std::pow(-std::log(mobservable_pJAK2_au)/M_LN10 + std::log(observable_pJAK2_au)/M_LN10, 2)/std::pow(sigma_observable_pJAK2_au, 3);
            break;
        case 16:
            dJydsigma[16] = 1.0/sigma_observable_pSTAT5B_rel - 1.0*std::pow(-mobservable_pSTAT5B_rel + observable_pSTAT5B_rel, 2)/std::pow(sigma_observable_pSTAT5B_rel, 3);
            break;
        case 17:
            dJydsigma[17] = 1.0/sigma_observable_pSTAT5_au - 1.0*std::pow(-std::log(mobservable_pSTAT5_au)/M_LN10 + std::log(observable_pSTAT5_au)/M_LN10, 2)/std::pow(sigma_observable_pSTAT5_au, 3);
            break;
        case 18:
            dJydsigma[18] = 1.0/sigma_observable_tSHP1_au - 1.0*std::pow(-std::log(mobservable_tSHP1_au)/M_LN10 + std::log(observable_tSHP1_au)/M_LN10, 2)/std::pow(sigma_observable_tSHP1_au, 3);
            break;
        case 19:
            dJydsigma[19] = 1.0/sigma_observable_tSTAT5_au - 1.0*std::pow(-std::log(mobservable_tSTAT5_au)/M_LN10 + std::log(observable_tSTAT5_au)/M_LN10, 2)/std::pow(sigma_observable_tSTAT5_au, 3);
            break;
    }
}

} // namespace model_Bachmann_MSB2011
} // namespace amici
