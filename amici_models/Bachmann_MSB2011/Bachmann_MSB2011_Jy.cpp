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

void Jy_Bachmann_MSB2011(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_CISRNA_foldA, 2)*std::pow(sigma_observable_CISRNA_foldA, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(-std::log(mobservable_CISRNA_foldA)/M_LN10 + std::log(observable_CISRNA_foldA)/M_LN10, 2)/std::pow(sigma_observable_CISRNA_foldA, 2);
            break;
        case 1:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_CISRNA_foldB, 2)*std::pow(sigma_observable_CISRNA_foldB, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(-std::log(mobservable_CISRNA_foldB)/M_LN10 + std::log(observable_CISRNA_foldB)/M_LN10, 2)/std::pow(sigma_observable_CISRNA_foldB, 2);
            break;
        case 2:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_CISRNA_foldC, 2)*std::pow(sigma_observable_CISRNA_foldC, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(-std::log(mobservable_CISRNA_foldC)/M_LN10 + std::log(observable_CISRNA_foldC)/M_LN10, 2)/std::pow(sigma_observable_CISRNA_foldC, 2);
            break;
        case 3:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_CIS_abs, 2)*std::pow(sigma_observable_CIS_abs, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(-std::log(mobservable_CIS_abs)/M_LN10 + std::log(observable_CIS_abs)/M_LN10, 2)/std::pow(sigma_observable_CIS_abs, 2);
            break;
        case 4:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_CIS_au, 2)*std::pow(sigma_observable_CIS_au, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(-std::log(mobservable_CIS_au)/M_LN10 + std::log(observable_CIS_au)/M_LN10, 2)/std::pow(sigma_observable_CIS_au, 2);
            break;
        case 5:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_CIS_au1, 2)*std::pow(sigma_observable_CIS_au1, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(-std::log(mobservable_CIS_au1)/M_LN10 + std::log(observable_CIS_au1)/M_LN10, 2)/std::pow(sigma_observable_CIS_au1, 2);
            break;
        case 6:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_CIS_au2, 2)*std::pow(sigma_observable_CIS_au2, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(-std::log(mobservable_CIS_au2)/M_LN10 + std::log(observable_CIS_au2)/M_LN10, 2)/std::pow(sigma_observable_CIS_au2, 2);
            break;
        case 7:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_SHP1_abs, 2)*std::pow(sigma_observable_SHP1_abs, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(-std::log(mobservable_SHP1_abs)/M_LN10 + std::log(observable_SHP1_abs)/M_LN10, 2)/std::pow(sigma_observable_SHP1_abs, 2);
            break;
        case 8:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_SOCS3RNA_foldA, 2)*std::pow(sigma_observable_SOCS3RNA_foldA, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(-std::log(mobservable_SOCS3RNA_foldA)/M_LN10 + std::log(observable_SOCS3RNA_foldA)/M_LN10, 2)/std::pow(sigma_observable_SOCS3RNA_foldA, 2);
            break;
        case 9:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_SOCS3RNA_foldB, 2)*std::pow(sigma_observable_SOCS3RNA_foldB, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(-std::log(mobservable_SOCS3RNA_foldB)/M_LN10 + std::log(observable_SOCS3RNA_foldB)/M_LN10, 2)/std::pow(sigma_observable_SOCS3RNA_foldB, 2);
            break;
        case 10:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_SOCS3RNA_foldC, 2)*std::pow(sigma_observable_SOCS3RNA_foldC, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(-std::log(mobservable_SOCS3RNA_foldC)/M_LN10 + std::log(observable_SOCS3RNA_foldC)/M_LN10, 2)/std::pow(sigma_observable_SOCS3RNA_foldC, 2);
            break;
        case 11:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_SOCS3_abs, 2)*std::pow(sigma_observable_SOCS3_abs, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(-std::log(mobservable_SOCS3_abs)/M_LN10 + std::log(observable_SOCS3_abs)/M_LN10, 2)/std::pow(sigma_observable_SOCS3_abs, 2);
            break;
        case 12:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_SOCS3_au, 2)*std::pow(sigma_observable_SOCS3_au, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(-std::log(mobservable_SOCS3_au)/M_LN10 + std::log(observable_SOCS3_au)/M_LN10, 2)/std::pow(sigma_observable_SOCS3_au, 2);
            break;
        case 13:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_STAT5_abs, 2)*std::pow(sigma_observable_STAT5_abs, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(-std::log(mobservable_STAT5_abs)/M_LN10 + std::log(observable_STAT5_abs)/M_LN10, 2)/std::pow(sigma_observable_STAT5_abs, 2);
            break;
        case 14:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_pEpoR_au, 2)*std::pow(sigma_observable_pEpoR_au, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(-std::log(mobservable_pEpoR_au)/M_LN10 + std::log(observable_pEpoR_au)/M_LN10, 2)/std::pow(sigma_observable_pEpoR_au, 2);
            break;
        case 15:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_pJAK2_au, 2)*std::pow(sigma_observable_pJAK2_au, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(-std::log(mobservable_pJAK2_au)/M_LN10 + std::log(observable_pJAK2_au)/M_LN10, 2)/std::pow(sigma_observable_pJAK2_au, 2);
            break;
        case 16:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_pSTAT5B_rel, 2)) + 0.5*std::pow(-mobservable_pSTAT5B_rel + observable_pSTAT5B_rel, 2)/std::pow(sigma_observable_pSTAT5B_rel, 2);
            break;
        case 17:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_pSTAT5_au, 2)*std::pow(sigma_observable_pSTAT5_au, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(-std::log(mobservable_pSTAT5_au)/M_LN10 + std::log(observable_pSTAT5_au)/M_LN10, 2)/std::pow(sigma_observable_pSTAT5_au, 2);
            break;
        case 18:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_tSHP1_au, 2)*std::pow(sigma_observable_tSHP1_au, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(-std::log(mobservable_tSHP1_au)/M_LN10 + std::log(observable_tSHP1_au)/M_LN10, 2)/std::pow(sigma_observable_tSHP1_au, 2);
            break;
        case 19:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_tSTAT5_au, 2)*std::pow(sigma_observable_tSTAT5_au, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(-std::log(mobservable_tSTAT5_au)/M_LN10 + std::log(observable_tSTAT5_au)/M_LN10, 2)/std::pow(sigma_observable_tSTAT5_au, 2);
            break;
    }
}

} // namespace model_Bachmann_MSB2011
} // namespace amici
