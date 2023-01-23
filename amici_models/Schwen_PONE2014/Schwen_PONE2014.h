#ifndef _amici_TPL_MODELNAME_h
#define _amici_TPL_MODELNAME_h
#include <cmath>
#include <memory>
#include <gsl/gsl-lite.hpp>

#include "amici/model_ode.h"
#include "amici/solver_cvodes.h"

#include "sundials/sundials_types.h"

namespace amici {

class Solver;

namespace model_Schwen_PONE2014 {

extern std::array<const char*, 27> parameterNames;
extern std::array<const char*, 1> fixedParameterNames;
extern std::array<const char*, 11> stateNames;
extern std::array<const char*, 4> observableNames;
extern std::array<const ObservableScaling, 4> observableScalings;
extern std::array<const char*, 14> expressionNames;
extern std::array<const char*, 27> parameterIds;
extern std::array<const char*, 1> fixedParameterIds;
extern std::array<const char*, 11> stateIds;
extern std::array<const char*, 4> observableIds;
extern std::array<const char*, 14> expressionIds;
extern std::array<int, 11> stateIdxsSolver;
extern std::array<bool, 0> rootInitialValues;

extern void Jy_Schwen_PONE2014(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydsigma_Schwen_PONE2014(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydy_Schwen_PONE2014(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydy_colptrs_Schwen_PONE2014(SUNMatrixWrapper &colptrs, int index);
extern void dJydy_rowvals_Schwen_PONE2014(SUNMatrixWrapper &rowvals, int index);

extern void dwdp_Schwen_PONE2014(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp);
extern void dwdp_colptrs_Schwen_PONE2014(SUNMatrixWrapper &colptrs);
extern void dwdp_rowvals_Schwen_PONE2014(SUNMatrixWrapper &rowvals);
extern void dwdx_Schwen_PONE2014(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl);
extern void dwdx_colptrs_Schwen_PONE2014(SUNMatrixWrapper &colptrs);
extern void dwdx_rowvals_Schwen_PONE2014(SUNMatrixWrapper &rowvals);



extern void dxdotdw_Schwen_PONE2014(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void dxdotdw_colptrs_Schwen_PONE2014(SUNMatrixWrapper &colptrs);
extern void dxdotdw_rowvals_Schwen_PONE2014(SUNMatrixWrapper &rowvals);






extern void dydx_Schwen_PONE2014(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx);
extern void dydp_Schwen_PONE2014(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *tcl, const realtype *dtcldp);
extern void sigmay_Schwen_PONE2014(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y);
extern void dsigmaydp_Schwen_PONE2014(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip);

extern void w_Schwen_PONE2014(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl);
extern void x0_Schwen_PONE2014(realtype *x0, const realtype t, const realtype *p, const realtype *k);
extern void x0_fixedParameters_Schwen_PONE2014(realtype *x0_fixedParameters, const realtype t, const realtype *p, const realtype *k, gsl::span<const int> reinitialization_state_idxs);
extern void sx0_Schwen_PONE2014(realtype *sx0, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip);
extern void sx0_fixedParameters_Schwen_PONE2014(realtype *sx0_fixedParameters, const realtype t, const realtype *x0, const realtype *p, const realtype *k, const int ip, gsl::span<const int> reinitialization_state_idxs);
extern void xdot_Schwen_PONE2014(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void y_Schwen_PONE2014(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);




extern void x_solver_Schwen_PONE2014(realtype *x_solver, const realtype *x_rdata);












/**
 * @brief AMICI-generated model subclass.
 */
class Model_Schwen_PONE2014 : public amici::Model_ODE {
  public:
    /**
     * @brief Default constructor.
     */
    Model_Schwen_PONE2014()
        : amici::Model_ODE(
              amici::ModelDimensions(
                  11,                            // nx_rdata
                  11,                        // nxtrue_rdata
                  11,                           // nx_solver
                  11,                       // nxtrue_solver
                  0,                    // nx_solver_reinit
                  27,                                  // np
                  1,                                  // nk
                  4,                                  // ny
                  4,                              // nytrue
                  0,                                  // nz
                  0,                              // nztrue
                  0,                              // nevent
                  1,                          // nobjective
                  14,                                  // nw
                  20,                               // ndwdx
                  20,                               // ndwdp
                  0,                               // ndwdw
                  32,                            // ndxdotdw
                  std::vector<int>{1,1,1,1},                              // ndjydy
                  0,                    // ndxrdatadxsolver
                  0,                        // ndxrdatadtcl
                  0,                        // ndtotal_cldx_rdata
                  0,                                       // nnz
                  11,                                 // ubw
                  11                                  // lbw
              ),
              amici::SimulationParameters(
                  std::vector<realtype>{0.0}, // fixedParameters
                  std::vector<realtype>{47.237017282009603, 16.457631927604002, 0.0093798043666388307, 2.0907692381484, 6.7226916916103399, 9.6185010765549297, 3.75654890101743, 0.54530402971450898, 9.6457882009227003, 19.941427249128001, 324.72583814527798, 0.052907956097648701, 0.010042166937868901, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}        // dynamic parameters
              ),
              amici::SecondOrderMode::none,                                  // o2mode
              std::vector<realtype>(11, 0.0),   // idlist
              std::vector<int>{},                          // z2event
              true,                                        // pythonGenerated
              0,                       // ndxdotdp_explicit
              0,                       // ndxdotdx_explicit
              0                        // w_recursion_depth
          ) {
                 root_initial_values_ = std::vector<bool>(
                     rootInitialValues.begin(), rootInitialValues.end()
                 );
          }

    /**
     * @brief Clone this model instance.
     * @return A deep copy of this instance.
     */
    amici::Model *clone() const override {
        return new Model_Schwen_PONE2014(*this);
    }

    /**
     * @brief model specific implementation of fJrz
     * @param nllh regularization for event measurements z
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param z model event output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     */
    void fJrz(realtype *nllh, const int iz, const realtype *p,
              const realtype *k, const realtype *rz,
              const realtype *sigmaz) override {}

    void fJy(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my) override {
        Jy_Schwen_PONE2014(Jy, iy, p, k, y, sigmay, my);
    }


    /**
     * @brief model specific implementation of fJz
     * @param nllh negative log-likelihood for event measurements z
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param z model event output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     * @param mz event measurements at timepoint
     */
    void fJz(realtype *nllh, const int iz, const realtype *p,
             const realtype *k, const realtype *z,
             const realtype *sigmaz, const realtype *mz) override {}

    /**
     * @brief model specific implementation of fdJrzdsigma
     * @param dJrzdsigma Sensitivity of event penalization Jrz w.r.t.
     * standard deviation sigmaz
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param rz model root output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     */
    void fdJrzdsigma(realtype *dJrzdsigma, const int iz,
                     const realtype *p, const realtype *k,
                     const realtype *rz,
                     const realtype *sigmaz) override {}

    /**
     * @brief model specific implementation of fdJrzdz
     * @param dJrzdz partial derivative of event penalization Jrz
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param rz model root output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     */
    void fdJrzdz(realtype *dJrzdz, const int iz, const realtype *p,
                 const realtype *k, const realtype *rz,
                 const realtype *sigmaz) override {}

    void fdJydsigma(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my) override {
        dJydsigma_Schwen_PONE2014(dJydsigma, iy, p, k, y, sigmay, my);
    }


    /**
     * @brief model specific implementation of fdJzdsigma
     * @param dJzdsigma Sensitivity of event measurement
     * negative log-likelihood Jz w.r.t. standard deviation sigmaz
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param z model event output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     * @param mz event measurement at timepoint
     */
    void fdJzdsigma(realtype *dJzdsigma, const int iz,
                    const realtype *p, const realtype *k,
                    const realtype *z, const realtype *sigmaz,
                    const realtype *mz) override {}

    /**
     * @brief model specific implementation of fdJzdz
     * @param dJzdz partial derivative of event measurement negative
     * log-likelihood Jz
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param z model event output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     * @param mz event measurement at timepoint
     */
    void fdJzdz(realtype *dJzdz, const int iz, const realtype *p,
                const realtype *k, const realtype *z,
                const realtype *sigmaz, const realtype *mz) override {}

    /**
     * @brief model specific implementation of fdeltasx
     * @param deltaqB sensitivity update
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     * @param ip sensitivity index
     * @param ie event index
     * @param xdot new model right hand side
     * @param xdot_old previous model right hand side
     * @param xB adjoint state
     */
    void fdeltaqB(realtype *deltaqB, const realtype t,
                  const realtype *x, const realtype *p,
                  const realtype *k, const realtype *h, const int ip,
                  const int ie, const realtype *xdot,
                  const realtype *xdot_old,
                  const realtype *xB) override {}

    void fdeltasx(realtype *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *sx, const realtype *stau, const realtype *tcl) override {}


    void fdeltax(double *deltax, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_old) override {}


    /**
     * @brief model specific implementation of fdeltaxB
     * @param deltaxB adjoint state update
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     * @param ie event index
     * @param xdot new model right hand side
     * @param xdot_old previous model right hand side
     * @param xB current adjoint state
     */
    void fdeltaxB(realtype *deltaxB, const realtype t,
                  const realtype *x, const realtype *p,
                  const realtype *k, const realtype *h, const int ie,
                  const realtype *xdot, const realtype *xdot_old,
                  const realtype *xB) override {}

    /**
     * @brief model specific implementation of fdrzdp
     * @param drzdp partial derivative of root output rz w.r.t. model parameters
     * p
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     * @param ip parameter index w.r.t. which the derivative is requested
     */
    void fdrzdp(realtype *drzdp, const int ie, const realtype t,
                const realtype *x, const realtype *p, const realtype *k,
                const realtype *h, const int ip) override {}

    /**
     * @brief model specific implementation of fdrzdx
     * @param drzdx partial derivative of root output rz w.r.t. model states x
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     */
    void fdrzdx(realtype *drzdx, const int ie, const realtype t,
                const realtype *x, const realtype *p, const realtype *k,
                const realtype *h) override {}

    void fdsigmaydp(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip) override {
        dsigmaydp_Schwen_PONE2014(dsigmaydp, t, p, k, y, ip);
    }


    void fdsigmaydy(realtype *dsigmaydy, const realtype t, const realtype *p, const realtype *k, const realtype *y) override {}


    /**
     * @brief model specific implementation of fsigmaz
     * @param dsigmazdp partial derivative of standard deviation of event
     * measurements
     * @param t current time
     * @param p parameter vector
     * @param k constant vector
     * @param ip sensitivity index
     */
    void fdsigmazdp(realtype *dsigmazdp, const realtype t,
                    const realtype *p, const realtype *k,
                    const int ip) override {}

    void fdJydy(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my) override {
        dJydy_Schwen_PONE2014(dJydy, iy, p, k, y, sigmay, my);
    }

    void fdJydy_colptrs(SUNMatrixWrapper &colptrs, int index) override {        dJydy_colptrs_Schwen_PONE2014(colptrs, index);
    }

    void fdJydy_rowvals(SUNMatrixWrapper &rowvals, int index) override {        dJydy_rowvals_Schwen_PONE2014(rowvals, index);
    }


    void fdwdp(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp) override {
        dwdp_Schwen_PONE2014(dwdp, t, x, p, k, h, w, tcl, dtcldp);
    }

    void fdwdp_colptrs(SUNMatrixWrapper &colptrs) override {        dwdp_colptrs_Schwen_PONE2014(colptrs);
    }

    void fdwdp_rowvals(SUNMatrixWrapper &rowvals) override {        dwdp_rowvals_Schwen_PONE2014(rowvals);
    }


    void fdwdx(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl) override {
        dwdx_Schwen_PONE2014(dwdx, t, x, p, k, h, w, tcl);
    }

    void fdwdx_colptrs(SUNMatrixWrapper &colptrs) override {        dwdx_colptrs_Schwen_PONE2014(colptrs);
    }

    void fdwdx_rowvals(SUNMatrixWrapper &rowvals) override {        dwdx_rowvals_Schwen_PONE2014(rowvals);
    }


    void fdwdw(realtype *dwdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl) override {}

    void fdwdw_colptrs(SUNMatrixWrapper &colptrs) override {}

    void fdwdw_rowvals(SUNMatrixWrapper &rowvals) override {}


    void fdxdotdw(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        dxdotdw_Schwen_PONE2014(dxdotdw, t, x, p, k, h, w);
    }

    void fdxdotdw_colptrs(SUNMatrixWrapper &colptrs) override {        dxdotdw_colptrs_Schwen_PONE2014(colptrs);
    }

    void fdxdotdw_rowvals(SUNMatrixWrapper &rowvals) override {        dxdotdw_rowvals_Schwen_PONE2014(rowvals);
    }


    void fdxdotdp_explicit(realtype *dxdotdp_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {}

    void fdxdotdp_explicit_colptrs(SUNMatrixWrapper &colptrs) override {}

    void fdxdotdp_explicit_rowvals(SUNMatrixWrapper &rowvals) override {}


    void fdxdotdx_explicit(realtype *dxdotdx_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {}

    void fdxdotdx_explicit_colptrs(SUNMatrixWrapper &colptrs) override {}

    void fdxdotdx_explicit_rowvals(SUNMatrixWrapper &rowvals) override {}


    void fdydx(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) override {
        dydx_Schwen_PONE2014(dydx, t, x, p, k, h, w, dwdx);
    }


    void fdydp(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *tcl, const realtype *dtcldp) override {
        dydp_Schwen_PONE2014(dydp, t, x, p, k, h, ip, w, tcl, dtcldp);
    }


    /**
     * @brief model specific implementation of fdzdp
     * @param dzdp partial derivative of event-resolved output z w.r.t. model
     * parameters p
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     * @param ip parameter index w.r.t. which the derivative is requested
     */
    void fdzdp(realtype *dzdp, const int ie, const realtype t,
               const realtype *x, const realtype *p, const realtype *k,
               const realtype *h, const int ip) override {}

    /**
     * @brief model specific implementation of fdzdx
     * @param dzdx partial derivative of event-resolved output z w.r.t. model
     * states x
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     */
    void fdzdx(realtype *dzdx, const int ie, const realtype t,
               const realtype *x, const realtype *p, const realtype *k,
               const realtype *h) override {}

    void froot(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl) override {}


    /**
     * @brief model specific implementation of frz
     * @param rz value of root function at current timepoint (non-output events
     * not included)
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     */
    void frz(realtype *rz, const int ie, const realtype t,
             const realtype *x, const realtype *p, const realtype *k,
             const realtype *h) override {}

    void fsigmay(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y) override {
        sigmay_Schwen_PONE2014(sigmay, t, p, k, y);
    }


    /**
     * @brief model specific implementation of fsigmaz
     * @param sigmaz standard deviation of event measurements
     * @param t current time
     * @param p parameter vector
     * @param k constant vector
     */
    void fsigmaz(realtype *sigmaz, const realtype t, const realtype *p,
                 const realtype *k) override {}

    /**
     * @brief model specific implementation of fsrz
     * @param srz Sensitivity of rz, total derivative
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param sx current state sensitivity
     * @param h heaviside vector
     * @param ip sensitivity index
     */
    void fsrz(realtype *srz, const int ie, const realtype t,
              const realtype *x, const realtype *p, const realtype *k,
              const realtype *h, const realtype *sx,
              const int ip) override {}

    void fstau(realtype *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *sx, const int ip, const int ie) override {}

    void fsx0(realtype *sx0, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip) override {
        sx0_Schwen_PONE2014(sx0, t, x, p, k, ip);
    }

    void fsx0_fixedParameters(realtype *sx0_fixedParameters, const realtype t, const realtype *x0, const realtype *p, const realtype *k, const int ip, gsl::span<const int> reinitialization_state_idxs) override {
        sx0_fixedParameters_Schwen_PONE2014(sx0_fixedParameters, t, x0, p, k, ip,  reinitialization_state_idxs);
    }


    /**
     * @brief model specific implementation of fsz
     * @param sz Sensitivity of rz, total derivative
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     * @param sx current state sensitivity
     * @param ip sensitivity index
     */
    void fsz(realtype *sz, const int ie, const realtype t,
             const realtype *x, const realtype *p, const realtype *k,
             const realtype *h, const realtype *sx,
             const int ip) override {}

    void fw(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl) override {
        w_Schwen_PONE2014(w, t, x, p, k, h, tcl);
    }


    void fx0(realtype *x0, const realtype t, const realtype *p, const realtype *k) override {
        x0_Schwen_PONE2014(x0, t, p, k);
    }


    void fx0_fixedParameters(realtype *x0_fixedParameters, const realtype t, const realtype *p, const realtype *k, gsl::span<const int> reinitialization_state_idxs) override {
        x0_fixedParameters_Schwen_PONE2014(x0_fixedParameters, t, p, k,  reinitialization_state_idxs);
    }


    void fxdot(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        xdot_Schwen_PONE2014(xdot, t, x, p, k, h, w);
    }


    void fy(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        y_Schwen_PONE2014(y, t, x, p, k, h, w);
    }


    /**
     * @brief model specific implementation of fz
     * @param z value of event output
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     */
    void fz(realtype *z, const int ie, const realtype t,
            const realtype *x, const realtype *p, const realtype *k,
            const realtype *h) override {}

    

    void fx_solver(realtype *x_solver, const realtype *x_rdata) override {
        x_solver_Schwen_PONE2014(x_solver, x_rdata);
    }


    void ftotal_cl(realtype *total_cl, const realtype *x_rdata, const realtype *p, const realtype *k) override {}


    void fdx_rdatadx_solver(realtype *dx_rdatadx_solver, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k) override {}

    void fdx_rdatadx_solver_colptrs(SUNMatrixWrapper &colptrs) override {}

    void fdx_rdatadx_solver_rowvals(SUNMatrixWrapper &rowvals) override {}


    void fdx_rdatadp(realtype *dx_rdatadp, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k, const int ip) override {}


    void fdx_rdatadtcl(realtype *dx_rdatadtcl, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k) override {}

    void fdx_rdatadtcl_colptrs(SUNMatrixWrapper &colptrs) override {}

    void fdx_rdatadtcl_rowvals(SUNMatrixWrapper &rowvals) override {}


    void fdtotal_cldp(realtype *dtotal_cldp, const realtype *x_rdata, const realtype *p, const realtype *k, const int ip) override {}


    void fdtotal_cldx_rdata(realtype *dtotal_cldx_rdata, const realtype *x_rdata, const realtype *p, const realtype *k, const realtype *tcl) override {}

    void fdtotal_cldx_rdata_colptrs(SUNMatrixWrapper &colptrs) override {}

    void fdtotal_cldx_rdata_rowvals(SUNMatrixWrapper &rowvals) override {}


    std::string getName() const override {
        return "Schwen_PONE2014";
    }

    /**
     * @brief Get names of the model parameters
     * @return the names
     */
    std::vector<std::string> getParameterNames() const override {
        return std::vector<std::string>(parameterNames.begin(),
                                        parameterNames.end());
    }

    /**
     * @brief Get names of the model states
     * @return the names
     */
    std::vector<std::string> getStateNames() const override {
        return std::vector<std::string>(stateNames.begin(), stateNames.end());
    }

    /**
     * @brief Get names of the solver states
     * @return the names
     */
    std::vector<std::string> getStateNamesSolver() const override {
        std::vector<std::string> result;
        result.reserve(stateIdxsSolver.size());
        for(auto idx: stateIdxsSolver) {
            result.push_back(stateNames[idx]);
        }
        return result;
    }

    /**
     * @brief Get names of the fixed model parameters
     * @return the names
     */
    std::vector<std::string> getFixedParameterNames() const override {
        return std::vector<std::string>(fixedParameterNames.begin(),
                                        fixedParameterNames.end());
    }

    /**
     * @brief Get names of the observables
     * @return the names
     */
    std::vector<std::string> getObservableNames() const override {
        return std::vector<std::string>(observableNames.begin(),
                                        observableNames.end());
    }

    /**
     * @brief Get names of model expressions
     * @return Expression names
     */
    std::vector<std::string> getExpressionNames() const override {
        return std::vector<std::string>(expressionNames.begin(),
                                        expressionNames.end());
    }

    /**
     * @brief Get ids of the model parameters
     * @return the ids
     */
    std::vector<std::string> getParameterIds() const override {
        return std::vector<std::string>(parameterIds.begin(),
                                        parameterIds.end());
    }

    /**
     * @brief Get ids of the model states
     * @return the ids
     */
    std::vector<std::string> getStateIds() const override {
        return std::vector<std::string>(stateIds.begin(), stateIds.end());
    }

    /**
     * @brief Get ids of the solver states
     * @return the ids
     */
    std::vector<std::string> getStateIdsSolver() const override {
        std::vector<std::string> result;
        result.reserve(stateIdxsSolver.size());
        for(auto idx: stateIdxsSolver) {
            result.push_back(stateIds[idx]);
        }
        return result;
    }

    /**
     * @brief Get ids of the fixed model parameters
     * @return the ids
     */
    std::vector<std::string> getFixedParameterIds() const override {
        return std::vector<std::string>(fixedParameterIds.begin(),
                                        fixedParameterIds.end());
    }

    /**
     * @brief Get ids of the observables
     * @return the ids
     */
    std::vector<std::string> getObservableIds() const override {
        return std::vector<std::string>(observableIds.begin(),
                                        observableIds.end());
    }

    /**
     * @brief Get IDs of model expressions
     * @return Expression IDs
     */
    std::vector<std::string> getExpressionIds() const override {
        return std::vector<std::string>(expressionIds.begin(),
                                        expressionIds.end());
    }

    /**
     * @brief function indicating whether reinitialization of states depending
     * on fixed parameters is permissible
     * @return flag indicating whether reinitialization of states depending on
     * fixed parameters is permissible
     */
    bool isFixedParameterStateReinitializationAllowed() const override {
        return true;
    }

    /**
     * @brief returns the AMICI version that was used to generate the model
     * @return AMICI version string
     */
    std::string getAmiciVersion() const override {
        return "0.11.32";
    }

    /**
     * @brief returns the amici version that was used to generate the model
     * @return AMICI git commit hash
     */
    std::string getAmiciCommit() const override {
        return "unknown";
    }

    bool hasQuadraticLLH() const override {
        return true;
    }

    ObservableScaling getObservableScaling(int iy) const override {
        return observableScalings.at(iy);
    }
};


} // namespace model_Schwen_PONE2014

} // namespace amici

#endif /* _amici_TPL_MODELNAME_h */
