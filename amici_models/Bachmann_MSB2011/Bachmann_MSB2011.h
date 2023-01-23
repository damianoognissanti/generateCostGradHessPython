#ifndef _amici_TPL_MODELNAME_h
#define _amici_TPL_MODELNAME_h
#include <cmath>
#include <memory>
#include <gsl/gsl-lite.hpp>

#include "amici/model_ode.h"

namespace amici {

class Solver;

namespace model_Bachmann_MSB2011 {

extern std::array<const char*, 69> parameterNames;
extern std::array<const char*, 10> fixedParameterNames;
extern std::array<const char*, 25> stateNames;
extern std::array<const char*, 20> observableNames;
extern std::array<const ObservableScaling, 20> observableScalings;
extern std::array<const char*, 36> expressionNames;
extern std::array<const char*, 69> parameterIds;
extern std::array<const char*, 10> fixedParameterIds;
extern std::array<const char*, 25> stateIds;
extern std::array<const char*, 20> observableIds;
extern std::array<const char*, 36> expressionIds;
extern std::array<int, 25> stateIdxsSolver;
extern std::array<bool, 0> rootInitialValues;

extern void Jy_Bachmann_MSB2011(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydsigma_Bachmann_MSB2011(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydy_Bachmann_MSB2011(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydy_colptrs_Bachmann_MSB2011(SUNMatrixWrapper &colptrs, int index);
extern void dJydy_rowvals_Bachmann_MSB2011(SUNMatrixWrapper &rowvals, int index);







extern void dwdp_Bachmann_MSB2011(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp);
extern void dwdp_colptrs_Bachmann_MSB2011(SUNMatrixWrapper &colptrs);
extern void dwdp_rowvals_Bachmann_MSB2011(SUNMatrixWrapper &rowvals);
extern void dwdx_Bachmann_MSB2011(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl);
extern void dwdx_colptrs_Bachmann_MSB2011(SUNMatrixWrapper &colptrs);
extern void dwdx_rowvals_Bachmann_MSB2011(SUNMatrixWrapper &rowvals);



extern void dxdotdw_Bachmann_MSB2011(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void dxdotdw_colptrs_Bachmann_MSB2011(SUNMatrixWrapper &colptrs);
extern void dxdotdw_rowvals_Bachmann_MSB2011(SUNMatrixWrapper &rowvals);






extern void dydx_Bachmann_MSB2011(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx);
extern void dydp_Bachmann_MSB2011(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *tcl, const realtype *dtcldp);




extern void sigmay_Bachmann_MSB2011(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y);

extern void dsigmaydp_Bachmann_MSB2011(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip);


extern void w_Bachmann_MSB2011(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl);
extern void x0_Bachmann_MSB2011(realtype *x0, const realtype t, const realtype *p, const realtype *k);
extern void x0_fixedParameters_Bachmann_MSB2011(realtype *x0_fixedParameters, const realtype t, const realtype *p, const realtype *k, gsl::span<const int> reinitialization_state_idxs);
extern void sx0_Bachmann_MSB2011(realtype *sx0, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip);
extern void sx0_fixedParameters_Bachmann_MSB2011(realtype *sx0_fixedParameters, const realtype t, const realtype *x0, const realtype *p, const realtype *k, const int ip, gsl::span<const int> reinitialization_state_idxs);
extern void xdot_Bachmann_MSB2011(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void y_Bachmann_MSB2011(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);






extern void x_solver_Bachmann_MSB2011(realtype *x_solver, const realtype *x_rdata);












/**
 * @brief AMICI-generated model subclass.
 */
class Model_Bachmann_MSB2011 : public amici::Model_ODE {
  public:
    /**
     * @brief Default constructor.
     */
    Model_Bachmann_MSB2011()
        : amici::Model_ODE(
              amici::ModelDimensions(
                  25,                            // nx_rdata
                  25,                        // nxtrue_rdata
                  25,                           // nx_solver
                  25,                       // nxtrue_solver
                  0,                    // nx_solver_reinit
                  69,                                  // np
                  10,                                  // nk
                  20,                                  // ny
                  20,                              // nytrue
                  0,                                  // nz
                  0,                              // nztrue
                  0,                              // nevent
                  1,                          // nobjective
                  36,                                  // nw
                  60,                               // ndwdx
                  70,                               // ndwdp
                  0,                               // ndwdw
                  61,                            // ndxdotdw
                  std::vector<int>{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},                              // ndjydy
                  0,                    // ndxrdatadxsolver
                  0,                        // ndxrdatadtcl
                  0,                        // ndtotal_cldx_rdata
                  0,                                       // nnz
                  25,                                 // ubw
                  25                                  // lbw
              ),
              amici::SimulationParameters(
                  std::vector<realtype>{1.0, 1.0, 1.2499999999999999e-7, 1.2499999999999999e-7, 1.2499999999999999e-7, 1.2499999999999999e-7, 1.0, 0.0, 0.0, 0.0}, // fixedParameters
                  std::vector<realtype>{432.86041343491303, 785269991.45049596, 0.14477775532111001, 999.999999999946, 0.0083988695167017003, 0.53026444711960896, 0.26730484933305798, 999999.99999991199, 5.4298069390344796, 633167.43060080602, 142.72332309737999, 0.0010000000000000601, 0.0081622049095037407, 2.8256815341155499, 0.679165515556864, 173.64470023135999, 10.4078649133666, 1.06458446742251, 0.0083091764312036905, 9999.9999999991196, 38.995799107394802, 0.078106885579546703, 0.074515081901642299, 0.026886508382968501, 3.97622369384192, 26.725116427710901, 79.753639937710005, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}        // dynamic parameters
              ),
              amici::SecondOrderMode::none,                                  // o2mode
              std::vector<realtype>(25, 0.0),   // idlist
              std::vector<int>{},               // z2events
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
        return new Model_Bachmann_MSB2011(*this);
    }

    void fJrz(realtype *Jrz, const int iz, const realtype *p, const realtype *k, const realtype *rz, const realtype *sigmaz) override {}


    void fJy(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my) override {
        Jy_Bachmann_MSB2011(Jy, iy, p, k, y, sigmay, my);    
}


    void fJz(realtype *Jz, const int iz, const realtype *p, const realtype *k, const realtype *z, const realtype *sigmaz, const realtype *mz) override {}


    void fdJrzdsigma(realtype *dJrzdsigma, const int iz, const realtype *p, const realtype *k, const realtype *rz, const realtype *sigmaz) override {}


    void fdJrzdz(realtype *dJrzdz, const int iz, const realtype *p, const realtype *k, const realtype *rz, const realtype *sigmaz) override {}


    void fdJydsigma(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my) override {
        dJydsigma_Bachmann_MSB2011(dJydsigma, iy, p, k, y, sigmay, my);    
}


    void fdJzdsigma(realtype *dJzdsigma, const int iz, const realtype *p, const realtype *k, const realtype *z, const realtype *sigmaz, const realtype *mz) override {}


    void fdJzdz(realtype *dJzdz, const int iz, const realtype *p, const realtype *k, const realtype *z, const realtype *sigmaz, const double *mz) override {}


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

    void fdrzdp(realtype *drzdp, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip) override {}


    void fdrzdx(realtype *drzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {}


    void fdsigmaydp(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip) override {
        dsigmaydp_Bachmann_MSB2011(dsigmaydp, t, p, k, y, ip);    
}


    void fdsigmaydy(realtype *dsigmaydy, const realtype t, const realtype *p, const realtype *k, const realtype *y) override {}


    void fdsigmazdp(realtype *dsigmazdp, const realtype t, const realtype *p, const realtype *k, const int ip) override {}


    void fdJydy(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my) override {
        dJydy_Bachmann_MSB2011(dJydy, iy, p, k, y, sigmay, my);    
}

    void fdJydy_colptrs(SUNMatrixWrapper &colptrs, int index) override {        dJydy_colptrs_Bachmann_MSB2011(colptrs, index);
    }

    void fdJydy_rowvals(SUNMatrixWrapper &rowvals, int index) override {        dJydy_rowvals_Bachmann_MSB2011(rowvals, index);
    }


    void fdwdp(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp) override {
        dwdp_Bachmann_MSB2011(dwdp, t, x, p, k, h, w, tcl, dtcldp);    
}

    void fdwdp_colptrs(SUNMatrixWrapper &colptrs) override {        dwdp_colptrs_Bachmann_MSB2011(colptrs);
    }

    void fdwdp_rowvals(SUNMatrixWrapper &rowvals) override {        dwdp_rowvals_Bachmann_MSB2011(rowvals);
    }


    void fdwdx(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl) override {
        dwdx_Bachmann_MSB2011(dwdx, t, x, p, k, h, w, tcl);    
}

    void fdwdx_colptrs(SUNMatrixWrapper &colptrs) override {        dwdx_colptrs_Bachmann_MSB2011(colptrs);
    }

    void fdwdx_rowvals(SUNMatrixWrapper &rowvals) override {        dwdx_rowvals_Bachmann_MSB2011(rowvals);
    }


    void fdwdw(realtype *dwdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl) override {}

    void fdwdw_colptrs(SUNMatrixWrapper &colptrs) override {}

    void fdwdw_rowvals(SUNMatrixWrapper &rowvals) override {}


    void fdxdotdw(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        dxdotdw_Bachmann_MSB2011(dxdotdw, t, x, p, k, h, w);    
}

    void fdxdotdw_colptrs(SUNMatrixWrapper &colptrs) override {        dxdotdw_colptrs_Bachmann_MSB2011(colptrs);
    }

    void fdxdotdw_rowvals(SUNMatrixWrapper &rowvals) override {        dxdotdw_rowvals_Bachmann_MSB2011(rowvals);
    }


    void fdxdotdp_explicit(realtype *dxdotdp_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {}

    void fdxdotdp_explicit_colptrs(SUNMatrixWrapper &colptrs) override {}

    void fdxdotdp_explicit_rowvals(SUNMatrixWrapper &rowvals) override {}


    void fdxdotdx_explicit(realtype *dxdotdx_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {}

    void fdxdotdx_explicit_colptrs(SUNMatrixWrapper &colptrs) override {}

    void fdxdotdx_explicit_rowvals(SUNMatrixWrapper &rowvals) override {}


    void fdydx(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) override {
        dydx_Bachmann_MSB2011(dydx, t, x, p, k, h, w, dwdx);    
}


    void fdydp(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *tcl, const realtype *dtcldp) override {
        dydp_Bachmann_MSB2011(dydp, t, x, p, k, h, ip, w, tcl, dtcldp);    
}


    void fdzdp(realtype *dzdp, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip) override {}


    void fdzdx(realtype *dzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {}


    void froot(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl) override {}


    void frz(realtype *rz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {}


    void fsigmay(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y) override {
        sigmay_Bachmann_MSB2011(sigmay, t, p, k, y);    
}


    void fsigmaz(realtype *sigmaz, const realtype t, const realtype *p, const realtype *k) override {}


    void fstau(realtype *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *sx, const int ip, const int ie) override {}

    void fsx0(realtype *sx0, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip) override {
        sx0_Bachmann_MSB2011(sx0, t, x, p, k, ip);    
}

    void fsx0_fixedParameters(realtype *sx0_fixedParameters, const realtype t, const realtype *x0, const realtype *p, const realtype *k, const int ip, gsl::span<const int> reinitialization_state_idxs) override {
        sx0_fixedParameters_Bachmann_MSB2011(sx0_fixedParameters, t, x0, p, k, ip,  reinitialization_state_idxs);    
}


    void fw(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl) override {
        w_Bachmann_MSB2011(w, t, x, p, k, h, tcl);    
}


    void fx0(realtype *x0, const realtype t, const realtype *p, const realtype *k) override {
        x0_Bachmann_MSB2011(x0, t, p, k);    
}


    void fx0_fixedParameters(realtype *x0_fixedParameters, const realtype t, const realtype *p, const realtype *k, gsl::span<const int> reinitialization_state_idxs) override {
        x0_fixedParameters_Bachmann_MSB2011(x0_fixedParameters, t, p, k,  reinitialization_state_idxs);    
}


    void fxdot(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        xdot_Bachmann_MSB2011(xdot, t, x, p, k, h, w);    
}


    void fy(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        y_Bachmann_MSB2011(y, t, x, p, k, h, w);    
}


    void fz(realtype *z, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {}


    

    void fx_solver(realtype *x_solver, const realtype *x_rdata) override {
        x_solver_Bachmann_MSB2011(x_solver, x_rdata);    
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
        return "Bachmann_MSB2011";
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
        return "0.15.0";
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


} // namespace model_Bachmann_MSB2011

} // namespace amici

#endif /* _amici_TPL_MODELNAME_h */
