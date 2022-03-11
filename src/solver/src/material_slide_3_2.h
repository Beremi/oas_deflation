#ifndef _MATERIAL_SLIDE_3_2
#define _MATERIAL_SLIDE_3_2

#include "material.h"

/******************************************************************************
*                       Code generated with sympy 1.4                        *
*                                                                            *
*              See http://www.sympy.org/ for more information.               *
*                                                                            *
*                       This file is part of 'project'                       *
******************************************************************************/

// NOTE JK: for material params, it would be better to acess them directly "const double &var" instead of making a copy of each one (is it just the cpp syntax?)
// NOTE 2 JK - in case of state variables, why copies: I do not want to change them, just evaluate their copy and change thame after the converged step
// TODO out variable is not a double, but array (vector) of doubles, or even a matrix
// TODO does it work with my double[8] arrays, how to send them by pointer?
void get_Sig(double E_s, double E_w, double K_s, double alpha_x, double alpha_y, double gamma_s, double omega_s, double omega_w, double s_pi_x, double s_pi_y, double s_x, double s_y, double w, double w_pi, double z, double *out_2915105096036016098);
void get_dSig_dEps(double E_s, double E_w, double K_s, double gamma_s, double omega_s, double omega_w, double s_pi_x, double s_pi_y, double s_x, double s_y, double w, double w_pi, double *out_2558847501224171655);
double get_f(double X_x, double X_y, double Z, double bartau, double sigma_pi, double tau_pi_x, double tau_pi_y, double f_c, double f_c0, double f_t, double m, double omega_s, double omega_w);
void get_df_dSig(double X_x, double X_y, double Z, double bartau, double sigma_pi, double tau_pi_x, double tau_pi_y, double f_c, double f_c0, double f_t, double m, double omega_s, double omega_w, double *out_6697882517837858595);
void get_df_dEps(double X_x, double X_y, double Z, double bartau, double sigma_pi, double tau_pi_x, double tau_pi_y, double f_c, double f_c0, double f_t, double m, double omega_s, double omega_w, double *out_5723301903165105508);
void get_Phi(double S_s, double S_w, double X_x, double X_y, double Y_s, double Y_w, double Z, double bartau, double sigma_pi, double tau_pi_x, double tau_pi_y, double c_s, double c_w, double eta, double f_c, double f_c0, double f_t, double m, double omega_s, double omega_w, double *out_5145644395765063632);

/******************************************************************************
******************************************************************************/

class Slide32Material;
class Slide32MaterialStatus : public DisMechMaterialStatus
{
private:
    // state variables
    // TODO make this array of state variables (double[8])
    double strain_slip_multiplier;

    double *Eps;
    double *temp_Eps;
    // double s_pi_x, s_pi_y, w_pi, z, alpha_x, alpha_y, omega_s, omega_w;

    // thermodynamic forces
    double *Sig;
    double *temp_Sig;
    // double tau_pi_x, tau_pi_y, sig_pi, Z, X_x, X_y, Y_s, Y_w;

    unsigned max_iter; ///> maximum iteration for return mapping

    // threshold fn:
    void get_f_df(const double &s_x_n1, const double &s_y_n1, const double &w_n1, double &f_k, double &df_k);

    void get_Eps_k1(const double &s_x_n1, const double &s_y_n1, const double &w_n1, const double &lam_k);

    void print_state_vars(const bool &temp = true);
    void print_thrmdyn_forces(const bool &temp = true);

    void check_state_variable_ranges();

public:
    Slide32MaterialStatus(Slide32Material *m, Element *e, unsigned ipnum);
    virtual ~Slide32MaterialStatus() {};
    void init();
    virtual void update();
    virtual MyMatrix giveStiffnessTensor(std :: string type, unsigned dim) const;
    virtual MyVector giveStress(const MyVector &strain, double timeStep);
    virtual void giveValues(std :: string code, MyVector &result) const;
};


class Slide32Material : public DisMechMaterial
{
private:
    // tangential material params:
    // double E_s = E0 * alpha; ///> nonnegative
    double bartau; ///> nonnegative  [Pa]
    double gamma_s; ///> nonnegative  [Pa]
    double K_s; ///>   [Pa]
    double S_s; ///>   [Pa]
    double r_s; ///>   [-]
    double c_s; ///>   [-]
    // normal material params:
    // double E_w = E0; ///> nonnegative
    double S_w; ///> [Pa]
    double r_w; ///>   [-]
    double c_w; ///>   [-]
    double f_t; ///> nonnegative
    double f_c; ///> nonnegative
    double f_c0; ///> nonnegative
    double m; ///> nonnegative compression shear slope [-]

    double eta; ///> nonnegative [-]

    bool use_discrete_values; // if true calcualte in slips (or strains if false)
    unsigned max_iter; // maximum naumber of iteration in return mapping
public:
    Slide32Material() { name = "Slide32 material"; };
    ~Slide32Material() {};
    void readFromLine(std :: istringstream &iss);
    MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);

    double giveE_s() const { return E0 * alpha; }
    double giveGamma_s() const { return gamma_s; }
    double giveK_s() const { return K_s; }
    double giveS_s() const { return S_s; }
    double giveR_s() const { return r_s; }
    double giveC_s() const { return c_s; }
    double giveBartau() const { return bartau; }
    // normal material params:
    double giveE_w() const { return E0; }
    double giveS_w() const { return S_w; }
    double giveR_w() const { return r_w; }
    double giveC_w() const { return c_w; }
    double giveF_t() const { return f_t; }
    double giveF_c() const { return f_c; }
    double giveF_c0() const { return f_c0; }
    double giveM() const { return m; }

    double giveEta() const { return eta; }

    bool calc_slips() const { return use_discrete_values; }
    unsigned giveMax_iter() const { return max_iter; }

    virtual void init() {};
};

#endif /* _MATERIAL_SLIDE_3_2 */
