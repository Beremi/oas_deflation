#include "material_slide_3_2.h"
#include "element.h"
#include <algorithm>

/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                       This file is part of 'project'                       *
 ******************************************************************************/


void get_Sig(double E_s, double E_w, double K_s, double alpha_x, double alpha_y, double gamma_s, double omega_s, double omega_w, double s_pi_x, double s_pi_y, double s_x, double s_y, double w, double w_pi, double z, double *out_2915105096036016098) {
   out_2915105096036016098[0] = -1.0/2.0*E_s*(1 - omega_s)*(2*s_pi_x - 2*s_x);
   out_2915105096036016098[1] = -1.0/2.0*E_s*(1 - omega_s)*(2*s_pi_y - 2*s_y);
   out_2915105096036016098[2] = -E_w*(1.0/2.0 - 1.0/2.0*omega_w)*(-2*w + 2*w_pi);
   out_2915105096036016098[3] = K_s*z;
   out_2915105096036016098[4] = alpha_x*gamma_s;
   out_2915105096036016098[5] = alpha_y*gamma_s;
   out_2915105096036016098[6] = (1.0/2.0)*E_s*pow(-s_pi_x + s_x, 2) + (1.0/2.0)*E_s*pow(-s_pi_y + s_y, 2);
   out_2915105096036016098[7] = (1.0/2.0)*E_w*pow(w - w_pi, 2);
}
void get_dSig_dEps(double E_s, double E_w, double K_s, double gamma_s, double omega_s, double omega_w, double s_pi_x, double s_pi_y, double s_x, double s_y, double w, double w_pi, double *out_2558847501224171655) {
   out_2558847501224171655[0] = -E_s*(1 - omega_s);
   out_2558847501224171655[1] = 0;
   out_2558847501224171655[2] = 0;
   out_2558847501224171655[3] = 0;
   out_2558847501224171655[4] = 0;
   out_2558847501224171655[5] = 0;
   out_2558847501224171655[6] = (1.0/2.0)*E_s*(2*s_pi_x - 2*s_x);
   out_2558847501224171655[7] = 0;
   out_2558847501224171655[8] = 0;
   out_2558847501224171655[9] = -E_s*(1 - omega_s);
   out_2558847501224171655[10] = 0;
   out_2558847501224171655[11] = 0;
   out_2558847501224171655[12] = 0;
   out_2558847501224171655[13] = 0;
   out_2558847501224171655[14] = (1.0/2.0)*E_s*(2*s_pi_y - 2*s_y);
   out_2558847501224171655[15] = 0;
   out_2558847501224171655[16] = 0;
   out_2558847501224171655[17] = 0;
   out_2558847501224171655[18] = -2*E_w*(1.0/2.0 - 1.0/2.0*omega_w);
   out_2558847501224171655[19] = 0;
   out_2558847501224171655[20] = 0;
   out_2558847501224171655[21] = 0;
   out_2558847501224171655[22] = 0;
   out_2558847501224171655[23] = (1.0/2.0)*E_w*(-2*w + 2*w_pi);
   out_2558847501224171655[24] = 0;
   out_2558847501224171655[25] = 0;
   out_2558847501224171655[26] = 0;
   out_2558847501224171655[27] = K_s;
   out_2558847501224171655[28] = 0;
   out_2558847501224171655[29] = 0;
   out_2558847501224171655[30] = 0;
   out_2558847501224171655[31] = 0;
   out_2558847501224171655[32] = 0;
   out_2558847501224171655[33] = 0;
   out_2558847501224171655[34] = 0;
   out_2558847501224171655[35] = 0;
   out_2558847501224171655[36] = gamma_s;
   out_2558847501224171655[37] = 0;
   out_2558847501224171655[38] = 0;
   out_2558847501224171655[39] = 0;
   out_2558847501224171655[40] = 0;
   out_2558847501224171655[41] = 0;
   out_2558847501224171655[42] = 0;
   out_2558847501224171655[43] = 0;
   out_2558847501224171655[44] = 0;
   out_2558847501224171655[45] = gamma_s;
   out_2558847501224171655[46] = 0;
   out_2558847501224171655[47] = 0;
   out_2558847501224171655[48] = (1.0/2.0)*E_s*(2*s_pi_x - 2*s_x);
   out_2558847501224171655[49] = (1.0/2.0)*E_s*(2*s_pi_y - 2*s_y);
   out_2558847501224171655[50] = 0;
   out_2558847501224171655[51] = 0;
   out_2558847501224171655[52] = 0;
   out_2558847501224171655[53] = 0;
   out_2558847501224171655[54] = 0;
   out_2558847501224171655[55] = 0;
   out_2558847501224171655[56] = 0;
   out_2558847501224171655[57] = 0;
   out_2558847501224171655[58] = (1.0/2.0)*E_w*(-2*w + 2*w_pi);
   out_2558847501224171655[59] = 0;
   out_2558847501224171655[60] = 0;
   out_2558847501224171655[61] = 0;
   out_2558847501224171655[62] = 0;
   out_2558847501224171655[63] = 0;
}
double get_f(double X_x, double X_y, double Z, double bartau, double sigma_pi, double tau_pi_x, double tau_pi_y, double f_c, double f_c0, double f_t, double m, double omega_s, double omega_w) {
   double get_f_result;
   if ((sqrt(pow(X_x + tau_pi_x/(omega_s - 1), 2) + pow(X_y + tau_pi_y/(omega_s - 1), 2)) - (Z + bartau)*(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m))*(Z + bartau - 2*f_t*m)/(pow(f_t, 2)*m))*(((f_t) > 0) - ((f_t) < 0))*(((m) > 0) - ((m) < 0)) < 0) {
      get_f_result = sqrt((pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/((Z + bartau)*pow(-Z - bartau + f_t*m, 2)*(Z + bartau - 2*f_t*m)) + pow(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m), 2)*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/(pow(f_t, 2)*pow(Z + bartau - f_t*m, 2))) - (pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2))/(Z + bartau - 2*f_t*m);
   }
   else if ((sqrt(pow(X_x + tau_pi_x/(omega_s - 1), 2) + pow(X_y + tau_pi_y/(omega_s - 1), 2)) - (Z + bartau + f_c0*m)*(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m))*(Z + bartau + 2*f_c*m - f_c0*m)/(m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))))*(((m) > 0) - ((m) < 0))*(((-f_c + f_c0) > 0) - ((-f_c + f_c0) < 0)) < 0) {
      get_f_result = sqrt((pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(-Z - bartau - f_c*m, 2)*(Z + bartau + f_c0*m)*(Z + bartau + 2*f_c*m - f_c0*m)) + pow(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m), 2)*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(f_c - f_c0, 2)*pow(Z + bartau + f_c*m, 2))) - (pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2))/(Z + bartau + 2*f_c*m - f_c0*m);
   }
   else {
      get_f_result = -Z - bartau + sigma_pi*m/(1 - omega_w) + sqrt(pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2));
   }
   return get_f_result;
}
void get_df_dSig(double X_x, double X_y, double Z, double bartau, double sigma_pi, double tau_pi_x, double tau_pi_y, double f_c, double f_c0, double f_t, double m, double omega_s, double omega_w, double *out_6697882517837858595) {
   if ((sqrt(pow(X_x + tau_pi_x/(omega_s - 1), 2) + pow(X_y + tau_pi_y/(omega_s - 1), 2)) - (Z + bartau)*(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m))*(Z + bartau - 2*f_t*m)/(pow(f_t, 2)*m))*(((f_t) > 0) - ((f_t) < 0))*(((m) > 0) - ((m) < 0)) < 0) {
      out_6697882517837858595[0] = (-X_x + tau_pi_x/(1 - omega_s))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/((1 - omega_s)*(Z + bartau)*sqrt((pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/((Z + bartau)*pow(-Z - bartau + f_t*m, 2)*(Z + bartau - 2*f_t*m)) + pow(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m), 2)*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/(pow(f_t, 2)*pow(Z + bartau - f_t*m, 2)))*pow(-Z - bartau + f_t*m, 2)*(Z + bartau - 2*f_t*m));
   }
   else if ((sqrt(pow(X_x + tau_pi_x/(omega_s - 1), 2) + pow(X_y + tau_pi_y/(omega_s - 1), 2)) - (Z + bartau + f_c0*m)*(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m))*(Z + bartau + 2*f_c*m - f_c0*m)/(m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))))*(((m) > 0) - ((m) < 0))*(((-f_c + f_c0) > 0) - ((-f_c + f_c0) < 0)) < 0) {
      out_6697882517837858595[0] = (-X_x + tau_pi_x/(1 - omega_s))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/((1 - omega_s)*sqrt((pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(-Z - bartau - f_c*m, 2)*(Z + bartau + f_c0*m)*(Z + bartau + 2*f_c*m - f_c0*m)) + pow(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m), 2)*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(f_c - f_c0, 2)*pow(Z + bartau + f_c*m, 2)))*pow(-Z - bartau - f_c*m, 2)*(Z + bartau + f_c0*m)*(Z + bartau + 2*f_c*m - f_c0*m));
   }
   else {
      out_6697882517837858595[0] = (-X_x + tau_pi_x/(1 - omega_s))/((1 - omega_s)*sqrt(pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2)));
   }
   if ((sqrt(pow(X_x + tau_pi_x/(omega_s - 1), 2) + pow(X_y + tau_pi_y/(omega_s - 1), 2)) - (Z + bartau)*(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m))*(Z + bartau - 2*f_t*m)/(pow(f_t, 2)*m))*(((f_t) > 0) - ((f_t) < 0))*(((m) > 0) - ((m) < 0)) < 0) {
      out_6697882517837858595[1] = (-X_y + tau_pi_y/(1 - omega_s))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/((1 - omega_s)*(Z + bartau)*sqrt((pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/((Z + bartau)*pow(-Z - bartau + f_t*m, 2)*(Z + bartau - 2*f_t*m)) + pow(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m), 2)*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/(pow(f_t, 2)*pow(Z + bartau - f_t*m, 2)))*pow(-Z - bartau + f_t*m, 2)*(Z + bartau - 2*f_t*m));
   }
   else if ((sqrt(pow(X_x + tau_pi_x/(omega_s - 1), 2) + pow(X_y + tau_pi_y/(omega_s - 1), 2)) - (Z + bartau + f_c0*m)*(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m))*(Z + bartau + 2*f_c*m - f_c0*m)/(m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))))*(((m) > 0) - ((m) < 0))*(((-f_c + f_c0) > 0) - ((-f_c + f_c0) < 0)) < 0) {
      out_6697882517837858595[1] = (-X_y + tau_pi_y/(1 - omega_s))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/((1 - omega_s)*sqrt((pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(-Z - bartau - f_c*m, 2)*(Z + bartau + f_c0*m)*(Z + bartau + 2*f_c*m - f_c0*m)) + pow(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m), 2)*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(f_c - f_c0, 2)*pow(Z + bartau + f_c*m, 2)))*pow(-Z - bartau - f_c*m, 2)*(Z + bartau + f_c0*m)*(Z + bartau + 2*f_c*m - f_c0*m));
   }
   else {
      out_6697882517837858595[1] = (-X_y + tau_pi_y/(1 - omega_s))/((1 - omega_s)*sqrt(pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2)));
   }
   if ((sqrt(pow(X_x + tau_pi_x/(omega_s - 1), 2) + pow(X_y + tau_pi_y/(omega_s - 1), 2)) - (Z + bartau)*(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m))*(Z + bartau - 2*f_t*m)/(pow(f_t, 2)*m))*(((f_t) > 0) - ((f_t) < 0))*(((m) > 0) - ((m) < 0)) < 0) {
      out_6697882517837858595[2] = (sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/(pow(f_t, 2)*(1 - omega_w)*sqrt((pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/((Z + bartau)*pow(-Z - bartau + f_t*m, 2)*(Z + bartau - 2*f_t*m)) + pow(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m), 2)*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/(pow(f_t, 2)*pow(Z + bartau - f_t*m, 2)))*pow(Z + bartau - f_t*m, 2));
   }
   else if ((sqrt(pow(X_x + tau_pi_x/(omega_s - 1), 2) + pow(X_y + tau_pi_y/(omega_s - 1), 2)) - (Z + bartau + f_c0*m)*(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m))*(Z + bartau + 2*f_c*m - f_c0*m)/(m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))))*(((m) > 0) - ((m) < 0))*(((-f_c + f_c0) > 0) - ((-f_c + f_c0) < 0)) < 0) {
      out_6697882517837858595[2] = (sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/((1 - omega_w)*pow(f_c - f_c0, 2)*sqrt((pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(-Z - bartau - f_c*m, 2)*(Z + bartau + f_c0*m)*(Z + bartau + 2*f_c*m - f_c0*m)) + pow(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m), 2)*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(f_c - f_c0, 2)*pow(Z + bartau + f_c*m, 2)))*pow(Z + bartau + f_c*m, 2));
   }
   else {
      out_6697882517837858595[2] = m/(1 - omega_w);
   }
   if ((sqrt(pow(X_x + tau_pi_x/(omega_s - 1), 2) + pow(X_y + tau_pi_y/(omega_s - 1), 2)) - (Z + bartau)*(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m))*(Z + bartau - 2*f_t*m)/(pow(f_t, 2)*m))*(((f_t) > 0) - ((f_t) < 0))*(((m) > 0) - ((m) < 0)) < 0) {
      out_6697882517837858595[3] = -(2*Z + 2*bartau - 2*f_t*m)/(Z + bartau - 2*f_t*m) + (pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2))/pow(Z + bartau - 2*f_t*m, 2) + (-m*(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/(pow(Z + bartau - 2*f_t*m, 2)*pow(Z + bartau - f_t*m, 2)) + (1.0/2.0)*(pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*(4*Z + 4*bartau - 4*f_t*m)*(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2))/((Z + bartau)*pow(-Z - bartau + f_t*m, 2)*(Z + bartau - 2*f_t*m)) - 1.0/2.0*(pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/((Z + bartau)*pow(-Z - bartau + f_t*m, 2)*pow(Z + bartau - 2*f_t*m, 2)) + (pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/((Z + bartau)*pow(-Z - bartau + f_t*m, 3)*(Z + bartau - 2*f_t*m)) - 1.0/2.0*(pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/(pow(Z + bartau, 2)*pow(-Z - bartau + f_t*m, 2)*(Z + bartau - 2*f_t*m)) + (1.0/2.0)*pow(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m), 2)*(4*Z + 4*bartau - 4*f_t*m)*(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2))/(pow(f_t, 2)*pow(Z + bartau - f_t*m, 2)) - pow(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m), 2)*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/(pow(f_t, 2)*pow(Z + bartau - f_t*m, 3)))/sqrt((pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/((Z + bartau)*pow(-Z - bartau + f_t*m, 2)*(Z + bartau - 2*f_t*m)) + pow(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m), 2)*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/(pow(f_t, 2)*pow(Z + bartau - f_t*m, 2)));
   }
   else if ((sqrt(pow(X_x + tau_pi_x/(omega_s - 1), 2) + pow(X_y + tau_pi_y/(omega_s - 1), 2)) - (Z + bartau + f_c0*m)*(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m))*(Z + bartau + 2*f_c*m - f_c0*m)/(m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))))*(((m) > 0) - ((m) < 0))*(((-f_c + f_c0) > 0) - ((-f_c + f_c0) < 0)) < 0) {
      out_6697882517837858595[3] = -(2*Z + 2*bartau + 2*f_c*m)/(Z + bartau + 2*f_c*m - f_c0*m) + (pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2))/pow(Z + bartau + 2*f_c*m - f_c0*m, 2) + (-m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))*(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(f_c - f_c0, 2)*pow(Z + bartau + f_c*m, 2)*pow(Z + bartau + 2*f_c*m - f_c0*m, 2)) + (1.0/2.0)*(pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*(4*Z + 4*bartau + 4*f_c*m)*(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2))/(pow(-Z - bartau - f_c*m, 2)*(Z + bartau + f_c0*m)*(Z + bartau + 2*f_c*m - f_c0*m)) - 1.0/2.0*(pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(-Z - bartau - f_c*m, 2)*(Z + bartau + f_c0*m)*pow(Z + bartau + 2*f_c*m - f_c0*m, 2)) - 1.0/2.0*(pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(-Z - bartau - f_c*m, 2)*pow(Z + bartau + f_c0*m, 2)*(Z + bartau + 2*f_c*m - f_c0*m)) + (pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(-Z - bartau - f_c*m, 3)*(Z + bartau + f_c0*m)*(Z + bartau + 2*f_c*m - f_c0*m)) + (1.0/2.0)*(4*Z + 4*bartau + 4*f_c*m)*pow(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m), 2)*(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2))/(pow(f_c - f_c0, 2)*pow(Z + bartau + f_c*m, 2)) - pow(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m), 2)*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(f_c - f_c0, 2)*pow(Z + bartau + f_c*m, 3)))/sqrt((pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(-Z - bartau - f_c*m, 2)*(Z + bartau + f_c0*m)*(Z + bartau + 2*f_c*m - f_c0*m)) + pow(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m), 2)*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(f_c - f_c0, 2)*pow(Z + bartau + f_c*m, 2)));
   }
   else {
      out_6697882517837858595[3] = -1;
   }
   if ((sqrt(pow(X_x + tau_pi_x/(omega_s - 1), 2) + pow(X_y + tau_pi_y/(omega_s - 1), 2)) - (Z + bartau)*(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m))*(Z + bartau - 2*f_t*m)/(pow(f_t, 2)*m))*(((f_t) > 0) - ((f_t) < 0))*(((m) > 0) - ((m) < 0)) < 0) {
      out_6697882517837858595[4] = (1.0/2.0)*(2*X_x - 2*tau_pi_x/(1 - omega_s))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/((Z + bartau)*sqrt((pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/((Z + bartau)*pow(-Z - bartau + f_t*m, 2)*(Z + bartau - 2*f_t*m)) + pow(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m), 2)*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/(pow(f_t, 2)*pow(Z + bartau - f_t*m, 2)))*pow(-Z - bartau + f_t*m, 2)*(Z + bartau - 2*f_t*m));
   }
   else if ((sqrt(pow(X_x + tau_pi_x/(omega_s - 1), 2) + pow(X_y + tau_pi_y/(omega_s - 1), 2)) - (Z + bartau + f_c0*m)*(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m))*(Z + bartau + 2*f_c*m - f_c0*m)/(m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))))*(((m) > 0) - ((m) < 0))*(((-f_c + f_c0) > 0) - ((-f_c + f_c0) < 0)) < 0) {
      out_6697882517837858595[4] = (1.0/2.0)*(2*X_x - 2*tau_pi_x/(1 - omega_s))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(sqrt((pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(-Z - bartau - f_c*m, 2)*(Z + bartau + f_c0*m)*(Z + bartau + 2*f_c*m - f_c0*m)) + pow(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m), 2)*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(f_c - f_c0, 2)*pow(Z + bartau + f_c*m, 2)))*pow(-Z - bartau - f_c*m, 2)*(Z + bartau + f_c0*m)*(Z + bartau + 2*f_c*m - f_c0*m));
   }
   else {
      out_6697882517837858595[4] = (X_x - tau_pi_x/(1 - omega_s))/sqrt(pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2));
   }
   if ((sqrt(pow(X_x + tau_pi_x/(omega_s - 1), 2) + pow(X_y + tau_pi_y/(omega_s - 1), 2)) - (Z + bartau)*(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m))*(Z + bartau - 2*f_t*m)/(pow(f_t, 2)*m))*(((f_t) > 0) - ((f_t) < 0))*(((m) > 0) - ((m) < 0)) < 0) {
      out_6697882517837858595[5] = (1.0/2.0)*(2*X_y - 2*tau_pi_y/(1 - omega_s))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/((Z + bartau)*sqrt((pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/((Z + bartau)*pow(-Z - bartau + f_t*m, 2)*(Z + bartau - 2*f_t*m)) + pow(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m), 2)*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/(pow(f_t, 2)*pow(Z + bartau - f_t*m, 2)))*pow(-Z - bartau + f_t*m, 2)*(Z + bartau - 2*f_t*m));
   }
   else if ((sqrt(pow(X_x + tau_pi_x/(omega_s - 1), 2) + pow(X_y + tau_pi_y/(omega_s - 1), 2)) - (Z + bartau + f_c0*m)*(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m))*(Z + bartau + 2*f_c*m - f_c0*m)/(m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))))*(((m) > 0) - ((m) < 0))*(((-f_c + f_c0) > 0) - ((-f_c + f_c0) < 0)) < 0) {
      out_6697882517837858595[5] = (1.0/2.0)*(2*X_y - 2*tau_pi_y/(1 - omega_s))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(sqrt((pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(-Z - bartau - f_c*m, 2)*(Z + bartau + f_c0*m)*(Z + bartau + 2*f_c*m - f_c0*m)) + pow(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m), 2)*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(f_c - f_c0, 2)*pow(Z + bartau + f_c*m, 2)))*pow(-Z - bartau - f_c*m, 2)*(Z + bartau + f_c0*m)*(Z + bartau + 2*f_c*m - f_c0*m));
   }
   else {
      out_6697882517837858595[5] = (X_y - tau_pi_y/(1 - omega_s))/sqrt(pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2));
   }
   out_6697882517837858595[6] = 0;
   out_6697882517837858595[7] = 0;
}
void get_df_dEps(double X_x, double X_y, double Z, double bartau, double sigma_pi, double tau_pi_x, double tau_pi_y, double f_c, double f_c0, double f_t, double m, double omega_s, double omega_w, double *out_5723301903165105508) {
   out_5723301903165105508[0] = 0;
   out_5723301903165105508[1] = 0;
   out_5723301903165105508[2] = 0;
   out_5723301903165105508[3] = 0;
   out_5723301903165105508[4] = 0;
   out_5723301903165105508[5] = 0;
   if ((sqrt(pow(X_x + tau_pi_x/(omega_s - 1), 2) + pow(X_y + tau_pi_y/(omega_s - 1), 2)) - (Z + bartau)*(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m))*(Z + bartau - 2*f_t*m)/(pow(f_t, 2)*m))*(((f_t) > 0) - ((f_t) < 0))*(((m) > 0) - ((m) < 0)) < 0) {
      out_5723301903165105508[6] = (1.0/2.0)*(2*tau_pi_x*(-X_x + tau_pi_x/(1 - omega_s))/pow(1 - omega_s, 2) + 2*tau_pi_y*(-X_y + tau_pi_y/(1 - omega_s))/pow(1 - omega_s, 2))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/((Z + bartau)*sqrt((pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/((Z + bartau)*pow(-Z - bartau + f_t*m, 2)*(Z + bartau - 2*f_t*m)) + pow(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m), 2)*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/(pow(f_t, 2)*pow(Z + bartau - f_t*m, 2)))*pow(-Z - bartau + f_t*m, 2)*(Z + bartau - 2*f_t*m));
   }
   else if ((sqrt(pow(X_x + tau_pi_x/(omega_s - 1), 2) + pow(X_y + tau_pi_y/(omega_s - 1), 2)) - (Z + bartau + f_c0*m)*(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m))*(Z + bartau + 2*f_c*m - f_c0*m)/(m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))))*(((m) > 0) - ((m) < 0))*(((-f_c + f_c0) > 0) - ((-f_c + f_c0) < 0)) < 0) {
      out_5723301903165105508[6] = (1.0/2.0)*(2*tau_pi_x*(-X_x + tau_pi_x/(1 - omega_s))/pow(1 - omega_s, 2) + 2*tau_pi_y*(-X_y + tau_pi_y/(1 - omega_s))/pow(1 - omega_s, 2))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(sqrt((pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(-Z - bartau - f_c*m, 2)*(Z + bartau + f_c0*m)*(Z + bartau + 2*f_c*m - f_c0*m)) + pow(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m), 2)*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(f_c - f_c0, 2)*pow(Z + bartau + f_c*m, 2)))*pow(-Z - bartau - f_c*m, 2)*(Z + bartau + f_c0*m)*(Z + bartau + 2*f_c*m - f_c0*m));
   }
   else {
      out_5723301903165105508[6] = (tau_pi_x*(-X_x + tau_pi_x/(1 - omega_s))/pow(1 - omega_s, 2) + tau_pi_y*(-X_y + tau_pi_y/(1 - omega_s))/pow(1 - omega_s, 2))/sqrt(pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2));
   }
   if ((sqrt(pow(X_x + tau_pi_x/(omega_s - 1), 2) + pow(X_y + tau_pi_y/(omega_s - 1), 2)) - (Z + bartau)*(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m))*(Z + bartau - 2*f_t*m)/(pow(f_t, 2)*m))*(((f_t) > 0) - ((f_t) < 0))*(((m) > 0) - ((m) < 0)) < 0) {
      out_5723301903165105508[7] = sigma_pi*(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/(pow(f_t, 2)*pow(1 - omega_w, 2)*sqrt((pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/((Z + bartau)*pow(-Z - bartau + f_t*m, 2)*(Z + bartau - 2*f_t*m)) + pow(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m), 2)*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/(pow(f_t, 2)*pow(Z + bartau - f_t*m, 2)))*pow(Z + bartau - f_t*m, 2));
   }
   else if ((sqrt(pow(X_x + tau_pi_x/(omega_s - 1), 2) + pow(X_y + tau_pi_y/(omega_s - 1), 2)) - (Z + bartau + f_c0*m)*(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m))*(Z + bartau + 2*f_c*m - f_c0*m)/(m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))))*(((m) > 0) - ((m) < 0))*(((-f_c + f_c0) > 0) - ((-f_c + f_c0) < 0)) < 0) {
      out_5723301903165105508[7] = sigma_pi*(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(1 - omega_w, 2)*pow(f_c - f_c0, 2)*sqrt((pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(-Z - bartau - f_c*m, 2)*(Z + bartau + f_c0*m)*(Z + bartau + 2*f_c*m - f_c0*m)) + pow(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m), 2)*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(f_c - f_c0, 2)*pow(Z + bartau + f_c*m, 2)))*pow(Z + bartau + f_c*m, 2));
   }
   else {
      out_5723301903165105508[7] = sigma_pi*m/pow(1 - omega_w, 2);
   }
}
void get_Phi(double S_s, double S_w, double X_x, double X_y, double Y_s, double Y_w, double Z, double bartau, double sigma_pi, double tau_pi_x, double tau_pi_y, double c_s, double c_w, double eta, double f_c, double f_c0, double f_t, double m, double omega_s, double omega_w, double *out_5145644395765063632) {
   if ((sqrt(pow(X_x + tau_pi_x/(omega_s - 1), 2) + pow(X_y + tau_pi_y/(omega_s - 1), 2)) - (Z + bartau)*(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m))*(Z + bartau - 2*f_t*m)/(pow(f_t, 2)*m))*(((f_t) > 0) - ((f_t) < 0))*(((m) > 0) - ((m) < 0)) < 0) {
      out_5145644395765063632[0] = (-X_x + tau_pi_x/(1 - omega_s))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/((1 - omega_s)*(Z + bartau)*sqrt((pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/((Z + bartau)*pow(-Z - bartau + f_t*m, 2)*(Z + bartau - 2*f_t*m)) + pow(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m), 2)*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/(pow(f_t, 2)*pow(Z + bartau - f_t*m, 2)))*pow(-Z - bartau + f_t*m, 2)*(Z + bartau - 2*f_t*m));
   }
   else if ((sqrt(pow(X_x + tau_pi_x/(omega_s - 1), 2) + pow(X_y + tau_pi_y/(omega_s - 1), 2)) - (Z + bartau + f_c0*m)*(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m))*(Z + bartau + 2*f_c*m - f_c0*m)/(m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))))*(((m) > 0) - ((m) < 0))*(((-f_c + f_c0) > 0) - ((-f_c + f_c0) < 0)) < 0) {
      out_5145644395765063632[0] = (-X_x + tau_pi_x/(1 - omega_s))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/((1 - omega_s)*sqrt((pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(-Z - bartau - f_c*m, 2)*(Z + bartau + f_c0*m)*(Z + bartau + 2*f_c*m - f_c0*m)) + pow(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m), 2)*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(f_c - f_c0, 2)*pow(Z + bartau + f_c*m, 2)))*pow(-Z - bartau - f_c*m, 2)*(Z + bartau + f_c0*m)*(Z + bartau + 2*f_c*m - f_c0*m));
   }
   else {
      out_5145644395765063632[0] = (-X_x + tau_pi_x/(1 - omega_s))/((1 - omega_s)*sqrt(pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2)));
   }
   if ((sqrt(pow(X_x + tau_pi_x/(omega_s - 1), 2) + pow(X_y + tau_pi_y/(omega_s - 1), 2)) - (Z + bartau)*(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m))*(Z + bartau - 2*f_t*m)/(pow(f_t, 2)*m))*(((f_t) > 0) - ((f_t) < 0))*(((m) > 0) - ((m) < 0)) < 0) {
      out_5145644395765063632[1] = (-X_y + tau_pi_y/(1 - omega_s))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/((1 - omega_s)*(Z + bartau)*sqrt((pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/((Z + bartau)*pow(-Z - bartau + f_t*m, 2)*(Z + bartau - 2*f_t*m)) + pow(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m), 2)*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/(pow(f_t, 2)*pow(Z + bartau - f_t*m, 2)))*pow(-Z - bartau + f_t*m, 2)*(Z + bartau - 2*f_t*m));
   }
   else if ((sqrt(pow(X_x + tau_pi_x/(omega_s - 1), 2) + pow(X_y + tau_pi_y/(omega_s - 1), 2)) - (Z + bartau + f_c0*m)*(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m))*(Z + bartau + 2*f_c*m - f_c0*m)/(m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))))*(((m) > 0) - ((m) < 0))*(((-f_c + f_c0) > 0) - ((-f_c + f_c0) < 0)) < 0) {
      out_5145644395765063632[1] = (-X_y + tau_pi_y/(1 - omega_s))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/((1 - omega_s)*sqrt((pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(-Z - bartau - f_c*m, 2)*(Z + bartau + f_c0*m)*(Z + bartau + 2*f_c*m - f_c0*m)) + pow(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m), 2)*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(f_c - f_c0, 2)*pow(Z + bartau + f_c*m, 2)))*pow(-Z - bartau - f_c*m, 2)*(Z + bartau + f_c0*m)*(Z + bartau + 2*f_c*m - f_c0*m));
   }
   else {
      out_5145644395765063632[1] = (-X_y + tau_pi_y/(1 - omega_s))/((1 - omega_s)*sqrt(pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2)));
   }
   if ((sqrt(pow(X_x + tau_pi_x/(omega_s - 1), 2) + pow(X_y + tau_pi_y/(omega_s - 1), 2)) - (Z + bartau)*(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m))*(Z + bartau - 2*f_t*m)/(pow(f_t, 2)*m))*(((f_t) > 0) - ((f_t) < 0))*(((m) > 0) - ((m) < 0)) < 0) {
      out_5145644395765063632[2] = (sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/(pow(f_t, 2)*(1 - omega_w)*sqrt((pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/((Z + bartau)*pow(-Z - bartau + f_t*m, 2)*(Z + bartau - 2*f_t*m)) + pow(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m), 2)*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/(pow(f_t, 2)*pow(Z + bartau - f_t*m, 2)))*pow(Z + bartau - f_t*m, 2));
   }
   else if ((sqrt(pow(X_x + tau_pi_x/(omega_s - 1), 2) + pow(X_y + tau_pi_y/(omega_s - 1), 2)) - (Z + bartau + f_c0*m)*(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m))*(Z + bartau + 2*f_c*m - f_c0*m)/(m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))))*(((m) > 0) - ((m) < 0))*(((-f_c + f_c0) > 0) - ((-f_c + f_c0) < 0)) < 0) {
      out_5145644395765063632[2] = (sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/((1 - omega_w)*pow(f_c - f_c0, 2)*sqrt((pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(-Z - bartau - f_c*m, 2)*(Z + bartau + f_c0*m)*(Z + bartau + 2*f_c*m - f_c0*m)) + pow(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m), 2)*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(f_c - f_c0, 2)*pow(Z + bartau + f_c*m, 2)))*pow(Z + bartau + f_c*m, 2));
   }
   else {
      out_5145644395765063632[2] = m/(1 - omega_w);
   }
   out_5145644395765063632[3] = -(((sqrt(pow(X_x + tau_pi_x/(omega_s - 1), 2) + pow(X_y + tau_pi_y/(omega_s - 1), 2)) - (Z + bartau)*(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m))*(Z + bartau - 2*f_t*m)/(pow(f_t, 2)*m))*(((f_t) > 0) - ((f_t) < 0))*(((m) > 0) - ((m) < 0)) < 0) ? (
      -(2*Z + 2*bartau - 2*f_t*m)/(Z + bartau - 2*f_t*m) + (pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2))/pow(Z + bartau - 2*f_t*m, 2) + (-m*(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/(pow(Z + bartau - 2*f_t*m, 2)*pow(Z + bartau - f_t*m, 2)) + (1.0/2.0)*(pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*(4*Z + 4*bartau - 4*f_t*m)*(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2))/((Z + bartau)*pow(-Z - bartau + f_t*m, 2)*(Z + bartau - 2*f_t*m)) - 1.0/2.0*(pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/((Z + bartau)*pow(-Z - bartau + f_t*m, 2)*pow(Z + bartau - 2*f_t*m, 2)) + (pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/((Z + bartau)*pow(-Z - bartau + f_t*m, 3)*(Z + bartau - 2*f_t*m)) - 1.0/2.0*(pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/(pow(Z + bartau, 2)*pow(-Z - bartau + f_t*m, 2)*(Z + bartau - 2*f_t*m)) + (1.0/2.0)*pow(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m), 2)*(4*Z + 4*bartau - 4*f_t*m)*(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2))/(pow(f_t, 2)*pow(Z + bartau - f_t*m, 2)) - pow(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m), 2)*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/(pow(f_t, 2)*pow(Z + bartau - f_t*m, 3)))/sqrt((pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/((Z + bartau)*pow(-Z - bartau + f_t*m, 2)*(Z + bartau - 2*f_t*m)) + pow(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m), 2)*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/(pow(f_t, 2)*pow(Z + bartau - f_t*m, 2)))
   )
   : (((sqrt(pow(X_x + tau_pi_x/(omega_s - 1), 2) + pow(X_y + tau_pi_y/(omega_s - 1), 2)) - (Z + bartau + f_c0*m)*(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m))*(Z + bartau + 2*f_c*m - f_c0*m)/(m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))))*(((m) > 0) - ((m) < 0))*(((-f_c + f_c0) > 0) - ((-f_c + f_c0) < 0)) < 0) ? (
      -(2*Z + 2*bartau + 2*f_c*m)/(Z + bartau + 2*f_c*m - f_c0*m) + (pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2))/pow(Z + bartau + 2*f_c*m - f_c0*m, 2) + (-m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))*(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(f_c - f_c0, 2)*pow(Z + bartau + f_c*m, 2)*pow(Z + bartau + 2*f_c*m - f_c0*m, 2)) + (1.0/2.0)*(pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*(4*Z + 4*bartau + 4*f_c*m)*(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2))/(pow(-Z - bartau - f_c*m, 2)*(Z + bartau + f_c0*m)*(Z + bartau + 2*f_c*m - f_c0*m)) - 1.0/2.0*(pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(-Z - bartau - f_c*m, 2)*(Z + bartau + f_c0*m)*pow(Z + bartau + 2*f_c*m - f_c0*m, 2)) - 1.0/2.0*(pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(-Z - bartau - f_c*m, 2)*pow(Z + bartau + f_c0*m, 2)*(Z + bartau + 2*f_c*m - f_c0*m)) + (pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(-Z - bartau - f_c*m, 3)*(Z + bartau + f_c0*m)*(Z + bartau + 2*f_c*m - f_c0*m)) + (1.0/2.0)*(4*Z + 4*bartau + 4*f_c*m)*pow(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m), 2)*(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2))/(pow(f_c - f_c0, 2)*pow(Z + bartau + f_c*m, 2)) - pow(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m), 2)*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(f_c - f_c0, 2)*pow(Z + bartau + f_c*m, 3)))/sqrt((pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(-Z - bartau - f_c*m, 2)*(Z + bartau + f_c0*m)*(Z + bartau + 2*f_c*m - f_c0*m)) + pow(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m), 2)*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(f_c - f_c0, 2)*pow(Z + bartau + f_c*m, 2)))
   )
   : (
      -1
   )));
   out_5145644395765063632[4] = -(((sqrt(pow(X_x + tau_pi_x/(omega_s - 1), 2) + pow(X_y + tau_pi_y/(omega_s - 1), 2)) - (Z + bartau)*(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m))*(Z + bartau - 2*f_t*m)/(pow(f_t, 2)*m))*(((f_t) > 0) - ((f_t) < 0))*(((m) > 0) - ((m) < 0)) < 0) ? (
      (1.0/2.0)*(2*X_x - 2*tau_pi_x/(1 - omega_s))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/((Z + bartau)*sqrt((pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/((Z + bartau)*pow(-Z - bartau + f_t*m, 2)*(Z + bartau - 2*f_t*m)) + pow(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m), 2)*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/(pow(f_t, 2)*pow(Z + bartau - f_t*m, 2)))*pow(-Z - bartau + f_t*m, 2)*(Z + bartau - 2*f_t*m))
   )
   : (((sqrt(pow(X_x + tau_pi_x/(omega_s - 1), 2) + pow(X_y + tau_pi_y/(omega_s - 1), 2)) - (Z + bartau + f_c0*m)*(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m))*(Z + bartau + 2*f_c*m - f_c0*m)/(m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))))*(((m) > 0) - ((m) < 0))*(((-f_c + f_c0) > 0) - ((-f_c + f_c0) < 0)) < 0) ? (
      (1.0/2.0)*(2*X_x - 2*tau_pi_x/(1 - omega_s))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(sqrt((pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(-Z - bartau - f_c*m, 2)*(Z + bartau + f_c0*m)*(Z + bartau + 2*f_c*m - f_c0*m)) + pow(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m), 2)*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(f_c - f_c0, 2)*pow(Z + bartau + f_c*m, 2)))*pow(-Z - bartau - f_c*m, 2)*(Z + bartau + f_c0*m)*(Z + bartau + 2*f_c*m - f_c0*m))
   )
   : (
      (X_x - tau_pi_x/(1 - omega_s))/sqrt(pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))
   )));
   out_5145644395765063632[5] = -(((sqrt(pow(X_x + tau_pi_x/(omega_s - 1), 2) + pow(X_y + tau_pi_y/(omega_s - 1), 2)) - (Z + bartau)*(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m))*(Z + bartau - 2*f_t*m)/(pow(f_t, 2)*m))*(((f_t) > 0) - ((f_t) < 0))*(((m) > 0) - ((m) < 0)) < 0) ? (
      (1.0/2.0)*(2*X_y - 2*tau_pi_y/(1 - omega_s))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/((Z + bartau)*sqrt((pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/((Z + bartau)*pow(-Z - bartau + f_t*m, 2)*(Z + bartau - 2*f_t*m)) + pow(sigma_pi/(1 - omega_w) + pow(f_t, 2)*m/(Z + bartau - 2*f_t*m), 2)*pow(pow(f_t, 2)*pow(m, 2) - 2*f_t*m*(Z + bartau) + pow(Z + bartau, 2), 2)/(pow(f_t, 2)*pow(Z + bartau - f_t*m, 2)))*pow(-Z - bartau + f_t*m, 2)*(Z + bartau - 2*f_t*m))
   )
   : (((sqrt(pow(X_x + tau_pi_x/(omega_s - 1), 2) + pow(X_y + tau_pi_y/(omega_s - 1), 2)) - (Z + bartau + f_c0*m)*(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m))*(Z + bartau + 2*f_c*m - f_c0*m)/(m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))))*(((m) > 0) - ((m) < 0))*(((-f_c + f_c0) > 0) - ((-f_c + f_c0) < 0)) < 0) ? (
      (1.0/2.0)*(2*X_y - 2*tau_pi_y/(1 - omega_s))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(sqrt((pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(-Z - bartau - f_c*m, 2)*(Z + bartau + f_c0*m)*(Z + bartau + 2*f_c*m - f_c0*m)) + pow(sigma_pi/(1 - omega_w) + f_c0 + m*(pow(f_c, 2) - 2*f_c*f_c0 + pow(f_c0, 2))/(Z + bartau + 2*f_c*m - f_c0*m), 2)*pow(pow(f_c, 2)*pow(m, 2) - 2*f_c*f_c0*pow(m, 2) + 2*f_c*m*(Z + bartau + f_c0*m) + pow(f_c0, 2)*pow(m, 2) - 2*f_c0*m*(Z + bartau + f_c0*m) + pow(Z + bartau + f_c0*m, 2), 2)/(pow(f_c - f_c0, 2)*pow(Z + bartau + f_c*m, 2)))*pow(-Z - bartau - f_c*m, 2)*(Z + bartau + f_c0*m)*(Z + bartau + 2*f_c*m - f_c0*m))
   )
   : (
      (X_y - tau_pi_y/(1 - omega_s))/sqrt(pow(-X_x + tau_pi_x/(1 - omega_s), 2) + pow(-X_y + tau_pi_y/(1 - omega_s), 2))
   )));
   out_5145644395765063632[6] = pow(1 - omega_s, c_s)*(2*Y_s/S_s + Y_w*eta/S_s) + Y_w*eta*pow(1 - omega_w, c_w)/S_w;
   out_5145644395765063632[7] = pow(1 - omega_w, c_w)*(Y_s*eta/S_w + 2*Y_w/S_w) + Y_s*eta*pow(1 - omega_s, c_s)/S_s;
}

/******************************************************************************
 ******************************************************************************/

 ////////////////////////////////////////////////////////////////////////////////
 // Slide material status
 ////////////////////////////////////////////////////////////////////////////////

Slide32MaterialStatus :: Slide32MaterialStatus(Slide32Material *m, Element *e) : DisMechMaterialStatus(m, e) {
   name = "slide 3.2 mat. status";
}


void Slide32MaterialStatus :: init() {
  // TODO
  Slide32Material *m = static_cast< Slide32Material * >( mat );
  if ( m->calc_slips() ){
    RigidBodyContact *rbc = static_cast< RigidBodyContact * >( element );
    strain_slip_multiplier = rbc->giveLength();
  } else {
    strain_slip_multiplier = 1.0;
  }
  max_iter = m->giveMax_iter();

  Eps = new double[8];
  temp_Eps = new double[8];

  Sig = new double[8];
  temp_Sig = new double[8];

  for (size_t i = 0; i < 8; i++) {
    Eps[i] = 0;
    Sig[i] = 0;
  }
}


void Slide32MaterialStatus :: update() {
  // TODO
  std :: copy(temp_Sig, temp_Sig+8, Sig);
  std :: copy(temp_Eps, temp_Eps+8, Eps);
}


void Slide32MaterialStatus :: print_state_vars(const bool &temp) {
  if ( temp ) {
    std::cout << "temp_Eps: " << '\t';
    for (size_t i = 0; i < 8; i++) {
      std::cout << temp_Eps[i] << "\t";
    }
    std::cout << '\n';
  } else {
    std::cout << "Eps: " << '\t';
    for (size_t i = 0; i < 8; i++) {
      std::cout << Eps[i] << "\t";
    }
    std::cout << '\n';
  }
}


void Slide32MaterialStatus :: print_thrmdyn_forces(const bool &temp) {
  if ( temp ) {
    std::cout << "temp_Sig: " << '\t';
    for (size_t i = 0; i < 8; i++) {
      std::cout << temp_Sig[i] << "\t";
    }
    std::cout << '\n';
  } else {
    std::cout << "Sig: " << '\t';
    for (size_t i = 0; i < 8; i++) {
      std::cout << Sig[i] << "\t";
    }
    std::cout << '\n';
  }
}


void print_array_of_doubles(const double *arr, const string &str="whatever", const size_t &range=8) {
    std::cout << str << ": " << '\t';
    for (size_t i = 0; i < range; i++) {
      std::cout << arr[i] << "\t";
    }
    std::cout << '\n';
}

void check_array_of_doubles_for_nans(double *arr, const size_t &range=8){
  for (size_t i = 0; i < range; i++) {
    if ( std :: isnan(arr[ i ])) {
      arr[ i ] = 0;
    }
  }
}

// FNS with temp variables:
///////////////////////////////////////////////////////////////////////
void Slide32MaterialStatus :: get_f_df(const double &s_x_n1, const double &s_y_n1, const double &w_n1, double &f_k, double &df_k) {
    Slide32Material *m = static_cast< Slide32Material * >( mat );
    // uses sympy generated fns

    // trial state variables
    // TODO move this outside this fn, move only those that are returned, delete other at the end
    get_Sig(m->giveE_s(), m->giveE0(), m->giveK_s(), temp_Eps[4], temp_Eps[5], m->giveGamma_s(), temp_Eps[6], temp_Eps[7], temp_Eps[0], temp_Eps[1], s_x_n1, s_y_n1, w_n1, temp_Eps[2], temp_Eps[3], temp_Sig);
    // check_array_of_doubles_for_nans(temp_Sig);


    double * dSig_dEps_k = new double[64];
    get_dSig_dEps(m->giveE_s(), m->giveE0(), m->giveK_s(), m->giveGamma_s(), temp_Eps[6], temp_Eps[7], temp_Eps[0], temp_Eps[1], s_x_n1, s_y_n1, w_n1, temp_Eps[2], dSig_dEps_k);
    // check_array_of_doubles_for_nans(dSig_dEps_k, 64);


    f_k = get_f(temp_Sig[4], temp_Sig[5], temp_Sig[3], m->giveBartau(), temp_Sig[2], temp_Sig[0], temp_Sig[1], m->giveF_c(), m->giveF_c0(), m->giveF_t(), m->giveM(), temp_Eps[6], temp_Eps[7]);

    double * df_dSig_k = new double[8];
    get_df_dSig(temp_Sig[4], temp_Sig[5], temp_Sig[3], m->giveBartau(), temp_Sig[2], temp_Sig[0], temp_Sig[1], m->giveF_c(), m->giveF_c0(), m->giveF_t(), m->giveM(), temp_Eps[6], temp_Eps[7], df_dSig_k);
    // check_array_of_doubles_for_nans(df_dSig_k);

    double * ddf_dEps_k = new double[8];
    get_df_dEps(temp_Sig[4], temp_Sig[5], temp_Sig[3], m->giveBartau(), temp_Sig[2], temp_Sig[0], temp_Sig[1], m->giveF_c(), m->giveF_c0(), m->giveF_t(), m->giveM(), temp_Eps[6], temp_Eps[7], ddf_dEps_k);
    // check_array_of_doubles_for_nans(ddf_dEps_k);

    double * df_dEps_k = new double[8];
    for ( unsigned i = 0; i < 8; i++){
      df_dEps_k[ i ] = 0.0;
      for ( unsigned j = 0; j < 8; j++){
        df_dEps_k[ i ] += df_dSig_k[ j ] * dSig_dEps_k[ j + i*8 ];
      }
      df_dEps_k[ i ] += ddf_dEps_k[ i ];
    }
    // check_array_of_doubles_for_nans(df_dEps_k);

    double * dEps_dlambda_k = new double[8];
    get_Phi(m->giveS_s(), m->giveS_w(), temp_Sig[4], temp_Sig[5], temp_Sig[6], temp_Sig[7], temp_Sig[3], m->giveBartau(), temp_Sig[2], temp_Sig[0], temp_Sig[1], m->giveC_s(), m->giveC_w(), m->giveEta(), m->giveF_c(), m->giveF_c0(), m->giveF_t(), m->giveM(), temp_Eps[6], temp_Eps[7], dEps_dlambda_k);
    // check_array_of_doubles_for_nans(dEps_dlambda_k);

    // print_array_of_doubles(df_dEps_k, "df_dEps_k");
    // print_array_of_doubles(dEps_dlambda_k, "dEps_dlambda_k");
    df_k = 0;
    for ( unsigned i = 0; i < 8; i++){
      df_k += df_dEps_k[ i ] * dEps_dlambda_k[ i ];

    }
    delete [] dSig_dEps_k;
    delete [] df_dSig_k;
    delete [] ddf_dEps_k;
    delete [] df_dEps_k;
    delete [] dEps_dlambda_k;
}

void Slide32MaterialStatus :: get_Eps_k1(const double &s_x_n1, const double &s_y_n1, const double &w_n1, const double &lam_k
) {
    Slide32Material *m = static_cast< Slide32Material * >( mat );

    get_Sig(m->giveE_s(), m->giveE0(), m->giveK_s(), temp_Eps[4], temp_Eps[5], m->giveGamma_s(), temp_Eps[6], temp_Eps[7], temp_Eps[0], temp_Eps[1], s_x_n1, s_y_n1, w_n1, temp_Eps[2], temp_Eps[3], temp_Sig);
    // check_array_of_doubles_for_nans(temp_Sig);

    double * Phi_k = new double[8];
    get_Phi(m->giveS_s(), m->giveS_w(), temp_Sig[4], temp_Sig[5], temp_Sig[6], temp_Sig[7], temp_Sig[3], m->giveBartau(), temp_Sig[2], temp_Sig[0], temp_Sig[1], m->giveC_s(), m->giveC_w(), m->giveEta(), m->giveF_c(), m->giveF_c0(), m->giveF_t(), m->giveM(), temp_Eps[6], temp_Eps[7], Phi_k);
    // check_array_of_doubles_for_nans(Phi_k);

    for ( unsigned i = 0; i < 8; i++){
      temp_Eps[ i ] = Eps[ i ] + lam_k * Phi_k[ i ];
    }
    delete Phi_k;
    // check_state_variable_ranges();
}


void Slide32MaterialStatus :: check_state_variable_ranges(){
  // damage in <0, 1>
  // if ( temp_Eps[6] < 0 ) {
  //   temp_Eps[6] = 0.0;
  // } else if ( temp_Eps[6] > 1.0 ) {
  //   temp_Eps[6] = 1.0;
  // }
  // if ( temp_Eps[7] < 0 ) {
  //   temp_Eps[7] = 0.0;
  // } else if ( temp_Eps[7] > 1.0 ) {
  //   temp_Eps[7] = 1.0;
  // }
  // damage in <0, 1)
  if ( temp_Eps[6] < 0 ) {
    temp_Eps[6] = 0.0;
  } else if ( temp_Eps[6] >= 1.0 ) {
    temp_Eps[6] = 1.0 - 1e-10;
  }
  if ( temp_Eps[7] < 0 ) {
    temp_Eps[7] = 0.0;
  } else if ( temp_Eps[7] >= 1.0 ) {
    temp_Eps[7] = 1.0 - 1e-10;
  }
  // Z is nonzero
  if ( temp_Sig[3] < 0 ){
    temp_Sig[3] = 0.0;
  }
}


Vector Slide32MaterialStatus :: giveStress(const Vector &strain
){
    Vector stress(strain.size());

    double s_x_n1, s_y_n1, w_n1;
    w_n1 = strain[0] * strain_slip_multiplier;
    s_x_n1 = strain[1] * strain_slip_multiplier;
    if ( strain.size() > 2 ){
      s_y_n1 =  strain[2] * strain_slip_multiplier;
    } else {
      s_y_n1 = 0.0;
    }

    double f_k, df_k;

    std :: copy(Eps, Eps+8, temp_Eps);

    double dlam;
    double lam_k = 0;
    get_f_df(s_x_n1, s_y_n1, w_n1, f_k, df_k);
    double f_k_norm = abs(f_k);
    double f_k_trial = f_k;
    unsigned k = 0;
    while ( k < max_iter ){

        if ( f_k_trial < 0  || f_k_norm < 1e-6 ){
            break;
        }

        dlam = -f_k / df_k;
        if ( std :: isnan( dlam )) dlam = 0;
        lam_k += dlam;
        get_Eps_k1(s_x_n1, s_y_n1, w_n1, lam_k);
        get_f_df(s_x_n1, s_y_n1, w_n1, f_k, df_k);
        // f_k_trial = fk;
        f_k_norm = abs(f_k);
        k += 1;
    }
    if ( k >= max_iter ){
        std::cerr << "No convergence after " << k;
        std::cout << " num of iterations" << '\n';
        print_state_vars(true);
        print_thrmdyn_forces(true);
        // exit(1);
    }

    // check_state_variable_ranges();

    stress[ 0 ] = temp_Sig[2];
    stress[ 1 ] = temp_Sig[0];
    if ( stress.size() > 2 ){
      stress[ 2 ] = temp_Sig[ 1 ];
    }
    // std::cout << "num of iterations: " << k << '\n';
    return stress;
}




Matrix Slide32MaterialStatus :: giveStiffnessTensor(string type, unsigned dim) const {
  // TODO
  Matrix stiff = DisMechMaterialStatus :: giveStiffnessTensor(type, dim);
  if ( type.compare("elastic") == 0 ) {
      return stiff;
  } else if ( type.compare("secant") == 0 ) {
      return stiff;
  } else if ( type.compare("unloading") == 0 ) {
      return stiff;
  } else if ( type.compare("tangent") == 0 ) {
      return stiff; //not implemented, used unloading
  } else {
      cerr << "Error: Slide32MaterialStatus does not provide '" << type << "' stiffness";
      exit(1);
  };

}


double Slide32MaterialStatus :: giveValue(string code) const {
  // TODO
  return DisMechMaterialStatus :: giveValue(code);
}

////////////////////////////////////////////////////////////////////////////////
// Slide material
////////////////////////////////////////////////////////////////////////////////
void Slide32Material :: readFromLine(istringstream &iss) {
    DisMechMaterial :: readFromLine(iss); //read elastic parameters

    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream

    use_discrete_values = false;  // should be by default, but in windows, one never knows :D
    max_iter = 20;

    string param;
    bool bts, bks, bgs, bss, bcs, brs, bsw, brw, bcw, bft, bfc, bfc0, bm, bet;
    bts = bks = bgs = bss = bcs = brs = bsw = brw = bcw = bft = bfc = bfc0 = bm = bet = false;
    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("gamma_s") == 0 ||  param.compare("gamma") == 0 ) {
            bgs = true;
            iss >> gamma_s;
        } else if (  param.compare("K_s") == 0 || param.compare("Kin") == 0 ) {
            bks = true;
            iss >> K_s;
        } else if (  param.compare("bartau") == 0 || param.compare("tauBar") == 0 ) {
            bts = true;
            iss >> bartau;
        } else if ( param.compare("S_s") == 0 || param.compare("S") == 0 ) {
            bss = true;
            iss >> S_s;
        } else if ( param.compare("c_s") == 0 || param.compare("c") == 0 ) {
            bcs = true;
            iss >> c_s;
        } else if ( param.compare("r_s") == 0 || param.compare("r") == 0 ) {
            brs = true;
            iss >> r_s;
        //////////////////////////
        } else if (  param.compare("S_w") == 0) {
            bsw = true;
            iss >> S_w;
        } else if ( param.compare("r_w") == 0 ) {
            brw = true;
            iss >> r_w;
        } else if ( param.compare("c_w") == 0 ) {
            bcw = true;
            iss >> c_s;
        } else if (  param.compare("f_t") == 0 ) {
            bft = true;
            iss >> f_t;
        } else if (  param.compare("f_c") == 0 ) {
            bfc = true;
            iss >> f_c;
        } else if (  param.compare("f_c0") == 0 ) {
            bfc0 = true;
            iss >> f_c0;
        } else if ( param.compare("m") == 0 ) {
            bm = true;
            iss >> m;
        } else if ( param.compare("eta") == 0 ) {
            bet = true;
            iss >> eta;
        } else if ( param.compare("use_slips") == 0 ||  param.compare("use_slip") == 0  ) {
            use_discrete_values = true;
        } else if ( param.compare("max_iter") == 0 ) {
            iss >> max_iter;
        }
    }
    if ( !bts ) {
        cerr << name << ": material parameter 'bartau' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    if ( !bks ) {
        cout << name << ": material parameter 'K_s'was not specified, taking zero" << endl;
        K_s = 0;
    }
    if ( !bgs ) {
        cerr << name << ": material parameter 'gamma_s' was not specified, taking zero" << endl;
        gamma_s = 0;
    }
    if ( !bss ) {
        cerr << name << ": material parameter 'S_s' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    if ( !bcs ) {
        cout << name << ": material parameter 'c_s' was not specified, taking c = 1.0" << endl;
        c_s = 1.0;
    }
    if ( !brs ) {
        cout << name << ": material parameter 'r_s' was not specified, taking r = 1.0" << endl;
        r_s = 1.0;
    }


    if ( !bsw ) {
        cerr << name << ": material parameter 'S_w' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    if ( !bcw ) {
        cout << name << ": material parameter 'c_w' was not specified, taking c = 1.0" << endl;
        c_w = 1.0;
    }
    if ( !brw ) {
        cout << name << ": material parameter 'r_w' was not specified, taking r = 1.0" << endl;
        r_w = 1.0;
    }
    if ( !bft ) {
        cerr << name << ": material parameter 'f_t' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    if ( !bfc ) {
        cerr << name << ": material parameter 'f_c' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    if ( !bfc0 ) {
        cerr << name << ": material parameter 'f_c0' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    if ( f_c0 < f_c ) {
        std::cerr << "material param f_c0 = " << f_c0 << " must not be greater than parameter f_c = " << f_c << '\n';
        exit(EXIT_FAILURE);
        // TODO similar must be check of taubar vs f_t etc, some combinations does not make sense for the cap function
    }
    if ( !bm ) {
        cerr << name << ": material parameter 'm' was not specified" << endl;
        // cout << name << ": material parameter 'm' was not specified, taking r = 1.0" << endl;
        // m = 0;
        exit(EXIT_FAILURE);
    }
    if ( !bet ) {
        cerr << name << ": material parameter 'eta' was not specified" << endl;
        // cout << name << ": material parameter 'eta' was not specified, taking r = 1.0" << endl;
        // eta = 0;
        exit(EXIT_FAILURE);
    }
}


MaterialStatus *Slide32Material :: giveNewMaterialStatus(Element *e) {
    Slide32MaterialStatus *newStatus = new Slide32MaterialStatus(this, e); //needs to be deleted manually
    return newStatus;
};
