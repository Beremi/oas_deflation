#include "material.h"

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT MATERIAL

TrsprtMaterialStatus::TrsprtMaterialStatus(TrsprtMaterial *m){
    name="transport mat. status";
    mat=m;
}

//////////////////////////////////////////////////////////
double TrsprtMaterialStatus::giveConductivity() const{
    TrsprtMaterial *tmat = static_cast<TrsprtMaterial *>(mat);
    return tmat->giveConductivity();
} 

//////////////////////////////////////////////////////////
double TrsprtMaterialStatus::giveCapacity() const{
    TrsprtMaterial *tmat = static_cast<TrsprtMaterial *>(mat);
    return tmat->giveCapacity();
} 

//////////////////////////////////////////////////////////
void TrsprtMaterial::readFromLine(istringstream &iss){
    string param;
    bool bcapacity,bconductivity;
    bcapacity=bconductivity = false;

    while (!iss.eof()) {
        iss >> param;
        if (param.compare("capacity") == 0){
            bcapacity = true;
            iss >> capacity;
        }else if (param.compare("conductivity") == 0){
            bconductivity = true;
            iss >> conductivity;
        }        
    }    
    if (not bcapacity) {cerr << name << ": material parameter 'capacity' was not specified" << endl; exit(0);};
    if (not bconductivity) {cerr << name << ": material parameter 'conductivity' was not specified" << endl; exit(0);};
};

//////////////////////////////////////////////////////////
MaterialStatus* TrsprtMaterial::giveNewMaterialStatus(){
    TrsprtMaterialStatus* newStatus = new TrsprtMaterialStatus(this); //needs to be deleted manually
    return newStatus;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE MECHANICAL ELASTIC MATERIAL

DisMechMaterialStatus::DisMechMaterialStatus(DisMechMaterial *m){
    name="discrete mechanical mat. status";
    mat=m;
}

//////////////////////////////////////////////////////////
vector<double> DisMechMaterialStatus::giveNormalShearStiffness() const{
    DisMechMaterial *tmat = static_cast<DisMechMaterial *>(mat);
    return tmat->giveNormalShearStiffness();
} 

//////////////////////////////////////////////////////////
double DisMechMaterialStatus::giveDensity() const{
    DisMechMaterial *tmat = static_cast<DisMechMaterial *>(mat);
    return tmat->giveDensity();
} 

//////////////////////////////////////////////////////////
void DisMechMaterial::readFromLine(istringstream &iss){
    string param;

    bool bE0,balpha,bdensity;
    bE0=balpha=bdensity = false;

    while (!iss.eof()) {
        iss >> param;
        if (param.compare("E0") == 0){
            bE0 = true;
            iss >> E0;
        }else if (param.compare("alpha") == 0){
            balpha = true;
            iss >> alpha;
        }else if (param.compare("density") == 0){
            bdensity = true;
            iss >> density;
        }        
    }    
    if (not bE0) {cerr << name << ": material parameter 'E0' was not specified" << endl; exit(0);};
    if (not balpha) {cerr << name << ": material parameter 'alpha' was not specified" << endl; exit(0);};
    if (not bdensity) {cerr << name << ": material parameter 'density' was not specified" << endl; exit(0);};

};

//////////////////////////////////////////////////////////
MaterialStatus* DisMechMaterial::giveNewMaterialStatus(){
    DisMechMaterialStatus* newStatus = new DisMechMaterialStatus(this); //needs to be deleted manually
    return newStatus;
};

/*
//////////////////////////////////////////////////////////
// CUSATIS MATERIAL

MarsMaterialStatus::MarsMaterialStatus(MatsMaterial *m){
    name="MARS mat. status";
    mat=m;
    RAND_H=1.0;
}


void calculateDamage(eps);
    double Eps = sqrt(eps[0] * eps[0] + cLaw.alpha * (eps[1] * eps[1] + eps[2] * eps[2])); //equivalent strain
    double epsN = eps[0];
    double epsT = sqrt(eps[1] * eps[1] + eps[2] * eps[2]);

    double s, s2, c2;

    if (Eps > 0 && E.D < 1.0 && E.NL == true) {
        if (epsT > 0) omega = atan(epsN / (sqrt(cLaw.alpha) * epsT));
        else omega = 0.5 * M_PI * epsN / fabs(epsN);
        s = sin(omega);
        s2 = s*s;
        c2 = cos(omega) * cos(omega);

        double St = cLaw.St * E.rf;
        double Gt = cLaw.Gt * E.rf * E.rf;
        double Ss = cLaw.Ss * E.rf;
        double GS = cLaw.Gs * E.rf * E.rf;
        double sa;

        if (cLaw.type.compare("Cusatis") == 0) {
            if (omega < cLaw.omega0) {
                //compression
                S0 = cLaw.Sc / sqrt(s2 + (cLaw.alpha * c2) / cLaw.beta);
                chi = Eps;
                K0 = cLaw.Kc * E.stiffness * (1. - pow((omega + 0.5 * M_PI) / (cLaw.omega0 + 0.5 * M_PI), cLaw.nc));
            } else {
                //tension-shear
                epsmaxT = E.EpsMaxT;
                epsmaxN = E.EpsMaxN;
                if (E.EpsMaxN < epsN) epsmaxN = epsN;
                if (E.EpsMaxT < epsT) epsmaxT = epsT;
                emax = sqrt(epsmaxN*epsmaxN + cLaw.alpha*epsmaxT*epsmaxT);
                sa = .5 * St * (Ss2 / (cLaw.mu2 * St2) - 1.);
                S0 = (-(St + sa) * s + sqrt(pow((St + sa) * s, 2.)+(cLaw.alpha * (c2 / cLaw.mu2) - s2)*(St + 2*sa) * St)) / (cLaw.alpha * (c2 / cLaw.mu2) - s2);
                if (fabs(confiningStrain) / cLaw.lambda < 1. || confiningStrain > 0.) f = 1.;
                else f = 1. / (1. - confiningStrain / cLaw.lambda);
                K0 = -f * E.Kt0 * (1. - pow((omega - 0.5 * M_PI) / (cLaw.omega0 - 0.5 * M_PI), E.nt0));
                if (omega < 0.0) chi = Eps * omega / cLaw.omega0 + emax * (1. - omega / cLaw.omega0);
                else chi = emax;                
            }
            if (chi - S0 / E.stiffness > 0) S = S0 * exp(K0 / S0 * (chi - S0 / E.stiffness));
            else S = S0;
            D = 1. - S / (E.stiffness * Eps);
        }
        
        if (cLaw.type.compare("Bilinear") == 0) {
            //compression
            if (omega < cLaw.omega0) {
                S0 = cLaw.Sc / sqrt(s2 + (cLaw.alpha * c2) / cLaw.beta);
                chi = Eps;
                K0 = cLaw.Kc * E.stiffness * (1. - pow((omega + 0.5 * M_PI) / (cLaw.omega0 + 0.5 * M_PI), cLaw.nc));
                if (chi - S0 / E.stiffness > 0) S = S0 * exp(K0 / S0 * (chi - S0 / E.stiffness));
                else S = S0;
                D = 1. - S / (E.stiffness * Eps);
            } else {
            //tension-shear
                epsmaxT = E.EpsMaxT;
                epsmaxN = E.EpsMaxN;
                if (E.EpsMaxN < epsN) epsmaxN = epsN;
                if (E.EpsMaxT < epsT) epsmaxT = epsT;
                emax = sqrt(epsmaxN*epsmaxN + cLaw.alpha*epsmaxT*epsmaxT);
                sa = .5 * St * (Ss2 / (cLaw.mu2 * St2) - 1.);
                S0 = (-(St + sa) * s + sqrt(pow((St + sa) * s, 2.)+(cLaw.alpha * (c2 / cLaw.mu2) - s2)*(St + 2*sa) * St)) / (cLaw.alpha * (c2 / cLaw.mu2) - s2);
                if (fabs(confiningStrain) / cLaw.lambda < 1. || confiningStrain > 0.) f = 1.;
                else f = 1. / (1. - confiningStrain / cLaw.lambda);
                K0 = -f * E.Kt0 * (1. - pow((omega - 0.5 * M_PI) / (cLaw.omega0 - 0.5 * M_PI), E.nt0));
                K1 = -f * E.Kt1 * (1. - pow((omega - 0.5 * M_PI) / (cLaw.omega0 - 0.5 * M_PI), E.nt1));
                if (omega < 0.0) chi = Eps * omega / cLaw.omega0 + emax * (1. - omega / cLaw.omega0);
                else chi = emax;
                double e0 = S0 / E.stiffness;
                double e1 = e0 - S0*(1-cLaw.delta)/K0;
                double ec = e1 - S0*cLaw.delta/K1;
                if (chi <= e0) S = S0;
                else if (chi <= e1) S = S0+K0*(chi-e0);
                else if (chi <= ec) S = cLaw.delta*S0 + K1*(chi-e1);
                else S = 0;
                D = 1. - S / (E.stiffness * Eps);
            }
        }
    } else D = 0.0;

    // no recovery	
    if (D < E.D_t) D = E.D_t;    	
 

    // compression recovery
    //if (D < E.D_t && epsN >= 0) D = E.D_t;    	
    //if (D < E.D_c && epsN < 0) D = E.D_c;

    return D;
*/
//////////////////////////////////////////////////////////
// COUPLED MECHANICAL MATERIAL


