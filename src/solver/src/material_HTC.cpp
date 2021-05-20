#include "material_HTC.h"
#include "element_continuous.h"


//////////////////////////////////////////////////////////
// HTC MATERIAL STATUS

HTCMaterialStatus :: HTCMaterialStatus(HTCMaterial *m, Element *e, unsigned ipnum) : TrsprtMaterialStatus(m, e, ipnum) {
    name = "HTC mat. status";
}

//////////////////////////////////////////////////////////
double HTCMaterialStatus :: giveValue(string code) const {
    if ( code.compare("alpha_c") == 0 ) {
        return alphac;
    } else if ( code.compare("alpha_s") == 0 ) {
        return alphas;
    } else {
        return TrsprtMaterialStatus :: giveValue(code);
    }
}

//////////////////////////////////////////////////////////
void HTCMaterialStatus :: init() {
    TrsprtMaterialStatus :: init();
    HTCMaterial * htc = static_cast< HTCMaterial * > (mat);
    alphac = temp_alphac = htc->giveInitAlphac();
    alphas = temp_alphas = htc->giveInitAlphas();
    updateMaterialParameters(-1);
}

//////////////////////////////////////////////////////////
void HTCMaterialStatus :: update() {
    TrsprtMaterialStatus :: update();

    alphac = temp_alphac;
    alphas = temp_alphas;
}

//////////////////////////////////////////////////////////
Matrix HTCMaterialStatus :: giveStiffnessTensor(string type, unsigned dim) const {
    (void) type;
    HTCMaterial * htc = static_cast< HTCMaterial * > (mat);
    Matrix P = htc->givePermeabilityTensor();
    double kappa = htc->giveKappa();

    Matrix s(2*dim,2*dim);
    for(unsigned i=0; i<dim; i++){
        for(unsigned j=0; j<dim; j++){
            s[i][j] = -Dh*P[i][j];
            s[i+dim][j+dim] = -kappa*P[i][j];
        }
    }
    return s;    
}

//////////////////////////////////////////////////////////
void HTCMaterialStatus :: updateMaterialParameters(double timeStep) {

    HTCMaterial * htc = static_cast< HTCMaterial * >(mat);
    if(timeStep<0) timeStep=0;

    double betah = 1./( 1.+pow( htc->giveA() - htc->giveA()*h ,htc->giveB()) );
    double Ac, As, alphacdot, alphasdot, old_temp_alphac, old_temp_alphas;

    //central difference method
    temp_alphac = alphac;    
    double midalphac;
    double error = 1e9;
    unsigned it = 0;
    while ( error > 1e-12 && it<100){
        old_temp_alphac = temp_alphac;
        midalphac = (temp_alphac+alphac)/2.;
        Ac = htc->giveAc1() * ( htc->giveAc2()/htc->giveAlphacinf() + midalphac) * (htc->giveAlphacinf() - midalphac) * exp( -htc->giveEtac()*midalphac/htc->giveAlphacinf());  
        alphacdot = Ac*betah*exp(-htc->giveEacR()/T); 
        temp_alphac = alphac +  alphacdot*timeStep;
        error = abs(old_temp_alphac - temp_alphac);
        it ++;
    }
    temp_alphac = max(min(temp_alphac,1.),0.);
    midalphac = temp_alphac; //use value a the end of the time step

    //central difference method
    temp_alphas = alphas;    
    error = 1e9;
    it = 0;
    double midalphas;
    while ( error > 1e-12 && it<100){
        old_temp_alphas = temp_alphas;
        midalphas = (temp_alphas+alphas)/2.;
        As = htc->giveAs1() * ( htc->giveAs2()/htc->giveAlphasinf() + midalphas) * (htc->giveAlphasinf() - midalphas) * exp( -htc->giveEtas()*midalphas/htc->giveAlphasinf());  
        alphasdot = As*exp(-htc->giveEasR()/T); 
        temp_alphas = alphas +  alphasdot*timeStep;
        error = abs(old_temp_alphas - temp_alphas);
        it ++;
    }
    temp_alphas = max(min(temp_alphas,1.),0.);
    midalphas = temp_alphas; //use value a the end of the time step

    double psi = exp ( htc->giveEadR()/htc->giveT0() - htc->giveEadR()/T );
    double G1 = htc->giveKcvg()*htc->giveC()*midalphac + htc->giveKsvg()*htc->giveS()*midalphas;
    double K1 = (htc->giveW0() - 0.188*midalphac*htc->giveC() + 0.22*midalphas*htc->giveS() - G1*(1.-exp(-10*(htc->giveG1()*htc->giveAlphacinf()-midalphac)) ) )/ (exp( 10*(htc->giveG1()*htc->giveAlphacinf()-midalphac)) -1.);

    double G1mult = 1. - 1./exp(10.*(htc->giveG1()*htc->giveAlphacinf()-midalphac)*h);
    double K1mult = exp(10.*(htc->giveG1()*htc->giveAlphacinf()-midalphac)*h) -1.;
    //double we = G1*G1mult + K1*K1mult;
    //double wn = htc->giveKappac() * midalphac * htc->giveC();

    double dG1_dac = htc->giveKcvg()*htc->giveC();
    double dG1_das = htc->giveKsvg()*htc->giveS();
    double dK1_dac = (
        (exp( 10*(htc->giveG1()*htc->giveAlphacinf()-midalphac)) -1.) * (-0.188*htc->giveC() - dG1_dac*(1.-exp(-10*(htc->giveG1()*htc->giveAlphacinf()-midalphac)) )  - G1*(-10*exp(-10*(htc->giveG1()*htc->giveAlphacinf()-midalphac)) ) )
        -(-10.*exp( 10*(htc->giveG1()*htc->giveAlphacinf()-midalphac)) -1.) * (htc->giveW0() - 0.188*midalphac*htc->giveC() + 0.22*midalphas*htc->giveS() - G1*(1.-exp(-10*(htc->giveG1()*htc->giveAlphacinf()-midalphac)) ) )
        )
        / pow(exp( 10*(htc->giveG1()*htc->giveAlphacinf()-midalphac)) -1.,2);
    double dK1_das = (0.22*htc->giveS()  - dG1_das*(1.-exp(-10*(htc->giveG1()*htc->giveAlphacinf()-midalphac)) ))/ (exp( 10*(htc->giveG1()*htc->giveAlphacinf()-midalphac)) -1.);
    double dG1mult_dac = -10.* h *exp(-10.*(htc->giveG1()*htc->giveAlphacinf()-midalphac)*h);
    double dK1mult_dac = -10.* h *exp(10.*(htc->giveG1()*htc->giveAlphacinf()-midalphac)*h);

    double dwe_dac = dG1_dac * G1mult + dK1_dac * K1mult + G1 * dG1mult_dac + K1 * dK1mult_dac;
    double dwe_das = dG1_das * G1mult + dK1_das * K1mult;
    double dwn_dac = htc->giveKappac()*htc->giveC();

    Dh = psi*htc->giveD1()/ (1.+(htc->giveD1()/htc->giveD0()-1.)*pow(1.-h,htc->giveN()));
    qh = dwe_dac* alphacdot + dwe_das * alphasdot + dwn_dac*alphacdot;
    qT = alphacdot*htc->giveC()*htc->giveQcinf() + alphasdot*htc->giveS()*htc->giveQsinf();
    dwe_dh = G1*( 10.*(htc->giveG1()*htc->giveAlphacinf()-midalphac)) * exp(-10.*(htc->giveG1()*htc->giveAlphacinf()-midalphac)*h) + K1*(10.*(htc->giveG1()*htc->giveAlphacinf()-midalphac))*exp(10.*(htc->giveG1()*htc->giveAlphacinf()-midalphac)*h);   
}

//////////////////////////////////////////////////////////
Vector HTCMaterialStatus :: giveStress(const Vector &strain, double timeStep) {
    updateMaterialParameters(timeStep);
    return HTCMaterialStatus :: giveStressWithFrozenIntVars(strain, timeStep);
}

//////////////////////////////////////////////////////////
Vector HTCMaterialStatus :: giveStressWithFrozenIntVars(const Vector &strain, double timeStep) {
    temp_strain = strain;

    Vector hstrain(3);
    Vector tstrain(3);
    for(unsigned i=0; i<3; i++) {
        hstrain[i] = strain[i];
        tstrain[i] = strain[3+i];
    }

    HTCMaterial * htc = static_cast< HTCMaterial * > (mat);
    Matrix P = htc->givePermeabilityTensor();

    Vector hstress = -Dh*matrix_vector_multiply(P,hstrain);    
    Vector tstress = -htc->giveKappa()*matrix_vector_multiply(P,tstrain);    
    temp_stress.resize(6);
    for(unsigned i=0; i<3; i++) {
        temp_stress[i] = hstress[i];
        temp_stress[3+i] = tstress[i];
    }

    return temp_stress;
}

//////////////////////////////////////////////////////////
Vector HTCMaterialStatus :: giveInternalSource() const{
    Vector ints(2);
    ints[0] = -qh;
    ints[1] = qT;
    return ints;
}

//////////////////////////////////////////////////////////
Matrix HTCMaterialStatus :: giveDampingTensor() const{
    HTCMaterial * htc = static_cast< HTCMaterial * > (mat);
    Matrix S(2,2);
    S[0][0] = -dwe_dh;
    S[1][1] = -htc->giveRho()*htc->giveCt();
    return S;
}
//////////////////////////////////////////////////////////
Matrix HTCMaterialStatus :: giveMassTensor() const{
    return Matrix(2,2);
}

//////////////////////////////////////////////////////////
void HTCMaterialStatus :: setParameterValue(string code, double value){
    if (code.compare("humidity")==0) {
        h = value;
    } else if  (code.compare("temperature")==0) {
        T = value;
    } else TrsprtMaterialStatus :: setParameterValue(code, value);
}


//////////////////////////////////////////////////////////
// HTC MATERIAL

//////////////////////////////////////////////////////////
void HTCMaterial :: readFromLine(istringstream &iss) {
    string param;
    
    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("kappa") == 0 ) {            
            iss >> kappa;
        } else if ( param.compare("D1") == 0 ) {            
            iss >> D1;
        } else if ( param.compare("rho") == 0 ) {            
            iss >> rho;
        } else if ( param.compare("c") == 0 ) {            
            iss >> c;
        } else if ( param.compare("ct") == 0 ) {            
            iss >> ct;
        } else if ( param.compare("Qcinf") == 0 ) {            
            iss >> Qcinf;
        } else if ( param.compare("Qsinf") == 0 ) {            
            iss >> Qsinf;
        } else if ( param.compare("s") == 0 ) {            
            iss >> s;
        } else if ( param.compare("EacR") == 0 ) {            
            iss >> EacR;
        } else if ( param.compare("EasR") == 0 ) {            
            iss >> EasR;
        } else if ( param.compare("Ac1") == 0 ) {            
            iss >> Ac1;
        } else if ( param.compare("Ac2") == 0 ) {            
            iss >> Ac2;
        } else if ( param.compare("As1") == 0 ) {            
            iss >> As1;
        } else if ( param.compare("As2") == 0 ) {            
            iss >> As2;
        } else if ( param.compare("alphacinf") == 0 ) {            
            iss >> alphacinf;
        } else if ( param.compare("alphasinf") == 0 ) {            
            iss >> alphasinf;
        } else if ( param.compare("a") == 0 ) {            
            iss >> a;
        } else if ( param.compare("b") == 0 ) {            
            iss >> b;
        } else if ( param.compare("etas") == 0 ) {            
            iss >> etas;
        } else if ( param.compare("etac") == 0 ) {            
            iss >> etac;
        } else if ( param.compare("kcvg") == 0 ) {            
            iss >> kcvg;
        } else if ( param.compare("ksvg") == 0 ) {            
            iss >> ksvg;
        } else if ( param.compare("g1") == 0 ) {            
            iss >> g1;
        } else if ( param.compare("kappac") == 0 ) {            
            iss >> kappac;
        } else if ( param.compare("D0") == 0 ) {            
            iss >> D0;
        } else if ( param.compare("EadR") == 0 ) {            
            iss >> EadR;
        } else if ( param.compare("T0") == 0 ) {            
            iss >> T0;
        } else if ( param.compare("init_alphas") == 0 ) {            
            iss >> init_alphas;
        } else if ( param.compare("init_alphac") == 0 ) {            
            iss >> init_alphac;
        } else if ( param.compare("w0") == 0 ) {            
            iss >> w0;
        } else if ( param.compare("n") == 0 ) {            
            iss >> n;
        }
    }
};

//////////////////////////////////////////////////////////
MaterialStatus *HTCMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    HTCMaterialStatus *newStatus = new HTCMaterialStatus(this, e, ipnum); //needs to be deleted manually
    return newStatus;
};


//////////////////////////////////////////////////////////
void HTCMaterial :: init() {
    permeabilityTensor  = Matrix(3,3);
    
    permeabilityTensor[0][0] = permeabilityTensor[1][1] = permeabilityTensor[2][2] = 0.90020548;
};


