#include "fatigue_material.h"
#include "element.h"


template <typename T> int sgn(T &val) {
    // NOTE this returns 1 for val = 0 (this is an intention, do not repair it!!)
    return (T(0) <= val) - (val < T(0));
}


//////////////////////////////////////////////////////////
// MATERIAL ACCORDING TO ALGORITHM FROM https://doi.org/10.1016/j.ijfatigue.2018.04.020
// TODO:
// NORMAL STIFNESS IS KEPT ELASTIC, DAMAGE INFLUENCES ONLY SHEAR DIRECTION
//////////////////////////////////////////////////////////
// FATIGUE SHEAR MATERIAL STATUS

FatigueShearMaterialStatus :: FatigueShearMaterialStatus(FatigueShearMaterial *m, Element *e) : DisMechMaterialStatus(m, e) {
    name = "Fatigue Shear mat. status";
    RAND_H = 1.0;
}

//////////////////////////////////////////////////////////
double FatigueShearMaterialStatus :: giveValue(string code) const {
    if ( code.compare("damage") == 0 ) {
        return temp_damageShear;
    } else if ( code.compare("cumSlip") == 0 ) {
        return temp_sPi;
    } else if ( code.compare("slip") == 0 ) {
        return slip_cur;
    } else if ( code.compare("multip") == 0 ) {
        return temp_stiffMultip;
    } else {
        return DisMechMaterialStatus :: giveValue(code);
    }
}

//////////////////////////////////////////////////////////
void FatigueShearMaterialStatus :: init() {
    maxEpsT = temp_maxEpsT = 0;
    damageShear = temp_damageShear = 0;
    sPi = temp_sPi = 0;
    stiffMultip = temp_stiffMultip = 1;
    // tauPi = temp_tauPi = 0;
    tauPiTrial = tauTildaPiTrial = 0;
    alphaKin = temp_alphaKin = 0;
    zIso = temp_zIso = 0;
    EAlg = temp_EAlg = 0;
}


double FatigueShearMaterialStatus :: computeTrials(const double &stressN, const FatigueShearMaterial *m){
  tauTildaPiTrial = m->giveEb() * (slip_cur - sPi);
  tauPiTrial = (1 - damageShear) * tauTildaPiTrial;
  return abs(tauTildaPiTrial - m->giveGamma() * alphaKin) - (m->giveKin() * zIso) - m->giveTauBar() + (m->giveM() * stressN);
}

void FatigueShearMaterialStatus :: computeShearStifness(const Vector &strain){
  double f_trial, stressN;
  double mult = 1;
  RigidBodyContact *rbc = static_cast< RigidBodyContact * >( element );
  FatigueShearMaterial *m = static_cast< FatigueShearMaterial * >( mat );
  if ( strain.size() == 2 ) {
      // 1e3 is added to use the same parameters
      slip_cur = strain [ 1 ] * mult;  // 2D
  } else {
      // NOTE the sign of shear matters (mainly when crossing the zero value and continuing further into negative)
      slip_cur = sqrt(pow(strain [ 1 ], 2) + pow(strain [ 2 ], 2) ) * sgn(strain [ 1 ]) * sgn(strain [ 2 ]) * mult;  // 3D
  }
  // std::cout << "slip = " << slip << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<," << '\n';
  stressN = strain[0] * giveSecantNormalShearStiffness()[0];
  f_trial = computeTrials(stressN, m);
  if (f_trial <= 0){
    // std::cout << "elastic step" << '\n';
    // elastic step, variables stay the same, only tauPi is updated
    // temp_tauPi = tauPiTrial;
    temp_stiffMultip = stiffMultip;
  } else {
    // std::cout << "nonelastic step" << '\n';
    // nonelastic step, update variables
    double dLambda, Ynext, part1;
    int sgn1;

    dLambda = f_trial / ((m->giveEb()/(1 - damageShear)) + m->giveGamma() + m->giveKin());
    double hlep = tauTildaPiTrial - m->giveGamma() * temp_alphaKin;
    sgn1 = sgn(hlep);
    temp_sPi = sPi + dLambda * sgn1 / (1 - damageShear);

    Ynext = 0.5 * m->giveEb() * pow(slip_cur - temp_sPi, 2);

    part1 = pow(1 - damageShear, m->giveC()) * (m->giveTauBar()/(m->giveTauBar() - m->giveM() * stressN)) * pow(Ynext / m->giveS(), m->giveR());
    temp_damageShear = fmin(1-1e-10, damageShear + dLambda * part1);

    // temp_tauPi = (1 - temp_damageShear) * m->giveEb() * (slip_cur - temp_sPi);

    temp_zIso = zIso + dLambda;

    temp_alphaKin = alphaKin + dLambda * sgn1;

    if (abs(slip_cur - slip_prev) < 1e-10){
      temp_stiffMultip = stiffMultip;
    } else if (slip_cur == 0) {
      temp_stiffMultip = (1 - temp_damageShear) * (rbc->giveArea() / rbc->giveLength()) * mult;
    } else {
      temp_stiffMultip = (1 - temp_damageShear) * (1 - temp_sPi / slip_cur)  * (rbc->giveArea() / rbc->giveLength()) * mult; // normal stiffness stays elastic
    }

    // // calculate algorithmic (tangent) shear stifness
    double  partA, partB, partC;

    partA = (1 - temp_damageShear) * m->giveEb();

    partB = ((1 - temp_damageShear) * pow(m->giveEb(), 2)) / (m->giveEb() + (m->giveGamma() + m->giveKin()) * (1 - temp_damageShear));

    partC = (pow(m->giveEb(), 2) * (slip_cur - temp_sPi)  * part1 * sgn1) / ((m->giveEb() / (1 - temp_damageShear)) + m->giveGamma() + m->giveKin());

    temp_EAlg = fmax(0, partA - partB - partC);
    // std::cout << "Eb = " << m->giveEb() << ", Ealg = " << temp_EAlg << '\n';
    // std::cout << "damage = " << temp_damageShear << ", slip = " << slip << ", sPi = " << sPi << ", stiff = " << temp_tauPi/slip << '\n';

    // damAlg = (1 - temp_EAlg / m->giveEb());
    // std::cout << "Eb = " << m->giveEb() << ", tauBar = " << m->giveTauBar() << ", Kin = " << m->giveKin() << ", gamma = " << m->giveGamma() << ", S = " << m->giveS() << ", c = " << m->giveC() << ", r = " << m->giveR() << ", m = " << m->giveM() << '\n';
  }

}

//////////////////////////////////////////////////////////
void FatigueShearMaterialStatus :: update() {
  maxEpsT = temp_maxEpsT;
  damageShear = temp_damageShear;
  sPi = temp_sPi;
  stiffMultip = temp_stiffMultip;
  // tauPi = temp_tauPi;
  alphaKin = temp_alphaKin;
  zIso = temp_zIso;
  slip_prev = slip_cur;
  EAlg = temp_EAlg;
}

void FatigueShearMaterialStatus :: print() const {
  std::cout << "damageShear = " << damageShear << '\n';
  std::cout << "temp_damageShear = " << temp_damageShear << '\n';
  std::cout << "stiffMultip = " << stiffMultip << '\n';
  std::cout << "temp_stiffMultip = " << temp_stiffMultip << '\n';
  std::cout << "sPi = " << sPi << '\n';
  std::cout << "temp_sPi = " << temp_sPi << '\n';
  // std::cout << "tauPi = " << tauPi << '\n';
  // std::cout << "temp_tauPi = " << temp_tauPi << '\n';
  std::cout << "tauPiTrial = " << tauPiTrial << '\n';
  std::cout << "tauTildaPiTrial = " << tauTildaPiTrial << '\n';
  std::cout << "alphaKin = " << alphaKin << '\n';
  std::cout << "temp_alphaKin = " << temp_alphaKin << '\n';
  std::cout << "zIso = " << zIso << '\n';
  std::cout << "temp_zIso = " << temp_zIso << '\n';
  std::cout << "temp_EAlg = " << temp_EAlg << '\n';
  std::cout << "EAlg = " << EAlg << '\n';
}

//////////////////////////////////////////////////////////
Vector FatigueShearMaterialStatus :: giveSecantNormalShearStiffness() const {
    Vector stiff = giveElasticNormalShearStiffness();
    stiff[1] *= temp_stiffMultip; // normal stiffness stays elastic
    return stiff;
}

//////////////////////////////////////////////////////////
Vector FatigueShearMaterialStatus :: giveSecantNormalShearStiffness(const Vector &strain) {
  Vector stiff = giveElasticNormalShearStiffness();
  computeShearStifness(strain);
  stiff[1] *= temp_stiffMultip; // normal stiffness stays elastic
  return stiff;
}

Vector FatigueShearMaterialStatus :: giveTangentNormalShearStiffness() const{
  Vector stiff = giveElasticNormalShearStiffness();
  stiff[1] = temp_EAlg; // normal stiffness stays elastic
  return stiff;
}

Vector FatigueShearMaterialStatus :: giveTangentNormalShearStiffness(const Vector &strain){
  Vector stiff = giveElasticNormalShearStiffness();
  computeShearStifness(strain);
  stiff[1] = temp_EAlg;  // normal stiffness stays elastic
  return stiff;
}

//////////////////////////////////////////////////////////
Vector FatigueShearMaterialStatus :: giveStress(const Vector &strain) {
  computeShearStifness(strain);
  // this->print();
  Vector stiff = giveElasticNormalShearStiffness();
  stiff[1] *= temp_stiffMultip; // normal stiffness stays elastic
  Vector stress(strain.size() );
  stress [ 0 ] = stiff [ 0 ] * strain [ 0 ];
  for ( unsigned i = 1; i < strain.size(); i++ ) {
      stress [ i ] = stiff [ 1 ] * strain [ i ];
  }
  // std::cout << "damage = " << temp_damageShear << ", slip = " << slip_cur << ", slip_prev = " << slip_prev << ", sPi = " << sPi << ", E_tru = " << stiff[1] << '\n';
  return stress;
}
//
//
//////////////////////////////////////////////////////////
// FATIGUE SHEAR MATERIAL

//////////////////////////////////////////////////////////
void FatigueShearMaterial :: readFromLine(istringstream &iss) {
    DisMechMaterial :: readFromLine(iss); //read elastic parameters

    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream

    string param;
    bool btau, bkin, bgam, bs, bc, br, bm;
    btau = bkin = bgam = bs = bc = br = bm = false;
    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("tauBar") == 0 ) {
            btau = true;
            iss >> tauBar;
        } else if ( param.compare("Kin") == 0 )    {
            bkin = true;
            iss >> Kin;
        } else if ( param.compare("gamma") == 0 )    {
            bgam = true;
            iss >> gamma;
        } else if ( param.compare("S") == 0 )    {
            bs = true;
            iss >> S;
        } else if ( param.compare("c") == 0 )    {
            bc = true;
            iss >> c;
        } else if ( param.compare("r") == 0 )    {
            br = true;
            iss >> r;
        } else if ( param.compare("m") == 0 )    {
            bm = true;
            iss >> m;
        }
    }
    if ( !btau ) {
        cerr << name << ": material parameter 'tauBar' was not specified" << endl;
        exit(0);
    }
    ;
    if ( !bkin ) {
        cerr << name << ": material parameter 'Kin' was not specified" << endl;
        exit(0);
    }
    if ( !bgam ) {
        cerr << name << ": material parameter 'gamma' was not specified" << endl;
        exit(0);
    }
    if ( !bs ) {
        cerr << name << ": material parameter 'S' was not specified" << endl;
        exit(0);
    }
    if ( !bc ) {
        cout << name << ": material parameter 'c' was not specified, taking c = 1.0" << endl;
        c = 1.0;
    }
    if ( !br ) {
        cout << name << ": material parameter 'r' was not specified, taking r = 1.0" << endl;
        r = 1.0;
    }
    if ( !bm ) {
        cout << name << ": material parameter 'm' was not specified, taking m = 1.0" << endl;
        m = 1.0;
    }
    ;
};

//////////////////////////////////////////////////////////
MaterialStatus *FatigueShearMaterial :: giveNewMaterialStatus(Element *e) {
    FatigueShearMaterialStatus *newStatus = new FatigueShearMaterialStatus(this, e); //needs to be deleted manually
    return newStatus;
};


//////////////////////////////////////////////////////////
void FatigueShearMaterial :: init() {

};
