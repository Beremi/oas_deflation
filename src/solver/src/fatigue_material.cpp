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
}

//////////////////////////////////////////////////////////
double FatigueShearMaterialStatus :: giveValue(string code) const {
    if ( code.compare("damage") == 0 ) {
        return damageShear;
    } else if ( code.compare("SlipPi") == 0 ) {
        return sPi;
    } else if ( code.compare("slip") == 0 ) {
        return slip;
    } else {
        return DisMechMaterialStatus :: giveValue(code);
    }
}

//////////////////////////////////////////////////////////
void FatigueShearMaterialStatus :: init() {
    damageShear = temp_damageShear = 0;
    sPi = temp_sPi = 0;
    alphaKin = temp_alphaKin = 0;
    zIso = temp_zIso = 0;
    slip = temp_slip = 0;
}

//////////////////////////////////////////////////////////
Vector FatigueShearMaterialStatus :: giveStress(const Vector &strain) {

  Vector stiff = giveElasticNormalShearStiffness();
  Vector stress( strain.size());

  stress [ 0 ] = stiff [ 0 ] * strain [ 0 ]; //normal stress

  FatigueShearMaterial *m = static_cast< FatigueShearMaterial * >( mat );

  //kill element when excessive tension occur
  if (m->giveTauBar() - m->giveM() * stress[0] <=0){
    for(unsigned i=1; i<stress.size(); i++)  stress[ i ] = 0;
    temp_damageShear = 1;
    tang_stiff = 0;
    return stress;
  }


  double f_trial;

  if ( strain.size() == 2 ) {   // 2D
      //temp_slip = abs(strain [ 1 ]) * mult;  // 2D
      temp_slip = strain [ 1 ];
  } else {  // 3D
      temp_slip = sqrt(pow(strain [ 1 ], 2) + pow(strain [ 2 ], 2) ); //NOT CORRECT
  }

  //compute trials
  double tauTildaPiTrial = stiff [1] * (temp_slip - sPi);
  f_trial = abs(tauTildaPiTrial - m->giveGamma() * alphaKin) - (m->giveKin() * zIso) - (m->giveTauBar() - (m->giveM() * stress [ 0 ]));

  if (f_trial <= 0){
    // internal variables unchanged
    temp_zIso = zIso;
    temp_alphaKin = alphaKin;
    temp_damageShear = damageShear;
    temp_sPi = sPi;
    stress [ 1 ] = (1 - damageShear) * tauTildaPiTrial; //shear stress
  } else {
    // inelastic step, update variables
    double dLambda, Ynext, part1;
    int sgn1;

    dLambda = f_trial / ((stiff[1]/(1 - damageShear)) + m->giveGamma() + m->giveKin());
    double h = tauTildaPiTrial - m->giveGamma() * temp_alphaKin;
    sgn1 = sgn(h);
    temp_sPi = sPi + dLambda * sgn1 / (1 - damageShear);

    Ynext = 0.5 * stiff[1] * pow(temp_slip - temp_sPi, 2);

    part1 = pow(1 - damageShear, m->giveC()) * (m->giveTauBar()/(m->giveTauBar() - m->giveM() * stress[0])) * pow(Ynext / m->giveS(), m->giveR());
    temp_damageShear = fmax(1e-10,fmin(1-1e-10, damageShear + dLambda * part1)); //limited by <0 1>

    temp_zIso = zIso + dLambda;
    temp_alphaKin = alphaKin + dLambda * sgn1;

    stress [ 1 ] = (1 - temp_damageShear) * stiff[ 1 ] * (temp_slip - temp_sPi);   //shear stress

    // calculate algorithmic (tangent) shear stifness
    //computed here only for convenience
    double  partA, partB, partC;
    partA = (1 - temp_damageShear) * stiff [ 1 ];
    partB = ((1 - temp_damageShear) * pow(stiff [ 1 ], 2)) / (stiff [ 1 ] + (m->giveGamma() + m->giveKin()) * (1 - temp_damageShear));
    partC = (pow(stiff [ 1 ], 2) * (temp_slip - temp_sPi)  * part1 * sgn1) / ((stiff [ 1 ] / (1 - temp_damageShear)) + m->giveGamma() + m->giveKin());
    tang_stiff = fmax(0, partA - partB - partC);
  }

  return stress;
}

//////////////////////////////////////////////////////////
void FatigueShearMaterialStatus :: update() {
  damageShear = temp_damageShear;
  sPi = temp_sPi;
  alphaKin = temp_alphaKin;
  zIso = temp_zIso;
  slip = temp_slip;
}

//////////////////////////////////////////////////////////
void FatigueShearMaterialStatus :: print() const {
  std::cout << "damageShear = " << damageShear << '\n';
  std::cout << "temp_damageShear = " << temp_damageShear << '\n';
  std::cout << "sPi = " << sPi << '\n';
  std::cout << "temp_sPi = " << temp_sPi << '\n';
  std::cout << "alphaKin = " << alphaKin << '\n';
  std::cout << "temp_alphaKin = " << temp_alphaKin << '\n';
  std::cout << "zIso = " << zIso << '\n';
  std::cout << "temp_zIso = " << temp_zIso << '\n';
  std::cout << "slip = " << slip << '\n';
  std::cout << "temp_slip = " << temp_slip << '\n';
}


//////////////////////////////////////////////////////////
Vector FatigueShearMaterialStatus :: giveNormalShearStiffness(string type) const {
    Vector stiff = giveElasticNormalShearStiffness();
    if (type.compare("elastic")==0) return stiff;
    else if (type.compare("secant")==0){  //not implemented, used unloading
        stiff[1] *= (1-temp_damageShear); // normal stiffness stays elastic
        return  stiff;
    }
    else if (type.compare("unloading")==0){
        stiff[1] *= (1-temp_damageShear); // normal stiffness stays elastic
        return  stiff;
    }
    else if (type.compare("tangent")==0){
        stiff[1] *= tang_stiff; // normal stiffness stays elastic
        return  stiff;
    }
    else{
        cerr << "Error: FatigueShearMaterialStatus does not provide '"<< type << "' stiffness";
        exit(1);
    };
}

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
