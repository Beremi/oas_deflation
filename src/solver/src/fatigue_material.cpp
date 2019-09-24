#include "fatigue_material.h"
#include "element.h"


template <typename T> int sgn(T &val) {
    // NOTE this returns 1 for val = 0 (this is an intention, do not repair it!!)
    return (T(0) <= val) - (val < T(0));
}

//////////////////////////////////////////////////////////
// FATIGUE SHEAR MATERIAL STATUS

FatigueShearMaterialStatus :: FatigueShearMaterialStatus(FatigueShearMaterial *m, Element *e) : DisMechMaterialStatus(m, e) {
    name = "Fatigue Shear mat. status";
}

//////////////////////////////////////////////////////////
double FatigueShearMaterialStatus :: giveValue(string code) const {
    if ( (code.compare("damage") == 0) ||  (code.compare("damageT") == 0) ) {
        return damageShear;
    } else if ( (code.compare("strainPL") == 0) ||  (code.compare("strainPLT") == 0) ) {
        return sPi.norm();
    } else if ( (code.compare("strainT") == 0) ||  (code.compare("strain") == 0) ) {
        return slip.norm();
    } else if ( (code.compare("strainTY") == 0)) {
        return slip.getY();
    } else if ( (code.compare("strainTZ") == 0)) {
        return slip.getZ();
    } else if ( (code.compare("strainPLTY") == 0)) {
        return sPi.getY();
    } else if ( (code.compare("strainPLTZ") == 0)) {
        return sPi.getZ();
    } else if ( (code.compare("stressTY") == 0)) {
        return stressT.getY();
    } else if ( (code.compare("stressTZ") == 0)) {
        return stressT.getZ();
    } else {
        return DisMechMaterialStatus :: giveValue(code);
    }
}

//////////////////////////////////////////////////////////
void FatigueShearMaterialStatus :: init() {
    damageShear = temp_damageShear = 0;
    zIso = temp_zIso = 0;
    sPi = temp_sPi = Point();
    alphaKin = temp_alphaKin = Point();
    slip = temp_slip = Point();
}

//////////////////////////////////////////////////////////
Vector FatigueShearMaterialStatus :: giveStress(const Vector &strain) {
  // TENSORIAL FORM OF CONST LAW ACCORDING TO FRAMCOS PAPER BY ABEDULGADER BAKTHER et al doi.org/10.21012/FC10.233196
  ////////////////////////////////////////////////////////

  Vector stiff = giveElasticNormalShearStiffness();
  Vector stress( strain.size());
  // Point stressT;  // declared in class

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
  double x = 0;
  double y = 0;
  double z = 0;

  for (unsigned i = 1; i < strain.size(); i++){
    if (i == 1) y = strain[ i ];
    else if (i == 2) z = strain[ i ];
    else {
      std::cerr << "should never get here, exit" << '\n';
      exit(1);
    }
  }
  temp_slip = Point(x, y, z);

  //compute trials
  Point tauTildaPiTrial = (temp_slip - sPi) * stiff [1];
  f_trial = (tauTildaPiTrial - alphaKin * m->giveGamma()).norm() - (m->giveKin() * zIso) - (m->giveTauBar() - (m->giveM() * stress [ 0 ]));

  if (f_trial <= 0){
    // internal variables unchanged
    // it is necessary to asign them to temp, because temp values could have been changed in the previous iterration
    temp_zIso = zIso;
    temp_alphaKin = alphaKin;
    temp_damageShear = damageShear;
    temp_sPi = sPi;
    stressT =  tauTildaPiTrial * (1 - damageShear); //shear stress
  } else {
    // inelastic step, update variables
    double dLambda, Ynext, part1;
    Point sgn1;

    dLambda = f_trial / ((stiff[1]/(1 - damageShear)) + m->giveGamma() + m->giveKin());
    Point h = tauTildaPiTrial - temp_alphaKin * m->giveGamma();
    sgn1 = h / h.norm();
    temp_sPi = sPi + sgn1 * dLambda / (1 - damageShear);

    Ynext = 0.5 * stiff[1] * (temp_slip - temp_sPi).sqNorm(); // sqNorm = self dot product

    part1 = pow(1 - damageShear, m->giveC()) * (m->giveTauBar()/(m->giveTauBar() - m->giveM() * stress[0])) * pow(Ynext / m->giveS(), m->giveR());
    temp_damageShear = fmax(1e-10,fmin(1-1e-10, damageShear + dLambda * part1)); //limited by <0 1>

    temp_zIso = zIso + dLambda;
    temp_alphaKin = alphaKin + sgn1 * dLambda;

    stressT = (temp_slip - temp_sPi) * (1 - temp_damageShear) * stiff[ 1 ];   //shear stress

    // calculate algorithmic (tangent) shear stifness
    //computed here only for convenience
    double  partA, partB, partC;
    partA = (1 - temp_damageShear) * stiff [ 1 ];
    partB = ((1 - temp_damageShear) * pow(stiff [ 1 ], 2)) / (stiff [ 1 ] + (m->giveGamma() + m->giveKin()) * (1 - temp_damageShear));
    // in partC, norms are used to obtain scalar, not sure if it is correct
    partC = (pow(stiff [ 1 ], 2) * (temp_slip - temp_sPi).norm() * part1 * sgn1.norm()) / ((stiff [ 1 ] / (1 - temp_damageShear)) + m->giveGamma() + m->giveKin());
    tang_stiff = fmax(0, partA - partB - partC);
  }

  for (unsigned i = 1; i < strain.size(); i++){
    if (i == 1) stress[ i ] = stressT.getY();
    else if (i == 2) stress[ i ] = stressT.getZ();
    else {
      std::cerr << "should never get here, exit" << '\n';
      exit(1);
    }
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
  std::cout << "sPi = " << sPi.norm() << '\n';
  std::cout << "temp_sPi = " << temp_sPi.norm() << '\n';
  std::cout << "alphaKin = " << alphaKin.norm() << '\n';
  std::cout << "temp_alphaKin = " << temp_alphaKin.norm() << '\n';
  std::cout << "zIso = " << zIso << '\n';
  std::cout << "temp_zIso = " << temp_zIso << '\n';
  std::cout << "slip = " << slip.norm() << '\n';
  std::cout << "temp_slip = " << temp_slip.norm() << '\n';
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
        } else if ( param.compare("a") == 0 )    {
            // pressure sensitivity renamed due to collision with normal loading plasticity parameter
            // for coupled model (cumulative sliding + plasticity damage)
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
        cout << name << ": material parameter 'a' was not specified, taking a = 0.0" << endl;
        m = 0.0;
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


//////////////////////////////////////////////////////////
// DAMAGE PLASTIC MATERIAL STATUS

DamagePlasticMaterialStatus :: DamagePlasticMaterialStatus(DamagePlasticMaterial *m, Element *e) : DisMechMaterialStatus(m, e) {
    name = "Damage Plastic mat. status";
}

//////////////////////////////////////////////////////////
double DamagePlasticMaterialStatus :: giveValue(string code) const {
    if ( (code.compare("damage") == 0) || (code.compare("damageN") == 0) ) {
        return damage;
    } else if ( (code.compare("strainN") == 0) || (code.compare("strain") == 0) ) {
        return epsN;
    } else if ( (code.compare("stressN") == 0) || (code.compare("stress") == 0) ) {
        return stressN;
    } else if ( (code.compare("strainPL") == 0) || (code.compare("strainPLN") == 0) ) {
        return epsNP;
    } else {
        std::cout << "do not have such parameter " << code << '\n';
        return DisMechMaterialStatus :: giveValue(code);
    }
}

//////////////////////////////////////////////////////////
void DamagePlasticMaterialStatus :: init() {
    damage = temp_damage = 0;
    epsN = temp_epsN = 0;
    epsNP = temp_epsNP = 0;
    alphaN = temp_alphaN = 0;
    zN = temp_zN = 0;
    rN = temp_rN = 0;
}

//////////////////////////////////////////////////////////
Vector DamagePlasticMaterialStatus :: giveStress(const Vector &strain) {
  // TODO transition from compression to tension is very simply done here, should be improved
  Vector stress( strain.size() );
  Vector stiff = giveElasticNormalShearStiffness();
  for (size_t i = 1; i < stress.size(); i++){
    stress[ i ] = stiff [ i ] * strain [ i ];
  }
  temp_epsN = strain[ 0 ];

  DamagePlasticMaterial *m = static_cast< DamagePlasticMaterial * >( mat );

  ////////////////////////////////////////////////////////
  double f_trialT, f_trialC;
  double SigmaTilda;
  // distinguish between tension and compression

  int Heaviside;
  if (temp_epsN > 0){
    Heaviside = 1;
    // TODO following should work with plasticity, but does not
    // double Y_next = Heaviside * 0.5 * stiff [ 0 ] * pow( temp_epsN - epsNP, 2 );
    double Y_next = Heaviside * 0.5 * stiff [ 0 ] * pow( temp_epsN, 2 );
    double Y_n0 = 0.5 * stiff [ 0 ] * pow( m->giveElasticLimit(), 2);
    double Rn = ( 1 / m->giveAd() ) * ( - rN / ( 1 + rN ) );

    f_trialT = Y_next - ( Y_n0 + Rn );
    if ( f_trialT <= 0 ){
      temp_damage = damage;
      temp_rN = rN;
    } else {
      // update tensile internal variables
      temp_damage = fmin( 1-1e-10, fmax( 0, 1 - 1 / ( 1 + m->giveAd() * ( Y_next - Y_n0) ) ) );
      temp_rN = -temp_damage;
    }
    stress[ 0 ] = ( 1 - Heaviside * temp_damage ) * stiff [ 0 ] * ( temp_epsN );
    // apply damage also in shear direction
    for ( unsigned i = 1; i < strain.size(); i++){
      stress[ i ] = ( 1 - Heaviside * temp_damage ) * stiff [ i ] * ( strain[ i ] );
    }
    stressN = stress[ 0 ];
    return stress;
  } else {
    Heaviside = 0;
    SigmaTilda = ( 1 - Heaviside * damage ) * stiff [ 0 ] * ( temp_epsN - epsNP );
    double Xn = m->giveGammaN() * alphaN;
    double Zn = m->giveKinN() * zN;
    int pos_iso;
    if ( m->giveYieldStress() + Zn > 0 ){
      pos_iso = 1;
    } else {
      pos_iso = 0;
    }
    f_trialC =  fabs( SigmaTilda - Xn ) - pos_iso * ( m->giveYieldStress() + Zn );
    if ( f_trialC <= 0 ){
      temp_epsNP = epsNP;
      temp_alphaN = alphaN;
      temp_zN = zN;
    } else {
      // update compressive internal variables
      double h = SigmaTilda - Xn;
      int sgn1 = sgn( h );
      ////////////////////////////////////////////////////////////////
      // TODO the following expression is taken from article, but does not work yet
      // double dLambda = ( stiff [ 0 ] * ( temp_epsN - epsN ) * sgn1 / ( stiff [ 0 ] + m->giveKinN()) + m->giveGammaN() *
      //                         ( 1 + m->giveM() * Xn * sgn1 ) );
      ////////////////////////////////////////////////////////////////
      double dLambda = ( f_trialC / ( stiff [ 0 ] + m->giveKinN() + m->giveGammaN() ) );
      temp_alphaN = alphaN + dLambda * (sgn1 + m->giveM() * Xn );
      temp_zN = zN + dLambda;
      temp_epsNP = epsNP + dLambda * sgn1;

    }
    stress[ 0 ] = ( 1 - Heaviside * temp_damage ) * stiff [ 0 ] * ( temp_epsN - temp_epsNP );
    // apply the same in shear direction
    for ( unsigned i = 1; i < strain.size(); i++){
      // stress[ i ] = ( 1 - Heaviside * temp_damage ) * stiff [ i ] * ( strain[ i ] - (strain[ i ]/temp_epsN)*temp_epsNP );
    }
    if (stress[ 0 ] > 0){
      stress[ 0 ] = 0;
    }
    stressN = stress[ 0 ];
    return stress;
  }
}

//////////////////////////////////////////////////////////
void DamagePlasticMaterialStatus :: update() {
  damage = temp_damage;
  epsN = temp_epsN;
  epsNP = temp_epsNP;
  alphaN = temp_alphaN;
  zN = temp_zN;
  rN = temp_rN;
}

//////////////////////////////////////////////////////////
void DamagePlasticMaterialStatus :: print() const {
  std::cout << "damage = " << damage << '\n';
  std::cout << "temp_damage = " << temp_damage << '\n';
  std::cout << "strain = " << epsN << '\n';
  std::cout << "temp_strain = " << temp_epsN << '\n';
  std::cout << "strainP = " << epsNP << '\n';
  std::cout << "temp_strainP = " << temp_epsNP << '\n';
  std::cout << "alphaN = " << alphaN << '\n';
  std::cout << "temp_alphaN = " << temp_alphaN << '\n';
  std::cout << "zN = " << zN << '\n';
  std::cout << "temp_zN = " << temp_zN << '\n';
  std::cout << "rN = " << rN << '\n';
  std::cout << "temp_rN = " << temp_rN << '\n';
}


//////////////////////////////////////////////////////////
Vector DamagePlasticMaterialStatus :: giveNormalShearStiffness(string type) const {
    Vector stiff = giveElasticNormalShearStiffness();
    if (type.compare( "elastic" ) == 0 ){
      return stiff;
    } else if (type.compare( "secant" ) == 0 ){  //not implemented, used unloading
        if (temp_epsN >= 0 ){
          for (unsigned i = 0; i < stiff.size(); i++){
            stiff[ i ] *= ( 1 - temp_damage ); // this is just unloading one
          }
        }
        return  stiff;
    } else if (type.compare( "unloading" ) == 0 ){
      if (temp_epsN >= 0 ){
        stiff[ 1 ] *= ( 1 - temp_damage );
      }
      return  stiff;
    }
    // else if (type.compare("tangent")==0){
    //     stiff[1] *= tang_stiff; // normal stiffness stays elastic
    //     return  stiff;
    // }
    else{
        cerr << "Error: DamagePlasticMaterialStatus does not provide '"<< type << "' stiffness";
        exit(1);
    };
}

//////////////////////////////////////////////////////////
// DAMAGE PLASTIC MATERIAL
//////////////////////////////////////////////////////////
void DamagePlasticMaterial :: readFromLine(istringstream &iss) {
    DisMechMaterial :: readFromLine(iss); //read elastic parameters

    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream

    string param;
    bool bfc, bft, bgam, bkin, bAd, bm;
    bfc = bft = bgam = bkin = bAd = bm = false;
    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("fc") == 0 ) {
            bfc = true;
            iss >> fc;
        } else if ( param.compare("ft") == 0 ) {
            bft = true;
            iss >> ft;
        } else if ( param.compare("KinN") == 0 )    {
            bkin = true;
            iss >> KinN;
        } else if ( param.compare("gammaN") == 0 )    {
            bgam = true;
            iss >> gammaN;
        } else if ( param.compare("Ad") == 0 )    {
            bAd = true;
            iss >> Ad;
        } else if ( param.compare("m") == 0 )    {
            bm = true;
            iss >> m;
        }
    }
    if ( !bfc ) {
        cerr << name << ": material parameter 'fc' was not specified" << endl;
        exit(0);
    }
    if ( !bft ) {
        cerr << name << ": material parameter 'ft' was not specified" << endl;
        exit(0);
    }
    if ( !bkin ) {
        cerr << name << ": material parameter 'KinN' was not specified" << endl;
        exit(0);
    }
    if ( !bgam ) {
        cerr << name << ": material parameter 'gammaN' was not specified" << endl;
        exit(0);
    }
    if ( !bAd ) {
        cerr << name << ": material parameter 'Ad' was not specified" << endl;
        exit(0);
    }
    if ( !bm ) {
        cout << name << ": material parameter 'm' was not specified, taking m = 0.0" << endl;
        m = 0.0;
    }
    ;
};

//////////////////////////////////////////////////////////
MaterialStatus *DamagePlasticMaterial :: giveNewMaterialStatus(Element *e) {
    DamagePlasticMaterialStatus *newStatus = new DamagePlasticMaterialStatus(this, e); //needs to be deleted manually
    return newStatus;
};

//////////////////////////////////////////////////////////
void DamagePlasticMaterial :: init() {

};

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// normal (damagePlastic) and tangential (cumulative slip - fatigue shear) const. laws together

FatigueMaterialStatus :: FatigueMaterialStatus(FatigueMaterial *m, Element *e) : FatigueShearMaterialStatus(m, e), DamagePlasticMaterialStatus(m, e) {
    FatigueShearMaterialStatus :: name = "Fatigue mat. status";
    DamagePlasticMaterialStatus :: name = "Fatigue mat. status";
}

//////////////////////////////////////////////////////////
double FatigueMaterialStatus :: giveValue(string code) const {
    if ( code.compare("damageN") == 0 ) {
        return DamagePlasticMaterialStatus :: giveValue("damageN");
    } else if ( code.compare("strainN") == 0 ) {
        return DamagePlasticMaterialStatus :: giveValue("strainN");
    } else if ( code.compare("stressN") == 0 ) {
        return DamagePlasticMaterialStatus :: giveValue("stressN");
    } else if ( code.compare("strainPLN") == 0 ) {
        return DamagePlasticMaterialStatus :: giveValue("strainPLN");
    } else {
      return FatigueShearMaterialStatus :: giveValue(code);
      // return 0;
    // return DisMechMaterialStatus :: giveValue(code);
    }
}

//////////////////////////////////////////////////////////
void FatigueMaterialStatus :: init() {
    FatigueShearMaterialStatus :: init();
    DamagePlasticMaterialStatus :: init();
}

//////////////////////////////////////////////////////////
Vector FatigueMaterialStatus :: giveStress(const Vector &strain) {
  // TODO transition from compression to tension is very simply done here, should be improved
  Vector stress( strain.size() );
  Vector stiff = FatigueShearMaterialStatus :: giveElasticNormalShearStiffness();

  for (size_t i = 0; i < stress.size(); i++){
    if ( i == 0 ){
      stress[ i ] = DamagePlasticMaterialStatus :: giveStress(strain)[ i ];
    } else {
      stress[ i ] = FatigueShearMaterialStatus :: giveStress(strain)[ i ];
    }
  }

  return stress;
}

//////////////////////////////////////////////////////////
void FatigueMaterialStatus :: update() {
  FatigueShearMaterialStatus :: update();
  DamagePlasticMaterialStatus :: update();
}

//////////////////////////////////////////////////////////
Vector FatigueMaterialStatus :: giveNormalShearStiffness(string type) const {
    Vector stiff = FatigueShearMaterialStatus :: giveElasticNormalShearStiffness();
    if (type.compare( "elastic" ) == 0 ){
      return stiff;
    } else if (type.compare( "secant" ) == 0 || type.compare( "unloading" ) == 0 ){  //not implemented, used unloading
      for (unsigned i = 0; i < stiff.size(); i++){
        if ( i == 0 ){
          stiff[ i ] = DamagePlasticMaterialStatus :: giveNormalShearStiffness(type)[ i ];
        } else {
          stiff[ i ] = FatigueShearMaterialStatus :: giveNormalShearStiffness(type)[ i ];
        }
      }
      return  stiff;
    } else {
      cerr << "Error: FatigueMaterialStatus does not provide '"<< type << "' stiffness";
      exit(1);
    };
}

//////////////////////////////////////////////////////////
// fatigue material in normal and tangential direction
//////////////////////////////////////////////////////////
void FatigueMaterial :: readFromLine(istringstream &iss) {
    FatigueShearMaterial :: readFromLine(iss);
    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream
    DamagePlasticMaterial :: readFromLine(iss);
};

//////////////////////////////////////////////////////////
MaterialStatus *FatigueMaterial :: giveNewMaterialStatus(Element *e) {
    FatigueMaterialStatus *newStatus1 = new FatigueMaterialStatus(this, e); //needs to be deleted manually
    FatigueShearMaterialStatus *newStatus = (FatigueMaterialStatus *) newStatus1;
    return newStatus;
};

//////////////////////////////////////////////////////////
void FatigueMaterial :: init() {

};

//////////////////////////////////////////////////////////
// ALLICHE MATERIAL STATUS
//////////////////////////////////////////////////////////

AllicheMaterialStatus :: AllicheMaterialStatus(AllicheMaterial *m, Element *e) : DisMechMaterialStatus(m, e) {
    name = "Alliche mat. status";
}

//////////////////////////////////////////////////////////
double AllicheMaterialStatus :: giveValue(string code) const {
    if ( (code.compare("damage") == 0) ) {
        // damage is different in each local direction
        return damage.norm();
    } else if ( (code.compare("damageX") == 0) ) {
        return damage.getX();
    } else if ( (code.compare("damageY") == 0) ) {
        return damage.getY();
    } else if ( (code.compare("damageZ") == 0) ) {
        return damage.getZ();
    } else if ( (code.compare("stressX") == 0) ) {
        return sigma.getX();
    } else if ( (code.compare("stressY") == 0) ) {
        return sigma.getY();
    } else if ( (code.compare("stressZ") == 0) ) {
        return sigma.getZ();
    } else {
        return DisMechMaterialStatus :: giveValue(code);
    }
}

//////////////////////////////////////////////////////////
void AllicheMaterialStatus :: init() {
    damage = temp_damage = 0;
    eps = eps_plus = eps_plus_prev = Point();
    // shear_eps_cur = shear_eps_prev = Point();
    sigma = Point();
    Y = Y_plus = Y_prev = Point();
}

//////////////////////////////////////////////////////////
void AllicheMaterialStatus :: calculateDamage(const Vector &strain) {
  AllicheMaterial *m = static_cast< AllicheMaterial * >( mat );

  // std::cout << "strain inside = (" << strain[ 0 ] << ", " << strain[ 1 ] << ", " << strain[ 2 ] << ")" << '\n';

  for (unsigned i = 0; i < strain.size(); i++){
    if ( i == 0 ){
      eps.setX( strain[ i ] );
      if ( strain[ i ] > 0 ){
        eps_plus.setX( strain[ i ] );
      }
    } else if ( i == 1 ){
      eps.setY( strain[ i ] );
      // if ( strain[ i ] > 0 )
        eps_plus.setY( abs(strain[ i ]) );  // shear is considered always positive...?
    } else if ( i == 2 ){
      eps.setZ( strain[ i ] );
      // if ( strain[ i ] > 0 )
        eps_plus.setZ( abs(strain[ i ]) );  // shear is considered always positive...?
    } else {
      std::cerr << __func__ << " should never get here, exiting" << '\n';
    }
  }

  // I stands for identity matrix (coord vector in this case)
  Point I(1, 1, 1);

  Y = eps * (-1) * m->giveG() - eps * eps.sum() * m->giveAlphaDam() -
      Point(pow(eps.getX(), 2), pow(eps.getY(), 2), pow(eps.getZ(), 2)) * 2 * m->giveBetaDam();

  if ( Y.getX() < 0 ) Y.setX(0);  ///< only positive part needed
  if ( Y.getY() < 0 ) Y.setY(0);
  if ( Y.getZ() < 0 ) Y.setZ(0);

  double treshold_func;
  // treshold_func = Y.sqNorm() - Y_prev.sqNorm();
  treshold_func = (Y - Y_prev).sum();


  if ( treshold_func > 0 ){
    double yield_func = ( m->giveG() / sqrt(2) ) * eps_plus.norm() - (m->giveC0() - m->giveC1() * damage.sum());
    temp_damage = damage + eps_plus * (1 / sqrt(2 * eps_plus.sqNorm() ) ) *
                  ( dot(eps_plus, eps_plus - eps_plus_prev) / (m->giveC1() * eps_plus.sum() ) ) *
                  pow( yield_func / m->giveK(), m->giveN() );

    temp_damage.set(fmax(0, fmin(temp_damage.getX(), 1) ),
                    fmax(0, fmin(temp_damage.getY(), 1) ),
                    fmax(0, fmin(temp_damage.getZ(), 1) ) );
  } else {
    temp_damage = damage;
  }

}

Vector AllicheMaterialStatus :: giveStress(const Vector &strain){
  Vector stiff = giveElasticNormalShearStiffness();
  Vector stress( strain.size());
  calculateDamage(strain);

  // NOTE use of damage as Point makes its calclulation more simple, but here it becomes quite messy
  stress[ 0 ] = ( 1 - temp_damage.getX() ) * stiff[ 0 ] * strain[ 0 ];
  stress[ 1 ] = ( 1 - temp_damage.getY() ) * stiff[ 1 ] * strain[ 1 ];;
  if ( stress.size() > 2 ){
    stress[ 2 ] = ( 1 - temp_damage.getZ() ) * stiff[ 2 ] * strain[ 2 ];;
  }

  return stress;
}

//////////////////////////////////////////////////////////
Vector AllicheMaterialStatus :: giveNormalShearStiffness(string type) const {
    Vector stiff = giveElasticNormalShearStiffness();

    if (type.compare("elastic")==0) return stiff;
    else if (type.compare("secant")==0){  //not implemented, used unloading
        stiff[0] *= (1 - temp_damage.getX() );
        stiff[1] *= (1 -  0.5 * ( temp_damage.getY() + temp_damage.getZ() ) );
        return  stiff;
    } else if (type.compare("unloading")==0){
        stiff[0] *= (1 - temp_damage.getX() );
        stiff[1] *= (1 -  0.5 * ( temp_damage.getY() + temp_damage.getZ() ) );
        return  stiff;
    } else if (type.compare("tangent")==0){
        stiff[0] *= (1 - temp_damage.getX() );
        stiff[1] *= (1 -  0.5 * ( temp_damage.getY() + temp_damage.getZ() ) );
        return  stiff;
    } else {
        cerr << "Error: FatigueShearMaterialStatus does not provide '"<< type << "' stiffness";
        exit(1);
    };
}


//////////////////////////////////////////////////////////
void AllicheMaterialStatus :: update() {
  damage = temp_damage;
  Y_prev  = Y;
  eps_plus_prev = eps_plus;
  // shear_eps_prev = shear_eps_cur;
}


//////////////////////////////////////////////////////////
// ALLICHE MATERIAL
//////////////////////////////////////////////////////////
void AllicheMaterial :: readFromLine(istringstream &iss) {
    DisMechMaterial :: readFromLine(iss); //read elastic parameters

    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream

    string param;
    bool  bC0, bC1, bK, bG, bAl, bBe, bN;
    bC0 = bC1 = bK = bG = bAl = bBe = bN = false;
    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("C0") == 0 ) {
            bC0 = true;
            iss >> C0;
        } else if ( param.compare("C1") == 0 ) {
            bC1 = true;
            iss >> C1;
        } else if ( param.compare("K") == 0 )    {
            bK = true;
            iss >> K;
        } else if ( param.compare("g") == 0 )    {
            bG = true;
            iss >> g;
        } else if ( param.compare("alphaDam") == 0 )    {
            bAl = true;
            iss >> alphaDam;
        } else if ( param.compare("betaDam") == 0 )    {
            bBe = true;
            iss >> betaDam;
        } else if ( param.compare("n") == 0 )    {
            bN = true;
            iss >> n;
        }
    }
    if ( !bC0 ) {
        cerr << name << ": material parameter 'C0' was not specified, taking C0 = 0.0" << endl;
    }
    if ( !bC1 ) {
        cerr << name << ": material parameter 'C1' was not specified" << endl;
        exit(1);
    }
    if ( !bK ) {
        cerr << name << ": material parameter 'K' was not specified" << endl;
        exit(1);
    }
    if ( !bG ) {
        cerr << name << ": material parameter 'g' was not specified" << endl;
        exit(1);
    }
    if ( !bAl ) {
        cout << name << ": material parameter 'alphaDam' was not specified" << endl;
        exit(1);
    }
    if ( !bBe ) {
        cout << name << ": material parameter 'betaDam' was not specified" << endl;
        exit(1);
    }
    if ( !bN ) {
        cout << name << ": material parameter 'n' was not specified" << endl;
        exit(1);
    }
    // calculate Lamé constatnts from elastic parameter
    lambda = E0 * alpha / ( ( 1 + alpha ) * ( 1 - 2 * alpha ) );
    mu = E0 / ( 2 * ( 1 + alpha ) );
};

//////////////////////////////////////////////////////////
MaterialStatus *AllicheMaterial :: giveNewMaterialStatus(Element *e) {
    AllicheMaterialStatus *newStatus = new AllicheMaterialStatus(this, e); //needs to be deleted manually
    return newStatus;
};

//////////////////////////////////////////////////////////
void AllicheMaterial :: init() {

};


//////////////////////////////////////////////////////////
// DESMORAT MATERIAL STATUS
//////////////////////////////////////////////////////////
DesmoratMaterialStatus :: DesmoratMaterialStatus(DesmoratMaterial *m, Element *e) : DisMechMaterialStatus(m, e) {
  name = "Desmorat mat. status";
}

void DesmoratMaterialStatus :: init(){
  temp_sigma = sigma = Point(); ///< stress
  temp_Y = Y = 0; ///< energy release rate
  epsN = 0;
  epsT = Point(); ///< strain
  temp_epsPi = epsPi = Point(); ///< irreversible strain
  temp_damage = damage = 0; ///< damage
  temp_zIso = zIso = 0;
  temp_alphaKin = alphaKin = Point();

}

void DesmoratMaterialStatus :: update(){
  sigma = temp_sigma;
  Y = temp_Y;
  epsPi = temp_epsPi;
  damage = temp_damage;
  zIso = temp_zIso;
  alphaKin = temp_alphaKin;
}

Vector DesmoratMaterialStatus :: giveNormalShearStiffness(string type) const {
  Vector stiff = giveElasticNormalShearStiffness();
  if (type.compare("elastic")==0) return stiff;
  else if (type.compare("secant")==0) return  stiff * ( 1 - temp_damage );
  else if (type.compare("unloading")==0) return  stiff * ( 1 - temp_damage );
  else if (type.compare("tangent")==0) return  stiff * ( 1 - temp_damage ); //not implemented, used unloading
  else{
      cerr << "Error: DesmoratMatStatus does not provide '"<< type << "' stiffness";
      exit(1);
  };
}

Vector DesmoratMaterialStatus :: giveStress(const Vector &strain){
    Vector stress( strain.size());

    DesmoratMaterial *m = static_cast< DesmoratMaterial * >( mat );
    for (unsigned i = 0; i < strain.size(); i++ ){
      if ( i == 0 ){ epsN = strain[ i ]; }
      else if ( i == 1 ){ epsT.setY( strain[ i ] ); }
      else if ( i == 2 ){ epsT.setZ( strain[ i ] ); }
    }

    Point sigmaPiTilda = (epsT - epsPi) * m->giveE2();

    double treshold_func = (sigmaPiTilda - alphaKin * m->giveGamma()).norm() - m->giveK() * zIso - m->giveSigma0();

    if ( treshold_func > 0 ){
      double deltaPi = treshold_func / ( m->giveE2() + ( m->giveK() * m->giveGamma() ) * ( 1 - damage ) );
      Point deltaEps = ( ( sigmaPiTilda - alphaKin * m->giveGamma() ) / ( ( sigmaPiTilda - alphaKin * m->giveGamma() ).norm() ) ) * deltaPi;
      temp_epsPi = epsPi + deltaEps;
      temp_Y = 0.5 * m->giveE0() * pow( epsN, 2 ) + 0.5 * m->giveE2() * ( epsT - temp_epsPi ).sqNorm();
      temp_damage = fmax(0, fmin(damage + ( temp_Y / m->giveS() ) * deltaPi, 1 ) );
      temp_zIso = zIso + deltaPi * ( 1 - temp_damage );
      temp_alphaKin = alphaKin + deltaEps * ( 1 - temp_damage );
    } else {
      temp_damage = damage;
    }

    stress[ 0 ] = m->giveE0() * ( 1 - temp_damage ) * epsN;
    stress[ 1 ] = m->giveE2() * ( 1 - temp_damage ) * epsT.getY();
    if ( strain.size() > 1) stress[ 2 ] = m->giveE2() * ( 1 - temp_damage ) * epsT.getZ();

    return stress;
}

double DesmoratMaterialStatus :: giveValue(string code) const {
  if ( (code.compare("damage") == 0) ) {
      return damage;
  } else if ( (code.compare("stressX") == 0) ) {
      return sigma.getX();
  } else if ( (code.compare("stressY") == 0) ) {
      return sigma.getY();
  } else if ( (code.compare("stressZ") == 0) ) {
      return sigma.getZ();
  } else {
      return DisMechMaterialStatus :: giveValue(code);
  }
}


//////////////////////////////////////////////////////////
// DESMORAT MATERIAL
//////////////////////////////////////////////////////////
void DesmoratMaterial :: readFromLine(istringstream &iss) {
    DisMechMaterial :: readFromLine(iss); //read elastic parameters

    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream

    string param;
    bool  bS0, bK, bG, bS;
    bS0 = bK = bG = bS = false;
    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("Sigma0") == 0 ) {
            bS0 = true;
            iss >> Sigma0;
        } else if ( param.compare("K") == 0 )    {
            bK = true;
            iss >> K;
        } else if ( param.compare("gamma") == 0 )    {
            bG = true;
            iss >> gamma;
        } else if ( param.compare("S") == 0 )    {
            bS = true;
            iss >> S;
        }
    }
    if ( !bS0 ) {
        cerr << name << ": material parameter 'Sigma0' was not specified" << endl;
        exit(1);
    }
    if ( !bK ) {
        cerr << name << ": material parameter 'K' was not specified" << endl;
        exit(1);
    }
    if ( !bG ) {
        cerr << name << ": material parameter 'gamma' was not specified" << endl;
        exit(1);
    }
    if ( !bS ) {
        cout << name << ": material parameter 'S' was not specified" << endl;
        exit(1);
    }
};

//////////////////////////////////////////////////////////
MaterialStatus *DesmoratMaterial :: giveNewMaterialStatus(Element *e) {
    DesmoratMaterialStatus *newStatus = new DesmoratMaterialStatus(this, e); //needs to be deleted manually
    return newStatus;
};

//////////////////////////////////////////////////////////
void DesmoratMaterial :: init() {

};
