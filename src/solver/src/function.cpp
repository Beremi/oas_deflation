#include "function.h"

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// PIECE-WISE LINEAR FUNCTION
void PieceWiseLinearFunction :: readFromLine(istringstream &iss) {
    unsigned num;
    iss >> num;
    x.resize(num);
    y.resize(num);
    for ( unsigned i = 0; i < num; i++ ) {
        iss >> x [ i ];
    }
    for ( unsigned i = 0; i < num; i++ ) {
        iss >> y [ i ];
    }
}

//////////////////////////////////////////////////////////
double PieceWiseLinearFunction :: giveY(double t)  {

	if (x.size() <= 0) return 0;

	//indexovani ma typ podle typu velikosti pole
	for (size_t i = 0; i < x.size(); i++) {
		if (x[i] > t)	return y[i - 1] + (y[i] - y[i - 1]) / (x[i] - x[i - 1]) * (t - x[i - 1]);
	}

	//nenasel jsem vetsi prvek, takze vracim posledni prvek pole
	return y[x.size() - 1];


	/*  Tohle je spatne poradi podminek: index i vzdy presahne delku pole a vznikne neopravneny pristup do pameti
    unsigned i = 0;
    // while ( x [ i ] < t && i < x.size() ) {
    // ok, this way would also be sufficient
    while ( i < x.size() && x [ i ] < t ) {
        i++;
    }
    if ( i == 0 ) {
        return 0.;
    } else if ( i == x.size() ) {
        return y [ x.size() - 1 ];
    } else {
        return y [ i - 1 ] + ( y [ i ] - y [ i - 1 ] ) / ( x [ i ] - x [ i - 1 ] ) * ( t - x [ i - 1 ] );
    }
	*/
}

//////////////////////////////////////////////////////////
double PieceWiseLinearFunction :: giveNextEtreme(const double &t) const {

	if(x.size() <= 0) return INFINITY;

	//indexovani ma typ podle typu velikosti pole
	for (size_t i = 0; i < x.size(); i++) {
		if (x[i] > t) return x[i];
	}

	return INFINITY;

  /*  Tohle je spatne poradi podminek: index i vzdy presahne delku pole a vznikne neopravneny pristup do pameti
  unsigned i = 0;
  // while ( x[i] <= t && i < x.size()) {
  while ( i < x.size() && x[i] <= t) {
      i++;
  }
  if ( i == 0 || i == x.size() ) {
	  return INFINITY;
  } else {
    return x [ i ];
  }
  */

}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// SAW TOOTH FUNCTION WITH CONSTATNT MAX VALUE
void ConstSawToothFunction :: readFromLine(istringstream &iss) {
  iss.clear(); // clear string stream
  iss.seekg(0, iss.beg); //reset position in string stream
  string param;
  double timo, temp;
  int num_cycles;
  bool bup, blow, bper, btim, bcyc, sym;
  bup = blow = bper = btim = bcyc = sym = false;
  while ( !iss.eof() ) {
    iss >> param;
    if ( param.compare("lower") == 0 ) {
      iss >> lower;
      blow = true;
    } else if ( param.compare("value") == 0 ){
      iss >> upper;
      bup = true;
    } else if ( param.compare("period") == 0 ){
      iss >> period;
      bper = true;
    } else if ( param.compare("time") == 0 ){
      iss >> timo;
      btim = true;
    } else if ( param.compare("num_cycles") == 0 ){
      iss >> num_cycles;
      bcyc = true;
    } else if ( param.compare("sym") == 0 ){
      iss >> sym;
    }
  }

  if ( !bper ) {
    if (bcyc && btim){
      period = timo / num_cycles;
      // cout << " function parameter 'period' for function " << typeid(this).name() << " calculated from number of cycles per given time" << endl;
    } else {
      cerr << " function parameter 'period' for function " << typeid(this).name() << " was not specified" << endl;
      exit(0);
    }
  }
  if ( !bup ) {
    cerr << " function parameter 'value' for function " << typeid(this).name() << " was not specified" << endl;
    exit(0);
  }
  if ( !blow ) {
    if (sym){
      lower = -upper;
    } else {
      lower = 0;
    }
  }
  if (lower < 0 && upper < 0){
    multip = -1;
    lower *= multip;
    upper *= multip;
  }
  if (upper < lower){
    temp = lower;
    lower = upper;
    upper = temp;
  }
  time_shift = 0.5 * period * abs(lower) / ( upper - lower );
  if (lower > 0){
    time_shift *= -1;
  }
}

//////////////////////////////////////////////////////////
double ConstSawToothFunction :: giveY(double t)  {
  if (lower > 0 && t < abs(time_shift)){
    return multip * t * lower / ( abs( time_shift ) );
  } else {
    return ( (upper - lower) - ( ( ( upper - lower ) / ( 0.5 * period ) ) *
            abs( fmod( ( t + time_shift ), period ) - 0.5 * period ) ) +
            lower ) * multip;
  }
}

//////////////////////////////////////////////////////////
double ConstSawToothFunction :: giveNextEtreme(const double &t) const {
  if ( t < abs(time_shift) )
    return time_shift;
  else
    return ( ( int( ( ( t - abs(time_shift) ) / ( 0.5 * period ) ) + 1.0 ) ) *
            ( 0.5 * period ) ) + abs(time_shift);
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// SAW TOOTH FUNCTION WITH CONSTATNT MAX VALUE
void LinSawToothFunction :: readFromLine(istringstream &iss) {
  // read the same constatnts from constant function
  ConstSawToothFunction :: readFromLine(iss);
  iss.clear(); // clear string stream
  iss.seekg(0, iss.beg); //reset position in string stream
  string param;
  double timo, val;
  bool bmult, btim, bval;
  bmult = btim = bval = false;

  while ( !iss.eof() ) {
    iss >> param;
    if ( param.compare("multip") == 0 ) {
      iss >> time_multiplier;
      bmult = true;
    } else if ( param.compare("time") == 0 ) {
      iss >> timo;
      btim = true;
    } else if ( param.compare("value") == 0 ) {
      iss >> val;
      bval = true;
    }
  }
  if ( !bmult ) {
    if (btim && bval){
      time_multiplier = 1/timo;
      // std::cout << "timo = " << timo << ", value = " << val << ", time_multiplier = " << time_multiplier << '\n';
    } else {
      cerr << " function parameter 'time_multiplier' for function " << typeid(this).name() << " was not specified" << endl;
      exit(0);
    }
  }
}

//////////////////////////////////////////////////////////
double LinSawToothFunction :: giveY(double t)  {
  // only multiply result of the constant saw tooth function by linearly increasing time function SawToot(t) * (t_m * t)
  double value = ConstSawToothFunction :: giveY(t);
  return value * time_multiplier * t;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// SAW TOOTH FUNCTION WITH VARYING MAX VALUE
void VaryingSawToothFunction :: readFromLine(istringstream &iss) {
  // read the same constatnts from constant function
  PieceWiseLinearFunction :: readFromLine(iss);
  iss.clear(); // clear string stream
  iss.seekg(0, iss.beg); //reset position in string stream
  ConstSawToothFunction :: readFromLine(iss);
  iss.clear(); // clear string stream
  iss.seekg(0, iss.beg); //reset position in string stream
  string param;
  while ( !iss.eof() ) {
    iss >> param;
    if ( param.compare("shift") == 0 ) {
      iss >> shift;
    }
  }
}

//////////////////////////////////////////////////////////
double VaryingSawToothFunction :: giveY(double t)  {
  // only multiply result of the constant saw tooth function by linearly increasing time function SawToot(t) * (t_m * t)
  double value = ConstSawToothFunction :: giveY(t);
  double value2 = PieceWiseLinearFunction :: giveY(t);
  return (value * value2) + shift;
}


//////////////////////////////////////////////////////////
// JM: Rotation function using angle from a ConstSawToothFunction
void ConstSawToothRotationFunction :: readFromLine(istringstream &iss){
  // read ConstSawToothFunction from file
  ConstSawToothFunction :: readFromLine (iss);
  iss.clear(); // clear string stream
  iss.seekg(0, iss.beg); //reset position in string stream
  string param;
  //
  double angX;
  double angY;
  double angZ;
  double nodeX;
  double nodeY;
  double nodeZ;
  //
  currentTime = 0;
  previousTime = 0;
  //
  while ( !iss.eof() ){
    iss >> param;
    //cout << param << endl;
    if ( param.compare("rotAngles") == 0 ){
      iss >> angX >> angY >> angZ;
      rotationAngles = Point (angX, angY, angZ);
    }
    else if ( param.compare("initNodeCrds") == 0){
      iss >> nodeX >> nodeY >> nodeZ;
      initNodePosition = Point (nodeX, nodeY, nodeZ);
    }
    else if ( param.compare("displType") == 0){
      iss >> displacementType;
    }
  }

}

void ConstSawToothRotationFunction :: setCurrentTime(double t){
  if (t >= currentTime){
    previousTime = currentTime;
    currentTime = t;
  }
  if (t < currentTime){
    currentTime = t;
  }
}
// JM: Return value depending on displacement type
//////////////////////////////////////////////////////////
double ConstSawToothRotationFunction :: giveY(double t)  {
  //
  ConstSawToothRotationFunction :: setCurrentTime(t);
  //
  double currentAngleX  = rotationAngles.x * ConstSawToothFunction :: giveY(currentTime);
  double previousAngleX = 0;

  // delta Coordinate X
  if (displacementType == 0 ){
      double val = 0;
      return val;
  }
  // delta Coordinate Y
  else if ( displacementType == 1){
      double val = - (cos(currentAngleX) *initNodePosition.y - sin(currentAngleX) *initNodePosition.z)
                   + (cos(previousAngleX)*initNodePosition.y - sin(previousAngleX)*initNodePosition.z);
      return val;
  }
  // delta Coordinate Z
  else if ( displacementType == 2){
      double val = - (sin(currentAngleX) *initNodePosition.y + cos(currentAngleX) *initNodePosition.z)
                   + (sin(previousAngleX)*initNodePosition.y + cos(previousAngleX)*initNodePosition.z);
      return val;
  }
  // Angle around X
  /*
  else if ( displacementType == 3){
    return 0;
    //  return currentAngleX;
  }
  // Angle around Y
  else if ( displacementType == 4){
      return 0;
    //  return currentAngleY;
  }
  // Angle around Z
  else if ( displacementType == 5){
      return 0;
    //  return currentAngleZ;
  }
  //return 0;
  */
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// SAW TOOTH FUNCTION WITH CONSTATNT MAX VALUE
void SinusFunction :: readFromLine(istringstream &iss) {
  iss.clear(); // clear string stream
  iss.seekg(0, iss.beg); //reset position in string stream
  string param;
  bool bper, bamp, bshift;
  bper = bamp = bshift = false;
  while ( !iss.eof() ) {
    iss >> param;
    if ( param.compare("period") == 0 ) {
      iss >> period;
      bper = true;
    } else if ( param.compare("amplitude") == 0 ){
      iss >> amplitude;
      bamp = true;
    } else if ( param.compare("shift") == 0 ){
      iss >> shift;
      bshift = true;
    }
  }
  if ( !bper ) {
    cerr << " function parameter 'period' for function " << typeid(this).name() << " was not specified" << endl;
    exit(0);
  }
  if ( !bamp ) {
    cerr << " function parameter 'amplitude' for function " << typeid(this).name() << " was not specified" << endl;
    exit(0);
  }
  if ( !bshift ) {
    // cout << " function parameter 'shift' (optional) for function " << typeid(this).name() << " was not specified" << endl;
  }
}

//////////////////////////////////////////////////////////
double SinusFunction :: giveY(double t)  {
  return amplitude * sin(2 * M_PI * t / period) + shift;
}

//////////////////////////////////////////////////////////
double SinusFunction :: giveNextEtreme(const double &t) const {
  return ( ( int ( ( t ) / ( 0.5 * period ) ) + 1 ) * ( 0.5 * period ) );
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CONTAINER FOR FUNCTIONS
FunctionContainer :: ~FunctionContainer() {
    for ( vector< Function * > :: iterator f = functions.begin(); f != functions.end(); ++f ) {
        delete * f;
    }
}

//////////////////////////////////////////////////////////
void FunctionContainer :: readFromFile(const string filename) {
    size_t origsize = functions.size();
    string line, ftype;
    ifstream inputfile(filename.c_str() );
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std::ws, line) ) {
            if ( line.empty() ){
                continue;
            }
            if ( line.at(0) == '#' ) {
                continue;
            }
            istringstream iss(line);
            iss >> ftype;
            if ( !ftype.rfind("#", 0) == 0 ) {
                if ( ftype.compare("PWLFunction") == 0 ) {
                    PieceWiseLinearFunction *newf = new PieceWiseLinearFunction();
                    newf->readFromLine(iss);
                    functions.push_back(newf);
                } else if ( ftype.compare("ConstSawToothFn") == 0 ) {
                    ConstSawToothFunction *newf = new ConstSawToothFunction();
                    newf->readFromLine(iss);
                    functions.push_back(newf);
                } else if ( ftype.compare("LinSawToothFn") == 0 ) {
                    LinSawToothFunction *newf = new LinSawToothFunction();
                    newf->readFromLine(iss);
                    functions.push_back(newf);
                } else if ( ftype.compare("SinusFn") == 0 ) {
                    SinusFunction *newf = new SinusFunction();
                    newf->readFromLine(iss);
                    functions.push_back(newf);
                } else if ( ftype.compare("VaryingSawToothFn") == 0 ) {
                    VaryingSawToothFunction *newf = new VaryingSawToothFunction();
                    newf->readFromLine(iss);
                    // it is necessary to specify which of two parent classes will the tree go throug
                    functions.push_back((ConstSawToothFunction *) newf);
                }
                //JM
                else if ( ftype.compare("ConstSawToothRotationFunction") == 0){
                  ConstSawToothRotationFunction *newf = new ConstSawToothRotationFunction();
                  newf->readFromLine(iss);
                  functions.push_back(newf);
                }
                else  {
                    cerr << "Error: function '" <<  ftype <<  "' is not implemented yet." << endl;
                    exit(0);
                }
            }
        }
        inputfile.close();
        cout << "Input file '" <<  filename << "' succesfully loaded; " << functions.size() - origsize << " functions found" << endl;
    } else {
        cerr << "Error: unable to open input file '" <<  filename <<  "'" << endl;
        exit(0);
    }
}

//////////////////////////////////////////////////////////
double FunctionContainer :: giveY(unsigned f, double t)  {
    return functions [ f ]->giveY(t);
}

//////////////////////////////////////////////////////////
/*
  returns time of next extreme - important point on time scale
  peaks for periodic functions (SawTooth, Sinus ... )
  changes of slope for PieceWiseLinearFunction
*/
double FunctionContainer :: giveTimeOfNextExtreme(const double &t) const {
  double nextExtreme = INFINITY;
  double thisFnExtreme;
  for ( auto const &fn : functions ){
    thisFnExtreme = fn->giveNextEtreme( t );
    if ( fn->isActive() && thisFnExtreme < nextExtreme && thisFnExtreme > t ){
      nextExtreme = thisFnExtreme;
    }
  }
  return nextExtreme;
}
