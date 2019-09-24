#include "boundarycondition.h"
#include "nodecontainer.h"
#include "linearalgebra.h" 

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
    unsigned i = 0;
    while ( x [ i ] < t && i < x.size() ) {
        i++;
    }
    if ( i == 0 ) {
        return 0.;
    } else if ( i == x.size() )   {
        return y [ x.size() - 1 ];
    } else                                                       {
        return y [ i - 1 ] + ( y [ i ] - y [ i - 1 ] ) / ( x [ i ] - x [ i - 1 ] ) * ( t - x [ i - 1 ] );
    }
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
      cout << " function parameter 'period' for function " << typeid(this).name() << " calculated from number of cycles per given time" << endl;
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
}

//////////////////////////////////////////////////////////
double ConstSawToothFunction :: giveY(double t)  {
  double time_shift, up_low;
  up_low = upper - lower;
  time_shift = 0.5 * period * abs(lower) / up_low;
  if (lower > 0){
    time_shift *= -1;
  }
  if (lower > 0 && t < abs(time_shift)){
    return multip * t * lower / (abs(time_shift));
  } else {
    return (up_low - ((up_low/(0.5*period))*abs(fmod((t + time_shift), period) - 0.5 * period)) + lower) * multip;
  }
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
      std::cout << "timo = " << timo << ", value = " << val << ", time_multiplier = " << time_multiplier << '\n';
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
  if (t>0){
    rotationAngles.x * ConstSawToothFunction :: giveY(previousTime);
  }

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
    cout << " function parameter 'shift' (optional) for function " << typeid(this).name() << " was not specified" << endl;
  }
}

//////////////////////////////////////////////////////////
double SinusFunction :: giveY(double t)  {
  return amplitude * sin(2 * M_PI * t / period) + shift;
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
    unsigned origsize = functions.size();
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
//////////////////////////////////////////////////////////
// DIRRICHLET AND NEUMANN BOUNDARY CONDITION
void BoundaryCondition :: init() {
    blockedDoFNum = 0;
    loadedDoFNum = 0;
    for ( vector< int > :: const_iterator i = dirrichBC.begin(); i != dirrichBC.end(); ++i ) {
        if ( * i >= 0 ) {
            blockedDoFNum++;
        }
    }
    for ( vector< int > :: const_iterator i = neumannBC.begin(); i != neumannBC.end(); ++i ) {
        if ( * i >= 0 ) {
            loadedDoFNum++;
        }
    }
}

//////////////////////////////////////////////////////////
vector< unsigned >BoundaryCondition :: giveBlockedDoFs() const {
    vector< unsigned >blocked;
    blocked.resize(blockedDoFNum);
    unsigned k = 0;
    unsigned s = 0;
    for ( vector< int > :: const_iterator i = dirrichBC.begin(); i != dirrichBC.end(); ++i, k++ ) {
        if ( * i >= 0 ) {
            blocked [ s ] = node->giveStartingDoF() + k;
            s++;
        }
    }
    ;
    return blocked;
}

//////////////////////////////////////////////////////////
vector< unsigned >BoundaryCondition :: giveLoadedDoFs() const {
    vector< unsigned >loaded;
    loaded.resize(loadedDoFNum);
    unsigned k = 0;
    unsigned s = 0;
    for ( vector< int > :: const_iterator i = neumannBC.begin(); i != neumannBC.end(); ++i, k++ ) {
        if ( * i >= 0 ) {
            loaded [ s ] = node->giveStartingDoF() + k;
            s++;
        }
    }
    ;
    return loaded;
}

//////////////////////////////////////////////////////////
vector< unsigned >BoundaryCondition :: giveBlockedFunctions() const {
    vector< unsigned >blocked;
    blocked.resize(blockedDoFNum);
    unsigned k = 0;
    unsigned s = 0;
    for ( vector< int > :: const_iterator i = dirrichBC.begin(); i != dirrichBC.end(); ++i, k++ ) {
        if ( * i >= 0 ) {
            blocked [ s ] = dirrichBC [ k ];
            s++;
        }
    }
    ;
    return blocked;
}

//////////////////////////////////////////////////////////
vector< unsigned >BoundaryCondition :: giveLoadedFunctions() const {
    vector< unsigned >loaded;
    loaded.resize(loadedDoFNum);
    unsigned k = 0;
    unsigned s = 0;
    for ( vector< int > :: const_iterator i = neumannBC.begin(); i != neumannBC.end(); ++i, k++ ) {
        if ( * i >= 0 ) {
            loaded [ s ] = neumannBC [ k ];
            s++;
        }
    }
    ;
    return loaded;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CONTAINER FOR BOUNDARY CONDITIONS
BCContainer :: ~BCContainer() {
    for ( vector< BoundaryCondition * > :: iterator bc = BC.begin(); bc != BC.end(); ++bc ) {
        delete * bc;
    }
}

//////////////////////////////////////////////////////////
void BCContainer :: readFromFile(const string filename, NodeContainer *nodes) {
    unsigned origsize = BC.size();
    string line, aux;
    unsigned intnum, nDoFs;
    vector< int >dirrichBC, neumannBC;
    Node *node;
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
            iss >> intnum;
            node = nodes->giveNode(intnum);
            nDoFs = node->giveNumberOfDoFs();
            dirrichBC.resize(nDoFs);
            neumannBC.resize(nDoFs);
            for ( unsigned i = 0; i < nDoFs; i++ ) {
                iss >> dirrichBC [ i ];
            }
            for ( unsigned i = 0; i < nDoFs; i++ ) {
                iss >> neumannBC [ i ];
            }
            for ( unsigned i = 0; i < nDoFs; i++ ) {
                if ( neumannBC [ i ] >= 0 && dirrichBC [ i ] >= 0 ) {
                    cerr << "Error: Dirrichlet and Neumann boundary conditions assigned simulatneuosly" << endl;
                    cerr << line << endl;
                    exit(0);
                }
            }
            ;
            BoundaryCondition *newBC = new BoundaryCondition(node, dirrichBC, neumannBC);
            BC.push_back(newBC);
        }
        inputfile.close();
        for ( unsigned i = 0; i < BC.size() - origsize; i++ ) {
            BC [ i + origsize ]->giveNode()->setBC(BC [ i + origsize ]);
        }
        cout << "Input file '" <<  filename << "' succesfully loaded; " << BC.size() - origsize << " boundary conditions found" << endl;
    } else {
        cerr << "Error: unable to open input file '" <<  filename <<  "'" << endl;
        exit(0);
    }
}

//////////////////////////////////////////////////////////
void BCContainer :: init() {
    dirrichDoFs.resize(0);
    neumannDoFs.resize(0);
    for ( vector< BoundaryCondition * > :: iterator bc = BC.begin(); bc != BC.end(); ++bc ) {
        ( * bc )->init();
    }
}

//////////////////////////////////////////////////////////
void BCContainer :: calculateDoFfields() {
    vector< unsigned >help;
    for ( vector< BoundaryCondition * > :: iterator bc = BC.begin(); bc != BC.end(); ++bc ) {
        help = ( * bc )->giveBlockedDoFs();
        dirrichDoFs.insert(dirrichDoFs.end(), help.begin(), help.end() );
        help = ( * bc )->giveLoadedDoFs();
        neumannDoFs.insert(neumannDoFs.end(), help.begin(), help.end() );
        help = ( * bc )->giveBlockedFunctions();
        dirrichF.insert(dirrichF.end(), help.begin(), help.end() );
        help = ( * bc )->giveLoadedFunctions();
        neumannF.insert(neumannF.end(), help.begin(), help.end() );
    }
}

//////////////////////////////////////////////////////////
vector< double >BCContainer :: giveBlockedDoFValues(double t) const {
    vector< double >blocked(dirrichDoFs.size() );
    for ( unsigned h = 0; h < dirrichDoFs.size(); h++ ) {
        blocked [ h ] = functions->giveY(dirrichF [ h ], t);
    }
    return blocked;
}

//////////////////////////////////////////////////////////
vector< double >BCContainer :: giveLoadedDoFValues(double t) const {
    vector< double >loaded(neumannDoFs.size() );
    for ( unsigned h = 0; h < neumannDoFs.size(); h++ ) {
        loaded [ h ] = functions->giveY(neumannF [ h ], t);
    }
    return loaded;
}
