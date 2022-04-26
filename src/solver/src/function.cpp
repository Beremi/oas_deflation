#include "function.h"

using namespace std;

#ifdef __EXPRTK_MODULE
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// GENERAL SPATIAL FUNCTION
void GeneralSpatialFunction :: readFromLine(istringstream &iss) {
    iss >> expression_string;

    symbol_table_t symbols;
    parser_t parser;
    symbols.add_variable("x", x);
    symbols.add_variable("y", y);
    symbols.add_variable("z", z);
    symbols.add_constants();
    expression.register_symbol_table(symbols);
    parser.compile(expression_string, expression);
}

//////////////////////////////////////////////////////////
double GeneralSpatialFunction :: giveY(const Point *xyz) {
    x = xyz->x();
    y = xyz->y();
    z = xyz->z();

    return expression.value();
}

//////////////////////////////////////////////////////////
double GeneralSpatialFunction :: giveNextEtreme(const double &t) const {
    ( void ) t;
    return INFINITY;
}
#endif

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
double PieceWiseLinearFunction :: giveY(double t) const {
    if ( x.size() <= 0 ) {
        return 0;
    }

    //indexovani ma typ podle typu velikosti pole
    for ( size_t i = 0; i < x.size(); i++ ) {
        if ( x [ i ] > t ) {
            return y [ i - 1 ] + ( y [ i ] - y [ i - 1 ] ) / ( x [ i ] - x [ i - 1 ] ) * ( t - x [ i - 1 ] );
        }
    }

    //nenasel jsem vetsi prvek, takze vracim posledni prvek pole
    return y [ x.size() - 1 ];
}

//////////////////////////////////////////////////////////
double PieceWiseLinearFunction :: giveNextEtreme(const double &t) const {
    ( void ) t;
    return INFINITY;
}

//////////////////////////////////////////////////////////
double PieceWiseLinearFunctionWithExtremes :: giveNextEtreme(const double &t) const {
    if ( x.size() <= 0 ) {
        return INFINITY;
    }

    //indexovani ma typ podle typu velikosti pole
    for ( size_t i = 0; i < x.size(); i++ ) {
        if ( x [ i ] > t ) {
            return x [ i ];
        }
    }

    return INFINITY;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// SAW TOOTH FUNCTION WITH CONSTATNT MAX VALUE
void ConstSawToothFunction :: readFromLine(istringstream &iss) {
    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream
    string param;
    double timo = 0;
    int num_cycles = 1;
    double temp;
    bool bup, blow, bper, btim, bcyc, sym;
    bup = blow = bper = btim = bcyc = sym = false;
    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("lower") == 0 ) {
            iss >> lower;
            blow = true;
        } else if ( param.compare("value") == 0 ) {
            iss >> upper;
            bup = true;
        } else if ( param.compare("period") == 0 ) {
            iss >> period;
            bper = true;
        } else if ( param.compare("time") == 0 ) {
            iss >> timo;
            btim = true;
        } else if ( param.compare("num_cycles") == 0 ) {
            iss >> num_cycles;
            bcyc = true;
        } else if ( param.compare("sym") == 0 ) {
            iss >> sym;
        }
    }

    if ( !bper ) {
        if ( bcyc && btim ) {
            period = timo / num_cycles;
            // cout << " function parameter 'period' for function " << typeid(this).name() << " calculated from number of cycles per given time" << endl;
        } else {
            cerr << " function parameter 'period' for function " << typeid( this ).name() << " was not specified" << endl;
            exit(1);
        }
    }
    if ( !bup ) {
        cerr << " function parameter 'value' for function " << typeid( this ).name() << " was not specified" << endl;
        exit(1);
    }
    if ( !blow ) {
        if ( sym ) {
            lower = -upper;
        } else {
            lower = 0;
        }
    }
    if ( lower < 0 && upper < 0 ) {
        multip = -1;
        lower *= multip;
        upper *= multip;
    }
    if ( upper < lower ) {
        temp = lower;
        lower = upper;
        upper = temp;
    }
    time_shift = 0.5 * period * abs(lower) / ( upper - lower );
    if ( lower > 0 ) {
        time_shift *= -1;
    }
}

//////////////////////////////////////////////////////////
double ConstSawToothFunction :: giveY(double t) const {
    if ( lower > 0 && t < abs(time_shift) ) {
        return multip * t * lower / ( abs(time_shift) );
    } else {
        return ( ( upper - lower ) - ( ( ( upper - lower ) / ( 0.5 * period ) ) *
                                       abs(fmod( ( t + time_shift ), period) - 0.5 * period) ) +
                 lower ) * multip;
    }
}

//////////////////////////////////////////////////////////
double ConstSawToothFunction :: giveNextEtreme(const double &t) const {
    if ( t < abs(time_shift) ) {
        return time_shift;
    } else {
        return ( ( int( ( ( t - abs(time_shift) ) / ( 0.5 * period ) ) + 1.0 ) ) *
                 ( 0.5 * period ) ) + abs(time_shift);
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
            iss >> timo;  //TODO potentially uninitialized local variable 'timo' used
            btim = true;
        } else if ( param.compare("value") == 0 ) {
            iss >> val;
            bval = true;
        }
    }
    if ( !bmult ) {
        if ( btim && bval ) {
            time_multiplier = 1 / timo; //todo:  warning C4701: potentially uninitialized local variable 'timo' used
            // std::cout << "timo = " << timo << ", value = " << val << ", time_multiplier = " << time_multiplier << '\n';
        } else {
            cerr << " function parameter 'time_multiplier' for function " << typeid( this ).name() << " was not specified" << endl;
            exit(EXIT_FAILURE);
        }
    }
}

//////////////////////////////////////////////////////////
double LinSawToothFunction :: giveY(double t) const {
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
double VaryingSawToothFunction :: giveY(double t) const {
    // only multiply result of the constant saw tooth function by linearly increasing time function SawToot(t) * (t_m * t)
    double value = ConstSawToothFunction :: giveY(t);
    double value2 = PieceWiseLinearFunction :: giveY(t);
    return ( value * value2 ) + shift;
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
        } else if ( param.compare("amplitude") == 0 ) {
            iss >> amplitude;
            bamp = true;
        } else if ( param.compare("shift") == 0 ) {
            iss >> shift;
            bshift = true;
        }
    }

    if ( !bper ) {
        cerr << " function parameter 'period' for function " << typeid( this ).name() << " was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    if ( !bamp ) {
        cerr << " function parameter 'amplitude' for function " << typeid( this ).name() << " was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    if ( !bshift ) {
        // cout << " function parameter 'shift' (optional) for function " << typeid(this).name() << " was not specified" << endl;
    }
}

//////////////////////////////////////////////////////////
double SinusFunction :: giveY(double t) const {
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
        if ( * f != nullptr ) {
            delete * f;
        }
    }
}

//////////////////////////////////////////////////////////
void FunctionContainer :: readFromFile(const string filename) {
    size_t origsize = functions.size();
    string line, ftype;
    ifstream inputfile( filename.c_str() );
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std :: ws, line) ) {
            if ( line.empty() || (line.at(0) == '#') ) {
                continue;
            }
            istringstream iss(line);
            iss >> ftype;
            if ( !( ftype.rfind("#", 0) == 0 ) ) {
                if ( ftype.compare("GeneralSpatialFunction") == 0 ) {
#ifdef __EXPRTK_MODULE
                    GeneralSpatialFunction *newf = new GeneralSpatialFunction();
                    newf->readFromLine(iss);
                    functions.push_back(newf);
#else
                    cout << "This binary was not build with EXPRTK module." << endl;
                    exit(EXIT_FAILURE);
#endif
                } else if ( ftype.compare("PWLFunction") == 0 ) {
                    PieceWiseLinearFunction *newf = new PieceWiseLinearFunction();
                    newf->readFromLine(iss);
                    functions.push_back(newf);
                } else if ( ftype.compare("PWLFunctionWithExtremes") == 0 ) {
                    PieceWiseLinearFunctionWithExtremes *newf = new PieceWiseLinearFunctionWithExtremes();
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
                    functions.push_back( ( ConstSawToothFunction * ) newf);
                } else {
                    cerr << "Error: function '" <<  ftype <<  "' is not implemented yet." << endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
        inputfile.close();
        cout << "Input file '" <<  filename << "' succesfully loaded; " << functions.size() - origsize << " functions found" << endl;
    } else {
        cerr << "Error: unable to open input file '" <<  filename <<  "'" << endl;
        exit(EXIT_FAILURE);
    }
}

//////////////////////////////////////////////////////////
Function *FunctionContainer :: giveFunction(unsigned k) {
    if ( k >= functions.size() ) {
        cerr << "Function Container Error: requested function number " << k << " but the container contains only " << functions.size() << " functions" << endl;
        exit(1);
    } else {
        return functions [ k ];
    }
};

//////////////////////////////////////////////////////////
/*
 * returns time of next extreme - important point on time scale
 * peaks for periodic functions (SawTooth, Sinus ... )
 * changes of slope for PieceWiseLinearFunction
 */
double FunctionContainer :: giveTimeOfNextExtreme(const double &t) const {
    double nextExtreme = INFINITY;
    double thisFnExtreme;
    for ( auto const &fn : functions ) {
        thisFnExtreme = fn->giveNextEtreme(t);
        // if ( fn->isActive() && thisFnExtreme < nextExtreme && thisFnExtreme > t ){
        if ( thisFnExtreme < nextExtreme && thisFnExtreme > t ) {
            nextExtreme = thisFnExtreme;
        }
    }
    return nextExtreme;
}

//////////////////////////////////////////////////////////
void FunctionContainer :: removeFunction(unsigned i) {
    if ( i > functions.size() - 1 ) {
        cerr << "FunctionContainer Error: requester function number " << i << " out of " << functions.size() << endl;
        exit(1);
    }
    delete functions [ i ];
    functions [ i ] = nullptr;
    functions.erase(functions.begin() + i);
}
