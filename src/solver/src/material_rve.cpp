#include "material_rve.h"
#include "model.h"
#include "element_discrete.h"
#include "element_ldpm.h"
#include "element_continuous.h"
#include "periodic_bc.h"

using namespace std;

# define M_PI           3.14159265358979323846  /* pi */

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// general RVE status
RVEMaterialStatus :: RVEMaterialStatus(RVEMaterial *m, Element *e, unsigned ipnum, fs :: path masterfile, unsigned ndim) : MaterialStatus(m, e, ipnum) {
    name = "general RVE mat. status";
    inputfile = masterfile;
    RVE = new Model(false);

    axDirs = Matrix :: Zero(ndim, ndim);
    for ( unsigned i = 0; i < ndim; i++ ) {
        axDirs(i, i) = 1.;
    }
}

//////////////////////////////////////////////////////////
RVEMaterialStatus :: ~RVEMaterialStatus() {
    if ( RVE != nullptr ) {
        delete RVE;
    }
}

//////////////////////////////////////////////////////////
void RVEMaterialStatus :: init() {
    RVE->readFromFile( inputfile.string() );
    //here the structre of model initialization should be coppied
    //we needed to insert volumetric average generation after applying preprocessing block, otherwise constraints from preprocessing block would not be active
    //therefore, preprocessing blocks are now called in reader instead of model initialization

    //generateVolumetricAverageBC(); //not good idea because the stiffness matrix becomes full
    generateRandomFixedBC();

    stringstream appendname;
    appendname << "_" << std :: setfill('0') << std :: setw(4) << element->giveID() << "_" << std :: setw(2) << idx;
    RVE->giveExporters()->appendToAllNames( appendname.str() );

    RVE->init();
}

//////////////////////////////////////////////////////////
void RVEMaterialStatus :: update() {
    Solver *solver = RVE->giveSolver();
    Solver *masterSolver = masterModel->giveSolver();
    solver->runAfterEachStep();   //update material statuses
    RVE->giveExporters()->exportData( masterSolver->giveStepNumber(), -1, masterSolver->giveTime(), masterSolver->isTerminated() );
}

//////////////////////////////////////////////////////////
vector< bool >RVEMaterialStatus :: calculateElemDiscreteness() const {
    ElementContainer *elems = RVE->giveElements();
    Material *emat;
    Element *e;
    vector< bool >is_discrete( elems->giveSize() );
    for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
        e = elems->giveElement(i);
        emat = e->giveMaterial();
        if ( dynamic_cast< VectMechMaterial * >( emat ) ) {
            is_discrete [ i ] = true;
        } else if ( dynamic_cast< TensMechMaterial * >( emat ) ) {
            is_discrete [ i ] = false;
        }
    }
    return is_discrete;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// GENERAL RVE MATERIAL
MaterialStatus *RVEMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    RVEMaterialStatus *newstat = new RVEMaterialStatus(this, e, ipnum, inputfile, dim);
    return newstat;
}

//////////////////////////////////////////////////////////
void RVEMaterial :: readFromLine(istringstream &iss) {
    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream

    string param;
    bool bdim, bfile;
    bdim = bfile = false;

    while (  iss >> param ) {
        if ( param.compare("dim") == 0 ) {
            bdim = true;
            iss >> dim;
        } else if ( param.compare("RVEfolder") == 0 ) {
            bfile = true;
            iss >> inputfile;
        } else if ( param.compare("RVEfolder") == 0 ) {
            bfile = true;
            iss >> inputfile;
        } else if ( param.compare("enforce_linearity") == 0 ) {
            nonlinear = false;
        } else if ( param.compare("do_not_precompute") == 0 ) {
            start_from_precomputed = false;
        }
    }
    if ( !bdim ) {
        cerr << name << ": material parameter 'dim' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    if ( !bfile ) {
        cerr << name << ": material parameter 'RVEfolder' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    inputfile = GlobPaths :: BASEDIR  / inputfile;
}

//////////////////////////////////////////////////////////
void RVEMaterial :: init(MaterialContainer *matcont) {
    Material :: init(matcont);
    if ( !nonlinear ) {
        start_from_precomputed = true;
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE TRANSPORT RVE MATERIAL
//////////////////////////////////////////////////////////
DiscreteTransportRVEMaterialStatus :: DiscreteTransportRVEMaterialStatus(RVEMaterial *m, Element *e, unsigned ipnum, fs :: path masterfile, unsigned ndim) : RVEMaterialStatus(m, e, ipnum, masterfile, ndim) {
    name = "transport RVE mat. status";
    is_precomputed = true;
    is_master_status = false;
    temp_nonlin = 1;
    macro_pressure = 0;
}

/////////////////////////////////./////////////////////////
void DiscreteTransportRVEMaterialStatus :: applyEigenStrains() {
    //set eigenstrains and set nonlinearity
    DiscreteTransportRVEMaterial *macromat = static_cast< DiscreteTransportRVEMaterial * >( mat );
    unsigned ndim = macromat->giveDimension();
    ElementContainer *elems = RVE->giveElements();
    Point normal;
    Vector eigstr = Vector :: Zero(1);
    DiscreteTrsprtElem *e;
    VectTrsprtMaterial *m;
    VectTrsprtMaterialStatus *s;
    for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
        eigstr.setZero();
        e = static_cast< DiscreteTrsprtElem * >( elems->giveElement(i) );
        s = static_cast< VectTrsprtMaterialStatus * >( e->giveMatStatus(0) );
        m = static_cast< VectTrsprtMaterial * >( s->giveMaterial() );
        m->setPermeability(orig_mater_params [ 2 * i ]); //set back original permeability
        m->setPermeability( s->calculatePressureDependentPermeability(macro_pressure) );         //calculating pressure depedent conductivity
        m->setParamA(-1.); //switch of linearity
        normal = e->giveNormal();
        for ( unsigned v = 0; v < ndim; v++ ) {
            eigstr [ 0 ] -= local_strain [ v ] * normal(v);
        }
        for ( unsigned k = 0; k < e->giveNumIP(); k++ ) {
            e->giveMatStatus(k)->setEigenStrain(eigstr);
        }
    }
}



/////////////////////////////////./////////////////////////
void DiscreteTransportRVEMaterialStatus :: collectStresses() {
    DiscreteTransportRVEMaterial *macromat = static_cast< DiscreteTransportRVEMaterial * >( mat );
    unsigned ndim = macromat->giveDimension();
    ElementContainer *elems = RVE->giveElements();
    double volume = 0.;
    Vector factor;
    DiscreteTrsprtElem *e;
    VectTrsprtMaterial *m;
    Point normal;
    local_stress.setZero();
    for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
        e = static_cast< DiscreteTrsprtElem * >( elems->giveElement(i) );
        m = static_cast< VectTrsprtMaterial * >( e->giveMatStatus(0)->giveMaterial() );
        m->setParamA(orig_mater_params [ 2 * i + 1 ]); //set back a parameter
        normal = e->giveNormal();
        factor  = e->giveArea() * e->giveLength() * e->giveMatStatus(0)->giveTempStress();

        for ( unsigned v = 0; v < ndim; v++ ) {
            local_stress [ v ] += factor [ 0 ] * normal(v);
        }
        volume += e->giveVolume();
    }
    for ( unsigned v = 0; v < ndim; v++ ) {
        local_stress [ v ] /= volume;
    }
}

//////////////////////////////////////////////////////////
Vector DiscreteTransportRVEMaterialStatus :: giveStressPrecomputed(const Vector &strain, double timeStep) {
    ( void ) timeStep;

    temp_strain = strain;

    transformStrain();    //stiffness is already transformed

    DiscreteTransportRVEMaterial *macromaterial = static_cast< DiscreteTransportRVEMaterial * >( mat );
    temp_nonlin = macromaterial->giveMasterStatus()->calculatePressureDependentPermeability(macro_pressure) / macromaterial->giveMasterMaterial()->givePermeability();
    local_stress = ( macromaterial->givePrecomputedConductivity() * local_strain ) * temp_nonlin;

    transformStress();

    return temp_stress;
}

//////////////////////////////////////////////////////////
bool DiscreteTransportRVEMaterialStatus ::  giveValues(string code, Vector &result) const {
    if ( code.compare("permeability_tensor") == 0 ) {
        //BE CAREFULL, THIS PART OF CODE MIGHT MODIFY THE RESULTS
        DiscreteTransportRVEMaterial *macromat = static_cast< DiscreteTransportRVEMaterial * >( mat );
        unsigned ndim = macromat->giveDimension();
        DiscreteTransportRVEMaterialStatus *newThis = const_cast< DiscreteTransportRVEMaterialStatus * >( this );
        Matrix Keff1 = newThis->giveStiffnessTensorLocalExact("secant");
        Keff1 = ( transf.transpose() * Keff1 ) * transf;
        result = Vector :: Zero(3 * ( ndim - 1 ) );
        for ( unsigned i = 0; i < 3; i++ ) {
            result [ i ] = Keff1(i, i);
        }
        if ( ndim == 2 ) {
            result [ 2 ] = Keff1(1, 0);
        } else if ( ndim == 3 ) {
            result [ 3 ] = 0.5 * ( Keff1(1, 2) + Keff1(2, 1) );
            result [ 4 ] = 0.5 * ( Keff1(0, 2) + Keff1(2, 0) );
            result [ 5 ] = 0.5 * ( Keff1(0, 1) + Keff1(1, 0) );
        }
        return true;
    } else if ( code.compare("flux") == 0 || code.compare("stress") == 0 ) {
        result.resize( temp_stress.size() );
        for ( unsigned k = 0; k < temp_stress.size(); k++ ) {
            result [ k ] = temp_stress [ k ];
        }
        return true;
    } else {
        return RVEMaterialStatus :: giveValues(code, result);
    }
}


/////////////////////////////////./////////////////////////
Vector DiscreteTransportRVEMaterialStatus :: giveStress(const Vector &strain, double timeStep) {
    ( void ) timeStep;

    //precomputed material
    if ( is_precomputed ) {
        return giveStressPrecomputed(strain, timeStep);
    }

    cout << "Solving transport RVE" << endl;

    temp_strain  = strain;
    temp_strain = addEigenStrain(temp_strain); //macroscopic eigenstrain

    transformStrain();

    //set eigenstrains
    applyEigenStrains();

    //solve
    RVE->resetTime();
    RVE->giveSolver()->runBeforeEachStep();
    RVE->giveSolver()->solve();


    //collect results
    local_stress.resize(temp_strain.size() );
    local_stress.setZero();
    collectStresses();

    transformStress();

    return temp_stress;
}

/////////////////////////////////./////////////////////////
void DiscreteTransportRVEMaterialStatus :: setParameterValue(string code, double value) {
    if ( code.compare("pressure") == 0 ) {
        macro_pressure = value;
    } else {
        RVEMaterialStatus :: setParameterValue(code, value);
    }
}
/////////////////////////////////./////////////////////////
Vector DiscreteTransportRVEMaterialStatus :: giveStressWithFrozenIntVars(const Vector &strain, double timeStep) {
    ( void ) timeStep;

    temp_strain = strain;
    transformStrain();
    local_stress = giveStiffnessTensorLocal("elastic") * local_strain;
    transformStress();

    return temp_stress;
}

//////////////////////////////////////////////////////////
void DiscreteTransportRVEMaterialStatus :: update() {
    if ( !is_precomputed ) {
        RVEMaterialStatus :: update();
    }
}


//////////////////////////////////////////////////////////
void DiscreteTransportRVEMaterialStatus :: generateRandomFixedBC() {
    BCContainer *bconds = RVE->giveBC();
    FunctionContainer *funcs = RVE->giveFunctions();
    ConstraintContainer *constrs = RVE->giveConstraints();

    RVEMaterial *macromat = static_cast< RVEMaterial * >( mat );
    size_t ndim = macromat->giveDimension();

    //mechanics
    JointDoF *jd;
    Node *masternode;
    bool found = false;
    for ( unsigned j = 0; j < constrs->giveConstraintsSize(); j++ ) {
        if ( found ) {
            break;
        }
        jd = constrs->giveConstraint(j);
        for ( unsigned k = 0; k < jd->giveNumOfDoFMasters(); k++ ) {
            masternode = jd->giveMasterNode(0);
            if ( dynamic_cast< MechDoF * >( masternode ) == nullptr && dynamic_cast< TrsDoF * >( masternode ) == nullptr ) {
                BoundaryCondition *bc;
                vector< int >dBC, nBC;
                dBC.resize(masternode->giveNumberOfDoFs(), -1);     //todo: warning C4267: 'argument': conversion from 'size_t' to 'const _Ty', possible loss of data
                nBC.resize(masternode->giveNumberOfDoFs(), -1);
                for ( unsigned p = 0; p < min(dBC.size(), ndim); p++ ) {
                    dBC [ p ] = funcs->giveSize();                                         //only pressure and translations
                }
                bc = new BoundaryCondition(masternode, dBC, nBC);
                bconds->addBoundaryCondition(bc);

                found = true;
                break;
            }
        }
    }

    //add constant function
    vector< double >x, y;
    x.resize(1, 0);
    y.resize(1, 0);
    PieceWiseLinearFunction *newf = new PieceWiseLinearFunction(x, y);
    funcs->addFunction(newf);
}


//////////////////////////////////////////////////////////
void DiscreteTransportRVEMaterialStatus :: generateVolumetricAverageBC() {
    //this function applies volumetric constraint
    //it is not really wise to use it as it leads to relatively full stiffness matrix after transfer to constraint space

    DiscreteTransportRVEMaterial *macromat = static_cast< DiscreteTransportRVEMaterial * >( mat );
    unsigned ndim = macromat->giveDimension();

    NodeContainer *nodes = RVE->giveNodes();
    BCContainer *bconds = RVE->giveBC();
    ElementContainer *elems = RVE->giveElements();
    ConstraintContainer *constrs = RVE->giveConstraints();

    VolumetricAverage *va;
    unsigned funsize = RVE->giveFunctions()->giveSize();
    vector< double >x, y;
    x.resize(1, 0);
    y.resize(1, 0);
    volumAverFunc = new PieceWiseLinearFunction(x, y);
    RVE->giveFunctions()->addFunction(volumAverFunc);
    vector< unsigned >dirs;

    vector< Node * >vm;

    for ( unsigned n = 0; n < nodes->giveSize(); n++ ) {
        if ( nodes->giveNode(n)->doesMechanics() && ( dynamic_cast< MechDoF * >( nodes->giveNode(n) ) == nullptr && dynamic_cast< TrsDoF * >( nodes->giveNode(n) ) == nullptr ) ) {
            vm.push_back( nodes->giveNode(n) );
        }
    }
    if ( vm.size() > 0 ) {
        unsigned nDoFs = 3;
        if ( ndim == 3 ) {
            nDoFs = 6;
        }
        MechDoF *pn = new MechDoF(ndim, nDoFs);   //?? for transport
        nodes->addNode(pn);

        dirs.resize( vm.size() );

        for ( unsigned vi = 0; vi < nDoFs; vi++ ) {
            fill(dirs.begin(), dirs.end(), vi);
            va = new VolumetricAverage(vm, dirs, pn, vi, elems, constrs);
            constrs->addConstraint(va);
        }

        BoundaryCondition *bc;
        vector< int >dBC, nBC;
        dBC.resize(nDoFs, funsize);
        nBC.resize(nDoFs, -1);
        bc = new BoundaryCondition(pn, dBC, nBC);
        bconds->addBoundaryCondition(bc);
    }
}

//////////////////////////////////////////////////////////
Matrix DiscreteTransportRVEMaterialStatus :: giveStiffnessTensorLocal(string type) const {
    ( void ) type;
    if ( is_precomputed ) {
        return giveStiffnessTensorPrecomputedLocal(type);
    }

    unsigned strain_size = mat->giveStrainSize();
    Matrix Keff = Matrix :: Zero(strain_size, strain_size);

    ElementContainer *elems = RVE->giveElements();

    unsigned dimension = mat->giveDimension();
    double volume = 0;
    Point normal;
    Vector n = Vector :: Zero(dimension);
    DiscreteTrsprtElem *e;
    volume = 0;
    VectTrsprtMaterialStatus *tstat;

    for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
        e = static_cast< DiscreteTrsprtElem * >( elems->giveElement(i) );
        normal = e->giveNormal();
        for ( unsigned v = 0; v < dimension; v++ ) {
            n [ v ] = normal(v);
        }
        tstat = static_cast< VectTrsprtMaterialStatus * >( e->giveMatStatus(0) );
        Keff -= dyadicProduct(n, n) * ( e->giveLength() * e->giveArea() * tstat->giveEffectiveConductivity("tangent") );
        volume += e->giveVolume();
    }
    Keff /= volume;
    return Keff;
}

//////////////////////////////////////////////////////////
Matrix DiscreteTransportRVEMaterialStatus :: giveStiffnessTensorLocalExact(string type) {
    ( void ) type;
    if ( is_precomputed ) {
        return giveStiffnessTensorPrecomputedLocal(type);
    }

    //compute exactly
    cout << "Precomputing primary fields for transport RVE permeability" << endl;
    unsigned dimension = mat->giveDimension();
    Matrix Keff = Matrix :: Zero(dimension, dimension);
    Vector help_strain = Vector :: Zero(dimension);
    Vector help_stress;
    double factor = 1e-10;
    for ( unsigned i = 0; i < dimension; i++ ) {
        help_strain [ i ] = factor;
        help_stress = giveStress(help_strain, -1);
        help_strain [ i ] = 0.;
        for ( unsigned j = 0; j < dimension; j++ ) {
            Keff(i, j) = help_stress [ j ] / factor;
        }
    }
    return Keff;
}

//////////////////////////////////////////////////////////
Matrix DiscreteTransportRVEMaterialStatus :: giveStiffnessTensor(string type) const {
    Matrix Keff = giveStiffnessTensorLocal(type);
    return ( transf.transpose() * Keff ) * transf;
}


//////////////////////////////////////////////////////////
Matrix DiscreteTransportRVEMaterialStatus :: giveStiffnessTensorPrecomputedLocal(string type) const {
    ( void ) type;
    DiscreteTransportRVEMaterial *macromaterial = static_cast< DiscreteTransportRVEMaterial * >( mat );
    Matrix Keff = macromaterial->givePrecomputedConductivity() * temp_nonlin;
    return Keff;
}


//////////////////////////////////////////////////////////
Matrix DiscreteTransportRVEMaterialStatus :: giveDampingTensor() const {
    if ( is_precomputed ) {
        return giveDampingTensorPrecomputed();
    }

    ElementContainer *elems = RVE->giveElements();
    DiscreteTrsprtElem *e;
    double volume = 0;
    double mass = 0;
    for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
        e = dynamic_cast< DiscreteTrsprtElem * >( elems->giveElement(i) );
        if ( e ) {
            mass += e->giveVolume() * e->giveMatStatus(0)->giveDampingTensor()(0, 0);
            volume += e->giveVolume();
        }
    }
    Matrix m = Matrix :: Zero(1, 1);
    m(0, 0) = mass / volume;
    return m;
};

//////////////////////////////////////////////////////////
Matrix DiscreteTransportRVEMaterialStatus :: giveDampingTensorPrecomputed() const {
    DiscreteTransportRVEMaterial *macromaterial = static_cast< DiscreteTransportRVEMaterial * >( mat );
    Matrix M = Matrix :: Zero(1, 1);
    M(0, 0) = macromaterial->givePrecomputedCapacity();
    return M;
};

//////////////////////////////////////////////////////////
void DiscreteTransportRVEMaterialStatus :: setFromPrecomputedToFullModel() {
    if ( !is_precomputed ) {
        return;
    }
    is_precomputed = false;
    if ( is_master_status ) {
        return;                   //this was the master status already initialized
    }
    cout << "*** initialization of transport RVE ***" << endl;
    RVEMaterialStatus :: init();

    ElementContainer *elems = RVE->giveElements();
    orig_mater_params.resize(elems->giveSize() * 2); //initial permeability and a parameter
    DiscreteTrsprtElem *e;
    VectTrsprtMaterial *m;
    for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
        e = static_cast< DiscreteTrsprtElem * >( elems->giveElement(i) );
        m = static_cast< VectTrsprtMaterial * >( e->giveMatStatus(0)->giveMaterial() );
        orig_mater_params [ 2 * i  ] = m->givePermeability();
        orig_mater_params [ 2 * i + 1 ] = m->giveParamA();
    }
}

//////////////////////////////////////////////////////////
void DiscreteTransportRVEMaterialStatus :: init() {
    DiscreteTransportRVEMaterial *macromaterial = static_cast< DiscreteTransportRVEMaterial * >( mat );
    calculateTransformationMatrix();
    Matrix stiff = macromaterial->givePrecomputedConductivity();
    if ( macromaterial->givePrecomputedCapacity() < 0 ) {
        DiscreteTransportRVEMaterialStatus :: setFromPrecomputedToFullModel();
        Matrix c = DiscreteTransportRVEMaterialStatus :: giveDampingTensor();
        macromaterial->setPrecomputedCapacity(c(0, 0) );

        VectTrsprtMaterialStatus *status = static_cast< VectTrsprtMaterialStatus * >( RVE->giveElements()->giveElement(0)->giveMatStatus(0) );
        VectTrsprtMaterial *material = static_cast< VectTrsprtMaterial * >( RVE->giveElements()->giveElement(0)->giveMaterial() );
        macromaterial->setMasterMaterial(status, material);
        is_master_status = true;

        if ( macromaterial->shouldStartFromPrecomputed() && stiff.size() == 0 ) {
            unsigned ndim = macromaterial->giveDimension();
            Matrix Keff = Matrix :: Zero(ndim, ndim);
            if ( macromaterial->isElasticSolutionVoigt() ) {
                //use Voigt constraint
                Keff = giveStiffnessTensorLocal("elastic");
            } else {
                //compute exactly
                cout << "Precomputing primary fields for transport RVE" << endl;
                Vector help_strain = Vector :: Zero(ndim);
                Vector help_stress;
                double factor = 1e-10;
                for ( unsigned i = 0; i < ndim; i++ ) {
                    cout << "precomputing for pressure gradient component " << i << " out of " << ndim << endl;
                    help_strain [ i ] = factor;
                    help_stress = giveStress(help_strain, -1);
                    help_strain [ i ] = 0.;
                    for ( unsigned j = 0; j < ndim; j++ ) {
                        Keff(i, j) = help_stress [ j ] / factor;
                    }
                }
            }
            macromaterial->setPrecomputedConductivity(Keff);
            is_precomputed = true; //set back to precomputed to use it
        }
    }
    if ( !macromaterial->shouldStartFromPrecomputed() ) {
        DiscreteTransportRVEMaterialStatus :: setFromPrecomputedToFullModel();
    }
}

//////////////////////////////////////////////////////////
void DiscreteTransportRVEMaterialStatus :: transformStrain() {
    local_strain = transf * temp_strain;
}

//////////////////////////////////////////////////////////
void DiscreteTransportRVEMaterialStatus :: transformStress() {
    temp_stress = transf.transpose() * local_stress;
}

//////////////////////////////////////////////////////////
void DiscreteTransportRVEMaterialStatus :: calculateTransformationMatrix() {
    transf = axDirs;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE TRANSPORT RVE MATERIAL
//////////////////////////////////////////////////////////
DiscreteTransportRVEMaterial :: DiscreteTransportRVEMaterial(unsigned dimension) : RVEMaterial(dimension) {
    name = "transport RVE material";
    conductivity = Matrix(0, 0);
    start_from_precomputed = true;
    capacity = -1;
    strainsize = dim;
}

//////////////////////////////////////////////////////////
MaterialStatus *DiscreteTransportRVEMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    DiscreteTransportRVEMaterialStatus *newstat = new DiscreteTransportRVEMaterialStatus(this, e, ipnum, inputfile, dim);
    return newstat;
}

//////////////////////////////////////////////////////////
void DiscreteTransportRVEMaterial :: setPrecomputedCapacity(double c) {
    capacity = c;
}

//////////////////////////////////////////////////////////
void DiscreteTransportRVEMaterial :: setPrecomputedConductivity(Matrix lam) {
    conductivity = lam;

    /*
     * const char *buffer = "precomputed_conductivity.out";
     * ofstream outputfile( ( masterModel->giveResultDirectory() / buffer ).string() );
     * if ( outputfile.is_open() ) {
     *  outputfile << std :: scientific;
     *  outputfile.precision(10);
     *
     *  unsigned size = conductivity.rows();
     *  for ( unsigned i = 0; i < size; i++ ) {
     *      for ( unsigned j = 0; j < size; j++ ) {
     *          if ( j > 0 ) {
     *              outputfile << "\t";
     *          }
     *          outputfile << conductivity(i, j);
     *      }
     *      outputfile << endl;
     *  }
     *  outputfile.close();
     * }
     */
}

//////////////////////////////////////////////////////////
void DiscreteTransportRVEMaterial :: setMasterMaterial(VectTrsprtMaterialStatus *masterS, VectTrsprtMaterial *masterM) {
    masterStatus = masterS;
    masterMaterial = masterM;

    /*
     * const char *buffer = "precomputed_conductivity.out";
     * ofstream outputfile( ( masterModel->giveResultDirectory() / buffer ).string() );
     * if ( outputfile.is_open() ) {
     *  outputfile << std :: scientific;
     *  outputfile.precision(10);
     *
     *  unsigned size = conductivity.rows();
     *  for ( unsigned i = 0; i < size; i++ ) {
     *      for ( unsigned j = 0; j < size; j++ ) {
     *          if ( j > 0 ) {
     *              outputfile << "\t";
     *          }
     *          outputfile << conductivity(i, j);
     *      }
     *      outputfile << endl;
     *  }
     *  outputfile.close();
     * }
     */
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE MECHANICAL RVE MATERIAL
//////////////////////////////////////////////////////////
DiscreteMechanicalRVEMaterialStatus :: DiscreteMechanicalRVEMaterialStatus(RVEMaterial *m, Element *e, unsigned ipnum, fs :: path masterfile, unsigned ndim) : DiscreteTransportRVEMaterialStatus(m, e, ipnum, masterfile, ndim) {
    name = "mechanical RVE mat. status";
    is_precomputed = true;
    is_master_status = false;
    macro_volumetricStrain = 0;
}

/////////////////////////////////./////////////////////////
void DiscreteMechanicalRVEMaterialStatus :: applyEigenStrains() {
    //set eigenstrains
    DiscreteMechanicalRVEMaterial *macromat = static_cast< DiscreteMechanicalRVEMaterial * >( mat );
    vector< vector< Matrix > > *projectors = macromat->giveProjectors();
    ElementContainer *elems = RVE->giveElements();
    Element *e;
    for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
        if ( ( ( * projectors ) [ i ] ).size() == 0 ) {
            continue;
        }
        e = elems->giveElement(i);
        for ( unsigned ip = 0; ip < e->giveNumIP(); ip++ ) {
            Vector eigstr =  -( * projectors ) [ i ] [ ip ] * local_strain;
            e->giveMatStatus(ip)->setEigenStrain(eigstr);
        }
    }
}

/////////////////////////////////./////////////////////////
void DiscreteMechanicalRVEMaterialStatus :: collectStresses() {
    DiscreteMechanicalRVEMaterial *macromat = static_cast< DiscreteMechanicalRVEMaterial * >( mat );
    unsigned ndim = macromat->giveDimension();
    vector< vector< Matrix > > *projectors = macromat->giveProjectors();
    vector< bool > *is_elem_discrete = macromat->giveElemDiscreteness();
    ElementContainer *elems = RVE->giveElements();
    double volume = 0.;
    Element *e;
    local_stress *= 0;
    double multipl = 1;
    for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
        if ( ( ( * projectors ) [ i ] ).size() == 0 ) {
            continue;
        }
        e = elems->giveElement(i);
        if ( ( * is_elem_discrete ) [ i ] ) {
            multipl = ndim;
        } else {
            multipl = 1;
        }
        for ( unsigned ip = 0; ip < e->giveNumIP(); ip++ ) {
            local_stress += multipl * e->giveIPVolume(ip) * ( ( ( * projectors ) [ i ] [ ip ] ).transpose() * e->giveMatStatus(ip)->giveTempStress() );
        }
        volume += e->giveVolume();
    }
    local_stress /= volume;
}


/////////////////////////////////./////////////////////////
Vector DiscreteMechanicalRVEMaterialStatus :: giveStressPrecomputed(const Vector &strain, double timeStep) {
    ( void ) timeStep;
    DiscreteMechanicalRVEMaterial *macromaterial = static_cast< DiscreteMechanicalRVEMaterial * >( mat );
    temp_strain = addEigenStrain(macromaterial->strainToCosserat(strain) );  //macroscopic eigenstrain
    transformStrain();
    local_stress = macromaterial->givePrecomputedElasticTensor() * local_strain;
    transformStress();
    if ( macromaterial->isNonlinear() && checkOttosenCriterion() ) {
        setFromPrecomputedToFullModel();
        return macromaterial->stressToCauchy(giveStress(strain, timeStep) );
    }
    return macromaterial->stressToCauchy(temp_stress);
}


//////////////////////////////////////////////////////////
bool DiscreteMechanicalRVEMaterialStatus ::  giveValues(string code, Vector &result) const {
    if ( code.compare("stress") == 0 ) {
        result.resize( temp_stress.size() );
        for ( unsigned k = 0; k < temp_stress.size(); k++ ) {
            result [ k ] = temp_stress [ k ];
        }
        return true;
    } else {
        return RVEMaterialStatus :: giveValues(code, result);
    }
}

/////////////////////////////////./////////////////////////
Vector DiscreteMechanicalRVEMaterialStatus :: giveStressWithFrozenIntVars(const Vector &strain, double timeStep) {
    ( void ) timeStep;
    DiscreteMechanicalRVEMaterial *macromaterial = static_cast< DiscreteMechanicalRVEMaterial * >( mat );
    temp_strain = macromaterial->strainToCosserat(strain);
    transformStrain();
    local_stress = giveStiffnessTensorLocal("elastic") * local_strain;
    transformStress();
    return macromaterial->stressToCauchy(temp_stress);
}

/////////////////////////////////./////////////////////////
Vector DiscreteMechanicalRVEMaterialStatus :: giveStress(const Vector &strain, double timeStep) {
    //delete curvatutures according to an updated homogenization theory
    ( void ) timeStep;
    //precomputed material
    if ( is_precomputed ) {
        return giveStressPrecomputed(strain, timeStep);
    }
    cout << "Solving mechanical RVE" << endl;
    DiscreteMechanicalRVEMaterial *macromaterial = static_cast< DiscreteMechanicalRVEMaterial * >( mat );
    temp_strain = addEigenStrain(macromaterial->strainToCosserat(strain) );  //macroscopic eigenstrain

    transformStrain();

    //set eigenstrains
    applyEigenStrains();

    //solve
    RVE->resetTime();
    RVE->giveSolver()->runBeforeEachStep();
    RVE->giveSolver()->solve();


    //collect results
    local_stress.resize( temp_strain.size() );
    local_stress.setZero();
    collectStresses();

    transformStress();

    return macromaterial->stressToCauchy(temp_stress);
}

/////////////////////////////////./////////////////////////
bool DiscreteMechanicalRVEMaterialStatus :: checkOttosenCriterion() {

    if (mat->giveDimension() != 3){
        cerr << "Ottosen criterion implemented only for 3 dimensions" << endl;
        exit(1);
    }

    vector< Vector >eigvecs;
    Vector eignums;

    //from paper Adaptive multiscale homogenization of the lattice discrete particle model for the analysis of damage and fracture in concrete, Roozbeh Rezakhani and Xinwei Zhou and Gianluca Cusatis
    double fc =  20e6 / 2.;
    double ft = 2.35e6 / 2.;
    double c1 = 4.775 / fc;
    double c2 = 7.048 / fc;
    double c3 = 0.88 / pow(fc, 2);

    Vector stress_without_couple = Vector :: Zero(9);
    for ( unsigned p = 0; p < 9; p++ ) {
        stress_without_couple [ p ] = temp_stress [ p ];
    }
    LinalgEigenSolver(stress_without_couple, eignums, eigvecs);
    double I1 = eignums [ 0 ] + eignums [ 1 ] + eignums [ 2 ];
    if ( abs(I1) < 1e-10 ) {
        return false;
    }

    double J2 = ( pow(eignums [ 0 ] - eignums [ 1 ], 2) + pow(eignums [ 0 ] - eignums [ 2 ], 2) + pow(eignums [ 1 ] - eignums [ 2 ], 2) ) / 6;
    double J3 = ( eignums [ 0 ] - I1 / 3. ) * ( eignums [ 1 ] - I1 / 3. ) * ( eignums [ 2 ] - I1 / 3. );
    double xi = I1 / sqrt(3);
    double rho = sqrt(2 * J2);
    double cos3theta = ( 3. * sqrt(3.) * J3 ) / ( 2. * pow(J2, 3. / 2.) );
    double K = 1 - 6.8 * pow(ft / fc - 0.07, 2);
    double r;
    if ( cos3theta > 0 ) {
        r = cos(acos(K * cos3theta) / 3);
    } else {
        r = cos(M_PI / 3 - acos(-K * cos3theta) / 3);
    }

    double otto = c1 * xi + c2 * rho * r + c3 * pow(rho, 2) - 1;
    return otto > 0;
}

/////////////////////////////////./////////////////////////
double DiscreteMechanicalRVEMaterialStatus :: giveCrackVolume() const {
    if ( is_precomputed ) {
        return 0;
    }
    double crackVolume = 0.;
    RigidBodyContact *e;
    ElementContainer *elems = RVE->giveElements();
    Vector res;
    for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
        e = dynamic_cast< RigidBodyContact * >( elems->giveElement(i) );
        if ( e ) {
            e->giveIPValues("tempCrackOpening", 0, res);
            if ( res.size() > 0 ) {
                crackVolume += res [ 0 ] * e->giveArea();
            }
        }
    }

    return crackVolume;
}

//////////////////////////////////////////////////////////
Matrix DiscreteMechanicalRVEMaterialStatus :: giveStiffnessTensorLocal(string type) const {
    DiscreteMechanicalRVEMaterial *macromat = static_cast< DiscreteMechanicalRVEMaterial * >( mat );
    if ( is_precomputed || macromat->hasPrecomputedTensorsStored()) {
        return macromat->givePrecomputedElasticTensor();
    } else {
        unsigned strain_size = mat->giveStrainSize();
        unsigned dim = mat->giveDimension();
        Matrix Keff = Matrix :: Zero(strain_size, strain_size);
        vector< vector< Matrix > > *projectors = macromat->giveProjectors();
        ElementContainer *elems = RVE->giveElements();
        double volume = 0;
        Point normal;
        Vector n = Vector :: Zero(dim);
        Element *e;
        volume = 0;
        Matrix stiff;
        for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
            if ( ( ( * projectors ) [ i ] ).size() == 0 ) {
                continue;
            }
            e = elems->giveElement(i);
            for ( unsigned ip = 0; ip < e->giveNumIP(); ip++ ) {
                stiff  = e->giveMatStatus(ip)->giveStiffnessTensor(type);
                Keff += e->giveIPVolume(ip) *  ( ( * projectors ) [ i ] [ ip ] ).transpose() * stiff * ( * projectors ) [ i ] [ ip ];
            }
            volume += e->giveVolume();
        }
        Keff /= volume;
        return Keff;
    }
}

//////////////////////////////////////////////////////////
Matrix DiscreteMechanicalRVEMaterialStatus :: giveStiffnessTensor(string type) const {
    Matrix Keff = giveStiffnessTensorLocal(type);
    DiscreteMechanicalRVEMaterial *macromaterial = static_cast< DiscreteMechanicalRVEMaterial * >( mat );
    return macromaterial->matrixToCauchy( ( transf.transpose() * Keff ) * transf);
}

//////////////////////////////////////////////////////////
Matrix DiscreteMechanicalRVEMaterialStatus :: giveDampingTensor() const {
    DiscreteMechanicalRVEMaterial *macromat = static_cast< DiscreteMechanicalRVEMaterial * >( mat );
    if ( is_precomputed ) {
        return macromat->givePrecomputedDampingTensor();
    } else {
        /*
         * unsigned ndim = macromat->giveDimension();
         * unsigned dof_size = 3 * ( ndim - 1 );
         * Matrix Keff = Matrix :: Zero(dof_size, dof_size);
         * for ( unsigned i = 0; i < dof_size; i++ ) {
         *  Keff(i, i) = 1e-20;
         * }
         * return Keff;
         */
        return Matrix :: Zero(1, 1);
    }
}

//////////////////////////////////////////////////////////
Matrix DiscreteMechanicalRVEMaterialStatus :: giveInertiaTensor() const {
    DiscreteMechanicalRVEMaterial *macromat = static_cast< DiscreteMechanicalRVEMaterial * >( mat );
    if ( is_precomputed ) {
        return macromat->givePrecomputedInertiaTensor();
    } else {
        unsigned ndim = macromat->giveDimension();
        unsigned dof_size = 3 * ( ndim - 1 );
        Matrix Keff = Matrix :: Zero(dof_size, dof_size);
        return Keff;
    }
}
//////////////////////////////////////////////////////////
Point DiscreteMechanicalRVEMaterialStatus :: calculateCentroid() {  //only for mechanics
    Point centroid = Point(0., 0., 0.);

    //find centroid
    ElementContainer *elems = RVE->giveElements();

    //mechanics
    Element *e;
    VectMechMaterial *vms;
    TensMechMaterial *tms;
    double density = 0;
    double weight;
    double totalweight = 0;
    for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
        e = elems->giveElement(i);
        vms = dynamic_cast< VectMechMaterial * >( e->giveMaterial() );
        if ( vms ) {
            density = vms->giveDensity();
        } else {
            tms = dynamic_cast< TensMechMaterial * >( e->giveMaterial() );
            if ( tms ) {
                density = tms->giveDensity();
            } else {
                cout << "Mechanical RVE contains element with material not derived from MechanicalMaterial" << endl;
                exit(1);
            }
        }
        for ( unsigned k = 0; k < e->giveNumIP(); k++ ) {
            weight = e->giveIPVolume(k) * density;
            totalweight += weight;
            centroid += e->giveIPLoc(k) * weight;
        }
    }
    centroid /= totalweight;
    return centroid;
}

//////////////////////////////////////////////////////////
double DiscreteMechanicalRVEMaterialStatus :: computeAverageDensity() const {
    ElementContainer *elems = RVE->giveElements();

    //mechanics
    Element *e;
    VectMechMaterial *vms;
    TensMechMaterial *tms;
    double density = 0;
    double weight = 0;
    double volume = 0;

    for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
        e = elems->giveElement(i);
        vms = dynamic_cast< VectMechMaterial * >( e->giveMaterial() );
        if ( vms ) {
            density = vms->giveDensity();
        } else {
            tms = dynamic_cast< TensMechMaterial * >( e->giveMaterial() );
            if ( tms ) {
                density = tms->giveDensity();
            } else {
                cout << "Mechanical RVE contains element with material not derived from MechanicalMaterial" << endl;
                exit(1);
            }
        }
        for ( unsigned k = 0; k < e->giveNumIP(); k++ ) {
            weight += e->giveIPWeight(k) * density;
            volume += e->giveVolume();
        }
    }
    return weight / volume;
}


//////////////////////////////////////////////////////////
vector< vector< Matrix > >DiscreteMechanicalRVEMaterialStatus :: calculateProjectors(const Point centroid) {   //only for mechanics
    ElementContainer *elems = RVE->giveElements();
    Material *emat;
    Point normal;
    Element *e;
    vector< vector< Matrix > >projectors( elems->giveSize() );
    for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
        e = elems->giveElement(i);
        emat = e->giveMaterial();
        if ( dynamic_cast< VectMechMaterial * >( emat ) ) {
            projectors [ i ] = calculateVectProjector(e, centroid);
        } else if ( dynamic_cast< TensMechMaterial * >( emat ) ) {
            projectors [ i ] = calculateTensProjector(e, centroid);
        }
    }
    return projectors;
}

//////////////////////////////////////////////////////////
vector< Matrix >DiscreteMechanicalRVEMaterialStatus :: calculateVectProjector(const Element *e, const Point centroid) {    //only for mechanics
    vector< Matrix >projector;

    DiscreteMechanicalRVEMaterial *macromat = static_cast< DiscreteMechanicalRVEMaterial * >( mat );
    unsigned ndim = macromat->giveDimension();
    unsigned projNum = e->giveMaterial()->giveStrainSize();
    unsigned strain_size = macromat->giveStrainSize();
    projector.resize( e->giveNumIP() );


    const RigidBodyContact *rbc;
    const LDPMTetra *tet;
    Point xc;
    Point alphaVec = Point(0, 0, 0);
    Point normal, t1, t2;

    rbc = dynamic_cast< const RigidBodyContact * >( e );
    tet = dynamic_cast< const LDPMTetra * >( e );

    for ( unsigned ip = 0; ip < e->giveNumIP(); ip++ ) {
        if ( rbc ) {
            normal = rbc->giveNormal();
            t1 = rbc->giveT1();
            t2 = rbc->giveT2();
        } else if ( tet ) {
            normal = tet->giveNormal(ip);
            t1 = tet->giveT1(ip);
            t2 = tet->giveT2(ip);
        } else {
            continue;
        }


        projector [ ip ] = Matrix :: Zero(projNum, strain_size);
        xc =  e->giveIPLoc(ip) - centroid; //here we make corrections to have origin of the reference system at the centroid
        unsigned v = 0;
        //displacements
        for ( ; v < ndim; v++ ) {
            if ( v == 0 ) {
                alphaVec = normal;
            } else if ( v == 1 ) {
                alphaVec = t1;
            } else if ( v == 2 ) {
                alphaVec = t2;
            }
            if ( ndim == 2 ) {
                projector [ ip ](v, 0) =  normal.x() * alphaVec.x();
                projector [ ip ](v, 1) =  normal.y() * alphaVec.y();
                projector [ ip ](v, 2) =  normal.x() * alphaVec.y();
                projector [ ip ](v, 3) =  normal.y() * alphaVec.x();
                //in paper EliCus22, this projection is not active
                double factor2D = 0;
                if ( macromat->projectCurvature() ) {
                    factor2D = xc.x() * alphaVec.y() - xc.y() * alphaVec.x();
                    projector [ ip ](v, 4) =  factor2D * normal.x();
                    projector [ ip ](v, 5) =  factor2D * normal.y();
                }
            } else {
                projector [ ip ](v, 0) =  normal.x() * alphaVec.x();
                projector [ ip ](v, 1) =  normal.y() * alphaVec.y();
                projector [ ip ](v, 2) =  normal.z() * alphaVec.z();
                projector [ ip ](v, 3) =  normal.z() * alphaVec.y();
                projector [ ip ](v, 4) =  normal.y() * alphaVec.z();
                projector [ ip ](v, 5) =  normal.z() * alphaVec.x();
                projector [ ip ](v, 6) =  normal.x() * alphaVec.z();
                projector [ ip ](v, 7) =  normal.y() * alphaVec.x();
                projector [ ip ](v, 8) =  normal.x() * alphaVec.y();
                //in paper EliCus22, this projection is not active
                Point factor3D = Point(0, 0, 0);
                if ( macromat->projectCurvature() ) {
                    factor3D = xc.cross(alphaVec);
                    projector [ ip ](v, 9)  = normal.x() * factor3D.x();
                    projector [ ip ](v, 10) =  normal.y() * factor3D.y();
                    projector [ ip ](v, 11) =  normal.z() * factor3D.z();
                    projector [ ip ](v, 12) =  normal.z() * factor3D.y();
                    projector [ ip ](v, 13) =  normal.y() * factor3D.z();
                    projector [ ip ](v, 14) =  normal.z() * factor3D.x();
                    projector [ ip ](v, 15) =  normal.x() * factor3D.z();
                    projector [ ip ](v, 16) =  normal.y() * factor3D.x();
                    projector [ ip ](v, 17) =  normal.x() * factor3D.y();
                }
            }
        }
        //rotations
        for ( ; v < projNum; v++ ) {
            if ( v - ndim == 0 ) {
                alphaVec = normal;
            } else if ( v - ndim == 1 ) {
                alphaVec = t1;
            } else if ( v - ndim == 2 ) {
                alphaVec = t2;
            }

            if ( ndim == 2 ) {
                projector [ ip ](v, 4) =  normal.x();
                projector [ ip ](v, 5) =  normal.y();
            } else {
                projector [ ip ](v, 9) =  normal.x() * alphaVec.x();
                projector [ ip ](v, 10) =  normal.y() * alphaVec.y();
                projector [ ip ](v, 11) =  normal.z() * alphaVec.z();
                projector [ ip ](v, 12) =  normal.z() * alphaVec.y();
                projector [ ip ](v, 13) =  normal.y() * alphaVec.z();
                projector [ ip ](v, 14) =  normal.z() * alphaVec.x();
                projector [ ip ](v, 15) =  normal.x() * alphaVec.z();
                projector [ ip ](v, 16) =  normal.y() * alphaVec.x();
                projector [ ip ](v, 17) =  normal.x() * alphaVec.y();
            }
        }
    }
    return projector;
}

//////////////////////////////////////////////////////////
vector< Matrix >DiscreteMechanicalRVEMaterialStatus :: calculateTensProjector(const Element *e, const Point centroid) {    //only for mechanics
    vector< Matrix >projector;

    DiscreteMechanicalRVEMaterial *macromat = static_cast< DiscreteMechanicalRVEMaterial * >( mat );
    unsigned ndim = macromat->giveDimension();
    unsigned projNum = e->giveMaterial()->giveStrainSize();
    unsigned strain_size = macromat->giveStrainSize();
    projector.resize( e->giveNumIP() );

    Point xc;
    for ( unsigned ip = 0; ip < e->giveNumIP(); ip++ ) {
        projector [ ip ] = Matrix :: Zero(projNum, strain_size);
        xc =  e->giveIPLoc(ip) - centroid; //here we make corrections to have origin of the reference system at the centroid
        unsigned v = 0;
        //displacements
        for ( ; v < ndim * ndim; v++ ) {
            projector [ ip ](v, v) =  1.;
        }
        //in paper EliCus22, this projection is not active
        if ( macromat->projectCurvature() ) {
            if ( ndim == 2 ) {
                projector [ ip ](0, 4) = -xc.y();
                projector [ ip ](1, 5) = xc.x();
                //projector[ip](2,4 ) = xc.y();
                projector [ ip ](2, 5) = xc.y();
                projector [ ip ](3, 4) = -xc.x();
                //projector[ip](3,5 ) = xc.y();
            } else {
                projector [ ip ](0, 17) = xc.z();
                projector [ ip ](0, 15) = -xc.y();
                projector [ ip ](1, 13) = xc.x();
                projector [ ip ](1, 16) = -xc.z();
                projector [ ip ](2, 14) = xc.y();
                projector [ ip ](2, 12) = -xc.x();

                projector [ ip ](3, 11) = -xc.x();
                projector [ ip ](3, 14) =  xc.z();
                projector [ ip ](4, 16) = -xc.y();
                projector [ ip ](4, 10) =  xc.x();
                projector [ ip ](5, 12) = -xc.z();
                projector [ ip ](5, 11) =  xc.y();
                projector [ ip ](6, 9) = -xc.y();
                projector [ ip ](6, 17) =  xc.x();
                projector [ ip ](7, 10) = -xc.z();
                projector [ ip ](7, 13) =  xc.y();
                projector [ ip ](8, 15) = -xc.x();
                projector [ ip ](8, 9) =  xc.z();
            }
        }

        //rotations
        for ( ; v < projNum; v++ ) {
            projector [ ip ](v, v) =  1;
        }
    }
    return projector;
}


//////////////////////////////////////////////////////////
void DiscreteMechanicalRVEMaterialStatus :: setFromPrecomputedToFullModel() {
    if ( !is_precomputed ) {
        return;
    }
    is_precomputed = false;
    if ( is_master_status ) {
        return;                   //this was the master status already initialized
    }
    cout << "*** initialization of mechanical RVE ***" << endl;
    RVEMaterialStatus :: init();
}

//////////////////////////////////////////////////////////
void DiscreteMechanicalRVEMaterialStatus :: init() {
    DiscreteMechanicalRVEMaterial *macromaterial = static_cast< DiscreteMechanicalRVEMaterial * >( mat );
    calculateTransformationMatrix();
    Matrix D = macromaterial->givePrecomputedElasticTensor();
    Matrix I = macromaterial->givePrecomputedInertiaTensor();

    if ( I.size() == 0 ) {
        setFromPrecomputedToFullModel();
        Point centroid = calculateCentroid();
        vector< bool >is_discrete = calculateElemDiscreteness();
        macromaterial->setElemDiscreteness(is_discrete);
        vector< vector< Matrix > >projectors = calculateProjectors(centroid);
        macromaterial->setCentroidAndProjectors(centroid, projectors);

        macromaterial->setPrecomputedDampingTensor( giveDampingTensor() );
        macromaterial->setPrecomputedInertiaTensor( giveInertiaTensor() );

        if ( macromaterial->shouldStartFromPrecomputed() && D.size() == 0 ) {
            is_master_status = true;
            unsigned strain_size = macromaterial->giveStrainSize();
            Matrix Keff = Matrix :: Zero(strain_size, strain_size);
            if ( macromaterial->isElasticSolutionVoigt() ) {
                //use Voigt constraint
                Keff = giveStiffnessTensorLocal("secant");
            } else {
                //compute exactly
                cout << "Precomputing primary fields for mechanical RVE" << endl;
                bool previously_cauchy = macromaterial->isConvertingFromCauchy();
                macromaterial->setConversionFromCauchy(false);
                Vector help_strain = Vector :: Zero(strain_size);
                Vector help_stress;
                double factor = 1e-10;
                for ( unsigned i = 0; i < strain_size; i++ ) {
                    cout << "precomputing for strain component " << i << " out of " << strain_size << endl;
                    help_strain [ i ] = factor;     
                    help_stress = giveStress(help_strain, -1);
                    help_strain [ i ] = 0.;
                    for ( unsigned j = 0; j < strain_size; j++ ) {
                        Keff(i, j) = help_stress [ j ] / factor;
                    }
                }
                //cout << Keff << endl;
                //exit(1);
                macromaterial->setConversionFromCauchy(previously_cauchy);
            }
            macromaterial->setPrecomputedElasticTensor(Keff);
            macromaterial->setAverageDensity(computeAverageDensity() );
            is_precomputed = true; //set back to precomputed to use it
        }
    }
    if ( !macromaterial->shouldStartFromPrecomputed() ) {
        DiscreteMechanicalRVEMaterialStatus :: setFromPrecomputedToFullModel();
    }
}


//////////////////////////////////////////////////////////
void DiscreteMechanicalRVEMaterialStatus :: transformStrain() {
    local_strain = transf * temp_strain;
    //DiscreteMechanicalRVEMaterial *macromat = static_cast< DiscreteMechanicalRVEMaterial * >( mat );
    //unsigned ndim = macromat->giveDimension();
    //for ( unsigned i = ndim * ndim; i < local_strain.size(); i++ ) {
    //    local_strain [ i ] = 0.;
    //}
}

//////////////////////////////////////////////////////////
void DiscreteMechanicalRVEMaterialStatus :: transformStress() {
    temp_stress = transf.transpose() * local_stress;
}

//////////////////////////////////////////////////////////
void giveTensorCoords3D(unsigned pos, unsigned *i, unsigned *j) {
    switch ( pos ) {
    case 0:
        * i = 0;
        * j = 0;
        break;
    case 1:
        * i = 1;
        * j = 1;
        break;
    case 2:
        * i = 2;
        * j = 2;
        break;
    case 3:
        * i = 2;
        * j = 1;
        break;
    case 4:
        * i = 1;
        * j = 2;
        break;
    case 5:
        * i = 2;
        * j = 0;
        break;
    case 6:
        * i = 0;
        * j = 2;
        break;
    case 7:
        * i = 1;
        * j = 0;
        break;
    case 8:
        * i = 0;
        * j = 1;
        break;
    default:
        cerr << "giveTensorCoords should never go here" << endl;
        exit(1);
    }
}

//////////////////////////////////////////////////////////
void giveVectorCoords2D(unsigned *pos, unsigned i, unsigned j) {
    if ( i == 0 ) {
        if ( j == 0 ) {
            * pos = 0;
        } else if ( j == 1 ) {
            * pos = 3;
        } else {
            cerr << "giveTensorCoords should never go here" << endl;
            exit(1);
        }
    } else if ( i == 1 ) {
        if ( j == 0 ) {
            * pos = 2;
        } else if ( j == 1 ) {
            * pos = 1;
        } else {
            cerr << "giveTensorCoords should never go here" << endl;
            exit(1);
        }
    } else {
        cerr << "giveTensorCoords should never go here" << endl;
        exit(1);
    }
}

//////////////////////////////////////////////////////////
void giveVectorCoords3D(unsigned *pos, unsigned i, unsigned j) {
    if ( i == 0 ) {
        if ( j == 0 ) {
            * pos = 0;
        } else if ( j == 1 ) {
            * pos = 8;
        } else if ( j == 2 ) {
            * pos = 6;
        } else {
            cerr << "giveTensorCoords should never go here" << endl;
            exit(1);
        }
    } else if ( i == 1 ) {
        if ( j == 0 ) {
            * pos = 7;
        } else if ( j == 1 ) {
            * pos = 1;
        } else if ( j == 2 ) {
            * pos = 4;
        } else {
            cerr << "giveTensorCoords should never go here" << endl;
            exit(1);
        }
    } else if ( i == 2 ) {
        if ( j == 0 ) {
            * pos = 5;
        } else if ( j == 1 ) {
            * pos = 3;
        } else if ( j == 2 ) {
            * pos = 2;
        } else {
            cerr << "giveTensorCoords should never go here" << endl;
            exit(1);
        }
    } else {
        cerr << "giveTensorCoords should never go here" << endl;
        exit(1);
    }
}
//////////////////////////////////////////////////////////
void DiscreteMechanicalRVEMaterialStatus :: calculateTransformationMatrix() {
    DiscreteMechanicalRVEMaterial *macromat = static_cast< DiscreteMechanicalRVEMaterial * >( mat );
    unsigned ndim = macromat->giveDimension();

    if ( ndim == 2 ) {
        transf = Matrix :: Zero(6, 6);
        unsigned posA, posB;
        unsigned iA, jA, iB, jB;
        for ( iA = 0; iA < 2; iA++ ) {
            for ( jA = 0; jA < 2; jA++ ) {
                for ( iB = 0; iB < 2; iB++ ) {
                    for ( jB = 0; jB < 2; jB++ ) {
                        giveVectorCoords2D(& posA, iA, iB);
                        giveVectorCoords2D(& posB, jA, jB);
                        transf(posA, posB) = axDirs(iA, jA) * axDirs(iB, jB);
                    }
                }
                transf(iA + 4, jA + 4) = axDirs(iA, jA); //might be transposed, nt sure. Check!
            }
        }
    } else {
        transf = Matrix :: Zero(18, 18);
        unsigned posA, posB;
        unsigned iA, jA, iB, jB;
        for ( iA = 0; iA < 3; iA++ ) {
            for ( jA = 0; jA < 3; jA++ ) {
                for ( iB = 0; iB < 3; iB++ ) {
                    for ( jB = 0; jB < 3; jB++ ) {
                        giveVectorCoords3D(& posA, iA, iB);
                        giveVectorCoords3D(& posB, jA, jB);
                        transf(posA, posB) = transf(posA + 9, posB + 9) = axDirs(iA, jA) * axDirs(iB, jB);
                    }
                }
            }
        }
    }
}

//////////////////////////////////////////////////////////
Matrix DiscreteMechanicalRVEMaterialStatus :: giveMassTensor() const {
    DiscreteMechanicalRVEMaterial *macromat = static_cast< DiscreteMechanicalRVEMaterial * >( mat );
    Matrix M = Matrix :: Zero(macromat->giveStrainSize(), macromat->giveStrainSize() );
    for ( unsigned i = 0; i < macromat->giveDimension(); i++ ) {
        M(i, i) = macromat->giveAverageDensity();
    }
    return M;
};

//////////////////////////////////////////////////////////
void DiscreteMechanicalRVEMaterialStatus :: setParameterValue(string code, double value) {
    if ( code.compare("volumetricStrain") == 0 ) {
        macro_volumetricStrain = value;  //this is unfortunately not active at elements
    } else {
        MaterialStatus :: setParameterValue(code, value);
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE MECHANICAL RVE MATERIAL
//////////////////////////////////////////////////////////
DiscreteMechanicalRVEMaterial :: DiscreteMechanicalRVEMaterial(unsigned dimension) : DiscreteTransportRVEMaterial(dimension) {
    name = "mechanical RVE material";
    precompElastic = Matrix(0, 0);
    start_from_precomputed = true;
    strainsize = dim * dim;
    strainsize += ( dim == 3 ) ? dim * dim : dim; //in 2D only vector
    project_curvature = false;
    average_density = 0.;
    convert_from_cauchy = false;
    CauchyToCosseratMatrix = Matrix :: Zero(0, 0);
    storedPrecomputeTensors = false;
}

/////////////////////////////////./////////////////////////
Vector DiscreteMechanicalRVEMaterial :: strainToCosserat(Vector strain) {
    if ( convert_from_cauchy ) {
        return CauchyToCosseratMatrix * strain;
    } else {
        return strain;
    }
}

/////////////////////////////////./////////////////////////
Vector DiscreteMechanicalRVEMaterial :: stressToCauchy(Vector stress) {
    if ( convert_from_cauchy ) {
        return CauchyToCosseratMatrix.transpose() * stress;
    } else {
        return stress;
    }
}

/////////////////////////////////./////////////////////////
Matrix DiscreteMechanicalRVEMaterial :: matrixToCauchy(Matrix matrix) {
    if ( convert_from_cauchy ) {
        return CauchyToCosseratMatrix.transpose() * matrix * CauchyToCosseratMatrix;
    } else {
        return matrix;
    }
}

//////////////////////////////////////////////////////////
MaterialStatus *DiscreteMechanicalRVEMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    DiscreteMechanicalRVEMaterialStatus *newstat = new DiscreteMechanicalRVEMaterialStatus(this, e, ipnum, inputfile, dim);
    return newstat;
}

//////////////////////////////////////////////////////////
void DiscreteMechanicalRVEMaterial :: setCentroidAndProjectors(Point c, vector< vector< Matrix > >p) {
    centroid = c;
    projectors = p;
}

//////////////////////////////////////////////////////////
void DiscreteMechanicalRVEMaterial :: readFromLine(istringstream &iss) {
    RVEMaterial :: readFromLine(iss);

    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream

    string param;
    while (  iss >> param ) {
        if ( param.compare("project_curvatures") == 0 || param.compare("project_curvature") == 0 ) {
            //switch on projection of curvatures accordong to paper RezCus2017
            project_curvature = true;
        } else if ( param.compare("convert_from_cauchy") == 0 ) {
            //the element is Cauchy but the material is Cosserat
            convert_from_cauchy = true;
            if ( dim == 2 ) {
                CauchyToCosseratMatrix = Matrix :: Zero(6, 3);
                CauchyToCosseratMatrix(0, 0) = 1.;
                CauchyToCosseratMatrix(1, 1) = 1.;
                CauchyToCosseratMatrix(2, 2) = 0.5;
                CauchyToCosseratMatrix(3, 2) = 0.5;
            } else if ( dim == 3 ) {
                CauchyToCosseratMatrix = Matrix :: Zero(18, 6);
                CauchyToCosseratMatrix(0, 0) = 1.;
                CauchyToCosseratMatrix(1, 1) = 1.;
                CauchyToCosseratMatrix(2, 2) = 1.;
                CauchyToCosseratMatrix(3, 3) = 0.5;
                CauchyToCosseratMatrix(4, 3) = 0.5;
                CauchyToCosseratMatrix(5, 4) = 0.5;
                CauchyToCosseratMatrix(6, 4) = 0.5;
                CauchyToCosseratMatrix(7, 5) = 0.5;
                CauchyToCosseratMatrix(8, 5) = 0.5;
            }
        }
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE COUPLED RVE MATERIAL STATUS
//////////////////////////////////////////////////////////
DiscreteCoupledRVEMaterialStatus ::  DiscreteCoupledRVEMaterialStatus(RVEMaterial *m, Element *e, unsigned ipnum, fs :: path masterfileM, fs :: path masterfileT, unsigned ndim) : RVEMaterialStatus(m, e, ipnum, "none", ndim) {
    name = "coupled discrete RVE mat. status";
    DiscreteCoupledRVEMaterial *coupledm = dynamic_cast< DiscreteCoupledRVEMaterial * >( m );
    if ( coupledm ) {
        mechRVEstat = new DiscreteMechanicalRVEMaterialStatus(coupledm->giveMechanicalRVEmat(), e, idx, masterfileM, ndim);
        trspRVEstat = new DiscreteTransportRVEMaterialStatus(coupledm->giveTransportRVEmat(), e, idx, masterfileT, ndim);

        temp_volumetricStrain = volumetricStrain = 0;
        temp_pressure = pressure = 0;
        temp_crackVolume = crackVolume  = 0;
        volStrainRate = pressureRate = crackVolumeRate = 0;

        is_precomputed = true;
        is_master_status = false;
    } else {
        cerr << "DiscreteCoupledRVEMaterialStatus accepts only DiscreteCoupledRVEMaterial" << endl;
        exit(1);
    }
}

//////////////////////////////////////////////////////////
DiscreteCoupledRVEMaterialStatus ::  ~DiscreteCoupledRVEMaterialStatus() {
    if ( mechRVEstat != nullptr ) {
        delete mechRVEstat;
    }
    if ( trspRVEstat != nullptr ) {
        delete trspRVEstat;
    }
}

//////////////////////////////////////////////////////////
void DiscreteCoupledRVEMaterialStatus ::  setFromPrecomputedToFullModel() {
    if ( !is_precomputed ) {
        return;
    }
    is_precomputed = false;
    mechRVEstat->setFromPrecomputedToFullModel();
    trspRVEstat->setFromPrecomputedToFullModel();
    if ( !is_master_status ) {
        findFriends();
    }
}

//////////////////////////////////////////////////////////
void DiscreteCoupledRVEMaterialStatus ::  init() {
    DiscreteCoupledRVEMaterial *macromaterial = static_cast< DiscreteCoupledRVEMaterial * >( mat );
    mechRVEstat->init();
    trspRVEstat->init();
    if ( macromaterial->givePrecomputedDampingTensor().size() == 0 ) {
        setFromPrecomputedToFullModel();
        Matrix dam, ine;
        dam = giveDampingTensor();
        ine = giveInertiaTensor();
        macromaterial->setPrecomputedDampingAndInertiaTensors(dam, ine);

        TransportPeriodicBC *tp;
        double PUCVolume = 0.;
        PBlockContainer *pblocks = trspRVEstat->giveWholeRVE()->givePBlockContainer();
        for ( unsigned i = 0; i < pblocks->giveSize(); i++ ) {
            tp = dynamic_cast< TransportPeriodicBC * >( pblocks->givePBlock(i) );
            if ( tp ) {
                PUCVolume = tp->giveVolume();
                break;
            }
        }
        if ( PUCVolume == 0 ) {
            cerr << "Error in DiscreteCoupledRVEMaterialStatus: PUC volume was not determined" << endl;
            exit(1);
        } else {
            macromaterial->setPUCVolume(PUCVolume);
        }

        DiscreteTransportRVEMaterial *dtRVEmat = static_cast< DiscreteTransportRVEMaterial * >( trspRVEstat->giveMaterial() );
        VectTrsprtCoupledMaterial *dctm = dynamic_cast< VectTrsprtCoupledMaterial * >( dtRVEmat->giveMasterMaterial() );

        if ( !dctm ) {
            cerr << "Error in " << name << ": transport material in RVE is not VectTrsprtCoupledMaterial" << endl;
            exit(1);
        }
        macromaterial->setMasterMaterial(dctm);
        is_master_status = true;
        mechRVEstat->setToMasterStatus();
        trspRVEstat->setToMasterStatus();
        if ( macromaterial->shouldStartFromPrecomputed() && macromaterial->givePrecomputedElasticTensor().size() == 0 ) {
            //set to precomputed AFTER calculating inertia and damping
            setToPrecomputed();
            Matrix ela = giveStiffnessTensor("secant");
            macromaterial->setPrecomputedElasticTensor(ela);
        }
    }
    if ( !macromaterial->shouldStartFromPrecomputed() ) {
        setFromPrecomputedToFullModel();
    }
}

//////////////////////////////////////////////////////////
void DiscreteCoupledRVEMaterialStatus ::  update() {
    mechRVEstat->update();
    trspRVEstat->update();

    volumetricStrain = temp_volumetricStrain;
    pressure = temp_pressure;
    crackVolume = temp_crackVolume;
}

//////////////////////////////////////////////////////////
void DiscreteCoupledRVEMaterialStatus :: setParameterValue(string code, double value) {
    if ( code.compare("volumetricStrain") == 0 ) {
        temp_volumetricStrain = value;
    } else if ( code.compare("pressure") == 0 ) {
        trspRVEstat->setParameterValue(code, value);
        temp_pressure = value;
    } else {
        MaterialStatus :: setParameterValue(code, value);
    }
}

//////////////////////////////////////////////////////////
void DiscreteCoupledRVEMaterialStatus ::  updateRateVariables(double timeStep) {
    if ( timeStep > 0 ) {
        volStrainRate = ( temp_volumetricStrain - volumetricStrain ) / ( timeStep );
        pressureRate =  ( temp_pressure - pressure ) / ( timeStep );
        crackVolumeRate =  ( temp_crackVolume - crackVolume ) / ( timeStep );
    } else {
        volStrainRate = 0.;
        pressureRate = 0.;
        crackVolumeRate = 0.;
    }
}

//////////////////////////////////////////////////////////
Vector DiscreteCoupledRVEMaterialStatus ::  giveStress(const Vector &strain, double timeStep) {
    DiscreteCoupledRVEMaterial *dcm = static_cast< DiscreteCoupledRVEMaterial * >( mat );
    unsigned ndim = dcm->giveDimension();
    unsigned sizeM = mechRVEstat->giveMaterial()->giveStrainSize();
    unsigned sizeT = trspRVEstat->giveMaterial()->giveStrainSize();
    temp_strain.resize(sizeM + sizeT);
    Vector strainM = Vector :: Zero(sizeM);
    Vector strainT = Vector :: Zero(sizeT + 1); //adding macroscopic pressure
    unsigned i;
    for ( i = 0; i < sizeM; i++ ) {
        strainM [ i ] = temp_strain [ i ] = strain [ i ];
    }
    for ( i = 0; i < sizeT; i++ ) {
        strainT [ i ] = temp_strain [ i + sizeM ] = strain [ i + sizeM ];
    }
    strainT [ sizeT ] = strain [ sizeT + sizeM ]; //macroscopic pressure

    //solve mechanics
    Vector stressM = mechRVEstat->giveStress(strainM, timeStep);
    if ( !mechRVEstat->isPrecomputed() && trspRVEstat->isPrecomputed() ) {
        setFromPrecomputedToFullModel();                                                                 //when mechanical RVE changes to full version
    }
    double mechBiot = -temp_pressure *dcm->giveBiotCoefficient();  //take pressure strored in element for integration point IDX

    //solve transport
    temp_crackVolume = mechRVEstat->giveCrackVolume();
    updateRateVariables(timeStep);
    Vector stressT = trspRVEstat->giveStress(strainT, timeStep);


    temp_stress.resize(sizeM + sizeT);
    for ( i = 0; i < sizeM; i++ ) {
        temp_stress [ i ] = stressM [ i ];
    }

    for ( i = 0; i < ndim; i++ ) {
        temp_stress [ i ] += mechBiot;
    }

    for ( i = 0; i < sizeT; i++ ) {
        temp_stress [ i + sizeM ] = stressT [ i ];
    }
    return temp_stress;
}

//////////////////////////////////////////////////////////
Vector DiscreteCoupledRVEMaterialStatus :: giveInternalSource() const {
    DiscreteCoupledRVEMaterial *dcm = static_cast< DiscreteCoupledRVEMaterial * >( mat );
    unsigned ndim = dcm->giveDimension();


    DiscreteCoupledRVEMaterial *dcRVEmat = static_cast< DiscreteCoupledRVEMaterial * >( mat );
    VectTrsprtCoupledMaterial *dtcm = dcRVEmat->giveMasterMaterial();

    double PUCVolume = dcm->givePUCVolume();
    Vector intS = Vector :: Zero(3 * ( ndim - 1 ) + 1);
    unsigned TDoF = 3 * ( ndim - 1 );

    intS [ TDoF ] = -dcm->giveBiotCoefficient() * volStrainRate * 3.;
    intS [ TDoF ] -= temp_crackVolume * pressureRate / ( PUCVolume * dtcm->giveKw() );
    intS [ TDoF ] -= crackVolumeRate / PUCVolume * ( 1. - dcm->giveBiotCoefficient() + ( temp_pressure - dtcm->giveReferencePressure() ) / dtcm->giveKw() );
    intS [ TDoF ] *= dtcm->giveDensity();
    return intS;
}

//////////////////////////////////////////////////////////
Vector DiscreteCoupledRVEMaterialStatus ::  giveStressWithFrozenIntVars(const Vector &strain, double timeStep) {
    ( void ) timeStep;
    temp_strain = strain;
    DiscreteCoupledRVEMaterial *dcm = static_cast< DiscreteCoupledRVEMaterial * >( mat );
    unsigned ndim = dcm->giveDimension();
    temp_stress = giveStiffnessTensor("secant") * strain;
    updateRateVariables(timeStep);
    double mechBiot = -temp_pressure *dcm->giveBiotCoefficient();  //take pressure strored in element for integration point IDX
    for ( unsigned i = 0; i < ndim; i++ ) {
        temp_stress [ i ] += mechBiot;
    }
    return temp_stress;
}

//////////////////////////////////////////////////////////
bool DiscreteCoupledRVEMaterialStatus ::  giveValues(string code, Vector &result) const {
    if ( code.compare("crack_volume_density") == 0 ) {
        DiscreteCoupledRVEMaterial *dcm = static_cast< DiscreteCoupledRVEMaterial * >( mat );
        result.resize(1);
        result [ 0 ] = temp_crackVolume / dcm->givePUCVolume();
        return true;
    } else if ( code.compare("permeability_tensor") == 0 ) {
        return trspRVEstat->giveValues(code, result);
    } else if ( code.compare("stress") == 0 ) {
        return mechRVEstat->giveValues(code, result);
    } else if ( code.compare("flux") == 0 ) {
        return trspRVEstat->giveValues(code, result);
    } else {
        return RVEMaterialStatus :: giveValues(code, result);
    }
}

//////////////////////////////////////////////////////////
Matrix DiscreteCoupledRVEMaterialStatus :: giveStiffnessTensor(string type) const {
    Matrix M = mechRVEstat->giveStiffnessTensor(type);
    Matrix T = trspRVEstat->giveStiffnessTensor(type);
    unsigned sizeM = M.rows();
    unsigned sizeT = T.rows();

    Matrix D = Matrix :: Zero(sizeM + sizeT, sizeM + sizeT);
    for ( unsigned i = 0; i < sizeM; i++ ) {
        for ( unsigned j = 0; j < sizeM; j++ ) {
            D(i, j) = M(i,  j);
        }
    }
    for ( unsigned i = 0; i < sizeT; i++ ) {
        for ( unsigned j = 0; j < sizeT; j++ ) {
            D(i + sizeM, j + sizeM) = T(i, j);
        }
    }
    return D;
}

//////////////////////////////////////////////////////////
Matrix DiscreteCoupledRVEMaterialStatus :: giveDampingTensor() const {
    if ( is_precomputed ) {
        DiscreteCoupledRVEMaterial *macromat = static_cast< DiscreteCoupledRVEMaterial * >( mat );
        return macromat->givePrecomputedDampingTensor();
    } else {
        Matrix M = mechRVEstat->giveDampingTensor();
        Matrix T = trspRVEstat->giveDampingTensor();
        unsigned sizeM = M.rows();
        unsigned sizeT = T.rows();
        Matrix D = Matrix :: Zero(sizeM + sizeT, sizeM + sizeT);
        for ( unsigned i = 0; i < sizeM; i++ ) {
            for ( unsigned j = 0; j < sizeM; j++ ) {
                D(i, j) = M(i,  j);
            }
        }
        for ( unsigned i = 0; i < sizeT; i++ ) {
            for ( unsigned j = 0; j < sizeT; j++ ) {
                D(i + sizeM, j + sizeM) = T(i, j);
            }
        }
        return D;
    }
}

//////////////////////////////////////////////////////////
Matrix DiscreteCoupledRVEMaterialStatus :: giveInertiaTensor() const {
    if ( is_precomputed ) {
        DiscreteCoupledRVEMaterial *macromat = static_cast< DiscreteCoupledRVEMaterial * >( mat );

        return macromat->givePrecomputedInertiaTensor();
    } else {
        Matrix M = mechRVEstat->giveInertiaTensor();
        Matrix T = trspRVEstat->giveDampingTensor(); //to get correct size
        unsigned sizeM = M.rows();
        unsigned sizeT = T.rows();
        Matrix D = Matrix :: Zero(sizeM + sizeT, sizeM + sizeT);
        for ( unsigned i = 0; i < sizeM; i++ ) {
            for ( unsigned j = 0; j < sizeM; j++ ) {
                D(i, j) = M(i,  j);
            }
        }
        return D;
    }
}
//////////////////////////////////////////////////////////
void DiscreteCoupledRVEMaterialStatus :: findFriends() {
    ElementContainer *elemsM = mechRVEstat->giveWholeRVE()->giveElements();
    ElementContainer *elemsT = trspRVEstat->giveWholeRVE()->giveElements();
    NodeContainer *nodesM = mechRVEstat->giveWholeRVE()->giveNodes();
    NodeContainer *nodesT = trspRVEstat->giveWholeRVE()->giveNodes();

    vector< double >RVEsize;
    PBlockContainer *pblocksM =  mechRVEstat->giveWholeRVE()->givePBlockContainer();

    for ( unsigned i = 0; i < pblocksM->giveSize(); i++ ) {
        MechanicalPeriodicBC *mpb = dynamic_cast< MechanicalPeriodicBC * >( pblocksM->givePBlock(i) );
        if ( mpb ) {
            RVEsize = mpb->giveDimensions();
        }
    }

    DiscreteCoupledRVEMaterial *macromat = static_cast< DiscreteCoupledRVEMaterial * >( mat );
    unsigned ndim = macromat->giveDimension();

    //attach mech elems to node numbers
    vector< vector< RigidBodyContact * > >attachedRBC( nodesM->giveSize() );
    RigidBodyContact *rbc;
    Point insideP;
    Node *foundN;
    double coord, dist;
    bool is_inside;
    for ( unsigned k = 0; k < elemsM->giveSize(); k++ ) {
        rbc = static_cast< RigidBodyContact * >( elemsM->giveElement(k) );
        for ( unsigned p = 0; p < 2; p++ ) {
            is_inside = true;
            insideP = rbc->giveNode(p)->givePoint();
            for ( unsigned v = 0; v < ndim; v++ ) {
                coord = insideP(v);
                if ( coord < 0. ) {
                    insideP(v) = ( coord + RVEsize [ v ] );
                    is_inside = false;
                } else if ( coord > RVEsize [ v ] ) {
                    insideP(v) = ( coord - RVEsize [ v ] );
                    is_inside = false;
                }
            }
            if ( is_inside ) {
                attachedRBC [ nodesM->giveNodeNumber( rbc->giveNode(p) ) ].push_back(rbc);
            } else {
                foundN = nodesM->findClosestMechanicalNode(insideP, & dist);
                if ( dist > 1e-10 ) {
                    cerr << "DiscreteCoupledRVEMaterialStatus Error: searching for periodic image of mechanical node failed, distance to closest point is " << dist << endl;
                    exit(1);
                }
                attachedRBC [ nodesM->giveNodeNumber(foundN) ].push_back(rbc);
            }
        }
    }

    //attach mech elems to transport elements
    DiscreteTrsprtCoupledElem *trsp;
    vector< Node * >vertices;
    vector< unsigned >vnums;
    vector< RigidBodyContact * >mechelems, mechelemsold;
    unsigned pold;
    double weight = 0;
    for ( unsigned k = 0; k < elemsT->giveSize(); k++ ) {
        trsp = dynamic_cast< DiscreteTrsprtCoupledElem * >( elemsT->giveElement(k) );
        if ( !trsp ) {
            cerr << "Coupled transport elements are required" << endl;
            exit(1);
        }
        vertices = trsp->giveVertices();
        vnums.resize( vertices.size() );
        for ( unsigned p = 0; p < vertices.size(); p++ ) {
            is_inside = true;
            insideP = vertices [ p ]->givePoint();
            for ( unsigned v = 0; v < ndim; v++ ) {
                coord = insideP(v);
                if ( coord < 0. ) {
                    insideP(v) = ( coord + RVEsize [ v ] );
                    is_inside = false;
                } else if ( coord > RVEsize [ v ] ) {
                    insideP(v) = ( coord - RVEsize [ v ] );
                    is_inside = false;
                }
            }
            if ( is_inside ) {
                vnums [ p ] = nodesT->giveNodeNumber(vertices [ p ]);
            } else {
                foundN = nodesT->findClosestAuxiliaryNode(insideP, & dist);
                if ( dist > 1e-10 ) {
                    cerr << "DiscreteCoupledRVEMaterialStatus Error: searching for periodic image of transport node failed, distance to closest point is " << dist << endl;
                    exit(1);
                }
                vnums [ p ] = nodesT->giveNodeNumber(foundN);
            }
        }
        mechelemsold = attachedRBC [ vnums.back() ];
        pold = vnums.size() - 1;
        for ( unsigned p = 0; p < vnums.size(); p++ ) {
            mechelems = attachedRBC [ vnums [ p ] ];
            for ( auto me: mechelemsold ) {
                if ( find(mechelems.begin(), mechelems.end(), me) != mechelems.end() ) {
                    if ( ndim == 2 ) {
                        weight = me->giveArea();
                    } else if ( ndim == 3 ) {
                        weight = ( trsp->giveIPLoc(0) - ( vertices [ p ]->givePoint() + vertices [ pold ]->givePoint() ) / 2. ).norm();
                    }
                    trsp->addNewFriend(me, weight);
                    break;
                }
            }
            mechelemsold = mechelems;
            pold = p;
        }
    }
}

//////////////////////////////////////////////////////////
void DiscreteCoupledRVEMaterialStatus ::  setEigenStrain(Vector &x) {
    unsigned sizeM = mechRVEstat->giveMaterial()->giveStrainSize();
    unsigned sizeT = trspRVEstat->giveMaterial()->giveStrainSize();

    Vector eigenstrainM = Vector :: Zero(sizeM);
    Vector eigenstrainT = Vector :: Zero(sizeT);
    unsigned i;
    for ( i = 0; i < sizeM; i++ ) {
        eigenstrainM [ i ] = x [ i ];
    }
    for ( i = 0; i < sizeT; i++ ) {
        eigenstrainT [ i ] = x [ i + sizeM ];
    }
    mechRVEstat->setEigenStrain(eigenstrainM);
    trspRVEstat->setEigenStrain(eigenstrainT);
}

//////////////////////////////////////////////////////////
string DiscreteCoupledRVEMaterialStatus :: giveLineToSave() const {
    return mechRVEstat->giveLineToSave() + trspRVEstat->giveLineToSave();
}

//////////////////////////////////////////////////////////
void DiscreteCoupledRVEMaterialStatus ::  setReferenceSystemDirections(Matrix r) {
    mechRVEstat->setReferenceSystemDirections(r);
    trspRVEstat->setReferenceSystemDirections(r);
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE COUPLED RVE MATERIAL
//////////////////////////////////////////////////////////
DiscreteCoupledRVEMaterial :: DiscreteCoupledRVEMaterial(unsigned dimension) : RVEMaterial(dimension) {
    name = "coupled RVE material";
    precompElastic = Matrix(0, 0);
    produceInternalSources = true;
    start_from_precomputed = true;
    PUCVolume = 0;
    dampMatUpdate = true;
};


//////////////////////////////////////////////////////////
void DiscreteCoupledRVEMaterial :: init(MaterialContainer *matcont) {
    RVEMaterial :: init(matcont);
    if ( !nonlinear ) {
        mechRVEmat->setStartFromPrecomputed(true);
        trspRVEmat->setStartFromPrecomputed(true);
    }
}

//////////////////////////////////////////////////////////
DiscreteCoupledRVEMaterial :: ~DiscreteCoupledRVEMaterial() {
    if ( mechRVEmat != nullptr ) {
        delete mechRVEmat;
    }
    if ( trspRVEmat != nullptr ) {
        delete trspRVEmat;
    }
}

//////////////////////////////////////////////////////////
MaterialStatus *DiscreteCoupledRVEMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    DiscreteCoupledRVEMaterialStatus *newstat = new DiscreteCoupledRVEMaterialStatus(this, e, ipnum, mechRVEmat->givePathToInputFile(), trspRVEmat->givePathToInputFile(), dim);
    return newstat;
};

//////////////////////////////////////////////////////////
void DiscreteCoupledRVEMaterial :: readFromLine(istringstream &iss) {
    mechRVEmat = new DiscreteMechanicalRVEMaterial(dim);
    trspRVEmat = new DiscreteTransportRVEMaterial(dim);
    strainsize = mechRVEmat->giveStrainSize() + trspRVEmat->giveStrainSize();

    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream

    string param;
    bool bdim, bfilem, bfilet, bbiot;
    bdim = bfilem = bfilet = bbiot = false;
    string fileM, fileT;

    while (  iss >> param ) {
        if ( param.compare("dim") == 0 ) {
            bdim = true;
            iss >> dim;
        } else if ( param.compare("RVEfolderMech") == 0 ) {
            bfilem = true;
            iss >> fileM;
            mechRVEmat->setPathToInputFolder(fileM);
        } else if ( param.compare("RVEfolderTrsprt") == 0 ) {
            bfilet = true;
            iss >> fileT;
            trspRVEmat->setPathToInputFolder(fileT);
        } else if ( param.compare("biot_coeff") == 0 ) {
            bbiot = true;
            iss >> biotCoeff;
        } else if ( param.compare("enforce_linearity") == 0 ) {
            nonlinear = false;
            mechRVEmat->enforceLinearity();
            trspRVEmat->enforceLinearity();
        } else if ( param.compare("do_not_precompute") == 0 ) {
            start_from_precomputed = false;
            mechRVEmat->setStartFromPrecomputed(false);
            trspRVEmat->setStartFromPrecomputed(false);
        }
    }
    if ( !bdim ) {
        cerr << name << ": material parameter 'dim' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    if ( !bfilem ) {
        cerr << name << ": material parameter 'RVEfolderMech' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    if ( !bfilet ) {
        cerr << name << ": material parameter 'RVEfolderTrsprt' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    if ( !bbiot ) {
        cerr << name << ": material parameter 'biot_coeff' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
}

//////////////////////////////////////////////////////////
void DiscreteCoupledRVEMaterial :: setPrecomputedElasticTensor(Matrix ela) {
    precompElastic = ela;
}

//////////////////////////////////////////////////////////
void DiscreteCoupledRVEMaterial :: setPrecomputedDampingAndInertiaTensors(Matrix dam, Matrix ine) {
    precompDamping = dam;
    precompInertia = ine;
}
