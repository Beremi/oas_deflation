#include "material_RVE.h"
#include "model.h"
#include "element_discrete.h"
#include "element_continuous.h"

using namespace std;

# define M_PI           3.14159265358979323846  /* pi */

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// general RVE status
RVEMaterialStatus :: RVEMaterialStatus(RVEMaterial *m, Element *e, unsigned ipnum, fs :: path masterfile, unsigned ndim) : MaterialStatus(m, e, ipnum) {
    name = "general RVE mat. status";
    inputfile = masterfile;
    RVE = new Model(false);

    axDirs = MyMatrix :: Zero(ndim, ndim);
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
    RVE->readFromFile(inputfile.string() );
    //here the structre of model initialization should be coppied
    //we needed to insert volumetric average generation after applying preprocessing block, otherwise constraints from preprocessing block would not be active
    //therefore, preprocessing blocks are now called in reader instead of model initialization

    //generateVolumetricAverageBC(); //not good idea because the stiffness matrix becomes full
    generateRandomFixedBC();
    RVE->init();

    stringstream appendname;
    appendname << "_" << std :: setfill('0') << std :: setw(4) << element->giveID() << "_" << std :: setw(2) << idx;
    RVE->giveExporters()->appendToAllNames(appendname.str() );

    RVEMaterial *macromat = static_cast< RVEMaterial * >( mat );
    unsigned ndim = RVE->giveDimension();

    macromat->setNumOfDimensions(ndim);
}

//////////////////////////////////////////////////////////
void RVEMaterialStatus :: update() {
    Solver *solver = RVE->giveSolver();
    Solver *masterSolver = masterModel->giveSolver();
    solver->runAfterEachStep();   //update material statuses
    RVE->giveExporters()->exportData(masterSolver->giveStepNumber(), masterSolver->giveTime(), solver->giveDoFValues(), solver->giveNodalForces(), masterSolver->isTerminated() );
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// GENERAL RVE MATERIAL
MaterialStatus *RVEMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    RVEMaterialStatus *newstat = new RVEMaterialStatus(this, e, ipnum, inputfile, ndim);
    return newstat;
}

//////////////////////////////////////////////////////////
void RVEMaterial :: readFromLine(istringstream &iss) {
    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream

    string param;
    bool bdim, bfile;
    bdim = bfile = false;

    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("dim") == 0 ) {
            bdim = true;
            iss >> ndim;
        } else if ( param.compare("RVEfolder") == 0 ) {
            bfile = true;
            iss >> inputfile;
        } else if ( param.compare("RVEfolder") == 0 ) {
            bfile = true;
            iss >> inputfile;
        } else if ( param.compare("enforce_linearity") == 0 ) {
            nonlinear = false;
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
    unsigned ndim = macromat->giveNumOfDimensions();
    ElementContainer *elems = RVE->giveElements();
    Point normal;
    MyVector eigstr = MyVector :: Zero(1);
    Transp1D *e;
    TrsprtMaterial *m;
    TrsprtMaterialStatus *s;
    for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
        eigstr *= 0.;
        e = static_cast< Transp1D * >( elems->giveElement(i) );
        s = static_cast< TrsprtMaterialStatus * >( e->giveMatStatus(0) );
        m = static_cast< TrsprtMaterial * >( s->giveMaterial() );
        m->setPermeability(orig_mater_params [ 2 * i ]); //set back original permeability
        m->setPermeability(s->calculatePressureDependentPermeability(macro_pressure) );   //calculating pressure depedent conductivity
        m->setParamA(-1.); //switch of linearity
        normal = e->giveNormal();
        for ( unsigned v = 0; v < ndim; v++ ) {
            eigstr [ 0 ] -= local_strain [ v ] * normal(v);
        }
        e->giveMatStatus(0)->setEigenStrain(eigstr);
    }
}



/////////////////////////////////./////////////////////////
void DiscreteTransportRVEMaterialStatus :: collectStresses() {
    DiscreteTransportRVEMaterial *macromat = static_cast< DiscreteTransportRVEMaterial * >( mat );
    unsigned ndim = macromat->giveNumOfDimensions();
    ElementContainer *elems = RVE->giveElements();
    double volume = 0.;
    MyVector factor;
    Transp1D *e;
    TrsprtMaterial *m;
    Point normal;
    local_stress *= 0.;
    for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
        e = static_cast< Transp1D * >( elems->giveElement(i) );
        m = static_cast< TrsprtMaterial * >( e->giveMatStatus(0)->giveMaterial() );
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
MyVector DiscreteTransportRVEMaterialStatus :: giveStressPrecomputed(const MyVector &strain, double timeStep) {
    ( void ) timeStep;

    temp_strain = strain;

    transformStrain();    //stiffness is already transformed

    DiscreteTransportRVEMaterial *macromaterial = static_cast< DiscreteTransportRVEMaterial * >( mat );
    temp_nonlin = macromaterial->giveMasterStatus()->calculatePressureDependentPermeability(macro_pressure) / macromaterial->giveMasterMaterial()->givePermeability();
    local_stress = (macromaterial->givePrecomputedConductivity() * local_strain) * temp_nonlin;

    transformStress();

    return temp_stress;
}

/////////////////////////////////./////////////////////////
MyVector DiscreteTransportRVEMaterialStatus :: giveStress(const MyVector &strain, double timeStep) {
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
    local_stress.resize( temp_strain.size() );
    local_stress *= 0;
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
MyVector DiscreteTransportRVEMaterialStatus :: giveStressWithFrozenIntVars(const MyVector &strain, double timeStep) {
    ( void ) timeStep;

    temp_strain = strain;
    DiscreteTransportRVEMaterial *macromat = static_cast< DiscreteTransportRVEMaterial * >( mat );
    unsigned ndim = macromat->giveNumOfDimensions();
    transformStrain();
    local_stress = giveStiffnessTensorLocal("elastic", ndim) * local_strain;
    transformStress();

    return temp_stress;
}

/////////////////////////////////./////////////////////////
unsigned DiscreteTransportRVEMaterialStatus :: giveStrainSize(unsigned dimension) const {
    return dimension;
}


/////////////////////////////////./////////////////////////
unsigned DiscreteTransportRVEMaterialStatus :: giveStrainSize() const {
    RVEMaterial *macromat = static_cast< RVEMaterial * >( mat );
    return giveStrainSize( macromat->giveNumOfDimensions() );
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

    //mechanics
    JointDoF *jd;
    Node *masternode;
    bool found = false;
    for ( unsigned j = 0; j < constrs->giveSize(); j++ ) {
        if ( found ) {
            break;
        }
        jd = constrs->giveConstraint(j);
        for ( unsigned k = 0; k < jd->giveNumOfDoFMasters(); k++ ) {
            masternode = jd->giveMasterNode(0);
            if ( dynamic_cast< MechDoF * >( masternode ) == nullptr && dynamic_cast< TrsDoF * >( masternode ) == nullptr ) {
                BoundaryCondition *bc;
                vector< int >dBC, nBC;
                dBC.resize(masternode->giveNumberOfDoFs(), funcs->giveSize() );   //todo: warning C4267: 'argument': conversion from 'size_t' to 'const _Ty', possible loss of data
                nBC.resize(masternode->giveNumberOfDoFs(), -1);
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
    unsigned ndim = macromat->giveNumOfDimensions();

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
            vm.push_back(nodes->giveNode(n) );
        }
    }
    if ( vm.size() > 0 ) {
        unsigned nDoFs = 3;
        if ( ndim == 3 ) {
            nDoFs = 6;
        }
        MechDoF *pn = new MechDoF(ndim, nDoFs);   //?? for transport
        nodes->addNode(pn);

        dirs.resize(vm.size() );

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
MyMatrix DiscreteTransportRVEMaterialStatus :: giveStiffnessTensorLocal(string type, unsigned dimension) const {
    ( void ) type;

    if ( is_precomputed ) {
        return giveStiffnessTensorPrecomputedLocal(type, dimension);
    }

    unsigned strain_size = giveStrainSize(dimension);
    MyMatrix Keff = MyMatrix :: Zero(strain_size, strain_size);

    ElementContainer *elems = RVE->giveElements();

    double volume = 0;
    Point normal;
    MyVector n = MyVector :: Zero(dimension);
    Transp1D *e;
    volume = 0;
    TrsprtMaterialStatus *tstat;

    for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
        e = static_cast< Transp1D * >( elems->giveElement(i) );
        normal = e->giveNormal();
        for ( unsigned v = 0; v < dimension; v++ ) {
            n [ v ] = normal(v);
        }
        tstat = static_cast< TrsprtMaterialStatus * >( e->giveMatStatus(0) );
        Keff -= dyadicProduct(n, n) * ( e->giveLength() * e->giveArea() * tstat->giveEffectiveConductivity("tangent") );
        volume += e->giveVolume();
    }
    Keff /= volume;

    return Keff;
}

//////////////////////////////////////////////////////////
MyMatrix DiscreteTransportRVEMaterialStatus :: giveStiffnessTensor(string type, unsigned dimension) const {
    MyMatrix Keff = giveStiffnessTensorLocal(type, dimension);
    return ( transf.transpose() * Keff ) * transf;
}


//////////////////////////////////////////////////////////
MyMatrix DiscreteTransportRVEMaterialStatus :: giveStiffnessTensorPrecomputedLocal(string type, unsigned rdim) const {
    ( void ) rdim;
    ( void ) type;
    DiscreteTransportRVEMaterial *macromaterial = static_cast< DiscreteTransportRVEMaterial * >( mat );
    MyMatrix Keff = macromaterial->givePrecomputedConductivity() * temp_nonlin;
    return Keff;
}


//////////////////////////////////////////////////////////
MyMatrix DiscreteTransportRVEMaterialStatus :: giveDampingTensor() const {
    if ( is_precomputed ) {
        return giveDampingTensorPrecomputed();
    }

    ElementContainer *elems = RVE->giveElements();
    Transp1D *e;
    double volume = 0;
    double mass = 0;
    for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
        e = dynamic_cast< Transp1D * >( elems->giveElement(i) );
        if ( e ) {
            mass += e->giveVolume() * e->giveMatStatus(0)->giveDampingTensor()(0, 0);
            volume += e->giveVolume();
        }
    }
    MyMatrix m = MyMatrix :: Zero(1, 1);
    m(0, 0) = mass / volume;
    return m;
};

//////////////////////////////////////////////////////////
MyMatrix DiscreteTransportRVEMaterialStatus :: giveDampingTensorPrecomputed() const {
    DiscreteTransportRVEMaterial *macromaterial = static_cast< DiscreteTransportRVEMaterial * >( mat );
    MyMatrix M = MyMatrix :: Zero(1, 1);
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
    Transp1D *e;
    TrsprtMaterial *m;
    for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
        e = static_cast< Transp1D * >( elems->giveElement(i) );
        m = static_cast< TrsprtMaterial * >( e->giveMatStatus(0)->giveMaterial() );
        orig_mater_params [ 2 * i  ] = m->givePermeability();
        orig_mater_params [ 2 * i + 1 ] = m->giveParamA();
    }
}

//////////////////////////////////////////////////////////
void DiscreteTransportRVEMaterialStatus :: init() {
    DiscreteTransportRVEMaterial *macromaterial = static_cast< DiscreteTransportRVEMaterial * >( mat );
    MyMatrix stiff = macromaterial->givePrecomputedConductivity();
    if ( stiff.size() == 0 ) {
        DiscreteTransportRVEMaterialStatus :: setFromPrecomputedToFullModel();
        is_master_status = true;
        unsigned ndim = macromaterial->giveNumOfDimensions();
        calculateTransformationMatrix();

        MyMatrix c = DiscreteTransportRVEMaterialStatus :: giveDampingTensor();


        MyMatrix Keff = MyMatrix :: Zero(ndim, ndim);
        if ( macromaterial->isElasticSolutionVoigt() ) {
            //use Voigt constraint
            Keff = giveStiffnessTensorLocal("elastic", ndim);
        } else  {
            //compute exactly
            cout << "Precomputing primary fields for transport RVE" << endl;
            MyVector help_strain = MyVector :: Zero(ndim);
            MyVector help_stress;
            double factor = 1e-10;
            for ( unsigned i = 0; i < ndim; i++ ) {
                help_strain [ i ] = factor;
                help_stress = giveStress(help_strain, -1);
                help_strain [ i ] = 0.;
                for ( unsigned j = 0; j < ndim; j++ ) {
                    Keff(i, j) = help_stress [ j ] / factor;
                }
            }
        }

        TrsprtMaterialStatus *status = static_cast< TrsprtMaterialStatus * >( RVE->giveElements()->giveElement(0)->giveMatStatus(0) );
        TrsprtMaterial *material = static_cast< TrsprtMaterial * >( RVE->giveElements()->giveElement(0)->giveMaterial() );
        macromaterial->setPrecomputedConductivityAndCapacityAndMasterMaterial(Keff, c(0, 0), status, material);

        is_precomputed = true; //set back to precomputed to use it
    } else {
        calculateTransformationMatrix();
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
MaterialStatus *DiscreteTransportRVEMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    DiscreteTransportRVEMaterialStatus *newstat = new DiscreteTransportRVEMaterialStatus(this, e, ipnum, inputfile, ndim);
    return newstat;
}


//////////////////////////////////////////////////////////
void DiscreteTransportRVEMaterial :: setPrecomputedConductivityAndCapacityAndMasterMaterial(MyMatrix lam, double c, TrsprtMaterialStatus *masterS, TrsprtMaterial *masterM) {
    conductivity = lam;
    capacity = c;
    masterStatus = masterS;
    masterMaterial = masterM;

    const char *buffer = "precomputed_conductivity.out";
    ofstream outputfile( ( masterModel->giveResultDirectory() / buffer ).string() );
    if ( outputfile.is_open() ) {
        outputfile << std :: scientific;
        outputfile.precision(10);

        unsigned size = conductivity.rows();
        for ( unsigned i = 0; i < size; i++ ) {
            for ( unsigned j = 0; j < size; j++ ) {
                if ( j > 0 ) {
                    outputfile << "\t";
                }
                outputfile << conductivity(i, j);
            }
            outputfile << endl;
        }
        outputfile.close();
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE MECHANICAL RVE MATERIAL
//////////////////////////////////////////////////////////
DiscreteMechanicalRVEMaterialStatus :: DiscreteMechanicalRVEMaterialStatus(RVEMaterial *m, Element *e, unsigned ipnum, fs :: path masterfile, unsigned ndim) : DiscreteTransportRVEMaterialStatus(m, e, ipnum, masterfile, ndim) {
    name = "mechancial RVE mat. status";
    is_precomputed = true;
    is_master_status = false;
}

/////////////////////////////////./////////////////////////
void DiscreteMechanicalRVEMaterialStatus :: applyEigenStrains() {
    //set eigenstrains
    DiscreteMechanicalRVEMaterial *macromat = static_cast< DiscreteMechanicalRVEMaterial * >( mat );
    vector< vector< MyVector > > *projectors = macromat->giveProjectors();
    ElementContainer *elems = RVE->giveElements();
    unsigned ndim = macromat->giveNumOfDimensions();
    MyVector eigstr = MyVector :: Zero(ndim);
    unsigned num = 0;
    for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
        eigstr *= 0.;
        for ( unsigned v = 0; v < ndim; v++ ) {
            for ( unsigned r = 0; r < local_strain.size(); r++ ) {
                eigstr [ v ] -= local_strain [ r ] * ( * projectors ) [ v ] [ num ] [ r ];
            }
        }
        num++;
        elems->giveElement(i)->giveMatStatus(0)->setEigenStrain(eigstr);
    }
}

/////////////////////////////////./////////////////////////
void DiscreteMechanicalRVEMaterialStatus :: collectStresses() {
    DiscreteMechanicalRVEMaterial *macromat = static_cast< DiscreteMechanicalRVEMaterial * >( mat );
    vector< vector< MyVector > > *projectors = macromat->giveProjectors();
    ElementContainer *elems = RVE->giveElements();
    unsigned ndim = macromat->giveNumOfDimensions();
    double volume = 0.;
    MyVector factor;
    unsigned num = 0;
    RigidBodyContact *e;
    for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
        e = static_cast< RigidBodyContact * >( elems->giveElement(i) );
        factor  = e->giveArea() * e->giveLength() * e->giveMatStatus(0)->giveTempStress();
        for ( unsigned v = 0; v < ndim; v++ ) {
            for ( unsigned r = 0; r < local_strain.size(); r++ ) {
                local_stress [ r ] += factor [ v ] * ( * projectors ) [ v ] [ num ] [ r ];
            }
        }
        num++;
        volume += e->giveVolume();
    }
    for ( unsigned v = 0; v < local_strain.size(); v++ ) {
        local_stress [ v ] /= volume;
    }
}


/////////////////////////////////./////////////////////////
MyVector DiscreteMechanicalRVEMaterialStatus :: giveStressPrecomputed(const MyVector &strain, double timeStep) {
    ( void ) timeStep;
    temp_strain = addEigenStrain(strain); //macroscopic eigenstrain
    transformStrain();
    DiscreteMechanicalRVEMaterial *macromaterial = static_cast< DiscreteMechanicalRVEMaterial * >( mat );
    local_stress = macromaterial->givePrecomputedElasticTensor() * local_strain;
    transformStress();
    if ( macromaterial->isNonlinear() && checkOttosenCriterion() ) {
        setFromPrecomputedToFullModel();
        return giveStress(strain, timeStep);
    }
    return temp_stress;
}

/////////////////////////////////./////////////////////////
MyVector DiscreteMechanicalRVEMaterialStatus :: giveStress(const MyVector &strain, double timeStep) {
    ( void ) timeStep;

    //cout << "strain " << endl;
    //for (auto s: strain) cout << s << endl;

    //precomputed material
    if ( is_precomputed ) {
        return giveStressPrecomputed(strain, timeStep);
    }

    cout << "Solving mechanical RVE" << endl;

    temp_strain = addEigenStrain(strain); //macroscopic eigenstrain

    transformStrain();

    //set eigenstrains
    applyEigenStrains();

    //solve
    RVE->resetTime();
    RVE->giveSolver()->runBeforeEachStep();
    RVE->giveSolver()->solve();


    //collect results
    local_stress.resize(temp_strain.size() );
    local_stress *= 0;
    collectStresses();

    transformStress();

    /*
     * cout << "STRAIN MECH";
     * cout << endl;
     * for ( unsigned v = 0; v < 6; v++ ) {
     * cout << " " << local_strain [ v ];
     * }
     * cout << endl;
     * for ( unsigned v = 0; v < 6; v++ ) {
     * cout << " " << temp_strain [ v ];
     * }
     * cout << endl;
     *
     * cout << "STRESS MECH";
     * cout << endl;
     * for ( unsigned v = 0; v < 6; v++ ) {
     * cout << " " << local_stress [ v ];
     * }
     * cout << endl;
     * for ( unsigned v = 0; v < 6; v++ ) {
     * cout << " " << temp_stress [ v ];
     * }
     * cout << endl;
     */
    /*
     * cout << "STRAIN MECH";
     * for ( unsigned v = 0; v < temp_strain.size(); v++ ) {
     *  cout << " " << temp_strain [ v ];
     * }
     * cout << endl;
     *
     * cout << "STRESS MECH";
     * for ( unsigned v = 0; v < temp_strain.size(); v++ ) {
     *  cout << " " << temp_stress [ v ];
     * }
     * cout << endl;
     */

    return temp_stress;
}

/////////////////////////////////./////////////////////////
bool DiscreteMechanicalRVEMaterialStatus :: checkOttosenCriterion() {
    vector< MyVector >eigvecs;
    MyVector eignums;

    //from paper Adaptive multiscale homogenization of the lattice discrete particle model for the analysis of damage and fracture in concrete, Roozbeh Rezakhani and Xinwei Zhou and Gianluca Cusatis
    double fc =  20e6 / 2.;
    double ft = 2.35e6 / 2.;
    double c1 = 4.775 / fc;
    double c2 = 7.048 / fc;
    double c3 = 0.88 / pow(fc, 2);

    MyVector stress_without_couple = MyVector :: Zero(9);
    for(unsigned p=0; p<9; p++) stress_without_couple[p] = temp_stress[p];
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
MyVector DiscreteMechanicalRVEMaterialStatus :: giveStressWithFrozenIntVars(const MyVector &strain, double timeStep) {
    ( void ) timeStep;
    temp_strain = strain;
    transformStrain();
    DiscreteMechanicalRVEMaterial *macromat = static_cast< DiscreteMechanicalRVEMaterial * >( mat );
    unsigned ndim = macromat->giveNumOfDimensions();
    local_stress = giveStiffnessTensorLocal("elastic", ndim) * local_strain;
    transformStress();
    return temp_stress;
}

/////////////////////////////////./////////////////////////
unsigned DiscreteMechanicalRVEMaterialStatus :: giveStrainSize(unsigned rdim) const {
    unsigned strain_size = rdim * rdim;
    strain_size += ( rdim == 3 ) ? rdim * rdim : rdim; //in 2D only vector
    return strain_size;
}

/////////////////////////////////./////////////////////////
double DiscreteMechanicalRVEMaterialStatus :: giveCrackVolume() const {
    if ( is_precomputed ) {
        return 0;
    }
    double crackVolume = 0.;
    RigidBodyContact *e;
    ElementContainer *elems = RVE->giveElements();
    MyVector res;
    for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
        e = static_cast< RigidBodyContact * >( elems->giveElement(i) );
        e->giveIPValues("tempCrackOpening",0, res);
        if (res.size()>0) crackVolume += res[0] * e->giveArea();
    }

    return crackVolume;
}

/////////////////////////////////./////////////////////////
unsigned DiscreteMechanicalRVEMaterialStatus :: giveStrainSize() const {
    RVEMaterial *macromat = static_cast< RVEMaterial * >( mat );
    return giveStrainSize( macromat->giveNumOfDimensions() );
}

//////////////////////////////////////////////////////////
MyMatrix DiscreteMechanicalRVEMaterialStatus :: giveStiffnessTensorLocal(string type, unsigned rdim) const {
    DiscreteMechanicalRVEMaterial *macromat = static_cast< DiscreteMechanicalRVEMaterial * >( mat );
    if ( is_precomputed ) {
        return macromat->givePrecomputedElasticTensor();
    } else {
        unsigned strain_size = giveStrainSize(rdim);
        MyMatrix Keff = MyMatrix :: Zero(strain_size, strain_size);

        vector< vector< MyVector > > *projectors = macromat->giveProjectors();
        ElementContainer *elems = RVE->giveElements();
        double volume = 0;
        Point normal;
        MyVector n = MyVector :: Zero(rdim);

        RigidBodyContact *e;
        volume = 0;
        MyMatrix stiff;
        unsigned num = 0;
        for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
            e = static_cast< RigidBodyContact * >( elems->giveElement(i) );
            stiff  = e->giveMatStatus(0)->giveStiffnessTensor(type, rdim);
            for ( unsigned v = 0; v < rdim; v++ ) {
                Keff += dyadicProduct( ( * projectors ) [ v ] [ num ], ( * projectors ) [ v ] [ num ]) * e->giveLength() * e->giveArea() * stiff(v, v);
            }
            num++;
            volume += e->giveVolume();
        }
        Keff /= volume;

        return Keff;
    }
}

//////////////////////////////////////////////////////////
MyMatrix DiscreteMechanicalRVEMaterialStatus :: giveStiffnessTensor(string type, unsigned rdim) const {
    MyMatrix Keff = giveStiffnessTensorLocal(type, rdim);
    return ( transf.transpose() * Keff ) * transf;
}

//////////////////////////////////////////////////////////
MyMatrix DiscreteMechanicalRVEMaterialStatus :: giveDampingTensor() const {
    DiscreteMechanicalRVEMaterial *macromat = static_cast< DiscreteMechanicalRVEMaterial * >( mat );
    if ( is_precomputed ) {
        return macromat->givePrecomputedDampingTensor();
    } else {
        unsigned ndim = macromat->giveNumOfDimensions();
        unsigned dof_size = 3 * ( ndim - 1 );
        MyMatrix Keff = MyMatrix :: Zero(dof_size, dof_size);
        for ( unsigned i = 0; i < dof_size; i++ ) {
            Keff(i, i) = 1e-20;
        }
        return Keff;
    }
}

//////////////////////////////////////////////////////////
MyMatrix DiscreteMechanicalRVEMaterialStatus :: giveInertiaTensor() const {
    DiscreteMechanicalRVEMaterial *macromat = static_cast< DiscreteMechanicalRVEMaterial * >( mat );
    if ( is_precomputed ) {
        return macromat->givePrecomputedInertiaTensor();
    } else {
        unsigned ndim = macromat->giveNumOfDimensions();
        unsigned dof_size = 3 * ( ndim - 1 );
        MyMatrix Keff = MyMatrix :: Zero(dof_size, dof_size);
        return Keff;
    }
}
//////////////////////////////////////////////////////////
Point DiscreteMechanicalRVEMaterialStatus :: calculateCentroid() {  //only for mechanics
    Point centroid = Point(0., 0., 0.);

    //find centroid
    ElementContainer *elems = RVE->giveElements();

    //mechanics
    RigidBodyContact *e;
    DisMechMaterialStatus *ms;
    double weight;
    double totalweight = 0;
    for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
        e = static_cast< RigidBodyContact * >( elems->giveElement(i) );
        ms = static_cast< DisMechMaterialStatus * >( e->giveMatStatus(0) );
        weight = e->giveVolume() * ms->giveDensity();
        totalweight += weight;
        centroid += e->giveIPLoc(0) * weight;
    }
    centroid /= totalweight;
    return centroid;
}

//////////////////////////////////////////////////////////
vector< vector< MyVector > >DiscreteMechanicalRVEMaterialStatus :: calculateProjectors(const Point centroid) {   //only for mechanics
    vector< vector< MyVector > >projectors;

    DiscreteTransportRVEMaterial *macromat = static_cast< DiscreteTransportRVEMaterial * >( mat );
    unsigned ndim = macromat->giveNumOfDimensions();
    ElementContainer *elems = RVE->giveElements();
    projectors.resize(ndim);
    unsigned strain_size = giveStrainSize(ndim);

    RigidBodyContact *e;
    Point xc;
    Point alphaVec;
    Point normal;
    MyVector PQ = MyVector :: Zero(strain_size);
    for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
        e = dynamic_cast< RigidBodyContact * >( elems->giveElement(i) );
        if ( e ) {
            xc =  e->giveIPLoc(0) - centroid; //here we make corrections to have origin of the reference system at the centroid
            normal = e->giveNormal();
            for ( unsigned v = 0; v < ndim; v++ ) {
                if ( v == 0 ) {
                    alphaVec = normal;
                } else if ( v == 1 ) {
                    alphaVec = e->giveT1();
                } else if ( v == 2 ) {
                    alphaVec = e->giveT2();
                }

                if ( ndim == 2 ) {
                    PQ [ 0 ] = normal.x() * alphaVec.x();
                    PQ [ 1 ] = normal.y() * alphaVec.y();
                    PQ [ 2 ] = normal.x() * alphaVec.y();
                    PQ [ 3 ] = normal.y() * alphaVec.x();
                    double factor2D = xc.x() * alphaVec.y() - xc.y() * alphaVec.x();
                    PQ [ 4 ] = factor2D * normal.x();
                    PQ [ 5 ] = factor2D * normal.y();
                } else {
                    PQ [ 0 ] = normal.x() * alphaVec.x();
                    PQ [ 1 ] = normal.y() * alphaVec.y();
                    PQ [ 2 ] = normal.z() * alphaVec.z();
                    PQ [ 3 ] = normal.z() * alphaVec.y();
                    PQ [ 4 ] = normal.y() * alphaVec.z();
                    PQ [ 5 ] = normal.z() * alphaVec.x();
                    PQ [ 6 ] = normal.x() * alphaVec.z();
                    PQ [ 7 ] = normal.y() * alphaVec.x();
                    PQ [ 8 ] = normal.x() * alphaVec.y();
                    Point factor3D = xc.cross(alphaVec);
                    PQ [ 9 ]  = normal.x() * factor3D.x();
                    PQ [ 10 ] = normal.y() * factor3D.y();
                    PQ [ 11 ] = normal.z() * factor3D.z();
                    PQ [ 12 ] = normal.z() * factor3D.y();
                    PQ [ 13 ] = normal.y() * factor3D.z();
                    PQ [ 14 ] = normal.z() * factor3D.x();
                    PQ [ 15 ] = normal.x() * factor3D.z();
                    PQ [ 16 ] = normal.y() * factor3D.x();
                    PQ [ 17 ] = normal.x() * factor3D.y();
                }
                projectors [ v ].push_back(PQ);
            }
        }
    }
    return projectors;
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
    MyMatrix D = macromaterial->givePrecomputedElasticTensor();
    if ( D.size() == 0 ) {
        DiscreteMechanicalRVEMaterialStatus :: setFromPrecomputedToFullModel();
        is_master_status = true;
        unsigned ndim = RVE->giveDimension();
        calculateTransformationMatrix();

        Point centroid = calculateCentroid();
        vector< vector< MyVector > >projectors = calculateProjectors(centroid);
        macromaterial->setCentroidAndProjectors(centroid, projectors);


        unsigned strain_size = giveStrainSize(ndim);
        MyMatrix Keff = MyMatrix :: Zero(strain_size, strain_size);
        if ( macromaterial->isElasticSolutionVoigt() ) {
            //use Voigt constraint
            Keff = giveStiffnessTensorLocal("secant", ndim);
        } else  {
            //compute exactly
            cout << "Precomputing primary fields for mechanical RVE" << endl;
            MyVector help_strain = MyVector :: Zero(strain_size);
            MyVector help_stress;
            double factor = 1e-10;
            for ( unsigned i = 0; i < strain_size; i++ ) {
                help_strain [ i ] = factor;
                help_stress = giveStress(help_strain, -1);
                help_strain [ i ] = 0.;
                for ( unsigned j = 0; j < strain_size; j++ ) {
                    Keff(i, j) = help_stress [ j ] / factor;
                }
            }
        }

        macromaterial->setPrecomputedElasticTensor(Keff);
        macromaterial->setPrecomputedDampingTensor(giveDampingTensor() );
        macromaterial->setPrecomputedInertiaTensor(giveInertiaTensor() );

        is_precomputed = true; //set back to precomputed to use it
    } else {
        calculateTransformationMatrix();
    }
}


//////////////////////////////////////////////////////////
void DiscreteMechanicalRVEMaterialStatus :: transformStrain() {
    local_strain = transf * temp_strain;
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
void giveVectorCoords3D(unsigned *pos, unsigned i, unsigned j) {
    if ( i == 0 ) {
        if ( j == 0 ) {
            * pos = 0;
        } else if ( j == 1 ) {
            * pos = 8;
        } else if ( j == 2 )                                 {
            * pos = 6;
        } else                                                              {
            cerr << "giveTensorCoords should never go here" << endl;
            exit(1);
        }
    } else if ( i == 1 )       {
        if ( j == 0 ) {
            * pos = 7;
        } else if ( j == 1 ) {
            * pos = 1;
        } else if ( j == 2 )                                 {
            * pos = 4;
        } else                                                              {
            cerr << "giveTensorCoords should never go here" << endl;
            exit(1);
        }
    } else if ( i == 2 )      {
        if ( j == 0 ) {
            * pos = 5;
        } else if ( j == 1 ) {
            * pos = 3;
        } else if ( j == 2 )                                 {
            * pos = 2;
        } else                                                              {
            cerr << "giveTensorCoords should never go here" << endl;
            exit(1);
        }
    } else  {
        cerr << "giveTensorCoords should never go here" << endl;
        exit(1);
    }
}

//////////////////////////////////////////////////////////
void DiscreteMechanicalRVEMaterialStatus :: calculateTransformationMatrix() {
    DiscreteMechanicalRVEMaterial *macromat = static_cast< DiscreteMechanicalRVEMaterial * >( mat );
    unsigned ndim = macromat->giveNumOfDimensions();

    if ( ndim == 2 ) {
        cerr << "RVE transformation matrix not implemented yet" << endl;
        exit(1);
    } else {
        transf = MyMatrix :: Zero(18, 18);
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
//////////////////////////////////////////////////////////
// DISCRETE MECHANICAL RVE MATERIAL
//////////////////////////////////////////////////////////
MaterialStatus *DiscreteMechanicalRVEMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    DiscreteMechanicalRVEMaterialStatus *newstat = new DiscreteMechanicalRVEMaterialStatus(this, e, ipnum, inputfile, ndim);
    return newstat;
}

//////////////////////////////////////////////////////////
void DiscreteMechanicalRVEMaterial :: setCentroidAndProjectors(Point c, vector< vector< MyVector > >p) {
    centroid = c;
    projectors = p;
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
    if ( macromaterial->givePrecomputedElasticTensor().size() == 0 ) {
        unsigned ndim = trspRVEstat->giveWholeRVE()->giveDimension();
        macromaterial->setNumOfDimensions(ndim);

        setFromPrecomputedToFullModel();

        is_master_status = true;
        mechRVEstat->setToMasterStatus();
        trspRVEstat->setToMasterStatus();

        MyMatrix ela, dam, ine;
        dam = giveDampingTensor();
        ine = giveInertiaTensor();
        //set to precomputed AFTER calculating inertia and damping
        setToPrecomputed();
        ela = giveStiffnessTensor("secant", ndim);


        macromaterial->setPrecomputedElasticDampingAndInertiaTensors(ela, dam, ine);



        DiscreteTransportRVEMaterial *dtRVEmat = static_cast< DiscreteTransportRVEMaterial * >( trspRVEstat->giveMaterial() );
        DiscreteTrsprtCoupledMaterial *dctm = dynamic_cast< DiscreteTrsprtCoupledMaterial * >( dtRVEmat->giveMasterMaterial() );

        if ( !dctm ) {
            cerr << "Error in " << name << ": transport material in RVE is not DiscreteTrsprtCoupledMaterial" << endl;
            exit(1);
        }
        macromaterial->setMasterMaterial(dctm);

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
MyVector DiscreteCoupledRVEMaterialStatus ::  giveStress(const MyVector &strain, double timeStep) {
    DiscreteCoupledRVEMaterial *dcm = static_cast< DiscreteCoupledRVEMaterial * >( mat );

    unsigned ndim = dcm->giveNumOfDimensions();
    unsigned sizeM = mechRVEstat->giveStrainSize();
    unsigned sizeT = trspRVEstat->giveStrainSize();
    temp_strain.resize(sizeM + sizeT);
    MyVector strainM = MyVector :: Zero(sizeM);
    MyVector strainT = MyVector :: Zero(sizeT + 1); //adding macroscopic pressure
    unsigned i;
    for ( i = 0; i < sizeM; i++ ) {
        strainM [ i ] = temp_strain [ i ] = strain [ i ];
    }
    for ( i = 0; i < sizeT; i++ ) {
        strainT [ i ] = temp_strain [ i + sizeM ] = strain [ i + sizeM ];
    }
    strainT [ sizeT ] = strain [ sizeT + sizeM ]; //macroscopic pressure

    //solve mechanics
    MyVector stressM = mechRVEstat->giveStress(strainM, timeStep);
    if ( !mechRVEstat->isPrecomputed() && trspRVEstat->isPrecomputed() ) {
        setFromPrecomputedToFullModel();                                                                 //when mechanical RVE changes to full version
    }
    double mechBiot = -temp_pressure *dcm->giveBiotCoefficient();  //take pressure strored in element for integration point IDX

    //solve transport
    temp_crackVolume = mechRVEstat->giveCrackVolume();
    updateRateVariables(timeStep);
    MyVector stressT = trspRVEstat->giveStress(strainT, timeStep);


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
MyVector DiscreteCoupledRVEMaterialStatus :: giveInternalSource() const {
    DiscreteCoupledRVEMaterial *dcm = static_cast< DiscreteCoupledRVEMaterial * >( mat );
    unsigned ndim = dcm->giveNumOfDimensions();


    DiscreteCoupledRVEMaterial *dcRVEmat = static_cast< DiscreteCoupledRVEMaterial * >( mat );
    DiscreteTrsprtCoupledMaterial *dtcm = dcRVEmat->giveMasterMaterial();

    double PUCVolume = dcm->givePUCVolume();
    MyVector intS;
    intS.resize(3 * ( ndim - 1 ) + 1);
    unsigned TDoF = 3 * ( ndim - 1 );

    intS [ TDoF ] = -dcm->giveBiotCoefficient() * volStrainRate * 3.;
    intS [ TDoF ] -= temp_crackVolume * pressureRate / ( PUCVolume * dtcm->giveKw() );
    intS [ TDoF ] -= crackVolumeRate / PUCVolume * ( 1. - dcm->giveBiotCoefficient() + ( temp_pressure - dtcm->giveReferencePressure() ) / dtcm->giveKw() );
    intS [ TDoF ] *= dtcm->giveDensity();
    return intS;
}

//////////////////////////////////////////////////////////
MyVector DiscreteCoupledRVEMaterialStatus ::  giveStressWithFrozenIntVars(const MyVector &strain, double timeStep) {
    ( void ) timeStep;
    temp_strain = strain;
    DiscreteCoupledRVEMaterial *dcm = static_cast< DiscreteCoupledRVEMaterial * >( mat );
    unsigned ndim = dcm->giveNumOfDimensions();
    temp_stress = giveStiffnessTensor("secant", ndim) * strain;

    updateRateVariables(timeStep);

    double mechBiot = -temp_pressure *dcm->giveBiotCoefficient();  //take pressure strored in element for integration point IDX
    for ( unsigned i = 0; i < ndim; i++ ) {
        temp_stress [ i ] += mechBiot;
    }



    return temp_stress;
}

//////////////////////////////////////////////////////////
void DiscreteCoupledRVEMaterialStatus ::  giveValues(string code, MyVector & result) const {
    if ( code.compare("crack_volume_density") == 0 ) {
        DiscreteCoupledRVEMaterial *dcm = static_cast< DiscreteCoupledRVEMaterial * >( mat );
        result.resize(1);
        result[0] = temp_crackVolume / dcm->givePUCVolume();
    } else {
        RVEMaterialStatus :: giveValues(code, result);
    }
}

//////////////////////////////////////////////////////////
MyMatrix DiscreteCoupledRVEMaterialStatus :: giveStiffnessTensor(string type, unsigned dimension) const {
    MyMatrix M = mechRVEstat->giveStiffnessTensor(type, dimension);
    MyMatrix T = trspRVEstat->giveStiffnessTensor(type, dimension);
    unsigned sizeM = M.rows();
    unsigned sizeT = T.rows();
    MyMatrix D = MyMatrix :: Zero(sizeM + sizeT, sizeM + sizeT);
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
MyMatrix DiscreteCoupledRVEMaterialStatus :: giveDampingTensor() const {
    if ( is_precomputed ) {
        DiscreteCoupledRVEMaterial *macromat = static_cast< DiscreteCoupledRVEMaterial * >( mat );
        return macromat->givePrecomputedDampingTensor();
    } else {
        MyMatrix M = mechRVEstat->giveDampingTensor();
        MyMatrix T = trspRVEstat->giveDampingTensor();
        unsigned sizeM = M.rows();
        unsigned sizeT = T.rows();
        MyMatrix D = MyMatrix :: Zero(sizeM + sizeT, sizeM + sizeT);
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
MyMatrix DiscreteCoupledRVEMaterialStatus :: giveInertiaTensor() const {
    if ( is_precomputed ) {
        DiscreteCoupledRVEMaterial *macromat = static_cast< DiscreteCoupledRVEMaterial * >( mat );

        return macromat->givePrecomputedInertiaTensor();
    } else {
        MyMatrix M = mechRVEstat->giveInertiaTensor();
        MyMatrix T = trspRVEstat->giveDampingTensor(); //to get correct size
        unsigned sizeM = M.rows();
        unsigned sizeT = T.rows();
        MyMatrix D = MyMatrix :: Zero(sizeM + sizeT, sizeM + sizeT);
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
    cout << "searching for mechanical friends of transport elements across RVEs" << endl;

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
    unsigned ndim = macromat->giveNumOfDimensions();

    //attach mech elems to node numbers
    vector< vector< RigidBodyContact * > >attachedRBC(nodesM->giveSize() );
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
                    insideP(v) = (coord + RVEsize [ v ]);
                    is_inside = false;
                } else if ( coord > RVEsize [ v ] ) {
                    insideP(v) = (coord - RVEsize [ v ]);
                    is_inside = false;
                }
            }
            if ( is_inside ) {
                attachedRBC [ nodesM->giveNodeNumber(rbc->giveNode(p) ) ].push_back(rbc);
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
    Transp1DCoupled *trsp;
    vector< Node * >vertices;
    vector< unsigned >vnums;
    vector< RigidBodyContact * >mechelems, mechelemsold;
    unsigned pold;
    double weight = 0;
    for ( unsigned k = 0; k < elemsT->giveSize(); k++ ) {
        trsp = dynamic_cast< Transp1DCoupled * >( elemsT->giveElement(k) );
        if ( !trsp ) {
            cerr << "Coupled transport elements are required" << endl;
            exit(1);
        }
        vertices = trsp->giveVertices();
        vnums.resize(vertices.size() );
        for ( unsigned p = 0; p < vertices.size(); p++ ) {
            is_inside = true;
            insideP = vertices [ p ]->givePoint();
            for ( unsigned v = 0; v < ndim; v++ ) {
                coord = insideP(v);
                if ( coord < 0. ) {
                    insideP(v) = (coord + RVEsize [ v ]);
                    is_inside = false;
                } else if ( coord > RVEsize [ v ] ) {
                    insideP(v) = (coord - RVEsize [ v ]);
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
void DiscreteCoupledRVEMaterialStatus ::  setEigenStrain(MyVector &x) {
    unsigned sizeM = mechRVEstat->giveStrainSize();
    unsigned sizeT = trspRVEstat->giveStrainSize();

    MyVector eigenstrainM = MyVector :: Zero(sizeM);
    MyVector eigenstrainT = MyVector :: Zero(sizeT);
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
void DiscreteCoupledRVEMaterialStatus ::  setReferenceSystemDirections(MyMatrix r) {
    mechRVEstat->setReferenceSystemDirections(r);
    trspRVEstat->setReferenceSystemDirections(r);
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE COUPLED RVE MATERIAL
//////////////////////////////////////////////////////////

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
    DiscreteCoupledRVEMaterialStatus *newstat = new DiscreteCoupledRVEMaterialStatus(this, e, ipnum, mechRVEmat->givePathToInputFile(), trspRVEmat->givePathToInputFile(), ndim);
    return newstat;
};

//////////////////////////////////////////////////////////
void DiscreteCoupledRVEMaterial :: readFromLine(istringstream &iss) {
    mechRVEmat = new DiscreteMechanicalRVEMaterial();
    trspRVEmat = new DiscreteTransportRVEMaterial();

    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream

    string param;
    bool bdim, bfilem, bfilet, bbiot;
    bdim = bfilem = bfilet = bbiot = false;
    string fileM, fileT;

    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("dim") == 0 ) {
            bdim = true;
            iss >> ndim;
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
void DiscreteCoupledRVEMaterial :: setPrecomputedElasticDampingAndInertiaTensors(MyMatrix ela, MyMatrix dam, MyMatrix ine) {
    precompElastic = ela;
    precompDamping = dam;
    precompInertia = ine;
}
