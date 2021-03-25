#include "material_RVE.h"
#include "model.h"


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// general RVE status
RVEMaterialStatus :: RVEMaterialStatus(RVEMaterial *m, Element *e, fs :: path masterfile) : MaterialStatus(m, e) {
    name = "general RVE mat. status";
    inputfile = masterfile;
    RVE = new Model(false);
}

//////////////////////////////////////////////////////////
RVEMaterialStatus :: ~RVEMaterialStatus() {
    delete RVE;
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

    ndim = RVE->giveDimension();
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
MaterialStatus *RVEMaterial :: giveNewMaterialStatus(Element *e) {
    RVEMaterialStatus *newstat = new RVEMaterialStatus(this, e, inputfile);
    return newstat;
}

//////////////////////////////////////////////////////////
void RVEMaterial :: readFromLine(istringstream &iss) {
    string filename;
    iss >> filename;
    inputfile = GlobPaths :: BASEDIR  / filename;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE TRANSPORT RVE MATERIAL
//////////////////////////////////////////////////////////
DiscreteTransportRVEMaterialStatus :: DiscreteTransportRVEMaterialStatus(RVEMaterial *m, Element *e, fs :: path masterfile) : RVEMaterialStatus(m, e, masterfile) {
    name = "transport RVE mat. status";
}

/////////////////////////////////./////////////////////////
void DiscreteTransportRVEMaterialStatus :: applyEigenStrains() {
    //set eigenstrains and set nonlinearity
    ElementContainer *elems = RVE->giveElements();
    Point normal;
    Vector eigstr(1);
    Transp1D *e;
    TrsprtMaterial *m;
    TrsprtMaterialStatus *s;
    for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
        eigstr *= 0.;
        e = static_cast< Transp1D * >( elems->giveElement(i) );
        s = static_cast< TrsprtMaterialStatus * >( e->giveMatStatus(0) );
        m = static_cast< TrsprtMaterial * >( s->giveMaterial() );
        m->setPermeability(orig_mater_params [ 2 * i ]); //set back original permeability
        m->setPermeability(s->calculatePressureDependentPermeability(aux_params [ 0 ]) ); //calculating pressure depedent conductivity
        m->setParamA(-1.); //switch of linearity
        normal = e->giveNormal();
        for ( unsigned v = 0; v < ndim; v++ ) {
            eigstr [ 0 ] -= temp_strain [ v ] * normal.giveCoord(v);
        }        
        e->giveMatStatus(0)->setEigenStrain(eigstr);
    }
}

/////////////////////////////////./////////////////////////
void DiscreteTransportRVEMaterialStatus :: collectStresses() {
    ElementContainer *elems = RVE->giveElements();
    double volume = 0.;
    Vector factor;
    Transp1D *e;
    TrsprtMaterial *m;
    Point normal;
    for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
        e = static_cast< Transp1D * >( elems->giveElement(i) );
        m = static_cast< TrsprtMaterial * >( e->giveMatStatus(0)->giveMaterial() );
        m->setParamA(orig_mater_params [ 2 * i + 1 ]); //set back a parameter
        normal = e->giveNormal();
        factor  = e->giveArea() * e->giveLength() * e->giveMatStatus(0)->giveTempStress();
        for ( unsigned v = 0; v < ndim; v++ ) {
            temp_stress [ v ] += factor [ 0 ] * normal.giveCoord(v);
        }
        volume += e->giveVolume();
    }
    for ( unsigned v = 0; v < ndim; v++ ) {
        temp_stress [ v ] /= volume;
    }
}

/////////////////////////////////./////////////////////////
Vector DiscreteTransportRVEMaterialStatus :: giveStress(const Vector &strain) {

    temp_strain.resize(giveStrainSize(ndim) );
    for ( unsigned i = 0; i < temp_strain.size(); i++ ) {
        temp_strain [ i ] = strain [ i ];
    }
    temp_strain = addEigenStrain(temp_strain); //macroscopic eigenstrain

    aux_params.resize(strain.size() );
    for ( unsigned i = temp_strain.size(); i < strain.size(); i++ ) {
        aux_params [ i - temp_strain.size() ] = strain [ i ];
    }

    //set eigenstrains
    applyEigenStrains();

    //solve
    RVE->resetTime();
    RVE->giveSolver()->runBeforeEachStep();
    RVE->giveSolver()->solve();


    //collect results
    temp_stress.resize(temp_strain.size() );
    temp_stress *= 0;
    collectStresses();

    //*
    cout << "STRAIN ";
    for(unsigned v=0; v<temp_strain.size(); v++) cout << " " << temp_strain[v];
    cout << endl;
     
    cout << "STRESS ";
    for(unsigned v=0; v<temp_strain.size(); v++) cout << " " << temp_stress[v];
    cout << endl;
    //*/

    return temp_stress;
}

/////////////////////////////////./////////////////////////
Vector DiscreteTransportRVEMaterialStatus :: giveStressWithFrozenIntVars(const Vector &strain) {
    //TODO: this is WRONG, needs fix in future
    return giveStress(strain);
}

/////////////////////////////////./////////////////////////
unsigned DiscreteTransportRVEMaterialStatus :: giveStrainSize(unsigned ndim) const {
    return ndim;
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
        for ( unsigned k = 0; k < jd->giveNumOfMasters(); k++ ) {
            masternode = jd->giveMasterNode(0);
            if ( dynamic_cast< MechDoF * >( masternode ) == nullptr && dynamic_cast< TrsDoF * >( masternode ) == nullptr ) {
                BoundaryCondition *bc;
                vector< int >dBC, nBC;
                dBC.resize(masternode->giveNumberOfDoFs(), funcs->giveSize() ); //todo: warning C4267: 'argument': conversion from 'size_t' to 'const _Ty', possible loss of data
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
    NodeContainer *nodes = RVE->giveNodes();
    BCContainer *bconds = RVE->giveBC();
    ElementContainer *elems = RVE->giveElements();
    ConstraintContainer *constrs = RVE->giveConstraints();

    unsigned ndim = RVE->giveDimension();

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
        MechDoF *pn = new MechDoF(nDoFs);   //?? for transport
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
Matrix DiscreteTransportRVEMaterialStatus :: giveStiffnessTensor(string type, unsigned ndim) const {
    unsigned strain_size = giveStrainSize(ndim);
    Matrix Keff(strain_size, strain_size);

    ElementContainer *elems = RVE->giveElements();
    double volume = 0;
    Point normal;
    Vector n(ndim);
    Transp1D *e;
    volume = 0;
    TrsprtMaterialStatus *tstat;
    for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
        e = static_cast< Transp1D * >( elems->giveElement(i) );
        normal = e->giveNormal();
        for ( unsigned v = 0; v < ndim; v++ ) {
            n [ v ] = normal.giveCoord(v);
        }
        tstat = static_cast< TrsprtMaterialStatus * >( e->giveMatStatus(0) );
        Keff -= dyadicProduct(n, n) * ( e->giveLength() * e->giveArea() * tstat->giveEffectiveConductivity("tangent") );
        volume += e->giveVolume();
    }
    Keff /= volume;

    return Keff;
}

//////////////////////////////////////////////////////////
double DiscreteTransportRVEMaterialStatus :: giveDampingConstant() const {
    ElementContainer *elems = RVE->giveElements();
    Transp1D *e;
    double volume = 0;
    double mass = 0;
    for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
        e = dynamic_cast< Transp1D * >( elems->giveElement(i) );
        if ( e ) {
            mass += e->giveVolume() * e->giveMatStatus(0)->giveDampingConstant();
            volume += e->giveVolume();
        }
    }
    return mass / volume;
};

//////////////////////////////////////////////////////////
void DiscreteTransportRVEMaterialStatus :: init() {
    RVEMaterialStatus :: init();

    ElementContainer *elems = RVE->giveElements();
    aux_params.resize(1); // macro pressure
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
//////////////////////////////////////////////////////////
// DISCRETE TRANSPORT RVE MATERIAL
//////////////////////////////////////////////////////////
MaterialStatus *DiscreteTransportRVEMaterial :: giveNewMaterialStatus(Element *e) {
    DiscreteTransportRVEMaterialStatus *newstat = new DiscreteTransportRVEMaterialStatus(this, e, inputfile);
    return newstat;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE MECHANICAL RVE MATERIAL
//////////////////////////////////////////////////////////
DiscreteMechanicalRVEMaterialStatus :: DiscreteMechanicalRVEMaterialStatus(RVEMaterial *m, Element *e, fs :: path masterfile) : DiscreteTransportRVEMaterialStatus(m, e, masterfile) {
    name = "mechancial RVE mat. status";
}

/////////////////////////////////./////////////////////////
void DiscreteMechanicalRVEMaterialStatus :: applyEigenStrains() {
    //set eigenstrains
    ElementContainer *elems = RVE->giveElements();
    Vector eigstr(ndim);
    unsigned num = 0;
    for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
        eigstr *= 0.;
        for ( unsigned v = 0; v < ndim; v++ ) {
            for ( unsigned r = 0; r < temp_strain.size(); r++ ) {
                eigstr [ v ] -= temp_strain [ r ] * projectors [ v ] [ num ] [ r ];
            }
        }
        num++;
        elems->giveElement(i)->giveMatStatus(0)->setEigenStrain(eigstr);
    }
}

/////////////////////////////////./////////////////////////
void DiscreteMechanicalRVEMaterialStatus :: collectStresses() {
    vector< double >zeta(temp_strain.size() );
    Vector zeta2(9);

    ElementContainer *elems = RVE->giveElements();
    double volume = 0.;
    Vector factor;
    unsigned num = 0;
    RigidBodyContact *e;
    for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
        e = static_cast< RigidBodyContact * >( elems->giveElement(i) );
        factor  = e->giveArea() * e->giveLength() * e->giveMatStatus(0)->giveTempStress();
        for ( unsigned v = 0; v < ndim; v++ ) {
            for ( unsigned r = 0; r < temp_strain.size(); r++ ) {
                temp_stress [ r ] += factor [ v ] * projectors [ v ] [ num ] [ r ];
                if ( v == 0 ) {
                    zeta [ r ] += e->giveArea() * e->giveLength() * projectors [ v ] [ num ] [ r ];
                }
            }
        }
        

        Point normal = e->giveNormal();
        Point xc =  e->giveIPLoc(0) - centroid;
        Point q = cross(xc,normal)*e->giveArea() * e->giveLength() ;
        zeta2[0] += normal.getX() * q.getX();
        zeta2[1] += normal.getY() * q.getY();
        zeta2[2] += normal.getZ() * q.getZ();
        zeta2[3] += normal.getZ() * q.getY();
        zeta2[4] += normal.getY() * q.getZ();
        zeta2[5] += normal.getZ() * q.getX();
        zeta2[6] += normal.getX() * q.getZ();
        zeta2[7] += normal.getY() * q.getX();
        zeta2[8] += normal.getX() * q.getY();




        num++;
        volume += e->giveVolume();
    }
    for ( unsigned v = 0; v < temp_strain.size(); v++ ) {
        temp_stress [ v ] /= volume;
        zeta [ v ] /= volume;
    }


    //*
    cout << "ZETA1 ";
    for ( unsigned v = 9; v < temp_strain.size(); v++ ) {
        cout << " " << zeta [ v ];
    }
    cout << endl;
    cout << "ZETA2 ";
    for ( unsigned v = 0; v < 9; v++ ) {
        cout << " " << zeta2 [ v ]/volume;
    }
    cout << endl;
    cout << "Centroid " << centroid.getX() << " " << centroid.getY() << " " << centroid.getZ() << endl;
    //*//
}


/////////////////////////////////./////////////////////////
unsigned DiscreteMechanicalRVEMaterialStatus :: giveStrainSize(unsigned ndim) const {
    unsigned strain_size = ndim * ndim;
    strain_size += ( ndim == 3 ) ? ndim * ndim : ndim; //in 2D only vector
    return strain_size;
}

//////////////////////////////////////////////////////////
Matrix DiscreteMechanicalRVEMaterialStatus :: giveStiffnessTensor(string type, unsigned ndim) const {
    unsigned strain_size = giveStrainSize(ndim);
    Matrix Keff(strain_size, strain_size);

    ElementContainer *elems = RVE->giveElements();
    double volume = 0;
    Point normal;
    Vector n(ndim);

    RigidBodyContact *e;
    volume = 0;
    Matrix stiff;
    unsigned num = 0;
    for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
        e = static_cast< RigidBodyContact * >( elems->giveElement(i) );
        stiff  = e->giveMatStatus(0)->giveStiffnessTensor(type, ndim);
        for ( unsigned v = 0; v < ndim; v++ ) {
            Keff += dyadicProduct(projectors [ v ] [ num ], projectors [ v ] [ num ]) * e->giveLength() * e->giveArea() * stiff [ v ] [ v ];
        }
        num++;
        volume += e->giveVolume();
    }
    Keff /= volume;

    return Keff;
}

//////////////////////////////////////////////////////////
void DiscreteMechanicalRVEMaterialStatus :: calculateCentroid() {  //only for mechanics
    centroid = Point(0., 0., 0.);

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

    /* not needed, update is done in mechanical projectors
     * //move coordinates
     * NodeContainer *nodes = RVE->giveNodes();
     * for ( unsigned i = 0; i < nodes->giveSize(); i++ ) {
     *  nodes->giveNode(i)->subtructFromPoint(&centroid);
     * }
     + UPDATE
     */
}

//////////////////////////////////////////////////////////
void DiscreteMechanicalRVEMaterialStatus :: init() {
    RVEMaterialStatus :: init();
    calculateCentroid();

    ElementContainer *elems = RVE->giveElements();
    projectors.resize(ndim);
    unsigned strain_size = giveStrainSize(ndim);

    RigidBodyContact *e;
    Point xc;
    Point alphaVec;
    Point normal;
    Vector PQ(strain_size);
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
                    PQ [ 0 ] = normal.getX() * alphaVec.getX();
                    PQ [ 1 ] = normal.getY() * alphaVec.getY();
                    PQ [ 2 ] = normal.getX() * alphaVec.getY();
                    PQ [ 3 ] = normal.getY() * alphaVec.getX();
                    double factor2D = xc.getX() * alphaVec.getY() - xc.getY() * alphaVec.getX();
                    PQ [ 4 ] = factor2D * normal.getX();
                    PQ [ 5 ] = factor2D * normal.getY();
                } else {
                    PQ [ 0 ] = normal.getX() * alphaVec.getX();
                    PQ [ 1 ] = normal.getY() * alphaVec.getY();
                    PQ [ 2 ] = normal.getZ() * alphaVec.getZ();
                    PQ [ 3 ] = normal.getZ() * alphaVec.getY();
                    PQ [ 4 ] = normal.getY() * alphaVec.getZ();
                    PQ [ 5 ] = normal.getZ() * alphaVec.getX();
                    PQ [ 6 ] = normal.getX() * alphaVec.getZ();
                    PQ [ 7 ] = normal.getY() * alphaVec.getX();
                    PQ [ 8 ] = normal.getX() * alphaVec.getY();
                    Point factor3D = cross(xc, alphaVec);
                    PQ [ 9 ]  = normal.getX() * factor3D.getX();
                    PQ [ 10 ] = normal.getY() * factor3D.getY();
                    PQ [ 11 ] = normal.getZ() * factor3D.getZ();
                    PQ [ 12 ] = normal.getZ() * factor3D.getY();
                    PQ [ 13 ] = normal.getY() * factor3D.getZ();
                    PQ [ 14 ] = normal.getZ() * factor3D.getX();
                    PQ [ 15 ] = normal.getX() * factor3D.getZ();
                    PQ [ 16 ] = normal.getY() * factor3D.getX();
                    PQ [ 17 ] = normal.getX() * factor3D.getY();
                }
                projectors [ v ].push_back(PQ);
            }
        }
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE MECHANICAL RVE MATERIAL
//////////////////////////////////////////////////////////
MaterialStatus *DiscreteMechanicalRVEMaterial :: giveNewMaterialStatus(Element *e) {
    DiscreteMechanicalRVEMaterialStatus *newstat = new DiscreteMechanicalRVEMaterialStatus(this, e, inputfile);
    return newstat;
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE RVE TRANSPORT MATERIAL STATUS PRE-COMPUTED
//////////////////////////////////////////////////////////
DiscreteTransportRVEMaterialPrecomputedStatus :: DiscreteTransportRVEMaterialPrecomputedStatus(RVEMaterial *m, Element *e, fs :: path masterfile) : DiscreteTransportRVEMaterialStatus(m, e, masterfile) {
    name = "transport RVE precomputed mat. status";
    is_master_status = false;
    temp_nonlin = 1;
}

/////////////////////////////////./////////////////////////
Vector DiscreteTransportRVEMaterialPrecomputedStatus :: giveStress(const Vector &strain) {
    double macro_pressure = strain [ strain.size() - 1 ];
    temp_strain.resize(giveStrainSize(ndim) );
    for ( unsigned i = 0; i < temp_strain.size(); i++ ) {
        temp_strain [ i ] = strain [ i ];
    }

    DiscreteTransportRVEMaterialPrecomputed *macromaterial = static_cast< DiscreteTransportRVEMaterialPrecomputed * >( mat );
    temp_nonlin = macromaterial->giveMasterStatus()->calculatePressureDependentPermeability(macro_pressure) / macromaterial->giveMasterMaterial()->givePermeability();
    temp_stress = matrix_vector_multiply(macromaterial->givePrecomputedConductivity(), temp_strain) * temp_nonlin;

    return temp_stress;
}

/////////////////////////////////./////////////////////////
Vector DiscreteTransportRVEMaterialPrecomputedStatus :: giveStressWithFrozenIntVars(const Vector &strain) {
    return giveStress(strain);
}

//////////////////////////////////////////////////////////
Matrix DiscreteTransportRVEMaterialPrecomputedStatus :: giveStiffnessTensor(string type, unsigned ndim) const {
    DiscreteTransportRVEMaterialPrecomputed *macromaterial = static_cast< DiscreteTransportRVEMaterialPrecomputed * >( mat );
    return macromaterial->givePrecomputedConductivity() * temp_nonlin;
}

//////////////////////////////////////////////////////////
double DiscreteTransportRVEMaterialPrecomputedStatus :: giveDampingConstant() const {
    DiscreteTransportRVEMaterialPrecomputed *macromaterial = static_cast< DiscreteTransportRVEMaterialPrecomputed * >( mat );
    return macromaterial->givePrecomputedCapacity();
};

//////////////////////////////////////////////////////////
void DiscreteTransportRVEMaterialPrecomputedStatus :: init() {
    DiscreteTransportRVEMaterialPrecomputed *macromaterial = static_cast< DiscreteTransportRVEMaterialPrecomputed * >( mat );
    Matrix stiff = macromaterial->givePrecomputedConductivity();
    if ( stiff.size() == 0 ) {
        is_master_status = true;
        DiscreteTransportRVEMaterialStatus :: init();

        //TODO: if this is needed, one should call it only once and distribute the result among all identical RVEs
        cout << "Precomputing pressure fields on RVE" << endl;

        unsigned ndim = RVE->giveDimension();
        stiff = Matrix(ndim, ndim);
        Vector strainDoFs(ndim + 1);
        Vector stress;
        for ( unsigned i = 0; i < ndim; i++ ) {
            strainDoFs [ i ] = 1.;
            stress = DiscreteTransportRVEMaterialStatus :: giveStress(strainDoFs);
            strainDoFs [ i ] = 0.;
            for ( unsigned j = 0; j < ndim; j++ ) {
                stiff [ i ] [ j ] = stress [ j ];
            }
        }
        double c = DiscreteTransportRVEMaterialStatus :: giveDampingConstant();

        TrsprtMaterialStatus *status = static_cast< TrsprtMaterialStatus * >( RVE->giveElements()->giveElement(0)->giveMatStatus(0) );
        TrsprtMaterial *material = static_cast< TrsprtMaterial * >( RVE->giveElements()->giveElement(0)->giveMaterial() );
        macromaterial->setPrecomputedConductivityAndCapacityAndMasterMaterial(stiff, c, status, material);
    } else   {
        RVEMaterialStatus :: init();
    }
}

//////////////////////////////////////////////////////////
void DiscreteTransportRVEMaterialPrecomputedStatus :: update() {
    if ( is_master_status ) {
        RVEMaterialStatus :: update();
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE RVE TRANSPORT MATERIAL PRE-COMPUTED
//////////////////////////////////////////////////////////
MaterialStatus *DiscreteTransportRVEMaterialPrecomputed :: giveNewMaterialStatus(Element *e) {
    DiscreteTransportRVEMaterialPrecomputedStatus *newstat = new DiscreteTransportRVEMaterialPrecomputedStatus(this, e, inputfile);
    return newstat;
}

//////////////////////////////////////////////////////////
void DiscreteTransportRVEMaterialPrecomputed :: setPrecomputedConductivityAndCapacityAndMasterMaterial(Matrix lam, double c, TrsprtMaterialStatus *masterS,  TrsprtMaterial *masterM) {
    conductivity = lam;
    capacity = c;
    masterStatus = masterS;
    masterMaterial = masterM;

    const char *buffer = "precomputed_conductivity.out";
    ofstream outputfile( ( masterModel->giveResultDirectory() / buffer ).string() );
    if ( outputfile.is_open() ) {
        outputfile << std :: scientific;
        outputfile.precision(10);

        unsigned size = conductivity.numRows();
        for ( unsigned i = 0; i < size; i++ ) {
            for ( unsigned j = 0; j < size; j++ ) {
                if ( j > 0 ) {
                    outputfile << "\t";
                }
                outputfile << conductivity [ i ] [ j ];
            }
            outputfile << endl;
        }
        outputfile.close();
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE COUPLED RVE MATERIAL STATUS
//////////////////////////////////////////////////////////
DiscreteCoupledRVEMaterialStatus ::  DiscreteCoupledRVEMaterialStatus(Material *m, Element *e, fs :: path masterfileM, fs :: path masterfileT): MaterialStatus(m,e){
    name = "coupled discrete RVE mat. status";
    DiscreteCoupledRVEMaterial* coupledm = dynamic_cast<DiscreteCoupledRVEMaterial*>(m);
    if (coupledm){
        mechRVEstat = new DiscreteMechanicalRVEMaterialStatus(coupledm->giveMechanicalRVEmat(),e,masterfileM);
        trspRVEstat = new DiscreteTransportRVEMaterialStatus(coupledm->giveTransportRVEmat(),e,masterfileT);
    }else{
        cerr << "DiscreteCoupledRVEMaterialStatus accepts only DiscreteCoupledRVEMaterial" << endl;
        exit(1);
    }
}

//////////////////////////////////////////////////////////
DiscreteCoupledRVEMaterialStatus ::  ~DiscreteCoupledRVEMaterialStatus(){
    delete mechRVEstat;
    delete trspRVEstat;
}

//////////////////////////////////////////////////////////
void DiscreteCoupledRVEMaterialStatus ::  init(){
    mechRVEstat->init();
    trspRVEstat->init();

    findFriends();
}

//////////////////////////////////////////////////////////
void DiscreteCoupledRVEMaterialStatus ::  update(){
    mechRVEstat->update();
    trspRVEstat->update();
} 

//////////////////////////////////////////////////////////
Vector DiscreteCoupledRVEMaterialStatus ::  giveStress(const Vector &strain){
    unsigned sizeM = mechRVEstat->giveStrainSize();
    unsigned sizeT = trspRVEstat->giveStrainSize();
    temp_strain.resize(sizeM+sizeT);
    Vector strainM(sizeM);
    Vector strainT(sizeT+1); //adding macroscopic pressure
    unsigned i;
    for(i=0; i<sizeM; i++){
        strainM[i] = temp_strain[i] = strain[i];
    }
    for(i=0; i<sizeT; i++){
        strainT[i] = temp_strain[i+sizeM] = strain[i+sizeM];
    }
    strainT[sizeT] = strain[sizeT+sizeM]; //macroscopic pressure

    Vector stressM = mechRVEstat->giveStress(strainM);
    Vector stressT = trspRVEstat->giveStress(strainT); 

    temp_stress.resize(sizeM+sizeT);
    for(i=0; i<sizeM; i++){
        temp_stress[i] = stressM[i];
    }
    for(i=0; i<sizeT; i++){
        temp_stress[i+sizeM] = stressT[i];
    }    
    return temp_stress;
}

//////////////////////////////////////////////////////////
Vector DiscreteCoupledRVEMaterialStatus ::  giveStressWithFrozenIntVars(const Vector &strain){
    return giveStress(strain);
}

//////////////////////////////////////////////////////////
double DiscreteCoupledRVEMaterialStatus ::  giveValue(string code){
    return 0;
}

//////////////////////////////////////////////////////////
Matrix DiscreteCoupledRVEMaterialStatus :: giveStiffnessTensor(string type, unsigned dimension) const {

    Matrix M = mechRVEstat->giveStiffnessTensor(type,dimension);
    Matrix T = trspRVEstat->giveStiffnessTensor(type,dimension);
    unsigned sizeM = M.numRows();
    unsigned sizeT = T.numRows();
    Matrix D(sizeM+sizeT,sizeM+sizeT);
    for (unsigned i=0; i<sizeM; i++){
        for (unsigned j=0; j<sizeM; j++){ 
            D[i][j] = M[i][j];
        }
    }
    for (unsigned i=0; i<sizeT; i++){
        for (unsigned j=0; j<sizeT; j++){ 
            D[i+sizeM][j+sizeM] = T[i][j];
        }
    }
    return D;    
}

//////////////////////////////////////////////////////////
void DiscreteCoupledRVEMaterialStatus :: findFriends() {

    cout << "searching for mechanical friends of transport elements" << endl;

    unsigned ndim = mechRVEstat->giveWholeRVE()->giveDimension();
    ElementContainer *elemsM = mechRVEstat->giveWholeRVE()->giveElements();
    ElementContainer *elemsT = trspRVEstat->giveWholeRVE()->giveElements();
    NodeContainer *nodesM = mechRVEstat->giveWholeRVE()->giveNodes();
    NodeContainer *nodesT = trspRVEstat->giveWholeRVE()->giveNodes();

    vector<double> RVEsize;
    PBlockContainer* pblocksM =  mechRVEstat->giveWholeRVE()->givePBlockContainer();
    
    for(unsigned i=0; i<pblocksM->giveSize(); i++){
        MechanicalPeriodicBC* mpb = dynamic_cast< MechanicalPeriodicBC* >(pblocksM->givePBlock(i));
        if (mpb){
            RVEsize = mpb->giveDimensions();
        }
    }

    //attach mech elems to node numbers
    vector< vector < RigidBodyContact * > > attachedRBC(nodesM->giveSize());    
    RigidBodyContact *rbc;    
    Point insideP;
    Node* foundN;
    double coord, dist;
    bool is_inside;    
    for ( unsigned k = 0; k < elemsM->giveSize(); k++ ) {
        rbc = static_cast< RigidBodyContact * >( elemsM->giveElement(k) );
        for( unsigned p=0; p<2; p++){
            is_inside = true;
            insideP = rbc->giveNode(p)->givePoint();
            for(unsigned v=0; v<ndim; v++){
                coord = insideP.giveCoord(v);
                if(coord<0.) { 
                    insideP.setCoord(v, coord + RVEsize[v]);
                    is_inside = false;
                }
                else if(coord>RVEsize[v]) { 
                    insideP.setCoord(v, coord - RVEsize[v]);
                    is_inside = false;
                }
            }
            if (is_inside){
                attachedRBC[nodesM->giveNodeNumber(rbc->giveNode(p))].push_back(rbc);
            }else{
                foundN = nodesM->findClosestMechanicalNode(insideP,&dist);
                if(dist>1e-10) {
                    cerr << "DiscreteCoupledRVEMaterialStatus Error: searching for periodic image of mechanical node failed, distance to closest point is " << dist << endl;
                    exit(1);
                }
                attachedRBC[nodesM->giveNodeNumber(foundN)].push_back(rbc);
            }            
        }
    }

    //attach mech elems to transport elements
    Transp1DCoupled *trsp;    
    vector<Node*> vertices;
    vector<unsigned> vnums;
    vector< RigidBodyContact* > mechelems, mechelemsold;
    unsigned pold;
    double weight = 0;
    for ( unsigned k = 0; k < elemsT->giveSize(); k++ ) {
        trsp = static_cast< Transp1DCoupled * >( elemsT->giveElement(k) );
        vertices = trsp->giveVertices();
        vnums.resize(vertices.size());
        for( unsigned p=0; p<vertices.size(); p++){
            is_inside = true;
            insideP = vertices[p]->givePoint();
            for(unsigned v=0; v<ndim; v++){
                coord = insideP.giveCoord(v);
                if(coord<0.) { 
                    insideP.setCoord(v, coord + RVEsize[v]);
                    is_inside = false;
                }
                else if(coord>RVEsize[v]) { 
                    insideP.setCoord(v, coord - RVEsize[v]);
                    is_inside = false;
                }
            }
            if (is_inside){
                vnums[p] = nodesT->giveNodeNumber(vertices[p]);
            }else{
                foundN = nodesT->findClosestAuxiliaryNode(insideP,&dist);
                if(dist>1e-10) {
                    cerr << "DiscreteCoupledRVEMaterialStatus Error: searching for periodic image of transport node failed, distance to closest point is " << dist << endl;
                    exit(1);
                }
                vnums[p] = nodesT->giveNodeNumber(foundN);
            }  
        }          

        mechelemsold = attachedRBC[vnums.back()];
        pold = vnums.size()-1;
        for( unsigned p = 0; p < vnums.size(); p++){
            mechelems = attachedRBC[vnums[p]];
            for( auto me: mechelemsold){
                if( find(mechelems.begin(), mechelems.end(), me) != mechelems.end() ){
                    if (ndim==2)     weight = me->giveArea();
                    else if(ndim==3) weight = ( trsp->giveIPLoc(0) - ( vertices [ p ]->givePoint() + vertices [ pold ]->givePoint() ) / 2. ).norm();
                    trsp->addNewFriend(me, weight );
                    break;
                }
            }            
            mechelemsold = mechelems;
            pold = p;
        }
    }
}

//////////////////////////////////////////////////////////
void DiscreteCoupledRVEMaterialStatus ::  setEigenStrain(Vector &x){   
    unsigned sizeM = mechRVEstat->giveStrainSize();
    unsigned sizeT = trspRVEstat->giveStrainSize();

    Vector eigenstrainM(sizeM);
    Vector eigenstrainT(sizeT);
    unsigned i;
    for(i=0; i<sizeM; i++){
        eigenstrainM[i] = x[i];
    }
    for(i=0; i<sizeT; i++){
        eigenstrainT[i] = x[i+sizeM];
    }
    mechRVEstat -> setEigenStrain(eigenstrainM);
    trspRVEstat -> setEigenStrain(eigenstrainT);
}

//////////////////////////////////////////////////////////
string DiscreteCoupledRVEMaterialStatus :: giveLineToSave()const{
    return mechRVEstat->giveLineToSave() + trspRVEstat->giveLineToSave();
}

//////////////////////////////////////////////////////////
void DiscreteCoupledRVEMaterialStatus :: setID(unsigned i){
    MaterialStatus::setID(i);
    mechRVEstat->setID(i);
    trspRVEstat->setID(i);
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE COUPLED RVE MATERIAL
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
DiscreteCoupledRVEMaterial :: ~DiscreteCoupledRVEMaterial(){
    delete mechRVEmat;
    delete trspRVEmat;
}

//////////////////////////////////////////////////////////
MaterialStatus* DiscreteCoupledRVEMaterial ::giveNewMaterialStatus(Element *e){
    DiscreteCoupledRVEMaterialStatus *newstat = new DiscreteCoupledRVEMaterialStatus(this, e, mechRVEmat->givePathToInputFile(), trspRVEmat->givePathToInputFile());
    return newstat;
};

//////////////////////////////////////////////////////////
void DiscreteCoupledRVEMaterial :: readFromLine(istringstream &iss) {
    mechRVEmat = new DiscreteMechanicalRVEMaterial();
    trspRVEmat = new DiscreteTransportRVEMaterial();
    mechRVEmat->readFromLine(iss);
    trspRVEmat->readFromLine(iss);
}


/*
 * //// TO BE DELETED
 *
 *
 * //////////////////////////////////////////////////////////
 * //////////////////////////////////////////////////////////
 * // DISCRETE RVE MATERIAL
 * DiscreteRVEMaterialStatus :: DiscreteRVEMaterialStatus(DiscreteRVEMaterial *m, Element *e, fs :: path masterfile) : RVEMaterialStatus(m, e, masterfile) {
 *  name = "transport RVE mat. status";
 * }
 *
 * /////////////////////////////////./////////////////////////
 * Vector DiscreteRVEMaterialStatus :: giveStress(const Vector &strain) {
 *
 *  double macro_pressure = 0;
 *  if ( active_transport ) {
 *      temp_strain.resize(strain.size()-1);
 *      for (unsigned i = 0; i<temp_strain.size(); i++) temp_strain[i] = strain[i];
 *      macro_pressure = strain[strain.size()-1];
 *  } else {
 *      temp_strain = strain;
 *  }
 *
 *  unsigned ndim = RVE->giveDimension();
 *
 *  unsigned strain_size = ndim * ndim;
 *  strain_size += ( ndim == 3 ) ? ndim * ndim : ndim; //in 2D only vector
 *  if ( !active_mechanics ) {
 *      strain_size = 0;
 *  }
 *
 *  unsigned pres_size = ndim;
 *  if ( !active_transport ) {
 *      pres_size = 0;
 *  }
 *
 *  //set eigenstrains
 *  ElementContainer *elems = RVE->giveElements();
 *  Point normal;
 *
 *  if ( active_mechanics ) {
 *      RigidBodyContact *e;
 *      Vector eigstr(ndim);
 *      unsigned num = 0;
 *      for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
 *          e = dynamic_cast< RigidBodyContact * >( elems->giveElement(i) );
 *          if ( !e ) {
 *              continue;
 *          }
 *          eigstr *= 0.;
 *          for ( unsigned v = 0; v < ndim; v++ ) {
 *              for ( unsigned r = 0; r < strain_size; r++ ) {
 *                  eigstr [ v ] -= temp_strain [ r ] * projectors [ v ] [ num ] [ r ];
 *              }
 *          }
 *          num++;
 *          e->giveMatStatus(0)->setEigenStrain(eigstr);
 *      }
 *  }
 *
 *  vector < double >  initial_permeability, initial_a;
 *  if ( active_transport ) {
 *      Transp1D *e;
 *      TrsprtMaterialStatus *s;
 *      TrsprtMaterial *m;
 *      Vector eigstr(1);
 *      for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
 *          e = dynamic_cast< Transp1D * >( elems->giveElement(i) );
 *          if ( !e ) {
 *              continue;
 *          }
 *          s = static_cast< TrsprtMaterialStatus* >( e->giveMatStatus(0) );
 *          m = static_cast< TrsprtMaterial* >( s->giveMaterial() );
 *          initial_permeability.push_back( m->givePermeability());
 *          initial_a.push_back( m->giveParamA());
 *          m->setPermeability( s->calculatePressureDependentPermeability(macro_pressure) ); //calculating pressure depedent conductivity
 *          m->setParamA(-1.); //switch of linearity
 *          normal = e->giveNormal();
 *          eigstr *= 0.;
 *          for ( unsigned v = 0; v < ndim; v++ ) {
 *              eigstr [ 0 ] -= temp_strain [ strain_size + v ] * normal.giveCoord(v);
 *          }
 *          s->setEigenStrain(eigstr);
 *      }
 *  }
 *
 *  //solve
 *  RVE->resetTime();
 *  //RVE->giveSolver()->updateSteadyStateMatrix(); //needed when using linear solver - NO, internal forces automatically takes care, it can be switched on optionally
 *  RVE->giveSolver()->runBeforeEachStep();
 *  RVE->giveSolver()->solve();
 *
 *  //REMOVE MACROPRESSURE
 *
 *  //collect results
 *  temp_stress.resize(strain_size + pres_size);
 *  temp_stress *= 0;
 *
 *  vector<double> zeta(strain_size + pres_size);
 *
 *  if ( active_mechanics ) {
 *      double volume = 0.;
 *      RigidBodyContact *e;
 *      Vector factor;
 *      unsigned num = 0;
 *      for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
 *          e = dynamic_cast< RigidBodyContact * >( elems->giveElement(i) );
 *          if ( e ) {
 *              factor  = e->giveArea() * e->giveLength() * e->giveMatStatus(0)->giveTempStress();
 *              for ( unsigned v = 0; v < ndim; v++ ) {
 *                  for ( unsigned r = 0; r < strain_size; r++ ) {
 *                      temp_stress [ r ] += factor [ v ] * projectors [ v ] [ num ] [ r ];
 *                      if(v==0) zeta [ r ] += e->giveArea() * e->giveLength() * projectors [ v ] [ num ] [ r ];
 *                  }
 *
 *              }
 *              num++;
 *              volume += e->giveVolume();
 *          }
 *      }
 *      for ( unsigned v = 0; v < strain_size; v++ ) {
 *          temp_stress [ v ] /= volume;
 *          zeta [ v ] /= volume;
 *      }
 *  }
 *
 *
 *  if ( active_transport ) {
 *      unsigned k = 0;
 *      double volume = 0.;
 *      double factor;
 *      Transp1D *e;
 *      TrsprtMaterial *m;
 *      for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
 *          e = dynamic_cast< Transp1D * >( elems->giveElement(i) );
 *          if ( e ) {
 *              normal = e->giveNormal();
 *              factor = e->giveArea() * e->giveLength() * e->giveMatStatus(0)->giveTempStress() [ 0 ];
 *              for ( unsigned v = 0; v < ndim; v++ ) {
 *                  temp_stress [ strain_size + v ] += factor * normal.giveCoord(v);
 *              }
 *              volume += e->giveVolume();
 *
 *              m = static_cast< TrsprtMaterial* >( e->giveMatStatus(0)->giveMaterial() );
 *              m->setPermeability( initial_permeability[k] );
 *              m->setParamA( initial_a[k] ); //switch of nonlinearity
 *              k ++;
 *          }
 *      }
 *      for ( unsigned v = 0; v < ndim; v++ ) {
 *          temp_stress [ strain_size + v ] /= volume;
 *      }
 *  }
 *
 *  return temp_stress;
 * }
 *
 * /////////////////////////////////./////////////////////////
 * Vector DiscreteRVEMaterialStatus :: giveStressWithFrozenIntVars(const Vector &strain) {
 *  //TODO: this is WRONG, needs fix in future
 *  return giveStress(strain);
 * }
 *
 * //////////////////////////////////////////////////////////
 * Matrix DiscreteRVEMaterialStatus :: giveStiffnessTensor(string type, unsigned ndim) const {
 *  unsigned strain_size = ndim * ndim;
 *  strain_size += ( ndim == 3 ) ? ndim * ndim : ndim; //in 2D only vector
 *  if ( !active_mechanics ) {
 *      strain_size = 0;
 *  }
 *  Matrix Keff(strain_size, strain_size);
 *
 *  unsigned pres_size = ndim;
 *  if ( !active_transport ) {
 *      pres_size = 0;
 *  }
 *  Matrix Leff(pres_size, pres_size);
 *
 *
 *  ElementContainer *elems = RVE->giveElements();
 *  double volume = 0;
 *  Point normal;
 *  Vector n(ndim);
 *
 *  if ( active_mechanics ) {
 *      RigidBodyContact *e;
 *      volume = 0;
 *      Matrix stiff;
 *      unsigned num = 0;
 *      for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
 *          e = dynamic_cast< RigidBodyContact * >( elems->giveElement(i) );
 *          if ( e ) {
 *              stiff  = e->giveMatStatus(0)->giveStiffnessTensor(type, ndim);
 *              for ( unsigned v = 0; v < ndim; v++ ) {
 *                  Keff += dyadicProduct(projectors [ v ] [ num ], projectors [ v ] [ num ]) * e->giveLength() * e->giveArea() * stiff [ v ] [ v ];
 *              }
 *              num++;
 *              volume += e->giveVolume();
 *          }
 *      }
 *      Keff /= volume;
 *  }
 *
 *
 *  if ( active_transport ) {
 *      Transp1D *e;
 *      volume = 0;
 *      for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
 *          e = dynamic_cast< Transp1D * >( elems->giveElement(i) );
 *          if ( e ) {
 *              normal = e->giveNormal();
 *              n [ 0 ] = normal.getX();
 *              n [ 1 ] = normal.getY();
 *              if ( ndim == 3 ) {
 *                  n [ 2 ] = normal.getZ();
 *              }
 *              Leff += dyadicProduct(n, n) * ( e->giveLength() * e->giveArea() * e->giveMatStatus(0)->giveStiffnessTensor(type, ndim) [ 0 ] [ 0 ] );
 *              volume += e->giveVolume();
 *          }
 *      }
 *      Leff /= volume;
 *  }
 *
 *  //Keff.print();
 *
 *  if ( !active_transport ) {
 *      return Keff;
 *  }
 *  if ( !active_mechanics ) {
 *      return Leff;
 *  }
 *
 *  Matrix KL(strain_size + pres_size, strain_size + pres_size);
 *  for ( unsigned r = 0; r < strain_size; r++ ) {
 *      for ( unsigned s = 0; s < strain_size; s++ ) {
 *          KL [ r ] [ s ] = Keff [ r ] [ s ];
 *      }
 *  }
 *
 *  for ( unsigned r = 0; r < pres_size; r++ ) {
 *      for ( unsigned s = 0; s < pres_size; s++ ) {
 *          KL [ strain_size + r ] [ strain_size + s ] = Leff [ r ] [ s ];
 *      }
 *  }
 *
 *  return KL;
 * }
 *
 * //////////////////////////////////////////////////////////
 * double DiscreteRVEMaterialStatus :: giveDampingConstant() const {
 *  ElementContainer *elems = RVE->giveElements();
 *  Transp1D *e;
 *  double volume = 0;
 *  double mass = 0;
 *  for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
 *      e = dynamic_cast< Transp1D * >( elems->giveElement(i) );
 *      if ( e ) {
 *          mass += e->giveVolume() * e->giveMatStatus(0)->giveDampingConstant();
 *          volume += e->giveVolume();
 *      }
 *  }
 *  return mass / volume;
 * };
 *
 * //////////////////////////////////////////////////////////
 * void DiscreteRVEMaterialStatus :: generateRandomFixedBC() {
 *
 *  BCContainer *bconds = RVE->giveBC();
 *  FunctionContainer *funcs = RVE->giveFunctions();
 *  ConstraintContainer *constrs = RVE->giveConstraints();
 *
 *  //mechanics
 *  JointDoF *jd;
 *  Node * masternode;
 *  active_mechanics = false;
 *  for ( unsigned j = 0; j < constrs->giveSize(); j++ ) {
 *      if (active_mechanics) break;
 *      jd = constrs->giveConstraint(j);
 *      for (unsigned k=0; k<jd->giveNumOfMasters(); k++){
 *          masternode = jd->giveMasterNode(0);
 *          if ( masternode->doesMechanics() && ( dynamic_cast< MechDoF * >( masternode ) == nullptr ) ) {
 *              BoundaryCondition *bc;
 *              vector< int >dBC, nBC;
 *              dBC.resize(masternode->giveNumberOfDoFs(), funcs->giveSize() ); //todo: warning C4267: 'argument': conversion from 'size_t' to 'const _Ty', possible loss of data
 *              nBC.resize(masternode->giveNumberOfDoFs(), -1);
 *              bc = new BoundaryCondition(masternode, dBC, nBC);
 *              bconds->addBoundaryCondition(bc);
 *
 *              active_mechanics = true;
 *              break;
 *          }
 *      }
 *  }
 *
 *  //transport
 *  active_transport = false;
 *  for ( unsigned j = 0; j < constrs->giveSize(); j++ ) {
 *      if (active_transport) break;
 *      jd = constrs->giveConstraint(j);
 *      for (unsigned k=0; k<jd->giveNumOfMasters(); k++){
 *          masternode = jd->giveMasterNode(0);
 *          if ( masternode->doesTransport() && ( dynamic_cast< TrsDoF * >( masternode ) == nullptr ) ) {
 *              BoundaryCondition *bc;
 *              vector< int >dBC, nBC;
 *              dBC.resize(masternode->giveNumberOfDoFs(), funcs->giveSize() ); //todo: warning C4267: 'argument': conversion from 'size_t' to 'const _Ty', possible loss of data
 *              nBC.resize(masternode->giveNumberOfDoFs(), -1);
 *              bc = new BoundaryCondition(masternode, dBC, nBC);
 *              bconds->addBoundaryCondition(bc);
 *
 *              active_transport = true;
 *              break;
 *          }
 *      }
 *  }
 *  //active_transport = true;
 *
 *  //add constant function
 *  vector< double >x, y;
 *  x.resize(1, 0);
 *  y.resize(1, 0);
 *  PieceWiseLinearFunction *newf = new PieceWiseLinearFunction(x, y);
 *  funcs->addFunction(newf);
 * }
 *
 *
 *
 * //////////////////////////////////////////////////////////
 * void DiscreteRVEMaterialStatus :: generateVolumetricAverageBC() {
 *  //this function applies volumetric constraint
 *  //it is not really wise to use it as it leads to relatively full stiffness matrix after transfer to constraint space
 *  NodeContainer *nodes = RVE->giveNodes();
 *  BCContainer *bconds = RVE->giveBC();
 *  ElementContainer *elems = RVE->giveElements();
 *  ConstraintContainer *constrs = RVE->giveConstraints();
 *
 *  unsigned ndim = RVE->giveDimension();
 *
 *  VolumetricAverage *va;
 *  unsigned funsize = RVE->giveFunctions()->giveSize();
 *  vector< double >x, y;
 *  x.resize(1, 0);
 *  y.resize(1, 0);
 *  volumAverFunc = new PieceWiseLinearFunction(x, y);
 *  RVE->giveFunctions()->addFunction(volumAverFunc);
 *  vector< unsigned >dirs;
 *
 *  vector< Node * >vm;
 *
 *  //mechanics
 *  for ( unsigned n = 0; n < nodes->giveSize(); n++ ) {
 *      if ( nodes->giveNode(n)->doesMechanics() && ( dynamic_cast< MechDoF * >( nodes->giveNode(n) ) == nullptr ) ) {
 *          vm.push_back(nodes->giveNode(n) );
 *      }
 *  }
 *  if ( vm.size() > 0 ) {
 *      unsigned nDoFs = 3;
 *      if ( ndim == 3 ) {
 *          nDoFs = 6;
 *      }
 *      MechDoF *pn = new MechDoF(nDoFs);
 *      nodes->addNode(pn);
 *
 *      dirs.resize(vm.size() );
 *
 *      for ( unsigned vi = 0; vi < nDoFs; vi++ ) {
 *          fill(dirs.begin(), dirs.end(), vi);
 *          va = new VolumetricAverage(vm, dirs, pn, vi, elems, constrs);
 *          constrs->addConstraint(va);
 *      }
 *
 *      BoundaryCondition *bc;
 *      vector< int >dBC, nBC;
 *      dBC.resize(nDoFs, funsize);
 *      nBC.resize(nDoFs, -1);
 *      bc = new BoundaryCondition(pn, dBC, nBC);
 *      bconds->addBoundaryCondition(bc);
 *
 *      active_mechanics = true;
 *  } else   {
 *      active_mechanics = false;
 *  }
 *
 *  //transport
 *  vm.clear();
 *  for ( unsigned n = 0; n < nodes->giveSize(); n++ ) {
 *      if ( nodes->giveNode(n)->doesTransport() && ( dynamic_cast< TrsDoF * >( nodes->giveNode(n) ) == nullptr ) ) {
 *          vm.push_back(nodes->giveNode(n) );
 *      }
 *  }
 *  if ( vm.size() > 0 ) {
 *      TrsDoF *tn = new TrsDoF(1);
 *      nodes->addNode(tn);
 *
 *      dirs.resize(vm.size() );
 *      va = new VolumetricAverage(vm, dirs, tn, 0, elems, constrs);
 *      constrs->addConstraint(va);
 *
 *      BoundaryCondition *bc;
 *      vector< int >dBC, nBC;
 *      dBC.resize(1, funsize);
 *      nBC.resize(1, -1);
 *      bc = new BoundaryCondition(tn, dBC, nBC);
 *      bconds->addBoundaryCondition(bc);
 *
 *      active_transport = true;
 *  } else   {
 *      active_transport = false;
 *  }
 * }
 *
 * //////////////////////////////////////////////////////////
 * void DiscreteRVEMaterialStatus :: calculateCentroid() {
 *  centroid = Point(0., 0., 0.);
 *
 *  if ( !active_mechanics ) {
 *      return;
 *  }
 *
 *  //find centroid
 *  ElementContainer *elems = RVE->giveElements();
 *
 *  //mechanics
 *  RigidBodyContact *e;
 *  DisMechMaterialStatus *ms;
 *  double weight;
 *  double totalweight = 0;
 *  for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
 *      e = dynamic_cast< RigidBodyContact * >( elems->giveElement(i) );
 *      if ( e ) {
 *          ms = static_cast< DisMechMaterialStatus * >( e->giveMatStatus(0) );
 *          weight = e->giveLength() * e->giveArea() * ms->giveDensity();
 *          ;
 *          totalweight += weight;
 *          centroid += e->giveIPLoc(0) * weight;
 *      }
 *  }
 *  centroid /= totalweight;
 *
 * }
 *
 * //////////////////////////////////////////////////////////
 * void DiscreteRVEMaterialStatus :: init() {
 *  RVEMaterialStatus :: init();
 *  calculateCentroid();
 *
 *  unsigned ndim = RVE->giveDimension();
 *  ElementContainer *elems = RVE->giveElements();
 *  Point normal;
 *
 *  if ( active_mechanics ) {
 *      projectors.resize(ndim);
 *
 *      unsigned strain_size = ndim * ndim;
 *      strain_size += ( ndim == 3 ) ? ndim * ndim : ndim; //in 2D only vector
 *
 *      RigidBodyContact *e;
 *      Point xc;
 *      Point alphaVec;
 *      Vector PQ(strain_size);
 *      for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
 *          e = dynamic_cast< RigidBodyContact * >( elems->giveElement(i) );
 *          if ( e ) {
 *              xc =  e->giveIPLoc(0) - centroid; //here we make corrections to have origin of the reference system at the centroid
 *              normal = e->giveNormal();
 *              for ( unsigned v = 0; v < ndim; v++ ) {
 *                  if ( v == 0 ) {
 *                      alphaVec = normal;
 *                  } else if ( v == 1 ) {
 *                      alphaVec = e->giveT1();
 *                  } else if ( v == 2 )                                  {
 *                      alphaVec = e->giveT2();
 *                  }
 *
 *                  if ( ndim == 2 ) {
 *                      PQ [ 0 ] = normal.getX() * alphaVec.getX();
 *                      PQ [ 1 ] = normal.getY() * alphaVec.getY();
 *                      PQ [ 2 ] = normal.getX() * alphaVec.getY();
 *                      PQ [ 3 ] = normal.getY() * alphaVec.getX();
 *                      double factor2D = xc.getX() * alphaVec.getY() - xc.getY() * alphaVec.getX();
 *                      PQ [ 4 ] = factor2D * normal.getX();
 *                      PQ [ 5 ] = factor2D * normal.getY();
 *                  } else   {
 *                      PQ [ 0 ] = normal.getX() * alphaVec.getX();
 *                      PQ [ 1 ] = normal.getY() * alphaVec.getY();
 *                      PQ [ 2 ] = normal.getZ() * alphaVec.getZ();
 *                      PQ [ 3 ] = normal.getZ() * alphaVec.getY();
 *                      PQ [ 4 ] = normal.getY() * alphaVec.getZ();
 *                      PQ [ 5 ] = normal.getZ() * alphaVec.getX();
 *                      PQ [ 6 ] = normal.getX() * alphaVec.getZ();
 *                      PQ [ 7 ] = normal.getY() * alphaVec.getX();
 *                      PQ [ 8 ] = normal.getX() * alphaVec.getY();
 *                      Point factor3D = cross(xc, alphaVec);
 *                      PQ [ 9 ] = normal.getX() * factor3D.getX();
 *                      PQ [ 10 ] = normal.getY() * factor3D.getY();
 *                      PQ [ 11 ] = normal.getZ() * factor3D.getZ();
 *                      PQ [ 12 ] = normal.getY() * factor3D.getZ();
 *                      PQ [ 13 ] = normal.getZ() * factor3D.getY();
 *                      PQ [ 14 ] = normal.getX() * factor3D.getZ();
 *                      PQ [ 15 ] = normal.getZ() * factor3D.getX();
 *                      PQ [ 16 ] = normal.getX() * factor3D.getY();
 *                      PQ [ 17 ] = normal.getY() * factor3D.getZ();
 *                  }
 *                  projectors [ v ].push_back(PQ);
 *              }
 *          }
 *      }
 *  }
 * }
 *
 * //////////////////////////////////////////////////////////
 * //////////////////////////////////////////////////////////
 * // DISCRETE RVE MATERIAL
 * MaterialStatus *DiscreteRVEMaterial :: giveNewMaterialStatus(Element *e) {
 *  DiscreteRVEMaterialStatus *newstat = new DiscreteRVEMaterialStatus(this, e, inputfile);
 *  return newstat;
 * }
 *
 *
 */
