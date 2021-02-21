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

    //generateVolumetricAverageBC();
    generateRandomFixedBC();    
    RVE->init();

    stringstream appendname;
    appendname << "_" << std :: setfill('0') << std :: setw(4) << element->giveID() << "_" << std :: setw(2) << idx;
    RVE->giveExporters()->appendToAllNames(appendname.str() );
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
// general RVE material
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
// DISCRETE RVE MATERIAL
DiscreteRVEMaterialStatus :: DiscreteRVEMaterialStatus(DiscreteRVEMaterial *m, Element *e, fs :: path masterfile) : RVEMaterialStatus(m, e, masterfile) {
    name = "transport RVE mat. status";
}

/////////////////////////////////./////////////////////////
Vector DiscreteRVEMaterialStatus :: giveStress(const Vector &strain) {

    double macro_pressure = 0;
    if ( active_transport ) {
        temp_strain.resize(strain.size()-1);
        for (unsigned i = 0; i<temp_strain.size(); i++) temp_strain[i] = strain[i];
        macro_pressure = strain[strain.size()-1];      
    } else {
        temp_strain = strain;
    }

    unsigned ndim = RVE->giveDimension();
    
    unsigned stra_size = ndim * ndim;
    stra_size += ( ndim == 3 ) ? ndim * ndim : ndim; //in 2D only vector
    if ( !active_mechanics ) {
        stra_size = 0;
    }

    unsigned pres_size = ndim;
    if ( !active_transport ) {
        pres_size = 0;
    }

    //set eigenstrains
    ElementContainer *elems = RVE->giveElements();
    Point normal;

    if ( active_mechanics ) {
        RigidBodyContact *e;
        Vector eigstr(ndim);
        unsigned num = 0;
        for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
            e = dynamic_cast< RigidBodyContact * >( elems->giveElement(i) );
            if ( !e ) {
                continue;
            }
            eigstr *= 0.;
            for ( unsigned v = 0; v < ndim; v++ ) {
                for ( unsigned r = 0; r < stra_size; r++ ) {
                    eigstr [ v ] -= temp_strain [ r ] * mechProjectors [ v ] [ num ] [ r ];
                }
            }
            num++;
            e->giveMatStatus(0)->setEigenStrain(eigstr);
        }
    }
  
    vector < double >  initial_permeability, initial_a;
    if ( active_transport ) {
        Transp1D *e;
        TrsprtMaterialStatus *s;
        TrsprtMaterial *m;
        Vector eigstr(1);
        for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
            e = dynamic_cast< Transp1D * >( elems->giveElement(i) );
            if ( !e ) {
                continue;
            }
            s = static_cast< TrsprtMaterialStatus* >( e->giveMatStatus(0) );
            m = static_cast< TrsprtMaterial* >( s->giveMaterial() );   
            initial_permeability.push_back( m->givePermeability());
            initial_a.push_back( m->giveParamA());
            m->setPermeability( s->calculatePressureDependentPermeability(macro_pressure) ); //calculating pressure depedent conductivity
            m->setParamA(-1.); //switch of linearity    
            normal = e->giveNormal();            
            eigstr *= 0.;
            for ( unsigned v = 0; v < ndim; v++ ) {
                eigstr [ 0 ] -= temp_strain [ stra_size + v ] * normal.giveCoord(v);
            }
            s->setEigenStrain(eigstr);
        }
    }

    //solve
    RVE->resetTime();
    //RVE->giveSolver()->updateSteadyStateMatrix(); //needed when using linear solver - NO, internal forces automatically takes care, it can be switched on optionally
    RVE->giveSolver()->runBeforeEachStep();
    RVE->giveSolver()->solve();

    //REMOVE MACROPRESSURE

    //collect results
    temp_stress.resize(stra_size + pres_size);
    temp_stress *= 0;


    if ( active_mechanics ) {
        double volume = 0.;
        RigidBodyContact *e;
        Vector factor;
        unsigned num = 0;
        for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
            e = dynamic_cast< RigidBodyContact * >( elems->giveElement(i) );
            if ( e ) {
                factor  = e->giveArea() * e->giveLength() * e->giveMatStatus(0)->giveTempStress();
                for ( unsigned v = 0; v < ndim; v++ ) {
                    for ( unsigned r = 0; r < stra_size; r++ ) {
                        temp_stress [ r ] += factor [ v ] * mechProjectors [ v ] [ num ] [ r ];
                    }
                }
                num++;
                volume += e->giveVolume();
            }
        }
        for ( unsigned v = 0; v < stra_size; v++ ) {
            temp_stress [ v ] /= volume;
        }
    }


    if ( active_transport ) {   
        unsigned k = 0;
        double volume = 0.;
        double factor;
        Transp1D *e;
        TrsprtMaterial *m;
        for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
            e = dynamic_cast< Transp1D * >( elems->giveElement(i) );
            if ( e ) {
                normal = e->giveNormal();
                factor = e->giveArea() * e->giveLength() * e->giveMatStatus(0)->giveTempStress() [ 0 ];
                for ( unsigned v = 0; v < ndim; v++ ) {
                    temp_stress [ stra_size + v ] += factor * normal.giveCoord(v); 
                }
                volume += e->giveVolume();

                m = static_cast< TrsprtMaterial* >( e->giveMatStatus(0)->giveMaterial() );
                m->setPermeability( initial_permeability[k] );
                m->setParamA( initial_a[k] ); //switch of nonlinearity   
                k ++;
            }
        }
        for ( unsigned v = 0; v < ndim; v++ ) {
            temp_stress [ stra_size + v ] /= volume;
        }
    }

    /*
    cout << "STRAIN ";
    for(unsigned v=0; v<temp_strain.size(); v++) cout << " " << temp_strain[v];
    cout << endl;
    
    cout << "STRESS ";
    for(unsigned v=0; v<temp_strain.size(); v++) cout << " " << temp_stress[v];
    cout << endl;
    */

    return temp_stress;
}

/////////////////////////////////./////////////////////////
Vector DiscreteRVEMaterialStatus :: giveStressWithFrozenIntVars(const Vector &strain) {
    //TODO: this is WRONG, needs fix in future
    return giveStress(strain);
}

//////////////////////////////////////////////////////////
Matrix DiscreteRVEMaterialStatus :: giveStiffnessTensor(string type, unsigned ndim) const {
    unsigned stra_size = ndim * ndim;
    stra_size += ( ndim == 3 ) ? ndim * ndim : ndim; //in 2D only vector
    if ( !active_mechanics ) {
        stra_size = 0;
    }
    Matrix Keff(stra_size, stra_size);

    unsigned pres_size = ndim;
    if ( !active_transport ) {
        pres_size = 0;
    }
    Matrix Leff(pres_size, pres_size);


    ElementContainer *elems = RVE->giveElements();
    double volume = 0;
    Point normal;
    Vector n(ndim);

    if ( active_mechanics ) {
        RigidBodyContact *e;
        volume = 0;
        Matrix stiff;
        unsigned num = 0;
        for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
            e = dynamic_cast< RigidBodyContact * >( elems->giveElement(i) );
            if ( e ) {
                stiff  = e->giveMatStatus(0)->giveStiffnessTensor(type, ndim);
                for ( unsigned v = 0; v < ndim; v++ ) {
                    Keff += dyadicProduct(mechProjectors [ v ] [ num ], mechProjectors [ v ] [ num ]) * e->giveLength() * e->giveArea() * stiff [ v ] [ v ];
                }
                num++;
                volume += e->giveVolume();
            }
        }
        Keff /= volume;
    }


    if ( active_transport ) {
        Transp1D *e;
        volume = 0;
        for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
            e = dynamic_cast< Transp1D * >( elems->giveElement(i) );
            if ( e ) {
                normal = e->giveNormal();
                n [ 0 ] = normal.getX();
                n [ 1 ] = normal.getY();
                if ( ndim == 3 ) {
                    n [ 2 ] = normal.getZ();
                }
                Leff += dyadicProduct(n, n) * ( e->giveLength() * e->giveArea() * e->giveMatStatus(0)->giveStiffnessTensor(type, ndim) [ 0 ] [ 0 ] );
                volume += e->giveVolume();
            }
        }
        Leff /= volume;
    }

    //Keff.print();

    if ( !active_transport ) {
        return Keff;
    }
    if ( !active_mechanics ) {
        return Leff;
    }

    Matrix KL(stra_size + pres_size, stra_size + pres_size);
    for ( unsigned r = 0; r < stra_size; r++ ) {
        for ( unsigned s = 0; s < stra_size; s++ ) {
            KL [ r ] [ s ] = Keff [ r ] [ s ];
        }
    }

    for ( unsigned r = 0; r < pres_size; r++ ) {
        for ( unsigned s = 0; s < pres_size; s++ ) {
            KL [ stra_size + r ] [ stra_size + s ] = Leff [ r ] [ s ];
        }
    }

    return KL;
}

//////////////////////////////////////////////////////////
double DiscreteRVEMaterialStatus :: giveDampingConstant() const {
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
void DiscreteRVEMaterialStatus :: generateRandomFixedBC() {
    BCContainer *bconds = RVE->giveBC();
    FunctionContainer *funcs = RVE->giveFunctions();
    ConstraintContainer *constrs = RVE->giveConstraints();

    //mechanics
    JointDoF *jd;
    Node * masternode;    
    active_mechanics = false;
    for ( unsigned j = 0; j < constrs->giveSize(); j++ ) {
        if (active_mechanics) break; 
        jd = constrs->giveConstraint(j);
        for (unsigned k=0; k<jd->giveNumOfMasters(); k++){
            masternode = jd->giveMasterNode(0);
            if ( masternode->doesMechanics() && ( dynamic_cast< MechDoF * >( masternode ) == nullptr ) ) {      
                BoundaryCondition *bc;
                vector< int >dBC, nBC;
                dBC.resize(masternode->giveNumberOfDoFs(), funcs->giveSize() ); //todo: warning C4267: 'argument': conversion from 'size_t' to 'const _Ty', possible loss of data
                nBC.resize(masternode->giveNumberOfDoFs(), -1);
                bc = new BoundaryCondition(masternode, dBC, nBC);
                bconds->addBoundaryCondition(bc);                

                active_mechanics = true;              
                break;
            }  
        }           
    }

    //transport
    active_transport = false;
    for ( unsigned j = 0; j < constrs->giveSize(); j++ ) {
        if (active_transport) break;
        jd = constrs->giveConstraint(j);
        for (unsigned k=0; k<jd->giveNumOfMasters(); k++){
            masternode = jd->giveMasterNode(0);
            if ( masternode->doesTransport() && ( dynamic_cast< TrsDoF * >( masternode ) == nullptr ) ) {   
                BoundaryCondition *bc;
                vector< int >dBC, nBC;
                dBC.resize(masternode->giveNumberOfDoFs(), funcs->giveSize() ); //todo: warning C4267: 'argument': conversion from 'size_t' to 'const _Ty', possible loss of data
                nBC.resize(masternode->giveNumberOfDoFs(), -1);
                bc = new BoundaryCondition(masternode, dBC, nBC);
                bconds->addBoundaryCondition(bc);                

                active_transport = true;              
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
void DiscreteRVEMaterialStatus :: generateVolumetricAverageBC() {
    //this function applies volumetric constraint
    //it is not really wise to use it as it leads to relatively full stiffness matrix after transfer to constraint space
    NodeContainer *nodes = RVE->giveNodes();
    BCContainer *bconds = RVE->giveBC();
    FunctionContainer *funcs = RVE->giveFunctions();
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

    //mechanics
    for ( unsigned n = 0; n < nodes->giveSize(); n++ ) {
        if ( nodes->giveNode(n)->doesMechanics() && ( dynamic_cast< MechDoF * >( nodes->giveNode(n) ) == nullptr ) ) {
            vm.push_back(nodes->giveNode(n) );
        }
    }
    if ( vm.size() > 0 ) {
        unsigned nDoFs = 3;
        if ( ndim == 3 ) {
            nDoFs = 6;
        }
        MechDoF *pn = new MechDoF(nDoFs);
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

        active_mechanics = true;
    } else   {
        active_mechanics = false;
    }

    //transport
    vm.clear();
    for ( unsigned n = 0; n < nodes->giveSize(); n++ ) {
        if ( nodes->giveNode(n)->doesTransport() && ( dynamic_cast< TrsDoF * >( nodes->giveNode(n) ) == nullptr ) ) {
            vm.push_back(nodes->giveNode(n) );
        }
    }
    if ( vm.size() > 0 ) {
        TrsDoF *tn = new TrsDoF(1);
        nodes->addNode(tn);

        dirs.resize(vm.size() );
        va = new VolumetricAverage(vm, dirs, tn, 0, elems, constrs);
        constrs->addConstraint(va);

        BoundaryCondition *bc;
        vector< int >dBC, nBC;
        dBC.resize(1, funsize);
        nBC.resize(1, -1);
        bc = new BoundaryCondition(tn, dBC, nBC);
        bconds->addBoundaryCondition(bc);

        active_transport = true;
    } else   {
        active_transport = false;
    }
}

//////////////////////////////////////////////////////////
void DiscreteRVEMaterialStatus :: calculateCentroid() {
    centroid = Point(0., 0., 0.);

    if ( !active_mechanics ) {
        return;
    }

    //find centroid
    ElementContainer *elems = RVE->giveElements();
    unsigned ndim = RVE->giveDimension();

    //mechanics
    RigidBodyContact *e;
    DisMechMaterialStatus *ms;
    double weight;
    double totalweight = 0;
    for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
        e = dynamic_cast< RigidBodyContact * >( elems->giveElement(i) );
        if ( e ) {
            ms = static_cast< DisMechMaterialStatus * >( e->giveMatStatus(0) );
            weight = e->giveLength() * e->giveArea() * ms->giveDensity();
            ;
            totalweight += weight;
            centroid += e->giveIPLoc(0) * weight;
        }
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
void DiscreteRVEMaterialStatus :: init() {
    RVEMaterialStatus :: init();
    calculateCentroid();

    unsigned ndim = RVE->giveDimension();
    ElementContainer *elems = RVE->giveElements();
    Point normal;

    if ( active_mechanics ) {
        mechProjectors.resize(ndim);

        unsigned stra_size = ndim * ndim;
        stra_size += ( ndim == 3 ) ? ndim * ndim : ndim; //in 2D only vector

        RigidBodyContact *e;
        Point xc;
        Point alphaVec;
        Vector PQ(stra_size);
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
                    } else if ( v == 2 )                                  {
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
                    } else   {
                        PQ [ 0 ] = normal.getX() * alphaVec.getX();
                        PQ [ 1 ] = normal.getY() * alphaVec.getY();
                        PQ [ 2 ] = normal.getZ() * alphaVec.getZ();
                        PQ [ 3 ] = normal.getY() * alphaVec.getZ();
                        PQ [ 4 ] = normal.getZ() * alphaVec.getY();
                        PQ [ 5 ] = normal.getX() * alphaVec.getZ();
                        PQ [ 6 ] = normal.getZ() * alphaVec.getX();
                        PQ [ 7 ] = normal.getX() * alphaVec.getY();
                        PQ [ 8 ] = normal.getY() * alphaVec.getZ();
                        Point factor3D = cross(xc, alphaVec);
                        PQ [ 9 ] = normal.getX() * factor3D.getX();
                        PQ [ 10 ] = normal.getY() * factor3D.getY();
                        PQ [ 11 ] = normal.getZ() * factor3D.getZ();
                        PQ [ 12 ] = normal.getY() * factor3D.getZ();
                        PQ [ 13 ] = normal.getZ() * factor3D.getY();
                        PQ [ 14 ] = normal.getX() * factor3D.getZ();
                        PQ [ 15 ] = normal.getZ() * factor3D.getX();
                        PQ [ 16 ] = normal.getX() * factor3D.getY();
                        PQ [ 17 ] = normal.getY() * factor3D.getZ();
                    }
                    mechProjectors [ v ].push_back(PQ);
                }
            }
        }
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE RVE MATERIAL
MaterialStatus *DiscreteRVEMaterial :: giveNewMaterialStatus(Element *e) {
    DiscreteRVEMaterialStatus *newstat = new DiscreteRVEMaterialStatus(this, e, inputfile);
    return newstat;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE RVE MATERIAL STATUS PRE-COMPUTED
DiscreteRVEMaterialPrecomputedStatus :: DiscreteRVEMaterialPrecomputedStatus(DiscreteRVEMaterialPrecomputed *m, Element *e, fs :: path masterfile) : DiscreteRVEMaterialStatus(m, e, masterfile) {
    name = "transport RVE precomputed mat. status";
    is_master_status = false;
    temp_nonlin = 1;
}

/////////////////////////////////./////////////////////////
Vector DiscreteRVEMaterialPrecomputedStatus :: giveStress(const Vector &strain) {

    double macro_pressure = 0;
    if ( active_transport ) {
        temp_strain.resize(strain.size()-1);
        for (unsigned i = 0; i<temp_strain.size(); i++) temp_strain[i] = strain[i];
        macro_pressure = strain[strain.size()-1];      
    } else {
        temp_strain = strain;
    }

    DiscreteRVEMaterialPrecomputed* macromaterial = static_cast < DiscreteRVEMaterialPrecomputed* >(mat);
    temp_nonlin = macromaterial->giveMasterStatus()->calculatePressureDependentPermeability(macro_pressure)/macromaterial->giveMasterMaterial()->givePermeability();
    temp_stress = matrix_vector_multiply(macromaterial->givePrecomputedConductivity(),temp_strain)*temp_nonlin;

    return temp_stress;
}

/////////////////////////////////./////////////////////////
Vector DiscreteRVEMaterialPrecomputedStatus :: giveStressWithFrozenIntVars(const Vector &strain) {
    return giveStress(strain);
}

//////////////////////////////////////////////////////////
Matrix DiscreteRVEMaterialPrecomputedStatus :: giveStiffnessTensor(string type, unsigned ndim) const {
    DiscreteRVEMaterialPrecomputed* macromaterial = static_cast < DiscreteRVEMaterialPrecomputed* >(mat);
    return macromaterial->givePrecomputedConductivity() * temp_nonlin;

}

//////////////////////////////////////////////////////////
double DiscreteRVEMaterialPrecomputedStatus :: giveDampingConstant() const {
    DiscreteRVEMaterialPrecomputed* macromaterial = static_cast < DiscreteRVEMaterialPrecomputed* >(mat);
    return macromaterial->givePrecomputedCapacity();
};

//////////////////////////////////////////////////////////
void DiscreteRVEMaterialPrecomputedStatus :: init() {

    DiscreteRVEMaterialPrecomputed* macromaterial = static_cast < DiscreteRVEMaterialPrecomputed* >(mat);    
    Matrix stiff = macromaterial -> givePrecomputedConductivity();
    if (stiff.size()==0){

        is_master_status = true;
        DiscreteRVEMaterialStatus :: init();

        //TODO: if this is needed, one should call it only once and distribute the result among all identical RVEs
        cout << "Precomputing pressure fields on RVE" << endl;

        unsigned ndim = RVE->giveDimension();
        stiff = Matrix(ndim,ndim);
        Vector strainDoFs(ndim+1);
        Vector stress;
        for(unsigned i=0; i<ndim; i++){
            strainDoFs[i] = 1.;
            stress = DiscreteRVEMaterialStatus :: giveStress(strainDoFs);
            strainDoFs[i] = 0.;   
            for(unsigned j=0; j<ndim; j++) stiff[i][j] = stress[j];
        }
        double c = DiscreteRVEMaterialStatus :: giveDampingConstant();

        TrsprtMaterialStatus* status = static_cast < TrsprtMaterialStatus* >(RVE->giveElements()->giveElement(0)->giveMatStatus(0));
        TrsprtMaterial* material = static_cast < TrsprtMaterial* >(RVE->giveElements()->giveElement(0)->giveMaterial());
        macromaterial->setPrecomputedConductivityAndCapacityAndMasterMaterial(stiff,c, status, material);
    }   

    else{
        active_transport = true;
    }
}

//////////////////////////////////////////////////////////
void DiscreteRVEMaterialPrecomputed :: setPrecomputedConductivityAndCapacityAndMasterMaterial(Matrix lam, double c, TrsprtMaterialStatus* masterS,  TrsprtMaterial* masterM){
    conductivity=lam; 
    capacity = c;
    masterStatus = masterS;
    masterMaterial = masterM;

    char* buffer = "precomputed_conductivity.out";   
    ofstream outputfile( ( masterModel->giveResultDirectory() / buffer ).string() );
    if ( outputfile.is_open() ) {
        outputfile << std :: scientific;
        outputfile.precision(10);

        unsigned size = conductivity.numRows();        
        for ( unsigned i = 0; i < size; i++ ) {
            for ( unsigned j = 0; j < size; j++ ) {
                if (j>0) outputfile << "\t";
                outputfile << conductivity[i][j];
            }
            outputfile << endl;
        }
        outputfile.close();
    }
}    


//////////////////////////////////////////////////////////
void DiscreteRVEMaterialPrecomputedStatus :: update(){
    if (is_master_status) RVEMaterialStatus :: update();
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE RVE MATERIAL
MaterialStatus *DiscreteRVEMaterialPrecomputed :: giveNewMaterialStatus(Element *e) {
    DiscreteRVEMaterialPrecomputedStatus *newstat = new DiscreteRVEMaterialPrecomputedStatus(this, e, inputfile);
    return newstat;
}
