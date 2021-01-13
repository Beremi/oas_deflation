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
    generateVolumetricAverageBC();
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
// transport RVE status
DiscreteRVEMaterialStatus :: DiscreteRVEMaterialStatus(DiscreteRVEMaterial *m, Element *e, fs :: path masterfile) : RVEMaterialStatus(m, e, masterfile) {
    name = "transport RVE mat. status";
}

/////////////////////////////////./////////////////////////
Vector DiscreteRVEMaterialStatus :: giveStress(const Vector &strain) {
    temp_strain = strain;

    volumAverFunc->setYValue(0., 0);

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

    if ( active_transport ) {
        Transp1D *e;
        Vector eigstr(0);
        for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
            e = dynamic_cast< Transp1D * >( elems->giveElement(i) );
            if ( !e ) {
                continue;
            }
            normal = e->giveNormal();
            eigstr *= 0.;
            for ( unsigned v = 0; v < ndim; v++ ) {
                eigstr [ 0 ] = -strain [ stra_size + v ] * normal.giveCoord(v);
            }
            e->giveMatStatus(0)->setEigenStrain(eigstr);
        }
    }

    //solve
    RVE->resetTime();
    RVE->giveSolver()->runBeforeEachStep();
    RVE->giveSolver()->solve();

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
        double volume = 0.;
        double factor;
        Transp1D *e;
        for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
            e = dynamic_cast< Transp1D * >( elems->giveElement(i) );
            if ( e ) {
                normal = e->giveNormal();
                factor = e->giveArea() * e->giveLength() * e->giveMatStatus(0)->giveTempStress() [ 0 ];
                for ( unsigned v = 0; v < ndim; v++ ) {
                    temp_stress [ stra_size + v ] += factor * normal.giveCoord(v);
                }
                volume += e->giveVolume();
            }
        }
        for ( unsigned v = 0; v < ndim; v++ ) {
            temp_stress [ stra_size + v ] /= volume;
        }
    }

    /*
     * cout << "STRAIN ";
     * for(unsigned v=0; v<stra_size; v++) cout << " " << temp_strain[v];
     * cout << endl;
     *
     * cout << "STRESS ";
     * for(unsigned v=0; v<stra_size; v++) cout << " " << temp_stress[v];
     * cout << endl;
     */

    return temp_stress;
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
double DiscreteRVEMaterialStatus :: giveMassConstant() const {
    ElementContainer *elems = RVE->giveElements();
    Transp1D *e;
    double volume = 0;
    double mass = 0;
    for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
        e = dynamic_cast< Transp1D * >( elems->giveElement(i) );
        if ( e ) {
            mass += e->giveVolume() * e->giveMatStatus(0)->giveMassConstant();
            volume += e->giveVolume();
        }
    }
    return mass / volume;
};

//////////////////////////////////////////////////////////
void DiscreteRVEMaterialStatus :: generateVolumetricAverageBC() {
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
// dicrete RVE material
MaterialStatus *DiscreteRVEMaterial :: giveNewMaterialStatus(Element *e) {
    DiscreteRVEMaterialStatus *newstat = new DiscreteRVEMaterialStatus(this, e, inputfile);
    return newstat;
}
