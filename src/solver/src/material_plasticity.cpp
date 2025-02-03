#include "material_plasticity.h"


using namespace std;


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// VON MISES PLASTICITY
//////////////////////////////////////////////////////////
VonMisesPlasticMaterialStatus :: VonMisesPlasticMaterialStatus(VonMisesPlasticMaterial *m, Element *e, unsigned ipnum) : TensMechMaterialStatus(m, e, ipnum) {
    name = "tensorial mechanical von Mises plastic mat. status";
    unsigned strainsize = m->giveStrainSize();
    temp_plasticstrain = Vector :: Zero(strainsize);
    plasticstrain = Vector :: Zero(strainsize);
    temp_backstress = Vector :: Zero(strainsize);
    backstress = Vector :: Zero(strainsize);
    N = Vector :: Zero(6);
    sigmay = m->giveSigma0();
    temp_sigmay = sigmay;
    temp_outplane_plasticstrain = 0.;
    outplane_plasticstrain = 0.;
    temp_outplane_backstress = 0.;
    outplane_backstress = 0.;
    beta_t = 1.;
    f = -1;

}

//////////////////////////////////////////////////////////
Vector VonMisesPlasticMaterialStatus :: giveStress(const Vector &strain, double timeStep) {
    temp_strain = addEigenStrain(strain);

    VonMisesPlasticMaterial *vmpm = static_cast< VonMisesPlasticMaterial * >( mat );

    Matrix D;
    Vector temp_strain_full = Vector :: Zero(6);
    Vector updt_strain_full = Vector :: Zero(6);
    Vector plasticstrain_full = Vector :: Zero(6);
    Vector updt_stress_full = Vector :: Zero(6);
    Vector backstress_full = Vector :: Zero(6);

    unsigned dimension = vmpm->giveDimension();
    if ( dimension == 2 ) {
        if ( vmpm->isPlaneStress() ) { //Plane Stress - out of plane stress iteration necessary to implement
            // f = sqrt( pow(sigmaEff[0],2) - sigmaEff[0]*sigmaEff[1] + pow(sigmaEff[1],2) + 3*pow(sigmaEff[2],2))  - temp_sigmay;
            cerr << name << " error: " <<  "Plane Stress not implemented" << endl;
            exit(1);
        } else {   //Plane Strain
            D = giveElasticStiffnessTensor3D();

            temp_strain_full [ 0 ] = temp_strain [ 0 ];
            temp_strain_full [ 1 ] = temp_strain [ 1 ];
            temp_strain_full [ 5 ] = temp_strain [ 2 ];

            updt_strain_full [ 0 ] = updt_strain [ 0 ];
            updt_strain_full [ 1 ] = updt_strain [ 1 ];
            updt_strain_full [ 5 ] = updt_strain [ 2 ];

            plasticstrain_full [ 0 ] = plasticstrain [ 0 ];
            plasticstrain_full [ 1 ] = plasticstrain [ 1 ];
            plasticstrain_full [ 5 ] = plasticstrain [ 2 ];
            plasticstrain_full [ 2 ] = outplane_plasticstrain;

            updt_stress_full = D * ( updt_strain_full - plasticstrain_full );

            backstress_full [ 0 ] = backstress [ 0 ];
            backstress_full [ 1 ] = backstress [ 1 ];
            backstress_full [ 5 ] = backstress [ 2 ];
            backstress_full [ 2 ] = outplane_backstress;
        }
    }  else if ( dimension == 3 ) {
        D = giveStiffnessTensor("elastic");

        temp_strain_full = temp_strain;
        updt_strain_full = updt_strain;
        plasticstrain_full = plasticstrain;
        updt_stress_full = updt_stress;
        backstress_full = backstress;
    } else {
        std :: cout << "Dimension not supported for Mises plasticity" << "\n";
    }

    Vector sigmaTrial = updt_stress_full + D * ( temp_strain_full - updt_strain_full );
    Vector sigmaEff = sigmaTrial - backstress_full;
    Vector sigmaIso = Vector :: Zero(6);
    sigmaIso [ 0 ] = sigmaIso [ 1 ] = sigmaIso [ 2 ] = 1. / 3. * ( sigmaEff [ 0 ] + sigmaEff [ 1 ] + sigmaEff [ 2 ] );
    Vector n = sigmaEff - sigmaIso; // Deviatoric part

    f = sqrt(3. * ( 1. / 6. * ( pow( ( sigmaEff [ 0 ] - sigmaEff [ 1 ] ), 2) + pow( ( sigmaEff [ 0 ] - sigmaEff [ 2 ] ), 2) + pow( ( sigmaEff [ 1 ] - sigmaEff [ 2 ] ), 2) ) + pow(sigmaEff [ 3 ], 2) + pow(sigmaEff [ 4 ], 2) + pow(sigmaEff [ 5 ], 2) ) ) - temp_sigmay;

    if ( f <= 0.000001 ) {  // Elastic regime
        if ( dimension == 2 ) {
            temp_stress << sigmaTrial [ 0 ], sigmaTrial [ 1 ], sigmaTrial [ 5 ];
        } else {
            temp_stress = sigmaTrial;
        }
        temp_plasticstrain = plasticstrain;
        temp_backstress = backstress;
        temp_sigmay = sigmay;
        temp_outplane_plasticstrain = outplane_plasticstrain;
        temp_outplane_backstress = outplane_backstress;

    } else { // Plastic regime
        double n_norm = sqrt(n [ 0 ] * n [ 0 ] + n [ 1 ] * n [ 1 ] + n [ 2 ] * n [ 2 ] + 2 * n [ 3 ] * n [ 3 ] + 2 * n [ 4 ] * n [ 4 ] + 2 * n [ 5 ] * n [ 5 ]);
        N = n / n_norm;
        double G = vmpm->giveElasticModulus() / ( 2. * ( 1. + vmpm->givePoissonsRatio() ) );
        double lambda = 1. / ( 2. * G ) * 1. / ( 1. + ( vmpm->giveHardeningModulus() / ( 3. * G ) ) ) * ( n_norm - sigmay * sqrt(2. / 3.) );

        Vector N2 = N;
        N2 [ 3 ] = 2 * N [ 3 ];
        N2 [ 4 ] = 2 * N [ 4 ];
        N2 [ 5 ] = 2 * N [ 5 ];

        Vector temp_stress_full = Vector :: Zero(6);
        Vector temp_backstress_full = Vector :: Zero(6);
        Vector temp_plasticstrain_full = Vector :: Zero(6);

        temp_plasticstrain_full = plasticstrain_full + lambda * N2;   // N2 namisto N @ P viz Python

        temp_backstress_full = backstress_full + 2. / 3. * ( 1 - vmpm->giveBetaRatio() ) * vmpm->giveHardeningModulus() * lambda * N;
        temp_sigmay = sigmay + 2. / 3. * vmpm->giveBetaRatio() * vmpm->giveHardeningModulus() * lambda * sqrt(3. / 2.);
        temp_stress_full = sigmaTrial - 2. * G * lambda * N;

        beta_t = sqrt(2./3.) * (temp_sigmay + sqrt(2./3.) * lambda * ( 1 - vmpm->giveBetaRatio() ) * vmpm->giveHardeningModulus()) / n_norm;

        if ( dimension == 2 ) {
            temp_stress << temp_stress_full [ 0 ], temp_stress_full [ 1 ], temp_stress_full [ 5 ];

            temp_plasticstrain [ 0 ] = temp_plasticstrain_full [ 0 ];
            temp_plasticstrain [ 1 ] = temp_plasticstrain_full [ 1 ];
            temp_plasticstrain [ 2 ] = temp_plasticstrain_full [ 5 ];

            temp_backstress [ 0 ] = temp_backstress_full [ 0 ];
            temp_backstress [ 1 ] = temp_backstress_full [ 1 ];
            temp_backstress [ 2 ] = temp_backstress_full [ 5 ];

            temp_outplane_plasticstrain = temp_plasticstrain_full [ 2 ];
            temp_outplane_backstress = temp_backstress_full [ 2 ];
        } else {
            temp_stress = temp_stress_full;
            temp_plasticstrain = temp_plasticstrain_full;
            temp_backstress = temp_backstress_full;
        }
        // cout << "PLASTIC" << "\n";

    }
    return temp_stress;
};

//////////////////////////////////////////////////////////
Vector VonMisesPlasticMaterialStatus :: giveStressWithFrozenIntVars(const Vector &strain, double timeStep) {
    ( void ) timeStep;
    temp_strain = addEigenStrain(strain);

    VonMisesPlasticMaterial *vmpm = static_cast< VonMisesPlasticMaterial * >( mat );

    unsigned dimension = vmpm->giveDimension();
    if ( dimension == 2 ) {
        Vector temp_strain_full = Vector :: Zero(6);
        Vector temp_plasticstrain_full = Vector :: Zero(6);
        Vector temp_stress_full = Vector :: Zero(6);

        temp_strain_full [ 0 ] = temp_strain [ 0 ];
        temp_strain_full [ 1 ] = temp_strain [ 1 ];
        temp_strain_full [ 5 ] = temp_strain [ 2 ];

        temp_plasticstrain_full [ 0 ] = temp_plasticstrain [ 0 ];
        temp_plasticstrain_full [ 1 ] = temp_plasticstrain [ 1 ];
        temp_plasticstrain_full [ 5 ] = temp_plasticstrain [ 2 ];
        temp_plasticstrain_full [ 2 ] = temp_outplane_plasticstrain;

        temp_stress_full = giveElasticStiffnessTensor3D() * ( temp_strain_full - temp_plasticstrain_full );

        temp_stress << temp_stress_full [ 0 ], temp_stress_full [ 1 ], temp_stress_full [ 5 ];
    } else {
        temp_stress = giveStiffnessTensor("elastic") * ( temp_strain - temp_plasticstrain );
    }

    return temp_stress;
};

//////////////////////////////////////////////////////////
Matrix VonMisesPlasticMaterialStatus :: giveStiffnessTensor(string type) const {
    if ( type == "tangent" || type == "consistent") {
        ( void ) type; 
        VonMisesPlasticMaterial *vmpm = static_cast< VonMisesPlasticMaterial * >( mat );
        Matrix elastic_tensor = TensMechMaterialStatus :: giveStiffnessTensor(type);
        Matrix elastoplastic_tensor = elastic_tensor;
        Vector N_dim = Vector :: Zero(vmpm->giveStrainSize());
        Matrix I_min_13 = Matrix :: Zero(vmpm->giveStrainSize(), vmpm->giveStrainSize());

        if ( f <= 0.000001 ) {  // Elastic regime
            elastoplastic_tensor = elastic_tensor;

        } else { // Plastic regime

            double G = vmpm->giveElasticModulus() / ( 2. * ( 1. + vmpm->givePoissonsRatio() ) );
            
            double gamma = 1. / (1. + vmpm->giveHardeningModulus() / (3. * G));

            unsigned dimension = vmpm->giveDimension();
            if ( dimension == 2 ) {
                I_min_13(0,0) = I_min_13(1,1) = 2./3.;
                I_min_13(0,1) = I_min_13(1,0) = -1./3.;
                I_min_13(2,2) = 0.5;

                N_dim[0] = N[0];
                N_dim[1] = N[1];
                N_dim[2] = N[5];

            } else {
                I_min_13(0,0) = I_min_13(1,1) = I_min_13(2,2) = 2./3.;
                I_min_13(0,1) = I_min_13(0,2) = I_min_13(1,0) = I_min_13(1,2) = I_min_13(2,0) = I_min_13(2,1) = -1./3.;
                I_min_13(3,3) = I_min_13(4,4) = I_min_13(5,5) = 0.5;
                N_dim = N;
            }

            if ( type == "tangent") {
                elastoplastic_tensor = elastic_tensor - 2.*G*gamma* N_dim * N_dim.transpose(); // continuum elastoplastic tensor
            } else {
                elastoplastic_tensor = elastic_tensor - (1. - beta_t) * 2.*G* I_min_13  - 2.*G* (gamma - (1. - beta_t))* N_dim * N_dim.transpose(); // consistent elastoplastic tensor 
            }

        }

        return elastoplastic_tensor;  
        
    } else {
        ( void ) type;
        return TensMechMaterialStatus :: giveStiffnessTensor(type);
    }

};


//////////////////////////////////////////////////////////
void VonMisesPlasticMaterialStatus :: update() {
    TensMechMaterialStatus :: update();
    backstress = temp_backstress;
    plasticstrain = temp_plasticstrain;
    sigmay = temp_sigmay;
    outplane_plasticstrain = temp_outplane_plasticstrain;
    outplane_backstress = temp_outplane_backstress;
}

//////////////////////////////////////////////////////////
bool VonMisesPlasticMaterialStatus :: giveValues(string code, Vector &result) const {
    if ( code.compare("plastic_strain") == 0 || code.compare("plasticstrain") == 0 ) {
        unsigned size = temp_plasticstrain.size();
        result.resize(size);
        for ( unsigned p = 0; p < size; p++ ) {
            result [ p ] = temp_plasticstrain [ p ];
        }
        return true;
    } else if ( code.compare("backstress") == 0 || code.compare("back_stress") == 0 ) {
        unsigned size = temp_backstress.size();
        result.resize(size);
        for ( unsigned p = 0; p < size; p++ ) {
            result [ p ] = temp_backstress [ p ];
        }
        return true;
    } else if ( code.compare("sigmay") == 0 ) {
        result.resize(1);
        result [ 0 ] = temp_sigmay;
        return true;
    } else if ( code.compare("outplane_backstress") == 0 ) {
        result.resize(1);
        result [ 0 ] = temp_outplane_backstress;
        return true;
    } else {
        return TensMechMaterialStatus :: giveValues(code, result);
    }
}

//////////////////////////////////////////////////////////
VonMisesPlasticMaterial :: VonMisesPlasticMaterial(unsigned dimension) : TensMechMaterial(dimension) {
    name = "von Mises plastic material";
    strainsize = ( dim - 1 ) * 3;
    sigma0 = 0;
    H = 0;
    beta = 0;
};

//////////////////////////////////////////////////////////
void VonMisesPlasticMaterial :: readFromLine(istringstream &iss) {
    TensMechMaterial :: readFromLine(iss); //read elastic parameters

    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream

    string param;
    bool bs0 = false;

    while (  iss >> param ) {
        if ( param.compare("sigma0") == 0 ) {
            bs0 = true;
            iss >> sigma0;
        } else if ( param.compare("H") == 0 ) {
            iss >> H;
        } else if ( param.compare("beta") == 0 ) {
            iss >> beta;
        }
    }
    if ( !bs0 ) {
        cerr << name << ": material parameter 'sigma0' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
};

//////////////////////////////////////////////////////////
MaterialStatus *VonMisesPlasticMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    VonMisesPlasticMaterialStatus *newStatus = new VonMisesPlasticMaterialStatus(this, e, ipnum);
    return newStatus;
};
