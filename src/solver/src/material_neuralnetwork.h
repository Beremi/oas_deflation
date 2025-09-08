#ifndef _MATERIAL_NEURAL_H
#define _MATERIAL_NEURAL_H

#include "linalg.h"
#include "material.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <torch/script.h>

struct Layer {
    std :: string name = "Dense";
    int num_layers = 0;
    int hidden_size = 0;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// RECURRENT NETWORK TENSORIAL MATERIAL

class NeuralNetworkMaterial;
class NeuralNetworkMaterialStatus : public MaterialStatus
{
protected:
    Matrix giveElasticStiffnessTensor3D() const;

    std :: vector< std :: vector< torch :: Tensor > >hc_vectors; // h & c tensors for LSTM or h tensor for GRU layers
    std :: vector< std :: vector< torch :: Tensor > >temp_hc_vectors; // temp -  h & c tensors for LSTM or h tensor for GRU layers


public:
    NeuralNetworkMaterialStatus(NeuralNetworkMaterial *m, Element *e, unsigned ipnum);
    virtual ~NeuralNetworkMaterialStatus() {};
    virtual Matrix giveStiffnessTensor(std :: string type) const;
    virtual Vector giveStress(const Vector &strain, double timeStep);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    virtual double giveMassConstant() const;
    virtual bool giveValues(std :: string code, Vector &result) const;
    virtual void update();
    virtual Matrix giveMassTensor() const;
    virtual Matrix giveDampingTensor() const;

    // virtual torch :: Tensor giveHiddenState() const;
    virtual std :: vector< std :: vector< torch :: Tensor > >giveHiddenState() const;
};

//////////////////////////////////////////////////////////
class NeuralNetworkMaterial : public Material
{
protected:
    double E, nu, density;
    double E0, ft, Gt, RVEsize;
    bool planeStress;
    bool stiffness_from_matrix;

    fs :: path sm_path;
    fs :: path nm_path;
    Matrix readStiffMatrixFromFile() const;
    Matrix readDataNormalizationMatrix(int size, fs :: path matrix_path) const;
    Matrix stiffmat_elastic;
    Matrix norm;

    fs :: path ml_path;
    torch :: jit :: script :: Module network;
    std :: vector< Layer >layers; // layer info


public:
    NeuralNetworkMaterial(unsigned dimension) : Material(dimension) { name = "Neural network mechanical material"; planeStress = true; strainsize = ( dim - 1 ) * 3; };
    ~NeuralNetworkMaterial() {};
    virtual void readFromLine(std :: istringstream &iss);
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    double giveElasticModulus() const { return E; }
    double givePoissonsRatio() const { return nu; }
    Matrix giveStiffnessMatrix() const { return stiffmat_elastic; };
    double giveDensity() const { return density; };
    double giveE0() const { return E0; }
    double giveft() const { return ft; }
    double giveGt() const { return Gt; }
    double giveRVEsize() const { return RVEsize; }
    bool isPlaneStress() { return planeStress; };
    bool isStiffnessFromMatrix() { return stiffness_from_matrix; };
    virtual void init(MaterialContainer *matcont) { Material :: init(matcont);  };

    Matrix giveNormMatrix() const { return norm; };
    std :: vector< Layer >giveNetworkProps() const { return layers; };

    // std::vector<torch::jit::IValue> predictNetwork(std::vector<torch::jit::IValue> inputs); // I am not doing this because it requires type definition, but "auto" is more versatile
    torch :: jit :: script :: Module giveNetwork() const { return network; };
};



#endif /* _MATERIAL_NEURAL_H */
