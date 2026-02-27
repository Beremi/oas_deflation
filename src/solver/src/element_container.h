#ifndef _ELEMENT_CONTAINER_H
#define _ELEMENT_CONTAINER_H

#include "globals.h"
#include "element.h"
#include "boundary_condition.h"
#include "node_container.h"

class Model; //forward declaration
class ConstraintContainer; //forward declaration
class CrossSectionContainer; //forward declaration

//////////////////////////////////////////////////////////
class ElementContainer
{
private:
    std :: vector< Element * >elems;
    NodeContainer *nodes = nullptr;
    ;
    BCContainer *bconds = nullptr;
    ;
    MaterialContainer *materials = nullptr;
    Model *model = nullptr;
    ConstraintContainer *constcont = nullptr;
    CrossSectionContainer *crosssects = nullptr;
    ;
    unsigned max_sol_order = 0; //maximum number of successive rounds of internal force evaluations
    void prepareStructuralMatrix(CoordinateIndexedSparseMatrix &K, unsigned diffType, bool lumped, bool BC_applied = true) const;
    void updateStructuralMatrix(CoordinateIndexedSparseMatrix &K, unsigned diffType, std :: string matrixType, bool lumped, bool BC_applied = true, bool solver_numbering = true) const;
    void integrateDampingOrInertiaForces(const Vector &full_v, Vector &full_f, unsigned diffType) const;
    std :: vector< std :: string >file_to_load_from;

public:
    ElementContainer() { model = nullptr; };
    ~ElementContainer();
    void setModel(Model *mod);
    void setContainers(ConstraintContainer *c) { constcont = c; };
    void readFromFile(const std :: string filename, const unsigned ndim, MaterialContainer *matrs);
    // void saveToFile(const std :: string &filepath, std :: vector< unsigned > &elems_to_save) const;
    void saveElemStatsToFile(const std :: string &filepath, const std :: vector< unsigned > &elems_to_save, const double time_now = 0, const unsigned step = 0, const bool saveNodeIds = true, const double idc_time = 0, const double time_step = 1e-4) const;
    void readMatStatsFromFile(double &ini_time, unsigned &ini_step, double &ini_time_step, double &ini_idc_time, const bool get_time_from_file = true);
    void setFileToLoadStatsFrom(const std :: string &str);
    void init();
    void clear();
    size_t giveSize() const { return elems.size(); }
    void findElementFriends();
    void updateMaterialStatuses();
    void resetMaterialStatuses();
    void addElement(Element *newelem);
    CoordinateIndexedSparseMatrix prepareOutputStiffnessMatrix(bool BC_applied = true) const;
    CoordinateIndexedSparseMatrix updateOutputStiffnessMatrix(CoordinateIndexedSparseMatrix K_out, std :: string param, bool BC_applied = true, bool solver_numbering = true) const;
    void prepareStiffnessMatrix(CoordinateIndexedSparseMatrix &K) const;
    void updateStiffnessMatrix(CoordinateIndexedSparseMatrix &K, std :: string param) const;
    void prepareDampingMatrix(CoordinateIndexedSparseMatrix &C) const;
    void updateDampingMatrix(CoordinateIndexedSparseMatrix &C) const;
    void prepareMassMatrix(CoordinateIndexedSparseMatrix &M, bool lumped) const;
    void updateMassMatrix(CoordinateIndexedSparseMatrix &M, bool lumped) const;
    //void updateLumpedMassMatrix(Vector &M) const;
    double integrateKineticEnergy(const Vector &velocity) const;
    void integrateInternalForces(Vector &full_r, Vector &full_f, double time, double timeStep);                        ///< return internal forces with temporary update of internal variables
    void integrateInternalForcesWithFrozenIntVariables(Vector &full_r, Vector &full_f, double time, double timeStep);  ///< return internal forces based on current state of internal variables
    void integrateInternalForces(const Vector &full_r, Vector &full_f, bool frozen, double time, double timeStep);          ///< return internal forces with or without frozen internal variables
    void integrateDampingForces(const Vector &full_v, Vector &full_f) const;
    void integrateInertiaForces(const Vector &full_a, Vector &full_f) const;
    Element *giveElement(unsigned const num) const;
    Element *giveElementConnectingNodes(std :: vector< unsigned > &node_ids) const;
    unsigned giveElemId(const Element *elem) const;
    bool findElementOwningPoint(Element **elem, Point *xn, const Point *x) const;
    Element *findClosestElement(const Point *x) const;
    void extrapolateValuesFromIntegrationPointsToNodes(std :: string code, std :: vector< Vector > &results) const;
    void assignFibersToElems();
    void giveValues(std :: string code, Vector &result) const;
    void sumFromElements(std :: string code, Vector &result) const;
    std :: vector< Vector >computePrincipalStresses() const;
    void replaceTrueMassMatricesByLumpedOnes();
    double giveDissipatedEnergy() const;
    void resetEigenStrain(double time);

    std :: vector< Element * > :: iterator begin() { return elems.begin(); }
    std :: vector< Element * > :: iterator end() { return elems.end(); }
    std :: vector< Element * > :: const_iterator begin() const { return elems.begin(); }
    std :: vector< Element * > :: const_iterator end() const { return elems.end(); }
};



#endif  /* _ELEMENT_CONTAINER_H */
