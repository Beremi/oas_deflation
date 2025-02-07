#ifndef _SOLVER_EIGEN_H
#define _SOLVER_EIGEN_H

#include  "solver.h"
#include <Eigen/Core>
#include <Eigen/SparseLU>

//////////////////////////////////////////////////////////
class EigenvalueMechanicalSolver : public Solver
{
protected:
    Vector lumpedM;
    CoordinateIndexedSparseMatrix M;
    unsigned num_eigs;

private:
public:
    EigenvalueMechanicalSolver();
    virtual ~EigenvalueMechanicalSolver();
    virtual void init(std :: string init_r_file, std :: string init_v_file, const bool initial = true);
    virtual void solve();
    virtual void giveValues(std :: string code, Vector &result) const;
    virtual Solver *readFromFile(const std :: string filename);
    virtual void runBeforeEachStep();
};

#endif /* _SOLVER_EIGEN_H */
