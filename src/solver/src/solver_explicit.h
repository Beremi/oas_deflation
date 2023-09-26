#ifndef _SOLVER_EXPLICIT_H
#define _SOLVER_EXPLICIT_H

#include  "solver.h"
#include <Eigen/Core>
#include <Eigen/SparseLU>

//////////////////////////////////////////////////////////
class TransientCentralDifferenceMechanicalSolver : public Solver
{
protected:
    Vector v_red_old, v_red, a_red, v, a;    
    CoordinateIndexedSparseMatrix M, Minv;

    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;

    unsigned show_period;
private:
public:
    TransientCentralDifferenceMechanicalSolver();
    virtual ~TransientCentralDifferenceMechanicalSolver();  
    virtual void init(std :: string init_r_file, std :: string init_v_file, const bool initial = true);
    virtual void solve();
    virtual void giveValues(std :: string code, Vector &result) const;
    void computeAcceleration();
    virtual Solver *readFromFile(const std::string filename);
    virtual void runBeforeEachStep();
};

#endif /* _SOLVER_EXPLICIT_H */





