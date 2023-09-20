#ifndef _SOLVER_EXPLICIT_H
#define _SOLVER_EXPLICIT_H

#include  "solver.h"

//////////////////////////////////////////////////////////
class TransientCentralDifferenceMechanicalSolver : public Solver
{
protected:
    Vector v_old, v, accel;    
private:
public:
    TransientCentralDifferenceMechanicalSolver();
    virtual ~TransientCentralDifferenceMechanicalSolver();  
    virtual void init(std :: string init_r_file, std :: string init_v_file, const bool initial = true);
    virtual void solve();
    virtual void giveValues(std :: string code, Vector &result) const;
};

#endif /* _SOLVER_EXPLICIT_H */





