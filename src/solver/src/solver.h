#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include <iostream>
#include <fstream>
#include "linear_algebra.h"
#include "element_container.h"
#include "node_container.h"


class Solver{
    protected:
        string name;
        ElementContainer *elems;
        NodeContainer *nodes;
        double time, dt, termination_time;        
        CoordinateIndexedSparseMatrix K11,K12;
        Vector f, pbc, r, dr, full_r;
        unsigned freeDoFnum, fixedDoFnum;
        int step;
        bool terminated;
        virtual void runBeforeEachStep();
        virtual void runAfterEachStep();
        
    public:
        Solver(){name="basic solver"; terminated=false;};
        ~Solver(){};
        virtual void init();
        virtual Solver* readFromLine(istringstream &iss);
        virtual void solveStep(){};
        void setElementContainer(ElementContainer *e){elems = e;}
        void setNodeContainer(NodeContainer *n){nodes=n;};
        string giveName()const {return name;}
        bool isTerminated(){return terminated;}
        Vector giveFullDoFValues(){return full_r;}
        int giveStepNumber()const{return step;};
};

class SteadyStateLinearSolver: public Solver{
    protected:
        virtual void runBeforeEachStep();
        virtual void runAfterEachStep();
        virtual void solve();
    private:

    public:
        SteadyStateLinearSolver();        
        ~SteadyStateLinearSolver(); //destructor
        virtual void init();
        virtual void solveStep(){runBeforeEachStep(); solve(); runAfterEachStep();};
        virtual Solver* readFromLine(istringstream &iss);
};

#endif
