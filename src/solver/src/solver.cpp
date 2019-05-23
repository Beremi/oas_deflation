#include "solver.h"

//////////////////////////////////////////////////////////
Solver* Solver::readFromLine(istringstream &iss){
    string param;
    iss >> param;
    if (param.compare("SteadyStateLinearSolver") == 0){
        SteadyStateLinearSolver* newsolver = new SteadyStateLinearSolver();
        newsolver->readFromLine(iss);
        return newsolver;

    }else {cerr << "Error: Solver "<< param << " is not implemented" << endl; exit(0);};
}

//////////////////////////////////////////////////////////
void Solver::runBeforeEachStep(){
    step += 1; 
    time +=dt;  
}

//////////////////////////////////////////////////////////
void Solver::runAfterEachStep(){   
    if (abs(time-termination_time)<1e-14) terminated=true;
    if (time+dt>termination_time) dt -= time+dt-termination_time;
    nodes->giveFullDoFArray(r,pbc, full_r);
}

//////////////////////////////////////////////////////////
void Solver::init(){
    step = -1;
    time = 0.;
    freeDoFnum = nodes->giveNumFreeDoFs();
    fixedDoFnum = (nodes->giveTotalNumDoFs()-freeDoFnum);
    
    elems->prepareSteadyStateMatrices(K11, K12);
    elems->updateSteadyStateMatrices(K11, K12);
    f = Vector(freeDoFnum);
    r = Vector(freeDoFnum);
    dr = Vector(freeDoFnum);
    pbc = Vector(fixedDoFnum);
    full_r = Vector(freeDoFnum+fixedDoFnum);
}

//////////////////////////////////////////////////////////
SteadyStateLinearSolver::SteadyStateLinearSolver(){
    name = "SteadyStateLinearSolver";
}

//////////////////////////////////////////////////////////
SteadyStateLinearSolver::~SteadyStateLinearSolver(){
}

//////////////////////////////////////////////////////////
void SteadyStateLinearSolver::init(){
    Solver::init();
}

//////////////////////////////////////////////////////////
Solver* SteadyStateLinearSolver::readFromLine(istringstream &iss){
    string param;
    bool bdt,bttime;
    bdt=bttime = false;

    while (!iss.eof()) {
        iss >> param;
        if (param.compare("time_step") == 0){
            bdt = true;
            iss >> dt;
        }else if (param.compare("total_time") == 0){
            bttime = true;
            iss >> termination_time;
        }        
    }    
    if (not bdt) {cerr << name << ": solver parameter 'time_step' was not specified" << endl; exit(0);};
    if (not bttime) {cerr << name << ": solver parameter 'total_time' was not specified" << endl; exit(0);};
    cout << name << " succesfully loaded"<< endl;
    return this;
};

//////////////////////////////////////////////////////////
void SteadyStateLinearSolver::solve(){

    cout << "######### Solving step "<< step << " at time " << time << "; time step " << dt << " #########" << endl;
    f = f*0.;  //clear nodal load
    nodes->addRHS_nodalLoad(f, time);  //add nodal load
    nodes->updateDirrichletBC(pbc, time); //give prescribed DoFs

    if (ConjGrad(K11, dr, f-K12*pbc-K11*r, dr) == false) cerr << "Conjugate gradients did not converge" << endl;
    
    for(int i=0; i<freeDoFnum; i++) r[i] += dr[i];
    
}

//////////////////////////////////////////////////////////
void SteadyStateLinearSolver::runBeforeEachStep(){
    Solver::runBeforeEachStep();
}
//////////////////////////////////////////////////////////
void SteadyStateLinearSolver::runAfterEachStep(){ 
    Solver::runAfterEachStep();    
}
