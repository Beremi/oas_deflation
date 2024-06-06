#include "element_superelem.h"
#include "boundary_condition.h"

using namespace std;

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// ML ELEMENT
MLMechElement :: MLMechElement(unsigned dim) :Element(dim) {
    numOfNodes = 0;
    name = "MLElement";
    vtk_cell_type = 0;
    shafunc = new NullShapeF(ndim);
    inttype = new EmptyIntegration();
    physicalFields [ 0 ] = true; //mechanics
}

//////////////////////////////////////////////////////////
Matrix MLMechElement :: readStiffMatrixFromFile() const {
    vector<double> matrixEntries;
    ifstream matrixDataFile(sm_path);
    string matrixRowString;
    string matrixEntry;
    int matrixRowNumber = 0;
    while (getline(matrixDataFile, matrixRowString))
    {
        stringstream matrixRowStringStream(matrixRowString); 
        while (matrixRowStringStream >> matrixEntry)
        {
            matrixEntries.push_back(stod(matrixEntry));   //here we convert the string to double and fill in the row vector storing all the matrix entries
        }
        matrixRowNumber++; //update the column numbers
    } 

    if (matrixEntries.size()!=matrixRowNumber*matrixRowNumber) {
        cerr << "Loaded matrix is not a square matrix, found "<< matrixEntries.size() << ", expected "<< matrixRowNumber*matrixRowNumber << endl;
        exit(1);
    }
    Matrix A(matrixRowNumber,matrixRowNumber);
    for(int i=0; i<matrixRowNumber; i++){
        for(int j=0; j<matrixRowNumber; j++){
            A(i,j) = matrixEntries[i*matrixRowNumber+j];
        }   
    }

    return A;
}

//////////////////////////////////////////////////////////
void MLMechElement :: readFromLine(std :: istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs){
    (void) fullmatrs;

    unsigned num;
    iss >> numOfNodes;
    nodes.resize(numOfNodes);
    for ( unsigned k = 0; k < numOfNodes; k++ ) {
        iss >> num;
        nodes [ k ] = fullnodes->giveNode(num);
    }

    string param;    
    while (  iss >> param ) {
        if ( param.compare("poly_degree") == 0 ) {            
            iss >> poly_degree;                     
        } else if ( param.compare("stiff_mat") == 0 ) {        
            string filepath;    
            iss >> filepath;
            sm_path = GlobPaths :: BASEDIR  / filepath;
        }
    }

}

//////////////////////////////////////////////////////////
void MLMechElement :: init(){
    Element::init();

    //reduce stiffness matrix
    vector<unsigned>keep_ind; 
    //nodal displacement
    for(unsigned i=0; i<8; i++){
        keep_ind.push_back(i);
    }
    // i runs over sides
    for(unsigned i=0; i<8; i++){        
        // j runs over polynomial constants, linear basis means no additional dofs
        for(unsigned j=0; j<poly_degree-1; j++){       
            keep_ind.push_back(8+i*4+j);
        }
    } 


    Matrix A = readStiffMatrixFromFile();
    stiffmat = A(keep_ind,keep_ind);    

    if(outDoFs!=keep_ind.size()){
        cerr << "Error in MLMechElement: there are " << outDoFs << " DoFs on input of the element, but " << keep_ind.size() << " DoFs are required for polynom of degree " << poly_degree << endl;
        exit(1);
    }
}

//////////////////////////////////////////////////////////
Matrix MLMechElement :: giveStiffnessMatrix(std :: string matrixType) const{
    (void) matrixType;
    return stiffmat;
}

//////////////////////////////////////////////////////////
Vector MLMechElement :: giveInternalForces(const Vector &DoFs, bool frozen, double timeStep){
    (void) timeStep;
    if (frozen){
        return stiffmat*DoFs;
    } else {
        return stiffmat*DoFs;
    }
}

//////////////////////////////////////////////////////////
Vector MLMechElement :: integrateInternalSources(){
    return Vector::Zero(0);
}

//////////////////////////////////////////////////////////
void MLMechElement :: giveValues(std :: string code, Vector &result) const{
    (void) code;    (void) result;
}


/*
//////////////////////////////////////////////////////////
Matrix MLMechElement :: giveBMatrix(const Point *x) const {
    cout << "MLMechElement :: giveBMatrix should not be called" << endl;
    exit(1);
    return Matrix::Zeros(0,0);
}

//////////////////////////////////////////////////////////
Matrix MLMechElement :: giveHMatrix(const Point *x) const {
    cout << "MLMechElement :: giveHMatrix should not be called" << endl;
    exit(1);
    return Matrix::Zeros(0,0);
}
*/
//////////////////////////////////////////////////////////
Vector MLMechElement :: giveStrain(unsigned i, const Vector &DoFs) {
    (void) i; (void) DoFs;
    cout << "MLMechElement :: giveStrain should not be called" << endl;
    exit(1);
    return Vector::Zero(0);
}


//////////////////////////////////////////////////////////
Matrix MLMechElement :: giveMassMatrix(){
    cout << "MLMechElement :: giveBMatrix should not be called" << endl;
    exit(1);
    return Matrix::Zero(0,0);
}

//////////////////////////////////////////////////////////
Matrix MLMechElement :: giveDampingMatrix(){
    cout << "MLMechElement :: giveBMatrix should not be called" << endl;
    exit(1);
    return Matrix::Zero(0,0);
}

//////////////////////////////////////////////////////////
Vector MLMechElement :: giveLumpedMassMatrix(){
    cout << "MLMechElement :: giveBMatrix should not be called" << endl;
    exit(1);
    return Vector::Zero(0);
}

