//Imran Adham

#include<vector>
#ifndef QC
#define QC
using namespace std;

//Calculates log_{base}(x), rounded up
int log(int base, int x);


//Prints a vector of ints as a number
void printbitstring(vector<int> v);

// Class to hold complex numbers
// a is the real part
// b is the imaginary part
class complex
{
public:
	double a;
	double b;
	complex();
	complex(double real);
	complex(double real, double imaginary);

	// Define the = operator
	complex operator=(const complex& c);

	// Define the * operator
	complex operator*(const complex& c);

	// Define the + operator
	complex operator+(const complex& c);

	// Calculates the absolute value
	double abs();

	// Calculates the absolute value squared
	double abs_squared();

	// Prints the complex number
	void print();
};

// Creates a complex object with the given parameters
complex to_complex(double real);

// Creates a complex object with the given parameters
complex to_complex(double real, double imaginary);

// Calculates the tensor product of two vectors
vector<complex> tensor(vector<complex> a, vector<complex> b);

// Calculates the tensor product of two matrices
vector<vector<complex>> tensor(vector<vector<complex>> a, vector<vector<complex>> b);

// Multiply a matrix and a vector, A*x
vector<complex> multiply(vector< vector<complex> > A, vector<complex> x);

// Scale a vector by a complex number, a*v
vector<complex> scale(complex a, vector<complex> v);

// Scale a matrix by a complex number, a*A
vector<vector<complex>> scale(complex a, vector<vector<complex>> A);

// Scale a matrix by a real number, a*A
vector<vector<complex>> scale(double a, vector<vector<complex>> A);

// Turn U gate into CU
vector<vector<complex>> CU(vector<vector<complex>> U);

// Creates the nxn identity matrix
vector<vector<complex>> I_gate(int n);

// Creates the 2x2 identity matrix
vector<vector<complex>> I_gate();

// Creates the matrix for the H gate
vector<vector<complex>> H_gate();

// Creates the matrix for the Pauli-X gate (same as NOT)
vector<vector<complex>> X_gate();

// Creates the matrix for the NOT gate (same as Pauli-X)
vector<vector<complex>> NOT_gate();

// Creates the matrix for the Pauli-Y gate
vector<vector<complex>> Y_gate();

// Creates the matrix for the Pauli-Z gate
vector<vector<complex>> Z_gate();

// Creates the matrix for the sqrtNOT gate
vector<vector<complex>> sqrtNOT_gate();

// Creates the matrix for the phase shift gate
vector<vector<complex>> R_gate(double phi);

// Creates the matrix for the SWAP gate
vector<vector<complex>> SWAP_gate();

// Creates the matrix for the sqrtSWAP gate
vector<vector<complex>> sqrtSWAP_gate();

// Creates the matrix for the CNOT gate
vector<vector<complex>> CNOT_gate();

// Creates the matrix for the CX gate
vector<vector<complex>> CX_gate();

// Creates the matrix for the CY gate
vector<vector<complex>> CY_gate();

// Creates the matrix for the CZ gate
vector<vector<complex>> CZ_gate();

// Creates the matrix for the CCNOT gate
vector<vector<complex>> CCNOT_gate();

// Creates the matrix for the CSWAP gate
vector<vector<complex>> CSWAP_gate();

// Creates the matrix for the XX (Ising) gate
vector<vector<complex>> XX_gate(double phi);

// Calculates base^exp
int power(int base, int exp);

class q_state
{
private:
	vector<complex> state;
	int dimension;
	int bitlength;
public:

	// Defaults to |0>
	q_state();

	// Initializes to |init>
	q_state(int init);

	// Initializes to a custom bitstring
	q_state(vector<int> init);

	// Initializes to a custom state vector
	q_state(vector<complex> init);

	// Return dimension of state vector
	int get_dimension();

	// Return bitlength
	int get_bitlength();

	// Return state vector
	vector<complex> get_state();

	// Copy constructor
	//q_state(q_state init);

	~q_state();

	// Define the = operator
	void operator=(const q_state& c);

	// Define the * operator (tensor product)
	q_state operator*(const q_state& c);

	// Measure the quantum state
	vector<int> measure();

	// Prints the state vector and dimension
	void print();
};

// Calculates the tensor product of two q_state objects
q_state tensor(q_state a, q_state b);

// Apply a custom gate at the location bit, indexing from 0 (EX. type 1 to apply from 2nd bit, or 0 for 1st bit)
// gate is formatted as a vector of column vectors
q_state apply_gate(vector< vector<complex> > gate, q_state q, int location);

// Apply a custom gate
// gate is formatted as a vector of column vectors
q_state apply_gate(vector< vector<complex> > gate, q_state q);

// Apply Hadamard gate
q_state H(q_state q);

// Apply Hadamard gate
q_state H(q_state q, int location);

// Apply Pauli-X gate (same as NOT)
q_state X(q_state q);

// Apply Pauli-X gate (same as NOT)
q_state X(q_state q, int location);

// Apply NOT gate (same as Pauli-X)
q_state NOT(q_state q);

// Apply NOT gate (same as Pauli-X)
q_state NOT(q_state q, int location);

// Apply Pauli-Y gate
q_state Y(q_state q);

// Apply Pauli-Y gate
q_state Y(q_state q, int location);

// Apply Pauli-Z gate
q_state Z(q_state q);

// Apply Pauli-Z gate
q_state Z(q_state q, int location);

// Apply square root of NOT gate
q_state sqrtNOT(q_state q);

// Apply square root of NOT gate
q_state sqrtNOT(q_state q, int location);

// Apply phase shift gate by phi radians
q_state R(double phi, q_state q);

// Apply phase shift gate by phi radians
q_state R(double phi, q_state q, int location);

// Apply SWAP gate
q_state SWAP(q_state q);

// Apply SWAP gate
q_state SWAP(q_state q, int location);

// Apply square root of SWAP gate
q_state sqrtSWAP(q_state q);

// Apply square root of SWAP gate
q_state sqrtSWAP(q_state q, int location);

// Apply CNOT gate (same as cX)
q_state CNOT(q_state q);

// Apply CNOT gate (same as cX)
q_state CNOT(q_state q, int location);

// Apply cX gate (same as CNOT)
q_state CX(q_state q);

// Apply cX gate (same as CNOT)
q_state CX(q_state q, int location);

// Apply cY gate
q_state CY(q_state q);

// Apply cY gate
q_state CY(q_state q, int location);

// Apply cZ gate
q_state CZ(q_state q);

// Apply cZ gate
q_state CZ(q_state q, int location);

// Apply CCNOT gate (Toffoli gate)
q_state CCNOT(q_state q);

// Apply CCNOT gate (Toffoli gate)
q_state CCNOT(q_state q, int location);

// Apply CSWAP gate (Fredkin gate)
q_state CSWAP(q_state q);

// Apply CSWAP gate (Fredkin gate)
q_state CSWAP(q_state q, int location);

// Apply XX gate (Ising gate)
q_state XX(double phi, q_state q);

// Apply XX gate (Ising gate)
q_state XX(double phi, q_state q, int location);

#endif