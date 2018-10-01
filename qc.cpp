//Imran Adham

#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include "qc.h"
//using namespace std;


int log(int base, int x)
{
	int temp=base;
	int count=1;
	while(temp<x)
	{
		temp=temp*base;
		count++;
	}
	return count;
}

void printbitstring(vector<int> v)
{
	for (int i=0;i<v.size();i++)
	{
		cout<<v[i];
	}
}

complex::complex()
{
	a=0.0;
	b=0.0;
}
complex::complex(double real)
{
	a=real;
	b=0;
}
complex::complex(double real, double imaginary)
{
	a=real;
	b=imaginary;
}
complex complex::operator=(const complex& c)
{
	this->a = c.a;
	this->b = c.b;
}

// Define the * operator
complex complex::operator*(const complex& c)
{
	complex ans;
	ans.a = this->a * c.a - this->b * c.b;
	ans.b = this->a * c.b + this->b * c.a;
	return ans;
}

	// Define the + operator
complex complex::operator+(const complex& c)
{
	complex ans;
	ans.a = this->a + c.a;
	ans.b = this->b + c.b;
	return ans;
}

	// Calculates the absolute value
double complex::abs()
{
	return sqrt(a*a+b*b);
}

	// Calculates the absolute value squared
double complex::abs_squared()
{
	return a*a+b*b;
}

	// Prints the complex number
void complex::print()
{
	if(b==0.0)
	{
		cout<<a;
	}
	else if (a==0.0)
	{
		cout<<b<<"i";
	}
	else
	{
		cout<<a<<"+"<<b<<"i";
	}
}

// Creates a complex object with the given parameters
complex to_complex(double real)
{
	complex c(real);
	return c;
}

// Creates a complex object with the given parameters
complex to_complex(double real, double imaginary)
{
	complex c(real, imaginary);
	return c;
}

// Calculates the tensor product of two vectors
vector<complex> tensor(vector<complex> a, vector<complex> b)
{
	vector<complex> out;
	for(int i=0;i<a.size();i++)
	{
		for(int j=0;j<b.size();j++)
		{
			out.push_back(a[i]*b[j]);
		}
	}
	return out;
}

// Calculates the tensor product of two matrices
vector<vector<complex>> tensor(vector<vector<complex>> a, vector<vector<complex>> b)
{
	vector<vector<complex>> out;
	vector<complex> temp;
	/*for(int i=0;i<b.size();i++)
	{
		for(int j=0;j<a.size();j++)
		{
			for(int k=0;k<b.size();k++)
			{
				for(int l=0;l<a.size();l++)
				{
					temp.push_back(a[j][l]*b[i][k]);
				}
			}
			out.push_back(temp);
			temp.clear();
		}
	}*/

	for(int i=0;i<a.size();i++)
	{
		for(int j=0;j<b.size();j++)
		{
			for(int k=0;k<a.size();k++)
			{
				for (int l=0;l<b.size();l++)
				{
					temp.push_back(a[i][k]*b[j][l]);
				}
			}
			out.push_back(temp);
			temp.clear();
		}
	}

	return out;
}

// Multiply a matrix and a vector, A*x
vector<complex> multiply(vector< vector<complex> > A, vector<complex> x)
{
	vector<complex> out;
	for (int i=0;i<A[0].size();i++)
	{
		complex temp=to_complex(0.0);
		for (int j=0;j<A.size();j++)
		{
			complex temp2=A[j][i]*x[j];
			complex temp3=temp+temp2;
			temp=temp3;
		}
		out.push_back(temp);
	}
	return out;
}

// Scale a vector by a complex number, a*v
vector<complex> scale(complex a, vector<complex> v)
{
	vector<complex> out;
	for (int i=0;i<v.size();i++)
	{
		out.push_back(v[i]*a);
	}
	return out;
}

vector<vector<complex>> scale(complex a, vector<vector<complex>> A)
{
	vector<vector<complex>> out;
	vector<complex> temp;
	for (int i=0;i<A.size();i++)
	{
		for (int j=0;j<A.size();j++)
		{
			temp.push_back(A[i][j]*a);
		}
		out.push_back(temp);
		temp.clear();
	}
	return out;
}

vector<vector<complex>> scale(double a, vector<vector<complex>> A)
{
	vector<vector<complex>> out;
	vector<complex> temp;
	for (int i=0;i<A.size();i++)
	{
		for (int j=0;j<A.size();j++)
		{
			temp.push_back(A[i][j]*to_complex(a));
		}
		out.push_back(temp);
		temp.clear();
	}
	return out;
}

// Turn U gate into CU
vector<vector<complex>> CU(vector<vector<complex>> U)
{
	vector< vector<complex> > gate;
	vector <complex> temp;

	for (int i=0;i<U.size();i++)
	{
		for (int j=0;j<2*U.size();j++)
		{
			if (i==j)
			{
				temp.push_back(to_complex(1.0));
			}
			else
			{
				temp.push_back(to_complex(0.0));
			}
		}
		gate.push_back(temp);
		temp.clear();
	}

	for (int i=0;i<U.size();i++)
	{
		for (int j=0;j<U.size();j++)
		{
			temp.push_back(to_complex(0.0));
		}
		for (int j=0;j<U.size();j++)
		{
			temp.push_back(U[i][j]);
		}
		gate.push_back(temp);
		temp.clear();
	}

	return gate;
}

// Creates the nxn identity matrix
vector<vector<complex>> I_gate(int n)
{
	vector<vector<complex>> mat;
	vector<complex> temp;
	for (int i=0;i<n;i++)
	{
		for (int j=0;j<n;j++)
		{
			if(i==j)
			{
				temp.push_back(to_complex(1.0));
			}
			else
			{
				temp.push_back(to_complex(0.0));
			}
		}
		mat.push_back(temp);
		temp.clear();
	}
	return mat;
}

vector<vector<complex>> I_gate()
{
	return I_gate(2);
}

vector<vector<complex>> H_gate()
{
	vector< vector<complex> > gate;
	vector <complex> temp;
	temp.push_back(to_complex(1.0));
	temp.push_back(to_complex(1.0));
	gate.push_back(temp);
	vector <complex> temp2;
	temp2.push_back(to_complex(1.0));
	temp2.push_back(to_complex(-1.0));
	gate.push_back(temp2);
	gate=scale(0.70711,gate);
	return gate;
}

vector<vector<complex>> X_gate()
{
	vector< vector<complex> > gate;
	vector <complex> temp;
	temp.push_back(to_complex(0.0));
	temp.push_back(to_complex(1.0));
	gate.push_back(temp);
	vector <complex> temp2;
	temp2.push_back(to_complex(1.0));
	temp2.push_back(to_complex(0.0));
	gate.push_back(temp2);
	return gate;
}

vector<vector<complex>> NOT_gate()
{
	vector< vector<complex> > gate;
	vector <complex> temp;
	temp.push_back(to_complex(0.0));
	temp.push_back(to_complex(1.0));
	gate.push_back(temp);
	vector <complex> temp2;
	temp2.push_back(to_complex(1.0));
	temp2.push_back(to_complex(0.0));
	gate.push_back(temp2);
	return gate;
}

vector<vector<complex>> Y_gate()
{
	vector< vector<complex> > gate;
	vector <complex> temp;
	temp.push_back(to_complex(0.0));
	temp.push_back(to_complex(0.0, 1.0));
	gate.push_back(temp);
	vector <complex> temp2;
	temp2.push_back(to_complex(0.0, -1.0));
	temp2.push_back(to_complex(0.0));
	gate.push_back(temp2);
	return gate;
}

vector<vector<complex>> Z_gate()
{
	vector< vector<complex> > gate;
	vector <complex> temp;
	temp.push_back(to_complex(1.0));
	temp.push_back(to_complex(0.0));
	gate.push_back(temp);
	vector <complex> temp2;
	temp2.push_back(to_complex(0.0));
	temp2.push_back(to_complex(-1.0));
	gate.push_back(temp2);
	return gate;
}

vector<vector<complex>> sqrtNOT_gate()
{
	vector< vector<complex> > gate;
	vector <complex> temp;
	temp.push_back(to_complex(0.5, 0.5));
	temp.push_back(to_complex(0.5, -0.5));
	gate.push_back(temp);
	vector <complex> temp2;
	temp2.push_back(to_complex(0.5, -0.5));
	temp2.push_back(to_complex(0.5, 0.5));
	gate.push_back(temp2);
	return gate;
}

vector<vector<complex>> R_gate(double phi)
{
	vector< vector<complex> > gate;
	vector <complex> temp;
	temp.push_back(to_complex(1.0));
	temp.push_back(to_complex(0.0));
	gate.push_back(temp);
	vector <complex> temp2;
	temp2.push_back(to_complex(0.0));
	temp2.push_back(to_complex(cos(phi))+to_complex(0.0,1.0)*to_complex(sin(phi)));
	gate.push_back(temp2);
	return gate;
}

vector<vector<complex>> SWAP_gate()
{
	vector< vector<complex> > gate;
	vector <complex> temp;
	temp.push_back(to_complex(1.0));
	temp.push_back(to_complex(0.0));
	temp.push_back(to_complex(0.0));
	temp.push_back(to_complex(0.0));
	gate.push_back(temp);
	temp.clear();
	temp.push_back(to_complex(0.0));
	temp.push_back(to_complex(0.0));
	temp.push_back(to_complex(1.0));
	temp.push_back(to_complex(0.0));
	gate.push_back(temp);
	temp.clear();
	temp.push_back(to_complex(0.0));
	temp.push_back(to_complex(1.0));
	temp.push_back(to_complex(0.0));
	temp.push_back(to_complex(0.0));
	gate.push_back(temp);
	temp.clear();
	temp.push_back(to_complex(0.0));
	temp.push_back(to_complex(0.0));
	temp.push_back(to_complex(0.0));
	temp.push_back(to_complex(1.0));
	gate.push_back(temp);
	return gate;
}

vector<vector<complex>> sqrtSWAP_gate()
{
	vector< vector<complex> > gate;
	vector <complex> temp;
	temp.push_back(to_complex(1.0));
	temp.push_back(to_complex(0.0));
	temp.push_back(to_complex(0.0));
	temp.push_back(to_complex(0.0));
	gate.push_back(temp);
	temp.clear();
	temp.push_back(to_complex(0.0));
	temp.push_back(to_complex(0.5, 0.5));
	temp.push_back(to_complex(0.5, -0.5));
	temp.push_back(to_complex(0.0));
	gate.push_back(temp);
	temp.clear();
	temp.push_back(to_complex(0.0));
	temp.push_back(to_complex(0.5, -0.5));
	temp.push_back(to_complex(0.5, 0.5));
	temp.push_back(to_complex(0.0));
	gate.push_back(temp);
	temp.clear();
	temp.push_back(to_complex(0.0));
	temp.push_back(to_complex(0.0));
	temp.push_back(to_complex(0.0));
	temp.push_back(to_complex(1.0));
	gate.push_back(temp);
	return gate;
}

vector<vector<complex>> CNOT_gate()
{
	vector< vector<complex> > gate=X_gate();
	gate=CU(gate);
	return gate;
}

vector<vector<complex>> CX_gate()
{
	vector< vector<complex> > gate=X_gate();
	gate=CU(gate);
	return gate;
}

vector<vector<complex>> CY_gate()
{
	vector< vector<complex> > gate=Y_gate();
	gate=CU(gate);
	return gate;
}

vector<vector<complex>> CZ_gate()
{
	vector< vector<complex> > gate=Z_gate();
	gate=CU(gate);
	return gate;
}

vector<vector<complex>> CCNOT_gate()
{
	vector< vector<complex> > gate=NOT_gate();
	gate=CU(gate);
	gate=CU(gate);
	return gate;
}

vector<vector<complex>> CSWAP_gate()
{
	vector< vector<complex> > gate=SWAP_gate();
	gate=CU(gate);
	return gate;
}

vector<vector<complex>> XX_gate(double phi)
{
	vector< vector<complex> > gate;
	vector <complex> temp;
	temp.push_back(to_complex(1.0));
	temp.push_back(to_complex(0.0));
	temp.push_back(to_complex(0.0));
	temp.push_back(to_complex(cos(-1.0*phi-1.5708))+to_complex(0.0,1.0)*to_complex(sin(-1.0*phi-1.5708)));
	gate.push_back(temp);
	temp.clear();
	temp.push_back(to_complex(0.0));
	temp.push_back(to_complex(1.0));
	temp.push_back(to_complex(0.0, -1.0));
	temp.push_back(to_complex(0.0));
	gate.push_back(temp);
	temp.clear();
	temp.push_back(to_complex(0.0));
	temp.push_back(to_complex(0.0, -1.0));
	temp.push_back(to_complex(1.0));
	temp.push_back(to_complex(0.0));
	gate.push_back(temp);
	temp.clear();
	temp.push_back(to_complex(cos(phi-1.5708))+to_complex(0.0,1.0)*to_complex(sin(phi-1.5708)));
	temp.push_back(to_complex(0.0));
	temp.push_back(to_complex(0.0));
	temp.push_back(to_complex(1.0));
	gate.push_back(temp);
	gate=scale(0.70711,gate);
	return gate;
}

// Calculates base^exp
int power(int base, int exp)
{
	int temp=1;
	for(int i=0;i<exp;i++)
	{
		temp=temp*base;
	}
	return temp;
}

	// Defaults to |0>
q_state::q_state()
{
	state.push_back(to_complex(1.0));
	state.push_back(to_complex(0.0));
	dimension=state.size();
	bitlength=1;
}

	// Initializes to |init>
q_state::q_state(int init)
{
	if (init==0)
	{
		state.push_back(to_complex(1.0));
		state.push_back(to_complex(0.0));
		dimension=state.size();
	}
	else
	{
		state.push_back(to_complex(0.0));
		state.push_back(to_complex(1.0));
		dimension=state.size();
	}
	bitlength=1;
}

	// Initializes to a custom bitstring
q_state::q_state(vector<int> init)
{
	vector< vector<complex> > vec;
	for (int i=0;i<init.size();i++)
	{
		if (init[i]==0)
		{
			vector<complex> temp;
			temp.push_back(to_complex(1.0));
			temp.push_back(to_complex(0.0));
			vec.push_back(temp);
		}
		if (init[i]==1)
		{
			vector<complex> temp;
			temp.push_back(to_complex(0.0));
			temp.push_back(to_complex(1.0));
			vec.push_back(temp);
		}
	}	
	vector<complex> temp=vec[0];
	for (int i=0;i<vec.size()-1;i++){
		temp=tensor(temp,vec[i+1]);
	}
	state=temp;
	dimension=state.size();
	bitlength=init.size();
}

	// Initializes to a custom state vector
q_state::q_state(vector<complex> init)
{
	state=init;
	//cout<<state.size()<<endl;
	dimension=state.size();
	bitlength= (int) (log2(dimension)+0.2);
}
	// Return dimension
int q_state::get_dimension()
{
	return dimension;
}

	// Return bitlength
int q_state::get_bitlength()
{
	return bitlength;
}

void q_state::operator=(const q_state& c)
{
	this->state.clear();
	for(int i=0;i<c.dimension;i++)
	{
		this->state.push_back(c.state[i]);
	}
	this->dimension=c.dimension;
	this->bitlength=c.bitlength;
}

q_state q_state::operator*(const q_state& c)
{
	vector<complex> out;
	for(int i=0;i<this->dimension;i++)
	{
		for(int j=0;j<c.dimension;j++)
		{
			out.push_back(this->state[i]*c.state[j]);
		}
	}
	q_state temp(out);
	return temp;
}

q_state::~q_state()
{
	state.clear();
}

	// Return state vector
vector<complex> q_state::get_state()
{
	return state;
}

	// Copy constructor
/*q_state::q_state(q_state init)
{
	state=init.get_state();
	dimension=init.get_dimension();
	bitlength=init.get_bitlength();
}*/

	// Measure the quantum state
vector<int> q_state::measure()
{
	std::srand((unsigned)time(NULL));
	vector<int> v;
	int count=0;
	double num=(double) rand()/RAND_MAX;
	for (int i=0;i<dimension;i++)
	{
		if (num<=state[i].abs_squared())
		{
			break;
		}
		else
		{
			if (i==dimension-1)
			{
				break;
			}
			num-=state[i].abs_squared();
			count++;
		}
	}
	for (int i=0;i<dimension;i++)
	{
		if(i==count)
		{
			state[i]=to_complex(1.0);
		}
		else
		{
			state[i]=to_complex(0.0);
		}
	}
	vector<int> ans;
	for(int i=0;i<bitlength;i++)
	{
		if(count>=power(2,bitlength-i-1))
		{
			ans.push_back(1);
			count-=power(2,bitlength-i-1);
		}
		else
		{
			ans.push_back(0);
		}
	}
	return ans;
}

	// Prints the state vector and dimension
void q_state::print()
{
	cout<<"State: ";
	for (int i=0;i<dimension;i++)
	{
		state[i].print();
		cout<<" ";
	}
	cout<<endl<<"Dimension: "<<dimension<<"     Bits: "<<bitlength<<endl;
}


// Calculates the tensor product of two q_state objects
q_state tensor(q_state a, q_state b)
{
	vector<complex> out;
	for(int i=0;i<a.get_dimension();i++)
	{
		for(int j=0;j<b.get_dimension();j++)
		{
			out.push_back(a.get_state()[i]*b.get_state()[j]);
		}
	}
	q_state temp(out);
	return temp;
}

// Apply a custom gate at the location bit, indexing from 0 (EX. type 1 to apply from 2nd bit, or 0 for 1st bit)
// default is 0
// gate is formatted as a vector of column vectors
q_state apply_gate(vector< vector<complex> > gate, q_state q, int location)
{
	vector< vector<complex> > temp=gate;
	for (int i=0;i<location;i++)
	{
		temp=tensor(I_gate(),temp);
	}
	for (int i=0;i<q.get_bitlength()-( ((int) (log2(gate.size())+0.1)) + location);i++)
	{
		temp=tensor(temp,I_gate());
	}
	vector<complex> ans = multiply(temp,q.get_state());
	q_state out(ans);
	return out;
}

// Apply a custom gate at location=0
// gate is formatted as a vector of column vectors
q_state apply_gate(vector< vector<complex> > gate, q_state q)
{
	return apply_gate(gate, q, 0);
}

// Apply Hadamard gate
q_state H(q_state q)
{
	q_state out=apply_gate(H_gate(), q.get_state());
	return out;
}

// Apply Hadamard gate
q_state H(q_state q, int location)
{
	q_state out=apply_gate(H_gate(), q.get_state(), location);
	return out;
}

// Apply Pauli-X gate (same as NOT)
q_state X(q_state q)
{
	q_state out=apply_gate(X_gate(), q.get_state());
	return out;
}

// Apply Pauli-X gate
q_state X(q_state q, int location)
{
	q_state out=apply_gate(X_gate(), q.get_state(), location);
	return out;
}

// Apply NOT gate (same as Pauli-X)
q_state NOT(q_state q)
{
	q_state out=apply_gate(NOT_gate(), q.get_state());
	return out;
}

// Apply NOT gate (same as Pauli-X)
q_state NOT(q_state q, int location)
{
	q_state out=apply_gate(NOT_gate(), q.get_state(), location);
	return out;
}

// Apply Pauli-Y gate
q_state Y(q_state q)
{
	q_state out=apply_gate(Y_gate(), q.get_state());
	return out;
}

// Apply Pauli-Y gate
q_state Y(q_state q, int location)
{
	q_state out=apply_gate(Y_gate(), q.get_state(), location);
	return out;
}

// Apply Pauli-Z gate
q_state Z(q_state q)
{
	q_state out=apply_gate(Z_gate(), q.get_state());
	return out;
}

// Apply Pauli-Z gate
q_state Z(q_state q, int location)
{
	q_state out=apply_gate(Z_gate(), q.get_state(), location);
	return out;
}

// Apply square root of NOT gate
q_state sqrtNOT(q_state q)
{
	q_state out=apply_gate(sqrtNOT_gate(), q.get_state());
	return out;
}

// Apply square root of NOT gate
q_state sqrtNOT(q_state q, int location)
{
	q_state out=apply_gate(sqrtNOT_gate(), q.get_state(), location);
	return out;
}

// Apply phase shift gate by phi radians
q_state R(double phi, q_state q)
{
	q_state out=apply_gate(R_gate(phi), q.get_state());
	return out;
}

// Apply phase shift gate by phi radians
q_state R(double phi, q_state q, int location)
{
	q_state out=apply_gate(R_gate(phi), q.get_state(), location);
	return out;
}

// Apply SWAP gate
q_state SWAP(q_state q)
{
	q_state out=apply_gate(SWAP_gate(), q.get_state());
	return out;
}

// Apply SWAP gate
q_state SWAP(q_state q, int location)
{
	q_state out=apply_gate(SWAP_gate(), q.get_state(), location);
	return out;
}

// Apply square root of SWAP gate
q_state sqrtSWAP(q_state q)
{
	q_state out=apply_gate(sqrtSWAP_gate(), q.get_state());
	return out;
}

// Apply square root of SWAP gate
q_state sqrtSWAP(q_state q, int location)
{
	q_state out=apply_gate(sqrtSWAP_gate(), q.get_state(), location);
	return out;
}

// Apply CNOT gate (same as cX)
q_state CNOT(q_state q)
{
	q_state out=apply_gate(CNOT_gate(), q.get_state());
	return out;
}

// Apply CNOT gate (same as cX)
q_state CNOT(q_state q, int location)
{
	q_state out=apply_gate(CNOT_gate(), q.get_state(), location);
	return out;
}

// Apply cX gate (same as CNOT)
q_state CX(q_state q)
{
	q_state out=apply_gate(CX_gate(), q.get_state());
	return out;
}

// Apply cX gate (same as CNOT)
q_state CX(q_state q, int location)
{
	q_state out=apply_gate(CX_gate(), q.get_state(), location);
	return out;
}

// Apply cY gate
q_state CY(q_state q)
{
	q_state out=apply_gate(CY_gate(), q.get_state());
	return out;
}

// Apply cY gate
q_state CY(q_state q, int location)
{
	q_state out=apply_gate(CY_gate(), q.get_state(), location);
	return out;
}

// Apply cZ gate
q_state CZ(q_state q)
{
	q_state out=apply_gate(CZ_gate(), q.get_state());
	return out;
}

// Apply cZ gate
q_state CZ(q_state q, int location)
{
	q_state out=apply_gate(CZ_gate(), q.get_state(), location);
	return out;
}

// Apply CCNOT gate (Toffoli gate)
q_state CCNOT(q_state q)
{
	q_state out=apply_gate(CCNOT_gate(), q.get_state());
	return out;
}

// Apply CCNOT gate (Toffoli gate)
q_state CCNOT(q_state q, int location)
{
	q_state out=apply_gate(CCNOT_gate(), q.get_state(), location);
	return out;
}

// Apply CSWAP gate (Fredkin gate)
q_state CSWAP(q_state q)
{
	q_state out=apply_gate(CSWAP_gate(), q.get_state());
	return out;
}

// Apply CSWAP gate (Fredkin gate)
q_state CSWAP(q_state q, int location)
{
	q_state out=apply_gate(CSWAP_gate(), q.get_state(), location);
	return out;
}

// Apply XX gate (Ising gate)
q_state XX(double phi, q_state q)
{
	q_state out=apply_gate(XX_gate(phi), q.get_state());
	return out;
}

// Apply XX gate (Ising gate)
q_state XX(double phi, q_state q, int location)
{
	q_state out=apply_gate(XX_gate(phi), q.get_state(), location);
	return out;
}