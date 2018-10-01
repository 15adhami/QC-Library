//Just need to include qc.h

#include "qc.h"
#include <iostream>
#include <vector>

int main()
{
	//Quantum teleportation example
	//Bob will receive the information in the info qbit
	q_state alice(0);
	q_state bob(0);
	q_state info(1);
	q_state state;
	//Applying Hadamard gate to Alice's qbit
	alice=H(alice);
	//Getting the tensor product of Alice and Bob, storing it into a quantum state called ab
	q_state ab=alice*bob;
	//Applying controlled NOT gate from alice to bob, creating an entangled state
	ab=CNOT(ab);
	//Getting the tensor product of info and ab
	state=info*ab;
	//Applying CNOT to info and alice
	state=CNOT(state);
	//Applying Hadamard gate to info
	state=H(state);
	//Applying CNOT to alice and bob (1 means to apply to second qbit in the state)
	state=CNOT(state,1);
	//measure the state, bob should always measure what "info" originally had
	vector<int> out = state.measure();
	printbitstring(out);
	cout<<endl;

	return 0;
}