import pennylane as qml
import numpy as np
import scipy.optimize

class PennyLaneTSPSolverNaive:
    def __init__(self, distance_matrix, steps=1, ftol=1.0e-2, xtol=1.0e-2, use_constraints=False):
        self.distance_matrix = distance_matrix
        self.steps = steps
        self.ftol = ftol
        self.xtol = xtol
        self.use_constraints = use_constraints
        self.number_of_qubits = self.get_number_of_qubits()
        self.dev = qml.device("default.qubit", wires=self.number_of_qubits)
        self.circuit = qml.QNode(self.circuit_definition, self.dev)
        self.probability_circuit = qml.QNode(self.probability_circuit_definition, self.dev)
        self.sampling_results = {}
    def get_number_of_qubits(self):
        return len(self.distance_matrix)**2

    def create_cost_hamiltonian(self):
        cost_operators = []
        number_of_nodes = len(self.distance_matrix)
        
        for i in range(number_of_nodes):
            for j in range(i+1, number_of_nodes):
                for t in range(number_of_nodes - 1):
                    weight = self.distance_matrix[i][j]
                    qubit_1 = t * number_of_nodes + i
                    qubit_2 = (t + 1) * number_of_nodes + j
                    cost_operators.append((weight, [qubit_1, qubit_2], "ZZ"))

        return qml.Hamiltonian([op[0] for op in cost_operators], [qml.PauliZ(op[1][0]) @ qml.PauliZ(op[1][1]) for op in cost_operators])

    def circuit_definition(self, params):
        cost_hamiltonian = self.create_cost_hamiltonian()
        
        for i in range(self.number_of_qubits):
            qml.Hadamard(wires=i)
        
        for step in range(self.steps):
            qml.ApproxTimeEvolution(cost_hamiltonian, params[step], 1)
            for i in range(self.number_of_qubits):
                qml.RX(2 * params[self.steps + step], wires=i)
        
        return qml.expval(cost_hamiltonian)

    def probability_circuit_definition(self, params):
        cost_hamiltonian = self.create_cost_hamiltonian()
        
        for i in range(self.number_of_qubits):
            qml.Hadamard(wires=i)
        
        for step in range(self.steps):
            qml.ApproxTimeEvolution(cost_hamiltonian, params[step], 1)
            for i in range(self.number_of_qubits):
                qml.RX(2 * params[self.steps + step], wires=i)
        
        return qml.probs(wires=range(self.number_of_qubits))

    def solve_tsp(self):
        init_params = np.random.uniform(0, np.pi, 2 * self.steps)
        
        opt = scipy.optimize.minimize(
            self.circuit,
            init_params,
            method="Nelder-Mead",
            options={"fatol": self.ftol, "xatol": self.xtol}
        )
        
        optimal_params = opt.x
        
        probs = self.probability_circuit(optimal_params)
         #Store the sampling results
        self.sampling_results = {}
        for i, prob in enumerate(probs):
            if prob > 1e-6:  # Only store non-zero probabilities
                binary_state = format(i, f'0{self.number_of_qubits}b')
                self.sampling_results[binary_state] = prob
        
        most_probable_bitstring = max(self.sampling_results, key=self.sampling_results.get)
        solution = self.binary_to_points(list(map(int, most_probable_bitstring)))
        
        return solution, probs

    @staticmethod
    def points_to_binary(points):
        no_points = len(points)
        binary_conv = np.zeros((len(points))**2)
        for i in range(len(points)):
            p = points[i]
            binary_conv[(no_points) * (i) + p[1]] = 1
        return binary_conv

    @staticmethod
    def binary_to_points(binary):
        no_points = int(np.sqrt(len(binary)))
        points = []
        for i in range(no_points):
            for j in range(no_points):
                if binary[i*no_points + j] == 1:
                    points.append(j)
        return points
    
    def create_cost_hamiltonian(self):
        cost_operators = []
        number_of_cities = len(self.distance_matrix)

        for i in range(number_of_cities):
            for j in range(i, number_of_cities):
                for t in range(number_of_cities - 1):
                    weight = -self.distance_matrix[i][j] / 2
                    if self.distance_matrix[i][j] != 0:
                        qubit_1 = t * number_of_cities + i
                        qubit_2 = (t + 1) * number_of_cities + j
                        cost_operators.append(qml.Hamiltonian([weight, -weight], 
                                              [qml.Identity(0), qml.PauliZ(qubit_1) @ qml.PauliZ(qubit_2)]))
                        print(f"City {i} to {j} at t = {t} costs {weight}, Qubits: {qubit_1}, {qubit_2}")

        return qml.sum(cost_operators)