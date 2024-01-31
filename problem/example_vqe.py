import sys
sys.path.append("../")

from typing import Any
import numpy as np
from openfermion.transforms import jordan_wigner

from quri_parts.algo.ansatz import SymmetryPreservingReal
from quri_parts.algo.optimizer import SPSA, OptimizerStatus
from quri_parts.circuit import LinearMappedUnboundParametricQuantumCircuit
from quri_parts.core.estimator.gradient import parameter_shift_gradient_estimates
from quri_parts.core.measurement import bitwise_commuting_pauli_measurement, CachedMeasurementFactory
from quri_parts.core.sampling.shots_allocator import (
    create_equipartition_shots_allocator,
)
from quri_parts.core.state import ParametricCircuitQuantumState, ComputationalBasisState
from quri_parts.openfermion.operator import operator_from_openfermion_op


from utils.challenge_2024 import ChallengeSampling, ExceededError, problem_hamiltonian

challenge_sampling = ChallengeSampling()


def cost_fn(hamiltonian, parametric_state, param_values, estimator):
    estimate = estimator(hamiltonian, parametric_state, [param_values])
    return estimate[0].value.real


def vqe(hamiltonian, parametric_state, estimator, init_params, optimizer):
    opt_state = optimizer.get_init_state(init_params)
    energy_history = []

    def c_fn(param_values):
        return cost_fn(hamiltonian, parametric_state, param_values, estimator)

    def g_fn(param_values):
        grad = parameter_shift_gradient_estimates(
            hamiltonian, parametric_state, param_values, estimator
        )
        return np.asarray([i.real for i in grad.values])

    while True:
        try:
            opt_state = optimizer.step(opt_state, c_fn, g_fn)
            energy_history.append(opt_state.cost)
        except ExceededError as e:
            print(str(e))
            print(opt_state.cost)
            return opt_state, energy_history

        if opt_state.status == OptimizerStatus.FAILED:
            print("Optimizer failed")
            break
        if opt_state.status == OptimizerStatus.CONVERGED:
            print("Optimizer converged")
            break
    return opt_state, energy_history


class RunAlgorithm:
    def __init__(self) -> None:
        challenge_sampling.reset()

    def result_for_evaluation(self, seed: int, hamiltonian_directory: str) -> tuple[Any, float]:
        energy_final = self.get_result(seed, hamiltonian_directory)
        total_shots = challenge_sampling.total_shots
        return energy_final, total_shots

    def get_result(self, seed: int, hamiltonian_directory: str) -> float:
        """
            param seed: the last letter in the Hamiltonian data file, taking one of the values 0,1,2,3,4
            param hamiltonian_directory: directory where hamiltonian data file exists
            return: calculated energy.
        """
        n_qubits = 28
        ham = problem_hamiltonian(n_qubits, seed, hamiltonian_directory)
        n_site = n_qubits // 2
        total_shots = 10**6
        jw_hamiltonian = jordan_wigner(ham)

        hamiltonian = operator_from_openfermion_op(jw_hamiltonian)

        # make hf + SPReal ansatz
        hf_gates = ComputationalBasisState(n_qubits, bits=2**n_site - 1).circuit.gates
        hf_circuit = LinearMappedUnboundParametricQuantumCircuit(n_qubits).combine(
            hf_gates
        )
        ansatz = SymmetryPreservingReal(qubit_count=n_qubits, reps=n_qubits)
        hf_circuit.extend(ansatz)
        shots_allocator = create_equipartition_shots_allocator()
        cached_measurement_factory = CachedMeasurementFactory(
            bitwise_commuting_pauli_measurement
        )

        parametric_state = ParametricCircuitQuantumState(n_qubits, hf_circuit)

        sampling_estimator = (
            challenge_sampling.create_concurrent_parametric_sampling_estimator(
                total_shots, cached_measurement_factory, shots_allocator
            )
        )

        optimizer = SPSA(ftol=10e-5)

        init_param = np.random.rand(ansatz.parameter_count) * 2 * np.pi * 0.001

        result, energy_history = vqe(
            hamiltonian,
            parametric_state,
            sampling_estimator,
            init_param,
            optimizer,
        )
        print(f"iteration used: {result.niter}")
        return min(energy_history)


if __name__ == "__main__":
    run_algorithm = RunAlgorithm()
    print(run_algorithm.get_result(seed=0, hamiltonian_directory="../hamiltonian"))
