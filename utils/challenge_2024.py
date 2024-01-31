import juliacall
from juliacall import Main as jl

from typing import TYPE_CHECKING, Any, Optional, Sequence

from collections import Counter
from collections.abc import Collection, Iterable
from openfermion import load_operator
from quri_parts.circuit import QuantumCircuit, QuantumGate, NonParametricQuantumCircuit
from quri_parts.core.estimator.sampling.estimator import _ConstEstimate, _Estimate
from quri_parts.core.estimator import (
    ConcurrentParametricQuantumEstimator,
    ConcurrentQuantumEstimator,
    Estimatable,
    Estimate,
    ParametricQuantumEstimator,
    QuantumEstimator,
    create_parametric_estimator,
)
from quri_parts.core.operator import PAULI_IDENTITY, Operator
from quri_parts.core.measurement import (
    CommutablePauliSetMeasurementFactory,
)
from quri_parts.core.sampling import (
    Sampler,
    ConcurrentSampler,
    PauliSamplingShotsAllocator,
    MeasurementCounts,
)
from quri_parts.core.state import CircuitQuantumState, ParametricCircuitQuantumState
from quri_parts.core.utils.concurrent import execute_concurrently
from quri_parts.itensor.circuit import convert_circuit
from quri_parts.itensor.load_itensor import ensure_itensor_loaded
from quri_parts.itensor.sampler import create_itensor_mps_sampler
from time import time

if TYPE_CHECKING:
    from concurrent.futures import Executor

max_num_shots = 10**7
max_run_time = 6 * 10**5
maxdim, cutoff = 50, 0.0


class ChallengeSampling:
    def __init__(self) -> None:
        self.total_shots: int = 0
        self.total_jobs: int = 0
        self.init_time: float = time()
        # test sampling so that later sampling can be performed quickly
        test_circuit = QuantumCircuit(2)
        test_circuit.add_H_gate(1)
        test_circuit.add_CNOT_gate(0, 1)
        ts_test = time()
        sampler = create_itensor_mps_sampler()
        test_res = sampler(test_circuit, 10)
        print(f"test mps sampling took: {time() - ts_test, test_res}")

    def _sample(
        self, circuit: NonParametricQuantumCircuit, shots: int, **kwargs: Any
    ) -> MeasurementCounts:
        if len(circuit.gates) == 0:
            return Counter({0: shots})
        run_time = time() - self.init_time
        exp_total_shots = self.total_shots + shots
        if exp_total_shots > max_num_shots or run_time > max_run_time:
            raise ExceededError(self.total_shots, exp_total_shots, run_time)
        qubits = circuit.qubit_count
        s: juliacall.VectorValue = jl.siteinds("Qubit", qubits)
        psi_initial: juliacall.AnyValue = jl.init_state(s, qubits)
        circuit_ops = convert_circuit(circuit, s)
        psi = jl.apply(circuit_ops, psi_initial, **kwargs)
        psi = jl.normalize(psi)
        result: list[int] = jl.sampling(psi, shots)
        self.total_shots = exp_total_shots
        self.total_jobs += 1
        return Counter(result)

    def _sample_concurrently(
        self,
        circuit_shots_tuples: Iterable[tuple[NonParametricQuantumCircuit, int]],
        executor: Optional["Executor"],
        concurrency: int = 1,
        **kwargs: Any,
    ) -> Iterable[MeasurementCounts]:
        def _sample_sequentially(
            _: Any,
            _circuit_shots_tuples: Iterable[tuple[NonParametricQuantumCircuit, int]],
        ) -> Iterable[MeasurementCounts]:
            return [
                self._sample(circuit, shots, **kwargs)
                for circuit, shots in _circuit_shots_tuples
            ]

        return execute_concurrently(
            _sample_sequentially, None, circuit_shots_tuples, executor, concurrency
        )

    def create_sampler(self) -> Sampler:
        ensure_itensor_loaded()

        kwargs = {"maxdim": maxdim, "cutoff": cutoff}

        def sampler(
            circuit: NonParametricQuantumCircuit, n_shots: int
        ) -> MeasurementCounts:
            return self._sample(circuit, n_shots, **kwargs)

        return sampler

    def create_concurrent_sampler(
        self,
        executor: Optional["Executor"] = None,
        concurrency: int = 1,
    ) -> ConcurrentSampler:
        ensure_itensor_loaded()
        kwargs = {"maxdim": maxdim, "cutoff": cutoff}

        def sampler(
            circuit_shots_tuples: Iterable[tuple[NonParametricQuantumCircuit, int]]
        ) -> Iterable[MeasurementCounts]:
            return self._sample_concurrently(
                circuit_shots_tuples, executor, concurrency, **kwargs
            )

        return sampler

    # estimator

    def _sample_improve(
        self, circuit: Sequence[NonParametricQuantumCircuit], shots: int, s: juliacall.VectorValue, psi: juliacall.AnyValue, **kwargs: Any
    ) -> MeasurementCounts:
        if len(circuit[0].gates) == 0:
            return Counter({0: shots})
        run_time = time() - self.init_time
        exp_total_shots = self.total_shots + shots
        if exp_total_shots > max_num_shots or run_time > max_run_time:
            raise ExceededError(self.total_shots, exp_total_shots, run_time)
        measurement = circuit[1]
        meas_circ = QuantumCircuit(circuit[0].qubit_count, gates=measurement)
        if len(meas_circ.gates) == 0:
            psi_new = psi
        else:
            meas_circuit_mps = convert_circuit(meas_circ, s)
            psi_new = jl.apply(meas_circuit_mps, psi, **kwargs)
        if any(k in kwargs for k in ["mindim", "maxdim", "cutoff"]):
            psi_new = jl.normalize(psi_new)
        result: list[int] = jl.sampling(psi_new, shots)
        self.total_jobs += 1
        self.total_shots = exp_total_shots
        return Counter(result)

    def _sampling_estimate_concurrently(
        self,
        circuit_shots_tuples: Iterable[tuple[list[NonParametricQuantumCircuit, Sequence[QuantumGate]], int]],
        s: juliacall.VectorValue=None,
        psi: juliacall.AnyValue = None,
        **kwargs: Any,
    ) -> Iterable[MeasurementCounts]:
        def _sample_sequentially(
            _: Any,
            _circuit_shots_tuples: Iterable[tuple[list[NonParametricQuantumCircuit, Sequence[QuantumGate]], int]],
        ) -> Iterable[MeasurementCounts]:
            return [
                self._sample_improve(circuit, shots, s, psi, **kwargs)
                for circuit, shots in _circuit_shots_tuples
            ]

        return execute_concurrently(
            _sample_sequentially, None, circuit_shots_tuples, None, 1
        )

    def sampling_estimate(
        self,
        op: Estimatable,
        state: CircuitQuantumState,
        total_shots: int,
        measurement_factory: CommutablePauliSetMeasurementFactory,
        shots_allocator: PauliSamplingShotsAllocator,
    ) -> Estimate[complex]:
        kwargs = {"maxdim": maxdim, "cutoff": cutoff}

        if not isinstance(op, Operator):
            op = Operator({op: 1.0})

        if len(op) == 0:
            return _ConstEstimate(0.0)

        if len(op) == 1 and PAULI_IDENTITY in op:
            return _ConstEstimate(op[PAULI_IDENTITY])

        # If there is a standalone Identity group then eliminate, else set const 0.
        const: complex = 0.0
        measurements = []
        for m in measurement_factory(op):
            if m.pauli_set == {PAULI_IDENTITY}:
                const = op[PAULI_IDENTITY]
            else:
                measurements.append(m)

        pauli_sets = tuple(m.pauli_set for m in measurements)
        shot_allocs = shots_allocator(op, pauli_sets, total_shots)
        shots_map = {pauli_set: n_shots for pauli_set, n_shots in shot_allocs}

        # First convert the part excluding the measurements circuit to psi

        qubits = state.circuit.qubit_count
        s: juliacall.VectorValue = jl.siteinds("Qubit", qubits)
        psi_initial: juliacall.AnyValue = jl.init_state(s, qubits)
        circuit_ops = convert_circuit(state.circuit, s)
        psi = jl.apply(circuit_ops, psi_initial, **kwargs)

        # Eliminate pauli sets which are allocated no shots
        measurement_circuit_shots = [
            (m, [state.circuit, m.measurement_circuit], shots_map[m.pauli_set])
            for m in measurements
            if shots_map[m.pauli_set] > 0
        ]

        circuit_and_shots = [
            (circuit, shots) for _, circuit, shots in measurement_circuit_shots
        ]

        sampling_counts = self._sampling_estimate_concurrently(
            circuit_and_shots, s, psi, **kwargs
        )

        pauli_sets = tuple(m.pauli_set for m, _, _ in measurement_circuit_shots)
        pauli_recs = tuple(
            m.pauli_reconstructor_factory for m, _, _ in measurement_circuit_shots
        )
        return _Estimate(op, const, pauli_sets, pauli_recs, tuple(sampling_counts))

    def concurrent_sampling_estimate(
        self,
        operators: Collection[Estimatable],
        states: Collection[CircuitQuantumState],
        n_shots: int,
        measurement_factory: CommutablePauliSetMeasurementFactory,
        shots_allocator: PauliSamplingShotsAllocator,
    ) -> Iterable[Estimate[complex]]:
        num_ops = len(operators)
        num_states = len(states)

        if num_ops == 0:
            raise ValueError("No operator specified.")

        if num_states == 0:
            raise ValueError("No state specified.")

        if num_ops > 1 and num_states > 1 and num_ops != num_states:
            raise ValueError(
                f"Number of operators ({num_ops}) does not match"
                f"number of states ({num_states})."
            )

        if num_states == 1:
            states = [next(iter(states))] * num_ops
        if num_ops == 1:
            operators = [next(iter(operators))] * num_states
        return [
            self.sampling_estimate(
                op,
                state,
                n_shots,
                measurement_factory,
                shots_allocator,
            )
            for op, state in zip(operators, states)
        ]

    def create_sampling_estimator(
        self,
        n_shots: int,
        measurement_factory: CommutablePauliSetMeasurementFactory,
        shots_allocator: PauliSamplingShotsAllocator,
    ) -> QuantumEstimator[CircuitQuantumState]:
        def estimator(op: Estimatable, state: CircuitQuantumState) -> Estimate[complex]:
            return self.sampling_estimate(
                op,
                state,
                n_shots,
                measurement_factory,
                shots_allocator,
            )

        return estimator

    def create_sampling_concurrent_estimator(
        self,
        n_shots: int,
        measurement_factory: CommutablePauliSetMeasurementFactory,
        shots_allocator: PauliSamplingShotsAllocator,
    ) -> ConcurrentQuantumEstimator[CircuitQuantumState]:
        def estimator(
            operators: Collection[Estimatable],
            states: Collection[CircuitQuantumState],
        ) -> Iterable[Estimate[complex]]:
            return self.concurrent_sampling_estimate(
                operators,
                states,
                n_shots,
                measurement_factory,
                shots_allocator,
            )

        return estimator

    def create_parametric_sampling_estimator(
        self,
        n_shots: int,
        measurement_factory: CommutablePauliSetMeasurementFactory,
        shots_allocator: PauliSamplingShotsAllocator,
    ) -> ParametricQuantumEstimator[ParametricCircuitQuantumState]:
        sampling_estimator = self.create_sampling_estimator(
            n_shots,
            measurement_factory,
            shots_allocator,
        )
        return create_parametric_estimator(sampling_estimator)

    def create_concurrent_parametric_sampling_estimator(
        self,
        n_shots: int,
        measurement_factory: CommutablePauliSetMeasurementFactory,
        shots_allocator: PauliSamplingShotsAllocator,
    ) -> ConcurrentParametricQuantumEstimator[ParametricCircuitQuantumState]:
        def concurrent_parametric_sampling_estimater(
            operator: Operator,
            state: ParametricCircuitQuantumState,
            params: Sequence[Sequence[float]],
        ) -> Iterable[Estimate[complex]]:
            bind_states = [state.bind_parameters(param) for param in params]
            concurrent_estimator = self.concurrent_sampling_estimate(
                [operator],
                bind_states,
                n_shots,
                measurement_factory,
                shots_allocator,
            )
            return concurrent_estimator

        return concurrent_parametric_sampling_estimater

    def reset(self):
        self.total_shots = 0
        self.total_jobs = 0
        self.init_time = time()


def problem_hamiltonian(n_qubits: int, seed: int, hamiltonian_directory: str):
    print(f"{str(n_qubits).zfill(2)}qubits_{str(seed).zfill(2)}")
    ham = load_operator(
        file_name=f"{str(n_qubits).zfill(2)}qubits_{str(seed).zfill(2)}",
        data_directory=hamiltonian_directory,
        plain_text=False,
    )
    return ham


class ExceededError(Exception):
    def __init__(self, shots_used, shots_exp: float, run_time: float):
        self.shots_used = shots_used
        self.shots_exp = shots_exp
        self.run_time = run_time

    def __str__(self) -> str:
        if self.run_time > max_run_time:
            return (
                f"Reached maximum runtime: {max_run_time}. " f"Run time {self.run_time}"
            )
        else:
            return (
                f"Reached maximum number of shots: {max_num_shots}. "
                f"Number of shots used: {self.shots_used}. "
                f"Job terminated before reaching {self.shots_exp} shots. "
            )
