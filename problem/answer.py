import sys
from typing import Any

sys.path.append("../")
from utils.challenge_2024 import ChallengeSampling, problem_hamiltonian

challenge_sampling = ChallengeSampling()

"""
####################################
add codes here
####################################
"""


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
        """
        ####################################
        add codes here
        ####################################
        """

        return 0.


if __name__ == "__main__":
    run_algorithm = RunAlgorithm()
    print(run_algorithm.get_result(seed=0, hamiltonian_directory="../hamiltonian"))
