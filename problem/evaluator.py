import numpy as np
import sys
import traceback

from typing import Optional

from example_adaptqsci import RunAlgorithm

###################################### Edit here for Calculating Accuracy ######################################
seeds_list = [0]  # list of hamiltonian seeds to be used (e.g. [0, 1, 2, 3, 4], ...)
# For example, if you want to run the algorithm three times with a Hamiltonian of 0 to get a score, please set seeds_list = [0, 0, 0].
ref_value = -30.748822808
################################################################################################################


"""
Available Parameter seed
[0, 1, 2, 3, 4]

reference values (n_qubits: reference_value)
4: -4.00000000, 
12: -13.433353608,
20: -22.046059902,
28: -30.748822808,
"""


############################################ Do NOT Edit Below ############################################
class EvaluateResults:
    def __init__(self) -> None:
        self.n_shots_history: list[float] = []
        self.result_history: list[float] = []
        self.accuracy: Optional[float] = None
        self.num_run: int = 0

    def get_accuracy(
            self, reference_value=ref_value, seeds=None, hamiltonian_directory="../hamiltonian",
    ) -> Optional[float]:
        """
        :return: Get the accuracy of the algorithm results.
        """
        if seeds is None:
            seeds = seeds_list
        for seed in seeds:
            self.num_run += 1
            run_algorithm = RunAlgorithm()
            print(f"Running algorithm({self.num_run})..")
            try:
                energy, n_shots = run_algorithm.result_for_evaluation(seed, hamiltonian_directory)
                acc = abs(reference_value - energy)
                print("\n############## Result ##############")
                print(f"Resulting Energy = {energy}")
                print(f"Number of Shots = {n_shots}")
                print("####################################\n")
                self.n_shots_history.append(n_shots)
                self.result_history.append(acc)
            except Exception as e:
                type_, value_, traceback_ = sys.exc_info()
                print(f"{type(e)}: {str(e)}")
                traceback_message = traceback.format_exception(
                    type_, value_, traceback_
                )
                for l_ in traceback_message:
                    print(l_.strip("\n"))
                return None
        result_ave = np.average(self.result_history)
        self.accuracy = result_ave
        print("\n############## Final Result ##############")
        print(f"Average accuracy = {np.round(self.accuracy, 8)}")
        print(f"reference energy = {reference_value}")
        print("##########################################")
        return self.accuracy


if __name__ == "__main__":
    evaluator = EvaluateResults()
    evaluator.get_accuracy(seeds=seeds_list)
    print(f"Accuracy: {evaluator.accuracy}")
    print(f"Energy result history: {evaluator.result_history}")
    print("Finished")
