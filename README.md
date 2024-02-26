# Quantum Algorithm Grand Challenge

# Table of Contents
1. [Overview of Quantum Algorithm Grand Challenge](#Overview)
2. [Introduction](#Introduction)
    - [Background](#Introduction_1)
    - [Model description](#Introduction_2)
    - [Simulator](#Introduction_3)
3. [Problem description](#problem)
    - [Fermi-Hubbard Model](#problem_1)
    - [Problem statement](#problem_2)
4. [Evaluation](#Evaluation)
5. [Implementation](#Implementation)
6. [How to submit](#Submission)
7. [Description of the Provided Program](#program)
8. [Available Packages](#Packages)
9. [Notes and Prohibited Items](#forbidden)
    - [Notes on Evaluation](#forbidden_1)
    - [Prohibited Items](#forbidden_2)
10. [Terms](#Terms)


# Overview of Quantum Algorithm Grand Challenge <a id="Overview"></a>
Quantum Algorithm Grand Challenge (QAGC) is a global online contest for students, researchers, and others who learn quantum computation and quantum chemistry around the world.

Through this challenge, we aim to explore practical uses of Noisy Intermediate-Scale Quantum (NISQ) devices, visualize bottlenecks in NISQ device utilization, and create metrics to benchmark NISQ algorithms.

## Date
From February 1, 2024, 0:00 (JST) to June 30, 2024, 23:59 (JST)

## Web-site
 https://www.qagc.org/

## Awards

- 1st Place - $5,000
- 2nd Place - $2,500
- 3rd Place - $1,500

Currently, QunaSys is planning to apply to hold a workshop to be hosted by QunaSys at IEEE Quantum Week 2024. IEEE Quantum Week 2024 will be held as an in-person event with virtual participation on Sep 15–20 at Palais des Congrès Montréal, Québec, Canada. As in previous years, if applications are approved, the top three teams will be invited to present their algorithms at this workshop. We plan to share the details with participants as soon as they are finalized.

For more information on IEEE Quantum Week 2024, please visit the website https://qce.quantum.ieee.org/2024/.


# Introduction <a id="Introduction"></a>

As Quantum Computing technology evolves with qubit capacity regularly duplicating, we need to understand how to make better use of NISQ devices and create algorithms that will enable industrial applications. To identify how to shape the direction for promoting the NISQ algorithm for practical industrial application, it is important to clarify the evaluation criteria to compare the algorithm's performance and define the key factors to take into account. 

We hold a global online contest, the QAGC, to explore practical uses for NISQ devices, visualize bottlenecks in NISQ device utilization, and create a metric for benchmarking the NISQ algorithms.

## Background <a id="Introduction_1"></a>

The materials around us are constructed from molecules and the microscopic behavior is dominated by quantum mechanics. Quantum chemistry is widely used to understand the chemical behavior of these materials in not only academic studies but also material design in industries.

Quantum chemistry is considered one of the most promising fields for considering practical industrial applications of NISQ devices and algorithms.  However, although various NISQ algorithms have been proposed in recent years, it is still far from practical industrial applications. 

For practical industrial applications of NISQ algorithms, it is important to develop new useful NISQ algorithms and define evaluation criteria for accurately comparing the performance of various algorithms and define the key factors to take into account.

Based on these situations, the focuses of QAGC are on the industrial application and defining evaluation criteria for appropriate performance comparison of NISQ algorithms. We have prepared a simulator that reflects the features of NISQ devices and a suitable model for the problem to achieve these goals. Below, we will explain each of them.


## Model description <a id="Introduction_2"></a> 

The ground-state energy of a molecule is an important quantity for understanding its properties and behavior, and many quantum chemistry studies focus on the ground-state energy of individual atoms or molecules. 

In QAGC, the task of participants is to calculate the ground-state energy of a model (Hamiltonian) which we have prepared. From the focus of QAGC, the Hamiltonian should have some properties as follows:

- Has similar properties as the molecular Hamiltonian used in quantum chemistry.

  - The number of terms of the Hamiltonian is $O(N^4)$ , which is the same as the molecular Hamiltonian. Then we can compare the performance of grouping methods that reduce the number of measurements.
  - This Hamiltonian has the same order of operator norm as the molecular Hamiltonian. Therefore, the resulting ground-state energy scale is similar to the scale in quantum chemistry.


- The exact value of ground-state energy of this Hamiltonian can be calculated classically for the arbitrary size of the system.
  
  - Our aim through QAGC is also to create a common metric that can evaluate various NISQ algorithms. For evaluating algorithms in large qubit systems that cannot be simulated in classical computers, it will be necessary to know the exact value of the quantity to be measured as a reference value for evaluating NISQ algorithms. 

We have prepared a Hamiltonian that satisfies all of these properties. The detail of the Hamiltonian and the problem statement in QAGC is written in [Problem](#problem).

## Simulator <a id="Introduction_3"></a>

In order to explore the practical application of NISQ devices and to visualize the bottlenecks in their use, it is necessary to perform simulations that reflect the features of NISQ devices in systems with enough qubits to address practical problems.

In QAGC, the participants need to use a sampling simulator we have provided. In this sampling simulator, quantum circuits are internally converted to Matrix Product State (MPS) and sampled using the MPS simulator, which is faster than ordinary simulators. This enables sampling simulation at 28 qubits, which is difficult to achieve with ordinary simulators. This simulator is used by calling MPS simulator in [ITensor](https://www.scipost.org/SciPostPhysCodeb.4) through QURI Parts.

In addition to this, noise is an important feature of NISQ devices, and it is necessary to consider how to deal with noise when considering the practical use of NISQ devices. Therefore, in this simulator, the approximation error in converting the quantum circuit to MPS is regarded as noise, and a simulation with noise at 28 qubits is performed. It should be noted that the noise of this simulator is essentially different from the noise of the NISQ device due to the above-mentioned origin.

Also, when sampling within the algorithm, the total number of shots is limited to no more than  $10^7$.

# Problem description <a id="problem"></a>

## Fermi-Hubbard Model <a id="problem_1"></a>
Below is an explanation of the Hamiltonian used in QAGC. For further details, please refer to the following paper.
[https://arxiv.org/abs/2402.11869](https://arxiv.org/abs/2402.11869)

The Fermi-Hubbard model is a model used to describe the properties of strongly correlated electron systems, which are solids with strong electron correlation effects. It is used to explain important physical phenomena such as magnetism, Mott insulators, and high-temperature superconductors. 


In QAGC, we deal with a one-dimensional orbital rotated Fermi-Hubbard model with **periodic boundary conditions**. The Hamiltonian of one-dimensional Fermi-Hubbard model is as follows:

$$
    H = - t \sum_{i=0}^{N-1} \sum_{\sigma=\uparrow, \downarrow} (a^\dagger_{i, \sigma}  a_{i+1, \sigma} +  a^\dagger_{i+1, \sigma}  a_{i, \sigma})  - \mu \sum_{i=0}^{N-1} \sum_{\sigma=\uparrow, \downarrow}  a^\dagger_{i, \sigma} a_{i, \sigma} + U \sum_{i=0}^{N-1} a^\dagger_{i, \uparrow}  a_{i, \uparrow}  a^\dagger_{i, \downarrow} a_{i, \downarrow},
$$

where $t$ is the tunneling amplitude, $\mu$ is the chemical potential, and $U$ is the Coulomb potential. For the case of half-filling, i.e. the number of electrons is equal to the number of sites, the exact value of the ground-state energy for this Hamiltonian can be calculated by using Bethe Ansatz method. 

This time we consider the orbital rotated one-dimensional Fermi-Hubbard model. The orbital rotation means the linear transformation of the creation operator $a_i^\dagger$ and annihilation operator $a_i$ by using unitary matrices

$$
    \tilde c_i^\dagger = \sum_{k=0}^{2N-1} u_{ik} c_k^\dagger, \quad 
    \tilde c_i = \sum_{k=0}^{2N-1} u_{ik}^* c_k.
$$

where we label the creation operator $a_{i, \sigma}^\dagger$ as follows:

$$
    a_{i, \uparrow}^\dagger = c_{2i}^\dagger, \quad 
    a_{i, \downarrow}^\dagger = c_{2i + 1}^\dagger.
$$

The annihilator operator is labeled in the same way.

By performing this spin-involved orbital rotation, the 1D FH Hamiltonian can be expressed as

$$
    \tilde{H} = \sum_{p, q =0}^{2N-1}h_{pq} c^\dagger_p c_q + \frac{1}{2}\sum_{p, q, r, s =0}^{2N-1}h_{pqrs} c^\dagger_p c_qc^\dagger_r c_s.
$$

Here, the coefficients $h_{pq}$, $h_{pqrs}$ are the one-body and two-body coefficients, respectively, and they are expressed as follows:

$$
    h_{pq} = -\sum_{i =0}^{N-1}(u_{2i,p} u_{2i+2,q}^* + u_{2i+2,p} u_{2i,q}^* + u_{2i+1,p} u_{2i+3,q}^* + u_{2i+3,p} u_{2i+1,q}^* + \mu \delta_{pq}),
$$

$$
    h_{pqrs} = 2 U \sum_{i=0}^{N-1} u_{2i,p} u_{2i,q}^* u_{2i+1,r} u_{2i+1,s}^*.
$$

This form is the same as the second quantized molecular Hamiltonian in quantum chemistry and the number of terms is $\mathcal{O}(N^4)$ which is the same as the molecular Hamiltonian. 

After performing orbital rotation, the Hartree-Fock calculation can be performed similarly to the molecular Hamiltonian. The resulting Hartree-Fock state becomes:

$$
    |HF\rangle = |00\dots 0011\dots 11\rangle
$$

where electrons are filled from the bottom up for a number of sites.

For QAGC, we will provide the Hamiltonian obtained by orbital rotation using a real orthogonal matrix. By using a real orthogonal matrix, the Hamiltonian becomes real and it become more similar to the actual molecular Hamiltonian.

## Problem statement <a id="problem_2"></a>

Find the energy of the ground state of the one-dimensional orbital rotated Fermi-Hubbard model.

$$
    \tilde{H} = - t \sum_{i=0}^{2N-1}(\tilde c^\dagger_i \tilde c_{i+1} + \tilde c^\dagger_{i+1} \tilde c_i)  - \mu \sum_{i=0}^{2N-1}  \tilde c^\dagger_i \tilde c_i + U \sum_{i=0}^{N-1} \tilde c^\dagger_{2i} \tilde c_{2i} \tilde c^\dagger_{2i + 1} \tilde c_{2i + 1} 
$$

The value of each parameter is $N = 14,\ t=1, U=3,\ \mu=U/2 = 1.5$. 

For QAGC, we prepared an orbital rotated Hamiltonian with the random real orthogonal matrix $u$ and performed Hartree-Fock calculation. 

In the `hamiltonian` folder, in addition to the 28 qubit Hamiltonian, we also prepared 4, 12, and 20 qubit Hamiltonians in the format `.data`. Participants can use these Hamiltonians in implementing their algorithms.


# Evaluation<a id="Evaluation"></a>

First, the submitted answers are checked for compliance with the prohibited items. Then, we calculate the score based on the answers, and the ranking is determined. 

During the evaluation, Hamiltonians constructed using a real orthogonal matrix different from the one used to construct the Hamiltonian provided to the participants will be used.


## Score

The score $S$ is calculated as the average accuracy computed from each absolute error over 10 runs of the algorithm, rounded to the nearest $10^{-8}$. The smaller the score, the higher the ranking you can achieve.

For each run, the absolute error is calculated as follows:

$$
    e_i = \left|E_i - E_{\mathrm{exact}}\right|.
$$

Here $E_i$ is the output result of the $i$ th algorithm and $E_{\mathrm{exact}}$ is the exact value of the ground-state energy.

The score is determined as follows:

$$
    S = \frac{1}{10}\sum_{i=1}^{10}e_i
$$

## Limitation by the number of shots

Reducing the number of shots is crucial for considering the industrial application of NISQ algorithms. To reflect this, participants will be imposed a limit on the number of shots. 

For the QAGC, the total number of shots is limited to $10^7$.

## Limitation by Run Time During the Evaluation

During the evaluation period by the management, the submitted algorithm must be completed within one week  ( $6\times10^5$ sec) on a computer with the following specifications. If the evaluation period exceeds this limit and is not completed, it will be forcibly stopped and the score at that time will be the final score.

- PC: MacBook Pro (13-inch, M1, 2020)
- Processor: Apple M1 (8-core)
- Memory: 16 GB
- Storage: 256 GB

## Maximum Usable Qubits

In QAGC, the number of the maximum usable qubits has been set to 55 to prevent a substantial increase in the number of shots resulting from running two same 28-qubit circuits in parallel. Within this limit, you are free to add more qubits, but the calculation must be completed within the time limit.

# Implementation <a id="Implementation"></a>

Here, we will explain the necessary items for participants to implement their answer code.

Participants need to write their algorithms in `answer.py`.
- Participants should write their code in `get_result()` in `RunAlgorithm`. 

- It is also possible to add functions outside of RunAlgorithm as needed.

- The only codes that participants can modify are those in the problem folder. Do not modify the codes in the utils folder.

We have prepared an answer example in `example_***.py`, so please refer to it. 

Below, we will explain the sampling function and how to use the Hamiltonian of the problem.

-  ## Sampling Function

    In QAGC, all participants need to use the sampling function we have provided. Please refer to the `sampling.ipynb` in the `tutorials` folder for instructions on how to use it.

    This sampling function has the property that when the expected shot limit is reached, the error **ExceededError** will be output.

-  ## Hamiltonian

    The orbital rotated Fermi-Hubbard Hamiltonian is stored in the `hamiltonian` folder in `.data` format. To load it, use `problem_hamiltonian()` as follows:
    ``` python
    from utils.challenge_2024 import problem_hamiltonian
    
    n_qubits = 28

    ham = problem_hamiltonian(n_qubits, seed, hamiltonian_directory)
    ```

    - The parameter `seed` is a parameter that specifies which Hamiltonian to use by selecting from the different Hamiltonians in the `hamiltonian` folder. There are five Hamiltonians available for each qubit, and you can specify which Hamiltonian to use by selecting one seed from `[0, 1, 2, 3, 4]`.

    - The parameter `hamiltonian_directory` in `problem_hamiltonian` must be left unchanged.


    In the `hamiltonian` folder, in addition to the 28 qubit Hamiltonian, we also prepared 4, 12, and 20 qubit Hamiltonians in the format `.data`.  
    Participants can use these Hamiltonians in implementing their algorithms.
    
    The important point is that during the evaluation, we will use Hamiltonians constructed using different real orthogonal matrices than the one used to construct the Hamiltonian in the `hamiltonian` folder.

Participants can calculate the score by running `evaluator.py`.
  - **seeds_list**: This parameter specifies which Hamiltonian is used in the evaluation. For example, if you want to run the algorithm three times with a Hamiltonian of 0 to get a score, please set `seeds_list = [0, 0, 0]`.
  - **ref_value**: The reference value (exact value of the ground-state energy) for each Hamiltonian is listed. The score is evaluated based on this value.

Since we are dealing with a large qubits system such as 28 qubits, running evaluator.py using the code in example_***.py takes one or two days for a single execution.

# How to submit <a id="Submission"></a>

The participants's code will be submitted as an issue using this template summarizing your project. Specifically, this issue should contain:

1. Team name: Your team's name
2. Team members: Listup all member's name
3. Project Description: A brief description of your project (1-2 paragraphs).
4. Presentation: A link of the presentation with slides of your team’s hackathon project.
5. Source code: A link to the final source code for your team's hackathon project (e.g., a GitHub repo).

The score will be calculated by the management side, and the rankings will be determined and published on the [QAGC web site](https://www.qagc.org/).

- Participants can submit their work as many times as they want during the period of QAGC.

- Participants can form teams with members.

- Submitted code is evaluated once a week and rankings are presented along with scores.

Here are some points to note when submitting.

- The participants's code can be viewed by other participants.

- If for some reason you do not wish to make the code public, please send the reason and the code directly to the management (qagc@qunasys.com). If this reason is accepted, your score will be calculated and your ranking determined without publishing your code.

# Description of the Provided Program <a id="program"></a>

We have provided some codes for QAGC. The descriptions of each code are as follows.

  - `README.md`:
  
    This is the explanation of the QAGC.

  - `tutorials`:

    This contains some tutorials.

  - `hamiltonian`:
  
    The Hamiltonian to be used in the problem is stored in this folder in `.data` format.

The code in `problem` is structured as follows

  - `answer.py`:

    This is the file for implementing the participants's code. See [How to Submit an Algorithm](#Submission) for details.

  - `evaluator.py`:
    
    This is the code to evaluate the answer and calculate the score.

  - `example_***.py`:
    
    These are example codes of the answer prepared by QunaSys. 
    - `example_vqe.py`

      This is the code that implements [VQE](https://www.nature.com/articles/ncomms5213) as an example.
      

    - `example_adaptqsci.py`
    
      This is the code that implements [ADAPT-QSCI](https://arxiv.org/abs/2311.01105) as an example. The score of QunaSys in the scoreboard is calculated using this code.
      Please refer to the [Terms](#Terms) regarding the handling of this code.

The code in `utils` is structured as follows.

  - `challenge_2024.py`:
    
    This contains the sampling function used in QAGC. 


# Available Packages <a id="Packages"></a>

The following Python software library can be used in QAGC.

- [QURI Parts](https://quri-parts.qunasys.com/)

- [Qiskit](https://qiskit.org/)

- [Cirq](https://quantumai.google/cirq)

- [Amazon Braket Python SDK](https://amazon-braket-sdk-python.readthedocs.io/en/latest/#)

**QURI Parts** is an open-source quantum computing library that is modular, efficient, and platform-independent, developed by QunaSys.

- Platform-independent: Run one algorithm code on various simulators and platforms.

- Modularity and Scalability: Combine parts to create your own algorithm, and easily create and use your own parts.

- High-speed: Classical processing and simulator calls associated with quantum computing are efficient. It is the fastest platform-independent library using Qulacs.

- Open source: Released under Apache License 2.0.

All codes we have prepared are written by using **QURI Parts**.

In QURI Parts, there are codes to convert circuits of **Braket** and circuits and operators of **Qiskit** and **Cirq** to **QURI Parts** circuits and operators. Therefore, when implementing with **Braket**, **Qiskit**, and **Cirq**, these conversion codes can be used to take advantage of the sampling features provided.

Below is an example of code that converts the circuit and operator of **Cirq** to the circuit and operator of **QURI Parts**. We also have provided some example codes of these conversion codes.

```python
from quri_parts.cirq.circuit import circuit_from_cirq
from quri_parts.cirq.operator import operator_from_cirq_op

quri_parts_circuit = circuit_from_cirq(cirq_circuit)
quri_parts_operator = operator_from_cirq_op(cirq_operator)
```

## Version

The version of the main package used in the challenge for participants will be fixed as follows:

```
python >= 3.9.8, <=3.11

quri-parts == 0.16.1
qiskit == 0.41.1
cirq == 1.1.0
amazon-braket-sdk = 1.66.0
openfermion == 1.5.1
qulacs == 0.5.6
numpy == 1.23.5
```

If you use a version other than the specified one, or use other packages, please specify the name of that package and its version in the issue to be registered when submitting.

# Notes and Prohibited Items <a id="forbidden"></a>

## Notes on Evaluation <a id="forbidden_1"></a>

The validity of the final answer will be judged by the judge based on whether it falls under the prohibited answers below. If it is deemed valid, a score will be calculated. The final decision on the validity of the answer and the score will be made by the operator.

## Prohibited Items <a id="forbidden_2"></a>

- Answers that do not essentially use quantum computers.
  
  In QAGC, we prohibit answers that do not essentially use quantum computers, such as the following examples.

  - Example 1: Pushing the exponentially time-consuming parts onto classical computation
    - Calculate the wave function classically and only calculate the final expectation value with the quantum algorithm.
  - Example 2: Not using the quantum algorithm at all.
    - Algorithms that are all classically computable are regarded as essentially not using a quantum computer.

- Algorithms that use explicitly obtained values from classical computation.

  - Example: Algorithms that use an exact value of the ground-state energy of the Hamiltonian which is obtained by using the Bethe ansatz.

- The implemented algorithm must be scalable with respect to the number of qubits.
  
  - For example, if you consider diagonalizing a truncated Hamiltonian, the truncation must be selected to be scalable.

- Hard-coded *good* parameter sets.

    - The selection of initial parameters or hyper parameters must be done in a scalable manner. You cannot hard-code a *good* parameter set into the answer code.
      - Example 1: Non-trivial parameters determined without any basis in physical or mathematical theory.

      - Example 2: Knowing in advance the seed value that will give good results and specifying it in the algorithm.

- Answers that output values that are not calculated using the algorithm.

- Modifying code that is not allowed to be modified.

  - The only codes that participants can modify are those in the `problem` folder. Do not modify the codes in the utils folder.

- Codes with bugs are invalid.

# Terms <a id="Terms"></a>

Conditions of Participation in the Quantum Algorithm Grand Challenge

I, or our company (the “Participant”), agree to the following conditions of participation (the "Terms") and will participate in the Quantum Algorithm Grand Challenge (“QAGC”) conducted or operated by QunaSys Inc. (“QunaSys”). If any of our employees participate in the QAGC, they will also comply with these Terms.

1.	The purpose of the QAGC is to engage the Participant in practical problem-solving learning by collaborating with themselves or other participants and utilizing the challenges, programs, or data ("Challenge Data,") provided by QunaSys.

2.	The Participants are expected to analyze the Challenge Data, create responses to the challenges, and develop or modify programs.

3.	All intellectual property rights arising from the challenge data provided by QunaSys belong exclusively to QunaSys.

4.	The intellectual property rights to the results created or generated by the Participant using the Challenge Data (The "Results") belong to the respective the Participant. The Results include, but are not limited to, new ideas, responses to challenges, and programs.

5.	The Participants are required to submit the Results to QunaSys by the end of the QAGC.

6.	QunaSys will not use, exploit, or implement the Results beyond the scope of considering awards in the QAGC or the purpose of operating QAGC.

7.	Unless the Participants explicitly refuse in advance, the Results will be made publicly available via Github.

8.	The Participant agrees not to use the Challenge Data for ADAPT_QSCI provided by QunaSys, including but not limited to `example_adaptqsci.py`, for purposes other than those set forth herein. Additionally, the Challenge Data contains patents owned by QunaSys, and these patents are authorized for use by the Participants solely to achieve the purpose stated herein and are not permitted for any other purposes.

9.	QunaSys will award a prize to the Participants who have been selected as winners based on the evaluation of the Results.

10.	The Participant must comply with all laws, regulations, and public order and morals and must not infringe upon any third-party intellectual property rights or any other rights in participating in the QAGC.

11.	The Participant shall resolve any disputes arising from the QAGC independently and shall not seek compensation or indemnification from QunaSys.

12.	If the Participant violates any provisions of these Terms and causes damage to QunaSys or other participants, they shall be liable to compensate for such damages.

13.	If the Participant is a legal entity, the responsibility for any violations of these Terms by employees who actually participate in the project will be borne by that legal entity.