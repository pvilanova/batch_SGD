### Adaptive Batch Size Selection for SGD

#### Overview
This repository contains the MATLAB implementation of the methods discussed in the paper "On the Equivalence of Different Adaptive Batch Size Selection Strategies for Stochastic Gradient Descent Methods" by Luis Espath, Sebastian Krumscheid, Raul Tempone, and Pedro Vilanova. The paper presents a theoretical analysis and numerical experiments demonstrating the equivalence of two batch size selection strategies for Stochastic Gradient Descent (SGD) methods: the norm test and the inner product/orthogonality test.

#### Key Concepts
- **Statistical Error Decomposition**: The statistical error in estimating the expected gradient is decomposed, leading to two distinct batch size selection tests: the norm test and the inner/orthogonality test.
- **Adaptive Batch Size Strategies**: The paper shows that these two tests are equivalent in terms of convergence rates under specific conditions, particularly when the squared batch size parameters (\( \epsilon^2 \)) equal the sum of the squared errors in the parallel (\( \theta^2 \)) and orthogonal (\( \nu^2 \)) directions to the exact gradient.
- **Numerical Experiments**: The implementation includes numerical experiments to validate the theoretical analysis, comparing the efficiency and effectiveness of both batch size selection strategies in SGD optimization.

#### Implementation Details
The MATLAB code provided in this repository implements the SGD algorithm with adaptive batch size selection based on the norm and inner/orth tests. The key components of the code include:
- Calculation of the Hessian matrix and its adaptation over iterations.
- Definition of the objective function and its gradient.
- Adaptive batch size calculation using both the norm and inner/orth tests.
- Iterative optimization process to minimize the given function using the calculated batch sizes.

#### Usage
The code is structured to allow easy experimentation with different parameters and functions. Users can modify the objective function, gradient definitions, and batch size parameters to explore the behavior of the SGD algorithm under various conditions.

#### Conclusion
This implementation provides a practical exploration of the theoretical findings in the paper, offering insights into the efficiency of adaptive batch size selection strategies in stochastic optimization. By comparing the norm and inner/orth tests, users can gain a deeper understanding of their equivalence and applicability in different optimization scenarios.
