# Modular Gauss Sum Analysis Tools

This repository contains programs developed to investigate a conjecture about modular counter examples to the classical converse theorem of Gauss sums. While the classical theorem (over ℂ) uniquely determines multiplicative characters without counter examples, our work focuses on the modular setting (over the algebraic closure of Fℓ, Fℓ̄) where counter examples do occur. According to the conjecture by Bakeberg, Gerbelli-Gauthier, Goodson, Iyengar, Moss, and Zhang, such counter examples occur precisely when:

q = 2ℓ^i + 1 (for some integer i > 0)

## Tools

- **pair_generator.py**  
  Identifies candidate (ℓ, q) pairs that follow the conjectured form.

- **conjecture.sage**  
  Constructs Gauss sum tables over Fℓ̄ using candidate pairs from `pair_generator.py` and checks if they produce a theta group of size greater than 2 (indicating a counter example).

- **converse_theorem_info.sage**  
  Investigates specific (ℓ, q) pairs by computing Gauss sums modulo ℓ to verify whether they serve as counter examples, regardless of their form.

## Table of Contents
- [Tools](#tools)
- [Getting Started](#getting-started)
- [Usage](#usage)
- [Contributing](#contributing)
- [Contact](#contact)

## Getting Started

### Prerequisites

- **SageMath**: Required for running the `.sage` scripts.
- **Python 3**: Required for running `pair_generator.py`.
- Additional dependencies as needed.

### Installation

1. **Clone the repository**:
`git clone https://github.com/jamwevan/MATH-440.git`  
`cd MATH-440/Code`

2. **Set up your environment**:  
Ensure SageMath is installed and accessible via the command line. Optionally, create a Python virtual environment:
`python3 -m venv venv`  
`source venv/bin/activate`  
`pip install -r requirements.txt`

## Usage

1. **Run the Pair Generator**:  
`python pair_generator.py`  
This script scans for candidate (ℓ, q) pairs that follow the conjectured form and outputs the maximum candidate values (l_max and q_max). Note these values.

2. **Test Candidate Pairs with conjecture.sage**:  
`sage conjecture.sage`  
When prompted, enter the l_max and q_max values from the pair generator. The script builds Gauss sum tables over Fℓ̄ and checks if the candidate pairs yield a counter example (e.g., a theta group larger than 2).

3. **Test Specific (ℓ, q) Pairs with converse_theorem_info.sage**:  
`sage converse_theorem_info.sage`  
Follow the prompts to enter any (ℓ, q) values. This script computes Gauss sums modulo ℓ and verifies if the pair is a counter example.

## Contributors

James Evans, Xinning Ma, & Yanshun Zhang

## Contact

**James Evans**  
Email: [jamwevan@umich.edu](mailto:jamwevan@umich.edu)
