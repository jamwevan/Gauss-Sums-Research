# Import necessary functions from SageMath
from sage.all import is_prime, carmichael_lambda, GF

# -----------------------------------------------
# Class to generate a Gauss Sum Table
# This table computes sums of additive and multiplicative characters
# over a finite field, which are used in number theory and cryptography.
# -----------------------------------------------
class GaussSumTable:
    def __init__(self, q, additive_character_generator, multiplicative_character_generator):
        # Store field size parameter q
        self.q = q
        # Store the additive and multiplicative character generators
        self.additive_character_generator = additive_character_generator
        self.multiplicative_character_generator = multiplicative_character_generator
        # Define a finite field GF(q^2) (a quadratic extension of GF(q))
        self.finite_field = GF(q^2)
        # Get a generator for the field
        self.generator = self.finite_field.gen()
        # Get all elements of the finite field
        self.finite_field_elements = list(self.finite_field)
        # Extract the multiplicative group (excluding zero)
        self.finite_field_multiplicative_group = [x for x in self.finite_field_elements if x != 0]
        # Initialize a table of size (q^2 - 1) x (q - 1) filled with zeros
        self.table = [[0 for _ in range(q - 1)] for _ in range(q^2 - 1)]
        # Compute the Gauss sum table upon initialization
        self.compute_gauss_sum_table()

    # -----------------------------------------------
    # Compute the full Gauss sum table
    # Iterates over all valid (theta, alpha) pairs
    # -----------------------------------------------
    def compute_gauss_sum_table(self):
        for theta in range(self.q^2 - 1):
            for alpha in range(self.q - 1):
                self.table[theta][alpha] = self.compute_gauss_sum(theta, alpha)

    # -----------------------------------------------
    # Compute a single Gauss sum for given (theta, alpha)
    # -----------------------------------------------
    def compute_gauss_sum(self, theta, alpha):
        total = 0  # Initialize sum
        # Iterate over all nonzero elements in the finite field
        for x in self.finite_field_multiplicative_group:
            # Compute the additive character value
            additive_character_value = self.additive_character_generator^self.trace(x)
            # Compute the theta character value
            theta_character_value = self.multiplicative_character_generator^(theta * self.log(x))
            # Compute the alpha character value
            alpha_character_value = self.multiplicative_character_generator^(alpha * self.get_norm_log(x))
            # Multiply character values and add to the total sum
            total += additive_character_value * theta_character_value * alpha_character_value
        return total  # Return computed Gauss sum

    # -----------------------------------------------
    # Compute the norm logarithm of an element in GF(q^2)
    # Norm operation maps an element in GF(q^2) to GF(q)
    # -----------------------------------------------
    def get_norm_log(self, x):
        return (self.q + 1) * self.log(x)

    # -----------------------------------------------
    # Compute the logarithm of x in the finite field
    # Logarithm is computed relative to the field generator
    # -----------------------------------------------
    def log(self, x):
        return x.log(self.generator) if x != 0 else 0  # Log(0) is undefined, return 0

    # -----------------------------------------------
    # Compute the trace function in GF(q^2) to GF(q)
    # Trace sums the conjugates of an element in the field extension
    # -----------------------------------------------
    def trace(self, x):
        return x.trace()  # Built-in SageMath function

# -----------------------------------------------
# Compute the highest power of l that divides N
# This function is used to determine how N can be factored in a power of l
# -----------------------------------------------
def max_power(N, l):
    m = 0  # Initialize exponent count
    while N % l == 0:  # While N is divisible by l
        N //= l  # Divide N by l
        m += 1  # Increase exponent count
    return m  # Return highest power of l dividing N

# -----------------------------------------------
# Function to construct a Gauss sum table
# Takes in q (a prime power) and l (a prime)
# Computes necessary field parameters and constructs the table
# -----------------------------------------------
def fL_bar_gauss_sum_table(q, l):
    # Ensure l is prime
    if not is_prime(l):
        raise ValueError("l must be a prime number!")
    # Check if q is a prime power
    prime_power_result = q.is_prime_power(get_data=True)
    if prime_power_result[1] == 0:  # If q is not a prime power, raise error
        raise ValueError("Expected a prime power!")
    # Extract prime base p from the prime power decomposition of q
    p = prime_power_result[0]
    # Compute N, which is a multiple of q^2 - 1
    N = p * (q*q - 1)
    # Compute the highest power of l that divides N
    m = max_power(N, l)
    # Compute N' by factoring out the highest power of l from N
    N_prime = N // (l**m)
    # Compute the Carmichael function value of N_prime
    c = carmichael_lambda(N_prime)
    # Construct a finite field of order l^c
    F = GF(l**c)
    # Get a generator of the multiplicative group
    h = F.gen()
    # Construct and return the GaussSumTable object
    return GaussSumTable(q, h**((l**c - 1) // p), h**((p * (l**c - 1)) // N_prime))

# -----------------------------------------------
# Example function call to construct a Gauss sum table for q=9 and l=5
# -----------------------------------------------
gauss_sum_table_object = fL_bar_gauss_sum_table(9, 5)

# Store the computed table
table_of_gauss_sum = gauss_sum_table_object.table

# Print the computed Gauss sum table
print(table_of_gauss_sum)
