class GaussSumTable:
    def __init__(self, q, additive_character_generator, multiplicative_character_generator):
        self.q = q
        self.additive_character_generator = additive_character_generator
        self.multiplicative_character_generator = multiplicative_character_generator
        # TODO: Does this stay the same?
        self.finite_field = GF(q^2)
        self.generator = self.finite_field.gen()
        self.finite_field_elements = list(self.finite_field)
        self.finite_field_multiplicative_group = [x for x in self.finite_field_elements if x != 0]
        self.table = [[0 for _ in range(q - 1)] for _ in range(q^2 - 1)]
        self.compute_gauss_sum_table()

    def compute_gauss_sum_table(self):
        for theta in range(self.q^2 - 1):
            for alpha in range(self.q - 1):
                self.table[theta][alpha] = self.compute_gauss_sum(theta, alpha)

    def compute_gauss_sum(self, theta, alpha):
        total = 0
        for x in self.finite_field_multiplicative_group:
            additive_character_value = self.additive_character_generator^self.trace(x)
            theta_character_value = self.multiplicative_character_generator^(theta * self.log(x))
            alpha_character_value = self.multiplicative_character_generator^(alpha * self.get_norm_log(x))
            total += additive_character_value * theta_character_value * alpha_character_value
        return total

    def get_norm_log(self, x):
        return (self.q + 1) * self.log(x)

    def log(self, x):
        return x.log(self.generator) if x != 0 else 0  # Using the logarithm in finite fields

    def trace(self, x):
        return x.trace()  # Using the field trace function from GF(q^2) to GF(q)

# Change: Implementing helper function to calculate m
    def max_power(N, l):
        m = 0
        while N % l == 0:
            N //= l
            m += 1
        return m
    
# Change: Updated Function Arguments
def complex_gauss_sum_table(q, l):
    # Change: Checked if l is prime
    if not isprime(l):
        raise ValueError("l must be a prime number!")
    # Check if q is a prime power
    prime_power_result = q.is_prime_power(get_data=True)
    if  prime_power_result[1] == 0:
        raise ValueError("Expected a prime power!")
    p = prime_power_result[0]
    # Change: Implementing Extension Algo
    N = p*(q*q-1)
    m = max_power(n, l)
    N_prime = n/(l**m)
    c = carmichael_function(n_prime)
    F = GF(l**c)
    h = F.gen()
    # Change: Updated arguments in function call
    return GaussSumTable(q, h**((l**c-1)/p), h**((p*(l**c-1))/N_prime))

# TODO: Update funciton call by adding an l value here William thinks is good
gauss_sum_table_object = complex_gauss_sum_table(3)
table_of_gauss_sum = gauss_sum_table_object.table
