class GaussSumTable:
    def __init__(self, q, additive_character_generator, multiplicative_character_generator):
        self.q = q
        self.additive_character_generator = additive_character_generator
        self.multiplicative_character_generator = multiplicative_character_generator
        
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

def complex_gauss_sum_table(q):
    prime_power_result = q.is_prime_power(get_data=True)
    if  prime_power_result[1] == 0:
        raise ValueError("Expected a prime power!")
    p = prime_power_result[0]
        
    return GaussSumTable(q, e^(2 * pi * I / p), e^(2 * pi * I / (q*q - 1)))

gauss_sum_table_object = complex_gauss_sum_table(3)
table_of_gauss_sum = gauss_sum_table_object.table
# Print the table to see the output
print(table_of_gauss_sum)
