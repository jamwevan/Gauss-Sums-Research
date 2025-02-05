from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from sage.all import GF, e, pi, I, latex

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
    if prime_power_result[1] == 0:
        raise ValueError("Expected a prime power!")
    p = prime_power_result[0]
        
    return GaussSumTable(q, e^(2 * pi * I / p), e^(2 * pi * I / (q*q - 1)))

def save_gauss_sum_table_as_latex_pdf(table, filename="GaussSumTable_output.pdf"):
    """
    Saves the Gauss sum table in a structured LaTeX format inside a PDF.
    Each entry appears on a separate row and is sorted.
    """
    # Flatten the table into a single list and sort it
    sorted_entries = sorted([latex(entry) for row in table for entry in row])

    with PdfPages(filename) as pdf:
        fig, ax = plt.subplots(figsize=(8.5, 11))  # Letter-sized page
        ax.axis('off')  # Remove axes

        # Title
        ax.text(0.5, 0.95, "Gauss Sum Table", fontsize=14, ha='center')

        y_position = 0.85  # Start position for the equations

        for entry in sorted_entries:
            latex_equation = r"$" + entry + r"$"  # Convert to LaTeX equation format
            
            ax.text(0.1, y_position, latex_equation, fontsize=12, ha='left', family='monospace')
            y_position -= 0.05  # Move down for the next line
            
            if y_position < 0.1:  # Create a new page if out of space
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)
                fig, ax = plt.subplots(figsize=(8.5, 11))
                ax.axis('off')
                y_position = 0.85  # Reset for new page
        
        pdf.savefig(fig, bbox_inches='tight')  # Save the final page
        plt.close(fig)
    
    print(f"Saved Gauss Sum Table as {filename}")

# Generate the Gauss Sum Table
gauss_sum_table_object = complex_gauss_sum_table(3)
table_of_gauss_sum = gauss_sum_table_object.table

# Save as LaTeX-formatted PDF
save_gauss_sum_table_as_latex_pdf(table_of_gauss_sum)
