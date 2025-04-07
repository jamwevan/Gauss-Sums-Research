# Import functions to check primality and generate primes from sympy
from sympy import isprime, primerange

# Check whether n is a prime power (i.e., n = p^k for some prime p and integer k >= 1)
def is_prime_power(n):
    """Returns True if n is a prime power (p^k), False otherwise."""

    # Reject n < 2 immediately
    if n < 2:
        return False

    # Try all primes up to n
    for p in primerange(2, n + 1):
        k = 1

        # Raise p to increasing powers until it exceeds n
        while (power := p**k) <= n:

            # If p^k equals n, n is a prime power
            if power == n:
                return True
            k += 1

    # If no match was found, n is not a prime power
    return False

# Generate the first `count` many (ℓ, q) pairs such that q = 1 + 2 * ℓ^j and q is a prime power
def generate_lq_pairs(count):

    # Stores valid (ℓ, q) pairs
    results = []

    # Start testing from the first prime ℓ = 2
    l = 2

    # Continue until we have the desired number of pairs
    while len(results) < count:

        # Skip ℓ if it is not prime
        if not isprime(l):
            l += 1
            continue

        # Try increasing values of j for this ℓ
        j = 1
        while True:

            # Construct q in the required form
            q = 1 + 2 * (l**j)

            # If q is a prime power, store the pair
            if is_prime_power(q):
                results.append((l, q))

                # Stop once enough pairs are found
                if len(results) >= count:
                    break

            # Prevent generating excessively large q values
            elif q > 5000:
                break

            # Try next exponent j
            j += 1

        # Move on to next prime ℓ
        l += 1

    return results

# Main execution block
if __name__ == "__main__":

    # Ask user how many (ℓ, q) pairs to generate
    x = int(input("How many (ℓ, q) pairs do you want? "))

    # Generate the list of valid pairs
    pairs = generate_lq_pairs(x)

    # Print header message
    print("\nValid (ℓ, q) pairs such that q = 1 + 2 * ℓ^j and q is a prime power:\n")

    # Print each pair in nice format
    for l, q in pairs:
        print(f"(ℓ = {l}, q = {q})")

    # Find the largest ℓ value among all pairs
    max_l = max(pairs, key=lambda t: t[0])[0]

    # Find the largest q value among all pairs
    max_q = max(pairs, key=lambda t: t[1])[1]

    # Display the largest ℓ and q summary
    print(f"\n(ℓ_max, q_max) = ({max_l}, {max_q})")
