from sympy import symbols, I, expand, simplify, collect
from sympy.physics.quantum import TensorProduct, Dagger, Operator, Commutator, qapply
from sympy.physics.quantum.operator import Operator

# 1. Define Symbols and Operators
t, t1 = symbols('t t1', real=True)
nu = symbols('nu', real=True)   # nu(t-t1)
eta = symbols('eta', real=True) # eta(t-t1)

# Define abstract operators
# We use a class to ensure they behave like non-commuting quantum operators
class TimeOp(Operator):
    pass

A_t = TimeOp('A(t)')
B_t = TimeOp('B(t)')
A_t1 = TimeOp('A(t1)')
B_t1 = TimeOp('B(t1)')
RhoS = TimeOp('RhoS')
RhoB = TimeOp('RhoB')

# 2. Define Helper Functions (Mathematica Equivalents)

def ham_i(time):
    """Returns A(time) (x) B(time)"""
    if time == t:
        return TensorProduct(A_t, B_t)
    elif time == t1:
        return TensorProduct(A_t1, B_t1)

def comm(a, b):
    """Computes commutator [A, B]"""
    return a * b - b * a

def liouvillian(time, rho):
    """L(t, rho) = -i [H(t), rho]"""
    h = ham_i(time)
    return -I * comm(h, rho)

# 3. Implement the Logic: K2 Calculation
# Corresponds to: K2[t_, Rho_] := TrB[P[L[t, L[t1, P[Rho]]]]]
# Note: In your Mathematica code, P[Rho] projects to RhoS (x) RhoB.
# We assume the input 'rho' is already the projected state RhoS (x) RhoB for the derivation.

# Initial State: RhoS (x) RhoB
rho_initial = TensorProduct(RhoS, RhoB)

# Inner Liouvillian: L(t1, rho)
term1 = liouvillian(t1, rho_initial)

# Apply "P" (Projector) effectively resets the Bath to RhoB and traces the previous bath part.
# However, standard TCL2 derivation usually keeps the structure until the final trace.
# Looking closely at your code: K2 = TrB[ L(t, L(t1, Rho)) ] (if P is effectively identity on uncorrelated states initially)
# Let's compute the double commutator directly, then apply the trace rules.

# Expression: -i [H(t), -i [H(t1), RhoS (x) RhoB]]
expr = liouvillian(t, term1)

# Expand the expression to get individual terms
# qapply handles the TensorProduct algebra: (A (x) B) * (C (x) D) -> (AC) (x) (BD)
expanded_expr = expand(qapply(expr))

# 4. Custom Trace and Substitution Logic
# This replaces your 'TensorSimplification', 'PartialTraceRule', and 'CyclicTraceRule'

def get_bath_trace_val(bath_part):
    """
    Evaluates Tr_B[ bath_part ].
    Matches B(t), B(t1), and RhoB patterns to returns nu/eta.
    """
    # We expect terms like: B(t)*B(t1)*RhoB or B(t1)*B(t)*RhoB
    # Note: The cyclic property means Tr[X Y RhoB] = Tr[Y RhoB X] etc.
    
    # Check the order of operators in the product
    # We simply check the string representation or args to identify the sequence
    # Target: Tr[ B(t) B(t1) RhoB ] -> nu + i*eta
    # Target: Tr[ B(t1) B(t) RhoB ] -> nu - i*eta
    
    # Convert to string for easy pattern matching (robust enough for this specific derivation)
    s = str(bath_part)
    
    # Case 1: B(t) ... B(t1) ... RhoB (Normal time order)
    if "B(t)" in s and "B(t1)" in s and s.find("B(t)") < s.find("B(t1)"):
        return nu + I*eta
    
    # Case 2: B(t1) ... B(t) ... RhoB (Reverse time order)
    if "B(t1)" in s and "B(t)" in s and s.find("B(t1)") < s.find("B(t)"):
        return nu - I*eta
        
    return 0

final_terms = []

# Loop over the expanded sum
# If it's a sum, iterate args. If single term, make it a list.
terms = expanded_expr.args if expanded_expr.is_Add else [expanded_expr]

for term in terms:
    coeff = 1
    operators = term
    
    # Separate coefficient (like -1, I) from operators
    if term.is_Mul:
        coeff = term.as_coeff_Mul()[0]
        operators = term.as_coeff_Mul()[1]
        
    # 'operators' should be a TensorProduct (System (x) Bath)
    if isinstance(operators, TensorProduct):
        sys_part = operators.args[0]
        bath_part = operators.args[1]
        
        # Calculate Trace of the Bath part
        # This replaces TrB[...] -> scalar
        trace_val = get_bath_trace_val(bath_part)
        
        # Combine: coeff * sys_part * trace_val
        final_terms.append(coeff * sys_part * trace_val)

# Sum up the results
result = sum(final_terms)

# 5. Final Simplification
# Group by nu and eta, similar to Collect[..., {nu, eta}]
result = collect(expand(result), [nu, eta])

print("--- Final TCL2 Generator Output ---")
print(result)

# Optional: Verify against manually expected form
# Expected: [A(t), A(t1) RhoS] (nu + i eta) + [RhoS A(t1), A(t)] (nu - i eta) ... (grouped appropriately)
