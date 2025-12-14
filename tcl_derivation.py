from sympy import (
    Add,
    Function,
    I,
    Mul,
    collect,
    expand,
    pprint,
    simplify,
    symbols,
)
from sympy.physics.quantum import TensorProduct

# Non‑commutative operator-valued functions
A = Function("A", commutative=False)
B = Function("B", commutative=False)

# Partial trace over the bath
TrB = Function("TrB", commutative=False)

# Symbolic commutator / anticommutator placeholders (non-evaluating)
Comm = Function("Comm", commutative=False)
AntiComm = Function("AntiComm", commutative=False)

# Correlation functions (scalars)
nu = Function("nu")
eta = Function("eta")

# Basic operators and symbols
Rho = symbols("Rho", commutative=False)   # full density operator
RhoB = symbols("RhoB", commutative=False) # bath reference state
RhoS = symbols("RhoS", commutative=False) # system reduced state Tr_B[Rho]

t, t1 = symbols("t t1")


def commutator(a, b):
    return a * b - b * a


def anticommutator(a, b):
    return a * b + b * a


def HamI(tt):
    return TensorProduct(A(tt), B(tt))


def P(rho):
    return TensorProduct(TrB(rho), RhoB)


def Q(rho):
    return rho - P(rho)


def L(tt, rho):
    return -I * commutator(HamI(tt), rho)


def K2(tt, rho):
    return TrB(P(L(tt, L(t1, P(rho)))))


# --- Rewrite machinery (close translation of Mathematica rules) ---

def fixed_point(expr, transform, max_iters=10):
    """Apply transform repeatedly until convergence or max_iters."""
    for _ in range(max_iters):
        new_expr = transform(expr)
        if new_expr == expr:
            return new_expr
        expr = new_expr
    return expr


def combine_tensor_products_once(expr):
    """TensorProduct[A,B] * TensorProduct[C,D] -> TensorProduct[A*C, B*D]."""
    if not expr.is_Mul:
        return expr
    comm, nc = expr.args_cnc()
    changed = False
    merged = []
    i = 0
    while i < len(nc):
        if (
            i + 1 < len(nc)
            and isinstance(nc[i], TensorProduct)
            and isinstance(nc[i + 1], TensorProduct)
        ):
            tp1, tp2 = nc[i], nc[i + 1]
            merged.append(TensorProduct(tp1.args[0] * tp2.args[0], tp1.args[1] * tp2.args[1]))
            i += 2
            changed = True
        else:
            merged.append(nc[i])
            i += 1
    if changed:
        return Mul(*comm) * Mul(*merged, evaluate=False)
    return expr


def combine_tensor_products(expr):
    return fixed_point(expr, lambda e: e.replace(lambda x: x.is_Mul, combine_tensor_products_once))


def trb_linearity(expr):
    """Linearity of TrB: pull out scalars and split sums."""
    def _rule(e):
        if e.func == TrB:
            arg = e.args[0]
            if arg.is_Add:
                return Add(*[TrB(term) for term in arg.args])
            if arg.is_Mul:
                comm, nc = arg.args_cnc()
                if comm:
                    return Mul(*comm) * TrB(Mul(*nc, evaluate=False))
        return e

    return expr.replace(lambda e: hasattr(e, "func") and e.func == TrB, _rule)


def partial_trace_rules(expr):
    """TrB[TensorProduct[A,B]] -> TensorProduct[A, TrB[B]]; TrB[Rho] -> RhoS."""
    expr = expr.replace(
        lambda e: e.func == TrB and isinstance(e.args[0], TensorProduct),
        lambda e: TensorProduct(e.args[0].args[0], TrB(e.args[0].args[1])),
    )
    return expr.subs(TrB(Rho), RhoS)


def drop_trb_rhob(expr):
    """TensorProduct[X, TrB[RhoB]] -> X."""
    return expr.replace(
        lambda e: isinstance(e, TensorProduct) and e.args[1] == TrB(RhoB),
        lambda e: e.args[0],
    )


def cyclic_trace(expr):
    """Cyclic property anchored on RhoB inside TrB."""
    def _cycle(e):
        if e.func == TrB:
            arg = e.args[0]
            if arg.is_Mul:
                factors = list(arg.args)
                if RhoB in factors and factors[-1] != RhoB:
                    idx = factors.index(RhoB)
                    pre, post = factors[:idx], factors[idx + 1 :]
                    return TrB(Mul(*post, *pre, RhoB, evaluate=False))
        return e

    return expr.replace(lambda e: hasattr(e, "func") and e.func == TrB, _cycle)


def bath_correlation(expr):
    """TensorProduct[X, TrB[B(t) B(t1) RhoB]] -> (nu + i eta) X."""
    Bt, Bt1 = B(t), B(t1)

    def _rule(e):
        if isinstance(e, TensorProduct) and getattr(e.args[1], "func", None) == TrB:
            arg = e.args[1].args[0]
            if arg == Bt * Bt1 * RhoB:
                return (nu(t - t1) + I * eta(t - t1)) * e.args[0]
            if arg == Bt1 * Bt * RhoB:
                return (nu(t1 - t) + I * eta(t1 - t)) * e.args[0]
        return e

    return expr.replace(
        lambda e: isinstance(e, TensorProduct) and getattr(e.args[1], "func", None) == TrB,
        _rule,
    )


def sign_rules(expr):
    """Use stationarity: eta(-τ)=-eta(τ), nu(-τ)=nu(τ)."""
    expr = expr.subs(eta(t1 - t), -eta(t - t1))
    expr = expr.subs(nu(t1 - t), nu(t - t1))
    expr = expr.subs(eta(-t + t1), -eta(t - t1))
    expr = expr.subs(nu(-t + t1), nu(t - t1))
    return expr


def tensor_simplify(expr):
    """Bundle of expansion + tensor-product merging + TrB linearity."""
    return trb_linearity(combine_tensor_products(expand(expr)))


def derive_tcl2():
    """Close translation of PerturbativeTCLDerivation.wl into SymPy."""
    r0 = K2(t, Rho)
    r1 = tensor_simplify(r0)
    r2 = partial_trace_rules(r1)
    r3 = drop_trb_rhob(r2)
    r3 = combine_tensor_products(r3)
    r4 = cyclic_trace(r3)
    r4 = tensor_simplify(r4)
    r5 = partial_trace_rules(r4)
    r5 = drop_trb_rhob(r5)
    r5 = cyclic_trace(r5)
    r6 = bath_correlation(r5)
    r6 = sign_rules(r6)
    r7 = simplify(collect(r6, [nu(t - t1), eta(t - t1)]))
    return r7


# Expose the final symbolic generator
tcl2_generator = derive_tcl2()


def group_bath_terms(expr):
    """Rewrite as nu(τ)*X + eta(τ)*Y, isolating bath scalars."""
    tau = t - t1
    nu_coeff = expr.expand().coeff(nu(tau))
    eta_coeff = expr.expand().coeff(eta(tau))
    return nu(tau) * nu_coeff + eta(tau) * eta_coeff, nu_coeff, eta_coeff


tcl2_grouped, tcl2_nu_coeff, tcl2_eta_coeff = group_bath_terms(tcl2_generator)


def split_by_A_ends(coeff):
    """
    Split a coefficient into terms where A(t) sits on the far left vs far right.
    Left group: noncommutative part starts with A(t).
    Right group: noncommutative part ends with A(t).
    """
    left_terms = []
    right_terms = []
    expr = coeff.expand()
    terms = expr.as_ordered_terms()
    for term in terms:
        if term.is_Mul:
            comm, nc = term.args_cnc()
            nc = list(nc)
            if nc:
                if nc[0] == A(t):
                    left_terms.append(term)
                    continue
                if nc[-1] == A(t):
                    right_terms.append(term)
                    continue
        # fallback: if no noncommutative factors, drop into left (unlikely)
        left_terms.append(term)
    left_sum = Add(*left_terms) if left_terms else 0
    right_sum = Add(*right_terms) if right_terms else 0
    return left_sum, right_sum


def factor_out_left(expr, head):
    """Factor `head` from the extreme left of every term (if present)."""
    if expr == 0:
        return 0
    terms = expr.as_ordered_terms()
    stripped = []
    for term in terms:
        if term.is_Mul:
            comm, nc = term.args_cnc()
            nc = list(nc)
            if nc and nc[0] == head:
                rest_nc = nc[1:]
                rest_mul = Mul(*comm) * (Mul(*rest_nc, evaluate=False) if rest_nc else 1)
                stripped.append(rest_mul)
                continue
        # fallback: cannot strip; return original expression un-factored
        return expr
    return head * Add(*stripped)


def factor_out_right(expr, tail):
    """Factor `tail` from the extreme right of every term (if present)."""
    if expr == 0:
        return 0
    terms = expr.as_ordered_terms()
    stripped = []
    for term in terms:
        if term.is_Mul:
            comm, nc = term.args_cnc()
            nc = list(nc)
            if nc and nc[-1] == tail:
                rest_nc = nc[:-1]
                rest_mul = Mul(*comm) * (Mul(*rest_nc, evaluate=False) if rest_nc else 1)
                stripped.append(rest_mul)
                continue
        # fallback: cannot strip; return original expression un-factored
        return expr
    return Add(*stripped) * tail


nu_left, nu_right = split_by_A_ends(tcl2_nu_coeff)
eta_left, eta_right = split_by_A_ends(tcl2_eta_coeff)

nu_left_factored = factor_out_left(nu_left, A(t))
nu_right_factored = factor_out_right(nu_right, A(t))
eta_left_factored = factor_out_left(eta_left, A(t))
eta_right_factored = factor_out_right(eta_right, A(t))


def to_comm_or_anticomm(expr):
    """
    Recognize two-term sums/differences of the form X*Y ± Y*X and rewrite
    as (anti)commutators. Assumes overall scalars already factored out.
    """
    from sympy import Mul

    if not expr.is_Add or len(expr.args) != 2:
        return expr

    a, b = expr.args

    def decompose(term):
        coeff, factors = term.as_coeff_mul()
        factors = list(factors)
        if len(factors) != 2:
            return None
        return coeff, factors[0], factors[1]

    da = decompose(a)
    db = decompose(b)
    if not da or not db:
        return expr

    ca, xa, ya = da
    cb, xb, yb = db

    reversed_pair = xa == yb and ya == xb
    if not reversed_pair:
        return expr

    if ca == 1 and cb == 1:
        return AntiComm(xa, ya)
    if ca == 1 and cb == -1:
        return Comm(xa, ya)
    if ca == -1 and cb == 1:
        return -Comm(xa, ya)
    if ca == -1 and cb == -1:
        return -AntiComm(xa, ya)
    return expr


def compress_side(expr, side):
    """
    Pull scalars next to A(t) and rewrite inner sum/difference as comm/anticomm.
    side: 'left' means A(t) sits on the left; 'right' means on the right.
    """
    if expr == 0:
        return 0

    if not expr.is_Mul:
        return expr

    comm, nc = expr.args_cnc()
    scalar = Mul(*comm) if comm else 1

    if side == "left":
        if not nc or nc[0] != A(t):
            return expr
        rest_nc = nc[1:]
        rest = Mul(*rest_nc, evaluate=False) if rest_nc else 1
    else:  # right
        if not nc or nc[-1] != A(t):
            return expr
        rest_nc = nc[:-1]
        rest = Mul(*rest_nc, evaluate=False) if rest_nc else 1

    # Factor out additional scalars from rest
    from sympy import factor_terms

    rest_fact = factor_terms(rest)
    rcomm, rnc = rest_fact.args_cnc() if rest_fact.is_Mul else ([], [rest_fact])
    scalar *= Mul(*rcomm) if rcomm else 1
    rest_core = Mul(*rnc, evaluate=False) if rnc else 1

    rest_core = to_comm_or_anticomm(rest_core.expand())

    if side == "left":
        return scalar * A(t) * rest_core
    else:
        return rest_core * scalar * A(t)


nu_left_compact = compress_side(nu_left_factored, "left")
nu_right_compact = compress_side(nu_right_factored, "right")
eta_left_compact = compress_side(eta_left_factored, "left")
eta_right_compact = compress_side(eta_right_factored, "right")


def combine_bidirectional(expr, head):
    """
    If expr = s1*head*X + s2*X*head, map to:
      - s1*Comm(head, X) when s2 == -s1
      - s1*AntiComm(head, X) when s2 == s1
    Otherwise return expr unchanged.
    """
    from sympy import Wild

    expr_exp = expr.expand()
    a = Wild("a", properties=[lambda s: s.is_commutative])
    b = Wild("b", properties=[lambda s: s.is_commutative])
    x = Wild("x", commutative=False)
    y = Wild("y", commutative=False)

    m = expr_exp.match(a * head * x + b * y * head)
    if not m:
        return expr

    a_val = m.get(a, 1)
    b_val = m.get(b, 1)
    X = m[x]
    Y = m[y]

    def canon_comm(expr_in):
        if expr_in.func == Comm and len(expr_in.args) == 2:
            p, q = expr_in.args
            if p.sort_key() > q.sort_key():
                return -Comm(q, p)
        if expr_in.func == AntiComm and len(expr_in.args) == 2:
            p, q = expr_in.args
            if p.sort_key() > q.sort_key():
                return AntiComm(q, p)
        return expr_in

    Xc = canon_comm(X)
    Yc = canon_comm(Y)

    if Yc == Xc:
        same = True
        neg = False
    elif Yc == -Xc or simplify(Yc + Xc) == 0:
        same = False
        neg = True
    else:
        return expr

    if same and b_val == -a_val:
        return a_val * Comm(head, Xc)
    if same and b_val == a_val:
        return a_val * AntiComm(head, Xc)
    if neg and b_val == a_val:
        return a_val * Comm(head, Xc)
    if neg and b_val == -a_val:
        return a_val * AntiComm(head, Xc)
    return expr


nu_combined = combine_bidirectional(nu_left_compact + nu_right_compact, A(t))
eta_combined = combine_bidirectional(eta_left_compact + eta_right_compact, A(t))


if __name__ == "__main__":
    print("Full TCL2 generator:")
    pprint(tcl2_generator)
    print("\nGrouped by bath correlations:")
    pprint(tcl2_grouped)
    print("\nCoefficient of nu(t - t1):")
    pprint(tcl2_nu_coeff)
    print("\nCoefficient of eta(t - t1):")
    pprint(tcl2_eta_coeff)
    print("\nnu left factored A(t)*[...] :")
    pprint(nu_left_factored)
    print("\nnu right factored [...]*A(t) :")
    pprint(nu_right_factored)
    print("\neta left factored A(t)*[...] :")
    pprint(eta_left_factored)
    print("\neta right factored [...]*A(t) :")
    pprint(eta_right_factored)
    print("\nnu left compact (comm/anticomm):")
    pprint(nu_left_compact)
    print("\nnu right compact (comm/anticomm):")
    pprint(nu_right_compact)
    print("\neta left compact (comm/anticomm):")
    pprint(eta_left_compact)
    print("\neta right compact (comm/anticomm):")
    pprint(eta_right_compact)
    print("\nnu final compact:")
    pprint(nu_combined)
    print("\neta final compact:")
    pprint(eta_combined)
    print("\nFinal TCL2 generator (compact, grouped by bath correlations):")
    final_compact = nu(t - t1) * nu_combined + eta(t - t1) * eta_combined
    pprint(final_compact)
