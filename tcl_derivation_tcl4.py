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
    factor_terms,
)
from tcl_derivation import (
    A,
    B,
    TrB,
    RhoB,
    RhoS,
    commutator,
    HamI,
    L,
    P,
    tensor_simplify,
    partial_trace_rules,
    drop_trb_rhob,
    combine_tensor_products,
    trb_linearity,
    cyclic_trace,
    bath_correlation,
    sign_rules,
    nu,
    eta,
    t,
    t1,
    Comm,
    AntiComm,
)

# Additional time symbols for TCL4
t2, t3 = symbols("t2 t3")


def pairings(indices):
    """Generate all pairings (perfect matchings) of a list of indices, preserving order semantics."""
    if not indices:
        return [[]]
    first = indices[0]
    rest = indices[1:]
    results = []
    for i, partner in enumerate(rest):
        remaining = rest[:i] + rest[i + 1 :]
        for sub in pairings(remaining):
            results.append([(first, partner)] + sub)
    return results


def corr_pair(time_i, time_j):
    """Bath two-point correlation with ordering."""
    return nu(time_i - time_j) + I * eta(time_i - time_j)


def wick_decompose(expr):
    """
    Replace TensorProduct[X, TrB[B(ti)...B(tj) RhoB]] with
    sum over Wick pairings of products of two-point correlations times X.
    """

    from sympy.physics.quantum import TensorProduct

    def target(e):
        return isinstance(e, TensorProduct) and len(e.args) == 2 and getattr(e.args[1], "func", None) == TrB

    def repl(e):
        trb_arg = e.args[1].args[0]
        if not trb_arg.is_Mul:
            return e
        factors = list(trb_arg.args)
        if RhoB not in factors:
            return e
        idx = factors.index(RhoB)
        bath_ops = factors[:idx]
        if factors[idx + 1 :]:
            return e  # extra operators after RhoB, skip
        times = []
        for op in bath_ops:
            if getattr(op, "func", None) == B and len(op.args) == 1:
                times.append(op.args[0])
            else:
                return e
        n = len(times)
        if n < 2 or n % 2 == 1:
            return e

        total = 0
        for pairing in pairings(list(range(n))):
            term = 1
            for i, j in pairing:
                term *= corr_pair(times[i], times[j])
            total += term
        return total * e.args[0]

    return expr.replace(target, repl)


def K4(rho):
    """
    TCL4 kernel (unsymmetrised) acting on rho, before trace over bath:
    P L(t) L(t1) L(t2) L(t3) P
      - P L(t) L(t1) P L(t2) L(t3) P
      - P L(t) L(t2) P L(t1) L(t3) P
      - P L(t) L(t3) P L(t1) L(t2) P
    """
    term1 = P(L(t, L(t1, L(t2, L(t3, P(rho))))))
    term2 = P(L(t, L(t1, P(L(t2, L(t3, P(rho)))))))
    term3 = P(L(t, L(t2, P(L(t1, L(t3, P(rho)))))))
    term4 = P(L(t, L(t3, P(L(t1, L(t2, P(rho)))))))
    return TrB(term1 - term2 - term3 - term4)


def derive_tcl4():
    expr = K4(symbols("Rho", commutative=False))
    # Mirror the TCL2 simplification pipeline
    expr = tensor_simplify(expr)
    expr = partial_trace_rules(expr)
    expr = drop_trb_rhob(expr)
    expr = combine_tensor_products(expr)
    expr = cyclic_trace(expr)
    expr = tensor_simplify(expr)
    expr = partial_trace_rules(expr)
    expr = drop_trb_rhob(expr)
    expr = cyclic_trace(expr)
    # Wick decompose multi-time bath correlations into two-point functions
    expr = wick_decompose(expr)
    # Catch any remaining two-point correlators
    expr = bath_correlation(expr)
    expr = sign_rules(expr)

    # Collect over all pairwise bath correlators that may appear
    taus = [
        t - t1,
        t - t2,
        t - t3,
        t1 - t2,
        t1 - t3,
        t2 - t3,
    ]
    collect_terms = []
    for tau in taus:
        collect_terms.append(nu(tau))
        collect_terms.append(eta(tau))

    expr = simplify(collect(expr, collect_terms))
    return expr


def group_by_correlators(expr):
    """Group terms sharing identical products of nu/eta factors."""
    taus = [t - t1, t - t2, t - t3, t1 - t2, t1 - t3, t2 - t3]
    corr_syms = [nu(tau) for tau in taus] + [eta(tau) for tau in taus]
    terms = expand(expr).as_ordered_terms()
    grouped = {}
    for term in terms:
        comm, nc = term.args_cnc()
        comm_expr = simplify(Mul(*comm))
        key = comm_expr
        val = Mul(*nc)
        grouped[key] = grouped.get(key, 0) + val
    grouped_expr = Add(*[k * v for k, v in grouped.items()])
    return grouped_expr, grouped


def startswith_head(term, head):
    comm, nc = term.args_cnc()
    nc = list(nc)
    if nc and nc[0] == head:
        coeff = Mul(*comm) if comm else 1
        rest = Mul(*nc[1:], evaluate=False) if len(nc) > 1 else 1
        return coeff, rest
    return None


def endswith_head(term, head):
    comm, nc = term.args_cnc()
    nc = list(nc)
    if nc and nc[-1] == head:
        coeff = Mul(*comm) if comm else 1
        rest = Mul(*nc[:-1], evaluate=False) if len(nc) > 1 else 1
        return coeff, rest
    return None


def sort_by_head(expr, head):
    """Group terms with head on left/right; pair by identical remainders -> comm/anticomm."""
    left_map = {}
    right_map = {}
    others = 0
    for term in expand(expr).as_ordered_terms():
        ls = startswith_head(term, head)
        if ls:
            coeff, rest = ls
            left_map[rest] = left_map.get(rest, 0) + coeff
            continue
        rs = endswith_head(term, head)
        if rs:
            coeff, rest = rs
            right_map[rest] = right_map.get(rest, 0) + coeff
            continue
        others += term

    combined = 0
    remainders = set(left_map) | set(right_map)
    for rem in remainders:
        cL = left_map.get(rem, 0)
        cR = right_map.get(rem, 0)
        if cL != 0 and cR != 0:
            if simplify(cL - cR) == 0:
                combined += cL * AntiComm(head, rem)
            elif simplify(cL + cR) == 0:
                combined += cL * Comm(head, rem)
            else:
                combined += head * cL * rem + cR * rem * head
        elif cL != 0:
            combined += head * cL * rem
        else:
            combined += cR * rem * head

    return simplify(combined + others)


def nest_heads(expr, heads):
    out = expr
    for h in heads:
        out = sort_by_head(out, h)
    return out


if __name__ == "__main__":
    tcl4_generator = derive_tcl4()
    grouped_expr, grouped = group_by_correlators(tcl4_generator)
    heads_order = [A(t), A(t1), A(t2), A(t3)]
    nested_terms = []
    for pref, op in grouped.items():
        nested_op = nest_heads(op, heads_order)
        nested_terms.append(pref * nested_op)
    nested_expr = Add(*nested_terms)
    print("TCL4 generator, compacted into nested comm/anticomm with respect to A(t),A(t1),A(t2),A(t3):")
    pprint(nested_expr)
