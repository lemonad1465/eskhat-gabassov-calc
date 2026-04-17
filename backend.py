"""
4D Numbers Algebra System CAS Calculator v3
- Arithmetic: multiplication table (exact algebra)
- Functions: spectral method (exact analytical solutions, no truncation errors)
- Supports 6 anisotropic 4D spaces (M2~M7)
- Exact spectrum/inverse formulas for both isotropic & anisotropic spaces
"""

from sympy import (
    sqrt, simplify, sympify,
    cos, sin, tan, exp, log, cosh, sinh, tanh,
    pi, E, I, re, im, Abs, conjugate, SympifyError
)
import sympy as sp
from typing import Dict, Any, List, Tuple, Optional

# ══════════════════════════════════════════════════════════════════
#  Space Definitions
# ══════════════════════════════════════════════════════════════════

SPACES: Dict[str, Dict[str, Any]] = {
    "M2": {"signs": {"beta": -1, "gamma": +1, "delta": +1},
           "desc": "β<0, γ>0, δ>0  (J1,J2 real; J3,J4 imaginary)"},
    "M3": {"signs": {"beta": +1, "gamma": -1, "delta": +1},
           "desc": "β>0, γ<0, δ>0  (J1,J3 real; J2,J4 imaginary)"},
    "M4": {"signs": {"beta": +1, "gamma": +1, "delta": -1},
           "desc": "β>0, γ>0, δ<0  (J1,J4 real; J2,J3 imaginary)"},
    "M5": {"signs": {"beta": -1, "gamma": -1, "delta": +1},
           "desc": "β<0, γ<0, δ>0"},
    "M6": {"signs": {"beta": -1, "gamma": +1, "delta": -1},
           "desc": "β<0, γ>0, δ<0"},
    "M7": {"signs": {"beta": +1, "gamma": -1, "delta": -1},
           "desc": "β>0, γ<0, δ<0  (J1,J2 and J3,J4 symmetric conjugate)"},
}

# ══════════════════════════════════════════════════════════════════
#  Expression Parsing
# ══════════════════════════════════════════════════════════════════

def parse_expr(raw: str):
    raw = raw.strip().replace("^", "**")
    try:
        return sympify(raw, locals={
            "pi": pi, "e": E, "E": E,
            "cos": cos, "sin": sin, "tan": tan,
            "cosh": cosh, "sinh": sinh, "tanh": tanh,
            "exp": exp, "log": log, "ln": log,
            "sqrt": sqrt, "Abs": Abs, "I": I,
        })
    except SympifyError as err:
        raise ValueError(f"Cannot parse '{raw}': {err}")

def input_expr(prompt: str):
    while True:
        raw = input(prompt)
        try:
            return parse_expr(raw)
        except ValueError as e:
            print(f"  ⚠  {e}, please try again.")

def input_quaternion(label: str):
    print(f"\n  Input 4D number {label} (supports expressions, e.g., cos(pi/3), sqrt(2)):")
    return tuple(input_expr(f"    {label}_{i} = ") for i in range(1, 5))

# ══════════════════════════════════════════════════════════════════
#  4D Number Class (Multiplication Table Operations)
# ══════════════════════════════════════════════════════════════════

class QuadNum:
    def __init__(self, comps, params, space_key: str, is_isotropic: bool):
        self.c = [simplify(c) for c in comps]
        self.p = params
        self.space_key = space_key
        self.is_isotropic = is_isotropic

    def __add__(self, o):
        return QuadNum([self.c[i] + o.c[i] for i in range(4)], self.p, self.space_key, self.is_isotropic)

    def __sub__(self, o):
        return QuadNum([self.c[i] - o.c[i] for i in range(4)], self.p, self.space_key, self.is_isotropic)

    def __mul__(self, o):
        x, y = self.c, o.c
        a, b, g, d = self.p["alpha"], self.p["beta"], self.p["gamma"], self.p["delta"]
        return QuadNum([
            a*x[0]*y[0] + (g*d/a)*x[1]*y[1] + (b*d/a)*x[2]*y[2] + (b*g/a)*x[3]*y[3],
            a*x[1]*y[0] + a*x[0]*y[1] + b*x[3]*y[2] + b*x[2]*y[3],
            a*x[2]*y[0] + g*x[3]*y[1] + a*x[0]*y[2] + g*x[1]*y[3],
            a*x[3]*y[0] + d*x[2]*y[1] + d*x[1]*y[2] + a*x[0]*y[3],
        ], self.p, self.space_key, self.is_isotropic)

    def symp_conj(self):
        return QuadNum([self.c[0], -self.c[1], -self.c[2], -self.c[3]], self.p, self.space_key, self.is_isotropic)

    def symp_norm_sq(self):
        return simplify((self * self.symp_conj()).c[0])

    def inv(self):
        mu1, mu2, mu3, mu4 = compute_spectrum(self.c, self.p, self.space_key, self.is_isotropic)
        if simplify(Abs(mu1)) == 0 or simplify(Abs(mu3)) == 0:
            raise ValueError("Zero Divisor Error: Spectrum contains 0. Inverse does not exist.")
        lam1 = simplify(1 / mu1)
        lam3 = simplify(1 / mu3)
        return spectrum_to_quad(lam1, lam3, self.p, self.space_key, self.is_isotropic)

    def __truediv__(self, o): return self * o.inv()

    def __pow__(self, n: int):
        if n == 0:  return QuadNum([sp.Integer(1), 0, 0, 0], self.p, self.space_key, self.is_isotropic)
        if n < 0:   return self.inv() ** (-n)
        
        mu1, mu2, mu3, mu4 = compute_spectrum(self.c, self.p, self.space_key, self.is_isotropic)
        lam1 = simplify(mu1 ** n)
        lam3 = simplify(mu3 ** n)
        
        return spectrum_to_quad(lam1, lam3, self.p, self.space_key, self.is_isotropic)

    def scale(self, s):
        return QuadNum([simplify(s * v) for v in self.c], self.p, self.space_key, self.is_isotropic)

# ══════════════════════════════════════════════════════════════════
#  Spectral Method: Spectrum Formula 
# ══════════════════════════════════════════════════════════════════

def compute_spectrum(comps, params, space_key: str, is_isotropic: bool):
    """
    Returns (mu1, mu2, mu3, mu4).
    mu2 = conj(mu1), mu4 = conj(mu3) always hold.
    """
    x1, x2, x3, x4 = comps
    a, b, g, d = params["alpha"], params["beta"], params["gamma"], params["delta"]

    sGD = sqrt(Abs(g * d))
    sBD = sqrt(Abs(b * d))
    sBG = sqrt(Abs(b * g))

    if space_key == "M2":
        mu1 = (a*x1 - sGD*x2) + I*(sBD*x3 - sBG*x4)
        mu3 = (a*x1 + sGD*x2) + I*(sBD*x3 + sBG*x4)

    elif space_key == "M3":
        mu1 = (a*x1 - sBD*x3) + I*(sGD*x2 - sBG*x4)
        mu3 = (a*x1 + sBD*x3) + I*(sGD*x2 + sBG*x4)

    elif space_key == "M4":
        mu1 = (a*x1 - sBG*x4) + I*(sGD*x2 - sBD*x3)
        mu3 = (a*x1 + sBG*x4) + I*(sGD*x2 + sBD*x3)

    elif space_key == "M5":
        mu1 = (a*x1 - sBG*x4) + I*(sGD*x2 + sBD*x3)
        mu3 = (a*x1 + sBG*x4) + I*(sGD*x2 - sBD*x3)

    elif space_key == "M6":
        mu1 = (a*x1 - sBD*x3) + I*(sGD*x2 + sBG*x4)
        mu3 = (a*x1 + sBD*x3) + I*(sGD*x2 - sBG*x4)

    elif space_key == "M7":
        mu1 = (a*x1 - sGD*x2) + I*(sBD*x3 + sBG*x4)
        mu3 = (a*x1 + sGD*x2) + I*(sBD*x3 - sBG*x4)

    else:
        raise NotImplementedError(f"Space {space_key} is not implemented yet.")

    mu2 = conjugate(mu1)
    mu4 = conjugate(mu3)
    return tuple(simplify(m) for m in (mu1, mu2, mu3, mu4))

# ══════════════════════════════════════════════════════════════════
#  Spectral Method: Inverse Formula (New Spectrum to 4D Coordinates)
# ══════════════════════════════════════════════════════════════════

def spectrum_to_quad(lam1, lam3, params, space_key: str, is_isotropic: bool):
    a_sym = params["alpha"]
    b_sym = params["beta"]
    g_sym = params["gamma"]
    d_sym = params["delta"]

    A = simplify(re(lam1) + re(lam3))
    B = simplify(re(lam3) - re(lam1))
    C = simplify(im(lam1) + im(lam3))
    D = simplify(im(lam3) - im(lam1))

    sGD = sqrt(Abs(g_sym * d_sym))
    sBD = sqrt(Abs(b_sym * d_sym))
    sBG = sqrt(Abs(b_sym * g_sym))

    if space_key == "M2":
        w1 = A / (2 * a_sym)
        w2 = B / (2 * sGD)
        w3 = C / (2 * sBD)
        w4 = D / (2 * sBG)

    elif space_key == "M3":
        w1 = A / (2 * a_sym)
        w2 = C / (2 * sGD)
        w3 = B / (2 * sBD)
        w4 = D / (2 * sBG)

    elif space_key == "M4":
        w1 = A / (2 * a_sym)
        w2 = C / (2 * sGD)
        w3 = D / (2 * sBD)
        w4 = B / (2 * sBG)

    elif space_key == "M5":
        w1 = A / (2 * a_sym)
        w2 = C / (2 * sGD)
        w3 = -D / (2 * sBD)
        w4 = B / (2 * sBG)

    elif space_key == "M6":
        w1 = A / (2 * a_sym)
        w2 = C / (2 * sGD)
        w3 = B / (2 * sBD)
        w4 = -D / (2 * sBG)

    elif space_key == "M7":
        w1 = A / (2 * a_sym)
        w2 = B / (2 * sGD)
        w3 = C / (2 * sBD)
        w4 = -D / (2 * sBG)

    else:
        raise NotImplementedError(f"Space {space_key} inverse formula is not implemented yet.")

    return QuadNum([simplify(w) for w in (w1, w2, w3, w4)], params, space_key, is_isotropic)

# ══════════════════════════════════════════════════════════════════
#  Spectral Method Framework
# ══════════════════════════════════════════════════════════════════

def apply_spectral_func(X: QuadNum, func_c, space_key: str, is_isotropic: bool):
    """
    Mapping -> Complex Function -> Inverse Mapping
    """
    mu1, mu2, mu3, mu4 = compute_spectrum(X.c, X.p, X.space_key, X.is_isotropic)
    print("  (Computing spectrum...)")
    lam1 = simplify(func_c(mu1))
    lam3 = simplify(func_c(mu3))
    print("  (Restoring 4D coordinates from spectrum...)")
    return spectrum_to_quad(lam1, lam3, X.p, X.space_key, X.is_isotropic)

def quad_exp(X, sk, iso):  return apply_spectral_func(X, exp,  sk, iso)
def quad_cos(X, sk, iso):  return apply_spectral_func(X, cos,  sk, iso)
def quad_sin(X, sk, iso):  return apply_spectral_func(X, sin,  sk, iso)
def quad_cosh(X, sk, iso): return apply_spectral_func(X, cosh, sk, iso)
def quad_sinh(X, sk, iso): return apply_spectral_func(X, sinh, sk, iso)
def quad_tan(X, sk, iso):  return apply_spectral_func(X, tan,  sk, iso)
def quad_ln(X, sk, iso):   return apply_spectral_func(X, log,  sk, iso)
def quad_sqrt(X, sk, iso): return apply_spectral_func(X, sqrt, sk, iso)

# ══════════════════════════════════════════════════════════════════
#  Print Outputs
# ══════════════════════════════════════════════════════════════════

def print_result(label: str, Q: QuadNum):
    print(f"\n  {'═'*10} Result for {label} {'═'*10}")
    for i, v in enumerate(Q.c):
        print(f"    z_{i+1} = {simplify(v)}")
    print()

def print_spectrum_result(space_key: str, mus, is_isotropic: bool):
    iso_str = "Isotropic" if is_isotropic else "Anisotropic"
    print(f"\n  {'═'*10} Spectrum Λ(X)  [{space_key} · {iso_str}] {'═'*10}")
    for i, mu in enumerate(mus, 1):
        print(f"    μ_{i} = {simplify(mu)}")
    print()

# ══════════════════════════════════════════════════════════════════
#  Menu
# ══════════════════════════════════════════════════════════════════

MENU: List[Tuple[str, str]] = [
    # (name,              category)
    ("X + Y",             "dual"),
    ("X - Y",             "dual"),
    ("X * Y",             "dual"),
    ("X / Y",             "dual"),
    ("X ^ n (integer power)", "pow"),
    ("X ^ (m/n) (fractional power)", "pow"),
    ("exp(X)",            "spec"),
    ("cos(X)",            "spec"),
    ("sin(X)",            "spec"),
    ("cosh(X)",           "spec"),
    ("sinh(X)",           "spec"),
    ("tan(X)",            "spec"),
    ("ln(X)",             "spec"),
    ("sqrt(X)",           "spec"),
    ("View spectrum Λ(X)", "spectrum"),
    ("Symplectic conjugate conj(X)", "unary"),
    ("Symplectic norm squared |X|²", "unary"),
    ("Inverse X⁻¹",       "unary"),
    ("Switch space/parameters", "switch"),
]

SECTION: Dict[str, str] = {
    "dual":     "─ Dual Operations ────────────",
    "pow":      "─ Power Operations ───────────",
    "spec":     "─ Functions (Spectral) ───────",
    "spectrum": "─ Spectrum Analysis ──────────",
    "unary":    "─ Algebraic Structure ────────",
    "switch":   "─ System ─────────────────────",
}

def print_menu(space_key, is_isotropic, params):
    iso = "Isotropic" if is_isotropic else "Anisotropic"
    a, b, g, d = params["alpha"], params["beta"], params["gamma"], params["delta"]
    print(f"\n  Current space: {space_key} ({iso})  α={a} β={b} γ={g} δ={d}")
    print("┌" + "─"*50 + "┐")
    print(f"│   0. Exit                                        │")
    prev: Optional[str] = None
    for i, (name, cat) in enumerate(MENU, 1):
        if cat != prev and cat in SECTION:
            print(f"│  {SECTION[str(cat)]:<48}│")
            prev = cat
        print(f"│   {i:2d}. {name:<45}│")
    print("└" + "─"*50 + "┘")

# ══════════════════════════════════════════════════════════════════
#  Startup: Select Space & Parameters
# ══════════════════════════════════════════════════════════════════

def choose_space() -> str:
    print("\n╔" + "═"*52 + "╗")
    print("║  4D Numbers Algebra CAS Calculator v3 (Spectral)   ║")
    print("╚" + "═"*52 + "╝")
    print("\n  Please select a 4D space:")
    keys = list(SPACES.keys())
    for i, k in enumerate(keys, 1):
        print(f"    {i}. {k}  —  {SPACES[k]['desc']}")
    while True:
        c = input("\n  Number [1-6]: ").strip()
        if c.isdigit() and 1 <= int(c) <= 6:
            return keys[int(c) - 1]
        print("  ⚠  Invalid number.")
    return ""

def choose_isotropy(space_key: str) -> Tuple[Dict[str, Any], bool]:
    s = SPACES[space_key]["signs"]
    print(f"\n  Please select space type:")
    print(f"    1. Isotropic  (α=1, β={s['beta']}, γ={s['gamma']}, δ={s['delta']})")
    print(f"    2. Anisotropic (manual input α, β, γ, δ)")
    while True:
        c = input("  Number [1/2]: ").strip()
        if c == "1":
            params = {
                "alpha": sp.Integer(1),
                "beta":  sp.Integer(s["beta"]),
                "gamma": sp.Integer(s["gamma"]),
                "delta": sp.Integer(s["delta"]),
            }
            print(f"  ✅  α=1  β={s['beta']}  γ={s['gamma']}  δ={s['delta']}")
            return params, True
        elif c == "2":
            print(f"  Hint: β{'<0' if s['beta']<0 else '>0'}  "
                  f"γ{'<0' if s['gamma']<0 else '>0'}  "
                  f"δ{'<0' if s['delta']<0 else '>0'}")
            params = {
                "alpha": input_expr("    α = "),
                "beta":  input_expr("    β = "),
                "gamma": input_expr("    γ = "),
                "delta": input_expr("    δ = "),
            }
            if simplify(params["alpha"] * params["beta"] * params["gamma"] * params["delta"]) == 0:
                print("  ⚠  Invalid parameters: α, β, γ, δ cannot be 0 (their product must not be 0).")
                continue
            return params, False
        print("  ⚠  Please enter 1 or 2.")
    return {}, False

# ══════════════════════════════════════════════════════════════════
#  Main Loop
# ══════════════════════════════════════════════════════════════════

def main_loop(params, space_key: str, is_isotropic: bool):
    while True:
        print_menu(space_key, is_isotropic, params)
        c = input("  Enter number: ").strip()
        
        if c == "0":
            print("\n  Goodbye!\n")
            break
            
        try:
            idx = int(c) - 1
            if not (0 <= idx < len(MENU)): raise ValueError
        except ValueError:
            print("  ⚠  Invalid number.")
            input("\n  Press Enter to continue...")
            continue

        op, cat = MENU[idx]

        if cat == "switch":
            space_key = choose_space()
            params, is_isotropic = choose_isotropy(space_key)
            iso_str = "Isotropic" if is_isotropic else "Anisotropic"
            print(f"\n  ✅  Switched to {space_key} ({iso_str})\n")

        elif cat == "dual":
            xc = input_quaternion("X")
            yc = input_quaternion("Y")
            X, Y = QuadNum(xc, params, space_key, is_isotropic), QuadNum(yc, params, space_key, is_isotropic)
            try:
                ops = {"X + Y": X+Y, "X - Y": X-Y, "X * Y": X*Y, "X / Y": X/Y}
                print_result(op, ops[op])
            except Exception as e:
                print(f"\n  ⚠  Calculation error: {e}")

        elif cat == "pow":
            if op == "X ^ n (integer power)":
                xc = input_quaternion("X")
                X = QuadNum(xc, params, space_key, is_isotropic)
                try: 
                    n_val = int(input("  Integer n = ").strip())
                except ValueError: 
                    print("  ⚠  n must be an integer.")
                    input("\n  Press Enter to continue...")
                    continue
                print_result(f"X^{n_val}", X ** n_val)
            
            elif op == "X ^ (m/n) (fractional power)":
                xc = input_quaternion("X")
                X = QuadNum(xc, params, space_key, is_isotropic)
                try: 
                    p_val = int(input("  Numerator m = ").strip())
                    q_val = int(input("  Denominator n (n > 0) = ").strip())
                    if q_val <= 0:
                        raise ValueError("n must be positive")
                except ValueError: 
                    print("  ⚠  Invalid input for m/n (must be integers, n > 0).")
                    input("\n  Press Enter to continue...")
                    continue
                
                print(f"\n  Computing all {q_val * q_val} branches of X^({p_val}/{q_val})...")
                mu1, mu2, mu3, mu4 = compute_spectrum(X.c, X.p, space_key, is_isotropic)
                
                r1, theta1 = simplify(Abs(mu1)), simplify(sp.arg(mu1))
                r3, theta3 = simplify(Abs(mu3)), simplify(sp.arg(mu3))
                
                is_mu1_zero = (simplify(r1) == 0)
                is_mu3_zero = (simplify(r3) == 0)
                
                solutions = []
                for k in range(q_val):
                    for m_idx in range(q_val):
                        if is_mu1_zero:
                            lam1 = sp.S.Zero
                        else:
                            lam1 = (r1 ** sp.Rational(p_val, q_val)) * exp(I * p_val * (theta1 + 2 * pi * k) / q_val)
                            
                        if is_mu3_zero:
                            lam3 = sp.S.Zero
                        else:
                            lam3 = (r3 ** sp.Rational(p_val, q_val)) * exp(I * p_val * (theta3 + 2 * pi * m_idx) / q_val)
                            
                        lam1 = simplify(lam1)
                        lam3 = simplify(lam3)
                        W = spectrum_to_quad(lam1, lam3, X.p, space_key, is_isotropic)
                        solutions.append({'quad': W, 'k': k, 'm': m_idx})
                
                print("\n  " + "=" * 70)
                for idx, sol in enumerate(solutions):
                    k = sol['k']
                    m_idx = sol['m']
                    marker = "*" if idx == 0 else " "
                    print(f"  {marker} Branch [k = {k}, m = {m_idx}]")
                    for i, v in enumerate(sol['quad'].c):
                        print(f"    z_{i+1} = {simplify(v)}")
                    print()
                print("  " + "-" * 70)
                print(f"  {len(solutions)} solutions in total.")
                print("  " + "=" * 70)

        elif cat == "spec":
            xc = input_quaternion("X")
            X = QuadNum(xc, params, space_key, is_isotropic)
            print(f"\n  Computing {op} using spectral method...")
            try:
                fn = {"exp(X)": quad_exp, "cos(X)": quad_cos, "sin(X)": quad_sin,
                      "cosh(X)": quad_cosh, "sinh(X)": quad_sinh, "tan(X)": quad_tan,
                      "ln(X)": quad_ln, "sqrt(X)": quad_sqrt}
                print_result(op, fn[op](X, space_key, is_isotropic))
            except Exception as e:
                print(f"\n  ⚠  Calculation error: {e}")

        elif cat == "spectrum":
            xc = input_quaternion("X")
            X = QuadNum(xc, params, space_key, is_isotropic)
            mus = compute_spectrum(X.c, X.p, space_key, is_isotropic)
            print_spectrum_result(space_key, mus, is_isotropic)

        elif cat == "unary":
            xc = input_quaternion("X")
            X = QuadNum(xc, params, space_key, is_isotropic)
            try:
                if op == "Symplectic conjugate conj(X)":  
                    print_result("conj(X)", X.symp_conj())
                elif op == "Symplectic norm squared |X|²":
                    print(f"\n  {'═'*10} |X|² {'═'*10}\n    |X|² = {simplify(X.symp_norm_sq())}\n")
                elif op == "Inverse X⁻¹":         
                    print_result("X⁻¹", X.inv())
            except Exception as e:
                print(f"\n  ⚠  Calculation error: {e}")

        input("\n  Press Enter to continue...")

# ══════════════════════════════════════════════════════════════════
#  Entry Point
# ══════════════════════════════════════════════════════════════════

def main():
    space_key = choose_space()
    params, is_isotropic = choose_isotropy(space_key)
    iso_str = "Isotropic" if is_isotropic else "Anisotropic"
    print(f"\n  ✅  Entered {space_key} ({iso_str})\n")
    main_loop(params, space_key, is_isotropic)

if __name__ == "__main__":
    main()