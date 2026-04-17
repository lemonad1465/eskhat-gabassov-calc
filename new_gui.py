"""
4D Numbers Algebra System CAS — Streamlit Edition
Requires: streamlit, sympy
Run: streamlit run app.py
"""

import streamlit as st
import sympy as sp
from sympy import simplify, I, pi, Abs, exp



hide_style = """
    <style>
    header[data-testid="stHeader"] { display: none !important; }
    .stAppDeployButton { display: none !important; }
    
    div[data-testid="stStatusWidget"] { display: none !important; }

    footer { visibility: hidden !important; }
    #MainMenu { visibility: hidden !important; }
    
    .block-container { 
        padding-top: 1rem !important; 
        padding-bottom: 0rem !important; 
    }
    </style>
"""

# Применяем стили 
st.markdown(hide_style, unsafe_allow_html=True)
# --- КОНЕЦ БЛОКА СКРЫТИЯ ИНТЕРФЕЙСА ---

# Импорт алгоритмического ядра
from backend import (
    SPACES, QuadNum, parse_expr,
    compute_spectrum, spectrum_to_quad,
    quad_exp, quad_cos, quad_sin, quad_cosh, quad_sinh, quad_tan, quad_ln
)

# ═══════════════════════════════════════════════════════════════════
# Инициализация состояния
# ═══════════════════════════════════════════════════════════════════
st.set_page_config(
    page_title="4D Numbers CAS",
    layout="wide",
)


CSS = """
<style>
[data-testid="stAppViewContainer"] { background: #ffffff; color: #111111; }
h1, h2, h3 { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; font-weight: 400; color: #222222 !important; }

.output-box {
    border-radius: 4px;
    padding: 1rem 1.4rem;
    font-family: 'Consolas', 'Courier New', monospace;
    font-size: 0.9rem;
    line-height: 1.6;
    min-height: 120px;
    white-space: pre-wrap;
    word-break: break-word;
    background: #f4f5f7; 
    border: 1px solid #dcdde1; 
    color: #2f3640;
}

.output-box .ok   { color: #0059b3; }
.output-box .err  { color: #c23616; }
.output-box .head { color: #27ae60; font-weight: bold; }

div[data-testid="stButton"] > button {
    border-radius: 4px;
    font-weight: 600;
    transition: all 0.2s;
    background: #f0f2f6; 
    border: 1px solid #dcdde1; 
    color: #2f3640;
}

div[data-testid="stButton"] > button:hover { 
    background: #e1e4e8; 
    border-color: #b2bec3; 
}

input, select, textarea { 
    background: #ffffff !important; 
    color: #111111 !important; 
    border: 1px solid #cccccc !important; 
}

label { color: #333333 !important; font-weight: 500 !important; }
hr { border-color: #eeeeee; }
</style>
"""

st.markdown(CSS, unsafe_allow_html=True)

# ═══════════════════════════════════════════════════════════════════
# Интерфейс
# ═══════════════════════════════════════════════════════════════════
st.title("4D Number Calculator")
st.caption("Spectral decomposition and transcendental functions for four-dimensional number spaces")

# ═══════════════════════════════════════════════════════════════════
# Вспомогательные функции
# ═══════════════════════════════════════════════════════════════════

def get_params(space_key: str) -> dict:
    """Возвращает параметры пространства на основе выбранного ключа. Коэффициенты строго фиксированы."""
    signs = SPACES[space_key]["signs"]
    return {
        "alpha": parse_expr("1"),
        "beta": parse_expr(str(signs["beta"])),
        "gamma": parse_expr(str(signs["gamma"])),
        "delta": parse_expr(str(signs["delta"]))
    }

def get_X(p, sk):
    comps = [parse_expr(st.session_state.get(f"x{i}", "0") or "0") for i in range(4)]
    return QuadNum(comps, p, sk, True) # iso жестко задан как True

def get_Y(p, sk):
    comps = [parse_expr(st.session_state.get(f"y{i}", "0") or "0") for i in range(4)]
    return QuadNum(comps, p, sk, True)

# ═══════════════════════════════════════════════════════════════════
#  Конфигурация
# ═══════════════════════════════════════════════════════════════════
with st.expander("Step 1: Space Configuration", expanded=True):
    space_key = st.selectbox("Select Space", list(SPACES.keys()), key="space_key")

# ═══════════════════════════════════════════════════════════════════
# Выбор операции
# ═══════════════════════════════════════════════════════════════════
OPERATIONS = [
    "X + Y", "X - Y", "X * Y", "X / Y",
    "X ^ n", "X ^ (m/n)",
    "exp(X)", "ln(X)", "cos(X)", "sin(X)",
    "cosh(X)", "sinh(X)", "tan(X)",
    "Spectrum S(X)", "conj(X)", "|X|² (Norm)", "X⁻¹ (Inv)"
]

with st.expander("Step 2: Select Operation", expanded=True):
    selected_op = st.radio(
        "Operation", OPERATIONS, key="selected_op", horizontal=True, label_visibility="collapsed"
    )

# ═══════════════════════════════════════════════════════════════════
# Ввод переменных
# ═══════════════════════════════════════════════════════════════════
with st.expander("Step 3: Variables Input", expanded=True):
    st.markdown("**X Component Vector** (z₁, z₂, z₃, z₄)")
    x_cols = st.columns(4)
    for i in range(4):
        with x_cols[i]:
            st.text_input(f"x_{i+1}", value="0", key=f"x{i}")

    DUAL_OPS = ("X + Y", "X - Y", "X * Y", "X / Y")
    if selected_op in DUAL_OPS:
        st.markdown("**Y Component Vector** (z₁, z₂, z₃, z₄)")
        y_cols = st.columns(4)
        for i in range(4):
            with y_cols[i]:
                st.text_input(f"y_{i+1}", value="0", key=f"y{i}")

    elif selected_op == "X ^ n":
        st.text_input("Integer Exponent (n)", value="2", key="n_int")

    elif selected_op == "X ^ (m/n)":
        pm_col, pn_col = st.columns(2)
        with pm_col:
            st.text_input("Numerator (m)", value="1", key="fp_m")
        with pn_col:
            st.text_input("Denominator (n)", value="2", key="fp_n")

# ═══════════════════════════════════════════════════════════════════
# Вычисления и Логика
# ═══════════════════════════════════════════════════════════════════
col_run, col_clear = st.columns([5, 1])
with col_run:
    run = st.button("Execute Calculation", use_container_width=True)
with col_clear:
    if st.button("Clear Output", use_container_width=True):
        st.session_state["output"] = ""

if "output" not in st.session_state:
    st.session_state["output"] = ""

def append(text: str, tag: str = "ok"):
    st.session_state["output"] += f'<span class="{tag}">{text}</span>\n'

def append_head(text: str):
    append(text, "head")

def append_err(text: str):
    append(text, "err")

FN_MAP = {
    "exp(X)":  quad_exp,
    "cos(X)":  quad_cos,
    "sin(X)":  quad_sin,
    "cosh(X)": quad_cosh,
    "sinh(X)": quad_sinh,
    "tan(X)":  quad_tan,
    "ln(X)":   quad_ln
}

if run:
    op = st.session_state["selected_op"]
    sk = st.session_state["space_key"]
    iso = True

    try:
        p = get_params(sk)
        append_head(f"Operation: [{op}] | Space: {sk}")
        
        X = get_X(p, sk)

        if op in DUAL_OPS:
            Y = get_Y(p, sk)
            ans = (X + Y) if op == "X + Y" else \
                  (X - Y) if op == "X - Y" else \
                  (X * Y) if op == "X * Y" else (X / Y)
            for i, c in enumerate(ans.c):
                append(f"  z_{i+1} = {simplify(c)}")

        elif op == "X ^ n":
            n_val = int(st.session_state.get("n_int", "2"))
            ans = X ** n_val
            for i, c in enumerate(ans.c):
                append(f"  z_{i+1} = {simplify(c)}")

        elif op == "X ^ (m/n)":
            pm = int(st.session_state.get("fp_m", "1"))
            qn = int(st.session_state.get("fp_n", "2"))
            if qn <= 0:
                raise ValueError("Знаменатель n должен быть больше нуля.")

            mu1, mu2, mu3, mu4 = compute_spectrum(X.c, X.p, sk, iso)
            r1, theta1 = simplify(Abs(mu1)), simplify(sp.arg(mu1))
            r3, theta3 = simplify(Abs(mu3)), simplify(sp.arg(mu3))
            is_m1z = simplify(r1) == 0
            is_m3z = simplify(r3) == 0

            count = 0
            append(f"  Computing {qn*qn} branches...")
            for k in range(qn):
                for m_idx in range(qn):
                    lam1 = sp.S.Zero if is_m1z else \
                           (r1 ** sp.Rational(pm, qn)) * exp(I * pm * (theta1 + 2 * pi * k) / qn)
                    lam3 = sp.S.Zero if is_m3z else \
                           (r3 ** sp.Rational(pm, qn)) * exp(I * pm * (theta3 + 2 * pi * m_idx) / qn)
                    W = spectrum_to_quad(simplify(lam1), simplify(lam3), X.p, sk, iso)
                    append(f"\n  * Branch [k={k}, m={m_idx}]", "head")
                    for i, v in enumerate(W.c):
                        append(f"    z_{i+1} = {simplify(v)}")
                    count += 1
            append(f"  Total solutions: {count}")

        elif op in FN_MAP:
            ans = FN_MAP[op](X, sk, iso)
            for i, c in enumerate(ans.c):
                append(f"  z_{i+1} = {simplify(c)}")

        elif op == "Spectrum S(X)":
            mus = compute_spectrum(X.c, X.p, sk, iso)
            for i, mu in enumerate(mus, 1):
                append(f"  μ_{i} = {mu}")

        elif op == "conj(X)":
            ans = X.symp_conj()
            for i, c in enumerate(ans.c):
                append(f"  z_{i+1} = {simplify(c)}")

        elif op == "|X|² (Norm)":
            append(f"  |X|² = {simplify(X.symp_norm_sq())}")

        elif op == "X⁻¹ (Inv)":
            ans = X.inv()
            for i, c in enumerate(ans.c):
                append(f"  z_{i+1} = {simplify(c)}")

        else:
            append_err("Operation not supported.")

    except Exception as exc:
        append_err(f"  Error: {exc}")

# ═══════════════════════════════════════════════════════════════════
# Вывод результатов
# ═══════════════════════════════════════════════════════════════════
st.markdown("### Output Log")
output_content = st.session_state["output"] if st.session_state["output"] else "<span style='color:#777777'>Awaiting execution...</span>"
st.markdown(f'<div class="output-box">{output_content}</div>', unsafe_allow_html=True)

st.divider()
st.caption("Author: Aset Eskhat;Stat - 201")