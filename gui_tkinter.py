"""
4D Numbers Algebra System CAS — Tkinter Edition
Requires: sympy (pip install sympy)
Run: python app_tkinter.py
"""

import tkinter as tk
from tkinter import ttk, font as tkfont
import sympy as sp
from sympy import simplify, I, pi, Abs, exp

# Импорт алгоритмического ядра (тот же модуль, что и в Streamlit-версии)
from backend import (
    SPACES, QuadNum, parse_expr,
    compute_spectrum, spectrum_to_quad,
    quad_exp, quad_cos, quad_sin, quad_cosh, quad_sinh, quad_tan, quad_ln
)

# ═══════════════════════════════════════════════════════════════════
# Цветовая палитра (точно из CSS Streamlit-версии)
# ═══════════════════════════════════════════════════════════════════
CLR = {
    "bg":           "#ffffff",
    "fg":           "#111111",
    "h_fg":         "#222222",
    "output_bg":    "#f4f5f7",
    "output_fg":    "#2f3640",
    "output_bd":    "#dcdde1",
    "btn_bg":       "#f0f2f6",
    "btn_bd":       "#dcdde1",
    "btn_fg":       "#2f3640",
    "btn_hover":    "#e1e4e8",
    "input_bd":     "#cccccc",
    "lf_bd":        "#dcdde1",
    "caption":      "#777777",
    "sep":          "#eeeeee",
    # теги вывода
    "ok":           "#0059b3",
    "err":          "#c23616",
    "head":         "#27ae60",
}

OPERATIONS = [
    "X + Y", "X - Y", "X * Y", "X / Y",
    "X ^ n", "X ^ (m/n)",
    "exp(X)", "ln(X)", "cos(X)", "sin(X)",
    "cosh(X)", "sinh(X)", "tan(X)",
    "Spectrum S(X)", "conj(X)", "|X|² (Norm)", "X⁻¹ (Inv)",
]

DUAL_OPS = ("X + Y", "X - Y", "X * Y", "X / Y")

FN_MAP = {
    "exp(X)":  quad_exp,
    "cos(X)":  quad_cos,
    "sin(X)":  quad_sin,
    "cosh(X)": quad_cosh,
    "sinh(X)": quad_sinh,
    "tan(X)":  quad_tan,
    "ln(X)":   quad_ln,
}


# ═══════════════════════════════════════════════════════════════════
# Главный класс приложения
# ═══════════════════════════════════════════════════════════════════
class CASApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("4D Number Calculator")
        self.configure(bg=CLR["bg"])
        self.resizable(True, True)
        self.minsize(860, 680)

        self._build_fonts()
        self._build_style()
        self._build_ui()

        # Инициализируем видимость секции переменных
        self._on_op_change()

    # ── Шрифты ──────────────────────────────────────────────────
    def _build_fonts(self):
        self.f_title   = tkfont.Font(family="Segoe UI", size=18, weight="normal")
        self.f_caption = tkfont.Font(family="Segoe UI", size=9)
        self.f_label   = tkfont.Font(family="Segoe UI", size=10, weight="bold")
        self.f_normal  = tkfont.Font(family="Segoe UI", size=10)
        self.f_mono    = tkfont.Font(family="Consolas",  size=10)
        self.f_section = tkfont.Font(family="Segoe UI", size=10, weight="bold")
        self.f_btn     = tkfont.Font(family="Segoe UI", size=10, weight="bold")

    # ── ttk-стиль ───────────────────────────────────────────────
    def _build_style(self):
        s = ttk.Style(self)
        s.theme_use("clam")

        s.configure("TFrame",       background=CLR["bg"])
        s.configure("TLabel",       background=CLR["bg"], foreground=CLR["fg"],
                    font=("Segoe UI", 10))
        s.configure("Caption.TLabel", foreground=CLR["caption"],
                    background=CLR["bg"], font=("Segoe UI", 9))
        s.configure("Head.TLabel",  foreground=CLR["h_fg"],
                    background=CLR["bg"], font=("Segoe UI", 10, "bold"))
        s.configure("Section.TLabel", foreground=CLR["h_fg"],
                    background=CLR["bg"], font=("Segoe UI", 10, "bold"))

        s.configure("TLabelframe", background=CLR["bg"],
                    bordercolor=CLR["lf_bd"], relief="solid", borderwidth=1)
        s.configure("TLabelframe.Label", background=CLR["bg"],
                    foreground=CLR["h_fg"], font=("Segoe UI", 10, "bold"))

        s.configure("TEntry",
                    fieldbackground="#ffffff",
                    foreground=CLR["fg"],
                    bordercolor=CLR["input_bd"],
                    relief="flat")
        s.map("TEntry", bordercolor=[("focus", "#888888")])

        s.configure("TCombobox",
                    fieldbackground="#ffffff",
                    foreground=CLR["fg"],
                    background=CLR["btn_bg"],
                    selectbackground="#ffffff",
                    selectforeground=CLR["fg"])

        s.configure("TRadiobutton",
                    background=CLR["bg"],
                    foreground=CLR["fg"],
                    font=("Segoe UI", 10))
        s.map("TRadiobutton", background=[("active", CLR["bg"])])

        # Кнопки — плоский вид как в Streamlit
        s.configure("Flat.TButton",
                    background=CLR["btn_bg"],
                    foreground=CLR["btn_fg"],
                    bordercolor=CLR["btn_bd"],
                    relief="solid",
                    borderwidth=1,
                    font=("Segoe UI", 10, "bold"),
                    padding=(10, 6))
        s.map("Flat.TButton",
              background=[("active", CLR["btn_hover"])],
              bordercolor=[("active", "#b2bec3")])

        s.configure("Run.TButton",
                    background=CLR["btn_bg"],
                    foreground=CLR["btn_fg"],
                    bordercolor=CLR["btn_bd"],
                    relief="solid",
                    borderwidth=1,
                    font=("Segoe UI", 10, "bold"),
                    padding=(10, 6))
        s.map("Run.TButton",
              background=[("active", CLR["btn_hover"])])

        s.configure("TSeparator", background=CLR["sep"])

    # ── Построение интерфейса ────────────────────────────────────
    def _build_ui(self):
        # Скроллируемый канвас — вся страница
        outer = ttk.Frame(self)
        outer.pack(fill="both", expand=True)

        canvas = tk.Canvas(outer, bg=CLR["bg"], highlightthickness=0)
        vbar   = ttk.Scrollbar(outer, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=vbar.set)
        vbar.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)

        self.page = ttk.Frame(canvas)
        self._canvas_window = canvas.create_window((0, 0), window=self.page, anchor="nw")

        def _on_frame_configure(e):
            canvas.configure(scrollregion=canvas.bbox("all"))
        def _on_canvas_resize(e):
            canvas.itemconfig(self._canvas_window, width=e.width)

        self.page.bind("<Configure>", _on_frame_configure)
        canvas.bind("<Configure>", _on_canvas_resize)
        canvas.bind_all("<MouseWheel>",
                        lambda e: canvas.yview_scroll(-1 * (e.delta // 120), "units"))

        pad = {"padx": 20, "pady": (6, 4)}

        # ── Заголовок ──────────────────────────────────────────
        ttk.Label(self.page, text="4D Number Calculator",
                  font=self.f_title, foreground=CLR["h_fg"],
                  background=CLR["bg"]).pack(anchor="w", **pad)
        ttk.Label(self.page,
                  text="Spectral decomposition and transcendental functions "
                       "for four-dimensional number spaces",
                  style="Caption.TLabel").pack(anchor="w", padx=20, pady=(0, 10))

        ttk.Separator(self.page, orient="horizontal").pack(fill="x", padx=20, pady=4)

        # ── Step 1: Space Configuration ────────────────────────
        lf1 = ttk.LabelFrame(self.page, text="Step 1: Space Configuration", padding=12)
        lf1.pack(fill="x", padx=20, pady=6)

        ttk.Label(lf1, text="Select Space", style="Head.TLabel").pack(anchor="w")
        self.space_var = tk.StringVar(value=list(SPACES.keys())[0])
        cb = ttk.Combobox(lf1, textvariable=self.space_var,
                          values=list(SPACES.keys()),
                          state="readonly", width=40, font=("Consolas", 10))
        cb.pack(anchor="w", pady=(4, 0))

        # ── Step 2: Select Operation ───────────────────────────
        lf2 = ttk.LabelFrame(self.page, text="Step 2: Select Operation", padding=12)
        lf2.pack(fill="x", padx=20, pady=6)

        self.op_var = tk.StringVar(value=OPERATIONS[0])
        self._build_op_radios(lf2)

        # ── Step 3: Variables Input ────────────────────────────
        self.lf3 = ttk.LabelFrame(self.page, text="Step 3: Variables Input", padding=12)
        self.lf3.pack(fill="x", padx=20, pady=6)
        self._build_var_section(self.lf3)

        # ── Кнопки ────────────────────────────────────────────
        btn_row = ttk.Frame(self.page)
        btn_row.pack(fill="x", padx=20, pady=(8, 4))

        ttk.Button(btn_row, text="Execute Calculation",
                   style="Run.TButton",
                   command=self._run).pack(side="left", fill="x", expand=True, padx=(0, 6))
        ttk.Button(btn_row, text="Clear Output",
                   style="Flat.TButton",
                   command=self._clear_output).pack(side="left")

        # ── Output Log ────────────────────────────────────────
        ttk.Label(self.page, text="Output Log",
                  font=("Segoe UI", 13), foreground=CLR["h_fg"],
                  background=CLR["bg"]).pack(anchor="w", padx=20, pady=(10, 2))

        out_frame = ttk.Frame(self.page)
        out_frame.pack(fill="x", padx=20, pady=(0, 6))

        self.output = tk.Text(
            out_frame,
            bg=CLR["output_bg"],
            fg=CLR["output_fg"],
            relief="flat",
            bd=0,
            highlightthickness=1,
            highlightbackground=CLR["output_bd"],
            font=self.f_mono,
            wrap="word",
            height=14,
            state="disabled",
            padx=14, pady=10,
            cursor="arrow",
        )
        out_sb = ttk.Scrollbar(out_frame, orient="vertical",
                               command=self.output.yview)
        self.output.configure(yscrollcommand=out_sb.set)
        out_sb.pack(side="right", fill="y")
        self.output.pack(side="left", fill="both", expand=True)

        # Теги цвета вывода
        self.output.tag_configure("ok",   foreground=CLR["ok"])
        self.output.tag_configure("err",  foreground=CLR["err"])
        self.output.tag_configure("head", foreground=CLR["head"], font=self.f_mono)

        # Начальная подсказка
        self._show_placeholder()

        ttk.Separator(self.page, orient="horizontal").pack(fill="x", padx=20, pady=8)
        ttk.Label(self.page,
                  text="Author: Aset Eskhat · Stat-201",
                  style="Caption.TLabel").pack(anchor="w", padx=20, pady=(0, 12))

    def _build_op_radios(self, parent):
        """Рисует радиокнопки операций в несколько строк."""
        # 4 столбца, как horizontal radio в Streamlit
        cols_per_row = 4
        rows_frame = ttk.Frame(parent)
        rows_frame.pack(fill="x")
        for idx, op in enumerate(OPERATIONS):
            r, c = divmod(idx, cols_per_row)
            rb = ttk.Radiobutton(
                rows_frame, text=op,
                variable=self.op_var, value=op,
                command=self._on_op_change,
            )
            rb.grid(row=r, column=c, sticky="w", padx=8, pady=2)
        for c in range(cols_per_row):
            rows_frame.columnconfigure(c, weight=1)

    def _build_var_section(self, parent):
        """Строит секцию ввода переменных (динамическая)."""
        # X вектор — постоянный
        self.x_vars = [tk.StringVar(value="0") for _ in range(4)]
        self.y_vars = [tk.StringVar(value="0") for _ in range(4)]
        self.n_var  = tk.StringVar(value="2")
        self.fm_var = tk.StringVar(value="1")
        self.fn_var = tk.StringVar(value="2")

        # Блок X
        self.x_frame = ttk.Frame(parent)
        self.x_frame.pack(fill="x")
        ttk.Label(self.x_frame, text="X Component Vector  (z₁, z₂, z₃, z₄)",
                  style="Head.TLabel").pack(anchor="w", pady=(0, 4))
        x_cols_frame = ttk.Frame(self.x_frame)
        x_cols_frame.pack(fill="x")
        self.x_entries = []
        for i in range(4):
            f = ttk.Frame(x_cols_frame)
            f.grid(row=0, column=i, sticky="ew", padx=(0, 8))
            ttk.Label(f, text=f"x_{i+1}").pack(anchor="w")
            e = ttk.Entry(f, textvariable=self.x_vars[i], width=16,
                          font=("Consolas", 10))
            e.pack(fill="x")
            self.x_entries.append(e)
        for c in range(4):
            x_cols_frame.columnconfigure(c, weight=1)

        # Блок Y
        self.y_frame = ttk.Frame(parent)
        ttk.Label(self.y_frame, text="Y Component Vector  (z₁, z₂, z₃, z₄)",
                  style="Head.TLabel").pack(anchor="w", pady=(8, 4))
        y_cols_frame = ttk.Frame(self.y_frame)
        y_cols_frame.pack(fill="x")
        for i in range(4):
            f = ttk.Frame(y_cols_frame)
            f.grid(row=0, column=i, sticky="ew", padx=(0, 8))
            ttk.Label(f, text=f"y_{i+1}").pack(anchor="w")
            ttk.Entry(f, textvariable=self.y_vars[i], width=16,
                      font=("Consolas", 10)).pack(fill="x")
        for c in range(4):
            y_cols_frame.columnconfigure(c, weight=1)

        # Блок целого показателя
        self.n_frame = ttk.Frame(parent)
        ttk.Label(self.n_frame, text="Integer Exponent (n)",
                  style="Head.TLabel").pack(anchor="w", pady=(8, 4))
        ttk.Entry(self.n_frame, textvariable=self.n_var, width=14,
                  font=("Consolas", 10)).pack(anchor="w")

        # Блок дробного показателя
        self.fp_frame = ttk.Frame(parent)
        ttk.Label(self.fp_frame, text="Fractional Exponent  (m/n)",
                  style="Head.TLabel").pack(anchor="w", pady=(8, 4))
        fp_cols = ttk.Frame(self.fp_frame)
        fp_cols.pack(anchor="w")
        ttk.Label(fp_cols, text="Numerator (m)").grid(row=0, column=0, sticky="w", padx=(0, 20))
        ttk.Label(fp_cols, text="Denominator (n)").grid(row=0, column=1, sticky="w")
        ttk.Entry(fp_cols, textvariable=self.fm_var, width=10,
                  font=("Consolas", 10)).grid(row=1, column=0, padx=(0, 20))
        ttk.Entry(fp_cols, textvariable=self.fn_var, width=10,
                  font=("Consolas", 10)).grid(row=1, column=1)

    def _on_op_change(self, *_):
        op = self.op_var.get()
        # Скрываем все дополнительные блоки
        for w in (self.y_frame, self.n_frame, self.fp_frame):
            w.pack_forget()
        # Показываем нужный
        if op in DUAL_OPS:
            self.y_frame.pack(fill="x")
        elif op == "X ^ n":
            self.n_frame.pack(fill="x")
        elif op == "X ^ (m/n)":
            self.fp_frame.pack(fill="x")

    # ── Вспомогательные методы вывода ──────────────────────────
    def _show_placeholder(self):
        self.output.configure(state="normal")
        self.output.delete("1.0", "end")
        self.output.insert("end", "Awaiting execution...", "placeholder")
        self.output.tag_configure("placeholder", foreground=CLR["caption"])
        self.output.configure(state="disabled")

    def _clear_output(self):
        self._show_placeholder()

    def _append(self, text: str, tag: str = "ok"):
        self.output.configure(state="normal")
        self.output.insert("end", text + "\n", tag)
        self.output.see("end")
        self.output.configure(state="disabled")

    def _append_head(self, text: str):
        self._append(text, "head")

    def _append_err(self, text: str):
        self._append(text, "err")

    def _reset_output(self):
        self.output.configure(state="normal")
        self.output.delete("1.0", "end")
        self.output.configure(state="disabled")

    # ── Логика вычислений (порт из Streamlit-версии) ────────────
    def _get_params(self):
        sk = self.space_var.get()
        signs = SPACES[sk]["signs"]
        return {
            "alpha": parse_expr("1"),
            "beta":  parse_expr(str(signs["beta"])),
            "gamma": parse_expr(str(signs["gamma"])),
            "delta": parse_expr(str(signs["delta"])),
        }

    def _get_X(self, p, sk):
        comps = [parse_expr(v.get() or "0") for v in self.x_vars]
        return QuadNum(comps, p, sk, True)

    def _get_Y(self, p, sk):
        comps = [parse_expr(v.get() or "0") for v in self.y_vars]
        return QuadNum(comps, p, sk, True)

    def _run(self):
        op = self.op_var.get()
        sk = self.space_var.get()
        iso = True

        self._reset_output()

        try:
            p = self._get_params()
            self._append_head(f"Operation: [{op}] | Space: {sk}")

            X = self._get_X(p, sk)

            if op in DUAL_OPS:
                Y = self._get_Y(p, sk)
                ans = (X + Y) if op == "X + Y" else \
                      (X - Y) if op == "X - Y" else \
                      (X * Y) if op == "X * Y" else (X / Y)
                for i, c in enumerate(ans.c):
                    self._append(f"  z_{i+1} = {simplify(c)}")

            elif op == "X ^ n":
                n_val = int(self.n_var.get())
                ans = X ** n_val
                for i, c in enumerate(ans.c):
                    self._append(f"  z_{i+1} = {simplify(c)}")

            elif op == "X ^ (m/n)":
                pm = int(self.fm_var.get())
                qn = int(self.fn_var.get())
                if qn <= 0:
                    raise ValueError("Denominator n must be > 0.")

                mu1, mu2, mu3, mu4 = compute_spectrum(X.c, X.p, sk, iso)
                r1, theta1 = simplify(Abs(mu1)), simplify(sp.arg(mu1))
                r3, theta3 = simplify(Abs(mu3)), simplify(sp.arg(mu3))
                is_m1z = simplify(r1) == 0
                is_m3z = simplify(r3) == 0

                count = 0
                self._append(f"  Computing {qn*qn} branches...")
                for k in range(qn):
                    for m_idx in range(qn):
                        lam1 = sp.S.Zero if is_m1z else \
                               (r1 ** sp.Rational(pm, qn)) * \
                               exp(I * pm * (theta1 + 2 * pi * k) / qn)
                        lam3 = sp.S.Zero if is_m3z else \
                               (r3 ** sp.Rational(pm, qn)) * \
                               exp(I * pm * (theta3 + 2 * pi * m_idx) / qn)
                        W = spectrum_to_quad(simplify(lam1), simplify(lam3),
                                             X.p, sk, iso)
                        self._append_head(f"\n  * Branch [k={k}, m={m_idx}]")
                        for i, v in enumerate(W.c):
                            self._append(f"    z_{i+1} = {simplify(v)}")
                        count += 1
                self._append(f"  Total solutions: {count}")

            elif op in FN_MAP:
                ans = FN_MAP[op](X, sk, iso)
                for i, c in enumerate(ans.c):
                    self._append(f"  z_{i+1} = {simplify(c)}")

            elif op == "Spectrum S(X)":
                mus = compute_spectrum(X.c, X.p, sk, iso)
                for i, mu in enumerate(mus, 1):
                    self._append(f"  μ_{i} = {mu}")

            elif op == "conj(X)":
                ans = X.symp_conj()
                for i, c in enumerate(ans.c):
                    self._append(f"  z_{i+1} = {simplify(c)}")

            elif op == "|X|² (Norm)":
                self._append(f"  |X|² = {simplify(X.symp_norm_sq())}")

            elif op == "X⁻¹ (Inv)":
                ans = X.inv()
                for i, c in enumerate(ans.c):
                    self._append(f"  z_{i+1} = {simplify(c)}")

            else:
                self._append_err("Operation not supported.")

        except Exception as exc:
            self._append_err(f"  Error: {exc}")


# ═══════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    app = CASApp()
    app.mainloop()