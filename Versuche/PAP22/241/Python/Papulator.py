'''
In diesem Dokument werden wichtige Methoden für das Auswerten von Versucehn gesammelt und Dokumentiert.
'''

# ==================================================
#   Libaries
# ==================================================
import numpy as np
from decimal import Decimal, ROUND_HALF_UP
import sympy as sp
from IPython.display import display, Math, Latex, HTML
import string
from enum import Enum


# ==================================================
#   Default basic Methods
# ==================================================
# --------------------------------------------------
#   Rounding erros to Sig Position
# --------------------------------------------------
def round_err_sig(err, forceRoundDig: int = None):
    """
    Rounds a number according to custom significant-digit rules and returns a formatted string representation.

    
    Parameters :
    -----------
    **err** : float
        The number to be rounded.

    **forceRoundDig** : int, optional
        If provided, forces the number of decimal places used in rounding, overriding the automatic significant-digit logic.

    
    Returns :
    --------
    **results** : str
        The rounded number as a string, formatted to preserve or suppress trailing zeros.

    **rounded** : float
        The rounded number as a float for further calculations if needed. 

    **decimal** : int
        The rounding position. This can be used to round other values to the same length.
    """

    err = Decimal(str(err))  # Avoid float artifacts

    if err == 0:
        return "0", Decimal("0"), 0

    sign = "-" if err < 0 else ""
    x_abs = abs(err)

    # equivalent to floor(log10(x))
    order = x_abs.adjusted()

    # first significant digit
    first_digit = int(x_abs.scaleb(-order).to_integral_value(rounding=ROUND_HALF_UP))

    # determine significant digits
    if first_digit in (1, 2):
        sig_digits = 2
    else:
        sig_digits = 1

    decimals = sig_digits - order - 1

    # override if forced
    if forceRoundDig is not None:
        decimals = forceRoundDig

    # create quantization target
    quant = Decimal("1e{}".format(-decimals))

    rounded = x_abs.quantize(quant, rounding=ROUND_HALF_UP)

    # detect carry
    carry_happened = (
        x_abs < 1 and
        rounded != 0 and
        rounded.adjusted() > x_abs.adjusted()
    )

    # formatting
    if decimals > 0:
        result = format(rounded, f".{decimals}f")
    else:
        result = str(rounded.quantize(Decimal("1"), rounding=ROUND_HALF_UP))

    # force trailing zero
    if carry_happened and decimals > 0:
        result = format(rounded, f".{decimals}f")

    return result, rounded, decimals


# --------------------------------------------------
#   Rounding erros and a Value to Sig Position
# --------------------------------------------------
def round_sig_digs(err: float, val: float = None):
    '''
    Function to round a value and its uncertainty to the significant digits of the error.

    Parameters :
    -----------
    **err** : float
        Error to be rounded to significant digits

    **val** : float
        Value corresponding to the error. Will be rounded to the same decimal precision as the error.

    Returns :
    --------
    **rounded_val** : float
        Rounded value (numeric)

    **rounded_err** : float
        Rounded error (numeric)

    **rounded_combined_str** : str
        Formatted string "value ± error" with matching decimal places
    '''

    err_str, rounded_err, decimals = round_err_sig(err)

    # round value numerically
    rounded_val = round(val, decimals)

    # format value to match error precision
    if decimals > 0:
        val_str = f"{rounded_val:.{decimals}f}"
    else:
        val_str = str(int(rounded_val))

    rounded_combined_str = f"{val_str} \\pm {err_str}"

    return rounded_val, rounded_err, rounded_combined_str


# --------------------------------------------------
#   Gaußische Fehlerfortpfanzung
# --------------------------------------------------
def gff(func, errPronePar):
    """
    Kann die Fehlerformel einer gegebenen Gleichung nach gaußischer Fehlerfortpflanzung bestimmen.

    Parameters
    ----------
    **func** : sympy function
        Funktion dessen Fehler bestimmt werden soll.

    **errPronePar** : Array of sp.symbols
        Liste (Array) aller fehlerbehafteten Größen der Gleichung. Elemente des Arrays müssen sp.symbols sein.
        Diese Werte konnen bspw als x_sym, y_sym, z_sym etc. bezeichnet werden und sind ungleich den Werten für x, y, z.
        Für die Werte wird daher die Bezeichnung x_val, y_val, z_val etc. genutzt und für deren Fehler err_x, err_y, err_z etc.

    Return
    ----------
    **absolut_err** : sympy function
        Gibt die Fehlergleichung des absoluten Fehlers wieder. 

    **relativ_err** : sympy function
        Gibt die Fehlergleichung des relativen Fehlers wieder. 
        
    **errProneParamters** : array
        Liste aller Fehlerbehafteten Größen
    """ 

    error = 0
    errProneParamters = []
    for errPar in errPronePar:
        # Für das darstellen der Fehler wird \Delta in Latex später hier eingesetzt
        delta = sp.symbols('thisWillTurnDelta' + errPar.name)

        # Die Funktion wird nach der fehlerbehafteten Variable abgeleitet
        partial = sp.diff(func, errPar) * delta  

        # Fehler werden quadratisch aufsummiert
        error = error + partial**2
        errProneParamters.append((errPar,delta))

    absolut_err=sp.simplify(sp.sqrt(error),rational = True)             
    relativ_err=sp.simplify(sp.sqrt(error/func**2),rational = True)

    return absolut_err, relativ_err, errProneParamters


# --------------------------------------------------
#   Standard Abweichung Sigma
# --------------------------------------------------
def std_abw(p1, p2, err_p1 = 0.0, err_p2 = 0.0):
    """
    Funktion zum berechnen der Sigma-Abweichugn von zwei Messwerten, oder einem Messwert und einem Literaturwert.
    """
    return round(abs(p1 - p2)/(np.sqrt(err_p1**2 + err_p2**2)), 2)

# --------------------------------------------------
#   Berechnung von Wert und Fehler
# --------------------------------------------------
def calc_with_err(func, errFunc, values):
    """
    Methode zum berechnen von Werten und deren Fehler.

    Parameter
    ----------
    **func** : sympy function
         Funktionen mit Parametern. 

    **errFunc** : sympy function
        Die zu func gehörende Fehlerformel nach gff.

    **values** : Array
        Werte, die in die Funktionen eingesetzt werden.
        Als array von Tupeln der Form [(a,da),(b,db),...] oder als array/liste [a,da,b,db,...] 
        (Reihenfolge muss die sein, in der die Argumente in der Funktion genommen werden)

    Return
    ----------

    """
    #Falls der Input in mehrere Tupel aufgeteilt ist, werden diese zu einem Array zusammengefügt 
    if (np.ndim(values) != 1):                    
        values = np.concatenate(values)

    value = func(*values[::2])
    uncertainty = errFunc(*values)

    return value, uncertainty

# --------------------------------------------------
#   Formel zu Latex Code
# --------------------------------------------------
def tex_func(func):
    """
    Zeigt die Formel als gerenderte Math-Darstellung und darunter
    den Latex-Quelltext, der per Button kopiert werden kann. (Für leichtere Benutzung als HTML).
    
    Zudem wird die Latexformel als Variable in Latex gespeichert.

    Parameters
    ----------
    **f** : sympy function
    """
    display(Math(sp.latex(func, long_frac_ratio=2).replace('thisWillTurnDelta', r'\Delta ')))

    latex_str = sp.latex(func, long_frac_ratio=2).replace('thisWillTurnDelta', r'\Delta ')

    html = f"""
    <div style="margin-top:0.5em;">
        <code id="latex-code-{id(func)}" style="display:none;">
            {latex_str}
        </code>

        <button onclick="
            const tex_as_txt = document.getElementById('latex-code-{id(func)}').innerText;
            navigator.clipboard.writeText(tex_as_txt);
        " style="
            padding:4px 8px;
            cursor:pointer;
        ">
            LaTeX kopieren
        </button>
    </div>
    """
    display(HTML(html))

# --------------------------------------------------
#   Papulation (Clac everything)
# --------------------------------------------------
def do_it(function, params, data, params_without_error=None, print_formula: bool = True):
    """
    Calculates values and uncertainties for a given SymPy function.

    Automatically handles both:
    - 1D input (single measurement set)
    - 2D input (multiple measurements)

    Parameters
    ----------
    function : sympy expression
        Function to evaluate

    params : list of sympy symbols
        Parameters of the function

    data : array-like
        Measurement data (1D or 2D). Format:
        [x1, dx1, x2, dx2, ...] or rows of such

    params_without_error : list, optional
        Parameters without uncertainty

    print_formula : bool
        If True, prints LaTeX formulas

    Returns
    -------
    np.ndarray
        Array of [value, uncertainty] for each row
    """

    if params_without_error is None:
        params_without_error = []

    data = np.array(data)

    # Normalize input
    if data.ndim == 1:
        data = data.reshape(1, -1)

    # Expand data to include missing uncertainty columns
    exp_data = np.zeros((data.shape[0], data.shape[1] + len(params_without_error)))

    i = j = z = 0

    while i < len(params):
        if params[i] in params_without_error:
            exp_data[:, j] = data[:, z]
            i += 1
            j += 2
            z += 1
        else:
            exp_data[:, j] = data[:, z]
            exp_data[:, j + 1] = data[:, z + 1]
            i += 1
            j += 2
            z += 2

    # Build parameter list without errors
    params_with_error = [
        p for p in params if p not in params_without_error
    ]

    # Lambdify functions
    f = sp.lambdify(params, function, "numpy")

    absolut_err, relativ_err, parameters = gff(function, params)
    err_abs = sp.lambdify(np.concatenate(parameters), absolut_err, "numpy")

    # Compute results
    results = np.zeros((data.shape[0], 2))

    for n in range(data.shape[0]):
        results[n, :] = calc_with_err(f, err_abs, exp_data[n, :])

    # Remove error contributions for params_without_error
    for p in params_without_error:
        absolut_err = absolut_err.subs('thisWillTurnDelta' + p.name, 0)
        relativ_err = relativ_err.subs('thisWillTurnDelta' + p.name, 0)

    # LaTeX output
    function = sp.simplify(function, rational=True)
    function = sp.separatevars(function)

    if print_formula:
        print("gegebene Funktion:")
        tex_func(function)
        print("-----------------------")

        print("Formel des absoluten Fehlers:")
        tex_func(absolut_err)
        print("-----------------------")

        print("Formel des relativen Fehlers:")
        tex_func(relativ_err)
        print("-----------------------")

    return results






# ==================================================
#   Konstanten (in Heidelberg)
# ==================================================
class const:
    # Ortsgravitationsbeschleunigung in Heidelberg [m/s^2]
    g_hd = 9.80984
    err_g_hd = 0.00002

    # Lichtgeschwindigkeit im Vakuum [m/s]
    c = 299792458
    err_c = 0  # exakt definiert

    # Planck-Konstante [J s]
    h = 6.62607015e-34
    err_h = 0  # exakt definiert

    # reduzierte Planck-Konstante [J s]
    h_bar = 1.054571817e-34
    err_h_bar = 0

    # Elementarladung [C]
    e = 1.602176634e-19
    err_e = 0  # exakt definiert

    # Elektronenmasse [kg]
    m_e = 9.1093837015e-31
    err_m_e = 2.8e-40

    # Protonenmasse [kg]
    m_p = 1.67262192369e-27
    err_m_p = 5.1e-37

    # Neutronenmasse [kg]
    m_n = 1.67492749804e-27
    err_m_n = 9.5e-37

    # Avogadro-Konstante [1/mol]
    N_A = 6.02214076e23
    err_N_A = 0  # exakt definiert

    # Boltzmann-Konstante [J/K]
    k_B = 1.380649e-23
    err_k_B = 0  # exakt definiert

    # Universelle Gaskonstante [J/(mol K)]
    R = 8.314462618
    err_R = 0

    # Gravitationskonstante [m^3/(kg s^2)]
    G = 6.67430e-11
    err_G = 1.5e-15

    # elektrische Feldkonstante (Vakuumpermittivität) [F/m]
    epsilon_0 = 8.8541878128e-12
    err_epsilon_0 = 1.3e-21

    # magnetische Feldkonstante [H/m]
    mu_0 = 1.25663706212e-6
    err_mu_0 = 1.9e-16

    # Feinstrukturkonstante [-]
    alpha = 7.2973525693e-3
    err_alpha = 1.1e-12

    # Rydberg-Konstante [1/m]
    R_inf = 10973731.568160
    err_R_inf = 2.1e-5

# ==================================================
#   Sympy default symbols
# ==================================================
class Sympy_Symbols:
    '''
    Class, welche die Standart symbolde für Sympy functions beinhält. 
    '''
    # Latin lowercase
    for _c in string.ascii_lowercase:
        locals()[_c] = sp.symbols(rf"{_c}")

    # Latin uppercase
    for _c in string.ascii_uppercase:
        locals()[_c] = sp.symbols(rf"{_c}")

    # Greek lowercase
    _greek_lower = [
        "alpha","beta","gamma","delta","epsilon","zeta","eta","theta",
        "iota","kappa","lambda","mu","nu","xi","omicron","pi","rho",
        "sigma","tau","upsilon","phi","chi","psi","omega"
    ]

    for _g in _greek_lower:
        locals()[_g] = sp.symbols(rf"\{_g}")

    # Greek uppercase (only those commonly distinct in math)
    _greek_upper = [
        "Gamma","Delta","Theta","Lambda","Xi","Pi","Sigma","Upsilon","Phi","Psi","Omega"
    ]

    for _g in _greek_upper:
        locals()[_g] = sp.symbols(rf"\{_g}")

    # cleanup helper variables
    del _c, _g, _greek_lower, _greek_upper

# ==================================================
#   Default colors
# ==================================================
class Colors(str, Enum):
    LIGHT_GREEN = "#1ba97f"

    ORANGE = "#D56A11"
    BLUE = "#073984"
    PINK = "#B42B6F"
    MINT = "#67DAB8"
    TEAL = "#106E7D"
    YELLOW = "#FFC107"
    DARK_PURPLE = "#4B2D83"
    PEACH = "#FEA75F"
    LIGHT_BLUE = "#3485FD"
    PURPLE = "#6610f2"

    RED = "#b92d3a"
    WINE_RED = "#721c24"
    GREEN = "#1b722f"


    def rgb(self):
        h = self.value.lstrip("#")
        return tuple(int(h[i:i+2], 16) for i in (0, 2, 4))
    
# ==================================================
#   Default Plotting system
# ==================================================
# Default Sytling
def plot_me(title:str=None, xlabel:str=None, ylabel:str=None, safe_as:str=None):
    from cycler import cycler
    import matplotlib.pyplot as plt

    # Apply style
    plt.rcParams["axes.prop_cycle"] = cycler(
        color=[c.value for c in Colors]
    )
    plt.rcParams["grid.linestyle"] = "--"
    plt.rcParams["grid.alpha"] = 0.6
    plt.rcParams["axes.grid"] = True

    # Apply plot-specific settings
    ax = plt.gca()
    ax.ticklabel_format(axis='both', style='sci', scilimits=(0,3))


    if title:
        ax.set_title(title)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)

    # Always show legend (only works if labeled plots exist)
    ax.legend()

    if safe_as != None:
        plt.savefig(f'{safe_as}.pdf', dpi=300)

    plt.show()

    return ax

# Vertical Errorbarlines
def plot_errorlines(
    ax,
    x_val,
    err,
    y_min,
    y_max,
    names:str = None,
    color_line=Colors.RED,
    color_err=Colors.RED,
    add_label:bool=False,
    label:str='Messreihe'
):
    import matplotlib.patches as patches

    # Error rectangles
    for x in x_val:
        rect = patches.Rectangle(
            (x - err, y_min),
            2 * err,
            y_max - y_min,
            color=color_err,
            alpha=0.5
        )
        ax.add_patch(rect)

    # Vertical lines
    for i, x in enumerate(x_val):
        ax.axvline(
            x,
            color=color_line,
            linestyle=':',
            alpha=0.9,
            label=label if (add_label and i == 0) else None
        )

    # Labels
    height_scale = [0.95, 0.85, 0.75, 0.65]


    if names != None:
        for i, (x, name) in enumerate(zip(x_val, names)):
            ax.text(
                x,
                y_max * height_scale[i % len(height_scale)],
                name,
                color='black',
                ha='center',
                va='top',
                fontsize=10,
                bbox=dict(
                    boxstyle='round',
                    facecolor='white',
                    edgecolor=color_line,
                    alpha=0.9,
                    pad=0.45
                )
            )

    return ax