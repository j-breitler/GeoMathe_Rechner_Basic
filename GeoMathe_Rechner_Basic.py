#-------------------- GeoMathe_Rechner_Advanced.py -----------#
#-------------------- Author: Ai & Julian Breitler -------------#
#------- Co-Author & Helferlein: Fenja Runfors, Klara Gössler, Christopher Stering---#
#------------------ Viel Spaß beim Verwenden! ----------------#

#-------------------------------------------------------------#
#------------------Benützung auf eigene Gefahr !!! -----------#
#-------------------------------------------------------------#



import math
import tkinter as tk
from tkinter import messagebox, ttk

# -----------------------------------------------------------
# --- GLOBAL CONSTANTS & PURE CALCULATION FUNCTIONS (HA_combined) ---
# -----------------------------------------------------------

# Define the Gon constant: 400 gon = 2 * pi radians
GON_PER_RADIAN = 200 / math.pi
RADIAN_PER_GON = math.pi / 200

def normalize_azimuth(azimuth_gon):
    """
    Normalizes an azimuth value to the standard geodetic range [0, 400) gon 
    and returns the normalized azimuth and a message if normalization occurred.
    """
    original_az = azimuth_gon
    message = ""
    
    normalized_az = azimuth_gon % 400
    
    # Check if normalization actually occurred for the message
    if normalized_az != original_az:
        # Simplified message generation based on original logic
        if original_az < 0:
            message = f"Modulo operation: {original_az:.4f} + {math.ceil(-original_az / 400) * 400} gon"
        elif original_az >= 400:
            message = f"Modulo operation: {original_az:.4f} - {math.floor(original_az / 400) * 400} gon"
        
    if normalized_az == 400:
        normalized_az = 0.0
        
    return normalized_az, message

def calculate_geodetic_inverse(x1, y1, x2, y2):
    """
    Calculates the distance (S12) and forward/backward azimuths (Az12, Az21).
    """
    dx = x2 - x1
    dy = y2 - y1
    
    if dx == 0 and dy == 0:
        return 0.0, 0.0, 0.0, "", ""
        
    S12 = math.sqrt(dx**2 + dy**2)
    
    Az12_rad = math.atan2(dy, dx)
    Az12_gon_unnormalized = Az12_rad * GON_PER_RADIAN
    
    Az12_gon, Az12_msg = normalize_azimuth(Az12_gon_unnormalized)
        
    Az21_unnormalized = Az12_gon + 200
    Az21_gon, Az21_msg = normalize_azimuth(Az21_unnormalized)

    return S12, Az12_gon, Az21_gon, Az12_msg, Az21_msg

def calculate_geodetic_forward(x1, y1, S12, Az12_gon):
    """
    Calculates the coordinates of P2 (x2, y2) given P1 (x1, y1), 
    distance S12, and forward azimuth Az12.
    """
    Az12_gon_norm, Az12_msg = normalize_azimuth(Az12_gon)
    Az12_rad = Az12_gon_norm * RADIAN_PER_GON
    
    dx = S12 * math.cos(Az12_rad)
    dy = S12 * math.sin(Az12_rad)
    
    x2 = x1 + dx
    y2 = y1 + dy
    
    return x2, y2, Az12_msg

# -----------------------------------------------------------
# --- PURE CALCULATION FUNCTIONS (HWS_berechnung_alle_Winkel) ---
# -----------------------------------------------------------

def calculate_all_angles(a: float, b: float, c: float) -> dict:
    """
    Calculates all three angles (alpha, beta, gamma) of a triangle 
    using the Half-Angle Theorem (Halbwinkelsatz), returning the results in Gon.
    """
    
    # --- Input Validation (Triangle Inequality) ---
    if a + b <= c or a + c <= b or b + c <= a:
        raise ValueError("Ungültiges Dreieck (Dreiecksungleichung verletzt)")

    # --- Calculation of Semi-perimeter and differences ---
    s = (a + b + c) / 2
    s_minus_a = s - a
    s_minus_b = s - b
    s_minus_c = s - c

    # --- General Calculation Helper ---
    def calculate_angle(s_opp_num, s_adj1_num, s_adj2_num):
        """Calculates an angle based on the side differences."""
        zaehler = s_adj1_num * s_adj2_num
        nenner = s * s_opp_num
        
        # Check for non-physical/edge cases before sqrt/atan
        if nenner == 0 or zaehler < 0 or nenner < 0:
             return math.nan, math.nan 

        tan_half = math.sqrt(zaehler / nenner)
        angle_rad = 2 * math.atan(tan_half)
        angle_gon = angle_rad * GON_PER_RADIAN
        return angle_rad, angle_gon
    
    # --- Calculate Angles ---
    alpha_rad, alpha_gon = calculate_angle(s_minus_a, s_minus_b, s_minus_c)
    beta_rad, beta_gon = calculate_angle(s_minus_b, s_minus_a, s_minus_c)
    gamma_rad, gamma_gon = calculate_angle(s_minus_c, s_minus_a, s_minus_b)

    # --- Return Results ---
    return {
        "Semiperimeter (s)": s,
        "s - a": s_minus_a,
        "s - b": s_minus_b,
        "s - c": s_minus_c,
        
        "--- Ergebnisse im Gon-System ---": "", # Separator
        
        "Alpha (α) [Gon]": alpha_gon,
        "Beta (β) [Gon]": beta_gon,
        "Gamma (γ) [Gon]": gamma_gon,
        
        "--- Überprüfung ---": "", # Separator
        
        "Summe (α+β+γ) [Gon]": alpha_gon + beta_gon + gamma_gon
    }

# -----------------------------------------------------------
# --- PURE CALCULATION FUNCTIONS (Numerisch_stabiler_Algorithmus) ---
# -----------------------------------------------------------

def compute_N_from_LMR_gon(xL, yL, xM, yM, xR, yR, alpha_gon, beta_gon):
    """
    Calculates the coordinates of point N using the numerically stable 
    resection algorithm (Pothenot-Snellius).
    """
    
    # 1) Distances s_ML, s_MR from coordinates
    dx_ML = xL - xM
    dy_ML = yL - yM
    dx_MR = xR - xM
    dy_MR = yR - yM

    s_ML = math.hypot(dx_ML, dy_ML)
    s_MR = math.hypot(dx_MR, dy_MR)

    # 2) Convert angles to radians
    alpha = alpha_gon * RADIAN_PER_GON
    beta = beta_gon * RADIAN_PER_GON

    # 3) Azimuths ν_ML, ν_MR (in radians)
    nu_ML = math.atan2(dy_ML, dx_ML)
    nu_MR = math.atan2(dy_MR, dx_MR)

    # 4) a, b
    a = math.sin(alpha) / s_ML
    b = -math.sin(beta) / s_MR

    # 5) Denominator
    den = math.sin(nu_ML - nu_MR + alpha + beta)
    if abs(den) < 1e-12:
        raise ValueError("Denominator ~ 0 (Geometrieproblem / Singularität).")

    # 6) λ and µ
    lam = (a * math.cos(nu_MR - beta) - b * math.cos(nu_ML + alpha)) / den
    mu = (a * math.sin(nu_MR - beta) - b * math.sin(nu_ML + alpha)) / den

    # 7) s_mn^2 (Square of the distance M to N)
    s_mn_sq = 1.0 / (lam**2 + mu**2)

    # 8) Δx_mn, Δy_mn (Coordinate differences M to N)
    dx_mn = lam * s_mn_sq
    dy_mn = mu * s_mn_sq

    # 9) Point N
    x_n = xM + dx_mn
    y_n = yM + dy_mn

    return {
        "s_ML": s_ML,
        "s_MR": s_MR,
        "alpha_gon": alpha_gon,
        "beta_gon": beta_gon,
        "nu_ML_gon": nu_ML * GON_PER_RADIAN,
        "nu_MR_gon": nu_MR * GON_PER_RADIAN,
        "a": a,
        "b": b,
        "Nenner von Formel": den,
        "lambda": lam,
        "mu": mu,
        "s_mn^2": s_mn_sq,
        "dx_mn": dx_mn,
        "dy_mn": dy_mn,
        "x_n": x_n,
        "y_n": y_n,
    }

# -----------------------------------------------------------
# --- MAIN APPLICATION CLASS ---
# -----------------------------------------------------------

class GeodesyApp:
    def __init__(self, root):
        self.root = root
        self.root.title("GeoMathe Rechner")
        
        # Start with an EMPTY coordinate dictionary ---
        self.COORDINATE_DATA = {}
        self.POINT_IDS = sorted(self.COORDINATE_DATA.keys())
        # ---------------------------------------------------------
        
        # Dictionary to hold entry/combobox references for easy access
        self.widgets = {}

        # Create a tabbed interface (Notebook)
        self.notebook = ttk.Notebook(root)
        self.notebook.pack(pady=10, padx=10, expand=True, fill="both")

        # Setup all tabs
        self._setup_point_management()
        self._setup_inverse_tab()
        self._setup_forward_tab()
        self._setup_half_angle_tab()
        self._setup_n_point_tab()
        
        # Initial GUI synchronization
        self._load_points_to_treeview()
        self._refresh_comboboxes()
        self.notebook.bind("<<NotebookTabChanged>>", self._on_tab_change)

    # --- Utility Methods ---

    def _on_tab_change(self, event):
        """Called when the user switches tabs to refresh related data."""
        current_tab = self.notebook.tab(self.notebook.select(), "text")
        
        # Refresh coordinate labels if switching to a calculation tab
        if current_tab == '1. Geodaetische Hauptaufgabe':
            self._update_coords_fwd(None)
        elif current_tab == '2. Geodaetische Hauptaufgabe':
            self._update_coords_inv(None)

    def _refresh_comboboxes(self):
        """Updates the list of available points in all relevant Comboboxes."""
        self.POINT_IDS = sorted(self.COORDINATE_DATA.keys())
        
        for key in ['combo_p1_inv', 'combo_p2_inv', 'combo_p1_fwd']:
            if key in self.widgets:
                combo = self.widgets[key]
                combo['values'] = self.POINT_IDS
                # Maintain selection if possible, otherwise set to empty string
                if combo.get() not in self.POINT_IDS:
                    if self.POINT_IDS:
                        combo.set(self.POINT_IDS[0]) # Default to first point if list is not empty
                    else:
                        combo.set("") # Set to empty if list is empty

    # -----------------------------------------------------------
    # --- Tab 1: Point Management Setup and Logic ---
    # -----------------------------------------------------------

    def _setup_point_management(self):
        tab_manage = ttk.Frame(self.notebook)
        self.notebook.add(tab_manage, text='Point Management')

        # Input Frame for adding/editing points
        input_frame_manage = tk.LabelFrame(tab_manage, text="Add/Update Point Coordinates", padx=10, pady=10)
        input_frame_manage.pack(padx=10, pady=10, fill="x")

        tk.Label(input_frame_manage, text="Point ID:").grid(row=0, column=0, padx=5, pady=5, sticky="w")
        self.widgets['entry_p_id'] = tk.Entry(input_frame_manage, width=15)
        self.widgets['entry_p_id'].grid(row=0, column=1, padx=5, pady=5)

        tk.Label(input_frame_manage, text="X Coordinate:").grid(row=1, column=0, padx=5, pady=5, sticky="w")
        self.widgets['entry_x_val'] = tk.Entry(input_frame_manage, width=15)
        self.widgets['entry_x_val'].grid(row=1, column=1, padx=5, pady=5)

        tk.Label(input_frame_manage, text="Y Coordinate:").grid(row=1, column=2, padx=5, pady=5, sticky="w")
        self.widgets['entry_y_val'] = tk.Entry(input_frame_manage, width=15)
        self.widgets['entry_y_val'].grid(row=1, column=3, padx=5, pady=5)

        btn_add_update = tk.Button(input_frame_manage, text="Add/Update Point", command=self._add_update_point, bg="#008CBA", fg="white")
        btn_add_update.grid(row=0, column=2, padx=5, pady=5, sticky="ew")

        btn_delete = tk.Button(input_frame_manage, text="Delete Selected Point", command=self._delete_selected_point, bg="#f44336", fg="white")
        btn_delete.grid(row=0, column=3, padx=5, pady=5, sticky="ew")

        # Display Frame for points (Treeview)
        display_frame_manage = tk.LabelFrame(tab_manage, text="Current Points (Select to Edit)", padx=10, pady=10)
        display_frame_manage.pack(padx=10, pady=10, fill="both", expand=True)

        self.widgets['tree_points'] = ttk.Treeview(display_frame_manage, columns=("ID", "X", "Y"), show="headings")
        self.widgets['tree_points'].heading("ID", text="Point ID")
        self.widgets['tree_points'].heading("X", text="X (Northing)")
        self.widgets['tree_points'].heading("Y", text="Y (Easting)")
        self.widgets['tree_points'].column("ID", width=100, anchor="center")
        self.widgets['tree_points'].column("X", width=150, anchor="e")
        self.widgets['tree_points'].column("Y", width=150, anchor="e")
        self.widgets['tree_points'].pack(fill="both", expand=True)
        self.widgets['tree_points'].bind("<<TreeviewSelect>>", self._select_point_from_tree)

    def _load_points_to_treeview(self):
        """Clears and reloads all points from COORDINATE_DATA into the Treeview."""
        tree_points = self.widgets['tree_points']
        for i in tree_points.get_children():
            tree_points.delete(i)
            
        for name, (x, y) in self.COORDINATE_DATA.items():
            tree_points.insert("", "end", iid=name, values=(name, f"{x:.4f}", f"{y:.4f}"))

    def _add_update_point(self):
        """Reads input fields and adds/updates a point in COORDINATE_DATA."""
        try:
            p_id = self.widgets['entry_p_id'].get().strip()
            x_val = float(self.widgets['entry_x_val'].get())
            y_val = float(self.widgets['entry_y_val'].get())
            
            if not p_id:
                messagebox.showerror("Input Error", "Point ID cannot be empty.")
                return

            self.COORDINATE_DATA[p_id] = (x_val, y_val)
            
            self._load_points_to_treeview()
            self._refresh_comboboxes()
            
        except ValueError:
            messagebox.showerror("Input Error", "Please enter valid numbers for X and Y coordinates.")
        except Exception as e:
            messagebox.showerror("Error", f"An unexpected error occurred: {e}")

    def _select_point_from_tree(self, event):
        """Loads the selected point from the Treeview back into the input fields for editing."""
        tree_points = self.widgets['tree_points']
        selected_item = tree_points.focus()
        if selected_item:
            values = tree_points.item(selected_item, 'values')
            if values:
                self.widgets['entry_p_id'].delete(0, tk.END)
                self.widgets['entry_x_val'].delete(0, tk.END)
                self.widgets['entry_y_val'].delete(0, tk.END)
                
                self.widgets['entry_p_id'].insert(0, values[0])
                self.widgets['entry_x_val'].insert(0, values[1])
                self.widgets['entry_y_val'].insert(0, values[2])
            
    def _delete_selected_point(self):
        """Deletes the selected point from the Treeview and the data structure."""
        tree_points = self.widgets['tree_points']
        selected_item = tree_points.focus()
        if selected_item:
            p_id = tree_points.item(selected_item, 'values')[0]
            
            if messagebox.askyesno("Confirm Delete", f"Are you sure you want to delete point '{p_id}'?"):
                del self.COORDINATE_DATA[p_id]
                
                self._load_points_to_treeview()
                self._refresh_comboboxes()
                
                self.widgets['entry_p_id'].delete(0, tk.END)
                self.widgets['entry_x_val'].delete(0, tk.END)
                self.widgets['entry_y_val'].delete(0, tk.END)
        else:
            messagebox.showwarning("Selection Error", "Please select a point to delete.")

    # -----------------------------------------------------------
    # --- Tab 2: Inverse Geodetic Problem (2.HA) Setup and Logic ---
    # -----------------------------------------------------------

    def _setup_inverse_tab(self):
        tab_inverse = ttk.Frame(self.notebook)
        self.notebook.add(tab_inverse, text='2.Geodätische Hauptaufgabe') 

        input_frame_inv = tk.LabelFrame(tab_inverse, text="Input Point IDs (P1 | P2)", padx=10, pady=10)
        input_frame_inv.pack(padx=10, pady=10, fill="x")

        # Point 1 (P1 ID)
        tk.Label(input_frame_inv, text="Point ID P₁:").grid(row=0, column=0, padx=5, pady=5, sticky="w")
        self.widgets['combo_p1_inv'] = ttk.Combobox(input_frame_inv, values=self.POINT_IDS, state="readonly", width=12)
        self.widgets['combo_p1_inv'].grid(row=0, column=1, padx=5, pady=5)
        self.widgets['combo_p1_inv'].set("") 
        self.widgets['coord_x1_inv'] = tk.Label(input_frame_inv, text=f"X₁: --.--", anchor="w", fg="gray")
        self.widgets['coord_x1_inv'].grid(row=1, column=0, padx=5, pady=0, sticky="w")
        self.widgets['coord_y1_inv'] = tk.Label(input_frame_inv, text=f"Y₁: --.--", anchor="w", fg="gray")
        self.widgets['coord_y1_inv'].grid(row=1, column=1, padx=5, pady=0, sticky="w")

        # Point 2 (P2 ID)
        tk.Label(input_frame_inv, text="Point ID P₂:").grid(row=0, column=2, padx=5, pady=5, sticky="w")
        self.widgets['combo_p2_inv'] = ttk.Combobox(input_frame_inv, values=self.POINT_IDS, state="readonly", width=12)
        self.widgets['combo_p2_inv'].grid(row=0, column=3, padx=5, pady=5)
        self.widgets['combo_p2_inv'].set("")
        self.widgets['coord_x2_inv'] = tk.Label(input_frame_inv, text=f"X₂: --.--", anchor="w", fg="gray")
        self.widgets['coord_x2_inv'].grid(row=1, column=2, padx=5, pady=0, sticky="w")
        self.widgets['coord_y2_inv'] = tk.Label(input_frame_inv, text=f"Y₂: --.--", anchor="w", fg="gray")
        self.widgets['coord_y2_inv'].grid(row=1, column=3, padx=5, pady=0, sticky="w")

        self.widgets['combo_p1_inv'].bind("<<ComboboxSelected>>", self._update_coords_inv)
        self.widgets['combo_p2_inv'].bind("<<ComboboxSelected>>", self._update_coords_inv)

        # Calculation Button
        calc_button_inv = tk.Button(tab_inverse, text="Berechne 2.HA (S, Nu)", command=self._run_inverse_calculation, bg="#4CAF50", fg="white", font=('Arial', 10, 'bold'))
        calc_button_inv.pack(padx=10, pady=10, fill="x")

        # Result Frame
        result_frame_inv = tk.LabelFrame(tab_inverse, text="Results", padx=10, pady=10)
        result_frame_inv.pack(padx=10, pady=10, fill="x")

        tk.Label(result_frame_inv, text="Distance Sₙₘ:").grid(row=0, column=0, padx=5, pady=5, sticky="w")
        self.widgets['result_S12'] = tk.Label(result_frame_inv, text="--.--", width=15, anchor="w", fg="#004D40", font=('Arial', 10, 'bold'))
        self.widgets['result_S12'].grid(row=0, column=1, padx=5, pady=5, sticky="w")

        tk.Label(result_frame_inv, text="Orientierte Richtung νₙₘ:").grid(row=1, column=0, padx=5, pady=5, sticky="w")
        self.widgets['result_Az12'] = tk.Label(result_frame_inv, text="--.-- gon", width=15, anchor="w", fg="#004D40", font=('Arial', 10, 'bold'))
        self.widgets['result_Az12'].grid(row=1, column=1, padx=5, pady=5, sticky="w")
        self.widgets['msg_Az12'] = tk.Label(result_frame_inv, text="", fg="red")
        self.widgets['msg_Az12'].grid(row=1, column=2, padx=5, pady=5, sticky="w")

        tk.Label(result_frame_inv, text="Orientierte Richtung νₘₙ:").grid(row=2, column=0, padx=5, pady=5, sticky="w")
        self.widgets['result_Az21'] = tk.Label(result_frame_inv, text="--.-- gon", width=15, anchor="w", fg="#004D40", font=('Arial', 10, 'bold'))
        self.widgets['result_Az21'].grid(row=2, column=1, padx=5, pady=5, sticky="w")
        self.widgets['msg_Az21'] = tk.Label(result_frame_inv, text="", fg="red")
        self.widgets['msg_Az21'].grid(row=2, column=2, padx=5, pady=5, sticky="w")

    def _update_coords_inv(self, event):
        """Updates the displayed coordinates based on selected P1 and P2 IDs."""
        p1_id = self.widgets['combo_p1_inv'].get()
        p2_id = self.widgets['combo_p2_inv'].get()
        
        if p1_id in self.COORDINATE_DATA:
            x1, y1 = self.COORDINATE_DATA[p1_id]
            self.widgets['coord_x1_inv'].config(text=f"X₁: {x1:.4f}")
            self.widgets['coord_y1_inv'].config(text=f"Y₁: {y1:.4f}")
        else:
            self.widgets['coord_x1_inv'].config(text="X₁: --.--")
            self.widgets['coord_y1_inv'].config(text="Y₁: --.--")

        if p2_id in self.COORDINATE_DATA:
            x2, y2 = self.COORDINATE_DATA[p2_id]
            self.widgets['coord_x2_inv'].config(text=f"X₂: {x2:.4f}")
            self.widgets['coord_y2_inv'].config(text=f"Y₂: {y2:.4f}")
        else:
            self.widgets['coord_x2_inv'].config(text="X₂: --.--")
            self.widgets['coord_y2_inv'].config(text="Y₂: --.--")

    def _run_inverse_calculation(self):
        """Runs the Inverse Problem calculation."""
        try:
            p1_id = self.widgets['combo_p1_inv'].get()
            p2_id = self.widgets['combo_p2_inv'].get()
            
            if p1_id not in self.COORDINATE_DATA or p2_id not in self.COORDINATE_DATA:
                messagebox.showerror("Input Error", "Please select valid point IDs that exist in Point Management.")
                return
                
            x1, y1 = self.COORDINATE_DATA[p1_id]
            x2, y2 = self.COORDINATE_DATA[p2_id]
            
            S12, Az12, Az21, Az12_msg, Az21_msg = calculate_geodetic_inverse(x1, y1, x2, y2)
            
            self.widgets['result_S12'].config(text=f"{S12:.4f}")
            self.widgets['result_Az12'].config(text=f"{Az12:.4f} gon")
            self.widgets['result_Az21'].config(text=f"{Az21:.4f} gon")

            self.widgets['msg_Az12'].config(text=f"ⓘ {Az12_msg}" if Az12_msg else "", fg="red")
            self.widgets['msg_Az21'].config(text=f"ⓘ {Az21_msg}" if Az21_msg else "", fg="red")
            
        except Exception as e:
            messagebox.showerror("Calculation Error", f"An unexpected error occurred: {e}")

    # -----------------------------------------------------------
    # --- Tab 3: Forward Geodetic Problem (1.HA) Setup and Logic ---
    # -----------------------------------------------------------
    
    def _setup_forward_tab(self):
        tab_forward = ttk.Frame(self.notebook)
        self.notebook.add(tab_forward, text='1. Geodätische Hauptaufgabe')

        input_frame_fwd = tk.LabelFrame(tab_forward, text="Input (P1 ID | Sₙₘ | νₙₘ)", padx=10, pady=10)
        input_frame_fwd.pack(padx=10, pady=10, fill="x")

        # Point 1 (P1 ID)
        tk.Label(input_frame_fwd, text="Point ID P₁:").grid(row=0, column=0, padx=5, pady=5, sticky="w")
        self.widgets['combo_p1_fwd'] = ttk.Combobox(input_frame_fwd, values=self.POINT_IDS, state="readonly", width=12)
        self.widgets['combo_p1_fwd'].grid(row=0, column=1, padx=5, pady=5)
        self.widgets['combo_p1_fwd'].set("") 
        self.widgets['coord_x1_fwd'] = tk.Label(input_frame_fwd, text=f"X₁: --.--", anchor="w", fg="gray")
        self.widgets['coord_x1_fwd'].grid(row=1, column=0, padx=5, pady=0, sticky="w")
        self.widgets['coord_y1_fwd'] = tk.Label(input_frame_fwd, text=f"Y₁: --.--", anchor="w", fg="gray")
        self.widgets['coord_y1_fwd'].grid(row=1, column=1, padx=5, pady=0, sticky="w")

        self.widgets['combo_p1_fwd'].bind("<<ComboboxSelected>>", self._update_coords_fwd)

        # Distance (S12) and Azimuth (Az12)
        tk.Label(input_frame_fwd, text="Distance Sₙₘ:").grid(row=0, column=2, padx=5, pady=5, sticky="w")
        self.widgets['entry_S12_fwd'] = tk.Entry(input_frame_fwd, width=15)
        self.widgets['entry_S12_fwd'].grid(row=0, column=3, padx=5, pady=5)
        self.widgets['entry_S12_fwd'].insert(0, "0") 

        tk.Label(input_frame_fwd, text="Azimuth νₙₘ:").grid(row=1, column=2, padx=5, pady=5, sticky="w")
        self.widgets['entry_Az12_fwd'] = tk.Entry(input_frame_fwd, width=15)
        self.widgets['entry_Az12_fwd'].grid(row=1, column=3, padx=5, pady=5)
        self.widgets['entry_Az12_fwd'].insert(0, "0") 

        # Calculation Button
        calc_button_fwd = tk.Button(tab_forward, text="Berechne 1.HA (Xₙ, Yₙ)", command=self._run_forward_calculation, bg="#008CBA", fg="white", font=('Arial', 10, 'bold'))
        calc_button_fwd.pack(padx=10, pady=10, fill="x")

        # Result Frame
        result_frame_fwd = tk.LabelFrame(tab_forward, text="Results (P: Xₙ, Yₙ)", padx=10, pady=10)
        result_frame_fwd.pack(padx=10, pady=10, fill="x")

        tk.Label(result_frame_fwd, text="Xₙ (Northing):").grid(row=0, column=0, padx=5, pady=5, sticky="w")
        self.widgets['result_x2'] = tk.Label(result_frame_fwd, text="--.--", width=15, anchor="w", fg="#004D40", font=('Arial', 10, 'bold'))
        self.widgets['result_x2'].grid(row=0, column=1, padx=5, pady=5, sticky="w")

        tk.Label(result_frame_fwd, text="Yₙ (Easting):").grid(row=1, column=0, padx=5, pady=5, sticky="w")
        self.widgets['result_y2'] = tk.Label(result_frame_fwd, text="--.--", width=15, anchor="w", fg="#004D40", font=('Arial', 10, 'bold'))
        self.widgets['result_y2'].grid(row=1, column=1, padx=5, pady=5, sticky="w")

        self.widgets['msg_Az12_fwd'] = tk.Label(result_frame_fwd, text="", fg="red")
        self.widgets['msg_Az12_fwd'].grid(row=2, column=0, columnspan=2, padx=5, pady=5, sticky="w")

    def _update_coords_fwd(self, event):
        """Updates the displayed P1 coordinates based on selected P1 ID."""
        p1_id = self.widgets['combo_p1_fwd'].get()
        
        if p1_id in self.COORDINATE_DATA:
            x1, y1 = self.COORDINATE_DATA[p1_id]
            self.widgets['coord_x1_fwd'].config(text=f"X₁: {x1:.4f}")
            self.widgets['coord_y1_fwd'].config(text=f"Y₁: {y1:.4f}")
        else:
            self.widgets['coord_x1_fwd'].config(text="X₁: --.--")
            self.widgets['coord_y1_fwd'].config(text="Y₁: --.--")

    def _run_forward_calculation(self):
        """Runs the Forward Problem calculation."""
        try:
            p1_id = self.widgets['combo_p1_fwd'].get()
            
            S12 = float(self.widgets['entry_S12_fwd'].get())
            Az12_gon = float(self.widgets['entry_Az12_fwd'].get())
            
            if p1_id not in self.COORDINATE_DATA:
                messagebox.showerror("Input Error", "Please select a valid point ID for P1 that exists in Point Management.")
                return

            x1, y1 = self.COORDINATE_DATA[p1_id]
            
            if S12 < 0:
                raise ValueError("Distance (S12) must be non-negative.")
                
            x2, y2, Az12_msg = calculate_geodetic_forward(x1, y1, S12, Az12_gon)
            
            self.widgets['result_x2'].config(text=f"{x2:.4f}")
            self.widgets['result_y2'].config(text=f"{y2:.4f}")

            self.widgets['msg_Az12_fwd'].config(text=f"ⓘ Normalized Az: {Az12_msg}" if Az12_msg else "", fg="red")
            
        except ValueError as ve:
            messagebox.showerror("Input Error", f"Please enter valid numbers/selections: {ve}")
        except Exception as e:
            messagebox.showerror("Calculation Error", f"An unexpected error occurred: {e}")


    # -----------------------------------------------------------
    # --- Tab 4: Half-Angle Theorem (HWS) Setup and Logic ---
    # -----------------------------------------------------------

    def _setup_half_angle_tab(self):
        tab_hws = ttk.Frame(self.notebook)
        self.notebook.add(tab_hws, text='Halbwinkelsatz (HWS)')

        frame = ttk.Frame(tab_hws, padding=10)
        frame.pack(fill="both", expand=True)

        ttk.Label(frame, text="For comma, please use dot '.' ").grid(
            row=0, column=0, columnspan=2, sticky="w", pady=(0, 10)
        )

        # Input fields
        ttk.Label(frame, text="Seite a:").grid(row=1, column=0, sticky="w")
        self.widgets['entry_a_hws'] = ttk.Entry(frame)
        self.widgets['entry_a_hws'].grid(row=1, column=1)
        self.widgets['entry_a_hws'].insert(0, "0")

        ttk.Label(frame, text="Seite b:").grid(row=2, column=0, sticky="w")
        self.widgets['entry_b_hws'] = ttk.Entry(frame)
        self.widgets['entry_b_hws'].grid(row=2, column=1)
        self.widgets['entry_b_hws'].insert(0, "0")

        ttk.Label(frame, text="Seite c:").grid(row=3, column=0, sticky="w")
        self.widgets['entry_c_hws'] = ttk.Entry(frame)
        self.widgets['entry_c_hws'].grid(row=3, column=1)
        self.widgets['entry_c_hws'].insert(0, "0")

        # Button
        calc_btn = ttk.Button(frame, text="Berechnen", command=self._run_half_angle_calculation)
        calc_btn.grid(row=4, column=0, columnspan=2, pady=10)

        # Output box
        self.widgets['output_hws'] = tk.Text(frame, height=20, width=50)
        self.widgets['output_hws'].grid(row=5, column=0, columnspan=2)

    def _run_half_angle_calculation(self):
        try:
            a = float(self.widgets['entry_a_hws'].get())
            b = float(self.widgets['entry_b_hws'].get())
            c = float(self.widgets['entry_c_hws'].get())
            
            if a <= 0 or b <= 0 or c <= 0:
                raise ValueError("Side lengths must be positive.")

            result = calculate_all_angles(a, b, c)

            self.widgets['output_hws'].delete(1.0, tk.END)
            for key, value in result.items():
                if key.startswith("---"):
                    self.widgets['output_hws'].insert(tk.END, f"\n{key}\n")
                else:
                    self.widgets['output_hws'].insert(tk.END, f"{key}: {value:.4f}\n")

        except ValueError as e:
            messagebox.showerror("Fehler", str(e))
        except Exception:
            messagebox.showerror("Fehler", "Ungültige Eingabe oder unerwarteter Fehler!")

    # -----------------------------------------------------------
    # --- Tab 5: N-Point Calculation Setup and Logic ---
    # -----------------------------------------------------------

    def _setup_n_point_tab(self):
        tab_n_point = ttk.Frame(self.notebook)
        self.notebook.add(tab_n_point, text='4. Numerissch-stabiler Algorithmus (Rückwärtsschnitt)')

        main_frame = ttk.Frame(tab_n_point, padding="10")
        main_frame.pack(fill="both", expand=True)
        
        row = 0
        
        # Helper for creating input fields
        def create_input_field(parent, key, label, r, default_value):
            ttk.Label(parent, text=label).grid(row=r, column=0, sticky=tk.W, pady=2)
            entry = ttk.Entry(parent, width=20)
            entry.grid(row=r, column=1, sticky=(tk.W, tk.E), pady=2, padx=(10,0))
            entry.insert(0, default_value)
            self.widgets[key] = entry
            return r + 1

        # Punkt L
        ttk.Label(main_frame, text="Punkt L:", font=('Arial', 10, 'bold')).grid(row=row, column=0, columnspan=2, sticky=tk.W, pady=(10,5))
        row += 1
        row = create_input_field(main_frame, "xL", "x-Koordinate L:", row, "0")
        row = create_input_field(main_frame, "yL", "y-Koordinate L:", row, "0")
        
        # Punkt M
        ttk.Label(main_frame, text="Punkt M:", font=('Arial', 10, 'bold')).grid(row=row, column=0, columnspan=2, sticky=tk.W, pady=(10,5))
        row += 1
        row = create_input_field(main_frame, "xM", "x-Koordinate M:", row, "0")
        row = create_input_field(main_frame, "yM", "y-Koordinate M:", row, "0")
        
        # Punkt R
        ttk.Label(main_frame, text="Punkt R:", font=('Arial', 10, 'bold')).grid(row=row, column=0, columnspan=2, sticky=tk.W, pady=(10,5))
        row += 1
        row = create_input_field(main_frame, "xR", "x-Koordinate R:", row, "0")
        row = create_input_field(main_frame, "yR", "y-Koordinate R:", row, "0")
        
        # Winkel
        ttk.Label(main_frame, text="Winkel:", font=('Arial', 10, 'bold')).grid(row=row, column=0, columnspan=2, sticky=tk.W, pady=(10,5))
        row += 1
        row = create_input_field(main_frame, "alpha_gon", "Alpha (Gon):", row, "0")
        row = create_input_field(main_frame, "beta_gon", "Beta (Gon):", row, "0")
        
        # Calculate button
        calc_button = ttk.Button(main_frame, text="Berechne Punkt N", command=self._run_n_point_calculation)
        calc_button.grid(row=row, column=0, columnspan=2, pady=20)
        row += 1
        
        # Results frame
        results_frame = ttk.LabelFrame(main_frame, text="Ergebnisse", padding="10")
        results_frame.grid(row=row, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=10)
        
        # Results text widget with scrollbar
        self.widgets['results_text_n_point'] = tk.Text(results_frame, height=15, width=60, font=('Courier', 9))
        scrollbar = ttk.Scrollbar(results_frame, orient="vertical", command=self.widgets['results_text_n_point'].yview)
        self.widgets['results_text_n_point'].configure(yscrollcommand=scrollbar.set)
        
        self.widgets['results_text_n_point'].grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))

    def _run_n_point_calculation(self):
        """Get values from inputs and perform calculation for N-Point."""
        try:
            # Read values from entries
            inputs = {}
            for key in ["xL", "yL", "xM", "yM", "xR", "yR", "alpha_gon", "beta_gon"]:
                inputs[key] = float(self.widgets[key].get())
            
            # Perform calculation
            result = compute_N_from_LMR_gon(
                inputs["xL"], inputs["yL"], inputs["xM"], inputs["yM"], 
                inputs["xR"], inputs["yR"], inputs["alpha_gon"], inputs["beta_gon"]
            )
            
            # Display results
            self._display_n_point_results(result)
            
        except ValueError as e:
            if "could not convert" in str(e):
                messagebox.showerror("Fehler", "Bitte geben Sie gültige Zahlenwerte ein!")
            else:
                messagebox.showerror("Fehler", str(e))
        except Exception as e:
            messagebox.showerror("Fehler", f"Ein Fehler ist aufgetreten:\n{str(e)}")
    
    def _display_n_point_results(self, result):
        """Display calculation results in the text widget."""
        self.widgets['results_text_n_point'].delete(1.0, tk.END)
        
        output = "=" * 50 + "\n"
        output += "BERECHNUNGSERGEBNISSE\n"
        output += "=" * 50 + "\n\n"
        
        output += "Koordinaten von Punkt N:\n"
        output += f"  x_n = {result['x_n']:.4f}\n"
        output += f"  y_n = {result['y_n']:.4f}\n\n"
        
        output += "-" * 50 + "\n"
        output += "Zwischenwerte:\n"
        output += "-" * 50 + "\n"
        
        for key, value in result.items():
            if key not in ['x_n', 'y_n']:
                if isinstance(value, float):
                    output += f"{key:15s} = {value:12.6f}\n"
                else:
                    output += f"{key:15s} = {value}\n"
        
        self.widgets['results_text_n_point'].insert(1.0, output)


if __name__ == "__main__":
    root = tk.Tk()
    app = GeodesyApp(root)
    root.mainloop()