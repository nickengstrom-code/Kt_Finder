import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import shutil
import warnings
import time
from math import ceil
from scipy.stats import linregress
from scipy.signal import find_peaks, savgol_filter, peak_widths, peak_prominences, filtfilt
from scipy.signal.windows import gaussian
from scipy.ndimage import grey_dilation, grey_erosion
from scipy.optimize import minimize
from pathlib import Path
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

#Hiding warnings from the output box
warnings.filterwarnings('ignore')

#Defining Functions Used by the Program
def fit_circle_kasa(x, y, half='bottom'):
    """
    Fit a circle to (x, y) data using the Kåsa least squares method.
    Returns circle center, radius, and the coordinates of the fitted half-circle.

    Parameters
    ----------
    x, y : array_like
        Input data points (valley or arc).
    half : {'bottom', 'top', 'full'}
        Which portion of the fitted circle to generate for visualization.

    Returns
    -------
    a, b, r : float
        Center coordinates (a, b) and radius r of fitted circle.
    circle_x, circle_y : ndarray
        Coordinates of fitted circle arc.
    """

    x = np.asarray(x)
    y = np.asarray(y)

    #Build and solve the linear system
    A = np.c_[x, y, np.ones_like(x)]
    b = -(x**2 + y**2)
    coeff, *_ = np.linalg.lstsq(A, b, rcond=None)
    A_, B_, C_ = coeff

    #Convert to circle parameters
    a = -A_ / 2
    b = -B_ / 2
    r = np.sqrt((A_**2 + B_**2) / 4 - C_)

    #Generate circle coordinates ---
    if half == 'bottom':
        theta = np.linspace(np.pi/2, 3*np.pi/2, 300)
    elif half == 'top':
        theta = np.linspace(-np.pi/2, np.pi/2, 300)
    else:
        theta = np.linspace(0, 2*np.pi, 400)

    circle_x = a + r * np.cos(theta)
    circle_y = b + r * np.sin(theta)

    return a, b, r, circle_x, circle_y

#Checking fit of estimated circle radius to actual data
def fit_valley_circle(x, y_temp, idx1, radius_of_curvature):
    """
    Fit a circular arc (lower half) around a local valley region centered at x[idx1].

    Parameters
    ----------
    x : np.ndarray
        The full x-coordinate array.
    y_temp : np.ndarray
        The smoothed y data array used for fitting.
    idx1 : int
        Index of the valley minimum in x/y_temp.
    radius_of_curvature : float
        Estimated radius of curvature of the valley.

    Returns
    -------
    result : dict
        Dictionary containing:
            'x_fit', 'y_fit'               - Local data window
            'x_c_opt', 'y_c_opt'           - Optimized circle center
            'circle_x', 'circle_y'         - Circle arc coordinates
            'r_squared'                    - Goodness of fit
            'x_fit_valid', 'y_fit_valid'   - Data points used in R² calculation
            'y_model_valid'                - Corresponding circle predictions
    """
    #Select valley region around minimum
    window_half = int(radius_of_curvature * 3)
    mask = (x >= x[idx1] - window_half) & (x <= x[idx1] + window_half)
    x_fit = x[mask]
    y_fit = y_temp[mask]

    

    #Define circle equation (lower half only)
    def circle_y_lower(x_vals, x_c, y_c, R):
        dx = x_vals - x_c
        inside = R**2 - dx**2
        return y_c - np.sqrt(np.clip(inside, 0, None))  # lower semicircle

    #Define weighted objective function
    def objective(params):
        x_c, y_c = params
        y_model = circle_y_lower(x_fit, x_c, y_c, radius_of_curvature)

        # Gaussian weights centered at valley minimum
        weights = np.exp(-((x_fit - x[idx1])**2) / (2 * (radius_of_curvature / 2)**2))
        weights /= np.max(weights)

        return np.sum(weights * (y_fit - y_model) ** 2)

    #Optimize circle center position
    x0 = [x[idx1], y_temp[idx1] + radius_of_curvature]
    res = minimize(objective, x0, method='Nelder-Mead')
    x_c_opt, y_c_opt = res.x

    #Generate 30° bottom arc for fitting onto the curve
    theta = np.linspace(4 * np.pi / 3, 5 * np.pi / 3, 300)
    circle_x = x_c_opt + radius_of_curvature * np.cos(theta)
    circle_y = y_c_opt + radius_of_curvature * np.sin(theta)

    # Restrict valley data to same x-range as the arc
    arc_mask = (x_fit >= np.min(circle_x)) & (x_fit <= np.max(circle_x))
    x_fit_arc = x_fit[arc_mask]
    y_fit_arc = y_fit[arc_mask]

    # Predict y-values for that x-range
    y_model_arc = circle_y_lower(x_fit_arc, x_c_opt, y_c_opt, radius_of_curvature)
    valid_mask = np.isfinite(y_model_arc)
    x_fit_valid = x_fit_arc[valid_mask]
    y_fit_valid = y_fit_arc[valid_mask]
    y_model_valid = y_model_arc[valid_mask]

    #Compute R²
    if len(y_fit_valid) > 1:
        ss_res = np.sum((y_fit_valid - y_model_valid) ** 2)
        ss_tot = np.sum((y_fit_valid - np.mean(y_fit_valid)) ** 2)
        r_squared = 1 - ss_res / ss_tot if ss_tot != 0 else 0.0
    else:
        r_squared = 0.0

    r_squared = float(max(0.0, min(1.0, r_squared)))

    return {
        "x_fit": x_fit,
        "y_fit": y_fit,
        "x_c_opt": x_c_opt,
        "y_c_opt": y_c_opt,
        "circle_x": circle_x,
        "circle_y": circle_y,
        "r_squared": r_squared,
        "x_fit_valid": x_fit_valid,
        "y_fit_valid": y_fit_valid,
        "y_model_valid": y_model_valid
    }

#Closing profile function
def closing_profile_filter(x, y, minima_indices, r_ball, top_n=5):
    """
    Apply morphological closing to identify sharp valleys and compute effective void depths.

    Parameters:
        x (array): X-coordinates of the profile
        y (array): Smoothed Y profile
        minima_indices (array): Indices of detected minima
        r_ball (float): Ball radius in same units as x
        top_n (int): Number of deepest voids to return

    Returns:
        sharp_minima (array): Indices of sharp valleys
        y_closed (array): Closing profile
        void_depths (dict): Mapping index -> effective void depth
        top_voids (list of tuples): [(x_position, depth), ...] sorted by depth descending
    """
     # Structuring element size in samples (diameter of ball)
    # Compute sampling interval
    dx = np.mean(np.diff(x))
    n = int(np.ceil(r_ball / dx))  # number of points within radius
    
    # Construct circular (ball) footprint in height units
    x_ball = np.linspace(-n, n, 2*n + 1) * dx
    y_ball = np.sqrt(np.clip(r_ball**2 - x_ball**2, 0, None))
    struct_elem = r_ball - y_ball  # "ball" inverted for surface morphology
    
    # Apply dilation then erosion (Equivalent to closing)
    y_dilated = grey_dilation(y, footprint=None, structure=struct_elem)
    y_closed = grey_erosion(y_dilated, footprint=None, structure=struct_elem)
    
    # Effective void depth = closing profile - original profile at minima
    depths = y_closed[minima_indices] - y[minima_indices]
    
    # Keep only sharp valleys (positive depth above threshold)
    threshold = 0.05 * np.max(depths)  # adjustable
    sharp_mask = depths > threshold
    sharp_minima = minima_indices[sharp_mask]
    sharp_depths = depths[sharp_mask]
    
    #Save the void depths to a list
    void_depths = np.array(-sharp_depths)
    orig_depths = y_smoothed[sharp_minima]
    
    func_depths = np.maximum(void_depths, orig_depths)

    # Rank and get top-N deepest voids
    #ranked = sorted(void_depths.items(),key=lambda t: t[1], reverse=True)
    ranked = sorted([(i, d) for i, d in zip(sharp_minima, void_depths)],
                    key=lambda t: t[1], reverse=False)
    top_voids = ranked[:top_n]
     

    return sharp_minima, y_closed, void_depths, top_voids, ranked, func_depths

#Automatically Calculates the required roughness paramters from the data
def calculate_roughness(x, y, num_segments=5):
    """
    Calculate Ra, Ry, and Rz from a surface profile.

    Parameters:
        x (array): X-coordinates of the profile.
        y (array): Y-coordinates (height values) of the profile.
        num_segments (int): Number of equal segments for Rz calculation (default=5).

    Returns:
        tuple: (Ra, Ry, Rz)
    """
  

    #Ra: arithmetic mean roughness
    #5 Section Approach to Determining Ra
    n = len(y)
    section_length = n // num_segments
    ra_values = []

    for i in range(num_segments):
        start = i * section_length
        # Handle the last section to include all points
        end = (i + 1) * section_length if i < num_segments - 1 else n
        section = y[start:end]
        mean_section = np.mean(section)
        ra_section = np.mean(np.abs(section - mean_section))
        ra_values.append(ra_section)
    Ra = np.mean(ra_values)   
    
    #Ry: max peak-to-valley across entire profile
    Ry = np.max(y) - np.min(y)

    #Rz: average of max peak-to-valley per segment
    segment_length = len(y) // num_segments
    rz_values = []
    for i in range(num_segments):
        start = i * segment_length
        end = start + segment_length
        segment = y[start:end]

        if len(segment) > 0:
            Pmax = np.max(segment)
            Vmin = np.min(segment)
            rz_values.append(Pmax - Vmin)

    Rz = np.mean(rz_values) if rz_values else np.nan

    return Ra, Ry, Rz

#Significant Curvature Region finder
def find_positive_region(y_data, center_idx, mode='strict', tol=0.05, sign_change_sensitivity=0.0):
    """
    Finds contiguous region around a center index where curvature (y_data)
    remains positive (strict) or near-positive (range), and stops if curvature
    changes sign sharply (inflection point).

    Parameters
    ----------
    y_data : array_like
        1D array of curvature or similar values (e.g. second derivative).
    center_idx : int
        Index of the central point (valley bottom or curvature maximum).
    mode : {'strict', 'range'}, optional
        - 'strict': includes only y_data > 0
        - 'range': includes values >= -tol (allows slightly negative)
    tol : float, optional
        Allowed deviation below zero for 'range' mode. Default 0.05.
    sign_change_sensitivity : float, optional
        Minimum absolute jump in curvature sign (Δy) considered an inflection.
        Set to 0 to disable sign-based stopping entirely.

    Returns
    -------
    left, right : int
        Indices marking the left and right boundaries of the significant curvature region.
    """

    y_data = np.asarray(y_data)
    n = len(y_data)
    if not (0 <= center_idx < n):
        raise IndexError("center_idx out of bounds.")

    if mode not in ("strict", "range"):
        raise ValueError("mode must be 'strict' or 'range'")

    # Define the inclusion condition
    if mode == "strict":
        condition = lambda val: val > 0
    else:  # 'range'
        condition = lambda val: val > -tol

    #Initialize
    left = center_idx
    right = center_idx

    # Helper function to check for inflection between two values
    def is_inflection(y1, y2):
        if np.sign(y1) != np.sign(y2):
            return abs(y2 - y1) > sign_change_sensitivity
        return False

    # Searching left
    while left > 0 and condition(y_data[left - 1]):
        if sign_change_sensitivity > 0 and is_inflection(y_data[left], y_data[left - 1]):
            break
        left -= 1

    #Searching right
    while right < n - 1 and condition(y_data[right + 1]):
        if sign_change_sensitivity > 0 and is_inflection(y_data[right], y_data[right + 1]):
            break
        right += 1

    return left, right

#Function to Check the Name of the File       
def get_user_parameters():
    """Opens a Tkinter GUI to collect user parameters and returns them as a dictionary."""

    def browse_files():
        files = filedialog.askopenfilenames(
            title="Select One or More CSV Files",
            filetypes=(("CSV Files", "*.csv"), ("All Files", "*.*"))
        )
        if files:
            file_names_var.set(", ".join(files))  # display comma-separated names in box
            selected_files.clear()
            selected_files.extend(list(files))  # store as list

    def submit():
        nonlocal params

        if not selected_files:
            messagebox.showerror("Error", "Please select at least one CSV file.")
            return

        try:
            params = {
                "file_names": selected_files.copy(),  # list of files
                "thresh": float(thresh_var.get()),
                "level_method": level_var.get(),
                "max_curvature": float(max_curvature_var.get()),
                "max_radius": float(max_radius_var.get()),
                "width": [0, int(width_var.get())],
                "r_ball": int(r_ball_var.get()),
                "n_top": int(n_top_var.get()),
                "graphs" : graphs_var.get()
            }

            root.destroy()  # close GUI

        except ValueError:
            messagebox.showerror("Error", "Please enter valid numeric values.")

    # --- GUI Setup ---
    root = tk.Tk()
    root.title("Parameter Input GUI")
    root.geometry("420x700")
    root.resizable(False, False)

    selected_files = []
    params = None

    # Tkinter Variables (with defaults)
    file_names_var = tk.StringVar()
    thresh_var = tk.StringVar(value="0.2")
    level_var = tk.StringVar(value="Gaussian")
    graphs_var = tk.StringVar(value="No")
    max_curvature_var = tk.StringVar(value=str(2))
    max_radius_var = tk.StringVar(value="100")
    width_var = tk.StringVar(value="500")
    r_ball_var = tk.StringVar(value="100")
    n_top_var = tk.StringVar(value="25")

    padding = {'padx': 10, 'pady': 5}

    # File selection
    tk.Label(root, text="Select CSV File(s):").pack(**padding)
    tk.Entry(root, textvariable=file_names_var, width=50, state="readonly").pack(**padding)
    tk.Button(root, text="Browse Files", command=browse_files).pack(**padding)

    # Numeric entries
    tk.Label(root, text="Threshold (thresh):").pack(**padding)
    tk.Entry(root, textvariable=thresh_var).pack(**padding)

    # Level method dropdown
    tk.Label(root, text="Select Leveling Method:").pack(**padding)
    level_dropdown = ttk.Combobox(
        root,
        textvariable=level_var,
        values=["Gaussian", "Linear", "No-Leveling"],
        state="readonly"
    )
    level_dropdown.pack(**padding)
    level_dropdown.current(0)
    
    tk.Label(root, text="Generate Circle Fit Graphs?:").pack(**padding)
    output_graphs = ttk.Combobox(
        root,
        textvariable=graphs_var,
        values=["No", "Yes"],
        state="readonly"
    )
    output_graphs.pack(**padding)
    output_graphs.current(0)

    tk.Label(root, text="Profilometer Tip Radius:").pack(**padding)
    tk.Entry(root, textvariable=max_curvature_var).pack(**padding)

    tk.Label(root, text="Largest Radius of Consideration:").pack(**padding)
    tk.Entry(root, textvariable=max_radius_var).pack(**padding)

    tk.Label(root, text="Width Max:").pack(**padding)
    tk.Entry(root, textvariable=width_var).pack(**padding)

    tk.Label(root, text="R Ball:").pack(**padding)
    tk.Entry(root, textvariable=r_ball_var).pack(**padding)
    
    tk.Label(root, text="Number of Deepest Valleys to Highlight:").pack(**padding)
    tk.Entry(root, textvariable=n_top_var).pack(**padding)

    tk.Button(root, text="Run", command=submit, bg="#4CAF50", fg="white").pack(pady=15)

    root.mainloop()

    return params

params = get_user_parameters()

#Pulling Data Out of Params
thresh = params["thresh"]
level = params["level_method"]

#Curvature Paramaters
max_curvature = 1 / params["max_curvature"]
min_curvature= 1/ params["max_radius"]


#Peak Finding Variables
prominence_value = max_curvature #Ensures that the minimum prominence will be larger that the radius of the profilometer tip
width = params["width"]

#Running the Closing profile over the identified valleys
r_ball = params["r_ball"]
n_top = params["n_top"]

#Get the paths of all inputted files
file_paths = params["file_names"] 

#Data Aquisition

for file in file_paths:
    #Start Time Measurement
    start = time.process_time()
    
    #Read in the CSV File to a Panda Datafram
    data = pd.read_csv(file)
        
    #Pull out the path
    file_path = Path(file)
    
    # Define the folder location and file name
    folder_location = str(file_path.parent)
    
   #Creating New Folder Path to save results
    file_name = file_path.stem
    print(f'\n-------------------------------------\nFile Name: {file_name}')
    
    New_Folder_Path = file_path.parent.joinpath( file_name + " Kt Data")
    
    #Create New folder or empty the folder if it already exists
    if New_Folder_Path.exists() and New_Folder_Path.is_dir():
        shutil.rmtree(New_Folder_Path)
    New_Folder_Path.mkdir(parents=True, exist_ok=True)             
    
    # Extract the X and Y data from the provided CSV Files
    x = data.iloc[:, 0].values
    y = data.iloc[:, 1].values
    
    
    #Taking in Profilometer Data
    try:
        ra=data.iloc[0 , 3]
        ry=data.iloc[0 , 4]
        rz=data.iloc[0 , 5]
    except IndexError:
        print("\nRoughness Data not Found, Calculating Internally ...")
        ra, ry, rz = calculate_roughness(x, y, num_segments=5) #Internally calculate roughness data if none provided
    
    if level in ['Linear']:
        # Perform linear regression to get the slope and intercept of the trend line
        slope, intercept, _, _, _ = linregress(x, y)
    
        # Compute the trend line (y_trend) based on the slope and intercept
        y_trend = slope * x + intercept
    
        # Subtract the trend from the original
        y_corrected = y - y_trend
        compare=0
    elif level in ['Gaussian']:
        # Perfrom Gaussian Leveling
        dx = x[1]-x[0]       # sampling interval in um 
        lc = 800             # cutoff wavelength (um)
        
        #Gaussian filter 
        sigma = (lc / (2 * np.pi)) / dx   # standard deviation
        N = int(np.ceil(6 * sigma))       # kernel length (~±3σ)
        if N % 2 == 0:
            N += 1                        # make it odd for symmetry
    
        # Build and normalize Gaussian kernel
        g = gaussian(N, sigma)
        g /= np.sum(g)
    
        # Apply zero-phase Gaussian filter
        y_trend = filtfilt(g, 1, y)
    
        # Subtract the trend from the original
        y_corrected = y - y_trend
        compare = 1
    else:
        #No Leveling Will be performed for raw data
        y_corrected = y
        compare = 2
            
    #Very slight Savitzky-Golay pre-smoothing to reduce noise
    y_smoothed = savgol_filter(y_corrected, window_length=10, polyorder=3)
    
    # Identify local minima using selective criteria 
    minima_indices, prominence_data = find_peaks(-y_smoothed[100:len(y_smoothed)-100], prominence=prominence_value, height=0, distance=40, width=width) #
    minima_indices=minima_indices+100 #Buffer to exclude any irregularities at start of data
    for key in ['left_bases','left_ips','right_bases', 'right_ips']:
         prominence_data[key] += 100 #Fix the prominence data
    
    #Run the closing profile to identify deep voids
    sharp_minima, closing_profile, void_depths, top_voids, ranked, func_depths= closing_profile_filter(
        x, y_smoothed, minima_indices, r_ball, top_n = n_top)
    
    #Prepare arrays to store valley information [x, y, curvature, radius]
    curvature_table = [] 
    out_of_range_radii = []
    asymmetric_minima = []
    r_table=[]   
    
    #Setting the minima list to iterate through
    minima_type = sharp_minima
    
    for count, idx in enumerate(minima_type):
            
            #Pull in previously calculated valley depths and use them to calculate valley width
            valley_data = peak_prominences(-y_smoothed, [minima_type[count]])
            prominence = func_depths[count]
            valley_data[0][0]=abs(prominence)
            valley_width_data = peak_widths(-y_smoothed, [minima_type[count]], rel_height = 0.6, prominence_data= valley_data)
            valley_width = valley_width_data[0][0]    

            if valley_width >  1/ (2 * min_curvature):
                half_valley = ceil(valley_width*0.4)
            else:
                half_valley = ceil(valley_width*0.3)
            start_window_length = max(3, half_valley)
            if start_window_length % 2 == 0:
                start_window_length += 1
                
            # Extract data within ±50 µm around the minimum
            curve_window_mask = (x >= x[idx] - abs(valley_width//2 +4) ) & (x <= x[idx] + (valley_width//2 + 4) ) 
            x_window = x[curve_window_mask]
            y_window = y_smoothed[curve_window_mask]
            
            
            window_length = start_window_length
            
                     
            #conditional Smoothing of the Valley                        
            y_window_smoothed = savgol_filter(y_window, window_length=window_length, polyorder=2)
            #y_window_smoothed = y_window
                         
            # Compute derivatives using local smoothed data
            dy_dx = np.gradient(y_window_smoothed, x_window)
            d2y_dx2 = np.gradient(dy_dx, x_window)
    
            # Curvature calculation as before...
            curvatures = d2y_dx2 / (1 + dy_dx**2)**(3 / 2)
            
            
            local_center_idx = np.argmin(np.abs(x_window - x[idx]))
            pos_left, pos_right = find_positive_region(curvatures, local_center_idx)
            
            positive_region = curvatures[pos_left:pos_right + 1]
            min_idx_in_region = np.argmax(positive_region)
            min_idx = pos_left + min_idx_in_region
            
            valley_region = y_window_smoothed[pos_left:pos_right + 1]
            min_idx_yregion = np.argmin(valley_region)
            min_idx_y = pos_left + min_idx_yregion
            
            
            x_min = x_window[min_idx_y]  # actual x-coordinate of smoothed minimum
            idx1= np.argmin(np.abs(x - x_min))
            
            #Filtering Curvatures
            threshold = thresh * np.max(positive_region)
            
            curvature = float(np.mean([entry for entry in positive_region if entry >= threshold]))
            
            # Estimate the radius of curvature accounting for outliers
            if curvature >= min_curvature and curvature <= max_curvature:
                radius_of_curvature = 1 / abs(curvature)
            else:
                radius_of_curvature = 1 / abs(curvature)
                out_of_range_radii.append([x[idx1], y_window_smoothed[min_idx_y], curvature, radius_of_curvature])
                continue    
            
                  
            #Update the plot y data with the individually smoothed data
            y_window_smoothed=y_window_smoothed[0:np.count_nonzero(curve_window_mask)]
            
            y_temp = np.copy(y_smoothed)
            y_temp[curve_window_mask]=y_window_smoothed
            
            valley_fit_data = fit_valley_circle(x, y_temp, idx1, radius_of_curvature)
            r_squared = valley_fit_data["r_squared"]
            
            if r_squared <= 0.2:                
                asymmetric_minima.append([x[idx1], y_window_smoothed[min_idx_y], curvature, radius_of_curvature, r_squared])
                left, right = find_positive_region(curvatures, min_idx, mode = 'range' , tol = 0.01, sign_change_sensitivity = 0.05)
                if (right - left) <= len(valley_region):
                    left, right = pos_left, pos_right
                a, b, radius_of_curvature, circle_x, circle_y = fit_circle_kasa(x_window[left:right], y_window_smoothed[left:right] , half = 'full')
                valley_fit_data = fit_valley_circle(x, y_smoothed, idx1, radius_of_curvature)
                #if valley_fit_data['r_squared'] <= 0:
                    #continue
            else:                               
                theta = np.linspace(0, 2*np.pi, 300)  # only bottom arc
                circle_x = valley_fit_data["x_c_opt"] + radius_of_curvature * np.cos(theta)
                circle_y = valley_fit_data["y_c_opt"] + radius_of_curvature * np.sin(theta)
            
            # Store results in the curvature table
            curvature_table.append([x[idx1], y_window_smoothed[min_idx_y], curvature, radius_of_curvature, valley_fit_data["r_squared"]])
            running_average = float(np.mean([entry[3] for entry in curvature_table if entry[3] != np.inf]))
            kt_tension=1+2*(ra/running_average)*(ry/rz)
            kt_shear=1+(ra/running_average)*(ry/rz)
            curvature_table[len(curvature_table)-1].append(kt_tension)
            curvature_table[len(curvature_table)-1].append(kt_shear)    
            
            if params["graphs"] in ["Yes"]:
                plt.figure(figsize=(8, 8))
                plt.plot(valley_fit_data["x_fit"], valley_fit_data["y_fit"], label='Window Data', color='blue')        
                plt.plot(circle_x, circle_y, linestyle='--', color='orange', label='Estimated Circle of Best Fit')
                plt.plot(x[idx1], y_temp[idx1], 'rx', markersize=10, label='Local Minimum')
                plt.plot([],[],'', label = f'Confidence of fit = {valley_fit_data["r_squared"]:.3f}' )
                plt.axis([x[idx1]-radius_of_curvature*3,x[idx1]+radius_of_curvature*3,y_window_smoothed[min_idx]-radius_of_curvature*2,y_window_smoothed[min_idx]+4*radius_of_curvature])
                plt.title(f'Local Minimum at x = {x[idx1]},{y_temp[idx1]}')
                plt.xlabel('Scan Distance(µm)')
                plt.ylabel('Height(µm)')
                plt.axhline(0, color='gray', linewidth=0.5, linestyle='--')  # Reference line at y=0
                plt.legend(loc = 'upper right')
                plt.grid(True)
                plt.savefig(file_path.joinpath(New_Folder_Path, f"Min at {x[idx1]},{y_smoothed[idx1]}.png"), dpi=300, bbox_inches="tight")# Save the plot as a png
                plt.show()  # Show the individual plot for the minimum
    elapsed = (time.process_time()-start)              
    
            
    if compare == 2:
        # Plot the original Data Only
        plt.figure(figsize=(15, 10))
        
        # Original data plot
        plt.subplot(2, 1, 1)
        plt.plot(x, y_corrected, label='Raw Data', linestyle='-', color='green')
        plt.axhline(0, color='red', linewidth=0.5, linestyle='--')  # Reference line at y=0
        plt.xlabel('Traverse Distance (µm)', weight = 'bold')
        plt.ylabel('Depth (µm) (Corrected)', weight = 'bold')
        plt.title('Raw Data', weight = 'bold')
        plt.legend(bbox_to_anchor = (1, 1), loc='upper left')
        plt.grid(False)
        
        plt.subplot(2, 1, 2)
        plt.plot(x, y_smoothed, label='Raw Data', color='green')
        plt.plot(x, closing_profile, label='Closing Profile', color='black', linestyle='--')

        for xv, depth in top_voids:
            plt.axvline(x=x[xv], color='red', linestyle=':', alpha=0.6)   # vertical line
            plt.scatter(x[xv], y_smoothed[xv], 
                        color='red', s=50, label= f'Deepest {n_top} Valleys ' if xv == top_voids[0][0] else "")
    
        plt.xlabel('Traverse Distance (µm)', weight = 'bold')
        plt.ylabel('Depth (µm) (Corrected)', weight = 'bold')
        plt.title('Closing Profile with Deepest Valleys', weight = 'bold')
        plt.legend(bbox_to_anchor = (1, 1), loc='upper left')
        plt.grid(False)
        
    else:  
        #Plot original and various corrected data     
        plt.figure(figsize=(12, 10))
    
        #Original data with trend line
        plt.subplot(3, 1, 1)
        plt.plot(x, y, label='Original Data', linestyle='-', color='blue')
        plt.plot(x, y_trend, label='Trend Line', color='red', linestyle='--')
        plt.xlabel('Traverse Distance (µm)', weight = 'bold')
        plt.ylabel('Depth (µm) ', weight = 'bold')
        plt.title('Original Data with Trend Line', weight = 'bold')
        plt.legend(bbox_to_anchor = (1, 1), loc='upper left')
        plt.grid(False)
    
        #Corrected data
        plt.subplot(3, 1, 2)
        plt.plot(x, y_corrected, label='Leveled Data', linestyle='-', color='green')
        plt.axhline(0, color='red', linestyle='--', label = 'Trend Line')
        plt.xlabel('Traverse Distance (µm)', weight = 'bold')
        plt.ylabel('Depth (µm) (Corrected)', weight = 'bold')
        plt.title('Leveled Data', weight = 'bold')
        plt.legend(bbox_to_anchor = (1, 1), loc='upper left')
        plt.grid(False)
    
        #Corrected data + closing profile + sharp voids
        plt.subplot(3, 1, 3)
        plt.plot(x, y_smoothed, label='Leveled Data', color='green')
        plt.plot(x, closing_profile, label='Closing Profile', color='black', linestyle='--')
    
        # Highlight deepest sharp voids
        for xv, depth in top_voids:
            plt.axvline(x=x[xv], color='red', linestyle=':', alpha=0.6)   # vertical line
            plt.scatter(x[xv], y_smoothed[xv], 
                        color='red', s=50, label=f'Deepest {n_top} Valleys ' if xv == top_voids[0][0] else "")
     
        plt.xlabel('Traverse Distance (µm)', weight = 'bold')
        plt.ylabel('Depth (µm) (Corrected)', weight = 'bold')
        plt.title('Closing Profile with Deepest Valleys', weight = 'bold')
        plt.legend(bbox_to_anchor = (1, 1), loc='upper left')
        plt.grid(False)
    
    # Adjust layout and show the plots
    plt.tight_layout()
    plt.savefig(file_path.joinpath(New_Folder_Path, "Profilometer Plot.png"), dpi=300, bbox_inches="tight")# Save the plot as a png
    plt.show()
    
    #Finding the Max and Min Radius
    radii=[entry[3] for entry in curvature_table if entry[3] != np.inf or 0]
    
    if len(radii) == 0:
        print("\nNo Radii Within Bounds Found, Unable to Calculate Kt")
    else:    
        #Calculating p and Kt values for all valleys
        radius_std = float(np.std([entry[3] for entry in curvature_table]))
        average_radius = float(np.mean([entry[3] for entry in curvature_table]))  # Exclude infinite values
        lower_bound = max(1/max_curvature, average_radius - 2*radius_std)
        upper_bound = min(1/min_curvature, average_radius + 2*radius_std)
        filtered_radii = [entry[3] for entry in curvature_table if lower_bound <= entry[3] <= upper_bound]
        new_average_radius = float(np.mean([entry[3] for entry in curvature_table if lower_bound <= entry[3] <= upper_bound]))
        average_r_square = float(np.mean([entry[4] for entry in curvature_table]))
        print(f"\nAverage Radius of Curvature: {new_average_radius:.2f}")
    
        #Calculating Kt for the average radius
        kt_tension=1+2*(ra/new_average_radius)*(ry/rz)
        kt_shear=1+(ra/new_average_radius)*(ry/rz)
        print(f"\nKt Tension(Overall): {kt_tension:.2f}")
        print(f"\nKt Shear(Overall): {kt_shear:.2f}")
    
        if len(radii) < 6:
            print("\nWARNING: Number of Radii Found is Low, kt value may be inaccurate")
    
        if len(radii) <= len(out_of_range_radii):
            print("\nWARNING: Number of Radii outside bounds high, kt value may be inaccurate")
    
        #Creating an Output Excel File
        length = len(radii)
        
        # Create lists where the first element is the scalar value, and the rest are None
        kt_tension_list = [kt_tension] + [None] * (length - 1)
        kt_shear_list = [kt_shear] + [None] * (length - 1)
        
        # Combine into a dictionary for creating the DataFrame
        results = {
            'Radius': radii,
            'Kt Tension': kt_tension_list,
            'Kt Shear': kt_shear_list,
        }
        
        # Create and save the DataFrame
        results_df = pd.DataFrame(results)
        results_df.to_excel(file_path.joinpath(New_Folder_Path, file_name+" Results.xlsx"))
        #Display the curvature table as a pandas DataFrame
        curvature_df = pd.DataFrame(curvature_table, columns=['X', 'Y', 'Curvature', 'Radius', 'r_squared','kt_tension' , 'kt_shear'])
        #Save The dataframe as an excel file
        curvature_df.to_excel(file_path.joinpath(New_Folder_Path, "Tabulated Values.xlsx"))
        # Plot the curvature table using matplotlib
        fig, ax = plt.subplots(figsize=(8, 4))  # Create a new figure for the table
        ax.axis('off')  # Turn off the axis
        # Render the DataFrame as a table in the plot
        table = ax.table(cellText=curvature_df.values, colLabels=curvature_df.columns, loc='center', cellLoc='center')
        plt.savefig(file_path.joinpath(New_Folder_Path, "Tabulated Values.png"), dpi=300, bbox_inches="tight")# Save the plot as a png
        plt.show()
    print (f'\nProcessing Time : {elapsed:.2f} seconds')