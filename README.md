# Surface Stress Concentration Factor (Kt) Finder

## Background

This code was created as part of a research project conducted at the **University of Washington** under **Dr. Dwayne Arola**.  
The goal of this project is to provide a **programmatic implementation** of the **Arola–Ramulu method** for estimating *surface stress concentration factors (Kt)* in **additively manufactured (AM)** components.  

The program allows for consistent and accurate determination of Kt for both *as-printed* and *treated* surfaces, helping evaluate the effectiveness of surface treatments.

For more information about the Arola–Ramulu method, see references [1][2].

---

## High-Level Overview

The program:

1. **Reads profilometer data**  
   Accepts `x, y` surface data and roughness parameters (`Ra`, `Ry`, `Rz`) from profilometer software.  
   - If roughness data is missing, these can be internally calculated.  
   - Note: Calculated values highly depend on the selected leveling type[3].

2. **Levels the data**  
   Based on user input:
   - Gaussian leveling  
   - Linear leveling  
   - No leveling

3. **User parameters**  
   Users can input or modify:
   - Profilometer tip radius  
   - Maximum valley width  
   - Number of valleys to analyze  
   - Maximum valley radius  
   - Leveling type  

   *(Defaults are provided based on optimal values found by the authors.)*

4. **Valley detection**  
   Identifies valley locations in the profilometer scan and uses a *closing profile* approach [4] to find the **critical valleys** that most contribute to stress concentration and potential failure points.

5. **Data smoothing and circle fitting**  
   For each valley:
   - Applies a conditional **Savitzky–Golay filter** based on the valley width to remove noise without eliminating small valleys.  
   - Fits a circle to the valley to estimate the **valley radius**.

6. **Kt Calculation**  
   Combines the input (or calculated) roughness values (`Ra`, `Ry`, `Rz`) with estimated valley radii to compute **Kt** using the **Arola–Ramulu equation**(Eq. 1).
```math
   \bar{K}_t=1+n (\frac{R_a}{ρ ̅ })(\frac{R_y}{R_z} )                   \qquad \qquad \qquad (1)
```
---

## How to Use

### 1. Prepare Input Data

Use the provided **Kt Program Data Template** to format your profilometer data correctly.

**Requirements:**
- `X` and `Y` in **micrometers (µm)**  
- File format: **CSV**
- If you lack roughness values, leave `Ra`, `Ry`, and `Rz` blank.  
  *(Recommended: use Gaussian smoothing if you have the program estimate Ra, Ry, Rz.)*
![Excel Template](https://github.com/nickengstrom-code/Kt_Finder/blob/main/Excel%20Format.png)
---

### 2. Run the Program

It is recommended to use **Jupyter Notebook** or **Spyder**, though any Python environment should work.

When executed, the program:
- Opens a **GUI** to select input files and specify parameters.  
- Parameters are **prefilled with optimal defaults** (based on LPBF surfaces) but can be modified.

![GUI View](https://github.com/nickengstrom-code/Kt_Finder/blob/main/Gui%20Example.png)
|Variable	|Recommended Range|
|---|---|
|Threshold 	|0-1|
|Profilometer Tip Radius	|Relative to Equipment|
|Largest Radius of Consideration	|50-250|
|Width Max	|500-1000|
|R Ball	|50-250|
|Number of Deepest Valleys to Highlight	|6-25|


---

### 3. Output

The program generates:

- Original and leveled surface data  
- Closing profile plots with **highlighted valleys**  
- Estimated **Kt values** for **shear** and **tension** conditions  

It can also optionally produce **circle fit graphs** for each analyzed valley — to visually inspect fit quality.  
Enable setting **Generate Circle Fit Graphs** to **Yes**



---

## Example Output

> Multiple files can be processed consecutively.  
> Output includes leveled data visualizations, critical valley identification, and estimated Kt results for shear and tension conditions.
> 
![Example Output](https://github.com/nickengstrom-code/Kt_Finder/blob/main/Example%20Output.png)
---

## References

1. D. Arola and C. Williams, “Estimating the fatigue stress concentration factor of machined surfaces,” *International Journal of Fatigue*, vol. 24, no. 9, pp. 923–930, Sept. 2002.  
   DOI: [10.1016/S0142-1123(02)00012-9](https://doi.org/10.1016/S0142-1123(02)00012-9)

2. D. Arola and M. Ramulu, “An Examination of the Effects from Surface Texture on the Strength of Fiber Reinforced Plastics,” *Journal of Composite Materials*, vol. 33, no. 2, pp. 102–123, Jan. 1999.  
   DOI: [10.1177/002199839903300201](https://doi.org/10.1177/002199839903300201)

3. D. Nečas, M. Valtr, and P. Klapetek, “How levelling and scan line corrections ruin roughness measurement and how to prevent it,” *Scientific Reports*, vol. 10, no. 1, p. 15294, Sept. 2020.  
   DOI: [10.1038/s41598-020-72171-8](https://doi.org/10.1038/s41598-020-72171-8)

4. M. Zecchino and M. Malburg, “Predicting fatigue failure: a better texture parameter.” *Digital Metrology*, Feb. 10, 2025.  
   [Online Article](https://digitalmetrology.com/tutorials/predicting-fatigue-failure-better-texture-parameter/)
