#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 14 16:14:46 2024

@author: alejandrosoto
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import scienceplots
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.interpolate import interp1d
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
import scienceplots
import matplotlib.colors as mcolors
import matplotlib.cm as cm


plt.style.use(['science', 'no-latex', 'bright'])

# ==============
# 1) Dictionary for custom names
# ==============
custom_names = {
    "L-4": "LS16",
    "L0": "LS12",
    "L1": "LS11",
    "L2": "LS10",
    "L3": "LS9",
    "L4": "LS8",
    "L7": "LS5",
    "L11": "LS1",
    "Y+L11": "Y+LS1",
    "Y+L0": "Y+LS12",
    "Y": "Y"

}

def get_custom_label(label):
    """
    Helper function:
    1) Removes underscores (and any other unwanted characters).
    2) Looks up the cleaned label in `custom_names`.
    3) Returns the mapped name if found; otherwise returns the original cleaned label.
    """
    # Example transformation: "L_-4" -> "L-4", "L_0" -> "L0", etc.
    cleaned = label.replace("_", "")
    # Now find a match in custom_names
    return custom_names.get(cleaned, cleaned)

# Ruta base donde se encuentran los archivos
base_paths = {
    "L_-4": "./L_-4/hb",
    "L_0": "./L_0/hb",
    "L_1": "./L_1/hb",
    "L_2": "./L_2/hb",
    "L_3": "./L_3/hb",
    "L_4": "./L_4/hb",
    "L_7": "./L_7/hb",
    "L_11": "./L_11/hb",
    #"Y+L11": "./Y+L11/hb",
    #"Y+L0": "./Y+L0/hb",
    "Y": "./Y_0/hb"
}

temps = [
    14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48,
    50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84,
    86, 88, 90, 92
]

norm = 25 * (28 + 14)  # Normalización

def shuffled_sterr(_in):
    _chunk = 20
    N = _in.size // _chunk
    if N == 0:
        return 0
    out = 0
    for i in range(N):
        _ids = np.random.randint(0, high=_in.size, size=_chunk)
        out += np.std(_in[_ids])
    return np.sqrt(out / N) / np.sqrt(N)

def extrapolate_bulk_array(x):
    temp = x / (1 - x)
    return (1 + 1 / (2 * temp)) - np.sqrt((1 + 1. / (2 * temp)) ** 2 - 1)

def find_temp_for_fraction(temps, fractions, target_fraction=0.5):
    """Interpolar linealmente para encontrar la temperatura en la cual la fracción alcanza un valor dado."""
    f = interp1d(fractions, temps, kind='linear', fill_value="extrapolate")
    return f(target_fraction)

# Configurar estilo y fuentes
plt.rcParams.update({
    'font.size': 10,
    'font.family': 'serif',
    'axes.labelsize': 10,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 8,
})

marker_size = 10

fig = plt.figure(figsize=(9, 9))
grid = plt.GridSpec(3, 2, hspace=0.4, wspace=0.3)

colors = [
"#313695", "#416aaf", "#74add1", "#a3d3e6", "#e0f3f8", "#fffbb9", "#fee090", "#fca55d", "#eb5a3a", "#ce2827", "#c6171b", "#7f0d0b"
][::-1]
temp_at_05_values = []
temp_at_05_errors = []
labels = []
norm = None





ax1 = fig.add_subplot(grid[0, 0])
# ==============
# 2) Plot 1: Fraction unbounded vs. T
# ==============
for i, (label, base_path) in enumerate(base_paths.items()):
    ave_data = []
    std_data = []

    # Leer el archivo para la temperatura 14 y encontrar valor máximo para normalizar
    hb_file_14 = os.path.join(base_path, "hb_list.14.dat")
    if os.path.exists(hb_file_14):
        data_14 = np.loadtxt(hb_file_14)
        norm = np.max(data_14)  # Asignar norm como el máximo del primer archivo
    else:
        print(f"Advertencia: El archivo {hb_file_14} no existe.")

    for temp in temps:
        hb_file = os.path.join(base_path, f"hb_list.{temp}.dat")
        
        if os.path.exists(hb_file):
            data = np.loadtxt(hb_file)
            data = np.nan_to_num(data, nan=0.0)
            
            store_data = data[-data.size // 10:] / norm if norm != 0 else data[-data.size // 10:]
            
            ave_data.append(np.average(store_data))
            std_data.append(shuffled_sterr(store_data))
        else:
            ave_data.append(np.nan)
            std_data.append(np.nan)

    ave_data = np.nan_to_num(np.array(ave_data), nan=0.0)
    std_data = np.nan_to_num(np.array(std_data), nan=0.0)

    extr = extrapolate_bulk_array(ave_data)
    color = colors[i % len(colors)+1]
    
    # ============== 
    # Use custom function to rename label 
    # ==============
    custom_label = get_custom_label(label)  
    clean_label = label.replace("_", "")
    
    # Plot
    ax1.errorbar(
        temps, extr, yerr=std_data,
        fmt='o-', color=color, markeredgecolor='k', capsize=5,
        label=clean_label
    )

    # Encontrar T a la cual la fracción alcanza 0.5
    temp_at_05 = find_temp_for_fraction(temps, extr, target_fraction=0.5)
    temp_at_05_values.append(temp_at_05)
    temp_at_05_errors.append(np.std(std_data))  # Asumir std_data como error
    labels.append(label)
    


ax1.set_ylabel("Fraction unbounded", fontsize=10)
ax1.set_xlabel("Temperature ($^{o}$C)", fontsize=10)
ax1.set_ylim(-0.1, 1.1)
ax1.set_xticks(np.arange(10, 110, 10))
ax1.set_yticks(np.arange(0, 1.2, 0.2))
ax1.text(-0.04, 1.07, 'a)', transform=ax1.transAxes,
              fontsize=12, fontweight='bold', va='top', ha='right')
ax1.set_xlim(10, 100)
ax1.text(75, 1.05, '100 mM NaCl',  
              fontsize=10, color='black', ha='left', va='top', 
              bbox=dict(boxstyle='square,pad=0.1', facecolor='white',
                        edgecolor='none', alpha=0.1))
ax1.text(75, 1.01, r'$1\,\mu\mathrm{M}$',  
              fontsize=10, color='black', ha='left', va='top', 
              bbox=dict(boxstyle='square,pad=0.1', facecolor='white',
                        edgecolor='none', alpha=0.1))

# Agregar la línea horizontal
ax1.axhline(y=0.5, color='black', alpha=0.4, linestyle='--', linewidth=1, xmin=0.4, xmax=1)



# En vez de meter "labels=..." a legend, hazlo así:
# 1) Obtenemos los handles y las viejas etiquetas
handles, old_labels = ax1.get_legend_handles_labels()

# 2) Creamos las nuevas etiquetas (con Tm, etc.)
new_labels = []
for i, old_lbl in enumerate(old_labels):
    # Aquí supón que la i-ésima serie corresponde al i-ésimo 'temp_at_05_values'
    # y que tienes un método para mapear esa 'old_lbl' a tu 'custom_label'.
    custom_lbl = get_custom_label(labels[i])
    new_labels.append(
        f'{custom_lbl} ($T_{{m}}={temp_at_05_values[i]:.2f} ^{{\circ}}C$)'
    )

# 3) Reconstruimos la leyenda con los mismos 'handles' pero textos nuevos
ax1.legend(handles, new_labels, loc='lower left')


# ==============
# Zoom inset on first subplot
# ==============
axins1 = fig.add_axes([0.38, 0.8, 0.08, 0.08])
for i, (label, base_path) in enumerate(base_paths.items()):
    ave_data = []
    for temp in temps:
        hb_file = os.path.join(base_path, f"hb_list.{temp}.dat")
        
        if os.path.exists(hb_file):
            data = np.loadtxt(hb_file)
            if data.size > 0:
                store_data = data[-data.size // 10:] / norm if norm != 0 else data[-data.size // 10:]
                ave_data.append(np.average(store_data))
            else:
                ave_data.append(np.nan)
        else:
            ave_data.append(np.nan)

    ave_data = np.nan_to_num(np.array(ave_data), nan=0.0)
    axins1.plot(temps, ave_data, label=label,
                color=colors[i % len(colors)])

axins1.set_xlim(55, 75)
axins1.set_ylim(0.4, 0.6)
axins1.set_xticks(np.arange(55, 85, 10))
axins1.set_yticks(np.linspace(0.4, 0.6, num=3))
axins1.grid(True)
ax1.indicate_inset_zoom(axins1)

# ==============
# 3) Subplot of Tm vs. sticky-end length
#   (based on custom_values mapping)
# ==============
custom_values = {
    "L_-4": 16,
    "L_0": 12,
    "L_1": 11,
    "L_2": 10,
    "L_3": 9,
    "L_4": 8,
    "L_7": 5,
    "L_11": 1,
    "Y": 12,
    "Y+L11": 11,
    "Y+L0": 0
}

ax2 = fig.add_subplot(grid[0, 1])

for i, (label, color) in enumerate(zip(labels, colors)):
    x_value = custom_values.get(label, 0)  # fallback if not in dict
    ax2.errorbar(
        x_value,
        temp_at_05_values[i],
        yerr=temp_at_05_errors[i],
        fmt='o',
        color=color,
        markeredgecolor='k',
        capsize=5
    )
    


ax2.set_ylabel('Melting temperature $(^{o}C)$', fontsize=10)
ax2.set_xlabel("Length of the sticky end (unpaired base)", fontsize=10)
ax2.set_xticks(list(custom_values.values()))
ax2.set_ylim(50, 70)
ax2.text(-0.04, 1.07, 'b)', transform=ax2.transAxes,
              fontsize=12, fontweight='bold', va='top', ha='right')
"""
# ==============
# Create a colorbar for the first set
# ==============
plt.subplots_adjust(wspace=0.4, left=0.1, right=0.95, bottom=0.15, top=0.9)

cleaned_labels = [
    get_custom_label(label)  # use your custom mapping
    for label in labels
]

cmap = mcolors.ListedColormap(colors)
norm_cbar = mcolors.BoundaryNorm(boundaries=np.arange(len(cleaned_labels) + 1) - 0.5,
                                 ncolors=len(cleaned_labels))
sm = cm.ScalarMappable(cmap=cmap, norm=norm_cbar)
sm.set_array([])

cbar_ax = fig.add_axes([0.59, 0.6, 0.38, 0.03])  # [left, bottom, width, height]
cbar = fig.colorbar(sm, cax=cbar_ax, orientation='horizontal',
                    ticks=np.arange(len(cleaned_labels)))
cbar.ax.set_xticklabels(cleaned_labels)
cbar.ax.tick_params(labelsize=7)
cbar.set_label('Molecules')
"""
# ==============
# 4) Experimental Data and Comparison
# ==============
colors_1 = [
    "#313695", "#588cc0", "#a3d3e6", "#e9f6e8", "#fee99d",
    "#fca55d", "#e34933", "#d73027", "#7f0d0b"
][::-1]

#linkers = ['L0', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6','Y']
linkers = ['LS12', 'LS11', 'LS10', 'LS9', 'LS8', 'LS7', 'LS6','Y']

sticky_end_lengths = {
    'LS12': 12,
    'LS11': 11,
    'LS10': 10,
    'LS9': 9,
    'LS8': 8,
    'LS7': 7,
    'LS6': 6,
    'Y': 12
}

df_L0_L1 = pd.read_csv("L0_L1.csv")
df_L2_L3 = pd.read_csv("L2_L3.csv")
df_YL1_L4_L5_L6 = pd.read_csv("YL1_L4_L5_L6.csv")

# Convert to numeric
df_L0_L1 = df_L0_L1.apply(pd.to_numeric, errors='coerce').dropna().reset_index(drop=True)
df_L2_L3 = df_L2_L3.apply(pd.to_numeric, errors='coerce').dropna().reset_index(drop=True)
df_YL1_L4_L5_L6 = df_YL1_L4_L5_L6.apply(pd.to_numeric, errors='coerce').dropna().reset_index(drop=True)

columns_of_interest = {
    'LS12': {'h': (0, 1), 'c': (2, 3)},
    'LS11': {'h': (12, 13), 'c': (14, 15)},
    'LS10': {'h': (0, 1), 'c': (2, 3)},
    'LS9': {'h': (12, 13), 'c': (14, 15)},
    'LS8': {'h': (8, 9), 'c': (10, 11)},
    'LS7': {'h': (20, 21), 'c': (22, 23)},
    'LS6': {'h': (24, 25), 'c': (26, 27)},
    'Y':  {'h': (0, 1), 'c': (2, 3)},
}

fusion_temperatures = []
error_bars = []
ax3 = fig.add_subplot(grid[1, 0])
# ==============
# Plot 2 (bottom-left, ax3)
# ==============
for i, (linker, col_indices) in enumerate(columns_of_interest.items()):
    # (1) Calentamiento
    if linker in ['LS12', 'LS11']:
        df_selected = df_L0_L1[[df_L0_L1.columns[col_indices['h'][0]],
                                df_L0_L1.columns[col_indices['h'][1]]]]
    elif linker in ['LS10', 'LS9']:
        df_selected = df_L2_L3[[df_L2_L3.columns[col_indices['h'][0]],
                                df_L2_L3.columns[col_indices['h'][1]]]]
    else:
        df_selected = df_YL1_L4_L5_L6[[df_YL1_L4_L5_L6.columns[col_indices['h'][0]],
                                       df_YL1_L4_L5_L6.columns[col_indices['h'][1]]]]

    temp_h_col, abs_h_col = df_selected.columns[0], df_selected.columns[1]
    df_h = df_selected.apply(pd.to_numeric, errors='coerce')
    
    # (2) Enfriamiento
    if linker in ['LS12', 'LS11']:
        df_selected_c = df_L0_L1[[df_L0_L1.columns[col_indices['c'][0]],
                                  df_L0_L1.columns[col_indices['c'][1]]]]
    elif linker in ['LS10', 'LS9']:
        df_selected_c = df_L2_L3[[df_L2_L3.columns[col_indices['c'][0]],
                                  df_L2_L3.columns[col_indices['c'][1]]]]
    else:
        df_selected_c = df_YL1_L4_L5_L6[[df_YL1_L4_L5_L6.columns[col_indices['c'][0]],
                                         df_YL1_L4_L5_L6.columns[col_indices['c'][1]]]]

    temp_c_col, abs_c_col = df_selected_c.columns[0], df_selected_c.columns[1]
    df_c = df_selected_c.apply(pd.to_numeric, errors='coerce')
    df_c = df_c.iloc[::-1].reset_index(drop=True)

    # Normalización
    min_h = df_h[abs_h_col].min()
    max_h = df_h[abs_h_col].max()
    min_c = df_c[abs_c_col].min()
    max_c = df_c[abs_c_col].max()

    df_h[f'Absorbance_Norm_{linker}_h'] = 1 - (df_h[abs_h_col] - min_h) / (max_h - min_h)
    df_c[f'Absorbance_Norm_{linker}_c'] = 1 - (df_c[abs_c_col] - min_c) / (max_c - min_c)

    # Promedio entre calentamiento y enfriamiento
    min_len = min(len(df_h[temp_h_col]), len(df_c[temp_c_col]))
    temp_range = df_h[temp_h_col].iloc[:min_len]
    avg_absorbance = (df_h[f'Absorbance_Norm_{linker}_h'].iloc[:min_len] +
                      df_c[f'Absorbance_Norm_{linker}_c'].iloc[:min_len]) / 2

    # Tm => fracción 0.5
    interp_func = interp1d(avg_absorbance, temp_range, bounds_error=False, fill_value='extrapolate')
    temperature_at_05 = interp_func(0.5)
    fusion_temperatures.append((sticky_end_lengths[linker], temperature_at_05))

    # Error = diff en calentamiento vs enfriamiento
    temp_h_05 = interp_func(0.5)
    temp_c_05 = interp1d(df_c[f'Absorbance_Norm_{linker}_c'], df_c[temp_c_col],
                         bounds_error=False, fill_value='extrapolate')(0.5)
    error_bars.append(abs(temp_h_05 - temp_c_05))
    


    # Plot experimental average
    ax3.plot(
        temp_range, avg_absorbance, '-',
        color=colors_1[i],
        label=f'{linker}: $T_{{m}}={temperature_at_05:.2f} ^{{\circ}}C$'
    )
    ax3.fill_between(
        df_h[temp_h_col],
        df_h[f'Absorbance_Norm_{linker}_h'],
        df_c[f'Absorbance_Norm_{linker}_c'],
        color=colors_1[i], alpha=0.2
    )

ax3.legend(loc='lower left')
ax3.set_xlabel('Temperature $(^{o}C)$')
ax3.set_ylabel('Normalized absorbance @ 260 nm')
ax3.text(75, 1, '100 mM NaCl',  
              fontsize=10, color='black', ha='left', va='top', 
              bbox=dict(boxstyle='square,pad=0.1', facecolor='white',
                        edgecolor='none', alpha=0.1))
ax3.text(75, 0.97, r'$\approx 1 \ \mu\mathrm{M}$',  
              fontsize=10, color='black', ha='left', va='top', 
              bbox=dict(boxstyle='square,pad=0.1', facecolor='white',
                        edgecolor='none', alpha=0.1))
ax3.set_xlim(10,100)
ax3.set_xticks(np.arange(10, 110, 10))
ax3.axhline(y=0.5, color='black', alpha=0.4, linestyle='--', linewidth=1)
ax3.text(-0.04, 1.07, 'c)', transform=ax3.transAxes,
              fontsize=12, fontweight='bold', va='top', ha='right')

# Zoom inset on the experimental data
axins3 = ax3.inset_axes([0.72, 0.63, 0.24, 0.24])
for i, (linker, col_indices) in enumerate(columns_of_interest.items()):
    if linker in ['LS12', 'LS11']:
        df_selected = df_L0_L1[[df_L0_L1.columns[col_indices['h'][0]],
                                df_L0_L1.columns[col_indices['h'][1]]]]
    elif linker in ['LS10', 'LS9']:
        df_selected = df_L2_L3[[df_L2_L3.columns[col_indices['h'][0]],
                                df_L2_L3.columns[col_indices['h'][1]]]]
    else:
        df_selected = df_YL1_L4_L5_L6[[df_YL1_L4_L5_L6.columns[col_indices['h'][0]],
                                       df_YL1_L4_L5_L6.columns[col_indices['h'][1]]]]

    temp_h_col, abs_h_col = df_selected.columns[0], df_selected.columns[1]
    df_h = df_selected.apply(pd.to_numeric, errors='coerce')

    if linker in ['LS12', 'LS11']:
        df_selected_c = df_L0_L1[[df_L0_L1.columns[col_indices['c'][0]],
                                  df_L0_L1.columns[col_indices['c'][1]]]]
    elif linker in ['LS10', 'LS9']:
        df_selected_c = df_L2_L3[[df_L2_L3.columns[col_indices['c'][0]],
                                  df_L2_L3.columns[col_indices['c'][1]]]]
    else:
        df_selected_c = df_YL1_L4_L5_L6[[df_YL1_L4_L5_L6.columns[col_indices['c'][0]],
                                         df_YL1_L4_L5_L6.columns[col_indices['c'][1]]]]

    temp_c_col, abs_c_col = df_selected_c.columns[0], df_selected_c.columns[1]
    df_c = df_selected_c.apply(pd.to_numeric, errors='coerce')
    df_c = df_c.iloc[::-1].reset_index(drop=True)

    min_h = df_h[abs_h_col].min()
    max_h = df_h[abs_h_col].max()
    min_c = df_c[abs_c_col].min()
    max_c = df_c[abs_c_col].max()

    df_h[f'Absorbance_Norm_{linker}_h'] = (df_h[abs_h_col] - min_h) / (max_h - min_h)
    df_c[f'Absorbance_Norm_{linker}_c'] = (df_c[abs_c_col] - min_c) / (max_c - min_c)

    min_len = min(len(df_h[temp_h_col]), len(df_c[temp_c_col]))
    temp_range = df_h[temp_h_col].iloc[:min_len]
    avg_absorbance = (df_h[f'Absorbance_Norm_{linker}_h'].iloc[:min_len] +
                      df_c[f'Absorbance_Norm_{linker}_c'].iloc[:min_len]) / 2

    axins3.plot(temp_range, avg_absorbance, '-', color=colors_1[i], alpha=0.7)

#ax3.indicate_inset_zoom(axins)
axins3.set_xlim(50, 68)
axins3.set_xticks(np.arange(50, 72, 8))
axins3.set_ylim(0.4, 0.6)
axins3.set_yticks(np.arange(0.4, 0.7, 0.1))
axins3.grid(True)



ax4 = fig.add_subplot(grid[1, 1])

# ==============
# 5) Tm vs Sticky end length (Experimental) in bottom-right subplot

# ==============
for i, (length, temp) in enumerate(fusion_temperatures):
    ax4.errorbar(
        length, temp, yerr=error_bars[i],
        fmt='o', color=colors_1[i],
        markeredgecolor='black', capsize=5
    )

ax4.set_xlabel('Sticky end length (nucleotide)')
ax4.set_ylabel('Melting temperature $(^{o}C)$')
ax4.set_ylim(50,70)
ax4.set_xlim(-1,17)
ax4.text(-0.04, 1.07, 'd)', transform=ax4.transAxes,
              fontsize=12, fontweight='bold', va='top', ha='right')
"""
# Extra colorbar for the second set
cmap_1 = mcolors.ListedColormap(colors_1)
norm_1 = mcolors.BoundaryNorm(boundaries=np.arange(len(linkers) + 1) - 0.5,
                              ncolors=len(linkers))
sm_1 = cm.ScalarMappable(cmap=cmap_1, norm=norm_1)
sm_1.set_array([])

cbar_ax_1 = fig.add_axes([0.60, 0.214, 0.3, 0.03])
cbar_1 = plt.colorbar(sm_1, cax=cbar_ax_1, orientation='horizontal',
                      ticks=np.arange(len(linkers)))

cbar_1.ax.set_xticklabels(linkers, fontsize=9)  # Change 10 to desired size
# Set the size of the colorbar label
cbar_1.set_label('Molecules', fontsize=10)
"""


ax5 = fig.add_subplot(grid[2, 0])

def shuffled_sterr(_in):
    _chunk = 20
    N = _in.size // _chunk
    if N == 0:
        return 0
    out = 0
    for i in range(N):
        _ids = np.random.randint(0, high=_in.size, size=_chunk)
        out += np.std(_in[_ids])
    return np.sqrt(out / N) / np.sqrt(N)

def extrapolate_bulk_array(x):
    temp = x / (1 - x)
    return (1 + 1 / (2 * temp)) - np.sqrt((1 + 1. / (2 * temp)) ** 2 - 1)

def find_temp_for_fraction(temps, fractions, target_fraction=0.5):
    """Interpolar linealmente para encontrar la temperatura en la cual la fracción alcanza un valor dado."""
    f = interp1d(fractions, temps, kind='linear', fill_value="extrapolate")
    return f(target_fraction)
"""
# Configurar estilo y fuentes
plt.rcParams.update({
    'font.size': 10,
    'font.family': 'serif',
    'axes.labelsize': 10,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 8,
})

# Tamaño de los marcadores
marker_size = 10
"""
# Colores especificados
colors = [
    "#313695", "#3c59a6", "#588cc0", "#74add1", "#a3d3e6", "#c1e4ef",
    "#e9f6e8", "#fbfdc7", "#fee99d", "#fed283", "#fca55d", "#f67f4b",
    "#e34933", "#ce2827", "#d73027", "#a50f15", "#7f0d0b"
][::-1]



def find_melting_temperature(df):
    # Encontrar el índice de la temperatura más cercana a la fracción 0.5
    idx1 = np.abs(df.iloc[:, 1] - 0.5).idxmin()
    
    # Encontrar los índices de los puntos más cercanos a ambos lados del punto 0.5
    idx2 = idx1 + 1 if df.iloc[idx1, 1] < 0.5 else idx1 - 1
    
    # Realizar interpolación lineal entre los puntos más cercanos
    x1, y1 = df.iloc[idx1, 0], df.iloc[idx1, 1]
    x2, y2 = df.iloc[idx2, 0], df.iloc[idx2, 1]
    tm = x1 + ((0.5 - y1) * (x2 - x1) / (y2 - y1))
    
    return tm

# Directorios donde se encuentran los archivos CSV
directory_bases_profiles_1 = "./0.1M"
directory_nupack_1 = "./0.1M"

# Plot de los perfiles de desaparición de bases
archivos_bases_profiles_1 = [archivo for archivo in os.listdir(directory_bases_profiles_1) if archivo.endswith(".csv")]

# Lista para almacenar los nombres de los archivos y los valores de Tm
tm_info_1 = []

# Iterar sobre los archivos CSV y encontrar Tm para cada uno
for archivo in archivos_bases_profiles_1:
    df = pd.read_csv(os.path.join(directory_bases_profiles_1, archivo))
    tm = find_melting_temperature(df)
    tm_info_1.append((archivo, tm))

# Imprimir los nombres de los archivos y los valores de Tm
print("Información de Tm:")
for nombre_archivo, tm in tm_info_1:
    print(f"Archivo: {nombre_archivo}, Tm: {tm:.2f} °C")

# Crear un DataFrame con la información de Tm
    df_tm_info_1 = pd.DataFrame(tm_info_1, columns=['molecule', 'Tm'])

# Separar y procesar el nombre de las moléculas
    df_tm_info_1[['molecule_clean', 'type_of_file']] = df_tm_info_1['molecule'].str.split('.c', expand=True)
    df_tm_info_1.drop(columns=['molecule', 'type_of_file'], inplace=True)

# Procesar los archivos que comienzan con "L"
    df_tm_info_1[['type_of_mol', 'molecule_clean_id']] = df_tm_info_1['molecule_clean'].str.extract(r'([A-Za-z]+)([0-9\-]+)')

# Asignar valores personalizados para "Y0", "Y1" y "Y12"
    custom_nucleotide_values = {'Y0': 12, 'Y1': 11, 'Y12': 0}
    df_tm_info_1['nucleotides'] = df_tm_info_1['molecule_clean'].map(custom_nucleotide_values)

# Para los "L", asignar correctamente los sticky ends: L0 → 12, L1 → 11, ..., L12 → 0
    mask_l = df_tm_info_1['type_of_mol'] == 'L'
    df_tm_info_1.loc[mask_l, 'molecule_clean_id'] = df_tm_info_1.loc[mask_l, 'molecule_clean_id'].astype(int)
    df_tm_info_1.loc[mask_l, 'nucleotides'] = 12 - df_tm_info_1.loc[mask_l, 'molecule_clean_id']

# Resetear índice y mostrar el DataFrame final
    df_tm_info_1 = df_tm_info_1.sort_values(by='nucleotides', ascending=False).reset_index(drop=True)

    
    # Eliminar las columnas intermedias que ya no necesitamos
    #df_tm_info_1.drop(columns=['Concentracion_ext', 'ext'], inplace=True)

    print(df_tm_info_1)
# Graficar cada perfil de desaparición de bases por separado
# Graficar la primera curva de desaparición de bases
# Graficar la primera curva de desaparición de bases

colors_1 = [
"#313695", "#3c59a6", "#588cc0", "#74add1", "#a3d3e6", "#c1e4ef", "#e9f6e8", "#fbfdc7", "#fee99d", "#fed283", "#fca55d", "#f67f4b", "#e34933", "#ce2827", "#d73027", "#a50f15", "#7f0d0b"
]

custom_names = {
    "Y0": "Y",
    "L-4": "LS16",
    "L0": "LS12",
    "L1": "LS11",
    "L2": "LS10",
    "L3": "LS9",
    "L4": "LS8",
    "L5": "LS7",
    "L6": "LS6",
    "L7": "LS5", 
    "L8":"LS4",
    "L9":"LS3",
    "L10":"LS2",
    "L11":"LS1",
    "L12":"LS0",
}

name_original = df_tm_info_1['molecule_clean'].iloc[4]
name_plot = custom_names.get(name_original, name_original)

#ax5 = fig.add_subplot(grid[2, 0])

df1 = pd.read_csv(os.path.join(directory_bases_profiles_1, archivos_bases_profiles_1[3]))
ax5.plot(df1["Temperature(C)"], 1 - df1["Fraction of bases unpaired at equilibrium"], color=colors_1[0], label=custom_names["Y0"],linestyle="--")

df4 = pd.read_csv(os.path.join(directory_bases_profiles_1, archivos_bases_profiles_1[0]))
ax5.plot(df4["Temperature(C)"], 1 - df4["Fraction of bases unpaired at equilibrium"], color=colors_1[1], label=custom_names["L-4"])

df5 = pd.read_csv(os.path.join(directory_bases_profiles_1, archivos_bases_profiles_1[13]))
ax5.plot(df5["Temperature(C)"], 1 - df5["Fraction of bases unpaired at equilibrium"], color=colors_1[2], label=custom_names["L0"])

df6 = pd.read_csv(os.path.join(directory_bases_profiles_1, archivos_bases_profiles_1[16]))
ax5.plot(df6["Temperature(C)"], 1 - df6["Fraction of bases unpaired at equilibrium"], color=colors_1[3], label=custom_names["L1"])

df7 = pd.read_csv(os.path.join(directory_bases_profiles_1, archivos_bases_profiles_1[12]))
ax5.plot(df7["Temperature(C)"], 1 - df7["Fraction of bases unpaired at equilibrium"], color=colors_1[4], label=custom_names["L2"])

df8 = pd.read_csv(os.path.join(directory_bases_profiles_1, archivos_bases_profiles_1[11]))
ax5.plot(df8["Temperature(C)"], 1 - df8["Fraction of bases unpaired at equilibrium"], color=colors_1[5], label=custom_names["L3"])

df9 = pd.read_csv(os.path.join(directory_bases_profiles_1, archivos_bases_profiles_1[8]))
ax5.plot(df9["Temperature(C)"], 1 - df9["Fraction of bases unpaired at equilibrium"], color=colors_1[6], label=custom_names["L4"])

df10 = pd.read_csv(os.path.join(directory_bases_profiles_1, archivos_bases_profiles_1[7]))
ax5.plot(df10["Temperature(C)"], 1 - df10["Fraction of bases unpaired at equilibrium"], color=colors_1[7], label=custom_names["L5"])

df11 = pd.read_csv(os.path.join(directory_bases_profiles_1, archivos_bases_profiles_1[9]))
ax5.plot(df11["Temperature(C)"], 1 - df11["Fraction of bases unpaired at equilibrium"], color=colors_1[8], label=custom_names["L6"])

df12 = pd.read_csv(os.path.join(directory_bases_profiles_1, archivos_bases_profiles_1[10]))
ax5.plot(df12["Temperature(C)"], 1 - df12["Fraction of bases unpaired at equilibrium"], color=colors_1[9], label=custom_names["L7"])

df13 = pd.read_csv(os.path.join(directory_bases_profiles_1, archivos_bases_profiles_1[2]))
ax5.plot(df13["Temperature(C)"], 1 - df13["Fraction of bases unpaired at equilibrium"], color=colors_1[10], label=custom_names["L8"])

df14 = pd.read_csv(os.path.join(directory_bases_profiles_1, archivos_bases_profiles_1[1]))
ax5.plot(df14["Temperature(C)"], 1 - df14["Fraction of bases unpaired at equilibrium"], color=colors_1[11], label=custom_names["L9"])

df15 = pd.read_csv(os.path.join(directory_bases_profiles_1, archivos_bases_profiles_1[5]))
ax5.plot(df15["Temperature(C)"], 1 - df15["Fraction of bases unpaired at equilibrium"], color=colors_1[12], label=custom_names["L10"])

df16 = pd.read_csv(os.path.join(directory_bases_profiles_1, archivos_bases_profiles_1[6]))
ax5.plot(df16["Temperature(C)"], 1 - df16["Fraction of bases unpaired at equilibrium"], color=colors_1[13], label=custom_names["L11"])

df17 = pd.read_csv(os.path.join(directory_bases_profiles_1, archivos_bases_profiles_1[4]))
ax5.plot(df17["Temperature(C)"], 1 - df17["Fraction of bases unpaired at equilibrium"], color=colors_1[14], label=custom_names["L12"])



# Añadir leyenda con dos columnas
ax5.legend(loc='lower left', fontsize=10, ncol=2)

ax5.text(-0.04, 1.07, 'a)', transform=ax5.transAxes, fontsize=12, fontweight='bold', va='top', ha='right')

# Configuración del zoom
axins = ax5.inset_axes([0.7, 0.58, 0.25, 0.25])  # (x0, y0, width, height) del panel secundario
# Aplicar los colores a los plots
axins.plot(df1["Temperature(C)"], df1["Fraction of bases unpaired at equilibrium"], color=colors_1[0],linestyle="--")
#axins.plot(df2["Temperature(C)"], df2["Fraction of bases unpaired at equilibrium"], color=colors_1[1])
#axins.plot(df3["Temperature(C)"], df3["Fraction of bases unpaired at equilibrium"], color=colors_1[2])
axins.plot(df4["Temperature(C)"], df4["Fraction of bases unpaired at equilibrium"], color=colors_1[1])
axins.plot(df5["Temperature(C)"], df5["Fraction of bases unpaired at equilibrium"], color=colors_1[2])
axins.plot(df6["Temperature(C)"], df6["Fraction of bases unpaired at equilibrium"], color=colors_1[3])
axins.plot(df7["Temperature(C)"], df7["Fraction of bases unpaired at equilibrium"], color=colors_1[4])
axins.plot(df8["Temperature(C)"], df8["Fraction of bases unpaired at equilibrium"], color=colors_1[5])
axins.plot(df9["Temperature(C)"], df9["Fraction of bases unpaired at equilibrium"], color=colors_1[6])
axins.plot(df10["Temperature(C)"], df10["Fraction of bases unpaired at equilibrium"], color=colors_1[7])
axins.plot(df11["Temperature(C)"], df11["Fraction of bases unpaired at equilibrium"], color=colors_1[8])
axins.plot(df12["Temperature(C)"], df12["Fraction of bases unpaired at equilibrium"], color=colors_1[9])
axins.plot(df13["Temperature(C)"], df13["Fraction of bases unpaired at equilibrium"], color=colors_1[10])
axins.plot(df14["Temperature(C)"], df14["Fraction of bases unpaired at equilibrium"], color=colors_1[11])
axins.plot(df15["Temperature(C)"], df15["Fraction of bases unpaired at equilibrium"], color=colors_1[12])
axins.plot(df16["Temperature(C)"], df16["Fraction of bases unpaired at equilibrium"], color=colors_1[13])
axins.plot(df17["Temperature(C)"], df17["Fraction of bases unpaired at equilibrium"], color=colors_1[14])


axins.set_xlim(52, 65)

axins.set_ylim(0.46, 0.54)
axins.set_xticks(np.arange(52, 73, 8))
axins.set_yticks(np.arange(0.46, 0.56, 0.04))
ax5.set_ylabel("$P\,(\mathrm{bound \,state})$", fontsize=12)
axins.grid(True)

# Cambiar el tamaño de la fuente de los ticks del gráfico de zoom
axins.tick_params(axis='both', which='major', labelsize=10)
# Configuración del tamaño de fuente de los ticks del primer subplot
ax5.tick_params(axis='both', which='major', labelsize=10)  # Aquí ajusta el tamaño de fuente según tus preferencias

# Añadir guías al gráfico principal para mostrar el área de zoom
ax5.indicate_inset_zoom(axins, edgecolor="black")
#ax.plot(df_tm_info_1['con_wo_units'], df_tm_info_1['Tm'], marker='o', linestyle='-')
ax5.set_ylim(0, 1)
ax5.set_yticks(np.arange(0, 1.2, 0.2))
ax5.set_xlabel('Temperature ($^{o}$C)', fontsize=12)
ax5.axhline(y=0.5, color='black', alpha=0.4, linestyle='--', linewidth=1, xmin=0.47, xmax=1)
#ax5.tick_params(axis='y', which='both', left=True, right=True, labelleft=False, labelright=False)

#ax5.set_ylabel('Fraction unbounded', fontsize=12)
ax5.text(80, 0.97, '100 mM NaCl',  
               fontsize=10, color='black',  
               ha='left', va='top', 
               bbox=dict(boxstyle='square,pad=0.1', facecolor='white', edgecolor='none', alpha=0.1))

# Anotación en la segunda línea
ax5.text(80, 0.93, '1 $\mu$M',  
               fontsize=10, color='black',  
               ha='left', va='top',
               bbox=dict(boxstyle='square,pad=0.1', facecolor='white', edgecolor='none', alpha=0.1))


# Ordenar el DataFrame por la columna Tm
df_tm_info_1 = df_tm_info_1.sort_values(by='Tm', ascending=True).reset_index(drop=True)

# Función para excluir filas específicas por índice o valores
def excluir_filas(df, excluir_indices=None, excluir_valores=None, columna=None):
    if excluir_indices:
        df = df.drop(excluir_indices, axis=0).reset_index(drop=True)
    if excluir_valores and columna:
        df = df[~df[columna].isin(excluir_valores)].reset_index(drop=True)
    return df

# Definir los índices o valores a excluir
indices_a_excluir = [0, 1]  # Índices a excluir
#valores_a_excluir = [68.5]  # Valores a excluir en la columna Tm

# Aplicar exclusión de filas
df_tm_info_1 = excluir_filas(df_tm_info_1, excluir_indices=indices_a_excluir)
#df_tm_info_1 = excluir_filas(df_tm_info_1, excluir_valores=valores_a_excluir, columna='Tm')

# Definir colores (gradiente inverso)
colors = [
    "#313695", "#3c59a6", "#588cc0", "#74add1", "#a3d3e6", "#c1e4ef",
    "#e9f6e8", "#fbfdc7", "#fee99d", "#fed283", "#fca55d", "#f67f4b",
    "#e34933", "#ce2827", "#d73027", "#a50f15", "#7f0d0b"
][::1]

# Asignar colores a los puntos según el índice de orden
assigned_colors = [colors[i % len(colors)] for i in range(len(df_tm_info_1))]

ax6 = fig.add_subplot(grid[2, 1])

for i, (length, temp) in enumerate(fusion_temperatures):
    ax6.errorbar(
        length, temp, yerr=error_bars[i],
        fmt='o', color=colors_1[i],
        markeredgecolor='black', capsize=5
    )

#ax2 = fig.add_subplot(grid[0, 1])

for i, (label, color) in enumerate(zip(labels, colors)):
    x_value = custom_values.get(label, 0)  # fallback if not in dict
    ax6.errorbar(
        x_value,
        temp_at_05_values[i],
        yerr=temp_at_05_errors[i],
        fmt='o',
        color=color,
        markeredgecolor='k',
        capsize=5
    )
    


#ax2.set_ylabel('Melting temperature $(^{o}C)$', fontsize=10)
#ax2.set_xlabel("Length of the sticky end (unpaired base)", fontsize=10)
#ax2.set_xticks(list(custom_values.values()))
#ax2.set_ylim(50, 70)
#ax2.text(-0.04, 1.07, 'b)', transform=ax2.transAxes,
              #fontsize=12, fontweight='bold', va='top', ha='right')

# ==============
# 4) Experimental Data and Comparison
# ==============
colors_1 = [
    "#313695", "#588cc0", "#a3d3e6", "#e9f6e8", "#fee99d",
    "#fca55d", "#e34933", "#d73027", "#7f0d0b"
]

#linkers = ['L0', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6','Y']
linkers = ['LS12', 'LS11', 'LS10', 'LS9', 'LS8', 'LS7', 'LS6','Y']

sticky_end_lengths = {
    'LS12': 12,
    'LS11': 11,
    'LS10': 10,
    'LS9': 9,
    'LS8': 8,
    'LS7': 7,
    'LS6': 6,
    'Y': 12
}

df_L0_L1 = pd.read_csv("L0_L1.csv")
df_L2_L3 = pd.read_csv("L2_L3.csv")
df_YL1_L4_L5_L6 = pd.read_csv("YL1_L4_L5_L6.csv")

# Convert to numeric
df_L0_L1 = df_L0_L1.apply(pd.to_numeric, errors='coerce').dropna().reset_index(drop=True)
df_L2_L3 = df_L2_L3.apply(pd.to_numeric, errors='coerce').dropna().reset_index(drop=True)
df_YL1_L4_L5_L6 = df_YL1_L4_L5_L6.apply(pd.to_numeric, errors='coerce').dropna().reset_index(drop=True)

columns_of_interest = {
    'LS12': {'h': (0, 1), 'c': (2, 3)},
    'LS11': {'h': (12, 13), 'c': (14, 15)},
    'LS10': {'h': (0, 1), 'c': (2, 3)},
    'LS9': {'h': (12, 13), 'c': (14, 15)},
    'LS8': {'h': (8, 9), 'c': (10, 11)},
    'LS7': {'h': (20, 21), 'c': (22, 23)},
    'LS6': {'h': (24, 25), 'c': (26, 27)},
    'Y':  {'h': (0, 1), 'c': (2, 3)},
}

fusion_temperatures = []
error_bars = []


#aqui acaba los datos MD

#ax4 = fig.add_subplot(grid[1, 1])

# ==============
# 5) Tm vs Sticky end length (Experimental) in bottom-right subplot

# ==============


#ax4.set_xlabel('Sticky end length (nucleotide)')
#ax4.set_ylabel('Melting temperature $(^{o}C)$')
#ax4.set_ylim(50,70)
#ax4.set_xlim(-1,17)
#ax4.text(-0.04, 1.07, 'd)', transform=ax4.transAxes,
              #fontsize=12, fontweight='bold', va='top', ha='right')
"""
# Extra colorbar for the second set
cmap_1 = mcolors.ListedColormap(colors_1)
norm_1 = mcolors.BoundaryNorm(boundaries=np.arange(len(linkers) + 1) - 0.5,
                              ncolors=len(linkers))
sm_1 = cm.ScalarMappable(cmap=cmap_1, norm=norm_1)
sm_1.set_array([])

cbar_ax_1 = fig.add_axes([0.60, 0.214, 0.3, 0.03])
cbar_1 = plt.colorbar(sm_1, cax=cbar_ax_1, orientation='horizontal',
                      ticks=np.arange(len(linkers)))

cbar_1.ax.set_xticklabels(linkers, fontsize=9)  # Change 10 to desired size
# Set the size of the colorbar label
cbar_1.set_label('Molecules', fontsize=10)
"""


# Asignar colores a los puntos según el índice de orden
assigned_colors = [colors[i % len(colors)] for i in range(len(df_tm_info_1))]

# Primer punto con triángulo
ax6.scatter(
    df_tm_info_1['nucleotides'].iloc[0], df_tm_info_1['Tm'].iloc[0],
    marker='^', color=assigned_colors[0], s=30, edgecolor='black'
)

# Resto de los puntos con círculos
ax6.scatter(
    df_tm_info_1['nucleotides'].iloc[1:], df_tm_info_1['Tm'].iloc[1:],
    marker='o', color=assigned_colors[1:], s=30, edgecolor='black'
)

# Configurar los ejes y apariencia
ax6.set_xlabel('$n_{\mathrm{tail}}$', fontsize=12)
ax6.set_xticks(np.arange(0, 17, 1))  # Especificar los ticks del eje x
ax6.tick_params(axis='both', which='major', labelsize=10)
ax6.set_ylim(50, 70)
ax6.set_ylabel('Melting temperature $T_{m}$ $(^{o}C)$', fontsize=12)
ax6.text(-0.04, 1.07, 'b)', transform=ax6.transAxes, fontsize=12, fontweight='bold', va='top', ha='right')
labels = list(custom_names.values())





"""
# Ruta base
img_dir = '/Users/alejandrosoto/Documents/GitHub/numerical_study_tm_KIT/paperJCTC_Soto/Plots/Plots/Figure 2/'

# Imagen 1: linker_1.png con rotación
img1 = mpimg.imread(img_dir + 'linker_1.png')
img1_rotated = rotate(img1, angle=15, reshape=True)  # Rota 15 grados
imagebox1 = OffsetImage(img1_rotated, zoom=0.20)
ab1 = AnnotationBbox(imagebox1, (11,68), frameon=False, box_alignment=(0.5, 0.5), xycoords='data')
ax6.add_artist(ab1)

# Imagen 2: yshape_1.png sin rotación
img2 = mpimg.imread(img_dir + 'yshape_1.png')
imagebox2 = OffsetImage(img2, zoom=0.3)
ab2 = AnnotationBbox(imagebox2, (5.5,56), frameon=False, box_alignment=(0.5, 0.5), xycoords='data')
ax6.add_artist(ab2)
"""
#plt.tight_layout()
plt.savefig('figure_4.png', dpi=400)
plt.show()
