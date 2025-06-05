import pandas as pd
import matplotlib.pyplot as plt
import scienceplots
import numpy as np

plt.style.use(['science', 'no-latex', 'bright'])

# Datos de ejemplo
data = {
    'nucleotides': [0, 1, 2, 3, 4, 5, 6],
    'DM':           [-2.62, 0, -0.38, -1.36, -1.93, -2.24, -2.44],
    'oxDNA':        [-2.5,  0, -0.5,  -1.5,  -1,    -1.5,  -2   ],
    'US':           [-3.4,  0, -2.7,  None,  -3.8,  -4.3,  -4.8 ],
    'NUPACK':       [-0.46004402, 0, -2.035865541, -4.222682693, -6.494932418, -9.251695084, -12.67174151]
}
df = pd.DataFrame(data)

# Creamos la figura y el eje con estilo de revista
fig, ax = plt.subplots(figsize=(6, 6))

marker_size = 9  # Ajusta este valor según prefieras

# Trazamos cada serie en negro (línea y marcador)
ax.plot(
    df['nucleotides'],
    df['DM'],
    color='black',
    marker='D',
    markerfacecolor='black',
    markeredgecolor='black',
    linestyle='-',
    markersize=marker_size,
    label=r'$\mathit{Di\ Michele\ et.\ al}$'
)

ax.plot(
    df['nucleotides'],
    df['oxDNA'],
    color='red',
    marker='o',
    markerfacecolor='red',
    markeredgecolor='red',
    linestyle='-',
    markersize=marker_size,
    label='oxDNA unbiased MD'
)

mask_us = df['US'].notna()
ax.plot(
    df.loc[mask_us, 'nucleotides'],
    df.loc[mask_us, 'US'],
    color='black',
    marker='p',
    markerfacecolor='black',
    markeredgecolor='black',
    linestyle='-',
    markersize=marker_size,
    label='oxDNA Umbrella Sampling'
)

ax.plot(
    df['nucleotides'],
    df['NUPACK'],
    color='black',
    marker='^',
    markerfacecolor='black',
    markeredgecolor='black',
    linestyle='-',
    markersize=marker_size,
    label='NUPACK'
)

# Etiquetas de ejes y título (sin sintaxis LaTeX)
ax.set_xlabel('$n_{tail}$', fontsize=14, color='black')
ax.set_ylabel('$\delta T_{m}$ ($^{o}$C)', fontsize=14, color='black')
ax.set_xticks(np.arange(0, 7, 1))
ax.tick_params(axis='both', which='major', labelsize=14, colors='black')
ax.tick_params(axis='both', which='minor', colors='black')

# Leyenda sin recuadro y ubicada abajo a la izquierda (texto negro)
ax.legend(frameon=False, loc='lower left', fontsize=12)

# Ticks hacia adentro, con ticks menores activados
ax.tick_params(axis='both', which='major', direction='in', length=6)
ax.tick_params(axis='both', which='minor', direction='in', length=3)
ax.minorticks_on()

# Ajustar márgenes interiores
plt.tight_layout()
plt.savefig('figure_val.png', dpi=400)
# Mostrar el gráfico
plt.show()
