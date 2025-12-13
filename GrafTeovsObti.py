import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Load the datasets
# Assuming standard CSV format based on the user's upload
try:
    df_theo = pd.read_csv('CSF_Lab.xlsx - Teórico.csv')
    df_obt = pd.read_csv('CSF_Lab.xlsx - Obtido.csv')
except:
    # Fallback if parsing fails initially (sometimes CSVs from Excel have skipping rows)
    df_theo = pd.read_csv('CSF_Lab.xlsx - Teórico.csv', header=1) # Adjusting header row if needed
    df_obt = pd.read_csv('CSF_Lab.xlsx - Obtido.csv', header=2)   # Adjusting header row if needed

# Cleaning and Identifying columns (Heuristic approach based on file names)
# We look for numeric columns. Usually: Col 0 = Distance, Col 1 = Freq1, Col 2 = Freq2

def clean_data(df):
    # Select only numeric columns roughly
    df = df.select_dtypes(include=[np.number])
    # Drop rows with NaN
    df = df.dropna()
    return df

df_theo_clean = clean_data(df_theo)
df_obt_clean = clean_data(df_obt)

# Assigning axes based on typical structure: Dist, 800MHz, 3.4GHz
# Using .iloc to access by position regardless of column names
t_dist = df_theo_clean.iloc[:, 0]
t_800 = df_theo_clean.iloc[:, 1]
t_3400 = df_theo_clean.iloc[:, 2]

o_dist = df_obt_clean.iloc[:, 0]
o_800 = df_obt_clean.iloc[:, 1]
o_3400 = df_obt_clean.iloc[:, 2]

# Plotting
plt.figure(figsize=(10, 6))

# Theoretical Lines
plt.plot(t_dist, t_800, label='Teórico 800 MHz (LFSPL)', color='blue', linestyle='--', linewidth=2)
plt.plot(t_dist, t_3400, label='Teórico 3.4 GHz (LFSPL)', color='red', linestyle='--', linewidth=2)

# Obtained Points
plt.scatter(o_dist, o_800, label='Medido (Obtido) 800 MHz', color='blue', marker='o', s=50, alpha=0.7)
plt.scatter(o_dist, o_3400, label='Medido (Obtido) 3.4 GHz', color='red', marker='x', s=50, alpha=0.7)

# Formatting
plt.title('Comparação: Perda de Caminho Teórica vs. Medida', fontsize=14)
plt.xlabel('Distância (metros)', fontsize=12)
plt.ylabel('Atenuação / Perda de Caminho (dB)', fontsize=12)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()

# Show the plot
plt.tight_layout()
plt.show()