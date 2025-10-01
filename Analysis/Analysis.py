import numpy as np
import matplotlib.pyplot as plt

# Fixed configuration for the diffent action's choices
N = 20
a = 0.5
t_vals = np.arange(N - 1) * a


def load_data(file):
    #[avgE, sdevE] was saved as np.savez
    data = np.load(file)
    return data["avgE"], data["sdevE"]

file_standard = "d_E4_standard.npz"
file_improved = "d_E5_improved.npz"
file_noghost  = "d_E6_noghost.npz"

#Load data
avgE_std, err_std = load_data(file_standard)
avgE_imp, err_imp = load_data(file_improved)
avgE_ngh, err_ngh = load_data(file_noghost)


plt.figure(figsize=(8, 5), dpi=120)

plt.errorbar(t_vals, avgE_std, yerr=err_std, fmt='o', label='Standard action', color='blue', capsize=3)
plt.errorbar(t_vals, avgE_imp, yerr=err_imp, fmt='s', label='Improved action', color='green', capsize=3)
plt.errorbar(t_vals, avgE_ngh, yerr=err_ngh, fmt='^', label='No-ghost action', color='red', capsize=3)

plt.axhline(y=1.0, color='black', linestyle='--', linewidth=1.0, label='Exact $\\omega = 1$')
plt.title('Compared different actions (harmonic oscillator)')
plt.xlabel('t')
plt.ylabel(r'$\Delta E(t)$')
plt.xlim(0, 3.2)
plt.ylim(0, 2)
plt.grid(True)
plt.legend()
plt.tight_layout()

plt.savefig("Comparison_differentActions.png")
plt.close()
print("Plot saved as 'Comparison_differentActions'")
