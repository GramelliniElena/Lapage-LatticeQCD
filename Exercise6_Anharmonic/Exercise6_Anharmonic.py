#Come è stato trovato il tratteggio 1.933?
import numpy as np
from numba import njit
import matplotlib.pyplot as plt

N = 20  #number of steps
a = 0.25  #lattice spacing
eps = 1.4  #Metropolis' step
N_cor = 40  #decorrelation steps
N_cf_list = [25, 100, 1000, 10000]  #number of configurations
c = 2  #anharmonic potential coefficient

@njit
def anharmonic_potential(x):
    return (x**2 / 2) * (1 + c * x**2) + (a**2 / 24) * ((x + 2 * c * x**3)**2) - ((a * c * x**2) / 4) + ((a**3) / 2) * (c * x**2 / 4)**2

#Definition of the action without ghost modes
@njit
def action(j, x):
    jp = (j + 1) % N
    jm = (j - 1) % N

    kinetic = x[j] * (x[j] - x[jp] - x[jm]) / a

    potential = a * anharmonic_potential(x[j]) *(1 + a**2 / 12)
    return kinetic + potential

#Metropolis update of the path
@njit
def update(x):
    for j in range(N):
        old_x = x[j]
        old_action = action(j, x)
        x[j] += np.random.uniform(-eps, eps)  #add a random number in (-eps, eps) to x[j]
        d_action = action(j, x) - old_action
        if d_action > 0 and np.exp(-d_action) < np.random.uniform(0, 1):  #condition to keep the old x[j]
            x[j] = old_x

@njit
def compute_G(x, n):
    g = 0
    for j in range(N):
        g += x[j] * x[(j + n) % N]
    return g / N

def MCaverage(x, G, N_config): 
    for j in range(N): 
        x[j] = 0  #initialization
    for _ in range(5 * N_cor):  #thermalization
        update(x)
    for alpha in range(N_config):
        for _ in range(N_cor):  #erasing correlations
            update(x)
        for n in range(N):
            G[alpha][n] = compute_G(x, n)
    return G

def bootstrap(G, nbstrap=100):
    N_cf, N = G.shape
    deltaEs = np.zeros((nbstrap, N - 1))   #initialization of nbstrap rows and N-1 columns

    for b in range(nbstrap):
        idx = np.random.choice(N_cf, N_cf, replace=True)  #resampling
        sample = G[idx]
        avg = np.mean(sample, axis=0)
        for t in range(N - 1):
            if avg[t] > 0 and avg[t + 1] > 0:
                deltaEs[b, t] = np.log(avg[t] / avg[t + 1]) / a  #condition to have a valid logarithm
            else:
                deltaEs[b, t] = 0.0

    avgE = np.mean(deltaEs, axis=0)
    sdevE = np.std(deltaEs, axis=0, ddof=1)
    return avgE, sdevE

def analysis(N_cf_list):
    for N_config in N_cf_list:
        print(f"Simulazione con N_config = {N_config}")
        G = np.zeros((N_config, N))
        x = np.zeros(N)
        G = MCaverage(x, G, N_config)
        avgE, sdevE = bootstrap(G, nbstrap=100)

        t_vals = [a * q for q in range(len(avgE))]

        plt.figure(figsize=(8, 5), dpi=120)
        plt.errorbar(t_vals, avgE, yerr=sdevE, fmt='o', ecolor="blue", elinewidth=0.8,
                     capsize=4, markersize=4, label='Calculated ΔE')
        plt.scatter(t_vals, avgE, color='black', s=8)

        plt.axhline(y = 1.933, color='black', linestyle='--')

        plt.title(f'NO-Ghost Action (with anharmonic potential): N={N}, a={a}, N_config={N_config}')        
        plt.xlabel('t')
        plt.ylabel(r'$\Delta E(t)$')
        plt.xlim(0, 1.4)
        plt.ylim(-3, 3)
        plt.grid(True)
        plt.legend()
        plt.tight_layout()

        filename =f"E6_Anharmonic_Ghost_Nconfig{N_config}.png"
        plt.savefig(filename)
        plt.close()
        print(f"Plot saved as '{filename}'\n")
if __name__ == "__main__":
    analysis(N_cf_list)
