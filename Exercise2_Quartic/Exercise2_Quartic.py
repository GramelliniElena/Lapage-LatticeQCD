import numpy as np
from numba import njit
import matplotlib.pyplot as plt

N = 20  #number of steps
a = 0.5  #lattice spacing
eps = 1.4  #Metropolis' step
N_cor = 20  #decorrelation steps
N_cf_list = [25, 100, 1000, 10000]  #number of configurations


@njit
def quartic_potential(x):
    return 0.5 * x**4

@njit
def action(j, x):
    jp = (j + 1) % N
    jm = (j - 1) % N
    return a * quartic_potential(x[j]) + x[j] * (x[j] - x[jp] - x[jm]) / a

#Metropolis update of the path: randomizing x[j] at the jth site
@njit
def update(x):
    for j in range(N):
        old_x = x[j]
        old_action = action(j, x)
        x[j] += np.random.uniform(-eps, eps)  #add a random number in (-eps, eps) to x[j]
        d_action = action(j, x) - old_action
        if d_action > 0 and np.exp(-d_action) < np.random.uniform(0, 1):  #conditions to keep the old x[j]
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
    deltaEs = np.zeros((nbstrap, N - 1))  #initialization of nbstrap rows and N-1 columns

    for b in range(nbstrap):
        idx = np.random.choice(N_cf, N_cf, replace=True)  #resampling
        sample = G[idx]
        avg = np.mean(sample, axis=0)
        for t in range(N - 1):
            if avg[t] > 0 and avg[t + 1] > 0:  #condition to have a valid logarithm
                deltaEs[b, t] = np.log(avg[t] / avg[t + 1]) / a
            else:
                deltaEs[b, t] = 0.0

    avgE = np.mean(deltaEs, axis=0)
    sdevE = np.std(deltaEs, axis=0, ddof=1)
    return avgE, sdevE

def analysis(N_cf_list):
    for N_config in N_cf_list:
        print(f"Simulazione con N_config = {N_config}...")
        G = np.zeros((N_config, N))
        x = np.zeros(N)
        G = MCaverage(x, G, N_config)
        avgE, sdevE = bootstrap(G, nbstrap=100)

        t_vals = [a * q for q in range(len(avgE))]

        plt.figure(figsize=(8, 5), dpi=120)
        #plt.axhline(y=1.0, color='black', linestyle='--', linewidth=1.0, label='Riferimento (1.0)')
        plt.errorbar(t_vals, avgE, yerr=sdevE, fmt='o', ecolor="blue", elinewidth=0.8,
                     capsize=4, markersize=4, label='Calculated ΔE')
        plt.scatter(t_vals, avgE, color='black', s=8)

        plt.title(f'Quartic Potential: N={N}, a={a}, N_config={N_config}, N_cor={N_cor}')
        plt.xlabel('t')
        plt.ylabel(r'$\Delta E(t)$')
        plt.xlim(0, min(9.5, N * a))
        plt.ylim(-3, 3)
        plt.grid(True)
        plt.legend()
        plt.tight_layout()

        filename = f"E2_Quartic_Nconfig{N_config}.png"
        plt.savefig(filename)
        plt.close()
        print(f"Plot salvato come '{filename}'\n")

if __name__ == "__main__":
    analysis(N_cf_list)
