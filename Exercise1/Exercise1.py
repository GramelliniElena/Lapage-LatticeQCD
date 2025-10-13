import numpy as np
import matplotlib.pyplot as plt
import vegas

a = 0.5  #lattice spacing
N = 8  #number of steps
m = 1.0
T = a * N  #total time
E0_analytical = 0.5  #exact harmonic ground-state energy

def harmonic_potential(x):
    return 0.5 * m * x**2

def quartic_potential(x):
    return 0.5 * x**4

#Exact ground-state wave function
def psi0(x):
    return np.exp(- x**2 / 2) / (np.pi) ** 0.25

def action(x_path, x0, potential_choice):
    x_full = np.concatenate(([x0], x_path, [x0]))
    kinetic = np.sum((m / (2 * a)) * (x_full[1:] - x_full[:-1])**2)
    potential = a * np.sum(potential_choice(x_full[:-1])) 
    return kinetic + potential

def integrand(x_path, x0, potential_choice):
    prefactor = (m / (2 * np.pi * a)) ** (N / 2)
    return prefactor * np.exp(- action(np.array(x_path), x0, potential_choice))

#Numerical integration of the propagator with fixed endpoints (N-1 variables)
def evaluate_propagator(x0, potential_choice):
    integrator = vegas.Integrator([[-5, 5]] * (N - 1))  #interval for the integrations
    def f(x_path):
        return integrand(x_path, x0, potential_choice)
    result = integrator(f, nitn=10, neval=100000) 
    return result.mean

#Integrand for the trace (no fixed endpoints)
def integrand_trace(x_path, V, m, a, N):
    x_path = np.array(x_path)  # <-- convert to numpy array
    kinetic = (m / (2 * a)) * np.sum((x_path[1:] - x_path[:-1])**2)
    kinetic += (m / (2 * a)) * (x_path[0] - x_path[-1])**2
    potential = a * np.sum(V(x_path))
    prefactor = (m / (2 * np.pi * a)) ** (N / 2)
    return prefactor * np.exp(- (kinetic + potential))

def evaluate_trace(potential_choice):
    integrator = vegas.Integrator([[-5, 5]] * N)  
    def f(x_path):
        return integrand_trace(x_path, potential_choice, m, a, N)
    result = integrator(f, nitn=10, neval=100000)
    return result.mean

#Estimate E0 from the trace integral
def E0_estimate(potential_choice):
    Z = evaluate_trace(potential_choice)
    E0 = -np.log(Z) / T
    return E0

def run_analysis(potential_choice, potential_name="Potential", show_analytical=True):
    x0_values = np.linspace(0, 2, 15)  #generates 15 equidistant values between 0 and 2

    #Estimate ground-state energy from trace integral in N dimensions
    E0_numerical = E0_estimate(potential_choice)

    #Evaluate propagators with fixed endpoints
    propagators = [evaluate_propagator(x, potential_choice) for x in x0_values]
    
    #Analytical expressions (only for harmonic case)
    if show_analytical:
        propagator_analytical = [psi0(x)**2 * np.exp(-E0_analytical * T) for x in x0_values]

    plt.figure(figsize=(8, 4))
    plt.plot(x0_values, propagators, 'o', label="Numerical, a = 0.5")
    if show_analytical:
        plt.plot(x0_values, propagator_analytical, 'k--', label="Analytical (harmonic)")
    plt.title(f"Closed Euclidean Propagator — {potential_name}")
    plt.xlabel("$x_0$")
    plt.ylabel(r"$\langle x_0 | e^{-HT} | x_0 \rangle$")
    plt.ylim(-0.02, 0.1)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"E1_{potential_name.lower()}.png")

    print(f"Estimated ground state energy E0 (numerical) for {potential_name} = {E0_numerical:.5f}")

run_analysis(harmonic_potential, "Harmonic", show_analytical=True)
run_analysis(quartic_potential, "Quartic", show_analytical=False)
