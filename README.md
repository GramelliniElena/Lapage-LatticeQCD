# Lapage-LatticeQCD
These are my solutions to some of the exercises proposed by Lepage (https://arxiv.org/abs/hep-lat/0506036). The numbering matches the order in the article, except for Analysis.py which was not one of the suggested exercises.

# Code description
**Exercise1**: Solving path integral, both for harmonic and quartic potentials, using Vegas routine to integrate over closed paths; the calculation is performed considering different boundaries. The estimation of the ground state energy for the harmonic case is performed by Vegas integrating over all boundaries, as described in TNAN_LatticeQCD.pdf. My numerical results are compared with the exact ones.

**Exercise2_Harmonic**: Evaluation of the energy difference between first excited state and ground state, for the harmonic oscillator case, considering the correlator <x(t+a)x(t)>. The path is randomized with a Metropolis algorithm and errors are estimated via statistical bootstrap. My numerical results are compared with the exact ones.

**Excercise2_Quartic**: Evaluation of the energy difference between first excited state and ground state, for the quartic oscillator case, considering the correlator <x(t+a)x(t)>. The path is randomized with a Metropolis algorithm and errors are estimated via statistical bootstrap.

**Exercise3_Harmonic**: Evaluation of the energy difference between first excited state and ground state, for the harmonic oscillator case, considering the correlator <x^3(t+a)x^3(t)>. The path is randomized with a Metropolis algorithm and errors are estimated via statistical bootstrap. My numerical results are compared with the exact ones.

**Excercise3_Quartic**: Evaluation of the energy difference between first excited state and ground state, for the quartic oscillator case, considering the correlator <x^3(t+a)x^3(t)>. The path is randomized with a Metropolis algorithm and errors are estimated via statistical bootstrap.

**Exercise4_Harmonic**: Improvement of the error estimations in Exercise2_Harmonic by implementing a binning procedure. My numerical results are compared with the exact ones.

**Exercise4_Quartic**: Improvement of the error estimations in Exercise2_Quartic by implementing a binning procedure.

**Exercise5_Harmonic**: Substitution of an improved form for the action in Exercise4_Harmonic; details of the form of the new action can be found in TNAN_LatticeQCD.pdf. My numerical results are compared with the exact ones.

**Exercise5_Quartic**: Substitution of an improved form for the action in Exercise4_Quartic; details of the form of the new action can be found in TNAN_LatticeQCD.pdf.

**Exercise6_Harmonic**: Change of the harmonic potential form, by means of a change of variables that avoids the so-called "ghost modes". My numerical results are compared with the exact ones.

**Exercise6_Anharmonic**: Change of the anharmonic potential form, by means of a change of variables that avoids the so-called "ghost modes".

**Analysis**: Comparison between the different estimations of DeltaE depending on the choices of action discretization ('Standard', 'Improved' and 'No-ghost'). These results are also compared with the exact ones.
