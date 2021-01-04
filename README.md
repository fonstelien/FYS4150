# FYS4150 Computational Physics
This repository contains my reports and source code for each of the five projects in the course FYS4150 Computational Physics given at University of Oslo (UiO) in the fall of 2020. From [UiO's description of the course](https://www.uio.no/studier/emner/matnat/fys/FYS4150/index-eng.html):

*An introduction to numerical methods which are used in solving problems in physics and chemistry, including solutions of differential equations, matrix operations and eigenvalue problems, interpolation and numerical integration, modeling of data and Monte Carlo methods.*

The final grade in the course was given based on the last three projects:
* [Project 3](https://github.com/fonstelien/FYS4150/blob/master/project3/report/project3.pdf): *Celestial Mechanics with Numerical Methods* (score 97/100)

    The classical $n$-body problem. We model the Solar System with the Velocity Verlet method and investigate how the Solar System would look like if there were small changes in Newton's laws, as well as reproducing the [perihelion precession of Mercury](https://en.wikipedia.org/wiki/Tests_of_general_relativity#Perihelion_precession_of_Mercury).
* [Project 4](https://github.com/fonstelien/FYS4150/blob/master/project4/report/project4.pdf): *The 2-Dimensional Ising-Model and Monte Carlo Simulations* (score 94/100)

    An exercise in statistical physics. We implement the 2D [Ising Model](https://en.wikipedia.org/wiki/Ising_model) for magnetic orientation of the spins in a square lattice as a Markov process driven forward with Monte Carlo simulations and the Metropolis algorithm as acceptance criteria.

* [Project 5](https://github.com/fonstelien/FYS4150/blob/master/project5/report/project5.pdf): *Solutions to the SIRS Disease Model with Runge-Kutta and Monte Carlo Methods* (score 99/100)
    
    Susceptible &#8594; Infectous &#8594; Recovering &#8594; Susceptible ([SIRS](https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#The_SIR_model)). Discrete modelling suitable for small populations with Monte Carlo methods; and continuous modelling for large populations with Runge-Kutta. We investigate the effect of seasonality in transmission and vaccination rates, as well as the effect if government-imposed social lockdowns -- Highly relevant these days!