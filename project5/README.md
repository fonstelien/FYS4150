# FYS4150 H20 Project 5 - Disease Modeling - SIRS Model

For this project, I have chosen to use Python for code development. Native Python running on a single thread is not ideal for scientific computing, but for this project, the computational demands are not that heavy, while the code is rather complex -- this makes Python a good choice, in my opinion. In the end, running the whole notebook that I have made for the project takes about 10 minutes, so I have gotten by just fine. Efficiency in code development also matters -- Horses for courses!

You will find the report and supporting files in the `report` folder and source code in the `source` folder. The source code consists of two files `sirs.py` and `utils.py`. `sirs.py` contains two classes `Sirs` and `SirsSolver`. `Sirs` implements the SIRS model, and `SirsSolver` solves the SIRS model. To run a simulation, you start by instantiating a `Sirs` object by sending in the main disease parameters `N, I0, a, b, c` to its constructor. Then you instantiate a `SirsSolver` object, which takes a `Sirs` instance as constructor argument. Then do a `SirsSolver.run_mc()` or `SirsSolver.run_rk4()`, with optional arguments for time, and steps in the RK4 case. You retrieve the results by doing `SirsSolver.get_fractions()`, which returns the fractions over time as a tuple `(S, I, R, N)`, which are all `numpy.ndarray`s. The `utils.py` file contains some helper functions for data preparation, as well as bulding blocks for creating custom functions for the seasonally varying factors `a` and `f`. I have been careful to write extensive docstrings which I hope will make the code easy to use. 

I recommend testing the code by running the examples I have prepared in `test-runs.ipynb`, which you will find in the `source` folder. You may also want to look at `project5.ipynb`, which contains all the plots from the report. You will find that in the `source` folder as well. I have included a copy of my `conda` environment in `fys4155.yml`. You install it by doing
```
$ conda env create -f fys4150.yml
```
... and get rid of it again with
```
$ conda env remove --name fys4150.yml
```

Enjoy! :whale: