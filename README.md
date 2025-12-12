# Bayesian Causal Dynamic Borrowing

This code was developed by David Jesse as part of a PhD project at the University of Göttingen and F. Hoffmann-La Roche AG.
It accompanies a publication titled "Bayesian Methods Integrating Causal Inference Approaches for Borrowing Historical Control Data in RCTs: A Neutral Comparison Study" that is currently in review.
This work was co-funded by the European Union’s Horizon Europe Framework programme under grant agreement 101136365 (INVENTS), co-funded by the Swiss State Secretariat for Education, Research and Innovation (SERI) and co-funded by the UKRI Innovative UK under their Horizon Europe Guarantee scheme.
The code is the responsibility of the authors and not necessarily endorsed by the institutions mentioned.
The code is provided "as is", without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and noninfringement.


## Methods

We compare the following three methods:

- Propensity score integrated commensurate prior
  - Paper by [Wang et al. (2022)](https://doi.org/10.1080/10543406.2021.2011743)
- Propensity score weighted multi-source exchangeability models (PS-MEM)
  - Paper by [Wei et al. (2024)](https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.10158)
- Bayesian additive regression trees (BART)
  - Paper by [Zhou and Ji (2021)](https://onlinelibrary.wiley.com/doi/full/10.1002/sim.9191)


## Structure of the repository

The repository has a normal R project structure (using [renv](https://rstudio.github.io/renv/articles/renv.html)) but it is not a package.
To nicely organize the code and handle dependencies, we use the [box package](https://klmr.me/box/).
Here is a description of where to find what:

- `/R`: any user-written R functions, in particular...
  - `/R/borrow`: R functions for the different borrowing methods
  - `/R/simfuns`: R functions for the simulation study (e.g. data-generating models, utilities, ...)
  - ...
- `/data`: (dummy) example data sets
- `/playground`: place for experimenting and testing out newly developed functions or other stuff
    