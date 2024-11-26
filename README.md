This repository contains the code and data for the paper:

**Examining potential biases through prior predictive checks: Prior mis-specifications and their impact on Bayesian stock assessments**

*Kyuhan Kim, and Philipp Neubauer*

# Abstract

Bayesian stock assessment models are widely used to evaluate fish
stock status and inform management decisions. The Bayesian approach
enables the incorporation of prior information into assessment models,
improving estimate precision and constraining over-parameterised
models toward solutions that align with prior expectations. Understanding
the interaction between joint priors across model parameters and the model likelihood is
crucial for robust Bayesian inference, yet this aspect is seldom
addressed in applied Bayesian stock assessments. Depending on the model,
structural assumptions, and parameterisations, various prior
mis-specification issues may emerge, leading to well-documented
problems. In this study, we present two common prior
mis-specifications using three stock assessment models of increasing
complexity applied to South Atlantic albacore tuna data. The first
mis-specification stems from an inconsistency between the prior and
likelihood function, where the model structure implicitly modifies the
prior to avoid negative biomass calculations. The second
mis-specification involves uniform priors with constrained support on
parameters like harvest rates or population scale, which, despite
their intended non-informativeness, are highly informative for derived
quantities from the model. Simulations show that failing to address these
issues in the likelihood function can result in misleading
inference. We demonstrate that such issues can be mitigated by
carefully encoding prior information with the help from prior predictive distributions of
quantities directly linked to the likelihood function. To prevent
potentially misinformed inference in Bayesian stock assessments, we
recommend routinely conducting prior predictive checks to identify and
correct unintended interactions between the likelihood and implicit
information in joint priors.


# How to run the code

1. go to the folder `Analysis`

2. run the scripts by typing `makefile` in the terminal after setting the working directory to the subfolders, `SSPM`, and `ASPM`.

3. the order is as follows:

   i. `SSPM`
   ii. `ASPM`

*Running this entire process on a single laptop might take several days to complete due to the hundreds of simulations and estimation processes involved in the self-consistency check*
