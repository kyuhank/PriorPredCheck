# PriorPredCheck

This repository contains the code and data for the paper:

**Examining potential biases through prior predictive checks: Prior mis-specifications and their impact on Bayesian stock assessments**

*Kyuhan Kim, and Philipp Neubauer*

# Abstract

Bayesian stock assessment models are widely used to assess the status
of fish stocks and provide management advice. The Bayesian approach
allows for incorporating prior information into assessment models to
improve the precision of estimates and to constrain otherwise
over-parameterised models towards a solution that is aligned with
expectations for specific parameters. A prior distribution is specified
based on specific prior information or expert opinion, and a likelihood
function is used to update the prior to yield a posterior distribution
based on observed data. Understanding how the joint prior, defined
over estimated model parameters, interacts with
the model likelihood is essential for robust Bayesian
inference. However, Bayesian stock assessments often neglect this
important aspect of Bayesian inference. Depending on the model
structural assumptions and parameterisations, different prior
mis-specification problems can arise, which give rise to
well-documented problems. In this study, we show two common
prior mis-specifications, using three different illustrative stock
assessment models of increasing complexity applied to the South Atlantic albacore tuna
stock. The first mis-specification is an inconsistency between the prior and
likelihood function associated with an implicit change of the
prior, brought about by the model structure to avoid negative biomass calculations. The
second mis-specififcation arises from the use of uniform priors with
constrained support on harvest rates or
population scale parameters. Despite their intended lack of
information, these priors
are strongly informative for output space of a
model, regardless of model parameterisation. Our simulation studies demonstrate that failing to account for
the implications of these problems in a likelihood function can lead
to misleading inference. We also show that such misleading inference
can be avoided if prior information is encoded based on the prior
predictive distributions of quantities that are directly linked to the
likelihood function. To avoid potential prior
mis-informed inference in Bayesian stock assessments, we recommend
that prior predictive checks should be routinely used to understand
(and correct) potentially unintended interactions between the likelihood and the implicit information contained in joint priors for stock assessments.


# How to run the code

1. go to the folder `Analysis`

2. run the scripts by typing `makefile` in the terminal after setting the working directory to the subfolders, `SSPM`, and `ASPM`.

3. the order is as follows:

   i. `SSPM`
   ii. `ASPM`





