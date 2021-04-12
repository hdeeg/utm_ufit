This directory contains codes that generate models that can be used for fitting and which UFIT may invoke instead of UTM. They are much simpler than UTM, but they cannot be run as standalone modellers. In order use them to for the generation of a model without fitting, use UFIT,'setup.utm',/nofit

For example setups to run these modellers, see /examples/alt_modellers

- ufit_polynome.pro : modeller for polynomes of any degree
- ufit_lin_sine.pro : modeller for a superposition of a linear and a sinusoidal function

