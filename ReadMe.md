The codes in this repository can be used to reproduce the results in the paper

# A two-stage drop-the-losers design for time-to-event outcome using a historical control arm

by R. Abbas, J. Wason, S. Michiels, and G. Le Teuff (2020)

Phase II clinical trials are simulated according to a drop-the-losers and to a fixed design both using historical control arm.

# Drop-the-losers design with historical control arm

DTLHC_functions.r 

This file loads functions required to run the DTLHC_main.r simulation program.

DTLHC_main.r 

This file runs the two-stage drop-the-losers design for time-to-event outcome using a historical control arm and gives power or type I error rate as shown in tables 2-4.

# Drop-the-losers design with misspecified historical control arm

DTLHC_misspecified_functions.r 

This file loads functions required to run the DTLHC_misspecified_main.r simulation program.

DTLHC_misspecified_main.r 

This file runs the two-stage drop-the-losers design for time-to-event outcome in case of misspecification of the historical control arm and gives power or type I error rate as shown in tables 5-6.

# Fixed design (Multi-Arm design with No Interim analysis) with or without correction for multiplicity

MANI_functions.r 

This file loads functions required to run the MANI_main.r simulation program.

MANI_main.r 

This file runs the multi-arm fixed design and gives power or type I error rate as shown in tables A1-A2 in appendix.

# Simulation results: Power

In a multi-arm setting the concept of power is complex due to multiplicity of hypotheses. The probability to reject all false null hypothese is called the conjunctive power. The probability to reject at least one false null hypothesis is called the disjunctive power. Usually, it is harder to obtain a high conjunctive power. Here the drop-the-losers design do not allow to reject more than one null hypothesis, we report disjunctive power.
 
# Simulation results: Type 1 error rate or family-wise error rate

By design, the family-wise error rate is controlled at the global null hypothesis. Here we report the type 1 error rate per scenario. The scenario 1, considering only experimental arms with null effect (also called the global null hypothesis) reflects the empirical family-wise error rate. 
