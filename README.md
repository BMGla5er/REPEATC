# REPEATC - Rapid Exploration of Phase Eequilibria Across Temperature and Composition
- SARTRE -single and ready to equal
- TC-SCOPE - ThermoCalc Single-axis Calculation Over Phase Equilibria
- ThermALPS - Thermodynamic Alloy Landscape via Phase Simulation
- ThermoPACT - Phase Analysis by Composition and Temperature
- REPEATC - Rapid Exploration of Phase Equilibria Across Temperature and Compositions
- "Mapping Phase Evolution one Axis at a time"
Generic single axis equilibrium evaluator using ThermoCalc's TC-Python package.

When exploring a new alloy system, valuable exploration is phase development and transformation as a known set of conditions changes. When developing alloys for high-temperature, validating that strengthening phases are stable beyond the service temperature is a must. In our group, this has been the "first step" of any new metal, from diule Al-alloys to multi-principle element alloys (MPEAs), where we evaluate phases and distribution versus temperature, introduction of new elements, and changes to composition (such as along the ranqge of compositions).

This code is intended to accelerate these explorations by offering ag eneral script which, using a local ThermoCalc and TC-Python installation, will take an alloy, database, and set of input conditions (currently configured to be temepratures, "T", or fraction of a non-balancing element), then perform single equilibrium calculatiuons at each condition and automatically generate plots for phase fraction, composition of each phsae, and distribution of each element in phases.

The architecture of this is a custom class, "alloy_sys" with the following properties selections:
- dependent (balancing element, default: 'Al')
- composition (Dictionary containing alloying elements and initial composition in weight or mole fraction, default: {"Si":0.1,"Mg":0.0035})
- param (Axis for simulation, currently accepts temperature ("T") or an element in the composition ("Si","Mg"), default: 'T')
- param_range (Values for parameter to be varied on, equilibrium is performed at each condition, default: range(300,1300,100))
- unit (Character that describes if composition is in weight fraction ("W") or mole fraction ("X"), default: 'X')
- custom (Boolean value for if custom database or preloaded ThermoCalc database is used, default: False)
- database (If using ThermoCalc databse (custom is False) write string for databse selected, otherwise enter path to custom TDB file, default: 'TCAL9')
- suspended (List of phases to be suspended from calculation, default: None)

The class has four built-in functions:
- do_perform_single_axis_split
- phase_distribution() which plots the amount of phase at each composition
- composition_distribution() which creates subplots for each phase and shows composition trends
- element_distribution() which creates subplots for each element and shows which phases they are in
