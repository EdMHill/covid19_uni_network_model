# covid19_uni_network_model

This repository contains files for performing computational simulations of a network model framework to explore SARS-CoV-2 transmission amongst a university population.

The code was developed for the analysis presented in the scientific paper "Modelling SARS-CoV-2 transmission in a UK university setting" by Edward M. Hill, Benjamin D. Atkins, Matt J. Keeling, Michael J. Tildesley and Louise Dyson.

Preprint details: Hill et al. (2020) Modelling SARS-CoV-2 transmission in a UK university setting. *medRxiv* doi: 10.1101/2020.10.15.20208454. URL: https://doi.org/10.1101/2020.10.15.20208454.

Model simulations are performed using the programming language Julia.
Julia makes use of environments, allowing bespoke package lists for separate projects. Documentation on working with environments and installing packages in the same state that is given by the project manifest: https://julialang.github.io/Pkg.jl/v1.5/environments/#Using-someone-else's-project-1

Please find below an explainer of the directory structure within this repository.

## Data
Directory containing files that may be read in during model runs.

## Results
Directory to store simulation outputs

## src

**university_model.jl**  
The main model run file. Pass requested configuration to uni_network_run (in include_files_uni_model/main_function.jl) and saves outputs into an MAT file.

**uni_model_demo.ipynb**  
Jupyter notebook giving example run of the model.

**include_files_uni_model**  
Houses function files to be used when running the uni model.

- **main_function.jl**
    Outline of the code structure:  
    * Unpack required variables
    * Set the random number generator
    * Initialise students into cohorts, split by department and into first year undergrad, non-first year undergrad, postgrad (supporting functions in network_generation_fns.jl)
    * Assign society/sports club membership  (supporting functions in network_generation_fns.jl)
    * Generate study/cohort contacts, society contacts. (supporting functions in network_generation_fns.jl)
    * Assign households and household contacts (supporting functions in network_generation_fns.jl)
    * Generate social dynamic contacts (supporting functions in network_generation_fns.jl)
    * Generate on-campus accomodation dynamic contacts (supporting functions in network_generation_fns.jl)
    * Generate household specific transmission (supporting functions in network_generation_fns.jl)
    * Initialise variables
    * Iterate over replicates:
        - Reinitialisation phase (supporting functions in additional_fns.jl)
        - Set course of infection times (supporting functions in additional_fns.jl)
        - Set initial seed infections & those already recovered (supporting functions in seed_initial_states_fn.jl)
        - Update output time series with initial conditions
        - Reset contact tracing variables (supporting functions in additional_fns.jl)
        - Iterate over time
            - Reinitialise variables at start of timestep
            - Assign outputs
            - Increment counters (supporting functions in additional_fns.jl)
            - Increment infection process (if come to the end of latent time, move to infectious time etc.). Includes household infection loop. (supporting functions in additional_fns.jl)
            - Record whether individuals are in isolation. Set study and society attendance status for timestep.
            - Transmit infections (supporting functions in additional_fns.jl)
            - Perform contact tracing (supporting functions in contact_tracing_fns.jl)
            - Assign prevalence & isolation outputs
            - Reactive inactivation of in person classes/society meet ups
            - Run interventions (supporting functions inintervention_condition_affect_fns.jl)
            - Perform mass testing, if applicable (supporting functions in mass_testing_fns.jl)
    * Return outputs

- **additional_fns.jl**   
    Stash of functions to run the university model. Includes:  
    * Time in state increment fn (increment_counters!)
    * load configs for sensitivity runs, generates (load_configs)
    * find_network_parameters (Load relevant network params based on number of cohorts requested)
    * transmission functions
    * functions to set up transmission rates within household for each individual
    * Functions to reinitialise states at start of each simulation replicate
    * miscellaneous functions

- **seed_initial_states_fn.jl**   
    Function specifying how the initial disease states will be assigned each simulation replicate.

- **contact_tracing_fns.jl**  
    Functions that are used with the university network model for performing contact tracing. Includes:
    * Get portion of dynamic contacts to be recallable (recallable_dynamic_contacts)
    * Check contacts made in study setting (get_study_contacts)
    * Check contacts made in society setting (get_society_contacts)
    * Perform forward CT from an identified infector (forwardCT_from_infector! and trace_node!)

- **intervention_condition_affect_fns.jl**  
    Functions to implement/rescind trigger based interventions

- **mass_testing_fns.jl**  
    Mass testing module

- **network_generation_fns.jl**  
    Produces the network layers

- **parametertypes.jl**  
    Defines the parameter types used with the network model for universities. Fields accessible with dot notation. Example using type student_params, with a variable named student_info. student_info.cohort_ID accesses the value in the cohort_ID field.
