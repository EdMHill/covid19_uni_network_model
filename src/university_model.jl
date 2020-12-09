"""
Purpose:
A network model to explore transmission amongst a university population
"""

"""
Set paths & load environment
"""

#Set paths
cd(dirname(@__FILE__))

#Load environment
using Pkg
Pkg.activate("../")

"""
Load packages
"""
#Required packages
using MAT, Distributions
using LinearAlgebra, Random, DelimitedFiles, Parameters

"""
Load supporting files
"""
# Files containing other functions needed to run the model
include("include_files_uni_model/parametertypes.jl")
include("include_files_uni_model/contact_tracing_fns.jl")
include("include_files_uni_model/network_generation_fns.jl")
include("include_files_uni_model/additional_fns.jl")
include("include_files_uni_model/intervention_condition_affect_fns.jl")
include("include_files_uni_model/mass_testing_fns.jl")
include("include_files_uni_model/seed_initial_states_fns.jl")
include("include_files_uni_model/student_travel_fns.jl")
include("include_files_uni_model/main_function.jl")

"""
Set variables from ARGS
"""
args = ARGS

# If running locally from REPL, will not have any input ARGS
# Can set values for input parameters to main run function here
if length(ARGS)==0
    args = [ "1", "100", "20", "77", "25000", "7155", "84", "seed_states_using_ODEmodel", "[\"run_one_run_baseline_nointerv\"]"]
end

# To run from command line, example:
# julia university_model.jl 1 100 100 77 25000 7155 84 seed_states_using_ODEmodel '["runsets"]'

# Runset options, defined in file "include_files_uni_model/additional_fns.jl"
#       '["run_one_run_baseline_nointerv","run_one_run_baseline"]'
        #["sweep_adherence", "sweep_adherence_with_rehouse",
        #"alter_test_result_delay",
        #"change_teaching_space_transmission","amount_backwards_CT",
        #"class_CT_threshold", "CT_engagement",
        #"asymp_scen", "change_scaling",
        #"vary_household_allocation_strat",
        #"vary_accom_level_lockdown",
        #"run_mass_testing", "run_mass_testing_6configs",
        #"run_mass_testing_alter_adherence","run_mass_testing_frequent",
        #"run_baseline_mass_testing",
        #"run_no_interventions", "run_one_run_baseline_nointerv",
        #"change_recov_propn",
        #"rehouse_strat_check",
        #"trans_risk_scale_no_interv",
        #"probasymp_scale_no_interv","transasymp_scale_no_interv",
        #"suscep_scale_no_interv"]

# Set identifier for job
job_ID = parse(Int64, args[1])

# Set RNG seed
RNGseed_base = parse(Int64, args[2])
RNGseed = job_ID*RNGseed_base

# Set simulation run params
countfinal = parse(Int64, args[3])  # Number of simulations to be performed per scenario
endtime = parse(Int64, args[4]) # Timesteps for each individual simulation
n_students = parse(Int64, args[5]) # Number of students (7155 on campus, 27278 overall in 2018/2019.)
                    # Estimate 25,000 students in local region (~24,500 based on 2019/2020 records)
n_students_on_campus = parse(Int64, args[6])
n_cohorts = parse(Int64, args[7]) # Defining total number of study cohorts

# Specify function for setting proportion/number of nodes in each initial disease state
# Relevant files in "include_files_uni_model/seed_initial_states_fn.jl"
# e.g. seed_states_using_ODEmodel
s = Symbol(args[8])
seed_initial_states_fn = getfield(Main, s) #Make Symbol a callable function

# Specify scenario to be run
runsets = eval(Meta.parse(args[9]))

"""
Set up incubation & infectivity distributions
"""

# If needed, set up a different incubation period distribution
# Default: Erlang(6,0.88)
# d_incub_alt = LogNormal(log(5.2), log(1.4))
    # If using alternative then need to pass d_incub_alt into infection_params creation:
    #e.g. infection_parameters = infection_params(transrisk=0.75*ones(n_cohorts),
    #                                                d_incub = d_incub_alt)

# If needed, set a different distribution of infectivity from the default
# Will then need to pass dist_infectivity into infection_params creation:
#  e.g. infection_parameters = infection_params(transrisk=0.75*ones(n_cohorts),
#                                                  dist_infectivity = dist_infectivity)
# dist_infectivity = ones(10) # if we don't want to have a distribution of infectiousness


"""
Set up adherence delay distribution
"""
# Set a different distribution for delay in adhering in guidance
# Will need to pass delay_adherence_pmf_alt into infection_params creation:
#  e.g. infection_parameters = infection_params(...,
#                                                  delay_adherence_pmf = delay_adherence_pmf_alt)
# delay_adherence_pmf_alt = [0.5,0.25,0.25]
#
#  # Check delay_adherence_pmf_alt sums to 1
# if sum(delay_adherence_pmf_alt) != 1
#     error("delay_adherence_pmf_alt must sum to 1. Currently sums to $(sum(delay_adherence_pmf_alt))")
# end

"""
Set up testing related distributions
"""

# If needed, set a different probability mass function for delay until test result received.
# Will then need to pass CT_delay_until_test_result_pmf_alt into CT_params creation:
#  e.g. CT_parameters = CT_params(...,
#                                  CT_delay_until_test_result_pmf = CT_delay_until_test_result_pmf_alt)
# CT_delay_until_test_result_pmf_alt = [0.,0.,0.5,0.,0.,0.,0.5]
#
#  # Check CT_delay_until_test_result_pmf_alt sums to 1
# if sum(CT_delay_until_test_result_pmf_alt) != 1
#     error("CT_delay_until_test_result_pmf_alt must sum to 1. Currently sums to $(sum(CT_delay_until_test_result_pmf_alt))")
# end

# If needed, set an alternative false negative test prob. , wrt days since infected
# Will need to pass test_false_negative_vec_alt into CT_params creation:
#  e.g. CT_parameters = CT_params(...,
#                                 test_false_negative_vec = test_false_negative_vec_alt)
# test_false_negative_vec_alt = 0.2*ones(20)

"""
Specify use of any additional, trigger interventions
"""
# Have as a vector input.
# -> Entry per intervention
# -> Within intervention fns, can call independently made condition fns
# -> Will perform affect if condition is satisfied

# intervention_fns_alt = [society_incidence_check!,accom_lockdown_block_level!]
# intervention_fns_alt = [affect_close_example!;
#                         affect_open_example!]

    # Will need to pass as optional input to worker_pattern_network_run:
    #  e.g. ...=  worker_pattern_network_run(...,
    #                                 intervention_fns = intervention_fns_alt)



"""
Specify function to allocate students to department/cohort & classes
"""
# If needed, set a different class assignment function from the default
generate_classes_fn_alt = generate_classes_all_students # generate_classes_campus_only
#         # Will need to pass as optional input to worker_pattern_network_run:
#         #  e.g. ...=  worker_pattern_network_run(...,
#                                 generate_classes_fn = generate_classes_fn_alt)

"""
Specify function to assign households and construct household contact network
"""
# If needed, set a different household assignment function from the default
# generate_student_households_fn_alt = assign_households_campus_only #assign_households_all_students
# generate_student_households_fn_alt = assign_households_by_cohort # assign_households_all_students

#         # Will need to pass as optional input to worker_pattern_network_run:
#         #  e.g. ...=  worker_pattern_network_run(...,
#                                 generate_student_households_fn! = generate_student_households_fn_alt)


"""
Specify household transmission risk by group
"""
# # If needed, set a different distribution of transrisk_household_group from the default
# # Will then need to pass transrisk_household_group_alt into infection_params creation:
# #  e.g. infection_parameters = infection_params(...,
# #                                                  transrisk_household_group = transrisk_household_group_alt)
# transrisk_household_group_alt = [1.,0.5,0.2]
#
# # Specify function to perform household group allocation
assign_household_transrisk_fn_alt = assign_household_transmit_household_size!
#         # Will need to pass as optional input to worker_pattern_network_run:
#         #  e.g. ...=  worker_pattern_network_run(...,
#         #                                 assign_household_transrisk_fn = assign_household_transrisk_fn_alt)

"""
Specify society assignment function
"""
# If needed, set a different assign_societies_fn from the default
assign_societies_fn_alt = assign_societies_from_aggregated_data

#         # Will need to pass as optional input to worker_pattern_network_run:
#         #  e.g. ...=  worker_pattern_network_run(...,
#
      #                                 assign_societies_fn = assign_societies_fn_alt)

## Iterate over the various configurations
for run_it = 1:length(runsets)

    # Set if contact tracing is active or not (Bool type variable)
    contact_tracing_active = false

    # Set if class closures is active or not (Bool type variable)
    work_study_group_closure_active = false

    # set if backwards contact tracing is active or not (Bool type variable)
    perform_CT_from_infector = false

    runset = runsets[run_it]

    # Specify if mass testing module is in use
    if (runset == "run_mass_testing") || (runset == "run_mass_testing_6configs") ||
         (runset == "run_mass_testing_alter_adherence") || (runset == "run_mass_testing_frequent")
        mass_testing_active = true
    else
        mass_testing_active = false
    end

    adherence_config, sameday_config, ton_config, toff_config, work_perc_config,
    num_config, work_or_study_group_CT_threshold_config, prob_backwards_CT_config,
    infector_engage_with_CT_prob_config, CT_engagement_config, probasymp_config,
    transasymp_config, scaling_config, suscep_config,
    n_students_config, recov_propn_config,
    rehouse_strat_config, # RNGseed_config,
    CT_delay_until_test_result_config,
    intervention_fn_config, assign_household_fn_config,
    mass_testing_config = load_configs(runset,countfinal,n_cohorts,n_students)

    # Initialise output arrays. Store counts for each network configuration
    # put countfinal in the 2nd place so that we can concatenate variables later
    numlat_save = zeros(endtime+1,countfinal,num_config)
    numinf_save = zeros(endtime+1,countfinal,num_config)
    numrep_save = zeros(endtime+1,countfinal,num_config)
    prevlat_save = zeros(endtime+1,countfinal,num_config)
    prevsymp_save = zeros(endtime+1,countfinal,num_config)
    prevasymp_save = zeros(endtime+1,countfinal,num_config)
    prevpresymp_save = zeros(endtime+1,countfinal,num_config)
    prevrec_save = zeros(endtime+1,countfinal,num_config)
    newinf_save = zeros(endtime+1,countfinal,num_config)
    atworkinf_save = zeros(endtime+1,countfinal,num_config)
    atworkasymp_save = zeros(endtime+1,countfinal,num_config)
    newasymp_save = zeros(endtime+1,countfinal,num_config)
    num_CT_save = zeros(endtime+1,countfinal,num_config)             # Number contact traced on given day
    num_infected_save = zeros(Int64,n_students,countfinal,num_config)
    social_dynamic_infection_count_save = zeros(Int64,endtime+1,countfinal,num_config)
    accom_dynamic_infection_count_save = zeros(Int64,endtime+1,countfinal,num_config)
    household_infection_count_save = zeros(Int64,endtime+1,countfinal,num_config)
    var_num_infected_save = zeros(Float64,1,countfinal,num_config) # so that the countfinal is in the 2nd place
    mean_init_generation_time_save = zeros(Float64,1,countfinal,num_config) # so that the countfinal is in the 2nd place
    num_isolating_CTcause_save = zeros(endtime+1,countfinal,num_config)
    Rt_save = zeros(endtime+1,countfinal,num_config)
    new_rehoused_save = zeros(endtime+1,countfinal,num_config)
    current_rehoused_save = zeros(endtime+1,countfinal,num_config)
    num_init_infected_save = Array{Array{Array{Int64,1}}}(undef,num_config)

    # Initialise isolation outputs
    num_isolating_save = zeros(endtime+1,countfinal,num_config)      # Number isolating on given day
    num_isolating_oncampus_save = zeros(endtime+1,countfinal,num_config)      # Number oncampus resident students isolating on given day
    num_isolating_offcampus_save = zeros(endtime+1,countfinal,num_config)      # Number offcampus resident students isolating on given day
    num_symp_isolating_save = zeros(endtime+1,countfinal,num_config) # Number isolating due to having symptoms
    num_asymp_isolating_save = zeros(endtime+1,countfinal,num_config) # Number isolating due to finding asymptomatic infection (via testing)
    num_household_isolating_save = zeros(endtime+1,countfinal,num_config) # Number isolting due to household infection
    num_accom_lockdown_isolating_save = zeros(endtime+1,countfinal,num_config)

    # Initialise testing output arrays
    tests_performed_save = zeros(endtime+1,countfinal,num_config)
    test_outcomes_save = zeros(endtime+1,countfinal,4,num_config)

    # Initialise cohort 3D output arrays
    cohort_infection_count_save = zeros(Int64,endtime+1,countfinal,n_cohorts,num_config)
    society_infection_count_save = zeros(Int64,0,0,0)

    # Initialise 1D infection count arrays
    n_oncampus_inf_save = zeros(countfinal,num_config)
    n_offcampus_inf_save = zeros(countfinal,num_config)

    # Initialise 1D isolation count arrays
    n_isol_adhering_save = zeros(countfinal,num_config)
    n_isol_adhering_oncampus_save = zeros(countfinal,num_config)
    n_isol_adhering_offcampus_save = zeros(countfinal,num_config)

    # If applicable, initialise mass testing output variables
    if mass_testing_active == true
        mass_testing_n_tests_save = Array{Array{Array{Int64,1},1},1}(undef,num_config)
        mass_testing_n_positive_save = Array{Array{Array{Int64,1},1},1}(undef,num_config)
        mass_testing_n_all_isolating_save = Array{Array{Array{Int64,1},1},1}(undef,num_config)
        mass_testing_n_hh_isolating_save = Array{Array{Array{Int64,1},1},1}(undef,num_config)
        mass_testing_n_CT_isolating_save = Array{Array{Array{Int64,1},1},1}(undef,num_config)
        mass_testing_n_prev_infected_tested_save = Array{Array{Array{Int64,1},1},1}(undef,num_config)
    end

    for it = 1:num_config

        adherence = adherence_config[it]
        sameday = sameday_config[it]
        ton = ton_config[it]
        toff = toff_config[it]
        attendence_propns = work_perc_config[it,:,:]
        work_or_study_group_CT_threshold = work_or_study_group_CT_threshold_config[it]
        prob_backwards_CT = prob_backwards_CT_config[it]
        infector_engage_with_CT_prob = infector_engage_with_CT_prob_config[it]
        CT_engagement = CT_engagement_config[it]
        prob_asymp = probasymp_config[it]
        asymp_trans_scaling = transasymp_config[it]
        scaling = scaling_config[it]
        suscep_scaling = suscep_config[it]
        local n_students = n_students_config[it]
        recov_propn = recov_propn_config[it]
        rehouse_strat_active = rehouse_strat_config[it]
        # RNGseed = RNGseed_config[it]
        CT_delay_until_test_result_pmf_alt = CT_delay_until_test_result_config[it,:]
        intervention_fns_alt = intervention_fn_config[it,:]
        generate_student_households_fn_alt = assign_household_fn_config[it]
        mass_testing_parameters_alt = mass_testing_config[it]

        if runset=="amount_backwards_CT"
            perform_CT_from_infector = true
            contact_tracing_active = true
        end

        if runset=="class_CT_threshold"
            work_study_group_closure_active = true
        end

        if (runset=="CT_engagement")||(runset=="amount_backwards_CT")||
            (runset=="sweep_adherence")||(runset=="sweep_adherence_with_rehouse")||
            (runset=="alter_test_result_delay")||
            (runset=="change_teaching_space_transmission")||
            (runset=="vary_accom_level_lockdown")||
            (runset=="run_one_run_accom_lockdown_baseline") ||
            (runset=="vary_household_allocation_strat") ||
            (runset=="run_mass_testing") ||
            (runset=="run_mass_testing_6configs") ||
            (runset=="run_mass_testing_alter_adherence") ||
            (runset=="run_mass_testing_frequent") ||
            (runset=="run_baseline_mass_testing") ||
            (runset=="run_one_run_baseline") ||
            (runset=="rehouse_strat_check")
            contact_tracing_active = true
        end

        # Set society & sports club transmission risk.
        if (runset=="run_no_interventions") || (runset=="run_one_run_baseline_nointerv") ||
            (runset=="trans_risk_scale_no_interv") || (runset=="probasymp_scale_no_interv") ||
            (runset=="transasymp_scale_no_interv")
            # No interventions scenario.
            # Society risk matches social contacts. Sports club matches household risk.
            society_sports_transrisk_mean = [0.2414,0.34]
        else
            # Transmission risk with interventions in place
            society_sports_transrisk_mean = [0.12, 0.2414]
        end

        # Set up network parameters
        # Uses requested number of cohort types
        # For options, the function resides in "incldue_files_uni_model/additional_fns.jl"
        network_parameters, team_generation_parameters = find_network_parameters(n_cohorts,
                                                                                attendence_propns=attendence_propns,
                                                                                n_students = n_students,
                                                                                n_students_on_campus = n_students_on_campus
                                                                                )

        CT_parameters = CT_params(prob_backwards_CT = prob_backwards_CT,
                                    perform_CT_from_infector = perform_CT_from_infector,
                                    infector_engage_with_CT_prob = infector_engage_with_CT_prob,
                                    work_or_study_group_CT_threshold = work_or_study_group_CT_threshold,
                                    CT_engagement = CT_engagement,
                                    CT_delay_until_test_result_pmf = CT_delay_until_test_result_pmf_alt
                                    #test_false_negative_vec = test_false_negative_vec_alt
                                    )

        infection_parameters = infection_params(n_cohorts = n_cohorts,
                                                transrisk_household_18_34_group_mean=scaling,
                                                suscep_scaling = suscep_scaling,
                                                probasymp_dist = prob_asymp,
                                                asymp_trans_scaling_dist = asymp_trans_scaling,
                                                adherence = adherence,
                                                recov_propn = recov_propn,
                                                society_sports_transrisk_mean = society_sports_transrisk_mean,
                                                #delay_adherence_pmf = delay_adherence_pmf_alt
                                                #transrisk_household_group = transrisk_household_group_alt
                                                )

        # Establish society generation parameters
        society_generation_parameters = society_generation_params()

        # Call function (located in include_files_network_model\main_function.jl)
        if contact_tracing_active==true
                @time  output = uni_network_run(RNGseed,
                                                n_students,
                                                ton,toff,
                                                infection_parameters,
                                                sameday,
                                                seed_initial_states_fn,
                                                countfinal,
                                                endtime,
                                                contact_tracing_active,
                                                CT_parameters,
                                                network_parameters,
                                                team_generation_parameters,
                                                society_generation_parameters,
                                                work_study_group_closure_active,
                                                mass_testing_active,
                                                rehouse_strat_active,
                                                mass_testing_parameters = mass_testing_parameters_alt,
                                                intervention_fns = intervention_fns_alt,
                                                generate_student_households_fn = generate_student_households_fn_alt,
                                                assign_societies_fn = assign_societies_fn_alt,
                                                generate_classes_fn = generate_classes_fn_alt,
                                                assign_household_transrisk_fn = assign_household_transrisk_fn_alt
                                                )
                @unpack numlat, numinf, numrep,
                        prevlat, prevsymp, prevasymp, prevpresymp, prevrec,
                        newinf, newasymp, atworkinf, atworkasymp, infected_by,
                        num_isolating, num_isolating_oncampus, num_isolating_offcampus,
                        num_household_isolating, num_symp_isolating, num_asymp_isolating,
                        num_isolating_CTcause,
                        num_accom_lockdown_isolating,
                        num_CT, num_infected, social_dynamic_infection_count,
                        accommodation_dynamic_infection_count, household_infection_count,
                        var_num_infected, num_init_infected, mean_init_generation_time,
                        Rt,
                        new_rehoused, current_rehoused,
                        tests_performed, test_outcomes,
                        cohort_infection_count, society_infection_count,
                        n_oncampus_inf, n_offcampus_inf,
                        n_isol_adhering, n_isol_adhering_oncampus, n_isol_adhering_offcampus = output
           else
            # Profile.clear_malloc_data()
               @time  output = uni_network_run(RNGseed,
                                               n_students,
                                               ton,toff,
                                               infection_parameters,
                                               sameday,
                                               seed_initial_states_fn,
                                               countfinal,
                                               endtime,
                                               contact_tracing_active,
                                               CT_parameters,
                                               network_parameters,
                                               team_generation_parameters,
                                               society_generation_parameters,
                                               work_study_group_closure_active,
                                               mass_testing_active,
                                               rehouse_strat_active,
                                               mass_testing_parameters = mass_testing_parameters_alt,
                                               intervention_fns = intervention_fns_alt,
                                               generate_student_households_fn = generate_student_households_fn_alt,
                                               assign_societies_fn = assign_societies_fn_alt,
                                               generate_classes_fn = generate_classes_fn_alt,
                                               assign_household_transrisk_fn = assign_household_transrisk_fn_alt
                                               )
               @unpack  numlat, numinf, numrep,
                        prevlat, prevsymp, prevasymp, prevpresymp, prevrec,
                        newinf, newasymp, atworkinf, atworkasymp, infected_by,
                        num_isolating, num_isolating_oncampus, num_isolating_offcampus,
                        num_household_isolating, num_symp_isolating, num_asymp_isolating,
                       num_accom_lockdown_isolating,
                       num_CT, num_infected,
                       social_dynamic_infection_count,
                       accommodation_dynamic_infection_count, household_infection_count,
                       var_num_infected,
                       num_init_infected, mean_init_generation_time, Rt,
                       new_rehoused, current_rehoused,
                       tests_performed, test_outcomes,
                       cohort_infection_count,
                       society_infection_count,
                       n_oncampus_inf, n_offcampus_inf,
                       n_isol_adhering, n_isol_adhering_oncampus, n_isol_adhering_offcampus = output
                   end

            # For this iteration's network config, write results to output storage arrays
            numlat_save[:,:,it] = numlat
            numinf_save[:,:,it] = numinf
            numrep_save[:,:,it] = numrep
            prevlat_save[:,:,it] = prevlat
            prevsymp_save[:,:,it] = prevsymp
            prevasymp_save[:,:,it] = prevasymp
            prevpresymp_save[:,:,it] = prevpresymp
            prevrec_save[:,:,it] = prevrec
            newinf_save[:,:,it] = newinf
            atworkinf_save[:,:,it] = atworkinf
            atworkasymp_save[:,:,it] = atworkasymp
            newasymp_save[:,:,it] = newasymp
            num_CT_save[:,:,it] = num_CT
            if runset!="change_n_students"
                num_infected_save[:,:,it] = num_infected
            end
            social_dynamic_infection_count_save[:,:,it] = social_dynamic_infection_count
            accom_dynamic_infection_count_save[:,:,it] = accommodation_dynamic_infection_count
            household_infection_count_save[:,:,it] = household_infection_count
            var_num_infected_save[1,:,it] = var_num_infected
            num_init_infected_save[it] = num_init_infected
            mean_init_generation_time_save[1,:,it] = mean_init_generation_time
            Rt_save[:,:,it] = Rt
            new_rehoused_save[:,:,it] = new_rehoused
            current_rehoused_save[:,:,it] = current_rehoused
            tests_performed_save[:,:,it] = tests_performed
            test_outcomes_save[:,:,:,it] = test_outcomes

            # Initialise isolation outputs
            num_isolating_save[:,:,it] = num_isolating
            num_isolating_oncampus_save[:,:,it] = num_isolating_oncampus
            num_isolating_offcampus_save[:,:,it] = num_isolating_offcampus
            num_symp_isolating_save[:,:,it] = num_symp_isolating
            num_asymp_isolating_save[:,:,it] = num_asymp_isolating
            num_household_isolating_save[:,:,it] = num_household_isolating
            num_accom_lockdown_isolating_save[:,:,it] = num_accom_lockdown_isolating

            # Allocate infection 1D outputs to save variables
            n_oncampus_inf_save[:,it] = n_oncampus_inf
            n_offcampus_inf_save[:,it] = n_offcampus_inf

            # Allocate isolation 1D outputs to save variables
            n_isol_adhering_save[:,it] = n_isol_adhering
            n_isol_adhering_oncampus_save[:,it] = n_isol_adhering_oncampus
            n_isol_adhering_offcampus_save[:,it] = n_isol_adhering_offcampus

            # Allocate 3D outputs to save variables
            # If on first replicate, initialise the society infection count save variable
            if it == 1
                n_societies = length(network_parameters.society_info)
                society_infection_count_save = zeros(Int64,endtime+1,countfinal,n_societies,num_config)
            end
            cohort_infection_count_save[:,:,:,it] = cohort_infection_count
            society_infection_count_save[:,:,:,it] = society_infection_count

            # Outputs saved if flag condition satisfied
            if contact_tracing_active==true
                num_isolating_CTcause_save[:,:,it] = num_isolating_CTcause
            end

            if mass_testing_active == true
                # If mass testing in use, write values to file from the mass testing parameter structure
                mass_testing_n_tests_save[it] = deepcopy(mass_testing_parameters_alt.n_tests_performed)
                mass_testing_n_positive_save[it] = deepcopy(mass_testing_parameters_alt.n_tests_positive)
                mass_testing_n_all_isolating_save[it] = deepcopy(mass_testing_parameters_alt.n_all_isolations_caused)
                mass_testing_n_hh_isolating_save[it] = deepcopy(mass_testing_parameters_alt.n_hh_isolations_caused)
                mass_testing_n_CT_isolating_save[it] = deepcopy(mass_testing_parameters_alt.n_CT_isolations_caused)
                mass_testing_n_prev_infected_tested_save[it] = deepcopy(mass_testing_parameters_alt.n_prev_infected_tested)
            end
    end
    # Save outputs to file
    if contact_tracing_active==true
        output_file = "../Results/uni_model_infection_output_$(runset)_withCT_#$(job_ID).mat"
        file = matopen(output_file, "w")
        write(file, "num_CT", num_CT_save)
        write(file, "num_isolating_CTcause", num_isolating_CTcause_save)
    else
        output_file = "../Results/uni_model_infection_output_$(runset)_#$(job_ID).mat"
        file = matopen(output_file, "w")
    end
    write(file,"numlat",numlat_save)
    write(file,"numinf",numinf_save)
    write(file,"numrep",numrep_save)
    write(file,"prevlat",prevlat_save)
    write(file,"prevsymp",prevsymp_save)
    write(file,"prevasymp",prevasymp_save)
    write(file,"prevpresymp",prevpresymp_save)
    write(file,"prevrec",prevrec_save)
    write(file,"newinf",newinf_save)
    write(file,"atworkinf",atworkinf_save)
    write(file,"atworkasymp",atworkasymp_save)
    write(file,"newasymp",newasymp_save)
    write(file, "num_infected", num_infected_save)
    write(file, "social_dynamic_infection_count", social_dynamic_infection_count_save)
    write(file, "accom_dynamic_infection_count", accom_dynamic_infection_count_save)
    write(file, "household_infection_count", household_infection_count_save)
    write(file, "var_num_infected_save", var_num_infected_save)
    write(file, "num_init_infected_save", num_init_infected_save)
    write(file, "mean_init_generation_time_save", mean_init_generation_time_save)
    write(file, "Rt_save", Rt_save)
    write(file,"new_rehoused",new_rehoused_save)
    write(file,"current_rehoused",current_rehoused_save)
    write(file, "cohort_infection_count", cohort_infection_count_save)
    write(file, "society_infection_count", society_infection_count_save)
    write(file, "tests_performed", tests_performed_save)
    write(file, "test_outcomes", test_outcomes_save)

    # Initialise isolation outputs
    write(file, "num_isolating", num_isolating_save)
    write(file, "num_isolating_oncampus", num_isolating_oncampus_save)
    write(file, "num_isolating_offcampus", num_isolating_offcampus_save)
    write(file, "num_symp_isolating", num_symp_isolating_save)
    write(file, "num_asymp_isolating", num_asymp_isolating_save)
    write(file, "num_household_isolating", num_household_isolating_save)
    write(file, "num_accom_lockdown_isolating", num_accom_lockdown_isolating_save)

    # Allocate infection 1D outputs to save variables
    write(file,"n_oncampus_inf",n_oncampus_inf_save)
    write(file,"n_offcampus_inf",n_offcampus_inf_save)

    # Allocate isolation 1D outputs to save variables
    write(file,"n_isol_adhering",n_isol_adhering_save)
    write(file,"n_isol_adhering_oncampus",n_isol_adhering_oncampus_save)
    write(file,"n_isol_adhering_offcampus",n_isol_adhering_offcampus_save)

    if runset=="amount_backwards_CT"
        write(file, "prob_backwards_CT", prob_backwards_CT_config)
    end
    if mass_testing_active == true
        write(file,"mass_testing_n_tests_used",mass_testing_n_tests_save)
        write(file,"mass_testing_n_positive",mass_testing_n_positive_save)
        write(file,"mass_testing_n_all_isolating",mass_testing_n_all_isolating_save)
        write(file,"mass_testing_n_hh_isolating",mass_testing_n_hh_isolating_save)
        write(file,"mass_testing_n_CT_isolating",mass_testing_n_CT_isolating_save)
        write(file,"mass_testing_n_prev_infected_tested",mass_testing_n_prev_infected_tested_save)
    end
    close(file)
end
