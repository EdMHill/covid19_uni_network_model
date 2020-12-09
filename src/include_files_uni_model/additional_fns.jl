"""
Purpose:
Functions required to run university network model

The functions contained within this file include:
- populate_inclass!                 (arrays specifying when students are in class)
- populate_society_schedule!        (arrays specifying when societies are running events)
- increment_counters!               (Used each timestep to increase the value of time in state variables)
- load_configs                      (For running batches of scenarios, selected by string value)
- find_network_parameters           (Load relevant network params based on number of cohorts requested)

Transmission functions contained within this file include:
- transmit_over!                    (infection event check over list of class contacts, though not used for other cohort contacts)
- transmit_over_societies!          (infection event check over society contacts)

Functions to set up transmission rates within household for each individual
- assign_household_transmit_onegroup!  (use if everyone has the same household transmission risk)
- assign_household_transmit_halls_risk! (assign different transmission risk based on halls of residence. Off campus households get medium risk.)

Functions to set up transmission rates cohort and societies for each individual
- assign_cohort_transmit!
- assign_society_sports_transmit!
- assign_dynamic_social_transmit!

Functions to reinitialise states at start of each run
- reinitialise_node_states!         (multiply time series vectors by 0)
- reinitialise_daily_record_arrays! (reset arrays tracking status of attending certain settings, being in isolation etc)
- reinitialise_student_params!      (reinitialise the student parameter fields)
- reinitialise_class_params!         (reinitialise the class related quantities)
- reinitialise_society_params!       (reinitialise the society parameter type variables)
- reinitialise_CT_vars!             (reinitialise the contact tracing variables)

Miscellaneous functions contained within this file include:
- draw_sample_from_pmf!
- set_infection_related_times!      (set times to infection etc.: returns inftime, symptime, lattime, hh_isolation and delay_adherence)
- flattenall                        (Collapse nested structures into a single vector)
"""

function populate_inclass!(inclass::Array{Int64,2},
                            sameday::Int64,
                            ton::Int64,
                            toff::Int64,
                            rng::MersenneTwister)


# Initialise array to store indicator values of whether individuals are in class/at work or not each day
lmul!(0,inclass)

n_students,endtime = size(inclass)

    # If sameday=0, all those in the team type are at their work setting on the same set of consecutive days.
    # If sameday=1, all those in the team type are at their work setting on a random set of consecutive days.
    # If sameday=2, all those in the team type are at their work setting on the same number of days, but scattered randomly throughout the week.
    if sameday==0
        # iterate over each node
        for student_itr=1:n_students
            num_reps = floor(Int64,endtime/7) # how many weeks in the simulation
            # iterate over each week
            for reps_itr = 1:num_reps
                for days=ton:toff # all works work days ton to toff
                    inclass[student_itr,(reps_itr-1)*7+(days+1)] = 1
                end
            end
            # and put in the last bit, if there aren't an exact number of repetitions
            for days=ton:toff # all works work days ton to toff
                if (num_reps*7+(days+1))<=endtime
                    inclass[student_itr,num_reps*7+(days+1)] = 1
                end
            end
        end
    elseif sameday==1 # all those in the team type are at their work setting on a random set of consecutive days

        # initialise pap
        pap = zeros(7)
        # iterate over each node
        for student_itr=1:n_students
            num_reps = floor(Int64,endtime/7) # how many weeks in the simulation
            dayon = ceil(Int64,7*rand(rng)) # which day do they start work
            lmul!(0,pap) # reinitalise pap
            # stay at work for toff days
            for t_it = 0:toff

                # Get day number of week.
                # uses mod1, rather than mod, so mod1(7,7) returns 7 rather than 0
                day_of_week = mod1(dayon+t_it,7)

                # Give value to pap
                pap[day_of_week] = 1
            end

            # iterate over each week
            for reps_itr = 1:num_reps
                # and put in the days at work in pap
                for day_itr = 1:7
                    inclass[student_itr,(reps_itr-1)*7+day_itr] = pap[day_itr]
                end
            end

            # and put in the last bit, if there aren't an exact number of repetitions
            # and put in the days at work in pap
            for day_itr = 1:7
                if (num_reps*7+day_itr)<=endtime
                    inclass[student_itr,num_reps*7+day_itr] = pap[day_itr]
                end
            end
        end
    elseif sameday==2 # workdays randomly placed throughout the week (but repeat each week)
        # iterate over each node
        pap = zeros(Int64,7)
        for student_itr=1:n_students
            randperm!(rng,pap)  # Get permutation of 1:7
            num_reps = floor(Int64,endtime/7) # how many weeks in the simulation
            # iterate over each week
            for reps_itr = 1:num_reps
                # and put in the first toff days in the random permuation, pap
                for pap_itr = 1:(toff+1)
                    inclass[student_itr,(reps_itr-1)*7+pap[pap_itr]] = 1
                end
            end
            # and put in the last bit, if there aren't an exact number of repetitions
            # and put in the days at work in pap
            for pap_itr = 1:(toff+1)
                if (num_reps*7+pap[pap_itr])<=endtime
                    inclass[student_itr,num_reps*7+pap[pap_itr]] = 1
                end
            end
        end
    elseif sameday==3 # do ton weeks on followed by toff weeks off, all simultaneous
        # iterate over each node
        for student_itr=1:n_students
            num_reps = floor(Int64,endtime/7) # how many weeks in the simulation
            # iterate over each week
            for reps_itr = 1:num_reps
                if mod1(reps_itr,ton+toff) <= ton # using mod1, as mod1(1,1) = 1, whereas mod(1,1) = 0
                    for days=1:5 # only work weekdays
                        inclass[student_itr,(reps_itr-1)*7+days] = 1
                    end
                end
            end
        end
    elseif sameday==4 # do ton weeks on followed by toff weeks off, starting at a random week each
        # initialise pap
        pap = zeros(ton+toff)
        # iterate over each node
        for student_itr=1:n_students
            # which week will they start
            weekon = ceil(Int64,(ton+toff)*rand(rng))
            lmul!(0,pap) # reinitalise pap
            # stay at work for ton weeks
            for t_it = 1:ton
                pap[mod1(weekon+(t_it-1),ton+toff)] = 1
            end
            num_reps = floor(Int64,endtime/7) # how many weeks in the simulation
            # iterate over each week
            for reps_itr = 1:num_reps
                if pap[mod1(reps_itr,ton+toff)] == 1
                    for days=1:5 # only work weekdays
                        inclass[student_itr,(reps_itr-1)*7+days] = 1
                    end
                end
            end
        end
    end

    return nothing
end

function reinitialise_society_params!(society_info::Array{society_params,1},
                                    endtime::Int64,
                                    rng::MersenneTwister)

"""
Outline of steps, iterating over each society:
- Initialise the schedule array
- Assign the schedule array. Based on 3 days randomly scattered throughout week.
- reinitialise variables as required
"""


    # Get number of societies in the system
    n_societies = length(society_info)

    # Iterate over each society
    for society_itr = 1:n_societies

        # Initialise the schedule array
        lmul!(0,society_info[society_itr].schedule)

        # Assign the schedule array.
        # Based on 3 days randomly scattered throughout week.
        # iterate over each node
        pap = zeros(Int64,7)
        randperm!(rng,pap)  # Get permutation of 1:7
        num_reps = floor(Int64,endtime/7) # how many weeks in the simulation

        # iterate over each week
        for reps_itr = 1:num_reps
            # and put in the first toff days in the random permuation, pap
            for pap_itr = 1:3
                society_info[society_itr].schedule[(reps_itr-1)*7+pap[pap_itr]] = 1
            end
        end

        # and put in the last bit, if there aren't an exact number of repetitions
        # and put in the days at work in pap
        for pap_itr = 1:3
            if (num_reps*7+pap[pap_itr])<=endtime
                society_info[society_itr].schedule[num_reps*7+pap[pap_itr]] = 1
            end
        end

        # Reinitialise other fields in the society parameter structure type
        society_info[society_itr].f2f_activity = true
        society_info[society_itr].society_inactivation_time = 0
    end

    return nothing
end

function increment_counters!(states::student_states,
                            household_isoltime::Int64,
                            symp_isoltime::Int64,
                            asymp_isoltime::Int64,
                            contact_tracing_active::Bool;
                            timeisol_CTcause::Array{Int64,1}=zeros(Int64,1),
                            CT_caused_isol_limit::Int64=0)

@unpack timelat, timeinf, timesymp, timeisol, symp_timeisol, asymp_timeisol,
        lattime, inftime, symptime, asymp, timeisol, symp_timeisol = states

# Increments time counters
    n_students = length(timelat)
    for student_itr = 1:n_students
        if timelat[student_itr]>0
            timelat[student_itr] += 1
        end

        if timeinf[student_itr]>0
            timeinf[student_itr] += 1
        end

        if timesymp[student_itr]>0
            timesymp[student_itr] += 1
        end

        if timeisol[student_itr]>0
            timeisol[student_itr] += 1
        end

        if timeisol[student_itr]>household_isoltime
            timeisol[student_itr] = 0
        end

        # Isolation due to symptomatic infection
        if symp_timeisol[student_itr]>0
            symp_timeisol[student_itr] += 1
        end

        # Exit isolation at end of symptomatic infection
        # Also resets other isolation counters to zero
        if symp_timeisol[student_itr]>symp_isoltime
            symp_timeisol[student_itr] = 0
            timeisol[student_itr] = 0
            if contact_tracing_active == true
                timeisol_CTcause[student_itr] = 0
            end
        end

        # Isolation due to discovering asymptomatic infection
        if asymp_timeisol[student_itr]>0
            asymp_timeisol[student_itr] += 1
        end

        if asymp_timeisol[student_itr]>asymp_isoltime
            asymp_timeisol[student_itr] = 0
            timeisol[student_itr] = 0
            if contact_tracing_active == true
                timeisol_CTcause[student_itr] = 0
            end
        end

        if contact_tracing_active == true

            # Increment time in self-isolation, caused by Contact Tracing
            if timeisol_CTcause[student_itr]>0
                timeisol_CTcause[student_itr] += 1
            end

            # Reset contact tracing counter if limit is exceeded
            if timeisol_CTcause[student_itr]>CT_caused_isol_limit
                timeisol_CTcause[student_itr] = 0
            end
        end
    end
    return nothing
end


function load_configs(runset::String,
                        n_replicates::Int64,
                        n_cohorts::Int64,
                        n_students::Int64)

# first set default parameters (then overwrite those that need to be changed)
sameday = 3
ton = 1
toff = 0
attendance_propn = 1.

# classes close when 50% report infection
# only used if workplace_closure_active = true
class_CT_threshold = 0.5

# 10% of people correctly identify infector
# only used if perform_CT_from_infector = true
prob_backwards_CT = 0.
infector_engage_with_CT_prob = 1.0

# 70% of those adhering to self isolation engage with contact tracing
# only used if contact_tracing_active = true
CT_engagement = 0.7

# Set default adherence
adherence = 0.7

# probability of being asymptomatic
probasymp_dist = Uniform(0.5,0.8)  # 0.9

# Amend infection potential of asymptomatics
asymp_trans_scaling_dist::Uniform{Float64} = Uniform(0.3,0.7)

# change scaling of infection parameter for non-household contacts
scaling = 0.34

# Scale the infectiousness of all contacts
suscep_scaling = 0.8

# At beginning of simulation, the proportion of student population that has
# been infected previously
recov_propn = 0.1

# If sameday=0, all workers are at work on the same set of consecutive days.
# If sameday=1, workers go to work on a random set of consecutive days.
# If sameday=2, workers go to work on the same number of days, but scattered randomly throughout the week.
# If sameday=3, workers go to work ton weeks on followed by toff weeks off, all in synchrony
# If sameday=4, workers go to work ton weeks on followed by toff weeks off, beginning at a random time
if runset=="sweep_adherence"
    CT_engagement = 1.  # Assume those that are symptomatic will also engage with contact tracing
    adherence_config = [0:0.1:1;]
    n_configs = length(adherence_config)
elseif runset=="sweep_adherence_with_rehouse"
        CT_engagement = 1.  # Assume those that are symptomatic will also engage with contact tracing
        adherence_config = [0:0.1:1;]
        n_configs = length(adherence_config)

        # For this config, have that rehousing strategy IS in use
        rehouse_strat_config = repeat([true],n_configs)
elseif runset=="alter_test_result_delay"

    # Assume those that are symptomatic will report and also engage with contact tracing
    adherence = 1.
    CT_engagement = 1.

    # Array with row per config. pmf with column for delay length. Col 1 for 0 days, Col 2 for 1 days etc
    CT_delay_until_test_result_config = [1. 0. 0. 0. 0. 0.;
                                        0. 1. 0. 0. 0. 0.;
                                        0. 0. 1. 0. 0. 0.;
                                        0. 0. 0. 1. 0. 0.;
                                        0. 0. 0. 0. 1. 0.;
                                        0. 0. 0. 0. 0. 1.]
    n_configs::Int64 = size(CT_delay_until_test_result_config,1)
elseif runset=="trans_risk_scale_no_interv"
    # everyone works mon-fri
    sameday = 3
    ton = 1
    toff = 0
    attendance_propn = 1.0

    # No interventions in use
    adherence = 0.
    CT_engagement = 0.

    # Assume all nodes in system begin susceptible
    recov_propn = 0.

    # change the infection scaling
    scaling_config = [0.1:0.1:0.4;]
    n_configs = length(scaling_config)
elseif runset=="suscep_scale_no_interv"
    # everyone works mon-fri
    sameday = 3
    ton = 1
    toff = 0
    attendance_propn = 1.0

    # No interventions in use
    adherence = 0.
    CT_engagement = 0.

    # Assume all nodes in system begin susceptible
    recov_propn = 0.

    # change the infection scaling
    suscep_config = [0.1:0.1:1.0;]
    n_configs = length(suscep_config)
elseif runset=="probasymp_scale_no_interv"
    # everyone works mon-fri
    sameday = 3
    ton = 1
    toff = 0
    attendance_propn = 1.0

    # No interventions in use
    adherence = 0.
    CT_engagement = 0.

    # Assume all nodes in system begin susceptible
    recov_propn = 0.

    # change the infection scaling
    probasymp_config = [0.5:0.1:0.9;]
    n_configs = length(probasymp_config)
elseif runset=="transasymp_scale_no_interv"
    # everyone works mon-fri
    sameday = 3
    ton = 1
    toff = 0
    attendance_propn = 1.0

    # No interventions in use
    adherence = 0.
    CT_engagement = 0.

    # Assume all nodes in system begin susceptible
    recov_propn = 0.

    # change the infection scaling
    transasymp_config = collect(0.1:0.1:0.6) # Scaling for transmission potential of asymp
    n_configs = length(transasymp_config)
elseif runset=="amount_backwards_CT"
    prob_backwards_CT_config = [0:0.05:0.5;]
    n_configs = length(prob_backwards_CT_config)
elseif runset=="run_one_run"
    n_configs = 1

    # Assume those that are symtomatic will report and also engage with contact tracing
    adherence = 1.
    CT_engagement = 1.
# elseif runset=="run_no_interventions"
#     RNGseed_config = [100,200,300]
#     n_configs = length(RNGseed_config)
#
#     # Assume those that are symtomatic will report and also engage with contact tracing
#     adherence = 0.
#     CT_engagement = 0.
#     prob_backwards_CT = 0.
#
#     # Assume all nodes in system begin susceptible
#     recov_propn = 0.
elseif runset=="class_CT_threshold"
    work_or_study_group_CT_threshold_config = [0:0.05:0.5;]
    n_configs = length(work_or_study_group_CT_threshold_config)
elseif runset=="CT_engagement"
    CT_engagement_config = [0:0.1:1;]
    n_configs = length(CT_engagement_config)
elseif runset=="asymp_scen"
    # Set possible values
    probasymp_vals = collect(0.5:0.1:0.9)  # probability of being asymptomatic
    transasymp_vals = collect(0.1:0.1:1) # Scaling for transmission potential of asymp

    # Construct configurations
    n_configs = length(probasymp_vals)*length(transasymp_vals)
    probasymp_config = zeros(n_configs)
    transasymp_config = zeros(n_configs)
    for probasymp_itr = 1:length(probasymp_vals)
        for transasymp_itr = 1:length(transasymp_vals)
            config_idx = transasymp_itr + (probasymp_itr-1)*length(transasymp_vals)
            probasymp_config[config_idx] = probasymp_vals[probasymp_itr]
            transasymp_config[config_idx] = transasymp_vals[transasymp_itr]
        end
    end
elseif runset=="change_n_students"
    n_students_config = [10000:10000:50000;]
    n_configs = length(n_students_config)
elseif runset=="change_recov_propn"
    # Assume those that are symptomatic will report and also engage with contact tracing
    adherence = 1.
    CT_engagement = 1.

    # Set configuration
    recov_propn_config = [0., 0.05, 0.1, 0.15, 0.2]
    n_configs = length(recov_propn_config)
elseif runset=="rehouse_strat_check"
    # Assume those that are symptomatic will report and also engage with contact tracing
    adherence = 1.
    CT_engagement = 1.

    # Set configuration
    rehouse_strat_config = [false,true]
    n_configs = length(rehouse_strat_config)
elseif runset=="vary_accom_level_lockdown"

    # Assume those that are symptomatic will report and also engage with contact tracing
    adherence = 1.
    CT_engagement = 1.

    # Array with row per config.
    # List of intervention fns to be enacted
    # Have to pass dummy function as the main run fn accepts an array of functions
    intervention_fn_config::Array{Function,1} = [accom_lockdown_household_level!;
                                                  accom_lockdown_floor_level!;
                                                  accom_lockdown_block_level!;
                                                  accom_lockdown_hall_level!
                                                   ]
    n_configs = size(intervention_fn_config,1)

    # For this config, have that rehousing strategy is not in use
    rehouse_strat_config = repeat([false],n_configs)
elseif runset=="run_one_run_baseline"
    n_configs = 1

    # Assume those that are symptomatic will report and also engage with contact tracing
    adherence = 1.
    CT_engagement = 1.

    # For this config, have that rehousing strategy is not in use
    rehouse_strat_config = repeat([false],n_configs)
elseif runset=="run_one_run_baseline_nointerv"
    n_configs = 1

    # Assume those that are symtomatic will report and also engage with contact tracing
    adherence = 0.
    CT_engagement = 0.
    prob_backwards_CT = 0.

    # Assume all nodes in system begin susceptible
    recov_propn = 0.

    # For this config, have that rehousing strategy is not in use
    rehouse_strat_config = repeat([false],n_configs)
elseif runset == "vary_household_allocation_strat"
    # Assume those that are symptomatic will report and also engage with contact tracing
    adherence = 1.
    CT_engagement = 1.


    # Array with row per configuration
    assign_household_fn_config::Array{String,1} = ["assign_households_all_students";
                                                        "assign_households_by_cohort"
                                                        ]

    n_configs = size(assign_household_fn_config,1)
elseif runset == "run_mass_testing"

    # Assume those that are symptomatic will report and also engage with contact tracing
    adherence = 0.7
    CT_engagement = 1.

    # Specify the possible values for the mass testing parameters
    mass_test_desig_days = [14,35,56]
    on_campus_cov = [1.,0.5,0.]
    off_campus_cov = [1.,0.5,0.]

    # Get number of configs in use.
    # Subtract 1 as will not be using the on_campus_cov, off_campus_cov both zero option
    n_configs = length(mass_test_desig_days)*((length(on_campus_cov)*length(off_campus_cov)) - 1)

    # Initialise mass testing configuration variables
    mass_testing_config = Array{mass_testing_params,1}(undef,n_configs)

    # For this config, have that rehousing strategy IS in use
    rehouse_strat_config = repeat([true],n_configs)

    # Initialise mass testing parameter fields
    n_tests_performed_vecs = Array{Array{Int64,1},1}(undef,n_replicates)
    n_tests_positive_vecs = Array{Array{Int64,1},1}(undef,n_replicates)
    n_all_isolations_caused_vecs = Array{Array{Int64,1},1}(undef,n_replicates)
    n_hh_isolations_caused_vecs = Array{Array{Int64,1},1}(undef,n_replicates)
    n_CT_isolations_caused_vecs = Array{Array{Int64,1},1}(undef,n_replicates)
    n_prev_infected_tested_vecs = Array{Array{Int64,1},1}(undef,n_replicates)

    config_itr = 1 # Initialise allocation index
    for test_day_itr = 1:length(mass_test_desig_days)
        selected_test_day = mass_test_desig_days[test_day_itr]
        for on_campus_cov_itr = 1:length(on_campus_cov)
            selected_on_campus_cov = on_campus_cov[on_campus_cov_itr]
            for off_campus_cov_itr = 1:length(off_campus_cov)
                selected_off_campus_cov = off_campus_cov[off_campus_cov_itr]
                if (selected_on_campus_cov != 0.) || (selected_off_campus_cov != 0.)

                    # Initialise mass testing parameter fields for current configuration
                    for replicate_itr = 1:n_replicates
                        n_tests_performed_vecs[replicate_itr] = [0]
                        n_tests_positive_vecs[replicate_itr] = [0]
                        n_all_isolations_caused_vecs[replicate_itr] = [0]
                        n_hh_isolations_caused_vecs[replicate_itr] = [0]
                        n_CT_isolations_caused_vecs[replicate_itr] = [0]
                        n_prev_infected_tested_vecs[replicate_itr] = [0]
                    end

                    # Only set up configuration if a location has non-zero testing coverage
                    mass_testing_config[config_itr] = mass_testing_params(designated_test_times = [selected_test_day],
                                                                            on_campus_coverage_propn = [selected_on_campus_cov],
                                                                            off_campus_coverage_propn = [selected_off_campus_cov],
                                                                            n_tests_performed = deepcopy(n_tests_performed_vecs),
                                                                            n_tests_positive = deepcopy(n_tests_positive_vecs),
                                                                            n_all_isolations_caused = deepcopy(n_all_isolations_caused_vecs),
                                                                            n_hh_isolations_caused = deepcopy(n_hh_isolations_caused_vecs),
                                                                            n_CT_isolations_caused = deepcopy(n_CT_isolations_caused_vecs),
                                                                            n_prev_infected_tested = deepcopy(n_prev_infected_tested_vecs))

                    # Increment config_itr
                    config_itr += 1
                end
            end
        end
    end
elseif runset == "run_mass_testing_6configs"

    # Assume those that are symptomatic will report and also engage with contact tracing
    adherence = 0.7
    CT_engagement = 1.

    # Specify the possible values for the mass testing parameters
    mass_test_desig_days = [14,35]
    on_campus_cov = [1., 0.]
    off_campus_cov = [1., 0.]

    # Get number of configs in use.
    # Subtract 1 as will not be using the on_campus_cov, off_campus_cov both zero option
    n_configs = length(mass_test_desig_days)*((length(on_campus_cov)*length(off_campus_cov)) - 1)

    # Initialise mass testing configuration variables
    mass_testing_config = Array{mass_testing_params,1}(undef,n_configs)

    # For this config, have that rehousing strategy IS in use
    rehouse_strat_config = repeat([true],n_configs)

    # Initialise mass testing parameter fields
    n_tests_performed_vecs = Array{Array{Int64,1},1}(undef,n_replicates)
    n_tests_positive_vecs = Array{Array{Int64,1},1}(undef,n_replicates)
    n_all_isolations_caused_vecs = Array{Array{Int64,1},1}(undef,n_replicates)
    n_hh_isolations_caused_vecs = Array{Array{Int64,1},1}(undef,n_replicates)
    n_CT_isolations_caused_vecs = Array{Array{Int64,1},1}(undef,n_replicates)
    n_prev_infected_tested_vecs = Array{Array{Int64,1},1}(undef,n_replicates)

    config_itr = 1 # Initialise allocation index
    for test_day_itr = 1:length(mass_test_desig_days)
        selected_test_day = mass_test_desig_days[test_day_itr]
        for on_campus_cov_itr = 1:length(on_campus_cov)
            selected_on_campus_cov = on_campus_cov[on_campus_cov_itr]
            for off_campus_cov_itr = 1:length(off_campus_cov)
                selected_off_campus_cov = off_campus_cov[off_campus_cov_itr]
                if (selected_on_campus_cov != 0.) || (selected_off_campus_cov != 0.)

                    # Initialise mass testing parameter fields for current configuration
                    for replicate_itr = 1:n_replicates
                        n_tests_performed_vecs[replicate_itr] = [0]
                        n_tests_positive_vecs[replicate_itr] = [0]
                        n_all_isolations_caused_vecs[replicate_itr] = [0]
                        n_hh_isolations_caused_vecs[replicate_itr] = [0]
                        n_CT_isolations_caused_vecs[replicate_itr] = [0]
                        n_prev_infected_tested_vecs[replicate_itr] = [0]
                    end

                    # Only set up configuration if a location has non-zero testing coverage
                    mass_testing_config[config_itr] = mass_testing_params(designated_test_times = [selected_test_day],
                                                                            on_campus_coverage_propn = [selected_on_campus_cov],
                                                                            off_campus_coverage_propn = [selected_off_campus_cov],
                                                                            n_tests_performed = deepcopy(n_tests_performed_vecs),
                                                                            n_tests_positive = deepcopy(n_tests_positive_vecs),
                                                                            n_all_isolations_caused = deepcopy(n_all_isolations_caused_vecs),
                                                                            n_hh_isolations_caused = deepcopy(n_hh_isolations_caused_vecs),
                                                                            n_CT_isolations_caused = deepcopy(n_CT_isolations_caused_vecs),
                                                                            n_prev_infected_tested = deepcopy(n_prev_infected_tested_vecs))

                    # Increment config_itr
                    config_itr += 1
                end
            end
        end
    end
elseif runset == "run_mass_testing_alter_adherence"

    # Assume those that are symptomatic will report and also engage with contact tracing
    adherence_vals = [0.2,0.5,0.8]
    CT_engagement = 1.

    # Specify the possible values for the mass testing parameters
    mass_test_desig_days = [21,42,63]
    on_campus_cov = [1., 0.]
    off_campus_cov = [1., 0.]

    # Get number of configs in use.
    # Subtract 1 as will not be using the on_campus_cov, off_campus_cov both zero option
    n_configs_per_adherence_val = length(mass_test_desig_days)*((length(on_campus_cov)*length(off_campus_cov)) - 1)
    n_configs = n_configs_per_adherence_val*length(adherence_vals)

    # Initiliase adherence_config
    adherence_config = zeros(n_configs)

    # Initialise mass testing configuration variables
    mass_testing_config = Array{mass_testing_params,1}(undef,n_configs)

    # For this config, have that rehousing strategy IS in use
    rehouse_strat_config = repeat([true],n_configs)

    # Initialise mass testing parameter fields
    n_tests_performed_vecs = Array{Array{Int64,1},1}(undef,n_replicates)
    n_tests_positive_vecs = Array{Array{Int64,1},1}(undef,n_replicates)
    n_all_isolations_caused_vecs = Array{Array{Int64,1},1}(undef,n_replicates)
    n_hh_isolations_caused_vecs = Array{Array{Int64,1},1}(undef,n_replicates)
    n_CT_isolations_caused_vecs = Array{Array{Int64,1},1}(undef,n_replicates)
    n_prev_infected_tested_vecs = Array{Array{Int64,1},1}(undef,n_replicates)

    config_itr = 1 # Initialise allocation index
    for adherence_itr = 1:length(adherence_vals)
        start_adherence_idx = ((adherence_itr-1)*n_configs_per_adherence_val) + 1
        end_adherence_idx = adherence_itr*n_configs_per_adherence_val
        adherence_config[start_adherence_idx:end_adherence_idx] .= adherence_vals[adherence_itr]
        for test_day_itr = 1:length(mass_test_desig_days)
            selected_test_days = [mass_test_desig_days[test_day_itr]] # Vector to match input type of Array{Int64,1}
            n_test_days = length(selected_test_days)
            for on_campus_cov_itr = 1:length(on_campus_cov)
                selected_on_campus_cov = on_campus_cov[on_campus_cov_itr]
                for off_campus_cov_itr = 1:length(off_campus_cov)
                    selected_off_campus_cov = off_campus_cov[off_campus_cov_itr]
                    if (selected_on_campus_cov != 0.) || (selected_off_campus_cov != 0.)
                        # Initialise mass testing parameter fields for current configuration
                        for replicate_itr = 1:n_replicates
                            n_tests_performed_vecs[replicate_itr] = zeros(Int64,n_test_days)
                            n_tests_positive_vecs[replicate_itr] = zeros(Int64,n_test_days)
                            n_all_isolations_caused_vecs[replicate_itr] = zeros(Int64,n_test_days)
                            n_hh_isolations_caused_vecs[replicate_itr] = zeros(Int64,n_test_days)
                            n_CT_isolations_caused_vecs[replicate_itr] = zeros(Int64,n_test_days)
                            n_prev_infected_tested_vecs[replicate_itr] = zeros(Int64,n_test_days)
                        end

                        # Only set up configuration if a location has non-zero testing coverage
                        mass_testing_config[config_itr] = mass_testing_params(designated_test_times = selected_test_days,
                                                                                on_campus_coverage_propn = selected_on_campus_cov*ones(Int64,n_test_days),
                                                                                off_campus_coverage_propn = selected_off_campus_cov*ones(Int64,n_test_days),
                                                                                n_mass_tests_performed = zeros(Int64,n_replicates),
                                                                                n_tests_performed = deepcopy(n_tests_performed_vecs),
                                                                                n_tests_positive = deepcopy(n_tests_positive_vecs),
                                                                                n_all_isolations_caused = deepcopy(n_all_isolations_caused_vecs),
                                                                                n_hh_isolations_caused = deepcopy(n_hh_isolations_caused_vecs),
                                                                                n_CT_isolations_caused = deepcopy(n_CT_isolations_caused_vecs),
                                                                                n_prev_infected_tested = deepcopy(n_prev_infected_tested_vecs))
                        # Increment config_itr
                        config_itr += 1
                    end
                end
            end
        end
    end
elseif runset == "run_mass_testing_frequent"

    # Scenario: Fortnightly mass testing or weekly mass testing

    # Assume those that are symptomatic will report and also engage with contact tracing
    adherence_vals = [0.2,0.5,0.8]
    CT_engagement = 1.

    # Specify the possible values for the mass testing parameters
    on_campus_cov = [1., 0.]
    off_campus_cov = [1., 0.]

    # Get number of configs in use.
    # Subtract 1 as will not be using the on_campus_cov, off_campus_cov both zero option
    # Multiply by 2 to account for fornightly and weekly test freqency scenarios
    n_configs_per_adherence_val = 2*((length(on_campus_cov)*length(off_campus_cov)) - 1)
    n_configs = n_configs_per_adherence_val*length(adherence_vals)

    # Initiliase adherence_config
    adherence_config = zeros(n_configs)

    # Initialise mass testing configuration variables
    mass_testing_config = Array{mass_testing_params,1}(undef,n_configs)

    # For this config, have that rehousing strategy IS in use
    rehouse_strat_config = repeat([true],n_configs)

    # Initialise mass testing parameter fields
    n_tests_performed_vecs = Array{Array{Int64,1},1}(undef,n_replicates)
    n_tests_positive_vecs = Array{Array{Int64,1},1}(undef,n_replicates)
    n_all_isolations_caused_vecs = Array{Array{Int64,1},1}(undef,n_replicates)
    n_hh_isolations_caused_vecs = Array{Array{Int64,1},1}(undef,n_replicates)
    n_CT_isolations_caused_vecs = Array{Array{Int64,1},1}(undef,n_replicates)
    n_prev_infected_tested_vecs = Array{Array{Int64,1},1}(undef,n_replicates)

    config_itr = 1 # Initialise allocation index
    for adherence_itr = 1:length(adherence_vals)
        start_adherence_idx = ((adherence_itr-1)*n_configs_per_adherence_val) + 1
        end_adherence_idx = adherence_itr*n_configs_per_adherence_val
        adherence_config[start_adherence_idx:end_adherence_idx] .= adherence_vals[adherence_itr]
        for test_day_itr = 1:2
            if test_day_itr == 1
                # Fortnightly testing
                selected_test_days = [1,14,28,42,56,70]
            else
                selected_test_days = [1,7,14,21,28,35,42,49,56,63,70]
            end
            n_test_days = length(selected_test_days)

            for on_campus_cov_itr = 1:length(on_campus_cov)
                selected_on_campus_cov = on_campus_cov[on_campus_cov_itr]
                for off_campus_cov_itr = 1:length(off_campus_cov)
                    selected_off_campus_cov = off_campus_cov[off_campus_cov_itr]
                    if (selected_on_campus_cov != 0.) || (selected_off_campus_cov != 0.)

                        # Initialise mass testing parameter fields for current configuration
                        for replicate_itr = 1:n_replicates
                            n_tests_performed_vecs[replicate_itr] = zeros(Int64,n_test_days)
                            n_tests_positive_vecs[replicate_itr] = zeros(Int64,n_test_days)
                            n_all_isolations_caused_vecs[replicate_itr] = zeros(Int64,n_test_days)
                            n_hh_isolations_caused_vecs[replicate_itr] = zeros(Int64,n_test_days)
                            n_CT_isolations_caused_vecs[replicate_itr] = zeros(Int64,n_test_days)
                            n_prev_infected_tested_vecs[replicate_itr] = zeros(Int64,n_test_days)
                        end

                        # Only set up configuration if a location has non-zero testing coverage
                        mass_testing_config[config_itr] = mass_testing_params(designated_test_times = selected_test_days,
                                                                                on_campus_coverage_propn = selected_on_campus_cov*ones(Int64,n_test_days),
                                                                                off_campus_coverage_propn = selected_off_campus_cov*ones(Int64,n_test_days),
                                                                                n_mass_tests_performed = zeros(Int64,n_replicates),
                                                                                n_tests_performed = deepcopy(n_tests_performed_vecs),
                                                                                n_tests_positive = deepcopy(n_tests_positive_vecs),
                                                                                n_all_isolations_caused = deepcopy(n_all_isolations_caused_vecs),
                                                                                n_hh_isolations_caused = deepcopy(n_hh_isolations_caused_vecs),
                                                                                n_CT_isolations_caused = deepcopy(n_CT_isolations_caused_vecs),
                                                                                n_prev_infected_tested = deepcopy(n_prev_infected_tested_vecs))

                        # Increment config_itr
                        config_itr += 1
                    end
                end
            end
        end
    end
elseif runset=="run_baseline_mass_testing"

    # Assume those that are symptomatic will report and also engage with contact tracing
    adherence_config= [0.2,0.5,0.8]
    CT_engagement = 1.

    # Set number of configurations in use
    n_configs = length(adherence_config)

    # For this config, have that rehousing strategy is not in use
    rehouse_strat_config = repeat([true],n_configs)
end

# extend scalars to vectors where needed
if @isdefined(adherence_config)==false
    adherence_config = adherence*ones(Int64,n_configs)
end
if @isdefined(sameday_config)==false
    sameday_config = sameday*ones(Int64,n_configs)
end
if @isdefined(ton_config)==false
    ton_config = ton*ones(Int64,n_configs)
end
if @isdefined(toff_config)==false
    toff_config = toff*ones(Int64,n_configs)
end
if @isdefined(work_percent_config)==false
    n_departments = n_cohorts ÷ 3
    work_percent_config = attendance_propn*ones(n_configs,n_departments,3)
end
if @isdefined(work_or_study_group_CT_threshold_config)==false
    work_or_study_group_CT_threshold_config = class_CT_threshold*ones(n_configs)
end
if @isdefined(prob_backwards_CT_config)==false
    prob_backwards_CT_config = prob_backwards_CT*ones(n_configs)
end
if @isdefined(infector_engage_with_CT_prob_config)==false
    infector_engage_with_CT_prob_config = infector_engage_with_CT_prob*ones(n_configs)
end
if @isdefined(CT_engagement_config)==false
    CT_engagement_config = CT_engagement*ones(n_configs)
end
if @isdefined(probasymp_config)==false
    # probasymp_config = probasymp*ones(n_configs)
    probasymp_config = repeat([probasymp_dist],n_configs)
end
if @isdefined(transasymp_config)==false
    # transasymp_config = asymp_trans_scaling*ones(n_configs)
    transasymp_config = repeat([asymp_trans_scaling_dist],n_configs)
end
if @isdefined(scaling_config)==false
    scaling_config = scaling*ones(n_configs)
end
if @isdefined(suscep_config)==false
    suscep_config = suscep_scaling*ones(n_configs)
end
if @isdefined(n_students_config)==false
    n_students_config = n_students*ones(Int64,n_configs)
end
if @isdefined(recov_propn_config)==false
    recov_propn_config = recov_propn*ones(Int64,n_configs)
end
if @isdefined(rehouse_strat_config)==false
    rehouse_strat_config = repeat([false],n_configs)
end

# # Set up the RNG seed
# if @isdefined(RNGseed_config)==false
#     RNGseed_config = RNGseed*ones(Int64,n_configs)
# end

# Set up the contact tracing delay array
if @isdefined(CT_delay_until_test_result_config)==false
    # Set default of 2 day delay for test result. Repeat for number of configs being tested
    CT_delay_until_test_result_pmf = [0. 0. 1. 0. 0. 0.]
    CT_delay_until_test_result_config = repeat(CT_delay_until_test_result_pmf,n_configs,1)
end

# Set up the intervention fn vector
if @isdefined(intervention_fn_config)==false
    # Set default of no additional intervention fns being in use
    intervention_fn_vec::Array{Function,1} = [dummy_example!]
    intervention_fn_config = repeat(intervention_fn_vec,n_configs)
end

# Set up the household assignment fn
if @isdefined(assign_household_fn_config)==false
    assign_household_fn_vec::Array{String,1} = ["assign_households_all_students"]
    assign_household_fn_config = repeat(assign_household_fn_vec,n_configs)
end

# Set up mass testing fn
if @isdefined(mass_testing_config)==false
    mass_testing_vec = [mass_testing_params()]
    mass_testing_config = repeat(mass_testing_vec,n_configs)
end
    # Return configurations
    return adherence_config, sameday_config, ton_config, toff_config, work_percent_config,
    n_configs, work_or_study_group_CT_threshold_config, prob_backwards_CT_config,
    infector_engage_with_CT_prob_config, CT_engagement_config, probasymp_config,
    transasymp_config, scaling_config, suscep_config,
    n_students_config, recov_propn_config,
    rehouse_strat_config,
    CT_delay_until_test_result_config::Array{Float64,2},
    intervention_fn_config, assign_household_fn_config,
    mass_testing_config::Array{mass_testing_params,1}
end

function find_network_parameters(n_cohorts::Int64;
                                    attendence_propns::Array{Float64,2}=Array{Float64,2}[],
                                    n_students::Int64 = 0,
                                    n_students_on_campus::Int64 = 0)

    if n_cohorts==12
        # Faculty and undergrad/postgrad breakdown
        # # Arts (UG/PG), Social sciences (UG/PG), Science & Engineering (UG/PG), Medicine (UG/PG).

        if isempty(attendence_propns)
            attendence_propns = ones(4,3)  # proportion of those in each each cohort that will attend f2f classes if they are running
        end

        """
        Define network associated parameters
        """
        # Set up class_degree_distribution
        # Distinct values for undergrad & postgrad
        class_degree_distribution_array = repeat([Distributions.LogNormal(1.749,1.331) Distributions.LogNormal(1.749,1.331) Distributions.LogNormal(1.223,1.125)],
                                                            4, 1)

        # Set up dynamic_social_contact_degree_distribution
        # Distinct values for undergrad & postgrad
        dynamic_social_contact_degree_distribution_array = repeat([Distributions.LogNormal(1.646,1.211) Distributions.LogNormal(1.646,1.211) Distributions.LogNormal(1.590,1.128)],
                                                            4, 1)

        # Assign network parameters
        network_parameters = network_params(n_students = n_students,
            n_students_on_campus = n_students_on_campus,
            prob_worktype_contact = ones(4,3).*(10000/n_students),    # Scale relative to level used for 10,000 nodes
            class_degree_distribution = class_degree_distribution_array,  # If single entry, will be replicated across all cohorts
            between_class_contact_probs = [1.0],  # If single entry, will be replicated across all cohorts
            dynamic_social_contact_degree_distribution = dynamic_social_contact_degree_distribution_array)

        """
        Define study associated parameters
        """
        # Set up class size means and stan. dev.
        class_size_mean_array = repeat([25. 25. 5.], 4, 1)
        class_size_sd_array = repeat([0. 0. 1.], 4, 1)

        # Assign class parameters
        class_generation_parameters = class_generation_params(n_cohorts=12,
                    attendence_propns=attendence_propns,
                    classtype_proportion = [0.1054*0.3 0.1054*0.7 0.0186;
                                                0.2408*0.3 0.2408*0.7 0.2051;
                                                0.2401*0.3 0.2401*0.7 0.1327;
                                                0.0309*0.3 0.0309*0.7 0.0264],
                    class_size_mean = class_size_mean_array,
                    class_size_sd= class_size_sd_array)
    elseif n_cohorts == 84
        # 28 departments

        if isempty(attendence_propns)
            # Two columns. Undergrad/postgrad distinction
            attendence_propns = ones(28,3)  # proportion of those in each each cohort that will attend f2f classes if they are running
        end

        """
        Define network associated parameters
        """
        # Set up class_degree_distribution
        # Distinct values for undergrad & postgrad
        class_degree_distribution_array = repeat([Distributions.LogNormal(1.749,1.331) Distributions.LogNormal(1.749,1.331) Distributions.LogNormal(1.223,1.125)],
                                                            28, 1)

        # Set up dynamic_social_contact_degree_distribution
        # Distinct values for undergrad & postgrad
        dynamic_social_contact_degree_distribution_array = repeat([Distributions.LogNormal(1.646,1.211) Distributions.LogNormal(1.646,1.211) Distributions.LogNormal(1.590,1.128)],
                                                            28, 1)

        # Assign network parameters
        network_parameters = network_params(n_students = n_students,
            n_students_on_campus = n_students_on_campus,
            prob_worktype_contact = ones(28,3).*(10000/n_students),    # Scale relative to level used for 10,000 nodes
            class_degree_distribution = class_degree_distribution_array,  # If single entry, will be replicated across all cohorts
            between_class_contact_probs = [1.0],  # If single entry, will be replicated across all cohorts
            dynamic_social_contact_degree_distribution = dynamic_social_contact_degree_distribution_array)

        """
        Define study associated parameters
        """
        # Set up class size means and stan. dev.
        class_size_mean_array = repeat([25. 25. 5.], 28, 1)
        class_size_sd_array = repeat([0. 0. 1.], 28, 1)

        # Assign class parameters
        class_generation_parameters = class_generation_params(n_cohorts=84,
                attendence_propns=attendence_propns,
                classtype_proportion = zeros(28,3),
                class_size_mean = class_size_mean_array,
                class_size_sd= class_size_sd_array)
    end

    return network_parameters, class_generation_parameters
end

"""
Transmission related fns
"""

function transmit_over!(student_info::Array{student_params,1},
                    transmission_risk::Float64,
                    infected_by::Array{Int64,1},
                    output::sim_outputs,
                    states::student_states,
                    probasymp::Float64,
                    rng::MersenneTwister,
                    time::Int64,
                    count::Int64;
                    infecting_by::Int64=0,
                    contacts_to_check::Array{Int64,1}=Array{Int64,1}(),
                    dynamic_contact::Int64=0,
                    inisol::Array{Int64,2} = Array{Int64,2}(undef,0,0),
                    attendance_record::Array{Int64,2} = Array{Int64,2}(undef,0,0),
                    cohort_ID::Int64 = 0,
                    household_check_flag::Bool = false)
    # Inputs:
    # student_info - Fields with individual level information
    # transmission_risk::Float64 - Probability of transmission given contact made
    # infected_by::Array{Int64,1} - Record of ID of infectors of each node
    # output::sim_outputs - record of outputs from the simulations
    # states::student_states - status of each node
    # probasymp::Float64 - Asymptomatic probability
    # rng::MersenneTwister - The random number generator
    # time::Int64 - Current timestep
    # count::Int64 - Replicate ID
    # infecting_by::Int64 - Node ID of node transmitting infection
    # contacts_to_check::Int64 - list of nodes IDs to receive infection
    # dynamic_contact::Int64 - Flag for if this infection is over a dynamic contact link
    # inisol::Array{Int64,2} - Flag if nodes are in isolation
    # attendance_record::Array{Int64,2} - Flag if nodes are at work
    # cohort_ID::Int64 - If dealing with a study based transmission event, the relevant cohort identifier
    # household_check_flag::Bool - If true, then checking transmission with household contacts. Otherwise, has value false.

    for contact_itr = 1:length(contacts_to_check) # Iterate over each contact
        infecting_to = contacts_to_check[contact_itr]

        # Check if infecting_to is susceptible, not isolating and at work
        # only checks isolation or at work if given the arguments inisol / attendance_record
        if (states.timelat[infecting_to]==0) &&
            (isassigned(inisol) == false || inisol[time,infecting_to] == false) &&
            (isassigned(attendance_record) == false || attendance_record[time,infecting_to] == true)

            if rand(rng) < transmission_risk
                states.timelat[infecting_to] = 1
                infected_by[infecting_to] = infecting_by
                output.num_infected[infecting_by,count] += 1
                states.acquired_infection[infecting_to] = time

                # Update cohort infection event counter, if applicable
                # Will have retained value 0 if not in use
                if cohort_ID > 0
                    output.cohort_infection_count[time+1,count,cohort_ID] += 1
                        # Row 1 corresponds to timestep 0. So need to save to row "time + 1" for saving at value "time"
                elseif household_check_flag == true
                    output.household_infection_count[time+1,count] += 1
                        # Row 1 corresponds to timestep 0. So need to save to row "time + 1" for saving at value "time"
                end

                # adjust Rt(t) = mean number of infections generated by nodes that were infected at time t
                if states.acquired_infection[infecting_by]>0
                    # Offset time to array indexing. Row 1 of output.Rt is for day 0, Row 2 is for day 1 etc
                    output.Rt[(states.acquired_infection[infecting_by]+1),count] += 1
                end

                # if this was from an initial infection, modify the generation time
                # this sum will be divided by the total number of secondary initial infections
                if infected_by[infecting_by]==-1
                    output.mean_init_generation_time[count] += time + states.lattime[infecting_by]
                end

                # Check if infection will be asymptomatic
                if rand(rng) < probasymp
                    states.asymp[infecting_to] = 1
                end

                # Update latent event counter
                output.numlat[time+1,count] += 1

                # Update infection location counters
                if (student_info[infecting_to].household_info.on_campus_accom == true)
                    output.n_oncampus_inf[count] += 1
                elseif (student_info[infecting_to].household_info.on_campus_accom == false)
                    output.n_offcampus_inf[count] += 1
                else
                    error("Transmit infection check. on_campus_accom status invalid")
                end

                # if this was a dynamic contact, update counter
                if dynamic_contact==1 # Social dynamic contact
                    output.social_dynamic_infection_count[time+1,count] += 1
                elseif dynamic_contact==2 # accommodation dynamic contact
                    output.accommodation_dynamic_infection_count[time+1,count] += 1
                end
            end
        end
    end
end

function transmit_over_societies!(student_info::Array{student_params,1},
                                    society_info::Array{society_params,1},
                                    infection_parameters::infection_params,
                                    infectiousness::Float64,
                                    asymp_trans_scaling::Float64,
                                    infected_by::Array{Int64,1},
                                    output::sim_outputs,
                                    states::student_states,
                                    probasymp::Float64,
                                    rng::MersenneTwister,
                                    time::Int64,
                                    count::Int64;
                                    infecting_by::Int64=0,
                                    society_contacts::Array{Array{Int64,1},2} = Array{Array{Int64,1},2}(undef,0,0),
                                    inisol::Array{Int64,2} = Array{Int64,2}(undef,0,0),
                                    society_attendance_record::Array{Int64,3} = Array{Int64,3}(undef,0,0,0))
    # Inputs:
    # student_info - Fields with individual level information
    # society_info::society_params - Parameter structure, entry for each society.
    # infection_parameters::infection_params - Structure with fields relating to transmission related variables
    # infectiousness::Float64 - relative infectivity on given day of infection
    # asymp_trans_scaling::Float64 - The relative infectiousness of an asymptomatic (sampled earlier in main_function.jl for the simn replicate)
    # infected_by::Array{Int64,1} - Record of ID of infectors of each node
    # output::sim_outputs - record of outputs from the simulations
    # states::student_states - status of each node
    # probasymp::Float64 - Asymptomatic probability
    # rng::MersenneTwister - The random number generator
    # time::Int64 - Current timestep
    # count::Int64 - Replicate ID
    # infecting_by::Int64 - Node ID of node transmitting infection
    # society_contacts::Int64 - list of nodes IDs contacted within societies
    # inisol::Array{Int64,2} - Flag if nodes are in isolation
    # society_attendance_record::Array{Int64,2} - Flag if nodes partcipate in society events

    @unpack suscep_scaling = infection_parameters

    if states.asymp[infecting_by]>0 # Asymptomatic
        base_transmission = infectiousness*suscep_scaling*asymp_trans_scaling
    else
        base_transmission = infectiousness*suscep_scaling
    end

    # Get the relevant society IDs for the person of interest
    infecting_student = student_info[infecting_by]
    society_ID_list = infecting_student.society_IDs

    # Iterate over each society.
    for society_itr = 1:length(society_ID_list)
        current_society_ID = society_ID_list[society_itr]
        #  Update attendance status if society is
        # active (i.e. running face-to-face events)
        if society_attendance_record[time,current_society_ID,infecting_by] == true

            # Get transmission risk based on if group is a society or sports club
            # society_type_val: 1 for society, 2 for sports club.
            society_type_val = society_info[current_society_ID].society_type
            transmission_risk = base_transmission*infecting_student.transrisk_society_sports[society_type_val]

            # Attended society event on current timestep
            # From all possible contacts, check which contacts are made
            contacts_to_check = society_contacts[infecting_by,current_society_ID]

            for contact_itr = 1:length(contacts_to_check) # Iterate over each contact
                infecting_to = contacts_to_check[contact_itr]

                # Check if infecting_to is susceptible, not isolating and attended the society
                # only checks isolation or at work if given the arguments inisol / attendance_record
                if (states.timelat[infecting_to]==0) &&
                    (inisol[time,infecting_to] == false) &&
                    (society_attendance_record[time,current_society_ID,infecting_to] == true)

                    if rand(rng) < transmission_risk
                        states.timelat[infecting_to] = 1
                        infected_by[infecting_to] = infecting_by
                        output.num_infected[infecting_by,count] += 1
                        states.acquired_infection[infecting_to] = time

                        # Update society infection event counter
                        output.society_infection_count[time+1,count,current_society_ID] += 1
                            # Row 1 corresponds to timestep 0. So need to save to row "time + 1" for saving at value "time"

                        # adjust Rt(t) = mean number of infections generated by nodes that were infected at time t
                        if states.acquired_infection[infecting_by]>0
                            # Offset time to array indexing. Row 1 of output.Rt is for day 0, Row 2 is for day 1 etc
                            output.Rt[(states.acquired_infection[infecting_by]+1),count] += 1
                        end

                        # if this was from an initial infection, modify the generation time
                        # this sum will be divided by the total number of secondary initial infections
                        if infected_by[infecting_by]==-1
                            output.mean_init_generation_time[count] += time + states.lattime[infecting_by]
                        end

                        # Check if infection will be asymptomatic
                        if rand(rng) < probasymp
                            states.asymp[infecting_to] = 1
                        end

                        # Update latent event counter
                        output.numlat[time+1,count] += 1

                        # Update infection location counters
                        if (student_info[infecting_to].household_info.on_campus_accom == true)
                            output.n_oncampus_inf[count] += 1
                        elseif (student_info[infecting_to].household_info.on_campus_accom == false)
                            output.n_offcampus_inf[count] += 1
                        else
                            error("Transmit infection check. on_campus_accom status invalid")
                        end
                    end
                end
            end
        end
    end
end

"""
Functions to set up transmission rates within household for each individual
"""
# Single household risk type
function assign_household_transmit_onegroup!(RNGseed::Int64,
                                    network_parameters::network_params,
                                    household_contacts_per_node::Array{Int64,1},
                                    transrisk_household_group_mean::Array{Float64,1},
                                    transrisk_household_group_sd::Array{Float64,1})
# Inputs:
# RNGseed::Int64 - Set the number to seed the random number generator
# network_parameters::network_params - Quantities to construct the contacts & stores the node properties
# household_contacts_per_node::Array{Int64,1} - As described
# transrisk_household_group_mean/sd::Array{Float64,1} -  probability of transmission within a household, can differ per household group

    # Unpack parameters
    @unpack n_students, student_info = network_parameters

    # Set the random number generator
    rng = MersenneTwister(RNGseed)

    # Get number of household groups in use
    n_transrisk_household_group_mean = length(transrisk_household_group_mean)
    n_transrisk_household_group_sd = length(transrisk_household_group_sd)

    # Throw error if more than one group
    if (n_transrisk_household_group_mean != 1)
        error("Should only be a single household group SAR mean estimate, but have found there to be $(n_transrisk_household_group_mean). Please rectify.")
    end

    # Throw error if more than one group
    if (n_transrisk_household_group_sd != 1)
        error("Should only be a single household group SAR standard deviation, but have found there to be $(n_transrisk_household_group_sd). Please rectify.")
    end

    # Construct normal distribution to sample from
    mean_val = transrisk_household_group_mean[1]
    sd_val = transrisk_household_group_sd[1]
    norm_dist = Normal(mean_val,sd_val)

    # Iterate over each individual
    for node_itr = 1:n_students
        student_info[node_itr].transrisk_household = rand(rng,norm_dist)
    end

    return nothing
end

# Household risk based on household size
function assign_household_transmit_household_size!(RNGseed::Int64,
                                    network_parameters::network_params,
                                    household_contacts_per_node::Array{Int64,1},
                                    transrisk_household_group_mean::Array{Float64,1},
                                    transrisk_household_group_sd::Array{Float64,1})
# Inputs:
# RNGseed::Int64 - Set the number to seed the random number generator
# network_parameters::network_params - Quantities to construct the contacts & stores the node properties
# household_contacts_per_node::Array{Int64,1} - As described
# transrisk_household_group_mean/sd::Array{Float64,1} -  probability of transmission within a household, can differ per household group

    # Unpack parameters
    @unpack n_students, student_info = network_parameters

    # Set the random number generator
    rng = MersenneTwister(RNGseed)

    # Get number of household groups in use
    n_transrisk_household_group_mean = length(transrisk_household_group_mean)
    n_transrisk_household_group_sd = length(transrisk_household_group_sd)

    # Throw error if there is not four groups
    if (n_transrisk_household_group_mean != 4)
        error("Should be four group SAR mean estimates, but have found there to be $(n_transrisk_household_group_mean). Please rectify.")
    end

    # Throw error if more than one group
    if (n_transrisk_household_group_sd != 4)
        error("Should be four group SAR standard deviations, but have found there to be $(n_transrisk_household_group_sd). Please rectify.")
    end

    # Construct normal distribution to sample from
    norm_dists = Normal.(transrisk_household_group_mean,transrisk_household_group_sd)

    # Iterate over each individual
    for node_itr = 1:n_students

        # Check number of people in the household
        hh_size = household_contacts_per_node[node_itr] + 1
        if hh_size == 1
            # Sole person in household. No household contacts, set transmission risk to zero.
            student_info[node_itr].transrisk_household = 0
        elseif hh_size == 2
            # Household size two
            student_info[node_itr].transrisk_household = rand(rng,norm_dists[1])
        elseif hh_size == 3
            # Household size three
            student_info[node_itr].transrisk_household = rand(rng,norm_dists[2])
        elseif hh_size == 4
            # Household size four
            student_info[node_itr].transrisk_household = rand(rng,norm_dists[3])
        else
            # Household size five or more
            student_info[node_itr].transrisk_household = rand(rng,norm_dists[4])
        end
    end

    return nothing
end

# For on-campus accommodation, assign household risk based on halls
# For off-campus accommodation, assign the "medium" risk
function assign_household_transmit_halls_risk!(RNGseed::Int64,
                                    network_parameters::network_params,
                                    household_contacts_per_node::Array{Int64,1},
                                    transrisk_household_group_mean::Array{Float64,1},
                                    transrisk_household_group_sd::Array{Float64,1})
# Inputs:
# RNGseed::Int64 - Set the number to seed the random number generator
# network_parameters::network_params - Quantities to construct the contacts & stores the node properties
# household_contacts_per_node::Array{Int64,1} - As described
# transrisk_household_group_mean/sd::Array{Float64,1} -  probability of transmission within a household, can differ per household group

    # Unpack parameters
    @unpack n_students, student_info = network_parameters

    # Set the random number generator
    rng = MersenneTwister(RNGseed)

    # Get number of household groups in use
    # For this function, should be three groups. Risk status of "Low","Medium","High"
    n_transrisk_household_group_mean = length(transrisk_household_group_mean)
    n_transrisk_household_group_sd = length(transrisk_household_group_sd)

    # Throw error if more than one group
    if (n_transrisk_household_group_mean != 3)
        error("Should be three household group SAR mean estimate, but have found there to be $(n_transrisk_household_group_mean). Please rectify.")
    end

    # Throw error if more than one group
    if (n_transrisk_household_group_sd != 3)
        error("Should be three household group SAR standard deviation, but have found there to be $(n_transrisk_household_group_sd). Please rectify.")
    end


    # Construct normal distribution to sample from
    norm_dists = Normal.(transrisk_household_group_mean,transrisk_household_group_sd)

    # Iterate over each individual
    for node_itr = 1:n_students

        # Check student on campus accommodation status
        on_campus_accom_status = student_info[node_itr].household_info.on_campus_accom
        if  on_campus_accom_status == false
             # If living off-campus, assign the medium risk
             student_info[node_itr].transrisk_household = rand(rng,norm_dists[2])
         else
             # If living on-campus, check the hall ID
             hall_ID_for_student = student_info[node_itr].household_info.hall_ID

             # Assign to the appropriate ID based on that status.
             # SPECIFIC HALL IDs TO BE ADDED!
             if hall_ID_for_student <= 4
                 # Low risk halls
                 student_info[node_itr].transrisk_household = rand(rng,norm_dists[1])
             elseif hall_ID_for_student <= 8
                 # Medium risk halls
                 student_info[node_itr].transrisk_household = rand(rng,norm_dists[2])
             else
                 # High risk halls
                 student_info[node_itr].transrisk_household = rand(rng,norm_dists[3])
             end
         end
    end

    return nothing
end

function assign_household_transmit_multigrouptest!(RNGseed::Int64,
                                    network_parameters::network_params,
                                    household_contacts_per_node::Array{Int64,1},
                                    transrisk_household_group_mean::Array{Float64,1},
                                    transrisk_household_group_sd::Array{Float64,1})
# Inputs:
# RNGseed::Int64 - Set the number to seed the random number generator
# network_parameters::network_params - Quantities to construct the contacts & stores the node properties
# household_contacts_per_node::Array{Int64,1} - As described
# transrisk_household_group_mean/sd::Array{Float64,1} -  probability of transmission within a household, can differ per household group

    # Unpack parameters
    @unpack n_students, student_info = network_parameters

    # Set the random number generator
    rng = MersenneTwister(RNGseed)

    # Get number of household groups in use
    n_transrisk_household_group_mean = length(transrisk_household_group_mean)
    n_transrisk_household_group_sd = length(transrisk_household_group_sd)

    # Throw error if more than one group
    if (n_transrisk_household_group_mean < 1)
        error("Should be multiple household group SAR mean estimate, but have found there to be $(n_transrisk_household_group_mean). Please rectify.")
    end

    # Throw error if more than one group
    if (n_transrisk_household_group_sd < 1)
        error("Should be multiple household group SAR standard deviation, but have found there to be $(n_transrisk_household_group_sd). Please rectify.")
    end


    # Construct normal distribution to sample from
    norm_dists = Normal.(transrisk_household_group_mean,transrisk_household_group_sd)

    # Iterate over each individual
    for node_itr = 1:n_students
        if household_contacts_per_node[node_itr] == 0
            student_info[node_itr].transrisk_household = rand(rng,norm_dists[1])
        elseif household_contacts_per_node[node_itr] == 1
            student_info[node_itr].transrisk_household = rand(rng,norm_dists[2])
        else
            student_info[node_itr].transrisk_household = rand(rng,norm_dists[3])
        end
    end

    return nothing
end


"""
Functions to set up transmission rates within cohorts and societies for each individual
"""
# Cohort transmission risk
function assign_cohort_transmit!(RNGseed::Int64,
                                    network_parameters::network_params,
                                    class_generation_parameters::class_generation_params,
                                    transrisk_cohort_mean::Array{Float64,1},
                                    transrisk_cohort_sd::Array{Float64,1})
# Inputs:
# RNGseed::Int64 - Set the number to seed the random number generator
# network_parameters::network_params - Quantities to construct the contacts & stores the node properties
# class_generation_parameters::class_generation_params - As described
# transrisk_cohort_mean/sd::Array{Float64,1} -  probability of transmission with a cohort contact. Can differ by cohort.

    # Unpack parameters
    @unpack n_students, student_info = network_parameters
    @unpack n_cohorts = class_generation_parameters

    # Set the random number generator
    rng = MersenneTwister(RNGseed)

    # Get number of cohorts in use
    n_transrisk_mean = length(transrisk_cohort_mean)
    n_transrisk_sd = length(transrisk_cohort_sd)

    # Throw error if dimensions don't match
    if (n_transrisk_mean != n_cohorts) ||
         (n_transrisk_sd != n_cohorts)
         error("Dimension mismatch in transrisk_cohort_mean/sd parameter arrays. Please rectify.")
    end

    # Iterate over each individual
     for node_itr = 1:n_students
         student_cohort_ID = student_info[node_itr].cohort_ID
         mean_val = transrisk_cohort_mean[student_cohort_ID]
         sd_val = transrisk_cohort_sd[student_cohort_ID]
         student_info[node_itr].transrisk_cohort = rand(rng,Normal(mean_val,sd_val))
     end

    return nothing
end


# Society & sports group transmission risk
function assign_society_sports_transmit!(RNGseed::Int64,
                                    network_parameters::network_params,
                                    society_sports_transrisk_mean::Array{Float64,1},
                                    society_sports_transrisk_sd::Array{Float64,1})
# Inputs:
# RNGseed::Int64 - Set the number to seed the random number generator
# network_parameters::network_params - Quantities to construct the contacts & stores the node properties
# class_generation_parameters::class_generation_params - As described
# society_sports_transrisk_mean/sd::Array{Float64,1} -  probability of transmission with a contact in [society,sports club].

    # Unpack parameters
    @unpack n_students, student_info = network_parameters

    # Set the random number generator
    rng = MersenneTwister(RNGseed)

    # Iterate over each individual. Assign a transmission risk for societies and
    # sports clubs
     for node_itr = 1:n_students
         for society_sports_itr = 1:length(society_sports_transrisk_mean)
             mean_val = society_sports_transrisk_mean[society_sports_itr]
             sd_val = society_sports_transrisk_sd[society_sports_itr]
             student_info[node_itr].transrisk_society_sports[society_sports_itr] = rand(rng,Normal(mean_val,sd_val))
         end
     end

    return nothing
end

function assign_dynamic_social_transmit!(RNGseed::Int64,
                                    network_parameters::network_params,
                                    transrisk_dynamic_social_mean::Float64,
                                    transrisk_dynamic_social_sd::Float64)
# Inputs:
# RNGseed::Int64 - Set the number to seed the random number generator
# network_parameters::network_params - Quantities to construct the contacts & stores the node properties
# class_generation_parameters::class_generation_params - As described
# transrisk_dynamic_social_mean/sd::Array{Float64,1} -  transmission risk across a dynamic social contact.

    # Unpack parameters
    @unpack n_students, student_info = network_parameters

    # Set the random number generator
    rng = MersenneTwister(RNGseed)

     # Get number of household groups in use
     n_transrisk_dynamic_social_mean = length(transrisk_dynamic_social_mean)
     n_transrisk_dynamic_social_sd = length(transrisk_dynamic_social_sd)

     # Throw error if more than one group
     if (n_transrisk_dynamic_social_mean != 1)
         error("Should only be a single dynamic social transrisk mean estimate, but have found there to be $(n_transrisk_dynamic_social_mean). Please rectify.")
     end

     # Throw error if more than one group
     if (n_transrisk_dynamic_social_sd != 1)
         error("Should only be a single dynamic social transrisk sd estimate, but have found there to be $(n_transrisk_dynamic_social_sd). Please rectify.")
     end

     # Construct normal distribution to sample from
     dynamic_social_norm_dist = Normal(transrisk_dynamic_social_mean,transrisk_dynamic_social_sd)

     # Iterate over each individual
     for node_itr = 1:n_students
         student_info[node_itr].transrisk_dynamic_social = rand(rng,dynamic_social_norm_dist)
     end

    return nothing
end

"""
Functions to reinitialise states at start of each run
"""
# Node states, household inf delay & CT vars
function reinitialise_node_states!(states::student_states)
    lmul!(0,states.timelat)
    lmul!(0,states.timeinf)
    lmul!(0,states.timesymp)
    lmul!(0,states.asymp)
    lmul!(0,states.lattime)
    lmul!(0,states.timeisol)
    lmul!(0,states.symp_timeisol)
    lmul!(0,states.asymp_timeisol)
    lmul!(0,states.timeisol_CTcause)
    lmul!(0,states.hh_isolation)
    lmul!(0,states.delay_adherence)
    lmul!(0,states.acquired_infection)
end

# Reinitialise daily record arrays
function reinitialise_daily_record_arrays!(contacts::contacts_struct)
    lmul!(0,contacts.daily_record_inclass)
    lmul!(0,contacts.daily_record_inisol)
    lmul!(0,contacts.daily_record_atsociety)
end

function reinitialise_student_params!(n_students::Int64,
                                            student_info::Array{student_params,1})
# Inputs:
# n_students - Number of students in the system
# nodes - Parameter type to have information on each indivdiual

    # Iterate over each individual and set delays to default values
    for student_itr = 1:n_students
        student_info[student_itr].time_of_reporting_infection = 0
        student_info[student_itr].no_contacts_status = false

        # Remove any household lockdown
        student_info[student_itr].household_info.lockdown_status = false
    end



    return nothing
end

function reinitialise_class_params!(class_info::Array{Array{class_params,1},1})

    # Get number of team types in use
    n_team_types = length(class_info)

    # Iterate over each team type.
    # Within each, iterate over each student class/workforce group
    # Update field values
    for team_type_itr = 1:n_team_types
        n_teams = length(class_info[team_type_itr])
        for team_itr = 1:n_teams
            class_info[team_type_itr][team_itr].f2f_activity = true
            class_info[team_type_itr][team_itr].class_inactivation_time = 0
        end
    end

    return nothing
end

function reinitialise_CT_vars!(CT_vars::contact_tracing_vars,n_students::Int64, rng::MersenneTwister,
    CT_parameters::CT_params, delay_adherence::Array{Int64,1},
    csum_test_result_delay::Array{Float64,1},max_test_result_delay::Int64)

@unpack CT_days_before_symptom_included, CT_engagement = CT_parameters

    # Reset vector tracking symptomatic cases (positive confirmed or untested)
    lmul!(0,CT_vars.Symp_cases_per_household_pos_or_unknown)

    # Variables for waiting for test results
    lmul!(0,CT_vars.Time_to_test_result)
    CT_vars.Time_to_test_result .-= 1 # Reset so all values are -1

    # Repopulate Boolean vector stating whether a false negative test result would be returned
    # and the number of days relevant for contact tracing
    lmul!(0,CT_vars.relevant_prev_days_for_CT)
    for ii = 1:n_students

      # For each worker, initialise CT_vars.Test_result_false_negative as false
      CT_vars.Test_result_false_negative[ii] = false

      # For each worker, check if they engage with contact tracing
      engage_with_CT_rand = rand(rng)
      if engage_with_CT_rand < CT_engagement # engage with contact tracing
          CT_vars.Engage_with_CT[ii] = true
      else # do not engage with contact tracing
          CT_vars.Engage_with_CT[ii] = false
      end

      # Get amount of days to be looked back over
      # Have upper bound of 7 days post symp
      # if we put in reporting delay, needs to be above this
      CT_vars.relevant_prev_days_for_CT[ii] = min(CT_days_before_symptom_included + delay_adherence[ii],
                                          CT_days_before_symptom_included + 7)
    end

    # Repopulate time until test result received for each individual
    lmul!(0,CT_vars.CT_delay_until_test_result)
    for student_itr = 1:n_students
                      CT_vars.CT_delay_until_test_result[student_itr] = draw_sample_from_pmf(csum_test_result_delay,
                                                                                rng;
                                                                                idx_offset = 1)
    end

    # Set up vector of vectors for storing IDs of those to be contacted in CT
    CT_vars.Inds_to_be_contacted = Array{Array{Int64,1},1}(undef,n_students)

    # Initialise array to keep track of whether an infected recalls their infector
    lmul!(0,CT_vars.Recall_infector)
end

"""
Misc. fns
"""

function draw_sample_from_pmf(csum_pmf::Array{Float64,1},
                                rng::MersenneTwister;
                                idx_offset::Int64 = 0)
# Inputs:
# val_to_update::Int64 - Entry sampled value will be assigned to
# csum_pmf::Array{Float64,1} - Cumulative summed probability mass function. Used to draw value from.
# rng::MersenneTwister - The random number generator
# idx_offset::Int64 = 0 - Links bin index to the quantity value

    # Get number of elements in the pmf
    n_bins = length(csum_pmf)

    # Initialise output value
    val_to_update = 0

    # Draw random number
    # Set delay in adherence/symptoms becoming known to household
    # Find interval random number resides in
    r = rand(rng)
    allocated_flag = false # Intialise allocation flag. Switch to true when allocation done.
    bin_idx = 1   # Current interval being checked
    while (allocated_flag == false)
        if r <= csum_pmf[bin_idx]
            # Assign selected value
            # Subtract idx_offset
            val_to_update = bin_idx - idx_offset

            # Update allocation flag
            allocated_flag = true
        else
            # r does not reside in this interval. Update bin index.
            bin_idx += 1

            # Error check, if not assigned value after checked final bin value
            if bin_idx > n_bins
                error("bin_idx is now $bin_idx. The pmf only has $n_bins bins. Terminating programme.")
            end
        end
    end

    return val_to_update::Int64
end

function set_infection_related_times!(time_to_symps::Array{Int64,1},states::student_states,
    isolation::Int64,adherence::Float64,csum_delay_adherence::Array{Float64,1},
    d_incub::Distribution,n_students::Int64,rng::MersenneTwister)

    time_to_symps .= ceil.(rand(rng,d_incub,n_students)) # time to symptoms
    # (for asymptomatics, the same from a silent start of "symptoms")

    ## uncomment below to change default times
    # states.inftime = 2 # inftime is the time from infectiousness to symptoms (the same for everyone)
    # states.symptime = 7 # symptime is the time from symptoms to not infectious

    # iterate over nodes to set lattime and hh_isolation
    for student_itr = 1:n_students
        # lattime is the time from infection to infectiousness
        if time_to_symps[student_itr]-states.inftime<1
            states.lattime[student_itr] = 1            # Infectiousness can begin the day after becoming infected
        else
            states.lattime[student_itr] = time_to_symps[student_itr] - states.inftime
        end


        if isolation==1
            p1 = rand(rng)
            if p1 < adherence # those who adhere will isolate when they get symptoms
                states.hh_isolation[student_itr] = 1 # adherence to household isolation = 1 if adherent, 0 if not.
            end

            # Draw random number
            # Set delay in adherence/symptoms becoming known to household
            # Find interval random number resides in
            states.delay_adherence[student_itr] = draw_sample_from_pmf(csum_delay_adherence,
                                                                    rng;
                                                                    idx_offset = 1)
        end

    end
end


# Collapse nested structures into a single vector
function flattenall(a::AbstractArray)
    while any(x->typeof(x)<:AbstractArray, a)
        a = collect(Iterators.flatten(a))
    end
    return a
end
