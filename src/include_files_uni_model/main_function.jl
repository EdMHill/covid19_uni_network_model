"""
Main function
"""

# Run outbreak on network
function uni_network_run(RNGseed::Int64,
                                        n_students::Int64,
                                        ton::Int64,
                                        toff::Int64,
                                        infection_parameters::infection_params,
                                        sameday::Int64,
                                        seed_initial_states_fn::Function,
                                        countfinal::Int64,
                                        endtime::Int64,
                                        contact_tracing_active::Bool,
                                        CT_parameters::CT_params,
                                        network_parameters::network_params,
                                        class_generation_parameters::class_generation_params,
                                        society_generation_parameters::society_generation_params,
                                        work_study_group_closure_active::Bool,
                                        mass_testing_active::Bool,
                                        rehouse_strat_active::Bool;
                                        mass_testing_parameters::mass_testing_params = mass_testing_params(),
                                        intervention_fns::Array{Function,1} = Array{Function,1}(undef,0),
                                        generate_classes_fn::Function = generate_classes_default,
                                        generate_student_households_fn::String = "assign_households_no_hierarchy",
                                        assign_household_transrisk_fn::Function = assign_household_transmit_onegroup!,
                                        assign_societies_fn::Function = assign_societies_one_per_student,
                                        assign_cohort_transrisk_fn::Function = assign_cohort_transmit!,
                                        assign_society_sports_transrisk_fn::Function = assign_society_sports_transmit!,
                                        assign_dynamic_social_transrisk_fn::Function = assign_dynamic_social_transmit!
                                        )

# Inputs:
# RNGseed - Sets the random number generator
# n_students::Int64 - Number of students
# ton::Int64, toff::Int64 - Work pattern vars. Days on and days off, for sameday==3 ton weeks on followed by toff weeks off
# infection_parameters::infection_params - Parameter structure for infection params
# sameday::Int64 - Flag.  If sameday=0, all workers are at work on the same set of consecutive days. If sameday=1, workers go to work on a random set of consecutive days. If sameday=3, workers go to work on the same number of days, but scattered randomly throughout the week.
# initial_states_fn::Function - Sets amount of nodes to be seeded as latent/asymp/symp/recovered.
# countfinal::Int64, endtime::Int64 - % Simulation replicates & timesteps for each replicate
# endtime::Int64 - Number of timesteps each simulation replicate to be performed for
# contact_tracing_active::Bool - Set if contact tracing is active or not (Bool type variable)
# CT_parameters::CT_params - Parameter structure for contact tracing params
# network_parameters::network_params - Parameter structure for network params
# class_generation_parameters::class_generation_params - Parameter structure for class group generation params
# work_study_group_closure_active::Bool - Whether in person group interactions can be made inactive
# mass_testing_active::Bool - Whether large scale testing is planned to be used
# rehouse_strat_active::Bool - Intervention where those who are symptomatic are rehoused/completely isolate with no contacts.
# mass_testing_parameters::mass_testing_params - Specify information relevant to mass testing strategy
# intervention_fns - Specify use of any additional, trigger interventions
# generate_classes_fn - Specify function to allocate students to department/cohort & classes
# generate_student_households_fn! - Specify function to assign individuals to households and construct household contact network
# assign_household_transrisk_fn - Specify assignment of individuals to household groups with differing household transmission risk
# assign_societies_fn - Specify how students will be assigned to societies
# assign_cohort_transrisk_fn::Function - Specify how transmission risk from each inidividual to cohort contacts will be assigned
# assign_society_sports_transrisk_fn::Function - Specify how transmission risk from each inidividual to society & sports club contacts will be assigned
# assign_dynamic_social_transrisk_fn::Function - Specify how transmission risk from each inidividual to dynamic social contacts will be assigned

# For parameters within the parameter structures, see "include_files_network_model/parametertypes.jl"

##  OUTLINE OF THE CODE STRUCTURE
#         - Unpack required variables
#         - Set the random number generator
#         - Initialise students into cohorts, split by department and into first year undergrad, non-first year undergrad, postgrad (supporting functions in network_generation_fns.jl)
#         - Assign society/sports club membership  (supporting functions in network_generation_fns.jl)
#         - Generate study/cohort contacts, society contacts. (supporting functions in network_generation_fns.jl)
#         - Assign households and household contacts (supporting functions in network_generation_fns.jl)
#         - Generate social dynamic contacts (supporting functions in network_generation_fns.jl)
#         - Generate on-campus accomodation dynamic contacts (supporting functions in network_generation_fns.jl)
#         - Generate household specific transmission (supporting functions in network_generation_fns.jl)
#         - Initialise variables
#         - Iterate over replicates
#             - Reinitialisation phase (supporting functions in additional_fns.jl)
#             - Set course of infection times (supporting functions in additional_fns.jl)
#             - Update output time series with initial conditions
#             - Reset contact tracing variables (supporting functions in additional_fns.jl)
#             - Iterate over time
#                 - Reinitialise variables at start of timestep
#                 - Assign outputs
#                 - Increment counters (supporting functions in additional_fns.jl)
#                 - Increment infection process (if come to the end of latent time, move to infectious time etc.). Includes household infection loop. (supporting functions in additional_fns.jl)
#                 - Record whether individuals are in isolation. Set study and society attendance status for timestep.
#                 - Transmit infections (supporting functions in additional_fns.jl)
#                 - Perform contact tracing (supporting functions in contact_tracing_fns.jl)
#                 - Assign prevalence & isolation outputs
#                 - Reactive inactivation of in person classes/society meet ups
#                 - Run interventions (supporting functions inintervention_condition_affect_fns.jl)
#                 - Perform mass testing, if applicable (supporting functions in mass_testing_fns.jl)
#         - Return outputs

##

     """
     Unpack required variables
     """
    if contact_tracing_active==true
        @unpack CT_engagement, CT_delay_until_test_result_pmf, CT_days_before_symptom_included, test_false_negative_vec,
            CT_caused_isol_limit, dynamic_contacts_recalled_propn, social_contacts_recalled_propn, prob_backwards_CT,
            perform_CT_from_infector, infector_engage_with_CT_prob = CT_parameters
    end
    if work_study_group_closure_active==true
        @unpack work_or_study_group_CT_memory, work_or_study_group_CT_threshold, time_WC = CT_parameters
    end
    @unpack iso_trans_scaling, asymp_trans_scaling_dist,
            transrisk_household_group_mean, transrisk_household_group_sd,
            transrisk_cohort_mean, transrisk_cohort_sd,
            society_sports_transrisk_mean, society_sports_transrisk_sd,
            transrisk_dynamic_social_mean, transrisk_dynamic_social_sd,
            CS_scale_transrisk, suscep_scaling,
            probasymp_dist, isolation, symp_isoltime, asymp_isoltime, household_isoltime, adherence,
            n_social_mean_workday,n_social_mean_nonworkday,
            d_latent, dist_infectivity, delay_adherence_pmf,
            delay_household_infection_pmf,
            recov_propn = infection_parameters
    @unpack dynamic_conts_mean, dynamic_conts_sd, max_contacts_social_dynamic = network_parameters
    @unpack n_cohorts = class_generation_parameters
    @unpack society_types = society_generation_parameters

    """
    Set the random number generator
    """
    rng = MersenneTwister(RNGseed)

    """
    Initialise study group
    """

    # generate_classes in network_generation_fns.jl
    student_info::Array{student_params,1},
    class_sizes::Array{Array{Int64,1},1},
    class_info::Array{Array{class_params,1},1},
    nodes_by_class::Array{Array{Array{Int64,1},1},1} = generate_classes_fn(n_students, class_generation_parameters,RNGseed)
    network_parameters.student_info = student_info
    network_parameters.class_sizes = class_sizes
    network_parameters.class_info = class_info

    # Initialise sector open/closure status
    cohort_f2f_study_active = Array{Bool,1}(undef,n_cohorts)
    for sector_itr = 1:n_cohorts
        cohort_f2f_study_active[sector_itr] = true
    end
    network_parameters.cohort_f2f_study_active = cohort_f2f_study_active

    """
    Set up society membership
    """
    society_info::Array{society_params,1} = assign_societies_fn(student_info,
                                                                society_generation_parameters,
                                                                endtime,
                                                                RNGseed)

    # Assign society details to network parameter object
    network_parameters.society_info = society_info

    # Assign number of societies to variable
    n_societies = length(society_info)

    # Initialise society type activity status
    society_f2f_active = Array{Bool,1}(undef,society_types)
    for society_type_itr = 1:society_types
        society_f2f_active[society_type_itr] = true
    end
    network_parameters.society_f2f_active = society_f2f_active

    """
    Generate contacts within work settings, societies & households
    """
    # Generate contact network in work settings & societies.
    # Find relevant function in "include_files_uni_model/network_generation_fns.jl"
        contacts::contacts_struct,
        work_or_study_group_contacts_per_node::Array{Int64,1},
        cohort_contacts_per_node::Array{Int64,1},
        society_contacts_per_node::Array{Int64,2} = @time generate_contacts_uni(n_students,
                                                                            endtime,
                                                                            network_parameters,
                                                                            class_generation_parameters,
                                                                            nodes_by_class,
                                                                            RNGseed)

    # Assign households and generate contact network in that setting
    # Also modifies contacts.household_contacts
    # Find relevant function in "include_files_uni_model/network_generation_fns.jl"
    household_contacts_per_node::Array{Int64,1},
    n_households::Int64 = @time assign_households_fn!(generate_student_households_fn,
                                                                     n_students,
                                                                     contacts,
                                                                     network_parameters,
                                                                     RNGseed)

    # Update daily_record_atsociety in contacts, to use number of societies
    contacts.daily_record_atsociety = zeros(Int64,endtime,length(society_info),n_students)

    """
    Generate social dynamic contacts
    """
    # Set up dynamic worker contacts, in network_generation_fns.jl
    if network_parameters.network_generation_method == "ER"
        dynamic_social_contacts = generate_dynamic_student_contacts(RNGseed,
                                                            n_students,
                                                            endtime,
                                                            student_info,
                                                            dynamic_conts_mean,
                                                            dynamic_conts_sd)
    elseif network_parameters.network_generation_method == "configuration"
        dynamic_social_contacts = @time generate_dynamic_student_contacts(RNGseed,
                                                            n_students,
                                                            endtime,
                                                            student_info,
                                                            network_parameters.dynamic_social_contact_degree_distribution,
                                                            max_contacts_social_dynamic)
    else
        println("Invalid network generation method!")
    end

    contacts.dynamic_social_contacts = dynamic_social_contacts

    """
    Generate on-campus accomodation dynamic contacts
    """

    dynamic_accomodation_contacts = @time generate_dynamic_campus_accom_contacts(RNGseed,
                                                                                    n_students,
                                                                                    endtime,
                                                                                    student_info,
                                                                                    network_parameters,
                                                                                    contacts)

    # Assign to contacts parameter structure
    contacts.dynamic_accomodation_contacts = dynamic_accomodation_contacts


    """
    Generate transmission risks
    """
    # Relevant functions are listed in "include_files_network_model/additional_fns.jl"
    assign_household_transrisk_fn(RNGseed,
                                        network_parameters,
                                        household_contacts_per_node,
                                        transrisk_household_group_mean,
                                        transrisk_household_group_sd)

    assign_cohort_transrisk_fn(RNGseed,
                            network_parameters,
                            class_generation_parameters,
                            transrisk_cohort_mean,
                            transrisk_cohort_sd)

    assign_society_sports_transrisk_fn(RNGseed,
                                    network_parameters,
                                    society_sports_transrisk_mean,
                                    society_sports_transrisk_sd)

    assign_dynamic_social_transrisk_fn(RNGseed,
                                    network_parameters,
                                    transrisk_dynamic_social_mean,
                                    transrisk_dynamic_social_sd)

    """
    Initialise storage arrays
    """

    # Initialise output structure
    output = sim_outputs(endtime=endtime,countfinal=countfinal,n_students=n_students,
                            cohort_infection_count = zeros(Int64,endtime+1,countfinal,n_cohorts),
                            society_infection_count = zeros(Int64,endtime+1,countfinal,n_societies))

    # Initialise vectors used in each replicate
    dayon = zeros(Int64,n_students)
    ccount = zeros(Int64,n_students)
    pflag = zeros(Int64,n_students)

    # Initialise node states
    states = student_states(n_students=n_students)
    time_to_symps = zeros(Int64,n_students)
    infected_by = zeros(Int64,n_students)

    # delay_adherence = zeros(Int64,n_students)     # Individial may not report symptoms immediately.
    max_delay_adherence = length(delay_adherence_pmf) # Vars to be used when allocating reporting delay
    csum_delay_adherence = cumsum(delay_adherence_pmf)

    # Initialise undefined array
    # Used transmit_over! function
    undefined_array = Array{Int64,2}(undef,0,0)

    # Initialise variables used to allocate household infection delay
    # & perform check the pmf sums to 1 (or near 1)
    if (sum(delay_household_infection_pmf) > 1+eps()) || (sum(delay_household_infection_pmf) < 1-eps())
        error("delay_household_infection_pmf does not sum to 1. Check values.")
    end
    csum_delay_household_inf = cumsum(delay_household_infection_pmf)

    """
    If required, initialise contract tracing related variables
    """
    # If contact tracing in use, create variables
    if contact_tracing_active == true

        CT_vars = contact_tracing_vars(n_students=n_students,endtime=endtime,n_households=n_households)

        max_test_result_delay = length(CT_delay_until_test_result_pmf) # Vars to be used when allocating test delay times
        csum_test_result_delay = cumsum(CT_delay_until_test_result_pmf)

        # Array to track amount of ppl isolating as a result of contact tracing.
        # Row per timestep, column per replicate.
        num_isolating_CTcause = zeros(endtime,countfinal)
    else
        # otherwise we'll make a zero version of the variable
        CT_vars = contact_tracing_vars()
    end

    """
    If required, initialise work/study group inactivation variables
    """
    work_or_study_group_memory = Array{Array{Int64,2},1}(undef, n_cohorts) # initialise the memory for each work/study group

    if work_study_group_closure_active
        group_inactivation_time = Array{Array{Int64,1},1}(undef, n_cohorts) # initialise the timer for each group inactivation
        work_or_study_group_thresholds = Array{Array{Int64,1},1}(undef, n_cohorts) # initialise the threshold for group inactivation
        # work_or_study_group_memory = Array{Array{Int64,2},1}(undef, n_cohorts) # initialise the memory for each workplace
        num_teams = zeros(Int64,n_cohorts)
        for worktypeID = 1:n_cohorts
            num_teams[worktypeID] = length(class_sizes[worktypeID])
            group_inactivation_time[worktypeID] = zeros(num_teams[worktypeID])
            work_or_study_group_memory[worktypeID] = zeros(num_teams[worktypeID],work_or_study_group_CT_memory)
            work_or_study_group_thresholds[worktypeID] = zeros(num_teams[worktypeID])
            for class_ID = 1:num_teams[worktypeID]
                work_or_study_group_thresholds[worktypeID][class_ID] = ceil(Int64,class_sizes[worktypeID][class_ID]*work_or_study_group_CT_threshold)
            end
        end
    end

    """
    If required, initialise triggered intervention variables
    """
    # Check if any intervetion were specified
    if isassigned(intervention_fns)
        # Number of intervention sets provided is number of elements of intervention_fns
        n_intervention_fns = length(intervention_fns)
    end

    """
    Run replicates
    """
    # Perform countfinal number of replicates
    for count=1:countfinal

        """
        Set the RNG
        """
        rng = MersenneTwister(RNGseed+count)

        """
        Reinitialisation phase
        """

        # Re-initiailise node related arrays
        lmul!(0,dayon)
        lmul!(0,ccount)
        lmul!(0,pflag)

        # Construct & populate arrays signifiying when in class/in work setting
        populate_inclass!(states.inclass,sameday,ton,toff,rng)

        # Construct society schedule signifiying when society events run
        # Also reinitialise other fields in the society parameter structure
        reinitialise_society_params!(society_info,endtime,rng)

        # Reinitialise time series vectors
        reinitialise_node_states!(states)

        # Reinitialise daily record arrays
        reinitialise_daily_record_arrays!(contacts)

        # Reinitialise worker parameter fields
        reinitialise_student_params!(n_students,student_info)

        # Reinitialise workplace_params
        reinitialise_class_params!(class_info)

        # Reinitialise infected_by & time_to_symps
        lmul!(0,infected_by)
        lmul!(0,time_to_symps)

        """
        Set course of infection times
        """
        # set times to infection etc.: returns inftime, symptime, lattime, hh_isolation and delay_adherence
        set_infection_related_times!(time_to_symps,states,isolation,adherence,csum_delay_adherence,d_latent,n_students,rng)

        """
        Seed non-susceptible disease states
        """

        # Draw asymptomatic probability for current replicate
        probasymp = rand(rng,probasymp_dist)

        # Draw relative infectiousness of an asymptomatic for current replicate
        asymp_trans_scaling = rand(rng,asymp_trans_scaling_dist)

        # Sets latent, asymptomatic, symptomatic, recovered nodes
        n_initial_latent::Int64,
        n_initial_asymp::Int64,
        n_initial_symp::Int64,
        n_initial_recovereds::Int64 = seed_initial_states_fn(rng,
                                                        n_students,
                                                        count,
                                                        student_info,
                                                        states,
                                                        probasymp,
                                                        infected_by,
                                                        recov_propn,
                                                        output)

       """
       Update output time series with initial conditions
       """
       # Update time series for latent & infecteds after assigning initial
       # infecteds
       output.numlat[1,count] = n_initial_latent
       output.numinf[1,count] = n_initial_asymp + n_initial_symp

       # Update prevalences
       output.prevlat[1,count] = n_initial_latent
       output.prevasymp[1,count] = n_initial_asymp
       output.prevsymp[1,count] = n_initial_symp
       output.prevrec[1,count] = n_initial_recovereds


       """
       Reset contact tracing variables
       """
        # If required, set up and/or reinitialise contact tracing related variables
        if contact_tracing_active == true
            reinitialise_CT_vars!(CT_vars, n_students, rng, CT_parameters, states.delay_adherence,csum_test_result_delay,max_test_result_delay)
        end

       # Initialise counter for being able to identify the infector of an infectee
       recall_infector_count = 0

       # Counter for identified infector engaging in contact tracing
       # (having not participated in CT before)
       infector_trace_count = [0]

        """
        Run single replicate
        """
        for time=1:endtime

            # Initial timepoint is for initial conditions
            # Set row to accessed in output arrays for this timestep
            output_time_idx = time + 1

             """
             Reinitialise variables at start of timestep
             """
            # Reinitialise timestep specific values
            lmul!(0,states.rep_inf_this_timestep)

            # reinitalise the current work_or_study_group_memory slot
            if work_study_group_closure_active == true
                WP_memory_slot = mod1(time,work_or_study_group_CT_memory)
                for worktypeID = 1:n_cohorts
                    # Iterate over each work/study team for current work group type
                    n_teams = num_teams[worktypeID]
                    for team_itr = 1:n_teams
                        work_or_study_group_memory[worktypeID][team_itr,WP_memory_slot] = 0
                    end
                end
            end


            """
            Assign outputs
            """
            # Assign counts in each disease state to array
            output.numlat[output_time_idx,count] = output.numlat[output_time_idx-1,count]
            output.numinf[output_time_idx,count] = output.numinf[output_time_idx-1,count]
            output.numrep[output_time_idx,count] = output.numrep[output_time_idx-1,count]

            """
            Increment counters
            """
            # Increment counters if node is currently in that state.
            if contact_tracing_active==true
                increment_counters!(states,
                    household_isoltime,symp_isoltime,asymp_isoltime,
                    contact_tracing_active,
                    timeisol_CTcause=states.timeisol_CTcause,
                    CT_caused_isol_limit=CT_caused_isol_limit) # in additional_fns.jl
            else
                increment_counters!(states,
                    household_isoltime,symp_isoltime,asymp_isoltime,contact_tracing_active) # in additional_fns.jl
            end

            """
            Increment infection process
            """
           # If come to the end of latent time, move to infectious time etc
           for student_itr = 1:n_students
                # if the node has reached the end of latent infection
                if states.timelat[student_itr]>states.lattime[student_itr]
                    # move to being infectious
                    states.timelat[student_itr] = -1
                    states.timeinf[student_itr] = 1

                    # Increment time series counts
                    output.numinf[output_time_idx,count] += 1
                    output.newinf[output_time_idx,count] += 1

                    # check if new infected will be asymptomatic
                    # Also set household transmission risk
                    if states.asymp[student_itr] > 0
                        output.newasymp[output_time_idx,count] += 1
                    end

                    # check if it someone that would be attending in person class/
                    # working at the workplace that is newly infected
                    if student_info[student_itr].would_attend_f2f_classes==1
                        output.atworkinf[output_time_idx,count] += 1

                        # Count those that are asymptomatically infected
                        if states.asymp[student_itr] > 0
                            output.atworkasymp[output_time_idx,count] += 1
                        elseif work_study_group_closure_active==true
                            # if the newly infected person is symptomatic, add to
                            # the study group/work group memory
                            cohort_ID = student_info[student_itr].cohort_ID
                            class_ID = student_info[student_itr].class_ID
                            work_or_study_group_memory[cohort_ID][class_ID,WP_memory_slot] += 1
                        end
                    end
                end

                # Update node disease state time vectors
                if states.timeinf[student_itr]>states.inftime
                    # the node becomes symptomatic (if they develop symptoms)
                    states.timeinf[student_itr] = -1
                    states.timesymp[student_itr] = 1

                    # Increment time series counts
                    output.numrep[output_time_idx,count] += 1

                    # Check if index case are symptomatic & would have zero adherence delay
                    if (states.asymp[student_itr] == 0) && (states.delay_adherence[student_itr]==0)
                        # Check if infected will isolate
                        if (states.hh_isolation[student_itr]==1)
                            states.symp_timeisol[student_itr] = 1

                            # Set that the unit has reported infection this timestep
                            states.rep_inf_this_timestep[student_itr] = 1

                            # Assign time of reporting to field in student parameter type
                            student_info[student_itr].time_of_reporting_infection = time

                            # If an available option, no contacts state entered if student is living on-campus
                            # Will last until end of symptoms (irrespective of test result)
                            if (rehouse_strat_active == true) &&
                                    (student_info[student_itr].household_info.on_campus_accom == true)
                                student_info[student_itr].no_contacts_status = true

                                # Update rehoused output variable if living in communal
                                # bathroom type accomodation
                                if (student_info[student_itr].household_info.ensuite_flag == false)
                                    output.new_rehoused[output_time_idx,count] += 1
                                end
                            end
                        end

                        # Irrespective of whether index case self-isolates,
                        # adherent members of their household may also isolate.
                        for hh = 1:household_contacts_per_node[student_itr]
                            contact_ID = contacts.household_contacts[student_itr][hh]
                            if (states.hh_isolation[contact_ID]==1) &&
                                (states.symp_timeisol[contact_ID]==0) # Individual not already symptomatic themselves
                                states.timeisol[contact_ID] = 1
                            end
                        end
                    end

                    # If contact tracing active, increase number of symptomatic infections
                    # in household by one
                    if contact_tracing_active == true
                        if (states.asymp[student_itr] == 0) # Check case is symptomatic
                            current_node_household_ID = student_info[student_itr].household_info.household_ID
                            CT_vars.Symp_cases_per_household_pos_or_unknown[current_node_household_ID] += 1
                        end
                    end
                end

                # Check if node, if having a delayed adherence, begins adherence on current day
                if (states.timesymp[student_itr] > 1)&&((states.timesymp[student_itr]-1)==states.delay_adherence[student_itr]) # Condition for node beginning adherence on current day & has been symptomatic for at least one day
                    if states.asymp[student_itr] == 0 # Check node is symptomatic and will adhere
                        if states.hh_isolation[student_itr]==1 # Check node will adhere
                            states.symp_timeisol[student_itr] = 1 + states.delay_adherence[student_itr]

                            # Set that the unit has reported infection this timestep
                            states.rep_inf_this_timestep[student_itr] = 1

                            # Assign time of reporting to field in student parameter type
                            student_info[student_itr].time_of_reporting_infection = time

                            # If an available option, no contacts state entered
                            # Will last until end of symptoms (irrespective of test result)
                            if rehouse_strat_active == true
                                student_info[student_itr].no_contacts_status = true
                            end
                        end

                        # Household members now aware of index case having symptoms.
                        # Adherent members of their household may also now isolate,
                        # assuming infected displays symptoms
                        # Note they are delayed in isolating, in line with delay
                        # of the index case
                        for hh = 1:household_contacts_per_node[student_itr]
                            contact_ID = contacts.household_contacts[student_itr][hh]
                            if (states.hh_isolation[contact_ID]==1) && (states.symp_timeisol[contact_ID]==0) # Individual not already symptomatic themselves
                                states.timeisol[contact_ID] = 1 + states.delay_adherence[student_itr]
                                    # Individual shortens isolation by length of time since
                                    # unwell individual began displaying symptoms
                            end
                        end
                    end
                end

                # Check if node has reached end of symptom period
                if states.timesymp[student_itr]>states.symptime
                    states.timesymp[student_itr] = -1

                    # If an active intervention option, reset no contacts state
                    if rehouse_strat_active == true
                        student_info[student_itr].no_contacts_status = false
                    end

                    # If contact tracing active and case was symptomatic,
                    # decrease number of symptomatic infections in household by one
                    if contact_tracing_active == true
                        # Check case is symptomatic & not returned a false negative (if false negative, has already been subtracted)
                        if (states.asymp[student_itr] == 0) && (CT_vars.Test_result_false_negative[student_itr] == false)
                            current_node_household_ID = student_info[student_itr].household_info.household_ID
                            CT_vars.Symp_cases_per_household_pos_or_unknown[current_node_household_ID] -= 1

                            # Error check
                            if CT_vars.Symp_cases_per_household_pos_or_unknown[current_node_household_ID] < 0
                                error("CT_vars.Symp_cases_per_household_pos_or_unknown contains a negative entry. Terminate programme.")
                            end

                        end
                    end
                end
            end

            """
            Record whether individuals are in isolation. Set study and society attendance status for timestep.
            """
            # record whether individuals are in isolation,
            # & whether attending in person classes/staff are at workplace
            # & whether they attend a society that day
            for student_itr = 1:n_students
                # Record whether the individual is in isolation on current day
                if (isolation>0) &&
                    ((states.timeisol[student_itr]>0) || (states.symp_timeisol[student_itr]>0) ||
                       (states.timeisol_CTcause[student_itr]>0) || (student_info[student_itr].household_info.lockdown_status == true))
                    # Default value is 0. So only update if node is in isolation for any reason
                    contacts.daily_record_inisol[time,student_itr] = 1
                end

                # Record whether the person is attending a class in person/at work on current day
                # Needs to not be WFH, and a day where does attend class/go to work
                # AND not in isolation
                # AND class/work setting has in person activity
                cohort_ID = student_info[student_itr].cohort_ID
                class_ID = student_info[student_itr].class_ID
                node_team_info = class_info[cohort_ID][class_ID]
                if (student_info[student_itr].would_attend_f2f_classes==1) &&
                    (states.inclass[student_itr,time] == true) &&
                    (contacts.daily_record_inisol[time,student_itr] == false) &&
                    (node_team_info.f2f_activity == true)

                    # Default value is 0. So only update if person is attending
                    # an in-person class/is attending a work setting
                    contacts.daily_record_inclass[time,student_itr] = 1
                end

                # Record whether the person attends the societies there are members
                # of on the current timestep
                n_societies_joined = length(student_info[student_itr].society_IDs)
                if (n_societies_joined > 0) && (contacts.daily_record_inisol[time,student_itr] == false)
                    # Pass check that they are members of a society and not currently isolating

                    # Get the relevant society IDs for the person of interest
                    society_ID_list = student_info[student_itr].society_IDs

                    # Iterate over each society.
                    for society_itr = 1:n_societies_joined
                        current_society_ID = society_ID_list[society_itr]

                        # Update attendance status if society is active (i.e. running face-to-face events)
                        # & society event occurred on current timestep
                        if (society_info[current_society_ID].f2f_activity == true) &&
                             (society_info[current_society_ID].schedule[time] == true)

                            # Update if person is not isolating and society is active
                            contacts.daily_record_atsociety[time,current_society_ID,student_itr] = 1
                        end
                    end
                end
            end


            """
            Transmit infections

            Structure:
            - Household
            - In class/at workplace setting
            - Society contacts
            - Dynamic contacts
            """
            # Iterate over nodes that may be able to transmit infection
            for student_itr = 1:n_students
                if ((states.timeinf[student_itr]>0) | (states.timesymp[student_itr]>0)) &&
                        (student_info[student_itr].no_contacts_status == false)
                        # Only enter loop if node is capable of transmitting infection
                        # & not in rehousing/room isolation.

                    # find the total time infectious
                    if states.timeinf[student_itr]>0
                        tot_time_inf = states.timeinf[student_itr]
                    else
                        tot_time_inf = states.timesymp[student_itr]+states.inftime
                    end

                    # find the infectiousness
                    infectiousness = dist_infectivity[tot_time_inf]
                    current_student = student_info[student_itr]
                    current_study_work_type_ID = student_info[student_itr].cohort_ID
                    if states.asymp[student_itr]>0 # Asymptomatic
                        transtemp_household = current_student.transrisk_household*infectiousness*asymp_trans_scaling*suscep_scaling
                        transtemp_cohort = current_student.transrisk_cohort*infectiousness*asymp_trans_scaling*suscep_scaling
                        transtemp_dynamic_social = current_student.transrisk_dynamic_social*infectiousness*asymp_trans_scaling*suscep_scaling
                    else
                        if states.timesymp[student_itr]>0 # symptomatic & less infectious due to cautionary behaviour
                            transtemp_household = current_student.transrisk_household*infectiousness*iso_trans_scaling*suscep_scaling
                            transtemp_cohort = current_student.transrisk_cohort*infectiousness*iso_trans_scaling*suscep_scaling
                            transtemp_dynamic_social = current_student.transrisk_dynamic_social*infectiousness*iso_trans_scaling*suscep_scaling
                        else  # infected, unscaled infectiousness
                            transtemp_household = current_student.transrisk_household*infectiousness*suscep_scaling
                            transtemp_cohort = current_student.transrisk_cohort*infectiousness*suscep_scaling
                            transtemp_dynamic_social = current_student.transrisk_dynamic_social*infectiousness*suscep_scaling
                        end
                    end

                    # Infection check for other household members
                    # Transmit over household_contacts[student_itr]
                    # checking that contacts are susceptible
                    n_hh_contacts = length(contacts.household_contacts[student_itr])
                    if n_hh_contacts > 0
                           transmit_over!(student_info,transtemp_household,infected_by,output,states,probasymp,rng,time,count,
                                                infecting_by = student_itr,
                                                contacts_to_check = contacts.household_contacts[student_itr],
                                                household_check_flag = true)
                    end

                    # Check individual is not in household isolation.
                    # If not, can see if transmission across cohort, society,
                    # dynamic household contacts occured.
                    if (contacts.daily_record_inisol[time,student_itr] == false)

                            # Check if node is at in-person classes/at the workplace
                            if contacts.daily_record_inclass[time,student_itr] == true

                                # These are the contacts with those in the same class
                                n_class_atwork_contacts = work_or_study_group_contacts_per_node[student_itr]
                                if n_class_atwork_contacts > 0
                                    # Satisfied condition that student_itr has any class/at work links
                                    # transmit over contacts.class_contacts[student_itr]
                                    # checking that contacts are susceptible, not isolating and at work
                                    transmit_over!(student_info,transtemp_cohort,infected_by,output,states,probasymp,rng,time,count,
                                            infecting_by = student_itr,
                                            contacts_to_check = contacts.class_contacts[student_itr],
                                            inisol = contacts.daily_record_inisol,
                                            attendance_record = contacts.daily_record_inclass,
                                            cohort_ID = current_study_work_type_ID)
                                end

                                # Contacts made with those in the same cohort, but a different class
                                n_cohort_contacts = cohort_contacts_per_node[student_itr]
                                if n_cohort_contacts > 0
                                    # Satisfied condition that student_itr has any cohort links
                                    # transmit over contacts.cohort_contacts[student_itr]
                                    # checking that contacts are susceptible, not isolating and at work
                                    transmit_over!(student_info,transtemp_cohort,infected_by,output,states,probasymp,rng,time,count,
                                            infecting_by = student_itr,
                                            contacts_to_check = contacts.cohort_contacts[student_itr],
                                            inisol = contacts.daily_record_inisol,
                                            attendance_record = contacts.daily_record_inclass,
                                            cohort_ID = current_study_work_type_ID)
                                end
                            end

                            # Check society contacts
                            n_societies_joined = length(student_info[student_itr].society_IDs)
                            if (n_societies_joined > 0)

                                # Satisfied condition that student_itr has society links
                                # transmit over contacts.society_contacts[student_itr]
                                # checking that contacts are susceptible, not isolating and at work
                                transmit_over_societies!(student_info,
                                                            society_info,
                                                            infection_parameters,
                                                            infectiousness,
                                                            asymp_trans_scaling,
                                                            infected_by,output,states,probasymp,rng,time,count,
                                                            infecting_by = student_itr,
                                                            society_contacts = contacts.society_contacts,
                                                            inisol = contacts.daily_record_inisol,
                                                            society_attendance_record = contacts.daily_record_atsociety)
                            end

                            # Check if any social dynamic contacts made by this student on this timestep
                            transmit_over!(student_info,transtemp_dynamic_social,infected_by,output,states,
                                        probasymp,rng,time,count,
                                        infecting_by = student_itr,
                                        contacts_to_check = contacts.dynamic_social_contacts[time,student_itr],
                                        inisol = contacts.daily_record_inisol,
                                        attendance_record = undefined_array,
                                        dynamic_contact = 1)

                            # Check if any accomodation dynamic contacts made by this student on this timestep
                            # Only check if student lives on campus
                            if student_info[student_itr].household_info.on_campus_accom == true
                                for accom_level = 1:3
                                    # Iterate over hall, block, floor contacts
                                    # If there may be dynamic accomodation contacts assigned, check if transmission occurs
                                    if isassigned(contacts.dynamic_accomodation_contacts,time,accom_level,student_itr)
                                        transmit_over!(student_info,transtemp_household,infected_by,output,states,
                                                                probasymp,rng,time,count,
                                                                infecting_by = student_itr,
                                                                contacts_to_check = contacts.dynamic_accomodation_contacts[time,accom_level,student_itr],
                                                                inisol = contacts.daily_record_inisol,
                                                                attendance_record = undefined_array,
                                                                dynamic_contact = 2)
                                    end
                                end
                            end
                    end
                end
            end

            """
            Perform contact tracing
            """
            # If in use, enact contact tracing from index cases reported today
            if contact_tracing_active == true

                    # Store contacts made during day.
                    # For those reporting symptoms, start delay to test result (if needed)
                    for student_itr = 1:n_students

                        # Increment time to test result if currently waiting for that to occur
                        if CT_vars.Time_to_test_result[student_itr]>=0
                            CT_vars.Time_to_test_result[student_itr] += 1
                        end

                        # For current person, check if would be leaving pre-symptomatic phase
                        # and displaying symptoms
                        # If so, and will not return a false negative test result, gather traceable contacts
                        if (states.rep_inf_this_timestep[student_itr] == 1)

                            # Initialise CT_vars.Time_to_test_result value
                            CT_vars.Time_to_test_result[student_itr] = 0

                            # Increment test counter
                            output.tests_performed[output_time_idx,count] += 1

                            # Determine whether test result will return a negative outcome
                            # - Get time since student_itr became infected
                            # - Given time since infected, look up probability case will return negative test result
                            if states.timeinf[student_itr]>0
                                tot_time_inf = states.timeinf[student_itr]
                            else
                                tot_time_inf = states.timesymp[student_itr]+states.inftime
                            end
                            test_false_negative_prob = test_false_negative_vec[tot_time_inf]

                            # Bernoulli trial to determine if false negative returned
                            if rand(rng) < test_false_negative_prob
                                CT_vars.Test_result_false_negative[student_itr] = true
                            end

                            # Check if the individual is destined to return a false negative result
                            # If false negative to be returned, do not need to work out who traceable contacts are
                            # Otherwise, gather traceable contacts
                            # Also, if no isolation in use, no need to gather contacts.
                            CT_vars.Inds_to_be_contacted[student_itr] = Int64[] # Initialise vector to store contacts
                            if (CT_vars.Test_result_false_negative[student_itr] == false) &&
                                    (isolation>0) && (CT_vars.Engage_with_CT[student_itr] == true)

                                trace_node!(student_itr,time,CT_vars,contacts,CT_parameters,network_parameters,rng)

                                # if we are doing "backward contact tracing"
                                # some small chance the infector is included in this
                                # don't try to backwards trace the initial infections
                                if ( (rand(rng)<prob_backwards_CT) && (infected_by[student_itr]!=-1) &&
                                        (infected_by[student_itr]!=-2) && (infected_by[student_itr]!=-3) )
                                    append!(CT_vars.Inds_to_be_contacted[student_itr],infected_by[student_itr])
                                        CT_vars.Recall_infector[student_itr] = 1
                                end
                            end

                            # Remove duplicates in CT_vars.Inds_to_be_contacted[student_itr]
                            unique!(CT_vars.Inds_to_be_contacted[student_itr])
                        end

                        # Check if delay to test result reached
                        if CT_vars.Time_to_test_result[student_itr] >= CT_vars.CT_delay_until_test_result[student_itr]

                            # Reset the CT_vars.Time_to_test_result counter
                            CT_vars.Time_to_test_result[student_itr] = -1

                            # If delay time passed, check if test returned a false negative.
                            if CT_vars.Test_result_false_negative[student_itr] == true

                                # Increment false negative counter
                                output.test_outcomes[output_time_idx,count,2] += 1

                                # Release index case from symptomatic caused isolation
                                # If yes, they cannot be released
                                # If no, they can be released.
                                states.symp_timeisol[student_itr] = 0

                                # Amend tracker of symptomatic cases, unknown test result
                                current_node_household_ID = student_info[student_itr].household_info.household_ID
                                CT_vars.Symp_cases_per_household_pos_or_unknown[current_node_household_ID] -= 1

                                # Error check
                                if CT_vars.Symp_cases_per_household_pos_or_unknown[current_node_household_ID] < 0
                                    error("CT_vars.Symp_cases_per_household_pos_or_unknown contains a negative entry. Terminate programme.")
                                end

                                # If returning a false negative, can release other household members from isolation
                                # Only release if no other household members are symptomatic
                                if CT_vars.Symp_cases_per_household_pos_or_unknown[current_node_household_ID] == 0
                                    n_household_contacts = length(contacts.household_contacts[student_itr])
                                    for household_contact_itr = 1:n_household_contacts
                                        household_contact_ID = contacts.household_contacts[student_itr][household_contact_itr]
                                        if states.timeisol_CTcause[household_contact_ID] == 0 # Check not also self isolating due to contact tracing
                                            states.timeisol[household_contact_ID] = 0
                                        end
                                    end
                                end

                            else

                                # Increment true positive counter
                                output.test_outcomes[output_time_idx,count,1] += 1

                                # If test is positive, contacts told to isolate
                                # Get number of recallable contacts
                                n_recallable_contacts = length(CT_vars.Inds_to_be_contacted[student_itr])
                                output.num_CT[output_time_idx,count] += n_recallable_contacts

                                # Do contacts isolate or not?
                                # Isolation based on adherence to household isolation of that individual
                                if isolation>0
                                    for recallable_contact_idx = 1:n_recallable_contacts
                                        recallable_contact_ID = CT_vars.Inds_to_be_contacted[student_itr][recallable_contact_idx]

                                        # Check if individual will adhere to guidance
                                        # If so, they self-isolate
                                        if states.hh_isolation[recallable_contact_ID] == 1
                                            states.timeisol_CTcause[recallable_contact_ID] = 1
                                        end
                                    end

                                    # Perform forwards CT from infector, if infector has been identified
                                    # and such a policy is active
                                    if perform_CT_from_infector == true
                                        if CT_vars.Recall_infector[student_itr] == 1
                                            recall_infector_count += 1

                                            forwardCT_from_infector!(infected_by[student_itr],
                                                                        CT_vars,
                                                                        contacts,
                                                                        CT_parameters,
                                                                        CT_vars.Engage_with_CT,
                                                                        states.inclass,
                                                                        time,
                                                                        count,
                                                                        output.num_CT,
                                                                        states.hh_isolation,
                                                                        states.timeisol_CTcause,
                                                                        network_parameters,
                                                                        rng,
                                                                        infector_trace_count)
                                        end
                                    end
                                end
                            end
                        end
                    end
            end


            """
            Assign prevalence & isolation outputs
            """
            # For this timestep, get number isolating
            # and if they are latent infected or infectious on current timestep
            for student_itr=1:n_students

                # Isolating due to housemate symptoms
                isolating_for_any_reason = false
                if states.timeisol[student_itr]>0
                    output.num_household_isolating[output_time_idx,count] += 1
                    isolating_for_any_reason = true
                end

                # Isolating due to symptoms
                if states.symp_timeisol[student_itr]>0
                    output.num_symp_isolating[output_time_idx,count] += 1
                    isolating_for_any_reason = true
                end

                # Isolating due to finding asymptomatic infection (due to testing)
                if states.asymp_timeisol[student_itr]>0
                    output.num_asymp_isolating[output_time_idx,count] += 1
                    isolating_for_any_reason = true
                end

                # Isolating due to accomodation lockdown
                if (student_info[student_itr].household_info.lockdown_status == true)
                    output.num_accom_lockdown_isolating[output_time_idx,count] += 1
                    isolating_for_any_reason = true
                end

                # Isolating as close contact of positive test symptomtic
                if contact_tracing_active == true
                    if states.timeisol_CTcause[student_itr]>0
                        output.num_isolating_CTcause[output_time_idx,count] += 1
                        isolating_for_any_reason = true
                    end
                end

                # Isolating for any reason
                if isolating_for_any_reason == true
                    output.num_isolating[output_time_idx,count] += 1

                    # Check if on-campus or off-campus resident. Update relevant output counter
                    if (student_info[student_itr].household_info.on_campus_accom == true)
                         # On-campus resident
                        output.num_isolating_oncampus[output_time_idx,count] += 1
                    else
                        # Off-campus resident
                        output.num_isolating_offcampus[output_time_idx,count] += 1
                    end
                end

                # Check if latently infected
                if states.timelat[student_itr]>0
                    output.prevlat[output_time_idx,count] += 1

                    # Check if on-campus or off-campus resident. Update relevant output counter
                    if (student_info[student_itr].household_info.on_campus_accom == true)
                         # On-campus resident
                        output.prevlat_oncampus[output_time_idx,count] += 1
                    else
                        # Off-campus resident
                        output.prevlat_offcampus[output_time_idx,count] += 1
                    end
                end

                # In presymptomatic infectious period.
                if (states.timeinf[student_itr]>0)
                    if states.asymp[student_itr] > 0 # asymptomatic
                        output.prevasymp[output_time_idx,count] += 1

                        # Check if on-campus or off-campus resident. Update relevant output counter
                        if (student_info[student_itr].household_info.on_campus_accom == true)
                             # On-campus resident
                            output.prevasymp_oncampus[output_time_idx,count] += 1
                        else
                            # Off-campus resident
                            output.prevasymp_offcampus[output_time_idx,count] += 1
                        end
                    else # will be symptomatic
                        output.prevpresymp[output_time_idx,count] += 1

                        # Check if on-campus or off-campus resident. Update relevant output counter
                        if (student_info[student_itr].household_info.on_campus_accom == true)
                             # On-campus resident
                            output.prevpresymp_oncampus[output_time_idx,count] += 1
                        else
                            # Off-campus resident
                            output.prevpresymp_offcampus[output_time_idx,count] += 1
                        end
                    end
                end

                # After presymp period, check if symptomatic or asymptomatic
                if (states.timesymp[student_itr]>0)
                    if states.asymp[student_itr] > 0 # asymptomatic
                        output.prevasymp[output_time_idx,count] += 1

                        # Check if on-campus or off-campus resident. Update relevant output counter
                        if (student_info[student_itr].household_info.on_campus_accom == true)
                             # On-campus resident
                            output.prevasymp_oncampus[output_time_idx,count] += 1
                        else
                            # Off-campus resident
                            output.prevasymp_offcampus[output_time_idx,count] += 1
                        end
                    else # symptomatic
                        output.prevsymp[output_time_idx,count] += 1

                        # Check if on-campus or off-campus resident. Update relevant output counter
                        if (student_info[student_itr].household_info.on_campus_accom == true)
                             # On-campus resident
                            output.prevsymp_oncampus[output_time_idx,count] += 1
                        else
                            # Off-campus resident
                            output.prevsymp_offcampus[output_time_idx,count] += 1
                        end
                    end
                end

                # Check if recovered
                if states.timesymp[student_itr] == -1
                    output.prevrec[output_time_idx,count] += 1

                    # Check if on-campus or off-campus resident. Update relevant output counter
                    if (student_info[student_itr].household_info.on_campus_accom == true)
                         # On-campus resident
                        output.prevrec_oncampus[output_time_idx,count] += 1
                    else
                        # Off-campus resident
                        output.prevrec_offcampus[output_time_idx,count] += 1
                    end
                end

                # Check if student is currently rehoused
                if (student_info[student_itr].no_contacts_status == true) &&
                        (student_info[student_itr].household_info.ensuite_flag == false)
                    output.current_rehoused[output_time_idx,count] += 1
                end
            end

            """
            Reactive inactivation of in person classes/society meet ups
            Check if needed
            """
            # close workplaces with too many infections
            if work_study_group_closure_active==true
                for worktypeID = 1:n_cohorts
                    for class_ID = 1:num_teams[worktypeID]
                        if group_inactivation_time[worktypeID][class_ID]>0
                            # if the activity is inactive, move on the time counter
                            group_inactivation_time[worktypeID][class_ID] += 1
                        else
                            # otherwise decide if the activity should be allowed to continue
                            total_infections = 0
                            for ii=1:work_or_study_group_CT_memory
                                total_infections += work_or_study_group_memory[worktypeID][class_ID,ii]
                            end

                            # Number of infections in time window exceeds threshold
                            # Set workplace to be closed
                            if total_infections>work_or_study_group_thresholds[worktypeID][class_ID]
                                group_inactivation_time[worktypeID][class_ID] = 1
                                network_parameters.class_info[worktypeID][class_ID].f2f_activity = false
                            end
                        end

                        if group_inactivation_time[worktypeID][class_ID]>time_WC
                            # if the group has been inactive for long enough, restart it
                            group_inactivation_time[worktypeID][class_ID] = 0

                            # Check f2f study is permitted for the cohort.
                            # If it is permitted, the individual class can be set to be active.
                            if (network_parameters.cohort_f2f_study_active[worktypeID] == true)
                                network_parameters.class_info[worktypeID][class_ID].f2f_activity = true
                            end
                        end
                    end
                end

            end

            """
            Run interventions
            """
            # Check if any interventions are triggered
            # Update statuses as needed
            if isassigned(intervention_fns) # Check if any intervetion were specified

                # Package health outcome measures that may be used in decision
                # to enact an intervention
                intervention_trigger_input_data = intervention_data_feeds(rep_inf_this_timestep = states.rep_inf_this_timestep,
                                                                            output = output,
                                                                            network_params = network_parameters,
                                                                            time = time,
                                                                            output_time_idx = output_time_idx,
                                                                            replicate_id = count
                                                                            )

                for interv_trig_itr = 1:n_intervention_fns
                    chosen_affect_fn = intervention_fns[interv_trig_itr]
                    chosen_affect_fn(network_parameters,
                                        intervention_trigger_input_data,
                                        contacts,
                                        time)
                end
            end


            """
            Perform mass testing (if applicable)
            """
            if mass_testing_active == true
                perform_mass_test!(mass_testing_parameters,
                                                      time,
                                                      count,
                                                      states,
                                                      CT_vars,
                                                      CT_parameters,
                                                      network_parameters,
                                                      contacts,
                                                      output,
                                                      household_contacts_per_node,
                                                      rng,
                                                      rehouse_strat_active)
            end
        end

        # Find how many nodes were infected by the initial nodes
        initial_nodes = findall(infected_by.==-1)
        sum_infections = 0
        output.num_init_infected[count] = zeros(Int64,length(initial_nodes)) # Initialise output vector
        for initial_node_it = 1:length(initial_nodes)
            output.num_init_infected[count][initial_node_it] = output.num_infected[initial_nodes[initial_node_it],count]
            sum_infections+=output.num_init_infected[count][initial_node_it]
        end

        # Find how many nodes were set to adhere to isolation guidance for this replicate
        output.n_isol_adhering[count] = sum(states.hh_isolation)
        for student_itr=1:n_students # Iterate over each student. Update oncampus and offcampus specific variables
            # Check if student would have been isolating
            if states.hh_isolation[student_itr] == 1
                # If yes, check if on-campus or off-campus resident. Increment counters
                if (student_info[student_itr].household_info.on_campus_accom == true)
                        output.n_isol_adhering_oncampus[count] += 1
                else
                        output.n_isol_adhering_offcampus[count] += 1
                end
            end
        end


        # find mean generation time
        if sum_infections>0
            output.mean_init_generation_time[count] = output.mean_init_generation_time[count]/sum_infections
        end

        # divide number of infections by number of infectors to find Rt
        for time=1:(endtime+1)
            # divide by the number of nodes that were infected at time
            output.Rt[time,count] = output.Rt[time,count] / output.newinf[time,count]
        end

        # Print to screen info on run just completed
        println("Run $count complete.")
    end

    # Compute variance in number of infected per node
    var_num_infected = zeros(countfinal)
    for count=1:countfinal
        var_num_infected[count] = var(output.num_infected[:,count])
    end

    # Specify what is output from the function
    return output
end
