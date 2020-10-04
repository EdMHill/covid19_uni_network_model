"""
Purpose:
Functions to produce the network layers

Functions for assigning students to classes, without data
- generate_classes_default     (Does not use data)

Functions for assigning students to classes, using data
- generate_classes_campus_only (Oncampus assignments based on provided data)
- generate_classes_all_students (Includes offcampus assignment from appropriate distribution.
                                    Oncampus assignments based on provided data)
   Uses the following sub-functions:
    - get_oncampus_cohort_data
    - create_class_variables
    - create_class_param_type!

Functions for assigning societies
- assign_societies

Functions for producing the contact networks
- generate_contacts_uni                     (Generate the networks)
- configuration_model!                     (Network construction using configuration model)
- ER_model!                                (Network construction using erdos-renyi)
- generate_dynamic_student_contacts         (Premake dynamic social contacts,
                                            to be loaded in ahead of simulation)
- generate_dynamic_campus_accom_contacts   (Premake dynamic contacts within halls of residence,
                                            to be loaded in ahead of simulation)
- generate_contacts_ER_class_and_society_groups (Network construction using ER)

Functions for assigning students to societies
- assign_societies_one_per_student
- assign_societies_from_aggregated_data
- assign_societies_from_individual_data

Batch of functions for generating student households
- get_on_campus_households_from_data, get_on_campus_households_by_cohort,
   populate_oncampus_contact_lists! & assign_offcampus_households! (used in assign_households_fn!)
- assign_households_fn!


Also a batch of functions to regenerate layers of the network if needed (can be made for ER and configuration model variants)
- regenerate_class_contacts                (static worker contact layer)
- regenerate_society_contacts               (society contact layer)
- regenerate_dynamic_student_contacts        (dynamic worker contact layer)
"""

"""
Functions for allocating students to classes, not using data
"""

# General structure:
# Inputs:
# n_students::Int64 - Number of individuals in the system
# class_generation_parameters::class_generation_params - Parameter structure with information relevant to class structure
# RNGseed::Int64 - Used to seed the random number generator

# Outputs:
# student_info::Array{student_params,1} - Parameter type, with relevant fields of info for each student
# class_sizes::Array{Array{Int64,1},1} - List of class sizes, with separate vector per cohort
# class_info::Array{Array{class_params,1},1} - Parameter stucture with class fields. Separate vector per cohort.
# nodes_by_class::Array{Array{Array{Int64,1},1},1} - List of student IDs per class (vector of vectors per cohort))

# Initialise vectors determining how many people there are per class
# and which study group each person belongs to
function generate_classes_default(n_students::Int64,
                                        class_generation_parameters::class_generation_params,
                                        RNGseed::Int64)

    # Set the RNG
    rng = MersenneTwister(RNGseed)

    # Get team params
    @unpack n_cohorts, attendence_propns, classtype_proportion,
            class_size_mean, class_size_sd = class_generation_parameters

    # total number of people in each type of study/work group
    student_numbers = [round(Int64,n_students*classtype_proportion[i]) for i=1:length(classtype_proportion)]

    # correct for rounding
    if sum(student_numbers) != n_students
        diff = sum(student_numbers) - n_students
        student_numbers[end] -= diff
    end

    # Initialise places to store study/work team sizes per work group type
    # & team information
    class_sizes = Array{Array{Int64,1},1}(undef, n_cohorts)
    class_info = Array{Array{class_params,1},1}(undef, n_cohorts)

    # Initialise ID vars
    class_ids = Array{Array{Int64,1},1}(undef, n_cohorts)  # representing team ID per study/work group
    cohort_ids = Int64[]  # representing worktype ID per person

    # Initialise array to store ids of students within each class
    nodes_by_class = Array{Array{Array{Int64,1},1},1}(undef,n_cohorts)

    ## Now generate classes within each study/work group type
    for worktype_idx = 1:n_cohorts
        # initialise arrays for this study/work group type
        class_info[worktype_idx] =  class_params[]
        class_sizes[worktype_idx] = Int64[]
        class_ids[worktype_idx] = Int64[]

        # Initialise array to store ids of workers within each workplace
        nodes_by_class[worktype_idx] = Array[]

        # Set up tracking counters
        team_count = 0   # track number of teams in this study/work group type
        person_count = 0        # track number of people assigned to each team

        # Add groups until person_count exceeds target number of people in that study/work group
        while person_count < student_numbers[worktype_idx]

            # Update counter
            team_count += 1

            # generate working groups of random size and add to array
            push!(class_sizes[worktype_idx], ceil(Int64,abs(randn(rng)*class_size_sd[worktype_idx]+class_size_mean[worktype_idx])))
            append!(class_ids[worktype_idx], repeat([team_count], class_sizes[worktype_idx][end]))

            # Add empty class to nodes_by_class array
            push!(nodes_by_class[worktype_idx], Int64[])

            # Update tracking counter values
            person_count += class_sizes[worktype_idx][end]
        end

        # remove extra people
        class_sizes[worktype_idx][end] = student_numbers[worktype_idx] - sum(class_sizes[worktype_idx][1:end-1])
        class_ids[worktype_idx] = class_ids[worktype_idx][1:sum(class_sizes[worktype_idx])]

        # randomly shuffle people between study/work groups
        shuffle!(rng,class_ids[worktype_idx])

        append!(cohort_ids, repeat([worktype_idx],sum(class_sizes[worktype_idx])))

        # Create team parameter type for this study/work group type
        class_info[worktype_idx] =  Array{class_params,1}(undef,team_count)
        for team_itr = 1:team_count
            class_info[worktype_idx][team_itr] = class_params()
        end
    end

    # randomly shuffle n_cohorts
    shuffle!(rng,cohort_ids)

    student_info = Array{student_params,1}(undef,n_students)
    for ii = 1:n_students
        # decide if returning to work
        would_attend_f2f_classes = Int64(rand(rng)<attendence_propns[cohort_ids[ii]])

        # Construct node parameter type for current individual
        student_info[ii] = student_params(would_attend_f2f_classes = would_attend_f2f_classes,
                                            cohort_ID = cohort_ids[ii],
                                            class_ID = class_ids[cohort_ids[ii]][1])

        # Add node to class if would attend in person
        if would_attend_f2f_classes == 1
            push!(nodes_by_class[cohort_ids[ii]][class_ids[cohort_ids[ii]][1]], ii)
        end
        deleteat!(class_ids[cohort_ids[ii]],1) # remove worker from top of list
    end

    return student_info::Array{student_params,1},
            class_sizes::Array{Array{Int64,1},1},
            class_info::Array{Array{class_params,1},1},
            nodes_by_class::Array{Array{Array{Int64,1},1},1}
end


"""
Functions for allocating students to classes, using data
"""

# Get cohort data from file
function get_oncampus_cohort_data(n_students::Int64, n_depts::Int64)

    # Load data from file
    # Two columns. Column 1 for department ID. Column 2 for undergrad (=1) or postgrad (=2)
    oncampus_dept_data_raw = readdlm("../Data/department_student_data.txt",Int64)

    # Number of rows of input data gives size of campus population
    n_students_living_oncampus::Int64 = size(oncampus_dept_data_raw,1)

    # Assign the students to relevant cohort
    # Track student numbers in each cohort (will be used for class generation)
    cohort_ids = zeros(Int64,n_students)
    student_numbers_per_cohort = zeros(Int64,n_cohorts)
    for student_itr = 1:n_students_living_oncampus
        assigned_dept = oncampus_dept_data_raw[student_itr,1]
        study_level = oncampus_dept_data_raw[student_itr,2]

        # Get current student cohort ID based on UG/PG status
        if study_level == 1 # First year undergrad study
            student_cohort_ID = assigned_dept
        elseif study_level == 2 # Non-first year undergrad study
            student_cohort_ID = assigned_dept + n_depts
        elseif study_level == 3 # Postgrad study. In 2nd half of list of cohort IDs.
            student_cohort_ID = assigned_dept + (2*n_depts)
        else
            # Sanity check. Throw error if invalid value
            error("study_level has value $study_level. Must take value 1, 2 or 3.")
        end

        # Update cohort ID vector & counts per cohort
        cohort_ids[student_itr] = student_cohort_ID
        student_numbers_per_cohort[student_cohort_ID] += 1
    end

    return cohort_ids::Array{Int64,1},
            student_numbers_per_cohort::Array{Int64,1}
end

function create_class_variables(n_cohorts::Int64,
                                student_numbers_per_cohort::Array{Int64,1},
                                class_generation_parameters::class_generation_params,
                                rng::MersenneTwister)
# Inputs:
# As names describe

# Outputs:
# class_sizes::Array{Array{Int64,1},1}
# class_info::Array{Array{class_params,1},1}
# class_ids::Array{Array{Int64,1},1}
# nodes_by_class::Array{Array{Array{Int64,1},1},1}

    # Get class params
    @unpack class_size_mean, class_size_sd = class_generation_parameters

    # Initialise places to store class sizes per cohort type
    # & team information
    class_sizes = Array{Array{Int64,1},1}(undef, n_cohorts)
    class_info = Array{Array{class_params,1},1}(undef, n_cohorts)

    # Initialise ID vars
    class_ids = Array{Array{Int64,1},1}(undef, n_cohorts)  # representing team ID per study/work group

    # Initialise array to store ids of students within each class
    nodes_by_class = Array{Array{Array{Int64,1},1},1}(undef,n_cohorts)

    ## Now generate classes within each study/work group type
    for cohort_idx = 1:n_cohorts
        # initialise arrays for this study/work group type
        class_info[cohort_idx] =  class_params[]
        class_sizes[cohort_idx] = Int64[]
        class_ids[cohort_idx] = Int64[]

        # Initialise array to store ids of workers within each workplace
        nodes_by_class[cohort_idx] = Array[]

        # Set up tracking counters
        team_count = 0   # track number of teams in this study/work group type
        person_count = 0        # track number of people assigned to each team

        # Add groups until person_count exceeds target number of people in that study/work group
        while person_count < student_numbers_per_cohort[cohort_idx]

            # Update counter
            team_count += 1

            # generate working groups of random size and add to array
            push!(class_sizes[cohort_idx], ceil(Int64,abs(randn(rng)*class_size_sd[cohort_idx]+class_size_mean[cohort_idx])))
            append!(class_ids[cohort_idx], repeat([team_count], class_sizes[cohort_idx][end]))

            # Add empty class to nodes_by_class array
            push!(nodes_by_class[cohort_idx], Int64[])

            # Update tracking counter values
            person_count += class_sizes[cohort_idx][end]
        end

        # For non-empty cohorts, cleanup class sizes and create class parameter type variables
        if student_numbers_per_cohort[cohort_idx] > 0
            # remove extra people
            class_sizes[cohort_idx][end] = student_numbers_per_cohort[cohort_idx] - sum(class_sizes[cohort_idx][1:end-1])
            class_ids[cohort_idx] = class_ids[cohort_idx][1:sum(class_sizes[cohort_idx])]

            # randomly shuffle people between class groups
            shuffle!(rng,class_ids[cohort_idx])

            # Create team parameter type for this study/work group type
            class_info[cohort_idx] = Array{class_params,1}(undef,team_count)
            for team_itr = 1:team_count
                class_info[cohort_idx][team_itr] = class_params()
            end
        end
    end

    return class_sizes::Array{Array{Int64,1},1},
            class_info::Array{Array{class_params,1},1},
            class_ids::Array{Array{Int64,1},1},
            nodes_by_class::Array{Array{Array{Int64,1},1},1}
end

function  create_class_param_type!(student_info::Array{student_params,1},
                                        n_students::Int64,
                                        attendence_propns::Array{Float64,2},
                                        rng::MersenneTwister,
                                        nodes_by_class::Array{Array{Array{Int64,1},1},1},
                                        class_ids::Array{Array{Int64,1},1},
                                        cohort_ids::Array{Int64,1})

    for ii = 1:n_students
        # decide if returning to work
        would_attend_f2f_classes = Int64(rand(rng)<attendence_propns[cohort_ids[ii]])

        # Construct node parameter type for current individual
        student_info[ii] = student_params(would_attend_f2f_classes = would_attend_f2f_classes,
                                            cohort_ID = cohort_ids[ii],
                                            class_ID = class_ids[cohort_ids[ii]][1])

        # Add node to class if would attend in person
        if would_attend_f2f_classes == 1
            push!(nodes_by_class[cohort_ids[ii]][class_ids[cohort_ids[ii]][1]], ii)
        end
        deleteat!(class_ids[cohort_ids[ii]],1) # remove worker from top of list
    end

    return nothing
end

# Load departmental data from file and then allocate classes
function generate_classes_campus_only(n_students::Int64,
                                        class_generation_parameters::class_generation_params,
                                        RNGseed::Int64)

    # Set the RNG
    rng = MersenneTwister(RNGseed)

    # Get team params
    @unpack n_cohorts, attendence_propns, classtype_proportion,
            class_size_mean, class_size_sd = class_generation_parameters

    # Check correct number of cohorts is in use
    if n_cohorts != 84
        error("The value of n_cohorts, $n_cohorts, is not valid for use with generate_classes_campus_only. Must be set at 84.")
    end

    # Get number of departments in use
    n_depts = n_cohorts ÷ 3

    # Load in oncampus relevant data
    cohort_ids,
    student_numbers_per_cohort = get_oncampus_cohort_data(n_students,n_depts)

    # Sanity check. Number of records loaded should match n_students
    n_students_living_oncampus = length(cohort_ids)
    if n_students != n_students_living_oncampus
        error("Mismatch between the specified student population size (set as $n_students, but should be 7155) and number of data records loaded")
    end

    # Initialise and then populate places to store class sizes per cohort type
    # & class information
    class_sizes::Array{Array{Int64,1},1},
            class_info::Array{Array{class_params,1},1},
            class_ids::Array{Array{Int64,1},1},
            nodes_by_class::Array{Array{Array{Int64,1},1},1} = create_class_variables(n_cohorts,
                                                                                        student_numbers_per_cohort,
                                                                                        class_generation_parameters,
                                                                                        rng)

    # Create class parameter type variables
    student_info = Array{student_params,1}(undef,n_students)
    create_class_param_type!(student_info,
                                n_students,
                                attendence_propns,
                                rng,
                                nodes_by_class,
                                class_ids,
                                cohort_ids)

    return student_info::Array{student_params,1},
            class_sizes::Array{Array{Int64,1},1},
            class_info::Array{Array{class_params,1},1},
            nodes_by_class::Array{Array{Array{Int64,1},1},1}
end


# Load departmental data from file and then allocate classes
function generate_classes_all_students(n_students::Int64,
                                        class_generation_parameters::class_generation_params,
                                        RNGseed::Int64)

    # Set the RNG
    rng = MersenneTwister(RNGseed)

    # Get team params
    @unpack n_cohorts, attendence_propns, classtype_proportion,
            class_size_mean, class_size_sd = class_generation_parameters

    # Check correct number of cohorts is in use
    if n_cohorts != 84
        error("The value of n_cohorts, $n_cohorts, is not valid for use with generate_classes_campus_only. Must be set at 84.")
    end

    # Get number of departments in use
    n_depts = n_cohorts ÷ 3

    # Load in oncampus relevant data
    cohort_ids,
    student_numbers_per_cohort = get_oncampus_cohort_data(n_students,n_depts)

    # Assign cohort status to offcampus students
    n_students_living_oncampus = sum(student_numbers_per_cohort)
    offcampus_dept_dist = readdlm("../Data/offcampus_departmental_split_dist.csv",',')[:]
    csum_offcampus_dept_dist = cumsum(offcampus_dept_dist)  # Get cumulative sum, to be used for sampling from
    for student_itr = (n_students_living_oncampus+1):n_students

        # Draw from relevant distribution
        student_cohort_ID = draw_sample_from_pmf(csum_offcampus_dept_dist,
                                                    rng)

        # Update cohort ID vector & counts per cohort
        cohort_ids[student_itr] = student_cohort_ID
        student_numbers_per_cohort[student_cohort_ID] += 1
    end

    # Initialise and then populate places to store class sizes per cohort type
    # & class information
    class_sizes::Array{Array{Int64,1},1},
            class_info::Array{Array{class_params,1},1},
            class_ids::Array{Array{Int64,1},1},
            nodes_by_class::Array{Array{Array{Int64,1},1},1} = create_class_variables(n_cohorts,
                                                                                        student_numbers_per_cohort,
                                                                                        class_generation_parameters,
                                                                                        rng)

    # Create class parameter type variables
    student_info = Array{student_params,1}(undef,n_students)
    create_class_param_type!(student_info,
                                n_students,
                                attendence_propns,
                                rng,
                                nodes_by_class,
                                class_ids,
                                cohort_ids)

    return student_info::Array{student_params,1},
            class_sizes::Array{Array{Int64,1},1},
            class_info::Array{Array{class_params,1},1},
            nodes_by_class::Array{Array{Array{Int64,1},1},1}
end

"""
Functions to produce the initial network layers
"""
# ASSUMPTIONS:
# Student contacts are split into 4 settings: class, social (fixed), household & dynamic social
# If a class is inactive (no face-to-face classes), there are NO class setting contacts
# Social contacts (fixed) are for society/sports groups contacts - instead take a subset each day (can be a larger subset on weekends)
# Household contacts are within a house and occur everyday
# Dynamic contacts can occur with any other student
function generate_contacts_uni(n_students::Int64,
                                endtime::Int64,
                                network_parameters::network_params,
                                class_generation_params::class_generation_params,
                                nodes_by_class::Array{Array{Array{Int64,1},1},1},
                                RNGseed::Int64)
    # Inputs:
    # n_students - Number of workers in the system
    # endtime - Number of timesteps per simulation
    # network_parameters - Fields relating to the network generation
    # class_generation_params - parameters relating to student class generation
    # nodes_by_class - array storing ids of students within each class
    # RNGseed:Int64 - Value to seed the random number generator

# Outputs:
    # class_contacts, work_contacts_other_workplace, household_contacts, social_contacts
    #           - Vector of vectors with IDs of contacts
    # work_or_study_group_contacts_per_node, work_contacts_other_workplace_per_node,
    #       household_contacts_per_node, social_contacts_per_node - Total number of regular contacts within each setting
    # n_households - Total number of households in the system

    @unpack student_info, class_sizes,
        prob_social_contact,dd_within_class,
        class_degree_distribution, between_class_contact_probs, dynamic_social_contact_degree_distribution,
        household_size_distribution,
        class_info, society_info, social_group_size_distribution,
        friend_of_friend_prob, max_contacts_social, network_generation_method = network_parameters

    # Set the RNG
    rng = MersenneTwister(RNGseed)

    # Initialise vector of vectors storing IDs of contacts for each node
    # in their work or study group (non-covid guidance setting)
    class_contacts = Array{Array{Int64,1},1}(undef,n_students)
    cohort_contacts = Array{Array{Int64,1},1}(undef,n_students)

    # Initialise vector of vectors storing IDs of contacts for each node in social
    # and household settings
    n_societies = length(society_info)
    society_contacts = Array{Array{Int64,1},2}(undef,n_students,n_societies)

    for ii = 1:n_students
        class_contacts[ii] = Int64[]
        cohort_contacts[ii] = Int64[]

        for jj = 1:n_societies
            society_contacts[ii,jj] = Int64[]
        end
    end

    # Initialise vectors giving total contacts made by node
    work_or_study_group_contacts_per_node = zeros(Int64,n_students)
    cohort_contacts_per_node = zeros(Int64,n_students)
    society_contacts_per_node = zeros(Int64,n_students,n_societies)

    # Construct network links for study setting contacts
    if network_generation_method == "configuration" # Use configuration model

        # Construct links for classes
        if length(class_degree_distribution)==1
            class_degree_distribution= repeat(class_degree_distribution, class_generation_params.n_cohorts)
        end
        if length(between_class_contact_probs)==1
            between_class_contact_probs = repeat(between_class_contact_probs, class_generation_params.n_cohorts)
        end
        if length(dynamic_social_contact_degree_distribution)==1
            dynamic_social_contact_degree_distribution = repeat(dynamic_social_contact_degree_distribution, class_generation_params.n_cohorts)
        end

        # Cycle through cohorts and classes
        for cohort_idx = 1:class_generation_params.n_cohorts

            nodes_outside_cluster = collect(Iterators.flatten(nodes_by_class[cohort_idx]))

            degree_distribution = class_degree_distribution[cohort_idx]
            external_contact_prob = between_class_contact_probs[cohort_idx]

            for class_idx = 1:length(class_sizes[cohort_idx])

                nodes_within_cluster = nodes_by_class[cohort_idx][class_idx]

                if length(nodes_within_cluster) > 0
                    configuration_model!(student_info,
                                        external_contact_prob,
                                        degree_distribution,
                                        nodes_within_cluster,
                                        nodes_outside_cluster,
                                        class_contacts,
                                        cohort_contacts,
                                        work_or_study_group_contacts_per_node,
                                        cohort_contacts_per_node,
                                        RNGseed::Int64)
                end
            end
        end

        # # Get social group contacts
        # degree_distribution = social_group_size_distribution
        # configuration_model!(n_students,
        #                         friend_of_friend_prob,
        #                         degree_distribution,
        #                         social_contacts,
        #                         social_contacts_per_node,
        #                         max_contacts_social,
        #                         RNGseed)

    elseif network_generation_method == "ER" # Construct network using Erdos-Renyi model

        ER_model_uni!(student_info,
                        n_students,
                        class_sizes,
                        dd_within_class,
                        class_contacts,
                        cohort_contacts,
                        work_or_study_group_contacts_per_node,
                        cohort_contacts_per_node,
                        RNGseed)
    end

    # Construct network links for society contacts
    for ii = 1:(n_students-1)
        person_ii_society_idxs::Array{Int64,1} = student_info[ii].society_IDs

        for jj = (ii+1):n_students
            ### Check society contacts ###
            person_jj_society_idxs::Array{Int64,1} = student_info[jj].society_IDs
            if (length(person_ii_society_idxs) > 0) &&
                (length(person_jj_society_idxs) > 0) # No checks needed if either person not a member of any society

                # Iterate over each society person ii is a member of
                for society_itr = 1:length(person_ii_society_idxs)
                    person_ii_current_society_ID = person_ii_society_idxs[society_itr]

                    # Check if person jj is also a member of that society
                    # If so, apply uniform possibility of contacts with all others in society.
                    # More people gives more contacts.
                    for person_jj_society_itr = 1:length(person_jj_society_idxs)
                        person_jj_current_society_ID = person_jj_society_idxs[person_jj_society_itr]
                        if (person_jj_current_society_ID == person_ii_current_society_ID)
                            # Get relevant society contact probability
                            society_type_val = society_info[person_ii_current_society_ID].society_type
                            if rand(rng) < prob_social_contact[society_type_val]
                                # Assign IDs of contact to tuple vectors
                                push!(society_contacts[ii,person_ii_current_society_ID],jj)
                                push!(society_contacts[jj,person_jj_current_society_ID],ii)

                                # Increment number of contacts for nodes ii & jj
                                society_contacts_per_node[ii,person_ii_current_society_ID] += 1
                                society_contacts_per_node[jj,person_jj_current_society_ID] += 1
                            end
                        end
                    end
                end
            end
        end
    end

    # Define what is returned from the function
    contacts = contacts_struct(n_students=n_students,endtime=endtime,
        class_contacts = class_contacts,
        cohort_contacts = cohort_contacts,
        society_contacts = society_contacts)


    # Return outputs
    return contacts::contacts_struct,
            work_or_study_group_contacts_per_node::Array{Int64,1},
            cohort_contacts_per_node::Array{Int64,1},
            society_contacts_per_node::Array{Int64,2}
end

"""
Functions for assigning students to societies
"""

function assign_societies_one_per_student(student_info::Array{student_params,1},
                            society_generation_parameters::society_generation_params,
                            endtime::Int64,
                            RNGseed::Int64)

    # Set the RNG
    rng = MersenneTwister(RNGseed)

    # Create societies
    # This is arbitrary at this point
    # Will have 250 societies with up to 40 members in each (equating to 10,000 individuals)
    n_societies = 1
    society_member_limit = 40
    n_members_per_society = zeros(Int64,1000)  # Store number of members in society.
                                                #  As total number of societies not known, set to be a large number


    # Initialise vector of vectors to store IDs of members in that society
    society_member_lists = Array{Array{Int64,1},1}(undef,1000)
    society_member_lists[1] = Int64[] # Initialise vector for first society

    # Iterate over each student
    n_students = length(student_info)
    shuffled_student_order = Vector(1:n_students)
    shuffle!(rng, shuffled_student_order)
    n_society_members = 1 # Initialise counter
    for person_itr = 1:n_students

        # Get student ID
        student_ID = shuffled_student_order[person_itr]

        # Assign society ID to student record
        student_info[student_ID].society_IDs = [n_societies]

        # Add student ID to society membership list
        push!(society_member_lists[n_societies],student_ID)

        # Generate the next society if needed
        if n_society_members == society_member_limit
            # Assign number of members of society to storage vector
            n_members_per_society[n_societies] = n_society_members

            # Increment number of societies
            n_societies += 1
            society_member_lists[n_societies] = Int64[] # Initialise vector

            # Reset counter
            n_society_members = 1
        else # Society not yet full. Increment member counter
            n_society_members += 1
        end
    end

    # Create society parameter type
    society_info = Array{society_params,1}(undef,n_societies)
    for society_itr = 1:n_societies
        society_info[society_itr] = society_params(schedule = zeros(Int64,endtime),
                                                                society_type = 1,
                                                                n_members = n_members_per_society[society_itr],
                                                                member_list = society_member_lists[society_itr]
                                                                )
            # Initialises the schedule field & the society type field
    end

    # Return society parameter type
    return society_info::Array{society_params,1}
end

function assign_societies_from_aggregated_data(student_info::Array{student_params,1},
                            society_generation_parameters::society_generation_params,
                            endtime::Int64,
                            RNGseed::Int64)

    # Set the RNG
    rng = MersenneTwister(RNGseed)

    # Get the society data to be loaded
    # society_size_data = readdlm("../Data/society_size_data.txt")
        # Envisage list of IDs per student
    # societies_joined_per_person_dist = readdlm("../Data/societies_joined_per_person_dist.txt")
        # For each group, whether it comes under society or sports club
    # society_type = readdlm("../Data/type_by_society.txt")
            # For each group, whether it comes under society or sports club

    # Dummy values in place of data files
    society_size_data = rand(rng,10:10:100,335) # 65 sports clubs, 270 societies
    @unpack societies_joined_per_person_dist = society_generation_parameters

    # Set the total number of societies
    n_societies = length(society_size_data)

    # Initialise society type & then shuffle list of society types
    society_type_vec = ones(Int64,n_societies)
    for society_itr = 271:n_societies
        society_type_vec[society_itr] = 2
    end
    shuffle!(rng,society_type_vec)

    # Initialise vector of vectors to store IDs of members in that society
    society_member_lists = Array{Array{Int64,1},1}(undef,n_societies)
    for society_itr = 1:n_societies
        society_member_lists[society_itr] = Int64[] # Initialise vector for first society
    end

    # Draw the number of societies each student is a member of
    n_students = length(student_info)
    n_societies_joinable = zeros(Int64,n_students)
    student_assignable_status = zeros(Int64,n_students)
    csum_societies_joined_per_person_dist = cumsum(societies_joined_per_person_dist)
    for student_itr = 1:n_students
        # Draw sample. Idx offset of 1 as entry one of csum_societies_joined_per_person_dist
        # corresponds to zero societies joined
        n_societies_joinable[student_itr] = draw_sample_from_pmf(csum_societies_joined_per_person_dist,
                                        rng,
                                        idx_offset = 1)

        # Alter assignable status to 1 if assigned a positive value for n_societies_joinable
        if n_societies_joinable[student_itr] > 0
            student_assignable_status[student_itr] = 1
        end
    end

    # Initialise selectable_student_IDs
    n_selectable_students = sum(student_assignable_status)
    selectable_student_IDs = zeros(Int64,n_selectable_students)
    selectable_student_idx = 1  # Initialise counter for allocating selectable students
    for student_itr = 1:n_students
        if student_assignable_status[student_itr] == 1
            # Student can be assigned to society. Get the student ID.
            selectable_student_IDs[selectable_student_idx] = student_itr

            # Increment the allocation counter
            selectable_student_idx += 1
        end
    end


    # Populate each society. Decrease number of societies each student is able to join
    for society_itr = 1:n_societies
        # Get size of this society
        society_size = society_size_data[society_itr]

        # Initialise vector to store member IDs
        society_member_vec = zeros(Int64,society_size)

        # Draw IDs of students to be assigned to society
        # Assign those students to the society
        # n_selectable_students = length(selectable_student_IDs)
        for selectable_student_itr = 1:society_size

            # Select student from those that are assignable
            if n_selectable_students > 0
                r = rand(rng,1:n_selectable_students)
                selected_student_ID = selectable_student_IDs[r]

                # Assign student ID to society membership vector
                society_member_vec[selectable_student_itr] = selected_student_ID

                # Push society ID to the student info parameter
                push!(student_info[selected_student_ID].society_IDs,society_itr)

                # Decrease the number of societies student can be allocated to
                # If joinable count hits zero, update assignable status vector
                n_societies_joinable[selected_student_ID] -= 1
                if n_societies_joinable[selected_student_ID] == 0
                    student_assignable_status[selected_student_ID] = 0

                    # Remove student from selectable list
                    deleteat!(selectable_student_IDs,r)
                    n_selectable_students -= 1  # Reduce number of selectable students by one
                end
            end
        end

        # Populate society membership list
        society_member_lists[society_itr] = copy(society_member_vec)
    end

    # Create society parameter type
    society_info = Array{society_params,1}(undef,n_societies)
    for society_itr = 1:n_societies
        society_info[society_itr] = society_params(schedule = zeros(Int64,endtime),
                                                                society_type = society_type_vec[society_itr],
                                                                n_members = length(society_member_lists[society_itr]),
                                                                member_list = society_member_lists[society_itr]
                                                                )
            # Initialises the schedule field & the society type field
    end

    # Return society parameter type
    return society_info::Array{society_params,1}
end

function assign_societies_from_individual_data(student_info::Array{student_params,1},
                            society_generation_parameters::society_generation_params,
                            endtime::Int64,
                            RNGseed::Int64)

    # Get the society data to be loaded
    society_data_raw = readdlm("../Data/society_data.txt")
        # Envisage list of IDs per student
    society_type = readdlm("../Data/type_by_society.txt")
        # For each group, whether it comes under society or sports club

    # Get total number of societies
    n_societies = maximum(society_data_raw)

    # Get maximum number of entries in the society record data
    # Given by the number of columns of society_data_raw
    max_n_societies_joined = size(society_data_raw,2)
        # Envisage data such as: 1 2 3
        #                        2 10 14 24
        #                        3 ...

        # If loaded using readdlm, "empty" entries in the array loaded as ""

    # Initialise vector of vectors to store IDs of members in that society
    society_member_lists = Array{Array{Int64,1},1}(undef,n_societies)
    for society_itr = 1:n_societies
        society_member_lists[society_itr] = Int64[] # Initialise vector for first society
    end

    # Get number of students. Check it matches number of records in the society data
    n_students = length(student_info)
    n_records_in_data::Int64 = size(society_data_raw,1)
    if n_students != n_records_in_data
        error("Number of students included in simulation, $n_students, does not match number of records in the society data, $n_records_in_data.")
    end

    # Iterate over each student
    for person_itr = 1:n_students

        # Iterate over relevant line from the input data
        non_empty_entries = 0
        for entry_itr = 1:max_n_societies_joined
            if society_data_raw[person_itr,entry_itr] != ""  # Check not an empty entry

                # Get the relevant society ID
                society_ID = society_data_raw[person_itr,entry_itr]

                # Add student ID to society membership list
                push!(society_member_lists[society_ID],person_itr)

                # Increment counter for number of non-empty entries
                non_empty_entries += 1
            end
        end

        # Construct vector to store IDs of vectors student is a member of
        if non_empty_entries > 0
            society_joined_vec = zeros(Int64,non_empty_entries)
            for index_itr = 1:non_empty_entries
                society_joined_vec[index_itr] = society_data_raw[person_itr,index_itr]
            end

            # Assign to society field in student parameter type
            student_info[person_itr].society_IDs = copy(society_joined_vec)
        end
    end

    # Create society parameter type
    society_info = Array{society_params,1}(undef,n_societies)
    for society_itr = 1:n_societies
        society_info[society_itr] = society_params(schedule = zeros(Int64,endtime),
                                                                society_type = society_type[society_itr],
                                                                n_members = length(society_member_lists[society_itr]),
                                                                member_list = society_member_lists[society_itr]
                                                                )
            # Initialises the schedule field & the society type field
    end

    # Return society parameter type
    return society_info::Array{society_params,1}
end

"""
Functions for generating househould contacts
"""

# General structure:
# Inputs:
# n_students::Int64 - Number of individuals in the system
# contacts::contacts_struct - Parameter structure for contacts. Will be modifying the household_contacts field.
# network_parameters::network_params - Parameter structure with information relevant to network structure
# RNGseed::Int64 - Used to seed the random number generator

# Outputs:
# household_contacts_per_node::Array{Int64,1} - Per individual, the number of household contacts they have
# n_households::Int64 - Number of households in the system

### Fn to read in campus population from data.
# Use hierarchy (hall > block > floor > household)
function get_on_campus_households_from_data!(network_parameters::network_params)

    @unpack lowest_campus_denomination, student_info = network_parameters

    # Get the campus data to be loaded
    campus_data_raw = readdlm("../Data/campus_data.txt",Int64) # Hall > Block > Floor > Household > Kitchen

    # Initialise to ensure Int64 type
    # Number of rows of input data gives size of campus population
    n_students_living_oncampus::Int64 = size(campus_data_raw,1)
    campus_data = Array{Int64,2}(undef, n_students_living_oncampus, 4)

    # Fill n_students_on_campus x 4 matrix according to lowest denomination
    if lowest_campus_denomination == "hall"
        campus_data[:,:] = campus_data_raw[:,[1,1,1,1]]
    elseif lowest_campus_denomination == "block"
        campus_data[:,:] = campus_data_raw[:,[1,2,2,2]]
    elseif lowest_campus_denomination == "floor"
        campus_data[:,:] = campus_data_raw[:,[1,2,3,3]]
    elseif lowest_campus_denomination == "household"
        campus_data[:,:] = campus_data_raw[:,[1,2,3,4]]
    elseif lowest_campus_denomination == "kitchen"
        # Need to account for multiple kitchens per household
        campus_data[:,:] = campus_data_raw[:,[1,2,3,4]]
        n_kitchens = 1
        campus_data[1,4] = n_kitchens
        for ii = 2:length(campus_data[:,1])
            # If original kitchen or household ID changes, add 1 to new kitchen ID
            if (campus_data_raw[ii-1,5] != campus_data_raw[ii,5]) || (campus_data_raw[ii-1,4] != campus_data_raw[ii,4])
                n_kitchens += 1
            end
            # If floor ID changes, reset new kitchen ID
            if campus_data_raw[ii-1,3] != campus_data_raw[ii,3]
                n_kitchens = 1
            end
            # Assign new kitchen ID
            campus_data[ii,4] = n_kitchens
        end
    end

    # Initialise array to store household member list (Hall > Block > Floor > Household)
    household_member_list = Array{Array{Array{Array{Array{Int64,1},1},1},1},1}(undef,maximum(Int64, campus_data[:,1]))

    # Initialise array to store member list in each hierarchy level
    hall_member_list = Array{Array{Int64,1},1}(undef, maximum(Int64, campus_data[:,1]))
    block_member_list = Array{Array{Array{Int64,1},1},1}(undef, maximum(Int64, campus_data[:,1]))
    floor_member_list = Array{Array{Array{Array{Int64,1},1},1},1}(undef, maximum(Int64, campus_data[:,1]))

    # Assign the information on hall, block, floor, household
    for student_itr = 1:n_students_living_oncampus

        # Assign location to node ii
        student_info[student_itr].household_info.location::String = "Campus"   # name of location where student is living

        # Read each row of data from campus_data
        student_info[student_itr].household_info.hall_ID::Int64 = campus_data[student_itr,1]
        student_info[student_itr].household_info.block_ID::Int64 = campus_data[student_itr,2]
        student_info[student_itr].household_info.floor_ID::Int64 = campus_data[student_itr,3]
        student_info[student_itr].household_info.household_ID_within_block::Int64 = campus_data[student_itr,4]

        # Assign communal/ensuite bathroom information to variable
        # Bathroom key:
        # 1 - Communal bathroom
    	# 2 - Ensuite single
    	# 3 - Ensuite studio
    	# 4 - Twin bathroom
    	# 5 - Twin bathroom: studio
        bathroom_type_val = campus_data_raw[student_itr,6]
        if (bathroom_type_val == 2) || (bathroom_type_val == 3)
            # Ensuite type accomodation
            student_info[student_itr].household_info.ensuite_flag = true
        else
            # Communal bathroom type accomodation
            student_info[student_itr].household_info.ensuite_flag = false
        end

        # Assign blocks to hall
        if !isassigned(household_member_list, campus_data[student_itr,1])
            household_member_list[campus_data[student_itr,1]] = Array{Array{Array{Array{Int64,1},1},1},1}(undef,maximum(Int64, campus_data[campus_data[:,1].==campus_data[student_itr,1],2]))
            block_member_list[campus_data[student_itr,1]] = Array{Array{Int64,1},1}(undef,maximum(Int64, campus_data[campus_data[:,1].==campus_data[student_itr,1],2]))
            floor_member_list[campus_data[student_itr,1]] = Array{Array{Array{Int64,1},1},1}(undef,maximum(Int64, campus_data[campus_data[:,1].==campus_data[student_itr,1],2]))
        end
        # Assign floors to block
        if !isassigned(household_member_list[campus_data[student_itr,1]], campus_data[student_itr,2])
            household_member_list[campus_data[student_itr,1]][campus_data[student_itr,2]] = Array{Array{Array{Int64,1},1},1}(undef,maximum(Int64, campus_data[(campus_data[:,1].==campus_data[student_itr,1]).&(campus_data[:,2].==campus_data[student_itr,2]),3]))
            floor_member_list[campus_data[student_itr,1]][campus_data[student_itr,2]] = Array{Array{Int64,1},1}(undef,maximum(Int64, campus_data[(campus_data[:,1].==campus_data[student_itr,1]).&(campus_data[:,2].==campus_data[student_itr,2]),3]))
        end
        # Assign households to floor
        if !isassigned(household_member_list[campus_data[student_itr,1]][campus_data[student_itr,2]], campus_data[student_itr,3])
            household_member_list[campus_data[student_itr,1]][campus_data[student_itr,2]][campus_data[student_itr,3]] = Array{Array{Int64,1},1}(undef,maximum(Int64, campus_data[(campus_data[:,1].==campus_data[student_itr,1]).&(campus_data[:,2].==campus_data[student_itr,2]).&(campus_data[:,3].==campus_data[student_itr,3]),4]))
        end
        # Initialise households
        if !isassigned(household_member_list[campus_data[student_itr,1]][campus_data[student_itr,2]][campus_data[student_itr,3]], campus_data[student_itr,4])
            household_member_list[campus_data[student_itr,1]][campus_data[student_itr,2]][campus_data[student_itr,3]][campus_data[student_itr,4]] = Int64[]
        end

        # Add students to household
        push!(household_member_list[campus_data[student_itr,1]][campus_data[student_itr,2]][campus_data[student_itr,3]][campus_data[student_itr,4]], student_itr)
    end

    return household_member_list::Array{Array{Array{Array{Array{Int64,1},1},1},1},1},
            hall_member_list::Array{Array{Int64,1},1},
            block_member_list::Array{Array{Array{Int64,1},1},1},
            floor_member_list::Array{Array{Array{Array{Int64,1},1},1},1},
            n_students_living_oncampus::Int64
end


### Fn to assign oncampus households by cohort
function get_on_campus_households_by_cohort!(network_parameters::network_params)

    @unpack lowest_campus_denomination, student_info = network_parameters

    """
    Use the campus accom structure data file
    Also take the deparmental allocation. Sort ascending. Allocate in order.
    """

    # Use the campus accom structure data file
    # Get the campus data to be loaded
    campus_data_raw = readdlm("../Data/campus_data.txt",Int64) # Hall > Block > Floor > Household > Kitchen

    # Initialise to ensure Int64 type
    # Number of rows of input data gives size of campus population
    n_students_living_oncampus::Int64 = size(campus_data_raw,1)
    campus_data = Array{Int64,2}(undef, n_students_living_oncampus, 4)

    # Fill n_students_on_campus x 4 matrix according to lowest denomination
    if lowest_campus_denomination == "hall"
        campus_data[:,:] = campus_data_raw[:,[1,1,1,1]]
    elseif lowest_campus_denomination == "block"
        campus_data[:,:] = campus_data_raw[:,[1,2,2,2]]
    elseif lowest_campus_denomination == "floor"
        campus_data[:,:] = campus_data_raw[:,[1,2,3,3]]
    elseif lowest_campus_denomination == "household"
        campus_data[:,:] = campus_data_raw[:,[1,2,3,4]]
    elseif lowest_campus_denomination == "kitchen"
        # Need to account for multiple kitchens per household
        campus_data[:,:] = campus_data_raw[:,[1,2,3,4]]
        n_kitchens = 1
        campus_data[1,4] = n_kitchens
        for ii = 2:length(campus_data[:,1])
            # If original kitchen or household ID changes, add 1 to new kitchen ID
            if (campus_data_raw[ii-1,5] != campus_data_raw[ii,5]) || (campus_data_raw[ii-1,4] != campus_data_raw[ii,4])
                n_kitchens += 1
            end
            # If floor ID changes, reset new kitchen ID
            if campus_data_raw[ii-1,3] != campus_data_raw[ii,3]
                n_kitchens = 1
            end
            # Assign new kitchen ID
            campus_data[ii,4] = n_kitchens
        end
    end

    # Initialise array to store household member list (Hall > Block > Floor > Household)
    household_member_list = Array{Array{Array{Array{Array{Int64,1},1},1},1},1}(undef,maximum(Int64, campus_data[:,1]))

    # Initialise array to store member list in each hierarchy level
    hall_member_list = Array{Array{Int64,1},1}(undef, maximum(Int64, campus_data[:,1]))
    block_member_list = Array{Array{Array{Int64,1},1},1}(undef, maximum(Int64, campus_data[:,1]))
    floor_member_list = Array{Array{Array{Array{Int64,1},1},1},1}(undef, maximum(Int64, campus_data[:,1]))


    # Also take the deparmental allocation.
    # Map study type (undergrad or postgrad value) to cohort ID.
    # Then sort ascending.
    dept_data_oncampus = readdlm("../Data/department_student_data.txt",Int64)
    n_depts = 28  # Number of departments in use with this dataset
    student_cohort_ID = zeros(Int64,n_students_living_oncampus)
    for student_itr = 1:n_students_living_oncampus
        # Map to cohort ID value. Determined by whether it study type is
        # undergraduate or postgraduate
        if dept_data_oncampus[student_itr,2] == 1 # First year undergrad
            student_cohort_ID[student_itr] = dept_data_oncampus[student_itr,1]
        elseif dept_data_oncampus[student_itr,2] == 2 # Non-first year undergrad
            student_cohort_ID[student_itr] = dept_data_oncampus[student_itr,1] + n_depts
        elseif dept_data_oncampus[student_itr,2] == 3 # Postgrad
            student_cohort_ID[student_itr] = dept_data_oncampus[student_itr,1] + (2*n_depts)
        end
    end
    sorted_order_student_cohort_ID = sortperm(student_cohort_ID) # Student IDs, sorted by ascending student_cohort_ID

    # Assign the information on hall, block, floor, household
    for student_itr = 1:n_students_living_oncampus

        # Get student ID from sorted list
        current_itr_student_ID = sorted_order_student_cohort_ID[student_itr]

        # Assign location to node ii
        student_info[current_itr_student_ID].household_info.location::String = "Campus"   # name of location where student is living

        # Read each row of data from campus_data
        student_info[current_itr_student_ID].household_info.hall_ID::Int64 = campus_data[student_itr,1]
        student_info[current_itr_student_ID].household_info.block_ID::Int64 = campus_data[student_itr,2]
        student_info[current_itr_student_ID].household_info.floor_ID::Int64 = campus_data[student_itr,3]
        student_info[current_itr_student_ID].household_info.household_ID_within_block::Int64 = campus_data[student_itr,4]

        # Assign communal/ensuite bathroom information to variable
        # Bathroom key:
        # 1 - Communal bathroom
    	# 2 - Ensuite single
    	# 3 - Ensuite studio
    	# 4 - Twin bathroom
    	# 5 - Twin bathroom: studio
        bathroom_type_val = campus_data_raw[student_itr,6]
        if (bathroom_type_val == 2) || (bathroom_type_val == 3)
            # Ensuite type accomodation
            student_info[student_itr].household_info.ensuite_flag = true
        else
            # Communal bathroom type accomodation
            student_info[student_itr].household_info.ensuite_flag = false
        end

        # Assign blocks to hall
        if !isassigned(household_member_list, campus_data[student_itr,1])
            household_member_list[campus_data[student_itr,1]] = Array{Array{Array{Array{Int64,1},1},1},1}(undef,maximum(Int64, campus_data[campus_data[:,1].==campus_data[student_itr,1],2]))
            block_member_list[campus_data[student_itr,1]] = Array{Array{Int64,1},1}(undef,maximum(Int64, campus_data[campus_data[:,1].==campus_data[student_itr,1],2]))
            floor_member_list[campus_data[student_itr,1]] = Array{Array{Array{Int64,1},1},1}(undef,maximum(Int64, campus_data[campus_data[:,1].==campus_data[student_itr,1],2]))
        end
        # Assign floors to block
        if !isassigned(household_member_list[campus_data[student_itr,1]], campus_data[student_itr,2])
            household_member_list[campus_data[student_itr,1]][campus_data[student_itr,2]] = Array{Array{Array{Int64,1},1},1}(undef,maximum(Int64, campus_data[(campus_data[:,1].==campus_data[student_itr,1]).&(campus_data[:,2].==campus_data[student_itr,2]),3]))
            floor_member_list[campus_data[student_itr,1]][campus_data[student_itr,2]] = Array{Array{Int64,1},1}(undef,maximum(Int64, campus_data[(campus_data[:,1].==campus_data[student_itr,1]).&(campus_data[:,2].==campus_data[student_itr,2]),3]))
        end
        # Assign households to floor
        if !isassigned(household_member_list[campus_data[student_itr,1]][campus_data[student_itr,2]], campus_data[student_itr,3])
            household_member_list[campus_data[student_itr,1]][campus_data[student_itr,2]][campus_data[student_itr,3]] = Array{Array{Int64,1},1}(undef,maximum(Int64, campus_data[(campus_data[:,1].==campus_data[student_itr,1]).&(campus_data[:,2].==campus_data[student_itr,2]).&(campus_data[:,3].==campus_data[student_itr,3]),4]))
        end
        # Initialise households
        if !isassigned(household_member_list[campus_data[student_itr,1]][campus_data[student_itr,2]][campus_data[student_itr,3]], campus_data[student_itr,4])
            household_member_list[campus_data[student_itr,1]][campus_data[student_itr,2]][campus_data[student_itr,3]][campus_data[student_itr,4]] = Int64[]
        end

        # Add students to household
        push!(household_member_list[campus_data[student_itr,1]][campus_data[student_itr,2]][campus_data[student_itr,3]][campus_data[student_itr,4]], current_itr_student_ID)
    end

    return household_member_list::Array{Array{Array{Array{Array{Int64,1},1},1},1},1},
            hall_member_list::Array{Array{Int64,1},1},
            block_member_list::Array{Array{Array{Int64,1},1},1},
            floor_member_list::Array{Array{Array{Array{Int64,1},1},1},1},
            n_students_living_oncampus::Int64
end

# Populate hall, block, floor, household lists
function populate_oncampus_contact_lists!(global_household_id::Int64,
                                            hall_member_list::Array{Array{Int64,1},1},
                                            block_member_list::Array{Array{Array{Int64,1},1},1},
                                            floor_member_list::Array{Array{Array{Array{Int64,1},1},1},1},
                                            household_member_list::Array{Array{Array{Array{Array{Int64,1},1},1},1},1},
                                            household_contacts::Array{Array{Int64,1},1},
                                            household_contacts_per_node::Array{Int64,1},
                                            student_info::Array{student_params,1}
                                            )

    for hall_id = 1:length(household_member_list)
        hall_member_list[hall_id] = flattenall(household_member_list[hall_id])
        for block_id = 1:length(household_member_list[hall_id])
            block_member_list[hall_id][block_id] = flattenall(household_member_list[hall_id][block_id])
            for floor_id = 1:length(household_member_list[hall_id][block_id])
                floor_member_list[hall_id][block_id][floor_id] = flattenall(household_member_list[hall_id][block_id][floor_id])
                for household_id = 1:length(household_member_list[hall_id][block_id][floor_id])

                    # Increment overall number of households
                    global_household_id += 1

                    for ii = 1:length(household_member_list[hall_id][block_id][floor_id][household_id])
                        # ID of student in household
                        student_id = household_member_list[hall_id][block_id][floor_id][household_id][ii]
                        # Update student's global household ID
                        student_info[student_id].household_info.household_ID::Int64 = global_household_id
                        # Add links to everyone else in household, except self
                        household_contacts[student_id] = household_member_list[hall_id][block_id][floor_id][household_id][1:end .!= ii]
                        household_contacts_per_node[student_id] = length(household_member_list[hall_id][block_id][floor_id][household_id]) - 1
                    end
                end
            end
        end
    end

    #Return number of households generated so far
    n_households_oncampus = global_household_id
    return n_households_oncampus::Int64
end

### Fn to perform household assignment for students resident off-campus
function assign_off_campus_households!(offcampus_student_household_size_distribution::Array{Float64,1},
                                        offcampus_student_household_location_dist::Array{Float64,1},
                                        offcampus_student_location_names::Array{String,1},
                                        n_students_living_oncampus::Int64,
                                        n_students::Int64,
                                        global_household_id::Int64,
                                        rng::MersenneTwister,
                                        household_contacts::Array{Array{Int64,1},1},
                                        household_contacts_per_node::Array{Int64,1},
                                        student_info::Array{student_params,1};
                                        assign_strat::String = "None"
                                        )
     # Get cumulative sum of offcampus household size & offcampus household location,
     # used in each repeat of while loop
    csum_offcampus_household_size_distribution = cumsum(offcampus_student_household_size_distribution)
    csum_offcampus_student_household_location_dist = cumsum(offcampus_student_household_location_dist)

    # Get starting index number for off campus students
    # Initialise variables to be incremented in while loop
    person_id = n_students_living_oncampus + 1
    global_household_id += 1 # Move on household_id to the next value
    # Iterate over off campus students
    while person_id <= n_students
        # generate random household size
        household_size = findfirst(x->x>rand(rng), csum_offcampus_household_size_distribution)

        # check we don't exceed remaining workers left
        household_size = min(household_size, n_students-person_id+1)

        # Draw offcampus location from distribution
        household_location_ID = draw_sample_from_pmf(csum_offcampus_student_household_location_dist,
                                                        rng)

        # create fully-connected household
        for ii = 0:(household_size-1)

            # Get index to access array at. Depends on assignment strategy
            if assign_strat == "None"
                person_ii_ID = person_id+ii
            elseif assign_strat == "Cohort"
                person_ii_ID = students_ordered_by_cohort[person_id+ii]
            else
                error("assign_strat has been passed an invalid string, $assign_strat. Please check.")
            end

            for jj = 0:(household_size-1)

                if assign_strat == "None"
                    person_jj_ID = person_id+jj
                elseif assign_strat == "Cohort"
                    person_jj_ID = students_ordered_by_cohort[person_id+jj]
                end

                # Get household contacts
                if ii != jj
                    push!(household_contacts[person_ii_ID], person_jj_ID)
                    household_contacts_per_node[person_ii_ID] += 1
                end
            end

            # Assign global household ID to node ii
            student_info[person_ii_ID].household_info.household_ID = global_household_id

            # Amend other household info fields as needed
            # Leaves on-campus fields (hall, block, floor, household within block) at default values of 0
            student_info[person_ii_ID].household_info.on_campus_accom = false   # Whether accomodation is on campus (true) or off campus (false)
            student_info[person_ii_ID].household_info.location = offcampus_student_location_names[household_location_ID] # name of location where student is living
        end

        # Increment counter variables
        person_id += household_size
        global_household_id += 1
    end

    # Return total number of househoulds
    return global_household_id::Int64
end


### Overall function for assigning househoulds
# Based on subfunction name, certain sections of the function are called
function assign_households_fn!(generate_student_households_fn_name::String,
                                n_students::Int64,
                                contacts::contacts_struct,
                                network_parameters::network_params,
                                RNGseed::Int64)

    @unpack n_students_on_campus, lowest_campus_denomination, offcampus_student_household_size_distribution,
            offcampus_student_household_location_dist, offcampus_student_location_names,
                student_info, household_size_distribution = network_parameters

    # Set the RNG
    rng = MersenneTwister(RNGseed)

    # Initialise vector of vectors storing IDs of contacts for each node in social
    # and household settings
    household_contacts = Array{Array{Int64,1},1}(undef,n_students)
    for ii = 1:n_students
        household_contacts[ii] = Int64[]
    end

    # Initialise vectors giving total contacts made by node
    household_contacts_per_node = zeros(Int64,n_students)

     if (generate_student_households_fn_name == "assign_households_no_hierarchy")
         # Create household contacts & assign a household ID
         csum_household_size_distribution = cumsum(household_size_distribution) # Get cumulative sum, used in each repeat of while loop
         person_id = 1     # Initialise variables to be incremented in while loop
         household_id = 1
         while person_id <= n_students
             household_size = findfirst(x->x>rand(rng), csum_household_size_distribution)  # generate random household size

             # check we don't exceed remaining workers left
             household_size = min(household_size, n_students-person_id+1)

             # create fully-connected household
             for ii = 0:(household_size-1)
                 for jj = 0:(household_size-1)
                     # Get household contacts
                     if ii != jj
                         push!(household_contacts[person_id+ii], person_id+jj)
                         household_contacts_per_node[person_id+ii] += 1
                     end
                 end

                 # Assign household ID to node ii
                 student_info[person_id+ii].household_info.household_ID = household_id
             end

             person_id += household_size
             household_id += 1
         end

         # Assign number of households in total to output variable
         n_households = household_id

         # Update fields in contacts parameter type
         contacts.household_contacts = household_contacts
     else
        # Assignment to households uses accomodation hierarchy

        """
        On campus household assignment
        """
        if (generate_student_households_fn_name == "assign_households_campus_only") ||
            (generate_student_households_fn_name == "assign_households_all_students")
            household_member_list,
            hall_member_list,
            block_member_list,
            floor_member_list,
            n_students_living_oncampus = get_on_campus_households_from_data!(network_parameters)
        elseif (generate_student_households_fn_name == "assign_households_by_cohort")
            household_member_list,
            hall_member_list,
            block_member_list,
            floor_member_list,
            n_students_living_oncampus = get_on_campus_households_by_cohort!(network_parameters)
                # Use the campus accom structure data file
                # Also take the deparmental allocation. Sort ascending. Allocate in order.
        end

        # Create fully-connected household & assign global household_ID
        # Use household_member_list
        global_household_id = 0
        n_oncampus_households = populate_oncampus_contact_lists!(global_household_id,
                                                                hall_member_list,
                                                                block_member_list,
                                                                floor_member_list,
                                                                household_member_list,
                                                                household_contacts,
                                                                household_contacts_per_node,
                                                                student_info
                                                                )

        """
        Off campus household assignment
        """
        n_households = 0
        if  (generate_student_households_fn_name == "assign_households_all_students") ||
                (generate_student_households_fn_name == "assign_households_by_cohort")
            # Create household contacts & assign a household ID
            # Assign number of households in total to output variable
            n_households = assign_off_campus_households!(offcampus_student_household_size_distribution,
                                            offcampus_student_household_location_dist,
                                            offcampus_student_location_names,
                                            n_students_living_oncampus,
                                            n_students,
                                            n_oncampus_households,
                                            rng,
                                            household_contacts,
                                            household_contacts_per_node,
                                            student_info
                                            )
        elseif (generate_student_households_fn_name == "assign_households_campus_only")
            n_households = n_oncampus_households
        end

        # Update fields in contacts parameter type
        contacts.household_contacts = household_contacts
        contacts.household_member_list = household_member_list
        contacts.hall_member_list = hall_member_list
        contacts.block_member_list = block_member_list
        contacts.floor_member_list = floor_member_list


        # Get IDs of students resident on-campus and off-campus assigned to vectors
        on_campus_student_IDs = zeros(Int64,n_students_on_campus)
        n_students_off_campus = n_students-n_students_on_campus
        off_campus_student_IDs = zeros(Int64,n_students_off_campus)
        on_campus_idx = 0; off_campus_idx = 0; # Set up values to access ID vectors
        for student_itr = 1:n_students
            # Check whether on-campus or off-campus resident
            if student_info[student_itr].household_info.on_campus_accom == true
                # Increment index counter and sanity check
                on_campus_idx += 1
                if on_campus_idx > n_students_on_campus
                    error("Error setting up on_campus_student_IDs. Exceeded length of allocated vector.")
                end

                # Allocate student ID to vector
                on_campus_student_IDs[on_campus_idx] = student_itr
            elseif student_info[student_itr].household_info.on_campus_accom == false
                # Increment index counter and sanity check
                off_campus_idx += 1
                if off_campus_idx > n_students_off_campus # Sanity check
                    error("Error setting up off_campus_student_IDs. Exceeded length of allocated vector.")
                end

                # Allocate student ID to vector
                off_campus_student_IDs[off_campus_idx] = student_itr
            else
                error("student_info[$student_itr].household_info.on_campus_accom has invalid value, $(student_info[student_itr].household_info.on_campus_accom).")
            end
        end

        # Update relevant field in network_params type
        network_parameters.on_campus_student_IDs = on_campus_student_IDs
        network_parameters.off_campus_student_IDs = off_campus_student_IDs
    end
    return household_contacts_per_node::Array{Int64,1},
                n_households::Int64
end

"""
Functions for use with configuration model
"""

#### Class/workplace version ####
function configuration_model!(student_info::Array{student_params,1},
                                external_contact_prob::Float64,
                                degree_distribution::Distribution,
                                nodes_within_cluster::Array{Int64,1},
                                nodes_outside_cluster::Array{Int64,1},
                                class_contacts::Array{Array{Int64,1},1},
                                cohort_contacts::Array{Array{Int64,1},1},
                                class_contacts_per_node::Array{Int64,1},
                                cohort_contacts_per_node::Array{Int64,1},
                                RNGseed::Int64)

    # Set the RNG
    rng = MersenneTwister(RNGseed)

    n_nodes = length(nodes_within_cluster)
    n_nodes_external = length(nodes_outside_cluster)

    edges_per_node = Distributions.rand(rng, degree_distribution, n_nodes)

    # Round degree to nearest whole number
    # Check if each node has already made contacts from other clusters
    # Limit maximum number of contacts to clustersize - 1
    for node_id = 1:n_nodes
        edges_per_node[node_id] = round(Int64, edges_per_node[node_id])
        edges_per_node[node_id] -= length(cohort_contacts[nodes_within_cluster[node_id]])
        edges_per_node[node_id] = min((n_nodes - 1), edges_per_node[node_id])
    end

    half_edges = cumsum(edges_per_node)

    n_stubs = half_edges[end]

    edges_within_group = zeros(Int64, n_nodes, n_nodes)
    [edges_within_group[i,i] = 1 for i=1:n_nodes]

    while half_edges[end] > 1

        stub_id1 = round(Int64, rand(rng)*(half_edges[end]-1) + 1)
        node_id1 = findfirst(x -> x >= stub_id1, half_edges)

        # This ensures no self-links or repeated edges within cluster
        nodes_remaining = findall(edges_within_group[node_id1,:] .== 0)
        half_edges_remaining = cumsum(edges_per_node[nodes_remaining])
        nodes_remaining = nodes_remaining[half_edges_remaining.>0]
        half_edges_remaining = half_edges_remaining[half_edges_remaining.>0]

        # Contact made external to cluster
        if (length(nodes_remaining) < 1) || (rand(rng) < external_contact_prob)

            node_id2 = round(Int64, rand(rng)*(n_nodes_external-1) + 1)

            # Don't allow links within cluster
            while node_id2 in nodes_within_cluster
                node_id2 = round(Int64, rand(rng)*(n_nodes_external-1) + 1)
            end

            if (student_info[nodes_within_cluster[node_id1]].would_attend_f2f_classes == 1) &&
                (student_info[nodes_outside_cluster[node_id2]].would_attend_f2f_classes == 1)
                push!(cohort_contacts[nodes_within_cluster[node_id1]], nodes_outside_cluster[node_id2])
                push!(cohort_contacts[nodes_outside_cluster[node_id2]], nodes_within_cluster[node_id1])

                cohort_contacts_per_node[nodes_within_cluster[node_id1]] += 1
                cohort_contacts_per_node[nodes_outside_cluster[node_id2]] += 1
            end

            half_edges[node_id1:end] .-= 1

        # Contact made within cluster
        else
            stub_id2 = round(Int64, rand(rng)*(half_edges_remaining[end]-1) + 1)
            node_id2 = nodes_remaining[findfirst(x -> x >= stub_id2, half_edges_remaining)]

            if (student_info[nodes_within_cluster[node_id1]].would_attend_f2f_classes == 1) &&
                (student_info[nodes_within_cluster[node_id2]].would_attend_f2f_classes == 1)
                push!(class_contacts[nodes_within_cluster[node_id1]], nodes_within_cluster[node_id2])
                push!(class_contacts[nodes_within_cluster[node_id2]], nodes_within_cluster[node_id1])

                class_contacts_per_node[nodes_within_cluster[node_id1]] += 1
                class_contacts_per_node[nodes_within_cluster[node_id2]] += 1
            end

            edges_within_group[node_id1,node_id2] += 1

            half_edges[node_id1:end] .-= 1
            half_edges[node_id2:end] .-= 1
        end
    end
end

#### Social version (friends-of-friends model) ####
function configuration_model!(n_students::Int64,
                                friend_of_friend_prob::Float64,
                                degree_distribution::Distribution,
                                social_contacts::Array{Array{Int64,1},1},
                                social_contacts_per_node::Array{Int64,1},
                                max_contacts_social::Int64,
                                RNGseed::Int64)

    # Set the RNG
    rng = MersenneTwister(RNGseed)

    n_nodes = n_students

    edges_per_node = Distributions.rand(rng, degree_distribution, n_nodes)

    # Round degree to nearest whole number
    # Limit maximum number of contacts to max_contacts_social
    for node_id = 1:n_nodes
        edges_per_node[node_id] = round(Int64, edges_per_node[node_id])
        edges_per_node[node_id] = min(max_contacts_social, edges_per_node[node_id])
    end

    half_edges = cumsum(edges_per_node)

    if half_edges[end] % 2 != 0
        half_edges[end] += 1
    end

    n_stubs = half_edges[end]

    edges_remaining_per_node = copy(edges_per_node)

    while half_edges[end] > 1

        stub_id1 = round(Int64, rand(rng)*(half_edges[end]-1) + 1)
        node_id1 = findfirst(x -> x >= stub_id1, half_edges)

        # Find friends of friends with edges remaining
        fof = collect(Iterators.flatten(social_contacts[social_contacts[node_id1]]))
        fof = fof[edges_remaining_per_node[fof].>0]

        # Remove self and anyone already friends with
        fof = fof[fof .!= node_id1]
        fof = fof[[fof[i] ∉ social_contacts[node_id1] for i=1:length(fof)]]

        # Contact made with friend of friend
        if (length(fof) > 0) & (rand(rng) < friend_of_friend_prob)

            node_id2 = fof[ceil(Int64, rand(rng)*length(fof))]

        # Contact made with not friend of friend
        else

            stub_id2 = round(Int64, rand(rng)*(half_edges[end]-1) + 1)
            node_id2 = findfirst(x -> x >= stub_id2, half_edges)

        end

        push!(social_contacts[node_id1], node_id2)
        push!(social_contacts[node_id2], node_id1)

        social_contacts_per_node[node_id1] += 1
        social_contacts_per_node[node_id2] += 1

        edges_remaining_per_node[node_id1] -= 1
        edges_remaining_per_node[node_id2] -= 1

        half_edges[node_id1:end] .-= 1
        half_edges[node_id2:end] .-= 1
    end
end

# Premake dynamic worker contacts,
# to be loaded in ahead of simulation
function generate_dynamic_student_contacts(RNGseed::Int64,
                                            n_students::Int64,
                                            endtime::Int64,
                                            student_info::Array{student_params,1},
                                            dynamic_social_contact_degree_distribution::Array{Distribution,2},
                                            max_contacts_social_dynamic::Int64)
# Inputs:
# RNGseed - Seed the random number generator
# n_students - Number of students in the system
# endtime - Number of timesteps simulation will run for
# student_info - entry for each person.
# dynamic_social_contact_degree_distribution - Distribution properties for students dynamic social contacts
# max_contacts_social_dynamic - Maxmimum permitted number of daily dynamic social contacts

# Outputs:
# dynamic_social_contacts - Per node, a record of dynamic social contacts made on each day

    """
    Set the RNG
    """
    rng = MersenneTwister(RNGseed)

    """
    Initialise vector of vectors storing IDs of contacts for each node
    """
    dynamic_social_contacts = Array{Array{Int64,1},2}(undef,endtime,n_students)

    """
    Iterate over all nodes
    Assign dynamic social contacts for each timestep
    """
    for student_itr = 1:n_students
        # Get dynamic group type for student_itr
        node_dynamic_grp_ID = student_info[student_itr].cohort_ID

        for time_itr = 1:endtime
            # Generate number of dynamic contacts from appropriate distribution
            gg = round(Int64, Distributions.rand(rng, dynamic_social_contact_degree_distribution[node_dynamic_grp_ID]))

            # Limit to specified maximum
            gg = min(gg, max_contacts_social_dynamic)

            # If dynamic worker contacts made, assign to output variable
            if gg > 0
                dynamic_social_contacts[time_itr,student_itr] = zeros(gg)

                # Generate required number of contacts
                for contact_itr = 1:gg
                    gg1 = ceil(Int64,rand(rng)*n_students) # Get IDs of nodes connected by dynamic link on current timestep
                    while gg1 == student_itr # Redraw if returned the index node
                        gg1 = ceil(Int64,rand(rng)*n_students) # Get IDs of nodes connected by dynamic link on current timestep
                    end
                    dynamic_social_contacts[time_itr,student_itr][contact_itr] = gg1
                end
            else # No dynamic worker contacts made on given timestep
                dynamic_social_contacts[time_itr,student_itr] = Int64[]
            end
        end
    end

    return dynamic_social_contacts::Array{Array{Int64,1},2}
end

"""
Functions used to generate dynamic accomodation contacts
"""
# Function to get dynamic accomodation contacts for each layer of accomodation hierarchy
function get_layer_dynamic_accom_contacts!(dynamic_accomodation_contacts::Array{Array{Int64,1},3},
                                            RNGseed::Int64,
                                            contact_prob::Float64,
                                            persons_to_check::Array{Int64,1},
                                            accom_hierarchy_name::String,
                                            endtime::Int64,
                                            student_ID::Int64
                                            )
# Inputs:
# dynamic_accomodation_contacts::Array{Array{Int64,1},3} - Per node, a record of dynamic accomodation contacts made on each day
#                                                           Relevant to those in on-campus accomodation only.
#                                                           Three dimensions of vectors. Row per timestep, column per accom level, slice per student
# RNGseed - Seed the random number generator
# contact_prob::Float64 - Chance of a contact occurring with each other person per timestep
# persons_to_check::Array{Int64,1} - Persons from the accomodation unit that may be contacted
# accom_hierarchy_string::Int64 - "Hall", "Block" or "Floor", dependent on where dynamic contacts are taking place
# endtime::Int64 - Number of timesteps performed in the simulation
# student_ID::Int64 - As decribed

# Outputs:
# Directly alters dynamic_accomodation_contacts

    """
    Set the RNG
    """
    rng = MersenneTwister(RNGseed)

    """
    Assign column entry of dynamic_accomodation_contacts based on accom_hierarchy_name
    """
    if accom_hierarchy_name == "Hall"
        accom_hierarchy_level = 1
    elseif accom_hierarchy_name == "Block"
        accom_hierarchy_level = 2
    elseif accom_hierarchy_name == "Floor"
        accom_hierarchy_level = 3
    else
        error("Illegal string provided for accom_hierarchy_name. Must be Hall, Block or Floor.")
    end

    """
    Check contacts in that accomodation unit. Push to vector if contact probability accepted
    """
    n_persons_to_check = length(persons_to_check)
    track_contact_status = zeros(Int64,n_persons_to_check)
    for time_itr = 1:endtime
        lmul!(0,track_contact_status)  # Reinitialise contacts made tracker
        for persons_to_check_itr = 1:n_persons_to_check
            if rand(rng) < contact_prob
                track_contact_status[persons_to_check_itr] = 1
            end
        end

        """
        Intialise & then populate the contacts vector
        """
        n_contacts_made = sum(track_contact_status)
        allocation_idx = 1  # Counter to allocate the ID of contact to vector in specified index position
        dynamic_accomodation_contacts[time_itr,accom_hierarchy_level,student_ID] = zeros(Int64,n_contacts_made)
        for persons_to_check_itr = 1:n_persons_to_check
            if track_contact_status[persons_to_check_itr] == 1
                person_contacted_ID = persons_to_check[persons_to_check_itr]
                dynamic_accomodation_contacts[time_itr,accom_hierarchy_level,student_ID][allocation_idx] = person_contacted_ID
                allocation_idx += 1 # Increment allocation index
            end
        end
    end
end

# Premake dynamic contacts within halls of residence,
# to be loaded in ahead of simulation
function generate_dynamic_campus_accom_contacts(RNGseed::Int64,
                                                n_students::Int64,
                                                endtime::Int64,
                                                student_info::Array{student_params,1},
                                                network_parameters::network_params,
                                                contacts::contacts_struct)
# Inputs:
# RNGseed - Seed the random number generator
# n_students - Number of students in the system
# endtime - Number of timesteps simulation will run for
# student_info - entry for each person.
# network_parameters - Items used to generate the contact network
# contacts - Record of who is in contact with whom in each layer

# Outputs:
# dynamic_accom_contacts - Per node, a record of dynamic accomodation (on campus) contacts made on each day


    """
    Unpack required variables
    """
    @unpack contact_prob_floor_level, contact_prob_block_level, contact_prob_hall_level = network_parameters
    @unpack household_member_list, floor_member_list, block_member_list = contacts


    """
    Initialise vector of vectors storing IDs of contacts for each node
    """
    dynamic_accomodation_contacts = Array{Array{Int64,1},3}(undef,endtime,3,n_students)

    """
    Iterate over all nodes
    Check if that student is living on campus
    If so, assign contacts with others from floor, block, hall
    """
    for student_itr = 1:n_students

        if student_info[student_itr].household_info.on_campus_accom==true # Add dynamic links, if living on campus

            # Get information on household, floor, block, hall
            student_household_ID_within_block = student_info[student_itr].household_info.household_ID_within_block
            student_floor_ID = student_info[student_itr].household_info.floor_ID
            student_block_ID = student_info[student_itr].household_info.block_ID
            student_hall_ID = student_info[student_itr].household_info.hall_ID


            # Get number of blocks, floors, households in the hall
            n_blocks_in_hall = length(block_member_list[student_hall_ID])
            n_floors_in_block = length(floor_member_list[student_hall_ID][student_block_ID])
            n_households_on_floor = length(household_member_list[student_hall_ID][student_block_ID][student_floor_ID])

            # Assign the contact probability based on level of accomodation hierarchy
            for block_itr = 1:n_blocks_in_hall
                # Checkings persons in different block. Use hall level contact prob.
                if (block_itr != student_block_ID)
                    contact_prob = contact_prob_hall_level

                    # Check contacts in that hall from a different block.
                    # Push to vector if contact probability accepted
                    persons_to_check = block_member_list[student_hall_ID][block_itr]
                    get_layer_dynamic_accom_contacts!(dynamic_accomodation_contacts,
                                                                RNGseed,
                                                                contact_prob,
                                                                persons_to_check,
                                                                "Hall",
                                                                endtime,
                                                                student_itr
                                                                )
                end
            end

            for floor_itr = 1:n_floors_in_block

                    # Checkings persons in same block, different floor.
                    # Use block level contact prob.
                    if  (floor_itr != student_floor_ID)
                        contact_prob = contact_prob_block_level

                        # Check contacts in that block, but on a different floor.
                        # Push to vector if contact probability accepted
                        persons_to_check = floor_member_list[student_hall_ID][student_block_ID][floor_itr]
                        get_layer_dynamic_accom_contacts!(dynamic_accomodation_contacts,
                                                                    RNGseed,
                                                                    contact_prob,
                                                                    persons_to_check,
                                                                    "Block",
                                                                    endtime,
                                                                    student_itr
                                                                    )
                    end
            end

            for household_itr = 1:n_households_on_floor
                if household_itr != student_household_ID_within_block # Only want to deal with other households on same floor
                    # Checkings persons in same floor, different household.
                    # Use floor level contact prob.
                    contact_prob = contact_prob_floor_level


                    # Check contacts on that floor, but in a different household
                    # Push to vector if contact probability accepted
                    persons_to_check = household_member_list[student_hall_ID][student_block_ID][student_floor_ID][household_itr]
                    get_layer_dynamic_accom_contacts!(dynamic_accomodation_contacts,
                                                                RNGseed,
                                                                contact_prob,
                                                                persons_to_check,
                                                                "Floor",
                                                                endtime,
                                                                student_itr
                                                                )
                end
            end
        end
    end

    return dynamic_accomodation_contacts::Array{Array{Int64,1},3}
end

"""
Functions for use with Erdos-Renyi network construction
"""

function ER_model_uni!(student_info::Array{student_params,1},
                    n_students::Int64,
                    class_sizes::Array{Array{Int64,1},1},
                    dd_within_class::Array{Float64,1},
                    class_contacts::Array{Array{Int64,1},1},
                    cohort_contacts::Array{Array{Int64,1},1},
                    work_or_study_group_contacts_per_node::Array{Int64,1},
                    cohort_contacts_per_node::Array{Int64,1},
                    RNGseed::Int64)
# Inputs:
# student_info - entry for each person.
# n_students - Number of students in the system
# class_sizes - Number of students in each class.
# dd_within_class - Class size degree distribution
# class_contacts, cohort_contacts - Log of contacts each person makes in study classes & with other cohort members.
# work_or_study_group_contacts_per_node - Number of contacts each individual has in work/study setting
# RNGseed - Seed the random number generator


    # Set the RNG
    rng = MersenneTwister(RNGseed)

    # Construct network links for study setting contacts
    for ii = 1:(n_students-1)

        would_attend_f2f_classes::Int64 = student_info[ii].would_attend_f2f_classes
        worktype_idx::Int64 = student_info[ii].cohort_ID
        team_idx::Int64 = student_info[ii].class_ID

        for jj = (ii+1):n_students

            ### On days when in study/work setting, increased contacts with others in that setting ###
            if (would_attend_f2f_classes == 1) & (student_info[jj].would_attend_f2f_classes==1)  # both returned to class
                ## For others in the same work/study setting, edges form according to ER graph
                if (student_info[jj].cohort_ID == worktype_idx) & (student_info[jj].class_ID == team_idx)
                    if rand(rng) < dd_within_class[worktype_idx]/(class_sizes[worktype_idx][team_idx] - 1) ## ER component
                        # Assign IDs of contact to tuple vectors
                        push!(class_contacts[ii],jj)
                        push!(class_contacts[jj],ii)

                        # Increment number of contacts for nodes ii & jj
                        work_or_study_group_contacts_per_node[ii] += 1
                        work_or_study_group_contacts_per_node[jj] += 1
                    end
                elseif (student_info[jj].cohort_ID == worktype_idx) & (student_info[jj].class_ID != team_idx)
                    # Student in same cohort but a different class

                    if rand(rng) < prob_worktype_contact[worktype_idx]
                        # Assign IDs of contact to tuple vectors
                        push!(cohort_contacts[ii],jj)
                        push!(cohort_contacts[jj],ii)

                        # Increment number of contacts for nodes ii & jj
                        cohort_contacts_per_node[ii] += 1
                        cohort_contact_per_node[jj] += 1
                    end
                end
            end
        end
    end
end

# Premake dynamic social contacts,
# to be loaded in ahead of simulation
function generate_dynamic_student_contacts(RNGseed::Int64,
                                            n_students::Int64,
                                            endtime::Int64,
                                            student_info::Array{student_params,1},
                                            dynamic_conts_mean::Array{Float64,1},
                                            dynamic_conts_sd::Array{Float64,1})
# Inputs:
# RNGseed - Seed the random number generator
# n_students - Number of students in the system
# endtime - Number of timesteps simulation will run for
# student_info - entry for each person.
# dynamic_conts_mean, dynamic_conts_sd - Distribution properties for students dynamic social contacts

# Outputs:
# dynamic_social_contacts - Per node, a record of dynamic social contacts made on each day

    """
    Set the RNG
    """
    rng = MersenneTwister(RNGseed)

    """
    Initialise vector of vectors storing IDs of contacts for each node
    """
    dynamic_social_contacts = Array{Array{Int64,1},2}(undef,endtime,n_students)

    """
    Iterate over all nodes
    Assign dynamic social contacts for each timestep
    """
    for student_itr = 1:n_students
            # Get dynamic  group type for student_itr
            node_dynamic_grp_ID = student_info[student_itr].cohort_ID

            for time_itr = 1:endtime
                # Generate number of dynamic contacts from appropriate distribution
                gg = round(Int64,abs(randn(rng)*dynamic_conts_sd[node_dynamic_grp_ID]+dynamic_conts_mean[node_dynamic_grp_ID]))

                # If dynamic worker contacts made, assign to output variable
                if gg > 0
                    dynamic_social_contacts[time_itr,student_itr] = zeros(gg)

                    # Generate required number of contacts
                    for contact_itr = 1:gg
                        gg1 = ceil(Int64,rand(rng)*n_students) # Get IDs of nodes connected by dynamic link on current timestep
                        while gg1 == student_itr # Redraw if returned the index node
                            gg1 = ceil(Int64,rand(rng)*n_students) # Get IDs of nodes connected by dynamic link on current timestep
                        end
                        dynamic_social_contacts[time_itr,student_itr][contact_itr] = gg1
                    end
                else # No dynamic worker contacts made on given timestep
                    dynamic_social_contacts[time_itr,student_itr] = Int64[]
                end
            end
    end

    return dynamic_social_contacts::Array{Array{Int64,1},2}
end

"""
Functions to regenerate specific layers of the network
"""
# Were worker proportion to change from initial value to a non-zero value
# this function allows regeneration of would_attend_f2f_classes for each node &
# reconstruct work networks.
function regenerate_class_contacts(n_students::Int64,network_parameters::network_params,RNGseed::Int64)
# Inputs:
# n_students - Number of workers in the system
# network_parameters
# RNGseed:Int64 - Value to seed the random number generator

# Outputs:
# work_contacts - Vector of vectors with IDs of contacts
# work_contacts_per_node - Total number of regular contacts within work setting

@unpack student_info, class_sizes,prob_workertype_contact,
    dd_within_class, household_size_distribution,
    class_info = network_parameters

    # Set the RNG
    rng = MersenneTwister(RNGseed)

    # Redo the return to work stetting values
    for ii = 1:n_students
        would_attend_f2f_classes = Int64(rand(rng)<attendence_propns[cohort_ids[ii]])   # decide if returning to work
        student_info[ii].would_attend_f2f_classes = would_attend_f2f_classes  # add worker info to node
    end

    # Initialise vector of vectors storing IDs of contacts for each node
    # at study/work setting
    class_contacts = Array{Array{Int64,1},1}(undef,n_students)

    for ii = 1:n_students
        class_contacts[ii] = Int64[]
    end

    # Initialise vectors giving total contacts made by node
    work_or_study_group_contacts_per_node = zeros(Int64,n_students)

    # Construct network links
    for ii = 1:(n_students-1)

        would_attend_f2f_classes::Int64 = student_info[ii].would_attend_f2f_classes
        worktype_idx::Int64 = student_info[ii].cohort_ID
        team_idx::Int64 = student_info[ii].class_ID

        for jj = (ii+1):n_students

            ### On days when at work, increased contacts with others at work  ###
            if (would_attend_f2f_classes == 1) & (student_info[jj].would_attend_f2f_classes==1)  # both returned to work

                ## For workers in the same group and workplace, edges form according to ER graph
                if (student_info[jj].cohort_ID == worktype_idx) & (student_info[jj].class_ID == team_idx)
                    if rand(rng) < dd_within_class[worktype_idx]/(class_sizes[worktype_idx][team_idx] - 1) ## ER component
                        # Assign IDs of contact to tuple vectors
                        push!(class_contacts[ii],jj)
                        push!(class_contacts[jj],ii)

                        # Increment number of contacts for nodes ii & jj
                        work_or_study_group_contacts_per_node[ii] += 1
                        work_or_study_group_contacts_per_node[jj] += 1
                    end
                end
            end
        end
    end

    # Define what is returned from the function
    return  class_contacts::Array{Array{Int64,1},1},
                work_or_study_group_contacts_per_node::Array{Int64,1}
end

function regenerate_society_contacts(network_parameters::network_params,
                                            RNGseed::Int64,
                                            n_students::Int64,
                                            starttime::Int64,
                                            endtime::Int64)
# Inputs:
# network_parameters - Fields relating to the network generation incl. student_info & prob_social_contact
# RNGseed - Seed the random number generator
# n_students - Number of nodes in the system
# starttime - Timestep to begin regenerating contacts from
# endtime - Number of timesteps simulation will run for

# Outputs:
# social_contacts_per_node - Total amount of social contacts each individual has (entry per individual)
# society_contacts
#       - A now amended record of society contacts

    @unpack student_info, prob_social_contact, society_info = network_parameters

    """
    Set the RNG
    """
    rng = MersenneTwister(RNGseed)

    """
    Draw new collection of possible social contacts
    """
    # Initialise vector of vectors storing IDs of contacts for each node in social
    # and household settings
    n_societies = length(society_info)
    society_contacts = Array{Array{Int64,1},2}(undef,n_students,n_societies)

    # Initialise vectors giving total contacts made by node
    society_contacts_per_node = zeros(Int64,n_students,n_societies)

    # Reconstruct network links for society contacts
    for ii = 1:(n_students-1)
        person_ii_society_idxs::Int64 = student_info[ii].society_IDs

        for jj = (ii+1):n_students
            person_jj_society_idxs::Int64 = student_info[jj].society_IDs
            if (length(person_ii_society_idxs) > 0) &&
                (length(person_jj_society_idxs) > 0) # No checks needed if either person not a member of any society

                # Iterate over each society person ii is a member of
                for society_itr = 1:length(person_ii_society_idxs)
                    person_ii_current_society_ID = person_ii_society_idxs[society_itr]

                    # Check if person jj is also a member of that society
                    # If so, apply uniform possibility of contacts with all others in society.
                    # More people gives more contacts.
                    for person_jj_society_itr = 1:length(person_jj_society_idxs)
                        person_jj_current_society_ID = person_jj_society_idxs[person_jj_society_itr]
                        if (person_jj_current_society_ID == person_ii_current_society_ID)
                            # Get relevant society contact probability
                            society_type_val = society_info[person_ii_current_society_ID].society_type
                            if rand(rng) < prob_social_contact[society_type_val]
                                # Assign IDs of contact to tuple vectors
                                push!(society_contacts[ii,person_ii_current_society_ID],jj)
                                push!(society_contacts[jj,person_jj_current_society_ID],ii)

                                # Increment number of contacts for nodes ii & jj
                                society_contacts_per_node[ii,person_ii_current_society_ID] += 1
                                society_contacts_per_node[jj,person_jj_current_society_ID] += 1
                            end
                        end
                    end
                end
            end
        end
    end

    return society_contacts_per_node::Array{Int64,2},
            society_contacts::Array{Array{Int64,1},2}
end


# Redraw dynamic social contacts from specified start time to end time
# Allocate to pre-existing array
function regenerate_dynamic_student_contacts(dynamic_social_contacts::Array{Array{Int64,1},2},
                                            RNGseed::Int64,
                                            n_students::Int64,
                                            starttime::Int64,
                                            endtime::Int64,
                                            student_info::Array{student_params,1},
                                            dynamic_social_contact_degree_distribution::Array{Distribution,1})

# Inputs:
# dynamic_social_contacts - Per node, a record of dynamic worker contacts made on each day
# RNGseed - Seed the random number generator
# n_students - Number of nodes in the system
# starttime - Timestep to begin regenerating contacts from
# endtime - Number of timesteps simulation will run for
# student_info - array with entry per worker
# dynamic_social_contact_degree_distribution- Distribution properties for students dynamic social contacts

# Outputs:
# dynamic_social_contacts - The now amended record of dynamic worker contacts made on each day

    """
    Set the RNG
    """
    rng = MersenneTwister(RNGseed)

    """
    Iterate over all nodes
    For those returning to work in role with dynamic contacts,
    assign dynamic contacts for each timestep
    """
    for student_itr = 1:n_students
        if student_info[student_itr].would_attend_f2f_classes==1 # Add dynamic links, if returned to work
            # Get dynamic worker group type for student_itr
            node_dynamic_grp_ID = student_info[student_itr].cohort_ID

            for time_itr = starttime:endtime
                # Generate number of dynamic contacts from appropriate distribution
                gg = round(Int64, Distributions.rand(rng, dynamic_social_contact_degree_distribution[node_dynamic_grp_ID]))

                # If dynamic worker contacts made, assign to output variable
                if gg > 0
                    dynamic_social_contacts[time_itr,student_itr] = zeros(gg)

                    # Generate required number of contacts
                    for contact_itr = 1:gg
                        gg1 = ceil(Int64,rand(rng)*n_students) # Get IDs of nodes connected by dynamic link on current timestep
                        while gg1 == student_itr # Redraw if returned the index node
                            gg1 = ceil(Int64,rand(rng)*n_students) # Get IDs of nodes connected by dynamic link on current timestep
                        end
                        dynamic_social_contacts[time_itr,student_itr][contact_itr] = gg1
                    end
                else # No dynamic worker contacts made on given timestep
                    dynamic_social_contacts[time_itr,student_itr] = Int64[]
                end
            end
        end
    end

    return dynamic_social_contacts::Array{Array{Int64,1},2}
end
