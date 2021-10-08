#=
Purpose:
Stash functions that are used with the university network model
for performing contact tracing
=#
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Get portion of dynamic contacts to be recallable
#-------------------------------------------------------------------------------

# In time d days before symptoms up to current day.
function recallable_dynamic_contacts(student_ID::Int64,
                                            all_dynamic_contacts::Array{Int64,1},
                                            dynamic_contacts_recalled_propn::Array{Float64,1},
                                            daily_record_inisol::Array{Int64,2},
                                            time_to_check::Int64,
                                            prev_day_val::Int64,
                                            rng::MersenneTwister
                                            )
# Inputs:
# student_ID
# all_dynamic_contacts - IDs of those dynamic contacts by the student for the current timestep
# dynamic_contacts_recalled_propn::Array{Float64,1} - Set up recall of dynamic contacts probability (Proportion of contacts remembered x days ago)
# daily_record_inisol - For each node, records per timestep whether in isolation
# time_to_check - Timestep of the simulation to be checked.
# prev_day_val - Number of days previous to be looked at
# rng - The random number generator in use

# Outputs:
# recallable_contacts - IDs of individuals (contacted by dynamic contact) to be retained

    # Get proportion of dynamic contacts from prev_day_val days ago to be retained
    if prev_day_val > length(dynamic_contacts_recalled_propn)
        # No contacts can possible be retained.
        # Return an empty vector
        recallable_contacts = Int64[]
    else
        # Get proportion of contacts expected to be retained from input distribution
        threshold = dynamic_contacts_recalled_propn[prev_day_val]

        # Check if any dynamic contacts were isolating.
        # If so, will not be a recallable contact
        n_possible_dynamic_contacts = length(all_dynamic_contacts)
        recallable_contact_check = zeros(Int64,n_possible_dynamic_contacts)
        for contact_itr = 1:n_possible_dynamic_contacts
            # If not isolating, then contact did occur
            contact_ID = all_dynamic_contacts[contact_itr]
            if daily_record_inisol[time_to_check,contact_ID]==false
                # Contact occurred
                # Check if contact will be remembered
                r = rand(rng)
                if r<threshold
                    recallable_contact_check[contact_itr] = 1
                end
            end
        end

        # Construct vector of IDs of contacts that did actually occur
        n_recallable_contacts = sum(recallable_contact_check)
        recallable_contacts = zeros(Int64,n_recallable_contacts) # Initialise vector to store IDs of recallable contacts
        recallable_contact_itr = 1 # Initialise counter for vector assignment index
        for contact_itr = 1:n_possible_dynamic_contacts
            if recallable_contact_check[contact_itr] == true
                contact_ID = all_dynamic_contacts[contact_itr]
                recallable_contacts[recallable_contact_itr] = contact_ID
                recallable_contact_itr += 1 # Increment vector assignment index
            end
        end
    end

    return recallable_contacts::Array{Int64,1}
end

#-------------------------------------------------------------------------------
# Check contacts made in class setting
#-------------------------------------------------------------------------------

# In time d days before symptoms up to current day.
# Go over list of usual class setting contacts & cohort contacts. Check if they actually occurred,
# or if contact was absent due to isolation/other closure.
function get_study_contacts(possible_study_contacts::Array{Int64,1},
                                daily_record_inisol::Array{Int64,2},
                                daily_record_inclass::Array{Int64,2},
                                time_to_check::Int64,
                                prev_day_val::Int64,
                                rng::MersenneTwister
                                )
# Inputs:
# possible_study_contacts - IDs of study contacts. Will check if actually occurred
# daily_record_inisol - For each node, records per timestep whether in isolation
# daily_record_inclass - For each node, records per timestep whether at workplace
# time_to_check - Timestep of the simulation to be checked.
# prev_day_val - Number of days previous to be looked at
# rng - The random number generator in use

# Outputs:
# class_contacts - IDs of individuals to be retained for CT

    # Check if any contacts were isolating and/or not in the class or work setting that day
    # If so, will not be a contact
    n_possible_study_contacts = length(possible_study_contacts)
    contact_occur_check = zeros(Int64,n_possible_study_contacts)
    for contact_itr = 1:n_possible_study_contacts

        # If contact not isolating and in the class or work setting, then contact can occur.
        contact_ID = possible_study_contacts[contact_itr]
        if (daily_record_inisol[time_to_check,contact_ID]==false) &&
            (daily_record_inclass[time_to_check,contact_ID]==true)

            # Contact in same workplace
            contact_occur_check[contact_itr] = 1
        end
    end

    # Construct vector of IDs of contacts that did actually occur
    n_occur_study_contacts = sum(contact_occur_check)
    study_contacts = zeros(Int64,n_occur_study_contacts) # Initialise vector to store IDs of recallable contacts
    contact_occur_itr = 1 # Initialise counter for vector assignment index
    for contact_itr = 1:n_possible_study_contacts
        if contact_occur_check[contact_itr] == 1
            contact_ID = possible_study_contacts[contact_itr]
            study_contacts[contact_occur_itr] = contact_ID
            contact_occur_itr += 1 # Increment vector assignment index
        end
    end

    return study_contacts::Array{Int64,1}
end

#-------------------------------------------------------------------------------
# Check contacts made in society setting
#-------------------------------------------------------------------------------

# In time d days before symptoms up to current day.
# Go over list of usual society contacts. Check if they actually occurred,
# or if contact was absent due to isolation/other closure.
function get_society_contacts(possible_society_contacts::Array{Int64,1},
                                society_ID::Int64,
                                daily_record_inisol::Array{Int64,2},
                                daily_record_atsociety::Array{Int64,3},
                                time_to_check::Int64,
                                prev_day_val::Int64,
                                rng::MersenneTwister
                                )
# Inputs:
# possible_society_contacts - IDs of usual society contacts. Will check if actually occurred
# society_ID -
# daily_record_inisol - For each node, records per timestep whether in isolation
# daily_record_atsociety - For each node, records per timestep whether attended society
# time_to_check - Timestep of the simulation to be checked.
# prev_day_val - Number of days previous to be looked at
# rng - The random number generator in use

# Outputs:
# recallable_society_contacts - IDs of individuals contacted in a society to be retained for CT

    # Check if any contacts were isolating and/or not in the class or work setting that day
    # If so, will not be a contact
    n_possible_society_contacts = length(possible_society_contacts)
    contact_occur_check = zeros(Int64,n_possible_society_contacts)
    for contact_itr = 1:n_possible_society_contacts

        # If contact not isolating and in the class or work setting, then contact can occur.
        contact_ID = possible_society_contacts[contact_itr]
        if (daily_record_inisol[time_to_check,contact_ID]==false) &&
            (daily_record_atsociety[time_to_check,society_ID,contact_ID]==true)

            # Contact in same workplace
            contact_occur_check[contact_itr] = 1
        end
    end

    # Construct vector of IDs of contacts that did actually occur
    n_occur_society_contacts = sum(contact_occur_check)
    recallable_society_contacts = zeros(Int64,n_occur_society_contacts) # Initialise vector to store IDs of recallable contacts
    contact_occur_itr = 1 # Initialise counter for vector assignment index
    for contact_itr = 1:n_possible_society_contacts
        if contact_occur_check[contact_itr] == 1
            contact_ID = possible_society_contacts[contact_itr]
            recallable_society_contacts[contact_occur_itr] = contact_ID
            contact_occur_itr += 1 # Increment vector assignment index
        end
    end

    return recallable_society_contacts::Array{Int64,1}
end

#-------------------------------------------------------------------------------
# Perform forward CT from an identified infector
#-------------------------------------------------------------------------------

# Comment break
function forwardCT_from_infector!(infector_ID::Int64,
                                CT_vars::contact_tracing_vars,
                                contacts::contacts_struct,
                                CT_parameters::CT_params,
                                engage_with_CT::Array{Bool,1},
                                inclass::Array{Int64,2},
                                time::Int64,
                                count::Int64,
                                num_CT::Array{Int64,2},
                                hh_isolation::Array{Int64,1},
                                timeisol_CTcause::Array{Int64,1},
                                network_parameters::network_params,
                                rng::MersenneTwister,
                                infector_trace_count::Array{Int64,1})

@unpack dynamic_contacts_recalled_propn, social_contacts_recalled_propn, infector_engage_with_CT_prob = CT_parameters
@unpack Inds_to_be_contacted, Test_result_false_negative = CT_vars

    # Already contact traced? If not, do it now
    if (isassigned(Inds_to_be_contacted,infector_ID)) # Checks if infector_ID reported infection.
                                                      # If not, no forward contact tracing will be done from infector
        if !(length(Inds_to_be_contacted[infector_ID])>0) # Inds_to_be_contacted[infector_ID] being empty signifies
                                                         # infector reported symptoms, but returned false neg and/or
                                                         # did not engage in contact tracing
            # check if the infector will engage with CT
            if (((engage_with_CT[infector_ID] == false) && (rand(rng)<infector_engage_with_CT_prob))
                || (engage_with_CT[infector_ID] == true)) &&
                (Test_result_false_negative[infector_ID] == false)
                        # didn't engage before but does now or already willing to engage
                        # and then checks infector has not previously tested negative

                infector_trace_count[1] += 1

                trace_node!(infector_ID,time,CT_vars,contacts,CT_parameters,network_parameters,rng)
            end
        end
    end
end

function trace_node!(student_itr::Int64,time::Int64,CT_vars::contact_tracing_vars,
    contacts::contacts_struct,
    CT_parameters::CT_params,network_parameters::network_params,rng)
    @unpack student_info, class_info = network_parameters

    # Preload study group/work setting contacts
    possible_class_contacts = contacts.class_contacts[student_itr]
    possible_cohort_contacts = contacts.cohort_contacts[student_itr]

    # Preload student info
    cohort_ID = student_info[student_itr].cohort_ID
    class_ID = student_info[student_itr].class_ID
    society_IDs = student_info[student_itr].society_IDs

    # Get contacts that will be contacted
    for time_itr = 1:CT_vars.relevant_prev_days_for_CT[student_itr]
        if time-time_itr > 0 # Can't look back before the simulation started

            # Get previous time being checked
            time_to_check = time-time_itr

            # If in class/in work setting, get contacts made that day
            if (contacts.daily_record_inclass[time_to_check,student_itr]==true)

                # Contacts made in teaching and work settings
                class_contacts = get_study_contacts(possible_class_contacts,
                                                            contacts.daily_record_inisol,
                                                            contacts.daily_record_inclass,
                                                            time_to_check,
                                                            time_itr,
                                                            rng)

                # If there are class/work setting contacts,
                # add to vector tracking traceable contacts
                if !isempty(class_contacts)
                    append!(CT_vars.Inds_to_be_contacted[student_itr],class_contacts)
                end

                # Contacts made with other cohort members, but in different classes
                # (on study days)
                cohort_contacts = get_study_contacts(possible_cohort_contacts,
                                                            contacts.daily_record_inisol,
                                                            contacts.daily_record_inclass,
                                                            time_to_check,
                                                            time_itr,
                                                            rng)

                # If there are class/work setting contacts,
                # add to vector tracking traceable contacts
                if !isempty(cohort_contacts)
                    append!(CT_vars.Inds_to_be_contacted[student_itr],cohort_contacts)
                end
            end

            # Iterate over each society the student is a member of
            # Get society contacts
            for society_itr = 1:length(society_IDs)
                current_society_ID = society_IDs[society_itr]
                if (contacts.daily_record_atsociety[time_to_check,current_society_ID,student_itr]==true)

                    # Contacts made in teaching and work settings
                    recallable_society_contacts = get_society_contacts(contacts.society_contacts[student_itr,current_society_ID], # possible society contacts
                                                                current_society_ID,
                                                                contacts.daily_record_inisol,
                                                                contacts.daily_record_atsociety,
                                                                time_to_check,
                                                                time_itr,
                                                                rng)

                    # If there are class/work setting contacts,
                    # add to vector tracking traceable contacts
                    if !isempty(recallable_society_contacts)
                        append!(CT_vars.Inds_to_be_contacted[student_itr],recallable_society_contacts)
                    end
                end
            end

            # Get recallable dynamic social contacts
            all_dynamic_contacts = contacts.dynamic_social_contacts[time_to_check,student_itr]  # Get IDs of those dynamic contacts on required day
            dynamic_recallable_contacts = recallable_dynamic_contacts(student_itr,
                                                                                all_dynamic_contacts,
                                                                                CT_parameters.dynamic_contacts_recalled_propn,
                                                                                contacts.daily_record_inisol,
                                                                                time_to_check,
                                                                                time_itr,
                                                                                rng) # in contact_tracing_fns.jl

            # If there are recallable dynamic contacts, add to vector tracking traceable contacts
            if !isempty(dynamic_recallable_contacts)
                append!(CT_vars.Inds_to_be_contacted[student_itr],dynamic_recallable_contacts)
            end

            # Get recallable dynamic accommodation contacts
            for accom_level = 1:3
                if isassigned(contacts.dynamic_accommodation_contacts,time,accom_level,student_itr)
                    accom_dynamic_contacts = contacts.dynamic_accommodation_contacts[time_to_check,accom_level,student_itr]  # Get IDs of those dynamic contacts on required day
                    dynamic_accom_recallable_contacts = recallable_dynamic_contacts(student_itr,
                                                                                    accom_dynamic_contacts,
                                                                                    CT_parameters.accom_dynamic_contacts_recalled_propn,
                                                                                    contacts.daily_record_inisol,
                                                                                    time_to_check,
                                                                                    time_itr,
                                                                                    rng) # in contact_tracing_fns.jl

                    # If there are recallable dynamic contacts, add to vector tracking traceable contacts
                    if !isempty(dynamic_accom_recallable_contacts)
                        append!(CT_vars.Inds_to_be_contacted[student_itr],dynamic_accom_recallable_contacts)
                    end
                end
            end
        end
    end
end
