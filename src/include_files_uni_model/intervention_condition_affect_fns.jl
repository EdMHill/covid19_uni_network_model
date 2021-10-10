#=
Purpose:
Define intervention functions

Conditions functions (used to check if an affect should be enacted):
 - condition_close_example
 - condition_open_example

The rest are affect functions.
May include conditions within it or make calls to specific condition fns

Cohort level functions
 - affect_close_example!, affect_open_example!

Society related functions
 - society_check

accommodation functions
- accom_lockdown_unit_check!  (Find those accommodation units whose lockdown measures need updating)
- accom_lockdown_unit_alter_status! (Update those students in accommodation units lockdown measures are being enforced)
- accom_lockdown_hall_level!, accom_lockdown_block_level!, accom_lockdown_floor_level!,
   accom_lockdown_household_level!  (lockdown accommodation, with conditions at the specified level)
=#
#-------------------------------------------------------------------------------

"""
    condition_close_example(intervention_trigger_input_data::intervention_data_feeds,
                             time::Int64,
                             network_parameters::network_params)

Example of a triggered intervention based on outbreak situation. if incidence of symptomatic infection exceeds given level, apply restriction.

Inputs:
- `intervention_trigger_input_data`: intervention_data_feeds structure
- `time`: Elapsed time
- `network_parameters`: network_params structure

Outputs: `output_bool`: Denotes whether condition outcome was true or false. \n
Location: intervention\\_condition\\_affect\\_fns.jl
"""
function condition_close_example(intervention_trigger_input_data::intervention_data_feeds,
                                        time::Int64,
                                        network_parameters::network_params)
   @unpack n_students = network_parameters
   @unpack rep_inf_this_timestep = intervention_trigger_input_data

   # Check if daily incidence condition surpassed
   daily_incidence = sum(rep_inf_this_timestep)/n_students
   if (daily_incidence > 0.01)
      output_bool = true
   else
      output_bool = false
   end

   return output_bool::Bool
end

"""
    condition_open_example(intervention_trigger_input_data::intervention_data_feeds,
                             time::Int64,
                             network_parameters::network_params)

Example of a triggered intervention based on outbreak situation. if incidence of symptomatic infection drops below given level, revert restriction.

Inputs:
- `intervention_trigger_input_data`: intervention_data_feeds structure
- `time`: Elapsed time
- `network_parameters`: network_params structure

Outputs: `output_bool`: Denotes whether condition outcome was true or false. \n
Location: intervention\\_condition\\_affect\\_fns.jl
"""
function condition_open_example(intervention_trigger_input_data::intervention_data_feeds,
                                        time::Int64,
                                        network_parameters::network_params)
 # Output:
#   output_bool - Denotes whether condition outcome was true or false.

   @unpack n_students  = network_parameters
   @unpack rep_inf_this_timestep = intervention_trigger_input_data

   # Check if daily incidence condition declined to low enough level
   daily_incidence = sum(rep_inf_this_timestep)/n_students
   if daily_incidence < 0.001
      output_bool = true
   else
      output_bool = false
   end

   return output_bool::Bool
end



#=
Cohort level affect functions.
=#
"""
    dummy_example!(network_parameters::network_params,
                            intervention_trigger_input_data::intervention_data_feeds,
                            contacts::contacts_struct,
                            time::Int64)

Dummy example for teaching cohort intervention.

Inputs:
- `network_parameters`: network_params structure
- `intervention_trigger_input_data`: intervention_data_feeds structure
- `contacts`: contacts_struct structure
- `time`: Elapsed time

Outputs: None \n
Location: intervention\\_condition\\_affect\\_fns.jl
"""
function dummy_example!(network_parameters::network_params,
                           intervention_trigger_input_data::intervention_data_feeds,
                           contacts::contacts_struct,
                           time::Int64)
   return nothing
end

"""
    affect_close_example!(network_parameters::network_params,
                            intervention_trigger_input_data::intervention_data_feeds,
                            contacts::contacts_struct,
                            time::Int64)

Stop face-to-face teaching for cohorts.

Inputs:
- `network_parameters`: network_params structure
- `intervention_trigger_input_data`: intervention_data_feeds structure
- `contacts`: contacts_struct structure
- `time`: Elapsed time

Outputs: None \n
Location: intervention\\_condition\\_affect\\_fns.jl
"""
function affect_close_example!(network_parameters::network_params,
                                 intervention_trigger_input_data::intervention_data_feeds,
                                 contacts::contacts_struct,
                                 time::Int64)

   @unpack class_info, cohort_f2f_study_active = network_parameters

   # Check if daily incidence condition surpassed
   output_bool = condition_close_example(intervention_trigger_input_data,
                                           time,
                                           network_parameters)
   # If condition satisfied, apply affect
   if output_bool == true

      # Get number of cohorts in use
      n_cohort_types = length(class_info)

      # # The cohorts to be closed
      # close_cohort_IDs = [1,2]
      # n_cohorts_closed = length(close_cohort_IDs)

      # Iterate over each group type to become inactive
      for close_cohort_itr = 1:n_cohort_types

         # Get group type of interest
         close_cohort_ID = close_cohort_itr

         if (cohort_f2f_study_active[close_cohort_ID] == true)
            # Get number of teams in that group type
            n_classes = length(class_info[close_cohort_ID])

            # Set status of each team in that group type to be inactive
            for class_itr = 1:n_classes
               class_info[close_cohort_ID][class_itr].f2f_class = false
            end

            # Update cohort_f2f_study_active field
            cohort_f2f_study_active[close_cohort_ID] = false
         end
      end
   end
end

"""
    affect_open_example!(network_parameters::network_params,
                            intervention_trigger_input_data::intervention_data_feeds,
                            contacts::contacts_struct,
                            time::Int64)

Restart face-to-face teaching for cohorts.

Inputs:
- `network_parameters`: network_params structure
- `intervention_trigger_input_data`: intervention_data_feeds structure
- `contacts`: contacts_struct structure
- `time`: Elapsed time

Outputs: None \n
Location: intervention\\_condition\\_affect\\_fns.jl
"""
function affect_open_example!(network_parameters::network_params,
                                 intervention_trigger_input_data::intervention_data_feeds,
                                 contacts::contacts_struct,
                                 time::Int64)

   @unpack class_info, cohort_f2f_study_active = network_parameters

   # Check if daily incidence condition surpassed
   output_bool = condition_open_example(intervention_trigger_input_data,
                                           time,
                                           network_parameters)
   # If condition satisfied, apply affect
   if output_bool == true

      # Get number of cohorts in use
      n_cohort_types = length(class_info)

      # # The group types to become active
      # open_cohort_IDs = [1,2]
      # n_cohorts_open = length(open_cohort_IDs)

      # Iterate over each cohort to be permitted to have f2f classes
      for open_cohort_itr = 1:n_cohort_types

         # Get group type of interest
         open_cohort_ID = open_cohort_itr

         if (cohort_f2f_study_active[open_cohort_ID] == false)
            # Get number of teams in that group type
            n_classes = length(class_info[open_cohort_ID])

            # Set status of each team in that group type to be active
            for class_itr = 1:n_classes
               class_info[open_cohort_ID][class_itr].f2f_class = true
            end

            # Update cohort_f2f_study_active field
            cohort_f2f_study_active[open_cohort_ID] = true
         end
      end
   end
end

#=
Society related functions
=#
"""
    society_incidence_check!(network_parameters::network_params,
                            intervention_trigger_input_data::intervention_data_feeds,
                            contacts::contacts_struct,
                            time::Int64)

Check if society can remain active based on recent reported cases.

Inputs:
- `network_parameters`: network_params structure
- `intervention_trigger_input_data`: intervention_data_feeds structure
- `contacts`: contacts_struct structure
- `time`: Elapsed time

Outputs: None \n
Location: intervention\\_condition\\_affect\\_fns.jl
"""
function society_incidence_check!(network_parameters::network_params,
                                 intervention_trigger_input_data::intervention_data_feeds,
                                 contacts::contacts_struct,
                                 time::Int64)

   @unpack society_info, student_info, society_f2f_active = network_parameters

    # Initialise intervention_params
    intervention_parameters = intervention_params(time_horizon = 14,
                                                   inactivation_length = 14,
                                                   absolute_case_threshold = 3,
                                                   propn_case_threshold = 0.1)

    # Want to iterate over each society
    n_societies = length(society_info)

    # Iterate over each society
    for society_itr = 1:n_societies

       # Use the number of members per society
       n_society_members = society_info[society_itr].n_members

       # Get students in that society.
       society_member_IDs = society_info[society_itr].member_list

       if society_info[society_itr].society_inactivation_time > 0
          # If previously made inactive,
          # and if activation time check not yet reached, increment the counter
          society_info[society_itr].society_inactivation_time += 1

          # if the group has been inactive for long enough, can check if can be restarted
          if society_info[society_itr].society_inactivation_time > intervention_parameters.inactivation_length

             # Check infection count is low enough to restart activity
             recent_reported_cases = 0 # Initialise counter of recent reported cases
             for student_member_itr = 1:n_society_members
                # For student, get ID and check if infected within the time horizon
                student_ID = society_member_IDs[student_member_itr]
                cutoff_date = min(1,time - intervention_parameters.time_horizon)
                if (student_info[student_ID].time_of_reporting_infection >= cutoff_date)
                  recent_reported_cases += 1  # Student reported infection within time horizon. Increment counter
                end
             end
             propn_members_reporting = recent_reported_cases/n_society_members

             if (recent_reported_cases < intervention_parameters.absolute_case_threshold) &&
                  (propn_members_reporting < intervention_parameters.propn_case_threshold)
                # If deemed okay to restart, reset the inactivation timer
                society_info[society_itr].society_inactivation_time = 0

                # Check entire society type is permitted to have activity
                # If it is permitted, the individual class can be set to be active.
                society_type = society_info[society_itr].society_type
                if (society_f2f_active[society_type] == true)
                     society_info[society_itr].f2f_activity = true
                end
             end
          end

       else
          # Currently active. Check if permitted to stay active

          # Get count & propn. of members that have reported infection in last X days
          # Uses time_of_reporting_infection field in student_info
          recent_reported_cases = 0 # Initialise counter of recent reported cases
          for student_member_itr = 1:n_society_members
             # For student, get ID and check if infected within the time horizon
             student_ID = society_member_IDs[student_member_itr]
             cutoff_date = min(1,time - intervention_parameters.time_horizon)
             if (student_info[student_ID].time_of_reporting_infection >= cutoff_date)
                recent_reported_cases += 1  # Student reported infection within time horizon. Increment counter
             end
          end
          propn_members_reporting = recent_reported_cases/n_society_members

          # If rate is above a certain level & a certain number of cases, stop f2f activity
          if (recent_reported_cases >= intervention_parameters.absolute_case_threshold) &&
                (propn_members_reporting > intervention_parameters.propn_case_threshold)
             society_info[society_itr].f2f_activity = false
             society_info[society_itr].society_inactivation_time = 1 # Begin inactivation counter
          end
       end
    end
end


#=
accommodation related functions
=#

"""
    accom_lockdown_unit_check(unit_resident_list::Array{Int64,1},
                               intervention_parameters::intervention_params,
                               student_info::Array{student_params,1},
                               time::Int64)

Check if accomodation unit enters lockdown.

Inputs:
- `unit_resident_list::Array{Int64,1}`: The IDs of students in the accommodation unit being checked
- `intervention_trigger_input_data::intervention_params`: Values used to check if measures should be enforced/relaxed
- `student_info::student_params `: Information of each student in the system
- `time`: Current time of the simulation

Outputs:
- `recent_reported_cases::Int64`: Number of cases within relevant time period in accommodation unit of interest
- `propn_unit_reporting::Float64`: Porportion of popn impacted within relevant time period in accommodation unit of interest

Location: intervention\\_condition\\_affect\\_fns.jl
"""
function accom_lockdown_unit_check(unit_resident_list::Array{Int64,1},
                                       intervention_parameters::intervention_params,
                                       student_info::Array{student_params,1},
                                       time::Int64
                                       )

   @unpack time_horizon, inactivation_length, absolute_case_threshold,
            propn_case_threshold, release_condition_time_horizon,
            release_condition_absolute_cases, release_condition_propn_cases  = intervention_parameters

   # Get number of students in the accommodation unit of interest
   n_students_in_unit = length(unit_resident_list)

   # Check if accommodation block is in lockdown already
   # Use first student, as status should be the same for all students in that accommodation
   first_student_ID = unit_resident_list[1]
   current_accom_lockdown_status = student_info[first_student_ID].household_info.lockdown_status

   # Get count & propn. of household that have reported infection during relevant time horizon
   # Uses time_of_reporting_infection field in student_info
   recent_reported_cases = 0 # Initialise counter of recent reported cases
   for student_itr = 1:n_students_in_unit
      # For student, get ID and check if infected within the time horizon
      student_ID = unit_resident_list[student_itr]
      cutoff_date = max(1,time - time_horizon)  # Ensure cutoff_date is not 0 or below
      if (student_info[student_ID].time_of_reporting_infection >= cutoff_date)
         recent_reported_cases += 1  # Student reported infection within time horizon. Increment counter
      end
   end
   propn_unit_reporting = recent_reported_cases/n_students_in_unit

   condition_values = [recent_reported_cases;propn_unit_reporting]

   return condition_values::Array{Float64,1}
end

"""
    accom_lockdown_unit_alter_status!(unit_resident_list::Array{Int64,1},
                                       intervention_parameters::intervention_params,
                                       student_info::Array{student_params,1},
                                       time::Int64,
                                       condition_values::Array{Float64,1})

Modify accomodation unit lockdown status.

Inputs:
- `unit_resident_list::Array{Int64,1}`: The IDs of students in the accommodation unit being checked
- `intervention_parameters::intervention_params`: Valeus used to check if measures should be enforced/relaxed
- `student_info::student_params`: Information of each student in the system
- `time::Int64`: Current time of the simulation
- `condition_values::Array{Float64,1}``: Criteria to check to determine of accomodation lockdown status requires amendment

Outputs: None \n
Location: intervention\\_condition\\_affect\\_fns.jl
"""
function accom_lockdown_unit_alter_status!(unit_resident_list::Array{Int64,1},
                                       intervention_parameters::intervention_params,
                                       student_info::Array{student_params,1},
                                       time::Int64,
                                       condition_values::Array{Float64,1}
                                       )

   @unpack time_horizon, inactivation_length, absolute_case_threshold,
            propn_case_threshold, release_condition_time_horizon,
            release_condition_absolute_cases, release_condition_propn_cases  = intervention_parameters

   # Initialise status for amending lockdown status
   alter_accom_lockdown_status = "no"

   # Disaggregate condition_values
   #  condition_values = [recent_reported_cases;propn_unit_reporting]
   recent_reported_cases::Int64 = convert(Int64,condition_values[1])
   propn_unit_reporting::Float64 = condition_values[2]

   # Get number of students in the accommodation unit of interest
   n_students_in_unit = length(unit_resident_list)

   # Check if accommodation unit is in lockdown already
   # Use first student, as status should be the same for all students in that accommodation
   first_student_ID = unit_resident_list[1]
   current_accom_lockdown_status = student_info[first_student_ID].household_info.lockdown_status

   if current_accom_lockdown_status == true
      # accommodation in lockdown. See if it can be released.

      # If rate is below a certain level & a certain number of cases
      # release lockdown
      if (recent_reported_cases < release_condition_absolute_cases) &&
            (propn_unit_reporting < release_condition_propn_cases)
         alter_accom_lockdown_status = "release"
      end

   else
      # accommodation not in lockdown. Check if it should be

      # If rate is above a certain level & a certain number of cases
      # enforce lockdown
      if (recent_reported_cases >= absolute_case_threshold) &&
            (propn_unit_reporting >= propn_case_threshold)
         alter_accom_lockdown_status = "enforce"
      end
   end

   # If lockdown activated/released, alter status for each student in that unit
   if alter_accom_lockdown_status != "no"
      for student_in_unit_itr = 1:n_students_in_unit
         # Get the ID of the student
         current_student_ID = unit_resident_list[student_in_unit_itr]

         # Update the status of accommodation lockdown for the student
         if alter_accom_lockdown_status == "enforce"
            student_info[current_student_ID].household_info.lockdown_status = true
         elseif alter_accom_lockdown_status == "release"
            student_info[current_student_ID].household_info.lockdown_status = false
         end
      end
   end
end

"""
    accom_lockdown_hall_level!(network_parameters::network_params,
                                    intervention_trigger_input_data::intervention_data_feeds,
                                    contacts::contacts_struct,
                                    time::Int64)

Check whether to "lock down" accommodation at hall level.

Applies to oncampus accommodation only. \n
Class/society activities no longer permitted for affected students. \n

Inputs:
- `network_parameters::network_params`: network_params structure
- `intervention_trigger_input_data::intervention_params`: Values used to check if measures should be enforced/relaxed
- `contacts::contacts_struct`: contacts_struct structure
- `time`: Current time of the simulation

Outputs: None \n
Location: intervention\\_condition\\_affect\\_fns.jl
"""
function accom_lockdown_hall_level!(network_parameters::network_params,
                                    intervention_trigger_input_data::intervention_data_feeds,
                                    contacts::contacts_struct,
                                    time::Int64)

   @unpack student_info = network_parameters
   @unpack hall_member_list = contacts

   # Initialise intervention_params
   intervention_parameters = intervention_params(time_horizon = 14,
                                                  absolute_case_threshold = 3,
                                                  propn_case_threshold = 0.1,
                                                  release_condition_time_horizon = 7,
                                                  release_condition_absolute_cases = 0,
                                                  release_condition_propn_cases = 0)

   # Want to iterate over each household
   # Get number of halls and iterate over each one
   n_halls = length(hall_member_list)
   for hall_itr = 1:n_halls

      # Is at this level we want to enact control
      # Get all students in unit
      unit_resident_list = hall_member_list[hall_itr]
      n_students_in_unit = length(unit_resident_list)

      # Call function to do condition check
      condition_values::Array{Float64,1} = accom_lockdown_unit_check(unit_resident_list,
                                                                     intervention_parameters,
                                                                     student_info,
                                                                     time
                                                                     )

      # Call function to modify statuses as appropriate
      accom_lockdown_unit_alter_status!(unit_resident_list,
                                       intervention_parameters,
                                       student_info,
                                       time,
                                       condition_values
                                       )
   end
end

"""
    accom_lockdown_block_level!(network_parameters::network_params,
                                    intervention_trigger_input_data::intervention_data_feeds,
                                    contacts::contacts_struct,
                                    time::Int64)

Check whether to "lock down" accommodation at block level.

Applies to oncampus accommodation only. \n
Class/society activities no longer permitted for affected students. \n

Inputs:
- `network_parameters::network_params`: network_params structure
- `intervention_trigger_input_data::intervention_params`: Values used to check if measures should be enforced/relaxed
- `contacts::contacts_struct`: contacts_struct structure
- `time`: Current time of the simulation

Outputs: None \n
Location: intervention\\_condition\\_affect\\_fns.jl
"""
function accom_lockdown_block_level!(network_parameters::network_params,
                                    intervention_trigger_input_data::intervention_data_feeds,
                                    contacts::contacts_struct,
                                    time::Int64)

   @unpack student_info = network_parameters
   @unpack block_member_list = contacts

   # Initialise intervention_params
   intervention_parameters = intervention_params(time_horizon = 14,
                                                  absolute_case_threshold = 3,
                                                  propn_case_threshold = 0.1,
                                                  release_condition_time_horizon = 14,
                                                  release_condition_absolute_cases = 0,
                                                  release_condition_propn_cases = 0)

   # Want to iterate over each household
   # Get number of halls and iterate over each one
   n_halls = length(block_member_list)
   for hall_itr = 1:n_halls

      # Get number of blocks & iterate over each one
      n_blocks = length(block_member_list[hall_itr])
      for block_itr = 1:n_blocks

         # Is at this level we want to enact control
         # Get all students in unit
         unit_resident_list = block_member_list[hall_itr][block_itr]
         n_students_in_unit = length(unit_resident_list)

         # Call function to do condition check
         condition_values::Array{Float64,1} = accom_lockdown_unit_check(unit_resident_list,
                                                                        intervention_parameters,
                                                                        student_info,
                                                                        time
                                                                        )

         # Call function to modify statuses as appropriate
         accom_lockdown_unit_alter_status!(unit_resident_list,
                                          intervention_parameters,
                                          student_info,
                                          time,
                                          condition_values
                                          )
      end
   end
end

"""
    accom_lockdown_floor_level!(network_parameters::network_params,
                                    intervention_trigger_input_data::intervention_data_feeds,
                                    contacts::contacts_struct,
                                    time::Int64)

Check whether to "lock down" accommodation at floor level.

Applies to oncampus accommodation only. \n
Class/society activities no longer permitted for affected students. \n

Inputs:
- `network_parameters::network_params`: network_params structure
- `intervention_trigger_input_data::intervention_params`: Values used to check if measures should be enforced/relaxed
- `contacts::contacts_struct`: contacts_struct structure
- `time`: Current time of the simulation

Outputs: None \n
Location: intervention\\_condition\\_affect\\_fns.jl
"""
function accom_lockdown_floor_level!(network_parameters::network_params,
                                    intervention_trigger_input_data::intervention_data_feeds,
                                    contacts::contacts_struct,
                                    time::Int64)

   @unpack student_info = network_parameters
   @unpack floor_member_list = contacts

   # Initialise intervention_params
   intervention_parameters = intervention_params(time_horizon = 14,
                                                  absolute_case_threshold = 3,
                                                  propn_case_threshold = 0.1,
                                                  release_condition_time_horizon = 14,
                                                  release_condition_absolute_cases = 0,
                                                  release_condition_propn_cases = 0)

   # Want to iterate over each household
   # Get number of halls and iterate over each one
   n_halls = length(floor_member_list)
   for hall_itr = 1:n_halls

      # Get number of blocks & iterate over each one
      n_blocks = length(floor_member_list[hall_itr])
      for block_itr = 1:n_blocks

         # Get number of floors & iterate over each one
         n_floors = length(floor_member_list[hall_itr][block_itr])
         for floor_itr = 1:n_floors

            # Is at this level we want to enact control
            # Get all students in unit
            unit_resident_list = floor_member_list[hall_itr][block_itr][floor_itr]
            n_students_in_unit = length(unit_resident_list)

            # Call function to do condition check
            condition_values::Array{Float64,1} = accom_lockdown_unit_check(unit_resident_list,
                                                                           intervention_parameters,
                                                                           student_info,
                                                                           time
                                                                           )

            # Call function to modify statuses as appropriate
            accom_lockdown_unit_alter_status!(unit_resident_list,
                                             intervention_parameters,
                                             student_info,
                                             time,
                                             condition_values
                                             )
         end
      end
   end
end

"""
    accom_lockdown_household_level!(network_parameters::network_params,
                                    intervention_trigger_input_data::intervention_data_feeds,
                                    contacts::contacts_struct,
                                    time::Int64)

Check whether to "lock down" accommodation at household level.

Applies to oncampus accommodation only. \n
Class/society activities no longer permitted for affected students. \n

Inputs:
- `network_parameters::network_params`: network_params structure
- `intervention_trigger_input_data::intervention_params`: Values used to check if measures should be enforced/relaxed
- `contacts::contacts_struct`: contacts_struct structure
- `time`: Current time of the simulation

Outputs: None \n
Location: intervention\\_condition\\_affect\\_fns.jl
"""
function accom_lockdown_household_level!(network_parameters::network_params,
                                          intervention_trigger_input_data::intervention_data_feeds,
                                          contacts::contacts_struct,
                                          time::Int64)

   @unpack student_info = network_parameters
   @unpack household_member_list = contacts

   # Initialise intervention_params
   intervention_parameters = intervention_params(time_horizon = 14,
                                                  absolute_case_threshold = 3,
                                                  propn_case_threshold = 0.1,
                                                  release_condition_time_horizon = 14,
                                                  release_condition_absolute_cases = 0,
                                                  release_condition_propn_cases = 0)

   # Want to iterate over each household
   # Get number of halls and iterate over each one
   n_halls = length(household_member_list)
   for hall_itr = 1:n_halls

      # Get number of blocks & iterate over each one
      n_blocks = length(household_member_list[hall_itr])
      for block_itr = 1:n_blocks

         # Get number of floors & iterate over each one
         n_floors = length(household_member_list[hall_itr][block_itr])
         for floor_itr = 1:n_floors

            # Get number of households & iterate over each one
            n_households = length(household_member_list[hall_itr][block_itr][floor_itr])
            for household_itr = 1:n_households

               # Is at this level we want to enact control
               # Get all students in unit
               unit_resident_list = household_member_list[hall_itr][block_itr][floor_itr][household_itr]
               n_students_in_unit = length(unit_resident_list)

               # Call function to do condition check
               condition_values::Array{Float64,1} = accom_lockdown_unit_check(unit_resident_list,
                                                                              intervention_parameters,
                                                                              student_info,
                                                                              time
                                                                              )

               # Call function to modify statuses as appropriate
               accom_lockdown_unit_alter_status!(unit_resident_list,
                                                intervention_parameters,
                                                student_info,
                                                time,
                                                condition_values
                                                )
            end
         end
      end
   end
end


# # Call function to do condition check
# # and modify statuses as appropraite
# function accom_lockdown_unit_check!(unit_resident_list::Array{Int64,1},
#                                        intervention_parameters::intervention_params,
#                                        student_info::Array{student_params,1},
#                                        time::Int64
#                                        )
# # Inputs:
# # unit_resident_list::Array{Int64,1} - The IDs of students in the accommodation unit being checked
# # intervention_parameters::intervention_params - Valeus used to check if measures should be enforced/relaxed
# # student_info::student_params - Information of each student in the system
# # time::Int64 - Current time of the simulation
#
# # Outputs:
# #  Direct modifications to student_info field on accommodation lockdown status
#
#    @unpack time_horizon, inactivation_length, absolute_case_threshold,
#             propn_case_threshold, release_condition_time_horizon,
#             release_condition_absolute_cases, release_condition_propn_cases  = intervention_parameters
#
#    # Initialise status for amending lockdown status
#    alter_accom_lockdown_status = "no"
#
#    # Get number of students in the accommodation unit of interest
#    n_students_in_unit = length(unit_resident_list)
#
#    # Check if accommodation block is in lockdown already
#    # Use first student, as status should be the same for all students in that accommodation
#    first_student_ID = unit_resident_list[1]
#    current_accom_lockdown_status = student_info[first_student_ID].household_info.lockdown_status
#
#    # Get count & propn. of household that have reported infection during relevant time horizon
#    # Uses time_of_reporting_infection field in student_info
#    recent_reported_cases = 0 # Initialise counter of recent reported cases
#    for student_itr = 1:n_students_in_unit
#       # For student, get ID and check if infected within the time horizon
#       student_ID = unit_resident_list[student_itr]
#       cutoff_date = max(1,time - time_horizon)  # Ensure cutoff_date is not 0 or below
#       if (student_info[student_ID].time_of_reporting_infection >= cutoff_date)
#          recent_reported_cases += 1  # Student reported infection within time horizon. Increment counter
#       end
#    end
#    propn_unit_reporting = recent_reported_cases/n_students_in_unit
#
#    if current_accom_lockdown_status == true
#       # accommodation in lockdown. See if it can be released.
#
#       # If rate is below a certain level & a certain number of cases
#       # release lockdown
#       if (recent_reported_cases < release_condition_absolute_cases) &&
#             (propn_unit_reporting < release_condition_propn_cases)
#          alter_accom_lockdown_status = "release"
#       end
#
#    else
#       # accommodation not in lockdown. Check if it should be
#
#       # If rate is above a certain level & a certain number of cases
#       # enforce lockdown
#       if (recent_reported_cases >= absolute_case_threshold) &&
#             (propn_unit_reporting >= propn_case_threshold)
#          alter_accom_lockdown_status = "enforce"
#       end
#    end
#
#    # If lockdown activated/released, alter status for each student in that unit
#    if alter_accom_lockdown_status != "no"
#       for student_in_unit_itr = 1:n_students_in_unit
#          # Get the ID of the student
#          current_student_ID = unit_resident_list[student_in_unit_itr]
#
#          # Update the status of accommodation lockdown for the student
#          if alter_accom_lockdown_status == "enforce"
#             student_info[current_student_ID].household_info.lockdown_status = true
#          elseif alter_accom_lockdown_status == "release"
#             student_info[current_student_ID].household_info.lockdown_status = false
#          end
#       end
#    end
# end
