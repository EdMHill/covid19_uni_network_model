"""
Purpose:
Define functions for use with the mass testing module

Supporting functions
- find_eligible_students

Main functions
- perform_mass_test!
"""


"""
Supporting functions
"""
# Get the IDs of nodes that would be eligable for mass testing
function find_eligible_nodes(student_info::Array{student_params,1},
                                 node_IDs_in_location::Array{Int64,1},
                                 CT_vars::contact_tracing_vars
                                 )

# Inputs:
# student_info::Array{student_params,1} - Fields of info for each student
# node_IDs_in_location::Array{Int64,1} - List of those residing in locale
# CT_vars::contact_tracing_vars - Variables used with contact tracing module

# Outputs:
# valid_node_IDs_vec - List of ID numbers for those that are eligible to be selected in mass testing

   # Get number of nodes in the locale of interest
   n_nodes_in_location = length(node_IDs_in_location)

   # Generate list of student IDs
   # Check if already reported infection previously. If so, should not include.
   # Find valid IDs
   total_valid_IDs = 0
   valid_node_IDs_flag = zeros(Int64,n_nodes_in_location)
   for node_itr = 1:n_nodes_in_location
      get_node_ID = node_IDs_in_location[node_itr]
      if (student_info[get_node_ID].time_of_reporting_infection == 0) ||
         ( (student_info[get_node_ID].time_of_reporting_infection > 0) && (CT_vars.Test_result_false_negative[get_node_ID] == true) )
              # Second condition gets those that reported infection and returned a negative test result

         # Update node validity vector
         valid_node_IDs_flag[node_itr] = 1
         total_valid_IDs += 1
      end
   end

   valid_ID_itr = 1
   valid_node_IDs_vec = zeros(Int64,total_valid_IDs)
   for node_itr = 1:n_nodes_in_location
      get_node_ID = node_IDs_in_location[node_itr]
      if valid_node_IDs_flag[node_itr] == 1
         valid_node_IDs_vec[valid_ID_itr] = get_node_ID

         # Increment index allocation
         valid_ID_itr += 1
      end
   end

   return valid_node_IDs_vec::Array{Int64,1}
end

# Get those isolating due to contact tracing
function apply_CT_isolation!(time::Int64,
                              Inds_to_be_contacted::Array{Int64,1},
                              replicate_ID::Int64,
                              output::sim_outputs,
                              states::student_states,
                              record_CT_isolated::Array{Int64,1})


   # Initial timepoint is for initial conditions
   # Set entry to be accessed in output arrays for this timestep
   output_time_idx = time + 1

   # Remove duplicates in list of nodes to be contact traced (with current_node_ID as index case)
   unique!(Inds_to_be_contacted)

   # Perform contact tracing immediately
   # Get number of recallable contacts
   n_recallable_contacts = length(Inds_to_be_contacted)
   output.num_CT[output_time_idx,replicate_ID] += n_recallable_contacts

   for recallable_contact_idx = 1:n_recallable_contacts
      recallable_contact_ID = Inds_to_be_contacted[recallable_contact_idx]

      # In node entry record as isolating due to CT
      record_CT_isolated[recallable_contact_ID] = 1

      # Check if individual will adhere to guidance
      # If so, they self-isolate
      if states.hh_isolation[recallable_contact_ID] == 1
          states.timeisol_CTcause[recallable_contact_ID] = 1
      end
   end

   return nothing
end



# Carry out tests
function carry_out_mass_tests_in_location(time::Int64,
                                             replicate_ID::Int64,
                                             test_coverage::Float64,
                                             n_valid_node_IDs::Int64,
                                             valid_node_IDs::Array{Int64,1},
                                             student_info::Array{student_params,1},
                                             CT_vars::contact_tracing_vars,
                                             CT_parameters::CT_params,
                                             network_parameters::network_params,
                                             states::student_states,
                                             contacts::contacts_struct,
                                             household_contacts_per_node::Array{Int64,1},
                                             output::sim_outputs,
                                             mass_testing_parameters::mass_testing_params,
                                             mass_test_ID::Int64,
                                             rng::MersenneTwister,
                                             record_test_postive::Array{Int64,1},
                                             record_hh_isolated::Array{Int64,1},
                                             record_CT_isolated::Array{Int64,1},
                                             rehouse_strat_active::Bool)
# Inputs:
# time::Int64 - Timestep of simulation being performed
# replicate_ID::Int64 - Identifies replicate of current simulation batch being performed
# test_coverage::Float64 - Proportion of eligible students, in given location, to request to be tested
# n_valid_node_IDs::Int64, valid_node_IDs::Array{Int64,1} - List of IDs that are selectable
# student_info::Array{student_params,1} - Fields of info for each student
# CT_vars::contact_tracing_vars - All parameter structures as described
# CT_parameters::CT_params, network_parameters::network_params, states::student_states, contacts::contacts_struct, output::sim_outputs
# mass_testing_parameters::mass_testing_params
# mass_test_ID::Int64 - ID number for the instance of mass testing that is taking place.
# record_test_postive, record_hh_isolated, record_CT_isolated - Node level vectors to track whether tested positive and/or isolated during this mass testing occurance
# rehouse_strat_active::Bool - Intervention where those who are symptomatic are rehoused/completely isolate with no contacts.

# Outputs:
# All changes made directly (in place) to the input variables

   # Unpack relevant variables
   @unpack perform_CT_from_infector, prob_backwards_CT = CT_parameters
   @unpack infected_by, new_rehoused = output
   @unpack designated_test_times, on_campus_coverage_propn,
            off_campus_coverage_propn, asymp_test_detection_prob_vec,
            n_mass_tests_performed,
            n_tests_performed, n_tests_positive,
            n_all_isolations_caused, n_hh_isolations_caused, n_CT_isolations_caused,
            n_prev_infected_tested = mass_testing_parameters


   # If coverage above zero, perform tests
   if (test_coverage > 0)

      # Get number tested and sample
      n_tested = ceil(Int64,test_coverage*n_valid_node_IDs)
      tested_nodes = rand(valid_node_IDs,n_tested)

      for location_node_itr = 1:n_tested

         # Get node ID
         current_node_ID = tested_nodes[location_node_itr]

         # Check if node would particpate in the test
         if (CT_vars.Engage_with_CT[current_node_ID] == true)
            # Increment the numbers of tests performed in this mass test instance by 1
            n_tests_performed[replicate_ID][mass_test_ID] += 1

            # Check if tested individual had been infected previously
            # If so, increment prev_infected_tested counter
            if (states.timesymp[current_node_ID] == -1)
               n_prev_infected_tested[replicate_ID][mass_test_ID] += 1
            end

            # Check if in pre-symptomatic phase, or asymptomatic in subsequent phase
            if (states.timelat[current_node_ID]>0) || (states.timeinf[current_node_ID]>0) ||
               ( (states.timesymp[current_node_ID]>0) && (states.asymp[current_node_ID]==1) )

               # Determine whether test result will return a negative outcome
               # - Get time since individual became infected
               # - Given time since infected, look up probability case will return negative test result
               tot_time_inf = time - states.acquired_infection[current_node_ID]

               # Get relevant sensitivity value based on time since infection
               # Cap total infection time at length of detection prob vector
               if tot_time_inf > length(asymp_test_detection_prob_vec)
                  tot_time_inf = length(asymp_test_detection_prob_vec)
               end

               # Get relevant test sensitivity value based on time since infection
               if tot_time_inf == 0
                  test_detection_prob = 0.
               else
                  test_detection_prob = asymp_test_detection_prob_vec[tot_time_inf]
               end

               # Bernoulli trial to determine if false negative returned
               if rand(rng) < test_detection_prob
                  # If not false negative, gather contacts
                  CT_vars.Inds_to_be_contacted[current_node_ID] = Int64[] # Initialise vector to store contacts
                  if (CT_vars.Engage_with_CT[current_node_ID] == true)

                     trace_node!(current_node_ID,time,CT_vars,contacts,CT_parameters,network_parameters,rng)

                     # if we are doing "backward contact tracing"
                     # some small chance the infector is included in this
                     # don't try to backwards trace the initial infections
                     if perform_CT_from_infector == true
                        if (rand(rng)<prob_backwards_CT) && (infected_by[current_node_ID]!=-1)
                             append!(CT_vars.Inds_to_be_contacted[current_node_ID],infected_by[current_node_ID])
                             CT_vars.Recall_infector[current_node_ID] = 1
                        end
                     end

                     # Apply isolation due to CT
                     apply_CT_isolation!(time,
                                          CT_vars.Inds_to_be_contacted[current_node_ID],
                                          replicate_ID,
                                          output,
                                          states,
                                          record_CT_isolated)

                     # Perform forwards CT from infector, if infector has been identified
                     # and such a policy is active
                     if perform_CT_from_infector == true
                        if CT_vars.Recall_infector[student_itr] == 1
                             recall_infector_count += 1

                             # Assign infector_ID to variable
                             infector_ID = infected_by[current_node_ID]

                             # Get recallable contacts of infector_ID
                             forwardCT_from_infector!(infector_ID,
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

                              # Apply isolation due to CT
                              apply_CT_isolation!(time,
                                             CT_vars.Inds_to_be_contacted[infector_ID],
                                             replicate_ID,
                                             output,
                                             states,
                                             record_CT_isolated)
                        end
                     end
                  end

                  # In node entry record as positive
                  record_test_postive[current_node_ID] = 1

                  # Check if enter isolation for designated period
                  if (states.hh_isolation[current_node_ID]==1)

                     # Enter isolation due to poisitive test (though non-symptomatic)
                     states.asymp_timeisol[current_node_ID] = 1

                     # Supercedes household isolation. Reset that counter
                     states.timeisol[current_node_ID] = 0

                      # Update number of persons in household confirmed infected
                      current_node_household_ID = student_info[current_node_ID].household_info.household_ID
                      CT_vars.Cases_per_household_pos_or_unknown[current_node_household_ID] += 1

                      # Assign time of reporting to field in student parameter type
                      student_info[current_node_ID].time_of_reporting_infection = time

                      # Increment the numbers of tests performed in this mass test instance by 1
                      n_tests_positive[replicate_ID][mass_test_ID] += 1

                      # If an available option, no contacts state entered if student is living on-campus
                     # Will last until end of symptoms (irrespective of test result)
                     if (rehouse_strat_active == true) &&
                             (student_info[current_node_ID].household_info.on_campus_accom == true)
                         student_info[current_node_ID].no_contacts_status = true

                         # Update rehoused output variable if living in communal
                         # bathroom type accommodation
                         if (student_info[current_node_ID].household_info.ensuite_flag == false)
                             altered_time_idx = time + 2  # +1 for first entry being time 0. +1 for rehousing effectively taking place the next timestep.
                             new_rehoused[altered_time_idx,replicate_ID] += 1
                         end
                     end
                  end

                  # Check if household contacts isolate
                  # Irrespective of whether index case self-isolates,
                  # adherent members of their household may also isolate.
                  for hh = 1:household_contacts_per_node[current_node_ID]
                     contact_ID = contacts.household_contacts[current_node_ID][hh]

                     # In node entry record as household isolated
                     record_hh_isolated[contact_ID] = 1

                     # If adheres to guidance, isolate
                     if (states.hh_isolation[contact_ID]==1) &&
                        (states.symp_timeisol[contact_ID]==0) && # Individual not already symptomatic themselves
                        (states.asymp_timeisol[contact_ID]==0) # Individual not already isolating as asymptomatic & tested positive
                          states.timeisol[contact_ID] = 1

                          # Update time of latest confirmed release from household isolation
                          CT_vars.Time_of_hh_isolation_release[contact_ID] = time + states.household_isoltime
                     end
                  end
               end
            end
         end
      end
   end

   return nothing
end

"""
Main functions
"""
# Run the primary mass testing function
function perform_mass_test!(mass_testing_parameters::mass_testing_params,
                                          time::Int64,
                                          replicate_ID::Int64,
                                          states::student_states,
                                          CT_vars::contact_tracing_vars,
                                          CT_params::CT_params,
                                          network_parameters::network_params,
                                          contacts::contacts_struct,
                                          output::sim_outputs,
                                          household_contacts_per_node::Array{Int64,1},
                                          rng::MersenneTwister,
                                          rehouse_strat_active::Bool)

   # Unpack student information and mass test parameter type
   @unpack n_students, n_students_on_campus,
            on_campus_student_IDs, off_campus_student_IDs,
            student_info = network_parameters

   # Unpack fields from mass testing parameters
   @unpack designated_test_times, on_campus_coverage_propn,
            off_campus_coverage_propn, asymp_test_detection_prob_vec,
            n_mass_tests_performed,
            n_tests_performed, n_tests_positive,
            n_all_isolations_caused, n_hh_isolations_caused, n_CT_isolations_caused,
            n_prev_infected_tested = mass_testing_parameters

   # Get number of off-campus students
   n_students_off_campus = length(off_campus_student_IDs)

   # Get index for accessing mass testing associated vectors.
   next_mass_test_count = n_mass_tests_performed[replicate_ID] + 1

   # Check if mass testing due to be carried out
   if next_mass_test_count <= length(designated_test_times)
      next_mass_test_time = designated_test_times[next_mass_test_count]

      if next_mass_test_time == time
         # Mass testing to be carried out

         # Generate list of student IDs
         valid_student_IDs_on_campus = find_eligible_nodes(student_info,
                                                               on_campus_student_IDs,
                                                               CT_vars)
         valid_student_IDs_off_campus = find_eligible_nodes(student_info,
                                                                  off_campus_student_IDs,
                                                                  CT_vars)

         # Get number of students eligible in each locale
         n_valid_student_IDs_on_campus = length(valid_student_IDs_on_campus)
         n_valid_student_IDs_off_campus = length(valid_student_IDs_off_campus)

         # Get coverage for this instance of mass testing
         on_campus_coverage = on_campus_coverage_propn[next_mass_test_count]
         off_campus_coverage = off_campus_coverage_propn[next_mass_test_count]

         # Set up variables to track whether each student tested positive and/or
         # should household isolate and/or isolate through CT
         record_test_postive = zeros(Int64,n_students)
         record_hh_isolated = zeros(Int64,n_students)
         record_CT_isolated = zeros(Int64,n_students)
         # Draw up lists of students to receive a test & check test result
         # On-campus resident check first. Then off-campus.
         carry_out_mass_tests_in_location(time,
                                             replicate_ID,
                                             on_campus_coverage,
                                             n_valid_student_IDs_on_campus,
                                             valid_student_IDs_on_campus,
                                             student_info,
                                             CT_vars,
                                             CT_params,
                                             network_parameters,
                                             states,
                                             contacts,
                                             household_contacts_per_node,
                                             output,
                                             mass_testing_parameters,
                                             next_mass_test_count,
                                             rng,
                                             record_test_postive,
                                             record_hh_isolated,
                                             record_CT_isolated,
                                             rehouse_strat_active)
         carry_out_mass_tests_in_location(time,
                                          replicate_ID,
                                          off_campus_coverage,
                                          n_valid_student_IDs_off_campus,
                                          valid_student_IDs_off_campus,
                                          student_info,
                                          CT_vars,
                                          CT_params,
                                          network_parameters,
                                          states,
                                          contacts,
                                          household_contacts_per_node,
                                          output,
                                          mass_testing_parameters,
                                          next_mass_test_count,
                                          rng,
                                          record_test_postive,
                                          record_hh_isolated,
                                          record_CT_isolated,
                                          rehouse_strat_active)
         # Compute the number of additional isolations that occurred due to mass testing
         # Increment the numbers of tests performed in this mass test instance by 1
         for student_itr = 1:n_students
            if (record_test_postive[student_itr] == 0)
               # Check not tested positive themselves first.
               # Then see if isolated due to household contact or CT of asymptomatic case
               if (record_hh_isolated[student_itr] == 1)
                  n_hh_isolations_caused[replicate_ID][next_mass_test_count] += 1
               end
               if (record_CT_isolated[student_itr] == 1)
                  n_CT_isolations_caused[replicate_ID][next_mass_test_count] += 1
               end
               isolation_reasons_sum = record_hh_isolated[student_itr] + record_CT_isolated[student_itr]
               if (isolation_reasons_sum > 0)
                  n_all_isolations_caused[replicate_ID][next_mass_test_count] += 1
               end
            end
         end

         # Update counter of mass tests performed
         n_mass_tests_performed[replicate_ID] += 1
      end
   end

   return nothing
end
