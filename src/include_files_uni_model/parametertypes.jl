#=
Purpose:
Parameter types to be used with the network model for universities

- class_params: Info on class size permitted and whether take place face-to-face.
- society_params: Info on society (if active, schedule, size, member list etc)
- household_params: Household parameter type to have information on the household the student is in
- student_params: Track individual level details.
- CTParams: Contact tracing related
- contacts_struct: Used for contact structures, that are also used outside of contact tracing
- contact_tracing_vars: Vars only needed for contact tracing
- infection_params: Parameters related to the transmission process
- network_params: Quantities to construct the contacts & stores the node properties
- class_generation_params: Used for generating a "population" of student classes & staff teams
- sim_outputs: Structures to hold data to be collected from simulation run
- student_states: Used for the status of each student
- intervention_data_feeds: The data streams that may be of use in implementing trigger-based interventions
- intervention_params: Used to implement the required intervention criteria.
- mass_testing_params: Fields associated with mass testing intervention
=#
#-------------------------------------------------------------------------------

"""
    Structure: class_params

Course parameter type to have information on the course being run with face-to-face classes.

Location: parametertypes.jl
"""
@with_kw mutable struct class_params
   f2f_activity::Bool = true      # For putting in place online learning for the course, can set flag to false
   class_inactivation_time::Int64 = 0  # Time class has not been permitted to have f2f teaching
end

"""
    Structure: society_params

Society parameter type to have information on the society being active.

Location: parametertypes.jl
"""
@with_kw mutable struct society_params
   f2f_activity::Bool = true
   schedule::Array{Int64,1} = Array{Int64,1}(undef,0)
   society_type::Int64 = 0
   society_inactivation_time::Int64 = 0
   n_members::Int64 = 0
   member_list::Array{Int64,1} = Array{Int64,1}(undef,0)
end

"""
    Structure: household_params

Household parameter type to have information on the household the student is in:
- Whether on-campus of off-campus, location
- Particularly relevant to on-campus accommodation: Halls >  Blocks > Floors > Households

Location: parametertypes.jl
"""
@with_kw mutable struct household_params
   household_ID::Int64 = 0        # Household the individual has been assigned to. ID for the global system
   on_campus_accom::Bool = true   # Whether accommodation is on campus (true) or off campus (false)
   location::String = "Campus"    # name of location where student is living
   hall_ID::Int64 = 0             # Relevant to on-campus accommodation (likewise for block_ID, floor_ID, household_ID_within_block)
   block_ID::Int64 = 0
   floor_ID::Int64 = 0
   household_ID_within_block::Int64 = 0
   lockdown_status::Bool = false  # Specifies if accommodation is on lockdown or not. If so, non-accom f2f contact not permitted.
   ensuite_flag::Bool = false    # Whether the accommodation is ensuite bathroom (true) or communal bathroom (false)
end

"""
    Structure: student_params

Track individual student details.

Location: parametertypes.jl
"""
@with_kw mutable struct student_params
   would_attend_f2f_classes::Int64        # Flag variable. If 1, individual would attend f2f if they are active. If 0, they will not.
   cohort_ID::Int64                       # Specify subjects and year groups. Also classifies staff members.
   class_ID::Int64                        # The ID of the class (within that cohort)
   society_IDs::Array{Int64,1} = Int64[]  # Social society group identifier. If empty, not in any societies.
   household_info::household_params = household_params()
   time_of_reporting_infection::Int64 = 0 # If student is symptomatically infected and it is reported, the time the case report occurs.
                                          # Asymptomatic infection also logged on day of positive test.
                                          # If never reporting infection, keeps value 0.
   transrisk_household::Float64 = 0.      # Secondary attack rate within household were that individual the index case
   transrisk_cohort::Float64 = 0.         # Transmission risk to cohort contacts
   transrisk_society_sports::Array{Float64,1} = [0.,0.] # Transmission risk to those in [society,sports club]
   transrisk_dynamic_social::Float64 = 0. # Transmission risk to dynamic social contacts
   no_contacts_status::Bool = false       # Defines if student is in rehoused/full quarantine state
end

"""
    Structure: CT_params

Parameter structure relating to contact tracing.

Location: parametertypes.jl
"""
@with_kw mutable struct CT_params
# used for parameters relating to contact tracing

   # For those adhering to self-isolating, engagement with contact tracing aspect
   # If 1, as well as self-isolating, all adhering individuals do give contacts
   CT_engagement::Float64 = 0.7

   # Set time delay. 0 corresponds to result on day of reporting.
   CT_delay_until_test_result_pmf::Array{Float64,1} = [0., 0., 1.,]

   # Set number of days before symptoms CT will attempt to capture
   CT_days_before_symptom_included::Int64 = 2

   # Propotions of tests that correctly detect infection (sensitivity)
   # Entry per day since infected
   test_detection_prob_vec::Array{Float64,1} = [0.,0.11,0.59,0.785,0.83,0.845,0.84,0.82,0.79,0.76,  # Days 1-10
                                                0.72,0.68,0.64,0.59,0.54,0.485,0.445,0.405,0.37,0.335, # Days 11-20
                                                0.30,0.27,0.24,0.22,0.20,0.18,0.16,0.15,0.14,0.13] # Days 21-30

   # Amount of time spent in isolation if contact traced
   CT_caused_isol_limit::Int64 = 14

   # Set up recall of dynamic contacts probability
   # Proportion of contacts remembered x days ago
   dynamic_contacts_recalled_propn::Array{Float64,1} = [0.5,0.4,0.3,0.2,0.1]
   accom_dynamic_contacts_recalled_propn::Array{Float64,1} = [1.,1.,1.,1.,1.]
   social_contacts_recalled_propn::Array{Float64,1} = [1.,1.,1.,1.,1.]

   # proportion of people that can identify their infector (set to 0 for no backwards CT)
   prob_backwards_CT::Float64 = 0.1

   # Parameters for performing forward contact tracing from identified infectors
   perform_CT_from_infector::Bool = false
   infector_engage_with_CT_prob::Float64 = 1.0  # For those not complying with other isolation guidance,
                                       # probability they do if idenfitied as possible infector
                                       # Note, those that are set to adhere to isolation guidance
                                       # are also assuemd fully compliant with CT from infector measures

   work_or_study_group_CT_memory::Int64 = 7 # how many days to look over when looking for clusters

   work_or_study_group_CT_threshold::Float64 = 0.5 # what proportion of employees need to be infected to trigger workplace closure

   time_WC::Int64 = 14 # how many days to close workplaces for
end

"""
    Structure: contacts_struct

Used for contact structures, that are also used outside of contact tracing.

Location: parametertypes.jl
"""
@with_kw mutable struct contacts_struct

   # sizes to initialise arrays
   n_students::Int64 = 0
   endtime::Int64 = 0

   # Stores work/study group contacts
   class_contacts::Array{Array{Int64,1},1} = Array{Array{Int64,1},1}(undef,n_students)

   # Store contacts with students in same cohort, but not same class
   cohort_contacts::Array{Array{Int64,1},1} = Array{Array{Int64,1},1}(undef,n_students)

   # Flag arrays to keep a log of whether a node is in class/at work
   # & in isolation each timestep
   # Used in each replicate
   daily_record_inclass::Array{Int64,2} = zeros(Int64,endtime,n_students)
   daily_record_inisol::Array{Int64,2} = zeros(Int64,endtime,n_students)

   # Flag arrays to keep a log of whether a person attends society group
   # Can particpate in more than one society, so use a 3D array
   # In array, row per timestep, column per society, slice per node.
   # Used in each replicate
   daily_record_atsociety::Array{Int64,3} = zeros(Int64,endtime,0,n_students)


   # # Total social contacts made by node
   # social_contacts_per_node::Array{Int64,1} = zeros(Int64,n_students)
   #
   # # Per node, a record of social contacts made on each day
   # workday_social_contacts_by_day::Array{Array{Int64,1},2} = Array{Array{Int64,1},2}(undef,0,0)
   # nonworkday_social_contacts_by_day::Array{Array{Int64,1},2} = Array{Array{Int64,1},2}(undef,0,0)

   # Per node, a record of contacts made in each society group
   society_contacts::Array{Array{Int64,1},2} = Array{Array{Int64,1},2}(undef,n_students,0)

   # Per node, a record of dynamic social contacts made on each day
   dynamic_social_contacts::Array{Array{Int64,1},2} = Array{Array{Int64,1},2}(undef,0,0)

   # Per node, household contacts
   household_contacts::Array{Array{Int64,1},1} = Array{Int64,1}[]

   # Per node, a record of dynamic accommodation contacts made on each day
   # Relevant to those in on-campus accommodation only.
   # Three dimensions of vectors. Row per timestep, column per accom level, slice per student
   dynamic_accommodation_contacts::Array{Array{Int64,1},3} = Array{Array{Int64,1},3}(undef,0,0,0)

   # Students in each household < floor < block < hall
   household_member_list::Array{Array{Array{Array{Array{Int64,1},1},1},1},1} = Array{Array{Array{Array{Array{Int64,1},1},1},1},1}(undef,0)
   hall_member_list::Array{Array{Int64,1},1} = Array{Array{Int64,1},1}(undef,0)
   block_member_list::Array{Array{Array{Int64,1},1},1} = Array{Array{Array{Int64,1},1},1}(undef,0)
   floor_member_list::Array{Array{Array{Array{Int64,1},1},1},1} = Array{Array{Array{Array{Int64,1},1},1},1}(undef,0)
end

"""
    Structure: contact_tracing_vars

Parameter structure for contact tracing variables.

Location: parametertypes.jl
"""
@with_kw mutable struct contact_tracing_vars
# only needed for contact tracing

   # sizes to initialise arrays
   n_students::Int64 = 0
   endtime::Int64 = 0
   n_households::Int64 = 0

   # The number of days prior to symptoms that each node remembers
   relevant_prev_days_for_CT::Array{Int64,1} = zeros(Int64,n_students)

   # Vector of vectors for storing IDs of those to be contacted in CT
   Inds_to_be_contacted::Array{Array{Int64,1},1} = Array{Array{Int64,1},1}(undef,n_students)

   # vector tracking symptomatic cases (positive confirmed or untested)
   Cases_per_household_pos_or_unknown::Array{Int64,1} = zeros(Int64,n_households)

   # Vector tracking the latest isolation release time due to household member cases
   hh_isolation_release_time::Array{Int64,1} = zeros(Int64,n_students)

   # Variables for waiting for test results, set at -1 until activated
   Time_to_test_result::Array{Int64,1} = -1*ones(Int64,n_students)

   # Boolean vector to store whether a false negative test result would be returned
   # Populated in each replicate
   Test_result_false_negative::Array{Bool,1} = Array{Bool,1}(undef,n_students)

   # Boolean vector to store if individual will engage with CT or not
   # Populated in each replicate
   Engage_with_CT::Array{Bool,1} = Array{Bool,1}(undef,n_students)

   # Array to keep track of whether an infected recalls their infector
   Recall_infector::Array{Int64,1} = zeros(Int64,n_students)

   # Delay before test result is returned
   CT_delay_until_test_result::Array{Int64,1} = zeros(Int64,n_students)
end

"""
    Structure: infection_params

Parameter structure with variables relating to the infection process, wih many of these fields reset in the configuration.

Location: parametertypes.jl
"""
@with_kw struct infection_params
   #Specify number of cohorts in use
   n_cohorts = 84

   # Scaling to infectiousness applied to those symptomatic
   # (if not isolating, value less than 1 could reflect being more cautious in public settings)
   iso_trans_scaling::Float64 = 1.

   # Also add in the possibility that individuals are asymptomatic.
   # In that case, their infection potential is lower but they dont isolate.
   # asymp_trans_scaling::Float64 = 0.165
   asymp_trans_scaling_dist::Uniform{Float64} = Uniform(0.3,0.7)

   # probability of transmission within a household, based on secondary attack rate
   # Can differ per household group. Values for household size of 2, 3, 4, 5 or more.
   transrisk_household_group_mean::Array{Float64,1} = [0.48,0.40,0.33,0.22]
   transrisk_household_group_sd::Array{Float64,1} = [0.06,0.06,0.05,0.05]
   transrisk_household_18_34_group_mean::Float64 = 0.34
   transrisk_household_18_34_group_sd::Float64 = 0.05

   # Baseline risk for cohort dependent transmission,
   # set relative to household transmission
   transrisk_cohort_mean::Array{Float64,1} = ((transrisk_household_18_34_group_mean/0.5)*0.355)*ones(n_cohorts)
   transrisk_cohort_sd::Array{Float64,1} = [transrisk_cohort_mean[ii]*(transrisk_household_18_34_group_sd/transrisk_household_18_34_group_mean) for ii=1:length(transrisk_cohort_mean)]

   # Society & sports club contacts transmission risk.
   # Pair of entries: [societies, sports clubs]
   society_sports_transrisk_mean::Array{Float64,1} = [0.2414,0.34]
   society_sports_transrisk_sd::Array{Float64,1} = society_sports_transrisk_mean.*(transrisk_household_18_34_group_sd/transrisk_household_18_34_group_mean)

   # Baseline risk for social transmission, set relative to household transmission
   transrisk_dynamic_social_mean::Float64 = (transrisk_household_18_34_group_mean/0.5)*0.355
   transrisk_dynamic_social_sd::Float64 = transrisk_dynamic_social_mean*(transrisk_household_18_34_group_sd/transrisk_household_18_34_group_mean)

   # For settings following covid guidance, scale the transmission risk
   # Zero would mean complete removal of transmission risk
   # One would mean no effect on the transmission risk
   CS_scale_transrisk::Array{Float64,1} = 0.2*ones(length(transrisk_cohort_mean))

   # Scale the infectiousness of all contacts
   suscep_scaling::Float64 = 1.

   # probability of being asymptomatic
   probasymp_dist::Uniform{Float64} = Uniform(0.5,0.8)

   # Flag variable indicating if self-isolation is an active measure
   isolation::Int64 = 1
   symp_isoltime::Int64 = 10
   asymp_isoltime::Int64 = 10
   household_isoltime::Int64 = 14

   # Set proportion of population who will adhere
   adherence::Float64 = 0.7

   # Average social contacts on a workday & nonworkday
   n_social_mean_workday::Int64 = 1
   n_social_mean_nonworkday::Int64 = 5

   # Distribution of incubation period
   d_incub::Distribution = Erlang(6,0.88)

   # Distribution of infectiousness
   dist_infectivity::Array{Float64,1} =  [0.0369, 0.0491, 0.0835, 0.1190, 0.1439, 0.1497, 0.1354, 0.1076, 0.0757, 0.0476, 0.0269, 0.0138, 0.0064, 0.0044] # Corrected He (4 day pre-symptomatic period) & 10 day symptomatic period
   # dist_infectivity::Array{Float64,1} = [0.0379, 0.0504, 0.0857, 0.1220, 0.1475, 0.1535, 0.1388, 0.1103, 0.0776, 0.0488, 0.0276] # Corrected He
   # dist_infectivity::Array{Float64,1} = [0.1847,0.2371,0.1968,0.1396,0.0911,0.0566,0.0339,0.0199,0.0114]/sum([0.1847,0.2371,0.1968,0.1396,0.0911,0.0566,0.0339,0.0199,0.0114])

   # Distribution of delay in reporting symptoms (for those that do/eventually adhere)
   # Short delays. E.g. Symptom onset while at work & isolate next day = 1 day delay
   delay_adherence_pmf::Array{Float64,1} = [1.,0.,0.] # Note first entry is 0 day delay,
                                                         # second entry 1 day delay etc

   # Distribution of delay in household infection (for those that will be infected)
   # First entry corresponds to 0 days, second entry to 1 day, etc.
   #delay_household_infection_pmf::Array{Float64,1} = [1.,0.,0.,0.,0.,0.]
   delay_household_infection_pmf::Array{Float64,1} = dist_infectivity./sum(dist_infectivity)

   # At beginning of simulation, the proportion of student population that has
   # been infected previously
   recov_propn = 0.1
end

"""
    Structure: network_params

Parameter structure with variables related to the network generation.

Location: parametertypes.jl
"""
@with_kw mutable struct network_params

   # Number of students in the system
   n_students::Int64 = 0

   # Number of students resident on campus
   n_students_on_campus::Int64 = 7155

   # Log of student IDs that are resident on-campus and off-campus
   on_campus_student_IDs::Array{Int64,1} = Int64[]
   off_campus_student_IDs::Array{Int64,1} = Int64[]

   # Information on each individual
   student_info::Array{student_params,1} = Array{student_params,1}(undef,0)

   # Vector of vectors for class/work team size
   # Vector per type of work (e.g. student class, lecturer, support staff)
   class_sizes::Array{Array{Int64,1},1} = Array{Int64,1}[]

   # Method used to generate contacts (can be "configuration" or "ER")
   network_generation_method::String = "configuration"

   # Lowest denomination for on-campus accommodation
   lowest_campus_denomination = "household"

   # Connection probabilities with others in society & sports clubs
   # First entry for societies. Second entry for sports clubs
   # Consistent for everyone
   prob_social_contact::Array{Float64,1} = [0.05, 0.1]

   # Household size distributions
   household_size_distribution::Array{Float64,1} = [0.314, 0.423, 0.075, 0.025, 0.006, 0.002]/sum([0.314, 0.423, 0.075, 0.025, 0.006, 0.002])

   # Offcampus household size distribution & locations
   offcampus_student_household_size_distribution::Array{Float64,1} = pdf.(Distributions.LogNormal(0.979,0.576),collect(1:20))
      # Gives distribution of household sizes from 1 through to 20
   offcampus_student_household_location_dist::Array{Float64,1} = [0.4, 0.1, 0.5]
   offcampus_student_location_names::Array{String,1} = ["Coventry","Kenilworth","Leamington"]

   # Work/Study type contact probability
   prob_worktype_contact::Array{Float64,2} = Array{Float64,1}(undef,0)

   # Degree distribution for contacts within class
   dd_within_class::Array{Float64,1} = Array{Float64,1}(undef,0)

   # Degree distribution of contacts in different classes
   class_degree_distribution::Array{Distribution,2} = [Distributions.LogNormal(1.896,1.233)]

   # Probability of making contact with other cohort members compared to class members
   between_class_contact_probs::Array{Float64,1} = [0.0]

   # Distribution of social group sizes
   social_group_size_distribution::Distribution = Distributions.LogNormal(log(6),1.056)

   # Distribution of social contacts per day
   social_workday_dd::Distribution = Distributions.LogNormal(1.397,1.27)
   social_nonworkday_dd::Distribution = Distributions.LogNormal(1.536,1.153)

   # Maximum contacts allowed for social group
   max_contacts_social::Int64 = 100
   max_contacts_social_dynamic::Int64  = 100

   # Probability of making contacts with f-o-f opposed to others
   friend_of_friend_prob::Float64 = 0.5

   # mean number & standard deviation of dynamic contacts for each study type
   dynamic_conts_mean::Array{Float64,1} = Array{Float64,1}(undef,0)
   dynamic_conts_sd::Array{Float64,1} = Array{Float64,1}(undef,0)

   # Degree distribution of dynamic social contacts
   dynamic_social_contact_degree_distribution::Array{Distribution,2} = [Distributions.LogNormal(1.646,1.211) Distributions.LogNormal(1.590,1.128)]

   # Associated with contacts within broader accommodation units
   contact_prob_floor_level::Float64 = 0.1
   contact_prob_block_level::Float64 = 0.05
   contact_prob_hall_level::Float64 = 0.

   # Info on each workplace/society
   class_info::Array{Array{class_params,1},1} = Array{class_params,1}[]
   society_info::Array{society_params,1} = Array{society_params,1}(undef,0)

   # Face-to-face activity flags
   cohort_f2f_study_active::Array{Bool,1} = Array{Bool,1}(undef,0)
   society_f2f_active::Array{Bool,1} = Array{Bool,1}(undef,0)

   # # Declares wether setting may have covid-secure (CS) status
   # # If false, all settings are non-CS
   # CS_active_flag::Bool = false
end

"""
    Structure: class_generation_params

Parameter structure with variables related to student class generation.

Location: parametertypes.jl
"""
@with_kw struct class_generation_params

# Total Number of Students 27,278
# Undergraduate	15,998
# Postgraduate	9,799
# Exchange/ Visiting/ Study Abroad/ IFP**/ Industry 1,481

# Undergraduate/postgraduate split: 62.01%/37.99%

# Faculty Populations (as % of total student numbers)
# Arts	85% undergraduates, 15% postgraduates	12.40%
#   As propn overall:  10.54% & 1.86%
# Social Sciences	54% undergraduates, 46% postgraduates	44.59%
#   As propn overall:  24.08% & 20.51%
# Medicine	54% undergraduates, 46% postgraduates	5.73%
#   As propn overall:  3.09% & 2.64%
# Science & Engineering (without Medicine) 37.28%
#   As propn overall:  24.01% & 13.27%
#        Note the breakdown for Science & Engineering (without Medicine) was calculated using
#        Science, Engineering and Medicine	63% undergraduates, 37% postgraduates	43.01%
#          As propn overall:  27.10% & 15.91%

   n_cohorts::Int64 = 3
   attendence_propns::Array{Float64,2} = ones(1,3)
   classtype_proportion::Array{Float64,2} = [0.5 0.5 0.5]
   class_size_mean::Array{Float64,2} = [25. 25. 5.]
   class_size_sd::Array{Float64,2} = [0. 0. 1.]
end

"""
    Structure: society_generation_params

Parameter structure with variables related to generation of societies.

Location: parametertypes.jl
"""
@with_kw struct society_generation_params
   society_types::Int64 = 2
   society_type_proportion::Array{Float64,1} = [0.5, 0.5]
   societies_joined_per_person_dist::Array{Float64,1} = [0.5,0.4,0.025,0.025,0.025,0.025]
end

"""
    Structure: sim_outputs

Parameter structure for outputs saved from the simulations.

Location: parametertypes.jl
"""
@with_kw mutable struct sim_outputs
   endtime::Int64 = 0
   countfinal::Int64 = 0
   n_students::Int64 = 0
   # n_initial_asymp::Int64 = 0
   # n_initial_symp::Int64 = 0

   # 2D outputs
   numlat::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # number latently infected
   numinf::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # number infectious
   numrep::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # number reporting infection
   prevlat::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # Prevalence for latently infected
   prevsymp::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # Prevalence for symptomatic infectious (post pre-symptomatic phase)
   prevasymp::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # Prevalence for asymptomatic infectious
   prevpresymp::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # Prevalence of pre-symptomatic symptomatic infectious
   prevrec::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # Prevalence for recovereds
   newinf::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # number of new infections
   newasymp::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # number of newly infected asymptomatics
   atworkinf::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # number of people that are newly infected while in study/work setting
   atworkasymp::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # number of peoplethat are newly infected asymptomatics while in study/work setting
   num_CT::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # total number of recallable contacts
   num_infected::Array{Int64,2} = zeros(Int64,n_students,countfinal) # number of infections caused by that node
   social_dynamic_infection_count::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # number infections over social dynamic links
   accommodation_dynamic_infection_count::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # number infections over dynamic on-campus accommodation links
   household_infection_count::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # number infections occuring within household
   num_init_infected::Array{Array{Int64,1}} = Array{Array{Int64,1}}(undef,countfinal) # number of infections caused by the initially infected nodes
   Rt::Array{Float64,2} = zeros(Float64,endtime+1,countfinal) # real time R value (number of secondary infections caused by nodes that were newly infected on that day)

   # 2D outputs. Residence location based prevalence. On-campus & off-campus
   prevlat_oncampus::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # Prevalence for latently infected
   prevsymp_oncampus::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # Prevalence for symptomatic infectious (post pre-symptomatic phase)
   prevasymp_oncampus::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # Prevalence for asymptomatic infectious
   prevpresymp_oncampus::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # Prevalence of pre-symptomatic symptomatic infectious
   prevrec_oncampus::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # Prevalence for recovereds

   prevlat_offcampus::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # Prevalence for latently infected
   prevsymp_offcampus::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # Prevalence for symptomatic infectious (post pre-symptomatic phase)
   prevasymp_offcampus::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # Prevalence for asymptomatic infectious
   prevpresymp_offcampus::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # Prevalence of pre-symptomatic symptomatic infectious
   prevrec_offcampus::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # Prevalence for recovereds

   # 2D outputs. Isolation related
   num_isolating::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # number of individuals isolating at each timepoint
   num_isolating_oncampus::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # number of oncampus resident individuals isolating at each timepoint
   num_isolating_offcampus::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # number of offcampus resident individuals isolating at each timepoint
   num_household_isolating::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # number of individuals isolating due to household members having symptoms at each timepoint
   num_symp_isolating::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # number of individuals isolating due to symptoms at each timepoint
   num_asymp_isolating::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # number of individuals isolating due to symptoms at each timepoint
   num_isolating_CTcause::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # number of individuals isolating as a contact of a positive test
   num_accom_lockdown_isolating::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # number of individuals isolating due to accommodation lockdown

   # 2D outputs. Rehousing related
   new_rehoused::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # new instances of students being rehoused
   current_rehoused::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # number of students in rehoused state each timestep

   #2D outputs. Testing related
   tests_performed::Array{Int64,2} = zeros(Int64,endtime+1,countfinal)

   #3D outputs. Testing related
    # Slices: 1 - True positive; 2 - false negative; 3 - True negative; 4 - false positive.
   test_outcomes::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,4)

   # 3D outputs (infection counts occuring in household, study, society setting from non-dynamic links)
   cohort_infection_count::Array{Int64,3} = zeros(Int64,0,0,0)
   society_infection_count::Array{Int64,3} = zeros(Int64,0,0,0)

   # 1D outputs - Infection counts
   n_oncampus_inf::Array{Int64,1} = zeros(Float64,countfinal) # Per replicate, number of oncampus residents infected individuals that are set to adhere to isolation guidance
   n_offcampus_inf::Array{Int64,1} = zeros(Float64,countfinal) # Per replicate, number of students resident off-campus infected

   # 1D outputs - Isolation
   n_isol_adhering::Array{Int64,1} = zeros(Float64,countfinal) # Per replicate, number of individuals that are set to adhere to isolation guidance
   n_isol_adhering_oncampus::Array{Int64,1} = zeros(Float64,countfinal) # Per replicate, number of oncampus resident individuals that are set to adhere to isolation guidance
   n_isol_adhering_offcampus::Array{Int64,1} = zeros(Float64,countfinal) # Per replicate, number of offcampus resident individuals that are set to adhere to isolation guidance

   # 1D outputs  - Misc
   infected_by::Array{Int64,1} = zeros(Int64,n_students)  # who each node was infected by
   var_num_infected::Array{Int64,1} = zeros(Int64,countfinal) # variance in the number of infections caused by each node
   mean_init_generation_time::Array{Float64,1} = zeros(Float64,countfinal) # mean initial generation time
end

"""
    Structure: student_states

Parameter structure with variables tracking status of each student.

Location: parametertypes.jl
"""
@with_kw mutable struct student_states
# used for the status of each node
   n_students::Int64 = 0
   timelat::Array{Int64,1} = zeros(Int64,n_students) # time currently spent in latent state
   timeinf::Array{Int64,1} = zeros(Int64,n_students) # time currently spent in infectious state
   timesymp::Array{Int64,1} = zeros(Int64,n_students) # time currently spent in symptomatic state
   lattime::Array{Int64,1} = zeros(Int64,n_students) # time to spend in latent state
   inftime::Int64 = 4 # time to spend in infectious state
   symptime::Int64 = 10 # time to spend in symptomatic state
   asymp::Array{Int64,1} = zeros(Int64,n_students) # whether the node will be asymptomatic
   timeisol::Array{Int64,1} = zeros(Int64,n_students) # time currently spent in isolation due to housemate symptoms
   symp_timeisol::Array{Int64,1} = zeros(Int64,n_students) # time currently spent isolating due to symptoms
   asymp_timeisol::Array{Int64,1} = zeros(Int64,n_students) # time currently spent isolating due to positive result in mass testing (finding asymptomatic infection)
   timeisol_CTcause::Array{Int64,1} = zeros(Int64,n_students) # time currently spent isolating due to being contact traced
   rep_inf_this_timestep::Array{Int64,1} = zeros(Int64,n_students) # whether the node reports symptoms this timestep
   inclass::Array{Int64,2} = zeros(Int64,n_students,endtime) # whether the student is in class
   hh_isolation::Array{Int64,1} = zeros(Int64,n_students) # Whether individual adheres to isolation guidance. (1) Yes. (0) No.
   delay_adherence::Array{Int64,1} = zeros(Int64,n_students) # Individual may not report symptoms immediately.
   acquired_infection::Array{Int64,1} = zeros(Int64,n_students) # time node acquired infection
end

"""
    Structure: intervention_data_feeds

Parameter structure with data inputs used for interventions.

Location: parametertypes.jl
"""
@with_kw mutable struct intervention_data_feeds
   rep_inf_this_timestep::Array{Int64,1} = Array{Float64,1}[] #Indicator of whether a node reports symptomatic infection during current timestep
   output::sim_outputs = sim_outputs()
   network_params::network_params = network_params()
   time::Int64 = 0
   output_time_idx::Int64 = 0
   replicate_id::Int64 = 0
end

"""
    Structure: intervention_params

Parameter structure to store intervention attributes.

Location: parametertypes.jl
"""
@with_kw mutable struct intervention_params
   time_horizon::Int64 = 0 # In determining whether to invoke an intervention, length of time to look back over to make decision
   inactivation_length::Int64 = 0 # Time for an activity to be inactive before checking if it can start again
   absolute_case_threshold::Int64 = 0  # Absolute number of group that needs to satisfy condition
   propn_case_threshold::Float64 = 0.  # Proportion of group that needs to satisfy condition
   release_condition_time_horizon::Int64 = 0   # Time period to look over when determining if conditions should be removed
   release_condition_absolute_cases::Int64 = 0 # For removing measures, level cases need to drop too
   release_condition_propn_cases::Float64 = 0. # For removing measures, condition related to propn of group that must be satisfied
end

"""
    Structure: mass_testing_params

Parameter structure to store mass testing intervention attributes.

Location: parametertypes.jl
"""
@with_kw mutable struct mass_testing_params
   designated_test_times::Array{Int64,1} = [0] # Timesteps on which mass testing will take place. Entry per mass testing event.
   on_campus_coverage_propn::Array{Float64,1} = [1.] # On campus residents, coverage.
   off_campus_coverage_propn::Array{Float64,1} = [1.] # Off campus residents, coverage.
   asymp_test_detection_prob_vec::Array{Float64,1} = [0.,0.11,0.59,0.785,0.83,0.845,0.84,0.82,0.79,0.76,  # Days 1-10
                                                            0.72,0.68,0.64,0.59,0.54,0.485,0.445,0.405,0.37,0.335, # Days 11-20
                                                            0.30,0.27,0.24,0.22,0.20,0.18,0.16,0.15,0.14,0.13] # Days 21-30
   n_mass_tests_performed::Array{Int64,1} = [0]         # Track the number of mass tests instances.  Entry per replicate.
   n_tests_performed::Array{Array{Int64,1},1} = Array{Int64,1}[]        # Track the number of tests performed in each mass testing instance. Vector per replicate.
   n_tests_positive::Array{Array{Int64,1},1} = Array{Int64,1}[]         # Track the number of positive tests in each mass testing instance. Vector per replicate.
   n_all_isolations_caused::Array{Array{Int64,1},1} = Array{Int64,1}[]  # Track the number of isolations that result from the mass testing. Vector per replicate.
   n_hh_isolations_caused::Array{Array{Int64,1},1} = Array{Int64,1}[]   # Track the number of isolations that result from the mass testing. Vector per replicate.
   n_CT_isolations_caused::Array{Array{Int64,1},1} = Array{Int64,1}[]   # Track the number of isolations that result from the mass testing. Vector per replicate.
   n_prev_infected_tested::Array{Array{Int64,1},1} = Array{Int64,1}[]   # Track the number of ppl that had been previously infected that received a test. Vector per replicate.
end
