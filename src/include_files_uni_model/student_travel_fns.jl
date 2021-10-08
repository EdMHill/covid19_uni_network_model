#=
Purpose:
Compute number of infected students travelling from and returning to each region
=#


#=
Compute number of latents, asymptomatics & symptomatics that would be returning to each region
=#
function get_return_students_infected(prev_latent::Int64,
                                       prev_asymp::Int64,
                                       prev_symp::Int64)
#Inputs:
# prev_latent, prev_asymp, prev_symp - Number of the student population in each infected state

#Outputs:
# return_students_infected_array::Array{Float64,2} - row by region, columns: latent,asymp,symp

 # Get student popn split by region domicile data
   # 1: North West
   # 2: North East,
   # 3: Yorkshire and The Humber,
   # 4: East Midlands,
   # 5: West Midlands,
   # 6: East of England,
   # 7: London,
   # 8: South East,
   # 9: South West,
   # 10: Wales,
   # 11: Scotland,
   # 12: Northern Ireland,
   # 13: Guernsey, Jersey and the Isle of Man,
   # 14: Non-UK
 student_travel_counts = [727,155,515,1158,3450,1592,3477,2594,992,365,125,80,35,9315]
 propn_of_students_per_region = student_travel_counts./sum(student_travel_counts)

 # Initialise output array.
    # row by region, columns: latent,asymp,symp
 n_regions = length(student_travel_counts)
 return_students_infected_array = zeros(n_regions,3)

 # Take counts of latent, asymp, symp at end of term
 # Scale counts to give counts returning to each region
 return_students_infected_array[:,1] = prev_latent.*propn_of_students_per_region
 return_students_infected_array[:,2] = prev_asymp.*propn_of_students_per_region
 return_students_infected_array[:,3] = prev_symp.*propn_of_students_per_region

 return return_students_infected_array::Array{Float64,2}
end

#=
Compute number of latents, asymptomatics & symptomatics that could be arriving from each region
=#
function get_arrive_students_infected_recovered(uni_age_latent_by_region::Array{Float64,1},
                                                uni_age_asymp_by_region::Array{Float64,1},
                                                uni_age_symp_by_region::Array{Float64,1},
                                                uni_age_rec_by_region::Array{Float64,1})
#Inputs:
# uni_age_XXX_by_region - For specified disease state, estimated number of those of uni age
#                       estimated to be in that state within that region

#Outputs:
# n_initial_XXX - For specified disease state, number of arriving student popn
#                       estimated to be in that state
# propn_rec - Fraction of uni population (domiciled in UK) that have been previously infected

    # Propn of popn in each region in age range 15-29
    region_popn_in_uni_ages = [1048194,      # 'East of England',
                                        1812370,   # 'London'
                                        2044955,   # 'Midlands'
                                        1567363,   # 'North East and Yorkshire'
                                        1372471,   # 'North West'
                                        1599254,   # 'South East'
                                        980546,   # 'South West'
                                        585874,   # 'Wales'
                                        1017528,   # 'Scotland'
                                        352445]   # 'Northern Ireland'

    # Use counts of students from each region, in order:
     #  'East of England','London','Midlands','North East and Yorkshire','North West','South East','South West','Wales','Scotland','Northern Ireland'};
    student_travel_counts_model_region_order = [1592,3477,4608,670,727,2594,992,365,125,80]

    # Get fraction of uni population that makes up that age range
    propn_uni_age_region_popn_arriving = student_travel_counts_model_region_order./region_popn_in_uni_ages

    # Get number of students being in each infected state AND attending the uni
    n_initial_latent = round(Int64,sum(uni_age_latent_by_region.*propn_uni_age_region_popn_arriving))
    n_initial_asymp = round(Int64,sum(uni_age_asymp_by_region.*propn_uni_age_region_popn_arriving))
    n_initial_symp = round(Int64,sum(uni_age_symp_by_region.*propn_uni_age_region_popn_arriving))

    # Get number and propn of students being in recovered state AND attending the uni
    n_initial_rec = round(Int64,sum(uni_age_rec_by_region.*propn_uni_age_region_popn_arriving))
    propn_rec = n_initial_rec/sum(student_travel_counts_model_region_order)

    # Return outputs
    return  n_initial_latent::Int64,
            n_initial_asymp::Int64,
            n_initial_symp::Int64,
            propn_rec::Float64
end
