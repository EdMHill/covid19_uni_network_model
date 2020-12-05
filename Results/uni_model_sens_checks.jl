"""
Purpose:
Plot outputs from university model sensitivity runs
"""

#Set paths
cd(dirname(@__FILE__))

#Load environment
using Pkg
Pkg.activate("../")

using Plots, MAT, LinearAlgebra, Statistics

# List of sensitivity configs to plot results for
svty_variables = [# "RNGseed",
                 # "mawidth",
                 # "trans_risk_no_interv",
                 # "probasymp_scale_no_interv",
                 # "transasymp_scale_no_interv",
                 # "suscep_scale_no_interv",
                  "No control simulations",
                 ]

# Values used in sensitivity configs
svty_variable_ops = Dict("RNGseed" => [100],  #[100, 200, 300],
                            "mawidth" => [1,3,5,7,14,28],
                            "trans_risk_no_interv" => [0.2,0.3,0.4,0.5],
                            "probasymp_scale_no_interv" => [0.5,0.6,0.7,0.8,0.9],
                            "transasymp_scale_no_interv" => [0.1,0.2,0.3,0.4,0.5,0.6],
                            "No control simulations" => [1000]
                            )

"""
Iterate over sensitivity configs
"""
# In turn, produce plots for each batch of runs
for svty_itr = 1:length(svty_variables)

    # Load configuration data for current sensitivity config
    variable_name = svty_variables[svty_itr]
    variable_ops = svty_variable_ops[variable_name]
    n_vars = length(variable_ops)

    # Load sensitivity output file
    if (variable_name == "RNGseed") || (variable_name == "mawidth")
        #file=matopen("../uni_model_infection_output_run_no_interventions.mat","r")
        file=matopen("../uni_model_infection_output_run_one_run_baseline_nointerv_#1.mat","r")
    elseif (variable_name == "trans_risk_no_interv")
        file=matopen("../uni_model_infection_output_trans_risk_scale_no_interv_#1.mat","r")
    elseif (variable_name == "probasymp_scale_no_interv")
        file=matopen("../uni_model_infection_output_probasymp_scale_no_interv_#1.mat","r")
    elseif (variable_name == "transasymp_scale_no_interv")
        file=matopen("../uni_model_infection_output_transasymp_scale_no_interv_#1.mat","r")
    elseif (variable_name =="suscep_scale_no_interv")
        file=matopen("../uni_model_infection_output_suscep_scale_no_interv_#1.mat","r")
    elseif (variable_name == "No control simulations")
        file=matopen("../uni_model_infection_output_run_one_run_baseline_nointerv_combined.mat","r")
    end

    if (variable_name == "No control simulations")
        Rt_unfiltered = read(file,"Rt_save_combined")
        GT_init = read(file,"mean_init_generation_time_save_combined")
        Rt_init = read(file, "num_init_infected_save_combined")
        numinf = read(file,"numinf_combined")
        secondary_infs = read(file,"num_infected_combined")

        # Get transmission setting count variables
        household_infection_count = read(file,"household_infection_count_combined")
        cohort_infection_count = read(file,"cohort_infection_count_combined")
        society_infection_count = read(file,"society_infection_count_combined")
        accom_dynamic_infection_count = read(file,"accom_dynamic_infection_count_combined")
        social_dynamic_infection_count = read(file,"social_dynamic_infection_count_combined")
    else
        Rt_unfiltered = read(file,"Rt_save")
        GT_init = read(file,"mean_init_generation_time_save")
        Rt_init = read(file, "num_init_infected_save")
        numinf = read(file,"numinf")
        secondary_infs = read(file,"num_infected")

        # Get transmission setting count variables
        household_infection_count = read(file,"household_infection_count")
        cohort_infection_count = read(file,"cohort_infection_count")
        society_infection_count = read(file,"society_infection_count")
        accom_dynamic_infection_count = read(file,"accom_dynamic_infection_count")
        social_dynamic_infection_count = read(file,"social_dynamic_infection_count")
    end



    # Close file
    close(file)

    """
    If required, get transmission setting data into combined array
    Setting order: Household, accom (dynamic), cohort, society, social (dynamic)
    """
    n_timesteps::Int64 = size(household_infection_count,1)
    if variable_name!="mawidth"
        # Initialise combined array
        n_replicates::Int64 = size(household_infection_count,2)
        n_settings = 5
        transmission_setting = zeros(Int64,n_timesteps,n_replicates,n_settings,n_vars)

        # Populate slices direct from loaded variables where possible
        transmission_setting[:,:,1,:] = household_infection_count
        transmission_setting[:,:,2,:] = accom_dynamic_infection_count
        transmission_setting[:,:,5,:] = social_dynamic_infection_count

        # Cohorts, need to sum across cohorts first
        total_cohort_infection_count = sum(cohort_infection_count,dims=3)
        transmission_setting[:,:,3,:] = total_cohort_infection_count

        # Cohorts, need to sum across cohorts first
        total_society_infection_count = sum(society_infection_count,dims=3)
        transmission_setting[:,:,4,:] = total_society_infection_count
    end


    # Have cutoff the simulation, so not all infections have been realised.
    # Consider values up to final timepoint less infectious period (so full infectious period has been seen)
    Rt_final_timepoint = n_timesteps

    """
    Construct plots for each configuration
    """
    for var_itr = 1:n_vars

        # Set default moving average width for Rt
        # Sets number of points either side of current timestep used to calculate moving average
        # e.g. ma_width = 3 would give 7 day moving average (3 timepoints before, current timepoint, 3 timepoints after)
        if variable_name=="mawidth"
            ma_width = variable_ops[var_itr]
        else
            ma_width = 7
        end

        # Set default population size
        if variable_name=="popsize"
            cmax = variable_ops[var_itr]
        else
            cmax = 25000
        end

        # Calculate MA for Rt
        if ma_width%2 == 1
            ma_length = Rt_final_timepoint - (ma_width - 1)
        else
            ma_length = Rt_final_timepoint - ma_width
        end

        # Set up ma_width based on if odd or even
        if ma_width%2 == 1
            # Odd
            ma_data_pts_lb_side = (ma_width - 1) ÷ 2
            ma_data_pts_ub_side = (ma_width - 1) ÷ 2
            ma_data_pts_plot_offset_lb =  (ma_width - 1) ÷ 2
            ma_data_pts_plot_offset_ub =  (ma_width - 1) ÷ 2
        else
            # Even
            ma_data_pts_lb_side = (ma_width÷2) - 1
            ma_data_pts_ub_side = (ma_width÷2)
            ma_data_pts_plot_offset_lb = (ma_width÷2)
            ma_data_pts_plot_offset_ub = (ma_width÷2)
        end

        # Get number of simn replicates in use
        if variable_name=="mawidth"
            countfinal = length(Rt_unfiltered[1,:,1,1])
        else
            countfinal = length(Rt_unfiltered[1,:,1,var_itr])
        end

        # Calculate moving average Rt
        Rt_ma = zeros(Float64, countfinal, ma_length)

        if variable_name=="mawidth"
            # Use first slice of data, var_itr = 1
            if ma_width%2 == 1
                # Central moving average using odd number of data points
                Rt_ma[1,:] = [mean(filter(!isnan,Rt_unfiltered[(lb-ma_data_pts_lb_side):(lb+ma_data_pts_ub_side),1,1])) for lb=(ma_data_pts_lb_side+1):(Rt_final_timepoint-ma_data_pts_ub_side)]
            else
                # Central moving average using even number of data points
                Rt_ma_temp = zeros(Float64, ma_length+1)
                Rt_ma_temp = [mean(filter(!isnan,Rt_unfiltered[(lb-ma_data_pts_lb_side):(lb+ma_data_pts_ub_side),1,1])) for lb=(ma_data_pts_lb_side+1):(Rt_final_timepoint-ma_data_pts_ub_side)]

                # Smooth the smoothed data points
                # E.g ma_width = 4. For timepoints 1-4, centralised average at 2.5.
                # For timepoints 2-5, centralised average at 3.5.
                # Smoothing these initial averages gives estimate at timepoint 3 (using ma_width = 4)
                for ma_itr = 1:ma_length
                    Rt_ma[1,ma_itr] = (Rt_ma_temp[ma_itr] + Rt_ma_temp[ma_itr+1])/2
                end
            end
        else
            if ma_width%2 == 1
                # Central moving average using odd number of data points
                Rt_ma[1,:] = [mean(filter(!isnan,Rt_unfiltered[(lb-ma_data_pts_lb_side):(lb+ma_data_pts_ub_side),1,1,var_itr])) for lb=(ma_data_pts_lb_side+1):(Rt_final_timepoint-ma_data_pts_ub_side)]
            else
                # Central moving average using even number of data points
                Rt_ma_temp = zeros(Float64, ma_length+1)
                Rt_ma_temp = [mean(filter(!isnan,Rt_unfiltered[(lb-ma_data_pts_lb_side):(lb+ma_data_pts_ub_side),1,1,var_itr])) for lb=(ma_data_pts_lb_side+1):(Rt_final_timepoint-ma_data_pts_ub_side)]

                # Smooth the smoothed data points
                # E.g ma_width = 4. For timepoints 1-4, centralised average at 2.5.
                # For timepoints 2-5, centralised average at 3.5.
                # Smoothing these initial averages gives estimate at timepoint 3 (using ma_width = 4)
                for ma_itr = 1:ma_length
                    Rt_ma[1,ma_itr] = (Rt_ma_temp[ma_itr] + Rt_ma_temp[ma_itr+1])/2
                end
            end
        end

        for count=2:countfinal
            if variable_name=="mawidth"
                if ma_width%2 == 1
                    # Central moving average using odd number of data points
                    Rt_ma[count,:] = [mean(filter(!isnan,Rt_unfiltered[(lb-ma_data_pts_lb_side):(lb+ma_data_pts_ub_side),count,1,1])) for lb=(ma_data_pts_lb_side+1):(Rt_final_timepoint-ma_data_pts_ub_side)]
                else
                    # Central moving average using even number of data points
                    Rt_ma_temp = zeros(Float64, ma_length+1)
                    Rt_ma_temp = [mean(filter(!isnan,Rt_unfiltered[(lb-ma_data_pts_lb_side):(lb+ma_data_pts_ub_side),count,1,1])) for lb=(ma_data_pts_lb_side+1):(Rt_final_timepoint-ma_data_pts_ub_side)]

                    # Smooth the smoothed data points
                    # E.g ma_width = 4. For timepoints 1-4, centralised average at 2.5.
                    # For timepoints 2-5, centralised average at 3.5.
                    # Smoothing these initial averages gives estimate at timepoint 3 (using ma_width = 4)
                    for ma_itr = 1:ma_length
                        Rt_ma[count,ma_itr] = (Rt_ma_temp[ma_itr] + Rt_ma_temp[ma_itr+1])/2
                    end
                end
            else
                if ma_width%2 == 1
                    # Central moving average using odd number of data points
                    Rt_ma[count,:] = [mean(filter(!isnan,Rt_unfiltered[(lb-ma_data_pts_lb_side):(lb+ma_data_pts_ub_side),count,1,var_itr])) for lb=(ma_data_pts_lb_side+1):(Rt_final_timepoint-ma_data_pts_ub_side)]
                    #plot!(p1,(ma_width+1):(n_timesteps-ma_width),Rt_ma[count,:],legend=false)
                else
                    # Central moving average using even number of data points
                    Rt_ma_temp = zeros(Float64, ma_length+1)
                    Rt_ma_temp = [mean(filter(!isnan,Rt_unfiltered[(lb-ma_data_pts_lb_side):(lb+ma_data_pts_ub_side),count,1,var_itr])) for lb=(ma_data_pts_lb_side+1):(Rt_final_timepoint-ma_data_pts_ub_side)]

                    # Smooth the smoothed data points
                    # E.g ma_width = 4. For timepoints 1-4, centralised average at 2.5.
                    # For timepoints 2-5, centralised average at 3.5.
                    # Smoothing these initial averages gives estimate at timepoint 3 (using ma_width = 4)
                    for ma_itr = 1:ma_length
                        Rt_ma[count,ma_itr] = (Rt_ma_temp[ma_itr] + Rt_ma_temp[ma_itr+1])/2
                    end
                end
            end
        end

        # Get quantiles of Rt_ma at each timepoint
        quantile_vals = [0.005, 0.05, 0.25, 0.50, 0.75, 0.95, 0.995]
        n_quantile_vals = length(quantile_vals)
        quantile_Rt_ma = zeros(ma_length,n_quantile_vals)
        for quantile_itr = 1:n_quantile_vals
            for tt=1:ma_length
                if !isempty(filter!(!isnan,Rt_ma[:,tt]))
                    quantile_Rt_ma[tt,quantile_itr] = quantile(filter!(!isnan,Rt_ma[:,tt]), quantile_vals[quantile_itr])
                end
            end
        end

        # Get average of Rt MA across all runs
        Rt_ma_av = [mean(filter!(!isnan,Rt_ma[:,t])) for t=1:ma_length]

        # Plot prediction intervals
        # 99%
        p1 = plot((ma_data_pts_plot_offset_lb):(Rt_final_timepoint-ma_data_pts_plot_offset_ub-1),[quantile_Rt_ma[:,4] quantile_Rt_ma[:,4]],
                        fillrange=[quantile_Rt_ma[:,1] quantile_Rt_ma[:,7]],
                        fillalpha=0.2,
                        c=:orange,
                        label = "")

        # 90%
        plot!(p1,(ma_data_pts_plot_offset_lb):(Rt_final_timepoint-ma_data_pts_plot_offset_ub-1),[quantile_Rt_ma[:,4] quantile_Rt_ma[:,4]],
                        fillrange=[quantile_Rt_ma[:,2] quantile_Rt_ma[:,6]],
                        fillalpha=0.5,
                        c=:orange,
                        label = "")


        # 50%
        plot!(p1,(ma_data_pts_plot_offset_lb):(Rt_final_timepoint-ma_data_pts_plot_offset_ub-1),[quantile_Rt_ma[:,4] quantile_Rt_ma[:,4]],
                        fillrange=[quantile_Rt_ma[:,3] quantile_Rt_ma[:,5]],
                        fillalpha=0.8,
                        c=:orange,
                        label = "")

        # Plot median
        plot!(p1,(ma_data_pts_plot_offset_lb):(Rt_final_timepoint-ma_data_pts_plot_offset_ub-1),quantile_Rt_ma[:,4],
                                                color = "black",
                                                # legend=true,
                                                widen=false,
                                                framestyle = :box,
                                                linewidth=1.5,
                                                label="")

        # Set up dummy legend
        plot!(p1,[],[],color = "black",linewidth=1.5,label="median",legend=(0.8,0.8), # Bottom left corner of legend is placed at (x,y).
                        legendfontsize=4)
        plot!(p1,[],[],
                    fillrange=[NaN NaN],
                    fillalpha=0.8,
                    linealpha=0.4,
                    c=:orange,
                    label = "50% prediction interval")
        plot!(p1,[],[],
                    fillrange=[NaN NaN],
                    fillalpha=0.5,
                    linealpha=0.1,
                    c=:orange,
                    label = "90% prediction interval")
        plot!(p1,[],[],
                fillrange=[NaN NaN],
                fillalpha=0.2,
                linealpha=0.05,
                c=:orange,
                label = "99% prediction interval")

        # Set up plot labels
        xlims!(p1,(0.,Rt_final_timepoint-1))
        ylabel!(p1, "Rt (MA: $(ma_width) days)")
        xlabel!(p1, "Time (days)")
        if variable_name!="mawidth"
            title!(p1, "$(variable_name): $(variable_ops[var_itr])")
        end

        # Other plots unless MA svty
        if variable_name!="mawidth"

        if (variable_name != "No control simulations")
            # R0. Get mean per replicate
            Rt_init_mean = zeros(n_replicates)
            for replicate_itr = 1:n_replicates
                Rt_init_mean[replicate_itr] = mean(Rt_init[var_itr][replicate_itr])
            end
            p2 = histogram(Rt_init_mean, #mean(Rt_init[:,:,var_itr],dims=1)[1,:],
                            legend=false,
                            xlabel="R0",
                            #normed=true,
                            normalize = :probability,
                            widen=false,
                            framestyle = :box,
                            #grid = false
                            )
        end

            # Mean generation time
            p3 = histogram(GT_init[1,:,var_itr],
                            legend=false,
                            xlabel="Mean generation time",
                            normalize = :probability,
                            widen=false,
                            framestyle = :box,
                            #grid = false
                            )

            # Total proportion infected
            p4 = histogram(numinf[end,:,var_itr]./cmax,
                            legend=false,
                            xlabel="Proportion infected",
                            normalize = :probability,
                            widen=false,
                            framestyle = :box,
                            #grid = false
                            )

            # # Outbreak duration
            # outbreak_duration = [findfirst(numinf[:,i,var_itr].==numinf[end,i,var_itr]) for i=1:countfinal]
            #
            # p5 = histogram(outbreak_duration,
            #                 legend=false,
            #                 xlabel="Duration",
            #                 normed=true)

            # New infections in each transmission setting over time
            p6 = plot(0:n_timesteps-1,mean(transmission_setting[:,:,:,var_itr],dims=2)[1:n_timesteps,1,:],
                        legend=false,
                        xlabel="Time (days)",
                        title="New infections",
                        widen=false,
                        framestyle = :box,
                        normalize = :probability,
                        #grid = false
                        )

            # Total proportion of infections in each transmission setting
            transmission_setting_props = sum(transmission_setting[:,:,:,var_itr],dims=1)[1,:,:]./repeat(numinf[end,:,var_itr],outer=[1,5])
            p7 = bar(["Household","Dynamic accom","Cohort","Society","Dynamic social"],
                    mean(transmission_setting_props,dims=1)[1,:],
                    xrotation=35,
                    xtickfontsize=font(8),
                    title="Infection setting",
                    color=palette(:default),
                    legend=false,
                    widen=false,
                    framestyle = :box,
                    #grid = false
                    )

            # Set up plot space based on variable name
            if (variable_name == "No control simulations")
                layout = @layout [a ; c ; e f]
                fig = plot(p1,p4,p7,p6, layout=layout, size=(450,600))
            else
                layout = @layout [a ; b ; c  d; e f]
                fig = plot(p1,p4,p2,p3,p7,p6, layout=layout, size=(450,600))
            end

        else
            fig = plot(p1, size=(450,200), dpi = 300)
        end


        savefig(fig, "Uni_model_sensitivity_$(variable_name)$(var_itr).pdf")
    end
end
