abstract type AdaptiveFlagsMarkingStrategy end;


# Implements Algorithm for calculating coarsening and refinement thresholds
# in Fig 5. of the following paper:
# Algorithms and data structures for massively parallel generic adaptive finite element codes
# W Bangerth, C Burstedde, T Heister, M Kronbichler
# ACM Transactions on Mathematical Software 38 (2)
struct FixedFractionAdaptiveFlagsMarkingStrategy{T<:Real} <: AdaptiveFlagsMarkingStrategy
  # Refine all cells s.t. #{e_i > \theta_r} \approx this%refinement_fraction * N_cells
  # Coarsen all cells s.t. #{e_i < \theta_c} \approx this%coarsening_fraction * N_cells 
  refinement_fraction::T
  coarsening_fraction::T
end

function compute_target_num_cells(num_global_cells,fraction)
  map(num_global_cells) do num_global_cells
    FT=eltype(fraction)
    NGCT=eltype(num_global_cells)
    round(NGCT, FT(num_global_cells)*fraction)
  end
end 

function _compute_thresholds(error_indicators, 
                             cell_partition,
                             refinement_fraction,
                             coarsening_fraction; 
                             verbose=false)

  num_owned_cells = map(cell_partition) do partition
    own_length(partition) 
  end 
  num_global_cells  = map(cell_partition) do partition
    global_length(partition)
  end

  target_num_cells_to_be_refined = 
       compute_target_num_cells(num_global_cells, 
                                refinement_fraction)

  target_num_cells_to_be_coarsened = 
       compute_target_num_cells(num_global_cells, 
                                coarsening_fraction)


  sq_error_indicators = map(error_indicators) do error_indicator
    error_indicator.*error_indicator
  end 

  ref_sq_min_estimate = map(sq_error_indicators, num_owned_cells) do sq_error_indicator, num_owned_cells
    minimum(view(sq_error_indicator,1:num_owned_cells))
  end

  ref_sq_max_estimate = map(sq_error_indicators, num_owned_cells) do sq_error_indicator, num_owned_cells
    maximum(view(sq_error_indicator,1:num_owned_cells))
  end

  ref_min_estimate=reduction(min,ref_sq_min_estimate,init=typemax(eltype(ref_sq_min_estimate)))
  ref_max_estimate=reduction(max,ref_sq_max_estimate,init=typemin(eltype(ref_sq_max_estimate)))

  # We compute refinement thresholds by bisection of the interval spanned by
  # the smallest and largest error indicator. this leads to a small problem:
  # if, for example, we want to refine zero per cent of the cells, then we
  # need to pick a threshold equal to the largest indicator, but of course
  # the bisection algorithm can never find a threshold equal to one of the
  # end points of the interval. So we slightly increase the interval before
  # we even start
  ref_min_estimate=map(ref_min_estimate) do ref_min_estimate
    if (ref_min_estimate>0.0)
      return ref_min_estimate*0.99
    else 
      return ref_min_estimate
    end  
  end

  ref_max_estimate=map(ref_max_estimate) do ref_max_estimate
    if (ref_max_estimate>0.0)
      return ref_max_estimate*1.01
    else 
      return ref_max_estimate
    end
  end

  coarsening_min_estimate = ref_min_estimate
  coarsening_max_estimate = ref_max_estimate
  num_iterations = 0
  refinement_converged = false 
  coarsening_converged = false 
  ref_split_estimate = nothing 
  coarsening_split_estimate = nothing
  current_num_cells_to_be_coarsened = nothing 
  current_num_cells_to_be_refined = nothing

  while (true)

    refinement_converged_old=refinement_converged

    map(ref_min_estimate,ref_max_estimate) do ref_min_estimate, ref_max_estimate 
      if (ref_min_estimate==ref_max_estimate) 
        refinement_converged = true
      end
    end

    if (verbose)
      if (!refinement_converged_old && refinement_converged)
        map_main(num_global_cells,current_num_cells_to_be_refined) do num_global_cells,current_num_cells_to_be_refined
          println("FixedFractionAdaptiveMarkingStrategy: $(current_num_cells_to_be_refined)/$(num_global_cells) to be refined")
        end
      end
    end 

    coarsening_converged_old=coarsening_converged
    map(coarsening_min_estimate,coarsening_max_estimate) do coarsening_min_estimate, coarsening_max_estimate 
      if (coarsening_min_estimate==coarsening_max_estimate) 
        coarsening_converged = true
      end
    end


    if (verbose)
      if (!coarsening_converged_old && coarsening_converged)
        map_main(num_global_cells,current_num_cells_to_be_coarsened) do num_global_cells,current_num_cells_to_be_coarsened
          println("FixedFractionAdaptiveMarkingStrategy: $(current_num_cells_to_be_coarsened)/$(num_global_cells) to be coarsened")
        end
      end
    end 

    if (refinement_converged && coarsening_converged)
      break
    end

    if (!refinement_converged)
      # Compute interval split point using the fact that the log of error estimators
      # is much better uniformly scattered than the error estimators themselves. This is required
      # in order to have faster convergence whenever the error estimators are scattered across very
      # different orders of magnitude
      # avg_estimate = exp(1/2*(log(min_estimate)+log(max_estimate))) = sqrt(min_estimate*max_estimate)
      ref_split_estimate=map(ref_min_estimate,ref_max_estimate) do ref_min_estimate, ref_max_estimate
        if (ref_min_estimate==0.0)
          ref_split_estimate = sqrt(1.0e-10*ref_max_estimate)
        else 
          ref_split_estimate = sqrt(ref_min_estimate*ref_max_estimate)
        end
        ref_split_estimate
      end
    end 
  
    if (!coarsening_converged)
      # Compute interval split point using the fact that the log of error estimators
      # is much better uniformly scattered than the error estimators themselves. This is required
      # in order to have faster convergence whenever the error estimators are scattered across very
      # different orders of magnitude
      # avg_estimate = exp(1/2*(log(min_estimate)+log(max_estimate))) = sqrt(min_estimate*max_estimate)
      coarsening_split_estimate=map(coarsening_min_estimate, coarsening_max_estimate) do coarsening_min_estimate,
                                                                                         coarsening_max_estimate
        if (coarsening_min_estimate==0.0)
          coarsening_split_estimate = sqrt(1.0e-10*coarsening_max_estimate)
        else 
          coarsening_split_estimate = sqrt(coarsening_min_estimate*coarsening_max_estimate)
        end
        coarsening_split_estimate
      end
    end

    # # Count how many cells have local error estimate larger or equal to avg_estimate
    # # count = #{ i: e_i >= avg_estimate }
    current_num_cells_to_be_refined,
       current_num_cells_to_be_coarsened=map(sq_error_indicators,
                                             ref_split_estimate,
                                             coarsening_split_estimate,
                                             num_owned_cells) do sq_error_indicators, 
                                                                 ref_split_estimate,
                                                                 coarsening_split_estimate, 
                                                                 num_owned_cells

      current_num_cells_to_be_refined = 0
      if (!refinement_converged)  
        for i=1:num_owned_cells
          if (sq_error_indicators[i]>=ref_split_estimate*ref_split_estimate)
            current_num_cells_to_be_refined += 1
          end
        end
      end

      current_num_cells_to_be_coarsened = 0
      if (!coarsening_converged)
        for i=1:num_owned_cells
          if (sq_error_indicators[i]<coarsening_split_estimate*coarsening_split_estimate)
            current_num_cells_to_be_coarsened += 1
          end
        end 
      end
      current_num_cells_to_be_refined, current_num_cells_to_be_coarsened
    end |> tuple_of_arrays

    if (!refinement_converged)  
       current_num_cells_to_be_refined=
           reduction(+,current_num_cells_to_be_refined,init=zero(eltype(current_num_cells_to_be_refined)))
       ref_min_estimate,ref_max_estimate=map(current_num_cells_to_be_refined, 
                                             target_num_cells_to_be_refined,
                                             ref_min_estimate,
                                             ref_max_estimate,
                                             ref_split_estimate) do current_num_cells_to_be_refined,
                                                                  target_num_cells_to_be_refined, 
                                                                  ref_min_estimate,
                                                                  ref_max_estimate, 
                                                                  ref_split_estimate
          if (current_num_cells_to_be_refined>target_num_cells_to_be_refined)
            ref_min_estimate = ref_split_estimate
          elseif (current_num_cells_to_be_refined<target_num_cells_to_be_refined)
            ref_max_estimate = ref_split_estimate
          else
            ref_min_estimate = ref_split_estimate
            ref_max_estimate = ref_split_estimate
          end
          ref_min_estimate,ref_max_estimate
        end |> tuple_of_arrays
    end

    if (!coarsening_converged)
        current_num_cells_to_be_coarsened=
            reduction(+,current_num_cells_to_be_coarsened,init=zero(eltype(current_num_cells_to_be_coarsened)))
        coarsening_min_estimate,coarsening_max_estimate = map(current_num_cells_to_be_coarsened,
                                                              target_num_cells_to_be_coarsened,         coarsening_min_estimate,
                                                              coarsening_max_estimate,
                                                              coarsening_split_estimate) do               current_num_cells_to_be_coarsened,
                                                                                  target_num_cells_to_be_coarsened,
                                                                                  coarsening_min_estimate,
                                                                                  coarsening_max_estimate,
                                                                                  coarsening_split_estimate
            if (current_num_cells_to_be_coarsened>target_num_cells_to_be_coarsened)
              coarsening_max_estimate = coarsening_split_estimate
            elseif (current_num_cells_to_be_coarsened<target_num_cells_to_be_coarsened)
              coarsening_min_estimate = coarsening_split_estimate
            else
              coarsening_min_estimate = coarsening_split_estimate
              coarsening_max_estimate = coarsening_split_estimate
            end
            coarsening_min_estimate,coarsening_max_estimate
          end |> tuple_of_arrays
    end    
    num_iterations += 1
    if ( num_iterations == 25 )
        ref_min_estimate  = ref_split_estimate
        ref_max_estimate  = ref_split_estimate
        coarsening_max_estimate = coarsening_split_estimate
        coarsening_min_estimate = coarsening_split_estimate
    end
  end


  θr2 = nothing
  map(ref_min_estimate) do ref_min_estimate 
    θr2=ref_min_estimate*ref_min_estimate
  end 
  
  θc2 = nothing
  map(coarsening_min_estimate) do coarsening_min_estimate
    θc2=coarsening_min_estimate*coarsening_min_estimate
  end

  if ( verbose )
    map_main(num_global_cells,current_num_cells_to_be_refined,current_num_cells_to_be_coarsened) do num_global_cells,current_num_cells_to_be_refined,current_num_cells_to_be_coarsened
      if (!refinement_converged)
        println("FixedFractionAdaptiveMarkingStrategy: $(current_num_cells_to_be_refined)/$(num_global_cells) to be refined")
      end 
      if (!coarsening_converged)
        println("FixedFractionAdaptiveMarkingStrategy: $(current_num_cells_to_be_coarsened)/$(num_global_cells) to be coarsened")
      end
      println("FixedFractionAdaptiveMarkingStrategy: Refinement threshold: $(sqrt(θr2))")
      println("FixedFractionAdaptiveMarkingStrategy: Coarsening threshold: $(sqrt(θc2))")
    end
  end
  return sqrt(θr2), sqrt(θc2)
end


function update_adaptivity_flags!(flags::AbstractArray{<:AbstractVector{<:Integer}},
                                  strategy::FixedFractionAdaptiveFlagsMarkingStrategy, 
                                  cell_partition,
                                  error_indicators::AbstractArray{<:AbstractVector{<:Real}};
                                  verbose=false)

  θr,θc=_compute_thresholds(error_indicators, 
                            cell_partition,
                            strategy.refinement_fraction,
                            strategy.coarsening_fraction; 
                            verbose=verbose)
  
  map(flags,cell_partition,error_indicators) do flags, partition, error_indicators
    flags .= nothing_flag
    for i=1:own_length(partition)
      if (error_indicators[i]>=θr)
        flags[i]=refine_flag
      elseif (error_indicators[i]<θc)
        flags[i]=coarsen_flag
      end
    end
  end
end







