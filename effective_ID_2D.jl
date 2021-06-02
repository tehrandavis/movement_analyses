#== effective_ID_2D.jl:
by: Tehran J. Davis, 2021

this script takes a scalar timeseries (vector) and calculates the effective ID
based on the average amplitude of the movements and the variability of the targeting.

effective_ID_2D(timeseries, peakrange), where:

timeseries = scalar 1D timeseries, e.g., [1,2,3,4,5,6,7,6,5,4,3,2,1]
peakrange = minimum number of timepoints between peaks, e.g. [3]

assume your data series is a vector named "test", and you desire a minimum
peakrange of 20:

effective_ID_2D(test, 20)

returns a Dictionary containing
 - mean amplitude of the oscillations
 - effective widths of the two targets (more positive v. more negative peaks)
 - resulting effective index of difficulty for the wo targets
==#

using Statistics, StatsBase, Peaks




function effective_ID_2D(timeseries, peakrange)

    pos_peaks_index = Peaks.maxima(timeseries, peakrange)
    neg_peaks_index = Peaks.minima(timeseries, peakrange)

    pos_peaks_values = timeseries[pos_peaks_index]
    neg_peaks_values = timeseries[neg_peaks_index]

    # remove the first and last values (transients)
    popfirst!(pos_peaks_values)
    pop!(pos_peaks_values)

    popfirst!(neg_peaks_values)
    pop!(neg_peaks_values)

    # ensure same number of peaks for both pos and neg
    end_peaks = minimum([length(neg_peaks_values), length(pos_peaks_values)])

    # get effective p1_eff_amplitude
    eff_amplitude = mean(pos_peaks_values[1:end_peaks]-neg_peaks_values[1:end_peaks])

    #pos_peaks_values_unit_transform = StatsBase.fit(UnitRangeTransform, pos_peaks_values; dims=1, unit=true)
    #neg_peaks_values_unit_transform = StatsBase.fit(UnitRangeTransform, neg_peaks_values; dims=1, unit=true)

    # get effective widths
    eff_width_pos_targets = 4.113 * std(pos_peaks_values[1:end_peaks])
    eff_width_neg_targets = 4.113 * std(neg_peaks_values[1:end_peaks])
    eff_width_both_CC = begin
        mean_pos = mean(pos_peaks_values[1:end_peaks])
        mean_neg = mean(neg_peaks_values[1:end_peaks])
        diff_mean = mean_pos-mean_neg
        neg_transformed_values = neg_peaks_values[1:end_peaks] .+ diff_mean
        4.133 * std([pos_peaks_values[1:end_peaks] neg_transformed_values])
    end
    eff_ID_pos_targets = log2(2*eff_amplitude/eff_width_pos_targets)
    eff_ID_neg_targets = log2(2*eff_amplitude/eff_width_neg_targets)
    eff_ID_both_CC = log2(2*eff_amplitude/eff_width_both_CC)


    Dict("amplitude" => eff_amplitude,
            "eff_width_pos" => eff_width_pos_targets,
            "eff_width_neg" => eff_width_neg_targets,
            "eff_width_both" => eff_width_both_CC,
            "eff_ID_pos" => eff_ID_pos_targets,
            "eff_ID_neg" => eff_ID_neg_targets,
            "eff_ID_both" => eff_ID_both_CC
            )
end
