module neutrinoDetector

    using CairoMakie

    include("neutrinoDetectorEffects/perfectDetector.jl")

    export PerfectDetector, PerfectDetectorEventTable
    export bin_edges_from_centers, bin_widths_from_centers
    export oscillation_channel_table, interacting_events
    export plot_interacting_events

end
