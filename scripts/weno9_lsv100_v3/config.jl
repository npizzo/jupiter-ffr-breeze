function base_config(equations,jet)
    merge(TOML.parsefile(joinpath(@__DIR__,"pilot.toml")),TOML.parsefile(joinpath(@__DIR__,"candidate.toml")),Dict(
    "equations"=>equations,"nx"=>128,"ny"=>256,"nz"=>80,"lx"=>4.,"ly"=>8.,"lz"=>1.,"initial_rms"=>.136,"perturbation"=>0.,
    "jet_speed"=>(jet ? .015625 : 0.),"jet_design"=>"cap","eddy_tilt_correction"=>true,"neutral_coriolis"=>true,
    "heat"=>1.,"bottom_factor"=>1.,"conduction_factor"=>1.,"top_factor"=>1.,"cap_factor"=>1.,"mean_restoration"=>true,
    "conductivity"=>.0004*.03125/(GRAV/2.5),"entropy_diffusion_factor"=>0.,"closure"=>"none","advection_order"=>9,
    "sponge"=>false,"horizontal_sponge"=>false,"initial_dt"=>.001,"max_dt"=>1/120,"stop_time"=>100.,"wall_seconds"=>6600.,
    "diagnostic_interval"=>1.,"movie_interval"=>1.,"checkpoint_interval"=>20.,"field_interval"=>50.,"save_3d"=>true,
    "mach_stop"=>2.,"outside_mach_stop"=>Inf,"diagnostic_extension"=>true,"temperature_stop"=>0.,"profile_drift_stop"=>.5,"theta_log_departure_stop"=>log(4.),
    "purpose"=>"Matched compressible/anelastic beta channel, jet/nojet; WENO9 with no closure; formation/persistence of large vortices"))
end
