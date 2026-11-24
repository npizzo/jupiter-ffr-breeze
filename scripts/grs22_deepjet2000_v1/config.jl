function base_config(equations,jet)
    merge(TOML.parsefile(joinpath(@__DIR__,"pilot.toml")),TOML.parsefile(joinpath(@__DIR__,"candidate.toml")),Dict(
    "equations"=>equations,"latitude"=>-22.,"beta"=>6cosd(22.)*250/71492.,"length_unit_km"=>250.,"nx"=>128,"ny"=>256,"nz"=>80,"lx"=>4.,"ly"=>8.,"lz"=>1.,"initial_rms"=>.136,"perturbation"=>0.,
    "jet_speed"=>(jet ? -.125 : 0.),"jet_design"=>"deep","eddy_tilt_correction"=>true,"neutral_coriolis"=>true,
    "heat"=>1.,"bottom_factor"=>1.,"conduction_factor"=>1.,"top_factor"=>1.,"cap_factor"=>1.,"mean_restoration"=>true,
    "conductivity"=>.0004*.03125/(GRAV/2.5),"entropy_diffusion_factor"=>0.,"closure"=>"none","advection_order"=>9,
    "sponge"=>false,"horizontal_sponge"=>false,"initial_dt"=>.001,"max_dt"=>1/120,"stop_time"=>2000.,"wall_seconds"=>26000.,
    "diagnostic_interval"=>2.,"movie_interval"=>2.,"checkpoint_interval"=>20.,"field_interval"=>100.,"save_3d"=>true,
    "mach_stop"=>.8,"outside_mach_stop"=>Inf,"diagnostic_extension"=>true,"temperature_stop"=>0.,"profile_drift_stop"=>.5,"theta_log_departure_stop"=>log(4.),
    "purpose"=>"Single fresh compressible 22S beta channel with full-depth anticyclonic shear; qualitative GRS-oriented formation/persistence test, not dimensional reproduction"))
end
