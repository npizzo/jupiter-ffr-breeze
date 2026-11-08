using Breeze, Oceananigans, Test, TOML
using Breeze.Thermodynamics: MoistureMassFractions, StaticEnergyState, mixture_heat_capacity,
    mixture_gas_constant, density, temperature, absolute_zero_latent_heat
using Breeze.Microphysics: adjust_thermodynamic_state

function validate_jovian_moist_thermodynamics(out)
    mkpath(out)
    c=ThermodynamicConstants(Float64;gravitational_acceleration=25.0,dry_air_molar_mass=0.0023,dry_air_heat_capacity=13000.0)
    Rd=c.molar_gas_constant/c.dry_air.molar_mass;Rv=c.molar_gas_constant/c.vapor.molar_mass
    L0=absolute_zero_latent_heat(c,c.liquid);I0=absolute_zero_latent_heat(c,c.ice)
    eq=MixedPhaseEquilibrium(Float64)
    micro=SaturationAdjustment(Float64;equilibrium=eq,solver=SecantSolver(Float64;abstol=1e-8,maxiter=80))
    rows=Any[]
    @testset "Water in hydrogen-helium: EOS, loading, saturation, latent heat" begin
        @test Rv<Rd
        dry=MoistureMassFractions(0.0);vapor=MoistureMassFractions(0.01);loaded=MoistureMassFractions(0.0,0.01)
        r0=density(300.0,5e5,dry,c);rv=density(300.0,5e5,vapor,c);rl=density(300.0,5e5,loaded,c)
        @test r0<rv<rl
        @test rv/r0 ≈ Rd/(0.99Rd+0.01Rv) rtol=1e-12
        # Independent triple-point saturation relation including the heavy molecular weight.
        qs=Breeze.adjustment_saturation_specific_humidity(273.16,5e5,0.01,c,WarmPhaseEquilibrium())
        @test qs ≈ 0.99*(Rd/Rv)*611.657/(5e5-611.657) rtol=1e-10
        for pressure in (1e5,5e5,1e6),T0 in (220.0,250.0,280.0,320.0),qt in (0.001,0.01,0.05)
            q0=MoistureMassFractions(qt)
            energy=mixture_heat_capacity(q0,c)*T0
            state=StaticEnergyState(energy,q0,0.0,pressure)
            adjusted=adjust_thermodynamic_state(state,micro,c)
            T1=temperature(adjusted,c);q=adjusted.moisture_mass_fractions
            reconstructed=mixture_heat_capacity(q,c)*T1-L0*q.liquid-I0*q.ice
            sat=Breeze.adjustment_saturation_specific_humidity(T1,pressure,qt,c,eq)
            @test all(isfinite,(T1,q.vapor,q.liquid,q.ice))
            @test minimum((q.vapor,q.liquid,q.ice)) >= -1e-12
            @test q.vapor+q.liquid+q.ice ≈ qt rtol=1e-11
            @test reconstructed ≈ energy rtol=1e-10
            @test q.vapor ≈ min(qt,sat) rtol=1e-5 atol=1e-10
            @test T1 >= T0-1e-6
            push!(rows,(pressure,T0,qt,T1,q.vapor,q.liquid,q.ice,(reconstructed-energy)/energy))
        end
    end
    open(joinpath(out,"equilibrium.csv"),"w") do io
        println(io,"pressure_Pa,initial_T_K,total_water,adjusted_T_K,vapor,liquid,ice,relative_energy_residual")
        for row in rows;println(io,join(row,','));end
    end
    # A prescribed hydrostatic pressure path, not a dynamically accelerating parcel.
    open(joinpath(out,"prescribed_ascent.csv"),"w") do io
        println(io,"z_m,pressure_Pa,total_water,T_K,vapor,liquid,ice,T_without_phase_change_K")
        @testset "Prescribed ascent: dry lapse and reversible phase partition" begin
            for qt in (0.0,0.02)
                q0=MoistureMassFractions(qt);cp0=mixture_heat_capacity(q0,c);energy=cp0*300.0
                lastcloud=0.0
                for z in 0.0:250.0:20000.0
                    pressure=5e5*(1-25z/(13000*300))^(13000/Rd)
                    a=adjust_thermodynamic_state(StaticEnergyState(energy,q0,z,pressure),micro,c)
                    T=temperature(a,c);q=a.moisture_mass_fractions
                    Tdry=300-25z/cp0
                    @test T>=Tdry-1e-6
                    @test q.vapor+q.liquid+q.ice ≈ qt atol=1e-12
                    qt==0 && (@test T≈Tdry atol=1e-8)
                    lastcloud=q.liquid+q.ice
                    println(io,join((z,pressure,qt,T,q.vapor,q.liquid,q.ice,Tdry),','))
                end
                qt>0 && (@test lastcloud>1e-5)
            end
        end
    end
    open(joinpath(out,"PASS.toml"),"w") do io
        TOML.print(io,Dict("status"=>"PASS","equilibrium_cases"=>length(rows),"carrier_molar_mass"=>0.0023,
            "carrier_heat_capacity"=>13000.0,"scope"=>"Thermodynamic implementation checks; not Jovian microphysics calibration"))
    end
end
