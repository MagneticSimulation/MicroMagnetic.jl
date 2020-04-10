mutable struct EnergyMinimization <: Driver
  gk::Array{Float64, 1}
  tau::Float64
  max_tau::Float64
  min_tau::Float64
  steps::Int64
end

mutable struct LLG <: Driver
  precession::Bool
  alpha::Float64
  gamma::Float64
  ode::Integrator
end

mutable struct LLG_STT <: Driver
  alpha::Float64
  beta::Float64
  gamma::Float64
  ode::Integrator
  tol::Float64
  ux::Array{Float64, 1}
  uy::Array{Float64, 1}
  uz::Array{Float64, 1}
  h_stt::Array{Float64, 1}
end

function create_driver(driver::String, integrator::String, nxyz::Int64)
    if driver=="SD" #Steepest Descent
        gk = zeros(Float64,3*nxyz)
        return EnergyMinimization(gk, 0.0, 1e-4, 1e-12, 0)
    elseif driver=="LLG"
        if integrator == "RungeKutta"
            rungekutta = RungeKutta(nxyz, llg_call_back, 5e-13)
            return LLG(true, 0.1, 2.21e5, rungekutta)
        elseif integrator == "DormandPrince"
            tol = 1e-6
            dopri5 = DormandPrince(nxyz, llg_call_back, tol)
            return LLG(true, 0.1, 2.21e5, dopri5)
        else
            tol = 1e-6
            dopri5 = DormandPrinceCayley(nxyz, llg_cay_call_back, tol)
            return LLG(true, 0.1, 2.21e5, dopri5)
        end
    elseif driver=="LLG_STT"
        tol = 1e-6
        if integrator == "DormandPrince"
            dopri5 = DormandPrince(nxyz, llg_stt_call_back, tol)
        else
            dopri5 = DormandPrinceCayley(nxyz, llg_stt_cay_call_back, tol)
        end
        ux = zeros(nxyz)
        uy = zeros(nxyz)
        uz = zeros(nxyz)
        hstt = zeros(3*nxyz)
        return LLG_STT(0.5, 0.0, 2.21e5, dopri5, tol, ux, uy, uz, hstt)

    else
       error("Supported drivers: SD, LLG, LLG_STT")
    end
    return nothing
end
