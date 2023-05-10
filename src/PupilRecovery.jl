module PupilRecovery

using Revise
using ZernikePolynomials
using FFTW
using Plots

export CreatePupil, CreatePSF, DMonPSF, AddingNoise

function CreatePupil(; random_phase=false, px=40, upsampling=4)
    x = LinRange(-1*upsampling, upsampling, px)
    mask = evaluateZernike(x, [0], [1.0], index=:OSA) 

    if random_phase
        matx = mask .* exp.(2*pi *  rand(Float64, (px, px)) * im)
        matx[abs.(matx) .< abs(eps())] .= zero(eltype(matx))

    else
        ϕ = evaluateZernike(x, [5], [1.0], index=:OSA) 
        matx = mask .* exp.(ϕ * im)
    end

end

function CreatePSF(Pupil)
    abs2.(fftshift(fft(ifftshift(Pupil))))
end

function AddingNoise(psf; scaling=0.1)
    return psf + scaling .*maximum(psf).*rand(size(psf)[1], size(psf)[2])
end


function DM(x, beta, gamma_S, gamma_M, M_in, S_in)

    function P_S(x, S_in)
        x_new = x .* S_in
        return x_new
    end
    
    function P_M(x, M_in)
        X = fft(x)
    
        X_new = M_in .* exp.(angle.(X) * im)
    
        x_new = ifft(X_new)
        return x_new
    end
    
    
    function R_M(x, gamma_M, M_in)
        return (1+gamma_M) .* P_M(x, M_in) .- gamma_M .* x
    end
    
    function R_S(x, gamma_S, S_in)
        return (1+gamma_S) .* P_S(x, S_in) .- gamma_S .* x
    end


    x_PMRS = P_M(R_S(x, gamma_S, S_in), M_in)
    x_PSRM = P_S(R_M(x, gamma_M, M_in), S_in)

    x_new = x + beta .* (x_PMRS - x_PSRM)

    return x_new, x_PSRM
end

function convolution_filter(x, kernel)
    return ifft(fft(x) .* kernel)
end


function DMonPSF(psf; beta=0.99, it_max=10000)
    gamma_M = -1.0/beta
    gamma_S = 1.0/beta
    psf = fftshift(psf)
    x = rand(ComplexF64, size(psf))
    x_sol = zeros(ComplexF64, size(psf))
    S_in = zeros(Float64, size(psf))
    M_in = zeros(Float64, size(psf))

    C_lp = zeros(size(psf))
    C_lp[20:size(psf)[1]-20, 20:size(psf)[2]-20] .= 1
    C_lp .= ifftshift(C_lp)

    supp = zeros(size(psf))
    supp[20:size(psf)[1]-20, 20:size(psf)[2]-20] .= 1
    supp .= supp

    S_in .= supp
    M_in .= psf


    for it in 1:it_max
        x, x_PS = DM(x, beta, gamma_S, gamma_M, M_in, S_in)
        
        x_PS[abs.(x_PS) .< abs(eps())] .= zero(eltype(x_PS))
    
        x_sol .= x_PS
        
        # Shrinkwrap
        if it % 10 ==9
            x_mod = convolution_filter(x_sol, C_lp)
            x_mod = abs.(x_mod)
            x_mod .= x_mod./ maximum(x_mod)
            supp .= x_mod .> 0.05
            S_in = supp
        end
    
        if it % (it_max/10) == 0
            Plots.display(plot(

                heatmap(angle.((x_sol)), title = "Rec phase, loop $(it)", 
                    aspect_ratio=1, c=:twilight, clim=(-1*pi, pi), legend = :none), 
                heatmap(abs.(fftshift(fft(x_sol))), title = "Rec PSF", aspect_ratio=1, legend = :none),
                heatmap(fftshift(psf), title = "True PSF", aspect_ratio=1, legend = :none),

                ))
            sleep(0.2)
        end
    end
    return x_sol, S_in
end


end # module PupilRecovery
