export DMonPSF

"""
    function DM(x, beta, gamma_S, gamma_M, M_in, S_in)

DM calculates a single update step of the Difference Map algorithm 

# Arguments
`x`: The object array as a vector in a vector space \n
`β`: "Slope" parameter or the speed of finding the best vector `x`
`η`: Multiplier to decrease the value of `β`
`γ_S`: Relaxation parameter for the support constraint
`γ_M`: Relaxation parameter for the Fourier space magnitude constraint
`S_in`: The support constraint in the pupil plane (in real space)
`M_in`: The Fourier space magnitude (or PSF values) as a constraint

# Returns
The updated `x` as the `x_new` and the estimated solution as `x_PSRM`
"""
function DM!(x::AbstractArray, x_PS::AbstractArray, β::AbstractFloat, η::AbstractFloat, γ_S::AbstractFloat, γ_M::AbstractFloat, S_in::AbstractArray, M_in::AbstractArray)

    # Performs the support constraint on the object vector in the pupil space
    function P_S(x, S_in)
        x_new = x .* S_in
        return x_new
    end
    
    # Performs the Fourier space magnitude constraint and setting the reference value
    function P_M(x, M_in)
        X = fft(x)

        ctr = size(x).÷2 .+ 1 #  (1,1)
        # enforce phase-only constraint
        # X_new = M_in .* exp.(angle.(X) * im)
        X_new = M_in .* cis.(angle.(X))
    
        x_new = ifft(X_new)
        x_new[ctr...] = real(x_new[ctr...])
        return x_new
    end
    
    # Performs the relaxed Fourier magnitude constraint
    function R_M(x, γ_M, M_in)
        return (1+γ_M) .* P_M(x, M_in) .- γ_M .* x
    end
    
    # Performs the relaxed support constraint 
    function R_S(x, γ_S, S_in)
        return (1+γ_S) .* P_S(x, S_in) .- γ_S .* x
    end


    x_PMRS = P_M(R_S(x, γ_S, S_in), M_in)
    x_PSRM = P_S(R_M(x, γ_M, M_in), S_in)

    # Updating the object vector 
    # x_new = x + η * β .* (x_PMRS - x_PSRM)
    x .+= η * β .* (x_PMRS - x_PSRM)
    x_PS .= x_PSRM
    # return x_new, x_PSRM
    return nothing
end

"""
    function convolution_filter(x, kernel)

Convolving a `kernel` with an array `x` 

"""
function convolution_filter(x, kernel)
    return ifft(fft(x) .* kernel)
end

"""
    function plotresult(pupil, x_sol, psf, mydiff, it, loss, rng)

Plots the results of each loop
"""
function plotresult(x_sol, psf, mydiff, it, loss, rng)
    try
        IJulia.clear_output(true) 
    catch
    end

    Plots.display(plot(

    heatmap(Float64.(angle.(ifftshift(x_sol))), title = "Rec phase, loop $(it), loss: $(loss)", 
        aspect_ratio=1, c=:twilight, clim=(-pi, pi),
        axis=([], false),
        legend = :none
        ), 


    heatmap(abs2.(fftshift(fft(x_sol))), title = "Rec PSF", aspect_ratio=1, legend = :none,
    axis=([], false),
    c=:grays,
    ),

    heatmap(mydiff, 
    title = "PSF diff, scaled", aspect_ratio=1,
    axis=([], false),
    c=:grays,
    ),#, clim=(0.0, 1.0)),

    heatmap(fftshift(psf), title = "Measured PSF", aspect_ratio=1, legend = :none,
    c=:grays,
    axis=([], false),
    )
    

    ))
    sleep(0.001)
    
end


"""
    function DMonPSF(psf; beta=0.7, eta=0.1, it_max=1000, tol=0.05, plotting=false, iterative_plotting=false)

# Arguments
`psf`: The Point Spread Function or the magnitude of the Fourier space (unshifted)
`β`: The β in the Difference Map algorithm, the speed of convergence
`η`: The multiplier of the `β` to decrease the speed of the algorithm and increase the accuracy
`it_max`: Maximum number of iterations 
`tol`: The tolerance level which ensures the convergence of the algorithm and breaks the iteration loop
`plotting`: If you need to plot to see the outputs
`iterative_plotting`: Plotting the iterations so that the output of each iteration is displayed

# Returns
`x_sol`: The estimated pupil function 
`loss_trace`: the trace of the values of the loss function in each iterations 

# Example
See the /examples folder for the example of using this function
"""
function DMonPSF(psf::AbstractArray; β=0.7, η=0.1, it_max=1000, tol=0.05)#, plotting=false, iterative_plotting=false)
    
    # for plotting purposes
    # default(size=(1200, 800))
    
    # the relaxation paramaters
    γ_M = -1.0/β
    γ_S =  1.0/β
    
    # upsampling the PSF
    # psf = upsample2(psf)
    
    psf = fftshift(psf)
    
    antialiasing_mask = fftshift(collect(disc(size(psf), size(psf)./4.0))) .+ 0.0
    
    start_pupil = fftshift(collect(disc(size(psf), size(psf)./12.0))) .+ 0.0
    x = complex(start_pupil) ./ sqrt(sum(abs2.(start_pupil)) * prod(size(psf))) .* sqrt(sum(psf))
    # @show sum(abs2.(fft(x))) 
    # @show sum(psf) 

    x_sol = zeros(ComplexF64, size(psf))
    x_PS = zeros(ComplexF64, size(psf))
    S_in = zeros(Float64, size(psf))
    M_in = zeros(Float64, size(psf))

    C_lp = zeros(size(psf))
    C_lp[10:size(psf)[1]-10, 10:size(psf)[2]-10] .= 1
    C_lp .= ifftshift(C_lp)


    S_in .= copy(start_pupil)
    M_in .= sqrt.(max.(0, psf))

    loss_trace = zeros(it_max)

    # antialiasing_mask = 1 # collect(box(size(x_sol), round.(size(x_sol)./1.5)))
    flag_converged = false
    rng = 1.0


    pupil = zeros(ComplexF64, size(psf))
    psf_meas = zeros(ComplexF64, size(psf))
    mydiff = zeros(Float64, size(psf))
    
    psf_sol = zeros(Float64, size(psf))
    x_mod = zeros(Float64, size(psf))
    loss = 0.0

    for it in 1:it_max

        # x, x_PS .= DM(x, β, η, γ_S, γ_M, S_in, M_in)
        DM!(x, x_PS, β, η, γ_S, γ_M, S_in, M_in)
        
        x_PS[abs.(x_PS) .< abs(eps())] .= zero(eltype(x_PS))
    
        x_sol .= x_PS
        
        # loss function calculation
        psf_sol .= abs2.(fftshift(fft(x_sol))) 
        psf_meas .= fftshift(psf)
        mydiff .= (psf_sol .- psf_meas) ./ maximum(psf)
        # pupil = ifftshift(Float64.(replace(angle.(x_sol), 0.0 => NaN)))
        pupil .= ifftshift(Float64.(angle.(x_sol)))        
        loss_trace[it] = sum(abs2.(mydiff))

        # shrinkwrap
        if it % 10 == 9 
            x_mod .= abs.(convolution_filter((x_sol), C_lp))
            # x_mod = abs.(x_mod)
            x_mod .= x_mod./ maximum(x_mod)
            S_in .= x_mod .> 0.03#8 # magic number! :D
            # S_in .= supp
            S_in .*= antialiasing_mask
        end


        if loss_trace[it]<tol
            # if plotting
            #     # plotresult(x_sol, psf, mydiff, it, loss_trace[it], rng)
            # end
            flag_converged = true

            println("\n\tThe algorithm is converged with these details: \n\n\t\titerations: $(it)\n\t\tloss value=$(loss_trace[it])")
            break
        else
            # if plotting && iterative_plotting
            #     # plotresult(x_sol, psf, mydiff, it, loss_trace[it], rng)
            # end
        end
    end

    if !flag_converged
        println("\n\tThe algorithm's loss function did not satisfy the tolerance threshold in the defined iterations.")

    end

    return x_sol, loss_trace
end