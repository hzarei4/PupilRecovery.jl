export DMonPSF

"""
    function DM(x, beta, gamma_S, gamma_M, M_in, S_in)

DM calculates a DifferenceMap
"""
function DM(x, beta, eta, gamma_S, gamma_M, M_in, S_in)

    function P_S(x, S_in)
        x_new = x .* S_in
        return x_new
    end
    
    function P_M(x, M_in)
        X = fft(x)

        ctr = (1,1) # size(x).÷2 .+ 1
        # enforce phase-only constraint
        # X_new = M_in .* exp.(angle.(X) * im)
        X_new = M_in .* cis.(angle.(X))
    
        x_new = ifft(X_new)
        x_new[ctr...] = real(x_new[ctr...])
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

    x_new = x + eta*beta .* (x_PMRS - x_PSRM)

    return x_new, x_PSRM
end

function convolution_filter(x, kernel)
    return ifft(fft(x) .* kernel)
end


function DMonPSF(psf; beta=0.7, eta=0.1, it_max=1000, plotting=false)
    psf = upsample2(psf)
    gamma_M = -1.0/beta
    gamma_S =  1.0/beta

    psf = fftshift(psf)
    
    antialiasing_mask = fftshift(collect(disc(size(psf), size(psf)./4.0))) .+ 0.0
    
    # x = rand(ComplexF64, size(psf))
    start_pupil = fftshift(collect(disc(size(psf), size(psf)./12.0))) .+ 0.0
    x = complex(start_pupil) ./ sqrt(sum(abs2.(start_pupil)) * prod(size(psf))) .* sqrt(sum(psf))
    # @show sum(abs2.(fft(x))) 
    # @show sum(psf) 
    # f_1 = zeros(Float64, size(psf)[1])
    x_sol = zeros(ComplexF64, size(psf))
    x_PS = zeros(ComplexF64, size(psf))
    S_in = zeros(Float64, size(psf))
    M_in = zeros(Float64, size(psf))

    C_lp = zeros(size(psf))
    C_lp[10:size(psf)[1]-10, 10:size(psf)[2]-10] .= 1
    C_lp .= ifftshift(C_lp)




    supp = copy(start_pupil)
    # supp .= (supp)

    S_in .= supp
    M_in .= sqrt.(max.(0,psf))

    loss_trace = zeros(it_max)
    beta_trace = zeros(it_max)
    convergernce_test = zeros(it_max)
    # antialiasing_mask = 1 # collect(box(size(x_sol), round.(size(x_sol)./1.5)))
    flag = false
    rng = 0.3

    for it in 1:it_max

        # if it >200 && it<400
        #     beta = 0.2
        # elseif it>=400 && it<600
        #     beta = 0.1
        # elseif it>800
        #     beta=0.05
        # end


        x, x_PS = DM(x, beta, eta, gamma_S, gamma_M, M_in, S_in)
        
        x_PS[abs.(x_PS) .< abs(eps())] .= zero(eltype(x_PS))
    
        x_sol .= x_PS
        


        psf_sol = abs2.(fftshift(fft(x_sol))) 
        psf_meas = fftshift(psf)
        mydiff = (psf_sol .- psf_meas) ./ maximum(psf)
        # pupil = ifftshift(Float64.(replace(angle.(x_sol), 0.0 => NaN)))
        pupil = ifftshift(Float64.(angle.(x_sol)))
        loss = sum(abs2.(mydiff))
        beta_trace[it] = beta
        
        if false #it > 100
            if loss > loss_trace[it-1]
                beta /= 1.01
                @show beta  
                loss_trace[it] = loss_trace[it-1]

                continue
            end


        end
        if false #it>150 && loss < 0.005

            if flag==false
                flag = true
                @show it 
            end
            
        end


        loss_trace[it] = loss

        # Shrinkwrap
        if it % 10 == 9 
            x_mod = convolution_filter((x_sol), C_lp)
            x_mod = abs.(x_mod)
            x_mod .= x_mod./ maximum(x_mod)
            supp .= x_mod .> 0.038
            S_in .= supp
            S_in .*= antialiasing_mask
        end

        if plotting &&  loss<0.005 #it>100 && plotting && abs(loss_trace[it-1]-loss)<0.0000001 #&&  it % it_max == 0
            # @show abs(loss-loss_trace[it-1])
            Plots.display(plot(

                heatmap(pupil, title = "Rec phase, loop $(it), loss: $(loss)", 
                    aspect_ratio=1, c=:twilight, clim=(-rng*pi, rng*pi), legend = :none), 
                heatmap(abs2.(fftshift(fft(x_sol))), title = "Rec PSF", aspect_ratio=1, legend = :none),
                heatmap(fftshift(psf), title = "True PSF, β=$(beta)", aspect_ratio=1, legend = :none),
                
                heatmap(mydiff, 
                    title = "PSF diff, scaled", aspect_ratio=1),#, clim=(0.0, 1.0)),
                ))
            sleep(0.001)
            @show it
            break
        end


    end
    return x_sol, S_in, loss_trace
end