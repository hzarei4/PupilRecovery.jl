export CreatePupil, CreatePSF, AddingNoise


function CreatePupil(;phase, px=40, upsampling=4)
    x = LinRange(-1*upsampling, upsampling, px)
    mask = evaluateZernike(x, [0], [1.0], index=:OSA) 

    if phase == "Random"
        matx = mask .* exp.(2*pi *  rand(Float64, (px, px)) * im)
        matx[abs.(matx) .< abs(eps())] .= zero(eltype(matx))

    elseif nameof(typeof(phase)) == :Array
        x = LinRange(-1*upsampling, upsampling, size(phase)[1])
        mask = evaluateZernike(x, [0], [1.0], index=:OSA) 
        
        ϕ = (2 .* pi .* phase)  .- pi
        matx = mask .* exp.(ϕ * im)
        matx[abs.(matx) .< abs(eps())] .= zero(eltype(matx))
    elseif phase == "zernike4"
        zer4 = evaluateZernike(x, [2, 7, 8, 9], [0.5, 0.2, 0.1, 0.2], index=:OSA)
        matx = mask .* exp.(zer4 * im)
        matx[abs.(matx) .< abs(eps())] .= zero(eltype(matx))
    end
    return matx
end

function CreatePSF(Pupil)
    abs2.(fftshift(fft((Pupil))))
end

function AddingNoise(psf; scaling=0.1)
    return psf .+ scaling .*maximum(psf).*rand(size(psf)[1], size(psf)[2])
end
