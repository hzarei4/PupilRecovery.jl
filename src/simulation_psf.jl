export CreatePupil, CreatePSF, AddingNoise


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
    abs.(fftshift(fft(ifftshift(Pupil))))
end

function AddingNoise(psf; scaling=0.1)
    return psf + scaling .*maximum(psf).*rand(size(psf)[1], size(psf)[2])
end
