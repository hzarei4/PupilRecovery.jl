module PupilRecovery

using Revise
using ZernikePolynomials
     

function CreatePupil(; px=120)

    x = LinRange(-2, 2, px)
    mask = evaluateZernike(x, [1], [1.0], index=:OSA) 
    matx = mask .* exp.(2*pi *  rand(Float64, (px, px)) * im)
    matx[abs.(matx) .< abs(eps())] .= zero(eltype(matx))
    return matx

end



end # module PupilRecovery
