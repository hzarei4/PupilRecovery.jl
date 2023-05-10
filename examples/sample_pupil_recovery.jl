using PupilRecovery
using Plots
using View5D


pupil = CreatePupil(px=60, upsampling=4)
@vt angle.(pupil)

psf = CreatePSF(pupil)
psf_n = AddingNoise(psf, scaling=0.05)

@vt psf, psf_n


default(size=(900, 800))


x_solution = DMonPSF(psf_n, it_max=200)

@vv angle.(x_solution)


