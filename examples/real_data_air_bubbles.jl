using Images
using FileIO
using Glob
using Plots
using View5D, FFTW
using PupilRecovery
using InverseModeling, PSFDistiller


# Path to your directory
path_to_images = "/home/hossein/Data/PSF/LowContrast/"

# Get all TIFF files in the directory
tiff_files = glob("*.tiff", path_to_images)

# Initialize an empty array to store images
images = []

# Loop over all TIFF files and read them
for file in tiff_files
    # Create the full file path
    file_path = joinpath(path_to_images, file)

    # Load the image
    img = load(file_path)

    # Append the image to the images array
    push!(images, img)
end

# Now, `images` is an array containing all the TIFF images in the directory
@vt Float64.(images[22])

main_img = Float64.(images[22])

psf =  main_img[247:327, 225:305]

@vt psf

default(size=(900, 800))


mypsf, rois, positions, selected, params, _ = distille_PSF(main_img, roi_size=Tuple(60 .* ones(Int, ndims(main_img))));

x_solution = DMonPSF(mypsf, it_max=1000);


@btime $x_solution = DMonPSF($mypsf, it_max=1000);