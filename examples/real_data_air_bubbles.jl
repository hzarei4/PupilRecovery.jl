using Images
using FileIO
using Glob

# Path to your directory
path_to_images = "/home/hossein/Data/PSF/20220815_beads/AirBubbleStack_LowConc unstacked"

# Get all TIFF files in the directory
tiff_files = glob("*.tif", path_to_images)

# Initialize an accumulator image with the same size as the first image
accumulator = load(joinpath(path_to_images, tiff_files[1]))
accumulator = Float64.(red.(accumulator))  # Convert to Float64 for normalization

# Normalize and add the first image to the accumulator
accumulator ./= maximum(accumulator)

# Loop over all TIFF files and read them
for i in 2:length(tiff_files)
    # Create the full file path
    file_path = joinpath(path_to_images, tiff_files[i])

    # Load the image
    img = load(file_path)
    img = Float64.(red.(img))  # Convert to Float64 for normalization

    # Normalize the image
    img ./= maximum(img)

    # Add the image to the accumulator
    accumulator .+= img
end

# Normalize the accumulator
accumulator ./= maximum(accumulator)