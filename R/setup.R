# --- setup.R ---

# 1. Define the complete list of required packages for the project
required_packages <- c(
    "RSpectra", "mvtnorm", "Matrix", "ggplot2",
    "astsa", "lpSolve", "lattice", "waveslim",
    "dplR", "complexplus", "MASS",
    "devtools", "gradfps", "multitaper", "LSPCA"
)

# 2. Function to check, load, and/or install packages
# Since we rely on the Dockerfile for installation, this focuses on loading and checking.
load_all_packages <- function(package_list) {
    message("Checking and loading required packages...")
    
    # Check if packages are installed (optional but good practice)
    missing_packages <- package_list[!(package_list %in% installed.packages()[,"Package"])]

    if(length(missing_packages) > 0) {
        stop(paste0("The following packages are missing (Did you run the Docker build?): ",
                    paste(missing_packages, collapse = ", ")))
    }
    
    # Load all packages
    sapply(package_list, library, character.only = TRUE, quietly = TRUE)
    
    message("All required packages loaded successfully!")
}

# 3. Execute the loading
load_all_packages(required_packages)
