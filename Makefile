# --- Docker Image and Container Configuration ---
# The name of the Docker image
IMAGE_NAME := lsdpca-project

# The name of the running container instance
CONTAINER_NAME := lsdpca

# The RStudio port mapping (HOST_PORT:CONTAINER_PORT)
HOST_PORT := 8787
CONTAINER_PORT := 8787

# The custom user and password, passed as environment variables to 'docker run'
R_USER := test
R_PASSWORD := 1234

# --- Project Volume Mapping ---
# Use 'pwd' (current directory) for the host side of the volume mount.
# Maps local project files to the container user's home directory.
VOLUME_MOUNT := -v $$(pwd):/home/$(R_USER)/project


.PHONY: build run stop clean all

all: build run
	@echo "--- Docker Image Built and Container is Running ---"
	@echo "Access RStudio Server at: http://localhost:$(HOST_PORT)"
	@echo "Username: $(R_USER), Password: $(R_PASSWORD)"
	@echo "To stop and remove the container, run: make stop"
	@echo "To remove the image, run: make clean"

# 1. Build the Docker Image
# Reads the Dockerfile and creates the lsdpca-project image.
build:
	@echo "--- Building Docker Image: $(IMAGE_NAME) ---"
	docker build -t $(IMAGE_NAME) .

# 2. Run the Docker Container
# Starts the container in detached mode (-d) and passes the necessary
# environment variables for the Rocker entrypoint to set the custom user/password.
run:
	@echo "--- Running Container: $(CONTAINER_NAME) ---"
	docker run --rm -d \
		-p $(HOST_PORT):$(CONTAINER_PORT) \
		-e USER=$(R_USER) \
		-e PASSWORD=$(R_PASSWORD) \
		-e RSTUDIO_USER=$(R_USER) \
		$(VOLUME_MOUNT) \
		--name $(CONTAINER_NAME) \
		$(IMAGE_NAME)

# 3. Stop the Container
# Stops the running container (gracefully shutting down RStudio Server).
stop:
	@echo "--- Stopping Container: $(CONTAINER_NAME) ---"
	-docker stop $(CONTAINER_NAME)

# 4. Clean/Kill (Removes the image and ensures the container is stopped)
# This removes the built image to save disk space.
clean: stop
	@echo "--- Cleaning up Docker Image: $(IMAGE_NAME) ---"
	-docker rmi $(IMAGE_NAME)
