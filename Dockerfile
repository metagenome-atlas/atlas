# Start with the Miniconda base image
FROM continuumio/miniconda3:24.9.2-0

# Set the working directory in the container
WORKDIR /main

# Copy the environment file and project code
COPY atlasenv.yml .

# Create a user with a specific UID and GID
RUN groupadd -g 1000 atlasgroup && \
    useradd -m -u 1000 -g atlasgroup -s /bin/bash atlasuser

# Set the HOME environment variable
ENV HOME=/home/atlasuser

# Change ownership of the home directory
RUN chown -R atlasuser:atlasgroup $HOME

# Switch to the new user
USER atlasuser

# Create and activate the environment
RUN conda env create -n atlas -f atlasenv.yml && \
    conda clean -afy && \
    echo "source activate atlas" > ~/.bashrc

# Set the working directory
WORKDIR /main


# Set the default command
CMD ["bash"]