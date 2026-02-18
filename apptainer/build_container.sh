
mkdir -p ../../OAS-build-container

# Build the Singularity/Apptainer container
echo "Building building container..."
echo "Rebuild container if there are changes in definition files."
apptainer build ../../OAS-build-container/OAS_build.sif ubuntu_24-dev.def

# Build OOFEM and OAS using the container
echo "Building OAS using the building container..."
cd ../../OAS-build-container
apptainer exec --bind ../OAS:/OAS,$(pwd):$(pwd) OAS_build.sif /bin/bash -c "
cmake -DCMAKE_BUILD_TYPE=Release -DUSE_CHOLMOD=ON /OAS

cmake --build . -- -j 10
"

cd ../OAS/apptainer

# Build minimal container first
echo "Building minimal release container..."
echo "Rebuild container if there are changes in definition files."
apptainer build ../../OAS-build-container/OAS_release.sif ubuntu_24-rel.def

# Build as writable sandbox to add files
echo "Converting to sandbox for file updates..."
apptainer build --sandbox ../../OAS-build-container/OAS_sandbox ../../OAS-build-container/OAS_release.sif

# Copy files into the sandbox
echo "Copying application files..."
cp ../../OAS-build-container/bin/OAS ../../OAS-build-container/OAS_sandbox/bin/OAS

# Rebuild the final SIF from sandbox
echo "Building final container from sandbox..."
apptainer build ../../OAS-build-container/OAS.sif ../../OAS-build-container/OAS_sandbox

# Clean up sandbox
rm -rf ../../OAS-build-container/OAS_sandbox

echo "Final container created: ../../OAS-build-container/OAS.sif"
