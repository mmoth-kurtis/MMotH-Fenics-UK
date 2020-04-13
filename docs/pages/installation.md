Currently, Docker is required to run FEniCS. Docker is a program that creates "containers" that allows code to be run in a controlled environment using the host computer's resources. A switch from Docker to Singularity may be made in the future to allow the code to be executed on a computing cluster.  
  * [Install Docker](#Install Docker)
  * [Load the MMotH-Fenics-UK](#Load Image)
  * [Create a Docker container with the loaded image](#Create Container)

## Install Docker
Install the latest version of [Docker](http://www.docker.com).

## Clone the MMoth-Fenics-UK Repository
All of the source code to run MMotH-Vent is located on a [GitHub repository](https://github.com/mmoth-kurtis/MMotH-Fenics-UK.git). Users experienced with Git can do this through a command line approach. Otherwise, a .zip file from the repository can be downloaded from the repository. Unzip the file in the desired directory.

## Load Image  
A Docker image is a copy of the environment used to execute the code. This allows standardization of the modules and their versions used by MMoth-Vent. The image that needs to be loaded by Docker is in the MMotH-Fenics-UK directory, saved as ```MMotH-Vent.tar```. Navigate to where this file is saved on your machine, and execute the following:  
```
docker load < MMotH-Vent.tar
```


## Create Container
