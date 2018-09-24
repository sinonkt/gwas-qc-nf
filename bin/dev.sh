DEV_HOME=/home/dev/Code
DEV_IMAGE=sinonkt/docker-centos7-singularity-nextflow:latest

docker run -it -u dev --volume $(pwd):$DEV_HOME --privileged $DEV_IMAGE /bin/bash