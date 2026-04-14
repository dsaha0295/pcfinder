#!/bin/bash

PROJECT_DIR=/storage1/fs1/christophermaher/Active/maherlab/saha.d/projects/prostate_ploidy

############################
# MOUNT PROJECT 
############################

export LSF_DOCKER_VOLUMES="\
${PROJECT_DIR}:${PROJECT_DIR} "


bsub \
  -J pcfinder \
  -G compute-christophermaher \
  -g /saha.d/max100 \
  -q general \
  -n 4 \
  -R 'select[mem>16000] rusage[mem=16000]' \
  -M 16000000 \
  -oo Logs/pcfinder.out \
  -eo Logs/pcfinder.err \
  -a 'docker(openjdk:11.0.11-jdk-slim)' /usr/local/openjdk-11/bin/java -Xms4g -Xmx12g \
       -Dconfig.file=cromwell_compute1.conf \
       -jar cromwell-86.jar \
       run "${PROJECT_DIR}"/R/crpc_analysis/pcfinder.wdl \
       --inputs pcfinder.json
