#!/bin/bash

PROJECT_DIR=/storage1/fs1/christophermaher/Active/maherlab/saha.d/projects/prostate_ploidy

############################
# MOUNT PROJECT 
############################

export LSF_DOCKER_VOLUMES="\
${PROJECT_DIR}:${PROJECT_DIR} "


bsub \
  -J sparc \
  -G compute-christophermaher \
  -g /saha.d/max100 \
  -q general \
  -n 4 \
  -R 'select[mem>16000] rusage[mem=16000]' \
  -M 16000000 \
  -oo Logs/sparc.out \
  -eo Logs/sparc.err \
  -a 'docker(openjdk:11.0.11-jdk-slim)' /usr/local/openjdk-11/bin/java -Xms4g -Xmx12g \
       -Dconfig.file=cromwell_compute1.conf \
       -jar cromwell-86.jar \
       run "${PROJECT_DIR}"/R/crpc_analysis/sparc.wdl \
       --inputs sparc.json
