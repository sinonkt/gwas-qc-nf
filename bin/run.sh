datetime=$(date -d "today" +"%Y%m%d%H%M")
NX_CONFIG=nextflow.config
LOG_DIR=logs/${datetime}
WORK_DIR=works/

mkdir -p $LOG_DIR

nextflow -log $LOG_DIR/.nextflow.log \
  -C $NX_CONFIG \
  run $1 \
  -w $WORK_DIR \
  -resume \
  -with-report ${LOG_DIR}/report.html \
  -with-dag ${LOG_DIR}/flowchart.svg \
  -with-timeline ${LOG_DIR}/timeline.html \
  --mapped .data/map-ped \
  --bedbimfam .data/bed-bim-fam \
  --output .data/ \
  "${@:2}"

