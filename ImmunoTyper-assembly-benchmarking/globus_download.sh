#!/bin/bash
#SBATCH --mail-type=END,FAIL

SOURCE_EBI_UUID=47772002-3e5b-4fd3-b97c-18cee38d6df2

globus transfer --no-verify-checksum --notify succeeded,failed,inactive --skip-source-errors \
  $SOURCE_EBI_UUID $TARGET_UUID --batch globus_urls.txt 
