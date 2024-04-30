import os
import subprocess
import sys

from hypothesis_testing.predictor import Predictor

ANALYSIS_FILE_PATH = "./hypothesis_testing/analysis.R"
QUERY_NAME = os.path.basename(sys.argv[1]).split(".")[0]
RT_RESULT_FILE = "./results/{}_results/{}_irt_results.csv".format(
    QUERY_NAME, QUERY_NAME
)
FRAG_RESULT_FILE = "./results/{}_results/{}_fragmentation_results.csv".format(
    QUERY_NAME, QUERY_NAME
)

# check if arguments are given correctly
if len(sys.argv) != 2:
    raise Exception("Invalid number of arguments!")

# run predictions on input file
predictor = Predictor()
#predictor.predict_retention_time(sys.argv[1])
#predictor.predict_fragmentation(sys.argv[1])

# run the analysis script on the prediction data
subprocess.run(f"Rscript {ANALYSIS_FILE_PATH} {RT_RESULT_FILE} {FRAG_RESULT_FILE}")
