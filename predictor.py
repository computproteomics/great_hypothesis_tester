import os
import csv
import pandas as pd
import numpy as np

from koinapy import Koina

"""
Global constant values
"""
KOINA_MODEL = 'Prosit_2019_irt'
KOINA_SERVER = 'koina.wilhelmlab.org:443'
DEEPLC_QUERY = 'deeplc.csv'
DEEPLC_RESULT = DEEPLC_QUERY.split(sep='.')[0] + '_deeplc_predictions.csv'


class Predictor:


    def predict_retention_time(self, input_file: str, result_file: str) -> None:

        """
        Predicts the retention time of peptide sequences using the deep learning models of 'Prosit' and 'DeepLC'.

        Parameters
        ----------
        input_file: str
            The path to the input .csv file which contains the sequences.

        result_file: str
            The path to the output .csv file where the predictions should be written to.

        Raises
        ------
        AttributeError
            If the given input file does not exist.
        """

        if not os.path.exists(input_file):
            raise AttributeError(f'File {input_file} does not exist!')
        
        # read csv file
        rt_input = pd.read_csv(input_file)

        # ----------------------------------------------------------------- #
        # use koina api for prediction                                      #
        # ----------------------------------------------------------------- #

        # map inputs
        koina_input = pd.DataFrame()
        koina_input['peptide_sequences'] = rt_input['seq'].values

        # build predictor
        koina_predictor = Koina(KOINA_MODEL, KOINA_SERVER)

        # make predictions
        koina_predictions = koina_predictor.predict(koina_input)['irt'].flatten()

        # ----------------------------------------------------------------- #
        # use deeplc for prediction                                         #
        # ----------------------------------------------------------------- #

        deeplc_field = ['seq', 'modifications']
        deeplc_rows = [[i, ''] for i in rt_input['seq'].values]

        # write to temporary csv file
        with open(DEEPLC_QUERY, 'w') as tmp:
            csvwriter = csv.writer(tmp)
            csvwriter.writerow(deeplc_field)
            csvwriter.writerows(deeplc_rows)

        # load predictions
        os.system(f'deeplc --file_pred {DEEPLC_QUERY}')

        # parse results
        deeplc_prediction = np.array(pd.read_csv(DEEPLC_RESULT).iloc[:, 1].values)

        os.remove(DEEPLC_QUERY)
        os.remove(DEEPLC_RESULT)

        rt_input['Prosit Prediction'] = koina_predictions
        rt_input['DeepLC Prediction'] = deeplc_prediction

        rt_input.to_csv(result_file)

