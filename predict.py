from dlomix.models.prosit import PrositRetentionTimePredictor
from dlomix.data import RetentionTimeDataset
from typing import Union
from koinapy import Koina
import numpy as np
import pandas as pd
import os

# models: https://figshare.com/projects/Prosit/35582

PEPTIDE_CHARGE = 2
FRAGMENTATION_TYPE = 'HCD'
COLLISION_ENERGY = 25

def predict(input_file: str, model_type: str = 'prosit', prediction_type: str = 'retention_time') -> Union[np.array, dict]:

    # check if input file exists
    if not os.path.exists(input_file):
        raise Exception('File not found!')
    
    # check the model and prediction type that should be used
    if model_type == 'prosit':

        if prediction_type == 'retention_time':
            
            prosit_retention_time_predictor = PrositRetentionTimePredictor()
            prosit_retention_time_predictor.load_weights('pretrained_models/prosit_retention_time/')
            
            extracted_data = RetentionTimeDataset(data_source=input_file,
                                                  seq_length=30,
                                                  batch_size=32,
                                                  test=True)
            
            prediction = prosit_retention_time_predictor.predict(extracted_data.test_data).ravel()
        
        elif prediction_type == 'fragmentation_pattern':

            sequences = pd.read_csv(input_file)['sequence'].to_list()

            inputs = pd.DataFrame()
            inputs['peptide_sequences'] = np.array(sequences)
            inputs['precursor_charges'] = np.array([PEPTIDE_CHARGE] * len(sequences))
            inputs['collision_energies'] = np.array([COLLISION_ENERGY] * len(sequences))

            model = Koina('Prosit_2020_intensity_HCD', 'koina.wilhelmlab.org:443')
            prediction = np.array(model.predict(inputs))

        else:
            raise Exception('Unkown prediction type!')
        
    elif model_type == 'deeplc':

        pass

    else:
        raise Exception('Unknown prediction model!')
    
    return prediction
    
pred = predict('test.csv', prediction_type='fragmentation_pattern')
print(pred)