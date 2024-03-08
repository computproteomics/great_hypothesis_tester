import os
import numpy as np

from dlomix.models.prosit import PrositRetentionTimePredictor
from dlomix.models.prosit import PrositIntensityPredictor

from dlomix.data import RetentionTimeDataset
from dlomix.data import IntensityDataset


class HypothesisTester():


    def __init__(self, seq_length: int = 30) -> None:
        
        self.seq_length = seq_length


    def predict_retention_time(self, input_file: str, model_type: str = 'prosit') -> np.array:

        # check if given input file exists
        if not os.path.exists(input_file):
            raise Exception('File ' + input_file + ' does not exist!')
        
        # evaluate model that should be used
        if model_type == 'prosit':

            prosit_retention_time_predictor = PrositRetentionTimePredictor()
            prosit_retention_time_predictor.load_weights('pretrained_models/prosit_retention_time/')

            extracted_data = RetentionTimeDataset(data_source=input_file,
                                                  seq_length=self.seq_length,
                                                  test=True)
            
            prediction = prosit_retention_time_predictor.predict(extracted_data.test_data).ravel()

        elif model_type == 'deepcl':

            # TODO
            ...

        else:
            raise Exception('Unknown model type!')

        return prediction


    def predict_intensities(self, input_file: str, model_type: str = 'prosit') -> np.array:
        
        # check if given input file exists
        if not os.path.exists(input_file):
            raise Exception('File ' + input_file + ' does not exist!')
        
        # evaluate model that should be used
        if model_type == 'prosit':

            prosit_intensity_predictor = PrositIntensityPredictor()
            prosit_intensity_predictor.load_weights('pretrained_models/prosit_fragmentation/')

            extracted_data = IntensityDataset(data_source=input_file,
                                              seq_length=self.seq_length,
                                              test=True)
            
            prediction = prosit_intensity_predictor.predict(extracted_data.test_data).ravel()

        elif model_type == 'deeplc':

            # TODO
            ...

        else:
            raise Exception('Unknown model type!')
        
        return prediction



tester = HypothesisTester()
result = tester.predict_intensities('intensity_test.csv')

print(result)