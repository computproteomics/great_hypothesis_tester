import os
import numpy as np
import pandas as pd

from dlomix.models.prosit import PrositRetentionTimePredictor
from dlomix.models.prosit import PrositIntensityPredictor
from dlomix.data import RetentionTimeDataset
from dlomix.data import IntensityDataset

import matplotlib.pyplot as plt

from dlomix.utils import convert_nested_list_to_numpy_array

class HypothesisTester():


    def __init__(self, collision_energy: float = 0.25, precursor_charge_one_hot: list = [0,0,1,0,0,0], fragmentation: float = 2.0, \
        sequence_length: int = 30)-> None:
                
        self.seq_length = sequence_length
        self.ce = collision_energy
        self.charge = precursor_charge_one_hot
        self.frag = fragmentation

        self.retention_time_stats = None
        self.intensity_stats = None


    def run_hypothesis(self, fasta_file: str, model_type: str = 'prosit') -> None:
        
        # run whole pipeline to get predictions for the hypothesis
        self.parse_fasta_sequences(fasta_file)
        self.predict_retention_time(model_type)
        self.predict_intensities(model_type)
        self.print_stats()
       

    def parse_fasta_sequences(self, file_path: str) -> None:

        # check if file path exists
        if not os.path.exists(file_path):
            raise FileNotFoundError(f'File {file_path} not found!')
        
        self.sequences = {}
        current_sequence_id = None
        current_sequence = ''

        # read file content
        with open(file_path, 'r') as content:

            for line in content:

                line = line.strip()

                # check if new sequence starts
                if line.startswith('>'):

                    # add sequence if new one starts
                    if current_sequence_id:
                        self.sequences[current_sequence_id] = current_sequence

                    # init new sequence
                    current_sequence_id = line[1:]
                    current_sequence = ''

                else:

                    # concat sequence
                    current_sequence += line

            # add last sequence
            if current_sequence_id:
                self.sequences[current_sequence_id] = current_sequence
            
        # save lengths of sequences
        self.sequence_lengths = [len(i) for i in self.sequences.values()]

    
    def predict_retention_time(self, model_type: str = 'prosit'):

        # evaluate model that should be used
        if model_type == 'prosit':

            # load pretrained predictor
            prosit_retention_time_predictor = PrositRetentionTimePredictor()
            prosit_retention_time_predictor.load_weights('pretrained_models/prosit_retention_time/')

            # load data from input file into dataset
            extracted_data = RetentionTimeDataset(data_source=np.array(list(self.sequences.values())),
                                                  seq_length=self.seq_length,
                                                  test=True)
                        
            # predict retention time
            prediction = prosit_retention_time_predictor.predict(extracted_data.test_data).ravel()

        elif model_type == 'deeplc':

            ...

        else:
            raise ValueError('Unknown model type!')
        
        # update stats
        self.__build_retention_time_stats__(prediction, model_type)
        
    
    def predict_intensities(self, model_type: str = 'prosit'):

        # evaluate model type that should be used
        if model_type == 'prosit':
            
            # load pretrained predictor
            prosit_intensity_predictor = PrositIntensityPredictor()
            prosit_intensity_predictor.load_weights('pretrained_models/prosit_fragmentation/')

            # map precursor charge to dlomix format
            formatted_charge = [self.charge for _ in range(len(self.sequences))]

            # map data to dlomix format
            formatted_data = tuple(
                [
                    np.array(list(self.sequences.values())),
                    np.array([self.ce] * len(self.sequences)),
                    convert_nested_list_to_numpy_array(formatted_charge, dtype=np.float64),
                    np.array([self.frag] * len(self.sequences))
                ]
            )

            # BIG TODO
            # check how to make without cvs

            # load data from input file into dataset
            extracted_data = IntensityDataset(data_source='intensity_test.csv',
                                              seq_length=self.seq_length,
                                              test=True)
            
            # predict intensities
            prediction = prosit_intensity_predictor.predict(extracted_data.test_data)

            # get dictionary containing tupel with y- and b-ions
            result = self.__postprocess_intensities__(prediction)

        elif model_type == 'deeplc':

            # TODO
            ...

        else:
            raise ValueError('Unknown model type!')
        
        # build intensity stats
        self.__build_intensity_stats__(result, model_type)


    def print_stats(self) -> None:

        # building string for retention time statistics printing
        retention_time_stats_string = \
                    '-------------------------------------------------------------------------------------------------------\n' + \
                    'Retention Time Statistics: \n\n' + \
                    (str(self.retention_time_stats) if self.retention_time_stats is not None else 'No statistics found for retention time!')
        
        # building string for intensity statistics printing
        intensity_stats_string = \
                    '-------------------------------------------------------------------------------------------------------\n' + \
                    'Intensity Statistics: \n\n' + \
                    (str(self.intensity_stats) if self.intensity_stats is not None else 'No statistics found for intensities!')
        
        # building string for the end of statistics
        ending_string = \
                    '-------------------------------------------------------------------------------------------------------\n'
        
        # print all components and plots
        print(retention_time_stats_string)

        self.retention_time_stats.plot(y='Retention Time', kind='bar')
        plt.show()

        print(intensity_stats_string)

        print(ending_string)
                       

    def __build_retention_time_stats__(self, prediction: np.array, model_type: str) -> None:

        # build DataFrame for retention time stats
        data = {'Sequence': list(self.sequences.values()),
                'Sequence Length': [len(i) for i in self.sequences.values()],
                'Model': [model_type for _ in self.sequences.keys()],
                'Retention Time': list(prediction)}

        self.retention_time_stats = pd.DataFrame(data=data,
                                                 index=list(self.sequences.keys()))


    def __build_intensity_stats__(self, prediction: dict, model_type: str) -> None:

        # build DataFrame
        data = {'Sequence': list(self.sequences.values()),
                'Sequence Length': [len(i) for i in self.sequences.values()],
                'Model': [model_type for _ in self.sequences.keys()]}
        
        self.intensity_stats = pd.DataFrame(data=data,
                                            index = list(self.sequences.keys()))

        # find max length of y- and b-ion list
        max_length = 0

        # iterate and compare each length of lists to find longest
        for value in prediction.values():
            max_length = max(max_length, len(value[0]), len(value[1]))

        # add columns to DataFrame
        for i in range(max_length):

            # extract ions
            extracted_y = [str(ions[0][i]) if i < len(ions[0]) else '-' for ions in prediction.values()]
            extracted_b = [str(ions[1][i]) if i < len(ions[1]) else '-' for ions in prediction.values()]

            self.intensity_stats[f'y_{i + 1}'] = extracted_y
            self.intensity_stats[f'b_{i + 1}'] = extracted_b


    def __postprocess_intensities__(self, prediction: np.array) -> dict:

        processed_intensities = { }

        for index, ions in enumerate(prediction):

            y_ions = []
            b_ions = []

            # cut prediction into chunks
            chunks = [ions[i * 6:(i+1) * 6] for i in range((len(ions) + 6 - 1 ) // 6)]

            # cut off ions that cannot be found
            filtered_chunks = chunks[:self.sequence_lengths[index]]

            # extract y and b ions from every chunk
            for chunk in filtered_chunks:

                y_ions.append(chunk[0])
                b_ions.append(chunk[3])

            # save tuple to dict
            processed_intensities[list(self.sequences.keys())[index]] = (y_ions, b_ions)

        return processed_intensities

