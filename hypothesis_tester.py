"""
Module imports
"""
import csv
import os

import numpy as np
import pandas as pd

from koinapy import Koina


"""
Global constants
"""
RT_QUERY_CSV  = 'retention_time_query.csv'
RT_RESULT_CSV = RT_QUERY_CSV.split(sep='.')[0] + '_deeplc_predictions.csv'
INTENS_QUERY_PEPREC = 'intensity_query.peprec'
INTENS_RESULT_CSV = INTENS_QUERY_PEPREC.split(sep='.')[0] + '_HCD_predictions.csv'
AA_ALPHABET   = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
CHARGE_ONE_HOTS = {2: [0,0,1,0,0,0]}


class HypothesisTester():

    """
    A class that is used to test a hypothesis by predicting the retention time and intensities of a peptide sequence by
    using deep learning models.

    Attributes
    ----------
    collision_energy: int
        The collision energy that theoretically was used (default is 25).
                
    precursor_charge: int
        The charge of the precursor (default is 2).

    max_seq_length: int
        The maximum length a peptide sequence can have in the FASTA file (default is 30).

    sequence_lengths: array
        The lengths of every sequence that was read in.

    intensity_stats: pd.DataFrame
        A DataFrame that contains the intensity results from the experiment.

    rt_stats: pd:DataFrame
        A DataFrame that contains the retention time results from the experiment.

    sequences: dict
        All sequences with their names that were specified in the FASTA file.

    Methods
    -------
    run_hypothesis(input_file: str, model_type: str = 'prosit')
        Runs the whole pipeline from reading the input to printing the output.

    parse_fasta_sequences(file_path: str)
        Reads and saves peptide sequences from FASTA file.

    predict_retention_time(model_type: str = 'prosit')        
        Predicts the retention time of peptide sequences using deep learning models.

    predict_intensities(model_type:str = 'prosit')
        Predicts the y- and b-ion intensities for peptide sequences using deep learning models.

    print_stats()
        Prints the experimental statistics in a formatted way.
    """

    def __init__(self, collision_energy: int = 25, precursor_charge = 2, max_seq_length = 30) -> None:

        """
        Parameters
        ----------
        collision_energy: int, optional
            The collision energy that theoretically was used (default is 25).
        
        precursor_charge: int, optional
            The charge of the precursor (default is 2).

        max_seq_length: int, optional
            The maximum length a peptide sequence can have in the FASTA file (default is 30).
        """
        
        self.collision_energy = collision_energy
        self.max_seq_length = max_seq_length
        self.charge = precursor_charge

        self.sequence_lengths = None    
        self.intensity_stats = None
        self.sequences = None
        self.rt_stats = None


    def run_hypothesis(self, input_file: str, model_type: str = 'prosit'):

        """
        Runs the whole pipeline from reading the input to printing the output.

        Parameters
        ----------
        input_file: str
            The path where the FASTA file is located.

        model_type: str, optional
            The deep learning model that should be used. Must either be 'prosit' or 'deeplc/ms2pip' (default is 'prosit').

        Raises
        ------
        FileNotFoundError
            If FASTA file is not found.

        ValueError
            If model type is not valid.
        """

        # check if input file exists
        if not os.path.exists(input_file):
            raise FileNotFoundError(f'File {input_file} not found!')
        
        # check if model type is valid
        if model_type != 'prosit' and model_type != 'deeplc/ms2pip':
            raise ValueError('Unknown model type!')
        
        # run whole pipeline from reading input to printing results
        self.parse_fasta_sequences(input_file)
        self.predict_retention_time(model_type=model_type if model_type == 'prosit' else 'deeplc')
        self.predict_intensities(model_type=model_type if model_type == 'prosit' else 'ms2pip')
        self.print_stats()


    def parse_fasta_sequences(self, file_path: str) -> None:

        """
        Reads and saves peptide sequences from FASTA file.

        Parameters
        ----------
        file_path: str
            The path where the FASTA file is located.

        Raises
        ------
        FileNotFoundError
            If given FASTA file was not found.
        ValueError
            If any letter is not in the one-letter-code of the amino acids.
        """

        # check if file path exists
        if not os.path.exists(file_path):
            raise FileNotFoundError(f'File {file_path} not found!')
        
        self.sequences = {}
        current_sequence_id = None
        current_sequence = ''

        # read file content
        with open(file_path, 'r') as content:

            for line in content:

                line = line.strip().upper()

                # check if new sequence starts
                if line.startswith('>'):

                    # add sequence if new one starts
                    if current_sequence_id:
                        self.sequences[current_sequence_id] = current_sequence

                    # init new sequence
                    current_sequence_id = line[1:]
                    current_sequence = ''

                else:

                    # check if all letters are in alphabet
                    if (all(aa in AA_ALPHABET for aa in line)):

                        # concat sequence
                        current_sequence += line

            # add last sequence
            if current_sequence_id:
                self.sequences[current_sequence_id] = current_sequence
            
        # save lengths of sequences
        self.sequence_lengths = [len(i) for i in self.sequences.values()]


    def predict_retention_time(self, model_type: str = 'prosit') -> None:

        """
        Predicts the retention time of peptide sequences using deep learning models.

        Parameters
        ----------
        model_type: str, optional
            The deep learning model that should be used. Must be either 'prosit' or 'deeplc' (default is 'prosit')

        Raises
        ------
        ValueError
            If given model type is not 'prosit' or 'deeplc'.
        AttributeError
            If sequences have not been read in yet.
        """

        # check if sequences have already been read in
        if self.sequences is None or len(self.sequences) == 0:
            raise AttributeError('No sequences found!')

        # evaluate model that should be used
        if model_type == 'prosit':

            # map inputs to koina format
            rt_inputs = pd.DataFrame()
            rt_inputs['peptide_sequences'] = np.array(list(self.sequences.values()))

            # build predictor
            rt_predictor = Koina('Prosit_2019_irt', 'koina.wilhelmlab.org:443')

            # make predictions
            rt_pred = rt_predictor.predict(rt_inputs)['irt'].flatten()

        elif model_type == 'deeplc':

            # write input data to temporary .csv file
            rt_fields = ['seq', 'modifications']

            rt_rows = [[i, ''] for i in self.sequences.values()]

            with open(RT_QUERY_CSV, 'w') as tmp:
                csvwriter = csv.writer(tmp)
                csvwriter.writerow(rt_fields)
                csvwriter.writerows(rt_rows)

            # load predictions
            os.system(f'deeplc --file_pred {RT_QUERY_CSV}')

            # parse results
            rt_pred = np.array(pd.read_csv(RT_RESULT_CSV).iloc[:, 1].values)

            #remove tmp files
            os.remove(RT_QUERY_CSV)
            os.remove(RT_RESULT_CSV)

        else:
            raise ValueError('Unknown model type!')
        
        # build statistics
        self.__build_retention_time_stats__(rt_pred, model_type)


    def predict_intensities(self, model_type:str = 'prosit') -> None:

        """
        Predicts the y- and b-ion intensities for peptide sequences using deep learning models.

        Parameters
        ----------
        model_type: str, optional
            The deep learning model that should be used. Must be either 'prosit' or 'ms2pip' (default is 'prosit')

        Raises
        ------
        ValueError
            If given model type is not 'prosit' or 'ms2pip'.
        AttributeError
            If sequences have not been read in yet.
        """

        # check if sequences have been read in
        if self.sequences is None or len(self.sequences) == 0:
            raise AttributeError('No sequences found!')
        
        # evaluate model that should be used
        if model_type == 'prosit':

            # map inputs to koina format
            intens_inputs = pd.DataFrame()
            intens_inputs['peptide_sequences'] = np.array(list(self.sequences.values()))
            intens_inputs['precursor_charges'] = np.array([self.charge for _ in self.sequences.keys()])
            intens_inputs['collision_energies'] = np.array([self.collision_energy for _ in self.sequences.keys()])

            # build predictor
            intens_predictor = Koina('Prosit_2020_intensity_HCD', 'koina.wilhelmlab.org:443')

            # make predictions
            intens_pred = intens_predictor.predict(intens_inputs)

            self.__build_intensity_stats__(intens_pred, model_type)

        elif model_type == 'ms2pip':

            # write input data to temporary .peprec file
            intens_fields = ['spec_id', 'modifications', 'peptide', 'charge']

            intens_rows = [[seq_name, '-', seq, self.charge] for seq_name, seq in self.sequences.items()]

            with open(INTENS_QUERY_PEPREC, 'w') as tmp:
                csvwriter = csv.writer(tmp, delimiter=' ')
                csvwriter.writerow(intens_fields)
                csvwriter.writerows(intens_rows)

            # load predictions
            os.system(f'ms2pip -c ms2pip_config.txt {INTENS_QUERY_PEPREC}')

            # parse results
            intens_pred = pd.read_csv(INTENS_RESULT_CSV)

            #remove tmp files
            os.remove(INTENS_QUERY_PEPREC)
            os.remove(INTENS_RESULT_CSV)

        else:
            raise ValueError('Unknown model type!')
        
        # build statistics
        self.__build_intensity_stats__(intens_pred, model_type)
        

    def print_stats(self) -> None:
        
        """
        Prints the experimental statistics in a formatted way.
        """

        # building string for rt stats
        rt_stats_output = \
            '-------------------------------------------------------------------------------------------------------\n' + \
            'Retention Time Statistics: \n\n' + \
            (str(self.rt_stats) if self.rt_stats is not None else 'No statistics found for retention time!') + \
            '\n'
        
        # building string for intensity stats
        intens_stats_output = \
            '-------------------------------------------------------------------------------------------------------\n' + \
            'Intensity Statistics: \n\n' + \
            (str(self.intensity_stats) if self.intensity_stats is not None else 'No statistics found for intensities!') + \
            '\n'
        
        # building string for the end of statistics
        ending_string = \
                    '-------------------------------------------------------------------------------------------------------\n'
        
        print(rt_stats_output, intens_stats_output, ending_string)


    def __build_retention_time_stats__(self, prediction: np.array, model_type: str) -> None:

        """
        Builds the statistics of the results from the retention time prediction.

        Parameters
        ----------
        prediction: np.array
            The results of the prediction as a numpy array.

        model_type: str
            The deep learning model that was used for the prediction.

        Raises
        ------
        ValueError
            If given prediction is empty or model type is not 'prosit' or 'deeplc'.
        """

        if len(prediction) <= 0:
            raise ValueError('No predictions found!')
        
        if model_type != 'prosit' and model_type != 'deeplc':
            raise ValueError('Unknown model type!')

        # parse sequences
        seq_list = list(self.sequences.values())

        # build DataFrame for stats
        data = {'Sequence': seq_list,
                'Sequence Length': [len(i) for i in seq_list],
                'Charge': [self.charge for _ in seq_list],
                'Model': [model_type for _ in seq_list],
                'Retention Time': list(prediction)}

        # set DataFrame
        self.rt_stats = pd.DataFrame(data=data,
                                     index=list(self.sequences.keys()))


    def __build_intensity_stats__(self, prediction, model_type) -> None:
        
        if model_type != 'prosit' and model_type != 'ms2pip':
            raise ValueError('Unknown model type!')

        seq = self.sequences.values()
        seq_lengths = [len(i) for i in seq]
        charge = [self.charge for _ in seq]
        model = [model_type for _ in seq]

        # construct first version of DataFrame
        data = {'Sequence': seq,
                'Sequence Length': seq_lengths,
                'Charge': charge,
                'Model': model}
        
        intens_df = pd.DataFrame(data=data,
                                 index=list(self.sequences.keys()))

        if (model_type == 'prosit'):

            if len(prediction) == 0:
                raise ValueError('No predictions found!')

            koina_df = pd.DataFrame(data={'annotations': np.array(prediction['annotation']).flatten(),
                                          'mz': np.array(prediction['mz']).flatten(),
                                          'intensities': np.array(prediction['intensities']).flatten()})
            
            koina_df = koina_df[koina_df['annotations'].apply(lambda x: x.endswith(b'+1'))]

            for i in range(1, self.max_seq_length):

                y_filtered_df = koina_df[koina_df['annotations'].str.decode('utf-8').str.contains(f'y{i}')]

                b_filtered_df = koina_df[koina_df['annotations'].str.decode('utf-8').str.contains(f'b{i}')]

                intens_df[f'mz_y{i}'] = y_filtered_df['mz'].values[0]
                intens_df[f'mz_b{i}'] = b_filtered_df['mz'].values[1]

                intens_df[f'intensity_y{i}'] = y_filtered_df['intensities'].values[0]
                intens_df[f'intensity_b{i}'] = b_filtered_df['intensities'].values[1]

        else:

            if prediction.empty:
                raise ValueError('No predictions found!')

            # add m/z and intensities
            for i in range(1, self.max_seq_length):

                # extract rows from prediction frame
                extracted_y = prediction.loc[(prediction['ion'] == 'Y') & (prediction['ionnumber'] == i)]
                extracted_b = prediction.loc[(prediction['ion'] == 'B') & (prediction['ionnumber'] == i)]

                y_mz = []
                y_intens = []

                b_mz = []
                b_intens = []
                
                for seq_name in self.sequences.keys():

                    y_mz.append(str(extracted_y.loc[extracted_y['spec_id'] == seq_name, 'mz'].iloc[0]) \
                                if seq_name in extracted_y['spec_id'].values \
                                else '-') 
                    
                    y_intens.append(str(extracted_y.loc[extracted_y['spec_id'] == seq_name, 'prediction'].iloc[0]) \
                                    if seq_name in extracted_y['spec_id'].values \
                                    else '-')

                    b_mz.append(str(extracted_b.loc[extracted_b['spec_id'] == seq_name, 'mz'].iloc[0]) \
                                if seq_name in extracted_b['spec_id'].values \
                                else '-')
                    
                    b_intens.append(str(extracted_b.loc[extracted_b['spec_id'] == seq_name, 'prediction'].iloc[0]) \
                                    if seq_name in extracted_b['spec_id'].values \
                                    else '-')
                    
                intens_df[f'y_{i}_mz'] = y_mz
                intens_df[f'y_{i}_intensity'] = y_intens
                intens_df[f'b_{i}_mz'] = b_mz
                intens_df[f'b_{i}_intensity'] = b_intens

        self.intensity_stats = intens_df


tester = HypothesisTester()
tester.run_hypothesis('sequences.fasta', model_type='deeplc/ms2pip')