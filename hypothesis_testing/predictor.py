import os

import numpy as np
import pandas as pd
from koinapy import Koina

"""
Global constant values
"""
KOINA_PROSIT_2019_RT_MODEL = "Prosit_2019_irt"
KOINA_DEEPLC_RT_MODEL = "Deeplc_hela_hf"
KOINA_ALPHAPEPT_RT_MODEL = "AlphaPept_rt_generic"

KOINA_PROSIT_2019_FRAG_MODEL = "Prosit_2019_intensity"
KOINA_DEEPLC_FRAG_MODEL = "ms2pip_2021_HCD"
KOINA_ALPHAPEPT_FRAG_MODEL = "AlphaPept_ms2_generic"

KOINA_SERVER = "koina.wilhelmlab.org:443"

PRECURSOR_CHARGE = 1
COLLISION_ENERGY = 25
INSTRUMENT_TYPE = "QE"

MAX_SEQ_LENGTH = 30


class Predictor:
    """
    A class that is able to request predictions from different models, provided by the
    Koina-API.

    Methods
    -------
    predict_retention_time(input_file: str, sequence_type: str = "random")
        predicts the retention time of all peptide sequences within the input file
    predict_fragmentation(input_file: str, sequence_type: str = "random")
        predicts the fragmentation pattern of all peptide sequences within the input
        file
    """

    def predict_retention_time(
        self, input_file: str, sequence_type: str = "random"
    ) -> None:
        """
        Predicts the retention time of all peptide sequences within the input file.

        Parameters
        ----------
        input_file: str
            The path of the .csv file, in which the peptide sequences are stored.
        sequences_type: str, optional
            A keyword for the type of sequences that are stored in the .csv file. Only
            important for the sake of understanding of the researchers (default is
            "random").

        Raises
        ------
        AttributeError
            If the given .csv file does not exist.
        """

        if not os.path.exists(input_file):
            raise AttributeError(f"File {input_file} does not exist!")

        print(
            "-----------------------------------------------------------------------------------------"
        )
        print("Starting Retention Time Prediction for " + input_file + " ...")

        # ----------------------------------------------------------------- #
        # prepare retention time input                                      #
        # ----------------------------------------------------------------- #

        # read .csv file
        input = pd.read_csv(input_file)

        rt_input = pd.DataFrame()
        rt_input["peptide_sequences"] = input["seq"].values

        # ----------------------------------------------------------------- #
        # use 2019 Prosit model for prediction                              #
        # ----------------------------------------------------------------- #

        print("Getting Prosit predictions...")

        # build predictor
        prosit_2019_predictor = Koina(KOINA_PROSIT_2019_RT_MODEL, KOINA_SERVER)

        # make predictions
        prosit_2019_prediction = prosit_2019_predictor.predict(rt_input)[
            "irt"
        ].flatten()

        # ----------------------------------------------------------------- #
        # use DeepLC model for prediction                                   #
        # ----------------------------------------------------------------- #

        print("Getting DeepLC predictions...")

        # build predictor
        deeplc_predictor = Koina(KOINA_DEEPLC_RT_MODEL, KOINA_SERVER)

        # make predictions
        deeplc_prediction = deeplc_predictor.predict(rt_input)["irt"].flatten()

        # ----------------------------------------------------------------- #
        # use AlphaPept model for prediction                                #
        # ----------------------------------------------------------------- #

        # build predictor
        # alphapept_predictor = Koina(KOINA_ALPHAPEPT_RT_MODEL, KOINA_SERVER)

        # make predictions
        # alphapept_prediction = alphapept_predictor.predict(rt_input)['irt'].flatten()

        # ----------------------------------------------------------------- #
        # write results to file                                             #
        # ----------------------------------------------------------------- #

        # build result dataframe
        rt_output = pd.DataFrame()
        rt_output["Sequence"] = rt_input["peptide_sequences"]
        rt_output["Sequence Length"] = rt_input["peptide_sequences"].str.len().values
        rt_output["Sequence Type"] = np.array(
            [sequence_type] * len(input["seq"].values)
        )
        rt_output["Prosit Prediction"] = prosit_2019_prediction
        rt_output["DeepLC Prediction"] = deeplc_prediction
        # rt_output['Alpha Pept Prediction'] = alphapept_prediction

        query_name = os.path.basename(input_file).split(".")[0]

        if not os.path.exists("./results/{}_results".format(query_name)):
            os.mkdir("./results/{}_results".format(query_name))

        # write result to .csv file
        result_filename = "./results/{}_results/{}_irt_results.csv".format(
            query_name, query_name
        )

        rt_output.to_csv(result_filename, index=False)

        print("Done! Results can be found in " + result_filename)
        print(
            "-----------------------------------------------------------------------------------------"
        )

    def predict_fragmentation(
        self, input_file: str, sequence_type: str = "random"
    ) -> None:
        """
        Predicts the fragmentation pattern of all peptide sequences within the input
        file.

        Parameters
        ----------
        input_file: str
            The path of the .csv file, in which the peptide sequences are stored.
        sequences_type: str, optional
            A keyword for the type of sequences that are stored in the .csv file. Only
            important for the sake of understanding of the researchers (default is
            "random").

        Raises
        ------
        AttributeError
            If the given .csv file does not exist.
        """

        if not os.path.exists(input_file):
            raise AttributeError(f"File {input_file} does not exist!")

        print(
            "-----------------------------------------------------------------------------------------"
        )
        print("Starting Fragmentation Prediction for " + input_file + " ...")

        # ----------------------------------------------------------------- #
        # prepare fragmentation input                                       #
        # ----------------------------------------------------------------- #

        # read .csv file
        input = pd.read_csv(input_file)

        frag_input = pd.DataFrame()
        frag_input["peptide_sequences"] = input["seq"].values
        frag_input["precursor_charges"] = np.array(
            [PRECURSOR_CHARGE] * len(input["seq"].values)
        )
        frag_input["collision_energies"] = np.array(
            [COLLISION_ENERGY] * len(input["seq"].values)
        )
        frag_input["instrument_types"] = np.array(
            [INSTRUMENT_TYPE] * len(input["seq"].values)
        )

        # ----------------------------------------------------------------- #
        # use Prosit model for prediction                                   #
        # ----------------------------------------------------------------- #

        print("Getting Prosit predictions...")

        # build predictor
        prosit_2019_predictor = Koina(KOINA_PROSIT_2019_FRAG_MODEL, KOINA_SERVER)

        # make predictions
        prosit_2019_prediction = prosit_2019_predictor.predict(frag_input)

        # ----------------------------------------------------------------- #
        # use DeepLC model for prediction                                   #
        # ----------------------------------------------------------------- #

        print("Getting DeepLC predictions...")

        # build predictor
        deeplc_predictor = Koina(KOINA_DEEPLC_FRAG_MODEL, KOINA_SERVER)

        # make predictions
        deeplc_prediction = deeplc_predictor.predict(frag_input)

        # ----------------------------------------------------------------- #
        # use AlphaPept model for prediction                                #
        # ----------------------------------------------------------------- #

        # print('Getting Alphapept predictions...')

        # build predictor
        # alphapept_predictor = Koina(KOINA_ALPHAPEPT_FRAG_MODEL, KOINA_SERVER)

        # alphapept_prediction = alphapept_predictor.predict(frag_input)

        # ----------------------------------------------------------------- #
        # map API responses                                                 #
        # ----------------------------------------------------------------- #

        print("Parsing results...")

        prosit_2019_result_df = self.__map_koina_prediction__(prosit_2019_prediction)
        deeplc_result_df = self.__map_koina_prediction__(deeplc_prediction, "deeplc")
        # alphapept_result_df = self.__map_koina_prediction__(alphapept_prediction, 'alphapept')

        # ----------------------------------------------------------------- #
        # write results to file                                             #
        # ----------------------------------------------------------------- #

        frag_output = pd.DataFrame()
        frag_output["Sequence"] = frag_input["peptide_sequences"]
        frag_output["Sequence Length"] = (
            frag_input["peptide_sequences"].str.len().values
        )
        frag_output["Sequence Type"] = np.array(
            [sequence_type] * len(input["seq"].values)
        )

        frag_output = pd.concat(
            [frag_output, prosit_2019_result_df, deeplc_result_df],
            axis=1,
            ignore_index=False,
        )

        query_name = os.path.basename(input_file).split(".")[0]

        if not os.path.exists("./results/{}_results".format(query_name)):
            os.mkdir("./results/{}_results".format(query_name))

        # write result to .csv file
        result_filename = "./results/{}_results/{}_fragmentation_results.csv".format(
            query_name, query_name
        )

        frag_output.to_csv(result_filename, index=False)

        print("Done! Results can be found in " + result_filename)
        print(
            "-----------------------------------------------------------------------------------------"
        )

    def __map_koina_prediction__(
        self, prediction: dict, predictor: str = "prosit"
    ) -> pd.DataFrame:
        """
        Maps the predictions of the Koina-API into the predefined format of the result
        file.

        Parameters
        ---------
        prediction: dict
            The predictions that were made by the deep learning model.
        predictor: str, optional
            The name of the model that was used for the predictions (default is
            "prosit").

        Raises
        ------
        ValueError
            If predictions are not valid.
        """

        if prediction is None or len(prediction) == 0:
            raise ValueError("No predictions found!")

        frames = []

        # build dataframe for each sequence
        for i in range(len(prediction["intensities"])):

            print("{} of {}".format(i + 1, len(prediction["intensities"])), end="\r")

            # extract information
            current_intensities = prediction["intensities"][i]
            current_mz = prediction["mz"][i]
            current_annotations = prediction["annotation"][i]

            # build df
            current_data = pd.DataFrame(
                data={
                    "seq": [i] * len(current_intensities),
                    "intensities": current_intensities,
                    "mz": current_mz,
                    "annotation": current_annotations,
                }
            )

            # filter ions with single charge
            current_data = current_data[
                current_data["annotation"].apply(lambda x: x.endswith(b"+1"))
            ]

            # pivot df
            current_data = current_data.pivot(
                index="seq", columns="annotation", values=["intensities", "mz"]
            )

            # append to temporary frames
            frames.append(current_data)

        # construct final df
        result_frame = pd.concat(frames, axis=0, ignore_index=True)

        # remove ions that were not found
        not_found_ions = result_frame.columns[result_frame.eq(-1).all()]
        result_frame.drop(columns=not_found_ions, inplace=True)

        # rename columns
        rename_column = (
            lambda column: f"{column[1].decode('utf-8').split('+')[0]}_{'intens' if column[0] == 'intensities' else 'mz'}_{predictor}"
        )

        result_frame.columns = [
            rename_column(column) for column in result_frame.columns
        ]

        return result_frame
