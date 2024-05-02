# great_hypothesis_tester

## Installation
The *great_hypothesis_tester* framework can easily be installed by simply cloning the Github Repository.

```
git clone https://github.com/computproteomics/great_hypothesis_tester
```

## Usage
The entire prediction pipeline can be run by executing the 'run_hypothesis.py' Python script and providing the path to the .csv file, in which the peptide sequences are stored.

```
python run_hypothesis.py <path_to_csv>
```

This command will request the predictions from all deep learning models and perform the basic analysis on the results. These raw result files can also be found inside the 'results' folder.

The .csv file containing the peptide sequences must be formatted as illustrated below to be processed.

|seq|
|---|
|FWANILGNVT|
|FLVAKIQMCV|
|EVGVQMIEYQ|
|CYCDEDGRSL|
|CIDMEPQHPD|
|PKHQKMWAES|

The length of the sequences can vary without effecting the functionality of the *great_hypothesis_tester*.


