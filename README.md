# G4review

### Run rG4-detector on G4detector data
#### Data preparation
From G4 detector project, download human data (pos and neg)
```
python data_gen.py
```

####  Run prediction
Clone rG4-detector project and run from the same directory 
```
python run_rg4_detector/predict_fasta.py
```

### Run rG4-detector on G4detector data
#### Data preparation
From G4 detector project, download human data (pos and neg)
```
python run_g4_detector/prepare_rg4_data/csv_data/extract_most_prominent_transcript.py
```

The output will be used to create the seq files:
```
python run_g4_detector/prepare_rg4_data/csv_data/csv2seq.py
```

####  Run prediction
Clone G4 project and run
```
python g4_inference.py
```

### Interpretability
To run models interpretability you'll need the following parameters:
"-rg4", dest="rg4_path", help="path to g4 data and rg4 predicts")
    parser.add_argument("-g4", dest="g4_path", help="path to rg4 data and g4 predicts"

```
python interpretability.py -rg4 <path to g4 data and rg4 predicts output> -g4 <path to rg4 data and g4 predicts>
```

for example:
python interpretability.py -rg4 'run_rg4_detector/' -g4 'run_g4_detector/'