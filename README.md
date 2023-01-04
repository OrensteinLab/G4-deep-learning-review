# G4review

### Run rG4-detector on G4detector data
#### Data preparation
From G4 detector project, download human data (pos and neg)
https://github.com/OrensteinLab/G4detector/tree/G4detector/benchmarks
This script generate test.fa and test.csv files that will be in used on the prediction
```
python data_gen.py -d <path to G4 data> -o <path to output file>
```

####  Run prediction
Clone rG4-detector project and run from the same directory 
```
python run_rg4_detector/predict_fasta.py -f <path to G4 fasta file> -m <path to rG4 model> -o <path to output>
```

### Run G4-detector on rG4detector data
#### Data preparation
From rG4 detector project, download human data
https://github.com/OrensteinLab/rG4detector/tree/main/data/human/csv_data
```
python code/prepare_rg4_data/csv_data/extract_most_prominent_transcript.py <pat to gff3 file> <output_dir>
```
The output is transcripts_locations.pkl file that will be the input to the next step to create the seq files:

Run under csv_data directory
```
python code/prepare_rg4_data/csv_data/csv2seq.py <chr_dict_path - transcripts_locations.pkl> <transcripts_fasta_path>
```
The output is test.fa file that wiil be used to run the prediction

####  Run prediction
Clone G4 project and run
```
python code/g4_inference.py -d <path to rG4 data fasta file> -m <path to G4 model> -o <path to output dir>
```

### Interpretability
To run models interpretability you'll need the following parameters:
"-rg4", dest="rg4_path", help="path to g4 data and rg4 predicts")
    parser.add_argument("-g4", dest="g4_path", help="path to rg4 data and g4 predicts"

```
python interpretability.py -rg4 <path to g4 data and rg4 predicts output> -g4 <path to rg4 data and g4 predicts>
```