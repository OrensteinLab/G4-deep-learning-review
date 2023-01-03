# G4review

### Interpretability
To run models interpretability you'll need the following parameters:
"-rg4", dest="rg4_path", help="path to g4 data and rg4 predicts")
    parser.add_argument("-g4", dest="g4_path", help="path to rg4 data and g4 predicts"

```
python interpretability.py -rg4 <path to g4 data and rg4 predicts output> -g4 <path to rg4 data and g4 predicts>
```

for example:
python interpretability.py -rg4 'run_rg4_detector/' -g4 'run_g4_detector/'