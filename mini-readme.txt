This folder contains the results of various predictions made on {%infile%}.
It has a file structure like this:

/{%infile%}
├── all-predictions.csv
├── all-preds.png
├── hpo-descriptions.csv
├── model-scores.csv
├── top-50-preds.png
├── variable-importances.csv
└── wmean-good-prediction.csv

all-predictions.csv contains information on each variant of the VCF,
starting with the Ensemble Gene ID and the gene name, and then after giving the
variant itself (in terms of the chromosome, location, ref and alt) each
subsequent "HP:XXXXXXX" column gives the predictions from the model build for
that HPO term.

Each model used has an associated self-evaluation score from the bootstrapped
test/train process, these scores are given in model-scores.csv. The variable
importances for each model are found in variable-importances.csv. Each model is
labelled according to its HPO id, descriptions for the HPO ids used can be found
in hpo-descriptions.csv.

The score-weighted average prediction across models with a score of at least 0.4
is given in wmean-good-prediction.csv.

top-50-preds.png is a heatmap of all predictions for the top 50 variants
overall, while all-preds.png is the heatmap for everything.
