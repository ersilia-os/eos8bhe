# imports
import os
import csv
import sys
# from rdkit import Chem
# from rdkit.Chem.Descriptors import MolWt
os.environ["TOKENIZERS_PARALLELISM"] = "false"
import safe as sf
import datamol as dm

# parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# current file directory
# root = os.path.dirname(os.path.abspath(__file__))


# read input csv file
# load pretrained
# loop smiles in list
#   convert smiles to safe 
#   convert safe to side chain
#   use side chain to generate new smiles
#   store result as a list for each smile

# convert SMILES to SAFE
def smiles_to_safe(smiles):
    try:
        safe = sf.encode(smiles)
        return safe
    except sf.EncoderError:
        print("Error in SMILES conversion")


#  convert core structure to side chain
def side_chains(core_structure):
    return core_structure


# generate new smiles
def generate_smiles(side_chains):
    generated_smiles = designer.scaffold_morphing(
        side_chains=side_chains,
        n_samples_per_trial=12,
        n_trials=1,
        sanitize=True,
        do_not_fragment_further=False,
        random_seed=100,)
    return generated_smiles

# my model
def my_model(side_chains_list):
    return [generate_smiles(side_chains) for side_chains in side_chains_list]


# read SMILES from .csv file, assuming one column with header
with open(input_file, "r") as f:
    reader = csv.reader(f)
    next(reader)  # skip header
    smiles_list = [r[0] for r in reader]

# Load pre-trained model
designer = sf.SAFEDesign.load_default(verbose=True)
designer.model

# convert to SAFE
safe_list = [smiles_to_safe(smi) for smi in smiles_list]

# convert molecule to side chains
side_chains_list = [side_chains(sf) for sf in safe_list]

# run model
outputs = my_model(smiles_list)

#check input and output have the same lenght
input_len = len(smiles_list)
output_len = len(outputs)
assert input_len == output_len

# write output in a .csv file
with open(output_file, "w") as f:
    writer = csv.writer(f)
    writer.writerow(["value"])  # header
    for o in outputs:
        writer.writerow([o])
