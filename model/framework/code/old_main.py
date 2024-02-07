# import libraries
import os
import sys
os.environ["TOKENIZERS_PARALLELISM"] = "false"
import safe as sf
# import datamol as dm
from safe.utils import compute_side_chains
import csv
from rdkit import Chem
from rdkit import Chem, rdBase
from rdkit.Chem import Draw
from rdkit.Chem.Scaffolds import rdScaffoldNetwork
from rdkit.Chem import Descriptors
# import matplotlib as mpl
from rdkit import Chem

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
    # except sf.EncoderError:
    except:
        print("Error in SMILES conversion")


def extract_core_structure(SMILES):
    # Define scaffold parameter network
    params = rdScaffoldNetwork.ScaffoldNetworkParams()
    # customize parameter attributes
    params.includeScaffoldsWithoutAttachments=False
    mol = Chem.MolFromSmiles(SMILES)
    net = rdScaffoldNetwork.CreateScaffoldNetwork([mol],params)
    nodemols = [Chem.MolFromSmiles(x) for x in net.nodes]

    filtered_list = []
    for mol in nodemols:
        # Check for the presence of attachment points and molecular weight range
        if "*" in Chem.MolToSmiles(mol) and 60 < Descriptors.MolWt(mol) < 100:
            filtered_list.append(mol)
    
    # If there are no scaffolds within the range, select the closest one
    if not filtered_list:
        closest_mol = min(nodemols, key=lambda x: abs(Descriptors.MolWt(x) - 80))
        filtered_list.append(closest_mol)
    
    # Sort the filtered list based on the number of heteroatoms (fewer carbons)
    filtered_list.sort(key=lambda x: x.GetNumHeavyAtoms())

    return filtered_list


# generate new smiles
def generate_smiles(side_chain):
    generated_smiles = designer.scaffold_morphing(
        side_chains=side_chain,
        n_samples_per_trial=12,
        n_trials=1,
        sanitize=True,
        do_not_fragment_further=False,
        random_seed=100,)
    return generated_smiles


def get_side_chain_pairs(side_chains):
    '''Function to break the side chains into pairs'''
    # side_chains_pairs = []
    side_chains_smiles = Chem.MolToSmiles(side_chains).split(".")

    # Generating all combinations of 2 and 3 elements
    combinations_2 = list(combinations(side_chains_smiles, 2))
    combinations_3 = list(combinations(side_chains_smiles, 3))
    all_combinations = combinations_2 + combinations_3

    # adjust the side  chain numbering
    modified_side_chain_pairs = []

    for i in all_combinations:
        modified_strings = []
        for index, j in enumerate(i):
            # Replace the second character with numbers starting from 1
            new_string = j[0] + str(index + 1) + j[2:]
            modified_strings.append(new_string)
        modified_side_chain_pairs.append(modified_strings)

    # join the side chains in each list
    joined_side_chain_pairs = []

    for i in modified_side_chain_pairs:
        joined_strings = ['.'.join(j for j in i)]
        joined_side_chain_pairs.append(joined_strings[0])

    return joined_side_chain_pairs



# my model
def my_model(SMILES):
    generated_smiles = []
    for i in SMILES:
        row = []
        # extract the core structures for each SMILES
        core_structures = extract_core_structure(i)
        # generate new molecules for each side chain of the smile
        for core in core_structures:
            # compute side chain
            side_chain = compute_side_chains(core=core, mol=i)
            # generate new molecules for each side chain of the smile
            output = generate_smiles(side_chain)
            row += output
        generated_smiles += [row]

    return generated_smiles 


# read SMILES from .csv file, assuming one column with header
with open("data/my_molecules.csv", "r") as f:
    reader = csv.reader(f)
    next(reader)  # skip header
    smiles_list = [r[0] for r in reader]

# Load pre-trained model
designer = sf.SAFEDesign.load_default(verbose=True)
# designer.model

# # convert to SAFE
safe_list = [smiles_to_safe(smi) for smi in smiles_list]

# run model
outputs = my_model(safe_list)
# remove duplicate smiles
outputs = list(set(outputs))

# #check input and output have the same lenght
# input_len = len(smiles_list)
# output_len = len(outputs)
# assert input_len == output_len

# write output in a .csv file
with open(output_file, "w") as f:
    writer = csv.writer(f)
    writer.writerow(["value"])  # header
    for o in outputs:
        writer.writerow([o])
