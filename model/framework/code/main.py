import mol_gen
import sys
import csv

# # parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

designer = mol_gen.MoleculeModel()

# read SMILES from .csv file, assuming one column with header
with open(input_file, "r") as f:
    reader = csv.reader(f)
    next(reader)  # skip header
    smiles_list = [r[0] for r in reader]

# Convert to SAFE
safe_list = [designer.smiles_to_safe(smi) for smi in smiles_list]

# run model
outputs = designer.run_model(safe_list)

#check input and output have the same length
input_len = len(smiles_list)
output_len = len(outputs)
assert input_len == output_len

N_COLS = 1000
HEADER = ["smiles_{0}".format(str(x).zfill(3)) for x in range(N_COLS)]

with open(output_file, "w") as fp:
    csv_writer = csv.writer(fp)
    # First Row: Header
    # Second Row onwards: Generated Smiles (Output)
    csv_writer.writerows([HEADER])
    for o in outputs:
        print(len(o))
        csv_writer.writerow(o)
