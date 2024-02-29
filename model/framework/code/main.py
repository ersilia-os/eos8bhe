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

# write output in a .csv file
N_COLS = 1000
with open(output_file, "w") as f:
    writer = csv.writer(f)
    h = ["smi_{}".format(i) for i in range(N_COLS)]
    writer.writerow(h)
    for o in outputs:
        o = o[:N_COLS]
        if len(o) < N_COLS:
            o = o + [None] * (N_COLS - len(o))
        writer.writerow(o)
