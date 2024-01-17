import safe

ibuprofen = "CC(Cc1ccc(cc1)C(C(=O)O)C)C"

# SMILES -> SAFE -> SMILES translation
try:
    ibuprofen_sf = safe.encode(ibuprofen)  # c12ccc3cc1.C3(C)C(=O)O.CC(C)C2
    ibuprofen_smi = safe.decode(ibuprofen_sf, canonical=True)  # CC(C)Cc1ccc(C(C)C(=O)O)cc1
except safe.EncoderError:
    pass
except safe.DecoderError:
    pass

ibuprofen_tokens = list(safe.split(ibuprofen_sf))
print(ibuprofen_tokens)