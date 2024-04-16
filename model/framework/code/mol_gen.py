# import libraries
import os
os.environ["TOKENIZERS_PARALLELISM"] = "false"
import safe as sf
# import datamol as dm
from safe.utils import compute_side_chains
from rdkit import Chem
from rdkit.Chem.Scaffolds import rdScaffoldNetwork
from rdkit.Chem import Descriptors
from rdkit import Chem
from itertools import combinations


class MoleculeModel:
    def __init__(self, n_trials=10, n_samples_per_trial=100, lower_molecular_weight=60, upper_molecular_weight=100):
        self.designer = sf.SAFEDesign.load_default(verbose=True)
        self.n_trials = n_trials
        self.n_samples_per_trial = n_samples_per_trial
        self.lower_molecular_weight = lower_molecular_weight
        self.upper_molecular_weight = upper_molecular_weight

    def smiles_to_safe(self, smiles):
        try:
            return sf.encode(smiles)
        except Exception as e:
            print(f" Error in SMILES conversion: {e}")
            return None

    def _extract_core_structure(self, safe):
        # Define scaffold parameter network
        params = rdScaffoldNetwork.ScaffoldNetworkParams()
        # customize parameter attributes
        params.includeScaffoldsWithoutAttachments=False
        if safe is not None:
            mol = Chem.MolFromSmiles(safe)
            net = rdScaffoldNetwork.CreateScaffoldNetwork([mol],params)
            nodemols = [Chem.MolFromSmiles(x) for x in net.nodes]

            filtered_list = []
            for mol in nodemols:
                # Check for the presence of attachment points and molecular weight range
                if "*" in Chem.MolToSmiles(mol) and self.lower_molecular_weight < Descriptors.MolWt(mol) < self.upper_molecular_weight:
                    filtered_list.append(mol)
            
            # If there are no scaffolds within the range, select the closest one
            if not filtered_list:
                closest_mol = min(nodemols, key=lambda x: abs(Descriptors.MolWt(x) - (self.lower_molecular_weight + self.upper_molecular_weight) / 2))
                filtered_list.append(closest_mol)
            
            # Sort the filtered list based on the number of heteroatoms (fewer carbons)
            filtered_list.sort(key=lambda x: x.GetNumHeavyAtoms())

            return filtered_list
        else:
            return None
        
    
    def _generate_smiles(self, side_chains):
        generated_smiles = self.designer.scaffold_morphing(
        side_chains=side_chains,
        n_samples_per_trial=self.n_samples_per_trial,
        n_trials=self.n_trials,
        sanitize=True,
        do_not_fragment_further=False,
        random_seed=100,)
        return generated_smiles


    def _get_side_chain_pairs(self, side_chains):
        '''Function to break the side chains into pairs'''
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


    def run_model(self, safe):
        generated_smiles = []
        for i in safe:
            row = []
            if i is not None:
                core_structures = self._extract_core_structure(i)
                for core in core_structures:
                    side_chain = compute_side_chains(core=core, mol=i)
                    side_chain_pairs = self._get_side_chain_pairs(side_chain)
                    for side_chain in side_chain_pairs:
                        output = self._generate_smiles(side_chain)
                        row += output
                generated_smiles += [row]
            else:
                generated_smiles += [row]
        return generated_smiles
