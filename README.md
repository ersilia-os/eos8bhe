# Scaffold morphing

The context discusses a novel notation system called Sequential Attachment-based Fragment Embedding (SAFE) that improves upon traditional molecular string representations like SMILES. SAFE reframes SMILES strings as an unordered sequence of interconnected fragment blocks while maintaining compatibility with existing SMILES parsers. This streamlines complex molecular design tasks by facilitating autoregressive generation under various constraints. The effectiveness of SAFE is demonstrated by training a GPT2-like model on a dataset of 1.1 billion SAFE representations that exhibited versatile and robust optimization performance for molecular design.

This model was incorporated on 2024-01-12.

## Information
### Identifiers
- **Ersilia Identifier:** `eos8bhe`
- **Slug:** `scaffold-morphing`

### Domain
- **Task:** `Sampling`
- **Subtask:** `Generation`
- **Biomedical Area:** `Any`
- **Target Organism:** `Not Applicable`
- **Tags:** `Compound generation`

### Input
- **Input:** `Compound`
- **Input Dimension:** `1`

### Output
- **Output Dimension:** `1000`
- **Output Consistency:** `Variable`
- **Interpretation:** Model generates new molecules from input molecule by replacing core structures of input molecule.

Below are the **Output Columns** of the model:
| Name | Type | Direction | Description |
|------|------|-----------|-------------|
| smiles_000 | string |  | Generated molecule index 0 using the SAFE scaffold morphing method |
| smiles_001 | string |  | Generated molecule index 1 using the SAFE scaffold morphing method |
| smiles_002 | string |  | Generated molecule index 2 using the SAFE scaffold morphing method |
| smiles_003 | string |  | Generated molecule index 3 using the SAFE scaffold morphing method |
| smiles_004 | string |  | Generated molecule index 4 using the SAFE scaffold morphing method |
| smiles_005 | string |  | Generated molecule index 5 using the SAFE scaffold morphing method |
| smiles_006 | string |  | Generated molecule index 6 using the SAFE scaffold morphing method |
| smiles_007 | string |  | Generated molecule index 7 using the SAFE scaffold morphing method |
| smiles_008 | string |  | Generated molecule index 8 using the SAFE scaffold morphing method |
| smiles_009 | string |  | Generated molecule index 9 using the SAFE scaffold morphing method |

_10 of 1000 columns are shown_
### Source and Deployment
- **Source:** `Local`
- **Source Type:** `External`
- **DockerHub**: [https://hub.docker.com/r/ersiliaos/eos8bhe](https://hub.docker.com/r/ersiliaos/eos8bhe)
- **Docker Architecture:** `AMD64`
- **S3 Storage**: [https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos8bhe.zip](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos8bhe.zip)

### Resource Consumption


### References
- **Source Code**: [https://github.com/datamol-io/safe/tree/main](https://github.com/datamol-io/safe/tree/main)
- **Publication**: [https://arxiv.org/pdf/2310.10773.pdf](https://arxiv.org/pdf/2310.10773.pdf)
- **Publication Type:** `Preprint`
- **Publication Year:** `2023`
- **Ersilia Contributor:** [Inyrkz](https://github.com/Inyrkz)

### License
This package is licensed under a [GPL-3.0](https://github.com/ersilia-os/ersilia/blob/master/LICENSE) license. The model contained within this package is licensed under a [Apache-2.0](LICENSE) license.

**Notice**: Ersilia grants access to models _as is_, directly from the original authors, please refer to the original code repository and/or publication if you use the model in your research.


## Use
To use this model locally, you need to have the [Ersilia CLI](https://github.com/ersilia-os/ersilia) installed.
The model can be **fetched** using the following command:
```bash
# fetch model from the Ersilia Model Hub
ersilia fetch eos8bhe
```
Then, you can **serve**, **run** and **close** the model as follows:
```bash
# serve the model
ersilia serve eos8bhe
# generate an example file
ersilia example -n 3 -f my_input.csv
# run the model
ersilia run -i my_input.csv -o my_output.csv
# close the model
ersilia close
```

## About Ersilia
The [Ersilia Open Source Initiative](https://ersilia.io) is a tech non-profit organization fueling sustainable research in the Global South.
Please [cite](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) the Ersilia Model Hub if you've found this model to be useful. Always [let us know](https://github.com/ersilia-os/ersilia/issues) if you experience any issues while trying to run it.
If you want to contribute to our mission, consider [donating](https://www.ersilia.io/donate) to Ersilia!
