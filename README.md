# Scaffold Morphing

The context discusses a novel notation system called Sequential Attachment-based Fragment Embedding (SAFE) that improves upon traditional molecular string representations like SMILES. SAFE reframes SMILES strings as an unordered sequence of interconnected fragment blocks while maintaining compatibility with existing SMILES parsers. This streamlines complex molecular design tasks by facilitating autoregressive generation under various constraints. The effectiveness of SAFE is demonstrated by training a GPT2-like model on a dataset of 1.1 billion SAFE representations that exhibited versatile and robust optimization performance for molecular design.

## Identifiers

* EOS model ID: `eos8bhe`
* Slug: `scaffold-morphing`

## Characteristics

* Input: `Compound`
* Input Shape: `Single`
* Task: `Generative`
* Output: `Compound`
* Output Type: `String`
* Output Shape: `List`
* Interpretation: Model generates new molecules from input molecule by replacing core structures of input molecule.

## References

* [Publication](https://arxiv.org/pdf/2310.10773.pdf)
* [Source Code](https://github.com/datamol-io/safe/tree/main)
* Ersilia contributor: [Inyrkz](https://github.com/Inyrkz)

## Ersilia model URLs
* [GitHub](https://github.com/ersilia-os/eos8bhe)
* [AWS S3](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos8bhe.zip)

## Citation

If you use this model, please cite the [original authors](https://arxiv.org/pdf/2310.10773.pdf) of the model and the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff).

## License

This package is licensed under a GPL-3.0 license. The model contained within this package is licensed under a CC license.

Notice: Ersilia grants access to these models 'as is' provided by the original authors, please refer to the original code repository and/or publication if you use the model in your research.

## About Us

The [Ersilia Open Source Initiative](https://ersilia.io) is a Non Profit Organization ([1192266](https://register-of-charities.charitycommission.gov.uk/charity-search/-/charity-details/5170657/full-print)) with the mission is to equip labs, universities and clinics in LMIC with AI/ML tools for infectious disease research.

[Help us](https://www.ersilia.io/donate) achieve our mission!