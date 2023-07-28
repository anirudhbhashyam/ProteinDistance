# Protein Distance
An implementation of the algorithm found in Chrysostomou, C., & Seker, H. (2013). Construction of protein distance matrix based on amino acid indices and Discrete Fourier Transform. Annual International Conference of the IEEE Engineering in Medicine and Biology Society. IEEE Engineering in Medicine and Biology Society. Annual International Conference, 2013, 4066–4069. https://doi.org/10.1109/EMBC.2013.6610438


# Usage
## Setup
It is recommended to use a virtual environment.
```bash
    pip install -r requirements.txt
```
## Run
```python
    python dc.py <path_to_fasta_file>
```
Currently the first two protein sequences from the file will be taken. And the distance will be output to `stdout`.

# TODO
- Add support for more than two sequences.
- *Proscale_4* feature is not present in the current implementation (which is used in the paper).
- Setup tests.