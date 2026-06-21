# ESM-2 Protein Embeddings and Similarity Analysis

This project demonstrates how to extract sequence-level and residue-level embeddings for a list of proteins using the ESM-2 pre-trained transformer model (`facebook/esm2_t6_8M_UR50D`), perform dimensionality reduction (UMAP), and analyze similarity patterns within and between protein families.

## Features
- **Embeddings Extraction**: Uses the smallest ESM-2 model (8M parameters) to embed protein sequences.
- **UMAP Dimensionality Reduction**: Projects high-dimensional protein embeddings (320 dimensions) into 2D space to visualize groupings.
- **Similarity Analysis**: Evaluates cosine similarity of sequence representations within same-family vs cross-family protein pairs, verified by a Mann-Whitney U test.

## Required Inputs and Outputs

### Inputs
1. `proteins.fasta`: FASTA file containing 45 human protein sequences across different families (kinases, GPCRs, immunoglobulins, oxygen binding proteins, etc.).
2. `protein_accessions.tsv`: Metadata file linking protein accession IDs to their respective protein family names.
3. Pre-trained Model: `facebook/esm2_t6_8M_UR50D` (downloaded dynamically from Hugging Face).

### Outputs
1. Interactive plots generated inline within the Jupyter Notebook:
   - **UMAP Scatter Plot**: Projects sequence embeddings colored by family.
   - **Cosine Similarity Distribution Plot**: Compares similarity of same-family pairs vs different-family pairs.
2. Cell output logs detailing loaded datasets, embedding shapes, average similarities, and Mann-Whitney U statistical test outcomes showing clear separation.

## Dependencies

The project relies on several key scientific and machine learning libraries:
- **PyTorch (`torch`)** - Deep learning framework.
- **Hugging Face Transformers (`transformers`)** - ESM-2 model and tokenizer interface.
- **Biopython (`biopython`)** - FASTA sequence parsing.
- **scikit-learn (`scikit-learn`)** - Pairwise cosine similarity calculations.
- **UMAP-learn (`umap-learn`)** - Non-linear dimensionality reduction.
- **SciPy (`scipy`)** - Mann-Whitney U statistical test.
- **Pandas (`pandas`)**, **NumPy (`numpy`)** - Data structures and numeric operations.
- **Matplotlib (`matplotlib`)**, **Seaborn (`seaborn`)** - Plotting and visualizations.
- **ipykernel** - Support for running cells in Jupyter notebook.

These dependencies are declared directly inside the [pyproject.toml](file:///c:/Users/leona/Documents/Summer%20School/bioinfo-school/Week_4_reproducible/pyproject.toml) configuration.

## How to Rerun

1. **Set up environment and install dependencies**:
   ```bash
   # Using uv (fastest)
   uv venv .venv
   # Windows:
   .venv\Scripts\activate.ps1
   # Linux/macOS:
   source .venv/bin/activate

   # Install the package and all its dependencies from pyproject.toml
   uv pip install .
   ```

   *Alternatively, using standard pip:*
   ```bash
   python -m venv .venv
   # Activate environment... (e.g., .venv\Scripts\activate)
   # Install the package and all its dependencies
   pip install .
   ```

2. **Open and Run the Jupyter Notebook**:
   - Start the Jupyter server:
     ```bash
     jupyter notebook
     ```
   - Open `B_protein_embeddings_esm2.ipynb`.
   - Select the `.venv` kernel (or the Python interpreter associated with the virtual environment).
   - Run all cells sequentially.

