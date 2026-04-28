# Carbonara

Carbonara bridges the gap between crystal-like and solution-state conformations by efficiently refining protein structures using experimental SAXS (Small Angle X-ray Scattering) data. Starting from AI-predicted models or crystallographic structures, Carbonara rapidly explores conformational space to identify physiologically relevant solution-state conformations.

The method can incorporate additional experimental constraints such as disulfide bonds, NMR distance measurements, contact predictions, or FRET data to further guide the refinement process.

![Method Overview](figures/method_overview_arrows.png)

Schematic representation of the Carbonara refinement pipeline. The workflow proceeds from an initial structure (a), identification of flexible regions (b), conformational sampling guided by SAXS and optional structural constraints (c), model selection based on optimal fit (d), and optional all-atom reconstruction (e) for downstream applications.

---

## Which Carbonara workflow should I use?

Carbonara can be used at several levels. The simplest route runs the core C++ fitting algorithm. The most complete route provides guided setup, real-time all-atom reconstruction, SAXS scoring, and analysis.

| Workflow | Best for | What it includes | What it does not include |
|---|---|---|---|
| **1. Core C++ engine** | Reproducing prepared runs or running existing `RunMe_*.sh` scripts | Fast Carbonara C++ fitting engine | Guided setup, all-atom reconstruction, real-time analysis |
| **2. Basic setup / one-shot run** | Users with a PDB/mmCIF and SAXS curve who want Carbonara to prepare a run using sensible defaults | Python setup script, automatic `RunMe_*.sh` creation, optional one-shot fitting | Real-time all-atom scoring/monitoring unless using the full workflow |
| **3. Full interactive/all-atom workflow** | Complex systems, exploratory fitting, custom flexibility/constraints, multimers, mixtures, and production analysis | Full Python setup, notebooks/front-end, real-time monitoring, all-atom backmapping, pyFoXS scoring, analysis tools | MODELLER and CG2ALL are optional external tools and may require separate installation |

**For most new users**, workflow 2 is the quickest way to start a Carbonara calculation.

**For users who want the full Carbonara experience**, workflow 3 is recommended. It requires more setup, but gives much better control over flexible regions, constraints, multimeric rotations, mixture fitting, real-time all-atom predictions, and live analysis. Carbonara's strength is this flexibility: the user can provide as little or as much structural information as they want.

---

## Clone the repository

```bash
git clone https://github.com/Prior-Lab-Durham-University/carbonara.git carbonara
cd carbonara
```

If you are using a specific branch, for example the current WAXSiS/pyFoXS workflow branch:

```bash
git checkout pseudoWaxsis
```

---

## 1. Core C++ engine

Use this route if you already have prepared Carbonara input files and a `RunMe_*.sh` script.

This is the lightest route. It runs the core Carbonara fitting algorithm but does not automatically perform all-atom reconstruction or real-time pyFoXS scoring.

### Build the C++ algorithm

Carbonara requires CMake and a C++ compiler.

```bash
mkdir build
cd build
cmake ..
make
cd ..
```

The expected executable is:

```text
build/bin/predictStructureQvary
```

### Run a prepared refinement

```bash
sh RunMe_humanSMARCAL1.sh
```

or:

```bash
sh RunMe_C239S.sh
```

---

## 2. Basic setup for new structures

Use this route if you have a starting structure and SAXS data and want Carbonara to generate the required input files and `RunMe_*.sh` script for you.

You need:

1. a starting structure, usually a PDB or mmCIF file from AlphaFold, crystallography, or another source;
2. SAXS data in Å units with columns for `q`, intensity, and experimental error;
3. the C++ Carbonara binary built with CMake.

### Setup a new run

```bash
python setup_carbonara.py \
    --pdb path/to/model.pdb \
    --saxs path/to/saxs.dat \
    --name ProteinName
```

Then run:

```bash
sh RunMe_ProteinName.sh
```

### One-shot default run

If you trust the default settings and want setup and fitting to run automatically:

```bash
python run_carbonara_oneshot.py \
    --pdb path/to/model.pdb \
    --saxs path/to/saxs.dat \
    --name ProteinName
```

This route is useful when you want to get a standard Carbonara run started quickly. It does not attempt to give the same level of interactivity, real-time all-atom reconstruction, or live analysis as the full workflow.

### Multimers and rigid-body rotations

For multimeric systems or cases where domains/chains may move as rigid bodies:

```bash
python setup_carbonara.py \
    --pdb path/to/model.pdb \
    --saxs path/to/saxs.dat \
    --name ProteinName \
    --rotation
```

or:

```bash
python run_carbonara_oneshot.py \
    --pdb path/to/model.pdb \
    --saxs path/to/saxs.dat \
    --name ProteinName \
    --rotation
```

### Using AlphaFold PAE to guide flexibility

If you have a PAE file and want Carbonara to use AlphaFold uncertainty to select flexible regions:

```bash
python setup_carbonara.py \
    -p path/to/model.pdb \
    -s path/to/saxs.dat \
    -f path/to/pae.json \
    --name ProteinName \
    --alphaFoldFlex \
    --rotation
```

### Mixture fitting

If the molecule may occupy multiple states in solution, or if you suspect substantial variation in radius of gyration, Carbonara can fit mixtures of conformational states.

```bash
python setup_carbonara.py \
    -p path/to/model.pdb \
    -s path/to/saxs.dat \
    -f path/to/pae.json \
    --name ProteinName \
    --alphaFoldFlex \
    --rotation \
    --mixture_n 2 \
    --max_mixture_combos 10
```

Useful options:

```text
--fit_n_times INT          Number of independent fitting runs, default 20
--min_q FLOAT              Minimum q-value, default 0.01
--max_q FLOAT              Maximum q-value, default 0.2
--max_q_start FLOAT        Maximum q-value used at the start of fitting
--max_fit_steps INT        Maximum number of fitting steps, default 10000
--pairedQ                  Use paired-distance/contact-style constraints
--rotation                 Apply affine rigid-body rotations
--alphaFoldFlex            Use PAE to specify flexible regions
--pae_flex_threshold       PAE threshold for selecting flexible regions, default 16
--mixture_n INT            Number of structures/species in a mixture refinement
--max_mixture_combos INT   Number of mixture combinations to test
```

---

## 3. Full interactive / all-atom workflow

Use this route for the most complete Carbonara workflow.

The full workflow is designed for users who want to interactively control the model and analyse results as they are produced. It supports:

- guided selection of flexible regions;
- user-defined flexible and fixed sections;
- multimers and rigid-body rotations;
- distance constraints such as disulfides, contact predictions, NMR-like distances, crosslinks, or FRET-style measurements;
- mixture fitting and ensemble-style comparisons;
- real-time monitoring of Carbonara predictions;
- real-time all-atom reconstruction;
- pyFoXS-based all-atom SAXS scoring;
- downstream analysis of χ², RMSD, TM-score, GDT-TS, radius of gyration, and model quality.

This is the recommended route for complex systems or production analysis.

### Install the Python workflow

For the full interactive/all-atom workflow, run:

```bash
bash setupPython.sh
```

The setup script installs the Python/Git-installable dependencies needed by the notebook/front-end tools. It also handles packages that need special installation order, such as `biobox`, and installs pyFoXS runtime dependencies.

We recommend running this inside a clean Python environment if possible, but this is not strictly required if you already manage Python packages another way.

### Optional: use a local virtual environment

If you want to avoid modifying your existing Python environment, create a virtual environment first:

```bash
python3 -m venv .venv
source .venv/bin/activate
bash setupPython.sh
```

For `tcsh`/`csh` shells:

```tcsh
python3 -m venv .venv
source .venv/bin/activate.csh
bash setupPython.sh
```

This is also a useful way to test that installation works from scratch.

### Build the C++ engine

If you have not already built the C++ engine:

```bash
mkdir build
cd build
cmake ..
make
cd ..
```

### Using notebooks

If using Jupyter, register your Python environment as a kernel:

```bash
python -m ipykernel install --user \
    --name carbonara \
    --display-name "Python (Carbonara)"
```

Inside notebooks, prefer:

```python
import sys
!{sys.executable} setup_carbonara_allAtom.py --help
```

rather than:

```python
!python setup_carbonara_allAtom.py --help
```

This ensures that notebook subprocesses use the same Python environment as the active kernel.

---

## All-atom reconstruction and optional external tools

The full workflow can generate all-atom models from Carbonara Cα predictions. This is optional and requires additional tooling.

### pyFoXS

`setupPython.sh` installs the Python dependencies needed by pyFoXS and creates a `pyfoxs` wrapper where possible. This enables all-atom SAXS scoring during or after a Carbonara run.

### MODELLER

MODELLER is optional and is **not** installed by `setupPython.sh`, because it requires a separate licence.

To check whether MODELLER is available in your current Python environment:

```bash
python -c "from modeller import environ; env=environ(); print('MODELLER OK')"
```

If MODELLER-based backmapping is selected but MODELLER is not installed or licensed, Carbonara should fail early with a clear message. This is expected behaviour.

The basic C++ Carbonara refinement does not require MODELLER.

### CG2ALL

CG2ALL can also be used for all-atom reconstruction, but it is usually installed separately in its own environment. It is not required for the basic Carbonara run.

---

## Constraints and flexibility

Carbonara is designed to allow as much or as little prior structural knowledge as the user wants to provide.

For simple cases, the setup scripts can choose reasonable default flexible regions. For complex cases, users can specify flexible regions, fixed regions, multimeric rigid-body rotations, mixture states, and structural restraints.

Examples of useful constraints include:

- disulfide bonds;
- contact predictions;
- NMR-style distance measurements;
- crosslinking data;
- FRET-style distances;
- known rigid domains;
- user-defined flexible linkers.

The one-shot route is designed to get a reasonable calculation running quickly. The full workflow is designed for careful system-specific modelling.

---

## Reproducing structures refined in the paper

To reproduce the refinement of the two structures presented in the paper, first ensure you are located in `/path/to/carbonara`, build the C++ code, and then run:

### human SMARCAL1

```bash
sh RunMe_humanSMARCAL1.sh
```

### ChiLob7/4 IgG2

```bash
sh RunMe_C239S.sh
```

---

## Colab implementation

Carbonara's key strength is its flexibility: users can specify as little or as much of the structure to vary, enforce rigid-body motions of subdomains, and apply a wide range of distance constraints.

The Colab implementation provides a guided graphical setup for choosing flexible regions and preparing Carbonara runs.

- Monomer version
- Multimer version

> ⚠️ These notebooks are shared in view-only mode. To use them, click **Copy to Drive** at the top of the Colab page. This creates your own editable copy in Google Drive.

---

## Troubleshooting

### `!python script.py` fails in a notebook, but imports work in the notebook

Use:

```python
import sys
!{sys.executable} script.py
```

This ensures that the subprocess uses the same Python environment as the notebook kernel.

### `biobox` fails to install

Use `setupPython.sh` rather than installing packages manually. `biobox` requires NumPy to be available during its build step, and the setup script handles this ordering.

### pyFoXS complains about missing `torch`

Use `setupPython.sh`. The full setup installs the pyFoXS runtime dependencies.

### MODELLER is missing

This is expected unless you installed and licensed MODELLER separately. MODELLER is only needed for MODELLER-based all-atom backmapping.

### The basic Carbonara engine runs but all-atom reconstruction does not

This usually means the core C++ run is correctly installed, but the optional all-atom backmapping route is missing one of its external tools. Check whether you selected MODELLER or CG2ALL, and confirm that the selected tool is installed and available from the same Python environment used by the notebook or script.

---

## Citation

If you use Carbonara in your research, please cite our preprint:

```bibtex
@article{carbonara2025,
  title={Carbonara: A Rapid Method for SAXS-Based Refinement of Protein Structures},
  author={McKeown, J. and Bale, A. and Brown, C. and Fisher, H. and Rambo, R. and Essex, J. and Degiacomi, M. and Prior, C.},
  journal={ResearchSquare},
  year={2025},
  doi={10.21203/rs.3.rs-6447099/v1},
  url={https://doi.org/10.21203/rs.3.rs-6447099/v1}
}
```

![CC BY-NC-SA 4.0](https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png)

This work is licensed under a [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License](https://creativecommons.org/licenses/by-nc-sa/4.0/).
