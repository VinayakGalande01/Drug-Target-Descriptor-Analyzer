# Drug-Target-Descriptor-Analyzer
 Analyze drug–protein interactions with ease! Paste or upload drug SMILES and protein FASTA, calculate key descriptors, check Lipinski’s drug-likeness, and get visual results—all in a user-friendly Streamlit web app. Perfect for drug discovery research.

A user-friendly tool for analyzing drug–target interactions and evaluating drug-likeness using Lipinski's Rule of Five. This project combines cheminformatics and bioinformatics to help researchers, students, and professionals quickly assess drug candidates and protein targets.

---

## Features
- **Paste or upload** drug SMILES and protein FASTA data
- **Automatic descriptor calculation** for drugs (molecular weight, LogP, H-bond donors/acceptors, rotatable bonds, TPSA)
- **Protein feature extraction** (molecular weight, aromaticity, instability, isoelectric point)
- **Lipinski's Rule of Five** evaluation with visual pass/fail indicators
- **Downloadable results** as CSV
- **View logs** for transparency and troubleshooting
- **Modern Streamlit web interface**

---

## Lipinski's Rule of Five (Drug-Likeness)
A drug is considered "drug-like" if it meets these criteria:
| Property                   | Threshold  |
|---------------------------|------------|
| Molecular Weight          | ≤ 500 Da   |
| LogP (lipophilicity)      | ≤ 5        |
| H-Bond Donors (HBD)       | ≤ 5        |
| H-Bond Acceptors (HBA)    | ≤ 10       |
| Rotatable Bonds           | ≤ 10       |

The app shows a green tick (✅) if a drug passes all rules, or a red cross (❌) if it fails any.

---

## Setup Instructions

### 1. Clone the repository and enter the project folder
```bash
# Example
cd "C:/Users/YourName/Desktop/Drug–Target Descriptor Analyzer"
```

### 2. Create and activate a virtual environment (recommended)
```bash
python -m venv drug_tar
# Activate (Windows)
drug_tar\Scripts\Activate.ps1
```

### 3. Install dependencies
```bash
pip install -r requirements.txt
```

---

## Usage

### **Run the Streamlit App**
```bash
streamlit run streamlit_app.py
```

### **Web Interface**
- **Paste or upload** your `drug.csv` (with columns `Name,SMILES`) and `proteins.fasta` (FASTA format)
- Click **Run Analysis**
- View/download results and see Lipinski's pass/fail for each drug
- View logs for troubleshooting

#### **Example drug.csv**
```
Name,SMILES
Aspirin,CC(=O)OC1=CC=CC=C1C(=O)O
Ibuprofen,CC(C)CC1=CC=C(C=C1)C(C)C(=O)O
```

#### **Example proteins.fasta**
```
>Protein_1
MVKVYAPASSANMSVGFDVLGAAVTPVDGALLGDVVTVEAAETFSLNNLGQKL
>Protein_2
GAVLGGEEAKKVAAWTLKAAAGGQFPLGADQLHFDVK
```

---

## Output
- Results are shown in the app and saved as `output.csv`
- The `Lipinski_Pass` column is `True`/`False`, and the next column shows a green tick (✅) or red cross (❌)
- All steps and warnings are logged in `analyser.log`

---

## Troubleshooting
- **ModuleNotFoundError:** Make sure you activate your virtual environment before running Streamlit.
- **KeyError: 'Name':** Ensure your `drug.csv` has the correct header: `Name,SMILES`.
- **Other issues:** Check the `analyser.log` file for details.

---

## License
This project is for educational and research use. See `LICENSE` for details. 
