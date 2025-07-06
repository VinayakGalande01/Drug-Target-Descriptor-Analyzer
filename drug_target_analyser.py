from rdkit import Chem
from rdkit.Chem import Descriptors
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import logging

# Set up logging
logging.basicConfig(
    filename='analyser.log',
    filemode='a',
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO
)

# 1. compute_drug_descriptors(smiles)
# Returns MW, logP, HBD, HBA, TPSA, RotBonds

def compute_drug_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return [None] * 6

    return [
        Descriptors.ExactMolWt(mol),
        Descriptors.MolLogP(mol),
        Descriptors.NumHDonors(mol),
        Descriptors.NumHAcceptors(mol),
        Descriptors.NumRotatableBonds(mol),
        Descriptors.TPSA(mol)
    ]

# 2. compute_protein_features(sequence)
# Returns MW, aromaticity, instability_index, isoelectric_point
def compute_protein_features(seq):
    analysis = ProteinAnalysis(seq)
    return [
        analysis.molecular_weight(),
        analysis.aromaticity(),
        analysis.instability_index(),
        analysis.isoelectric_point(),
    ]


# 3. read_drugs_csv(path)
# Reads drug names + SMILES

def read_drugs_csv(path):
    df = pd.read_csv(path)
    return df['Name'].tolist(), df['SMILES'].tolist()

# 4. read_proteins_fasta(path)
# Reads protein names (IDs)+ sequences

def read_fasta_proteins(path):
    protein_names = []
    protein_seqs = []
    for record in SeqIO.parse(path, "fasta"):
        protein_names.append(record.id)
        protein_seqs.append(str(record.seq))
    return protein_names, protein_seqs


#5. compute_interaction_score(mw, pI)

def compute_interaction_score(mw, pI):
    if mw and pI:
        return round((mw * pI) / 1000, 2)
    return None

# Lipinski's Rule of Five (with rotatable bonds extension)
def lipinski_pass(mw, logp, hbd, hba, rot_bonds):
    return (
        mw is not None and mw <= 500 and
        logp is not None and logp <= 5 and
        hbd is not None and hbd <= 5 and
        hba is not None and hba <= 10 and
        rot_bonds is not None and rot_bonds <= 10
    )

def main():
    logging.info('Starting drug-target analysis')
    drugs_names, drugs_smiles = read_drugs_csv("drug.csv")
    proteins_names, proteins_seqs = read_fasta_proteins("proteins.fasta")
    
    results = []
    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
    
    for name, smiles in zip(drugs_names, drugs_smiles):
        drug_feats = compute_drug_descriptors(smiles)
        lipinski = lipinski_pass(*drug_feats[:5])
        
        for protein_id, seq in zip(proteins_names, proteins_seqs):
            # Skip empty or invalid sequences
            if not seq or not set(seq.upper()).issubset(valid_aa):
                logging.warning(f"Skipping protein {protein_id} due to invalid or empty sequence.")
                continue
            prot_feats = compute_protein_features(seq)
            score = compute_interaction_score(drug_feats[0], prot_feats[3])  # MW × pI

            results.append([
                name, smiles, *drug_feats, protein_id, len(seq), *prot_feats, score, lipinski
            ])

    columns = [
        "Drug_Name", "SMILES", "MolWeight", "LogP", "HBD", "HBA", "RotBonds", "TPSA",
        "Protein_ID", "Protein_Length", "Protein_MW", "Aromaticity", "Instability", "pI", "Score", "Lipinski_Pass"
    ]
    
    df = pd.DataFrame(results, columns=columns)
    df.to_csv("output.csv", index=False)
    logging.info("Output saved to output.csv")
    print("Output saved to output.csv")

if __name__ == "__main__":
    main()


