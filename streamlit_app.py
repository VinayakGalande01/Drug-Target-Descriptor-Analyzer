import streamlit as st
import subprocess
import os
import pandas as pd

st.title('Drug–Target Descriptor Analyzer')

# Paste or upload drug.csv
drug_text = st.text_area('Paste drug.csv content here (Name,SMILES)', height=150)
drug_file = st.file_uploader('Or upload drug.csv', type='csv')

# Paste or upload proteins.fasta
protein_text = st.text_area('Paste proteins.fasta content here (FASTA format)', height=150)
protein_file = st.file_uploader('Or upload proteins.fasta', type=['fasta', 'fa'])

# Save uploaded or pasted files
def save_uploadedfile(uploadedfile, filename):
    with open(filename, 'wb') as f:
        f.write(uploadedfile.getbuffer())

def save_textfile(text, filename):
    with open(filename, 'w', encoding='utf-8') as f:
        f.write(text)

# Determine which input to use and save
drug_ready = False
protein_ready = False

if drug_text.strip():
    save_textfile(drug_text, 'drug.csv')
    drug_ready = True
    st.success('Drug data pasted and saved!')
elif drug_file:
    save_uploadedfile(drug_file, 'drug.csv')
    drug_ready = True
    st.success('Drug file uploaded and saved!')

if protein_text.strip():
    save_textfile(protein_text, 'proteins.fasta')
    protein_ready = True
    st.success('Protein data pasted and saved!')
elif protein_file:
    save_uploadedfile(protein_file, 'proteins.fasta')
    protein_ready = True
    st.success('Protein file uploaded and saved!')

if drug_ready and protein_ready:
    if st.button('Run Analysis'):
        with st.spinner('Running analysis...'):
            result = subprocess.run(['python', 'drug_target_analyser.py'], capture_output=True, text=True)
            if result.returncode == 0:
                st.success('Analysis complete!')
                if os.path.exists('output.csv'):
                    df = pd.read_csv('output.csv')
                    # Add visual indicator for Lipinski_Pass
                    def lipinski_icon(val):
                        if val is True or val == True or str(val).lower() == 'true':
                            return '✅'
                        else:
                            return '❌'
                    df['Lipinski_Visual'] = df['Lipinski_Pass'].apply(lipinski_icon)
                    # Move visual column next to Lipinski_Pass
                    cols = list(df.columns)
                    cols.insert(cols.index('Lipinski_Pass') + 1, cols.pop(cols.index('Lipinski_Visual')))
                    df = df[cols]
                    st.dataframe(df)
                    st.download_button('Download output.csv', df.to_csv(index=False), 'output.csv')
                else:
                    st.error('output.csv not found!')
            else:
                st.error('Error running analysis!')
                st.text(result.stderr)

# Show log file
if os.path.exists('analyser.log'):
    st.subheader('Log File')
    with open('analyser.log', 'r') as logf:
        st.text(logf.read())
