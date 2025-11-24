import os
base_dir='/home/nova/Desktop/alen-belen/results'
barcode_dirs=os.path.join(base_dir, 'genomad_bakta')
peptides_dir=os.path.join(base_dir,'bakta')

def extract_ids_from_bed (bed_file):
    ids=set()
    with open(bed_file, 'r') as bed:
        for line in bed:
            columns = line.strip().split('\t')
            if len(columns) >= 4:
                ids.add(columns[3])
    return ids

def extract_sequences_from_fasta(fasta_file,ids):
    sequences={}
    with open(fasta_file, 'r') as fasta:
        current_id= None
        current_seq=[]
        for line in fasta:
            if line.startswith('>'):
                if current_id and current_id in ids:
                    sequences[current_id]=''.join(current_seq)
                parts=line.split()
                current_id=parts[0][1:]
                current_seq=[]
            else:
                current_seq.append(line.strip()) 
        if current_seq and current_id in ids:
            sequences[current_id]=''.join(current_seq)
    return sequences

for barcode in os.listdir(barcode_dirs):
    barcode_path = os.path.join(barcode_dirs, barcode)
    if os.path.isdir(barcode_path):
        bed_file = os.path.join(barcode_path, 'virus', f'{barcode}_virus_genes_bakta.bed')
        fasta_file = os.path.join(peptides_dir, barcode, 'consensus.faa')

        if os.path.isfile(bed_file) and os.path.isfile(fasta_file):
            # Extract IDs from BED
            ids = extract_ids_from_bed(bed_file)
            # Extract matching sequences from FASTA
            sequences = extract_sequences_from_fasta(fasta_file, ids)

            # Write the extracted sequences to a new FASTA file
            output_file = os.path.join(barcode_path,'virus', f'{barcode}_extracted_sequences.faa')
            with open(output_file, 'w') as out_fasta:
                for seq_id, seq in sequences.items():
                    out_fasta.write(f'>{seq_id}\n{seq}\n')

            print(f'Extracted sequences for {barcode} written to {output_file}')



