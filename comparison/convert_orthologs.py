import pandas as pd
import mygene
import os

def convert_mouse_to_human_csv(input_file, output_file):
    print(f"Reading file: {input_file}")
    # Read CSV file
    mouse_degs_df = pd.read_csv(input_file)
    
    # Extract 'gene' column, remove duplicates and NaNs
    mouse_gene_list = mouse_degs_df['gene'].dropna().unique().tolist()
    print(f"Querying human orthologs for {len(mouse_gene_list)} mouse genes...")
    
    mg = mygene.MyGeneInfo()
    # Query human orthologs using mygene
    human_orthologs_info = mg.querymany(
        mouse_gene_list,
        scopes='symbol',
        species='mouse',
        fields='homologene',
        verbose=False
    )
    
    print("Query complete, processing results...")
    ortholog_map = {}
    human_entrez_ids = []
    
    # Map mouse genes to human Entrez IDs
    for gene_info in human_orthologs_info:
        query_gene = gene_info['query']
        # Ensure result contains homologene info
        if 'homologene' in gene_info and 'genes' in gene_info['homologene']:
            # Iterate through homologene list (TaxID, EntrezID)
            for tax_id, entrez_id in gene_info['homologene']['genes']:
                if tax_id == 9606:  # Homo sapiens TaxID
                    ortholog_map[query_gene] = str(entrez_id)
                    human_entrez_ids.append(str(entrez_id))
                    break  # Stop once human ortholog is found
    
    print(f"Found {len(human_entrez_ids)} ortholog Entrez IDs, converting back to symbols...")
    
    # Convert Entrez IDs back to gene symbols
    if human_entrez_ids:
        human_gene_symbols_info = mg.getgenes(list(set(human_entrez_ids)), fields='symbol', as_dataframe=True)
        id_to_symbol_map = human_gene_symbols_info['symbol'].to_dict()
    else:
        id_to_symbol_map = {}

    # Final map: Mouse Symbol -> Human Symbol
    final_ortholog_map = {
        mouse_gene: id_to_symbol_map.get(human_id)
        for mouse_gene, human_id in ortholog_map.items()
        if id_to_symbol_map.get(human_id) is not None
    }
    
    # Create new column 'human_ortholog'
    mouse_degs_df['human_ortholog'] = mouse_degs_df['gene'].map(final_ortholog_map)
    
    # Remove rows without human orthologs
    converted_df = mouse_degs_df.dropna(subset=['human_ortholog'])
    
    # Save results
    converted_df.to_csv(output_file, index=False)
    print(f"Successfully converted {len(converted_df)} rows. Saved to: {output_file}")

# --- Execution ---
# Define input and output filenames
# Note: Ensure these paths exist on your machine
input_csv = os.path.expanduser("~/volcano_402.csv")  ##volcano_403.csv volcano_404.csv volcano_405.csv
output_csv = os.path.expanduser("~/402.csv") ##403.csv 404.csv 405.csv 

# Run conversion
convert_mouse_to_human_csv(input_csv, output_csv)
