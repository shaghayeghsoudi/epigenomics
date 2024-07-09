

## Assuming you have a file with methylation calls for each cytosine, here’s a script to calculate conversion efficiency:
pip install pandas



import pandas as pd

# Load methylation data
# Example format: Chromosome, Position, Context, Strand, Methylated, Total
# Context could be "CHG", "CHH", "CpG"
methylation_data = pd.read_csv('methylation_calls.txt', sep='\t')

# Filter for non-CpG sites
non_cpg_data = methylation_data[(methylation_data['Context'] == 'CHG') | (methylation_data['Context'] == 'CHH')]

# Calculate total non-CpG sites and non-converted (methylated) non-CpG sites
total_non_cpg = non_cpg_data['Total'].sum()
non_converted_non_cpg = non_cpg_data['Methylated'].sum()

# Bisulfite conversion efficiency
conversion_efficiency = 1 - (non_converted_non_cpg / total_non_cpg)

print(f'Bisulfite conversion efficiency: {conversion_efficiency:.4f}')