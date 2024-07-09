import pandas as pd

# Load methylation data
# Example format: Chromosome, Position, Context, Strand, Methylated, Total
# Context could be "CHG", "CHH", "CpG"
methylation_data = pd.read_csv('methylation_calls.txt', sep='\t')

# Assume we have a column 'Sample_Type' that identifies spike-in control
# For simplicity, let's assume 'Sample_Type' is 'spike_in' for spike-in controls

# Filter for spike-in control data
spike_in_data = methylation_data[methylation_data['Sample_Type'] == 'spike_in']

# Filter for non-CpG sites (CHG and CHH contexts)
non_cpg_spike_in_data = spike_in_data[(spike_in_data['Context'] == 'CHG') | (spike_in_data['Context'] == 'CHH')]

# Calculate total non-CpG sites and non-converted (methylated) non-CpG sites
total_non_cpg_spike_in = non_cpg_spike_in_data['Total'].sum()
non_converted_non_cpg_spike_in = non_cpg_spike_in_data['Methylated'].sum()

# Bisulfite conversion efficiency
conversion_efficiency = 1 - (non_converted_non_cpg_spike_in / total_non_cpg_spike_in)

print(f'Bisulfite conversion efficiency: {conversion_efficiency:.4f}')
