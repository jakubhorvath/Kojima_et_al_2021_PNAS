import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(color_codes=True)
sns.set_style("white")

# Load dataset
hg_data = pd.read_csv("human_kmer_freq.csv", index_col=0)
virus_data = pd.read_csv("virus_kmer_freq.csv", index_col=0)
eve_data = pd.read_csv("EVE_kmer_freq.csv", index_col=0)

# Store all_data with merged data
all_data = pd.concat([hg_data, virus_data, eve_data])
all_data_label = pd.Series(np.concatenate([np.zeros(len(hg_data)),
                                           np.ones(len(virus_data)),
                                           np.full(len(eve_data), 2)],
                                          axis=0)
                           ).astype(int)

# Hierarchical clustering and plot the results.
fig = plt.figure(figsize=(50, 50))
lut = dict(zip(all_data_label.unique(), "bry"))
row_colors = all_data_label.map(lut)
row_colors.index = all_data.index
g = sns.clustermap(all_data,
                   figsize=(50, 50),
                   method="ward",
                   metric="euclidean",
                   row_colors=row_colors,
                   row_cluster=True,
                   col_cluster=True,
                   cmap="Oranges")

for tick_label, color in zip(g.ax_heatmap.axes.get_yticklabels(), row_colors.__iter__()):
    tick_label.set_color(color)

# Save the plot
fig = plt.gcf()
fig.savefig("plot_hierarchical_clustering.pdf", format="pdf")
plt.close()
