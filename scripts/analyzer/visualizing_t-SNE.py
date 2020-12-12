import pandas as pd
import numpy as np
from sklearn.manifold import TSNE
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

# Perform t-SNE
fit_tsne = TSNE(n_components=2, random_state=0).fit_transform(all_data)

# Plot the result of t-SNE
fig = plt.figure()
for label, name in zip(all_data_label.unique(), ["hg", "virus", "EVE"]):
    x = fit_tsne[all_data_label == label, 0]
    y = fit_tsne[all_data_label == label, 1]

    plt.scatter(x, y, s=5, label=name)
plt.xlabel("t-SNE1")
plt.ylabel("t-SNE2")

# Save the plot
fig = plt.gcf()
fig.savefig("plot_t-SNE.pdf", format="pdf")
plt.close()
