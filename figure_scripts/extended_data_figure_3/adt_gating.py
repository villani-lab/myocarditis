import os
import pegasus as pg
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.patches as patches
import matplotlib.pyplot as plt


def plot_scatterplot(ax, df, x_gene, y_gene, patches=[]):
    x = np.squeeze(np.asarray(df[:, x_gene].copy().get_matrix('X').todense()))
    y = np.squeeze(np.asarray(df[:, y_gene].copy().get_matrix('X').todense()))
    col = df.obs['lineage'].values
    plot_df = pd.DataFrame({'x': x + 1, 'y': y + 1, 'lineage': col})

    sns.scatterplot(plot_df, x='x', y='y', hue='lineage', edgecolor="none", s=5, ax=ax)
    ax.set_xlim([0.8, None])
    ax.set_ylim([0.8, None])
    ax.set_xlabel(x_gene)
    ax.set_ylabel(y_gene)
    ax.get_legend().remove()
    for patch in patches:
        ax.add_patch(patch)


def plot_hexbin(ax, df, x_gene, y_gene, patches=[], grid_size=80):
    x = np.squeeze(np.asarray(df[:, x_gene].copy().get_matrix('X').todense()))
    y = np.squeeze(np.asarray(df[:, y_gene].copy().get_matrix('X').todense()))
    plot_df = pd.DataFrame({'x': x + 1, 'y': y + 1})

    ax.hexbin(plot_df['x'],
              plot_df['y'],
              mincnt=1,
              gridsize=grid_size,
              edgecolor=None,
              bins='log')
    ax.set_xlabel(x_gene)
    ax.set_ylabel(y_gene)
    ax.set_xlim([0.8, None])
    ax.set_ylim([0.8, None])
    for patch in patches:
        ax.add_patch(patch)


def filter_df(df, x_gene, y_gene, x_line=None, y_line=None, x_gr=True, y_gr=True):
    # make dataframe with gene counts
    x = np.squeeze(np.asarray(df[:, x_gene].copy().get_matrix('X').todense()))
    y = np.squeeze(np.asarray(df[:, y_gene].copy().get_matrix('X').todense()))
    plot_df = pd.DataFrame({'x': x + 1, 'y': y + 1})

    # get info whether passed threshold
    x_mask = plot_df['x'] > x_line if x_gr else plot_df['x'] < x_line
    y_mask = plot_df['y'] > y_line if y_gr else plot_df['y'] < y_line

    return df[x_mask & y_mask].copy()


# read in data
os.chdir('/projects/home/ikernin/projects/myocarditis/updated_datasets')
adt_df = pg.read_input("adt_dataset.zarr")

# get all filtered dataframes to plot all at once
b_df = filter_df(adt_df,
                 x_gene='cite_CD3E',
                 y_gene='cite_CD19',
                 x_line=1.8,
                 y_line=1.7,
                 x_gr=False,
                 y_gr=True)
t_df = filter_df(adt_df,
                 x_gene='cite_CD3E',
                 y_gene='cite_CD19',
                 x_line=2,
                 y_line=1.7,
                 x_gr=True,
                 y_gr=False)
cd8_df = filter_df(t_df,
                   x_gene='cite_CD4',
                   y_gene='cite_CD8A',
                   x_line=2,
                   y_line=1.9,
                   x_gr=False,
                   y_gr=True)
cd4_df = filter_df(t_df,
                   x_gene='cite_CD4',
                   y_gene='cite_CD8A',
                   x_line=2,
                   y_line=1.7,
                   x_gr=True,
                   y_gr=False)
gate_b_df = filter_df(adt_df,
                      x_gene='cite_CD3E',
                      y_gene='cite_CD19',
                      x_line=1.4,
                      y_line=1.5,
                      x_gr=False,
                      y_gr=False)
nk_df = filter_df(gate_b_df,
                  x_gene='cite_CD14',
                  y_gene='cite_NCAM1',
                  x_line=1.4,
                  y_line=1.6,
                  x_gr=False,
                  y_gr=True)
gate_d_df = filter_df(gate_b_df,
                      x_gene='cite_CD14',
                      y_gene='cite_NCAM1',
                      x_line=2.8,
                      y_line=1.5,
                      x_gr=False,
                      y_gr=False)
pdc_df = filter_df(gate_d_df,
                   x_gene='cite_HLA-DRB1',
                   y_gene='cite_IL3RA',
                   x_line=1.5,
                   y_line=2.8,
                   x_gr=True,
                   y_gr=True)
gate_e_df = filter_df(gate_d_df,
                      x_gene='cite_HLA-DRB1',
                      y_gene='cite_IL3RA',
                      x_line=5,
                      x_gr=False,
                      y_line=2.0,
                      y_gr=False
                      )
e_df = filter_df(gate_e_df,
                 x_gene='cite_ITGAX',
                 y_gene='cite_CD1C',
                 x_line=7,  # one very large outlier distorting scale
                 x_gr=False,
                 y_line=0,
                 y_gr=True
                 )
dc2_df = filter_df(gate_e_df,
                   x_gene='cite_ITGAX',
                   y_gene='cite_CD1C',
                   x_line=2,
                   y_line=1.7,
                   x_gr=True,
                   y_gr=True)
gate_f_df = filter_df(gate_e_df,
                      x_gene='cite_ITGAX',
                      y_gene='cite_CD1C',
                      x_line=5,
                      x_gr=False,
                      y_line=1.6,
                      y_gr=False
                      )
mnp_df = filter_df(gate_f_df,
                   x_gene='cite_ITGAX',
                   y_gene='cite_CD94',
                   x_line=1.7,
                   y_line=1.6,
                   x_gr=True,
                   y_gr=False)

# plot scatterplot
fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(14, 8))
plot_scatterplot(
    ax=axes[0, 0],
    df=adt_df,
    x_gene='cite_CD3E',
    y_gene='cite_CD19',
    patches=
    [
        patches.Rectangle((2, 0), 2.8, 1.7, linewidth=2, facecolor='none', edgecolor='black'),
        # t (2, 1.7)
        patches.Rectangle((0, 0), 1.4, 1.5, linewidth=2, facecolor='none', edgecolor='black'),
        # other (1.4, 1.5)
        patches.Rectangle((0, 1.7), 1.8, 3.3, linewidth=2, facecolor='none', edgecolor='black')
        # b (1.8, 1.7)
    ]
)
plot_scatterplot(
    ax=axes[0, 1],
    df=t_df,
    x_gene='cite_CD4',
    y_gene='cite_CD8A',
    patches=[
        patches.Rectangle((2, 0), 2.8, 1.7, linewidth=2, facecolor='none', edgecolor='black'),  # cd4 (2, 1.7)
        patches.Rectangle((0, 1.9), 2, 2.5, linewidth=2, facecolor='none', edgecolor='black'),  # cd8 (2, 1.9)
    ]
)
fig.delaxes(axes[0, 2])
fig.delaxes(axes[0, 3])
plot_scatterplot(
    ax=axes[1, 0],
    df=gate_b_df,
    x_gene='cite_CD14',
    y_gene='cite_NCAM1',
    patches=[
        patches.Rectangle((0, 1.6), 1.4, 1.8, linewidth=2, facecolor='none', edgecolor='black'),
        # nk (1.4, y = 1.6)
        patches.Rectangle((0, 0), 2.8, 1.5, linewidth=2, facecolor='none', edgecolor='black')
    ]
)
plot_scatterplot(
    ax=axes[1, 1],
    df=gate_d_df,
    x_gene='cite_HLA-DRB1',
    y_gene='cite_IL3RA',
    patches=[
        patches.Rectangle((1.5, 2.8), 4.2, 2, linewidth=2, facecolor='none', edgecolor='black'),
        # pdc (1.5, 2.8)
        patches.Rectangle((0, 0), 5, 2, linewidth=2, facecolor='none', edgecolor='black')
    ]
)
plot_scatterplot(
    ax=axes[1, 2],
    df=gate_e_df,
    x_gene='cite_ITGAX',
    y_gene='cite_CD1C',
    patches=[
        patches.Rectangle((2, 1.7), 3, 3.3, linewidth=2, facecolor='none', edgecolor='black'),  # cdc (2, 1.7)
        patches.Rectangle((0, 0), 5, 1.6, linewidth=2, facecolor='none', edgecolor='black')
    ]
)
plot_scatterplot(
    ax=axes[1, 3],
    df=gate_f_df,
    x_gene='cite_ITGAX',
    y_gene='cite_CD94',
    patches=[
        patches.Rectangle((1.7, 0), 4, 1.6, linewidth=2, facecolor='none', edgecolor='black'),  # mnp (1.7, 1.6)
    ]
)
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels)
fig.suptitle('CLR Expression', x=0.1)
plt.tight_layout()
plt.savefig('/projects/home/ikernin/projects/myocarditis/updated_datasets/figures/adt_scatterplots.png')
plt.close('all')


# plot hexbin
fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(14, 8))
plot_hexbin(
    ax=axes[0, 0],
    df=adt_df,
    x_gene='cite_CD3E',
    y_gene='cite_CD19',
    patches=
    [
        patches.Rectangle((2, 0), 2.8, 1.7, linewidth=2, facecolor='none', edgecolor='black'),
        # t (2, 1.7)
        patches.Rectangle((0, 0), 1.4, 1.5, linewidth=2, facecolor='none', edgecolor='black'),
        # other (1.4, 1.5)
        patches.Rectangle((0, 1.7), 1.8, 3.3, linewidth=2, facecolor='none', edgecolor='black')
        # b (1.8, 1.7)
    ]
)
plot_hexbin(
    ax=axes[0, 1],
    df=t_df,
    x_gene='cite_CD4',
    y_gene='cite_CD8A',
    patches=[
        patches.Rectangle((2, 0), 2.8, 1.7, linewidth=2, facecolor='none', edgecolor='black'),  # cd4 (2, 1.7)
        patches.Rectangle((0, 1.9), 2, 2.5, linewidth=2, facecolor='none', edgecolor='black'),  # cd8 (2, 1.9)
    ]
)
fig.delaxes(axes[0, 2])
fig.delaxes(axes[0, 3])
plot_hexbin(
    ax=axes[1, 0],
    df=gate_b_df,
    x_gene='cite_CD14',
    y_gene='cite_NCAM1',
    patches=[
        patches.Rectangle((0, 1.6), 1.4, 1.8, linewidth=2, facecolor='none', edgecolor='black'),
        # nk (1.4, y = 1.6)
        patches.Rectangle((0, 0), 2.8, 1.5, linewidth=2, facecolor='none', edgecolor='black')
    ]
)
plot_hexbin(
    ax=axes[1, 1],
    df=gate_d_df,
    x_gene='cite_HLA-DRB1',
    y_gene='cite_IL3RA',
    patches=[
        patches.Rectangle((1.5, 2.8), 4.2, 2, linewidth=2, facecolor='none', edgecolor='black'),
        # pdc (1.5, 2.8)
        patches.Rectangle((0, 0), 5, 2, linewidth=2, facecolor='none', edgecolor='black')
    ]
)
plot_hexbin(
    ax=axes[1, 2],
    df=gate_e_df,
    x_gene='cite_ITGAX',
    y_gene='cite_CD1C',
    patches=[
        patches.Rectangle((2, 1.7), 3, 3.3, linewidth=2, facecolor='none', edgecolor='black'),  # cdc (2, 1.7)
        patches.Rectangle((0, 0), 5, 1.6, linewidth=2, facecolor='none', edgecolor='black')
    ]
)
plot_hexbin(
    ax=axes[1, 3],
    df=gate_f_df,
    x_gene='cite_ITGAX',
    y_gene='cite_CD94',
    patches=[
        patches.Rectangle((1.7, 0), 4, 1.6, linewidth=2, facecolor='none', edgecolor='black'),  # mnp (1.7, 1.6)
    ]
)
fig.suptitle('CLR Expression', x=0.1)
plt.tight_layout()
plt.savefig('/projects/home/ikernin/projects/myocarditis/updated_datasets/figures/adt_hexbins.pdf')
plt.close('all')
