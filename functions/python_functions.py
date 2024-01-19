import os
import pegasus as pg
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import matplotlib.patheffects as path_effects
from tqdm import tqdm


# colormap for hexplots
blues_cmap = clr.LinearSegmentedColormap.from_list('gene_cmap', ["#e0e0e1", '#4576b8', '#02024a'], N=200)

# umap cluster colors
tissue_global_pal = ['#ff0029', '#377eb8', '#66a61e', '#984ea3', '#00d2d5',
                     '#af8d00', '#7f80cd', '#b3e900', '#c42e60', '#ff7f00']
tissue_mnp_pal = ['#1b9e77', '#e7298a', '#a6761d', '#252525', '#f43600', '#356F83', '#eff26e', '#80b1d3']
tissue_t_pal = ['#cf8c00', '#ff4040', '#0097ff', '#00d067', '#bdbdbd', '#8a2be2', '#80b1d3']
tissue_nonimmune_pal = ['#a65628', '#f781bf', '#8dd3c7', '#bebada', '#fb8072', '#fdb462',
                        '#fccde5', '#bc80bd', '#ffed6f', '#c4eaff', '#d95f02', '#737373',
                        '#4ba93b', '#5779bb', '#927acc', '#bf3947', '#f48758', '#80b1d3']
tissue_fibroblast_pal = ['#a65628', '#f781bf', '#8dd3c7', '#bebada', '#fb8072', '#80b1d3', '#fdb462']
tissue_b_pal = ['#ff7373', '#f2b94f', '#80b1d3']
blood_global_pal = ['#1B9E77', '#D95F02', '#7570B3', '#66A61E', '#A6761D', '#E6AB02', '#E7298A']
blood_mnp_pal = ['#03396c', '#b3e900', '#cf8c00', '#f781bf', '#8dd3c7', '#fb8072', '#80b1d3',
                 '#7f80cd', '#000000', '#969696', '#cca3bf', '#4b7487', '#e7298a']
blood_cd8_pal = ['#edba00', '#8d86bb', '#f49ec4', '#90278f', '#8cd7f8', '#a12934', '#bf242a', '#cc702d',
                 '#000000', '#7b7b7b', '#00914c', '#2e3192', '#682326']
blood_cd4_pal = ['#266489', '#68B9C0', '#90D585', '#F3C151', '#F37F64', '#8F97A4', '#76846E', '#A65B69']
blood_b_pal = ['#00d067', '#927acc', '#9f5b00', '#f2b94f', '#e43872', '#4d5c82', '#bf3947', '#8caed6']


def plot_umap(df, title, color_palette,
              fig_size=(5, 5), wspace=0.8, marker_multiplier=20,
              cluster_font_size=12, axis_font_size=16, ncol_legend=1,
              save_fig=False, save_name=None):

    # get number of cells per cluster
    n_cells = pd.DataFrame(df.obs[['umap_number', 'umap_name']].value_counts()).reset_index()
    n_cells.columns = ['umap_number', 'umap_name', 'n_cells']
    n_cells = n_cells.sort_values('umap_number')
    dict_n_cells = dict(zip(n_cells['umap_name'], n_cells['n_cells']))

    # create labels for legend
    labels = []
    for clust, n in dict_n_cells.items():
        labels.append(f"{clust} (n={n:,})")

    # create base_figure
    fig = pg.scatter(df,
                     basis="umap",
                     attrs="umap_number",
                     legend_loc='on data',
                     palettes=','.join(color_palette),
                     panel_size=fig_size,
                     wspace=wspace,
                     return_fig=True)
    fig.patch.set_visible(True)
    ax = fig.axes[0]

    # add total number of cells
    ax.annotate(f'{df.shape[0]:,} cells',
                xy=(0.01, 0), xycoords='axes fraction',
                fontsize=cluster_font_size,
                horizontalalignment='left',
                verticalalignment='bottom')

    # edit cluster number font
    ax_is_text = [isinstance(i, matplotlib.text.Text) for i in ax.get_children()]
    ax_text_idx = np.where(ax_is_text)[0]
    for idx in ax_text_idx:
        if ax.get_children()[idx]._text in df.obs['umap_number'].cat.categories.astype('string'):
            text = ax.get_children()[idx]
            text.set_path_effects([path_effects.Stroke(linewidth=4, foreground='white'),
                                   path_effects.Normal()])
            text.set_fontsize(cluster_font_size)

    # add legend
    ax.legend(labels,
              bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0,
              markerscale=marker_multiplier, frameon=False,
              prop={'size': cluster_font_size},
              ncol=ncol_legend)

    # edit font sizes
    ax.xaxis.label.set_fontsize(axis_font_size)
    ax.yaxis.label.set_fontsize(axis_font_size)
    ax.set_title(title, fontsize=axis_font_size)

    # save and/or plot figure
    if save_fig:
        plt.savefig(save_name)
    plt.show()
    plt.close()


def make_gene_dotplot(df, cluster_order, gene_order, title, figsize=(8, 4), names='umap_name', save=False, save_name=None):
    # remove doublets if any
    df = df[['doublets' not in x.lower() for x in df.obs[names]], ].copy()
    df.obs[names] = df.obs[names].cat.remove_unused_categories()

    # set figure font
    plt.rcParams.update({'font.sans-serif': 'Arial'})
    plt.rcParams.update({'font.size': 19})

    # create figure
    fig, ax = plt.subplots(figsize=(12, 6))
    dotplot = sc.pl.dotplot(df, gene_order, names,
                            var_group_rotation=30, standard_scale='var', ax=ax,
                            categories_order=cluster_order,
                            figsize=figsize,
                            return_fig=True)
    dotplot.legend(width=2, colorbar_title='Scaled\nExpression')
    dotplot.make_figure()
    axes_dict = dotplot.get_axes()
    axes_dict["mainplot_ax"].set_axisbelow(True)
    axes_dict["mainplot_ax"].xaxis.tick_top()
    axes_dict["mainplot_ax"].set_title(f'{title}: Marker Genes')
    fig.tight_layout()
    if save:
        plt.savefig(save_name, dpi=300)
    plt.show()
    plt.close()




def hex_plot(df, title, n_genes=True, gridsize=200, cmap='YlOrRd',
             save=False, save_name=None):
    if n_genes:
        col = 'n_genes'
    else:
        col = 'percent_mito'

    # get umap coordinates
    umap_coords = pd.DataFrame(df.obsm['X_umap'], columns=['x', 'y'])
    x = umap_coords['x']
    y = umap_coords['y']

    # get data to plot
    data = df.obs[col]

    # create plot
    fig, ax = plt.subplots(figsize=(6, 4))
    hb = ax.hexbin(x=x,
                   y=y,
                   C=data,
                   cmap=cmap,
                   gridsize=gridsize,
                   edgecolors = "none")
    cb = fig.colorbar(hb, ax=ax, shrink=.75, aspect=10)
    cb.ax.set_title(col, loc='left', fontsize=14)
    ax.set_ylabel('UMAP2', fontsize=18)
    ax.set_xlabel('UMAP1', fontsize=18)
    ax.set_title(title, fontsize=18)
    ax.tick_params(left=False, labelleft=False,
                   bottom=False, labelbottom=False)

    if save:
        plt.savefig(save_name)
    plt.show()
    plt.close()


def hex_featureplot(df, gene, gridsize=200, cmap='YlOrRd',
                    save=False, save_name=None):
    # get umap coordinates for base
    umap_coords = pd.DataFrame(df.obsm['X_umap'], columns=['x', 'y'])
    x = umap_coords['x']
    y = umap_coords['y']

    # get gene counts
    norm_counts = df[:, gene].copy().get_matrix('X').todense().transpose()
    norm_counts = np.squeeze(np.asarray(norm_counts))

    if norm_counts.size == 0:  # ie all empty/zero in sparce matrix
        norm_counts = np.zeros(norm_counts.shape[1])

    # get ncells and perc cells
    ncells = (norm_counts != 0).sum()
    pcells = ncells / norm_counts.shape[0]

    # create plot
    fig, ax = plt.subplots(figsize=(6, 4))
    hb = ax.hexbin(x=x,
                   y=y,
                   C=norm_counts,
                   cmap=cmap,
                   gridsize=gridsize,
                   edgecolors = "none")
    ax.annotate(f'{ncells:,} ({pcells:.1%}) cells',
                xy=(0.01, 0), xycoords='axes fraction',
                fontsize=10,
                horizontalalignment='left',
                verticalalignment='bottom')
    cb = fig.colorbar(hb, ax=ax, shrink=.75, aspect=10)
    cb.ax.set_title('logCPM', loc='left', fontsize=14)
    ax.set_ylabel('UMAP2', fontsize=18)
    ax.set_xlabel('UMAP1', fontsize=18)
    ax.set_title(gene, fontsize=18)
    ax.tick_params(left=False, labelleft=False,
                   bottom=False, labelbottom=False)
    if save:
        plt.savefig(save_name)
    plt.show()
    plt.close()



def multi_hex_featureplot(df, genes, ncol, cmap='YlOrRd', gridsize=120, panel_size=(6, 4),
                          labels=True, save=False, save_name=None, dpi=300):
    # get umap coordinates for base
    umap_coords = pd.DataFrame(df.obsm['X_umap'], columns=['x', 'y'])
    x = umap_coords['x']
    y = umap_coords['y']

    # set up figure structure
    nrow = int(np.ceil(len(genes) / ncol))
    fig_size = (panel_size[0] * ncol, panel_size[1] * nrow)
    fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=fig_size,
                             sharex=True, sharey=True)
    ax = axes.ravel()

    # plot for each gene
    for num, gene in enumerate(tqdm(genes)):
        # get counts for each gene
        norm_counts = df[:, gene].copy().get_matrix('X').todense().transpose()
        norm_counts = np.squeeze(np.asarray(norm_counts))

        if norm_counts.size == 0:  # ie all empty/zero in sparce matrix
            norm_counts = np.zeros(norm_counts.shape[1])

        # get ncells and perc cells
        ncells = (norm_counts != 0).sum()
        pcells = ncells / norm_counts.shape[0]

        # create heatmap
        hb = ax[num].hexbin(x=x,
                            y=y,
                            C=norm_counts,
                            cmap=cmap,
                            gridsize=gridsize,
                            edgecolors = "none")

        # add percent expression
        ax[num].annotate(f'{ncells:,} ({pcells:.1%}) cells',
                xy=(0.01, 0), xycoords='axes fraction',
                fontsize=10,
                horizontalalignment='left',
                verticalalignment='bottom')

        # add axes and colorbar information
        if labels:
            cb = fig.colorbar(hb, ax=ax[num], shrink=.75, aspect=10)
            cb.ax.set_title('logCPM', loc='left', fontsize=14)
            ax[num].set_title(gene, fontsize=18)
            ax[num].tick_params(left=False, labelleft=False,
                                bottom=False, labelbottom=False)
            if (num + ncol) % ncol == 0:
                # start of row
                ax[num].set_ylabel('UMAP2', fontsize=18)
            if nrow == 1:
                # only one row
                ax[num].set_xlabel('UMAP1', fontsize=18)
            elif num > (len(genes) - ncol - 1):
                # last row if more than one row
                ax[num].set_xlabel('UMAP1', fontsize=18)
        else:
            ax[num].set_title(gene, fontsize=22, loc='left')
            ax[num].tick_params(left=False, labelleft=False,
                                bottom=False, labelbottom=False)
    for i in range(len(ax) - (len(ax) - len(genes)), len(ax)):
        ax[i].set_axis_off()
    fig.tight_layout()
    if save:
        plt.savefig(save_name, dpi=300)
    plt.show()
    plt.close()


def multi_hexfp_by_condition(df, genes, cmap='YlOrRd', gridsize=120, panel_size=(6, 4), 
                            save=False, save_name=None, dpi=300):
    # get umap coordinates for base
    ctrl_df = df[df.obs['condition'] == 'control'].copy()
    ctrl_coords = pd.DataFrame(ctrl_df.obsm['X_umap'], columns=['x', 'y'])
    ctrl_x = ctrl_coords['x']
    ctrl_y = ctrl_coords['y']

    myo_df = df[df.obs['condition'] == 'myocarditis'].copy()
    myo_coords = pd.DataFrame(myo_df.obsm['X_umap'], columns=['x', 'y'])
    myo_x = myo_coords['x']
    myo_y = myo_coords['y']

    x = [ctrl_x, myo_x]
    y = [ctrl_y, myo_y]

    # set up figure structure
    ncol = 2
    nrow = len(genes)
    fig_size = (panel_size[0] * ncol, panel_size[1] * nrow)
    fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=fig_size,
                             sharex=True, sharey=True)
    ax = axes.ravel()
    title = ['Control', 'Myocarditis']

    # plot for each gene
    for num, gene in enumerate(tqdm(genes)):
        # get gene counts
        ctrl_cts = ctrl_df[:, gene].copy().get_matrix('X').todense().transpose()
        ctrl_cts = np.squeeze(np.asarray(ctrl_cts))
        if ctrl_cts.size == 0:  # ie all empty/zero in sparce matrix
            ctrl_cts = np.zeros(ctrl_cts.shape[1])
        ctrl_n = (ctrl_cts != 0).sum()
        ctrl_perc = ctrl_n / ctrl_cts.shape[0]

        myo_cts = myo_df[:, gene].copy().get_matrix('X').todense().transpose()
        myo_cts = np.squeeze(np.asarray(myo_cts))
        if myo_cts.size == 0:  # ie all empty/zero in sparce matrix
            myo_cts = np.zeros(myo_cts.shape[1])
        myo_n = (myo_cts != 0).sum()
        myo_perc = (myo_cts != 0).sum() / myo_cts.shape[0]

        norm_counts = [ctrl_cts, myo_cts]
        cb_max = max(max(ctrl_cts), max(myo_cts))  # makes sure colobars are the same

        # create heatmap
        for i in range(2):
            hb = ax[num * 2 + i].hexbin(x=x[i],
                                        y=y[i],
                                        C=norm_counts[i],
                                        cmap=cmap,
                                        gridsize=gridsize,
                                        vmin=0,
                                        vmax=cb_max,
                                        edgecolors = "none")
            cb = fig.colorbar(hb, ax=ax[num * 2 + i], shrink=.75, aspect=10)
            cb.ax.set_title('logCPM', loc='left', fontsize=14)
            if i == 0:
                # ie control
                ax[num * 2 + i].set_ylabel('UMAP2', fontsize=18)
                ax[num * 2 + i].annotate(f'{ctrl_n:,} ({ctrl_perc:.1%}) cells',
                               xy=(0.01, 0), xycoords='axes fraction',
                               fontsize=10,
                               horizontalalignment='left',
                               verticalalignment='bottom')
            if num + 1 == len(genes):
                # ie control and only one row
                ax[num * 2 + i].set_xlabel('UMAP1', fontsize=18)
            if i == 1:
                # ie case
                ax[num * 2 + i].annotate(f'{myo_n:,} ({myo_perc:.1%}) cells',
                               xy=(0.01, 0), xycoords='axes fraction',
                               fontsize=10,
                               horizontalalignment='left',
                               verticalalignment='bottom')
            ax[num * 2 + i].set_title(f'{title[i]}: {gene}', fontsize=18)
            ax[num * 2 + i].tick_params(left=False, labelleft=False,
                                        bottom=False, labelbottom=False)
    fig.tight_layout()
    if save:
        plt.savefig(save_name, dpi=300)
    plt.show()
    plt.close()


def get_stacked_bar_df(df, lineage, cluster_col='umap_name'):
    print(f'Getting stacked bar info for: {lineage}')

    # remove post-steroid samples and CD45+
    post_samps = ['SIC_264_1047_1_heart', 'SIC_177_539_1_heart', 'SIC_177_539_2_heart', 'SIC_237_662_1_heart',
                  'SIC_175_551_1_heart']
    filtered_df = df[(~df.obs['Channel'].isin(post_samps)) &
                     (df.obs['source'] != 'CD45+') &
                     (df.obs[cluster_col] != 'Doublets/RBCs')]
    filtered_df.obs[cluster_col] = filtered_df.obs[cluster_col].cat.remove_unused_categories()
    filtered_df.obs['donor'] = filtered_df.obs['donor'].cat.remove_unused_categories()
    donor_cluster_df = filtered_df.obs[['donor', cluster_col]]

    # obtain information per donor
    perc_dict = {}
    for donor in donor_cluster_df['donor'].cat.categories:
        donor_df = donor_cluster_df[donor_cluster_df['donor'] == donor].copy()
        clust_perc_df = pd.DataFrame(donor_df[cluster_col].value_counts()).reset_index()
        clust_perc_df.columns = ['cluster_names', 'n_cells']
        clust_perc_df['perc'] = clust_perc_df['n_cells'] / sum(clust_perc_df['n_cells']) * 100
        clust_perc_df['donor'] = donor
        perc_dict[donor] = clust_perc_df

    # get data for plotting in R
    plot_df = pd.concat(perc_dict.values())
    return plot_df


def get_pseudobulk_info(df, save_name, sample_col='donor', cluster_col='umap_number', cols=['condition']):
    # remove post-steroid samples
    df = df[df.obs['on_steroids'] != 'True']

    # for each sample, get the cells in each cluster and sum them to form a pseudo-sample
    gene_sum_dict = {}
    meta_dict = {}
    for samp in set(df.obs[sample_col]):
        for clust in set(df.obs[cluster_col]):
            # get subset of data
            dat = df[(df.obs[sample_col] == samp) & (df.obs[cluster_col] == clust)].copy()
            if dat.shape[0] == 0:
                continue  # ie no cells from current sample in current cluster
            key = f'{samp}_c{clust}'
            # add to gene_sum_dict by summing counts
            count_sum = np.array(dat.get_matrix('raw.X').sum(axis=0)).flatten()
            gene_sum_dict[key] = count_sum
            # add to meta_dict by getting meta info
            samp_clust_dict = {'n_cells': dat.shape[0], 'cluster': clust, 'sample': samp}
            for col in cols:
                assert len(dat.obs[col].unique() == 1)
                samp_clust_dict[col] = dat.obs[col].iloc[0]
            meta_dict[key] = samp_clust_dict
    # put together pseudobulk counts matrix
    count_mtx = pd.DataFrame(gene_sum_dict, index=df.var_names)
    count_mtx.to_csv(f'{save_name}_pseudocounts.csv')
    # put together meta dataframe
    meta_df = pd.DataFrame.from_dict(meta_dict, orient='index', columns=['n_cells', 'cluster', 'sample'] + cols)
    meta_df.to_csv(f'{save_name}_metainfo.csv')


def blood_pb_info(adata_path, sample_label, cluster_label, cols):
    adata = pg.read_input(adata_path)

    adata.obs[cluster_label] = adata.obs[cluster_label].str.replace(' ', '_')

    # Iterate across samples and clusters
    gene_sum_dict = {}
    cell_num_dict = {}
    for samp in set(adata.obs[sample_label]):
        for clust in set(adata.obs[cluster_label]):

            dat = adata[(adata.obs[sample_label] == samp) & (adata.obs[cluster_label] == clust)].copy()
            if dat.shape[0] == 0:
                continue
            # Add info to my dictionaries
            key = f"{samp}_c{clust}"
            samp_clust_dict = {'n_cells': dat.shape[0], 'cluster': clust, 'sample': samp}
            for col in cols:
                assert len(dat.obs[col].unique() == 1), f"{sample_label} does not have unique {col}"
                samp_clust_dict[col] = dat.obs[col].iloc[0]
            cell_num_dict[key] = samp_clust_dict

            # Sum the counts
            count_sum = np.array(dat.get_matrix('raw.X').sum(axis=0)).flatten()
            gene_sum_dict[key] = count_sum

    count_mtx = pd.DataFrame(gene_sum_dict, index=adata.var_names)

    des_mtx = pd.DataFrame.from_dict(cell_num_dict, orient='index', columns=['n_cells', 'cluster', 'sample'] + cols)

    return count_mtx, des_mtx
