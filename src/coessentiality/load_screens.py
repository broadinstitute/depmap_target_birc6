import pandas as pd
from sklearn.decomposition import PCA
from taigapy import TaigaClient
tc = TaigaClient()

def load_screens():
    # Load screens
    screens = pd.read_csv('https://ndownloader.figshare.com/files/22543691').T
    screens.index = screens.index.str.split(' ').str.get(0)

    # Map Broad ID to CCLE name
    cell_lines = pd.read_csv('https://ndownloader.figshare.com/files/21522000')
    cell_lines=cell_lines.loc[:,['DepMap_ID','CCLE Name']]
    cell_lines=cell_lines.set_index('DepMap_ID')
    cell_lines=cell_lines.squeeze()


    ## Change names to CCLE IDs
    screens.columns = cell_lines[screens.columns].values

    # Lines with contamination
    lines_to_drop = screens.columns[screens.isna().sum() > 1000].values

    screens.drop(lines_to_drop, axis=1, inplace=True)
    #Drop genes from X and Y chromosome
    screens = screens.dropna()


    # Bias-correct using "molecular function:olfactory receptor activity" genes
    olfactory_genes = pd.read_csv('olfactory_genes.txt', header=None, squeeze=True)
    olfactory_data = screens.reindex(olfactory_genes).dropna()
    transformation = PCA(n_components=4)
    transformation.fit(olfactory_data)
    top_PC_effects = transformation.inverse_transform(transformation.transform(screens))
    screens -= top_PC_effects
    screens = screens.iloc[:, :-4]
    print(len(screens.index))
    return screens
