import pandas as pd
def make_MM(num_rxns,rxn_vol):
    '''
    Calculates the amount of each reagent to add to reach the desired master mix
    '''
    cutsmart = 1 * num_rxns
    atp = 1 * num_rxns
    vector = 0.25 * num_rxns
    ligase = 0.5 * num_rxns
    enzyme = 0.25 * num_rxns
    water = (rxn_vol - ((cutsmart + atp + vector + ligase + enzyme)/num_rxns)) * num_rxns
    master_mix = pd.DataFrame(
        {'Component':['H2O','Cutsmart','ATP','Vector','T4 Ligase','Restriction Enzyme','Total'],
        'Amount':[water,cutsmart,atp,vector,ligase,enzyme,rxn_vol*num_rxns]},
        columns=["Component","Amount"]
    )
    return master_mix
