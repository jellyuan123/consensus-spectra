import numpy as np
from sklearn.cluster import KMeans
import pandas as pd
import matplotlib.pyplot as plt

def initializeDF(df) :
        # print(df.columns)
        allSMILESGroups = df.groupby('smiles')
        # cluster_sizes = allSMILESGroups.size()

        # Filter groups by size (e.g., keep groups with size > 1)
        # filtered_clusters = cluster_sizes[cluster_sizes > 200]

        # Check the number of filtered groups
        # num_clusters = len(filtered_clusters)

        # Display the result
        # print(f"Number of groups with size > 200: {num_clusters}")
        # print(f"Avg Size: {cluster_sizes.mean()}")

        # largest_group_key = cluster_sizes.idxmax()  # Find the key of the largest group
        # getGroupFromSmiles(allSMILESGroups, largest_group_key)
        return allSMILESGroups

def getGroupFromSmiles(allSMILESGroups, smilesString):
        group = allSMILESGroups.get_group(smilesString)  # Retrieve the largest group

        # Print the largest group
        # print(f"SMILES of group: {smilesString}")
        # print(f"NUmber of spectra of group: {len(group)}")
        return group
              
        
def genListOfAllPeaks(group):

        all_mzs = []
        all_intensities = []
        seen_mzs = set()

        def addToAll(mz, intensity):
                seen_mzs.add(mz)
                all_mzs.append(mz)
                all_intensities.append(intensity)

        for _, row in group.iterrows():
                mzs = row['mzs'] 
                intensities = row['intensities']

                for x, y in enumerate(mzs):
                        mz = mzs[x]
                        intensity = intensities[x]
                        if mz in seen_mzs:
                                index = all_mzs.index(mz)
                                diff = abs(all_intensities[index] - intensity)
                                #if that exact peak has shown up and the intensities are within
                                # .01 of each other, average them - else add to list
                                if diff > .01:
                                        addToAll(mz, intensity)
                                else:
                                        all_intensities[index] = np.mean([all_intensities[index], intensity])  
                
                        else:
                                seen_mzs.add(mz)
                                all_mzs.append(mz)
                                all_intensities.append(intensity)
        
        # print(len(all_mzs), len(all_intensities))
        return all_mzs, all_intensities


def getAllVariation(allSMILESGroups):

        all_pcs = []
        all_adducts = []
        all_instruments = []
        all_collision_energy = []
        all_fold = []

        for smiles_value, group in allSMILESGroups:
                num_adducts = group['adduct'].nunique()
                all_adducts.append(num_adducts)

                num_instruments = group['instrument_type'].nunique()
                all_instruments.append(num_instruments)

                num_energy = group['collision_energy'].nunique()
                all_collision_energy.append(num_energy)

                num_fold = group['fold'].nunique()
                all_fold.append(num_fold)

                num_pcs = group['precursor_mz'].nunique()
                all_pcs.append(num_pcs)

        all_variation = {key: [] for key in ['precursor_mz', 'adduct', 'instrument_type', 'collision_energy', 'fold']}

        all_variation['precursor_mz'].append(np.mean(all_pcs))
        all_variation['precursor_mz'].append(np.std(all_pcs))

        all_variation['adduct'].append(np.mean(all_adducts))
        all_variation['adduct'].append(np.std(all_adducts))

        all_variation['instrument_type'].append(np.mean(all_instruments))
        all_variation['instrument_type'].append(np.std(all_instruments))

        all_variation['collision_energy'].append(np.mean(all_collision_energy))
        all_variation['collision_energy'].append(np.std(all_collision_energy))

        all_variation['fold'].append(np.mean(all_fold))
        all_variation['fold'].append(np.std(all_fold))

        return all_variation



def checkVariation(all_variation, group):

        keys = ['precursor_mz', 'adduct','instrument_type', 'collision_energy','fold']
        variation_flags = dict.fromkeys(keys, 0)
        variation_flags = {key: 0 for key in keys}


        for key in keys:
                group_num = group[key].nunique()
                all_variation_min = all_variation[key][0] - all_variation[key][1]
                all_variation_max = all_variation[key][0] + all_variation[key][1]

                if group_num > all_variation_max:
                        variation_flags[key] = 1
                
                if group_num < all_variation_min:
                        variation_flags[key] = -1
                
                else:
                        variation_flags[key] = 0
        
        return variation_flags

def createClusters(group, all_mzs, all_intensities):
        mean_length = group['mzs'].apply(len).mean()
        std_dev = group['mzs'].apply(len).std()

        # print(mean_length)
        # print(std_dev)

        data = {'mzs': all_mzs,
                'intensities': all_intensities}

        new_df = pd.DataFrame(data)
        num_clusters = int(mean_length)

        new_df = pd.DataFrame({'mzs': all_mzs, 'intensities': all_intensities})

        # KMeans Clustering
        kmeans = KMeans(n_clusters=num_clusters)
        new_df['Cluster'] = kmeans.fit_predict(new_df[['mzs', 'intensities']])

        return new_df



def getFinalPeaks(new_df):

        final_filtered_data = pd.DataFrame()
        group_clusters = new_df['Cluster'].unique()
        final_mzs = []
        final_intensities = []

        for cluster in group_clusters:
                cluster_data = new_df[new_df['Cluster'] == cluster].copy()
                # print(len(cluster_data))

                min_mz = cluster_data['mzs'].min()
                max_mz = cluster_data['mzs'].max()

                if (len(cluster_data) != 1):
                        bins = np.arange(min_mz, max_mz + 0.02, 0.02)

                        cluster_data['Bin'] = pd.cut(cluster_data['mzs'], bins, right=False)

                        most_populous_bin = cluster_data['Bin'].value_counts().idxmax()

                        filtered_cluster_data = cluster_data[cluster_data['Bin'] == most_populous_bin]

                        final_filtered_data = pd.concat([final_filtered_data, filtered_cluster_data], ignore_index=True)

                        avg_mzs = filtered_cluster_data['mzs'].mean()
                        final_mzs.append(avg_mzs)
                        avg_intensities = filtered_cluster_data['intensities'].mean()
                        final_intensities.append(avg_intensities)
                        final_filtered_data.drop(columns=['Bin'], inplace=True)

                else:
                        avg_mzs = cluster_data['mzs'].iloc[0]  # Get mz of the single entry
                        final_mzs.append(avg_mzs)

                        avg_intensities = cluster_data['intensities'].iloc[0]  # Get intensity of the single entry
                        final_intensities.append(avg_intensities)

                        final_filtered_data = pd.concat([final_filtered_data, cluster_data], ignore_index=True)

        return final_mzs, final_intensities




def createFinalSpectra(final_mzs, final_intensities, group):

        precursor_counts = group['precursor_mz'].value_counts()
        pc = precursor_counts.idxmax()

        adduct_counts = group['adduct'].value_counts()
        adduct = adduct_counts.idxmax()

        collision_energy_counts = group['collision_energy'].value_counts()
        if len(collision_energy_counts) == 0: 
                ce = 'NaN'
        else:
                ce = collision_energy_counts.idxmax()
            
        fold_counts = group['fold'].value_counts()
        f = fold_counts.idxmax()
        smile = group['smiles'].iloc[0]

        my_consensus_spectra = {
                "mzs": final_mzs,
                "intensities": final_intensities, 
                "smiles": smile,
                "precursor_mz": pc, 
                "adduct": adduct,
                "collision_energy": ce,
                "fold": f
        }

        return smile, my_consensus_spectra

def genAllSpectra(df):
        allSMILESGroups = initializeDF(df)
        smiles_list = [str(smile) for smile in allSMILESGroups['smiles'].unique()]
        print(f"{len(smiles_list)} spectra to generate")
        all_consensus = {}

        i = 0
        j = 0

        for sstring in smiles_list:
                cleaned_smiles = sstring[2:-2]
                try:
                        group = getGroupFromSmiles(allSMILESGroups, cleaned_smiles)
                        
                        all_group_mzs, all_group_intensities = genListOfAllPeaks(group)

                        all_variation = getAllVariation(allSMILESGroups)
                        group_variation = checkVariation(all_variation, group)

                        new_df = createClusters(group, all_group_mzs, all_group_intensities)
                        final_mzs, final_intensities = getFinalPeaks(new_df)
                        smiles, consensus_spectra = createFinalSpectra(final_mzs, final_intensities, group)
                        all_consensus[smiles] = [consensus_spectra, group_variation]
                        if (i % 20 == 0):
                                print(f"{j * 20 + 1} consensus spectra created ")
                                j += 1
                        i += 1
                except:
                        print(f"Error with {cleaned_smiles}")

        return all_consensus





        