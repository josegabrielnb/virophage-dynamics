#!/usr/bin/python3

#Parse simulation outputs and write a summary dataframe

import os
import re
import pandas as pd

my_dataframe = pd.DataFrame(columns=['cell_number', 'treatment', 'run', 'surviving_cells', 'max_virophage_amplitude', 'max_virus_amplitude', 't_extinction'])

os.chdir("/Users/user/Desktop/Julia2.0/simulation_results/simulations/")

my_folders = os.listdir()

for folder in my_folders:

    if re.match("^cells", folder):

        cell_number = re.match("^cells(.*)_", folder).group(1)

        #print(cell_number)

        if re.match(".*neutral.*", folder):

            treatment = 1

        elif re.match(".*inhibition.*", folder):

            treatment = 2

        elif re.match(".*pcd.*", folder):

            treatment = 3

        #print(cell_number,treatment)

        os.chdir(folder)

        runs = os.listdir()

        for run in runs:

            if re.match("^run", run):

                run_number = re.match("^run_(.*)", run).group(1)

                os.chdir(run)

                my_files = os.listdir()

                for file in my_files:

                    if file.endswith(".csv"):

                        #print(file)

                        df = pd.read_csv(file)

                        #print(df)

                        surviving_cells = df['cells'].iloc[-1]
                        max_virophage_amplitude = df['virophages'].max()
                        max_virus_amplitude = df['viruses'].max()
                        t_extinction = sum(df['viruses']>1)+1

                        my_dataframe_length = len(my_dataframe)

                        my_dataframe.loc[ my_dataframe_length] = [cell_number,treatment,run_number,surviving_cells,max_virophage_amplitude,max_virus_amplitude,t_extinction]

                        print(folder,run)

                os.chdir("..")

        os.chdir("..")

my_dataframe.to_csv("dynamics_summary_100replicates.csv",index=False)