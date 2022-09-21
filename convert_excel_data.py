# -*- coding: utf-8 -*-
"""
Said El-Barbari
EA-303
22.01.2021

"""
import os
import openpyxl
import numpy as np



path = os.path.dirname(__file__)

data_sheet_path = path +'/data/'
file_name = '2500A-3000A-Messung.xlsx'
excel_sheet_name=['030degC','125degC','177degC']


# Id Vds t
wb_col=[[0,1,4,5,8,9,12],[0,1,12],[0,1,8,9,12]]
wb_row=[5,50006]
#wb_row=[5,500]

#%%


wb_Rds_on = openpyxl.load_workbook(data_sheet_path + file_name)

data=[]
for sheets,col in zip(excel_sheet_name,wb_col):
    IU_data=np.zeros((wb_row[1]-wb_row[0],len(col)))
    ws_Rds_on = wb_Rds_on[sheets]
    for cn,curr_col in enumerate(col):
        for idx in range(wb_row[0],wb_row[1]): # get time current and voltage
            # MAP_data[idx,curr_col] = wb_Rds_on.cell_value(idx, curr_col) # Read the data in the current cell
            IU_data[idx-wb_row[0],cn] = ws_Rds_on.cell(idx+1,curr_col+1).value # Read the data in the current cell
    data.append(IU_data)
#%%
data_dic={
    
    'Id_30degC_2500A'  :data[0][:,0],
    'Vds_30degC_2500A' :data[0][:,1],
    'Id_30degC_3000A'  :data[0][:,2],
    'Vds_30degC_3000A' :data[0][:,3],
    'If_30degC_3000A'  :data[0][:,4],
    'Vsd_30degC_3000A' :data[0][:,5],
    't_30degC'         :data[0][:,6],
    
    'Id_125degC_2500A' :data[1][:,0],
    'Vds_125degC_2500A':data[1][:,1],
    't_125degC'        :data[1][:,2],
    
    'Id_177degC_2500A' :data[2][:,0],
    'Vds_177degC_2500A':data[2][:,1],
    'If_177degC_2500A' :data[2][:,2],
    'Vsd_177degC_2500A':data[2][:,3],
    't_177degC'        :data[2][:,4]}


np.savez('VI_curve_data.npz', **data_dic)


