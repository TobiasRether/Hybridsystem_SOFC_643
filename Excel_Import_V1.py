import pandas as pd
import os

def get_Excel():
    file_path = input('Enter Filepath:') # e.g. C:/Users/Guest/Desktop
    file_name = input('Enter File-Name: ') # File name including file type e.g. Test.xlsx
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    if '/' in file_path:
        file_path = file_path + '/' + file_name
    else:
        file_path = file_path + '\\' + file_name
        file_path = file_path.replace('\\','/')
    if '"' in file_path:
        file_path = file_path.replace('"','')
    if not os.path.isfile(file_path):
        print('File does not exist.')
        exit()
    else:
        data = pd.read_excel(file_path)
        df = pd.DataFrame(data)
        print(df)
    val1 = input('Enter LSV: ')
    val2 = input('Enter TEB: ')
    try:
        a=df[(df['LSV'] == int(val1)) & (df['TEB'] == int(val2))].index[0]
    except IndexError:
        print('LSV and TEB incorrect.')
        exit()
    data1 = [{'MV1': df.at[int(a),'MV1'],'MENEXT%': df.at[int(a),'MENEXT%'],'PIV': df.at[int(a),'PIV'],'ETAVI': df.at[int(a),'ETAVI']}]
    index = 'LSV'+': '+str(val1)+'; '+'TEB'+': '+str(val2)
    df1 = pd.DataFrame(data1, index=[index])
    print(df1)

get_Excel()