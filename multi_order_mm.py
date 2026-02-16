import sys
import subprocess
from pathlib import Path
import pandas as pd
import ast

k = 2
input_file_path = Path("projectData/chr4_200bp_bins.tsv")
tf_name = 'CTCF'

df_global = pd.DataFrame()

for m in range(1):
    command = [sys.executable, '-u', 'final_2.py', '--order', str(m), '--k', str(k), '--input', input_file_path, '--TF', tf_name]
    
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=None, text=True)
    output = result.stdout
    
    try:
        dict_start = output.rfind("{")
        dict_end = output.rfind("}") + 1
        
        if dict_start != -1 and dict_end > dict_start:
            dict_str = output[dict_start:dict_end]
            data_dict = ast.literal_eval(dict_str)
            df_m = pd.DataFrame(data_dict)
            df_global = pd.concat([df_global, df_m], ignore_index=True)
            print(f"Captured data for m={m}")
        else:
            print(f"No dictionary found for m={m}")
            
    except Exception as e:
        print(f"Error parsing m={m}: {e}")

print(df_global)

if not df_global.empty:
    df_global.to_csv(f"ultimate_log_k_{k}_{input_file_path.stem[:4]}_{tf_name}.csv", index=False)