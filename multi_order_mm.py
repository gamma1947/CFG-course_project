import sys
import subprocess
from pathlib import Path
import pandas as pd
<<<<<<< HEAD
import datetime # <--- REQUIRED: Import this so eval() understands the time objects
=======
import ast
>>>>>>> refs/remotes/origin/main

k = 2
input_file_path = Path("projectData/chr4_200bp_bins.tsv")
tf_name = 'CTCF'

df_global = pd.DataFrame()

for m in range(1):
<<<<<<< HEAD
    print(f"\n--- Starting subprocess for m={m} ---")
    
    command = [
        sys.executable, '-u', 'final_2.py', 
        '--order', str(m), 
        '--k', str(k), 
        '--input', str(input_file_path), 
        '--TF', tf_name
    ]
    
    output_buffer = []

    # Run the subprocess
    with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=None, text=True, bufsize=1) as proc:
        for line in proc.stdout:  # type: ignore
            print(line, end='') 
            output_buffer.append(line)
            
    full_output = "".join(output_buffer)
    
    try:
        dict_start = full_output.rfind("{")
        dict_end = full_output.rfind("}") + 1
        
        if dict_start != -1 and dict_end > dict_start:
            dict_str = full_output[dict_start:dict_end]
            
            # --- FIX IS HERE ---
            # ast.literal_eval fails on datetime objects. 
            # We use eval() instead, which executes the string as code.
            data_dict = eval(dict_str)
            
            df_m = pd.DataFrame(data_dict)
            df_global = pd.concat([df_global, df_m], ignore_index=True)
            print(f"\n[Success] DataFrame appended for m={m}")
        else:
            print(f"\n[Error] No dictionary structure found in output for m={m}")

    except Exception as e:
        print(f"\n[Exception] Parsing failed: {e}")

print("\nFinal Global DataFrame:")
print(df_global)

if not df_global.empty:
    output_filename = f"ultimate_log_k_{k}_{input_file_path.stem}_{tf_name}.txt"
    df_global.to_csv(output_filename, sep='\t', index=False)
    print(f"Saved to {output_filename}")
=======
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
>>>>>>> refs/remotes/origin/main
