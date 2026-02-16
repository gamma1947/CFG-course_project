import sys
import subprocess
from pathlib import Path
import pandas as pd
import datetime # <--- REQUIRED: Import this so eval() understands the time objects

k = 2
input_file_path = Path("projectData/chr4_200bp_bins.tsv")
tf_name = 'CTCF'

df_global = pd.DataFrame()

for m in range(1):
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
