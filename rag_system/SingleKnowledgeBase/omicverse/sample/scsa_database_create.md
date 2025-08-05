# scsa_database_create
*Converted from: omicverse/sample/scsa_database_create.ipynb*

```python#Run in pandas==1.5.3
import gzip
import pickle
import os
import pandas as pd

# Function to save variables as plain text files or CSV
def save_variable_as_txt(name, var):
    """
    Save any variable as a plain text file or CSV if it's a DataFrame.
    """
    os.makedirs('./db2plaintext', exist_ok=True)  # Create directory if it doesn't exist
    
    if isinstance(var, dict):  # For dictionaries, save as key-value pairs
        with open(f'./db2plaintext/{name}.txt', 'w') as f:
            for key, value in var.items():
                f.write(f"{key}: {value}\n")
    
    elif isinstance(var, list):  # For lists, save one item per line
        with open(f'./db2plaintext/{name}.txt', 'w') as f:
            for item in var:
                f.write(f"{item}\n")
    
    elif isinstance(var, pd.DataFrame):  # For DataFrames, save as CSV
        var.to_csv(f'./db2plaintext/{name}.csv', index=False)
    
    else:  # For other types, save as string representation
        with open(f'./db2plaintext/{name}.txt', 'w') as f:
            f.write(str(var))
    
    print(f"{name} saved.")

# Function to save DataFrames in a list as separate CSV files
def save_list_as_csv(list_name, data_list):
    """
    Save each DataFrame in the list as a separate CSV file.
    """
    os.makedirs('./db2plaintext', exist_ok=True)  # Create directory if it doesn't exist
    
    for i, df in enumerate(data_list):
        df.to_csv(f'./db2plaintext/{list_name}_{i + 1}.csv', index=False)
    
    print(f"{list_name} dataframes saved as separate CSV files.")

# Function to read the db file in Pandas 1.5.3 environment and save all variables
def read_and_save_db_pandas1_as_txt(db_file):
    """
    Read the db file and save all variables, including lists, in plain text (*.txt) format.
    """
    with gzip.open(db_file, 'rb') as f:
        gos = pickle.load(f)
        human_gofs = pickle.load(f)  # Expecting a list of 3 DataFrames
        mouse_gofs = pickle.load(f)  # Expecting a list of 3 DataFrames
        cmarkers = pickle.load(f)
        smarkers = pickle.load(f)
        snames = pickle.load(f)
        ensem_hgncs = pickle.load(f)
        ensem_mouse = pickle.load(f)
        if 'plus' in db_file:
            pmarkers = pickle.load(f)
        else:
            pmarkers = self.smarkers.copy()

    # Save each variable as a plain text file
    save_variable_as_txt('gos', gos)
    save_variable_as_txt('cmarkers', cmarkers)
    save_variable_as_txt('smarkers', smarkers)
    save_variable_as_txt('snames', snames)
    save_variable_as_txt('ensem_hgncs', ensem_hgncs)
    save_variable_as_txt('ensem_mouse', ensem_mouse)
    save_variable_as_txt('pmarkers', pmarkers)

    # Save human_gofs and mouse_gofs as separate CSV files
    save_list_as_csv('human_gofs', human_gofs)
    save_list_as_csv('mouse_gofs', mouse_gofs)

    print("All variables saved in './db2plaintext'.")

# Call the function
read_and_save_db_pandas1_as_txt('temp/pySCSA_2023_v2_plus.db')```
*Output:*
```gos saved.
cmarkers saved.
smarkers saved.
snames saved.
ensem_hgncs saved.
ensem_mouse saved.
pmarkers saved.
human_gofs dataframes saved as separate CSV files.
mouse_gofs dataframes saved as separate CSV files.
All variables saved in './db2plaintext'.
```

```python#Run in pandas>2
import pandas as pd
import gzip
import pickle
import os

# Function to load variables from the plain text (*.txt) files and CSV files
def load_from_txt_and_csv():
    """
    Load variables from plain text and CSV files in the './db2plaintext' directory.
    """
    # Load dictionary-like text files (key-value pairs)
    gos = {}
    with open('./db2plaintext/gos.txt', 'r') as f:
        for line in f:
            key, value = line.strip().split(': ', 1)
            gos[key] = value

    ensem_hgncs = {}
    with open('./db2plaintext/ensem_hgncs.txt', 'r') as f:
        for line in f:
            key, value = line.strip().split(': ', 1)
            ensem_hgncs[key] = value

    ensem_mouse = {}
    with open('./db2plaintext/ensem_mouse.txt', 'r') as f:
        for line in f:
            key, value = line.strip().split(': ', 1)
            ensem_mouse[key] = value

    snames = {}
    with open('./db2plaintext/snames.txt', 'r') as f:
        for line in f:
            key, value = line.strip().split(': ', 1)
            snames[key] = value

    # Load lists from CSV files
    human_gofs = [pd.read_csv(f'./db2plaintext/human_gofs_{i + 1}.csv') for i in range(3)]
    mouse_gofs = [pd.read_csv(f'./db2plaintext/mouse_gofs_{i + 1}.csv') for i in range(3)]

    # Load DataFrame-like text files (CSV formatted but saved as .txt)
    cmarkers = pd.read_csv('./db2plaintext/cmarkers.csv')
    smarkers = pd.read_csv('./db2plaintext/smarkers.csv')
    pmarkers = pd.read_csv('./db2plaintext/pmarkers.csv')

    return (gos, ensem_hgncs, ensem_mouse, snames, cmarkers, smarkers, pmarkers, human_gofs, mouse_gofs)

# Function to save the loaded variables to a new pickle file
def save_to_pickle(db_file, gos, ensem_hgncs, ensem_mouse, snames, cmarkers, smarkers,pmarkers, human_gofs, mouse_gofs):
    """
    Save the loaded variables to a new db file as a pickle.
    """
    with gzip.open(db_file, "wb") as handler:
        pickle.dump(gos, handler)
        pickle.dump(human_gofs, handler)
        pickle.dump(mouse_gofs, handler)
        pickle.dump(cmarkers, handler)
        pickle.dump(smarkers, handler)
        pickle.dump(snames, handler)
        pickle.dump(ensem_hgncs, handler)
        pickle.dump(ensem_mouse, handler)
        pickle.dump(pmarkers, handler)

    print(f"Database saved to {db_file} in pickle format.")

# Example usage
def main():
    # Load variables from txt and CSV files
    gos, ensem_hgncs, ensem_mouse, snames, cmarkers, smarkers, pmarkers, human_gofs, mouse_gofs = load_from_txt_and_csv()
    
    # Save to a new pickle file
    save_to_pickle('temp/pySCSA_2024_v1_plus.db', gos, ensem_hgncs, ensem_mouse, snames, cmarkers, smarkers, pmarkers,human_gofs, mouse_gofs)

# Run the main function
main()```
*Output:*
```Database saved to temp/pySCSA_2024_v1_plus.db in pickle format.
```

