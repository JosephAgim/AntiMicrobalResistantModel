import os
import time
import urllib.request
import urllib.error

# Ensure the directory exists
directory = 'Directory' 
if not os.path.exists(directory):
    os.makedirs(directory)

# Read PDB IDs from file
with open("Directory") as file:
    pdb_ids = file.read().strip().split(',')

print(f'Read {len(pdb_ids)} PDB IDs from file.')

# Download each PDB file
for pdb_id in pdb_ids:
    url = f'https://files.rcsb.org/download/{pdb_id}.pdb.gz'
    filename = f'{directory}/{pdb_id}.pdb.gz'  # Adjusted filename
    print(f'Trying to download PDB ID {pdb_id}...')
    
    try:
        urllib.request.urlretrieve(url, filename)
        print(f'Successfully downloaded {pdb_id}')
    except urllib.error.HTTPError as e:
        if e.code == 403:
            print(f'Failed to download {pdb_id}: Server refused the request (HTTP 403).')
        elif e.code == 404:
            print(f'Failed to download {pdb_id}: File not found on server (HTTP 404).')
        else:
            print(f'Failed to download {pdb_id}: HTTP error {e.code}.')
    except Exception as e:
        print(f'Failed to download {pdb_id}: {e}')

    time.sleep(0.5)  # wait for 0.5 second
