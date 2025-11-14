"""download ligands in pdbqt format from ZINC database

ZINC tranches are organized as:
- First 2 letters (WEIGHT+LOGP): first subdirectory, e.g., "BB", "CB"
- Next 4 letters (REACT+PURCHASE+PH+CHARGE): second subdirectory, e.g., "ABRP"
- All 6 letters: filename prefix, e.g., "CBABRP"
- Chunks: .xaa.pdbqt.gz, .xab.pdbqt.gz, etc.

This script:
1. Discovers all available chunks for the tranche
2. Downloads each chunk
3. Combines them into a single .pdbqt.gz file (without chunk indicator)
"""

import itertools
import os
import subprocess
import hashlib
import gzip
import sys

# Redirect all stdout/stderr to the log file
sys.stdout = open(snakemake.log[0], 'w', buffering=1)  # line buffering
sys.stderr = sys.stdout

# Use the snakemake output path directly
output_file = str(snakemake.output[0])

# Get ZINC mirror from config
zinc_mirror = snakemake.config.get("ZINC_MIRROR", "files.docking.org")
if not zinc_mirror.startswith("http://") and not zinc_mirror.startswith("https://"):
    zinc_mirror = "https://" + zinc_mirror
zinc_mirror = zinc_mirror.rstrip("/")

# Extract path components from the configuration file
zinc_params = snakemake.config.get("ZINC_INPUT", {})
print(f"ZINC input parameters from config: {zinc_params}")

# construct all download paths from the zinc_params
# The first two letters indicate the tranche sets (all permutations of WEIGHT+LOGP lists)
tranches = [
    "".join(_)
    for _ in itertools.product(
        zinc_params.get("WEIGHT", []), zinc_params.get("LOGP", [])
    )
]
# The following for letters are REACT+PURCHASE+PH+CHARGE
# We need to generate all combinations based on the provided parameters,
# keeping the order consistent with ZINC naming conventions.
subsets = []
for r in itertools.product(
    zinc_params.get("REACT", []),
    zinc_params.get("PURCHASE", []),
    zinc_params.get("PH", []),
    zinc_params.get("CHARGE", []),
):
    subsets.append("".join(r))

print(f"Generated tranches: {tranches}")
print(f"Generated subsets: {subsets}")

print(f"Output file: {output_file}")

# Ensure output directory exists
directory = os.path.dirname(output_file)
os.makedirs(directory, exist_ok=True)

# file to store hashes of all downloaded chunks in a format "<sha256>  <filename>"
hashfile = os.path.join(directory, "hashes.txt")

# memorize all downloaded chunks in the hashfile
downloaded_chunks = []
if os.path.exists(hashfile):
    with open(hashfile, "r") as hf:
        for line in hf:
            parts = line.strip().split(" ")
            if len(parts) == 2:
                downloaded_chunks.append(parts[1])

# Download all compbinations of tranches and subsets
for tranche in tranches:
    chunks = []
    for subset in subsets:
        full_tranche = tranche + subset
        print(f"Processing tranche: {full_tranche}")
        # there may (!) be more then one chunk per tranche
        for chunk in ("xaa", "xab", "xac"):
            chunk_url = (
                f"{zinc_mirror}/3D/{tranche}/{subset}/{full_tranche}.{chunk}.pdbqt.gz"
            )
            chunkfile = os.path.join(directory, f"{full_tranche}.{chunk}.pdbqt.gz")

            # skipt download if chunk file already exists and hashfile has an entry for this file
            if os.path.exists(chunkfile) and chunkfile in downloaded_chunks:
                chunks.append(chunkfile)
                print(f"  Found existing chunk file: {full_tranche}.{chunk}.pdbqt.gz")
                print("  Skipping download...")
                continue

            # Quick HEAD request to check if chunk exists using curl
            cmd = ["curl", "--remote-time", "--fail", "-o", chunkfile, chunk_url]
            proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            if proc.returncode == 0:
                chunks.append(chunkfile)
                # compute sha256 for the chunk and write to .sh2 file
                print(f"  Downloaded: {full_tranche}.{chunk}.pdbqt.gz")
                print("  Computing SHA-256 checksum...")
                hashobj = hashlib.sha256()
                with open(chunkfile, "rb") as f:
                    for c in iter(lambda: f.read(8192), b""):
                        hashobj.update(c)
                digest = hashobj.hexdigest()
                with open(hashfile, "w") as hf:
                    hf.write(f"{digest}  {chunkfile}\n")

                print(f"  Ready: {full_tranche}.{chunk}.pdbqt.gz")
            else:
                # No more chunks found - delete any partial file
                if os.path.exists(chunkfile):
                    os.remove(chunkfile)
                break
        if proc.returncode != 0 and len(chunks) > 0:
            continue
    if proc.returncode != 0 and len(chunks) > 0:
        continue

if not chunks:
    raise RuntimeError(f"No chunks found for tranche {tranche}/{subset}")

print(f"\nCombining {len(chunks)} chunk(s) into {output_file}")

# Combine all chunks into the final output file
# Decompress each chunk, concatenate, then recompress
with gzip.open(output_file, "wb") as outf:
    for chunk_path in chunks:
        print(f"  Reading {os.path.basename(chunk_path)}")
        with gzip.open(chunk_path, "rb") as inf:
            outf.write(inf.read())

        # Clean up chunk file
        os.remove(chunk_path)

print(f"Successfully created {output_file}")

# Compute sha256 and write to the .sh2 file (same format as `sha256sum`)
print("Computing SHA-256 checksum...")
hashed = output_file + ".sh2"
hashobj = hashlib.sha256()
with open(output_file, "rb") as f:
    for chunk in iter(lambda: f.read(8192), b""):
        hashobj.update(chunk)
digest = hashobj.hexdigest()
with open(hashed, "w") as hf:
    hf.write(f"{digest}  {output_file}\n")

print(f"Checksum written to {hashed}")
print("Done!")
