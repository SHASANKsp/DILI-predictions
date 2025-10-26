import pandas as pd
import time
from tqdm import tqdm
from chembl_webresource_client.new_client import new_client
import logging

# ---------------- CONFIG ----------------
INPUT_FILE = "DILIrank.csv"                 # Must have a column 'Name'
OUTPUT_META = "chembl_smiles_and_meta.csv"  # Output 1: SMILES and metadata
OUTPUT_TARGETS = "chembl_targets_mechanisms.csv"  # Output 2: Targets/mechanisms
SLEEP = 0.25     # seconds between queries to be polite to the API
MAX_HITS = 5     # maximum hits to inspect for each name
DEEP_CHECK = 3   # how many top hits to inspect deeply

# --------------- SETUP ------------------
molecule_client = new_client.molecule
mechanism_client = new_client.mechanism
logging.getLogger("chembl_webresource_client").setLevel(logging.WARNING)

# --------------- HELPERS ----------------
def pref_name_exact(name):
    """Try exact preferred-name match (fast, deterministic)."""
    try:
        res = molecule_client.filter(pref_name__iexact=name)
        for r in res:
            return r
    except Exception:
        pass
    return None

def limited_text_search(name):
    """Limited text search (avoid paging through hundreds)."""
    candidates = []
    try:
        res = molecule_client.search(name)
        for i, r in enumerate(res):
            if i >= MAX_HITS:
                break
            candidates.append(r)
    except Exception:
        pass
    return candidates

def choose_best_candidate(candidates, query_name):
    """Pick best chembl_id heuristically."""
    if not candidates:
        return None
    def score(c):
        s = 0
        pref = (c.get('pref_name') or '').lower().strip()
        q = query_name.lower().strip()
        if pref == q:
            s += 200
        if q in pref or pref in q:
            s += 20
        mp = c.get('max_phase')
        try:
            if mp is not None:
                s += int(mp) * 5
        except Exception:
            pass
        return s
    ranked = sorted(candidates, key=score, reverse=True)

    for c in ranked[:DEEP_CHECK]:
        chembl_id = c.get('molecule_chembl_id')
        try:
            full = molecule_client.get(chembl_id)
        except Exception:
            continue
        names = set()
        pref = (full.get('pref_name') or '').lower().strip()
        names.add(pref)
        for syn in full.get('molecule_synonyms') or []:
            if isinstance(syn, dict):
                s = syn.get('synonyms') or syn.get('molecule_synonym') or ""
            else:
                s = str(syn)
            if s:
                names.add(s.lower().strip())
        q = query_name.lower().strip()
        q_simplified = q.split()[0]
        if q in names or q_simplified in names:
            return chembl_id
        time.sleep(SLEEP)
    return ranked[0].get('molecule_chembl_id')

def fetch_molecule_data(chembl_id):
    """Fetch SMILES, InChIKey, molecule type."""
    try:
        mol = molecule_client.get(chembl_id)
    except Exception:
        return None, None, None
    mol_type = mol.get('molecule_type')
    struct = mol.get('molecule_structures') or {}
    smiles = struct.get('canonical_smiles')
    inchi_key = struct.get('standard_inchi_key')
    return smiles, inchi_key, mol_type

def fetch_targets(chembl_id):
    """Fetch mechanism/target data for a molecule."""
    results = []
    try:
        mechs = mechanism_client.filter(molecule_chembl_id=chembl_id)
        for m in mechs:
            results.append({
                "molecule_chembl_id": chembl_id,
                "target_chembl_id": m.get('target_chembl_id'),
                "target_pref_name": m.get('target_pref_name'),
                "target_organism": m.get('target_organism'),
                "action_type": m.get('action_type'),
                "mechanism_of_action": m.get('mechanism_of_action')
            })
    except Exception:
        pass
    return results

# --------------- MAIN -------------------
df = pd.read_csv(INPUT_FILE)
if "Name" not in df.columns:
    raise ValueError("Input CSV must contain a 'Name' column.")

names = df["Name"].dropna().unique().tolist()
meta_records = []
target_records = []

for name in tqdm(names, desc="Fetching from ChEMBL"):
    name = str(name).strip()
    if not name:
        continue

    # Try exact pref_name
    candidate = pref_name_exact(name)
    time.sleep(SLEEP)
    if candidate:
        candidates = [candidate]
    else:
        candidates = limited_text_search(name)

    if not candidates:
        meta_records.append({
            "query_name": name,
            "chembl_id": None,
            "smiles": None,
            "inchi_key": None,
            "molecule_type": None,
            "status": "no_hit"
        })
        continue

    chembl_id = choose_best_candidate(candidates, name)
    if not chembl_id:
        meta_records.append({
            "query_name": name,
            "chembl_id": None,
            "smiles": None,
            "inchi_key": None,
            "molecule_type": None,
            "status": "no_candidate"
        })
        continue

    smiles, inchi_key, mol_type = fetch_molecule_data(chembl_id)
    time.sleep(SLEEP)
    is_small = mol_type and mol_type.lower() == "small molecule" and smiles
    status = "ok" if is_small else "non_small_or_no_smiles"

    meta_records.append({
        "query_name": name,
        "chembl_id": chembl_id,
        "smiles": smiles,
        "inchi_key": inchi_key,
        "molecule_type": mol_type,
        "status": status
    })

    # Fetch targets only for small molecules
    if is_small:
        targets = fetch_targets(chembl_id)
        time.sleep(SLEEP)
        for t in targets:
            t["query_name"] = name
            target_records.append(t)

# --------------- SAVE OUTPUTS ---------------
meta_df = pd.DataFrame(meta_records)
target_df = pd.DataFrame(target_records)

meta_df.to_csv(OUTPUT_META, index=False)
target_df.to_csv(OUTPUT_TARGETS, index=False)

ok_count = (meta_df["status"] == "ok").sum()
skipped = len(meta_df) - ok_count
print(f"\n✅ Done. SMILES found for {ok_count} small molecules; skipped {skipped} non-small/no-hit entries.")
print(f"→ Saved metadata: {OUTPUT_META} ({len(meta_df)} rows)")
print(f"→ Saved targets:  {OUTPUT_TARGETS} ({len(target_df)} rows)")
