import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("TOKENIZERS_PARALLELISM", "false")
os.environ.setdefault("PYTHONFAULTHANDLER", "1")

import multiprocessing as mp
import safe as sf
from safe.utils import compute_side_chains
from rdkit import Chem
from rdkit.Chem.Scaffolds import rdScaffoldNetwork
from rdkit.Chem import Descriptors
from itertools import combinations

def _dbg(i, **kv):
  try:
    items = " ".join(f"{k}={v}" for k, v in kv.items())
  except Exception:
    items = ""
  print(f"[DBG] {i} {items}", flush=True)

def _torch_thread_guard():
  try:
    import torch
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
  except Exception:
    pass

def _morph_worker(side_chains, n_trials, n_samples_per_trial, q):
  os.environ.setdefault("OMP_NUM_THREADS", "1")
  os.environ.setdefault("MKL_NUM_THREADS", "1")
  os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
  os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
  os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
  os.environ.setdefault("TOKENIZERS_PARALLELISM", "false")
  os.environ.setdefault("PYTHONFAULTHANDLER", "1")
  _torch_thread_guard()
  import safe as sf_local
  designer = sf_local.SAFEDesign.load_default(verbose=False)
  try:
    out = designer.scaffold_morphing(
      side_chains=side_chains,
      n_samples_per_trial=n_samples_per_trial,
      n_trials=n_trials,
      sanitize=True,
      do_not_fragment_further=False,
      random_seed=100,
      num_beams=4,
      early_stopping=True,
    )
    q.put(("ok", out))
  except Exception as e:
    q.put(("err", repr(e)))

def safe_scaffold_morphing(side_chains, n_trials, n_samples_per_trial, retries=2, soft_fail=True):
  attempts = [
    (n_trials, n_samples_per_trial),
    (max(1, n_trials // 2), max(1, n_samples_per_trial // 2)),
    (1, max(1, n_samples_per_trial // 4)),
  ][: 1 + retries]

  for attempt_idx, (t, s) in enumerate(attempts):
    _dbg("morph:attempt", attempt=attempt_idx + 1, n_trials=t, n_samples=s)
    q = mp.Queue()
    p = mp.Process(target=_morph_worker, args=(side_chains, t, s, q))
    p.start()
    p.join()

    if p.exitcode == 0 and not q.empty():
      tag, payload = q.get()
      if tag == "ok":
        _dbg("morph:ok", n_out=len(payload) if payload is not None else None)
        return payload or []
      _dbg("morph:err", payload=payload)
    else:
      _dbg("morph:crash", exitcode=p.exitcode)

  if soft_fail:
    _dbg("morph:soft_fail_return_empty")
    return []

  raise RuntimeError("scaffold_morphing crashed repeatedly (likely native segfault/OOM).")

class MoleculeModel:
  def __init__(self, n_trials=10, n_samples_per_trial=100, lower_molecular_weight=60, upper_molecular_weight=100):
    _dbg("init:start", n_trials=n_trials, n_samples_per_trial=n_samples_per_trial)
    self.designer = sf.SAFEDesign.load_default(verbose=True)
    _dbg("init:designer_loaded")
    self.n_trials = n_trials
    self.n_samples_per_trial = n_samples_per_trial
    self.lower_molecular_weight = lower_molecular_weight
    self.upper_molecular_weight = upper_molecular_weight
    _dbg("init:done")

  def smiles_to_safe(self, smiles):
    _dbg("smiles_to_safe:in", smiles_len=len(smiles) if smiles is not None else None)
    try:
      out = sf.encode(smiles)
      _dbg("smiles_to_safe:out", out_len=len(out) if out is not None else None)
      return out
    except Exception as e:
      print(f"[ERR] smiles_to_safe {e}", flush=True)
      return None

  def _extract_core_structure(self, safe_str):
    _dbg("_extract_core_structure:in", safe_is_none=safe_str is None, safe_type=type(safe_str).__name__)
    params = rdScaffoldNetwork.ScaffoldNetworkParams()
    params.includeScaffoldsWithoutAttachments = False

    if safe_str is None:
      _dbg("_extract_core_structure:skip_none")
      return []

    _dbg("_extract_core_structure:mol_from_smiles:start")
    mol = Chem.MolFromSmiles(safe_str)
    _dbg("_extract_core_structure:mol_from_smiles:done", mol_is_none=mol is None)
    if mol is None:
      _dbg("_extract_core_structure:bad_smiles")
      return []

    _dbg("_extract_core_structure:create_network:start")
    net = rdScaffoldNetwork.CreateScaffoldNetwork([mol], params)
    _dbg("_extract_core_structure:create_network:done", n_nodes=len(net.nodes))

    _dbg("_extract_core_structure:nodemols:start")
    nodemols = [Chem.MolFromSmiles(x) for x in net.nodes]
    nodemols = [m for m in nodemols if m is not None]
    _dbg("_extract_core_structure:nodemols:done", n_nodemols=len(nodemols))
    if not nodemols:
      return []

    filtered_list = []
    _dbg("_extract_core_structure:filter:start")
    for idx, m in enumerate(nodemols):
      try:
        smi = Chem.MolToSmiles(m)
        wt = Descriptors.MolWt(m)
        if "*" in smi and self.lower_molecular_weight < wt < self.upper_molecular_weight:
          filtered_list.append(m)
      except Exception as e:
        print(f"[ERR] filter idx={idx} {e}", flush=True)
    _dbg("_extract_core_structure:filter:done", n_filtered=len(filtered_list))

    if not filtered_list:
      _dbg("_extract_core_structure:closest:start")
      target = (self.lower_molecular_weight + self.upper_molecular_weight) / 2
      closest_mol = min(nodemols, key=lambda x: abs(Descriptors.MolWt(x) - target))
      filtered_list.append(closest_mol)
      _dbg("_extract_core_structure:closest:done", target=target, closest_wt=Descriptors.MolWt(closest_mol))

    _dbg("_extract_core_structure:sort:start")
    filtered_list.sort(key=lambda x: x.GetNumHeavyAtoms())
    _dbg("_extract_core_structure:sort:done", first_wt=Descriptors.MolWt(filtered_list[0]) if filtered_list else None)

    return filtered_list

  def _get_side_chain_pairs(self, side_chains):
    _dbg("_get_side_chain_pairs:in", side_chains_is_none=side_chains is None, side_chains_type=type(side_chains).__name__)
    if side_chains is None:
      return []

    smi = Chem.MolToSmiles(side_chains)
    parts = smi.split(".") if smi else []
    _dbg("_get_side_chain_pairs:parts", n_parts=len(parts))

    comb2 = list(combinations(parts, 2))
    comb3 = list(combinations(parts, 3))
    allc = comb2 + comb3
    _dbg("_get_side_chain_pairs:combinations", n2=len(comb2), n3=len(comb3), n_all=len(allc))

    modified = []
    for tup in allc:
      ms = []
      for index, j in enumerate(tup):
        if len(j) >= 2:
          ms.append(j[0] + str(index + 1) + j[2:])
        elif len(j) == 1:
          ms.append(j[0] + str(index + 1))
        else:
          ms.append(str(index + 1))
      modified.append(ms)

    joined = [".".join(x) for x in modified]
    _dbg("_get_side_chain_pairs:out", n_joined=len(joined))
    return joined

  def _generate_smiles(self, side_chains):
    _dbg("_generate_smiles:in", side_chains=side_chains, n_trials=self.n_trials, n_samples=self.n_samples_per_trial)
    out = safe_scaffold_morphing(
      side_chains=side_chains,
      n_trials=self.n_trials,
      n_samples_per_trial=self.n_samples_per_trial,
      retries=2,
      soft_fail=True,
    )
    _dbg("_generate_smiles:out", n_out=len(out) if out is not None else None)
    return out

  def run_model(self, safe_list):
    _dbg("run_model:start", safe_type=type(safe_list).__name__, safe_len=len(safe_list) if safe_list is not None else None)
    generated_smiles = []

    for idx, i in enumerate(safe_list):
      _dbg("run_model:item", idx=idx, item_is_none=i is None, item_type=type(i).__name__)
      row = []

      if i is None:
        generated_smiles.append(row)
        continue

      _dbg("run_model:extract_core:start", idx=idx)
      core_structures = self._extract_core_structure(i)
      _dbg("run_model:extract_core:done", idx=idx, n_cores=len(core_structures))

      for cidx, core in enumerate(core_structures):
        _dbg("run_model:side_chains:start", idx=idx, core_idx=cidx)
        try:
          side_chain = compute_side_chains(core=core, mol=i)
          _dbg("run_model:side_chains:done", idx=idx, core_idx=cidx, side_chain_is_none=side_chain is None, side_chain_type=type(side_chain).__name__)
        except Exception as e:
          print(f"[ERR] compute_side_chains idx={idx} core_idx={cidx} {e}", flush=True)
          continue

        side_chain_pairs = self._get_side_chain_pairs(side_chain)
        _dbg("run_model:pairs", idx=idx, core_idx=cidx, n_pairs=len(side_chain_pairs))

        for pidx, sc in enumerate(side_chain_pairs):
          _dbg("run_model:generate:start", idx=idx, core_idx=cidx, pair_idx=pidx, pair=sc)
          out = self._generate_smiles(sc)
          _dbg("run_model:generate:done", idx=idx, core_idx=cidx, pair_idx=pidx, out_n=len(out) if out is not None else None)
          if out:
            row += out

      _dbg("run_model:row_done", idx=idx, row_len=len(row))
      generated_smiles.append(row)

    _dbg("run_model:done", n_rows=len(generated_smiles))
    return generated_smiles
