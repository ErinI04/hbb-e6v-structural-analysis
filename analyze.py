"""Reproducible HBB E6V sequence and ColabFold rank-1 comparison.
Run from project root: python src/analyze.py
"""
from pathlib import Path
import json, zipfile, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, Superimposer
from Bio.SeqUtils import seq1
from Bio.Align import PairwiseAligner, substitution_matrices

ROOT=Path(__file__).resolve().parents[1]; DATA=ROOT/'data'; OUT=ROOT/'results'; OUT.mkdir(exist_ok=True)
archives={'WT': DATA/'HBB_WildType_639b6_0.result.zip', 'E6V': DATA/'HBB_E6V_Mutant_b4a43.result.zip'}
structures={}; sequences={}; plddt={}; rank_records=[]
parser=PDBParser(QUIET=True)
for label,path in archives.items():
    with zipfile.ZipFile(path) as z:
        for rank in range(1,6):
            pdb_name=next(n for n in z.namelist() if f'_unrelaxed_rank_{rank:03d}_' in n and n.endswith('.pdb'))
            score_name=next(n for n in z.namelist() if f'_scores_rank_{rank:03d}_' in n and n.endswith('.json'))
            pdb_path=OUT/f'{label}_rank{rank}.pdb';pdb_path.write_bytes(z.read(pdb_name))
            model=parser.get_structure(f'{label}_{rank}',str(pdb_path))
            ca=[r['CA'] for r in model.get_residues() if r.id[0]==' ' and 'CA' in r]
            seq=''.join(seq1(a.get_parent().resname) for a in ca)
            scores=json.loads(z.read(score_name)); confidence=np.asarray(scores['plddt'],dtype=float)
            if len(seq)!=147 or len(confidence)!=147:raise ValueError(f'Unexpected length: {label} rank {rank}')
            rank_records.append({'variant':label,'rank':rank,'model':pdb_name.split('_rank_')[-1].removesuffix('.pdb'),'mean_plddt':float(confidence.mean()),'site7_plddt':float(confidence[6]),'ptm':scores.get('ptm')})
            if rank==1: structures[label]=ca; sequences[label]=seq; plddt[label]=confidence
assert len(sequences['WT'])==len(sequences['E6V'])==147
changes=[(i+1,a,b) for i,(a,b) in enumerate(zip(sequences['WT'],sequences['E6V'])) if a!=b]
assert changes==[(7,'E','V')],f'Unexpected mutation(s): {changes}'
aligner=PairwiseAligner();aligner.mode='global';aligner.substitution_matrix=substitution_matrices.load('BLOSUM62');aligner.open_gap_score=-10;aligner.extend_gap_score=-0.5
sup=Superimposer();sup.set_atoms(structures['WT'],structures['E6V'])
summary={'length':147,'substitutions':changes,'identical_residues':146,'identity_percent':146/147*100,'blosum62_global_score':float(aligner.score(sequences['WT'],sequences['E6V'])),'ca_rmsd_angstrom':float(sup.rms),'wt_mean_plddt':float(plddt['WT'].mean()),'mutant_mean_plddt':float(plddt['E6V'].mean()),'wt_site7_plddt':float(plddt['WT'][6]),'mutant_site7_plddt':float(plddt['E6V'][6])}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2))
with (OUT/'rank_confidence.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=rank_records[0].keys());writer.writeheader();writer.writerows(rank_records)
plt.rcParams.update({'font.size':11})
fig,ax=plt.subplots(figsize=(8,3.1));ax.barh(['WT vs E6V'],[146],label='Identical',color='#1294ac');ax.barh(['WT vs E6V'],[1],left=[146],label='E7V substitution',color='#d87576');ax.set_xlim(0,147);ax.set_xlabel('Aligned amino acid positions (n)');ax.legend(loc='lower center',bbox_to_anchor=(.5,-.65),ncol=2);ax.set_title('Global sequence comparison: 146/147 identical');fig.tight_layout();fig.savefig(OUT/'sequence_identity.png',dpi=200,bbox_inches='tight');plt.close(fig)
fig,ax=plt.subplots(figsize=(9,3.4));x=np.arange(1,148);ax.plot(x,plddt['WT'],label='Wild type',lw=1.5,color='#1294ac');ax.plot(x,plddt['E6V'],label='E6V mutant',lw=1.5,color='#d87576');ax.axvline(7,color='gray',ls=':',lw=1.2,label='Residue 7');ax.set(xlabel='Precursor residue number',ylabel='pLDDT (0–100)',ylim=(0,101),title='Rank-1 AlphaFold confidence by residue');ax.legend(ncol=3,fontsize=9);fig.tight_layout();fig.savefig(OUT/'plddt_profile.png',dpi=200);plt.close(fig)
fig,ax=plt.subplots(figsize=(8,3.5));for_plot=rank_records;xs=np.arange(1,6);ax.plot(xs,[r['mean_plddt'] for r in for_plot if r['variant']=='WT'],'o-',label='Wild type',color='#1294ac');ax.plot(xs,[r['mean_plddt'] for r in for_plot if r['variant']=='E6V'],'o-',label='E6V mutant',color='#d87576');ax.set(xlabel='ColabFold rank (1–5)',ylabel='Mean pLDDT (0–100)',xticks=xs,title='Confidence across all five ranked models');ax.legend();fig.tight_layout();fig.savefig(OUT/'five_model_confidence.png',dpi=200);plt.close(fig)
print(json.dumps(summary,indent=2))
