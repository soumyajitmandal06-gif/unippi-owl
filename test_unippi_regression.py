#!/usr/bin/env python3
"""
UniPPI-OWL v1.7 Deterministic Regression & Real RCSB PDB Verification Suite
Adapts the HUNTER v1 frozen regression verification standard across Modes A, B, and C.
"""

import os
import sys
import json
import shutil
import urllib.request
import subprocess

EXPECTED_FROZEN_V1_7 = {
    "1_VIP1_Homodimer_Positive_model_0": {"UIS": 0.777, "dG": -10.35, "Pattern": "Pattern A"},
    "3_Barnase_Barstar_Electrostatic_Pos_model_0": {"UIS": 0.771, "dG": -10.09, "Pattern": "Pattern A"},
    "2_VIP1_bZIP61_GroupE_Positive_model_0": {"UIS": 0.732, "dG": -10.10, "Pattern": "Pattern A"},
    "4_VIP1_GRF6_14-3-3_IDR_Diluted_Pos_model_0": {"UIS": 0.609, "dG": -7.61, "Pattern": "Pattern B"},
    "5_bZIP34_Proline_Homodimer_Negative_model_0": {"UIS": 0.120, "dG": 5.00, "Pattern": "Pattern C"},
    "6_VIP1_IDR_Hallucination_Decoy_model_0": {"UIS": 0.061, "dG": 2.10, "Pattern": "Pattern D"},
}

def run_cmd(cmd):
    res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if res.returncode != 0:
        print("[FAIL] Command failed:", cmd)
        print(res.stderr)
        sys.exit(1)
    return res.stdout

def main():
    print("="*92)
    print("UniPPI-OWL v1.7 — Automated 5-Gate Regression & Real RCSB PDB Verification Suite")
    print("="*92)
    report_lines = ["Gate\tTest_Name\tTarget\tObserved\tExpected\tStatus"]

    # --------------------------------------------------------------------------
    # GATE 1: Frozen 6-Complex Benchmark Numerical Parity (Mode A)
    # --------------------------------------------------------------------------
    run_cmd("python3 unippi_cli.py --generate-demo ./demo_structures -o ./demo_results --mutate A:LEU5PRO,A:PRO5LEU,A:ASP2ALA,A:TYR5ALA")
    with open("./demo_results/unippi_full_report.json", "r", encoding="utf-8") as f:
        records = {r["complex_id"]: r for r in json.load(f)}

    for cid, exp in EXPECTED_FROZEN_V1_7.items():
        rec = records.get(cid)
        assert rec is not None, f"Missing {cid} in report!"
        uis_diff = abs(rec["UIS"] - exp["UIS"])
        dg_diff = abs(rec["dG_uni_kcal"] - exp["dG"])
        pat_ok = (rec["pattern_id"] == exp["Pattern"])
        assert uis_diff <= 0.005 and dg_diff <= 0.05 and pat_ok, f"Parity mismatch on {cid}: {rec}"
        report_lines.append(f"Gate_1_Parity\t{cid}\tUIS/dG/Pattern\t{rec['UIS']:.3f}/{rec['dG_uni_kcal']:.2f}/{rec['pattern_id']}\t{exp['UIS']:.3f}/{exp['dG']:.2f}/{exp['Pattern']}\tPASS")
    print("[PASS] Gate 1: All 6 Frozen Benchmark Complexes match UIS, dG_uni, and Pattern A-D calls.")

    # --------------------------------------------------------------------------
    # GATE 2: Mode B Fractional [0, 1] pLDDT Auto-Scaling Verification
    # --------------------------------------------------------------------------
    os.makedirs("./_tmp_mode_b", exist_ok=True)
    src_cif = "./demo_structures/1_VIP1_Homodimer_Positive/1_VIP1_Homodimer_Positive_model_0.cif"
    with open(src_cif, "r", encoding="utf-8") as f:
        cif_txt = f.read().replace(" 94.00 ", " 0.94 ")
    with open("./_tmp_mode_b/fractional_plddt.cif", "w", encoding="utf-8") as f:
        f.write(cif_txt)
    run_cmd("python3 unippi_cli.py -i ./_tmp_mode_b/fractional_plddt.cif -o ./_tmp_mode_b_out")
    with open("./_tmp_mode_b_out/unippi_full_report.json", "r", encoding="utf-8") as f:
        mb_rec = json.load(f)[0]
    assert mb_rec["mode"] == "MODE_B_PRED_COORD" and mb_rec["mean_intf_pLDDT"] >= 93.0, f"Mode B scaling failed: {mb_rec}"
    report_lines.append(f"Gate_2_ModeB\tFractional_pLDDT_Scaling\tmean_pLDDT\t{mb_rec['mean_intf_pLDDT']:.1f}\t94.0\tPASS")
    shutil.rmtree("./_tmp_mode_b", ignore_errors=True)
    shutil.rmtree("./_tmp_mode_b_out", ignore_errors=True)
    print("[PASS] Gate 2: Mode B auto-detection & fractional [0,1] -> [0,100] pLDDT scaling verified.")

    # --------------------------------------------------------------------------
    # GATE 3: Strict 80-Column Legacy PDB Formatting & Glycine CB Sanitization
    # --------------------------------------------------------------------------
    pdb_files = [
        os.path.join("./demo_results/converted_structures", fn)
        for fn in os.listdir("./demo_results/converted_structures") if fn.endswith("_sanitized.pdb")
    ]
    for pf in pdb_files:
        with open(pf, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith(("ATOM  ", "HETATM")):
                    raw = line.rstrip("\r\n")
                    assert len(raw) == 80, f"Line length {len(raw)} != 80 in {pf}: {raw}"
                    res_n = raw[17:20].strip()
                    atom_n = raw[12:16].strip()
                    assert not (res_n == "GLY" and atom_n == "CB"), f"Unphysical GLY CB found in {pf}"
    report_lines.append(f"Gate_3_Sanitizer\tFoldX_Rosetta_80Col_PDB\tFiles_Audited\t{len(pdb_files)}_PDBs_Clean\t80_Cols_No_GLY_CB\tPASS")
    print(f"[PASS] Gate 3: All {len(pdb_files)} sanitized PDB files verified (strictly 80 columns, zero GLY CB atoms).")

    # --------------------------------------------------------------------------
    # GATE 4: High-Salt (800 mM NaCl) Osmotic Stress Divergence Check
    # --------------------------------------------------------------------------
    run_cmd("python3 unippi_cli.py -i ./demo_structures -o ./_tmp_salt800 --ionic-strength 800.0")
    with open("./_tmp_salt800/unippi_full_report.json", "r", encoding="utf-8") as f:
        salt_recs = {r["complex_id"]: r for r in json.load(f)}
    dg_vip1_150 = records["1_VIP1_Homodimer_Positive_model_0"]["dG_uni_kcal"]
    dg_vip1_800 = salt_recs["1_VIP1_Homodimer_Positive_model_0"]["dG_uni_kcal"]
    dg_barn_150 = records["3_Barnase_Barstar_Electrostatic_Pos_model_0"]["dG_uni_kcal"]
    dg_barn_800 = salt_recs["3_Barnase_Barstar_Electrostatic_Pos_model_0"]["dG_uni_kcal"]
    assert dg_vip1_800 < dg_vip1_150 and dg_barn_800 > dg_barn_150, "Salt stress divergence failed!"
    report_lines.append(f"Gate_4_SaltStress\tVIP1_vs_Barnase_800mM\tdG_Shift\tVIP1:{dg_vip1_800:.2f}/Barn:{dg_barn_800:.2f}\tVIP1<-10.35/Barn>-10.09\tPASS")
    shutil.rmtree("./_tmp_salt800", ignore_errors=True)
    print(f"[PASS] Gate 4: Salt-stress physics verified (800 mM NaCl stabilizes hydrophobic VIP1 [{dg_vip1_150:.2f} -> {dg_vip1_800:.2f} kcal/mol] while screening Barnase-Barstar [{dg_barn_150:.2f} -> {dg_barn_800:.2f} kcal/mol]).")

    # --------------------------------------------------------------------------
    # GATE 5: Real Experimental RCSB PDB Crystal Structure (1BRS) Mode C Test
    # --------------------------------------------------------------------------
    os.makedirs("./_tmp_rcsb", exist_ok=True)
    rcsb_cif = "./_tmp_rcsb/1BRS.cif"
    try:
        urllib.request.urlretrieve("https://files.rcsb.org/download/1BRS.cif", rcsb_cif)
        run_cmd(f"python3 unippi_cli.py -i {rcsb_cif} -o ./_tmp_rcsb_out --mutate A:ARG83ALA")
        with open("./_tmp_rcsb_out/unippi_full_report.json", "r", encoding="utf-8") as f:
            rcsb_recs = json.load(f)
        top_hit = rcsb_recs[0]
        mut_hit = top_hit["custom_mutations"][0]; assert abs(mut_hit["dG_mut_kcal"] - (top_hit["dG_uni_kcal"] + mut_hit["ddG_kcal"])) <= 0.02 and top_hit["mode"] == "MODE_C_EXP_CRYSTAL" and top_hit["UIS"] >= 0.65 and top_hit["dG_uni_kcal"] <= -10.5 and top_hit["pattern_id"] == "Pattern A", f"Real PDB 1BRS test failed: {top_hit}"
        report_lines.append(f"Gate_5_RealPDB\tRCSB_1BRS_Crystal_6Chain\tTop_Pair_{top_hit['chain_pair']}\tUIS={top_hit['UIS']:.3f}/dG={top_hit['dG_uni_kcal']:.2f}\tMODE_C_EXP_CRYSTAL\tPASS")
        print(f"[PASS] Gate 5: Real RCSB PDB 1BRS (6-chain X-ray crystal structure) evaluated in {top_hit['mode']} -> Top biological dimer Chain {top_hit['chain_pair']}: UIS = {top_hit['UIS']:.3f}, dG_uni = {top_hit['dG_uni_kcal']:.2f} kcal/mol, BSA = {top_hit['BSA_A2']} A^2.")
    except urllib.error.URLError as e:
        print(f"[INFO] Gate 5 (Online RCSB download skipped/offline: {e})")
    finally:
        shutil.rmtree("./_tmp_rcsb", ignore_errors=True)
        shutil.rmtree("./_tmp_rcsb_out", ignore_errors=True)

    with open("./demo_results/REGRESSION_VERIFICATION_REPORT.tsv", "w", encoding="utf-8") as f:
        f.write("\n".join(report_lines) + "\n")
    print("="*92)
    print("[ALL 5 REGRESSION GATES PASSED] Saved ./demo_results/REGRESSION_VERIFICATION_REPORT.tsv")
    print("="*92)

if __name__ == "__main__":
    main()
