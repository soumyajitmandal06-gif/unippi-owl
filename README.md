# UniPPI-OWL v1.7

**Universal Structure-to-Affinity, Interface Auditing, and Mutational Screening Resource**

* **Author:** Soumyajit Mandal (Department of Biochemistry, Central University of Rajasthan, India)
* **Model Version:** `UniPPI-v1.7` (Frozen Decision Threshold: `UIS = 0.5000`; High-Confidence Cutoff: `UIS = 0.6500`)
* **Architecture:** Dual-Layer — Zero-Dependency Python CLI (`unippi_cli.py`) + Client-Side Browser Application (`index.html`)
* **License:** Open-Source MIT / BSD 3-Clause

---

## 1. Overview and Scientific Objective

### What This Tool Does (In Plain Words)
Structure prediction models like AlphaFold3 (AF3), AlphaFold2, ESM3, and Boltz-1 predict the 3D shape of protein complexes (`.cif` or `.pdb` files) and output geometric confidence scores (`ipTM`, `pLDDT`, and `PAE` matrices). However, geometric confidence is not the same as physical binding strength. AlphaFold3 does not calculate thermodynamic binding free energy (`dG` in kcal/mol) or dissociation constants (`Kd` in nM), nor does it tell you how a point mutation (`ddG`) or a change in salt concentration (`[NaCl]`) and `pH` will affect binding. In addition, running external energy programs like FoldX or Rosetta requires manually converting AlphaFold3's `.cif` files into strict legacy 80-column `.pdb` files.

**UniPPI-OWL** solves these practical bottlenecks in a single workflow:
1. **Automated Format Sanitization:** Converts multi-chain AlphaFold3 `.cif` and `.pdb` structures into strict 80-column FoldX/Rosetta-ready `.pdb` files, extracts non-polymer small-molecule ligands into V2000 `.sdf` files, and outputs FoldX `individual_list.txt` mutation files.
2. **Structure-to-Affinity Estimation:** Calculates confidence-weighted atomic contacts, solvent-accessible surface burial (`wBSA`), binding free energy (`dG_uni`), and dissociation constants (`Kd`) in less than 0.1 seconds per complex.
3. **Disordered-Loop (IDR) Artifact Filtering:** Prevents low-confidence intrinsically disordered tails (common in plant and eukaryotic transcription factors) from generating false-positive binding scores, while rescuing true local interfaces whose global `ipTM` is diluted by long non-interacting tails.
4. **In-Silico Mutagenesis and Environmental Stress Modeling:** Evaluates single-site 20-amino-acid saturation mutagenesis (`ddG`), double-mutant epistatic cycles (`ddG_int`), and ionic strength (`10-1000 mM NaCl`) and `pH` (`4.0-9.5`) titration curves.

### Dual-Layer Architecture
To serve both computational bioinformaticians and wet-lab biologists, UniPPI-OWL is built in two synchronized layers that share identical mathematical equations and produce identical numerical outputs:
* **Layer 1 — Automated Python CLI (`unippi_cli.py`):** Written entirely in the Python 3 standard library (`math.dist` paired with 3D spatial voxel hashing) with zero external package dependencies (`NumPy`, `Biopython`, or C++ compilers are not required). Designed for batch folder screening and workflow managers such as Nextflow or Snakemake.
* **Layer 2 — Client-Side Browser Application (`index.html`):** A zero-server HTML5, Vanilla JavaScript, and WebGL (`3Dmol.js`) Single-Page Application (SPA) following the `-OWL` (*On the Web Lite*) standard. All calculations run locally in browser memory so unpublished structural coordinates never leave the user's computer.

---

## 2. Step-by-Step Methodology and Biophysical Logic

### Step 1: Universal Input Parsing and 3-Mode Confidence Auto-Detection
Biomolecular structures come from different sources that store different values in the coordinate B-factor column. In predicted structures (AF3/AF2/ESM), the B-factor column stores `pLDDT` (`0 to 100`, where higher values mean higher confidence). In experimental X-ray or Cryo-EM structures, the B-factor column stores Debye-Waller thermal displacement factors (`B = 8 * pi^2 * <u^2>`, where lower values mean higher certainty).

Our parser inspects file headers (`_exptl.method`, `_ma_qa_metric_local`, `EXPDTA`), B-factor distributions, and companion `.json` files to run automatically in one of three modes:
* **`MODE_A_PRED_FULL` (Predicted Structure + Companion JSON):** Reads 3D coordinates alongside AF3/AF2/ColabFold JSON files. Because AlphaFold3 indexes its `pae` matrix by *tokens* (splitting modified residues like phosphoserine `SEP` and ligands into one token per heavy atom), our parser maps `token_chain_ids` and `token_res_ids` back to exact residue keys `(chain_id, res_seq)` so the PAE matrix never shifts out of alignment.
* **`MODE_B_PRED_COORD` (Predicted Structure Only, No JSON):** Reads per-residue `pLDDT` directly from the coordinate B-factor column (automatically scaling fractional `[0, 1]` values by `100`).
* **`MODE_C_EXP_CRYSTAL` (Experimental Structure):** Applies inverse thermal normalization `c_i = exp(-max(0, B_i - B_min) / 60.0)` so experimental crystal structures are scored accurately without a JSON file.

### Step 2: Two-Stage 3D Voxel Spatial Hashing
Checking all-atom distances between two large chains takes `O(N_A * N_B)` operations. To make both the Python CLI and browser app run in milliseconds, we place residue centroids (`CA` for amino acids, `C4'` for DNA/RNA, geometric center for ligands) on Chain B into a 3D cubic lattice with voxel size `L_voxel = 14.0 A`:

    Voxel_Key(x, y, z) = ( floor(x / 14.0), floor(y / 14.0), floor(z / 14.0) )

For each residue on Chain A, we only evaluate heavy-atom distances against Chain B residues sitting in the `3 x 3 x 3 = 27` neighboring voxels whose centroids are within `16.0 A`.

### Step 3: Flat-Top Plateau Confidence Weighting (`w_ij`)
In plant and eukaryotic regulatory proteins, AlphaFold3 often wraps low-confidence disordered tails (`pLDDT < 40`, `PAE > 20 A`) around a partner protein. To silence these artifacts without penalizing normal thermal flexibility at well-folded interfaces, each residue `i` receives a plateau confidence weight `c_i` in `[0.02, 1.0]`:

    c_i = clip( (pLDDT_i - 30.0) / 55.0,  0.02,  1.00 )     [in Modes A and B]
    c_i = exp( -max(0.0, B_i - B_min) / 60.0 )              [in Mode C]

Any residue with `pLDDT >= 85.0` receives `c_i = 1.00` (crystal-grade). For an inter-chain residue pair `(i, j)` with mean cross-chain aligned error `PAE_ij = 0.5 * (PAE[i,j] + PAE[j,i])`, the joint confidence weight `w_ij` is:

    w_ij = c_i * c_j * exp( -max(0.0, PAE_ij - 3.5) / lambda_PAE )

where `lambda_PAE = 8.0 A`. A well-folded interface pair (`pLDDT = 92`, `PAE = 2.5 A`) gets `w_ij = 1.000`. A disordered loop hallucination (`pLDDT = 33`, `PAE = 27.5 A`) gets `w_ij = 0.00015`, removing `99.98%` of its false contact energy.

### Step 4: Seven-Channel Physicochemical Contact Classification
For every inter-chain residue pair with minimum heavy-atom distance `d_min <= 5.5 A`, we classify atomic interactions into seven physical channels:
1. **Salt Bridges (`d <= 4.0 A`):** Cationic side-chain nitrogens (`Arg NE,NH1,NH2`, `Lys NZ`, `His ND1,NE2`) to anionic oxygens (`Asp OD1,OD2`, `Glu OE1,OE2`, or nucleic acid phosphate `OP1,OP2,O1P,O2P`), capped at `2` per residue pair.
2. **Hydrogen Bonds (`d <= 3.5 A`):** Polar `N, O, S` donor-acceptor pairs not already counted as salt bridges, capped at `2` per residue pair.
3. **Hydrophobic Packing (`d <= 4.5 A`):** Side-chain apolar `C, S` atoms (excluding backbone `C, CA`) between hydrophobic residues (`Ala, Val, Leu, Ile, Met, Phe, Trp, Pro, Tyr` or ligands). To avoid coordinate jitter across different AF3 seeds, each apolar atom pair is weighted by a flat-top van der Waals plateau `phi_HP(d) = clip((4.75 - d) / 0.80, 0.35, 1.0)` that gives `1.0` for `d <= 3.95 A` and tapers to `0.35` at `4.50 A` (capped at `6.0` per residue pair).
4. **Aromatic Pi-Stacking (`d_ring <= 5.5 A`):** Distance between aromatic ring centroids (`Phe, Tyr, Trp, His`).
5. **Disulfide Bridges (`d <= 2.3 A`):** Covalent inter-chain `Cys SG` to `Cys SG` bonds.
6. **Steric Clashes (`d < 2.2 A`):** Unphysical non-bonded heavy-atom overlap (capped at `3` per residue pair).
7. **Intra-Helical Proline Kink Audit (`N_Pro_helix`):** Counts native Proline (`PRO`) residues located inside folded alpha-helical segments (`4.8 <= d_CA(i, i+4) <= 7.2 A` and `pLDDT >= 50`) within `+/- 2` residues of the interface.

### Step 5: Golden-Spiral Shrake-Rupley SASA and Confidence-Weighted Burial (`wBSA`)
We compute the Solvent-Accessible Surface Area (`SASA`) of the isolated chains (`SASA_A`, `SASA_B`) and bound complex (`SASA_AB`) using a 64-point Fibonacci Golden-Section Spiral around each atom with water probe radius `r_probe = 1.40 A`:

    BSA  = max( 0.0,  SASA_A + SASA_B - SASA_AB )
    wBSA = mean(w_ij) * BSA

Multiplying `BSA` by the mean interface confidence `mean(w_ij)` prevents disordered loop decoys from subtracting un-weighted geometric surface area energy. Surface-exposed non-interface residues (`SASA_complex >= 10.0 A^2`) define the Non-Interacting Surface fractions `%NIS_apolar` and `%NIS_charged`.

### Step 6: Hybrid Contact-Gated Free Energy (`dG_uni`) and Salt/pH Modulation
Standard linear contact models like raw PRODIGY omit `apolar-apolar` (`IC_aa`) contacts and use a fixed `-15.9433 kcal/mol` intercept, which under-scores hydrophobic leucine zippers (like *Arabidopsis* VIP1) while over-scoring sparse 1-contact decoys.

We combine the PRODIGY linear equation (`dG_PRODIGY`) with an explicit biophysical solvation and contact potential (`dG_phys`) that includes translational/rotational cratic entropy (`+2.10 kcal/mol`), explicit hydrophobic burial, electrostatic networks, polar desolvation, steric clashes, intra-helical Proline kink penalties (`+2.20 kcal/mol` per Proline), and environmental ionic strength (`I` in Molar) and `pH`:

    # Environmental modulators (at baseline I = 0.15 M NaCl and pH = 7.2, all equal 1.000):
    shield_f    = 1.0 - 0.65 * f_apolar
    debye_kappa = exp( -0.85 * shield_f * (sqrt(I)/3.04 - sqrt(0.15)/3.04) * 3.5 )
    ph_kappa    = [ (1 / (1 + 10^(4.2 - pH))) * (1 / (1 + 10^(pH - 10.4))) ] / 0.99837
    salting_out = 1.0 + (0.20 + 0.25 * f_apolar) * (I - 0.15)

    # Biophysical Solvation & Contact Potential:
    dG_phys = +2.10
              - salting_out * (0.110 * w_HP + 0.215 * IC_aa + 0.0026 * wBSA)
              - (0.765 * w_SB * debye_kappa * ph_kappa + 0.460 * w_HB + 0.615 * w_Pi)
              + 0.160 * IC_pp + 1.50 * w_clash + 2.20 * N_Pro_helix

    # Contact-density sigmoid gate (silences PRODIGY's -15.94 intercept on sparse decoys):
    gate   = 0.0                                         if IC_total < 3.0
    gate   = 1.0 / (1.0 + exp(-0.40 * (IC_total - 11.0))) otherwise
    w_phys = clip(0.55 + 0.35 * f_apolar + (0.15 if N_Pro_helix > 0 else 0.0), 0.55, 0.95)

    dG_uni = clip( gate * ((1 - w_phys)*dG_PRODIGY + w_phys*dG_phys) + (1 - gate)*dG_phys, -22.0, +5.0 )
    Kd     = exp( dG_uni / (R * T) )

### Step 7: Universal Interaction Score (`UIS`) and Mechanistic Patterns (`Patterns A-D`)
To rank candidates across batch screens, `UIS` in `[0, 1]` combines global chain topology (`ipTM`), local IDR-immune PAE coherence (`LIS`, evaluated strictly on cross-chain pairs with `PAE_ij <= 12.0 A`), size-normalized coordinate quality (`pDockQ_PAE`), and scaled affinity `sigma_aff = clip((-dG_uni - 5.0) / 8.5, 0.0, 1.0)`:

    UIS = 0.25 * ipTM + 0.30 * LIS + 0.25 * pDockQ_PAE + 0.20 * sigma_aff   [Mode A]
    UIS = 0.55 * pDockQ_coord + 0.45 * sigma_aff                            [Modes B & C]

Using the frozen decision threshold `UIS = 0.5000`, every complex is assigned to one of four mechanistic patterns:
* **Pattern A (Concordant Rigid / Coiled-Coil or Globular Interactor):** `UIS >= 0.5000`, `ipTM >= 0.60`, and `dG_uni <= -6.5 kcal/mol`. Global topology, local PAE coherence (`LIS`), and solvation thermodynamics all support high-affinity binding.
* **Pattern B (IDR-Diluted Local Interactor):** `UIS >= 0.5000` (or `LIS >= 0.50` and `dG_uni <= -7.0 kcal/mol`) but global `ipTM < 0.60`. Identifies full-length proteins where long non-interacting disordered tails lower global `ipTM` despite a well-folded local binding interface.
* **Pattern C (Thermodynamically Penalized / Kink-Disrupted Interface):** Chains are placed near each other in 3D space, but binding is opposed by internal helix-breaking Proline kinks (`N_Pro_helix > 0`), steric clashes, or unfavorable desolvation.
* **Pattern D (Below-Threshold Decoy or Disordered Hallucination):** `UIS < 0.5000` with low interface confidence `w_ij`. Spurious contacts from disordered loops are damped out, leaving `dG_uni` at the `+2.10 kcal/mol` cratic barrier.

---

## 3. Parameter Reference Table and Sensitivity Analysis

| Parameter | Value Used | Biophysical Reason | Effect of Changing the Parameter |
| :--- | :--- | :--- | :--- |
| `UIS_threshold` | `0.5000` | Frozen decision boundary separating true complexes from decoys (`High-Confidence >= 0.6500`). | Raising to `0.75` misses transient complexes (`Pattern B`); lowering to `0.25` admits weak decoys. |
| `voxel_size` | `14.0 A` | Maximum `CA-CA` distance capable of forming side-chain heavy-atom contacts (`<= 5.5 A`). | Lowering below `10.0 A` misses extended Arg-Glu salt bridges; raising above `25.0 A` slows runtime. |
| `WATER_PROBE` | `1.40 A` | Effective radius of a water molecule (`H2O`) in 64-point Golden-Spiral SASA integration. | Setting to `0.0 A` measures van der Waals surface rather than solvent-accessible surface. |
| `lambda_pae` | `8.0 A` | Exponentially damps uncertain inter-chain contacts above thermal noise (`PAE > 3.5 A`). | Lowering to `2.0 A` over-penalizes flexible loops; raising to `25.0 A` allows IDR decoys to leak through. |
| `LIS_cutoff` | `12.0 A` | Restricts local cross-chain PAE scoring to the binding domain, ignoring disordered tails. | Increasing to `30.0 A` dilutes local scores on full-length plant transcription factors. |
| `d_SB, d_HB, d_HP` | `4.0, 3.5, 4.5 A` | Standard heavy-atom distance limits for ion pairs, hydrogen bonds, and apolar methyl packing. | Using a flat `5-8 A` cutoff conflates non-bonded proximity with chemical bonds. |
| `dg_cratic` | `+2.10 kcal/mol` | Standard-state (`1 M`) translational and rotational entropy penalty of bimolecular association. | Setting to `0.0` causes non-contacting chains to register artificial negative free energy. |
| `dg_proline_kink` | `+2.20 kcal/mol` | Thermodynamic cost of losing the `i to i-4` backbone H-bond and kinking an alpha-helix. | Without this penalty, forced helical models of Group E bZIP homodimers (`AtbZIP34`) falsely score as binders. |

---

## 4. Six-Complex Validation Benchmark and Biological Interpretation

Both `unippi_cli.py` and `index.html` include a built-in 6-complex benchmark panel spanning *Arabidopsis thaliana* bZIP transcription factors (`VIP1 / AtbZIP51`, Group I; `AtbZIP61` and `AtbZIP34`, Group E), a cross-kingdom electrostatic globular standard (`Barnase-Barstar`), an IDR-diluted complex (`VIP1-14-3-3 GRF6`), and a disordered tail decoy:

| Rank | Complex ID | UIS | Margin | dG_uni (kcal/mol) | Kd (nM) | ipTM | LIS | wBSA (A^2) | SB | HB | HP | Pro Kinks | Pattern | Class |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **#1** | `1_VIP1_Homodimer_Positive` | **0.777** | `+0.277` | **-10.35** | **25.89** | 0.90 | 0.808 | 479.4 | 4 | 0 | 48 | 0 | **Pattern A** | `HIGH_CONFIDENCE` |
| **#2** | `3_Barnase_Barstar_Electrostatic_Pos` | **0.771** | `+0.271` | **-10.09** | **39.87** | 0.89 | 0.817 | 474.7 | 16 | 4 | 0 | 0 | **Pattern A** | `HIGH_CONFIDENCE` |
| **#3** | `2_VIP1_bZIP61_GroupE_Positive` | **0.732** | `+0.232` | **-10.10** | **39.68** | 0.83 | 0.742 | 406.8 | 4 | 0 | 48 | 0 | **Pattern A** | `HIGH_CONFIDENCE` |
| **#4** | `4_VIP1_GRF6_14-3-3_IDR_Diluted_Pos` | **0.609** | `+0.109` | **-7.61** | **2644.23** | 0.52 | 0.783 | 473.2 | 4 | 2 | 36 | 0 | **Pattern B** | `MODERATE` |
| **#5** | `5_bZIP34_Proline_Homodimer_Negative` | **0.120** | `-0.380` | **+5.00** | `>10mM (NB)` | 0.38 | 0.000 | 16.8 | 0 | 0 | 5 | 6 | **Pattern C** | `LOW_OR_DECOY` |
| **#6** | `6_VIP1_IDR_Hallucination_Decoy` | **0.061** | `-0.439` | **+2.10** | `>10mM (NB)` | 0.24 | 0.000 | 0.1 | 0 | 0 | 0 | 0 | **Pattern D** | `LOW_OR_DECOY` |

### Biological Takeaways from the Benchmark
1. **Selective Hydrophobic Zipper Dimerization (`#1` and `#3`) vs. Polar Globular Binding (`#2`):**
   *Arabidopsis* VIP1 (`AtbZIP51`) homodimerizes (`dG_uni = -10.35 kcal/mol`, `Kd = 25.89 nM`) and heterodimerizes with Group E `AtbZIP61` (`dG_uni = -10.10 kcal/mol`, `Kd = 39.68 nM`) through `48` hydrophobic side-chain contacts (`HP = 48`) along heptad `a` and `d` positions, stabilized by `4` flanking `g <-> e'` Glu-Lys salt bridges (`SB = 4`). At the same time, the engine accurately scores the purely polar/electrostatic `Barnase-Barstar` complex (`16 SB, 4 HB, 0 HP`, `dG_uni = -10.09 kcal/mol`), demonstrating universality across hydrophobic and electrostatic binding modes.
2. **Salt-Stress Resilience of VIP1 vs. Electrostatic Screening:**
   Increasing ionic strength from `150 mM` to `500-1000 mM NaCl` screens exposed salt bridges in `Barnase-Barstar` (weakening its binding affinity), whereas `VIP1` dimers maintain nanomolar stability through hydrophobic salting-out—consistent with VIP1's physiological function during plant osmotic and salt stress.
3. **Why Group E bZIPs (`AtbZIP34`) Fail to Homodimerize (`#5`, `Pattern C`):**
   *Arabidopsis* Group E bZIPs contain conserved internal Proline residues inside their zipper domain. UniPPI-OWL detects all `6` intra-helical Proline kinks in the `AtbZIP34` homodimer, penalizing `dG_uni` to `+5.00 kcal/mol` (`>10mM Non-Binding`). In-silico mutagenesis confirms that mutating `A:LEU5 -> PRO` in wild-type VIP1 destabilizes binding by `ddG = +4.20 kcal/mol` (`Kd` shifts from `25.89 nM` to `31,046.93 nM`), whereas mutating `PRO5 -> LEU` across both chains in `AtbZIP34` cooperatively rescues helical stability (`ddG_1,2 = -5.10 kcal/mol`).
4. **Rescuing IDR-Diluted Complexes (`#4`, `Pattern B`) and Silencing Loop Decoys (`#6`, `Pattern D`):**
   In `4_VIP1_GRF6_14-3-3_IDR_Diluted_Pos`, flanking disordered tails dilute global `ipTM` to `0.52`, which standard `ipTM >= 0.60` cutoffs reject; our local `LIS` (`0.783`) and `dG_uni` (`-7.61 kcal/mol`) recover the interaction as `Pattern B`. Conversely, in `6_VIP1_IDR_Hallucination_Decoy`, low confidence (`pLDDT = 33`, `PAE = 27.5 A`) reduces effective buried surface area to `wBSA = 0.1 A^2`, locking `dG_uni` at `+2.10 kcal/mol` (`Pattern D`).

---

## 5. Quick-Start Usage

### Option A: Zero-Install Browser Application (`index.html`)
1. Open `index.html` directly in any modern web browser (or visit the live GitHub Pages URL).
2. Click **Reset 6-Complex Benchmark Suite** to inspect the pre-loaded validation panel, or drag-and-drop any `.cif`, `.mmcif`, `.pdb`, and AF3 `.json` files into the upload box.
3. Click any row in the leaderboard to view its 3D interface in `3Dmol.js`, test single and double point mutations across all 20 amino acids, adjust temperature, `[NaCl]`, and `pH`, and download vector `.svg` plots, `.csv`/`.tsv` tables, sanitized `.pdb` files, `.sdf` ligands, and FoldX `individual_list.txt` files.

### Option B: Command-Line Interface (`unippi_cli.py`)
Requires only Python 3.8+ (zero external packages):

    # 1. Generate and run the 6-complex benchmark suite with custom point mutations:
    python3 unippi_cli.py --generate-demo ./demo_structures -o ./demo_results --mutate A:LEU5PRO,A:PRO5LEU,A:ASP2ALA,A:TYR5ALA

    # 2. Screen a local directory of AlphaFold3 / PDB outputs under 400 mM NaCl salt stress:
    python3 unippi_cli.py -i ./my_af3_outputs -o ./screening_results --temp 298.15 --ionic-strength 400.0 --ph 7.2