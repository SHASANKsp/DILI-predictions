**Mechanistic Interpretation of Descriptor and Substructure Features Predictive of Drug-Induced Liver Injury in SMILES-Driven Machine Learning Models**

---

## **Abstract**

Drug-induced liver injury (DILI) is a major cause of clinical attrition despite extensive in vitro toxicity screening. We applied interpretable machine learning to molecular descriptors and fingerprint bits derived from SMILES strings to identify structural and physicochemical features predictive of DILI risk. Six fingerprint bits, corresponding to distinct molecular substructures, emerged as highly informative. These bits, when interpreted in the context of top RDKit descriptors including aromaticity, heterocycle count, lipophilicity proxies, and topological features, reveal mechanistically plausible pathways to hepatotoxicity such as reactive metabolite formation, transporter inhibition, metabolic stress, clearance overload, and intracellular accumulation.

---

## **Introduction**

The liver’s central role in drug metabolism makes it especially susceptible to xenobiotic stress. DILI mechanisms span reactive metabolite formation, transporter inhibition leading to cholestasis, mitochondrial dysfunction, oxidative stress, and immune-mediated responses. Machine learning approaches can identify chemical features associated with DILI risk, but without mechanistic interpretation these remain opaque.

---

## **Methods Summary**

Molecular inputs were SMILES strings processed with RDKit to produce:

* Extended connectivity fingerprints (e.g., Morgan),
* 2D descriptors including aromatic counts, heterocycles, BCUT eigenvalues, SlogP-weighted surface areas, and topological indices (e.g., BalabanJ),
* Feature attribution for model interpretability to rank the most predictive bits.

Six fingerprint bits were isolated based on their SHAP contributions to positive DILI classification. Each bit was visualized with representative molecules from the dataset.

---

## **Results**

### **RDKit Descriptors Identified as Top Predictive Features**

The model highlighted the following descriptors:

* **fr_Ar_N**: Aromatic nitrogen count – linked to CYP interactions and bioactivation.
* **fr_NHO**: N-hydroxyl motifs – known precursors to reactive intermediates.
* **fr_aniline**: Aniline functional groups – classic bioactivation liability.
* **SlogP_VSA6 / SMR_VSA7**: Surface area weighted by lipophilicity/polarizability – proxies for hepatic accumulation.
* **NumAromaticHeterocycles / NumAromaticCarbocycles**: Counts of aromatic systems – linked with lipophilicity and metabolic persistence.
* **FpDensityMorgan2/3**: Fingerprint density – correlated with structural complexity and diverse metabolic pathways.
* **BCUT2D_CHGLO / BCUT2D_MRLOW**: Charge and polarizability eigenvalues – reflect electronic distribution relevant to bioactivation and protein interactions.

These descriptors align with known DILI risk factors such as high lipophilicity, aromatic complexity, and heteroatom-rich motifs prone to oxidative metabolism.

---

## **Fingerprint Bits and Mechanistic Interpretation**

Below are the six fingerprint bits, each interpreted mechanistically and in the context of relevant descriptors:

---

### **Bit 1 — Reactive Bioactivation Motif**

**Representative Substructure:** Phenolic / aniline / aromatic heterocycle core.

**Mechanistic Interpretation:**
Rich in aromatic nitrogens and hydroxyls, this motif undergoes CYP-mediated oxidation to generate reactive electrophiles, leading to covalent protein adducts and oxidative stress — a well-established DILI mechanism.
**Associated Descriptors:** fr_Ar_N, fr_NHO, fr_aniline, NumAromaticHeterocycles.

**Literature Support:** Bioactivation of aniline and phenolic motifs → quinone and nitroso intermediates implicated in DILI.

---

### **Bit 2 — Hepatic Accumulation / Transporter Interaction**

**Representative Substructure:** Bulky polyaromatic systems often with basic side chains.

**Mechanistic Interpretation:**
Cationic, lipophilic structures tend to accumulate in hepatocytes and inhibit canalicular transporters like BSEP, leading to cholestasis and bile acid accumulation.
**Associated Descriptors:** SlogP_VSA6, SMR_VSA7, BalabanJ, NumAromaticCarbocycles.

**Literature Support:** Structural features correlating with BSEP inhibition and cholestatic DILI.

---

### **Bit 3 — Metabolic Bottleneck Heterocycle**

**Representative Substructure:** Planar fused heterocycles with constrained pathways.

**Mechanistic Interpretation:**
Fused heterocyclic systems funnel metabolism into a limited set of oxidative pathways, increasing flux through stress-prone intermediates and generating oxidative stress.
**Associated Descriptors:** NumAromaticHeterocycles, BCUT2D_CHGLO.

**Literature Support:** Heterocycles with dense heteroatom patterns are frequently associated with hepatotoxic metabolites.

---

### **Bit 4 — Polar Clearance / Reactive Conjugates**

**Representative Substructure:** Aliphatic polar moieties including carboxylic acids and thioamides.

**Mechanistic Interpretation:**
These structures are avidly taken up by the liver and undergo extensive phase II metabolism (e.g., glucuronidation). Acyl-glucuronides and sulfur-derived metabolites can form protein adducts and stress detoxification pathways.
**Associated Descriptors:** fr_NHO, fingerprint density.

**Literature Support:** Carboxylic acid derivatives and thioamides are established as clearance-associated DILI risks.

---

### **Bit 5 — Mixed Polar/Lipophilic Accumulation**

**Representative Substructure:** Hybrid polar and aromatic frameworks.

**Mechanistic Interpretation:**
Strong lipophilicity encourages hepatic partitioning, whereas polar functional groups induce continued metabolism. This combination yields slow clearance with persistent intracellular exposure, increasing oxidative burden.
**Associated Descriptors:** SlogP_VSA6, SMR_VSA7, FpDensityMorgan3.

**Literature Support:** Mixed polarity and high lipophilicity correlate with chronic hepatic stress.

---

### **Bit 6 — Aliphatic Metabolic Stress Pattern**

**Representative Substructure:** Branched aliphatic chains with non-aromatic labile groups.

**Mechanistic Interpretation:**
While not strongly aromatic, these aliphatic motifs are substrates for high-flux metabolism, producing a broad array of metabolites and generating cumulative stress on hepatic detoxification processes.
**Associated Descriptors:** FpDensityMorgan2, BCUT2D_MRLOW.

**Literature Support:** High metabolic turnover and reactive intermediate formation from aliphatic substrates contribute to hepatocellular stress.

---

## **Discussion**

Together, the six bits reveal **distinct hepatic liability classes**:

1. **Reactive metabolite formation** — Bit 1 & Bit 3
2. **Transporter inhibition and accumulation** — Bit 2 & Bit 5
3. **Metabolic overload and conjugate reactivity** — Bit 4 & Bit 6

Their alignment with descriptive features such as aromaticity, heteroatom count, lipophilicity surface area, and topological complexity strengthens biological plausibility. These relationships mirror known mechanisms of DILI as reported in toxicology literature.

---

## **Conclusion**

The combined interpretation of molecular descriptors and fingerprint bits provides a **mechanistically grounded understanding** of structural features that contribute to prediction of DILI risk from SMILES alone. This integrative approach supports both improved interpretability and practical utility in early safety screening and medicinal chemistry optimization.

---

## **References**



