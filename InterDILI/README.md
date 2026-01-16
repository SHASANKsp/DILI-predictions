# **Mechanistic Interpretation of Descriptor and Substructure Features Predictive of Drug-Induced Liver Injury in SMILES-Driven Machine Learning Models : InterDILI**

## **Introduction**

The liver’s central role in drug metabolism makes it especially susceptible to xenobiotic stress. DILI mechanisms span reactive metabolite formation, transporter inhibition leading to cholestasis, mitochondrial dysfunction, oxidative stress, and immune-mediated responses. Machine learning approaches can identify chemical features associated with DILI risk, but without mechanistic interpretation these remain a black box.


## **Method**
The InterDILI data was downloaded and with molecular inputs as SMILES strings were processed with RDKit to produce:

* Extended connectivity fingerprints (Morgan fingerprint) and 3D and 2D descriptors 
* Feature attribution for model interpretability to rank the most predictive properties/bits.

Few descriptors along with six fingerprint bits were isolated based on their SHAP contributions to DILI classification. Each bit was visualized with representative molecules from the dataset.

## **Results**

### **RDKit Descriptors Identified as Top Predictive Features**

The InterDILI dataset model highlights the following descriptors as dominant contributors:

* **BCUT2D_LOGPLOW / BCUT2D_MRLOW / BCUT2D_MWLOW / BCUT2D_CHGLO:**
  Lowest eigenvalues of adjacency matrices weighted by lipophilicity, polarizability, atomic mass, and charge, respectively. These descriptors encode **electronic and physicochemical asymmetry** across the molecular graph.

* **SlogP_VSA8 / SlogP_VSA10**
  Surface area contributions from atoms within specific lipophilicity ranges, capturing **spatial clustering of lipophilic regions**.

* **PEOE_VSA3 / PEOE_VSA5 / PEOE_VSA13**
  Surface area weighted by partial charge, reflecting **charge segregation and membrane-interaction potential**.

* **VSA_EState10**
  Electrotopological state–weighted surface area, integrating electronic activation with steric accessibility.

* **NumAromaticRings / NumAromaticHeterocycles**
  Counts of aromatic and heteroaromatic rings, indicating **aromatic persistence and metabolic engagement**.

* **fr_Ar_N**
  Aromatic nitrogen-containing motifs associated with heterocycle metabolism.

* **fr_aniline / fr_NHO / fr_hydrazine**
  Classical structural alerts linked to **CYP-mediated bioactivation and reactive intermediate formation**.

Collectively, these descriptors emphasize *where* lipophilicity, charge, and reactivity are expressed within a molecule, rather than their absolute magnitudes.



### **Fingerprint Bits and Mechanistic Interpretation**

Below are the fingerprint bits identified as highly influential in this model, interpreted in the context of the above descriptors.

---

### **Bit 1 — Compact Heteroaromatic Carbonyl Cores**

**Representative Substructure:**
Small to medium heteroaromatic systems containing lactams, imides, or cyclic amides.

**Mechanistic Interpretation:**
These rigid, polar heterocycles funnel metabolism through a limited set of oxidative pathways, resulting in **metabolic focusing**. High flux through constrained CYP routes increases the likelihood of oxidative stress, mitochondrial perturbation, and hepatocellular injury despite apparently “clean” medicinal chemistry profiles.

**Associated Descriptors:**
BCUT2D_CHGLO, VSA_EState10, NumAromaticHeterocycles.

---

### **Bit 2 — Extended Conjugated and Macrocyclic Frameworks**

**Representative Substructure:**
Long conjugated systems, macrocycles, or heterocycles connected by flexible linkers.

**Mechanistic Interpretation:**
These structures exhibit **high hepatic exposure and prolonged intracellular residence** due to partial lipophilicity combined with continued metabolic engagement. Clearance is slow and incomplete, leading to chronic intracellular stress rather than acute toxicity.

**Associated Descriptors:**
SlogP_VSA8, SlogP_VSA10, PEOE_VSA*, BCUT2D_LOGPLOW.

---

### **Bit 3 — Fused Aromatic Heterocycles with Polar Exit Vectors**

**Representative Substructure:**
Fused bicyclic or tricyclic aromatic systems with one or two heteroatoms and a single polar substituent.

**Mechanistic Interpretation:**
Aromatic persistence promotes retention, while the polar substituent maintains metabolic turnover. This asymmetry results in **ongoing formation of low-level reactive metabolites** and chronic hepatocyte stress.

**Associated Descriptors:**
NumAromaticRings, NumAromaticHeterocycles, BCUT2D_MRLOW.

---

### **Bit 4 — Bulky Aromatic–Polar Hybrids**

**Representative Substructure:**
Multi-aromatic frameworks bearing phenols, amides, or carboxylic acids.

**Mechanistic Interpretation:**
These structures combine accumulation-prone aromatic bulk with metabolically active polar groups, producing a **clearance imbalance** where exposure persists but detoxification never fully resolves. This pattern is characteristic of delayed and idiosyncratic DILI phenotypes.

**Associated Descriptors:**
SlogP_VSA*, PEOE_VSA*, BCUT2D_LOGPLOW.



## **Conclusion**

The InterDILI model defines DILI risk through **structural organization rather than simple physicochemical thresholds**. Four dominant hepatic liability modes emerge:

1. **Electronic and lipophilic asymmetry driving metabolic hot spots**
   (BCUT2D_*LOW family)

2. **Surface-area partitioning promoting accumulation and transporter interaction**
   (SlogP_VSA*, PEOE_VSA*)

3. **Aromatic and heterocyclic persistence with continued metabolic engagement**
   (NumAromaticRings, NumAromaticHeterocycles, fr_Ar_N)

4. **Explicit bioactivation via classical structural alerts**
   (fr_aniline, fr_NHO, fr_hydrazine)

Compared to DILIrank model, this model shows **greater reliance on second-order descriptors**, indicating that increased data volume enables finer discrimination of mechanistic DILI drivers.


This InterDILI DILI model captures hepatotoxic risk as an emergent property of **electronic asymmetry, spatial distribution of lipophilicity and charge, aromatic persistence, and known bioactivation motifs**. Rather than relying on coarse descriptors such as molecular weight or total logP, the model encodes how chemical properties are organized within the molecular structure. This mechanistic grounding enhances confidence in both predictive performance and downstream use for safety-aware chemical exploration.


## **References**
1. David, Stefan, and James P. Hamilton. "Drug-induced liver injury." US gastroenterology & hepatology review 6 (2010): 73.
2. Andrade, Raul J., et al. "Drug-induced liver injury." Nature reviews Disease primers 5.1 (2019): 58.
3. Allison, Rebecca, et al. "Drug induced liver injury–a 2023 update." Journal of Toxicology and Environmental Health, Part B 26.8 (2023): 442-467.
4. Yuan, Liyun, and Neil Kaplowitz. "Mechanisms of drug induced liver injury." Clinics in liver disease 17.4 (2013): 507.
5. Saleh, Ahmed K., et al. "Exploring drug-induced liver injury: comprehensive insights into mechanisms and management of hepatotoxic agents." Future Journal of Pharmaceutical Sciences 11.1 (2025): 38.
6. Liu, Anika, et al. "Prediction and mechanistic analysis of drug-induced liver injury (DILI) based on chemical structure." Biology direct 16.1 (2021): 6.
7. Pizzo, Fabiola, et al. "A new structure-activity relationship (SAR) model for predicting drug-induced liver injury, based on statistical and expert-based structural alerts." Frontiers in Pharmacology 7 (2016): 442.
8. Lee, Taeyeub, and Joram M. Posma. "Improving drug-induced liver injury prediction using graph neural networks with augmented graph features from molecular optimisation." Journal of Cheminformatics 17.1 (2025): 124.
9. Xu, Youjun, et al. "Deep learning for drug-induced liver injury." Journal of chemical information and modeling 55.10 (2015): 2085-2093.
10. Kang, Myung-Gyun, and Nam Sook Kang. "Predictive model for drug-induced liver injury using deep neural networks based on substructure space." Molecules 26.24 (2021): 7548.
11. Seal, Srijit, et al. "Improved Detection of Drug-Induced Liver Injury by Integrating Predicted In Vivo and In Vitro Data." (2024).
