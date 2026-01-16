# **Mechanistic Interpretation of Descriptor and Substructure Features Predictive of Drug-Induced Liver Injury in SMILES-Driven Machine Learning Models : DILIrank**

## **Introduction**

The liver’s central role in drug metabolism makes it especially susceptible to xenobiotic stress. DILI mechanisms span reactive metabolite formation, transporter inhibition leading to cholestasis, mitochondrial dysfunction, oxidative stress, and immune-mediated responses. Machine learning approaches can identify chemical features associated with DILI risk, but without mechanistic interpretation these remain a black box.


## **Method**
The DILIrank2 data was downloaded from FDA website which had the drug name, severaty and the DILI classes. Their corresponding SMILES were downloaded using Chembl. Hence the molecular inputs were SMILES strings processed with RDKit to produce:

* Extended connectivity fingerprints (Morgan fingerprint) and 3D and 2D descriptors 
* Feature attribution for model interpretability to rank the most predictive properties/bits.

Few descriptors along with six fingerprint bits were isolated based on their SHAP contributions to DILI classification. Each bit was visualized with representative molecules from the dataset.  



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

### **Fingerprint Bits and Mechanistic Interpretation**

Below are the six fingerprint bits, each interpreted mechanistically and in the context of relevant descriptors:   
(Might be a bit bised as I tried to make scence of the sub-structures that shows up as high importance for the model)



### **Bit 1 — Reactive Bioactivation Motif** - Phenolic / aniline / aromatic heterocycle core.

**Interpretation:** Rich in aromatic nitrogens and hydroxyls, this motif undergoes CYP-mediated oxidation to generate reactive electrophiles, leading to covalent protein adducts and oxidative stress.  
**Associated Descriptors:** fr_Ar_N, fr_NHO, fr_aniline, NumAromaticHeterocycles.


### **Bit 2 — Hepatic Accumulation / Transporter Interaction** - Bulky polyaromatic systems often with basic side chains.

**Interpretation:**
Cationic, lipophilic structures tend to accumulate in hepatocytes and inhibit canalicular transporters like BSEP, leading to cholestasis and bile acid accumulation.  
**Associated Descriptors:** SlogP_VSA6, SMR_VSA7, BalabanJ, NumAromaticCarbocycles.


### **Bit 3 — Metabolic Bottleneck Heterocycle** - Planar fused heterocycles with constrained pathways.

**Interpretation:**
Fused heterocyclic systems funnel metabolism into a limited set of oxidative pathways, increasing flux through stress-prone intermediates and generating oxidative stress.  
**Associated Descriptors:** NumAromaticHeterocycles, BCUT2D_CHGLO.


### **Bit 4 — Mixed Polar/Lipophilic Accumulation** - Hybrid polar and aromatic frameworks.

**Interpretation:**
Strong lipophilicity encourages hepatic partitioning, whereas polar functional groups induce continued metabolism. This combination yields slow clearance with persistent intracellular exposure, increasing oxidative burden.  
**Associated Descriptors:** SlogP_VSA6, SMR_VSA7, FpDensityMorgan3.



## **Conclusion**

Together, all the bits reveal **distinct hepatic liability classes**:

1. **Reactive metabolite formation**
2. **Transporter inhibition and accumulation** 
3. **Metabolic overload and conjugate reactivity**

Their alignment with descriptive features such as aromaticity, heteroatom count, lipophilicity surface area, and topological complexity strengthens biological plausibility. These relationships mirror known mechanisms of DILI as reported in toxicology literature.


The combined interpretation of molecular descriptors and fingerprint bits provides a **mechanistically grounded understanding** of structural features that contribute to prediction of DILI risk from SMILES alone. This integrative approach supports both improved interpretability and practical utility in early safety screening and medicinal chemistry optimization. Although theres very limited performance that can be extracted from SMILES alone and including data from experiments will add biological context to the model


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



