# **Mechanistic Interpretation of Descriptor and Substructure Features Predictive of Drug-Induced Liver Injury in SMILES-Driven Machine Learning Models**

## **Abstract**
Drug-induced liver injury (DILI) is a major cause of clinical attrition despite extensive in vitro toxicity screening. I applied interpretable machine learning to molecular descriptors and fingerprint bits derived from SMILES strings to identify structural and physicochemical features predictive of DILI risk. These bits, when interpreted in the context of top RDKit descriptors including aromaticity, heterocycle count, lipophilicity proxies, and topological features, reveal mechanistically plausible pathways to hepatotoxicity such as reactive metabolite formation, transporter inhibition, metabolic stress, clearance overload, and intracellular accumulation. I build two model using the same approach with two dataset **DILIrank** and **InterDILI**.  
InterDILI: 1850 drugs  
DILIrank: 500 drugs

Both the SMILES-based machine learning models were developed to predict drug-induced liver injury (DILI) using RDKit descriptors and molecular fingerprint bits. Model-1 was trained on a smaller dataset (DILIrank), while Model-2 leveraged a substantially larger chemical space (InterDILI). Although both models identify mechanistically plausible DILI drivers, they differ in **resolution, emphasis, and abstraction level**. Below I have clarified how dataset scale influences the type of hepatotoxicity signals learned and how each model should be interpreted and used.

---

## **Overview of the Two Models**

| Aspect               | Model-1                         | Model-2                              |
| -------------------- | ------------------------------- | ------------------------------------ |
| Dataset size         | 500              | 1850                     |
| Dominant signal type | Explicit toxicophores           | Physicochemical asymmetry & exposure |
| Interpretability     | High, local                     | High, global                         |
| Best use             | Mechanistic explanation, alerts | Generalization, safety guidance      |


Explore each folder for in detail information about each model    

---

## **Descriptor-Level Comparison**

### **Model-1 Descriptor Emphasis**

Model-1 is dominated by descriptors and fragments that correspond to **classical, chemically explicit DILI mechanisms**, including:

* Direct bioactivation alerts (e.g., anilines, N-hydroxyls)
* Aromatic and heterocycle counts
* Lipophilicity surface area proxies
* Fingerprint density reflecting structural complexity

These features emphasize **what functional groups are present** and how they traditionally fail in the liver.

---

### **Model-2 Descriptor Emphasis**

Model-2 shifts decisively toward **second-order physicochemical organization**, highlighted by:

* The **BCUT2D family**, capturing electronic, lipophilic, mass, and charge asymmetry
* Multiple **VSA-weighted descriptors** partitioned by lipophilicity and charge
* Continued but more contextual use of aromaticity and heterocycles
* Retention of only the most severe classical alerts (aniline, NHO, hydrazine)

This indicates that with more data, the model becomes less reliant on “named toxicophores” and more sensitive to **how properties are distributed across the molecule**.

---

## **Fingerprint Bit Comparison**

### **Model-1 Bits: Mechanism-Centric Alerts**

Model-1 fingerprint bits cluster into **clear mechanistic classes**:

* Reactive metabolite precursors
* Transporter inhibition and accumulation
* Metabolic bottlenecks
* Clearance overload

Each bit can be readily explained by a dominant failure pathway, making Model-1 especially suitable for **post-hoc toxicological interpretation**.

---

### **Model-2 Bits: Structural Organization Classes**

Model-2 bits are broader and more architectural in nature:

* Compact heteroaromatic carbonyl cores (metabolic focusing)
* Extended conjugated or macrocyclic frameworks (residence time)
* Fused aromatic heterocycles with polar exit vectors (persistent stress)
* Bulky aromatic–polar hybrids (clearance imbalance)

Rather than pointing to a single toxic mechanism, these bits define **liability envelopes** — structural patterns that increase the probability of DILI through combined exposure, metabolism, and persistence.

---

## **Mechanistic Resolution: How the Models Differ**

| Dimension         | Model-1                 | Model-2                           |
| ----------------- | ----------------------- | --------------------------------- |
| Bioactivation     | Explicit and dominant   | Present but contextual            |
| Accumulation      | Clear via lipophilicity | Emergent via surface partitioning |
| Metabolic stress  | Fragment-driven         | Asymmetry-driven                  |
| Clearance failure | Identified via motifs   | Identified via organization       |
| DILI type         | Mechanism-specific      | Exposure-integrated               |

In short:

* **Model-1 asks:** *What chemical liabilities are present?*
* **Model-2 asks:** *How does this molecule behave in the liver as a system?*

---

## **Effect of Dataset Scale**

The transition from Model-1 to Model-2 reflects a well-known phenomenon in cheminformatics:

* Smaller datasets favor **explicit, high-signal features**
* Larger datasets allow learning of **distributed, second-order patterns**

Model-2’s reliance on BCUT and VSA descriptors strongly suggests that increased data volume enabled the model to resolve **subtle electronic and spatial effects** that are invisible in smaller datasets.

---

## **Complementarity of the Two Models**

Importantly, the models are **not redundant**.

* Model-1 provides **clear mechanistic flags** and aligns closely with toxicology intuition.
* Model-2 provides **generalizable safety guidance** across broader chemical space.

Used together:

* Model-1 explains *why* a molecule is risky.
* Model-2 estimates *how likely* that risk is to manifest given overall structure.

---

## **Practical Interpretation**

From a decision-making standpoint:

* Agreement between Model-1 and Model-2 indicates **high confidence DILI risk**
* Disagreement often highlights:

  * Borderline chemical space
  * Context-dependent risk
  * Optimization opportunities

This dual-model view mirrors real drug development, where both **structural alerts** and **exposure-driven liabilities** must be considered.

---

## **Conclusion**

The two DILI models represent successive layers of chemical understanding derived from SMILES alone. Model-1 captures well-established toxicological motifs and provides mechanistic clarity. Model-2, enabled by a larger dataset, abstracts these mechanisms into higher-order organizational principles governing hepatic exposure, metabolism, and stress.

## Future Work
Using the model for reinforcement learning in REINVENT4 for DILI aware molecule generation. 

