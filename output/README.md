## Cluster Tree (BuildTree Module)

### 1. Node Meaning
- Each colored circle represents a clone or subclone (e.g., Cluster 1, Cluster 2).
- The blue **“normal”** node at the top represents the **germline/root lineage**.
- Clones are inferred based on shared CCF (cancer cell fraction) patterns across samples.

### 2. Branch Labels (Numbers in Parentheses)
- Numbers such as `(2650)` or `(78)` represent the number of **private mutations** assigned to that branch.
- These mutations distinguish a child clone from its parent.
- ⚠️ These numbers **do not represent branch length, time, or evolutionary distance**.

### 3. Driver Annotations
- Gene names shown on branches (e.g., `SOX4`, `PIK3CA`) indicate **driver mutations** assigned to that evolutionary step.
- Only selected driver events are displayed; most mutations are not shown individually.

### 4. Tree Layout Does *Not* Represent Time
- Branch spacing, direction, and vertical position are purely for visualization.
- **There is no temporal information encoded in the tree layout.**
- Timing must be inferred from the **Timing module**, not from this plot.

---

## Timing Module (Somatic Event Ordering)

The Timing module produces a **directed acyclic graph (DAG)** to illustrate the **relative order** of somatic mutations and copy-number events.

### 1. Nodes
- Each node represents one or more genomic events whose timing estimates are **statistically indistinguishable**.

### 2. Edges (Arrows)
- An arrow between nodes indicates a **directional ordering** (e.g., Event Group A precedes Event Group B).
- Arrows are only drawn when the posterior probability supports the inferred direction.

### 3. No Quantitative Geometry
- Edge lengths, angles, and layout are **not proportional to time or evolutionary distance**.
- The graph is qualitative and only conveys ordering relationships.

### ✅ How to Interpret
- ✅ Use the **Cluster Tree** to understand clonal architecture and branching structure.
- ✅ Use the **Timing Module** to determine the relative order of somatic events.
- ❌ Do **not** interpret branch length, node height, or layout spacing as time.

## LeagueModel Output Figures

### 🟦 1. Log-Odds Plot (`<cohort>.log_odds.png` / `.pdf`)
This plot ranks events based on how **early vs. late** they tend to occur across the cohort.

**How to read:**
- **X-axis:**  
  Log-odds score comparing whether an event tends to occur earlier than others.  
  - More **negative** → occurs **earlier**  
  - More **positive** → occurs **later**

- **Right-side horizontal bar (one per event):**  
  - **Color:** Event class/type (e.g., SNV, WGD, arm-level CNV, focal CNV)  
  - **Length:** Event prevalence (number of samples where the event is observed)

✅ **Interpretation:**  
Events on the far **left** with long bars are **early and common**.  
Events on the **right** tend to be **later**.

---

### 🟩 2. Position (Rank) Plot (`<cohort>.positions.png` / `.pdf`)
This plot displays the **relative ordering (rank)** of events from earliest to latest.

**How to read:**
- **Event order (vertical or sorted layout):**  
  Events are arranged based on their **mean inferred rank** across samples.  
  - **Rank 1 = earliest**  
  - Higher rank = later

- **Right-side prevalence bar:**  
  Same as in the log-odds plot:  
  - **Color:** Event category  
  - **Length:** Event prevalence

✅ **Interpretation:**  
Use this plot to view the **relative timeline of events**, without relying on statistical magnitude.

---

### ✅ Summary of Differences

| Feature                | Log-Odds Plot                    | Position (Rank) Plot              |
|------------------------|----------------------------------|----------------------------------|
| X-axis meaning         | Log-odds score (early → late)    | Not numeric — rank-based ordering |
| Event ordering basis   | Strength of evidence for timing  | Relative mean inferred rank       |
| Bar on right           | Prevalence + type (color-coded)  | Same meaning                     |
| Best for               | Identifying statistically early/late events | Viewing overall ordering trends |

---

### Notes
- These figures describe **relative timing trends** across the entire cohort.