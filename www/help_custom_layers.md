
# How to Use the App

The **Explorer** tab allows you to visualize predicted species presence and interaction with human activities.
You can either explore the **default results** or run a **custom scenario** using your own data.

---

## 1. Explorer Modes

### Default mode
- Select a species and map type.
- Results are precomputed and displayed immediately.
- No user data is required.

### Custom scenario mode
This mode allows you to replace most environmental and fishing covariates with your own data.

**Important:** custom scenarios must follow the workflow below.

---

## 2. Custom Scenario Workflow (Required Order)

When using *Run custom scenario*:

1. **Upload your custom layer**  
   Upload an `.rds` file containing the variables you want to replace.

2. **Apply layer**  
   This validates your file and merges it with the prediction grid.

3. **Run prediction**  
   The species model is executed using the modified dataset.

Skipping or reordering steps will result in errors.

---

# Custom Covariate Layers – Data Requirements

## 3. Variables You Can Replace

You may replace **all variables except**:

- Port Proximity (`HW`)
- Urban Areas Proximity (`CW`)
- Rivers Proximity (`WR`)
- Depth (`Depth`)
- Distance from the weekend (`days_from_sunday`)

### Replaceable variables include:

**Environmental**
- `sst`, `so`, `chl`, `uo`, `vo`

**Human activity**
- `SD`, `FS`, `P`, `SS`

**Fishing effort**
- `trawlers`
- `tuna_purse_seines`
- `purse_seines`
- `other_purse_seines`
- `set_longlines`
- `fixed_gear`
- `set_gillnets`
- `pots_and_traps`

---

## 4. File Format (Strict)

Your file **must**:

- Be saved as an **RDS file**
- Contain **at least**:
  - `id` (integer)
  - `date` (Date)
- Use **exact variable names** as listed above
- Contain **only id/date combinations already present** in the dataset
- Use **numeric values** for all variables

Invalid IDs, dates, or column names will cause validation errors.

---

## 5. Recommended Preprocessing

To ensure compatibility with the models, we strongly recommend using the same units and transformations as the original dataset.

| Variable | Description | Unit |
|---------|------------|------|
| sst | Sea Surface Temperature | K |
| so | Salinity | PSU |
| chl | Chlorophyll-a | mg/m³ |
| uo | Eastward current | m/s |
| vo | Northward current | m/s |
| SD | Shipping density | – |
| FS | Fixed steel proximity | – |
| P | Floating platforms proximity | – |
| SS | Subsea steel proximity | – |
| Fishing variables | Apparent fishing hours | h |

### Additional recommendations

- Shipping density (`SD`) should be computed as a **weighted inverse squared distance** from shipping lanes.
- Installation proximity (`FS`, `P`, `SS`) should be computed as **weighted inverse distances** using importance indices.
- Missing fishing values should be replaced with `0`.
- Raw fishing effort should be **log10-transformed** before use.

