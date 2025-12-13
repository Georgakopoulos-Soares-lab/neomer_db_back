# NeomerDB Backend

This repository contains the backend service for **NeomerDB**, a system designed to analyze and serve genomic and exomic neomer (nullomer) data. The application is built using **Go (Golang)** and the **Gin** web framework, utilizing **DuckDB** as the high-performance analytical database engine.

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Configuration](#configuration)
3. [Installation & Usage](#installation--usage)
4. [API Reference](#api-reference)
   - [System & Metadata](#system--metadata)
   - [Genome Neomers](#genome-neomers)
   - [Exome Neomers](#exome-neomers)
   - [Patient & Analysis (Genome)](#patient--analysis-genome)
   - [Patient & Analysis (Exome)](#patient--analysis-exome)
   - [Statistical Distributions & Jaccard Indices](#statistical-distributions--jaccard-indices)

## Prerequisites

- **Go**: Version 1.18 or higher.
- **DuckDB**: The backend requires a valid DuckDB database file containing the neomer schemas (`neomers_{K}`, `exome_neomers_{K}`, `donor_data`, etc.).

## Configuration

The application relies on environment variables for configuration.

| Variable               | Description                                | Default Value          |
| :--------------------- | :----------------------------------------- | :--------------------- |
| `NEOMERS_DUCK_DB_FILE` | Absolute path to the DuckDB database file. | `/path_to/neomers.ddb` |

## Installation & Usage

1. **Clone the repository**

   ```bash
   git clone https://github.com/Georgakopoulos-Soares-lab/neomer_db_back
   cd neomer_db_back
   ```

2. **Install dependencies**

   ```bash
   go mod tidy
   ```

3. **Run the server**
   ```bash
   export NEOMERS_DUCK_DB_FILE="/path/to/your/database.ddb"
   go run main.go
   ```
   The server defaults to running on port `8080`.

## API Reference

All endpoints accept **GET** requests. The API supports Cross-Origin Resource Sharing (CORS) for all origins.

### System & Metadata

Basic endpoints for health checks and retrieving reference datasets.

#### `GET /healthcheck`

Returns the operational status of the API.

#### `GET /cancer_types`

Retrieves the full list of cancer types available in the database.

#### `GET /donor_data`

Retrieves metadata for all donors (genome dataset).

#### `GET /exomes_donor_data`

Retrieves metadata for all donors (exome dataset).

#### `GET /tcga_survival_data`

Retrieves survival analysis data associated with TCGA cohorts.

---

### Genome Neomers

Endpoints for querying neomers derived from whole-genome sequencing.

#### `GET /get_nullomers`

Retrieves a paginated list of neomers with associated metadata (GC content, frequencies, donor info).

**Parameters:**

- `length` (Required): Neomer length (e.g., `11` to `20`).
- `page`: Page number (0-indexed).
- `limit`: Rows per page.
- `filters`: SQL-like filter string (e.g., `AF < 0.01 AND gc_content > 30`).
- `specialFilters`: Specialized aggregation filters (e.g., `at_least_X_distinct_patients;3`).
- `column`, `filterType`, `value`: Used for range filtering (e.g., `filterType=between`, `value=0.1,0.5`).

#### `GET /get_suggestions`

Provides autocomplete suggestions for a specific column to assist UI filtering.

**Parameters:**

- `length` (Required): Neomer length.
- `column` (Required): The database column to search (e.g., `Hugo_Symbol`).
- `input`: The partial string to match.

#### `GET /get_nullomers_stats`

Returns aggregated statistics for neomers, supporting grouping and filtering.

**Parameters:**

- `length` (Required): Neomer length.
- `groupBy`: Comma-separated columns to group by.
- `topN`: Limit the number of returned groups (default: 10).
- `filters`: SQL-like filter string.

---

### Exome Neomers

Endpoints for querying neomers derived from exome sequencing.

#### `GET /get_exome_nullomers`

Retrieves a paginated list of exome-derived neomers.

**Parameters:**

- `length` (Required): Neomer length.
- `page`, `limit`: Pagination controls.
- `filters`: SQL-like filter string.
- `specialFilters`: Specialized aggregation filters.

#### `GET /get_exome_suggestions`

Provides autocomplete suggestions for exome data columns.

**Parameters:**

- `length` (Required): Neomer length.
- `column`: Column name.
- `input`: Search string.

#### `GET /get_exome_nullomers_stats`

Returns aggregated statistics for exome neomers.

**Parameters:**

- `length` (Required): Neomer length.
- `groupBy`: Columns to group by.
- `topN`: Limit results.

---

### Patient & Analysis (Genome)

Endpoints for analyzing specific patients or individual neomer sequences within the genome dataset.

#### `GET /patient_details`

Retrieves detailed clinical and metadata for a specific donor.

**Parameters:**

- `donor_id` (Required): The actual donor ID (e.g., ICGC Donor ID).

#### `GET /patient_neomers`

Retrieves all neomers found in a specific patient.

**Parameters:**

- `donor_id` (Required): The actual donor ID.
- `length` (Required): Neomer length.
- `top_n`: Limit results.
- `prefix`: Filter neomers by starting sequence.

#### `GET /analyze_neomer`

Performs a deep analysis of a single neomer sequence, returning prevalence across cancer types and organs.

**Parameters:**

- `neomer` (Required): The nucleotide sequence (e.g., `ACGT...`).

---

### Patient & Analysis (Exome)

Endpoints for analyzing specific patients or individual neomer sequences within the exome dataset.

#### `GET /exome_patient_details`

Retrieves detailed metadata for a specific exome donor.

**Parameters:**

- `donor_id` (Required): The BCR patient barcode.

#### `GET /exome_patient_neomers`

Retrieves neomers found in a specific exome donor.

**Parameters:**

- `donor_id` (Required): The BCR patient barcode.
- `length` (Required): Neomer length.
- `top_n`: Limit results.

#### `GET /exome_analyze_neomer`

Performs a deep analysis of a single exome neomer sequence, returning prevalence across cancer types and organs.

**Parameters:**

- `neomer` (Required): The nucleotide sequence.

---

### Statistical Distributions & Jaccard Indices

Endpoints for high-level statistical analysis of the dataset.

#### `GET /jaccard_index`

Calculates the Jaccard Similarity Index between pairwise **Cancer Types** based on shared neomers.

**Parameters:**

- `K` (Required): Neomer length.

#### `GET /jaccard_index_organs`

Calculates the Jaccard Similarity Index between pairwise **Organs** based on shared neomers.

**Parameters:**

- `K` (Required): Neomer length.

#### `GET /dataset_stats_cancer_types_varying_k`

Returns the count of neomers per cancer type for lengths K=11 through K=16.

#### `GET /distribution_neomer/:K/cancer_types`

Returns a list of distinct cancer types available for the distribution analysis at length K.

#### `GET /distribution_neomer/:K/organs`

Returns a list of distinct organs available for the distribution analysis at length K.

#### `GET /distribution_neomer/:K/data_by_cancer_type`

Returns distribution buckets (donor count vs. number of nullomers) for a specific cancer type.

**Parameters:**

- `cancerType` (Required): The cancer type identifier.

#### `GET /distribution_neomer/:K/data_by_organ`

Returns distribution buckets (donor count vs. number of nullomers) for a specific organ.

**Parameters:**

- `organ` (Required): The organ identifier.
