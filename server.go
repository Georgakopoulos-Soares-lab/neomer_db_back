package main

import (
    "database/sql"
    "fmt"
    "log"
    "net/http"
    "os"
    "strconv"
    "strings"
    "regexp"
    "github.com/gin-gonic/gin"
    _ "github.com/marcboeker/go-duckdb"
    "github.com/gin-contrib/cors"
    "time"
)

// Global database connection
var db *sql.DB

func main() {
    // Initialize global database connection once
    var err error
    db, err = sql.Open("duckdb", getDatabasePath())
    if err != nil {
        log.Fatalf("Failed to open database: %v", err)
    }
    defer db.Close()

    // Apply resource limits to prevent VM freezing
    _, err = db.Exec("SET memory_limit = '8GB'; SET threads = 3; SET temp_directory = '/tmp/duckdb_temp'; SET preserve_insertion_order = false;")
    if err != nil {
        log.Fatalf("Failed to set database resource limits: %v", err)
    }

    log.Println("Database initialized with resource limits: memory_limit=8GB, threads=3")
    router := gin.Default()

    // Disable CORS policy
        router.Use(cors.New(cors.Config{
                AllowOrigins:     []string{"*"},
                AllowMethods:     []string{"GET", "POST", "PUT", "PATCH", "DELETE", "OPTIONS"},
                AllowHeaders:     []string{"Origin", "Content-Type", "Accept", "Authorization"},
                ExposeHeaders:    []string{"Content-Length"},
                AllowCredentials: true,
                MaxAge:           12 * time.Hour,
        }))

    router.GET("/healthcheck", healthCheckHandler)
    router.GET("/cancer_types", makeHandler("SELECT * FROM cancer_types"))
    router.GET("/donor_data", makeHandler("SELECT * FROM donor_data"))
    router.GET("/exomes_donor_data", makeHandler("SELECT * FROM exome_donor_data"))

    router.GET("/tcga_survival_data", makeHandler("SELECT * FROM tcga_survival_data"))

    // Neomers
    router.GET("/get_nullomers", getNullomersHandler)
    router.GET("/get_suggestions", getSuggestionsHandler)
    router.GET("/get_nullomers_stats", getNullomersStatsHandler)

    // Exomes
    router.GET("/get_exome_nullomers",   getExomeNullomersHandler)
    router.GET("/get_exome_suggestions", getExomeSuggestionsHandler)
    router.GET("/get_exome_nullomers_stats", getExomeNullomersStatsHandler)

    // Patient Details
    router.GET("/patient_details", getPatientDetailsHandler)
    router.GET("/patient_neomers", getPatientNeomersHandler)
    router.GET("/analyze_neomer", analyzeNeomerHandler)
    


    router.GET("/jaccard_index", getJaccardIndexHandler)
    router.GET("/jaccard_index_organs", getJaccardIndexOrgansHandler)

    router.GET("/dataset_stats_cancer_types_varying_k", getDatasetStatsCancerTypesVaryingKHandler)
    
    router.GET("/distribution_neomer/:K/cancer_types", getDistNeomerKCancerTypes)
    router.GET("/distribution_neomer/:K/organs", getDistNeomerKOrgans)
    router.GET("/distribution_neomer/:K/data_by_cancer_type", getDistNeomerKDataByCancerType)
    router.GET("/distribution_neomer/:K/data_by_organ", getDistNeomerKDataByOrgan)
    
    router.GET("/exome_patient_details", getExomePatientDetailsHandler)
    router.GET("/exome_patient_neomers",  getExomePatientNeomersHandler)
    router.GET("/exome_analyze_neomer", getExomeAnalyzeNeomerHandler )


    
    if err := router.Run(); err != nil {
        log.Fatalf("Failed to run server: %v", err)
    }
}

// ------------------------------------------------------------------
// Basic Utility & Handlers
// ------------------------------------------------------------------

func corsMiddleware() gin.HandlerFunc {
    return func(c *gin.Context) {
        c.Writer.Header().Set("Access-Control-Allow-Origin", "*")
        c.Next()
    }
}

func healthCheckHandler(c *gin.Context) {
    c.JSON(http.StatusOK, gin.H{"status": "healthy"})
}

func makeHandler(query string) func(*gin.Context) {
    return func(c *gin.Context) {
        rows, err := db.Query(query)
        if err != nil {
            c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
            return
        }
        defer rows.Close()

        columns, err := rows.Columns()
        if err != nil {
            c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
            return
        }

        data := make([][]interface{}, 0)
        for rows.Next() {
            row := make([]interface{}, len(columns))
            rowPointers := make([]interface{}, len(columns))
            for i := range row {
                rowPointers[i] = &row[i]
            }

            if err := rows.Scan(rowPointers...); err != nil {
                c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
                return
            }

            // Convert []byte data to string, if applicable
            for i, val := range row {
                if b, ok := val.([]byte); ok {
                    row[i] = string(b)
                }
            }
            data = append(data, row)
        }

        result := map[string]interface{}{
            "headers": columns,
            "data":    data,
        }

        c.JSON(http.StatusOK, result)
    }
}

func getDatabasePath() string {
    if path := os.Getenv("NEOMERS_DUCK_DB_FILE"); path != "" {
        return path
    }
    // Default path if environment variable is not set
    return "/storage/group/izg5139/default/external/neo_database/staging.neomers.ddb"
}

// ------------------------------------------------------------------
// removeParentheses helper function
// ------------------------------------------------------------------
func removeParentheses(input string) string {
    re := regexp.MustCompile(`[()]`) // Matches ( and )
    return re.ReplaceAllString(input, "")
}


// ------------------------------------------------------------------
// getNullomersHandler
// ------------------------------------------------------------------
// Returns *all* columns from neomers_%s plus all columns from 
// cancer_type_details plus donor_data, plus an added column "gc_content."
//
// getNullomersHandler serves paginated nullomer records with support for “between” (AF-range) filters.
func getNullomersHandler(c *gin.Context) {
    length := c.Query("length")
    if length == "" {
        c.JSON(http.StatusBadRequest, gin.H{"error": "Missing required parameter 'length'"})
        return
    }

    // Pagination params
    pageStr := c.Query("page")
    limitStr := c.Query("limit")
    filters := c.Query("filters")             // e.g. "(gc_content > 10) AND (gc_content < 50)"
    specialFilters := c.Query("specialFilters") // e.g. "at_least_X_distinct_patients;3"

    // NEW: single “between” filter
    column := c.Query("column")       // e.g. "AF"
    filterType := c.Query("filterType") // should be "between"
    filterValue := c.Query("value")     // e.g. "0.10,0.50"

    // Default paging
    page := 0
    limit := 10000
    if p, err := strconv.Atoi(pageStr); err == nil && p >= 0 {
        page = p
    }
    if l, err := strconv.Atoi(limitStr); err == nil && l > 0 && l <= 10000 {
        limit = l
    }

    // 1) Base CTE
    baseQuery := fmt.Sprintf(`
        WITH base AS (
            SELECT
                n.* EXCLUDE (Donor_ID),
                c.*,
                d.*,
                di.Tumor_Sample_Barcode,
                di.Matched_Norm_Sample_Barcode,
                ROUND(
                    100.0 * (
                        LENGTH(n.nullomers_created)
                        - LENGTH(REPLACE(UPPER(n.nullomers_created), 'G', ''))
                        - LENGTH(REPLACE(UPPER(n.nullomers_created), 'C', ''))
                    ) / LENGTH(n.nullomers_created),
                    2
                ) * -1 AS gc_content
            FROM neomers_%[1]s n
            JOIN cancer_type_details c USING (Project_Code)
            LEFT JOIN donor_id_mapping di ON CAST(n."Donor_ID" AS INT) = di."Donor_ID"
            LEFT JOIN donor_data d ON di.Actual_Donor_ID = d.icgc_donor_id
        )
        SELECT * FROM base
    `, length)

    // Build WHERE clauses
    var whereClauses []string

    // a) “Between” filter for AF* columns
    if filterType == "between" && column != "" && filterValue != "" {
        parts := strings.Split(filterValue, ",")
        if len(parts) == 2 {
            minVal := parts[0]
            maxVal := parts[1]
            whereClauses = append(whereClauses,
                fmt.Sprintf(
                    `CAST("%s" AS FLOAT) >= %s AND CAST("%s" AS FLOAT) <= %s`,
                    column, minVal, column, maxVal,
                ),
            )
        }
    } else if filters != "" {
        // b) Legacy “filters” param
        for _, cond := range strings.Split(filters, " AND ") {
            parts := strings.Fields(cond)
            if len(parts) >= 3 {
                col := cleanColumnName(parts[0])
                if isNumericColumn(col) {
                    cond = fmt.Sprintf(`CAST("%s" AS FLOAT) %s %s`,
                        col, parts[1], removeParentheses(parts[2]))
                }
                whereClauses = append(whereClauses, cond)
            }
        }
    }

    // c) Special filters (e.g. at_least_X_distinct_patients)
    if specialFilters != "" {
        for _, part := range strings.Split(specialFilters, "|") {
            sf := strings.Split(part, ";")
            if sf[0] == "at_least_X_distinct_patients" && len(sf) == 2 {
                if n, err := strconv.Atoi(sf[1]); err == nil && n > 0 {
                    sub := fmt.Sprintf(`
                        nullomers_created IN (
                            SELECT nullomers_created
                            FROM neomers_%[1]s
                            JOIN cancer_type_details USING (Project_Code)
                            LEFT JOIN donor_id_mapping di ON CAST(n."Donor_ID" AS INT) = di."Donor_ID"
                            LEFT JOIN donor_data d ON di.Actual_Donor_ID = d.icgc_donor_id
                            GROUP BY nullomers_created
                            HAVING COUNT(DISTINCT donor_id) >= %d
                        )`, length, n)
                    whereClauses = append(whereClauses, sub)
                }
            }
        }
    }

    finalWhere := ""
    if len(whereClauses) > 0 {
        finalWhere = " WHERE " + strings.Join(whereClauses, " AND ")
    }
    fmt.Println("🔍  finalWhere", finalWhere)

    // 2) COUNT
    countQuery := fmt.Sprintf(`
        WITH base AS (
            SELECT
                n.*,
                c.*,
                d.*,
                ROUND(
                    100.0 * (
                        LENGTH(n.nullomers_created)
                        - LENGTH(REPLACE(UPPER(n.nullomers_created), 'G', ''))
                        - LENGTH(REPLACE(UPPER(n.nullomers_created), 'C', ''))
                    ) / LENGTH(n.nullomers_created),
                    2
                ) * -1 AS gc_content
            FROM neomers_%[1]s n
            JOIN cancer_type_details c USING (Project_Code)
            LEFT JOIN donor_id_mapping di ON CAST(n."Donor_ID" AS INT) = di."Donor_ID"
            LEFT JOIN donor_data d ON di.Actual_Donor_ID = d.icgc_donor_id
        )
        SELECT COUNT(*) FROM base
        %s
    `, length, finalWhere)

    var totalCount int
    if err := db.QueryRow(countQuery).Scan(&totalCount); err != nil {
        c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
        return
    }

    // 3) Data page
    offset := page * limit
    pageQuery := fmt.Sprintf("%s %s LIMIT %d OFFSET %d", baseQuery, finalWhere, limit, offset)
    rows, err := db.Query(pageQuery)
    if err != nil {
        c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
        return
    }
    defer rows.Close()

    cols, err := rows.Columns()
    if err != nil {
        c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
        return
    }

    var data [][]interface{}
    for rows.Next() {
        row := make([]interface{}, len(cols))
        ptrs := make([]interface{}, len(cols))
        for i := range row {
            ptrs[i] = &row[i]
        }
        if err := rows.Scan(ptrs...); err != nil {
            c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
            return
        }
        for i, v := range row {
            if b, ok := v.([]byte); ok {
                row[i] = string(b)
            }
            if v == nil {
                row[i] = nil
            }
        }
        data = append(data, row)
    }

    c.JSON(http.StatusOK, gin.H{
        "headers":    cols,
        "data":       data,
        "totalCount": totalCount,
    })
}
// ------------------------------------------------------------------
// cleanColumnName helper function
// ------------------------------------------------------------------
func cleanColumnName(column string) string {
    re := regexp.MustCompile(`[^a-zA-Z0-9_]`) // Allow only letters, numbers, and underscores
    cleaned := re.ReplaceAllString(column, "") // Remove unwanted characters
    return strings.TrimSpace(cleaned)          // Trim whitespace
}

// ------------------------------------------------------------------
// isNumericColumn hepler function to check if a column is numeric
// ------------------------------------------------------------------
func isNumericColumn(column string) bool {
    var columnTypes = map[string]string{
        "donor_age_at_diagnosis": "BIGINT",
        "gc_content": "FLOAT",
        "nullomers_created": "VARCHAR",
        "donor_id": "VARCHAR",
        "AF": "FLOAT",
        "AF_eas": "FLOAT",
        "AF_afr": "FLOAT",
        "AF_fin": "FLOAT",
        "AF_ami": "FLOAT",
        "AF_amr": "FLOAT",
        "AF_nfe": "FLOAT",
        "AF_sas": "FLOAT",
        "AF_asj": "FLOAT",
    }
    numericTypes := map[string]bool{
        "BIGINT": true,
        "INTEGER": true,
        "FLOAT": true,
        "DOUBLE": true,
    }
    column = cleanColumnName(column) // Ensure input is clean


    colType, exists := columnTypes[column]
    return exists && numericTypes[colType]
}

// ------------------------------------------------------------------
// getSuggestionsHandler
// ------------------------------------------------------------------
// ------------------------------------------------------------------
// getSuggestionsHandler
// ------------------------------------------------------------------
func getSuggestionsHandler(c *gin.Context) {
	column := c.Query("column")
	input := c.Query("input")
	length := c.Query("length")

	if length == "" {
		c.JSON(http.StatusBadRequest, gin.H{"error": "Missing required parameter 'length'"})
		return
	}
	if column == "" {
		c.JSON(http.StatusBadRequest, gin.H{"error": "Missing required parameter 'column'"})
		return
	}

	// It's not meaningful to get suggestions for a purely numeric column like this
	if column == "gc_content" {
		c.JSON(http.StatusOK, gin.H{"suggestions": []string{}})
		return
	}

	lowerInput := strings.ToLower(input)
	var whereClause string
	if lowerInput != "" {
		// Use CAST to VARCHAR to safely perform a LIKE search on any column type
		whereClause = fmt.Sprintf("WHERE LOWER(CAST(%s AS VARCHAR)) LIKE '%%%s%%'", column, lowerInput)
	}

	// This query now joins all relevant tables, allowing suggestions from any of them.
	query := fmt.Sprintf(`
        WITH base AS (
            SELECT DISTINCT %s
            FROM neomers_%s n
            JOIN cancer_type_details c USING (Project_Code)
            LEFT JOIN donor_id_mapping di ON CAST(n."Donor_ID" AS INT) = di."Donor_ID"
            LEFT JOIN donor_data d ON di.Actual_Donor_ID = d.icgc_donor_id
            %s
            LIMIT 10
        )
        SELECT %s
        FROM base
		WHERE %s IS NOT NULL
        ORDER BY LOWER(CAST(%s AS VARCHAR)) ASC
    `, column, length, whereClause, column, column, column)

	rows, err := db.Query(query)
	if err != nil {
		// Log the error and return empty suggestions to prevent frontend issues
		log.Printf("Suggestion query failed for column '%s': %v", column, err)
		c.JSON(http.StatusOK, gin.H{"suggestions": []string{}})
		return
	}
	defer rows.Close()

	var suggestions []string
	for rows.Next() {
		var val interface{}
		if err := rows.Scan(&val); err == nil {
			if v, ok := val.(string); ok {
				suggestions = append(suggestions, v)
			} else if val != nil {
				// Convert non-string types to their string representation for the response
				suggestions = append(suggestions, fmt.Sprintf("%v", val))
			}
		}
	}
	c.JSON(http.StatusOK, gin.H{"suggestions": suggestions})
}

// ------------------------------------------------------------------
// getNullomersStatsHandler
// ------------------------------------------------------------------
func getNullomersStatsHandler(c *gin.Context) {
    length := c.Query("length")
    if length == "" {
        c.JSON(http.StatusBadRequest, gin.H{"error": "Missing required parameter 'length'"})
        return
    }

    filters := c.Query("filters")
    specialFilters := c.Query("specialFilters")
    groupByStr := c.Query("groupBy")
    topNStr := c.Query("topN")

    // NEW: single “between” filter
    column := c.Query("column")
    filterType := c.Query("filterType")
    filterValue := c.Query("value")

    topN := 10
    if n, err := strconv.Atoi(topNStr); err == nil && n > 0 {
        topN = n
    }

    // Base CTE
    baseCTE := fmt.Sprintf(`
        WITH base AS (
            SELECT
                n.* EXCLUDE (Donor_ID),
                c.*,
                d.*,
                di.Tumor_Sample_Barcode,
                di.Matched_Norm_Sample_Barcode,
                ROUND(
                    100.0 * (
                        LENGTH(n.nullomers_created)
                        - LENGTH(REPLACE(UPPER(n.nullomers_created), 'G', ''))
                        - LENGTH(REPLACE(UPPER(n.nullomers_created), 'C', ''))
                    ) / LENGTH(n.nullomers_created),
                    2
                ) * -1 AS gc_content
            FROM neomers_%[1]s n
            JOIN cancer_type_details c USING (Project_Code)
            LEFT JOIN donor_id_mapping di ON CAST(n."Donor_ID" AS INT) = di."Donor_ID"
            LEFT JOIN donor_data d ON di.Actual_Donor_ID = d.icgc_donor_id
        )
    `, length)

    // Build WHERE clauses
    var whereClauses []string

    // a) “Between” filter
    if filterType == "between" && column != "" && filterValue != "" {
        parts := strings.Split(filterValue, ",")
        if len(parts) == 2 {
            minVal := parts[0]
            maxVal := parts[1]
            whereClauses = append(whereClauses,
                fmt.Sprintf(
                    `CAST("%s" AS FLOAT) >= %s AND CAST("%s" AS FLOAT) <= %s`,
                    column, minVal, column, maxVal,
                ),
            )
        }
    } else if filters != "" {
        // b) Legacy filters
        for _, cond := range strings.Split(filters, " AND ") {
            parts := strings.Fields(cond)
            if len(parts) >= 3 {
                col := cleanColumnName(parts[0])
                if isNumericColumn(col) {
                    cond = fmt.Sprintf(`CAST("%s" AS FLOAT) %s %s`,
                        col, parts[1], removeParentheses(parts[2]))
                }
                whereClauses = append(whereClauses, cond)
            }
        }
    }

    // c) Special filters
    if specialFilters != "" {
        for _, part := range strings.Split(specialFilters, "|") {
            sf := strings.Split(part, ";")
            if sf[0] == "at_least_X_distinct_patients" && len(sf) == 2 {
                if n, err := strconv.Atoi(sf[1]); err == nil && n > 0 {
                    sub := fmt.Sprintf(`
                        nullomers_created IN (
                            SELECT nullomers_created
                            FROM neomers_%[1]s
                            JOIN cancer_type_details USING (Project_Code)
                            LEFT JOIN donor_id_mapping di ON CAST(n."Donor_ID" AS INT) = di."Donor_ID"
                            LEFT JOIN donor_data d ON di.Actual_Donor_ID = d.icgc_donor_id
                            GROUP BY nullomers_created
                            HAVING COUNT(DISTINCT donor_id) >= %d
                        )`, length, n)
                    whereClauses = append(whereClauses, sub)
                }
            }
        }
    }

    finalWhere := ""
    if len(whereClauses) > 0 {
        finalWhere = " WHERE " + strings.Join(whereClauses, " AND ")
    }

    // Prepare GROUP BY
    groupByCols := []string{"nullomers_created"}
    selectCols := []string{"nullomers_created"}
    if groupByStr != "" {
        for _, col := range strings.Split(groupByStr, ",") {
            col = strings.TrimSpace(col)
            if col != "" {
                groupByCols = append(groupByCols, col)
                selectCols = append(selectCols, col)
            }
        }
    }
    selectClause := strings.Join(selectCols, ", ")
    groupByClause := strings.Join(groupByCols, ", ")

    // Final stats query
    query := fmt.Sprintf(`
        %s
        SELECT
            %s,
            COUNT(*) AS total_count
        FROM base
        %s
        GROUP BY %s
        ORDER BY total_count DESC
        LIMIT %d
    `, baseCTE, selectClause, finalWhere, groupByClause, topN)

    rows, err := db.Query(query)
    if err != nil {
        c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
        return
    }
    defer rows.Close()

    cols, err := rows.Columns()
    if err != nil {
        c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
        return
    }

    var data [][]interface{}
    for rows.Next() {
        row := make([]interface{}, len(cols))
        ptrs := make([]interface{}, len(cols))
        for i := range row {
            ptrs[i] = &row[i]
        }
        if err := rows.Scan(ptrs...); err != nil {
            c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
            return
        }
        for i, v := range row {
            if b, ok := v.([]byte); ok {
                row[i] = string(b)
            }
            if v == nil {
                row[i] = nil
            }
        }
        data = append(data, row)
    }

    c.JSON(http.StatusOK, gin.H{
        "headers": cols,
        "data":    data,
    })
}

// GET /exome_nullomers?length=<L>&page=<P>&limit=<N>&filters=…&specialFilters=…&column=…&filterType=between&value=…
func getExomeNullomersHandler(c *gin.Context) {
    length := c.Query("length")
    if length == "" {
        c.JSON(http.StatusBadRequest, gin.H{"error": "Missing required parameter 'length'"})
        return
    }
    page, _ := strconv.Atoi(c.DefaultQuery("page", "0"))
    limit, _ := strconv.Atoi(c.DefaultQuery("limit", "10000"))
    if limit < 1 || limit > 10000 {
        limit = 10000
    }

    filters := c.Query("filters")
    specialFilters := c.Query("specialFilters")

    // Build WHERE parts
    var whereParts []string
    if filters != "" {
        // For any numeric comparison, cast column to DOUBLE
        numericCols := []string{
            "gc_content",
            "AF", "AF_eas", "AF_afr", "AF_fin", "AF_ami",
            "AF_amr", "AF_nfe", "AF_sas", "AF_asj",
            "days_to_birth", "days_to_death", "days_to_last_followup",
        }
        pattern := `\b(` + strings.Join(numericCols, "|") + `)\s*(>=|<=|>|<|=)\s*([0-9]+(?:\.[0-9]+)?)`
        re := regexp.MustCompile(pattern)
        filtersCast := re.ReplaceAllString(filters, `CAST($1 AS DOUBLE) $2 $3`)
        whereParts = append(whereParts, filtersCast)
    }
    if specialFilters != "" {
        for _, part := range strings.Split(specialFilters, "|") {
            sf := strings.Split(part, ";")
            if len(sf) == 2 && sf[0] == "at_least_X_distinct_patients" {
                if n, err := strconv.Atoi(sf[1]); err == nil && n > 0 {
                    sub := fmt.Sprintf(
                        `"nullomers_created" IN (
                            SELECT "nullomers_created"
                            FROM exome_neomers_%s
                            GROUP BY "nullomers_created"
                            HAVING COUNT(DISTINCT "Donor_ID") >= %d
                         )`, length, n)
                    whereParts = append(whereParts, sub)
                }
            }
        }
    }

    finalWhere := ""
    if len(whereParts) > 0 {
        finalWhere = " WHERE " + strings.Join(whereParts, " AND ")
    }

    // 1) CTE definition
    baseCTE := fmt.Sprintf(`
WITH base AS (
    SELECT
        n.* EXCLUDE(Donor_ID),
        di.Tumor_Sample_Barcode,
        di.Matched_Norm_Sample_Barcode,
        d.*,
        ROUND(
            100.0 * (
                LENGTH(n.nullomers_created)
                - LENGTH(REPLACE(UPPER(n.nullomers_created), 'G', ''))
                - LENGTH(REPLACE(UPPER(n.nullomers_created), 'C', ''))
            ) / LENGTH(n.nullomers_created),
            2
        ) * -1 AS gc_content
    FROM exome_neomers_%s AS n
    JOIN exomes_donor_id_mapping AS di
        ON CAST(n."Donor_ID" AS INT) = di."Donor_ID"
    LEFT JOIN exome_donor_data AS d
        ON di.Actual_Donor_ID = d.bcr_patient_barcode
)
`, length)

    // 2) Total count
    countQ := baseCTE + "SELECT COUNT(*) FROM base" + finalWhere

    var totalCount int
    if err := db.QueryRow(countQ).Scan(&totalCount); err != nil {
        c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
        return
    }

    // 3) Page of data
    offset := page * limit
    pageQ := baseCTE +
        "SELECT * FROM base" + finalWhere +
        fmt.Sprintf(" LIMIT %d OFFSET %d", limit, offset)

    rows, err := db.Query(pageQ)
    if err != nil {
        c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
        return
    }
    defer rows.Close()

    cols, _ := rows.Columns()
    var data [][]interface{}
    for rows.Next() {
        vals := make([]interface{}, len(cols))
        ptrs := make([]interface{}, len(cols))
        for i := range vals {
            ptrs[i] = &vals[i]
        }
        if err := rows.Scan(ptrs...); err != nil {
            c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
            return
        }
        for i, v := range vals {
            if b, ok := v.([]byte); ok {
                vals[i] = string(b)
            }
            if v == nil {
                vals[i] = nil
            }
        }
        data = append(data, vals)
    }

    c.JSON(http.StatusOK, gin.H{
        "headers":    cols,
        "data":       data,
        "totalCount": totalCount,
    })
}


// getExomeSuggestionsHandler returns autocomplete suggestions for exome nullomer filters
func getExomeSuggestionsHandler(c *gin.Context) {
    column := c.Query("column")
    input := c.Query("input")
    length := c.Query("length")
    if length == "" {
        c.JSON(http.StatusBadRequest, gin.H{"error": "Missing length"})
        return
    }
    if column == "gc_content" {
        c.JSON(http.StatusOK, gin.H{"suggestions": []string{}})
        return
    }

    lowerInput := strings.ToLower(input)
    cond := ""
    if lowerInput != "" {
        cond = fmt.Sprintf("WHERE LOWER(%s) LIKE '%%%s%%'", column, lowerInput)
    }

    query := fmt.Sprintf(`
        WITH base AS (
            SELECT DISTINCT %s
            FROM exome_neomers_%s n
            JOIN exomes_donor_id_mapping di ON CAST(n."Donor_ID" AS INT) = di."Donor_ID"
            JOIN exome_donor_data d ON di.Actual_Donor_ID = d.bcr_patient_barcode
            %s
            LIMIT 10
        )
        SELECT %s
        FROM base
        ORDER BY LOWER(%s) ASC
    `, column, length, cond, column, column)

    rows, err := db.Query(query)
    if err != nil {
        c.JSON(http.StatusOK, gin.H{"suggestions": []string{}})
        return
    }
    defer rows.Close()

    var suggestions []string
    for rows.Next() {
        var val interface{}
        if err := rows.Scan(&val); err != nil {
            continue
        }
        if s, ok := val.(string); ok {
            suggestions = append(suggestions, s)
        } else if val != nil {
            suggestions = append(suggestions, fmt.Sprintf("%v", val))
        }
    }

    c.JSON(http.StatusOK, gin.H{"suggestions": suggestions})
}

// GET /exome_nullomers_stats?length=<L>&filters=…&specialFilters=…&groupBy=…&topN=…&column=…&filterType=between&value=…
func getExomeNullomersStatsHandler(c *gin.Context) {
    length := c.Query("length")
    if length == "" {
        c.JSON(http.StatusBadRequest, gin.H{"error": "Missing required parameter 'length'"})
        return
    }
    filters := c.Query("filters")
    specialFilters := c.Query("specialFilters")
    groupByStr := c.Query("groupBy")
    topNStr := c.Query("topN")

    topN := 10
    if n, err := strconv.Atoi(topNStr); err == nil && n > 0 {
        topN = n
    }

    // Build WHERE parts (with casting for numeric columns)...
    var whereParts []string
    if filters != "" {
        numericCols := []string{
            "gc_content",
            "AF", "AF_eas", "AF_afr", "AF_fin", "AF_ami",
            "AF_amr", "AF_nfe", "AF_sas", "AF_asj",
            "days_to_birth", "days_to_death", "days_to_last_followup",
        }
        pattern := `\b(` + strings.Join(numericCols, "|") + `)\s*(>=|<=|>|<|=)\s*([0-9]+(?:\.[0-9]+)?)`
        re := regexp.MustCompile(pattern)
        filtersCast := re.ReplaceAllString(filters, `CAST($1 AS DOUBLE) $2 $3`)
        whereParts = append(whereParts, filtersCast)
    }
    if specialFilters != "" {
        for _, part := range strings.Split(specialFilters, "|") {
            sf := strings.Split(part, ";")
            if len(sf) == 2 && sf[0] == "at_least_X_distinct_patients" {
                if n, err := strconv.Atoi(sf[1]); err == nil && n > 0 {
                    sub := fmt.Sprintf(
                        `"nullomers_created" IN (
                            SELECT "nullomers_created"
                            FROM exome_neomers_%s
                            GROUP BY "nullomers_created"
                            HAVING COUNT(DISTINCT "Donor_ID") >= %d
                         )`, length, n)
                    whereParts = append(whereParts, sub)
                }
            }
        }
    }
    finalWhere := ""
    if len(whereParts) > 0 {
        finalWhere = " WHERE " + strings.Join(whereParts, " AND ")
    }

    // 1) Base CTE: include ALL donor‐data columns (d.*) plus mapping cols
    baseCTE := fmt.Sprintf(`
WITH base AS (
    SELECT
        n.* EXCLUDE(Donor_ID),
        di.Tumor_Sample_Barcode,
        di.Matched_Norm_Sample_Barcode,
        d.*,  -- bring in every column from exome_donor_data
        ROUND(
            100.0 * (
                LENGTH(n.nullomers_created)
                - LENGTH(REPLACE(UPPER(n.nullomers_created), 'G', ''))
                - LENGTH(REPLACE(UPPER(n.nullomers_created), 'C', ''))
            ) / LENGTH(n.nullomers_created),
            2
        ) * -1 AS gc_content
    FROM exome_neomers_%s AS n
    JOIN exomes_donor_id_mapping AS di
        ON CAST(n."Donor_ID" AS INT) = di."Donor_ID"
    LEFT JOIN exome_donor_data AS d
        ON di.Actual_Donor_ID = d.bcr_patient_barcode
)
`, length)

    // 2) Build GROUP BY / SELECT
    groupByCols := []string{"nullomers_created"}
    selectCols := []string{"nullomers_created"}
    if groupByStr != "" {
        for _, col := range strings.Split(groupByStr, ",") {
            col = strings.TrimSpace(col)
            if col != "" {
                groupByCols = append(groupByCols, col)
                selectCols = append(selectCols, col)
            }
        }
    }
    selectClause := strings.Join(selectCols, ", ")
    groupByClause  := strings.Join(groupByCols, ", ")

    // 3) Stats query
    statsQ := baseCTE +
        fmt.Sprintf(
            "SELECT %s, COUNT(*) AS total_count FROM base%s GROUP BY %s ORDER BY total_count DESC LIMIT %d",
            selectClause, finalWhere, groupByClause, topN,
        )

    rows, err := db.Query(statsQ)
    if err != nil {
        c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
        return
    }
    defer rows.Close()

    cols, _ := rows.Columns()
    var data [][]interface{}
    for rows.Next() {
        vals := make([]interface{}, len(cols))
        ptrs := make([]interface{}, len(cols))
        for i := range vals {
            ptrs[i] = &vals[i]
        }
        if err := rows.Scan(ptrs...); err != nil {
            c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
            return
        }
        for i, v := range vals {
            if b, ok := v.([]byte); ok {
                vals[i] = string(b)
            }
            if v == nil {
                vals[i] = nil
            }
        }
        data = append(data, vals)
    }

    c.JSON(http.StatusOK, gin.H{
        "headers": cols,
        "data":    data,
    })
}





// ------------------------------------------------------------------
// getPatientDetailsHandler
// ------------------------------------------------------------------
func getPatientDetailsHandler(c *gin.Context) {
    // 1) Read the actual donor string
    actualID := c.Query("donor_id")
    if actualID == "" {
        c.JSON(http.StatusBadRequest, gin.H{"error": "Missing donor_id"})
        return
    }

    // 2) Map Actual_Donor_ID → internal Donor_ID
    var internalID int
    var err error
    err = db.QueryRow(`
        SELECT Donor_ID 
        FROM donor_id_mapping 
        WHERE Actual_Donor_ID = ?
    `, actualID).Scan(&internalID)
    if err != nil {
        // no mapping found → no patient
        c.JSON(http.StatusOK, gin.H{"patient": nil})
        return
    }

    // 3) Grab one Project_Code from neomers_15
    var projectCode string
    err = db.QueryRow(`
        SELECT Project_Code
        FROM neomers_15
        WHERE Donor_ID = ?
        LIMIT 1
    `, internalID).Scan(&projectCode)
    if err != nil {
        c.JSON(http.StatusOK, gin.H{"patient": nil})
        return
    }

    // 4) Now join donor_data + cancer_type_details on that Project_Code
    query := `
        SELECT d.*, c.Cancer_Type, c.Organ
        FROM donor_data d
        JOIN cancer_type_details c
          ON c.Project_Code = ?
        WHERE d.icgc_donor_id = ?
        LIMIT 1
    `
    row := db.QueryRow(query, projectCode, actualID)

    // --- rest is identical to before: dynamically fetch column names & scan ---

    columns, err := db.Query("SELECT * FROM donor_data LIMIT 0")
    if err != nil {
        c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
        return
    }
    colNames, _ := columns.Columns()
    columns.Close()

    // add our two extra fields
    colNames = append(colNames, "Cancer_Type", "Organ")

    vals := make([]interface{}, len(colNames))
    valPtrs := make([]interface{}, len(colNames))
    for i := range vals {
        valPtrs[i] = &vals[i]
    }

    if err := row.Scan(valPtrs...); err != nil {
        c.JSON(http.StatusOK, gin.H{"patient": nil})
        return
    }

    patientMap := make(map[string]interface{}, len(colNames))
    for i, name := range colNames {
        v := vals[i]
        if b, ok := v.([]byte); ok {
            patientMap[name] = string(b)
        } else {
            patientMap[name] = v
        }
    }

    c.JSON(http.StatusOK, gin.H{"patient": patientMap})
}


// ------------------------------------------------------------------
// getPatientNeomersHandler
// ------------------------------------------------------------------
func getPatientNeomersHandler(c *gin.Context) {
    donorID := c.Query("donor_id")
    lengthStr := c.Query("length")
    topNStr := c.Query("top_n")
    prefix := c.Query("prefix") // Optional search prefix

    if donorID == "" {
        c.JSON(http.StatusBadRequest, gin.H{"error": "Missing donor_id"})
        return
    }
    if lengthStr == "" {
        c.JSON(http.StatusBadRequest, gin.H{"error": "Missing length"})
        return
    }

    length, err := strconv.Atoi(lengthStr)
    if err != nil || length < 11 || length > 20 {
        c.JSON(http.StatusBadRequest, gin.H{"error": "Invalid length. Must be between 11 and 20."})
        return
    }

    topN := 10
    if topNStr != "" {
        if n, err := strconv.Atoi(topNStr); err == nil && n > 0 {
            topN = n
        }
    }


    tableName := fmt.Sprintf("neomers_%d", length)

    // ←––––– HERE’S THE ONLY CHANGE ––––––→
    baseQuery := fmt.Sprintf(`
        SELECT 
            n.nullomers_created AS neomer, 
            COUNT(*)            AS count
        FROM %s n
        JOIN cancer_type_details c         USING (Project_Code)
        JOIN donor_id_mapping   di        ON CAST(n."Donor_ID" AS INT) = di."Donor_ID"
        JOIN donor_data         d         ON di.Actual_Donor_ID = d.icgc_donor_id
        WHERE di.Actual_Donor_ID = ?
    `, tableName)

    var args []interface{}
    args = append(args, donorID)

    if prefix != "" {
        baseQuery += " AND n.nullomers_created LIKE ?"
        args = append(args, prefix+"%")
    }

    baseQuery += `
        GROUP BY neomer
        ORDER BY count DESC
        LIMIT ?
    `
    args = append(args, topN)

    rows, err := db.Query(baseQuery, args...)
    if err != nil {
        c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
        return
    }
    defer rows.Close()

    var result []map[string]interface{}
    for rows.Next() {
        var neomer string
        var count int
        if err := rows.Scan(&neomer, &count); err == nil {
            result = append(result, map[string]interface{}{
                "neomer": neomer,
                "count":  count,
            })
        }
    }
    c.JSON(http.StatusOK, gin.H{"neomers": result})
}



// analyzeNeomerHandler
// ------------------------------------------------------------------
//
// Endpoint: /analyze_neomer?neomer=ABCDEF
//
// Returns basic stats about the neomer across donors & cancer types,
// including a breakdown by Cancer_Type and Organ.
//
// ------------------------------------------------------------------

func analyzeNeomerHandler(c *gin.Context) {
	neomer := c.Query("neomer")
	if neomer == "" {
		c.JSON(http.StatusBadRequest, gin.H{"error": "Missing neomer parameter"})
		return
	}
	K := len(neomer)
	if K < 11 || K > 20 {
		c.JSON(http.StatusBadRequest, gin.H{"error": "Neomer length must be between 11 and 20"})
		return
	}

	tableName := fmt.Sprintf("neomers_%d", K)

	// —— First query: overall stats ——
	totalQuery := fmt.Sprintf(`
		SELECT
			COUNT(*)                          AS total_count,
			COUNT(DISTINCT di.Actual_Donor_ID) AS distinct_donors,
			COUNT(DISTINCT c.Cancer_Type)     AS distinct_cancer_types,
			COUNT(DISTINCT c.Organ)           AS distinct_organs
		FROM %s n
		JOIN cancer_type_details c         USING (Project_Code)
		JOIN donor_id_mapping   di        ON CAST(n."Donor_ID" AS INT) = di."Donor_ID"
		/* JOIN donor_data         d         ON di.Actual_Donor_ID = d.icgc_donor_id -- Not needed for these counts */
		WHERE n.nullomers_created = ?
	`, tableName)

	var totalCount, distinctDonors, distinctCancerTypes, distinctOrgans int
	if err := db.QueryRow(totalQuery, neomer).Scan(
		&totalCount, &distinctDonors, &distinctCancerTypes, &distinctOrgans,
	); err != nil {
		c.JSON(http.StatusInternalServerError, gin.H{"error": "Error fetching total stats: " + err.Error()})
		return
	}

	// —— Second query: breakdown ——
	breakdownQuery := fmt.Sprintf(`
		SELECT
			c.Cancer_Type,
			c.Organ,
			COUNT(*) AS count
		FROM %s n
		JOIN cancer_type_details c         USING (Project_Code)
		JOIN donor_id_mapping   di        ON CAST(n."Donor_ID" AS INT) = di."Donor_ID"
		/* JOIN donor_data         d         ON di.Actual_Donor_ID = d.icgc_donor_id -- Not needed for this breakdown */
		WHERE n.nullomers_created = ?
		GROUP BY c.Cancer_Type, c.Organ
	`, tableName)

	rows, err := db.Query(breakdownQuery, neomer)
	if err != nil {
		c.JSON(http.StatusInternalServerError, gin.H{"error": "Error fetching breakdown stats: " + err.Error()})
		return
	}
	defer rows.Close()

	type OrganCount struct {
		Organ string `json:"organ"`
		Count int    `json:"count"`
	}
	type CancerTypeCount struct {
		CancerType string       `json:"cancerType"`
		Count      int          `json:"count"`
		Organs     []OrganCount `json:"organs"`
	}

	cancerMap := make(map[string]*CancerTypeCount)
	for rows.Next() {
		var ct, organ string
		var cnt int
		if err := rows.Scan(&ct, &organ, &cnt); err != nil {
			c.JSON(http.StatusInternalServerError, gin.H{"error": "Error scanning breakdown row: " + err.Error()})
			return
		}
		if _, ok := cancerMap[ct]; !ok {
			cancerMap[ct] = &CancerTypeCount{CancerType: ct}
		}
		cancerMap[ct].Count += cnt
		cancerMap[ct].Organs = append(cancerMap[ct].Organs, OrganCount{Organ: organ, Count: cnt})
	}

	var breakdown []CancerTypeCount
	for _, v := range cancerMap {
		breakdown = append(breakdown, *v)
	}

	// —— Third query: distinct actual donor IDs ——
	distinctDonorIDsQuery := fmt.Sprintf(`
		SELECT DISTINCT di.Actual_Donor_ID
		FROM %s n
		JOIN donor_id_mapping di ON CAST(n."Donor_ID" AS INT) = di."Donor_ID"
		WHERE n.nullomers_created = ?
	`, tableName)

	donorRows, err := db.Query(distinctDonorIDsQuery, neomer)
	if err != nil {
		c.JSON(http.StatusInternalServerError, gin.H{"error": "Error fetching distinct donor IDs: " + err.Error()})
		return
	}
	defer donorRows.Close()

	var distinctDonorIDsList []string
	for donorRows.Next() {
		var actualDonorID string
		if err := donorRows.Scan(&actualDonorID); err != nil {
			c.JSON(http.StatusInternalServerError, gin.H{"error": "Error scanning distinct donor ID row: " + err.Error()})
			return
		}
		distinctDonorIDsList = append(distinctDonorIDsList, actualDonorID)
	}
	if err := donorRows.Err(); err != nil { // Check for errors during iteration
		c.JSON(http.StatusInternalServerError, gin.H{"error": "Error iterating distinct donor IDs: " + err.Error()})
		return
	}


	c.JSON(http.StatusOK, gin.H{
		"analysis": map[string]interface{}{
			"totalNeomers":        totalCount,
			"distinctDonors":      distinctDonors, // This is count of distinct Actual_Donor_ID
			"distinctCancerTypes": distinctCancerTypes,
			"distinctOrgans":      distinctOrgans,
			"cancerBreakdown":     breakdown,
			"distinctDonorIDs":    distinctDonorIDsList, // Added this list
		},
	})
}



// ------------------------------------------------------------------
// getJaccardIndexHandler
// ------------------------------------------------------------------
//
// Computes Jaccard Index for each pair of Cancer_Types based on 
// shared nullomers.
//
func getJaccardIndexHandler(c *gin.Context) {
    // Retrieve and validate the 'K' parameter
    K := c.Query("K")
    if K == "" {
        c.JSON(http.StatusBadRequest, gin.H{"error": "Missing parameter 'K'"})
        return
    }
    // Validate that K is a positive integer
    if _, err := strconv.Atoi(K); err != nil || K == "0" {
        c.JSON(http.StatusBadRequest, gin.H{"error": "Parameter 'K' must be a positive integer"})
        return
    }

    // Construct the table name safely
    tableName := fmt.Sprintf("neomers_%s", K)

    // Define the SQL query to compute Jaccard indices for all pairs,
    //  (though for Jaccard across cancer types,
    // the main references remain nullomers + cancer_type_details).
    //
    // The extra join ensures design consistency, but typically doesn't
    // affect the Jaccard logic for Cancer_Type. It's just "the same pattern."
    query := fmt.Sprintf(`
        WITH joined_data AS (
            SELECT n.nullomers_created, c.Cancer_Type
            FROM %s n
            JOIN cancer_type_details c USING (Project_Code)
        ),
        cancer_counts AS (
            SELECT Cancer_Type, COUNT(DISTINCT nullomers_created) AS count
            FROM joined_data
            GROUP BY Cancer_Type
        ),
        all_cancer_types AS (
            SELECT DISTINCT Cancer_Type
            FROM cancer_counts
        ),
        pairs AS (
            SELECT a.Cancer_Type AS Cancer_Type_A, b.Cancer_Type AS Cancer_Type_B
            FROM all_cancer_types a
            CROSS JOIN all_cancer_types b
        ),
        intersections AS (
            SELECT 
                jd1.Cancer_Type AS Cancer_Type_A, 
                jd2.Cancer_Type AS Cancer_Type_B, 
                COUNT(DISTINCT jd1.nullomers_created) AS intersection_count
            FROM joined_data jd1
            JOIN joined_data jd2 ON jd1.nullomers_created = jd2.nullomers_created
            GROUP BY jd1.Cancer_Type, jd2.Cancer_Type
        )
        SELECT 
            p.Cancer_Type_A, 
            p.Cancer_Type_B, 
            COALESCE(i.intersection_count, 0) AS intersection_count,
            (c1.count + c2.count - COALESCE(i.intersection_count, 0)) AS union_count,
            CASE 
                WHEN p.Cancer_Type_A = p.Cancer_Type_B THEN 1.0
                WHEN (c1.count + c2.count - COALESCE(i.intersection_count, 0)) = 0 THEN 0.0
                ELSE ROUND(
                    CAST(COALESCE(i.intersection_count, 0) AS DOUBLE) 
                    / (c1.count + c2.count - COALESCE(i.intersection_count, 0)), 
                    4
                )
            END AS jaccard_index
        FROM pairs p
        LEFT JOIN intersections i 
            ON p.Cancer_Type_A = i.Cancer_Type_A 
            AND p.Cancer_Type_B = i.Cancer_Type_B
        JOIN cancer_counts c1 
            ON p.Cancer_Type_A = c1.Cancer_Type
        JOIN cancer_counts c2 
            ON p.Cancer_Type_B = c2.Cancer_Type
        ORDER BY p.Cancer_Type_A, p.Cancer_Type_B;
    `, tableName)

    // Execute the query
    rows, err := db.Query(query)
    if err != nil {
        log.Printf("Error executing query: %v", err)
        c.JSON(http.StatusInternalServerError, gin.H{"error": "Failed to execute query"})
        return
    }
    defer rows.Close()

    // Define a struct to hold the results
    type JaccardResult struct {
        CancerTypeA  string  `json:"cancer_type_a"`
        CancerTypeB  string  `json:"cancer_type_b"`
        Intersection int     `json:"intersection_count"`
        Union        int     `json:"union_count"`
        JaccardIndex float64 `json:"jaccard_index"`
    }

    // Collect the results
    results := []JaccardResult{}
    for rows.Next() {
        var res JaccardResult
        if err := rows.Scan(
            &res.CancerTypeA,
            &res.CancerTypeB,
            &res.Intersection,
            &res.Union,
            &res.JaccardIndex,
        ); err != nil {
            log.Printf("Error scanning row: %v", err)
            c.JSON(http.StatusInternalServerError, gin.H{"error": "Failed to parse query results"})
            return
        }
        results = append(results, res)
    }

    // Check for errors from iterating over rows
    if err := rows.Err(); err != nil {
        log.Printf("Row iteration error: %v", err)
        c.JSON(http.StatusInternalServerError, gin.H{"error": "Error processing query results"})
        return
    }

    // Return the results as JSON
    c.JSON(http.StatusOK, gin.H{"jaccard_indices": results})
}


// ------------------------------------------------------------------
// getJaccardIndexOrgansHandler
// ------------------------------------------------------------------
//
// Computes Jaccard Index for each pair of Organs based on 
// shared nullomers.
//
func getJaccardIndexOrgansHandler(c *gin.Context) {
    // Retrieve and validate the 'K' parameter
    K := c.Query("K")
    if K == "" {
        c.JSON(http.StatusBadRequest, gin.H{"error": "Missing parameter 'K'"})
        return
    }
    // Validate that K is a positive integer
    if _, err := strconv.Atoi(K); err != nil || K == "0" {
        c.JSON(http.StatusBadRequest, gin.H{"error": "Parameter 'K' must be a positive integer"})
        return
    }

    // Construct the table name safely
    tableName := fmt.Sprintf("neomers_%s", K)

    // Define the SQL query to compute Jaccard indices for all pairs,
    // (though for Jaccard across organs
    // the main references remain nullomers + cancer_type_details).
    //
    // The extra join ensures design consistency, but typically doesn't
    // affect the Jaccard logic for Organs. It's just "the same pattern."
    query := fmt.Sprintf(`
        WITH joined_data AS (
            SELECT n.nullomers_created, c.Organ
            FROM %s n
            JOIN cancer_type_details c USING (Project_Code)
        ),
        organ_counts AS (
            SELECT Organ, COUNT(DISTINCT nullomers_created) AS count
            FROM joined_data
            GROUP BY Organ
        ),
        all_organs AS (
            SELECT DISTINCT Organ
            FROM organ_counts
        ),
        pairs AS (
            SELECT a.Organ AS Organ_A, b.Organ AS Organ_B
            FROM all_organs a
            CROSS JOIN all_organs b
        ),
        intersections AS (
            SELECT 
                jd1.Organ AS Organ_A, 
                jd2.Organ AS Organ_B, 
                COUNT(DISTINCT jd1.nullomers_created) AS intersection_count
            FROM joined_data jd1
            JOIN joined_data jd2 ON jd1.nullomers_created = jd2.nullomers_created
            GROUP BY jd1.Organ, jd2.Organ
        )
        SELECT 
            p.Organ_A, 
            p.Organ_B, 
            COALESCE(i.intersection_count, 0) AS intersection_count,
            (c1.count + c2.count - COALESCE(i.intersection_count, 0)) AS union_count,
            CASE 
                WHEN p.Organ_A = p.Organ_B THEN 1.0
                WHEN (c1.count + c2.count - COALESCE(i.intersection_count, 0)) = 0 THEN 0.0
                ELSE ROUND(
                    CAST(COALESCE(i.intersection_count, 0) AS DOUBLE) 
                    / (c1.count + c2.count - COALESCE(i.intersection_count, 0)), 
                    4
                )
            END AS jaccard_index
        FROM pairs p
        LEFT JOIN intersections i 
            ON p.Organ_A = i.Organ_A 
            AND p.Organ_B = i.Organ_B
        JOIN organ_counts c1 
            ON p.Organ_A = c1.Organ
        JOIN organ_counts c2 
            ON p.Organ_B = c2.Organ
        ORDER BY p.Organ_A, p.Organ_B;
    `, tableName)

    // Execute the query
    rows, err := db.Query(query)
    if err != nil {
        log.Printf("Error executing query: %v", err)
        c.JSON(http.StatusInternalServerError, gin.H{"error": "Failed to execute query"})
        return
    }
    defer rows.Close()

    // Define a struct to hold the results
    type JaccardResult struct {
        OrganA  string  `json:"organ_a"`
        OrganB  string  `json:"organ_b"`
        Intersection int     `json:"intersection_count"`
        Union        int     `json:"union_count"`
        JaccardIndex float64 `json:"jaccard_index"`
    }

    // Collect the results
    results := []JaccardResult{}
    for rows.Next() {
        var res JaccardResult
        if err := rows.Scan(
            &res.OrganA,
            &res.OrganB,
            &res.Intersection,
            &res.Union,
            &res.JaccardIndex,
        ); err != nil {
            log.Printf("Error scanning row: %v", err)
            c.JSON(http.StatusInternalServerError, gin.H{"error": "Failed to parse query results"})
            return
        }
        results = append(results, res)
    }

    // Check for errors from iterating over rows
    if err := rows.Err(); err != nil {
        log.Printf("Row iteration error: %v", err)
        c.JSON(http.StatusInternalServerError, gin.H{"error": "Error processing query results"})
        return
    }

    // Return the results as JSON
    c.JSON(http.StatusOK, gin.H{"jaccard_indices": results})
}


func getDatasetStatsCancerTypesVaryingKHandler(c *gin.Context){
    lower_k := 11
    upper_k := 16
    all_results := map[string]interface{}{}

    for i := lower_k; i <= upper_k; i++ {
        tableName := fmt.Sprintf("neomers_%d", i)
         // Build the query. Now we join with cancer_type_details and donor_data 
        // in the same pattern:
        baseQuery := fmt.Sprintf(`
        SELECT 
            c.Cancer_Type ,
            COUNT(nullomers_created) AS count_neomers
        FROM %s n
        JOIN cancer_type_details c USING (Project_Code)
        GROUP BY Cancer_Type
        `, tableName)


        rows, err := db.Query(baseQuery)
        if err != nil {
            c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
            return
        }
        defer rows.Close()

        result := []map[string]interface{}{}
        for rows.Next() {
            var cancer_type string
            var count string
            if err := rows.Scan(&cancer_type, &count); err == nil {
                result = append(result, map[string]interface{}{
                    "cancer_type": cancer_type,
                    "count":  count,
                })
            }
        }
        all_results[strconv.Itoa(i)] = result

    }
    c.JSON(http.StatusOK, gin.H{"stats": all_results})
  
}

func getDistNeomerKCancerTypes(c *gin.Context) {
    K := c.Param("K")
    if _, err := strconv.Atoi(K); err != nil {
        c.JSON(http.StatusBadRequest, gin.H{"error": "Parameter K must be an integer"})
        return
    }

    tableName := fmt.Sprintf("distribution_neomer_%s_per_cancer", K)
    query := fmt.Sprintf(`
        SELECT DISTINCT Cancer_Type
        FROM %s
        ORDER BY Cancer_Type
    `, tableName)

    rows, err := db.Query(query)
    if err != nil {
        // Handle cases where the table for K might not exist
        c.JSON(http.StatusInternalServerError, gin.H{"error": fmt.Sprintf("Query failed, possibly no data for K=%s: %v", K, err)})
        return
    }
    defer rows.Close()

    var types []string
    for rows.Next() {
        var t string
        if err := rows.Scan(&t); err != nil {
            c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
            return
        }
        types = append(types, t)
    }
    c.JSON(http.StatusOK, gin.H{"cancerTypes": types})
}

// GET /distribution_neomer/:K/organs
// Returns all distinct Organ from the distribution table for a given K.
func getDistNeomerKOrgans(c *gin.Context) {
    K := c.Param("K")
    if _, err := strconv.Atoi(K); err != nil {
        c.JSON(http.StatusBadRequest, gin.H{"error": "Parameter K must be an integer"})
        return
    }

    tableName := fmt.Sprintf("distribution_neomer_%s_per_organ", K)
    query := fmt.Sprintf(`
        SELECT DISTINCT Organ
        FROM %s
        ORDER BY Organ;
    `, tableName)

    rows, err := db.Query(query)
    if err != nil {
        c.JSON(http.StatusInternalServerError, gin.H{"error": fmt.Sprintf("Query failed, possibly no data for K=%s: %v", K, err)})
        return
    }
    defer rows.Close()

    var organs []string
    for rows.Next() {
        var o string
        if err := rows.Scan(&o); err != nil {
            c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
            return
        }
        organs = append(organs, o)
    }
    c.JSON(http.StatusOK, gin.H{"organs": organs})
}

// GET /distribution_neomer/:K/data_by_cancer_type?cancerType=...
// Returns [{ donor_count:int, num_nullomers:int }, …] for a specific cancer type and K.
func getDistNeomerKDataByCancerType(c *gin.Context) {
    K := c.Param("K")
    if _, err := strconv.Atoi(K); err != nil {
        c.JSON(http.StatusBadRequest, gin.H{"error": "Parameter K must be an integer"})
        return
    }
    ct := c.Query("cancerType")
    if ct == "" {
        c.JSON(http.StatusBadRequest, gin.H{"error": "Missing cancerType query parameter"})
        return
    }

    tableName := fmt.Sprintf("distribution_neomer_%s_per_cancer", K)
    stmt := fmt.Sprintf(`
      SELECT donor_count, num_nullomers
      FROM %s
      WHERE Cancer_Type = ?
      ORDER BY donor_count
    `, tableName)

    rows, err := db.Query(stmt, ct)
    if err != nil {
        c.JSON(http.StatusInternalServerError, gin.H{"error": fmt.Sprintf("Query failed, possibly no data for K=%s: %v", K, err)})
        return
    }
    defer rows.Close()

    type bucket struct {
        DonorCount   int64 `json:"donorCount"`
        NumNullomers int64 `json:"numNullomers"`
    }
    var data []bucket
    for rows.Next() {
        var b bucket
        if err := rows.Scan(&b.DonorCount, &b.NumNullomers); err != nil {
            c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
            return
        }
        data = append(data, b)
    }
    c.JSON(http.StatusOK, gin.H{
        "K":            K,
        "cancerType":   ct,
        "distribution": data,
    })
}

// GET /distribution_neomer/:K/data_by_organ?organ=...
// Returns [{ donor_count:int, num_nullomers:int }, …] for a specific organ and K.
func getDistNeomerKDataByOrgan(c *gin.Context) {
    K := c.Param("K")
    if _, err := strconv.Atoi(K); err != nil {
        c.JSON(http.StatusBadRequest, gin.H{"error": "Parameter K must be an integer"})
        return
    }
    organ := c.Query("organ")
    if organ == "" {
        c.JSON(http.StatusBadRequest, gin.H{"error": "Missing organ query parameter"})
        return
    }

    tableName := fmt.Sprintf("distribution_neomer_%s_per_organ", K)
    stmt := fmt.Sprintf(`
      SELECT donor_count, num_nullomers
      FROM %s
      WHERE Organ = ?
      ORDER BY donor_count
    `, tableName)

    rows, err := db.Query(stmt, organ)
    if err != nil {
        c.JSON(http.StatusInternalServerError, gin.H{"error": fmt.Sprintf("Query failed, possibly no data for K=%s: %v", K, err)})
        return
    }
    defer rows.Close()

    type bucket struct {
        DonorCount   int64 `json:"donorCount"`
        NumNullomers int64 `json:"numNullomers"`
    }
    var data []bucket
    for rows.Next() {
        var b bucket
        if err := rows.Scan(&b.DonorCount, &b.NumNullomers); err != nil {
            c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
            return
        }
        data = append(data, b)
    }
    c.JSON(http.StatusOK, gin.H{
        "K":            K,
        "organ":        organ,
        "distribution": data,
    })
}

// GET /exome_patient_details?donor_id=…
func getExomePatientDetailsHandler(c *gin.Context) {
    actualID := c.Query("donor_id")
    if actualID == "" {
        c.JSON(http.StatusBadRequest, gin.H{"error": "Missing donor_id"})
        return
    }

    // 1) Grab the column names dynamically
    cols, err := db.Query(`SELECT * FROM exome_donor_data LIMIT 0`)
    if err != nil {
        c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
        return
    }
    colNames, _ := cols.Columns()
    cols.Close()

    // 2) Fetch the one matching row (bcr_patient_barcode = actualID)
    row := db.QueryRow(`
        SELECT * 
        FROM exome_donor_data 
        WHERE bcr_patient_barcode = ? 
        LIMIT 1
    `, actualID)

    // 3) Scan into a slice of interfaces
    vals := make([]interface{}, len(colNames))
    ptrs := make([]interface{}, len(colNames))
    for i := range vals {
        ptrs[i] = &vals[i]
    }
    if err := row.Scan(ptrs...); err != nil {
        // no match
        c.JSON(http.StatusOK, gin.H{"patient": nil})
        return
    }

    // 4) Build a map[string]interface{} to JSON-ify
    patient := make(map[string]interface{}, len(colNames))
    for i, name := range colNames {
        v := vals[i]
        if b, ok := v.([]byte); ok {
            patient[name] = string(b)
        } else {
            patient[name] = v
        }
    }

    c.JSON(http.StatusOK, gin.H{"patient": patient})
}



// GET /exome_patient_neomers?donor_id=…&length=…&top_n=…[&prefix=…]
func getExomePatientNeomersHandler(c *gin.Context) {
    donorID := c.Query("donor_id")
    lengthStr := c.Query("length")
    topNStr := c.Query("top_n")
    prefix := c.Query("prefix")

    if donorID == "" {
        c.JSON(http.StatusBadRequest, gin.H{"error": "Missing donor_id"})
        return
    }
    if lengthStr == "" {
        c.JSON(http.StatusBadRequest, gin.H{"error": "Missing length"})
        return
    }
    length, err := strconv.Atoi(lengthStr)
    if err != nil || length < 11 || length > 20 {
        c.JSON(http.StatusBadRequest, gin.H{"error": "Invalid length. Must be between 11 and 20."})
        return
    }

    topN := 10
    if topNStr != "" {
        if n, err := strconv.Atoi(topNStr); err == nil && n > 0 {
            topN = n
        }
    }

    tableName := fmt.Sprintf("exome_neomers_%d", length)

    // build the base query
    baseQuery := fmt.Sprintf(`
        SELECT
          n.nullomers_created AS neomer,
          COUNT(*) AS count
        FROM %s n
        JOIN exomes_donor_id_mapping m
          ON n.Donor_ID = m.Donor_ID
        WHERE m.Actual_Donor_ID = ?
    `, tableName)

    args := []interface{}{donorID}
    if prefix != "" {
        baseQuery += " AND n.nullomers_created LIKE ?"
        args = append(args, prefix+"%")
    }
    baseQuery += `
        GROUP BY neomer
        ORDER BY count DESC
        LIMIT ?
    `
    args = append(args, topN)

    rows, err := db.Query(baseQuery, args...)
    if err != nil {
        c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
        return
    }
    defer rows.Close()

    var out []map[string]interface{}
    for rows.Next() {
        var nm string
        var ct int
        if err := rows.Scan(&nm, &ct); err == nil {
            out = append(out, map[string]interface{}{
                "neomer": nm,
                "count":  ct,
            })
        }
    }

    c.JSON(http.StatusOK, gin.H{"neomers": out})
}


func getExomeAnalyzeNeomerHandler(c *gin.Context) {
	neomer := c.Query("neomer")
	if neomer == "" {
		c.JSON(http.StatusBadRequest, gin.H{"error": "Missing neomer"})
		return
	}
	length := len([]rune(neomer))
	
	tableName := fmt.Sprintf("exome_neomers_%d", length)

	// 1) Cancer Type breakdown
	cancerQ := fmt.Sprintf(`
		SELECT d.Cancer_Type    AS cancerType,
			   COUNT(DISTINCT di.Actual_Donor_ID) AS count -- Count distinct actual donors for cancer type
		FROM %s AS n
		JOIN exomes_donor_id_mapping AS di
			ON n.Donor_ID = di.Donor_ID
		JOIN exome_donor_data AS d
			ON di.Actual_Donor_ID = d.bcr_patient_barcode
		WHERE n.nullomers_created = ?
		GROUP BY d.Cancer_Type
		ORDER BY count DESC
	`, tableName)

	cancerRows, err := db.Query(cancerQ, neomer)
	if err != nil {
		c.JSON(http.StatusInternalServerError, gin.H{"error": "Error fetching cancer breakdown: " + err.Error()})
		return
	}
	defer cancerRows.Close()

	var cancerBreakdown []map[string]interface{}
	for cancerRows.Next() {
		var ct string
		var cnt int
		if err := cancerRows.Scan(&ct, &cnt); err == nil {
			cancerBreakdown = append(cancerBreakdown, map[string]interface{}{
				"cancerType": ct,
				"count":      cnt,
			})
		} else {
			c.JSON(http.StatusInternalServerError, gin.H{"error": "Error scanning cancer breakdown row: " + err.Error()})
			return
		}
	}
	if err := cancerRows.Err(); err != nil {
		c.JSON(http.StatusInternalServerError, gin.H{"error": "Error iterating cancer breakdown: " + err.Error()})
		return
	}

	// 2) Organ breakdown
	organQ := fmt.Sprintf(`
		SELECT d.Organ         AS organ,
			   COUNT(DISTINCT di.Actual_Donor_ID) AS count -- Count distinct actual donors for organ
		FROM %s AS n
		JOIN exomes_donor_id_mapping AS di
			ON n.Donor_ID = di.Donor_ID
		JOIN exome_donor_data AS d
			ON di.Actual_Donor_ID = d.bcr_patient_barcode
		WHERE n.nullomers_created = ?
		GROUP BY d.Organ
		ORDER BY count DESC
	`, tableName)

	organRows, err := db.Query(organQ, neomer)
	if err != nil {
		c.JSON(http.StatusInternalServerError, gin.H{"error": "Error fetching organ breakdown: " + err.Error()})
		return
	}
	defer organRows.Close()

	var organBreakdown []map[string]interface{}
	for organRows.Next() {
		var org string
		var cnt int
		if err := organRows.Scan(&org, &cnt); err == nil {
			organBreakdown = append(organBreakdown, map[string]interface{}{
				"organ": org,
				"count": cnt,
			})
		} else {
			c.JSON(http.StatusInternalServerError, gin.H{"error": "Error scanning organ breakdown row: " + err.Error()})
			return
		}
	}
	if err := organRows.Err(); err != nil {
		c.JSON(http.StatusInternalServerError, gin.H{"error": "Error iterating organ breakdown: " + err.Error()})
		return
	}

	// 3) Distinct Actual Donor IDs List
	distinctDonorIDsQueryExome := fmt.Sprintf(`
		SELECT DISTINCT di.Actual_Donor_ID
		FROM %s AS n
		JOIN exomes_donor_id_mapping AS di ON n.Donor_ID = di.Donor_ID
		WHERE n.nullomers_created = ?
	`, tableName)

	exomeDonorRows, err := db.Query(distinctDonorIDsQueryExome, neomer)
	if err != nil {
		c.JSON(http.StatusInternalServerError, gin.H{"error": "Error fetching distinct exome donor IDs: " + err.Error()})
		return
	}
	defer exomeDonorRows.Close()

	var distinctExomeDonorIDsList []string
	for exomeDonorRows.Next() {
		var actualDonorID string
		if err := exomeDonorRows.Scan(&actualDonorID); err != nil {
			c.JSON(http.StatusInternalServerError, gin.H{"error": "Error scanning distinct exome donor ID row: " + err.Error()})
			return
		}
		distinctExomeDonorIDsList = append(distinctExomeDonorIDsList, actualDonorID)
	}
	if err := exomeDonorRows.Err(); err != nil {
		c.JSON(http.StatusInternalServerError, gin.H{"error": "Error iterating distinct exome donor IDs: " + err.Error()})
		return
	}

	c.JSON(http.StatusOK, gin.H{
		"analysis": gin.H{
			"cancerBreakdown":  cancerBreakdown,
			"organBreakdown":   organBreakdown,
			"distinctDonorIDs": distinctExomeDonorIDsList, // Added this list
		},
	})
}


  