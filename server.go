package main

import (
    "database/sql"
    "fmt"
    "log"
    "net/http"
    "os"
    "strconv"
    "strings"

    "github.com/gin-gonic/gin"
    _ "github.com/marcboeker/go-duckdb"
)

func main() {
    router := gin.Default()

    // Disable CORS policy
    router.Use(corsMiddleware())

    router.GET("/healthcheck", healthCheckHandler)
    router.GET("/cancer_types", makeHandler("SELECT * FROM cancer_types"))
    router.GET("/donor_data", makeHandler("SELECT * FROM donor_data"))
    router.GET("/tcga_survival_data", makeHandler("SELECT * FROM tcga_survival_data"))

    router.GET("/get_nullomers", getNullomersHandler)
    router.GET("/get_suggestions", getSuggestionsHandler)
    router.GET("/get_nullomers_stats", getNullomersStatsHandler)

    router.GET("/get_exomes", getExomesHandler)
    router.GET("/get_exomes_stats", getExomesStatsHandler)

    router.GET("/patient_details", getPatientDetailsHandler)
    router.GET("/patient_neomers", getPatientNeomersHandler)
    router.GET("/analyze_neomer", analyzeNeomerHandler)

    router.GET("/jaccard_index", getJaccardIndexHandler)
    router.GET("/jaccard_index_organs", getJaccardIndexOrgansHandler)

    router.GET("/dataset_stats_cancer_types_varying_k", getDatasetStatsCancerTypesVaryingKHandler)



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
        dbPath := getDatabasePath()
        db, err := sql.Open("duckdb", dbPath)
        if err != nil {
            c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
            return
        }
        defer db.Close()

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
// getNullomersHandler
// ------------------------------------------------------------------
//
// Returns *all* columns from neomers_%s plus all columns from 
// cancer_type_details plus donor_data, plus an added column "gc_content."
//
func getNullomersHandler(c *gin.Context) {
    length := c.Query("length")

    if length == "" {
        c.JSON(http.StatusBadRequest, gin.H{"error": "Missing required parameter 'length'"})
        return
    }

    // Pagination
    pageStr := c.Query("page")
    limitStr := c.Query("limit")
    filters := c.Query("filters")         // e.g. "(gc_content > 10) AND (gc_content < 50)"
    specialFilters := c.Query("specialFilters") // e.g. "at_least_X_distinct_patients;3"

    page := 0
    limit := 10000
    if pageStr != "" {
        if p, err := strconv.Atoi(pageStr); err == nil && p >= 0 {
            page = p
        }
    }
    if limitStr != "" {
        if l, err := strconv.Atoi(limitStr); err == nil && l > 0 && l <= 10000 {
            limit = l
        }
    }


    dbPath := getDatabasePath()

    db, err := sql.Open("duckdb", dbPath)
    if err != nil {
        c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
        return
    }
    defer db.Close()

    // 1) Base CTE returning *all* columns from the three tables, plus computed gc_content
    baseQuery := fmt.Sprintf(`
        WITH base AS (
            SELECT
                n.* EXCLUDE (Donor_ID),
                c.*,
                d.*,
                di.Tumor_Sample_Barcode, di.Matched_Norm_Sample_Barcode,
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
            LEFT JOIN donor_data d ON di.Actual_Donor_ID = d.icgc_donor_id        )
        SELECT * FROM base
    `, length)

         // Build WHERE Clause
         whereClauses := []string{}
         if filters != "" {
             filterConditions := strings.Split(filters, " AND ") // Split individual filter conditions
             for _, condition := range filterConditions {
                 parts := strings.Fields(condition) // Split by space
                 if len(parts) >= 3 {
                     column := cleanColumnName(parts[0]) // Ensure column name is cleaned properly
                     if isNumericColumn(column) {
                         condition = fmt.Sprintf(`CAST("%s" AS FLOAT) %s %s`, column, parts[1], removeParentheses(parts[2]))
                     }
         
                     whereClauses = append(whereClauses, condition)
                 }
             }
         }

    // Special filters (like at_least_X_distinct_patients)
    if specialFilters != "" {
        parts := strings.Split(specialFilters, "|")
        for _, part := range parts {
            sfPieces := strings.Split(part, ";")
            switch sfPieces[0] {
            case "at_least_X_distinct_patients":
                if len(sfPieces) == 2 {
                    distinctStr := sfPieces[1]
                    distinctCount, err := strconv.Atoi(distinctStr)
                    if err == nil && distinctCount > 0 {
                        subQuery := fmt.Sprintf(`
                            nullomers_created IN (
                                SELECT nullomers_created
                                FROM neomers_%[1]s
                                JOIN cancer_type_details USING (Project_Code)
                                LEFT JOIN donor_id_mapping di ON CAST(n."Donor_ID" AS INT) = di."Donor_ID"
                                LEFT JOIN donor_data d ON di.Actual_Donor_ID = d.icgc_donor_id)   
                                GROUP BY nullomers_created
                                HAVING COUNT(DISTINCT donor_id) >= %d
                            )
                        `, length, distinctCount)
                        whereClauses = append(whereClauses, subQuery)
                    }
                }
            }
        }
    }

    finalWhere := ""
    if len(whereClauses) > 0 {
        finalWhere = " WHERE " + strings.Join(whereClauses, " AND ")
    }    
    fmt.Println("🔍  finalWhere", finalWhere)


    // 2) COUNT query with same CTE
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
            LEFT JOIN donor_data d ON di.Actual_Donor_ID = d.icgc_donor_id        )        
        SELECT COUNT(*) FROM base
        %s
    `, length, finalWhere)

    var totalCount int
    err = db.QueryRow(countQuery).Scan(&totalCount)
    if err != nil {
        c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
        return
    }

    // 3) Final query with LIMIT/OFFSET
    offset := page * limit
    query := fmt.Sprintf("%s %s LIMIT %d OFFSET %d", baseQuery, finalWhere, limit, offset)

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

        // Convert []byte -> string, handle NULL
        for i, val := range row {
            if b, ok := val.([]byte); ok {
                row[i] = string(b)
            }
            if val == nil {
                row[i] = nil
            }
        }
        data = append(data, row)
    }

    result := map[string]interface{}{
        "headers":    columns,
        "data":       data,
        "totalCount": totalCount,
    }
    c.JSON(http.StatusOK, result)
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
func getSuggestionsHandler(c *gin.Context) {
    column := c.Query("column")
    input := c.Query("input")
    filterType := c.Query("filterType")
    length := c.Query("length")
    if length == "" {
        c.JSON(http.StatusBadRequest, gin.H{"error": "Missing length"})
        return
    }

    // If user tries to get suggestions for numeric columns like gc_content,
    // we can skip. Or return nothing, as below:
    if column == "gc_content" {
        c.JSON(http.StatusOK, gin.H{"suggestions": []string{}})
        return
    }

    dbPath := getDatabasePath()
    db, err := sql.Open("duckdb", dbPath)
    if err != nil {
        c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
        return
    }
    defer db.Close()

    lowerInput := strings.ToLower(input)
    var cond string
    if lowerInput == "" {
        cond = ""
    } else {
        switch filterType {
        case "contains":
            cond = fmt.Sprintf("WHERE LOWER(%s) LIKE '%%%s%%'", column, lowerInput)
        case "ends":
            cond = fmt.Sprintf("WHERE LOWER(%s) LIKE '%%%s'", column, lowerInput)
        default:
            // e.g. "starts", etc.
            cond = fmt.Sprintf("WHERE LOWER(%s) LIKE '%s%%'", column, lowerInput)
        }
    }

    // Join with donor_data using the same pattern:
    query := fmt.Sprintf(`
        WITH base AS (
            SELECT DISTINCT %s
            FROM neomers_%s n
            JOIN cancer_type_details c USING (Project_Code)
            LEFT JOIN donor_id_mapping di ON CAST(n."Donor_ID" AS INT) = di."Donor_ID"
            LEFT JOIN donor_data d ON di.Actual_Donor_ID = d.icgc_donor_id)   
            %s
            LIMIT 10
        )
        SELECT %s
        FROM base
        ORDER BY LOWER(%s) ASC
    `, column, length, cond, column, column)

    rows, err := db.Query(query)
    if err != nil {
        // Not necessarily an error, could be no suggestions
        c.JSON(http.StatusOK, gin.H{"suggestions": []string{}})
        return
    }
    defer rows.Close()

    suggestions := []string{}
    for rows.Next() {
        var val interface{}
        if err := rows.Scan(&val); err == nil {
            if v, ok := val.(string); ok {
                suggestions = append(suggestions, v)
            } else if val != nil {
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
    groupByStr := c.Query("groupBy")
    topNStr := c.Query("topN")
    specialFilters := c.Query("specialFilters")
    log.Printf("Received length: %s", length)


    topN := 10
    if topNStr != "" {
        if val, err := strconv.Atoi(topNStr); err == nil && val > 0 {
            topN = val
        }
    }

    dbPath := getDatabasePath()
    db, err := sql.Open("duckdb", dbPath)
    if err != nil {
        c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
        return
    }
    defer db.Close()

    // We'll keep all columns from all three tables, plus negative GC content
    baseCTE := fmt.Sprintf(`
        WITH base AS (
            SELECT
                n.* EXCLUDE (Donor_ID),
                c.*,
                d.*,
                di.Tumor_Sample_Barcode, di.Matched_Norm_Sample_Barcode,
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
            LEFT JOIN donor_data d ON di.Actual_Donor_ID = d.icgc_donor_id)           
    `, length)

      // Build WHERE Clause
         whereClauses := []string{}
         if filters != "" {
             filterConditions := strings.Split(filters, " AND ") // Split individual filter conditions
             for _, condition := range filterConditions {
                 parts := strings.Fields(condition) // Split by space
                 if len(parts) >= 3 {
                     column := cleanColumnName(parts[0]) // Ensure column name is cleaned properly
                     if isNumericColumn(column) {
                         condition = fmt.Sprintf(`CAST("%s" AS FLOAT) %s %s`, column, parts[1], removeParentheses(parts[2]))
                     }
         
                     whereClauses = append(whereClauses, condition)
                 }
             }
         }

    if specialFilters != "" {
        parts := strings.Split(specialFilters, "|")
        for _, part := range parts {
            sfPieces := strings.Split(part, ";")
            switch sfPieces[0] {
            case "at_least_X_distinct_patients":
                if len(sfPieces) == 2 {
                    distinctStr := sfPieces[1]
                    distinctCount, err := strconv.Atoi(distinctStr)
                    if err == nil && distinctCount > 0 {
                        subQuery := fmt.Sprintf(`
                            nullomers_created IN (
                                SELECT nullomers_created
                                FROM neomers_%[1]s
                                JOIN cancer_type_details USING (Project_Code)
                                LEFT JOIN exomes_donor_id_mapping di ON CAST(n."Donor_ID" AS INT) = di."Donor_ID"
                                LEFT JOIN donor_data d ON di.Actual_Donor_ID = d.icgc_donor_id
                                GROUP BY nullomers_created
                                HAVING COUNT(DISTINCT donor_id) >= %d
                            )
                        `, length, distinctCount)
                        whereClauses = append(whereClauses, subQuery)
                    }
                }
            }
        }
    }

    finalWhere := ""
    if len(whereClauses) > 0 {
        finalWhere = " WHERE " + strings.Join(whereClauses, " AND ")
    }

    // Figure out groupBy
    groupByCols := []string{"nullomers_created"}
    selectCols := []string{"nullomers_created"}

    if groupByStr != "" {
        additional := strings.Split(groupByStr, ",")
        for _, col := range additional {
            col = strings.TrimSpace(col)
            if col != "" {
                groupByCols = append(groupByCols, col)
                selectCols = append(selectCols, col)
            }
        }
    }

    selectClause := strings.Join(selectCols, ", ")
    groupByClause := strings.Join(groupByCols, ", ")

    // Build final stats query
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

        for i, val := range row {
            if b, ok := val.([]byte); ok {
                row[i] = string(b)
            }
            if val == nil {
                row[i] = nil
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

// ------------------------------------------------------------------
// getExomesHandler
// ------------------------------------------------------------------


func getExomesHandler(c *gin.Context) {
    length := c.Query("length")

    if length == "" {
        c.JSON(http.StatusBadRequest, gin.H{"error": "Missing required parameter 'length'"})
        return
    }

    // Pagination
    pageStr := c.Query("page")
    limitStr := c.Query("limit")
    filters := c.Query("filters")         // e.g. "(gc_content > 10) AND (gc_content < 50)"
    specialFilters := c.Query("specialFilters") // e.g. "at_least_X_distinct_patients;3"
    fmt.Println("🔍 Filters:", filters, "| Special Filters:", specialFilters)

    page := 0
    limit := 10000
    if pageStr != "" {
        if p, err := strconv.Atoi(pageStr); err == nil && p >= 0 {
            page = p
        }
    }
    if limitStr != "" {
        if l, err := strconv.Atoi(limitStr); err == nil && l > 0 && l <= 10000 {
            limit = l
        }
    }


    dbPath := getDatabasePath()

    db, err := sql.Open("duckdb", dbPath)
    if err != nil {
        c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
        return
    }
    defer db.Close()

    // 1) Base CTE returning *all* columns from the three tables, plus computed gc_content
    baseQuery := fmt.Sprintf(`
        WITH base AS (
            SELECT
            n.* EXCLUDE (Donor_ID), d.*, di.Tumor_Sample_Barcode, di.Matched_Norm_Sample_Barcode,
                ROUND(
                    100.0 * (
                        LENGTH(n.nullomers_created)
                        - LENGTH(REPLACE(UPPER(n.nullomers_created), 'G', ''))
                        - LENGTH(REPLACE(UPPER(n.nullomers_created), 'C', ''))
                    ) / LENGTH(n.nullomers_created),
                    2
                ) * -1 AS gc_content
            FROM exome_neomers_%[1]s n
            LEFT JOIN exomes_donor_id_mapping di ON CAST(n."Donor_ID" AS INT) = di."Donor_ID"
            LEFT JOIN donor_data d ON di.Actual_Donor_ID = d.icgc_donor_id
        )
        SELECT * FROM base
    `, length)

    // Build WHERE Clause
         whereClauses := []string{}
         if filters != "" {
             filterConditions := strings.Split(filters, " AND ") // Split individual filter conditions
             for _, condition := range filterConditions {
                 parts := strings.Fields(condition) // Split by space
                 if len(parts) >= 3 {
                     column := cleanColumnName(parts[0]) // Ensure column name is cleaned properly
                     if isNumericColumn(column) {
                         condition = fmt.Sprintf(`CAST("%s" AS FLOAT) %s %s`, column, parts[1], removeParentheses(parts[2]))
                     }
         
                     whereClauses = append(whereClauses, condition)
                 }
             }
         }
    
    // Special filters (like at_least_X_distinct_patients)
    if specialFilters != "" {
        parts := strings.Split(specialFilters, "|")
        for _, part := range parts {
            sfPieces := strings.Split(part, ";")
            switch sfPieces[0] {
            case "at_least_X_distinct_patients":
                if len(sfPieces) == 2 {
                    distinctStr := sfPieces[1]
                    distinctCount, err := strconv.Atoi(distinctStr)
                    if err == nil && distinctCount > 0 {
                        subQuery := fmt.Sprintf(`
                            nullomers_created IN (
                                SELECT nullomers_created
                                FROM exome_neomers_%[1]s n
                                LEFT JOIN exomes_donor_id_mapping di ON CAST(n."Donor_ID" AS INT) = di."Donor_ID"
                                LEFT JOIN donor_data d ON di.Actual_Donor_ID = d.icgc_donor_id
                                GROUP BY nullomers_created
                                HAVING COUNT(DISTINCT di.Actual_Donor_ID) >= %d
                            )
                        `, length, distinctCount)
                        whereClauses = append(whereClauses, subQuery)
                    }
                }
            }
        }
    }

    finalWhere := ""
    if len(whereClauses) > 0 {
        finalWhere = " WHERE " + strings.Join(whereClauses, " AND ")
    }

    // 2) COUNT query with same CTE
    countQuery := fmt.Sprintf(`
      WITH base AS (
            SELECT
                n.*, 
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
            FROM exome_neomers_%[1]s n
            LEFT JOIN exomes_donor_id_mapping di 
                ON CAST(n."Donor_ID" AS INT) = di."Donor_ID"
            LEFT JOIN donor_data d 
                ON di.Actual_Donor_ID = d.icgc_donor_id
        )
        SELECT COUNT(*) FROM base
        %s
    `, length, finalWhere)

    var totalCount int
    err = db.QueryRow(countQuery).Scan(&totalCount)
    if err != nil {
        c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
        return
    }

    // 3) Final query with LIMIT/OFFSET
    offset := page * limit
    query := fmt.Sprintf("%s %s LIMIT %d OFFSET %d", baseQuery, finalWhere, limit, offset)
    fmt.Println("Executing SQL Query:", query)


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

        // Convert []byte -> string, handle NULL
        for i, val := range row {
            if b, ok := val.([]byte); ok {
                row[i] = string(b)
            }
            if val == nil {
                row[i] = nil
            }
        }
        data = append(data, row)
    }

    result := map[string]interface{}{
        "headers":    columns,
        "data":       data,
        "totalCount": totalCount,
    }
    c.JSON(http.StatusOK, result)
}

// ------------------------------------------------------------------
// removeParentheses helper function
// ------------------------------------------------------------------
func removeParentheses(input string) string {
    re := regexp.MustCompile(`[()]`) // Matches ( and )
    return re.ReplaceAllString(input, "")
}

// ------------------------------------------------------------------
// getExomesStatsHandler
// ------------------------------------------------------------------

func getExomesStatsHandler(c *gin.Context) {
    length := c.Query("length")

    if length == "" {
        c.JSON(http.StatusBadRequest, gin.H{"error": "Missing required parameter 'length'"})
        return
    }

    filters := c.Query("filters")
    groupByStr := c.Query("groupBy")
    topNStr := c.Query("topN")
    specialFilters := c.Query("specialFilters")
    log.Printf("Received length: %s", length)


    topN := 10
    if topNStr != "" {
        if val, err := strconv.Atoi(topNStr); err == nil && val > 0 {
            topN = val
        }
    }

    dbPath := getDatabasePath()
    db, err := sql.Open("duckdb", dbPath)
    if err != nil {
        c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
        return
    }
    defer db.Close()

    // We'll keep all columns from all three tables, plus negative GC content
    baseCTE := fmt.Sprintf(`
       WITH base AS (
            SELECT
                n.* EXCLUDE (Donor_ID),
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
            FROM exome_neomers_%[1]s n
            LEFT JOIN exomes_donor_id_mapping di 
                ON CAST(n."Donor_ID" AS INT) = di."Donor_ID"
            LEFT JOIN donor_data d 
                ON di.Actual_Donor_ID = d.icgc_donor_id
        )
    `, length)

      // Build WHERE Clause
      whereClauses := []string{}
      if filters != "" {
          filterConditions := strings.Split(filters, " AND ") // Split individual filter conditions
          for _, condition := range filterConditions {
              parts := strings.Fields(condition) // Split by space
              if len(parts) >= 3 {
                  column := cleanColumnName(parts[0]) // Ensure column name is cleaned properly
                  if isNumericColumn(column) {
                      condition = fmt.Sprintf(`CAST("%s" AS FLOAT) %s %s`, column, parts[1], removeParentheses(parts[2]))
                  }
      
                  whereClauses = append(whereClauses, condition)
              }
          }
      }

    if specialFilters != "" {
        parts := strings.Split(specialFilters, "|")
        for _, part := range parts {
            sfPieces := strings.Split(part, ";")
            switch sfPieces[0] {
            case "at_least_X_distinct_patients":
                if len(sfPieces) == 2 {
                    distinctStr := sfPieces[1]
                    distinctCount, err := strconv.Atoi(distinctStr)
                    if err == nil && distinctCount > 0 {
                        subQuery := fmt.Sprintf(`
                             nullomers_created IN (
                                SELECT nullomers_created
                                FROM exome_neomers_%[1]s n
                                LEFT JOIN exomes_donor_id_mapping di ON CAST(n."Donor_ID" AS INT) = di."Donor_ID"
                                LEFT JOIN donor_data d ON di.Actual_Donor_ID = d.icgc_donor_id
                                GROUP BY nullomers_created
                                HAVING COUNT(DISTINCT di.Actual_Donor_ID) >= %d
                            )
                        `, length, distinctCount)
                        whereClauses = append(whereClauses, subQuery)
                    }
                }
            }
        }
    }

    finalWhere := ""
    if len(whereClauses) > 0 {
        finalWhere = " WHERE " + strings.Join(whereClauses, " AND ")
    }

    // Figure out groupBy
    groupByCols := []string{"nullomers_created"}
    selectCols := []string{"nullomers_created"}

    if groupByStr != "" {
        additional := strings.Split(groupByStr, ",")
        for _, col := range additional {
            col = strings.TrimSpace(col)
            if col != "" {
                groupByCols = append(groupByCols, col)
                selectCols = append(selectCols, col)
            }
        }
    }

    selectClause := strings.Join(selectCols, ", ")
    groupByClause := strings.Join(groupByCols, ", ")

    // Build final stats query
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

        for i, val := range row {
            if b, ok := val.([]byte); ok {
                row[i] = string(b)
            }
            if val == nil {
                row[i] = nil
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

// ------------------------------------------------------------------
// getPatientDetailsHandler
// ------------------------------------------------------------------
func getPatientDetailsHandler(c *gin.Context) {
    donorID := c.Query("donor_id")
    if donorID == "" {
        c.JSON(http.StatusBadRequest, gin.H{"error": "Missing donor_id"})
        return
    }

    dbPath := getDatabasePath()
    db, err := sql.Open("duckdb", dbPath)
    if err != nil {
        c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
        return
    }
    defer db.Close()

    // Example query (customize as needed).
    // We assume "Project_Code" is in "donor_data" so we can join with "cancer_type_details"
    query := `
        SELECT d.*, c.Cancer_Type, c.Organ
        FROM donor_data d, cancer_type_details c
        WHERE d.icgc_donor_id = ? 
          AND POSITION(c.Acronym IN d.project_code) > 0
        LIMIT 1
    `
    row := db.QueryRow(query, donorID)

    // Grab columns you want or do "SELECT ..." instead of "*".
    columns, err := db.Query("SELECT * FROM donor_data LIMIT 0") // just to get column names
    if err != nil {
        c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
        return
    }
    colNames, _ := columns.Columns()
    columns.Close()

    // Append the extra columns from cancer_type_details
    colNames = append(colNames, "Cancer_Type", "Organ")

    vals := make([]interface{}, len(colNames))
    valPtrs := make([]interface{}, len(colNames))
    for i := range vals {
        valPtrs[i] = &vals[i]
    }

    if err := row.Scan(valPtrs...); err != nil {
        // If not found, return "patient": nil instead of an error
        c.JSON(http.StatusOK, gin.H{"patient": nil})
        fmt.Println(err)
        return
    }

    // Convert to map
    patientMap := make(map[string]interface{})
    for i, colName := range colNames {
        patientMap[colName] = vals[i]
        if b, ok := vals[i].([]byte); ok {
            patientMap[colName] = string(b)
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

    topN := 10 // default
    if topNStr != "" {
        topNParsed, err := strconv.Atoi(topNStr)
        if err == nil && topNParsed > 0 {
            topN = topNParsed
        }
    }

    dbPath := getDatabasePath()
    db, err := sql.Open("duckdb", dbPath)
    if err != nil {
        c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
        return
    }
    defer db.Close()

    // Dynamically inject the length variable into the table name
    tableName := fmt.Sprintf("neomers_%d", length)

    // Build the query. Now we join with cancer_type_details and donor_data 
    // in the same pattern:
    baseQuery := fmt.Sprintf(`
        SELECT 
            n.nullomers_created AS neomer, 
            COUNT(*) AS count
        FROM %s n
        JOIN cancer_type_details c USING (Project_Code)
        JOIN donor_data d ON n.donor_id = d.icgc_donor_id
        WHERE n.donor_id = ?
    `, tableName)

    var args []interface{}
    args = append(args, donorID)

    // If prefix is provided, add a LIKE condition
    if prefix != "" {
        baseQuery += " AND n.nullomers_created LIKE ?"
        likePattern := prefix + "%"
        args = append(args, likePattern)
    }

    // Add GROUP BY, ORDER BY, and LIMIT clauses
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

    result := []map[string]interface{}{}
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

    neomerLength := len(neomer)
    if neomerLength < 11 || neomerLength > 20 {
        c.JSON(http.StatusBadRequest, gin.H{"error": "Neomer length must be between 11 and 20"})
        return
    }

    dbPath := getDatabasePath()
    db, err := sql.Open("duckdb", dbPath)
    if err != nil {
        c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
        return
    }
    defer db.Close()

    tableName := fmt.Sprintf("neomers_%d", neomerLength)

    // First Query: Basic Statistics
    totalQuery := fmt.Sprintf(`
        SELECT
            COUNT(*) AS total_count,
            COUNT(DISTINCT n.donor_id) AS distinct_donors,
            COUNT(DISTINCT c.Cancer_Type) AS distinct_cancer_types,
            COUNT(DISTINCT c.Organ) AS distinct_organs
        FROM %s n
        JOIN cancer_type_details c USING (Project_Code)
        JOIN donor_data d ON n.donor_id = d.icgc_donor_id
        WHERE n.nullomers_created = ?
    `, tableName)

    row := db.QueryRow(totalQuery, neomer)

    var totalCount, distinctDonors, distinctCancerTypes, distinctOrgans int
    err = row.Scan(&totalCount, &distinctDonors, &distinctCancerTypes, &distinctOrgans)
    if err != nil {
        c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
        return
    }

    // Second Query: Breakdown by Cancer_Type and Organ
    breakdownQuery := fmt.Sprintf(`
        SELECT
            c.Cancer_Type,
            c.Organ,
            COUNT(*) AS count
        FROM %s n
        JOIN cancer_type_details c USING (Project_Code)
        JOIN donor_data d ON n.donor_id = d.icgc_donor_id
        WHERE n.nullomers_created = ?
        GROUP BY c.Cancer_Type, c.Organ
    `, tableName)

    rows, err := db.Query(breakdownQuery, neomer)
    if err != nil {
        c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
        return
    }
    defer rows.Close()

    // Structures to hold the breakdown data
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
        var cancerType, organ string
        var count int
        err := rows.Scan(&cancerType, &organ, &count)
        if err != nil {
            c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
            return
        }

        if _, exists := cancerMap[cancerType]; !exists {
            cancerMap[cancerType] = &CancerTypeCount{
                CancerType: cancerType,
                Count:      0,
                Organs:     []OrganCount{},
            }
        }

        cancerMap[cancerType].Count += count
        cancerMap[cancerType].Organs = append(cancerMap[cancerType].Organs, OrganCount{
            Organ: organ,
            Count: count,
        })
    }

    // Convert map to slice for JSON serialization
    cancerTypes := make([]CancerTypeCount, 0, len(cancerMap))
    for _, ct := range cancerMap {
        cancerTypes = append(cancerTypes, *ct)
    }

    // Construct the final analysis response
    analysis := map[string]interface{}{
        "totalNeomers":        totalCount,
        "distinctDonors":      distinctDonors,
        "distinctCancerTypes": distinctCancerTypes,
        "distinctOrgans":      distinctOrgans,
        "cancerBreakdown":     cancerTypes,
    }

    c.JSON(http.StatusOK, gin.H{"analysis": analysis})
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

    // Open the database connection
    dbPath := getDatabasePath()
    db, err := sql.Open("duckdb", dbPath)
    if err != nil {
        log.Printf("Error opening database: %v", err)
        c.JSON(http.StatusInternalServerError, gin.H{"error": "Database connection failed"})
        return
    }
    defer db.Close()

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

    // Open the database connection
    dbPath := getDatabasePath()
    db, err := sql.Open("duckdb", dbPath)
    if err != nil {
        log.Printf("Error opening database: %v", err)
        c.JSON(http.StatusInternalServerError, gin.H{"error": "Database connection failed"})
        return
    }
    defer db.Close()

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
    
    dbPath := getDatabasePath()
    db, err := sql.Open("duckdb", dbPath)
    if err != nil {
        c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
        return
    }
    defer db.Close()

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