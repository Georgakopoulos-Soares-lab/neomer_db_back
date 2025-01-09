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

    router.GET("/patient_details", getPatientDetailsHandler)
    router.GET("/patient_neomers", getPatientNeomersHandler)
    router.GET("/analyze_neomer", analyzeNeomerHandler)

    router.GET("/jaccard_index", getJaccardIndexHandler)



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
    return "/storage/group/izg5139/default/external/neo_database/neomers.ddb"
}

// ------------------------------------------------------------------
// getNullomersHandler
// ------------------------------------------------------------------
//
// Returns *all* columns from nullomers_%s plus all columns from cancer_type_details,
// plus an added column "gc_content" (rounded 2 decimals, times -1).
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
    filters := c.Query("filters")             // e.g. "(gc_content > 10) AND (gc_content < 50)"
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

    // 1) Base CTE returning *all* columns from both tables, plus computed gc_content
    // We'll alias the tables to avoid collisions, e.g. n.*, c.*.
    baseQuery := fmt.Sprintf(`
        WITH base AS (
            SELECT
                n.*,
                c.*,
                ROUND(
                    100.0 * (
                        LENGTH(n.nullomers_created)
                        - LENGTH(REPLACE(UPPER(n.nullomers_created), 'G', ''))
                        - LENGTH(REPLACE(UPPER(n.nullomers_created), 'C', ''))
                    ) / LENGTH(n.nullomers_created),
                    2
                ) * -1 AS gc_content
            FROM nullomers_%[1]s n
            JOIN cancer_type_details c USING (Project_Code)
        )
        SELECT * FROM base
    `, length)

    // Build WHERE
    whereClauses := []string{}
    if filters != "" {
        whereClauses = append(whereClauses, filters)
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
                        // We'll do an IN referencing original table
                        subQuery := fmt.Sprintf(`
                            nullomers_created IN (
                                SELECT nullomers_created
                                FROM nullomers_%[1]s
                                JOIN cancer_type_details USING (Project_Code)
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

    // 2) COUNT query with same CTE
    countQuery := fmt.Sprintf(`
        WITH base AS (
            SELECT
                n.*,
                c.*,
                ROUND(
                    100.0 * (
                        LENGTH(n.nullomers_created)
                        - LENGTH(REPLACE(UPPER(n.nullomers_created), 'G', ''))
                        - LENGTH(REPLACE(UPPER(n.nullomers_created), 'C', ''))
                    ) / LENGTH(n.nullomers_created),
                    2
                ) * -1 AS gc_content
            FROM nullomers_%[1]s n
            JOIN cancer_type_details c USING (Project_Code)
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
    // we typically skip. But do as you wish:
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
            // equals, notEquals, starts, notStarts...
            cond = fmt.Sprintf("WHERE LOWER(%s) LIKE '%s%%'", column, lowerInput)
        }
    }

    query := fmt.Sprintf(`
        WITH base AS (
            SELECT DISTINCT %s
            FROM nullomers_%s
            JOIN cancer_type_details USING (Project_Code)
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

    // We'll keep all columns from both tables, plus negative GC content
    baseCTE := fmt.Sprintf(`
        WITH base AS (
            SELECT
                n.*,
                c.*,
                ROUND(
                    100.0 * (
                        LENGTH(n.nullomers_created)
                        - LENGTH(REPLACE(UPPER(n.nullomers_created), 'G', ''))
                        - LENGTH(REPLACE(UPPER(n.nullomers_created), 'C', ''))
                    ) / LENGTH(n.nullomers_created),
                    2
                ) * -1 AS gc_content
            FROM nullomers_%[1]s n
            JOIN cancer_type_details c USING (Project_Code)
        )
    `, length)

    // Build WHERE
    whereClauses := []string{}
    if filters != "" {
        whereClauses = append(whereClauses, filters)
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
                        // We'll do an IN referencing original table
                        subQuery := fmt.Sprintf(`
                            nullomers_created IN (
                                SELECT nullomers_created
                                FROM nullomers_%[1]s
                                JOIN cancer_type_details USING (Project_Code)
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

    // If you want to group by gc_content, you can do:
    // groupByCols = append(groupByCols, "gc_content")
    // selectCols = append(selectCols, "gc_content")

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

    // Example query (customize as needed)
    // We assume "Project_Code" is in "donor_data" so we can join with "cancer_type_details"
    query := `
        SELECT d.* , c.Cancer_Type , c.Organ
        FROM donor_data d, cancer_type_details c
        WHERE d.icgc_donor_id = ? AND POSITION(c.Acronym IN d.project_code) > 0
        LIMIT 1
    `
    row := db.QueryRow(query, donorID)

    // Grab columns you want or do "SELECT ..." instead of "*".
    // For brevity, we show a simplified approach:
    columns, err := db.Query("SELECT * FROM donor_data LIMIT 0") // to get columns
    if err != nil {
        c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
        return
    }
    colNames, _ := columns.Columns()
    columns.Close()

    colNames = append(colNames, "Cancer_Type", "Cancer_Organ")

    vals := make([]interface{}, len(colNames))
    valPtrs := make([]interface{}, len(colNames))
    for i := range vals {
        valPtrs[i] = &vals[i]
    }


    if err := row.Scan(valPtrs...); err != nil {
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

// getPatientNeomersHandler handles the /patient_neomers endpoint
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
    tableName := fmt.Sprintf("nullomers_%d", length)

    // Validate table name to prevent SQL injection (ensure length is between 11-20)
    if length < 11 || length > 20 {
        c.JSON(http.StatusBadRequest, gin.H{"error": "Invalid length. Must be between 11 and 20."})
        return
    }

    // Build the base query with mandatory donor_id condition
    baseQuery := fmt.Sprintf(`
        SELECT nullomers_created AS neomer, COUNT(*) AS count
        FROM %s
        WHERE donor_id = ?
    `, tableName)

    // If prefix is provided, add a LIKE condition
    if prefix != "" {
        baseQuery += " AND nullomers_created LIKE ?"
    }

    // Add GROUP BY, ORDER BY, and LIMIT clauses
    baseQuery += `
        GROUP BY neomer
        ORDER BY count DESC
        LIMIT ?
    `

    var rows *sql.Rows
    if prefix != "" {
        likePattern := prefix + "%"
        rows, err = db.Query(baseQuery, donorID, likePattern, topN)
    } else {
        rows, err = db.Query(baseQuery, donorID, topN)
    }

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
// Parameters:
// - neomer: string (required)
//
// Returns:
// {
//   "analysis": {
//     "totalNeomers": 100,
//     "distinctDonors": 20,
//     "distinctCancerTypes": 5,
//     "distinctOrgans": 3
//   }
// }

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

    tableName := fmt.Sprintf("nullomers_%d", neomerLength)

    // Validate table existence or sanitize the table name as needed
    // For simplicity, assuming table exists and name is safe since length is validated

    query := fmt.Sprintf(`
        SELECT
            COUNT(*) AS total_count,
            COUNT(DISTINCT donor_id) AS distinct_donors,
            COUNT(DISTINCT Cancer_Type) AS distinct_cancer_types,
            COUNT(DISTINCT Organ) AS distinct_organs
        FROM %s
        JOIN cancer_type_details USING (Project_Code)
        WHERE nullomers_created = ?
    `, tableName)

    row := db.QueryRow(query, neomer)

    var totalCount, distinctDonors, distinctCancerTypes, distinctOrgans int
    err = row.Scan(&totalCount, &distinctDonors, &distinctCancerTypes, &distinctOrgans)
    if err != nil {
        c.JSON(http.StatusInternalServerError, gin.H{"error": err.Error()})
        return
    }

    analysis := map[string]interface{}{
        "totalNeomers":        totalCount,
        "distinctDonors":      distinctDonors,
        "distinctCancerTypes": distinctCancerTypes,
        "distinctOrgans":      distinctOrgans,
    }

    c.JSON(http.StatusOK, gin.H{"analysis": analysis})
}


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
    tableName := fmt.Sprintf("nullomers_%s", K)

    // Open the database connection
    dbPath := getDatabasePath()
    db, err := sql.Open("duckdb", dbPath)
    if err != nil {
        log.Printf("Error opening database: %v", err)
        c.JSON(http.StatusInternalServerError, gin.H{"error": "Database connection failed"})
        return
    }
    defer db.Close()

   

    // Define the SQL query to compute Jaccard indices for all possible pairs
    query := fmt.Sprintf(`
        WITH joined_data AS (
            SELECT nt.nullomers_created, c.Cancer_Type
            FROM %s nt
            JOIN cancer_type_details c ON nt.Project_Code = c.Project_Code
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
                ELSE ROUND(CAST(COALESCE(i.intersection_count, 0) AS DOUBLE) / (c1.count + c2.count - COALESCE(i.intersection_count, 0)), 4)
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
        if err := rows.Scan(&res.CancerTypeA, &res.CancerTypeB, &res.Intersection, &res.Union, &res.JaccardIndex); err != nil {
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