CREATE INDEX cas_operon_id_idx ON cas (operon_id);
-- CREATE INDEX subtype_idx ON summary (subtype);
CREATE INDEX idx_cas_gene_evalue_operon ON cas (gene_name, evalue, operon_id);
CREATE INDEX idx_summary_subtype_operon ON summary (subtype, operon_id);






-------------------- Identifying & Managing Rarely Used Indexes ----------

-- Check when stats were last reset for the current database
SELECT
    datname,
    stats_reset
FROM pg_stat_database
WHERE datname = current_database();


-- Identify for rarely used indexes
SELECT
    schemaname,
    relname        AS table_name,
    indexrelname   AS index_name,
    idx_scan       AS times_used,
    pg_size_pretty(pg_relation_size(indexrelid)) AS index_size
FROM pg_stat_user_indexes
WHERE idx_scan < 10 -- Adjust threshold as needed
ORDER BY pg_relation_size(indexrelid) DESC;


-- Extend query to include index definition & uniqueness
SELECT
    s.schemaname,
    s.relname                          AS table_name,
    s.indexrelname                     AS index_name,
    s.idx_scan                         AS times_scanned,
    pg_size_pretty(pg_relation_size(s.indexrelid)) AS index_size,
    i.indisunique                      AS is_unique,
    i.indisprimary                     AS is_primary,
    pg_get_indexdef(i.indexrelid)      AS index_definition
FROM pg_stat_user_indexes s
JOIN pg_index i ON s.indexrelid = i.indexrelid
WHERE s.idx_scan < 10                  -- Adjust threshold as needed
  AND i.indisprimary = false           -- Exclude primary keys
  AND i.indisunique  = false           -- Exclude unique constraints
ORDER BY pg_relation_size(s.indexrelid) DESC;




-- Find Duplicate or Redundant Indexes
SELECT
    a.indrelid::regclass   AS table_name,
    a.indexrelid::regclass AS index_a,
    b.indexrelid::regclass AS index_b,
    pg_get_indexdef(a.indexrelid) AS def_a,
    pg_get_indexdef(b.indexrelid) AS def_b
FROM pg_index a
JOIN pg_index b
    ON a.indrelid  = b.indrelid
    AND a.indexrelid < b.indexrelid
    AND (a.indkey::int[] @> b.indkey::int[]
     OR  b.indkey::int[] @> a.indkey::int[])
-- Add this join and filter to exclude system catalogues
JOIN pg_class c ON c.oid = a.indrelid
JOIN pg_namespace n ON n.oid = c.relnamespace
WHERE n.nspname NOT IN ('pg_catalog', 'pg_toast', 'information_schema');



-- Once the query flags a pair, follow up with this to see actual usage:
SELECT
    s.indexrelname   AS index_name,
    s.idx_scan       AS times_used,
    pg_get_indexdef(s.indexrelid) AS definition
FROM pg_stat_user_indexes s
WHERE s.indexrelname IN ('cas_operon_id_idx', 'idx_cas_gene_evalue_operon');








-- Check the indexes
-- Lists all indexes defined in the database
SELECT schemaname, tablename, indexname, indexdef
FROM pg_indexes
WHERE schemaname = 'public'
ORDER BY tablename, indexname;

-- Monitor index usage
-- Track which indexes are being used
SELECT
    schemaname,
    relname        AS table_name,
    indexrelname   AS index_name,
    idx_scan       AS times_used,
    idx_tup_read   AS tuples_read,
    idx_tup_fetch  AS tuples_fetched
FROM pg_stat_user_indexes
ORDER BY idx_scan DESC;







-- Check If the Index Backs a Constraint
SELECT
    conname     AS constraint_name,
    contype     AS constraint_type,
    conrelid::regclass AS table_name
FROM pg_constraint
WHERE conindid = (
    SELECT oid FROM pg_class WHERE relname = 'your_index_name'
);

