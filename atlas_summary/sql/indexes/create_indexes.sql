-- Index gene_name in cas table
CREATE INDEX idx_cas_gene_name 
ON cas(gene_name varchar_pattern_ops);




-- All indexes in the database
SELECT tablename, indexname, indexdef
FROM pg_indexes
WHERE schemaname = 'public'
ORDER BY tablename;

-- For a specific table
SELECT tablename, indexname, indexdef
FROM pg_indexes
WHERE tablename = 'cas';