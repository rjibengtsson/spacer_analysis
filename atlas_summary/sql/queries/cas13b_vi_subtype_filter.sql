-- Query: Cas13b subtype VI operons with single CRISPR array
-- Author: Rebecca
-- Date: 2026-05
-- Filters: subtype VI, n_crispr = 1, evalue = 0

-- Get operon length distribution for Cas13b
SELECT 
    COUNT(*)                AS n,
    MIN(operon_length)      AS min,
    MAX(operon_length)      AS max,
    AVG(operon_length)      AS mean,
    STDDEV(operon_length)   AS std_dev,
    PERCENTILE_CONT(0.25) WITHIN GROUP (ORDER BY operon_length) AS q1,
    PERCENTILE_CONT(0.50) WITHIN GROUP (ORDER BY operon_length) AS median,
    PERCENTILE_CONT(0.75) WITHIN GROUP (ORDER BY operon_length) AS q3
FROM cas13b_vi
WHERE n_crispr = 1
AND evalue = 0;



-- Get spacers number distribution for Cas13b
SELECT 
    COUNT(*)                AS n,
    MIN(n_spacers)      AS min,
    MAX(n_spacers)      AS max,
    AVG(n_spacers)      AS mean,
    STDDEV(n_spacers)   AS std_dev,
    PERCENTILE_CONT(0.25) WITHIN GROUP (ORDER BY n_spacers) AS q1,
    PERCENTILE_CONT(0.50) WITHIN GROUP (ORDER BY n_spacers) AS median,
    PERCENTILE_CONT(0.75) WITHIN GROUP (ORDER BY n_spacers) AS q3
FROM cas13b_vi
WHERE n_crispr = 1
AND evalue = 0;


SELECT sp.spacer_length
FROM cas c
JOIN summary s 
  ON c.operon_id = s.operon_id
JOIN spacer sp
  ON sp.operon_id = s.operon_id
WHERE s.subtype = 'VI'
AND c.gene_name = 'Cas13b'
AND n_crispr = 1
AND evalue = 0;











SELECT n_crispr, COUNT(*) AS count
FROM summary
WHERE subtype = 'VI'
GROUP BY n_crispr
ORDER BY count DESC;


