-- This query retrieves all records from the summary table for operons that contain exactly one CRISPR array 
-- and are classified as subtype VI, which is associated with Cas13b.

CREATE OR REPLACE VIEW cas13b_vi AS
SELECT s.*, c.evalue, c.gene_name
FROM cas c
JOIN summary s 
  ON c.operon_id = s.operon_id
WHERE s.subtype = 'VI'
AND c.gene_name = 'Cas13b';



-- Retrieve spacers associated with Cas13b VI operons
