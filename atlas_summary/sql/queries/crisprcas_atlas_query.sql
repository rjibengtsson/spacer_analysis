SELECT gene_name, COUNT(*) AS count
FROM cas
GROUP BY gene_name
ORDER BY count DESC;


SELECT DISTINCT fp.phageid
FROM feature f
JOIN featureid_phageid fp ON f.featureid = fp.featureid
WHERE (
    f.product ILIKE '%terminase%'
    OR f.product ~* '\mterl\M'
);





SELECT fp.phageid, f.featureid, f.start, f.end, f.strand, f.product 
FROM feature f
JOIN featureid_phageid fp ON f.featureid = fp.featureid
WHERE fp.phageid = 'Han_2018_ERR1398082_NODE_13_length_96793_cov_101.246098'
AND (
f.product ILIKE '%terminase%'
OR f.product ~* '\mterl\M'
);




SELECT COUNT(DISTINCT c.operon_id) AS operon_count
FROM cas c
JOIN summary s
  ON c.operon_id = s.operon_id
WHERE c.evalue = 0
  AND c.gene_name IN ('Cas13f', 'Cas13i')
  AND s.subtype = 'VI';


SELECT *
FROM summary
WHERE subtype = 'VI';


SELECT operon_id, gene_name, hmm_name, evalue, score, truncated, length
FROM cas
WHERE gene_name ILIKE 'cas12%';


SELECT s.operon_id, s.spacer_length, s.spacer_sequence
FROM spacer s
JOIN cas c
  ON s.operon_id = c.operon_id
WHERE c.gene_name ILIKE 'cas13b%'
  AND c.evalue = 0
  AND (s.spacer_length >= 25 OR s.spacer_length <= 35);

  


-- -- DROP TABLE spacer;
-- -- DROP TABLE cas;

-- -- Index summary table

-- -- ALTER TABLE public.cas
-- -- ADD CONSTRAINT uq_cas_operon UNIQUE (operon_id, gene_name, hmm_name, evalue, score);