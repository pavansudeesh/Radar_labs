WITH filtered_tables AS (
    SELECT
        TABLE_NAME,
        REGEXP_REPLACE(TABLE_NAME, '_[[:alpha:]]{3}\d{2}$', '') AS TBLNAME,
        ROW_NUMBER() OVER (PARTITION BY REGEXP_REPLACE(TABLE_NAME, '_[[:alpha:]]{3}\d{2}$', '') ORDER BY TABLE_NAME) AS rn
    FROM
        all_tables
    WHERE
        owner = 'SLFR'
        AND table_name LIKE 'FAIR\_%\_OUT%'
        AND table_name NOT LIKE 'FAIR_OUT_CONTROL%'
)
SELECT TABLE_NAME, TBLNAME
FROM filtered_tables
WHERE rn = 1;
