# Radar_labs

DECLARE
    table_name VARCHAR2(100);
BEGIN
    FOR table_rec IN (
        SELECT table_name
        FROM all_tables
        WHERE owner = 'SLFR'
        AND table_name LIKE 'FAIR_OUT_%'
        AND TO_DATE(SUBSTR(table_name, -5), 'MONYY') < TO_DATE('DEC23', 'MONYY')
    )
    LOOP
        EXECUTE IMMEDIATE 'DROP TABLE ' || table_rec.table_name;
        DBMS_OUTPUT.PUT_LINE('Dropped table: ' || table_rec.table_name);
    END LOOP;
END;
/
