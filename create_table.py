import sqlite3
from sqlite3 import Error
def create_connection(db_file):
    """ create a database connection to the SQLite database
        specified by db_file
    :param db_file: database file
    :return: Connection object or None
    """
    conn = None
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except Error as e:
        print(e)

    return conn


def create_table(conn, create_table_sql):
    """ create a table from the create_table_sql statement
    :param conn: Connection object
    :param create_table_sql: a CREATE TABLE statement
    :return:
    """
    try:
        c = conn.cursor()
        c.execute(create_table_sql)
    except Error as e:
        print(e)


def main():
    database = r"database/pythonsqlite.db"

    sql_create_qc_results_table = """ CREATE TABLE IF NOT EXISTS qc_results (
                                        sample_id text PRIMARY KEY,
                                        run_id text NOT NULL,
                                        kraken2_top_species text NOT NULL,
                                        MASH_ID text,
                                        status text NOT NULL,
                                        QCtag text NOT NULL,
                                        cg_coverage double NOT NULL,
                                        quast_contigs integer NOT NULL,
                                        quast_length integer NOT NULL,
                                        quast_gc double,
                                        amr_genes text,
                                        assembly_file text,
                                        analysis_date date,
                                        armadillo_version text
                                    ); """

    sql_create_metadata_table = """CREATE TABLE IF NOT EXISTS metadata (
                                    sample_id text PRIMARY KEY,
                                    run_name text NOT NULL,
                                    workbook_id text NOT NULL,
                                    key text NOT NULL,
                                    source_site text,
                                    submitter text,
                                    source_city text,
                                    patient_age_years integer,
                                    patient_age_months integer,
                                    patient_age_days integer,
                                    patient_sex text,
                                    isolate_date date,
                                    received_date date,
                                    lab_id text,
                                    source_country text,
                                    source_state text,
                                    source_county text
                                );"""

    # create a database connection
    conn = create_connection(database)

    # create tables
    if conn is not None:
        # create projects table
        create_table(conn, sql_create_qc_results_table)

        # create tasks table
        create_table(conn, sql_create_metadata_table)
    else:
        print("Error! cannot create the database connection.")


if __name__ == '__main__':
    main()