import sqlite3
from sqlite3 import Error
import pandas as pd

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

def create_qc_result(conn, result):
    """
    Create a new result
    :param conn:
    :param result:
    :return:
    """

    sql = ''' REPLACE INTO qc_results(sample_id, run_id, kraken2_top_species, MASH_ID, status, QCtag, cg_coverage, quast_contigs, quast_length, quast_gc, amr_genes, assembly_file, analysis_date, armadillo_version)
              VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?) '''
    cur = conn.cursor()
    cur.execute(sql, result)
    conn.commit()

    return cur.lastrowid

def update_metadata(database, xlsx):
    conn = create_connection(database)
    db = sqlite3.connect(database)
    dfs = pd.read_excel(xlsx, sheet_name=None, engine='openpyxl')
    for table, df in dfs.items():
        df.to_sql("metadata", db, if_exists="replace")
    # cur = conn.cursor()
    # cur.execute("SELECT * FROM metadata")
    # rows = cur.fetchall()
    # for row in rows:
    #     print(row)
    conn.close()

if __name__ == '__main__':
    update_metadata("/home/dnalab/database/arln.sqlite", "/home/dnalab/database/ARLN_metadata.xlsx")