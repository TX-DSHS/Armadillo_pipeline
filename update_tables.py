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

def main():
    database = r"database/pythonsqlite.db"

    # create a database connection
    conn = create_connection(database)

    with conn:
     
        result = ('1TXAMD2100313','AR_220729_M05358','Proteus mirabilis','Proteus mirabilis','Complete','', 91.11, 56, 4371101, 39.5, 'Id,blaCMY-2,blaNDM-5,blaTEM-1,ble,catA,catA1,catB11,dfrA1,dfrA12,dfrA17,erm(42),merC,merP,merR,merT,qacE,qnrA1,sat2,sul1,sul2,terD,terZ,tet(J)', 's3://804609861260-bioinformatics-infectious-disease/cluster/Proteus_mirabilis/TXAMD2100313-TX-M05358-220729_contigs.fa','','')
        create_qc_result(conn, result)

if __name__ == '__main__':
    main()