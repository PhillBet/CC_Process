'''
    ARGUMENTS
'''
# standard
import argparse

# third party


# local


def get_arguments():
    '''
    Defines the arguments that the program will support.

    Returns
        arguments, argparse Object, defined arguments for the execution of the program.
    '''

    # ARGUMENTS DESCRIPTION

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="",
        epilog="You need to provided at least one entity argument")

    # GENERAL ARGUMENTS

    parser.add_argument(
        "-sum",
        "--summary",
        help="Path to read the GUs summary.",
        metavar="InputFiles/GensorUnits_Summary.json",
        default="InputFiles/GensorUnits_Summary.json",
    )

    parser.add_argument(
        "-grp",
        "--groups",
        help="Path to read the GUs groups.",
        metavar="InputFiles/GensorUnits_Groups.json",
        default="InputFiles/GensorUnits_Groups.json",
    )

    parser.add_argument(
        "-in",
        "--input",
        help="Path to read the origin file data.",
        metavar="../InputData/",
        default="../InputData/",
    )

    parser.add_argument(
        "-out",
        "--output",
        help="Path where the json files of the process will be stored.",
        metavar="../RawData/",
        default="../RawData/",
    )

    parser.add_argument(
        "-l",
        "--log",
        help="Path where the log of the process will be stored.",
        metavar="../logs/ht_etl_log/",
        default="../logs/ht_etl_log/",
    )

    parser.add_argument(
        "-r",
        "--report",
        help="Path where the report of the process will be stored.",
        metavar="../logs/ht_etl_log/",
        default="../logs/ht_etl_log/",
    )

    parser.add_argument(
        "-org",
        "--organism",
        help="Organism whose information is been downloaded.",
        default="ECOLI",
        metavar="ecoli",
    )

    parser.add_argument(
        "-v",
        "--version",
        help="Imput Data Verison.",
        default="0.0.1",
        metavar="0.0.1",
    )

    parser.add_argument(
        "-dstype",
        "--dataset-type",
        help="Dataset record source name.",
        choices=["TFBINDING", "GENE_EXPRESSION",
                 "TSS", "TUS", "TTS", "REGULONS", "GSELEX"]
    )

    parser.add_argument(
        "-email",
        "--email",
        help="User email address to connect to PUBMED database.",
        default="reguadm@ccg.unam.mx",
        metavar="reguadm@ccg.unam.mx",
    )

    parser.add_argument(
        "-u",
        "--url",
        help="URL to DB server.",
        default="mongodb://localhost",
        metavar="mongodb://localhost",
    )

    parser.add_argument(
        "-db",
        "--database",
        help="Name of the database where IDs are taken.",
        choices=["regulondbmultigenomic",
                 "regulondbht", "regulondbdatamarts"],
        default="regulondbmultigenomic",
        metavar="regulondbmultigenomic",
    )

    parser.add_argument(
        "-iddb",
        "--id-database",
        help="Name of the database where IDs are taken.",
        default="regulondbidentifiers",
        metavar="regulondbidentifiers",
    )

    parser.add_argument(
        "-tabname",
        "--table-name",
        help="table for testings.",
        default='all',
        metavar='all',
    )

    arguments = parser.parse_args()

    return arguments


def load_arguments():
    '''
    Load the arguments that the program will support.

    Returns
        arguments, argparse Object, loaded arguments for the execution of the program.
    '''

    arguments = get_arguments()
    return arguments
