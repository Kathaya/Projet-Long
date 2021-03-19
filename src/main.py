import os, argparse, sys;
from sys import exit

def get_arguments():
    """
    Lecture des arguments
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', help="Path to fasta file to analyse")
    parser.add_argument('-d', help="Path to database to search in")
    parser.add_argument('-t', type=int, default=30, help="Number of threads. Default: 30")
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    wdpath = sys.argv[0].replace("main.py", "")
    args = get_arguments()
    fastapath = args.i
    dbDir = args.d
    threads = args.t
    name  = fastapath.split("/")[-1].split(".")[0]
    if not os.path.isfile(fastapath) or not fastapath.endswith("fasta"):
        print("File doesn't exist or doesn't have fasta extensions.\nPlease correct this before running")
        exit()
   
    os.system("hhblits -i {} -d {} -cpu {} -oa3m {}".format(fastapath, dbDir, threads, fastapath.split(".")[0]+".a3m"))
    hhrpath = fastapath.split(".")[0]+".hhr"
    print(fastapath, hhrpath)
    os.system("python {}model_selection.py -p {} -f {}".format(wdpath, hhrpath, fastapath))
    """
    os.environ['threads'] = str(threads)
    os.environ['dataDir'] = dataDir
    os.environ['DB'] =  dbDir
    hh_time = time.time()
    #Effectue le hhblits sur les fichiers fastas
    for file in os.listdir(dataDir):
        if file.endswith(".fasta"):
            #print(file)
            os.environ['file'] = file
            os.system("hhblits -i $dataDir/$file -d $DB -cpu $threads -oa3m $dataDir/$file.a3m")
    hh_time = time.time() - hh_time
    print("hh_time={}".format(hh_time))
        
    """

        