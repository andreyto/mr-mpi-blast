////////////////////////////////////////////////////////////////////////
//
//  MPI-MR-BLAST: Parallelizing BLAST on MR-MPI
//
//  Author: Seung-Jin Sul
//          (ssul@jcvi.org)
//
//  Last updated: 02/08/2011
//
////////////////////////////////////////////////////////////////////////

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA 02110-1301, USA.
 */

/*
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 */

/// MPI AND MR-MPI
#include "mpi.h"
#include "mrmpi/mapreduce.h"
#include "mrmpi/keyvalue.h"

using namespace MAPREDUCE_NS;
using namespace std;

/// ----------------------------------------------------------------------------
/// NCBI C++ Toolkit
/// ----------------------------------------------------------------------------
#include <algo/blast/api/local_blast.hpp>
#include <algo/blast/api/objmgr_query_data.hpp>
#include <algo/blast/blastinput/blast_fasta_input.hpp>
 
/// Import search strategy
#include <objects/blast/Blast4_request.hpp>
#include <algo/blast/api/search_strategy.hpp>

USING_NCBI_SCOPE;
USING_SCOPE(blast);
/// ----------------------------------------------------------------------------


/// ----------------------------------------------------------------------------
/// Boost lib
/// ----------------------------------------------------------------------------
/// For tokenizer
#include <boost/algorithm/string.hpp>

/// For str -> int, int -> str
#include <boost/lexical_cast.hpp>

/// For processing command line arguments
#include <boost/program_options.hpp>
namespace po = boost::program_options;

/// For processing configuration file 
#include <boost/program_options/detail/config_file.hpp>
namespace pod = boost::program_options::detail;

/// For Boost memory mapped file
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/filesystem/operations.hpp>         /// for real file size
boost::iostreams::mapped_file_source MMAPFILE;     /// Read-only Boost mmap file
/// ----------------------------------------------------------------------------

/// For typedef unsigned long long int uint32_t
#include <stdint.h>

/// ----------------------------------------------------------------------------
/// Settings from mrblast.ini conf file
/// ----------------------------------------------------------------------------
int VERBOSITY;              /// log from mr-mpi lib
int TIMER;                  /// log elapsed time for each mapreduce call 
int MEMSIZE;                /// # the page size (in Mbytes)

bool EXCLUSIONORNOT;
int EXCLUSIONTHRESHOLD;     /// Exclusion threshold = 100bp
string EXCLUSIONHISTFILE;
 
/// DB options
string DBCHUNKLISTFILE;
string CONFFILE;
int NTOTALDBCHUNKS;

/// Log
#include <sys/time.h>
#include <sys/resource.h>
int LOGORNOT = 0;
int TIMING = 0;
int OPTDUMP = 0;            /// For dumping Blast opts out
string LOGFILENAME;
int CNT = 0;
/// ----------------------------------------------------------------------------

/// Log
string LOGMSG;
ofstream LOGFILE;
ostream* LOGSTREAM = NULL;
#define LOG (*LOGSTREAM)

/// For nucl or prot DB setting
bool ISPROTEIN = false;

/// Blast target DB setting
static CSearchDatabase *pTargetDb = 0;
string prevDbChunkName;

/// Import search strategy
string STRATEGYFILENAME;    /// Input blast search option file

/// Misc.
const int MAXSTR = 256;
char MPI_procName[MAXSTR];
const int MAXSTR2 = 256;    /// only for our simulated query sequence
const int SUBIDMAX = 20;    /// For subject ID of blast hits. possible KMV overflow

const int QUERY = 0;        /// To retireve query info from CSeq_align
const int SUBJECT = 1;      /// To retireve suject info from CSeq_align
int MYID;                   /// MPI rank
string OUTPREFIX;           /// Prefix string for output file names
string INDEXFILENAME;
string QUERYFILENAME;
int NITERATION;             /// iteration number
int MAPSTYLE = 0;
int NQUERYPERWORKITEM;      /// block size
#define NDEBUG 1

/// For syncronized timing
#ifndef MPI_WTIME_IS_GLOBAL
#define MPI_WTIME_IS_GLOBAL 1
#endif

/// To pass Blast hits following outfmt=6 format.
/// subject id, % identity, alignment length, mismatches, gap opens, 
/// q. start, q. end, s. start, s. end, evalue, bit score
typedef struct blastres {
    char subjectid[SUBIDMAX];
    double identity;
    uint32_t alignlen; 
    int mismatches;
    uint32_t gapopens;
    uint32_t qstart; 
    uint32_t qend;
    uint32_t sstart;
    uint32_t send;
    double evalue; 
    int bitscore;
} BLASTRES;
int BLASTRESSZ = sizeof(BLASTRES);

/// To sort Blast hits by evalue
typedef struct structEvalue {
    BLASTRES* p;
    double evalue;
} STRUCTEVALUE;
string WORKEROUTPUTFILENAME = "";
 
/// MR-MPI CALLS
void mr_run_blast(int itask, char *file, KeyValue *kv, void *ptr);                   
void mr_sort_multivalues_by_evalue(char *key, int keybytes, char *multivalue, 
                                   int nvalues, int *valuebytes, KeyValue *kv, 
                                   void *ptr);                         
inline bool mycompare(STRUCTEVALUE e1, STRUCTEVALUE e2);

/// Check hit exclusion
inline bool check_exclusion(string qGi, string sGi, int qCutLocStart, 
                            int qCutLocEnd, int sStart, int sEnd, int threshold);

 
/* -------------------------------------------------------------------------- */
int main(int argc, char **argv)
/* -------------------------------------------------------------------------- */
{     
    /// 
    /// nquery: number of queries for a work item. If there are q queries 
    /// in the input file, (q/nquery) blocks of queries are generated. In turn, 
    /// If there are p DB partitions, there are p*(q/nquery) work items.
    /// nquery=0 means there is only one block of queries and therefore p work
    /// items will be generated.
    ///
    /// NOTE: query bp size should be considered for nquery later.
    ///
    /// iteration: number of iterations. p*(q/nquery) work items are divided 
    /// into p*(q/nquery) / iteration files and each will be processed one by 
    /// one.
    ///
    /// If mapstyle = 2, to use the default master/slave mode in MR-MPI.
    /// If mapstyle = 3, to use the modified master/slave mode for controlling 
    /// DB partition locality.
    ///
    po::options_description generalDesc("General options");
    generalDesc.add_options() 
        ("help", "print help message")        
        ("query-file,i", po::value<string>(), "set input query file")
        ("index-file,d", po::value<string>(), "set input index file")
        ("import-search-strategy,s", po::value<string>(), 
            "set search strategy file")
        ("db-list,l", po::value<string>(), 
            "set DB partition name lst file")
        ("conf-file,c", po::value<string>(), 
            "set configuration file")
        /// should add dbchunks.txt and mrblast.ini
    ;
    
    po::options_description OptionalDesc("Optional");
    OptionalDesc.add_options() 
        ("block-size,b", 
            po::value<int>(&NQUERYPERWORKITEM)->default_value(0), 
            "set the number of queries per work item (default=0=all)")
        ("iteration,n", 
            po::value<int>(&NITERATION)->default_value(1), 
            "set the number of iterations")            
        ("output-prefix,o", 
            po::value<string>(&OUTPREFIX)->default_value("output"), 
            "set output prefix for output file names (default=output)")
        ("map-style,m", 
            po::value<int>(&MAPSTYLE)->default_value(2), 
            "set MR-MPI mapstyle: 2=master/slave, 3=new scheduler")
        ("is-protein,p", 
            po::value<bool>(&ISPROTEIN)->default_value(false), 
            "set the program as blastp (default=0)")
        ("self-exclusion,x", 
            po::value<bool>(&EXCLUSIONORNOT)->default_value(false), 
            "enable self exclusion")
    ;
    
    po::options_description allDesc("Allowed options");
    allDesc.add(generalDesc).add(OptionalDesc);

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, allDesc), vm);
    po::notify(vm); 
        
    if (argc < 2 || (!strcmp(argv[1], "-?") || !strcmp(argv[1], "--?") 
        || !strcmp(argv[1], "/?") || !strcmp(argv[1], "/h") 
        || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--h") 
        || !strcmp(argv[1], "--help") || !strcmp(argv[1], "/help") 
        || !strcmp(argv[1], "-help")  || !strcmp(argv[1], "help") )) {
        cout << "MR-MPI Blast\n" << "Author: Seung-Jin Sul (ssul@jcvi.org)\n\n"
             << allDesc;
        return 1;
    }
    else {
        if (vm.count("query-file"))  
            QUERYFILENAME = vm["query-file"].as<string>();
        else {
            cerr << "ERROR: query file was not set.\n\n"; 
            cout << allDesc;
            return 1;
        }
        if (vm.count("index-file"))  
            INDEXFILENAME = vm["index-file"].as<string>();
        else {
            cerr << "ERROR: index file was not set.\n\n"; 
            cout << allDesc;
            return 1;
        }
        if (vm.count("import-search-strategy")) 
            STRATEGYFILENAME = vm["import-search-strategy"].as<string>();
        else {
            cerr << "ERROR: option file was not set.\n\n"; 
            cout << allDesc;
            return 1;
        }
        if (vm.count("db-list")) 
            DBCHUNKLISTFILE = vm["db-list"].as<string>();
        else {
            cerr << "ERROR: DB name list file was not set.\n\n"; 
            cout << allDesc;
            return 1;
        }
        if (vm.count("conf-file")) 
            CONFFILE = vm["conf-file"].as<string>();
        else {
            cerr << "ERROR: configuration file was not set.\n\n"; 
            cout << allDesc;
            return 1;
        }
        
        /// 
        /// Read OPTIONAL arguments
        ///
        if (vm.count("output-prefix")) 
            OUTPREFIX = vm["output-prefix"].as<string>();
        if (vm.count("block-size")) 
            NQUERYPERWORKITEM = vm["block-size"].as<int>();
        if (vm.count("iteration")) 
            NITERATION = vm["iteration"].as<int>();
        if (vm.count("map-style")) 
            MAPSTYLE = vm["map-style"].as<int>();
        if (vm.count("self-exclusion")) 
            EXCLUSIONORNOT = vm["self-exclusion"].as<bool>();
    }
       
    ///    
    /// Read conf file, mrblast.ini and set parameters
    ///
    ifstream config(CONFFILE.c_str(), ios::in);
    if (!config) {
        cerr << "ERROR: configuration file not found" << endl;
        return 1;
    }
    
    set<string> options;
    map<string, string> parameters;
    options.insert("*");
    
    try {      
        for (pod::config_file_iterator i(config, options), e ; i != e; ++i) {
            parameters[i->string_key] = i->value[0];
        }
        //DBCHUNKLISTFILE = parameters["DBCHUNKLISTFILE"];
        EXCLUSIONHISTFILE = parameters["EXCLUSIONHISTFILE"];
        try { 
            VERBOSITY = boost::lexical_cast<int>(
                parameters["VERBOSITY"]);
            TIMER = boost::lexical_cast<int>(
                parameters["TIMER"]);
            MEMSIZE = boost::lexical_cast<int>(
                parameters["MEMSIZE"]);         
            EXCLUSIONTHRESHOLD = boost::lexical_cast<int>(
                parameters["EXCLUSIONTHRESHOLD"]);     
            LOGORNOT = boost::lexical_cast<int>(
                parameters["LOGORNOT"]);
            TIMING = boost::lexical_cast<int>(
                parameters["TIMING"]);
            OPTDUMP = boost::lexical_cast<int>(
                parameters["OPTDUMP"]);
            LOGFILENAME = boost::lexical_cast<string>(
                parameters["LOGFILENAME"]);
        }
        catch(const boost::bad_lexical_cast &) {
            cerr << "Exception: bad_lexical_cast" << endl;
        }
    }
    catch(exception& e) {
        cerr<< "Exception: " << e.what() << endl;
    }
    
    ///
    /// Creat memory-mapped file for feature vectors
    ///
    unsigned long int realFileSize = boost::filesystem::file_size(QUERYFILENAME);
    boost::iostreams::mapped_file_params params;
    params.path = QUERYFILENAME;
    params.length = realFileSize;
    /// mapped_file_source is actually read-only
    params.mode = std::ios_base::in;
    MMAPFILE.open(params);

    if (!MMAPFILE.is_open()) {
        cerr << "ERROR: failed to create mmap file\n";
        MPI_Finalize();
        exit(1);
    }
    
    ///
    /// MPI setup
    ///
    int MPI_nProcs, MPI_length;
    
    MPI_Init(&argc, &argv);    
    MPI_Comm_rank(MPI_COMM_WORLD, &MYID);
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_nProcs);
    MPI_Get_processor_name(MPI_procName, &MPI_length);
    
    /// Log file init
    if (LOGORNOT || TIMING) {
        LOGFILENAME = OUTPREFIX + "-" 
            + boost::lexical_cast<string>(MYID) + "-" + LOGFILENAME;
        LOGFILE.open(LOGFILENAME.c_str(), ios::out);
        LOGSTREAM = &LOGFILE;
    }
    LOGMSG = "Rank:" + boost::lexical_cast<string>(MYID) + " ";
    
    double profile_time;
    struct timeval totalStartTime;
    struct timeval totalEndTime;
    struct timeval totalStart_u_Time;
    struct timeval totalEnd_u_Time;
    struct timeval totalStart_s_Time;
    struct timeval totalEnd_s_Time;
    struct rusage ru_total;
        
    if (TIMING) {        
        /// MPI_Wtime 
        profile_time = MPI_Wtime();
        /// Wall-clock time
        gettimeofday(&totalStartTime, NULL);
        /// Process time
        getrusage(RUSAGE_SELF, &ru_total);
        totalStart_u_Time = ru_total.ru_utime;
        totalStart_s_Time = ru_total.ru_stime;
        LOG << "mr-mpi-blast starts," 
            << profile_time << ","
            << totalStartTime.tv_sec*1000000 + totalStartTime.tv_usec << ","
            << totalStart_u_Time.tv_sec*1000000 + totalStart_u_Time.tv_usec 
            << ","
            << totalStart_s_Time.tv_sec*1000000 + totalStart_s_Time.tv_usec 
            << endl;
        LOG.flush();
    }
    
    /// 
    /// MR-MPI init
    ///
    
    MapReduce *mr2 = new MapReduce(MPI_COMM_WORLD);
    
    /*
    * mapstyle = 0 (chunk) or 1 (stride) or 2 (master/slave)
    * all2all = 0 (irregular communication) or 1 (use MPI_Alltoallv)
    * verbosity = 0 (none) or 1 (summary) or 2 (histogrammed)
    * timer = 0 (none) or 1 (summary) or 2 (histogrammed)
    * memsize = N = number of Mbytes per page of memory
    * minpage = N = # of pages to pre-allocate per processor
    * maxpage = N = max # of pages allocatable per processor, 0 = no limit
    * keyalign = N = byte-alignment of keys
    * valuealign = N = byte-alignment of values
    * fpath = string 
    */
    mr2->verbosity = VERBOSITY;
    mr2->timer = TIMER;
    mr2->memsize = MEMSIZE;
    mr2->keyalign = sizeof(uint32_t); /// The key is a begin offset 
    mr2->mapstyle = MAPSTYLE;         /// master/slave mode=2, custom scheduler=3
    mr2->outofcore = -1;              /// disable out-of-core
    MPI_Barrier(MPI_COMM_WORLD);
        
    ///
    /// Make a file for map() which contains a list of file
    /// names. A file name is a form of "queryFile,dbChunkName"
    /// In map(), the file name is splitted into "queryFile" and
    /// "dbChunkName". We've got 109 DB chunks (Jul 2010).
    ///
    vector<string> vWorkItems;
    vector<string> vDbChunkNames;
    uint32_t nQueryBlocks = 0;
    uint32_t nWorkItems;    
    uint32_t nSubWorkItemFiles = 0;
    uint32_t nWorkItemsPerFile = 0;
    
    double master_init_time;
    if (LOGORNOT) master_init_time = MPI_Wtime();
           
    if (MYID == 0) {
        ///
        /// Read DB partition name list file
        ///
        string line;
        ifstream dbChunkNameFile(DBCHUNKLISTFILE.c_str(), ios::in);
        if (!dbChunkNameFile.is_open()) {
            cerr << "ERROR: dbchunks.txt open error.\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        while (!getline(dbChunkNameFile, line).eof()) {
            if (line.length() > 0)
                vDbChunkNames.push_back(line);
        }
        dbChunkNameFile.close();        
        NTOTALDBCHUNKS = vDbChunkNames.size();
        
        ///
        /// Read seqeunce line index file 
        /// Make work items from vDbChunkNames and vQueryFileNames
        /// a work item = <(begin_offset, end_offset), DBChunkName>
        /// *_offset mean the line index of the sequences in the 
        /// original query sequence file.
        ///
        ifstream indexFile(INDEXFILENAME.c_str(), ios::in);
        if (!indexFile.is_open()) {
            cerr << "ERROR: Index file open error.\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    
        string beginOffset;
        string endOffset;
        if (NQUERYPERWORKITEM) { /// If NQUERYPERWORKITEM > 0
            uint32_t n = 0;
            while (!getline(indexFile, line).eof()) {
                vector<string> vOffsets;
                boost::split(vOffsets, line, boost::is_any_of(","));
                            
                if (n == 0) beginOffset = vOffsets[0];
                endOffset = vOffsets[1];
                
                if (n == (unsigned)(NQUERYPERWORKITEM - 1)) {
                    endOffset = vOffsets[1];
                    for (size_t j = 0; j < (unsigned)NTOTALDBCHUNKS; ++j) { 
                        vWorkItems.push_back(beginOffset + "," + endOffset 
                            + "," + vDbChunkNames[j]);
                    }
                    nQueryBlocks++;
                    n = 0;
                }      
                else n++;
            }        
            if (n != 0) { /// If (total num seqs % NQUERYPERWORKITEM != 0)
                for (size_t j = 0; j < (unsigned)NTOTALDBCHUNKS; ++j) { 
                    vWorkItems.push_back(beginOffset + "," + endOffset + ","
                        + vDbChunkNames[j]);
                }
                nQueryBlocks++;
            }
        }
        else { /// If NQUERYPERWORKITEM = 0, Use all sequences as 1 block.
            getline(indexFile, line);
            vector<string> vOffsets;
            boost::split(vOffsets, line, boost::is_any_of(","));
            beginOffset = vOffsets[0];
            string temp;
            do { temp = line; } while (!getline(indexFile, line).eof());
            boost::split(vOffsets, temp, boost::is_any_of(","));   
            endOffset = vOffsets[1];   
            for (size_t j = 0; j < (unsigned)NTOTALDBCHUNKS; ++j) { 
                vWorkItems.push_back(beginOffset + "," + endOffset + ","
                    + vDbChunkNames[j]);
            }       
            nQueryBlocks = 1;      
        }        
        indexFile.close();
        nWorkItems = vWorkItems.size();
        vDbChunkNames.clear();
        if (LOGORNOT) {
            LOG << "[INFO] number of query blocks = " << nQueryBlocks << endl;     
            LOG << "[INFO] number of DB partitions = " << NTOTALDBCHUNKS << endl;     
            LOG << "[INFO] number of work items = " << nWorkItems << endl;
        }
        
        ///
        /// Split work item master file into sub master files and process the 
        /// work items iteratively while saving the results at each iteration.
        /// eq) if the total number of query files x and NITERATION = x/2, 
        /// the master file is divided into 2 sub master files and processed 
        /// individually.
        ///
        uint32_t nRemains = 0;
        if (NITERATION > 1) {
            nSubWorkItemFiles = NITERATION;
            nWorkItemsPerFile = (nQueryBlocks / NITERATION) * NTOTALDBCHUNKS;            
            nRemains = (nQueryBlocks % NITERATION) * NTOTALDBCHUNKS;              
        }
        else { /// If NITERATION = 1
            nWorkItemsPerFile = nWorkItems;
            nRemains = 0;
            nSubWorkItemFiles = 1;
        }
                
        size_t k = 0;
        size_t i = 0;
        
        for (; i < nSubWorkItemFiles; ++i) {
            string newFileName = OUTPREFIX + "-workitems-" 
                + boost::lexical_cast<string>(i) + ".txt";
            ofstream workItemFile(newFileName.c_str());
            if (LOGORNOT) LOG << "[INFO] Work item file name (" << i << ") = " 
                 << newFileName << endl;
            if (workItemFile.is_open()) {
                for (size_t j = 0; j < nWorkItemsPerFile; ++j, ++k) {
                    workItemFile << vWorkItems[k] << endl;
                }
                workItemFile.close();
            }
            else {
                cerr << "ERROR: work item file open error.\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
        
        if (nRemains) {
            string newFileName = OUTPREFIX + "-workitems-" 
                + boost::lexical_cast<string>(i) + ".txt";
            ofstream workItemFile(newFileName.c_str());
            if (LOGORNOT) LOG << "[INFO] Work item file name (remains) = " 
                 << newFileName << endl;
            if (workItemFile.is_open()) {
                for (size_t j = 0; j < nRemains; ++j, ++k) {
                    workItemFile << vWorkItems[k] << endl;
                }
                workItemFile.close();
                nSubWorkItemFiles++; /// for one remains file
            }
            else {
                cerr << "ERROR: work item file open error.\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
        if (LOGORNOT) LOG << "[INFO] Total number of sub work item files = " 
             << nSubWorkItemFiles << endl;
    } /// master
    MPI_Barrier(MPI_COMM_WORLD); 
        
    ///
    /// Iteratively call blast and save results for nSubWorkItemFiles
    /// times.
    ///

    /// Broadcast the number of sub work item files
    MPI_Bcast(&nSubWorkItemFiles, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    ///
    /// n iterations
    ///
    for (size_t n = 0; n < nSubWorkItemFiles; ++n) {
        if (LOGORNOT) LOG << "[INFO] iteration number = " << n << endl;   
        string subMasterFileName = OUTPREFIX + "-workitems-" 
            + boost::lexical_cast<string>(n) + ".txt";
   
        ///
        /// map, collate reduce
        ///
        uint32_t nvecRes;
        double map_time;
        if (LOGORNOT) {
            map_time = MPI_Wtime();
            LOG << LOGMSG 
                << "map() starts: " <<  map_time << endl;
        }
        
        ////////////////////////////////////////////////////
        ///
        /// For mrmpi 26AUG2010 version
        ///
        //nvecRes = mr2->map((char*)subMasterFileName.c_str(), 
                           //&mr_run_blast, (void*)NULL);
                           
        ///
        /// For mrmpi 11FEB2011 version
        ///
        //uint64_t MapReduce::map(
            //int nstr,         = 1
            //char **strings,   = subMasterFileName.c_str()
            //int self,         = 0
            //int recurse,      = 0
            //int readfile,     = 1 for using master file(s)
            //void (*mymap)(int, char *, KeyValue *, void *), void *ptr)
        char **masterFileLists;
        int nFiles = 1;
        masterFileLists = (char **) malloc(nFiles);
        masterFileLists[0] = (char *) malloc(subMasterFileName.length() + 1);
        strcpy(masterFileLists[0], subMasterFileName.c_str());

        nvecRes = mr2->map(1, masterFileLists, 0, 0, 1,
                           &mr_run_blast, (void*)NULL);
        ////////////////////////////////////////////////////
        
        if (LOGORNOT) LOG << LOGMSG 
            << "map() ends: " <<  MPI_Wtime() - map_time << endl;
                            
        WORKEROUTPUTFILENAME = OUTPREFIX + "-hits-" 
            + boost::lexical_cast<string>(n) + "-" 
            + boost::lexical_cast<string>(MYID) + ".txt";
        
        double collate_time;
        if (LOGORNOT) {
             collate_time = MPI_Wtime();
             LOG << LOGMSG << "collate starts: " << collate_time << endl;
        }
        
        ///////////////////
        mr2->collate(NULL);
        ///////////////////
        
        if (LOGORNOT) LOG << LOGMSG 
            << "collate ends: " << MPI_Wtime() - collate_time << endl;
        
        double reduce_time;
        if (LOGORNOT) {
             reduce_time= MPI_Wtime();
             LOG << LOGMSG << "reduce starts: " << reduce_time << endl;
        }
        
        //////////////////////////////////////////////////
        mr2->reduce(&mr_sort_multivalues_by_evalue, NULL);
        ////////////////////////////////////////////////////
        
        if (LOGORNOT) LOG << LOGMSG 
            << "reduce ends: " <<  MPI_Wtime() - reduce_time << endl;
        
        ///
        /// Save history
        ///
        if (MYID == 0) {                
            string histFileName = OUTPREFIX + "-history.txt";
            ofstream histFile(histFileName.c_str(), ios::out | ios::app);
            if (!histFile) {
                cerr << "ERROR: failed to open a history file" << endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            else {
                time_t timer;
                timer = time(NULL);
                histFile << subMasterFileName << "," 
                         << asctime(localtime(&timer));
            }
            histFile.close();
        }            
        MPI_Barrier(MPI_COMM_WORLD);           
    }
         
    if (TIMING) {
        /// Wall-clock time
        gettimeofday(&totalEndTime, NULL);
        double tS = totalStartTime.tv_sec*1000000 + (totalStartTime.tv_usec);
        double tE = totalEndTime.tv_sec*1000000  + (totalEndTime.tv_usec);
                
        /// process time
        getrusage(RUSAGE_SELF, &ru_total);
        totalEnd_u_Time = ru_total.ru_utime;
        totalEnd_s_Time = ru_total.ru_stime;
        double tS_user = totalStart_u_Time.tv_sec*1000000 
            + (totalStart_u_Time.tv_usec);
        double tE_user = totalEnd_u_Time.tv_sec*1000000  
            + (totalEnd_u_Time.tv_usec);
        double tS_sys = totalStart_s_Time.tv_sec*1000000 
            + (totalStart_s_Time.tv_usec);
        double tE_sys = totalEnd_s_Time.tv_sec*1000000  
            + (totalEnd_s_Time.tv_usec);
        
        double profile_etime = MPI_Wtime();
        LOG << "mr-mpi-blast ends,"
            << profile_etime << ","
            << totalEndTime.tv_sec*1000000 + totalEndTime.tv_usec << ","
            << totalEnd_u_Time.tv_sec*1000000 + totalEnd_u_Time.tv_usec << ","
            << totalEnd_s_Time.tv_sec*1000000 + totalEnd_s_Time.tv_usec << endl;
        
        /// MP_Wtime
        profile_time = profile_etime - profile_time;
        LOG << "Total MPI_Wtime time," << profile_time << endl;
        LOG << "Total wall-clock time," 
            << tE - tS << "," << (tE - tS) / 1000000 << endl;
        LOG << "Total process time (user)," 
            << tE_user - tS_user << "," << (tE_user - tS_user) / 1000000 << endl;
        LOG << "Total process time (sys)," 
            << tE_sys - tS_sys << "," << (tE_sys - tS_sys) / 1000000 << endl;
        LOG << "Total total process time (user+sys)," 
            << (tE_user - tS_user) + (tE_sys - tS_sys) << ","
            << ((tE_user - tS_user) + (tE_sys - tS_sys)) / 1000000 << endl;
        LOG.flush();
    }
    
    if (LOGORNOT || TIMING) LOGFILE.close();    
    
    if (MYID == 0) cout << "Done!" << endl;   
    
    delete mr2;
    delete pTargetDb;
    MMAPFILE.close();
    MPI_Finalize();
            
            
    return 0;
}

/** MR-MPI Map function - settup NCBI C++ Toolkit env and call blast
 * @param itask
 * @param file
 * @param kv
 * @param ptr
 */

void mr_run_blast(int itask, 
                  char *file, 
                  KeyValue *kv, 
                  void *ptr)
{
    struct timeval blastcallStartTime;
    struct timeval blastcallEndTime;
    struct timeval blastcallStart_u_Time;
    struct timeval blastcallEnd_u_Time;
    struct timeval blastcallStart_s_Time;
    struct timeval blastcallEnd_s_Time;    
    struct rusage ru_blastcall; 
    
    struct timeval qBuildStartTime;
    struct timeval qBuildStart_u_Time;
    struct timeval qBuildStart_s_Time; 
    struct rusage ru_qBuild; 
    
    /// 
    /// Set Blast opotions from file.
    ///
    ifstream strategyFile(STRATEGYFILENAME.c_str(), ios::in);
    if (!strategyFile.is_open()) {
        cerr << "ERROR: failed to open a search strategy file" << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    CNcbiIstream* in = &strategyFile;
    CRef<CBlast4_request> b4req;
    try { 
        b4req = ExtractBlast4Request(*in);
    } catch (const CSerialException&) {
        NCBI_THROW(CInputException, eInvalidInput, 
                   "Failed to read search strategy");
    }
    CImportStrategy strat(b4req);
    CRef<blast::CBlastOptionsHandle> opts = strat.GetOptionsHandle();    
    opts->Validate();

    if (MYID == 1 && OPTDUMP == 1) {
        OPTDUMP = 0;
        opts->GetOptions().DebugDumpText(cout, "opts", 1);
    }

    CRef<CObjectManager> objmgr = CObjectManager::GetInstance();
    if (!objmgr) {
        throw runtime_error("Could not initialize object manager");
    }

    SDataLoaderConfig dlconfig(ISPROTEIN);
    CBlastInputSourceConfig iconfig(dlconfig);
    
    ///
    /// Split work item: <beginOffset,endOffset,db partition name>
    ///
    string workItem(file);
    vector<string> vWorkItemTokens;
    boost::split(vWorkItemTokens, workItem, boost::is_any_of(","));
    assert(vWorkItemTokens.size() == 3);

    ///
    /// Target db name setting
    ///
    string dbChunkName = vWorkItemTokens[2];
    if(pTargetDb == 0 || dbChunkName != prevDbChunkName) {
        delete pTargetDb;
        if (ISPROTEIN) 
            pTargetDb = new CSearchDatabase(dbChunkName,
                                        CSearchDatabase::eBlastDbIsProtein);
        else
            pTargetDb = new CSearchDatabase(dbChunkName,
                                        CSearchDatabase::eBlastDbIsNucleotide);
    }
    prevDbChunkName = dbChunkName;

    ///
    /// Read sequence block and run blast
    ///
    uint32_t beginOffset = boost::lexical_cast<uint32_t>(vWorkItemTokens[0]);
    uint32_t endOffset   = boost::lexical_cast<uint32_t>(vWorkItemTokens[1]);
    
    /// query building timing
    double query_build_time;
    if (TIMING) {
        CNT++;
        query_build_time = MPI_Wtime();
        /// Wall-clock time
        gettimeofday(&qBuildStartTime, NULL);
        /// Process time
        getrusage(RUSAGE_SELF, &ru_qBuild);
        qBuildStart_u_Time = ru_qBuild.ru_utime;
        qBuildStart_s_Time = ru_qBuild.ru_stime;
        LOG << LOGMSG << "query build starts," 
            << query_build_time << ","
            << qBuildStartTime.tv_sec*1000000 
                + qBuildStartTime.tv_usec << ","
            << qBuildStart_u_Time.tv_sec*1000000 
                + qBuildStart_u_Time.tv_usec << ","
            << qBuildStart_s_Time.tv_sec*1000000 
                + qBuildStart_s_Time.tv_usec << ","
            << dbChunkName << "," << CNT << "," << MPI_procName << endl;
        LOG.flush();
    }
    
    
    //FILE *qf = fopen((char*)QUERYFILENAME.c_str(), "r");
    
    ///
    /// Read a block of sequeces from beginOffset
    ///
    const char* mmapQueryFile = (char*)MMAPFILE.data();
    string query(mmapQueryFile + beginOffset, endOffset - beginOffset);
    //cout << query << endl;

    //fseek(qf, beginOffset, SEEK_SET);    
    //string query;
    //vector<string> vHeaders;
    //vector<uint32_t> vBeginOffsets;
    //uint32_t nQuery = 0;
   
    //while ((uint32_t)ftell(qf) < endOffset) {
        //vBeginOffsets.push_back((uint32_t)ftell(qf));
        //fgets(buff2, MAXSTR2, qf);
        //if (buff2[0] == '>') {
            ///// Here I collect headers of fasta input seqs
            ///// for retrieving query sequence info.
            //if (EXCLUSIONORNOT) vHeaders.push_back(string(buff2));            
            //nQuery++;
        //}
        ///// Concatenate query sequces for Blast call
        //query += string(buff2);
    //}    
    //fclose(qf);       
 
    vector<string> vHeaders;        /// For collecting def lines of queries
    vector<uint32_t> vBeginOffsets; /// Fpr collecting beginoffsets of queries
    uint32_t nQuery = 0;
    char* c = (char*)(mmapQueryFile + beginOffset);
    uint32_t loc = beginOffset;
    char buff2[MAXSTR2];
    size_t buffIdx = 0;
    while (loc < endOffset) {
        if ((*c) == '>') {
            nQuery++;
            vBeginOffsets.push_back(loc);
            buffIdx = 0;
        }
        if (EXCLUSIONORNOT && (*c) != '\n') {
            if (buffIdx < MAXSTR2) buff2[buffIdx++] = (*c);
        }
        else if (EXCLUSIONORNOT && (*c) == '\n') 
            vHeaders.push_back(string(buff2));  
        c++;
        loc++;
    }
    if (LOGORNOT) LOG << LOGMSG 
        << "Number of queries for a Blast call = " << nQuery << endl;
 
    ///
    /// Set queries as fasta input
    ///
    CBlastFastaInputSource fasta_input(query, iconfig);
    CBlastInput blastInput(&fasta_input);
    CScope scope(*objmgr);
    TSeqLocVector queryLoc = blastInput.GetAllSeqLocs(scope);
    CRef<IQueryFactory> queryFactory(new CObjMgr_QueryFactory(queryLoc));    
    
    ///
    /// Run blast
    ///
    CLocalBlast blaster(queryFactory, opts, *pTargetDb);
         
    double blast_call_time;
    if (TIMING) {
        blast_call_time = MPI_Wtime();
        /// Wall-clock time
        gettimeofday(&blastcallStartTime, NULL);
        /// Process time
        getrusage(RUSAGE_SELF, &ru_blastcall);
        blastcallStart_u_Time = ru_blastcall.ru_utime;
        blastcallStart_s_Time = ru_blastcall.ru_stime;
        LOG << LOGMSG << "blast call starts," 
            << blast_call_time << ","
            << blastcallStartTime.tv_sec*1000000 
                + blastcallStartTime.tv_usec << ","
            << blastcallStart_u_Time.tv_sec*1000000 
                + blastcallStart_u_Time.tv_usec << ","
            << blastcallStart_s_Time.tv_sec*1000000 
                + blastcallStart_s_Time.tv_usec << ","
            << dbChunkName << "," << CNT << "," << MPI_procName << endl;
        LOG.flush();
    }
        
    //////////////////////////////////////////
    CSearchResultSet results = *blaster.Run();
    //////////////////////////////////////////
    
    if (TIMING) {
        double blast_call_etime = MPI_Wtime();
        
        /// Wall-clock time
        gettimeofday(&blastcallEndTime, NULL);
                
        /// process time
        getrusage(RUSAGE_SELF, &ru_blastcall);
        blastcallEnd_u_Time = ru_blastcall.ru_utime;
        blastcallEnd_s_Time = ru_blastcall.ru_stime;
                        
        LOG << LOGMSG << "blast call ends,"
            << blast_call_etime << ","
            << blastcallEndTime.tv_sec*1000000 
                + blastcallEndTime.tv_usec << ","
            << blastcallEnd_u_Time.tv_sec*1000000 
                + blastcallEnd_u_Time.tv_usec << ","
            << blastcallEnd_s_Time.tv_sec*1000000 
                + blastcallEnd_s_Time.tv_usec << "," 
            << dbChunkName << "," << CNT << "," << MPI_procName << endl;
        LOG.flush();
    }    
    
    ///
    /// Get warning messages
    ///
    for (size_t i = 0; i < results.GetNumResults(); ++i) {
        TQueryMessages messages = results[i].GetErrors(eBlastSevWarning);
        if (messages.size() > 0) {
            CConstRef<CSeq_id> seq_id = results[i].GetSeqId();
            if (seq_id.NotEmpty())
                cerr << "ID: " << seq_id->AsFastaString() << endl;
            else
                cerr << "ID: " << "Unknown" << endl;

            ITERATE(vector<CRef<CSearchMessage> >, it, messages) {
                cerr << (*it)->GetMessage() << endl;
            }
        }
    }

    ///
    /// Get the results
    ///   
    for (size_t i = 0; i < results.GetNumResults(); ++i) {

        CConstRef<CSeq_align_set> aln_set = results[i].GetSeqAlign();

        if (results[i].HasAlignments()) {
            ITERATE(CSeq_align_set::Tdata, itr, aln_set->Get()) {
                const CSeq_align& s = **itr;
                
                ///
                /// VERY IMPORTANT!!
                /// queryID is not unique. It's internal qid in fasta_input.
                ///
                string queryID = s.GetSeq_id(QUERY).GetSeqIdString(); 
                string subID = s.GetSeq_id(SUBJECT).GetSeqIdString(); /// GI
                
                ///
                /// Refer: http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/
                /// doxyhtml/classCSeq__align.html
                ///
                double pIdentity=0.0;
                s.GetNamedScore(CSeq_align::eScore_PercentIdentity, pIdentity);

                uint32_t alignLen=0, qStart=0, qEnd=0;
                uint32_t sStart=0, sEnd=0;
                uint32_t gapOpens=0;
                gapOpens    = s.GetNumGapOpenings();
                qStart      = s.GetSeqStart(QUERY);
                qEnd        = s.GetSeqStop(QUERY);
                sStart      = s.GetSeqStart(SUBJECT);
                sEnd        = s.GetSeqStop(SUBJECT);
                alignLen    = s.GetAlignLength();
                
                double eValue=0.0;
                int bitScore=0;
                int misMatches=0;
                
                s.GetNamedScore(CSeq_align::eScore_EValue, eValue);
                s.GetNamedScore(CSeq_align::eScore_BitScore, bitScore);
                s.GetNamedScore(CSeq_align::eScore_MismatchCount, misMatches);
                
                /// Get the saved beginOffset from vBeginOffsets
                uint32_t uniqueQID = 
                    vBeginOffsets[boost::lexical_cast<uint32_t>(queryID)-1];
                
                /// 
                /// Add a csv blast result to kv
                /// 
                if (EXCLUSIONORNOT) {
                    /// Tokenize query header
                    string qHeader = 
                        vHeaders[boost::lexical_cast<uint32_t>(queryID)-1];
                    vector<string> vQueryId;
                    boost::split(vQueryId, qHeader, boost::is_any_of("|"));
                                        
                    ///
                    /// The below infomation is only for our own simulated
                    /// sequence data sets which has 
                    /// - origin GI
                    /// - unique query ID                        
                    /// - original sequence length
                    /// - cut location start
                    /// - cut location end
                    ///
                    string qGi = vQueryId[1];   /// GI
                    //uint32_t origLen            /// length of the orig seq         
                        //= boost::lexical_cast<uint32_t>(vQueryId[3]);     
                    int qCutLocStart            /// cut coordinates - start
                        = boost::lexical_cast<int>(vQueryId[4]);      
                    int qCutLocEnd         /// cut coordinates - end
                        = boost::lexical_cast<int>(vQueryId[5]);  
                                           
                    if (!check_exclusion(qGi, subID, qCutLocStart, qCutLocEnd, 
                        sStart, sEnd, EXCLUSIONTHRESHOLD)) {
                                                
                        ///
                        /// To pass Blast hits using struct, outfmt=6
                        /// query id, subject id, % identity, alignment length, 
                        /// mismatches, gap opens, q. start, q. end, s. start, 
                        /// s. end, evalue, bit score
                        ///
                        BLASTRES res;
                        strncpy(res.subjectid, subID.c_str(), subID.length());
                        res.identity = pIdentity;
                        res.alignlen = alignLen;
                        res.mismatches = misMatches;
                        res.gapopens = gapOpens;
                        res.qstart = qStart;
                        res.qend = qEnd;
                        res.sstart = sStart;
                        res.send = sEnd;
                        res.evalue = eValue;
                        res.bitscore = bitScore;     
                                                
                        /// 
                        /// ADD <KEY = "QUERYID", VALUE="BLASTRESULT">
                        /// TO KV
                        ///
                        uint32_t newKey = uniqueQID;
                        kv->add((char*)&newKey, sizeof(uint32_t), (char*)&res, 
                        BLASTRESSZ);  
                    }
                    /// 
                    /// Found a self hits. Record the hits in a histroy file
                    ///
                    else {
                        string exFileName = OUTPREFIX + "-" 
                            + boost::lexical_cast<string>(MYID) + "-" 
                            + EXCLUSIONHISTFILE;
                        ofstream exFile(exFileName.c_str(), 
                                        ios::out | ios::app);
                            
                        if (!exFile) {
                            cerr << "ERROR: failed to open a exclusion "
                                 << "history file" << endl;
                            MPI_Abort(MPI_COMM_WORLD, 1);
                        }
                        else {
                            /// Format:
                            /// qGi sGi qCutLocStart qCutLocEnd sStart sEnd
                            exFile << qGi << "\t"
                                   << subID << "\t"
                                   << qCutLocStart << "\t"
                                   << qCutLocEnd << "\t" 
                                   << sStart << "\t"
                                   << sEnd << endl;
                        }
                        exFile.close();
                    }
                }
                else {
                    BLASTRES res;
                    strncpy(res.subjectid, subID.c_str(), subID.length());
                    res.identity    = pIdentity;
                    res.alignlen    = alignLen;
                    res.mismatches  = misMatches;
                    res.gapopens    = gapOpens;
                    res.qstart      = qStart;
                    res.qend        = qEnd;
                    res.sstart      = sStart;
                    res.send        = sEnd;
                    res.evalue      = eValue;
                    res.bitscore    = bitScore;     
                                            
                    uint32_t newKey = uniqueQID;
                    kv->add((char*)&newKey, sizeof(uint32_t), (char*)&res, 
                        BLASTRESSZ);                    
                }
            }
        }
    }
}


/** Sort function - Passed to MR-MPI sort_values() for sorting blast result
 * string by bit score.
 * @param e1
 * @param e2
 */
 
inline bool mycompare(STRUCTEVALUE e1, 
                      STRUCTEVALUE e2)
{
    return (e1.evalue < e2.evalue);
}

/** Sort by wvalue - Passed to MR-MPI reduce() for sorting KMVs by evalue.
 * @param key
 * @param keybytes
 * @param multivalue - collected blast result strings. There could be more than
 * nvalues results in it. The separator '>' should be used for splitting the 
 * resutls.
 * @param nvalues
 * @param valuebytes
 * @param kv
 * @param ptr
 */
 
void mr_sort_multivalues_by_evalue(char *key, 
                                   int keybytes, 
                                   char *multivalue,
                                   int nvalues, 
                                   int *valuebytes, 
                                   KeyValue *kv, 
                                   void *ptr) 
{   
    /// Check if there is KMV overflow
    assert(multivalue != NULL && nvalues != 0);
    
    ///
    /// Make STRUCTEVALUE = {BLASTRES* p; double evalue;}
    /// and sort by evalue
    ///
    vector<STRUCTEVALUE> vforsort;
    for (size_t n = 0; n < (unsigned)nvalues; n++) {                
        BLASTRES* res = (BLASTRES*)multivalue;
        STRUCTEVALUE structEvalue;
        structEvalue.p = res;
        structEvalue.evalue = res->evalue;
        vforsort.push_back(structEvalue);
        multivalue += BLASTRESSZ;         
    }
    sort(vforsort.begin(), vforsort.end(), mycompare);
    
    ///
    /// After sorting, each worker saves the set of results into a file.
    /// Note: The file open option is "a".
    /// Note: outfmt = 6 in Blast
    ///
    ofstream outputFile(WORKEROUTPUTFILENAME.c_str(), ios::out | ios::app);
    for (size_t n = 0; n < (unsigned)nvalues; n++) {
        BLASTRES* res = (BLASTRES*)(vforsort[n].p);
        outputFile << *(uint32_t *)key << "\t"
            << res->subjectid << "\t"
            << res->identity << "\t"
            << res->alignlen << "\t"
            << res->mismatches << "\t"
            << res->gapopens << "\t"
            << res->qstart << "\t"
            << res->qend << "\t"
            << res->sstart << "\t"
            << res->send << "\t"
            << res->evalue << "\t"
            << res->bitscore 
            << endl;
    }
    outputFile.close();
    vforsort.clear();
}   

/** Check exclusion - Based on the coordinates of query and subject, decide 
 * whether the result should be included in the final result or not.
 * @param qGi: query ID
 * @param sGi: subject ID
 * @param qCutLocStart: Cutting start location of the query from original 
 * refseq_genomic fasta input.
 * @param qCutLocEnd: Cutting end location of the query from original 
 * refseq_genomic fasta input.
 * @param sStart: subject alignment start location.
 * @param sEnd: subject alignemnt end location.
 * @param threshold: overlap threshold (default = 100bp).
 */
 
inline bool check_exclusion(string qGi, 
                            string sGi, 
                            int qCutLocStart,
                            int qCutLocEnd, 
                            int sStart, 
                            int sEnd,
                            int threshold)
{
    /// 
    /// To exclude Blast result from the original sequence from which
    /// the input query is originated (sampled). Basically if qGi == sGi
    /// and qCutLocStart is similar with sStart and qCutLocEnd is similar
    /// with sEnd in terms of coordinates, the result should be excluded.
    /// 
    /// Orig seq: ----------------XXXXXXXXXXXXXXXX----------------------
    ///                           |              |
    ///                     qCutLocStart      qCutLocEnd
    /// 
    /// Query:                    XXXXXXXXXXXXXXXX
    ///                              |          |
    ///                           qStart       qEnd
    /// 
    /// Subject:   ------------------XXXXXXXXXXXX-----------------------
    ///                              |          |
    ///                           sStart      sEnd
    /// 
    bool ret = false;

    if (qGi == sGi) {
        if (qCutLocStart < 0) {
            /// 
            /// In >gi|222299657|18|3605|-400|3604
            /// -400|3604 means query[-400:3604] in Python.
            /// 
            qCutLocStart = qCutLocEnd + 1 - qCutLocStart;
            qCutLocEnd += 1;
        }
        if ((qCutLocStart - threshold <= sStart 
            && sStart <= qCutLocStart + threshold) 
            && (qCutLocEnd - threshold <= sEnd 
            && sEnd <= qCutLocEnd + threshold))
            ret = true;
    }

    return ret;
}

/// EOF

