////////////////////////////////////////////////////////////////////////////////
//
//  MR-MPI-BLAST: Parallelizing BLAST with MR-MPI
//
//  Author: Seung-Jin Sul
//         (ssul@jcvi.org)
//
//  Last updated: 05/23/2011
//
////////////////////////////////////////////////////////////////////////////////

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

/// 
/// NCBI C++ Toolkit
/// 
#include <algo/blast/api/local_blast.hpp>
#include <algo/blast/api/objmgr_query_data.hpp>
#include <algo/blast/blastinput/blast_fasta_input.hpp>

/// Import search strategy
#include <objects/blast/Blast4_request.hpp>
#include <algo/blast/api/search_strategy.hpp>

/// Tabular
#include <algo/blast/format/blast_format.hpp>
#include <algo/blast/format/blastfmtutil.hpp> 
#include <objtools/align_format/tabular.hpp>

USING_NCBI_SCOPE;
USING_SCOPE(blast);

/// For typedef unsigned long long int uint32_t
#include <stdint.h>


/// 
/// Boost lib
/// 

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
#include <boost/filesystem/operations.hpp> /// for real file size
boost::iostreams::mapped_file_source g_memmapQueryFile; /// Read-only Boost mmap query file
uint32_t g_realFileSize;


/// 
/// Settings from mrblast.ini conf file
/// 

/// MR-MPI options
int g_verbosity;
int g_timer;                /// log elapsed time for each mapreduce call
int g_memSize;              /// page size (in Mbytes)
int g_outOfCore;

/// Self-Exclusion options
int g_exclEnabled;
int g_exclThreshold;        /// Exclusion threshold = 100bp
string g_exclHistFileName;

/// DB options
string g_dbFileName;
string g_configFileName;
int g_numDbFiles;
const int MAXSTR = 80;      /// For mpi proc name and query header

/// Filtering
double g_IDENT_CUTOFF = 0.5; /// Doug's identity for filtering
double g_COVER_CUTOFF = 0.9; /// Doug's coverage for filtering

/// 
/// Log
/// 
#include <sys/time.h>
#include <sys/resource.h>
//#define NDEBUG 1
int g_logEnabled = 0;
int g_timingEnabled = 0;
int g_optDumpEnabled = 0;   /// For dumping Blast pBlOpts out
int g_mapCallNo = 0;        /// Serial no of map calls in each rank
string g_logFileName;
string g_logMsg;
ofstream g_logFileStream;
static ostream *g_pLogStream = NULL;
#define LOG (*g_pLogStream)


/// 
/// New scheduler
/// 
typedef struct structIndex {
    uint32_t qStart;
    uint32_t qLength;
} structQueryIndex_t;
typedef struct structWorkItem {
    uint32_t dbNo;
    uint32_t blockBegin;
    uint32_t blockEnd;
} structWorkItem_t;
vector<string> g_vecDbFile;
vector<uint32_t> g_vecBlockBeginLoc;
vector<structQueryIndex_t> g_vecQueryIndex;
vector<structWorkItem_t> g_vecWorkItem;
uint32_t g_numQueries;
uint32_t g_numQueryBlocks;
uint32_t g_numWorkItems;
uint32_t g_blockSize; /// block size in base-pair
multimap<string, int> g_multimapProcNameRank; /// dict of rank by proc name
map<int, string> g_mapRankProcName;           /// dict of proc name by rank 
 
 
/// 
/// Misc.
/// 
const int SUB_ID_LEN = 20;  /// For KVing subject ID of blast hits. possible KMvec overflow
const int QUERY = 0;        /// To retireve query info from CSeq_align
const int SUBJECT = 1;      /// To retireve suject info from CSeq_align
string g_outFilePrefix;     /// Prefix string for output file names
string g_indexFileName;
string g_queryFileName;
int g_mapStyle;
int g_MPI_worldRank;         /// MPI rank
char g_MPI_procName[MAXSTR]; /// MPI procname
int g_MPI_numProcs;

/// For nucl or prot DB setting
bool g_bIsProtein = false;

/// Blast target DB setting
static CSearchDatabase *g_pTargetDb = 0;
string g_prevDbName = "";

/// Output file to store BLAST hits
string g_hitFileName = "";

/// Import search strategy
string g_strategyFileName;  /// Input blast search option file

/// For syncronized timing
#ifndef MPI_WTIME_IS_GLOBAL
#define MPI_WTIME_IS_GLOBAL 1
#endif


/// To pass Blast hits following outfmt=6 format.
/// subject id, % identity, alignment length, nMismatches, gap opens,
/// q. start, q. end, s. start, s. end, evalue, bit score
typedef struct structBlRes {
    char subjectId[SUB_ID_LEN];
    double identity;
    uint32_t alignLen;
    int misMatches;
    uint32_t gapOpens;
    uint32_t qStart;
    uint32_t qEnd;
    uint32_t sStart;
    uint32_t sEnd;
    double evalue;
    int bitScore;
    uint32_t cutStart;
    uint32_t cutEnd;
    double doug_identity;
    double doug_coverage;
} structBlRes_t;

 
/// To sort Blast hits by evalue
typedef struct structEvalue {
    structBlRes_t *pRec;
    char* pSubjectId;
    double evalue;
    int bitScore;
} structEValue_t;

/// Multiple iterations
int g_numIter = 1;   
uint32_t nSubWorkItemFiles = 0;
uint32_t nWorkItemsPerIter = 0;
uint32_t nRemains = 0;

/// To pass info to map() and reduce()
typedef struct structToPass {
    int rank;
    int iter;
} structToPass_t;

 
///
/// Function declarations
///
void        collect_mpi_node_name(int rank, int numProcs, MPI_Comm mpiComm);
void        run_mr_mpi_blast(MPI_Comm mpiComm, int rank); 
void        mr_run_blast(int itask, KeyValue *kv, void *ptr);
void        mr_sort_multivalues_by_evalue(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
inline bool evalue_compare(structEValue_t e1, structEValue_t e2);
inline bool check_exclusion(string qGi, string sGi, int qCutLocStart, int qCutLocEnd, int sStart, int sEnd, int threshold);                        

 
int main(int argc, char **argv)
{
    po::options_description generalDesc("General options");
    generalDesc.add_options()
    ("help,h", "print help message")
    ("query-file,i", po::value<string>(), "set input query file")
    ("index-file,d", po::value<string>(), "set input index file")
    ("import-search-strategy,s", po::value<string>(), "set search strategy file")
    ("db-list,l", po::value<string>(), "set DB partition name list file")
    ("conf-file,c", po::value<string>(), "set configuration file")
    ("block-size,b", po::value<int>(), "set the number of base-pairs per work item")
    ;

    po::options_description OptionalDesc("Optional options");
    OptionalDesc.add_options()
    ("output-prefix,o", po::value<string>(&g_outFilePrefix)->default_value("output"),
     "set output prefix for output file names (default=output)")
    ("map-style,m", po::value<int>(&g_mapStyle)->default_value(2),
     "set MR-MPI mapstyle: 2=master/slave, 3=new scheduler")
    ("is-protein,p", po::value<bool>(&g_bIsProtein)->default_value(false),
     "set db type (default=false)")
    ("iteration,n", po::value<int>(&g_numIter)->default_value(1), 
     "set the number of iterations")          
    ;

    po::options_description allDesc("Allowed options");
    allDesc.add(generalDesc).add(OptionalDesc);

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, allDesc), vm);
    po::notify(vm);

    if (
        argc < 2 || (!strcmp(argv[1], "-?") || !strcmp(argv[1], "--?")
                 || !strcmp(argv[1], "/?") || !strcmp(argv[1], "/h")
                 || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--h")
                 || !strcmp(argv[1], "--help") || !strcmp(argv[1], "/help")
                 || !strcmp(argv[1], "-help")  || !strcmp(argv[1], "help"))) {
        cout << "MR-MPI Blast\n" << "Author: Seung-Jin Sul (ssul@jcvi.org)\n\n"
             << allDesc;
        return 1;
    }
    else {
        if (vm.count("query-file"))
            g_queryFileName = vm["query-file"].as<string>();
        else {
            cerr << "ERROR: query file was not set.\n\n";
            cout << allDesc;
            return 1;
        }
        if (vm.count("index-file"))
            g_indexFileName = vm["index-file"].as<string>();
        else {
            cerr << "ERROR: index file was not set.\n\n";
            cout << allDesc;
            return 1;
        }
        if (vm.count("import-search-strategy"))
            g_strategyFileName = vm["import-search-strategy"].as<string>();
        else {
            cerr << "ERROR: option file was not set.\n\n";
            cout << allDesc;
            return 1;
        }
        if (vm.count("db-list"))
            g_dbFileName = vm["db-list"].as<string>();
        else {
            cerr << "ERROR: DB name list file was not set.\n\n";
            cout << allDesc;
            return 1;
        }
        if (vm.count("conf-file"))
            g_configFileName = vm["conf-file"].as<string>();
        else {
            cerr << "ERROR: configuration file was not set.\n\n";
            cout << allDesc;
            return 1;
        }
        if (vm.count("block-size"))
            g_blockSize = vm["block-size"].as<int>();
        else {
            cerr << "ERROR: block size was not set.\n\n";
            cout << allDesc;
            return 1;
        }

        ///
        /// Read OPTIONAL arguments
        ///
        if (vm.count("output-prefix"))
            g_outFilePrefix = vm["output-prefix"].as<string>();
        if (vm.count("map-style")) g_mapStyle = vm["map-style"].as<int>();
        if (vm.count("iteration")) 
            g_numIter = vm["iteration"].as<int>();
    }

    ///
    /// Read conf file, mrblast.ini and set parameters
    ///
    ifstream config(g_configFileName.c_str(), ios::in);
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
        g_exclHistFileName = parameters["EXCLHISTFNAME"];
        try {
            g_verbosity = boost::lexical_cast<int>(parameters["VERBOSITY"]);
            g_timer = boost::lexical_cast<int>(parameters["TIMER"]);
            g_memSize = boost::lexical_cast<int>(parameters["MEMSIZE"]);
            g_outOfCore = boost::lexical_cast<int>(parameters["OUTOFCORE"]);

            g_exclEnabled = boost::lexical_cast<int>(parameters["EXCLENABLED"]);
            g_exclThreshold = boost::lexical_cast<int>(parameters["EXCLTHRESHOLD"]);

            g_logEnabled = boost::lexical_cast<int>(parameters["LOGENABLED"]);
            g_timingEnabled = boost::lexical_cast<int>(parameters["TIMING"]);
            g_optDumpEnabled = boost::lexical_cast<int>(parameters["OPTDUMP"]);
            g_logFileName = parameters["LOGFNAME"];
            
            g_IDENT_CUTOFF = boost::lexical_cast<double>(parameters["IDENTCUTOFF"]);
            g_COVER_CUTOFF = boost::lexical_cast<double>(parameters["COVERCUTOFF"]);
        }
        catch (const boost::bad_lexical_cast &) {
            cerr << "Exception: bad_lexical_cast" << endl;
        }
    }
    catch (exception& e) {
        cerr << "Exception: " << e.what() << endl;
    }


    ///
    /// MPI setup
    ///
    int MPI_procNameLen;    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &g_MPI_worldRank);
    MPI_Comm_size(MPI_COMM_WORLD, &g_MPI_numProcs);
    MPI_Get_processor_name(g_MPI_procName, &MPI_procNameLen);

    ///
    /// Collect MPI node names and rank number for mapstyle=3 scheduler
    ///
    collect_mpi_node_name(g_MPI_worldRank, g_MPI_numProcs, MPI_COMM_WORLD);
    
    ///
    /// Creat memory-mapped file for query
    ///
    uint32_t g_realFileSize = boost::filesystem::file_size(g_queryFileName);
    g_memmapQueryFile.open(g_queryFileName, g_realFileSize, 0);
    if (!g_memmapQueryFile.is_open()) {
        cerr << "ERROR: failed to create mmap query file\n";
        MPI_Finalize();
        exit(1);
    }

    ///
    /// Create work items
    ///

    ///
    /// Load DB list 
    ///
    string line;
    ifstream dbListFile(g_dbFileName.c_str(), ios::in);
    if (!dbListFile.is_open()) {
        cerr << "ERROR: DB list file open error.\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    while (!getline(dbListFile, line).eof() && line.length() > 0)
        g_vecDbFile.push_back(line);
    dbListFile.close();
    g_numDbFiles = g_vecDbFile.size();

    ///
    /// Load index file
    ///
    ifstream indexFile(g_indexFileName.c_str(), ios::in);
    if (!indexFile.is_open()) {
        cerr << "ERROR: index file open error.\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    while (!getline(indexFile, line).eof() && line.length() > 0) {
        vector<string> vecIndexTokens;
        boost::split(vecIndexTokens, line, boost::is_any_of(","));
        structQueryIndex_t queryIndex;
        queryIndex.qStart  = boost::lexical_cast<uint32_t>(vecIndexTokens[0]);
        queryIndex.qLength = boost::lexical_cast<uint32_t>(vecIndexTokens[1]);
        g_vecQueryIndex.push_back(queryIndex);
    }
    indexFile.close();
    g_numQueries = g_vecQueryIndex.size();

    ///
    /// Fill g_vQueryStartLoc vector from the index file using g_blockSize
    ///
    uint32_t qSizeCurr = g_blockSize;
    for (size_t q = 0; q < g_numQueries; q++) {
        if (qSizeCurr >= (unsigned)g_blockSize) {
            g_vecBlockBeginLoc.push_back(g_vecQueryIndex[q].qStart);
            qSizeCurr = 0;
        }
        qSizeCurr += g_vecQueryIndex[q].qLength;
    }
    g_numQueryBlocks = g_vecBlockBeginLoc.size();

    ///
    /// Create work items
    /// Here we iterate query block ID with fixing each DB partition name
    ///
    for (size_t d = 0; d < (unsigned)g_numDbFiles; d++) {
        for (size_t b = 0; b < g_numQueryBlocks - 1; b++) {
            uint32_t blockBegin = g_vecBlockBeginLoc[b];
            uint32_t blockEnd = g_vecBlockBeginLoc[b + 1] - 1;
            structWorkItem_t aWorkItem;
            aWorkItem.blockBegin = blockBegin;
            aWorkItem.blockEnd   = blockEnd;
            aWorkItem.dbNo       = d;
            g_vecWorkItem.push_back(aWorkItem);
        }
        structWorkItem_t aWorkItem;
        aWorkItem.blockBegin = g_vecBlockBeginLoc[g_numQueryBlocks - 1];
        aWorkItem.blockEnd   = g_realFileSize;
        aWorkItem.dbNo       = d;
        g_vecWorkItem.push_back(aWorkItem);
    }
    g_numWorkItems = g_vecWorkItem.size();
    
    if (g_MPI_worldRank == 0) {
        cout << "Number of query blocks = " << g_numQueryBlocks << endl;
        cout << "Number of DB files = " << g_numDbFiles << endl;
        cout << "Number of total work items = " << g_numWorkItems << endl;
        cout.flush();
    }

    ///
    /// Calculate sub work item size for multiple iterations
    ///
    if (g_numIter != 1) {
        nSubWorkItemFiles = g_numIter;
        nWorkItemsPerIter = (g_numQueryBlocks / g_numIter) * g_numDbFiles;            
        nRemains = (g_numQueryBlocks % g_numIter) * g_numDbFiles;
    }
    else { 
        nWorkItemsPerIter = g_numWorkItems;
        nRemains = 0;
        nSubWorkItemFiles = 1;
    }
    /// for one more added file of remaining work items
    if (nRemains) nSubWorkItemFiles++; 

    if (g_MPI_worldRank == 0 && g_numIter != 1) {        
        cout << "Number of sub work item sets = " << nSubWorkItemFiles << endl;
        cout << "Number of work items per iteration  = " << nWorkItemsPerIter << endl;
        cout << "Number of work items remaining = " << nRemains << endl;
    }       
        
    MPI_Barrier(MPI_COMM_WORLD);

    ///
    /// Run mr-mpi calls
    ///
    run_mr_mpi_blast(MPI_COMM_WORLD, g_MPI_worldRank); 
    if (g_MPI_worldRank == 0) cout << "BLAST searching is done!" << endl;
    
    ///
    /// Clean up
    ///
    g_vecDbFile.clear();
    g_vecQueryIndex.clear();
    g_vecBlockBeginLoc.clear();
    g_vecWorkItem.clear();
    delete g_pTargetDb; 
    if (g_logEnabled || g_timingEnabled) g_logFileStream.close();
    g_memmapQueryFile.close(); /// close mmapped file    
    MPI_Finalize();
    
    return 0;
}
 
/** MapReduce fuction for BLAST search
 * @param mpiComm
 * @param rank
 */

void run_mr_mpi_blast(MPI_Comm mpiComm, int rank)
{
    ///
    /// Log file init
    ///
    if (g_logEnabled || g_timingEnabled) {
        g_logFileName = g_outFilePrefix + "-" + boost::lexical_cast<string>(rank) 
                        + "-" + g_logFileName;
        g_logFileStream.open(g_logFileName.c_str(), ios::out);
        g_pLogStream = &g_logFileStream;
    }
    g_logMsg = "Rank:" + boost::lexical_cast<string>(rank) + " ";

    double profileTime;
    struct timeval wallTimeStart;
    struct timeval wallTimeEnd;
    struct timeval userTimeStart;
    struct timeval userTimeEnd;
    struct timeval sysTimeStart;
    struct timeval sysTimeEnd;
    struct rusage ruTotal;

    if (g_timingEnabled) {
        profileTime = MPI_Wtime(); /// MPI_Wtime
        gettimeofday(&wallTimeStart, NULL); /// Wall-clock time
        getrusage(RUSAGE_SELF, &ruTotal); /// Process time
        userTimeStart = ruTotal.ru_utime;
        sysTimeStart = ruTotal.ru_stime;
        LOG << "mr-mpi-blast starts,"
            << profileTime << ","
            << wallTimeStart.tv_sec * 1000000 + wallTimeStart.tv_usec << ","
            << userTimeStart.tv_sec * 1000000 + userTimeStart.tv_usec << ","
            << sysTimeStart.tv_sec * 1000000 + sysTimeStart.tv_usec
            << endl;
        LOG.flush();
    }

    ///
    ///
    /// MR-MPI init (MAY2011 version)
    ///
    /// mapstyle = 0 (chunk) or 1 (stride) or 2 (master/slave)
    /// all2all = 0 (irregular communication) or 1 (use MPI_Alltoallv)
    /// verbosity = 0 (none) or 1 (summary) or 2 (histogrammed)
    /// timer = 0 (none) or 1 (summary) or 2 (histogrammed)
    /// memsize = N = number of Mbytes per page of memory
    /// minpage = N = # of pages to pre-allocate per processor
    /// maxpage = N = max # of pages allocatable per processor
    /// freepage = 1 if memory pages are freed in between operations, 0 if held
    /// outofcore = 1 if even 1-page data sets are forced to disk, 0 if not, -1
    ///               if cannot write to disk
    /// zeropage = 1 if zero out every allocated page, 0 if not
    /// keyalign = N = byte-alignment of keys
    /// valuealign = N = byte-alignment of values
    /// fpath = string
    ///
    MapReduce *pMr = new MapReduce(mpiComm);
    pMr->verbosity = g_verbosity;
    pMr->timer     = g_timer;
    pMr->memsize   = g_memSize;
    pMr->keyalign  = sizeof(uint32_t); /// The key is GI
    pMr->mapstyle  = g_mapStyle;       /// master/slave mode=2, custom scheduler=3
    pMr->outofcore = g_outOfCore;
    MPI_Barrier(mpiComm);

    for (int iter = 0; iter < nSubWorkItemFiles; iter++) {
        if (g_MPI_worldRank == 0) cout << "Iteration: " << iter << endl;
        if (g_logEnabled) LOG << g_logMsg << "Iteration: " << iter << endl;
        
        ///
        /// map, collate reduce
        ///
        double mapTime;
        if (g_logEnabled) {
            mapTime = MPI_Wtime();
            LOG << g_logMsg << "map() starts: " <<  mapTime << endl;
        }
    
        structToPass_t toPass;
        toPass.rank = rank;
        toPass.iter = iter;
        
        /////////////////////////////////////////////////////////////
        if (g_numIter != 1)
            if (nRemains != 0 && iter == nSubWorkItemFiles - 1)
                pMr->map(nRemains, &mr_run_blast, &toPass);
            else 
                pMr->map(nWorkItemsPerIter, &mr_run_blast, &toPass);
        else 
            pMr->map(g_numWorkItems, &mr_run_blast, &toPass);
        /////////////////////////////////////////////////////////////

        if (g_logEnabled)
            LOG << g_logMsg << "map() ends: " <<  MPI_Wtime() - mapTime << endl;

        g_hitFileName = g_outFilePrefix + "-hits-" 
                        + boost::lexical_cast<string>(iter) + "-"
                        + boost::lexical_cast<string>(rank) + ".txt";
                        
        double collateTime;
        if (g_logEnabled) {
            collateTime = MPI_Wtime();
            LOG << g_logMsg << "collate() starts: " << collateTime << endl;
        }

        ///////////////////
        pMr->collate(NULL);
        ///////////////////

        if (g_logEnabled)
            LOG << g_logMsg << "collate() ends: " << MPI_Wtime() - collateTime << endl;

        double reduceTime;
        if (g_logEnabled) {
            reduceTime = MPI_Wtime();
            LOG << g_logMsg << "reduce() starts: " << reduceTime << endl;
        }

        ///////////////////////////////////////////////////
        pMr->reduce(&mr_sort_multivalues_by_evalue, NULL);
        ///////////////////////////////////////////////////

        if (g_logEnabled)
            LOG << g_logMsg << "reduce() ends: " <<  MPI_Wtime() - reduceTime << endl;

        if (g_timingEnabled) {
            /// Wall-clock time
            gettimeofday(&wallTimeEnd, NULL);
            double tS = wallTimeStart.tv_sec * 1000000 + (wallTimeStart.tv_usec);
            double tE = wallTimeEnd.tv_sec * 1000000  + (wallTimeEnd.tv_usec);

            /// process time
            getrusage(RUSAGE_SELF, &ruTotal);
            userTimeEnd = ruTotal.ru_utime;
            sysTimeEnd = ruTotal.ru_stime;
            double tS_user = userTimeStart.tv_sec * 1000000
                             + (userTimeStart.tv_usec);
            double tE_user = userTimeEnd.tv_sec * 1000000
                             + (userTimeEnd.tv_usec);
            double tS_sys  = sysTimeStart.tv_sec * 1000000
                             + (sysTimeStart.tv_usec);
            double tE_sys  = sysTimeEnd.tv_sec * 1000000
                             + (sysTimeEnd.tv_usec);

            double profile_etime = MPI_Wtime();
            LOG << "mr-mpi-blast ends,"
                << profile_etime << ","
                << wallTimeEnd.tv_sec * 1000000 + wallTimeEnd.tv_usec << ","
                << userTimeEnd.tv_sec * 1000000 + userTimeEnd.tv_usec << ","
                << sysTimeEnd.tv_sec * 1000000 + sysTimeEnd.tv_usec << endl;

            /// MP_Wtime
            profileTime = profile_etime - profileTime;
            LOG << "Total MPI_Wtime time," << profileTime << endl;
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
    }
    
    delete pMr; /// MR-MPI
}


/** MR-MPI Map function - settup NCBI C++ Toolkit env and call blast
 * @param itask
 * @param kv
 * @param ptr
 */
 
void mr_run_blast(int itask,
                  KeyValue *kv,
                  void *ptr)
{
    structToPass_t *toPass = (structToPass_t *) ptr;
    int rank = toPass->rank;
    int iter = toPass->iter;
    
    uint32_t dbno, qBlockStart, qBlockEnd;
    
    if (g_numIter > 1) {
        dbno = g_vecWorkItem[itask + nWorkItemsPerIter * iter].dbNo;
        qBlockStart = g_vecWorkItem[itask + nWorkItemsPerIter * iter].blockBegin;
        qBlockEnd = g_vecWorkItem[itask + nWorkItemsPerIter * iter].blockEnd;
    }
    else {
        dbno = g_vecWorkItem[itask].dbNo;
        qBlockStart = g_vecWorkItem[itask].blockBegin;
        qBlockEnd = g_vecWorkItem[itask].blockEnd;
    }

    struct timeval blastcallStartTime;
    struct timeval blastcallEndTime;
    struct timeval blastcallStart_u_Time;
    struct timeval blastcallEnd_u_Time;
    struct timeval blastcallStart_s_Time;
    struct timeval blastcallEnd_s_Time;
    struct rusage  ru_blastcall;

    struct timeval qBuildStartTime;
    struct timeval qBuildStart_u_Time;
    struct timeval qBuildStart_s_Time;
    struct rusage  ru_qBuild;
    
    struct timeval dbLoadingStartTime;
    struct timeval dbLoadingStart_u_Time;
    struct timeval dbLoadingStart_s_Time;
    struct rusage  ru_dbLoading;

    ///
    /// Load Blast opotions from file
    ///
    ifstream strategyFile(g_strategyFileName.c_str(), ios::in);
    if (!strategyFile.is_open()) {
        cerr << "ERROR: failed to open a search strategy file" << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    CNcbiIstream *pStrategyFileStream = &strategyFile;
    CRef<CBlast4_request> pBl4Req;
    try {
        pBl4Req = ExtractBlast4Request(*pStrategyFileStream);
    }
    catch (const CSerialException&) {
        NCBI_THROW(CInputException, eInvalidInput, "Failed to read search strategy file");
    }
    CImportStrategy importStrategy(pBl4Req);
    CRef<blast::CBlastOptionsHandle> pBlOpts = importStrategy.GetOptionsHandle();
    pBlOpts->Validate();

    if (rank == 1 && g_optDumpEnabled == 1) {
        g_optDumpEnabled = 0;
        ofstream strategyOutFile("search_strategy.txt", ios::out);
        pBlOpts->GetOptions().DebugDumpText(strategyOutFile, "pBlOpts", 1);
        strategyOutFile.close();
    }

    string dbFileName = g_vecDbFile[dbno];

    ///
    /// Read a block of sequeces from qBlockStart
    ///
    double query_build_time;
    if (g_timingEnabled) {
        g_mapCallNo++;
        query_build_time = MPI_Wtime();
        gettimeofday(&qBuildStartTime, NULL); /// Wall-clock time
        getrusage(RUSAGE_SELF, &ru_qBuild);   /// Process time
        qBuildStart_u_Time = ru_qBuild.ru_utime;
        qBuildStart_s_Time = ru_qBuild.ru_stime;
        LOG << g_logMsg << "query build starts,"
            << query_build_time << ","
            << qBuildStartTime.tv_sec * 1000000
             + qBuildStartTime.tv_usec << ","
            << qBuildStart_u_Time.tv_sec * 1000000
             + qBuildStart_u_Time.tv_usec << ","
            << qBuildStart_s_Time.tv_sec * 1000000
             + qBuildStart_s_Time.tv_usec << ","
            << dbFileName << "," << g_mapCallNo << "," << g_MPI_procName
            << "," << qBlockStart << endl;
        LOG.flush();
    }

    const char *pMmapQueryFile = (char*)g_memmapQueryFile.data();
    ////////////////////////////////////////////////////////////////////
    string query(pMmapQueryFile + qBlockStart, qBlockEnd - qBlockStart);
    ////////////////////////////////////////////////////////////////////
    vector<string> vecHeader;        /// For collecting def lines of queries
    vector<uint32_t> vecBeginOffset; /// For collecting beginoffsets of queries
    uint32_t numQueries = 0;
    char *c = (char*)(pMmapQueryFile + qBlockStart);
    uint32_t loc = qBlockStart;
    char buff2[MAXSTR];
    size_t buffIdx = 0;

    while (loc < qBlockEnd) {
        if ((*c) == '>') {
            numQueries++;
            vecBeginOffset.push_back(loc);
            buffIdx = 0;
        }
        if (g_exclEnabled && (*c) != '\n' && buffIdx < (unsigned)MAXSTR)
            buff2[buffIdx++] = (*c);
        else if (g_exclEnabled && (*c) == '\n')
            vecHeader.push_back(string(buff2));
        c++;
        loc++;
    }

    if (g_logEnabled)
        LOG << g_logMsg << "Number of queries per a Blast call = " << numQueries << endl;

    ///
    /// Set queries as fasta input
    ///
    SDataLoaderConfig dataLoaderConf(g_bIsProtein);
    dataLoaderConf.OptimizeForWholeLargeSequenceRetrieval();
    CBlastInputSourceConfig blInputSourceConf(dataLoaderConf);
    /// Assign local query ID like "Query_1" with each query 
    //blInputSourceConf.SetQueryLocalIdMode();

    CBlastFastaInputSource blFastaInputSource(query, blInputSourceConf);
    CBlastInput blInput(&blFastaInputSource);
    CRef<CObjectManager> pObjmgr = CObjectManager::GetInstance();
    if (!pObjmgr) {
        throw runtime_error("Could not initialize object manager");
    }
    CScope scope(*pObjmgr);
    TSeqLocVector vecQueryLoc = blInput.GetAllSeqLocs(scope);
    CRef<IQueryFactory> pQueryFactory(new CObjMgr_QueryFactory(vecQueryLoc));

    ///
    /// Target db name setting
    ///    
    double db_loading_time;
    if (g_timingEnabled) {
        g_mapCallNo++;
        db_loading_time = MPI_Wtime();
        gettimeofday(&dbLoadingStartTime, NULL); /// Wall-clock time
        getrusage(RUSAGE_SELF, &ru_dbLoading);   /// Process time
        dbLoadingStart_u_Time = ru_dbLoading.ru_utime;
        dbLoadingStart_s_Time = ru_dbLoading.ru_stime;
        LOG << g_logMsg << "db_loading starts,"
            << db_loading_time << ","
            << dbLoadingStartTime.tv_sec * 1000000
             + dbLoadingStartTime.tv_usec << ","
            << dbLoadingStart_u_Time.tv_sec * 1000000
             + dbLoadingStart_u_Time.tv_usec << ","
            << dbLoadingStart_s_Time.tv_sec * 1000000
             + dbLoadingStart_s_Time.tv_usec << ","
            << dbFileName << "," << g_mapCallNo << "," << g_MPI_procName
            << "," << qBlockStart << endl;
        LOG.flush();
    }
    
    if (g_pTargetDb == 0 || dbFileName != g_prevDbName) {
        delete g_pTargetDb;        
        if (g_bIsProtein) 
            g_pTargetDb = new CSearchDatabase(dbFileName,
                                              CSearchDatabase::eBlastDbIsProtein);
        else
            g_pTargetDb = new CSearchDatabase(dbFileName,
                                              CSearchDatabase::eBlastDbIsNucleotide);
    }
    g_prevDbName = dbFileName;

    ///
    /// Use the CLocalBlast class to run a BLAST search 
    ///
    double blast_call_time;
    if (g_timingEnabled) {
        blast_call_time = MPI_Wtime();
        gettimeofday(&blastcallStartTime, NULL); /// Wall-clock time
        getrusage(RUSAGE_SELF, &ru_blastcall);   /// Process time
        blastcallStart_u_Time = ru_blastcall.ru_utime;
        blastcallStart_s_Time = ru_blastcall.ru_stime;
        LOG << g_logMsg << "blast call starts,"
            << blast_call_time << ","
            << blastcallStartTime.tv_sec * 1000000
             + blastcallStartTime.tv_usec << ","
            << blastcallStart_u_Time.tv_sec * 1000000
             + blastcallStart_u_Time.tv_usec << ","
            << blastcallStart_s_Time.tv_sec * 1000000
             + blastcallStart_s_Time.tv_usec << ","
            << dbFileName << "," << g_mapCallNo << "," << g_MPI_procName
            << "," << qBlockStart << endl;
        LOG.flush();
    }

    ////////////////////////////////////////////////////////////
    CSearchResultSet results;
    CLocalBlast lcl_blast(pQueryFactory, pBlOpts, *g_pTargetDb);
    results = *lcl_blast.Run(); 
    ////////////////////////////////////////////////////////////

    if (g_timingEnabled) {
        double blast_call_etime = MPI_Wtime();
        gettimeofday(&blastcallEndTime, NULL); /// Wall-clock time
        getrusage(RUSAGE_SELF, &ru_blastcall); /// process time
        blastcallEnd_u_Time = ru_blastcall.ru_utime;
        blastcallEnd_s_Time = ru_blastcall.ru_stime;
        LOG << g_logMsg << "blast call ends,"
            << blast_call_etime << ","
            << blastcallEndTime.tv_sec * 1000000
             + blastcallEndTime.tv_usec << ","
            << blastcallEnd_u_Time.tv_sec * 1000000
             + blastcallEnd_u_Time.tv_usec << ","
            << blastcallEnd_s_Time.tv_sec * 1000000
             + blastcallEnd_s_Time.tv_usec << ","
            << dbFileName << "," << g_mapCallNo << "," << g_MPI_procName
            << "," << qBlockStart << endl;
        LOG.flush();
    }

    ///
    /// Get warning messages
    ///
    for (size_t i = 0; i < results.GetNumResults(); ++i) {
        TQueryMessages messages = results[i].GetErrors(eBlastSevWarning);
        if (messages.size() > 0) {
            CConstRef<CSeq_id> pSeq_id = results[i].GetSeqId();
            if (pSeq_id.NotEmpty())
                cerr << "ID: " << pSeq_id->AsFastaString() << endl;
            else
                cerr << "ID: " << "Unknown" << endl;

            ITERATE(vector<CRef<CSearchMessage> >, it, messages) {
                cerr << (*it)->GetMessage() << endl;
            }
        }
    }

    ///
    /// Get the results and add to kv
    ///
    for (size_t i = 0; i < results.GetNumResults(); ++i) {

        CConstRef<CSeq_align_set> pAln_set = results[i].GetSeqAlign();
        if (results[i].HasAlignments()) {
            ITERATE(CSeq_align_set::Tdata, itr_res, pAln_set->Get()) {
                const CSeq_align& seqAlign = **itr_res;

                /// Note: queryID is not unique. It's internal qid in 
                /// blFastaInputSource.
                string queryID = seqAlign.GetSeq_id(QUERY).GetSeqIdString();
                string subID = seqAlign.GetSeq_id(SUBJECT).GetSeqIdString();

                ///
                /// Refer: 
                /// http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/doxyhtml/classCSeq__align.html
                ///
                /// eScore_Score: generic score, definable by any algorithm; not comparable across algorithms
                /// eScore_Blast
                /// eScore_BitScore: BLAST-specific bit score
                /// eScore_EValue: BLAST-specific e-value
                /// eScore_AlignLength: not a score per se, but a useful metric nonetheless. This is the sum of all aligned segments and all gaps; this excludes introns and discontinuities
                /// eScore_IdentityCount
                /// eScore_PositiveCount
                /// eScore_NegativeCount
                /// eScore_MismatchCount
                /// eScore_PercentIdentity_Gapped
                /// eScore_PercentIdentity_Ungapped
                /// eScore_PercentIdentity_GapOpeningOnly
                /// eScore_PercentCoverage
                /// eScore_SumEValue
                /// eScore_CompAdjMethod
                /// eScore_HighQualityPercentCoverage ?????
                /// eScore_PercentIdentity = eScore_PercentIdentity_Gapped 
                ///
                double pIdentity = 0.0;
                //double PercentIdentity_Gapped;
                //double PercentIdentity_Ungapped;
                //double PercentIdentity_GapOpeningOnly;
                //double percentCoverage = 0.0;
                //double HighQualityPercentCoverage;
                
                seqAlign.GetNamedScore(CSeq_align::eScore_PercentIdentity, pIdentity);
                //seqAlign.GetNamedScore(CSeq_align::eScore_PercentIdentity_Gapped, PercentIdentity_Gapped);
                //seqAlign.GetNamedScore(CSeq_align::eScore_PercentIdentity_Ungapped, PercentIdentity_Ungapped);
                //seqAlign.GetNamedScore(CSeq_align::eScore_PercentIdentity_GapOpeningOnly, PercentIdentity_GapOpeningOnly);
                //seqAlign.GetNamedScore(CSeq_align::eScore_PercentCoverage, percentCoverage);
                //seqAlign.GetNamedScore(CSeq_align::eScore_HighQualityPercentCoverage, HighQualityPercentCoverage);
                //cout << "pIdentity = " << double(pIdentity) << endl;
                //cout << "PercentIdentity_Gapped = " << PercentIdentity_Gapped << endl;
                //cout << "PercentIdentity_Ungapped = " << PercentIdentity_Ungapped << endl;
                //cout << "PercentIdentity_GapOpeningOnly = " << PercentIdentity_GapOpeningOnly << endl;
                //cout << "PercentCoverage = " << double(percentCoverage) << endl;
                //cout << "HighQualityPercentCoverage = " << HighQualityPercentCoverage << endl;

                uint32_t gapOpens    = seqAlign.GetNumGapOpenings();
                //cout << "gapOpens = " << gapOpens << endl;                
                //uint32_t totalGaps   = seqAlign.GetTotalGapCount();
                //cout << "totalGaps = " << totalGaps << endl;                
                uint32_t qStart      = seqAlign.GetSeqStart(QUERY);
                uint32_t qEnd        = seqAlign.GetSeqStop(QUERY);
                uint32_t sStart      = seqAlign.GetSeqStart(SUBJECT);
                uint32_t sEnd        = seqAlign.GetSeqStop(SUBJECT);
                uint32_t alignLen    = seqAlign.GetAlignLength();
                //cout << "alignLen = " << alignLen << endl;                

                double eValue = 0.0;
                int bitScore = 0, misMatches = 0;
                seqAlign.GetNamedScore(CSeq_align::eScore_EValue, eValue);
                seqAlign.GetNamedScore(CSeq_align::eScore_BitScore, bitScore);
                seqAlign.GetNamedScore(CSeq_align::eScore_MismatchCount, misMatches);
                //cout << "misMatches = " << misMatches << endl;  
                
                //double alignLenRatio = seqAlign.AlignLengthRatio();
                //TSeqRange r0 = seqAlign.GetSeqRange(QUERY);
                //TSeqRange r1 = seqAlign.GetSeqRange(SUBJECT);
                //double r = 0;
                //if (r0.GetLength()) {
                    //r = double(r1.GetLength()) / double(r0.GetLength());
                //}
                //cout << "alignLenRatio = " << r << " "<< r1.GetLength() << " / " << r0.GetLength() << endl;
                
                int identityCount = 0;
                //int positiveCount = 0, negativeCount = 0;
                seqAlign.GetNamedScore(CSeq_align::eScore_IdentityCount, identityCount);
                //seqAlign.GetNamedScore(CSeq_align::eScore_PositiveCount, positiveCount);
                //seqAlign.GetNamedScore(CSeq_align::eScore_NegativeCount, negativeCount);
                //cout << "identityCount = " << identityCount << endl;
                //cout << "positiveCount = " << positiveCount << endl;
                //cout << "negativeCount = " << negativeCount << endl;
                
                ///
                /// Retrieve GI and cutting locations    
                ///
                
                /// Get the saved qBlockStart from vBeginOffset
                uint32_t bOffset = vecBeginOffset[boost::lexical_cast<uint32_t>(queryID)-1];
                
                /// Get GI
                string line(pMmapQueryFile + bOffset, 80);
                vector<string> vecTokens;
                boost::split(vecTokens, line, boost::is_any_of("\n"));
                assert(vecTokens.size() == 2);
                string defLine = vecTokens[0];
                vecTokens.clear();
                boost::split(vecTokens, defLine, boost::is_any_of("|"));
                uint32_t gi = boost::lexical_cast<uint32_t>(vecTokens[1]);
                
                /// Get cutting coords               
                string coord = vecTokens[vecTokens.size()-1];
                vecTokens.clear();
                boost::split(vecTokens, coord, boost::is_any_of("_"));
                int cutStart = boost::lexical_cast<int>(vecTokens[vecTokens.size()-2]);
                int cutEnd   = boost::lexical_cast<int>(vecTokens[vecTokens.size()-1]);
                int length   = cutEnd - cutStart;
   
                ///
                /// Add a csv blast result to kv
                ///
                if (!g_exclEnabled) {
                    
                    ///
                    /// Doug's filtering
                    ///
                    /// identity = # of identical bases / length of read
                    /// coverage = (read end – read begin)/ length of read
                    ///
                    /// identity > 50%
                    /// coverage > 90%
                    ///
                    double ident = double(identityCount) / double(length);
                    double cover = double(qEnd - qStart) / double(length);
                    //cout << "length = " << length << endl;
                    //cout << "ident = " << ident << endl;
                    //cout << "cover = " << cover << endl;
 
                    if (ident >= g_IDENT_CUTOFF && cover >= g_COVER_CUTOFF) {
                        structBlRes_t res;
                        res.subjectId[0] = '\0';
                        strcpy(res.subjectId, subID.c_str());
                        res.identity      = pIdentity;
                        res.alignLen      = alignLen;
                        res.misMatches    = misMatches;
                        res.gapOpens      = gapOpens;
                        res.qStart        = qStart + cutStart;
                        res.qEnd          = qEnd + cutStart;
                        res.sStart        = sStart;
                        res.sEnd          = sEnd;
                        res.evalue        = eValue;
                        res.bitScore      = bitScore;
                        res.cutStart      = cutStart;
                        res.cutEnd        = cutEnd;
                        res.doug_identity = ident;
                        res.doug_coverage = cover;
                        
                        uint32_t newKey   = gi;
                        kv->add((char*)&newKey, sizeof(uint32_t), (char*)&res,
                                sizeof(structBlRes_t));     
                    }               
                }
                else {
                    ///
                    /// Tokenize query def line
                    ///
                    string qHeader =
                        vecHeader[boost::lexical_cast<uint32_t>(queryID)-1];
                    vector<string> vecQueryId;
                    boost::split(vecQueryId, qHeader, boost::is_any_of("|"));

                    ///
                    /// The below infomation is only for our own simulated
                    /// sequence data sets which has
                    /// - origin GI
                    /// - unique query ID
                    /// - original sequence length
                    /// - cut location start
                    /// - cut location end
                    ///
                    /// GI
                    string qGi = vecQueryId[1];   
                    /// length of the orig seq
                    //uint32_t origLen = boost::lexical_cast<uint32_t>(vecQueryId[3]);
                    /// cut coordinates - start
                    int qCutLocStart = boost::lexical_cast<int>(vecQueryId[4]);
                    /// cut coordinates - end
                    int qCutLocEnd =   boost::lexical_cast<int>(vecQueryId[5]);

                    if (!check_exclusion(qGi, subID, qCutLocStart, qCutLocEnd,
                                         sStart, sEnd, g_exclThreshold)) {

                        ///
                        /// To pass Blast hits using struct, outfmt=6
                        /// query id, subject id, % identity, alignment length,
                        /// misMatches, gap opens, q. start, q. end, s. start,
                        /// s. end, evalue, bit score
                        ///
                        structBlRes_t res;
                        res.subjectId[0] = '\0';
                        strcpy(res.subjectId, subID.c_str());
                        res.identity   = pIdentity;
                        res.alignLen   = alignLen;
                        res.misMatches = misMatches;
                        res.gapOpens   = gapOpens;
                        res.qStart     = qStart;
                        res.qEnd       = qEnd;
                        res.sStart     = sStart;
                        res.sEnd       = sEnd;
                        res.evalue     = eValue;
                        res.bitScore   = bitScore;
                        
                        ///
                        /// ADD <KEY = "QUERYID", VALUE="BLASTRESULT">
                        /// TO KV
                        ///
                        uint32_t newKey = boost::lexical_cast<uint32_t>(qGi);
                        kv->add((char*)&newKey, sizeof(uint32_t), (char*)&res,
                                sizeof(structBlRes_t));
                    }
                    ///
                    /// Found a self hits. Record the hits in a histroy file
                    ///
                    else {
                        string exFileName = g_outFilePrefix + "-"
                                            + boost::lexical_cast<string>(rank) 
                                            + "-" + g_exclHistFileName;
                        ofstream exFile(exFileName.c_str(), ios::out | ios::app);

                        if (!exFile) {
                            cerr << "ERROR: failed to open a exclusion "
                                 << "history file" << endl;
                            MPI_Abort(MPI_COMM_WORLD, 1);
                        }
                        else {
                            /// Format: qGi sGi qCutLocStart qCutLocEnd sStart sEnd
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
            }
        }
    }
    /// end of processing hits       
}


/** Sort function - Passed to MR-MPI sort_values() for sorting blast result
 * string by bit score.
 * @param e1
 * @param e2
 */

inline bool evalue_compare(structEValue_t e1,
                           structEValue_t e2)
{
    int ret = strcmp(e1.pSubjectId, e2.pSubjectId);
    /// ret > 0:  e1.pSubjectId > e2.pSubjectId
    /// ret == 0: e1.pSubjectId == e2.pSubjectId
    /// ret < 0: e1.pSubjectId < e2.pSubjectId
    if (ret == 0)
        return (e1.evalue != e2.evalue) ? 
            (e1.evalue < e2.evalue) : (e1.bitScore < e2.bitScore);
    else if (ret > 0) return true;
    else return false;
}

/** Sort by evalue - Passed to MR-MPI reduce() for sorting KMVs by evalue.
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
    /// Check if there is KMV overflow (which is not a failure)
    //assert(multivalue != NULL && nvalues != 0);
    //if (multivalue != NULL && nvalues != 0) {
        //cout << g_logMsg << "KMV overflow detected" << endl;
        //cout.flush();;
    //}

    ///
    /// Make structEValue_t = {structBlRes_t* pRec; double evalue;}
    /// and sort by evalue
    ///
    vector<structEValue_t> vecHit;
    for (size_t n = 0; n < (unsigned)nvalues; n++) {
        structBlRes_t* res = (structBlRes_t*)multivalue;
        structEValue_t structEvalue;
        structEvalue.pRec = res;
        structEvalue.pSubjectId = res->subjectId;
        structEvalue.evalue = res->evalue;
        structEvalue.bitScore = res->bitScore;
        vecHit.push_back(structEvalue);
        multivalue += sizeof(structBlRes_t);
    }
    
    /// ////////////////////////////////////////////////
    sort(vecHit.begin(), vecHit.end(), evalue_compare);
    /// ////////////////////////////////////////////////
    
    ///
    /// After sorting, each worker saves the set of results into a file.
    /// Note: The file open option is "a".
    /// Note: outfmt = 6 in Blast
    ///
    ofstream outputFile(g_hitFileName.c_str(), ios::out | ios::app);
    for (size_t n = 0; n < (unsigned)nvalues; n++) {
        structBlRes_t* res = (structBlRes_t*)(vecHit[n].pRec);
        outputFile << *(uint32_t*)key << "\t"
            << res->subjectId << "\t"
            << double(res->identity) << "\t"
            << res->alignLen << "\t"
            << res->misMatches << "\t"
            << res->gapOpens << "\t"
            << res->qStart << "\t"
            << res->qEnd << "\t"
            << res->sStart << "\t"
            << res->sEnd << "\t"
            << double(res->evalue) << "\t"
            << res->bitScore << "\t"
            << res->cutStart << "\t"
            << res->cutEnd << "\t"
            << res->doug_identity << "\t"
            << res->doug_coverage
            << endl;
    }
    outputFile.close();
    vecHit.clear();
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
    bool bRet = false;

    if (qGi == sGi) {
        if (qCutLocStart < 0) {
            ///
            /// In >gi|222299657|18|3605|-400|3604
            /// -400|3604 means query[-400:3604] in Python.
            ///
            qCutLocStart = qCutLocEnd + 1 - qCutLocStart;
            qCutLocEnd += 1;
        }
        if (
           (qCutLocStart - threshold <= sStart && sStart <= qCutLocStart + threshold)
           && 
           (qCutLocEnd - threshold <= sEnd && sEnd <= qCutLocEnd + threshold)
           ) 
        {
            bRet = true;
        }
    }

    return bRet;
}

/** Collect MPI node names and ranks
 * @param rank
 * @param numProcs
 * @param mpiComm
 */
 
void collect_mpi_node_name(int rank, int numProcs, MPI_Comm mpiComm)
{
    MPI_Status MPI_status;
    char procName[MAXSTR];
    int rankNo;
    int tag = 12345;
    if (rank != MPI_UNDEFINED) { /// if not -32766
        if (rank == 0) {
            for (size_t src = 1; src < (unsigned)numProcs; src++) {
                MPI_Recv(&rankNo, 1, MPI_INT, src, tag, mpiComm, &MPI_status);
                MPI_Recv(&procName, MAXSTR, MPI_CHAR, src, tag, mpiComm, &MPI_status);
                g_multimapProcNameRank.insert(pair<string, int>(procName, rankNo));
                g_mapRankProcName.insert(pair<int, string>(rankNo, procName));
            }
        }
        else {
            MPI_Send(&g_MPI_worldRank, 1, MPI_INT, 0, tag, mpiComm);
            MPI_Send(&g_MPI_procName, MAXSTR, MPI_CHAR, 0, tag, mpiComm);
        }
    }
}


 

/// EOF

