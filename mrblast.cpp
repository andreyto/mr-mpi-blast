////////////////////////////////////////////////////////////////////////////////
//
//  MR-MPI-BLAST: Parallelizing BLAST with MR-MPI
//
//  Author: Seung-Jin Sul
//         (ssul@jcvi.org)
//
//  Last updated: 06/23/2011
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
#include <algo/blast/api/blast_nucl_options.hpp>
#include <algo/blast/api/blast_prot_options.hpp>

/// Import search strategy
#include <objects/blast/Blast4_request.hpp>
#include <algo/blast/api/search_strategy.hpp>

/// Tabular
#include <algo/blast/format/blast_format.hpp>
#include <algo/blast/blastinput/blastn_args.hpp>
#include <algo/blast/blastinput/blastp_args.hpp>

/// sequence::GetTitle() for get def line
#include <objmgr/util/sequence.hpp>

USING_NCBI_SCOPE;
USING_SCOPE (blast);
//USING_SCOPE (objects);
USING_SCOPE (align_format);
    
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
uintmax_t g_realFileSize = 0;

/// For align info
#include "blast_app_util.hpp"

/// 
/// Settings from mrblast.ini conf file
/// 

/// MR-MPI options
int g_verbosity;
int g_timer;                /// log elapsed time for each mapreduce call
int g_memSize;              /// page size (in Mbytes)
int g_outOfCore;

/// DB options
string g_dbFileName;
string g_configFileName = "mrblast.ini";
int g_numDbFiles;
const int MAXSTR = 80;      /// For mpi proc name and query header

/// Filtering
float g_IDENT_CUTOFF = 0.0; /// Doug's identity for filtering
float g_COVER_CUTOFF = 0.0; /// Doug's coverage for filtering

/// 
/// Log
/// 
#include <sys/time.h>
#include <sys/resource.h>
//#define NDEBUG 1
int g_bLogEnabled = 0;
int g_timingEnabled = 0;
int g_optDumpEnabled = 0;   /// For dumping Blast opts_hndl out
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
    uintmax_t qStart;
    uint32_t qLength;
} structQueryIndex_t;

typedef struct structWorkItem {
    uint32_t dbNo;
    uintmax_t blockBegin;
    uintmax_t blockEnd;
} structWorkItem_t;

vector<string> g_vecDbFile;
vector<uintmax_t> g_vecBlockBeginLoc;
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
CRef<CSearchDatabase> g_searchDatabase;
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
    uint32_t    subjectId;    
    //float       identity;
    //uint32_t    alignLen;
    //int         misMatches;
    //uint32_t    gapOpens;
    uint32_t    qStart;
    uint32_t    qEnd;
    uint32_t    sStart;
    uint32_t    sEnd;
    double      eValue;
    float       bitScore;
    uint32_t    upperStart;
    uint32_t    upperEnd;
    float       doug_identity;
    float       doug_coverage;
} structBlRes_t;

/// For saving hits in ban format
typedef struct structBlRes2 {
    uint32_t    qId;
    uint32_t    subjectId;
    uint32_t    qStart;
    uint32_t    qEnd;
    uint32_t    sStart;
    uint32_t    sEnd;
    double      eValue;
    float       bitScore;
    uint32_t    upperStart;
    uint32_t    upperEnd;
    float       doug_identity;
    float       doug_coverage;
} structBlRes2_t;

/// To sort Blast hits by evalue
typedef struct structEvalue {
    structBlRes_t *pRec;
    uint32_t       subjectId;
    double         eValue;
    float          bitScore;
} structEValue_t;

/// Multiple iterations
int g_numIter = 1;   
uint32_t nSubWorkItemSets = 0;
uint32_t nWorkItemsPerIter = 0;
uint32_t nRemains = 0;

/// To pass info to map() and reduce()
typedef struct structToPass {
    int rank;
    int iter;
} structToPass_t;

/// HSP tiling
typedef struct subrange {
    uint32_t ss, se;
} subRange_t;
        
/// For NcbiApp
CRef<CBlastOptionsHandle> g_opts_hndl;
CRef<CBlastAppArgs> g_cmdLineArgs;

///
/// Function declarations
///
void        run_mr_mpi_blast(MPI_Comm mpiComm, int rank); 
void        mr_run_blast(int itask, KeyValue *kv, void *ptr);
inline void mpi_collect_node_name(int rank, int numProcs, MPI_Comm mpiComm);
inline void mr_sort_multivalues_by_evalue(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
inline bool compare_evalue(structEValue_t e1, structEValue_t e2);


class CMrMpiBlastApplication : public CNcbiApplication
{
private:
    virtual void Init(void);
    virtual int  Run(void); 
};

void CMrMpiBlastApplication::Init(void)
{
    ///
    /// Read conf file, mrblast.ini and set parameters
    ///
    ifstream config(g_configFileName.c_str(), ios::in);
    if (!config) {
        cerr << "ERROR: configuration file not found" << endl;
        exit(1);
    }

    set<string> options;
    map<string, string> parameters;
    options.insert("*");

    try {
        for (pod::config_file_iterator i(config, options), e ; i != e; ++i) {
            parameters[i->string_key] = i->value[0];
        }
        try {
            g_verbosity      = boost::lexical_cast<int>(parameters["VERBOSITY"]);
            g_timer          = boost::lexical_cast<int>(parameters["TIMER"]);
            g_memSize        = boost::lexical_cast<int>(parameters["MEMSIZE"]);
            g_outOfCore      = boost::lexical_cast<int>(parameters["OUTOFCORE"]);

            g_bLogEnabled    = boost::lexical_cast<bool>(parameters["LOGENABLED"]);
            g_timingEnabled  = boost::lexical_cast<int>(parameters["TIMING"]);
            g_optDumpEnabled = boost::lexical_cast<int>(parameters["OPTDUMP"]);
            g_logFileName    = parameters["LOGFNAME"];
            
            g_IDENT_CUTOFF   = boost::lexical_cast<float>(parameters["IDENTCUTOFF"]);
            g_COVER_CUTOFF   = boost::lexical_cast<float>(parameters["COVERCUTOFF"]);
            
            g_queryFileName  = parameters["QUERYFILENAME"];
            g_indexFileName  = parameters["INDEXFILENAME"];
            g_dbFileName     = parameters["DBLISTFILENAME"];
            g_outFilePrefix  = parameters["OUTFILEPREFIX"];
            g_blockSize      = boost::lexical_cast<int>(parameters["BLOCKSIZE"]);
            g_mapStyle       = boost::lexical_cast<int>(parameters["MAPSTYLE"]);
            g_numIter        = boost::lexical_cast<int>(parameters["NUMITER"]);
            g_bIsProtein     = boost::lexical_cast<bool>(parameters["ISPROTEIN"]);
        }
        catch (const boost::bad_lexical_cast &) {
            cerr << "Exception: bad_lexical_cast" << endl;
            cerr.flush();
        }
    }
    catch (exception& e) {
        cerr << "Exception: " << e.what() << endl;
        cerr.flush();
    }
     
    if (g_bIsProtein) g_cmdLineArgs.Reset(new CBlastpAppArgs());
    else g_cmdLineArgs.Reset(new CBlastnAppArgs());
    HideStdArgs(fHideLogfile | fHideConffile | fHideFullVersion | 
                fHideXmlHelp | fHideDryRun);
    SetupArgDescriptions(g_cmdLineArgs->SetCommandLine());  
}

int CMrMpiBlastApplication::Run(void) 
{     
    /// Allow the fasta reader to complain on invalid sequence input
    SetDiagPostLevel(eDiag_Warning);
    
    const CArgs& args = GetArgs();
    string allArgs;
    RecoverSearchStrategy(args, g_cmdLineArgs);
    CRef<CBlastOptionsHandle> opts_hndl(&*g_cmdLineArgs->SetOptions(args));
    g_opts_hndl = opts_hndl;
    
    ///
    /// Collect MPI node names and rank number for mapstyle=3 scheduler
    ///
    mpi_collect_node_name(g_MPI_worldRank, g_MPI_numProcs, MPI_COMM_WORLD);
    
    if (g_MPI_worldRank == 0) {
        for ( multimap<string, int>::const_iterator iter = g_multimapProcNameRank.begin();
              iter != g_multimapProcNameRank.end(); ++iter )
            cout << iter->first << '\t' << iter->second << '\n';
    }
   
    ///
    /// Creat memory-mapped file for query
    ///
    g_realFileSize = boost::filesystem::file_size(g_queryFileName);
    g_memmapQueryFile.open(g_queryFileName, g_realFileSize, 0);
    if (!g_memmapQueryFile.is_open()) {
        MPI_Finalize();
        cerr << "ERROR: failed to create mmap query file\n";
        exit(2);
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
        MPI_Finalize();
        cerr << "ERROR: DB list file open error.\n";
        exit(3);
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
        MPI_Finalize();
        cerr << "ERROR: index file open error.\n";
        exit(4);
    }
    while (!getline(indexFile, line).eof() && line.length() > 0) {
        vector<string> tok;
        boost::split(tok, line, boost::is_any_of(","));
        assert(tok.size() == 2);
        structQueryIndex_t queryIndex;
        queryIndex.qStart  = boost::lexical_cast<uintmax_t>(tok[0]);
        queryIndex.qLength = boost::lexical_cast<uint32_t>(tok[1]);
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
            structWorkItem_t aWorkItem;
            aWorkItem.blockBegin = g_vecBlockBeginLoc[b];
            aWorkItem.blockEnd   = g_vecBlockBeginLoc[b + 1] - 1;
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
    
    ///
    /// Save the list of workitem just in case
    ///
    if (g_MPI_worldRank == 0) {
        string workItemFName = g_outFilePrefix + "-workitems.txt";
        ofstream workItemFile(workItemFName.c_str(), ios::out);
        for (size_t i = 0; i < g_numWorkItems; i++)
            workItemFile << g_vecWorkItem[i].blockBegin << ","
                         << g_vecWorkItem[i].blockEnd << ","
                         << g_vecWorkItem[i].dbNo << endl;
        workItemFile.close();
    }        
    
    if (g_MPI_worldRank == 0) {
        cout << "Input query file = " << g_queryFileName << endl;
        cout << "Input query index file = " << g_indexFileName << endl;
        cout << "Database list file = " << g_dbFileName << endl;
        cout << "Block size (bp) = " << g_blockSize << endl;
        cout << "Number of query blocks = " << g_numQueryBlocks << endl;
        cout << "Number of DB files = " << g_numDbFiles << endl;
        cout << "Number of total work items = " << g_numWorkItems << endl;
        cout << "Number of iterations = " << g_numIter << endl;
        cout << "Identity cutoff = " << g_IDENT_CUTOFF << endl;
        cout << "Coverage cutoff = " << g_COVER_CUTOFF << endl;        
        cout.flush();
    }

    ///
    /// Calculate sub work item size for multiple iterations
    ///
    if (g_numIter != 1) {
        nSubWorkItemSets = g_numIter;
        nWorkItemsPerIter = (g_numQueryBlocks / g_numIter) * g_numDbFiles;            
        nRemains = (g_numQueryBlocks % g_numIter) * g_numDbFiles;
    }
    else { 
        nWorkItemsPerIter = g_numWorkItems;
        nRemains = 0;
        nSubWorkItemSets = 1;
    }
    /// for one more added file of remaining work items
    if (nRemains) nSubWorkItemSets++; 

    if (g_MPI_worldRank == 0 && g_numIter != 1) {        
        cout << "Number of sub work item sets = " << nSubWorkItemSets << endl;
        cout << "Number of work items per iteration  = " << nWorkItemsPerIter << endl;
        cout << "Number of work items remaining = " << nRemains << endl;
        cout.flush();
    }       
        
    MPI_Barrier(MPI_COMM_WORLD);

    ///
    /// Run mr-mpi calls
    ///
    run_mr_mpi_blast(MPI_COMM_WORLD, g_MPI_worldRank); 
    if (g_MPI_worldRank == 0) cout << "Done!" << endl;
    
    ///
    /// Clean up
    ///
    g_vecDbFile.clear();
    g_vecQueryIndex.clear();
    g_vecBlockBeginLoc.clear();
    g_vecWorkItem.clear();
    if (!g_searchDatabase.IsNull()) g_searchDatabase.Release();
    if (g_bLogEnabled || g_timingEnabled) g_logFileStream.close();
    if (g_memmapQueryFile.is_open()) g_memmapQueryFile.close(); /// close mmapped file  
    

    return 0;
}
 
 
int main(int argc, char** argv)
{
    ///
    /// MPI setup
    ///
    int MPI_procNameLen;    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &g_MPI_worldRank);
    MPI_Comm_size(MPI_COMM_WORLD, &g_MPI_numProcs);
    MPI_Get_processor_name(g_MPI_procName, &MPI_procNameLen);
        
    /// Execute main application function
    int ret = CMrMpiBlastApplication().AppMain(argc, argv);    
    
    MPI_Finalize();   
    
    return ret;
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
    if (g_bLogEnabled || g_timingEnabled) {
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
    /// mapstyle    = 0 (chunk) or 1 (stride) or 2 (master/slave)
    /// all2all     = 0 (irregular communication) or 1 (use MPI_Alltoallv)
    /// verbosity   = 0 (none) or 1 (summary) or 2 (histogrammed)
    /// timer       = 0 (none) or 1 (summary) or 2 (histogrammed)
    /// memsize     = N = number of Mbytes per page of memory
    /// minpage     = N = # of pages to pre-allocate per processor
    /// maxpage     = N = max # of pages allocatable per processor
    /// freepage    = 1 if memory pages are freed in between operations, 0 if held
    /// outofcore   = 1 if even 1-page data sets are forced to disk, 0 if not, -1
    ///               if cannot write to disk
    /// zeropage    = 1 if zero out every allocated page, 0 if not
    /// keyalign    = N = byte-alignment of keys
    /// valuealign  = N = byte-alignment of values
    /// fpath       = string
    ///
    MapReduce *pMr = new MapReduce(mpiComm);
    pMr->verbosity = g_verbosity;
    pMr->timer     = g_timer;
    pMr->memsize   = g_memSize;
    pMr->keyalign  = sizeof(uint32_t);  /// The key is query GI
    pMr->mapstyle  = g_mapStyle;        /// master/slave=2, custom scheduler=3
    pMr->outofcore = g_outOfCore;       /// "-1" to disable out-of-core operation
    MPI_Barrier(mpiComm);

    for (int iter = 0; iter < nSubWorkItemSets; iter++) {
        if (g_MPI_worldRank == 0) cout << "Iteration: " << iter << endl;
        if (g_bLogEnabled) LOG << g_logMsg << "Iteration: " << iter << endl;
        
        ///
        /// map, collate, reduce
        ///
        double mapTime;
        if (g_bLogEnabled) {
            mapTime = MPI_Wtime();
            LOG << g_logMsg << "map() starts: " <<  mapTime << endl;
        }
    
        structToPass_t toPass;
        toPass.rank = rank;
        toPass.iter = iter;
        
        /////////////////////////////////////////////////////////////
        if (g_numIter != 1)
            if (nRemains != 0 && iter == nSubWorkItemSets - 1)
                pMr->map(nRemains, &mr_run_blast, &toPass);
            else 
                pMr->map(nWorkItemsPerIter, &mr_run_blast, &toPass);
        else 
            pMr->map(g_numWorkItems, &mr_run_blast, &toPass);
        /////////////////////////////////////////////////////////////
                
        if (g_bLogEnabled)
            LOG << g_logMsg << "map() ends: " <<  MPI_Wtime() - mapTime << endl;

        double collateTime;
        if (g_bLogEnabled) {
            collateTime = MPI_Wtime();
            LOG << g_logMsg << "collate() starts: " << collateTime << endl;
        }

        ///////////////////
        pMr->collate(NULL);
        ///////////////////

        if (g_bLogEnabled)
            LOG << g_logMsg << "collate() ends: " << MPI_Wtime() - collateTime 
                << endl;

        double reduceTime;
        if (g_bLogEnabled) {
            reduceTime = MPI_Wtime();
            LOG << g_logMsg << "reduce() starts: " << reduceTime << endl;
        }
        
        g_hitFileName = g_outFilePrefix + "-hits-" 
                    + boost::lexical_cast<string>(iter) + "-"
                    + boost::lexical_cast<string>(rank) + ".txt";
                        
        ///////////////////////////////////////////////////
        pMr->reduce(&mr_sort_multivalues_by_evalue, NULL);
        ///////////////////////////////////////////////////

        if (g_bLogEnabled)
            LOG << g_logMsg << "reduce() ends: " <<  MPI_Wtime() - reduceTime 
                << endl;
       
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
                << tE_user - tS_user << "," << (tE_user - tS_user) / 1000000 
                << endl;
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
    if (g_bLogEnabled)
        LOG << g_logMsg << "Start processing itask = " << itask << endl;
            
    uint32_t dbno;
    uintmax_t qBlockStart, qBlockEnd;    
    if (g_numIter > 1) {
        dbno = g_vecWorkItem[itask + nWorkItemsPerIter * iter].dbNo;
        qBlockStart = g_vecWorkItem[itask + nWorkItemsPerIter * iter].blockBegin;
        qBlockEnd = g_vecWorkItem[itask + nWorkItemsPerIter * iter].blockEnd;
    }
    else {
        dbno = g_vecWorkItem[itask].dbNo;
        qBlockStart = g_vecWorkItem[itask].blockBegin;
        qBlockEnd   = g_vecWorkItem[itask].blockEnd;
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
     
    if (rank == 1 && g_optDumpEnabled == 1) {
        g_optDumpEnabled = 0;
        string strategyFName = g_outFilePrefix + "-search_strategy.txt";
        ofstream strategyOutFile(strategyFName.c_str(), ios::out);
        g_opts_hndl->GetOptions().DebugDumpText(strategyOutFile, "opts_hndl", 1);
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

    assert(g_memmapQueryFile.is_open());
    const char *pMmapQueryFile = (char*)g_memmapQueryFile.data();
    ////////////////////////////////////////////////////////////////////
    string query(pMmapQueryFile + qBlockStart, qBlockEnd - qBlockStart);    
    ////////////////////////////////////////////////////////////////////
    assert(query.length() > 0);
    assert(query[0] == '>');

    const CBlastOptions& opt = g_opts_hndl->GetOptions();
    
    /*** Get the query sequence(s) ***/
    CRef<CQueryOptionsArgs> query_opts = g_cmdLineArgs->GetQueryOptionsArgs();
    SDataLoaderConfig dlconfig(query_opts->QueryIsProtein());
    dlconfig.OptimizeForWholeLargeSequenceRetrieval();
    CBlastInputSourceConfig iconfig(dlconfig, query_opts->GetStrand(),
                                 query_opts->UseLowercaseMasks(),
                                 query_opts->GetParseDeflines(),
                                 query_opts->GetRange(),
                                 !g_cmdLineArgs->ExecuteRemotely());
    iconfig.SetLowercaseMask(true); /// Enforce lowercase mask option again
    CBlastFastaInputSource fasta(query, iconfig);
    CBlastInput input(&fasta);

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
    
    if (g_searchDatabase.IsNull() || dbFileName != g_prevDbName) {
        if (!g_searchDatabase.IsNull()) g_searchDatabase.Release();
        if (g_bIsProtein) {
            CRef<CSearchDatabase> searchDatabase(new CSearchDatabase(
                dbFileName, CSearchDatabase::eBlastDbIsProtein));
            g_searchDatabase = searchDatabase;
        }
        else {
            CRef<CSearchDatabase> searchDatabase(new CSearchDatabase(
                dbFileName, CSearchDatabase::eBlastDbIsNucleotide));
            g_searchDatabase = searchDatabase;
        }
    }
    g_prevDbName = dbFileName;
        
    /*** Initialize the database/subject ***/
    CRef<CBlastDatabaseArgs> db_args(g_cmdLineArgs->GetBlastDatabaseArgs());
    CRef<CLocalDbAdapter> db_adapter;
    CRef<CScope> scope;
    db_adapter.Reset();
    CRef<CSearchDatabase> search_db = g_searchDatabase;
    if (scope.Empty()) {
        scope.Reset(new CScope(*CObjectManager::GetInstance()));
    }
    _ASSERT(scope.NotEmpty());
    _ASSERT(search_db.NotEmpty());
    try { 
        /// Try to open the BLAST database even for remote searches, as if
        /// it is available locally, it will be better to fetch the
        /// sequence data for formatting from this (local) source
        CRef<CSeqDB> seqdb = search_db->GetSeqDb();
        db_adapter.Reset(new CLocalDbAdapter(*search_db));
        scope->AddDataLoader(RegisterOMDataLoader(seqdb));
    } catch (const CSeqDBException&) {
            /// The BLAST database couldn't be found, report this for local
            /// searches.
    }
    _ASSERT(db_adapter && scope);
        
    if (opt.GetUseIndex() && !g_cmdLineArgs->ExecuteRemotely()) {
        BlastSeqSrc* seqsrc = db_adapter->MakeSeqSrc();
        CRef<CBlastOptions> my_options(&(g_opts_hndl->SetOptions()));
        CSetupFactory::InitializeMegablastDbIndex(seqsrc, my_options);
    }
    
    /*** Get the formatting options ***/
    CRef<CFormattingArgs> fmt_args(g_cmdLineArgs->GetFormattingArgs());
    CBlastFormat formatter(opt, 
                           *db_adapter,
                           fmt_args->GetFormattedOutputChoice(),
                           query_opts->GetParseDeflines(),
                           g_cmdLineArgs->GetOutputStream(),
                           fmt_args->GetNumDescriptions(),
                           fmt_args->GetNumAlignments(),
                           *scope,
                           opt.GetMatrixName(),
                           fmt_args->ShowGis(),
                           fmt_args->DisplayHtmlOutput(),
                           opt.GetQueryGeneticCode(),
                           opt.GetDbGeneticCode(),
                           opt.GetSumStatisticsMode(),
                           g_cmdLineArgs->ExecuteRemotely(),
                           db_adapter->GetFilteringAlgorithm(),
                           fmt_args->GetCustomOutputFormatSpec(),
                           g_cmdLineArgs->GetTask() == "megablast",
                           opt.GetMBIndexLoaded());
            
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
    
    /*** Process the input ***/
    for (; !input.End(); formatter.ResetScopeHistory()) {
        CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
        CRef<IQueryFactory> queries(new CObjMgr_QueryFactory(*query_batch));
        CRef<CSearchResultSet> results;

        CLocalBlast lcl_blast(queries, g_opts_hndl, db_adapter);
        lcl_blast.SetNumberOfThreads(g_cmdLineArgs->GetNumThreads());
        results = lcl_blast.Run();
        
        /// Blastn original printing
        //ITERATE(CSearchResultSet, result, *results) {
            //formatter.PrintOneResultSet(**result, query_batch);
        //}
        
        for (size_t i = 0; i < results->GetNumResults(); ++i) {
            if ((*results)[i].HasAlignments()) {
                CConstRef<CSeq_align_set> aln_set = (*results)[i].GetSeqAlign();
                ITERATE(CSeq_align_set::Tdata, itr_res, aln_set->Get()) {  
                    const CSeq_align& s = **itr_res;                    
                    ///
                    /// Ref: http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/doxyhtml/classCSeq__align.html
                    ///
                    /// eScore_BitScore: BLAST-specific bit score
                    /// eScore_EValue: BLAST-specific e-value
                    /// eScore_AlignLength: not a score per se, but a useful metric nonetheless. This is the sum of all aligned segments and all gaps; this excludes introns and discontinuities
                    /// eScore_IdentityCount
                    /// eScore_MismatchCount
                    /// eScore_PercentIdentity_Gapped
                    /// eScore_PercentCoverage
                    /// eScore_PercentIdentity = eScore_PercentIdentity_Gapped 
                    ///                
                    string queryID      = s.GetSeq_id(QUERY).GetSeqIdString();
                    string subID        = s.GetSeq_id(SUBJECT).GetSeqIdString();
                    uint32_t qStart     = s.GetSeqStart(QUERY);
                    uint32_t qEnd       = s.GetSeqStop(QUERY);
                    uint32_t sStart     = s.GetSeqStart(SUBJECT);
                    uint32_t sEnd       = s.GetSeqStop(SUBJECT);
                    qEnd += 1;   
                    sEnd += 1;
                    
                    //uint32_t alignLenGap   = s.GetAlignLength();
                    //uint32_t alignLenUngap = s.GetAlignLength(false);                    
                    
                    ///
                    /// Retrieve def line 
                    ///
                    const CBioseq_Handle& query_bh2 = 
                        scope->GetBioseqHandle(s.GetSeq_id(QUERY));
                    string defLine = sequence::GetTitle(query_bh2);
                    
                    /// Get GI    
                    vector<string> vecTokens;
                    boost::split(vecTokens, defLine, boost::is_any_of("|"));
                    uint32_t gi = boost::lexical_cast<uint32_t>(vecTokens[1]);
                    
                    ///
                    /// Get cutting coords    
                    /// Format: orig_header_chunkID_x_y_v_w
                    ///
                    ///          lower    upper   lower case
                    ///         |xxxxxxx|XXXXXXX|xxxxxxx|
                    ///         x       v       w       y
                    ///
                    ///         x: cutStart,   y: cutEnd
                    ///         v: upperStart, w: upperEnd       
                    /// 
                    ///  chunkID Type-0:     XXX
                    ///          Type-2,4,6: xxxXXX
                    ///          Type-3,5,7: xxxXXXxxx
                    ///          Type-1:     XXXxxx
                    ///
                    string coord = vecTokens[vecTokens.size()-1];
                    vecTokens.clear();
                    boost::split(vecTokens, coord, boost::is_any_of("_"));
                    if (vecTokens.size() != 7) {
                        MPI_Finalize();
                        cerr << "ERROR: the number of tokens in coords string != 7\n";
                        cerr << "read line = " << defLine << endl;
                        cerr << "coord = " << coord << endl;
                        cerr.flush();                        
                        exit(5);
                    }                        
                    //assert(vecTokens.size() == 7);
                    //uint32_t cId 
                        //= boost::lexical_cast<uint32_t>(vecTokens[vecTokens.size()-6]);
                    //uint32_t cType           
                        //= boost::lexical_cast<uint32_t>(vecTokens[vecTokens.size()-5]);
                    uint32_t cutStart        
                        = boost::lexical_cast<uint32_t>(vecTokens[vecTokens.size()-4]);
                    //uint32_t cutEnd          
                        //= boost::lexical_cast<uint32_t>(vecTokens[vecTokens.size()-3]);
                    uint32_t upperStart      
                        = boost::lexical_cast<uint32_t>(vecTokens[vecTokens.size()-2]);
                    uint32_t upperEnd        
                        = boost::lexical_cast<uint32_t>(vecTokens[vecTokens.size()-1]);
                    
                    /// Just in case, remove out the hits in the lowercase parts
                    if ((qStart + cutStart < upperStart && qEnd + cutStart < upperStart) ||
                        (qStart + cutStart > upperEnd && qEnd + cutStart > upperEnd)) {
                        //cout << "Warning: A HSP is found in the flank areas!" << endl; 
                        ////cout << "cid, gi, sid, qstart, qend, sstart, send, cutstart, cutend, upperstart, upperend\n";
                        //cout << cId << ","
                             //<< gi << ","
                             //<< subID << ","
                             //<< qStart << ","
                             //<< qEnd << ","
                             //<< sStart << ","
                             //<< sEnd << ","
                             //<< cutStart << ","
                             //<< cutEnd << ","
                             //<< upperStart << ","
                             //<< upperEnd << endl;
                        //cout.flush();
                    }    
                    else {
                        ///
                        /// score, bit_score, evalue...
                        ///
                        int num_ident = -1;
                        int score = 0, sum_n = 0;
                        double bit_score = 0.0, evalue = 0.0;
                        list<int> use_this_gi;
                        CAlignFormatUtil::GetAlnScores(s, score, bit_score, evalue, 
                                                       sum_n, num_ident, use_this_gi);
                                                                
                        ///
                        /// Convert Std-seg and Dense-diag alignments to Dense-seg.
                        /// Ref: http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/doxyhtml/tabular_8cpp-source.html
                        ///      http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/doxyhtml/classCSeq__id.html
                        ///
                        bool query_is_na = false, subject_is_na = false;
                        const CBioseq_Handle& query_bh = 
                            scope->GetBioseqHandle(s.GetSeq_id(QUERY));
                        query_is_na = query_bh.IsNa();
                        const CBioseq_Handle& subject_bh = 
                            scope->GetBioseqHandle(s.GetSeq_id(SUBJECT));
                        subject_is_na = subject_bh.IsNa();
                                        
                        const bool kTranslated = s.GetSegs().IsStd();
                        CRef<CSeq_align> finalAln(0);                    
                        if (kTranslated) {
                            CRef<CSeq_align> densegAln = s.CreateDensegFromStdseg();
                            if (query_is_na && subject_is_na)
                                finalAln = densegAln->CreateTranslatedDensegFromNADenseg();
                            else
                                finalAln = densegAln;
                        }
                        else if (s.GetSegs().IsDendiag()) {
                            finalAln = CAlignFormatUtil::CreateDensegFromDendiag(s);
                        }   
                        const CDense_seg& ds = (finalAln ? finalAln->GetSegs().GetDenseg() :
                                                s.GetSegs().GetDenseg());
                                                
                        CRef<CAlnVec> alnVec;
                        if (!kTranslated && ds.IsSetStrands() && 
                            ds.GetStrands().front() == eNa_strand_minus) {
                            CRef<CDense_seg> reversed_ds(new CDense_seg);
                            reversed_ds->Assign(ds);
                            reversed_ds->Reverse();
                            alnVec.Reset(new CAlnVec(*reversed_ds, *scope));   
                        } else {
                            alnVec.Reset(new CAlnVec(ds, *scope));
                        }    
                        
                        ///
                        /// Ref: http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/doxyhtml/score__builder_8cpp-source.html
                        /// pct_identity = 100.0f * float(*identities) / count_aligned;
                        ///
                        //int align_length = 0, num_gaps = 0, num_gap_opens = 0;
                        //CAlignFormatUtil::GetAlignLengths(*alnVec, align_length, 
                                                        //num_gaps, num_gap_opens);
                        //float orig_perc_ident = (align_length > 0 ? 
                                       //((float)num_ident)/align_length * 100 : 0)
                        //int num_mismatches = align_length - num_ident - num_gaps;
                        
                        ///
                        /// Doug's filtering
                        ///
                        /// "Percent Identity” looks at every position in the aligned sequences 
                        /// and counts the number that have the same base or amino acid.  
                        /// This count is then divided by the aligned length and then multiplied by 100.
                        ///
                        /// identity = # of identical bases / length of read
                        /// coverage = (read end – read begin)/ length of read 
                        ///
                        /// identity > 50%
                        /// coverage > 90%
                        ///
                        /// Ref: http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/doxyhtml/score__builder_8cpp-source.html
                        /// pct_coverage = 100.0f * float(covered_bases) / float(seq_len);
                        ///
                        string querySeq = "";
                        string subjectSeq = "";
                        alnVec->SetGapChar('-');
                        alnVec->GetWholeAlnSeqString(QUERY, querySeq);
                        alnVec->GetWholeAlnSeqString(SUBJECT, subjectSeq);
                        assert(querySeq.size() > 0 && subjectSeq.size() > 0);
                        
                        /// Convert to the original coord in the query sequence
                        uint32_t qStartOrig = qStart + cutStart;
                        uint32_t qEndOrig = qEnd + cutStart;               
                        uint32_t scanStart = 0;
                                           
                        ///
                        ///     q.start                        q.end
                        ///        |xxxxxxxxxxxxxxxxxxxxxxxxxxxxx|
                        ///
                        ///    |xxxxxxxxxxx|XXXXXXXXXXXXXXX|xxxxxxxxxxx|
                        /// cutStart   upperStart      upperEnd      cutEnd
                        ///
                        ///         ------>|               |<-----
                        ///             scanStart       scanEnd
                        ///
                        
                        /// Why take min()? Ref: http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/doxyhtml/tabular_8cpp-source.html
                        uint32_t scanEnd = min(querySeq.size(), subjectSeq.size());
                        if (qStartOrig < upperStart) 
                            scanStart += (upperStart - qStartOrig);
                        if (qEndOrig > upperEnd) 
                            scanEnd -= (qEndOrig - upperEnd);
        
                        ///
                        /// Compare aligned query and subject sequences
                        /// Should consider the gaps in the align when counting 
                        /// num_ident_upperpart and the upperStart and upperEnd 
                        /// which represents the uppercase bases in the query
                        ///
                        int scanStartInAlign = 0;
                        int charCount = 0;
                        if (scanStart > 0) {
                            while (charCount < (upperStart - qStartOrig)) {
                                if (querySeq[scanStartInAlign] != '-') charCount++;
                                scanStartInAlign++;
                            }
                        }
                        charCount = 0;
                        int scanEndInAlign = min(querySeq.size(), subjectSeq.size());
                        if (scanEnd < min(querySeq.size(), subjectSeq.size())) {
                            while (charCount < (qEndOrig - upperEnd)) {
                                if (querySeq[scanEndInAlign-1] != '-') charCount++;
                                scanEndInAlign--;
                            }
                        }
                        int num_ident_upperpart = 0;                    
                        for (uint32_t i = scanStartInAlign; i < scanEndInAlign; ++i) {
                            if (querySeq[i] == subjectSeq[i]) ++num_ident_upperpart;
                        }         
                        
                        float doug_perc_ident = float(num_ident_upperpart) / 
                                                (upperEnd - upperStart) * 100;
                               
                        ///                 
                        /// Adjust qstart and qend for computing coverage
                        ///
                        int newQStart = qStartOrig, newQEnd = qEndOrig;
                        if (qStartOrig < upperStart) newQStart = upperStart;
                        if (qEndOrig > upperEnd) newQEnd = upperEnd;
                        float doug_perc_cover = float(newQEnd - newQStart) / 
                                                (upperEnd - upperStart) * 100;                   
                        
                        ///
                        /// Add a HSP to kv and emit
                        ///
                        if (doug_perc_ident >= g_IDENT_CUTOFF && 
                            doug_perc_cover >= g_COVER_CUTOFF) {
                            structBlRes_t res;
                            res.subjectId     = boost::lexical_cast<uint32_t>(subID);
                            //res.identity      = orig_perc_ident;
                            //res.alignLen      = align_length;
                            //res.misMatches    = num_mismatches;
                            //res.gapOpens      = num_gap_opens;
                            res.qStart        = qStart;
                            res.qEnd          = qEnd;
                            res.sStart        = sStart;
                            res.sEnd          = sEnd;
                            res.eValue        = evalue;
                            res.bitScore      = boost::lexical_cast<float>(bit_score);
                            res.upperStart    = upperStart;
                            res.upperEnd      = upperEnd;
                            res.doug_identity = doug_perc_ident;
                            res.doug_coverage = doug_perc_cover;
                            
                            uint32_t newKey   = gi;
                            kv->add((char*)&newKey, sizeof(uint32_t), (char*)&res,
                                    sizeof(structBlRes_t));     
                        }
                    } /// else
                }
            }            
        }    
    }
    
    ///
    /// BLAST call end time
    ///
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
}

 

/** Sort function - Passed to MR-MPI sort_values() for sorting blast result
 * string by bit score.
 * @param e1
 * @param e2
 */

inline bool compare_evalue(structEValue_t e1,
                           structEValue_t e2)
{
    return (e1.subjectId != e2.subjectId) ? (e1.subjectId < e2.subjectId) :
                ((e1.eValue != e2.eValue) ? (e1.eValue < e2.eValue) : 
                (e1.bitScore < e2.bitScore));
}


/** Sort by evalue - Passed to MR-MPI reduce() for sorting KMVs by evalue.
 * @param key
 * @param keybytes
 * @param multivalue - collected blast result strings.  
 * @param nvalues
 * @param valuebytes
 * @param kv
 * @param ptr
 */

inline void mr_sort_multivalues_by_evalue(char *key,
                                   int keybytes,
                                   char *multivalue,
                                   int nvalues,
                                   int *valuebytes,
                                   KeyValue *kv,
                                   void *ptr)
{
    ///
    /// Make structEValue_t = {structBlRes_t* pRec; float evalue;}
    /// and sort by evalue
    ///
    vector<structEValue_t> vecHit;
    for (size_t n = 0; n < (unsigned)nvalues; n++) {
        structBlRes_t* res = (structBlRes_t*)multivalue;
        structEValue_t strctEvalue;
        strctEvalue.pRec = res;
        strctEvalue.subjectId = res->subjectId;
        strctEvalue.eValue = res->eValue;
        strctEvalue.bitScore = res->bitScore;
        vecHit.push_back(strctEvalue);
        multivalue += sizeof(structBlRes_t);
    }
    
    ///////////////////////////////////////////////////
    sort(vecHit.begin(), vecHit.end(), compare_evalue);
    ///////////////////////////////////////////////////
    
    ///
    /// After sorting, each worker saves the set of results into a file.
    /// Note: The file open option is "a".
    /// Note: outfmt = 6 in Blast
    ///                       
    //ofstream outputFile(g_hitFileName.c_str(), ios::out | ios::app);
    //for (size_t n = 0; n < (unsigned)nvalues; n++) {
        //structBlRes_t* res = (structBlRes_t*)(vecHit[n].pRec);
        //outputFile 
            ////<< res->cId << "\t"
            //<< *(uint32_t*)key << "\t"
            //<< res->subjectId << "\t"
            ////<< res->identity << "\t"
            ////<< res->alignLen << "\t"
            ////<< res->misMatches << "\t"
            ////<< res->gapOpens << "\t"
            //<< res->qStart << "\t"
            //<< res->qEnd << "\t"
            //<< res->sStart << "\t"
            //<< res->sEnd << "\t"
            //<< res->eValue << "\t"
            //<< res->bitScore << "\t"
            //<< res->upperStart << "\t"
            //<< res->upperEnd << "\t"
            //<< res->doug_identity << "\t"
            //<< res->doug_coverage
            //<< endl;
    //}
    //outputFile.close();
    
    ///
    /// Write to file in binary format
    ///
    ofstream outputBinFile((g_hitFileName+".bin").c_str(), ios::binary | ios::app);
    for (size_t n = 0; n < (unsigned)nvalues; n++) {
        structBlRes_t* res = (structBlRes_t*)(vecHit[n].pRec);
        structBlRes2_t hit;
        hit.qId            = *(uint32_t*)key;
        hit.subjectId      = res->subjectId;
        hit.qStart         = res->qStart;
        hit.qEnd           = res->qEnd;
        hit.sStart         = res->sStart;
        hit.sEnd           = res->sEnd;
        hit.eValue         = res->eValue;
        hit.bitScore       = res->bitScore;
        hit.upperStart     = res->upperStart;
        hit.upperEnd       = res->upperEnd;
        hit.doug_identity  = res->doug_identity;
        hit.doug_coverage  = res->doug_coverage;
        
        outputBinFile.write((char*)&hit, sizeof(hit));
    }
    outputBinFile.close();
    vecHit.clear();
}


/** Collect MPI node names and ranks
 * @param rank
 * @param numProcs
 * @param mpiComm
 */
 
inline void mpi_collect_node_name(int rank, int numProcs, MPI_Comm mpiComm)
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

