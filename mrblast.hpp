//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the MGTAXA package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

/**    
    @file   mrblast.hpp
    @date   02/22/2012
    @author Seung-Jin Sul (ssul@jcvi.org)
    @brief  MR-MPI-BLAST: Parallelizing BLAST using MapReduce-MPI header file
      
            Version: 13.0.0
            This version of mr-mpi-blast supports both generic and classifier\n
            run modes. In the classifier run mode (ISCLASSIFIER=1 in
            "mrblast.ini"), two additional fields are added for classifier:
            percentage identity and percentage coverage for query. 
*/

#ifndef MRBLAST_HPP
#define MRBLAST_HPP

/// For typedef unsigned long long int uint32_t
#include <stdint.h>

/// For setfill and setw
#include <iomanip>
#include <sstream>

/// MPI AND MapReduce-MPI
#include "mpi.h"
#include "mrmpi/mapreduce.h"
#include "mrmpi/keyvalue.h"

/// For syncronized timing
#ifndef MPI_WTIME_IS_GLOBAL
#define MPI_WTIME_IS_GLOBAL 1
#endif


    
using namespace MAPREDUCE_NS;
using namespace std;

/// ----------------------------------------------------------------------------
/// NCBI C++ Toolkit
/// ----------------------------------------------------------------------------
#include <algo/blast/api/local_blast.hpp>
#include <algo/blast/api/objmgr_query_data.hpp>
#include <algo/blast/blastinput/blast_fasta_input.hpp>
#include <algo/blast/api/blast_nucl_options.hpp>
#include <algo/blast/api/blast_prot_options.hpp>

/// For importing search strategy
//#include <objects/blast/Blast4_request.hpp>
//#include <algo/blast/api/search_strategy.hpp>

/// For accessing fields in tabular output
//#include <algo/blast/format/blast_format.hpp>
#include <algo/blast/blastinput/blastn_args.hpp>
#include <algo/blast/blastinput/blastp_args.hpp>

/// sequence::GetTitle() for get def line
#include <objmgr/util/sequence.hpp>

/// For getting align info and RecoverSearchStrategy
#include "blast_app_util.hpp"

/// For blastdbcmd, setting the dbsize
#include <objtools/blast/seqdb_reader/seqdbexpert.hpp>

/// For EScoreType
#include <objects/seqalign/Seq_align.hpp>

/// For sequence::GetLength()
#include <objmgr/util/seq_loc_util.hpp>

/// Tabular
#include <algo/blast/format/blast_format.hpp>
#include <algo/blast/blastinput/blastn_args.hpp>
#include <algo/blast/blastinput/blastp_args.hpp>

USING_NCBI_SCOPE;
USING_SCOPE (blast);
USING_SCOPE (align_format);

const int QUERY = 0;        /**< To retireve query info from CSeq_align */
const int SUBJECT = 1;      /**< To retireve suject info from CSeq_align */
 
/// For NcbiApp
CRef<CBlastOptionsHandle> g_optsHndl;
CRef<CBlastAppArgs> g_cmdLineArgs;

/// Blast target DB setting
CRef<CSearchDatabase> g_searchDatabase;     /**< ref to database partition */
string g_prevDbName = "";                   /**< db partition name previously used */
    
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
/// for real file size
#include <boost/filesystem/operations.hpp> 
boost::iostreams::mapped_file_source g_memmapQueryFile; /**< Read-only Boost mmap query file */
uint64_t g_realFileSize = 0;

/// String trim
#include <boost/algorithm/string/trim.hpp>

/// ----------------------------------------------------------------------------
/// Settings from mrblast.ini file
/// ----------------------------------------------------------------------------
/// MapReduce-MPI options
int g_verbosity;
int g_timer;                /// Log elapsed time for each mapreduce call
int g_memSize;              /// Page size (in Mbytes)
int g_outOfCore;

/// DB options
int g_nDbFiles;                         /**< total number of db partitions */
string g_dbName;
string g_dbListFileName;                /**< db partition list file name */ 
string g_programName;
const string CONF_FILE_NAME = "mrblast.ini"; /**< configuration file name */
const int MAXPROCNAME = 80;             /**< For mpi proc name */
const int MAXSUBJID = 40;               /**< For subject id string */
const int MAXDBLENSTR = 40;             /**< For database length string */

/// ----------------------------------------------------------------------------
/// Log
/// ----------------------------------------------------------------------------
#include <sys/time.h>
#include <sys/resource.h>
//#define NDEBUG 1
int g_bLogEnabled = 0;
int g_timingEnabled = 0;
int g_optDumpEnabled = 0;               /**< For dumping Blast opts_hndl out */
int g_mapCallNo = 0;                    /**< Serial no of map calls in each rank */
string g_logFileName;
string g_logMsg;
ofstream g_logFileStream;
static ostream *g_pLogStream = NULL;
#define LOG (*g_pLogStream)

/// ----------------------------------------------------------------------------
/// custom scheduler
/// ----------------------------------------------------------------------------
typedef multimap<string, int> multimapSI_t;
typedef multimap<string, vector<int> > multimapSVI_t;
typedef pair<string, vector<int> > pairSVI_t;
typedef multimap<int, int> multimapII_t;
typedef map<int, string> mapIS_t;
typedef map<string, string> mapSS_t;

/**
 * Struct for loading query index file for constructing work items
 */
typedef struct structQueryIndex {
    uint64_t qStart;                    /**< query start loc in orig query file */
    uint32_t qLength;                   /**< quety len in bp */
    uint64_t uniqQueryId;               /**< unique query id (gi or serial num) */
} structQueryIndex_t;

/**
 * Struct for query block
 */
typedef struct structBlockBeginLoc {
    uint64_t blockBeginLoc;             /**< query start loc of a work item block */
    uint64_t qIdStart;                  /**< starting query id */
} structBlockBeginLoc_t;

/**
 * Struct for work item
 */
typedef struct structWorkItem {
    int      dbNo;                      /**< db ID */
    uint64_t blockBegin;                /**< query start offset of a work item block */
    uint64_t blockEnd;                  /**< query end offset of a work item block */
    uint64_t qIdStart;                  /**< starting query ID */
} structWorkItem_t;

vector<string> g_vecDbFile;
vector<structBlockBeginLoc_t> g_vecBlockBeginLoc; /// for generic version
vector<structQueryIndex_t> g_vecQueryIndex;
vector<structWorkItem_t> g_vecWorkItem;
uint32_t g_nQueries;                    /**< total number of queries */
uint32_t g_nQueryBlocks;                /**< total number of query blocks */
uint32_t g_nWorkItems;                  /**< total number of work items */
uint32_t g_blockSize;                   /**< block size in base-pair */
multimapSI_t g_multimapProcNameRank;    /**< dict of rank by proc name */
mapIS_t g_mapRankProcName;              /**< dictionary of proc name by rank */

/// ----------------------------------------------------------------------------
/// Misc.
/// ----------------------------------------------------------------------------
int g_mapStyle;                         /**< MapReduce-MPI scheduler */
int g_MPI_worldRank;                    /**< MPI rank */
int g_MPI_nProcs;                       /**< MPI number of procs */
char g_MPI_procName[MAXPROCNAME];       /**< MPI procname */
string g_outFilePrefix;                 /**< Prefix string for output file names */
string g_indexFileName;                 /**< index file name */
string g_queryFileName;                 /**< query file name */
bool g_bIsQueryProtein = false;         
bool g_bIsDbProtein = false;
bool g_bClassifier = false;             /**< Add classifier fields in final BLAST output */
string g_hitFileName;                   /**< Output file to store BLAST hits */
 
/**
 * Struct for emitting MapReduce KV for generic run mode (w/o percentage
 * identity and percentage coverage for query output)
 */
typedef struct structBlResGeneric {
    char        subjectId[MAXSUBJID];   /**< Subject Seq-id */
    double      identity;               /**< BLAST percentage identity (Percentage of identical matches) (%) */
    uint32_t    alignLen;               /**< Alignment length */
    uint32_t    nMismatches;            /**< Number of mismatches */
    uint32_t    nGaps;                  /**< Total number of gaps */
    uint32_t    qStart;                 /**< Start of alignment in query */
    uint32_t    qEnd;                   /**< End of alignment in query */
    uint32_t    sStart;                 /**< Start of alignment in subject */
    uint32_t    sEnd;                   /**< End of alignment in subject */
    double      eValue;                 /**< Expect value */
    double      bitScore;               /**< Bit score */
} structBlResGeneric_t;

/**
 * Struct for emitting MapReduce KV for classifier run mode (with percentage
 * identity and percentage coverage for query output)
 */
typedef struct structBlResClassifier {
    char        subjectId[MAXSUBJID];   /**< Subject Seq-id */
    double      identity;               /**< BLAST percentage identity (Percentage of identical matches) (%) */
    uint32_t    alignLen;               /**< Alignment length */
    uint32_t    nMismatches;            /**< Number of mismatches */
    uint32_t    nGaps;                  /**< Total number of gaps */
    uint32_t    qStart;                 /**< Start of alignment in query */
    uint32_t    qEnd;                   /**< End of alignment in query */
    uint32_t    sStart;                 /**< Start of alignment in subject */
    uint32_t    sEnd;                   /**< End of alignment in subject */
    double      eValue;                 /**< Expect value */
    double      bitScore;               /**< Bit score */
    double      percIdent;              /**< Percentage identity for query (identity count/query length * 100, computed if classifier option is enabled) */
    double      percCover;              /**< Percentage coverage for subject ((qend-qstart)/query length * 100, computed if classifier option is enabled) */
} structBlResClassifier_t;  

/** To sort Blast hits by evalue & bitscore in generic run mode
 */
typedef struct structEvalue {
    structBlResGeneric_t *pRec;         /**< pointer to BLAST hit result set */
    double                eValue;       /**< evalue to sort */
    double                bitScore;     /**< bit score to sort */
    double                identity;     /**< identity to sort */
} structEValue_t;

/** To sort Blast hits by evalue & bitscore in classifier run mode
 */
typedef struct structEvalueClassifier {
    structBlResClassifier_t *pRec;      /**< pointer to BLAST hit result set */
    double                   eValue;    /**< evalue to sort */
    double                   bitScore;  /**< bit score to sort */
    double                   identity;  /**< identity to sort */
} structEValueClassifier_t;

/** For saving hits in bin format in generic run mode
 */
typedef struct structBlResToSaveHits {
    uint64_t    queryId;                /**< Unique query ID */
    char        subjectId[MAXSUBJID];   /**< Subject Seq-id */
    double      identity;               /**< BLAST percentage identity (Percentage of identical matches) (%) */
    uint32_t    alignLen;               /**< Alignment length */
    uint32_t    nMismatches;            /**< Number of mismatches */
    uint32_t    nGaps;                  /**< Total number of gaps */
    uint32_t    qStart;                 /**< Start of alignment in query */
    uint32_t    qEnd;                   /**< End of alignment in query */
    uint32_t    sStart;                 /**< Start of alignment in subject */
    uint32_t    sEnd;                   /**< End of alignment in subject */
    double      eValue;                 /**< Expect value */
    double      bitScore;               /**< Bit score */
} structBlResToSaveHits_t;              

/** For saving hits in bin format in classifier run mode
 */
typedef struct structBlResToSaveHitsClassifier {
    uint64_t    queryId;                /**< Unique query ID */
    char        subjectId[MAXSUBJID];   /**< Subject Seq-id */
    double      identity;               /**< BLAST percentage identity (Percentage of identical matches) (%) */
    uint32_t    alignLen;               /**< Alignment length */
    uint32_t    nMismatches;            /**< Number of mismatches */
    uint32_t    nGaps;                  /**< Total number of gaps */
    uint32_t    qStart;                 /**< Start of alignment in query */
    uint32_t    qEnd;                   /**< End of alignment in query */
    uint32_t    sStart;                 /**< Start of alignment in subject */
    uint32_t    sEnd;                   /**< End of alignment in subject */
    double      eValue;                 /**< Expect value */
    double      bitScore;               /**< Bit score */
    double      percIdent;              /**< Percentage identity for query (identity count/query length * 100, computed if classifier option is enabled) */
    double      percCover;              /**< Percentage coverage for subject ((qend-qstart)/query length * 100, computed if classifier option is enabled) */
} structBlResToSaveHitsClassifier_t;

/// For multiple iterations
int g_nIter = 1;                        /**< number of iterations */
int g_currIter = 0;
uint32_t nSubWorkItemSets = 0;          /**< number of sub sets of work items */
uint32_t nWorkItemsPerIter = 0;         /**< number of work items per iteration */
uint32_t nRemains = 0;

/** For sorting hits by qid -
 * MPI_Allreduce is used for collecting the number of hits on each rank and
 * summing up, and broadcasting summed up array to all ranks.
 */
uint32_t *g_vecNumHitsPerQid;           /**< arrary for hit count per query */
uint32_t *g_vecNumHitsPerQid2;          /**< arrary for hit count per query */

/// prototypes
void mr_mpi_blast(); 
void mr_map_run_blast(int itask, KeyValue *kv, void *ptr);
inline void mr_reduce_sort_and_save_generic(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
inline void mr_reduce_sort_and_save_classifier(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
inline int  mr_myhash(char* key, int len);
inline void mpi_collect_node_name(int rank, int nProcs, MPI_Comm mpiComm);
inline bool compare_evalue_generic(structEValue_t e1, structEValue_t e2);
inline bool compare_evalue_classifier(structEValueClassifier_t e1, structEValueClassifier_t e2);

#endif
